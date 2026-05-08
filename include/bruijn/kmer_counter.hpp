#pragma once

#include "bruijn/bloom_filter.hpp"
#include "bruijn/compact_kmer_map.hpp"
#include "bruijn/dna.hpp"
#include "bruijn/read_buffer.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <thread>
#include <vector>

namespace bruijn {

struct CounterConfig {
    uint32_t min_count = 2;
    std::size_t bloom_counters = 1U << 24U;
    uint8_t bloom_hashes = 4;
    std::size_t thread_count = std::thread::hardware_concurrency();
    bool canonicalize = false;
};

template <std::size_t K>
class KmerCounter {
public:
    using CounterMap = CompactKmerMap<uint32_t>;

    explicit KmerCounter(CounterConfig config)
        : config_(config) {
        if (config_.thread_count == 0) {
            config_.thread_count = 1;
        }
        if (config_.bloom_counters == 0) {
            config_.bloom_counters = 1U << 20U;
        }
        if (config_.bloom_hashes == 0) {
            config_.bloom_hashes = 1;
        }
    }

    [[nodiscard]] CounterMap count(const ReadBuffer& reads) const {
        if (!config_.canonicalize) {
            if (config_.min_count <= 1) {
                return count_keys_without_bloom(reads, false);
            }
            return filter_by_min_count(count_keys_with_bloom(reads, false));
        }

        if (config_.min_count <= 1) {
            // Keep directional graph nodes when no abundance filter is applied.
            return count_keys_without_bloom(reads, false);
        }

        CounterMap canonical = filter_by_min_count(count_keys_with_bloom(reads, true));
        return count_raw_for_canonical_set(reads, canonical);
    }

private:
    struct Range {
        std::size_t begin = 0;
        std::size_t end = 0;
    };

    CounterConfig config_;

    [[nodiscard]] std::vector<Range> split_ranges(std::size_t total, std::size_t chunks) const {
        chunks = std::max<std::size_t>(1, std::min<std::size_t>(chunks, total == 0 ? 1 : total));
        std::vector<Range> ranges;
        ranges.reserve(chunks);

        const std::size_t base = total / chunks;
        const std::size_t extra = total % chunks;

        std::size_t begin = 0;
        for (std::size_t i = 0; i < chunks; ++i) {
            const std::size_t length = base + (i < extra ? 1U : 0U);
            ranges.push_back(Range{begin, begin + length});
            begin += length;
        }

        return ranges;
    }

    static void merge_counters(CounterMap& dst, const CounterMap& src) {
        src.for_each([&dst](uint64_t key, uint32_t value) {
            dst.increment(key, value);
        });
    }

    [[nodiscard]] static uint64_t normalize_key(uint64_t key, bool canonicalize_keys) {
        if (!canonicalize_keys) {
            return key;
        }
        return dna::canonical_kmer<K>(key);
    }

    [[nodiscard]] CounterMap filter_by_min_count(const CounterMap& counts) const {
        if (config_.min_count <= 1) {
            return counts;
        }

        CounterMap filtered(counts.size());
        counts.for_each([&](uint64_t key, uint32_t value) {
            if (value >= config_.min_count) {
                filtered.set(key, value);
            }
        });
        return filtered;
    }

    [[nodiscard]] CounterMap count_keys_without_bloom(const ReadBuffer& reads, bool canonicalize_keys) const {
        const std::size_t read_count = reads.read_count();
        if (read_count == 0) {
            return CounterMap{};
        }

        const std::size_t threads_to_use = std::min<std::size_t>(config_.thread_count, read_count);
        const auto ranges = split_ranges(read_count, threads_to_use);

        std::vector<CounterMap> local_maps;
        local_maps.reserve(ranges.size());
        for (const Range& range : ranges) {
            local_maps.emplace_back(std::max<std::size_t>(1, (range.end - range.begin) * 8U));
        }

        std::vector<std::thread> workers;
        workers.reserve(ranges.size());

        for (std::size_t i = 0; i < ranges.size(); ++i) {
            workers.emplace_back([&, i]() {
                CounterMap& map = local_maps[i];
                for (std::size_t r = ranges[i].begin; r < ranges[i].end; ++r) {
                    const std::string_view read = reads.read(r);
                    dna::for_each_kmer<K>(read, [&](uint64_t kmer, std::size_t) {
                        map.increment(normalize_key(kmer, canonicalize_keys));
                    });
                }
            });
        }

        for (std::thread& worker : workers) {
            worker.join();
        }

        CounterMap global(read_count * 16U);
        for (const CounterMap& local : local_maps) {
            merge_counters(global, local);
        }

        return global;
    }

    [[nodiscard]] CounterMap count_keys_with_bloom(const ReadBuffer& reads, bool canonicalize_keys) const {
        const std::size_t read_count = reads.read_count();
        if (read_count == 0) {
            return CounterMap{};
        }

        const std::size_t threads_to_use = std::min<std::size_t>(config_.thread_count, read_count);
        const auto ranges = split_ranges(read_count, threads_to_use);

        std::vector<CountingBloomFilter> local_blooms;
        local_blooms.reserve(ranges.size());
        for (std::size_t i = 0; i < ranges.size(); ++i) {
            local_blooms.emplace_back(config_.bloom_counters, config_.bloom_hashes);
        }

        std::vector<std::thread> workers;
        workers.reserve(ranges.size());

        for (std::size_t i = 0; i < ranges.size(); ++i) {
            workers.emplace_back([&, i]() {
                CountingBloomFilter& bloom = local_blooms[i];
                for (std::size_t r = ranges[i].begin; r < ranges[i].end; ++r) {
                    const std::string_view read = reads.read(r);
                    dna::for_each_kmer<K>(read, [&](uint64_t kmer, std::size_t) {
                        bloom.increment(normalize_key(kmer, canonicalize_keys));
                    });
                }
            });
        }

        for (std::thread& worker : workers) {
            worker.join();
        }

        CountingBloomFilter merged_bloom(config_.bloom_counters, config_.bloom_hashes);
        for (const CountingBloomFilter& local : local_blooms) {
            merged_bloom.merge_from(local);
        }

        std::vector<CounterMap> local_maps;
        local_maps.reserve(ranges.size());
        for (const Range& range : ranges) {
            local_maps.emplace_back(std::max<std::size_t>(1, (range.end - range.begin) * 4U));
        }

        workers.clear();
        for (std::size_t i = 0; i < ranges.size(); ++i) {
            workers.emplace_back([&, i]() {
                CounterMap& map = local_maps[i];
                for (std::size_t r = ranges[i].begin; r < ranges[i].end; ++r) {
                    const std::string_view read = reads.read(r);
                    dna::for_each_kmer<K>(read, [&](uint64_t kmer, std::size_t) {
                        const uint64_t key = normalize_key(kmer, canonicalize_keys);
                        if (merged_bloom.estimate(key) >= config_.min_count) {
                            map.increment(key);
                        }
                    });
                }
            });
        }

        for (std::thread& worker : workers) {
            worker.join();
        }

        CounterMap global(read_count * 8U);
        for (const CounterMap& local : local_maps) {
            merge_counters(global, local);
        }

        return global;
    }

    [[nodiscard]] CounterMap count_raw_for_canonical_set(const ReadBuffer& reads, const CounterMap& canonical_set) const {
        if (canonical_set.empty()) {
            return CounterMap{};
        }

        const std::size_t read_count = reads.read_count();
        if (read_count == 0) {
            return CounterMap{};
        }

        const std::size_t threads_to_use = std::min<std::size_t>(config_.thread_count, read_count);
        const auto ranges = split_ranges(read_count, threads_to_use);

        std::vector<CounterMap> local_maps;
        local_maps.reserve(ranges.size());
        for (const Range& range : ranges) {
            local_maps.emplace_back(std::max<std::size_t>(1, (range.end - range.begin) * 4U));
        }

        std::vector<std::thread> workers;
        workers.reserve(ranges.size());

        for (std::size_t i = 0; i < ranges.size(); ++i) {
            workers.emplace_back([&, i]() {
                CounterMap& map = local_maps[i];
                for (std::size_t r = ranges[i].begin; r < ranges[i].end; ++r) {
                    const std::string_view read = reads.read(r);
                    dna::for_each_kmer<K>(read, [&](uint64_t kmer, std::size_t) {
                        const uint64_t canonical = dna::canonical_kmer<K>(kmer);
                        if (canonical_set.find(canonical) != nullptr) {
                            map.increment(kmer);
                        }
                    });
                }
            });
        }

        for (std::thread& worker : workers) {
            worker.join();
        }

        CounterMap raw_counts(read_count * 8U);
        for (const CounterMap& local : local_maps) {
            merge_counters(raw_counts, local);
        }
        return raw_counts;
    }
};

}  // namespace bruijn
