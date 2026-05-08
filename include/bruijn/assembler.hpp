#pragma once

#include "bruijn/debruijn_graph.hpp"
#include "bruijn/kmer_counter.hpp"
#include "bruijn/read_buffer.hpp"

#include <cstddef>
#include <cstdint>
#include <thread>
#include <vector>

namespace bruijn {

struct AssemblerConfig {
    uint32_t min_count = 2;
    std::size_t thread_count = 0;
    std::size_t bloom_counters = 1U << 24U;
    uint8_t bloom_hashes = 4;
    bool canonical_kmers = false;
    std::size_t max_tip_length = 8;
    std::size_t bubble_depth = 64;
};

template <std::size_t K>
class DeBruijnAssembler {
public:
    struct Result {
        std::vector<typename DeBruijnGraph<K>::Contig> contigs;
        std::size_t retained_kmers = 0;
        std::size_t node_count = 0;
        std::size_t edge_count = 0;
    };

    explicit DeBruijnAssembler(AssemblerConfig config)
        : config_(config) {}

    [[nodiscard]] Result assemble(const ReadBuffer& reads) const {
        CounterConfig counter_config;
        counter_config.min_count = config_.min_count;
        counter_config.bloom_counters = config_.bloom_counters;
        counter_config.bloom_hashes = config_.bloom_hashes;
        counter_config.thread_count = (config_.thread_count == 0) ? std::thread::hardware_concurrency() : config_.thread_count;
        counter_config.canonicalize = config_.canonical_kmers;

        KmerCounter<K> counter(counter_config);
        auto counts = counter.count(reads);

        DeBruijnGraph<K> graph;
        graph.build(counts, reads);
        graph.remove_tips(config_.max_tip_length);
        graph.remove_bubbles(config_.bubble_depth);

        Result result;
        result.retained_kmers = counts.size();
        result.node_count = graph.node_count();
        result.edge_count = graph.edge_count();
        result.contigs = graph.generate_contigs();
        return result;
    }

private:
    AssemblerConfig config_;
};

}  // namespace bruijn
