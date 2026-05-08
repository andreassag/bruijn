#pragma once

#include "bruijn/dna.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <vector>

namespace bruijn {

class CountingBloomFilter {
public:
    CountingBloomFilter() = default;

    CountingBloomFilter(std::size_t counter_count, uint8_t hash_functions)
        : counter_count_(counter_count),
          hash_functions_(std::max<uint8_t>(hash_functions, 1)),
          counters_(counter_count, 0) {
        if (counter_count_ == 0) {
            throw std::invalid_argument("Bloom filter counter count must be > 0");
        }
    }

    [[nodiscard]] std::size_t counter_count() const noexcept { return counter_count_; }
    [[nodiscard]] uint8_t hash_functions() const noexcept { return hash_functions_; }
    [[nodiscard]] std::size_t memory_bytes() const noexcept { return counters_.size(); }

    void increment(uint64_t key) {
        for (uint8_t i = 0; i < hash_functions_; ++i) {
            const std::size_t idx = nth_index(key, i);
            uint8_t& value = counters_[idx];
            if (value != 0xFFU) {
                ++value;
            }
        }
    }

    [[nodiscard]] uint8_t estimate(uint64_t key) const {
        uint8_t min_value = 0xFFU;
        for (uint8_t i = 0; i < hash_functions_; ++i) {
            const std::size_t idx = nth_index(key, i);
            min_value = std::min(min_value, counters_[idx]);
        }
        return min_value;
    }

    void merge_from(const CountingBloomFilter& other) {
        if (counter_count_ != other.counter_count_ || hash_functions_ != other.hash_functions_) {
            throw std::invalid_argument("Bloom filters are incompatible for merge");
        }

        for (std::size_t i = 0; i < counters_.size(); ++i) {
            const uint16_t merged = static_cast<uint16_t>(counters_[i]) + static_cast<uint16_t>(other.counters_[i]);
            counters_[i] = static_cast<uint8_t>(std::min<uint16_t>(merged, 0xFFU));
        }
    }

private:
    std::size_t counter_count_ = 0;
    uint8_t hash_functions_ = 1;
    std::vector<uint8_t> counters_;

    [[nodiscard]] std::size_t nth_index(uint64_t key, uint8_t n) const noexcept {
        const uint64_t mixed = dna::splitmix64(key + 0x9e3779b97f4a7c15ULL * static_cast<uint64_t>(n + 1U));
        return static_cast<std::size_t>(mixed % counter_count_);
    }
};

}  // namespace bruijn
