#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <string_view>
#include <vector>

namespace bruijn::dna {

inline uint8_t normalize_base(char c) {
    switch (c) {
        case 'A':
        case 'a':
            return static_cast<uint8_t>('A');
        case 'C':
        case 'c':
            return static_cast<uint8_t>('C');
        case 'G':
        case 'g':
            return static_cast<uint8_t>('G');
        case 'T':
        case 't':
            return static_cast<uint8_t>('T');
        default:
            return static_cast<uint8_t>('N');
    }
}

inline bool base_to_bits(char c, uint8_t& bits) {
    switch (c) {
        case 'A':
        case 'a':
            bits = 0;
            return true;
        case 'C':
        case 'c':
            bits = 1;
            return true;
        case 'G':
        case 'g':
            bits = 2;
            return true;
        case 'T':
        case 't':
            bits = 3;
            return true;
        default:
            return false;
    }
}

inline char bits_to_base(uint8_t bits) {
    static constexpr char LUT[4] = {'A', 'C', 'G', 'T'};
    return LUT[bits & 0x3U];
}

inline constexpr uint8_t complement_bits(uint8_t bits) {
    return static_cast<uint8_t>(0x3U - (bits & 0x3U));
}

inline uint64_t splitmix64(uint64_t x) {
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30U)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27U)) * 0x94d049bb133111ebULL;
    return x ^ (x >> 31U);
}

template <std::size_t K>
constexpr uint64_t kmer_mask() {
    static_assert(K > 0 && K <= 31, "K must be in [1, 31] to fit in uint64_t");
    return (1ULL << (2ULL * K)) - 1ULL;
}

template <std::size_t K>
constexpr uint64_t reverse_complement(uint64_t kmer) {
    uint64_t rc = 0;
    for (std::size_t i = 0; i < K; ++i) {
        const uint8_t bits = static_cast<uint8_t>(kmer & 0x3U);
        rc = (rc << 2U) | static_cast<uint64_t>(complement_bits(bits));
        kmer >>= 2U;
    }
    return rc;
}

template <std::size_t K>
constexpr uint64_t canonical_kmer(uint64_t kmer) {
    const uint64_t rc = reverse_complement<K>(kmer);
    return std::min(kmer, rc);
}

template <std::size_t K, typename Fn>
void for_each_kmer(std::string_view sequence, Fn&& fn) {
    static_assert(K > 0 && K <= 31, "K must be in [1, 31]");

    if (sequence.size() < K) {
        return;
    }

    constexpr uint64_t MASK = kmer_mask<K>();
    uint64_t rolling = 0;
    std::size_t valid_run = 0;

    for (std::size_t i = 0; i < sequence.size(); ++i) {
        uint8_t bits = 0;
        if (!base_to_bits(sequence[i], bits)) {
            rolling = 0;
            valid_run = 0;
            continue;
        }

        rolling = ((rolling << 2U) | static_cast<uint64_t>(bits)) & MASK;
        ++valid_run;
        if (valid_run >= K) {
            const std::size_t start = i + 1U - K;
            fn(rolling, start);
        }
    }
}

template <std::size_t K>
void decode_kmer(uint64_t kmer, std::vector<char>& out) {
    out.assign(K, 'A');
    for (std::size_t i = 0; i < K; ++i) {
        const std::size_t pos = K - 1U - i;
        out[pos] = bits_to_base(static_cast<uint8_t>(kmer & 0x3U));
        kmer >>= 2U;
    }
}

template <std::size_t K>
constexpr uint64_t shift_left_append(uint64_t kmer, uint8_t base_bits) {
    return ((kmer << 2U) & kmer_mask<K>()) | static_cast<uint64_t>(base_bits & 0x3U);
}

template <std::size_t K>
constexpr uint64_t shift_right_prepend(uint64_t kmer, uint8_t base_bits) {
    return (kmer >> 2U) | (static_cast<uint64_t>(base_bits & 0x3U) << (2U * (K - 1U)));
}

}  // namespace bruijn::dna
