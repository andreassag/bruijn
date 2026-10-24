#include "bruijn/assembler.hpp"
#include "bruijn/debruijn_graph.hpp"
#include "bruijn/dna.hpp"
#include "bruijn/kmer_counter.hpp"
#include "bruijn/read_buffer.hpp"

#include <cstddef>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace {

template <std::size_t K>
std::string decode_to_string(uint64_t kmer) {
    std::vector<char> decoded;
    bruijn::dna::decode_kmer<K>(kmer, decoded);
    return std::string(decoded.begin(), decoded.end());
}

template <std::size_t K>
uint64_t encode_exact(std::string_view sequence) {
    uint64_t encoded = 0;
    bool found = false;
    bruijn::dna::for_each_kmer<K>(sequence, [&](uint64_t kmer, std::size_t start) {
        if (start == 0 && !found) {
            encoded = kmer;
            found = true;
        }
    });

    if (!found || sequence.size() != K) {
        throw std::runtime_error("encode_exact received invalid k-mer input");
    }
    return encoded;
}

#define CHECK_OR_FAIL(condition, message)         \
    do {                                           \
        if (!(condition)) {                        \
            std::cerr << "[FAIL] " << message     \
                      << "\n";                   \
            return false;                          \
        }                                          \
    } while (false)

bool test_kmer_extraction_and_canonicalization() {
    std::vector<std::string> observed;
    bruijn::dna::for_each_kmer<5>("AACGTNACGTA", [&](uint64_t kmer, std::size_t) {
        observed.push_back(decode_to_string<5>(kmer));
    });

    CHECK_OR_FAIL(observed.size() == 2, "expected two valid 5-mers after splitting on N");
    CHECK_OR_FAIL(observed[0] == "AACGT", "first extracted 5-mer mismatch");
    CHECK_OR_FAIL(observed[1] == "ACGTA", "second extracted 5-mer mismatch");

    const uint64_t forward = encode_exact<5>("ACGTT");
    const uint64_t reverse_complement = encode_exact<5>("AACGT");
    CHECK_OR_FAIL(
        bruijn::dna::canonical_kmer<5>(forward) == bruijn::dna::canonical_kmer<5>(reverse_complement),
        "canonical key should collapse reverse-complement pairs");

    bruijn::ReadBuffer reads;
    reads.add_read("ACGTT");
    reads.add_read("AACGT");

    bruijn::CounterConfig config;
    config.min_count = 2;
    config.thread_count = 1;
    config.bloom_counters = 1U << 12U;
    config.bloom_hashes = 3;
    config.canonicalize = true;

    bruijn::KmerCounter<5> counter(config);
    const auto counts = counter.count(reads);

    const uint32_t* forward_count = counts.find(forward);
    const uint32_t* reverse_count = counts.find(reverse_complement);
    CHECK_OR_FAIL(counts.size() == 2, "canonical filtering should retain both directional nodes");
    CHECK_OR_FAIL(forward_count != nullptr && reverse_count != nullptr, "expected directional keys missing");
    CHECK_OR_FAIL(*forward_count == 1 && *reverse_count == 1, "directional counts should be preserved after canonical filter");

    return true;
}

bool test_tip_clipping() {
    bruijn::ReadBuffer reads;
    for (std::size_t i = 0; i < 4; ++i) {
        reads.add_read("AAAACCCCGGGGAAAA");
    }
    reads.add_read("AAAACCCCTTTT");

    bruijn::CounterConfig config;
    config.min_count = 1;
    config.thread_count = 1;

    bruijn::KmerCounter<4> counter(config);
    const auto counts = counter.count(reads);

    bruijn::DeBruijnGraph<4> graph;
    graph.build(counts, reads);
    const std::size_t before = graph.node_count();

    graph.remove_tips(5);
    const std::size_t after = graph.node_count();
    CHECK_OR_FAIL(after < before, "tip clipping should remove low-support dead-end nodes");

    const auto contigs = graph.generate_contigs();
    bool has_tip_artifact = false;
    std::size_t longest = 0;
    for (const auto& contig : contigs) {
        const std::string sequence(contig.bases.begin(), contig.bases.end());
        longest = std::max(longest, sequence.size());
        if (sequence.find("TTTT") != std::string::npos) {
            has_tip_artifact = true;
        }
    }

    CHECK_OR_FAIL(!has_tip_artifact, "tip branch should be removed from output contigs");
    CHECK_OR_FAIL(longest >= 14, "main path should still assemble into a long contig after clipping");

    return true;
}

bool test_bubble_popping_multi_branch() {
    bruijn::ReadBuffer reads;
    for (std::size_t i = 0; i < 8; ++i) {
        reads.add_read("AACGAT");
    }
    reads.add_read("AACTGAT");
    reads.add_read("AACAAGAT");

    bruijn::CounterConfig config;
    config.min_count = 1;
    config.thread_count = 1;

    bruijn::KmerCounter<3> counter(config);
    const auto counts = counter.count(reads);

    bruijn::DeBruijnGraph<3> graph;
    graph.build(counts, reads);
    const std::size_t before = graph.node_count();

    graph.remove_bubbles(6);
    const std::size_t after = graph.node_count();
    CHECK_OR_FAIL(after < before, "bubble popping should remove lower-support branches");

    const auto contigs = graph.generate_contigs();
    bool has_main = false;
    bool has_alt_1 = false;
    bool has_alt_2 = false;

    for (const auto& contig : contigs) {
        const std::string sequence(contig.bases.begin(), contig.bases.end());
        if (sequence.find("AACGAT") != std::string::npos) {
            has_main = true;
        }
        if (sequence.find("AACTGAT") != std::string::npos) {
            has_alt_1 = true;
        }
        if (sequence.find("AACAAGAT") != std::string::npos) {
            has_alt_2 = true;
        }
    }

    CHECK_OR_FAIL(has_main, "highest-coverage bubble path should remain");
    CHECK_OR_FAIL(!has_alt_1 && !has_alt_2, "lower-support bubble paths should be removed");

    return true;
}

bool test_toy_assembly() {
    bruijn::ReadBuffer reads;
    // 17-bp reads with 16-bp overlap assembling into 22-bp target: TACGATCGATCGATCGAGCTAG
    reads.add_read("TACGATCGATCGATCGA");
    reads.add_read("ACGATCGATCGATCGAG");
    reads.add_read("CGATCGATCGATCGAGC");
    reads.add_read("GATCGATCGATCGAGCT");
    reads.add_read("ATCGATCGATCGAGCTA");
    reads.add_read("TCGATCGATCGAGCTAG");
    // Duplicate terminal reads to ensure all k-mers have count >= 2
    reads.add_read("TACGATCGATCGATCGA");
    reads.add_read("TCGATCGATCGAGCTAG");

    bruijn::AssemblerConfig config;
    config.min_count = 2;
    config.thread_count = 1;
    config.bloom_counters = 1U << 12U;
    config.bloom_hashes = 3;
    config.canonical_kmers = false;
    config.max_tip_length = 4;
    config.bubble_depth = 16;

    bruijn::DeBruijnAssembler<15> assembler(config);
    const auto result = assembler.assemble(reads);

    CHECK_OR_FAIL(result.retained_kmers == 8, "toy assembly should retain exactly 8 k-mers");
    CHECK_OR_FAIL(result.node_count == 8, "toy graph should contain exactly 8 nodes");
    CHECK_OR_FAIL(result.edge_count == 7, "toy graph should form a linear chain of 7 edges");
    CHECK_OR_FAIL(result.contigs.size() == 1, "toy assembly should yield exactly 1 contig");

    const std::string contig_seq(result.contigs[0].bases.begin(), result.contigs[0].bases.end());
    CHECK_OR_FAIL(contig_seq == "TACGATCGATCGATCGAGCTAG", "assembled contig sequence mismatch");
    CHECK_OR_FAIL(result.contigs[0].bases.size() == 22, "contig length should be exactly 22 bp");
    CHECK_OR_FAIL(result.contigs[0].mean_coverage == 3, "contig mean coverage should be exactly 3");

    // Test with canonical k-mers enabled
    bruijn::AssemblerConfig canonical_config = config;
    canonical_config.canonical_kmers = true;
    bruijn::DeBruijnAssembler<15> canonical_assembler(canonical_config);
    const auto canonical_result = canonical_assembler.assemble(reads);
    CHECK_OR_FAIL(canonical_result.contigs.size() == 1, "canonical assembly should yield 1 contig");
    CHECK_OR_FAIL(canonical_result.contigs[0].bases.size() == 22, "canonical contig length should be 22 bp");
    CHECK_OR_FAIL(canonical_result.contigs[0].mean_coverage == 3, "canonical mean coverage should be 3");
    const std::string canonical_seq(canonical_result.contigs[0].bases.begin(), canonical_result.contigs[0].bases.end());
    CHECK_OR_FAIL(canonical_seq == "TACGATCGATCGATCGAGCTAG", "canonical assembled contig sequence mismatch");

    // Test with multi-threading (2 threads)
    bruijn::AssemblerConfig mt_config = config;
    mt_config.thread_count = 2;
    bruijn::DeBruijnAssembler<15> mt_assembler(mt_config);
    const auto mt_result = mt_assembler.assemble(reads);
    CHECK_OR_FAIL(mt_result.contigs.size() == 1, "multi-threaded assembly should yield 1 contig");
    const std::string mt_seq(mt_result.contigs[0].bases.begin(), mt_result.contigs[0].bases.end());
    CHECK_OR_FAIL(mt_seq == "TACGATCGATCGATCGAGCTAG", "multi-threaded contig sequence mismatch");

    // Also test loading directly from FASTA file
    const auto file_reads = bruijn::ReadBuffer::from_fasta_fastq("data/toy_reads.fasta");
    const auto file_result = assembler.assemble(file_reads);
    CHECK_OR_FAIL(file_result.contigs.size() == 1, "assembly from FASTA file should yield 1 contig");
    const std::string file_contig_seq(file_result.contigs[0].bases.begin(), file_result.contigs[0].bases.end());
    CHECK_OR_FAIL(file_contig_seq == "TACGATCGATCGATCGAGCTAG", "FASTA assembled contig sequence mismatch");
    CHECK_OR_FAIL(file_result.contigs[0].mean_coverage == 3, "FASTA contig mean coverage should be 3");

    return true;
}

}  // namespace

int main() {
    const std::vector<std::pair<const char*, bool (*)()>> tests = {
        {"k-mer extraction and canonicalization", test_kmer_extraction_and_canonicalization},
        {"tip clipping", test_tip_clipping},
        {"multi-branch bubble popping", test_bubble_popping_multi_branch},
        {"toy dataset end-to-end assembly", test_toy_assembly},
    };

    bool all_ok = true;
    for (const auto& test : tests) {
        const bool ok = test.second();
        if (ok) {
            std::cout << "[PASS] " << test.first << "\n";
        } else {
            all_ok = false;
        }
    }

    return all_ok ? 0 : 1;
}
