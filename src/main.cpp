#include "bruijn/assembler.hpp"
#include "bruijn/read_buffer.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <thread>
#include <vector>

namespace {

struct CliOptions {
    std::filesystem::path input_path;
    std::filesystem::path output_path = "contigs.fasta";
    std::size_t k = 21;
    uint32_t min_count = 2;
    std::size_t threads = std::thread::hardware_concurrency();
    std::size_t bloom_counters = 1U << 24U;
    uint8_t bloom_hashes = 4;
    bool canonical_kmers = false;
    std::size_t tip_length = 8;
    std::size_t bubble_depth = 64;
};

void print_usage() {
    std::cout << "Usage:\n"
              << "  bruijn_assembler --input reads.fasta --k 21 [options]\n\n"
              << "Options:\n"
              << "  --output <path>          Output FASTA path (default: contigs.fasta)\n"
              << "  --k <odd integer>        K-mer size (supported: 15,17,19,21,23,25,27,29,31)\n"
              << "  --min-count <int>        K-mer count threshold (default: 2)\n"
              << "  --threads <int>          Worker threads (default: hardware concurrency)\n"
              << "  --bloom-counters <int>   Counting Bloom counters (default: 16777216)\n"
              << "  --bloom-hashes <int>     Bloom hash functions (default: 4)\n"
              << "  --canonical-kmers        Collapse reverse complements for k-mer filtering\n"
              << "  --tip-length <int>       Maximum tip length to clip (default: 8)\n"
              << "  --bubble-depth <int>     Max depth for bubble detection (default: 64)\n"
              << "  --help                   Show this message\n";
}

std::size_t parse_size(std::string_view value, const char* flag_name) {
    std::size_t parsed = 0;
    try {
        parsed = static_cast<std::size_t>(std::stoull(std::string(value)));
    } catch (const std::exception&) {
        throw std::runtime_error(std::string("Invalid numeric value for ") + flag_name + ": " + std::string(value));
    }
    return parsed;
}

CliOptions parse_args(int argc, char** argv) {
    CliOptions options;

    for (int i = 1; i < argc; ++i) {
        const std::string_view arg = argv[i];
        auto require_value = [&](const char* flag) -> std::string_view {
            if (i + 1 >= argc) {
                throw std::runtime_error(std::string("Missing value for ") + flag);
            }
            ++i;
            return argv[i];
        };

        if (arg == "--help") {
            print_usage();
            std::exit(0);
        }
        if (arg == "--input") {
            options.input_path = std::string(require_value("--input"));
            continue;
        }
        if (arg == "--output") {
            options.output_path = std::string(require_value("--output"));
            continue;
        }
        if (arg == "--k") {
            options.k = parse_size(require_value("--k"), "--k");
            continue;
        }
        if (arg == "--min-count") {
            options.min_count = static_cast<uint32_t>(parse_size(require_value("--min-count"), "--min-count"));
            continue;
        }
        if (arg == "--threads") {
            options.threads = parse_size(require_value("--threads"), "--threads");
            continue;
        }
        if (arg == "--bloom-counters") {
            options.bloom_counters = parse_size(require_value("--bloom-counters"), "--bloom-counters");
            continue;
        }
        if (arg == "--bloom-hashes") {
            options.bloom_hashes = static_cast<uint8_t>(parse_size(require_value("--bloom-hashes"), "--bloom-hashes"));
            continue;
        }
        if (arg == "--canonical-kmers") {
            options.canonical_kmers = true;
            continue;
        }
        if (arg == "--tip-length") {
            options.tip_length = parse_size(require_value("--tip-length"), "--tip-length");
            continue;
        }
        if (arg == "--bubble-depth") {
            options.bubble_depth = parse_size(require_value("--bubble-depth"), "--bubble-depth");
            continue;
        }

        throw std::runtime_error("Unknown argument: " + std::string(arg));
    }

    if (options.input_path.empty()) {
        throw std::runtime_error("Missing required --input argument");
    }

    if (options.k < 3 || options.k > 31 || options.k % 2 == 0) {
        throw std::runtime_error("K must be an odd integer between 3 and 31");
    }

    return options;
}

template <typename Contigs>
void write_contigs_fasta(const std::filesystem::path& output, const Contigs& contigs) {
    std::ofstream out(output);
    if (!out) {
        throw std::runtime_error("Failed to open output file: " + output.string());
    }

    for (std::size_t i = 0; i < contigs.size(); ++i) {
        const auto& contig = contigs[i];
        out << ">contig_" << (i + 1U)
            << " len=" << contig.bases.size()
            << " cov=" << contig.mean_coverage
            << "\n";

        constexpr std::size_t line_width = 80;
        for (std::size_t pos = 0; pos < contig.bases.size(); pos += line_width) {
            const std::size_t n = std::min(line_width, contig.bases.size() - pos);
            out.write(contig.bases.data() + pos, static_cast<std::streamsize>(n));
            out << '\n';
        }
    }
}

template <std::size_t K>
int run_assembler(const CliOptions& options, const bruijn::ReadBuffer& reads) {
    bruijn::AssemblerConfig config;
    config.min_count = options.min_count;
    config.thread_count = options.threads;
    config.bloom_counters = options.bloom_counters;
    config.bloom_hashes = options.bloom_hashes;
    config.canonical_kmers = options.canonical_kmers;
    config.max_tip_length = options.tip_length;
    config.bubble_depth = options.bubble_depth;

    bruijn::DeBruijnAssembler<K> assembler(config);
    auto result = assembler.assemble(reads);

    write_contigs_fasta(options.output_path, result.contigs);

    std::cout << "Assembled with k=" << K << '\n'
              << "Reads: " << reads.read_count() << '\n'
              << "Total bases: " << reads.total_bases() << '\n'
              << "Retained k-mers: " << result.retained_kmers << '\n'
              << "Nodes: " << result.node_count << '\n'
              << "Edges: " << result.edge_count << '\n'
              << "Contigs: " << result.contigs.size() << '\n'
              << "Output: " << options.output_path << '\n';

    return 0;
}

template <typename Fn>
bool dispatch_supported_k(std::size_t k, Fn&& fn) {
    switch (k) {
        case 15:
            fn.template operator()<15>();
            return true;
        case 17:
            fn.template operator()<17>();
            return true;
        case 19:
            fn.template operator()<19>();
            return true;
        case 21:
            fn.template operator()<21>();
            return true;
        case 23:
            fn.template operator()<23>();
            return true;
        case 25:
            fn.template operator()<25>();
            return true;
        case 27:
            fn.template operator()<27>();
            return true;
        case 29:
            fn.template operator()<29>();
            return true;
        case 31:
            fn.template operator()<31>();
            return true;
        default:
            return false;
    }
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const CliOptions options = parse_args(argc, argv);
        const bruijn::ReadBuffer reads = bruijn::ReadBuffer::from_fasta_fastq(options.input_path);

        int exit_code = 0;
        const bool supported = dispatch_supported_k(options.k, [&]<std::size_t K>() {
            exit_code = run_assembler<K>(options, reads);
        });

        if (!supported) {
            throw std::runtime_error("Unsupported k-mer size. Supported values: 15,17,19,21,23,25,27,29,31");
        }

        return exit_code;
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << '\n';
        return 1;
    }
}
