#include "bruijn/read_buffer.hpp"

#include "bruijn/dna.hpp"

#include <fstream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace bruijn {

void ReadBuffer::add_read(std::string_view read) {
    if (read.empty()) {
        return;
    }

    if (read.size() > static_cast<std::size_t>(std::numeric_limits<uint32_t>::max())) {
        throw std::runtime_error("Read is too long to fit into uint32_t span");
    }

    const std::size_t offset = bases_.size();
    bases_.reserve(bases_.size() + read.size());
    for (char c : read) {
        bases_.push_back(dna::normalize_base(c));
    }

    if (offset > static_cast<std::size_t>(std::numeric_limits<uint32_t>::max())) {
        throw std::runtime_error("Total sequence storage exceeds uint32_t addressable range");
    }

    spans_.push_back(ReadSpan{
        static_cast<uint32_t>(offset),
        static_cast<uint32_t>(read.size()),
    });
}

std::string_view ReadBuffer::read(std::size_t idx) const {
    if (idx >= spans_.size()) {
        throw std::out_of_range("Read index out of range");
    }

    const ReadSpan span = spans_[idx];
    const char* ptr = reinterpret_cast<const char*>(bases_.data() + span.offset);
    return std::string_view(ptr, span.length);
}

ReadBuffer ReadBuffer::from_fasta_fastq(const std::filesystem::path& path) {
    std::ifstream input(path);
    if (!input) {
        throw std::runtime_error("Failed to open input file: " + path.string());
    }

    enum class Mode { Unknown, Fasta, Fastq };

    ReadBuffer buffer;
    Mode mode = Mode::Unknown;
    std::vector<char> fasta_accumulator;
    std::string line;

    auto trim_cr = [](std::string& s) {
        if (!s.empty() && s.back() == '\r') {
            s.pop_back();
        }
    };

    auto flush_fasta = [&]() {
        if (!fasta_accumulator.empty()) {
            buffer.add_read(std::string_view(fasta_accumulator.data(), fasta_accumulator.size()));
            fasta_accumulator.clear();
        }
    };

    while (std::getline(input, line)) {
        trim_cr(line);
        if (line.empty()) {
            continue;
        }

        if (line[0] == '>') {
            if (mode != Mode::Fastq) {
                flush_fasta();
                mode = Mode::Fasta;
                continue;
            }
        }

        if (mode == Mode::Fasta) {
            fasta_accumulator.insert(fasta_accumulator.end(), line.begin(), line.end());
            continue;
        }

        if (line[0] == '@') {
            mode = Mode::Fastq;

            std::string sequence;
            std::string plus;
            if (!std::getline(input, sequence) || !std::getline(input, plus)) {
                throw std::runtime_error("Malformed FASTQ entry near header: " + line);
            }
            trim_cr(sequence);
            trim_cr(plus);
            if (plus.empty() || plus[0] != '+') {
                throw std::runtime_error("Malformed FASTQ plus line near header: " + line);
            }

            std::size_t qual_remaining = sequence.size();
            std::string quality;
            while (qual_remaining > 0 && std::getline(input, quality)) {
                trim_cr(quality);
                if (quality.size() >= qual_remaining) {
                    qual_remaining = 0;
                } else {
                    qual_remaining -= quality.size();
                }
            }

            if (qual_remaining != 0) {
                throw std::runtime_error("Malformed FASTQ quality section near header: " + line);
            }

            buffer.add_read(sequence);
            continue;
        }

        if (mode == Mode::Unknown) {
            buffer.add_read(line);
            continue;
        }

        buffer.add_read(line);
    }

    flush_fasta();
    return buffer;
}

}  // namespace bruijn
