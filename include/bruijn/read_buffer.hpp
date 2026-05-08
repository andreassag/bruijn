#pragma once

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <string_view>
#include <type_traits>
#include <vector>

namespace bruijn {

struct ReadSpan {
    uint32_t offset;
    uint32_t length;
};

static_assert(std::is_standard_layout_v<ReadSpan>);

class ReadBuffer {
public:
    void add_read(std::string_view read);

    [[nodiscard]] std::size_t read_count() const noexcept { return spans_.size(); }
    [[nodiscard]] std::size_t total_bases() const noexcept { return bases_.size(); }

    [[nodiscard]] std::string_view read(std::size_t idx) const;

    static ReadBuffer from_fasta_fastq(const std::filesystem::path& path);

private:
    std::vector<uint8_t> bases_;
    std::vector<ReadSpan> spans_;
};

}  // namespace bruijn
