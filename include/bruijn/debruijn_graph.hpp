#pragma once

#include "bruijn/compact_kmer_map.hpp"
#include "bruijn/dna.hpp"
#include "bruijn/read_buffer.hpp"

#include <algorithm>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <deque>
#include <optional>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace bruijn {

template <std::size_t K>
class DeBruijnGraph {
public:
    struct Node {
        uint64_t kmer = 0;
        uint32_t count = 0;
        uint8_t in_mask = 0;
        uint8_t out_mask = 0;
        uint8_t alive = 1;
    };

    struct Contig {
        std::vector<char> bases;
        uint32_t mean_coverage = 0;
    };

    void build(const CompactKmerMap<uint32_t>& counts, const ReadBuffer& reads) {
        nodes_.clear();
        index_ = CompactKmerMap<uint32_t>(counts.size() * 2U + 1U);

        const auto items = counts.items();
        nodes_.reserve(items.size());
        for (const auto& item : items) {
            const uint32_t idx = static_cast<uint32_t>(nodes_.size());
            nodes_.push_back(Node{item.first, item.second, 0, 0, 1});
            index_.set(item.first, idx);
        }

        build_edges_from_reads(reads);
        prune_dead_edges();
    }

    [[nodiscard]] std::size_t node_count() const {
        std::size_t count = 0;
        for (const Node& node : nodes_) {
            if (node.alive) {
                ++count;
            }
        }
        return count;
    }

    [[nodiscard]] std::size_t edge_count() const {
        std::size_t count = 0;
        for (const Node& node : nodes_) {
            if (node.alive) {
                count += out_degree(node);
            }
        }
        return count;
    }

    void remove_tips(std::size_t max_tip_len) {
        if (max_tip_len == 0) {
            return;
        }

        bool changed = true;
        while (changed) {
            changed = false;
            prune_dead_edges();

            for (uint32_t idx = 0; idx < nodes_.size(); ++idx) {
                if (!nodes_[idx].alive) {
                    continue;
                }
                const Node& node = nodes_[idx];
                if (in_degree(node) == 0 && out_degree(node) == 1) {
                    const auto tip = collect_forward_tip(idx, max_tip_len);
                    if (!tip.empty()) {
                        for (uint32_t tip_idx : tip) {
                            nodes_[tip_idx].alive = 0;
                        }
                        changed = true;
                    }
                }
            }

            for (uint32_t idx = 0; idx < nodes_.size(); ++idx) {
                if (!nodes_[idx].alive) {
                    continue;
                }
                const Node& node = nodes_[idx];
                if (out_degree(node) == 0 && in_degree(node) == 1) {
                    const auto tip = collect_backward_tip(idx, max_tip_len);
                    if (!tip.empty()) {
                        for (uint32_t tip_idx : tip) {
                            nodes_[tip_idx].alive = 0;
                        }
                        changed = true;
                    }
                }
            }
        }

        prune_dead_edges();
    }

    void remove_bubbles(std::size_t max_depth) {
        if (max_depth == 0) {
            return;
        }

        constexpr std::size_t MAX_PATHS_PER_BRANCH = 128;
        constexpr std::size_t MAX_STATES_PER_BRANCH = 2048;

        bool changed = true;
        while (changed) {
            changed = false;
            prune_dead_edges();

            for (uint32_t idx = 0; idx < nodes_.size(); ++idx) {
                if (!nodes_[idx].alive) {
                    continue;
                }

                const auto succ = successors(idx);
                if (succ.size() < 2) {
                    continue;
                }

                std::unordered_map<uint32_t, std::vector<BubblePath>> by_endpoint;
                for (std::size_t branch_idx = 0; branch_idx < succ.size(); ++branch_idx) {
                    auto paths = enumerate_bubble_paths(
                        succ[branch_idx],
                        max_depth,
                        MAX_STATES_PER_BRANCH,
                        MAX_PATHS_PER_BRANCH);
                    for (BubblePath& path : paths) {
                        if (path.truncated || path.nodes.empty() || path.endpoint == idx) {
                            continue;
                        }
                        path.branch_index = branch_idx;
                        by_endpoint[path.endpoint].push_back(std::move(path));
                    }
                }

                bool removed_in_this_branch = false;
                for (auto& [endpoint, endpoint_paths] : by_endpoint) {
                    std::vector<const BubblePath*> best_per_branch(succ.size(), nullptr);
                    for (const BubblePath& path : endpoint_paths) {
                        const BubblePath*& best = best_per_branch[path.branch_index];
                        if (best == nullptr || better_path(path, *best)) {
                            best = &path;
                        }
                    }

                    std::size_t branch_hits = 0;
                    for (const BubblePath* path : best_per_branch) {
                        if (path != nullptr) {
                            ++branch_hits;
                        }
                    }
                    if (branch_hits < 2) {
                        continue;
                    }

                    const BubblePath* winner = nullptr;
                    for (const BubblePath* path : best_per_branch) {
                        if (path == nullptr) {
                            continue;
                        }
                        if (winner == nullptr || better_path(*path, *winner)) {
                            winner = path;
                        }
                    }

                    if (winner == nullptr) {
                        continue;
                    }

                    std::unordered_set<uint32_t> keep_nodes;
                    keep_nodes.reserve(winner->nodes.size() + 2U);
                    keep_nodes.insert(idx);
                    keep_nodes.insert(endpoint);
                    for (uint32_t kept_node : winner->nodes) {
                        keep_nodes.insert(kept_node);
                    }

                    bool removed_any = false;
                    for (const BubblePath* candidate : best_per_branch) {
                        if (candidate == nullptr || candidate == winner) {
                            continue;
                        }

                        for (uint32_t doomed_idx : candidate->nodes) {
                            if (doomed_idx == idx || doomed_idx == endpoint) {
                                continue;
                            }
                            if (keep_nodes.find(doomed_idx) != keep_nodes.end()) {
                                continue;
                            }
                            if (!nodes_[doomed_idx].alive) {
                                continue;
                            }
                            if (in_degree(nodes_[doomed_idx]) <= 1 && out_degree(nodes_[doomed_idx]) <= 1) {
                                nodes_[doomed_idx].alive = 0;
                                removed_any = true;
                            }
                        }
                    }

                    if (removed_any) {
                        removed_in_this_branch = true;
                        break;
                    }
                }

                if (removed_in_this_branch) {
                    changed = true;
                    break;
                }
            }
        }

        prune_dead_edges();
    }

    [[nodiscard]] std::vector<Contig> generate_contigs() const {
        std::vector<Contig> contigs;
        if (nodes_.empty()) {
            return contigs;
        }

        std::vector<uint8_t> visited_out(nodes_.size(), 0);

        auto walk = [&](uint32_t start_idx, uint8_t start_base) {
            Contig contig;
            dna::decode_kmer<K>(nodes_[start_idx].kmer, contig.bases);

            uint64_t coverage_sum = nodes_[start_idx].count;
            uint64_t nodes_in_path = 1;

            uint32_t current_idx = start_idx;
            uint8_t edge_base = start_base;
            std::size_t safety = 0;
            const std::size_t safety_limit = nodes_.size() * 4U + K;

            while (safety++ < safety_limit) {
                const uint8_t edge_mask = static_cast<uint8_t>(1U << edge_base);
                if ((visited_out[current_idx] & edge_mask) != 0) {
                    break;
                }
                visited_out[current_idx] |= edge_mask;

                const std::optional<uint32_t> next_idx = successor_for_base(current_idx, edge_base);
                if (!next_idx.has_value()) {
                    break;
                }

                if (!nodes_[*next_idx].alive) {
                    break;
                }

                contig.bases.push_back(dna::bits_to_base(edge_base));
                coverage_sum += nodes_[*next_idx].count;
                ++nodes_in_path;

                const Node& next_node = nodes_[*next_idx];
                if (in_degree(next_node) != 1 || out_degree(next_node) != 1) {
                    break;
                }

                const auto next_base = single_out_base(next_node);
                if (!next_base.has_value()) {
                    break;
                }

                current_idx = *next_idx;
                edge_base = *next_base;
            }

            contig.mean_coverage = static_cast<uint32_t>(coverage_sum / std::max<uint64_t>(nodes_in_path, 1));
            return contig;
        };

        for (uint32_t idx = 0; idx < nodes_.size(); ++idx) {
            if (!nodes_[idx].alive) {
                continue;
            }

            const Node& node = nodes_[idx];
            if (out_degree(node) == 0) {
                continue;
            }

            if (in_degree(node) == 1 && out_degree(node) == 1) {
                continue;
            }

            for (uint8_t base = 0; base < 4U; ++base) {
                const uint8_t mask = static_cast<uint8_t>(1U << base);
                if ((node.out_mask & mask) == 0) {
                    continue;
                }
                if ((visited_out[idx] & mask) != 0) {
                    continue;
                }
                contigs.push_back(walk(idx, base));
            }
        }

        // Cover isolated cycles where every node has in=1 and out=1.
        for (uint32_t idx = 0; idx < nodes_.size(); ++idx) {
            if (!nodes_[idx].alive) {
                continue;
            }

            for (uint8_t base = 0; base < 4U; ++base) {
                const uint8_t mask = static_cast<uint8_t>(1U << base);
                if ((nodes_[idx].out_mask & mask) == 0) {
                    continue;
                }
                if ((visited_out[idx] & mask) != 0) {
                    continue;
                }
                contigs.push_back(walk(idx, base));
            }
        }

        std::sort(contigs.begin(), contigs.end(), [](const Contig& a, const Contig& b) {
            return a.bases.size() > b.bases.size();
        });

        return contigs;
    }

private:
    struct BubblePath {
        std::vector<uint32_t> nodes;
        uint32_t endpoint = 0;
        uint64_t coverage_sum = 0;
        std::size_t branch_index = 0;
        bool truncated = false;
    };

    struct BubbleSearchState {
        uint32_t node = 0;
        std::vector<uint32_t> path;
        uint64_t coverage_sum = 0;
        std::size_t depth = 0;
    };

    std::vector<Node> nodes_;
    CompactKmerMap<uint32_t> index_;

    [[nodiscard]] static uint8_t in_degree(const Node& node) {
        return static_cast<uint8_t>(std::popcount(static_cast<unsigned int>(node.in_mask)));
    }

    [[nodiscard]] static uint8_t out_degree(const Node& node) {
        return static_cast<uint8_t>(std::popcount(static_cast<unsigned int>(node.out_mask)));
    }

    [[nodiscard]] std::optional<uint32_t> node_index(uint64_t kmer) const {
        const uint32_t* idx = index_.find(kmer);
        if (idx == nullptr) {
            return std::nullopt;
        }
        return *idx;
    }

    [[nodiscard]] std::optional<uint8_t> single_out_base(const Node& node) const {
        if (out_degree(node) != 1) {
            return std::nullopt;
        }
        for (uint8_t base = 0; base < 4U; ++base) {
            if ((node.out_mask & static_cast<uint8_t>(1U << base)) != 0) {
                return base;
            }
        }
        return std::nullopt;
    }

    [[nodiscard]] std::optional<uint8_t> single_in_base(const Node& node) const {
        if (in_degree(node) != 1) {
            return std::nullopt;
        }
        for (uint8_t base = 0; base < 4U; ++base) {
            if ((node.in_mask & static_cast<uint8_t>(1U << base)) != 0) {
                return base;
            }
        }
        return std::nullopt;
    }

    [[nodiscard]] std::optional<uint32_t> successor_for_base(uint32_t idx, uint8_t base) const {
        const Node& node = nodes_[idx];
        if (!node.alive) {
            return std::nullopt;
        }
        if ((node.out_mask & static_cast<uint8_t>(1U << base)) == 0) {
            return std::nullopt;
        }

        const uint64_t next_kmer = dna::shift_left_append<K>(node.kmer, base);
        const auto next_idx = node_index(next_kmer);
        if (!next_idx.has_value() || !nodes_[*next_idx].alive) {
            return std::nullopt;
        }
        return *next_idx;
    }

    [[nodiscard]] std::optional<uint32_t> predecessor_for_base(uint32_t idx, uint8_t base) const {
        const Node& node = nodes_[idx];
        if (!node.alive) {
            return std::nullopt;
        }
        if ((node.in_mask & static_cast<uint8_t>(1U << base)) == 0) {
            return std::nullopt;
        }

        const uint64_t prev_kmer = dna::shift_right_prepend<K>(node.kmer, base);
        const auto prev_idx = node_index(prev_kmer);
        if (!prev_idx.has_value() || !nodes_[*prev_idx].alive) {
            return std::nullopt;
        }
        return *prev_idx;
    }

    [[nodiscard]] std::vector<uint32_t> successors(uint32_t idx) const {
        std::vector<uint32_t> out;
        out.reserve(4);
        for (uint8_t base = 0; base < 4U; ++base) {
            const auto next = successor_for_base(idx, base);
            if (next.has_value()) {
                out.push_back(*next);
            }
        }
        return out;
    }

    void add_edge(uint32_t from_idx, uint8_t outgoing_base, uint32_t to_idx, uint8_t incoming_base) {
        nodes_[from_idx].out_mask |= static_cast<uint8_t>(1U << outgoing_base);
        nodes_[to_idx].in_mask |= static_cast<uint8_t>(1U << incoming_base);
    }

    void build_edges_from_reads(const ReadBuffer& reads) {
        for (std::size_t r = 0; r < reads.read_count(); ++r) {
            std::optional<uint32_t> prev_idx;
            std::size_t prev_start = 0;

            const std::string_view read = reads.read(r);
            dna::for_each_kmer<K>(read, [&](uint64_t kmer, std::size_t start) {
                const auto curr_idx = node_index(kmer);
                if (!curr_idx.has_value()) {
                    prev_idx.reset();
                    prev_start = start;
                    return;
                }

                if (prev_idx.has_value() && start == prev_start + 1U) {
                    const uint8_t outgoing_base = static_cast<uint8_t>(kmer & 0x3U);
                    const uint8_t incoming_base = static_cast<uint8_t>((nodes_[*prev_idx].kmer >> (2U * (K - 1U))) & 0x3U);
                    add_edge(*prev_idx, outgoing_base, *curr_idx, incoming_base);
                }

                prev_idx = curr_idx;
                prev_start = start;
            });
        }
    }

    void prune_dead_edges() {
        for (Node& node : nodes_) {
            if (!node.alive) {
                node.in_mask = 0;
                node.out_mask = 0;
                continue;
            }

            uint8_t in_mask = 0;
            uint8_t out_mask = 0;

            for (uint8_t base = 0; base < 4U; ++base) {
                if ((node.out_mask & static_cast<uint8_t>(1U << base)) != 0) {
                    const uint64_t next_kmer = dna::shift_left_append<K>(node.kmer, base);
                    const auto next_idx = node_index(next_kmer);
                    if (next_idx.has_value() && nodes_[*next_idx].alive) {
                        out_mask |= static_cast<uint8_t>(1U << base);
                    }
                }

                if ((node.in_mask & static_cast<uint8_t>(1U << base)) != 0) {
                    const uint64_t prev_kmer = dna::shift_right_prepend<K>(node.kmer, base);
                    const auto prev_idx = node_index(prev_kmer);
                    if (prev_idx.has_value() && nodes_[*prev_idx].alive) {
                        in_mask |= static_cast<uint8_t>(1U << base);
                    }
                }
            }

            node.in_mask = in_mask;
            node.out_mask = out_mask;
        }
    }

    [[nodiscard]] std::vector<uint32_t> collect_forward_tip(uint32_t start, std::size_t max_tip_len) const {
        if (!nodes_[start].alive) {
            return {};
        }

        std::vector<uint32_t> path;
        path.reserve(max_tip_len);

        uint32_t current = start;
        for (std::size_t depth = 0; depth < max_tip_len; ++depth) {
            const Node& node = nodes_[current];
            path.push_back(current);

            if (out_degree(node) == 0) {
                return path;
            }
            if (out_degree(node) != 1) {
                return {};
            }

            const auto base = single_out_base(node);
            if (!base.has_value()) {
                return {};
            }

            const auto next = successor_for_base(current, *base);
            if (!next.has_value()) {
                return path;
            }

            const Node& next_node = nodes_[*next];
            if (in_degree(next_node) != 1) {
                return path;
            }

            current = *next;
            for (uint32_t seen : path) {
                if (seen == current) {
                    return {};
                }
            }
        }

        return {};
    }

    [[nodiscard]] std::vector<uint32_t> collect_backward_tip(uint32_t start, std::size_t max_tip_len) const {
        if (!nodes_[start].alive) {
            return {};
        }

        std::vector<uint32_t> path;
        path.reserve(max_tip_len);

        uint32_t current = start;
        for (std::size_t depth = 0; depth < max_tip_len; ++depth) {
            const Node& node = nodes_[current];
            path.push_back(current);

            if (in_degree(node) == 0) {
                return path;
            }
            if (in_degree(node) != 1) {
                return {};
            }

            const auto base = single_in_base(node);
            if (!base.has_value()) {
                return {};
            }

            const auto prev = predecessor_for_base(current, *base);
            if (!prev.has_value()) {
                return path;
            }

            const Node& prev_node = nodes_[*prev];
            if (out_degree(prev_node) != 1) {
                return path;
            }

            current = *prev;
            for (uint32_t seen : path) {
                if (seen == current) {
                    return {};
                }
            }
        }

        return {};
    }

    [[nodiscard]] static bool better_path(const BubblePath& lhs, const BubblePath& rhs) {
        const long double lhs_score = static_cast<long double>(lhs.coverage_sum) /
                                      static_cast<long double>(std::max<std::size_t>(lhs.nodes.size(), 1U));
        const long double rhs_score = static_cast<long double>(rhs.coverage_sum) /
                                      static_cast<long double>(std::max<std::size_t>(rhs.nodes.size(), 1U));
        if (lhs_score != rhs_score) {
            return lhs_score > rhs_score;
        }
        if (lhs.nodes.size() != rhs.nodes.size()) {
            return lhs.nodes.size() > rhs.nodes.size();
        }
        return lhs.coverage_sum > rhs.coverage_sum;
    }

    [[nodiscard]] std::vector<BubblePath> enumerate_bubble_paths(
        uint32_t start,
        std::size_t max_depth,
        std::size_t max_states,
        std::size_t max_paths) const {
        std::vector<BubblePath> paths;
        if (start >= nodes_.size() || !nodes_[start].alive) {
            return paths;
        }

        std::deque<BubbleSearchState> queue;
        BubbleSearchState initial;
        initial.node = start;
        initial.path.push_back(start);
        initial.coverage_sum = nodes_[start].count;
        queue.push_back(std::move(initial));

        std::size_t expanded = 0;
        while (!queue.empty() && expanded < max_states && paths.size() < max_paths) {
            BubbleSearchState state = std::move(queue.front());
            queue.pop_front();
            ++expanded;

            if (!nodes_[state.node].alive) {
                continue;
            }

            if (state.depth >= max_depth) {
                paths.push_back(BubblePath{std::move(state.path), state.node, state.coverage_sum, 0, true});
                continue;
            }

            const Node& node = nodes_[state.node];
            if (state.depth > 0 && in_degree(node) > 1) {
                paths.push_back(BubblePath{std::move(state.path), state.node, state.coverage_sum, 0, false});
                continue;
            }

            const auto next_nodes = successors(state.node);
            if (next_nodes.empty()) {
                paths.push_back(BubblePath{std::move(state.path), state.node, state.coverage_sum, 0, false});
                continue;
            }

            bool pushed_any = false;
            for (uint32_t next_idx : next_nodes) {
                if (!nodes_[next_idx].alive) {
                    continue;
                }
                if (std::find(state.path.begin(), state.path.end(), next_idx) != state.path.end()) {
                    continue;
                }

                BubbleSearchState next;
                next.node = next_idx;
                next.path = state.path;
                next.path.push_back(next_idx);
                next.coverage_sum = state.coverage_sum + nodes_[next_idx].count;
                next.depth = state.depth + 1U;
                queue.push_back(std::move(next));
                pushed_any = true;
            }

            if (!pushed_any) {
                paths.push_back(BubblePath{std::move(state.path), state.node, state.coverage_sum, 0, true});
            }
        }

        return paths;
    }
};

}  // namespace bruijn
