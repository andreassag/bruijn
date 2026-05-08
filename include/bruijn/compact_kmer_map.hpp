#pragma once

#include "bruijn/dna.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <type_traits>
#include <utility>
#include <vector>

namespace bruijn {

template <typename TValue>
class CompactKmerMap {
public:
    struct Slot {
        uint64_t key = 0;
        TValue value{};
        uint8_t used = 0;
    };

    CompactKmerMap() = default;

    explicit CompactKmerMap(std::size_t expected_items) {
        reserve(expected_items);
    }

    [[nodiscard]] std::size_t size() const noexcept { return size_; }
    [[nodiscard]] bool empty() const noexcept { return size_ == 0; }
    [[nodiscard]] std::size_t capacity() const noexcept { return slots_.size(); }

    void reserve(std::size_t expected_items) {
        if (expected_items == 0) {
            return;
        }

        std::size_t target = next_power_of_two(expected_items * 2U);
        if (target < 16U) {
            target = 16U;
        }

        if (target > slots_.size()) {
            rehash(target);
        }
    }

    TValue* find(uint64_t key) {
        if (slots_.empty()) {
            return nullptr;
        }

        const std::size_t m = mask();
        std::size_t idx = index_for(key);
        for (std::size_t probe = 0; probe <= m; ++probe) {
            Slot& slot = slots_[idx];
            if (!slot.used) {
                return nullptr;
            }
            if (slot.key == key) {
                return &slot.value;
            }
            idx = (idx + 1U) & m;
        }
        return nullptr;
    }

    const TValue* find(uint64_t key) const {
        return const_cast<CompactKmerMap*>(this)->find(key);
    }

    TValue& upsert(uint64_t key, const TValue& initial = TValue{}) {
        maybe_grow();
        const std::size_t m = mask();
        std::size_t idx = index_for(key);
        for (std::size_t probe = 0; probe <= m; ++probe) {
            Slot& slot = slots_[idx];
            if (!slot.used) {
                slot.used = 1;
                slot.key = key;
                slot.value = initial;
                ++size_;
                return slot.value;
            }
            if (slot.key == key) {
                return slot.value;
            }
            idx = (idx + 1U) & m;
        }
        // Should be unreachable because we grow before full.
        rehash(slots_.size() * 2U);
        return upsert(key, initial);
    }

    void set(uint64_t key, const TValue& value) {
        upsert(key) = value;
    }

    void increment(uint64_t key, TValue delta = static_cast<TValue>(1)) {
        TValue& value = upsert(key, static_cast<TValue>(0));
        value = static_cast<TValue>(value + delta);
    }

    [[nodiscard]] std::vector<std::pair<uint64_t, TValue>> items() const {
        std::vector<std::pair<uint64_t, TValue>> out;
        out.reserve(size_);
        for (const Slot& slot : slots_) {
            if (slot.used) {
                out.emplace_back(slot.key, slot.value);
            }
        }
        return out;
    }

    template <typename Fn>
    void for_each(Fn&& fn) {
        for (Slot& slot : slots_) {
            if (slot.used) {
                fn(slot.key, slot.value);
            }
        }
    }

    template <typename Fn>
    void for_each(Fn&& fn) const {
        for (const Slot& slot : slots_) {
            if (slot.used) {
                fn(slot.key, slot.value);
            }
        }
    }

private:
    std::vector<Slot> slots_;
    std::size_t size_ = 0;

    static std::size_t next_power_of_two(std::size_t v) {
        if (v <= 1U) {
            return 1U;
        }
        --v;
        for (std::size_t shift = 1U; shift < (sizeof(std::size_t) * 8U); shift <<= 1U) {
            v |= (v >> shift);
        }
        return v + 1U;
    }

    [[nodiscard]] std::size_t mask() const noexcept {
        return slots_.size() - 1U;
    }

    [[nodiscard]] std::size_t index_for(uint64_t key) const noexcept {
        return static_cast<std::size_t>(dna::splitmix64(key) & static_cast<uint64_t>(mask()));
    }

    void maybe_grow() {
        if (slots_.empty()) {
            rehash(16U);
            return;
        }

        const std::size_t load_times_100 = (size_ + 1U) * 100U / slots_.size();
        if (load_times_100 >= 70U) {
            rehash(slots_.size() * 2U);
        }
    }

    void rehash(std::size_t new_capacity) {
        if (new_capacity < 16U) {
            new_capacity = 16U;
        }
        new_capacity = next_power_of_two(new_capacity);

        std::vector<Slot> old_slots = std::move(slots_);
        slots_.assign(new_capacity, Slot{});
        size_ = 0;

        for (const Slot& slot : old_slots) {
            if (!slot.used) {
                continue;
            }
            upsert(slot.key, slot.value);
        }
    }
};

}  // namespace bruijn
