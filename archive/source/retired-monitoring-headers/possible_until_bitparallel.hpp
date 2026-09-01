#pragma once

#include <algorithm>
#include <array>
#include <bitset>
#include <cstdint>
#include <vector>

inline int untilTrailingZeroBits(std::uint64_t word)
{
    int count = 0;
    while ((word & 1) == 0)
    {
        word >>= 1;
        count++;
    }
    return count;
}

template <std::size_t N>
class UntilPackedGrid
{
public:
    UntilPackedGrid(int rows, int columns)
        : rows_(rows), columns_(columns), words_((columns + 63) / 64),
          data_(static_cast<std::size_t>(rows) * 4 * words_, 0)
    {
    }

    std::uint64_t *bits(int row, int flags)
    {
        return data_.data() +
               (static_cast<std::size_t>(row) * 4 + flags) * words_;
    }

    const std::uint64_t *bits(int row, int flags) const
    {
        return data_.data() +
               (static_cast<std::size_t>(row) * 4 + flags) * words_;
    }

    int rows() const
    {
        return rows_;
    }

    int columns() const
    {
        return columns_;
    }

    int words() const
    {
        return words_;
    }

    void clear()
    {
        std::fill(data_.begin(), data_.end(), 0);
    }

    void swap(UntilPackedGrid &other)
    {
        data_.swap(other.data_);
    }

private:
    int rows_;
    int columns_;
    int words_;
    std::vector<std::uint64_t> data_;
};

template <std::size_t N>
std::array<bool, 2> possibleUntilBitsBitParallel(
    const std::vector<std::bitset<N>> &lhs,
    const std::vector<std::bitset<N>> &rhs,
    const TimedRange &horizon,
    const TimedRange &window)
{
    if (timedEmpty(window))
    {
        return {true, false};
    }

    struct Cell
    {
        bool point;
        bool inWindow;
    };

    std::vector<__int128> cuts{horizon.left, horizon.right,
                               window.left, window.right};
    std::sort(cuts.begin(), cuts.end());
    cuts.erase(std::unique(cuts.begin(), cuts.end()), cuts.end());

    std::vector<Cell> cells;
    cells.reserve(7);
    for (std::size_t index = 0; index < cuts.size(); ++index)
    {
        if (timedContains(horizon, cuts[index]))
        {
            cells.push_back({true, timedContains(window, cuts[index])});
        }
        if (index + 1 < cuts.size())
        {
            TimedRange open{cuts[index], cuts[index + 1], false, false};
            if (!timedEmpty(timedIntersection(open, horizon)))
            {
                cells.push_back({false,
                    !timedEmpty(timedIntersection(open, window))});
            }
        }
    }

    constexpr std::uint64_t evenPositions = 0x5555555555555555ULL;
    constexpr std::uint64_t oddPositions = 0xAAAAAAAAAAAAAAAAULL;
    constexpr int maximumWords = static_cast<int>((N + 63) / 64);
    std::array<bool, 2> possible{false, false};

    for (int lhsFirst = 0; lhsFirst < 2; ++lhsFirst)
    {
        const int lhsMax = msb(lhs[lhsFirst]);
        if (lhsMax < 0)
        {
            continue;
        }

        for (int rhsFirst = 0; rhsFirst < 2; ++rhsFirst)
        {
            const int rhsMax = msb(rhs[rhsFirst]);
            if (rhsMax < 0)
            {
                continue;
            }

            const int lhsCount = lhsMax + 1;
            const int rhsCount = rhsMax + 1;
            UntilPackedGrid<N> current(lhsCount, rhsCount);
            UntilPackedGrid<N> next(lhsCount, rhsCount);
            const int words = current.words();
            const int tailBits = rhsCount & 63;
            const std::uint64_t tailMask = tailBits == 0
                ? ~std::uint64_t{0}
                : (std::uint64_t{1} << tailBits) - 1;
            const std::uint64_t rhsOneMask = rhsFirst
                ? evenPositions : oddPositions;
            const std::uint64_t rhsZeroMask = ~rhsOneMask;

            std::array<std::uint64_t, maximumWords> rhsAccepted{};
            for (int index = 0; index < rhsCount; ++index)
            {
                if (rhs[rhsFirst][index])
                {
                    rhsAccepted[index / 64] |=
                        std::uint64_t{1} << (index & 63);
                }
            }

            auto trim = [&](std::uint64_t *bits)
            {
                bits[words - 1] &= tailMask;
            };

            auto orShiftMasked = [&](std::uint64_t *target,
                                     const std::uint64_t *source,
                                     std::uint64_t mask,
                                     int shift)
            {
                if (shift == 0)
                {
                    for (int word = 0; word < words; ++word)
                    {
                        target[word] |= source[word] & mask;
                    }
                }
                else
                {
                    // Descending order also makes the operation safe when
                    // source and target alias.
                    for (int word = words - 1; word >= 0; --word)
                    {
                        const std::uint64_t value = source[word] & mask;
                        target[word] |= value << shift;
                        if (word + 1 < words)
                        {
                            target[word + 1] |= value >> (64 - shift);
                        }
                    }
                }
                trim(target);
            };

            auto orSuffixAfter = [&](std::uint64_t *target,
                                     const std::uint64_t *source1,
                                     const std::uint64_t *source2)
            {
                int firstWord = -1;
                std::uint64_t firstValue = 0;
                for (int word = 0; word < words; ++word)
                {
                    firstValue = source1[word] |
                                 (source2 == nullptr ? 0 : source2[word]);
                    if (firstValue != 0)
                    {
                        firstWord = word;
                        break;
                    }
                }
                if (firstWord < 0)
                {
                    return;
                }

                const int firstBit = untilTrailingZeroBits(firstValue);
                if (firstBit != 63)
                {
                    target[firstWord] |=
                        ~std::uint64_t{0} << (firstBit + 1);
                }
                for (int word = firstWord + 1; word < words; ++word)
                {
                    target[word] = ~std::uint64_t{0};
                }
                trim(target);
            };

            auto horizontalClosure = [&](int lhsIndex, bool lhsValue,
                                         bool inWindow)
            {
                std::uint64_t *state0 = current.bits(lhsIndex, 0);
                std::uint64_t *state1 = current.bits(lhsIndex, 1);
                std::uint64_t *state2 = current.bits(lhsIndex, 2);
                std::uint64_t *state3 = current.bits(lhsIndex, 3);

                if (!lhsValue)
                {
                    orSuffixAfter(state0, state0, state2);
                    orSuffixAfter(state1, state1, state3);
                }
                else if (!inWindow)
                {
                    orSuffixAfter(state0, state0, nullptr);
                    orSuffixAfter(state1, state1, nullptr);
                    orSuffixAfter(state2, state2, nullptr);
                    orSuffixAfter(state3, state3, nullptr);
                }
                else
                {
                    orSuffixAfter(state0, state0, nullptr);
                    orSuffixAfter(state1, state1, nullptr);
                    orSuffixAfter(state3, state3, state2);
                    orShiftMasked(state2, state2, rhsZeroMask, 1);
                }
            };

            auto evaluatedFlags = [](int flags, bool lhsValue,
                                     bool rhsValue, bool point,
                                     bool inWindow)
            {
                bool lhsGood = (flags & 2) != 0;
                bool witness = (flags & 1) != 0;
                if (point)
                {
                    witness = witness ||
                        (inWindow && lhsGood && rhsValue);
                }
                else
                {
                    witness = witness ||
                        (inWindow && lhsGood && lhsValue && rhsValue);
                }
                lhsGood = lhsGood && lhsValue;
                return 2 * static_cast<int>(lhsGood) + witness;
            };

            auto pointFlags = [](int flags, bool lhsPoint,
                                 bool rhsPoint, bool inWindow)
            {
                bool lhsGood = (flags & 2) != 0;
                bool witness = (flags & 1) != 0;
                witness = witness ||
                    (inWindow && lhsGood && rhsPoint);
                lhsGood = lhsGood && lhsPoint;
                return 2 * static_cast<int>(lhsGood) + witness;
            };

            auto hasIntersection = [&](const std::uint64_t *bits)
            {
                for (int word = 0; word < words; ++word)
                {
                    if ((bits[word] & rhsAccepted[word]) != 0)
                    {
                        return true;
                    }
                }
                return false;
            };

            current.bits(0, 2)[0] = 1;

            for (std::size_t cellIndex = 0;
                 cellIndex < cells.size(); ++cellIndex)
            {
                const Cell cell = cells[cellIndex];
                const bool lastCell = cellIndex + 1 == cells.size();

                for (int lhsIndex = 0; lhsIndex < lhsCount; ++lhsIndex)
                {
                    const bool lhsValue = lhsFirst ^ (lhsIndex & 1);
                    if (!cell.point)
                    {
                        horizontalClosure(lhsIndex, lhsValue, cell.inWindow);
                    }

                    std::array<std::array<std::uint64_t, maximumWords>, 4>
                        evaluated{};
                    for (int flags = 0; flags < 4; ++flags)
                    {
                        const std::uint64_t *source =
                            current.bits(lhsIndex, flags);
                        for (int rhsValue = 0; rhsValue < 2; ++rhsValue)
                        {
                            const int targetFlags = evaluatedFlags(
                                flags, lhsValue, rhsValue,
                                cell.point, cell.inWindow);
                            const std::uint64_t mask = rhsValue
                                ? rhsOneMask : rhsZeroMask;
                            for (int word = 0; word < words; ++word)
                            {
                                evaluated[targetFlags][word] |=
                                    source[word] & mask;
                            }
                        }
                    }

                    if (lastCell)
                    {
                        if (lhs[lhsFirst][lhsIndex])
                        {
                            possible[0] = possible[0] ||
                                hasIntersection(evaluated[0].data()) ||
                                hasIntersection(evaluated[2].data());
                            possible[1] = possible[1] ||
                                hasIntersection(evaluated[1].data()) ||
                                hasIntersection(evaluated[3].data());
                        }
                    }
                    else
                    {
                        for (int flags = 0; flags < 4; ++flags)
                        {
                            for (int lhsAdvance = 0;
                                 lhsAdvance <= 1; ++lhsAdvance)
                            {
                                if (lhsIndex + lhsAdvance >= lhsCount)
                                {
                                    continue;
                                }
                                for (int rhsAdvance = 0;
                                     rhsAdvance <= 1; ++rhsAdvance)
                                {
                                    orShiftMasked(
                                        next.bits(lhsIndex + lhsAdvance, flags),
                                        evaluated[flags].data(),
                                        ~std::uint64_t{0}, rhsAdvance);
                                }
                            }
                        }
                    }

                    if (!cell.point)
                    {
                        for (int lhsAdvance = 1;
                             lhsAdvance <= 2; ++lhsAdvance)
                        {
                            if (lhsIndex + lhsAdvance >= lhsCount)
                            {
                                continue;
                            }
                            const int lhsChoices =
                                lhsAdvance == 1 ? 2 : 1;
                            for (int rhsAdvance = 0;
                                 rhsAdvance <= 2; ++rhsAdvance)
                            {
                                const int rhsChoices =
                                    rhsAdvance == 1 ? 2 : 1;
                                for (int lhsChoice = 0;
                                     lhsChoice < lhsChoices; ++lhsChoice)
                                {
                                    const bool lhsPoint = lhsValue ^
                                        (lhsAdvance == 2 ||
                                         (lhsAdvance == 1 && lhsChoice == 1));
                                    for (int rhsChoice = 0;
                                         rhsChoice < rhsChoices; ++rhsChoice)
                                    {
                                        for (int flags = 0;
                                             flags < 4; ++flags)
                                        {
                                            for (int rhsValue = 0;
                                                 rhsValue < 2; ++rhsValue)
                                            {
                                                const bool rhsPoint = rhsValue ^
                                                    (rhsAdvance == 2 ||
                                                     (rhsAdvance == 1 &&
                                                      rhsChoice == 1));
                                                const int targetFlags =
                                                    pointFlags(
                                                        flags, lhsPoint,
                                                        rhsPoint,
                                                        cell.inWindow);
                                                orShiftMasked(
                                                    current.bits(
                                                        lhsIndex + lhsAdvance,
                                                        targetFlags),
                                                    evaluated[flags].data(),
                                                    rhsValue ? rhsOneMask
                                                             : rhsZeroMask,
                                                    rhsAdvance);
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }

                if (possible[0] && possible[1])
                {
                    return possible;
                }
                if (!lastCell)
                {
                    current.swap(next);
                    next.clear();
                }
            }
        }
    }

    return possible;
}
