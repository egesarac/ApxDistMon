// Optimized timed-operator kernels. This file is included from monitoring.hpp
// after the legacy differential-oracle implementations.

enum class TimedDirection
{
    Future,
    Past
};

inline TimedRange reverseTimedRange(const TimedRange &range)
{
    return {
        -range.right,
        -range.left,
        range.rightClosed,
        range.leftClosed};
}

template <std::size_t N>
vector<bitset<N>> bitsetConcat(
    const vector<bitset<N>> &left,
    const vector<bitset<N>> &right)
{
    static_assert(N > 0, "timed languages require positive bitset capacity");
    if (left[0].none() && left[1].none())
    {
        return right;
    }

    vector<bitset<N>> output(2);
    const array<int, 2> rightMaximum{msb(right[0]), msb(right[1])};
    constexpr int scalarSpanLimit = 8;

    for (int first = 0; first < 2; first++)
    {
        const int leftMaximum = msb(left[first]);
        for (int leftIndex = 0; leftIndex <= leftMaximum; leftIndex++)
        {
            if (!left[first][leftIndex])
            {
                continue;
            }

            const int last = first ^ (leftIndex % 2);
            for (int second = 0; second < 2; second++)
            {
                const int maximum = rightMaximum[second];
                if (maximum < 0)
                {
                    continue;
                }

                const int shift = leftIndex + (last != second);
                if (maximum + shift >= static_cast<int>(N))
                {
                    throw overflow_error(
                        "value expression exceeds the fixed bitset size");
                }

                // A word shift is substantially cheaper for a broad right
                // bucket. Preserve the scalar path only for short active
                // prefixes where scanning at most eight indices wins.
                if (maximum + 1 > scalarSpanLimit)
                {
                    output[first] |= right[second] << shift;
                }
                else
                {
                    for (int rightIndex = 0; rightIndex <= maximum;
                         rightIndex++)
                    {
                        if (right[second][rightIndex])
                        {
                            output[first][shift + rightIndex] = true;
                        }
                    }
                }
            }
        }
    }

    return output;
}

inline pair<size_t, size_t> timedIntersectingSegments(
    const vector<long long> &segmentation,
    const TimedRange &range)
{
    const size_t segmentCount = segmentation.empty()
        ? 0
        : segmentation.size() - 1;
    if (segmentCount == 0 || timedEmpty(range))
    {
        return {0, 0};
    }

    // First segment whose open right endpoint lies strictly after range.left.
    size_t low = 1;
    size_t high = segmentation.size();
    while (low < high)
    {
        const size_t middle = low + (high - low) / 2;
        if (2 * (__int128)(segmentation[middle]) <= range.left)
        {
            low = middle + 1;
        }
        else
        {
            high = middle;
        }
    }
    if (low == segmentation.size())
    {
        return {segmentCount, segmentCount};
    }
    const size_t first = low - 1;

    // A closed range endpoint includes a segment beginning at that point; an
    // open endpoint does not.
    low = 0;
    high = segmentCount;
    while (low < high)
    {
        const size_t middle = low + (high - low) / 2;
        const __int128 point = 2 * (__int128)(segmentation[middle]);
        const bool included = range.rightClosed
            ? point <= range.right
            : point < range.right;
        if (included)
        {
            low = middle + 1;
        }
        else
        {
            high = middle;
        }
    }
    const size_t last = max(first, low);
    return {first, min(last, segmentCount)};
}

template <std::size_t N>
vector<bitset<N>> timedRestrictedSegment(
    const vector<bitset<N>> &value,
    const TimedRange &segment,
    const TimedRange &part)
{
    if (part.left == segment.left && part.right == segment.right &&
        part.leftClosed == segment.leftClosed &&
        part.rightClosed == segment.rightClosed)
    {
        return value;
    }
    if (part.left == part.right && part.left == segment.left)
    {
        return segmentFirstBit(value);
    }
    if (part.left == segment.left)
    {
        vector<bitset<N>> restricted = part.right == segment.right
            ? value
            : segmentPrefix(value);
        if (!part.leftClosed)
        {
            restricted = segmentWithoutFirstPoint(restricted);
        }
        return restricted;
    }
    if (part.right == segment.right)
    {
        return segmentSuffix(value);
    }
    return segmentInfix(value);
}

template <std::size_t N>
vector<bitset<N>> timedProfile(
    const vector<vector<bitset<N>>> &values,
    const vector<long long> &segmentation,
    const TimedRange &horizon)
{
    vector<bitset<N>> profile(2);
    const auto [first, last] = timedIntersectingSegments(
        segmentation, horizon);
    for (size_t i = first; i < last; i++)
    {
        const TimedRange segment{
            2 * (__int128)(segmentation[i]),
            2 * (__int128)(segmentation[i + 1]), true, false};
        const TimedRange part = timedIntersection(segment, horizon);
        if (timedEmpty(part))
        {
            continue;
        }
        profile = bitsetConcat(
            profile,
            timedRestrictedSegment(values[i], segment, part));
    }
    return profile;
}

// Whole-range products are associative under the monitor invariant that every
// factor language is nonempty. bitsetTimedUntil checks that invariant and
// falls back to the legacy fold before constructing this index otherwise.
template <std::size_t N>
class TimedProfileIndex
{
public:
    TimedProfileIndex(
        const vector<vector<bitset<N>>> &values,
        const vector<long long> &segmentation)
        : segmentation_(segmentation)
    {
        while (base_ < values.size())
        {
            base_ *= 2;
        }
        products_.resize(2 * base_);
        for (size_t i = 0; i < values.size(); i++)
        {
            products_[base_ + i] = valueProduct(values[i]);
        }
        for (size_t i = base_ - 1; i > 0; i--)
        {
            products_[i] = combine(
                products_[2 * i], products_[2 * i + 1]);
        }
    }

    vector<bitset<N>> query(const TimedRange &horizon) const
    {
        vector<bitset<N>> profile(2);
        const auto [first, last] = timedIntersectingSegments(
            segmentation_, horizon);
        if (first == last)
        {
            return profile;
        }

        profile = bitsetConcat(
            profile, restricted(first, horizon));
        if (first + 1 < last - 1)
        {
            append(profile, wholeRange(first + 1, last - 1));
        }
        if (first + 1 < last)
        {
            profile = bitsetConcat(
                profile, restricted(last - 1, horizon));
        }
        return profile;
    }

private:
    enum class ProductState
    {
        identity,
        value,
        overflow
    };

    struct Product
    {
        ProductState state = ProductState::identity;
        vector<bitset<N>> language;
    };

    static Product valueProduct(const vector<bitset<N>> &language)
    {
        return {ProductState::value, language};
    }

    static Product combine(const Product &left, const Product &right)
    {
        if (left.state == ProductState::overflow ||
            right.state == ProductState::overflow)
        {
            return {ProductState::overflow, {}};
        }
        if (left.state == ProductState::identity)
        {
            return right;
        }
        if (right.state == ProductState::identity)
        {
            return left;
        }
        try
        {
            return {
                ProductState::value,
                bitsetConcat(left.language, right.language)};
        }
        catch (const overflow_error &)
        {
            // Eager cache construction must not make an unused wide range
            // fail. Propagate the overflow only if a query selects it.
            return {ProductState::overflow, {}};
        }
    }

    static void append(
        vector<bitset<N>> &profile,
        const Product &product)
    {
        if (product.state == ProductState::overflow)
        {
            throw overflow_error(
                "value expression exceeds the fixed bitset size");
        }
        if (product.state == ProductState::value)
        {
            profile = bitsetConcat(profile, product.language);
        }
    }

    Product wholeRange(size_t first, size_t last) const
    {
        Product left;
        Product right;
        first += base_;
        last += base_;
        while (first < last)
        {
            if (first & 1U)
            {
                left = combine(left, products_[first++]);
            }
            if (last & 1U)
            {
                right = combine(products_[--last], right);
            }
            first /= 2;
            last /= 2;
        }
        return combine(left, right);
    }

    vector<bitset<N>> restricted(
        size_t segmentIndex,
        const TimedRange &horizon) const
    {
        const TimedRange segment{
            2 * (__int128)(segmentation_[segmentIndex]),
            2 * (__int128)(segmentation_[segmentIndex + 1]),
            true,
            false};
        return timedRestrictedSegment(
            products_[base_ + segmentIndex].language,
            segment,
            timedIntersection(segment, horizon));
    }

    const vector<long long> &segmentation_;
    size_t base_ = 1;
    vector<Product> products_;
};

inline vector<vector<long long>> timedCriticalPoints(
    const vector<long long> &segmentation,
    const vector<long long> &offsets,
    bool includeSegmentStarts,
    TimedDirection direction)
{
    if (segmentation.size() < 2)
    {
        return {};
    }

    const size_t segmentCount = segmentation.size() - 1;
    vector<vector<long long>> buckets(segmentCount);
    if (includeSegmentStarts)
    {
        for (size_t i = 0; i < segmentCount; i++)
        {
            buckets[i].push_back(segmentation[i]);
        }
    }

    // For a fixed offset the shifted endpoints remain sorted, so a single
    // monotone segment pointer buckets all events in O(S). Across k distinct
    // offsets this is O(kS), rather than one binary search per event.
    for (const long long offset : offsets)
    {
        size_t segment = 0;
        for (const long long endpoint : segmentation)
        {
            const __int128 shifted = direction == TimedDirection::Future
                ? (__int128)(endpoint) - (__int128)(offset)
                : (__int128)(endpoint) + (__int128)(offset);
            if (shifted < segmentation.front())
            {
                continue;
            }
            if (shifted >= segmentation.back())
            {
                break;
            }
            while (segment + 1 < segmentation.size() &&
                   shifted >= segmentation[segment + 1])
            {
                segment++;
            }
            const long long point = static_cast<long long>(shifted);
            buckets[segment].push_back(point);
        }
    }

    for (auto &bucket : buckets)
    {
        sort(bucket.begin(), bucket.end());
        bucket.erase(unique(bucket.begin(), bucket.end()), bucket.end());
    }
    return buckets;
}

inline vector<vector<long long>> timedFutureCriticalPoints(
    const vector<long long> &segmentation,
    const vector<long long> &offsets,
    bool includeSegmentStarts)
{
    return timedCriticalPoints(
        segmentation, offsets, includeSegmentStarts,
        TimedDirection::Future);
}

template <std::size_t N>
class TimedChangeIndex
{
public:
    TimedChangeIndex(
        const vector<vector<bitset<N>>> &values,
        const vector<long long> &segmentation)
        : segmentation_(segmentation), prefix_(values.size() + 1, 0)
    {
        for (size_t i = 0; i < values.size(); i++)
        {
            prefix_[i + 1] = prefix_[i] +
                max(msb(values[i][0]), msb(values[i][1]));
        }
    }

    long long changes(
        long long left,
        long long right,
        long long offset,
        TimedDirection direction) const
    {
        const __int128 lower = direction == TimedDirection::Future
            ? (__int128)(left) + (__int128)(offset)
            : (__int128)(left) - (__int128)(offset);
        const __int128 upper = direction == TimedDirection::Future
            ? (__int128)(right) + (__int128)(offset)
            : (__int128)(right) - (__int128)(offset);
        if (lower >= upper)
        {
            return 0;
        }

        size_t firstEndpoint = 1;
        size_t high = segmentation_.size() - 1;
        while (firstEndpoint < high)
        {
            const size_t middle = firstEndpoint +
                (high - firstEndpoint) / 2;
            if ((__int128)(segmentation_[middle]) <= lower)
            {
                firstEndpoint = middle + 1;
            }
            else
            {
                high = middle;
            }
        }

        size_t lastEndpoint = 1;
        high = segmentation_.size() - 1;
        while (lastEndpoint < high)
        {
            const size_t middle = lastEndpoint +
                (high - lastEndpoint) / 2;
            if ((__int128)(segmentation_[middle]) < upper)
            {
                lastEndpoint = middle + 1;
            }
            else
            {
                high = middle;
            }
        }
        const long long endpointCount = lastEndpoint > firstEndpoint
            ? static_cast<long long>(lastEndpoint - firstEndpoint)
            : 0;

        const TimedRange range{
            2 * lower, 2 * upper, false, false};
        const auto [firstSegment, lastSegment] =
            timedIntersectingSegments(segmentation_, range);
        return endpointCount +
            prefix_[lastSegment] - prefix_[firstSegment];
    }

    long long changes(
        long long left,
        long long right,
        long long offset) const
    {
        return changes(
            left, right, offset, TimedDirection::Future);
    }

private:
    const vector<long long> &segmentation_;
    vector<long long> prefix_;
};

inline void validateTimedArguments(
    size_t firstSize,
    optional<size_t> secondSize,
    const vector<long long> &segmentation,
    long long a,
    long long b,
    bool upperInfinite,
    int s,
    int e)
{
    if (secondSize && *secondSize != firstSize)
    {
        throw invalid_argument("timed operands must have equal segment counts");
    }
    if (segmentation.size() != firstSize + 1 || segmentation.empty())
    {
        throw invalid_argument(
            "timed segmentation must contain one more endpoint than segments");
    }
    for (size_t i = 1; i < segmentation.size(); i++)
    {
        if (segmentation[i - 1] >= segmentation[i])
        {
            throw invalid_argument(
                "timed segmentation endpoints must be strictly increasing");
        }
    }
    if (a < 0 || (!upperInfinite && b < a))
    {
        throw invalid_argument("invalid timed interval bounds");
    }
    if (e == -1)
    {
        e = static_cast<int>(firstSize);
    }
    if (s < 0 || e < s || static_cast<size_t>(e) > firstSize)
    {
        throw out_of_range("invalid timed segment range");
    }
}

inline void validateTimedFutureArguments(
    size_t firstSize,
    optional<size_t> secondSize,
    const vector<long long> &segmentation,
    long long a,
    long long b,
    bool upperInfinite,
    int s,
    int e)
{
    validateTimedArguments(
        firstSize, secondSize, segmentation,
        a, b, upperInfinite, s, e);
}

#include "possible_until_bitparallel.hpp"

template <std::size_t N>
array<bool, 2> possibleUntilBits(
    const vector<bitset<N>> &lhs,
    const vector<bitset<N>> &rhs,
    const TimedRange &horizon,
    const TimedRange &window)
{
    return possibleUntilBitsBitParallel(lhs, rhs, horizon, window);
}

template <std::size_t N>
vector<vector<bitset<N>>> bitsetTimedUntil(
    const vector<vector<bitset<N>>> &lhs,
    const vector<vector<bitset<N>>> &rhs,
    const vector<long long> &segmentation,
    long long a,
    long long b,
    bool upperInfinite,
    bool leftClosed,
    bool rightClosed,
    int s,
    int e)
{
    validateTimedArguments(
        lhs.size(), rhs.size(), segmentation,
        a, b, upperInfinite, s, e);
    if (e == -1)
    {
        e = static_cast<int>(lhs.size());
    }

    for (size_t i = 0; i < lhs.size(); i++)
    {
        if (lhs[i].size() != 2 || rhs[i].size() != 2)
        {
            throw invalid_argument(
                "timed languages require exactly two buckets per segment");
        }
        if ((lhs[i][0].none() && lhs[i][1].none()) ||
            (rhs[i][0].none() && rhs[i][1].none()))
        {
            return bitsetTimedUntilLegacy(
                lhs, rhs, segmentation, a, b, upperInfinite,
                leftClosed, rightClosed, s, e);
        }
    }

    vector<vector<bitset<N>>> output(
        lhs.size(), vector<bitset<N>>(2));
    vector<long long> offsets{0, a};
    if (!upperInfinite)
    {
        offsets.push_back(b);
    }
    sort(offsets.begin(), offsets.end());
    offsets.erase(unique(offsets.begin(), offsets.end()), offsets.end());

    const vector<vector<long long>> criticalPoints =
        timedCriticalPoints(
            segmentation, offsets, false, TimedDirection::Future);
    const TimedProfileIndex<N> lhsProfiles(lhs, segmentation);
    const TimedProfileIndex<N> rhsProfiles(rhs, segmentation);
    const TimedChangeIndex<N> lhsChanges(lhs, segmentation);
    const TimedChangeIndex<N> rhsChanges(rhs, segmentation);
    const TimedRange domain{
        0, 2 * (__int128)(segmentation.back()), true, false};

    for (int i = s; i < e; i++)
    {
        vector<pair<long long, long long>> placements;
        placements.reserve(2 * criticalPoints[i].size() + 1);
        long long cursor = segmentation[i];
        for (const long long point : criticalPoints[i])
        {
            if (cursor < point)
            {
                placements.push_back({cursor, point});
            }
            placements.push_back({point, point});
            cursor = point;
        }
        if (cursor < segmentation[i + 1])
        {
            placements.push_back({cursor, segmentation[i + 1]});
        }

        vector<bitset<N>> result(2);
        for (const auto &placement : placements)
        {
            const bool point = placement.first == placement.second;
            const __int128 t = point
                ? 2 * (__int128)(placement.first)
                : (__int128)(placement.first) + placement.second;
            TimedRange horizon{
                t,
                upperInfinite
                    ? domain.right
                    : t + 2 * (__int128)(b),
                true,
                !upperInfinite};
            TimedRange window{
                t + 2 * (__int128)(a),
                upperInfinite
                    ? domain.right
                    : t + 2 * (__int128)(b),
                leftClosed,
                upperInfinite ? false : rightClosed};
            horizon = timedIntersection(horizon, domain);
            window = timedIntersection(window, domain);

            const vector<bitset<N>> lhsProfile =
                lhsProfiles.query(horizon);
            const vector<bitset<N>> rhsProfile =
                rhsProfiles.query(horizon);
            const array<bool, 2> bits = possibleUntilBits(
                lhsProfile, rhsProfile, horizon, window);

            vector<bitset<N>> region(2);
            if (bits[0] && bits[1])
            {
                long long bound = 1;
                if (!point)
                {
                    long long sum = 0;
                    for (const long long offset : offsets)
                    {
                        sum += lhsChanges.changes(
                            placement.first, placement.second, offset,
                            TimedDirection::Future);
                        sum += rhsChanges.changes(
                            placement.first, placement.second, offset,
                            TimedDirection::Future);
                    }
                    bound += 2 * sum;
                }
                if (bound > static_cast<long long>(N))
                {
                    throw overflow_error(
                        "timed until exceeds the fixed bitset size");
                }
                for (int length = 0; length < bound; length++)
                {
                    region[0][length] = true;
                    region[1][length] = true;
                }
            }
            else
            {
                region[0][0] = bits[0];
                region[1][0] = bits[1];
            }
            result = bitsetConcat(result, region);
        }
        output[i] = result;
    }
    return output;
}
