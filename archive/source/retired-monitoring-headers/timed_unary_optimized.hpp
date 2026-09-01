// Optimized timed Always/Eventually/Historically/Once kernel. This file is
// included from monitoring.hpp after bitsetTimedUnaryLegacy.

struct TimedTruthSummary
{
    bool allZero = true;
    bool allOne = true;
    bool hasZero = false;
    bool hasOne = false;
};

inline void combineTimedTruthSummary(
    TimedTruthSummary &left,
    const TimedTruthSummary &right)
{
    left.allZero = left.allZero && right.allZero;
    left.allOne = left.allOne && right.allOne;
    left.hasZero = left.hasZero || right.hasZero;
    left.hasOne = left.hasOne || right.hasOne;
}

template <std::size_t N>
TimedTruthSummary summarizeTimedTruth(
    const vector<bitset<N>> &language)
{
    const size_t falseSingleton = language[0][0] ? 1 : 0;
    const size_t trueSingleton = language[1][0] ? 1 : 0;
    return {
        language[0][0],
        language[1][0],
        language[0].any() || language[1].count() > trueSingleton,
        language[1].any() || language[0].count() > falseSingleton};
}

template <std::size_t N>
class TimedTruthSummaryIndex
{
public:
    TimedTruthSummaryIndex(
        const vector<vector<bitset<N>>> &values,
        const vector<long long> &segmentation)
        : values_(values), segmentation_(segmentation)
    {
        for (auto &prefix : prefixes_)
        {
            prefix.assign(values.size() + 1, 0);
        }
        for (size_t i = 0; i < values.size(); i++)
        {
            const TimedTruthSummary summary = summarizeTimedTruth(values[i]);
            prefixes_[0][i + 1] = prefixes_[0][i] + summary.allZero;
            prefixes_[1][i + 1] = prefixes_[1][i] + summary.allOne;
            prefixes_[2][i + 1] = prefixes_[2][i] + summary.hasZero;
            prefixes_[3][i + 1] = prefixes_[3][i] + summary.hasOne;
        }
    }

    TimedTruthSummary query(const TimedRange &window) const
    {
        const auto [first, last] = timedIntersectingSegments(
            segmentation_, window);
        TimedTruthSummary result;
        if (first == last)
        {
            // A geometrically nonempty window can still lie before the first
            // supplied segment when retained history starts after time zero.
            // Match the empty materialized profile: neither truth value is
            // represented. The default summary remains the combine identity.
            return {false, false, false, false};
        }

        if (first + 1 == last)
        {
            combineTimedTruthSummary(
                result, restrictedSummary(first, window));
            return result;
        }

        combineTimedTruthSummary(
            result, restrictedSummary(first, window));
        if (first + 1 < last - 1)
        {
            combineTimedTruthSummary(
                result, fullRange(first + 1, last - 1));
        }
        combineTimedTruthSummary(
            result, restrictedSummary(last - 1, window));
        return result;
    }

private:
    TimedTruthSummary fullRange(size_t first, size_t last) const
    {
        const int length = static_cast<int>(last - first);
        const auto count = [&](size_t property)
        {
            return prefixes_[property][last] -
                   prefixes_[property][first];
        };
        return {
            count(0) == length,
            count(1) == length,
            count(2) != 0,
            count(3) != 0};
    }

    TimedTruthSummary restrictedSummary(
        size_t segmentIndex,
        const TimedRange &window) const
    {
        const TimedRange segment{
            2 * (__int128)(segmentation_[segmentIndex]),
            2 * (__int128)(segmentation_[segmentIndex + 1]),
            true,
            false};
        const TimedRange part = timedIntersection(segment, window);
        return summarizeTimedTruth(timedRestrictedSegment(
            values_[segmentIndex], segment, part));
    }

    const vector<vector<bitset<N>>> &values_;
    const vector<long long> &segmentation_;
    array<vector<int>, 4> prefixes_;
};

template <std::size_t N>
vector<vector<bitset<N>>> bitsetTimedUnary(
    const vector<vector<bitset<N>>> &values,
    const vector<long long> &segmentation,
    const long long &a,
    const long long &b,
    const bool &upperInfinite,
    const bool &leftClosed,
    const bool &rightClosed,
    const bool &always,
    const bool &past,
    int s,
    int e)
{
    static_assert(N > 0, "timed languages require positive bitset capacity");

    validateTimedArguments(
        values.size(), nullopt, segmentation,
        a, b, upperInfinite, s, e);
    if (e == -1)
    {
        e = static_cast<int>(values.size());
    }
    for (const auto &segment : values)
    {
        if (segment.size() != 2)
        {
            throw invalid_argument(
                "timed languages require exactly two buckets per segment");
        }
        // The profile algebra assumes the valid monitor invariant that every
        // factor language is nonempty. Retain exact legacy behavior for test
        // or external callers that violate that invariant.
        if (segment[0].none() && segment[1].none())
        {
            return bitsetTimedUnaryLegacy(
                values, segmentation, a, b, upperInfinite,
                leftClosed, rightClosed, always, past, s, e);
        }
    }

    vector<vector<bitset<N>>> output(
        values.size(), vector<bitset<N>>(2));
    vector<long long> offsets{a};
    if (!upperInfinite)
    {
        offsets.push_back(b);
    }
    sort(offsets.begin(), offsets.end());
    offsets.erase(unique(offsets.begin(), offsets.end()), offsets.end());

    const TimedDirection direction = past
        ? TimedDirection::Past
        : TimedDirection::Future;
    const vector<vector<long long>> criticalPoints =
        timedCriticalPoints(segmentation, offsets, true, direction);
    const TimedTruthSummaryIndex<N> summaries(values, segmentation);
    const TimedChangeIndex<N> changeIndex(values, segmentation);
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
            TimedRange window;
            if (past)
            {
                window = upperInfinite
                    ? TimedRange{
                        domain.left,
                        t - 2 * (__int128)(a),
                        true,
                        leftClosed}
                    : TimedRange{
                        t - 2 * (__int128)(b),
                        t - 2 * (__int128)(a),
                        rightClosed,
                        leftClosed};
            }
            else
            {
                window = {
                    t + 2 * (__int128)(a),
                    upperInfinite
                        ? domain.right
                        : t + 2 * (__int128)(b),
                    leftClosed,
                    upperInfinite ? false : rightClosed};
            }
            window = timedIntersection(window, domain);

            array<bool, 2> bits{false, false};
            if (timedEmpty(window))
            {
                bits[always ? 1 : 0] = true;
            }
            else
            {
                TimedTruthSummary summary;
                if (window.left == window.right)
                {
                    summary = summarizeTimedTruth(segmentFirstBit(
                        timedProfile(values, segmentation, window)));
                }
                else
                {
                    summary = summaries.query(window);
                }
                bits = always
                    ? array<bool, 2>{summary.hasZero, summary.allOne}
                    : array<bool, 2>{summary.allZero, summary.hasOne};
            }

            vector<bitset<N>> region(2);
            if (bits[0] && bits[1])
            {
                __int128 bound = 1;
                if (!point)
                {
                    __int128 sum = 0;
                    for (const long long offset : offsets)
                    {
                        sum += changeIndex.changes(
                            placement.first, placement.second, offset,
                            direction);
                    }
                    bound += 2 * sum;
                }
                if (bound > static_cast<__int128>(N))
                {
                    throw overflow_error(always
                        ? "timed always exceeds the fixed bitset size"
                        : "timed eventually exceeds the fixed bitset size");
                }
                for (size_t length = 0;
                     length < static_cast<size_t>(bound);
                     length++)
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
