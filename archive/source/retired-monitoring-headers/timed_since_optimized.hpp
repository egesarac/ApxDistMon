// Optimized timed non-strict Since kernel. This file is included from
// monitoring.hpp after the frozen scalar Since implementation and the shared
// timed-future infrastructure.

template <std::size_t N>
vector<vector<bitset<N>>> bitsetTimedSince(
    const vector<vector<bitset<N>>> &lhs,
    const vector<vector<bitset<N>>> &rhs,
    const vector<long long> &segmentation,
    const long long &a,
    const long long &b,
    const bool &upperInfinite,
    const bool &leftClosed,
    const bool &rightClosed,
    int s,
    int e)
{
    static_assert(N > 0, "timed since requires positive bitset capacity");
    validateTimedArguments(
        lhs.size(), rhs.size(), segmentation,
        a, b, upperInfinite, s, e);
    if (e == -1)
    {
        e = static_cast<int>(lhs.size());
    }

    for (size_t segment = 0; segment < lhs.size(); segment++)
    {
        if (lhs[segment].size() != 2 || rhs[segment].size() != 2)
        {
            throw invalid_argument(
                "timed languages require exactly two buckets per segment");
        }
        // Cached concatenation products require nonempty factors. Preserve
        // the established behavior for valid callers outside that invariant.
        if ((lhs[segment][0].none() && lhs[segment][1].none()) ||
            (rhs[segment][0].none() && rhs[segment][1].none()))
        {
            return bitsetTimedSinceLegacy(
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

    const vector<vector<long long>> criticalPoints = timedCriticalPoints(
        segmentation, offsets, false, TimedDirection::Past);
    const TimedProfileIndex<N> lhsProfiles(lhs, segmentation);
    const TimedProfileIndex<N> rhsProfiles(rhs, segmentation);
    const TimedChangeIndex<N> lhsChanges(lhs, segmentation);
    const TimedChangeIndex<N> rhsChanges(rhs, segmentation);
    const TimedRange domain{
        0, 2 * (__int128)(segmentation.back()), true, false};

    const auto reverseProfile = [](const vector<bitset<N>> &profile)
    {
        const array<bitset<N>, 2> reversed =
            reverseAlternatingSegment(profile[0], profile[1]);
        return vector<bitset<N>>{reversed[0], reversed[1]};
    };

    for (int segment = s; segment < e; segment++)
    {
        vector<pair<long long, long long>> placements;
        placements.reserve(2 * criticalPoints[segment].size() + 1);
        long long cursor = segmentation[segment];
        for (const long long point : criticalPoints[segment])
        {
            if (cursor < point)
            {
                placements.push_back({cursor, point});
            }
            placements.push_back({point, point});
            cursor = point;
        }
        if (cursor < segmentation[segment + 1])
        {
            placements.push_back({cursor, segmentation[segment + 1]});
        }

        vector<bitset<N>> result(2);
        for (const auto &placement : placements)
        {
            const bool point = placement.first == placement.second;
            const __int128 t = point
                ? 2 * (__int128)(placement.first)
                : (__int128)(placement.first) + placement.second;

            TimedRange horizon{
                upperInfinite
                    ? domain.left
                    : t - 2 * (__int128)(b),
                t,
                true,
                true};
            TimedRange window{
                upperInfinite
                    ? domain.left
                    : t - 2 * (__int128)(b),
                t - 2 * (__int128)(a),
                upperInfinite ? true : rightClosed,
                leftClosed};
            horizon = timedIntersection(horizon, domain);
            window = timedIntersection(window, domain);

            const vector<bitset<N>> lhsProfile = reverseProfile(
                lhsProfiles.query(horizon));
            const vector<bitset<N>> rhsProfile = reverseProfile(
                rhsProfiles.query(horizon));
            const TimedRange reversedHorizon = reverseTimedRange(horizon);
            const TimedRange reversedWindow = reverseTimedRange(window);
            const array<bool, 2> bits = possibleUntilBits(
                lhsProfile, rhsProfile,
                reversedHorizon, reversedWindow);

            vector<bitset<N>> region(2);
            if (bits[0] && bits[1])
            {
                __int128 bound = 1;
                if (!point)
                {
                    __int128 changes = 0;
                    for (const long long offset : offsets)
                    {
                        changes += lhsChanges.changes(
                            placement.first, placement.second, offset,
                            TimedDirection::Past);
                        changes += rhsChanges.changes(
                            placement.first, placement.second, offset,
                            TimedDirection::Past);
                    }
                    bound += 2 * changes;
                }
                if (bound > static_cast<__int128>(N))
                {
                    throw overflow_error(
                        "timed since exceeds the fixed bitset size");
                }
                for (size_t length = 0;
                     length < static_cast<size_t>(bound); length++)
                {
                    region[0].set(length);
                    region[1].set(length);
                }
            }
            else
            {
                region[0][0] = bits[0];
                region[1][0] = bits[1];
            }
            result = bitsetConcat(result, region);
        }
        output[segment] = result;
    }
    return output;
}
