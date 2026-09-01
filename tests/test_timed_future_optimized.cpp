#include <cassert>
#include <iostream>
#include <random>
#include <vector>

#include "monitoring.hpp"

using BoolLanguage = vector<bitset<SIZE>>;
using SegmentedLanguage = vector<BoolLanguage>;

void initializeTimedMasks()
{
    evenMask = generateBitmask(0);
    oddMask = generateBitmask(1);
    allExceptFirstMask.set();
    allExceptFirstMask[0] = false;
}

BoolLanguage randomFactor(
    mt19937 &generator,
    int maximumIndex,
    bool allowEmpty = false)
{
    BoolLanguage factor(2);
    do
    {
        factor[0].reset();
        factor[1].reset();
        for (int bucket = 0; bucket < 2; bucket++)
        {
            for (int index = 0; index <= maximumIndex; index++)
            {
                factor[bucket][index] = (generator() & 3U) == 0;
            }
        }
    } while (!allowEmpty && factor[0].none() && factor[1].none());
    return factor;
}

BoolLanguage factorFromMask(unsigned int mask, int length)
{
    BoolLanguage factor(2);
    for (int bucket = 0; bucket < 2; bucket++)
    {
        for (int index = 0; index < length; index++)
        {
            factor[bucket][index] =
                (mask & (1U << (bucket * length + index))) != 0;
        }
    }
    return factor;
}

BoolLanguage singletonFactor(int firstBit, int length)
{
    assert(firstBit == 0 || firstBit == 1);
    assert(length >= 1 && length <= SIZE);
    BoolLanguage factor(2);
    factor[firstBit][length - 1] = true;
    return factor;
}

SegmentedLanguage randomTimedLanguage(
    mt19937 &generator,
    int segments,
    int maximumIndex)
{
    SegmentedLanguage language;
    language.reserve(segments);
    for (int i = 0; i < segments; i++)
    {
        language.push_back(randomFactor(generator, maximumIndex));
    }
    return language;
}

vector<long long> randomSegmentation(mt19937 &generator, int segments)
{
    vector<long long> segmentation{0};
    for (int i = 0; i < segments; i++)
    {
        segmentation.push_back(
            segmentation.back() + 1 + generator() % 5);
    }
    return segmentation;
}

long long legacyChanges(
    const SegmentedLanguage &values,
    const vector<long long> &segmentation,
    long long left,
    long long right,
    long long offset)
{
    const TimedRange range{
        2 * ((__int128)(left) + offset),
        2 * ((__int128)(right) + offset),
        false,
        false};
    long long result = 0;
    for (int i = 1; i + 1 < static_cast<int>(segmentation.size()); i++)
    {
        result += timedContains(
            range, 2 * (__int128)(segmentation[i]));
    }
    for (int i = 0; i < static_cast<int>(values.size()); i++)
    {
        const TimedRange segment{
            2 * (__int128)(segmentation[i]),
            2 * (__int128)(segmentation[i + 1]),
            true,
            false};
        if (!timedEmpty(timedIntersection(range, segment)))
        {
            result += max(msb(values[i][0]), msb(values[i][1]));
        }
    }
    return result;
}

vector<vector<long long>> legacyCriticalPoints(
    const vector<long long> &segmentation,
    const vector<long long> &offsets,
    bool includeSegmentStarts)
{
    vector<vector<long long>> result(segmentation.size() - 1);
    for (int segment = 0;
         segment + 1 < static_cast<int>(segmentation.size());
         segment++)
    {
        if (includeSegmentStarts)
        {
            result[segment].push_back(segmentation[segment]);
        }
        for (const long long endpoint : segmentation)
        {
            for (const long long offset : offsets)
            {
                const __int128 point = (__int128)(endpoint) - offset;
                if (point >= segmentation[segment] &&
                    point < segmentation[segment + 1])
                {
                    result[segment].push_back(
                        static_cast<long long>(point));
                }
            }
        }
        sort(result[segment].begin(), result[segment].end());
        result[segment].erase(
            unique(result[segment].begin(), result[segment].end()),
            result[segment].end());
    }
    return result;
}

void testHybridConcatenation()
{
    mt19937 generator(10101);
    for (int test = 0; test < 20000; test++)
    {
        const BoolLanguage left = randomFactor(generator, 80, true);
        const BoolLanguage right = randomFactor(generator, 80, true);
        assert(bitsetConcat(left, right) ==
               bitsetConcatLegacy(left, right));
    }

    BoolLanguage left(2);
    BoolLanguage right(2);
    left[0][SIZE - 1] = true;
    right[0][1] = true;
    bool legacyOverflow = false;
    bool optimizedOverflow = false;
    try
    {
        bitsetConcatLegacy(left, right);
    }
    catch (const overflow_error &)
    {
        legacyOverflow = true;
    }
    try
    {
        bitsetConcat(left, right);
    }
    catch (const overflow_error &)
    {
        optimizedOverflow = true;
    }
    assert(legacyOverflow && optimizedOverflow);
}

void testIndexedGeometry()
{
    mt19937 generator(20202);
    for (int test = 0; test < 20000; test++)
    {
        const int segments = 1 + generator() % 8;
        const vector<long long> segmentation =
            randomSegmentation(generator, segments);
        const long long duration = segmentation.back();
        const __int128 first =
            static_cast<long long>(generator() % (2 * duration + 7)) - 3;
        const __int128 second =
            static_cast<long long>(generator() % (2 * duration + 7)) - 3;
        const TimedRange range{
            min(first, second),
            max(first, second),
            (generator() & 1U) != 0,
            (generator() & 1U) != 0};
        const auto [indexedFirst, indexedLast] =
            timedIntersectingSegments(segmentation, range);

        vector<size_t> expected;
        for (int i = 0; i < segments; i++)
        {
            const TimedRange segment{
                2 * (__int128)(segmentation[i]),
                2 * (__int128)(segmentation[i + 1]),
                true,
                false};
            if (!timedEmpty(timedIntersection(segment, range)))
            {
                expected.push_back(static_cast<size_t>(i));
            }
        }
        vector<size_t> actual;
        for (size_t i = indexedFirst; i < indexedLast; i++)
        {
            actual.push_back(i);
        }
        assert(actual == expected);
    }
}

void testProfilesAndIndexes()
{
    mt19937 generator(30303);
    for (int test = 0; test < 12000; test++)
    {
        const int segments = 1 + generator() % 7;
        const vector<long long> segmentation =
            randomSegmentation(generator, segments);
        const SegmentedLanguage values =
            randomTimedLanguage(generator, segments, 7);
        const long long duration = segmentation.back();
        const __int128 first =
            static_cast<long long>(generator() % (2 * duration + 7)) - 3;
        const __int128 second =
            static_cast<long long>(generator() % (2 * duration + 7)) - 3;
        const TimedRange range{
            min(first, second),
            max(first, second),
            (generator() & 1U) != 0,
            (generator() & 1U) != 0};
        const BoolLanguage legacyProfile =
            timedProfileLegacy(values, segmentation, range);
        assert(timedProfile(values, segmentation, range) == legacyProfile);
        const TimedProfileIndex<SIZE> profiles(values, segmentation);
        assert(profiles.query(range) == legacyProfile);

        vector<long long> offsets{
            0,
            static_cast<long long>(generator() % (duration + 3)),
            static_cast<long long>(generator() % (duration + 3))};
        sort(offsets.begin(), offsets.end());
        offsets.erase(unique(offsets.begin(), offsets.end()), offsets.end());
        const bool includeStarts = (generator() & 1U) != 0;
        assert(timedFutureCriticalPoints(
                   segmentation, offsets, includeStarts) ==
               legacyCriticalPoints(
                   segmentation, offsets, includeStarts));

        const TimedChangeIndex<SIZE> index(values, segmentation);
        const int segment = generator() % segments;
        const long long left = segmentation[segment];
        const long long right = segmentation[segment + 1];
        for (const long long offset : offsets)
        {
            assert(index.changes(left, right, offset) ==
                   legacyChanges(
                       values, segmentation, left, right, offset));
        }
    }
}

void testPackedUntilSolver()
{
    const vector<pair<TimedRange, TimedRange>> exhaustiveRanges{
        {{0, 0, true, true}, {0, 0, true, true}},
        {{0, 0, true, true}, {1, 0, true, true}},
        {{0, 2, true, false}, {0, 2, true, false}},
        {{0, 2, true, false}, {0, 2, false, false}},
        {{0, 2, true, true}, {0, 2, true, true}},
        {{0, 2, true, true}, {0, 2, false, false}},
        {{0, 2, true, true}, {0, 2, false, true}},
        {{0, 2, true, true}, {0, 2, true, false}},
        {{0, 2, true, true}, {1, 0, true, true}},
        {{0, 8, false, false}, {2, 6, false, false}},
        {{0, 8, false, true}, {-2, 3, true, false}},
        {{0, 8, true, false}, {5, 10, false, true}}};
    constexpr int exhaustiveLength = 3;
    constexpr unsigned int subsetCount =
        1U << (2 * exhaustiveLength);
    for (const auto &ranges : exhaustiveRanges)
    {
        for (unsigned int lhsMask = 1; lhsMask < subsetCount; lhsMask++)
        {
            const BoolLanguage lhs =
                factorFromMask(lhsMask, exhaustiveLength);
            for (unsigned int rhsMask = 1;
                 rhsMask < subsetCount; rhsMask++)
            {
                const BoolLanguage rhs =
                    factorFromMask(rhsMask, exhaustiveLength);
                assert(possibleUntilBits(
                           lhs, rhs, ranges.first, ranges.second) ==
                       possibleUntilBitsLegacy(
                           lhs, rhs, ranges.first, ranges.second));
            }
        }
    }

    mt19937 generator(40404);
    for (int test = 0; test < 12000; test++)
    {
        const BoolLanguage lhs = randomFactor(generator, 8);
        const BoolLanguage rhs = randomFactor(generator, 8);
        const __int128 first = generator() % 5;
        const __int128 second = first + generator() % 5;
        const TimedRange horizon{
            first,
            second,
            (generator() & 1U) != 0,
            (generator() & 1U) != 0};
        if (timedEmpty(horizon))
        {
            continue;
        }
        const __int128 windowFirst = first + generator() % 5;
        const __int128 windowSecond = windowFirst + generator() % 5;
        const TimedRange window{
            windowFirst,
            windowSecond,
            (generator() & 1U) != 0,
            (generator() & 1U) != 0};
        assert(possibleUntilBits(lhs, rhs, horizon, window) ==
               possibleUntilBitsLegacy(lhs, rhs, horizon, window));
    }

    const array<int, 4> wordBoundaryMaximum{63, 64, 127, 128};
    for (int test = 0; test < 200; test++)
    {
        const BoolLanguage lhs = randomFactor(
            generator, wordBoundaryMaximum[test % 4]);
        const BoolLanguage rhs = randomFactor(
            generator, wordBoundaryMaximum[(test + 1) % 4]);
        const TimedRange horizon{0, 8, true, true};
        const TimedRange window{
            2, 6, (test & 1) != 0, (test & 2) != 0};
        assert(possibleUntilBits(lhs, rhs, horizon, window) ==
               possibleUntilBitsLegacy(lhs, rhs, horizon, window));
    }

    // Exercise exact accepted lengths on both sides of each packed-word
    // boundary. The random factors above usually contain many shorter words,
    // which can let the search terminate before reaching these states.
    const array<int, 4> wordBoundaryLengths{63, 64, 127, 128};
    for (const int length : wordBoundaryLengths)
    {
        for (int lhsFirst = 0; lhsFirst < 2; lhsFirst++)
        {
            for (int rhsFirst = 0; rhsFirst < 2; rhsFirst++)
            {
                const BoolLanguage lhs = singletonFactor(lhsFirst, length);
                const BoolLanguage rhs = singletonFactor(rhsFirst, length);
                const TimedRange horizon{0, 8, true, true};
                const TimedRange window{
                    2, 6, lhsFirst != 0, rhsFirst != 0};
                assert(possibleUntilBits(lhs, rhs, horizon, window) ==
                       possibleUntilBitsLegacy(lhs, rhs, horizon, window));
            }
        }
    }
}

void testDeferredProfileOverflow()
{
    const int segments = SIZE + 1;
    vector<long long> segmentation(segments + 1);
    SegmentedLanguage values(segments, BoolLanguage(2));
    for (int i = 0; i <= segments; i++)
    {
        segmentation[i] = i;
        if (i < segments)
        {
            values[i][i & 1][0] = true;
        }
    }

    const TimedProfileIndex<SIZE> profiles(values, segmentation);
    const TimedRange firstSegment{0, 2, true, false};
    assert(profiles.query(firstSegment) ==
           timedProfileLegacy(values, segmentation, firstSegment));

    const TimedRange fullDomain{
        0, 2 * (__int128)(segments), true, false};
    bool legacyOverflow = false;
    bool indexedOverflow = false;
    try
    {
        timedProfileLegacy(values, segmentation, fullDomain);
    }
    catch (const overflow_error &)
    {
        legacyOverflow = true;
    }
    try
    {
        profiles.query(fullDomain);
    }
    catch (const overflow_error &)
    {
        indexedOverflow = true;
    }
    assert(legacyOverflow && indexedOverflow);
}

void testCompleteFutureOperators()
{
    mt19937 generator(50505);
    for (int test = 0; test < 12000; test++)
    {
        const int segments = 1 + generator() % 5;
        const vector<long long> segmentation =
            randomSegmentation(generator, segments);
        const SegmentedLanguage lhs =
            randomTimedLanguage(generator, segments, 6);
        const SegmentedLanguage rhs =
            randomTimedLanguage(generator, segments, 6);
        const long long duration = segmentation.back();
        const long long a = generator() % (duration + 3);
        const bool upperInfinite = (generator() & 1U) != 0;
        const long long b = upperInfinite
            ? 0
            : a + generator() % (duration + 3);
        const bool leftClosed = (generator() & 1U) != 0;
        const bool rightClosed = (generator() & 1U) != 0;
        const int s = generator() % segments;
        const int e = s + 1 + generator() % (segments - s);

        for (const bool always : {false, true})
        {
            const auto legacy = bitsetTimedUnaryLegacy(
                lhs, segmentation, a, b, upperInfinite,
                leftClosed, rightClosed, always, false, s, e);
            const auto optimized = bitsetTimedUnary(
                lhs, segmentation, a, b, upperInfinite,
                leftClosed, rightClosed, always, false, s, e);
            assert(optimized == legacy);
        }

        const auto legacyUntil = bitsetTimedUntilLegacy(
            lhs, rhs, segmentation, a, b, upperInfinite,
            leftClosed, rightClosed, s, e);
        const auto optimizedUntil = bitsetTimedUntil(
            lhs, rhs, segmentation, a, b, upperInfinite,
            leftClosed, rightClosed, s, e);
        assert(optimizedUntil == legacyUntil);
    }
}

void testUntilEventuallyMetamorphism()
{
    mt19937 generator(60606);
    for (int test = 0; test < 2000; test++)
    {
        const int segments = 1 + generator() % 5;
        const vector<long long> segmentation =
            randomSegmentation(generator, segments);
        const SegmentedLanguage rhs =
            randomTimedLanguage(generator, segments, 4);
        SegmentedLanguage truth(segments, BoolLanguage(2));
        for (auto &factor : truth)
        {
            factor[1][0] = true;
        }

        const long long duration = segmentation.back();
        long long a = 0;
        long long b = 0;
        bool upperInfinite = false;
        switch (test % 6)
        {
        case 0:
            // Zero-lower-bound windows cover the compatibility case as well
            // as the open-endpoint refinement.
            b = duration;
            break;
        case 1:
            // A positive lower bound makes the witness window a proper
            // subset of the horizon at every unclipped placement.
            a = 1 + generator() % (duration + 2);
            b = a + 1 + generator() % (duration + 2);
            break;
        case 2:
            // Closed punctual windows exercise singleton profile domains.
            a = generator() % (duration + 2);
            b = a;
            break;
        case 3:
            // Degenerate open intervals are empty before clipping.
            a = generator() % (duration + 2);
            b = a;
            break;
        case 4:
            upperInfinite = true;
            break;
        default:
            upperInfinite = true;
            a = 1 + generator() % (duration + 2);
            break;
        }

        const int closureCase = (test / 6) % 4;
        bool leftClosed = (closureCase & 1) != 0;
        bool rightClosed = (closureCase & 2) != 0;
        if (test % 6 == 2)
        {
            leftClosed = true;
            rightClosed = true;
        }
        else if (test % 6 == 3)
        {
            leftClosed = false;
            rightClosed = (closureCase & 2) != 0;
        }

        const int s = generator() % segments;
        const int e = s + 1 + generator() % (segments - s);
        const auto eventuallyLegacy = bitsetTimedUnaryLegacy(
            rhs, segmentation, a, b, upperInfinite,
            leftClosed, rightClosed, false, false, s, e);
        const auto eventually = bitsetTimedUnary(
            rhs, segmentation, a, b, upperInfinite,
            leftClosed, rightClosed, false, false, s, e);
        assert(eventually == eventuallyLegacy);
        assert(bitsetTimedUntilLegacy(
                   truth, rhs, segmentation, a, b, upperInfinite,
                   leftClosed, rightClosed, s, e) == eventuallyLegacy);
        assert(bitsetTimedUntil(
                   truth, rhs, segmentation, a, b, upperInfinite,
                   leftClosed, rightClosed, s, e) == eventually);
    }
}

int main()
{
    initializeTimedMasks();
    testHybridConcatenation();
    testIndexedGeometry();
    testProfilesAndIndexes();
    testPackedUntilSolver();
    testDeferredProfileOverflow();
    testCompleteFutureOperators();
    testUntilEventuallyMetamorphism();
    cout << "optimized timed future differential tests passed" << endl;
    return 0;
}
