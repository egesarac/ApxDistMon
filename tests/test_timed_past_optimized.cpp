#include <algorithm>
#include <array>
#include <cassert>
#include <iostream>
#include <limits>
#include <optional>
#include <random>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

#include "monitoring.hpp"

using BoolLanguage = vector<bitset<SIZE>>;
using SegmentedLanguage = vector<BoolLanguage>;

struct IntervalSpec
{
    long long a;
    long long b;
    bool upperInfinite;
    bool leftClosed;
    bool rightClosed;
};

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

SegmentedLanguage randomTimedLanguage(
    mt19937 &generator,
    int segments,
    int maximumIndex)
{
    SegmentedLanguage result;
    result.reserve(segments);
    for (int segment = 0; segment < segments; segment++)
    {
        result.push_back(randomFactor(generator, maximumIndex));
    }
    return result;
}

vector<long long> randomSegmentation(
    mt19937 &generator,
    int segments,
    long long first = 0)
{
    vector<long long> result{first};
    for (int segment = 0; segment < segments; segment++)
    {
        result.push_back(result.back() + 1 + generator() % 5);
    }
    return result;
}

BoolLanguage constantFactor(bool value)
{
    BoolLanguage result(2);
    result[value][0] = true;
    return result;
}

SegmentedLanguage constantSignal(const vector<bool> &values)
{
    SegmentedLanguage result;
    result.reserve(values.size());
    for (const bool value : values)
    {
        result.push_back(constantFactor(value));
    }
    return result;
}

BoolLanguage makeLanguage(initializer_list<string> words)
{
    BoolLanguage result(2);
    for (const string &word : words)
    {
        assert(!word.empty());
        for (size_t index = 1; index < word.size(); index++)
        {
            assert(word[index - 1] != word[index]);
        }
        result[word.front() - '0'][word.size() - 1] = true;
    }
    return result;
}

void assertFactorWords(
    const SegmentedLanguage &values,
    int segment,
    const set<string> &expected)
{
    assert(segment >= 0);
    assert(static_cast<size_t>(segment) < values.size());
    assert(bitset2stringset(values[segment]) == expected);
}

template <typename Exception, typename Function>
void expectException(Function function)
{
    bool thrown = false;
    try
    {
        function();
    }
    catch (const Exception &)
    {
        thrown = true;
    }
    assert(thrown);
}

enum class OutcomeKind
{
    Value,
    InvalidArgument,
    OutOfRange,
    Overflow,
    Other
};

struct TimedOutcome
{
    OutcomeKind kind;
    optional<SegmentedLanguage> value;
};

template <typename Function>
TimedOutcome captureOutcome(Function function)
{
    try
    {
        return {OutcomeKind::Value, function()};
    }
    catch (const invalid_argument &)
    {
        return {OutcomeKind::InvalidArgument, nullopt};
    }
    catch (const out_of_range &)
    {
        return {OutcomeKind::OutOfRange, nullopt};
    }
    catch (const overflow_error &)
    {
        return {OutcomeKind::Overflow, nullopt};
    }
    catch (...)
    {
        return {OutcomeKind::Other, nullopt};
    }
}

template <typename Legacy, typename Optimized>
void assertSameOutcome(Legacy legacy, Optimized optimized)
{
    const TimedOutcome expected = captureOutcome(legacy);
    const TimedOutcome actual = captureOutcome(optimized);
    assert(actual.kind == expected.kind);
    if (actual.kind == OutcomeKind::Value)
    {
        assert(actual.value == expected.value);
    }
}

vector<vector<long long>> scalarPastCriticalPoints(
    const vector<long long> &segmentation,
    const vector<long long> &offsets,
    bool includeSegmentStarts)
{
    vector<vector<long long>> result(segmentation.size() - 1);
    for (size_t segment = 0; segment + 1 < segmentation.size(); segment++)
    {
        if (includeSegmentStarts)
        {
            result[segment].push_back(segmentation[segment]);
        }
        for (const long long endpoint : segmentation)
        {
            for (const long long offset : offsets)
            {
                const __int128 point =
                    (__int128)(endpoint) + (__int128)(offset);
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

long long scalarPastChanges(
    const SegmentedLanguage &values,
    const vector<long long> &segmentation,
    long long left,
    long long right,
    long long offset)
{
    const TimedRange range{
        2 * ((__int128)(left) - (__int128)(offset)),
        2 * ((__int128)(right) - (__int128)(offset)),
        false,
        false};
    long long result = 0;
    for (size_t endpoint = 1;
         endpoint + 1 < segmentation.size(); endpoint++)
    {
        result += timedContains(
            range, 2 * (__int128)(segmentation[endpoint]));
    }
    for (size_t segment = 0; segment < values.size(); segment++)
    {
        const TimedRange segmentRange{
            2 * (__int128)(segmentation[segment]),
            2 * (__int128)(segmentation[segment + 1]),
            true,
            false};
        if (!timedEmpty(timedIntersection(range, segmentRange)))
        {
            result += max(
                msb(values[segment][0]),
                msb(values[segment][1]));
        }
    }
    return result;
}

void testPastIndexedGeometry()
{
    mt19937 generator(330033);
    for (int test = 0; test < 20000; test++)
    {
        const int segments = 1 + generator() % 9;
        const long long first = generator() % 5;
        const vector<long long> segmentation =
            randomSegmentation(generator, segments, first);
        const SegmentedLanguage values =
            randomTimedLanguage(generator, segments, 12);
        vector<long long> offsets{
            0,
            static_cast<long long>(generator() % 8),
            static_cast<long long>(generator() % 8)};
        sort(offsets.begin(), offsets.end());
        offsets.erase(unique(offsets.begin(), offsets.end()), offsets.end());
        const bool includeStarts = (generator() & 1U) != 0;

        assert(timedCriticalPoints(
                   segmentation, offsets, includeStarts,
                   TimedDirection::Past) ==
               scalarPastCriticalPoints(
                   segmentation, offsets, includeStarts));

        const TimedChangeIndex<SIZE> index(values, segmentation);
        for (int query = 0; query < 5; query++)
        {
            const long long span = segmentation.back() - first + 9;
            const long long x = first - 4 + generator() % span;
            const long long y = first - 4 + generator() % span;
            const long long left = min(x, y);
            const long long right = max(x, y);
            for (const long long offset : offsets)
            {
                assert(index.changes(
                           left, right, offset,
                           TimedDirection::Past) ==
                       scalarPastChanges(
                           values, segmentation,
                           left, right, offset));
            }
        }
    }

    const long long maximum = numeric_limits<long long>::max();
    const vector<long long> nearLimit{
        maximum - 100,
        maximum - 70,
        maximum - 30,
        maximum - 5};
    const vector<long long> offsets{0, 10, 80, 200};
    assert(timedCriticalPoints(
               nearLimit, offsets, true, TimedDirection::Past) ==
           scalarPastCriticalPoints(nearLimit, offsets, true));
}

void testTimedRangeReversal()
{
    vector<TimedRange> ranges{
        {0, 0, true, true},
        {0, 4, true, false},
        {-7, 3, false, true},
        {5, 5, false, true},
        {9, 3, true, true},
        {(__int128)(numeric_limits<long long>::min()),
         (__int128)(numeric_limits<long long>::max()),
         false,
         true}};
    for (const TimedRange &range : ranges)
    {
        const TimedRange reversed = reverseTimedRange(range);
        const TimedRange restored = reverseTimedRange(reversed);
        assert(restored.left == range.left);
        assert(restored.right == range.right);
        assert(restored.leftClosed == range.leftClosed);
        assert(restored.rightClosed == range.rightClosed);
        for (int point = -12; point <= 12; point++)
        {
            assert(timedContains(range, point) ==
                   timedContains(reversed, -(__int128)(point)));
        }
        assert(timedEmpty(range) == timedEmpty(reversed));
    }
}

vector<IntervalSpec> timedSpecs()
{
    return {
        {0, 0, false, true, true},
        {0, 0, false, false, true},
        {0, 0, false, true, false},
        {0, 0, false, false, false},
        {0, 2, false, true, true},
        {0, 2, false, true, false},
        {0, 2, false, false, true},
        {0, 2, false, false, false},
        {1, 3, false, true, true},
        {1, 3, false, true, false},
        {1, 3, false, false, true},
        {1, 3, false, false, false},
        {0, 0, true, true, false},
        {1, 0, true, true, false},
        {1, 0, true, false, false},
        {7, 0, true, true, false}};
}

void compareTimedUnary(
    const SegmentedLanguage &values,
    const vector<long long> &segmentation,
    const IntervalSpec &spec,
    int s,
    int e)
{
    for (const bool always : {false, true})
    {
        assertSameOutcome(
            [&]()
            {
                return bitsetTimedUnaryLegacy(
                    values, segmentation,
                    spec.a, spec.b, spec.upperInfinite,
                    spec.leftClosed, spec.rightClosed,
                    always, true, s, e);
            },
            [&]()
            {
                return bitsetTimedUnary(
                    values, segmentation,
                    spec.a, spec.b, spec.upperInfinite,
                    spec.leftClosed, spec.rightClosed,
                    always, true, s, e);
            });
    }
}

void compareTimedSince(
    const SegmentedLanguage &lhs,
    const SegmentedLanguage &rhs,
    const vector<long long> &segmentation,
    const IntervalSpec &spec,
    int s,
    int e)
{
    assertSameOutcome(
        [&]()
        {
            return bitsetTimedSinceLegacy(
                lhs, rhs, segmentation,
                spec.a, spec.b, spec.upperInfinite,
                spec.leftClosed, spec.rightClosed, s, e);
        },
        [&]()
        {
            return bitsetTimedSince(
                lhs, rhs, segmentation,
                spec.a, spec.b, spec.upperInfinite,
                spec.leftClosed, spec.rightClosed, s, e);
        });
}

void assertSinceWords(
    const SegmentedLanguage &lhs,
    const SegmentedLanguage &rhs,
    const vector<long long> &segmentation,
    const IntervalSpec &spec,
    int segment,
    const set<string> &expected)
{
    const SegmentedLanguage legacy = bitsetTimedSinceLegacy(
        lhs, rhs, segmentation,
        spec.a, spec.b, spec.upperInfinite,
        spec.leftClosed, spec.rightClosed, segment, segment + 1);
    const SegmentedLanguage optimized = bitsetTimedSince(
        lhs, rhs, segmentation,
        spec.a, spec.b, spec.upperInfinite,
        spec.leftClosed, spec.rightClosed, segment, segment + 1);
    assert(optimized == legacy);
    assertFactorWords(legacy, segment, expected);
}

void assertTrueSinceMatchesPastEventually(
    const SegmentedLanguage &rhs,
    const vector<long long> &segmentation,
    const IntervalSpec &spec,
    int s,
    int e)
{
    const SegmentedLanguage lhs(
        rhs.size(), constantFactor(true));
    const SegmentedLanguage eventuallyLegacy = bitsetTimedUnaryLegacy(
        rhs, segmentation,
        spec.a, spec.b, spec.upperInfinite,
        spec.leftClosed, spec.rightClosed,
        false, true, s, e);
    const SegmentedLanguage eventually = bitsetTimedUnary(
        rhs, segmentation,
        spec.a, spec.b, spec.upperInfinite,
        spec.leftClosed, spec.rightClosed,
        false, true, s, e);
    const SegmentedLanguage sinceLegacy = bitsetTimedSinceLegacy(
        lhs, rhs, segmentation,
        spec.a, spec.b, spec.upperInfinite,
        spec.leftClosed, spec.rightClosed, s, e);
    const SegmentedLanguage since = bitsetTimedSince(
        lhs, rhs, segmentation,
        spec.a, spec.b, spec.upperInfinite,
        spec.leftClosed, spec.rightClosed, s, e);

    assert(eventually == eventuallyLegacy);
    assert(sinceLegacy == eventuallyLegacy);
    assert(since == eventually);
}

void testTimedUnaryDifferential()
{
    const vector<IntervalSpec> specs = timedSpecs();
    mt19937 generator(440044);
    for (int test = 0; test < 6000; test++)
    {
        const int segments = 1 + generator() % 5;
        const vector<long long> segmentation =
            randomSegmentation(generator, segments);
        const SegmentedLanguage values =
            randomTimedLanguage(generator, segments, 6);
        const SegmentedLanguage before = values;
        const int s = generator() % (segments + 1);
        const int e = s + generator() % (segments - s + 1);
        compareTimedUnary(
            values, segmentation, specs[test % specs.size()], s, e);
        if (test % 17 == 0)
        {
            compareTimedUnary(
                values, segmentation, specs[test % specs.size()], 0, -1);
        }
        assert(values == before);
    }

    SegmentedLanguage withEmpty{
        constantFactor(true), BoolLanguage(2), constantFactor(false)};
    const vector<long long> segmentation{0, 2, 5, 9};
    for (const IntervalSpec &spec : specs)
    {
        compareTimedUnary(withEmpty, segmentation, spec, 0, -1);
    }
}

void testTimedSinceDifferential()
{
    const vector<IntervalSpec> specs = timedSpecs();
    mt19937 generator(550055);
    for (int test = 0; test < 12000; test++)
    {
        const int segments = 1 + generator() % 5;
        const vector<long long> segmentation =
            randomSegmentation(generator, segments);
        const SegmentedLanguage lhs =
            randomTimedLanguage(generator, segments, 6);
        const SegmentedLanguage rhs =
            randomTimedLanguage(generator, segments, 6);
        const SegmentedLanguage lhsBefore = lhs;
        const SegmentedLanguage rhsBefore = rhs;
        const int s = generator() % (segments + 1);
        const int e = s + generator() % (segments - s + 1);
        compareTimedSince(
            lhs, rhs, segmentation,
            specs[(3 * test) % specs.size()], s, e);
        if (test % 17 == 0)
        {
            compareTimedSince(
                lhs, rhs, segmentation,
                specs[(3 * test) % specs.size()], 0, -1);
        }
        assert(lhs == lhsBefore);
        assert(rhs == rhsBefore);
    }

    SegmentedLanguage lhsWithEmpty{
        constantFactor(true), BoolLanguage(2), constantFactor(false)};
    SegmentedLanguage rhsWithEmpty{
        constantFactor(false), constantFactor(true), BoolLanguage(2)};
    const vector<long long> segmentation{0, 2, 5, 9};
    for (const IntervalSpec &spec : specs)
    {
        compareTimedSince(
            lhsWithEmpty, rhsWithEmpty, segmentation, spec, 0, -1);
    }
}

void testSinceWindowProfileDomain()
{
    const vector<long long> segmentation{0, 1, 2, 3, 4};
    const SegmentedLanguage lhs =
        constantSignal({true, true, true, true});
    const SegmentedLanguage rhs{
        makeLanguage({"0"}),
        makeLanguage({"0"}),
        makeLanguage({"1"}),
        makeLanguage({"0"})};

    // Every [t-2,t-1] witness window for t in [3,4) has profile 01.
    // The trailing zero belongs to the horizon, but not to the RHS profile.
    assertSinceWords(
        lhs, rhs, segmentation,
        {1, 2, false, true, true}, 3, {"1"});
}

void testSinceWindowOffsetLengthBound()
{
    const vector<long long> segmentation{0, 1, 2, 3, 4};
    const SegmentedLanguage lhs =
        constantSignal({true, true, true, true});
    const BoolLanguage ambiguous = makeLanguage({"0", "1"});
    const SegmentedLanguage rhsWithoutOutsideChanges{
        ambiguous, ambiguous, ambiguous, makeLanguage({"0"})};
    const SegmentedLanguage rhsWithOutsideChanges{
        ambiguous, ambiguous, ambiguous, makeLanguage({"010"})};
    const IntervalSpec spec{1, 2, false, true, true};
    const set<string> expected{"0", "01", "1", "10"};

    // Offset zero is needed for the LHS horizon, but it is not a witness-
    // window offset when a > 0. Changes in the final RHS segment therefore
    // cannot enlarge the result language.
    assertSinceWords(
        lhs, rhsWithoutOutsideChanges,
        segmentation, spec, 3, expected);
    assertSinceWords(
        lhs, rhsWithOutsideChanges,
        segmentation, spec, 3, expected);
}

void testSinceLhsHorizonOffsetLengthBound()
{
    const vector<long long> segmentation{0, 1, 2, 3};
    const SegmentedLanguage lhs{
        makeLanguage({"1"}),
        makeLanguage({"1"}),
        makeLanguage({"10"})};
    const BoolLanguage ambiguous = makeLanguage({"0", "1"});
    const SegmentedLanguage rhs{
        ambiguous, ambiguous, makeLanguage({"0"})};
    const IntervalSpec spec{1, 2, false, true, true};
    const set<string> expected{
        "0", "01", "010", "0101",
        "1", "10", "101", "1010"};

    // For J = [1,2], offset zero belongs to M_J but not C_J. The change in
    // lhs on [2,3) must therefore contribute to the open placement's bound.
    assertSinceWords(
        lhs, rhs, segmentation, spec, 2, expected);
}

void testTrueSinceMatchesPastEventually()
{
    const vector<IntervalSpec> specs = timedSpecs();
    mt19937 generator(660066);
    for (int test = 0; test < 480; test++)
    {
        constexpr int segments = 4;
        const vector<long long> segmentation =
            randomSegmentation(generator, segments);
        // Two-letter input words keep every intermediate language well below
        // the fixed capacity, so equality cannot be hidden by overflow.
        const SegmentedLanguage rhs =
            randomTimedLanguage(generator, segments, 1);

        int s = 0;
        int e = -1;
        if (test % 3 == 1)
        {
            e = 2;
        }
        else if (test % 3 == 2)
        {
            s = 1;
            e = 4;
        }
        assertTrueSinceMatchesPastEventually(
            rhs, segmentation, specs[test % specs.size()], s, e);
    }
}

void testEmptySinceWindowSkipsProfiles()
{
    const vector<long long> segmentation{0, 1, 2, 3, 4};
    BoolLanguage nearCapacity(2);
    nearCapacity[0].set(SIZE - 1);
    SegmentedLanguage lhs(4, nearCapacity);
    const SegmentedLanguage rhs(4, nearCapacity);
    const IntervalSpec outsideDomain{10, 11, false, true, true};

    // The past witness window is empty throughout [3,4). Constructing either
    // horizon profile would concatenate several maximum-length factors and
    // overflow, but an empty window is false without inspecting a profile.
    assertSinceWords(
        lhs, rhs, segmentation,
        outsideDomain, 3, {"0"});

    // An empty input factor selects the optimized kernel's legacy fallback.
    // The same empty-window short circuit must happen on that path as well.
    lhs[0] = BoolLanguage(2);
    assertSinceWords(
        lhs, rhs, segmentation,
        outsideDomain, 3, {"0"});
}

void testEndpointAndClosureCases()
{
    const vector<long long> segmentation{0, 2, 5, 9};
    const vector<SegmentedLanguage> signals{
        constantSignal({false, false, false}),
        constantSignal({true, true, true}),
        constantSignal({false, true, false}),
        constantSignal({true, false, true})};
    const vector<IntervalSpec> specs = timedSpecs();

    for (const SegmentedLanguage &lhs : signals)
    {
        for (const SegmentedLanguage &rhs : signals)
        {
            for (const IntervalSpec &spec : specs)
            {
                compareTimedUnary(lhs, segmentation, spec, 0, -1);
                compareTimedSince(
                    lhs, rhs, segmentation, spec, 0, -1);
                compareTimedUnary(lhs, segmentation, spec, 1, 2);
                compareTimedSince(
                    lhs, rhs, segmentation, spec, 1, 2);
            }
        }
    }

    // Retained history can start after time zero, while the monitor's timed
    // domain convention remains anchored at zero.
    const vector<long long> retainedSegmentation{4, 6, 9, 13};
    compareTimedUnary(
        signals[2], retainedSegmentation,
        {1, 4, false, false, true}, 0, -1);
    compareTimedSince(
        signals[3], signals[2], retainedSegmentation,
        {1, 0, true, true, false}, 1, 3);
}

void testWordBoundaries()
{
    const vector<long long> segmentation{0, 4};
    for (const int index : {63, 64, 127, 128})
    {
        SegmentedLanguage lhs(1, BoolLanguage(2));
        SegmentedLanguage rhs(1, BoolLanguage(2));
        lhs[0][0][0] = true;
        lhs[0][1][index] = true;
        rhs[0][1][0] = true;
        rhs[0][0][index] = true;
        const IntervalSpec punctual{0, 0, false, true, true};
        compareTimedUnary(lhs, segmentation, punctual, 0, 1);
        compareTimedSince(lhs, rhs, segmentation, punctual, 0, 1);
    }
}

void testWrappersAndPartialRanges()
{
    const vector<long long> segmentation{0, 2, 5, 9};
    const SegmentedLanguage lhs = constantSignal({true, true, false});
    const SegmentedLanguage rhs = constantSignal({false, true, false});

    assert(bitsetBoundedEventuallyPast(
               lhs, segmentation, 1, 3, false, true, 1, 3) ==
           bitsetTimedUnary(
               lhs, segmentation, 1, 3, false,
               false, true, false, true, 1, 3));
    assert(bitsetBoundedAlwaysPast(
               lhs, segmentation, 1, 3, true, false, 1, 3) ==
           bitsetTimedUnary(
               lhs, segmentation, 1, 3, false,
               true, false, true, true, 1, 3));
    assert(bitsetBoundedSince(
               lhs, rhs, segmentation, 1, 3,
               false, true, 1, 3) ==
           bitsetTimedSince(
               lhs, rhs, segmentation, 1, 3, false,
               false, true, 1, 3));
    assert(bitsetUnboundedEventuallyPast(
               lhs, segmentation, 1, false, 1, 3) ==
           bitsetTimedUnary(
               lhs, segmentation, 1, 0, true,
               false, false, false, true, 1, 3));
    assert(bitsetUnboundedAlwaysPast(
               lhs, segmentation, 1, true, 1, 3) ==
           bitsetTimedUnary(
               lhs, segmentation, 1, 0, true,
               true, false, true, true, 1, 3));
    assert(bitsetUnboundedSince(
               lhs, rhs, segmentation, 1, true, 1, 3) ==
           bitsetTimedSince(
               lhs, rhs, segmentation, 1, 0, true,
               true, false, 1, 3));
}

void testValidation()
{
    const SegmentedLanguage valid{constantFactor(true)};
    const SegmentedLanguage noSegments;
    const SegmentedLanguage malformed{BoolLanguage(1)};
    const vector<long long> validSegmentation{0, 2};
    const vector<long long> shortSegmentation{0};
    const vector<long long> repeatedSegmentation{0, 0};

    expectException<invalid_argument>([&]()
    {
        bitsetTimedUnary(
            malformed, validSegmentation,
            0, 1, false, true, true,
            false, true, 0, 1);
    });
    expectException<invalid_argument>([&]()
    {
        bitsetTimedUnary(
            valid, shortSegmentation,
            0, 1, false, true, true,
            false, true, 0, 1);
    });
    expectException<invalid_argument>([&]()
    {
        bitsetTimedUnary(
            valid, repeatedSegmentation,
            0, 1, false, true, true,
            false, true, 0, 1);
    });
    expectException<invalid_argument>([&]()
    {
        bitsetTimedUnary(
            valid, validSegmentation,
            -1, 1, false, true, true,
            false, true, 0, 1);
    });
    expectException<invalid_argument>([&]()
    {
        bitsetTimedUnary(
            valid, validSegmentation,
            2, 1, false, true, true,
            false, true, 0, 1);
    });
    expectException<out_of_range>([&]()
    {
        bitsetTimedUnary(
            valid, validSegmentation,
            0, 1, false, true, true,
            false, true, -1, 1);
    });

    expectException<invalid_argument>([&]()
    {
        bitsetTimedSince(
            valid, noSegments, validSegmentation,
            0, 1, false, true, true, 0, 1);
    });
    expectException<invalid_argument>([&]()
    {
        bitsetTimedSince(
            malformed, valid, validSegmentation,
            0, 1, false, true, true, 0, 1);
    });
    expectException<invalid_argument>([&]()
    {
        bitsetTimedSince(
            valid, valid, repeatedSegmentation,
            0, 1, false, true, true, 0, 1);
    });
    expectException<invalid_argument>([&]()
    {
        bitsetTimedSince(
            valid, valid, validSegmentation,
            -1, 1, false, true, true, 0, 1);
    });
    expectException<out_of_range>([&]()
    {
        bitsetTimedSince(
            valid, valid, validSegmentation,
            0, 1, false, true, true, 1, 0);
    });
}

int main()
{
    initializeTimedMasks();
    testPastIndexedGeometry();
    testTimedRangeReversal();
    testTimedUnaryDifferential();
    testTimedSinceDifferential();
    testSinceWindowProfileDomain();
    testSinceWindowOffsetLengthBound();
    testSinceLhsHorizonOffsetLengthBound();
    testTrueSinceMatchesPastEventually();
    testEmptySinceWindowSkipsProfiles();
    testEndpointAndClosureCases();
    testWordBoundaries();
    testWrappersAndPartialRanges();
    testValidation();
    cout << "optimized timed past differential tests passed" << endl;
    return 0;
}
