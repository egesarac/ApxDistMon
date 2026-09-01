#include <array>
#include <cassert>
#include <iostream>
#include <random>
#include <stdexcept>
#include <vector>

#include "monitoring.hpp"

using BoolLanguage = vector<bitset<SIZE>>;
using SegmentedLanguage = vector<BoolLanguage>;

void initializeMasks()
{
    evenMask = generateBitmask(0);
    oddMask = generateBitmask(1);
    allExceptFirstMask.set();
    allExceptFirstMask[0] = false;
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

BoolLanguage randomFactor(
    mt19937 &generator,
    int maximumIndex,
    bool allowEmpty = true)
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

SegmentedLanguage randomLanguage(
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

BoolLanguage constantFactor(bool value)
{
    BoolLanguage factor(2);
    factor[value][0] = true;
    return factor;
}

bool emptyFactor(const BoolLanguage &factor)
{
    return factor[0].none() && factor[1].none();
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

void testUnaryExhaustive()
{
    constexpr int length = 3;
    constexpr unsigned int languageCount = 1U << (2 * length);
    for (unsigned int mask = 0; mask < languageCount; mask++)
    {
        const SegmentedLanguage input{factorFromMask(mask, length)};
        for (const bool incomingFalse : {false, true})
        {
            for (const bool incomingTrue : {false, true})
            {
                assert(bitsetAlwaysPast(
                           input, 0, 1, incomingFalse, incomingTrue) ==
                       bitsetAlwaysPastLegacy(
                           input, 0, 1, incomingFalse, incomingTrue));
                assert(bitsetEventuallyPast(
                           input, 0, 1, incomingFalse, incomingTrue) ==
                       bitsetEventuallyPastLegacy(
                           input, 0, 1, incomingFalse, incomingTrue));
            }
        }
    }
}

void testSinceExhaustive()
{
    constexpr int length = 3;
    constexpr unsigned int languageCount = 1U << (2 * length);
    for (unsigned int lhsMask = 0; lhsMask < languageCount; lhsMask++)
    {
        const SegmentedLanguage lhs{factorFromMask(lhsMask, length)};
        for (unsigned int rhsMask = 0; rhsMask < languageCount; rhsMask++)
        {
            const SegmentedLanguage rhs{factorFromMask(rhsMask, length)};
            for (const bool incomingFalse : {false, true})
            {
                for (const bool incomingTrue : {false, true})
                {
                    assert(bitsetSinceNonStrict(
                               lhs, rhs, 0, 1,
                               incomingFalse, incomingTrue) ==
                           bitsetSinceNonStrictLegacy(
                               lhs, rhs, 0, 1,
                               incomingFalse, incomingTrue));
                }
            }
        }
    }
}

void testSinceTwoSegmentExhaustive()
{
    constexpr int length = 2;
    constexpr unsigned int languageCount = 1U << (2 * length);
    for (unsigned int lhsFirst = 0;
         lhsFirst < languageCount; lhsFirst++)
    {
        for (unsigned int lhsSecond = 0;
             lhsSecond < languageCount; lhsSecond++)
        {
            const SegmentedLanguage lhs{
                factorFromMask(lhsFirst, length),
                factorFromMask(lhsSecond, length)};
            for (unsigned int rhsFirst = 0;
                 rhsFirst < languageCount; rhsFirst++)
            {
                for (unsigned int rhsSecond = 0;
                     rhsSecond < languageCount; rhsSecond++)
                {
                    const SegmentedLanguage rhs{
                        factorFromMask(rhsFirst, length),
                        factorFromMask(rhsSecond, length)};
                    for (const bool incomingFalse : {false, true})
                    {
                        for (const bool incomingTrue : {false, true})
                        {
                            assert(bitsetSinceNonStrict(
                                       lhs, rhs, 0, 2,
                                       incomingFalse, incomingTrue) ==
                                   bitsetSinceNonStrictLegacy(
                                       lhs, rhs, 0, 2,
                                       incomingFalse, incomingTrue));
                        }
                    }
                }
            }
        }
    }
}

void testAllCarriesAndRanges()
{
    mt19937 generator(110011);
    for (int test = 0; test < 240; test++)
    {
        const int segments = 1 + generator() % 6;
        const SegmentedLanguage lhs =
            randomLanguage(generator, segments, 7);
        const SegmentedLanguage rhs =
            randomLanguage(generator, segments, 7);
        const SegmentedLanguage lhsBefore = lhs;
        const SegmentedLanguage rhsBefore = rhs;

        for (int s = 0; s <= segments; s++)
        {
            for (int e = s; e <= segments; e++)
            {
                for (const bool incomingFalse : {false, true})
                {
                    for (const bool incomingTrue : {false, true})
                    {
                        assert(bitsetAlwaysPast(
                                   lhs, s, e,
                                   incomingFalse, incomingTrue) ==
                               bitsetAlwaysPastLegacy(
                                   lhs, s, e,
                                   incomingFalse, incomingTrue));
                        assert(bitsetEventuallyPast(
                                   lhs, s, e,
                                   incomingFalse, incomingTrue) ==
                               bitsetEventuallyPastLegacy(
                                   lhs, s, e,
                                   incomingFalse, incomingTrue));
                        assert(bitsetSinceNonStrict(
                                   lhs, rhs, s, e,
                                   incomingFalse, incomingTrue) ==
                               bitsetSinceNonStrictLegacy(
                                   lhs, rhs, s, e,
                                   incomingFalse, incomingTrue));
                    }
                }
            }
        }

        assert(bitsetAlwaysPast(lhs) == bitsetAlwaysPastLegacy(lhs));
        assert(bitsetEventuallyPast(lhs) ==
               bitsetEventuallyPastLegacy(lhs));
        assert(bitsetSinceNonStrict(lhs, rhs) ==
               bitsetSinceNonStrictLegacy(lhs, rhs));
        assert(lhs == lhsBefore);
        assert(rhs == rhsBefore);
    }
}

void testReversal()
{
    constexpr int length = 5;
    constexpr unsigned int languageCount = 1U << (2 * length);
    for (unsigned int mask = 0; mask < languageCount; mask++)
    {
        const BoolLanguage language = factorFromMask(mask, length);
        const auto reversed = reverseAlternatingSegment(
            language[0], language[1]);
        const auto restored = reverseAlternatingSegment(
            reversed[0], reversed[1]);
        assert(restored[0] == language[0]);
        assert(restored[1] == language[1]);

        const AlternatingEndpointSummary originalSummary =
            summarizeAlternatingEndpoints(language[0], language[1]);
        const AlternatingEndpointSummary mappedSummary =
            reverseAlternatingSummary(originalSummary);
        const AlternatingEndpointSummary reversedSummary =
            summarizeAlternatingEndpoints(reversed[0], reversed[1]);
        assert(mappedSummary.singletonFalse ==
               reversedSummary.singletonFalse);
        assert(mappedSummary.singletonTrue ==
               reversedSummary.singletonTrue);
        assert(mappedSummary.nonSingletonLengths ==
               reversedSummary.nonSingletonLengths);
        assert(mappedSummary.nonSingletonLengths[0] ==
               originalSummary.nonSingletonLengths[0]);
        assert(mappedSummary.nonSingletonLengths[1] ==
               originalSummary.nonSingletonLengths[2]);
        assert(mappedSummary.nonSingletonLengths[2] ==
               originalSummary.nonSingletonLengths[1]);
        assert(mappedSummary.nonSingletonLengths[3] ==
               originalSummary.nonSingletonLengths[3]);
    }
}

void testChunkContinuation()
{
    mt19937 generator(220022);
    for (int test = 0; test < 500; test++)
    {
        const int segments = 2 + generator() % 7;
        const int split = 1 + generator() % (segments - 1);
        const SegmentedLanguage lhs =
            randomLanguage(generator, segments, 8);
        const SegmentedLanguage rhs =
            randomLanguage(generator, segments, 8);

        for (const bool initiallyFalse : {false, true})
        {
            for (const bool initiallyTrue : {false, true})
            {
                const SegmentedLanguage alwaysWhole = bitsetAlwaysPast(
                    lhs, 0, segments, initiallyFalse, initiallyTrue);
                SegmentedLanguage alwaysChunked = bitsetAlwaysPast(
                    lhs, 0, split, initiallyFalse, initiallyTrue);
                const bool alwaysFalse =
                    alwaysChunked[split - 1][0][0] ||
                    alwaysChunked[split - 1][1][1];
                const bool alwaysTrue =
                    alwaysChunked[split - 1][1][0];
                const SegmentedLanguage alwaysSuffix = bitsetAlwaysPast(
                    lhs, split, segments, alwaysFalse, alwaysTrue);
                for (int i = split; i < segments; i++)
                {
                    alwaysChunked[i] = alwaysSuffix[i];
                }
                assert(alwaysChunked == alwaysWhole);

                const SegmentedLanguage eventuallyWhole =
                    bitsetEventuallyPast(
                        lhs, 0, segments,
                        initiallyFalse, initiallyTrue);
                SegmentedLanguage eventuallyChunked =
                    bitsetEventuallyPast(
                        lhs, 0, split,
                        initiallyFalse, initiallyTrue);
                const bool eventuallyFalse =
                    eventuallyChunked[split - 1][0][0];
                const bool eventuallyTrue =
                    eventuallyChunked[split - 1][1][0] ||
                    eventuallyChunked[split - 1][0][1];
                const SegmentedLanguage eventuallySuffix =
                    bitsetEventuallyPast(
                        lhs, split, segments,
                        eventuallyFalse, eventuallyTrue);
                for (int i = split; i < segments; i++)
                {
                    eventuallyChunked[i] = eventuallySuffix[i];
                }
                assert(eventuallyChunked == eventuallyWhole);

                const SegmentedLanguage sinceWhole = bitsetSinceNonStrict(
                    lhs, rhs, 0, segments,
                    initiallyFalse, initiallyTrue);
                SegmentedLanguage sinceChunked = bitsetSinceNonStrict(
                    lhs, rhs, 0, split,
                    initiallyFalse, initiallyTrue);
                const auto [sinceLastFalse, sinceLastTrue] =
                    bitsetLastBits(sinceChunked[split - 1]);
                const auto [lhsLastFalse, lhsLastTrue] =
                    bitsetLastBits(lhs[split - 1]);
                const SegmentedLanguage sinceSuffix = bitsetSinceNonStrict(
                    lhs, rhs, split, segments,
                    sinceLastFalse || lhsLastFalse,
                    sinceLastTrue && lhsLastTrue);
                for (int i = split; i < segments; i++)
                {
                    sinceChunked[i] = sinceSuffix[i];
                }
                assert(sinceChunked == sinceWhole);
            }
        }
    }
}

void assertOnlyRangeIsPopulated(
    const SegmentedLanguage &language,
    int s,
    int e)
{
    for (int segment = 0;
         segment < static_cast<int>(language.size()); segment++)
    {
        if (segment < s || segment >= e)
        {
            assert(emptyFactor(language[segment]));
        }
    }
}

void testUnboundedHistoryPrefix()
{
    const vector<long long> segmentation{0, 1, 2, 3};

    const SegmentedLanguage onceInput{
        constantFactor(true),
        constantFactor(false),
        constantFactor(false)};
    const SegmentedLanguage onceFull = bitsetEventuallyPast(onceInput);
    const SegmentedLanguage oncePartial = bitsetUnboundedEventuallyPast(
        onceInput, segmentation, 0, true, 2, 3);
    const SegmentedLanguage onceClipped = bitsetEventuallyPast(
        onceInput, 2, 3);
    assert(oncePartial[2] == onceFull[2]);
    assert(oncePartial[2] != onceClipped[2]);
    assertOnlyRangeIsPopulated(oncePartial, 2, 3);

    const SegmentedLanguage historicallyInput{
        constantFactor(false),
        constantFactor(true),
        constantFactor(true)};
    const SegmentedLanguage historicallyFull =
        bitsetAlwaysPast(historicallyInput);
    const SegmentedLanguage historicallyPartial =
        bitsetUnboundedAlwaysPast(
            historicallyInput, segmentation, 0, true, 2, 3);
    const SegmentedLanguage historicallyClipped = bitsetAlwaysPast(
        historicallyInput, 2, 3);
    assert(historicallyPartial[2] == historicallyFull[2]);
    assert(historicallyPartial[2] != historicallyClipped[2]);
    assertOnlyRangeIsPopulated(historicallyPartial, 2, 3);

    const SegmentedLanguage sinceLhs{
        constantFactor(true),
        constantFactor(true),
        constantFactor(true)};
    const SegmentedLanguage sinceRhs{
        constantFactor(true),
        constantFactor(false),
        constantFactor(false)};
    const SegmentedLanguage sinceFull =
        bitsetSinceNonStrict(sinceLhs, sinceRhs);
    const SegmentedLanguage sincePartial = bitsetUnboundedSince(
        sinceLhs, sinceRhs, segmentation, 0, true, 2, 3);
    const SegmentedLanguage sinceClipped = bitsetSinceNonStrict(
        sinceLhs, sinceRhs, 2, 3);
    assert(sincePartial[2] == sinceFull[2]);
    assert(sincePartial[2] != sinceClipped[2]);
    assertOnlyRangeIsPopulated(sincePartial, 2, 3);
}

void testCapacityAndValidation()
{
    vector<vector<bitset<1>>> constants(
        3, vector<bitset<1>>(2));
    constants[0][0][0] = true;
    constants[1][1][0] = true;
    constants[2][1][0] = true;
    const auto once = bitsetEventuallyPast(constants);
    const auto historically = bitsetAlwaysPast(constants);
    const auto since = bitsetSinceNonStrict(constants, constants);
    assert(once.size() == constants.size());
    assert(historically.size() == constants.size());
    assert(since.size() == constants.size());

    static_assert(SIZE >= 5,
                  "past Since capacity fixtures require at least five bits");
    BoolLanguage untilBoundaryLhs(2);
    BoolLanguage untilBoundaryRhs(2);
    BoolLanguage untilOverflowLhs(2);
    BoolLanguage untilOverflowRhs(2);
    if (SIZE % 2 == 0)
    {
        untilBoundaryLhs[0][1] = true;
        untilBoundaryRhs[0][SIZE - 2] = true;
        untilOverflowLhs[1][2] = true;
        untilOverflowRhs[1][SIZE - 1] = true;
    }
    else
    {
        untilBoundaryLhs[1][2] = true;
        untilBoundaryRhs[1][SIZE - 2] = true;
        untilOverflowLhs[0][1] = true;
        untilOverflowRhs[0][SIZE - 1] = true;
    }
    const auto boundaryLhsArray = reverseAlternatingSegment(
        untilBoundaryLhs[0], untilBoundaryLhs[1]);
    const auto boundaryRhsArray = reverseAlternatingSegment(
        untilBoundaryRhs[0], untilBoundaryRhs[1]);
    const auto overflowLhsArray = reverseAlternatingSegment(
        untilOverflowLhs[0], untilOverflowLhs[1]);
    const auto overflowRhsArray = reverseAlternatingSegment(
        untilOverflowRhs[0], untilOverflowRhs[1]);
    const SegmentedLanguage boundaryLhs{
        BoolLanguage{boundaryLhsArray[0], boundaryLhsArray[1]}};
    const SegmentedLanguage boundaryRhs{
        BoolLanguage{boundaryRhsArray[0], boundaryRhsArray[1]}};
    const SegmentedLanguage overflowLhs{
        BoolLanguage{overflowLhsArray[0], overflowLhsArray[1]}};
    const SegmentedLanguage overflowRhs{
        BoolLanguage{overflowRhsArray[0], overflowRhsArray[1]}};
    const SegmentedLanguage boundary = bitsetSinceNonStrict(
        boundaryLhs, boundaryRhs, 0, 1, false, true);
    assert(boundary[0][0][SIZE - 1] ||
           boundary[0][1][SIZE - 1]);
    expectException<overflow_error>([&]()
    {
        bitsetSinceNonStrict(
            overflowLhs, overflowRhs, 0, 1, false, true);
    });

    const SegmentedLanguage valid{constantFactor(true)};
    const SegmentedLanguage noSegments;
    const SegmentedLanguage malformed{BoolLanguage(1)};
    expectException<invalid_argument>([&]()
    {
        bitsetAlwaysPast(malformed);
    });
    expectException<invalid_argument>([&]()
    {
        bitsetEventuallyPast(malformed);
    });
    expectException<invalid_argument>([&]()
    {
        bitsetSinceNonStrict(valid, noSegments);
    });
    expectException<invalid_argument>([&]()
    {
        bitsetSinceNonStrict(malformed, valid);
    });
    expectException<out_of_range>([&]()
    {
        bitsetAlwaysPast(valid, -1, 1);
    });
    expectException<out_of_range>([&]()
    {
        bitsetEventuallyPast(valid, 1, 0);
    });
    expectException<out_of_range>([&]()
    {
        bitsetSinceNonStrict(valid, valid, 0, 2);
    });

    const vector<long long> malformedSegmentation{0, 0};
    expectException<invalid_argument>([&]()
    {
        bitsetUnboundedEventuallyPast(
            valid, malformedSegmentation, 1, true);
    });
    expectException<invalid_argument>([&]()
    {
        bitsetUnboundedSince(
            valid, valid, malformedSegmentation, 1, true);
    });
}

int main()
{
    initializeMasks();
    testUnaryExhaustive();
    testSinceExhaustive();
    testSinceTwoSegmentExhaustive();
    testAllCarriesAndRanges();
    testReversal();
    testChunkContinuation();
    testUnboundedHistoryPrefix();
    testCapacityAndValidation();
    cout << "optimized untimed past differential tests passed" << endl;
    return 0;
}
