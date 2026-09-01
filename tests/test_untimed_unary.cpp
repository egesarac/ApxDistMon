#include <algorithm>
#include <bitset>
#include <cassert>
#include <iostream>
#include <random>
#include <stdexcept>
#include <vector>

#include "monitoring.hpp"

namespace
{

using Segment = std::vector<std::bitset<SIZE>>;
using Language = std::vector<Segment>;

Language referenceAlways(const Language &input)
{
    Language result(input.size(), Segment(2));
    bool firstBit0 = false;
    bool firstBit1 = true;

    for (int i = static_cast<int>(input.size()) - 1; i >= 0; --i)
    {
        if (firstBit1)
        {
            std::vector<int> lengths(4);
            lengths[0] = msb(input[i][0] & evenMask) + 1;
            lengths[1] = msb(input[i][0] & oddMask) + 1;
            lengths[2] = msb(input[i][1] & oddMask) + 1;
            lengths[3] = msb(input[i][1] & evenMask) + 1;

            result[i][0][0] = std::max(lengths[0], lengths[2]) > 0;
            result[i][0][1] = std::max(lengths[1], lengths[3] - 2) > 0;
            result[i][1][0] = input[i][1][0];
        }

        if (firstBit0)
        {
            result[i][0][0] = true;
        }

        firstBit0 = result[i][0][0] || result[i][0][1];
        firstBit1 = result[i][1][0];
    }

    return result;
}

Language referenceEventually(const Language &input)
{
    Language result(input.size(), Segment(2));
    bool firstBit0 = true;
    bool firstBit1 = false;

    for (int i = static_cast<int>(input.size()) - 1; i >= 0; --i)
    {
        if (firstBit0)
        {
            std::vector<int> lengths(4);
            lengths[0] = msb(input[i][0] & evenMask) + 1;
            lengths[1] = msb(input[i][0] & oddMask) + 1;
            lengths[2] = msb(input[i][1] & oddMask) + 1;
            lengths[3] = msb(input[i][1] & evenMask) + 1;

            result[i][0][0] = input[i][0][0];
            result[i][1][0] = std::max(lengths[1], lengths[3]) > 0;
            result[i][1][1] = std::max(lengths[0] - 2, lengths[2]) > 0;
        }

        if (firstBit1)
        {
            result[i][1][0] = true;
        }

        firstBit0 = result[i][0][0];
        firstBit1 = result[i][1][0] || result[i][1][1];
    }

    return result;
}

Segment segmentFromCode(unsigned int code)
{
    Segment segment(2);
    for (std::size_t bucket = 0; bucket < 2; ++bucket)
    {
        for (std::size_t bit = 0; bit < 3; ++bit)
        {
            segment[bucket][bit] = (code & (1U << (bucket * 3 + bit))) != 0;
        }
    }
    return segment;
}

void assertMatchesReference(const Language &input)
{
    assert(bitsetAlways(input) == referenceAlways(input));
    assert(bitsetEventually(input) == referenceEventually(input));
}

void testExhaustiveTwoSegmentLanguages()
{
    Language input(2, Segment(2));
    for (unsigned int first = 0; first < 64; ++first)
    {
        input[0] = segmentFromCode(first);
        for (unsigned int second = 0; second < 64; ++second)
        {
            input[1] = segmentFromCode(second);
            assertMatchesReference(input);
        }
    }
}

void testRandomFullWidthLanguages()
{
    std::mt19937 generator(0x51A17EEDU);
    std::uniform_int_distribution<int> segmentCount(1, 64);
    std::uniform_int_distribution<int> bucket(0, 1);
    std::uniform_int_distribution<int> bit(0, SIZE - 1);

    for (int trial = 0; trial < 2000; ++trial)
    {
        Language input(segmentCount(generator), Segment(2));
        const int draws = trial % 2 == 0 ? 4 : 80;

        for (Segment &segment : input)
        {
            for (int draw = 0; draw < draws; ++draw)
            {
                segment[bucket(generator)][bit(generator)] = true;
            }

            if (trial % 4 == 0 || trial % 4 == 1)
            {
                segment[1][0] = true;
            }
            if (trial % 4 == 0 || trial % 4 == 2)
            {
                segment[0][0] = true;
            }
        }

        assertMatchesReference(input);
    }
}

template <typename Exception, typename Action>
void expectException(Action action)
{
    bool thrown = false;
    try
    {
        action();
    }
    catch (const Exception &)
    {
        thrown = true;
    }
    assert(thrown);
}

Segment constantSegment(bool value)
{
    Segment segment(2);
    segment[value ? 1 : 0][0] = true;
    return segment;
}

void testUnboundedFutureWrappers()
{
    const std::vector<long long> segmentation{0, 1, 2, 3};

    const Language eventualInput{
        constantSegment(false),
        constantSegment(false),
        constantSegment(true)};
    const Language eventualInputBefore = eventualInput;
    const std::vector<long long> segmentationBefore = segmentation;
    const Language eventualFull = bitsetEventually(eventualInput);
    assert(bitsetUnboundedEventually(
        eventualInput, segmentation, 0, true) == eventualFull);

    Language eventualExpected = eventualFull;
    eventualExpected[0][0].reset();
    eventualExpected[0][1].reset();
    eventualExpected[2][0].reset();
    eventualExpected[2][1].reset();
    const Language eventualPartial = bitsetUnboundedEventually(
        eventualInput, segmentation, 0, true, 1, 2);
    assert(eventualPartial == eventualExpected);
    const Language eventualClippedInput{eventualInput[1]};
    assert(eventualPartial[1] != bitsetEventually(eventualClippedInput)[0]);

    const Language alwaysInput{
        constantSegment(true),
        constantSegment(true),
        constantSegment(false)};
    const Language alwaysInputBefore = alwaysInput;
    const Language alwaysFull = bitsetAlways(alwaysInput);
    assert(bitsetUnboundedAlways(
        alwaysInput, segmentation, 0, true) == alwaysFull);

    Language alwaysExpected = alwaysFull;
    alwaysExpected[0][0].reset();
    alwaysExpected[0][1].reset();
    alwaysExpected[2][0].reset();
    alwaysExpected[2][1].reset();
    const Language alwaysPartial = bitsetUnboundedAlways(
        alwaysInput, segmentation, 0, true, 1, 2);
    assert(alwaysPartial == alwaysExpected);
    const Language alwaysClippedInput{alwaysInput[1]};
    assert(alwaysPartial[1] != bitsetAlways(alwaysClippedInput)[0]);

    const Language eventualTimed = bitsetUnboundedEventually(
        eventualInput, segmentation, 1, true, 0, 2);
    assert(eventualTimed == bitsetTimedUnary(
        eventualInput, segmentation, 1, 0, true,
        true, false, false, false, 0, 2));
    const Language alwaysTimed = bitsetUnboundedAlways(
        alwaysInput, segmentation, 1, true, 0, 2);
    assert(alwaysTimed == bitsetTimedUnary(
        alwaysInput, segmentation, 1, 0, true,
        true, false, true, false, 0, 2));

    const std::vector<long long> mismatchedSegmentation{0, 1, 2};
    expectException<std::invalid_argument>([&]()
    {
        bitsetUnboundedEventually(
            eventualInput, mismatchedSegmentation, 0, true);
    });
    expectException<std::invalid_argument>([&]()
    {
        bitsetUnboundedAlways(alwaysInput, mismatchedSegmentation, 0, true);
    });

    const Language malformed{Segment(1)};
    const std::vector<long long> oneSegment{0, 1};
    expectException<std::invalid_argument>([&]()
    {
        bitsetUnboundedEventually(malformed, oneSegment, 0, true);
    });
    expectException<std::invalid_argument>([&]()
    {
        bitsetUnboundedAlways(malformed, oneSegment, 0, true);
    });

    expectException<std::out_of_range>([&]()
    {
        bitsetUnboundedEventually(
            eventualInput, segmentation, 0, true, -1, 2);
    });
    expectException<std::out_of_range>([&]()
    {
        bitsetUnboundedEventually(
            eventualInput, segmentation, 0, true, 0, 4);
    });
    expectException<std::out_of_range>([&]()
    {
        bitsetUnboundedAlways(alwaysInput, segmentation, 0, true, 2, 1);
    });
    expectException<std::out_of_range>([&]()
    {
        bitsetUnboundedAlways(alwaysInput, segmentation, 0, true, 0, 4);
    });

    assert(eventualInput == eventualInputBefore);
    assert(alwaysInput == alwaysInputBefore);
    assert(segmentation == segmentationBefore);
}

} // namespace

int main()
{
    evenMask = generateBitmask(0);
    oddMask = generateBitmask(1);

    // The optimized operators must not depend on driver initialization of this
    // unrelated global mask.
    allExceptFirstMask.reset();

    assertMatchesReference(Language{});
    testExhaustiveTwoSegmentLanguages();
    testRandomFullWidthLanguages();
    testUnboundedFutureWrappers();
    std::cout << "untimed unary tests passed\n";
    return 0;
}
