#include <algorithm>
#include <array>
#include <bitset>
#include <cassert>
#include <iostream>
#include <limits>
#include <random>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "monitoring_case_studies.hpp"

namespace
{

ScalarSignal scalar(std::initializer_list<double> values)
{
    ScalarSignal signal;
    long long timestamp = 0;
    for (const double value : values)
    {
        signal.push_back({timestamp++, value});
    }
    return signal;
}

PointSignal points(std::initializer_list<Point> values)
{
    PointSignal signal;
    long long timestamp = 0;
    for (const Point &value : values)
    {
        signal.push_back({timestamp++, value});
    }
    return signal;
}

void assertVerdict(const Evaluation &evaluation, int expected)
{
    assert(evaluation.verdict() == expected);
    if (expected == 0)
    {
        assert(evaluation.falsePossible());
        assert(!evaluation.truePossible());
    }
    else if (expected == 1)
    {
        assert(!evaluation.falsePossible());
        assert(evaluation.truePossible());
    }
    else
    {
        assert(expected == 2);
        assert(evaluation.falsePossible());
        assert(evaluation.truePossible());
    }
}

void testApproximateBitsetKernels()
{
    std::bitset<130> startsFalse;
    startsFalse.set(0);
    startsFalse.set(2);
    startsFalse.set(63);
    startsFalse.set(64);
    startsFalse.set(129);
    std::bitset<130> startsTrue;
    startsTrue.set(1);
    startsTrue.set(126);
    startsTrue.set(127);
    startsTrue.set(128);

    const auto summary = summarizeAlternatingEndpoints(
        startsFalse, startsTrue);
    assert(summary.singletonFalse);
    assert(!summary.singletonTrue);
    assert((summary.nonSingletonLengths ==
            std::array<int, 4>{65, 130, 128, 129}));

    const auto negated = summarizeAlternatingEndpoints(
        startsFalse, startsTrue, true);
    assert(!negated.singletonFalse);
    assert(negated.singletonTrue);
    assert((negated.nonSingletonLengths ==
            std::array<int, 4>{129, 128, 130, 65}));
    assert(wordOrientedMsb(startsFalse) == 129);
    assert(wordOrientedMsb(std::bitset<130>{}) == -1);

    using Language = std::vector<std::vector<std::bitset<SIZE>>>;
    auto smallLanguage = [](unsigned int mask)
    {
        Language result(1, std::vector<std::bitset<SIZE>>(2));
        for (int bit = 0; bit < 6; bit++)
        {
            if ((mask & (1U << bit)) != 0)
            {
                result[0][bit / 3].set(bit % 3);
            }
        }
        return result;
    };

    Language trueConstant = smallLanguage(0);
    trueConstant[0][1][0] = true;
    Language falseConstant = smallLanguage(0);
    falseConstant[0][0][0] = true;

    for (unsigned int leftMask = 0; leftMask < 64; leftMask++)
    {
        const Language left = smallLanguage(leftMask);
        assert(bitsetConjunction(left, trueConstant) == left);
        assert(bitsetDisjunction(left, falseConstant) == left);

        for (unsigned int rightMask = 0; rightMask < 64; rightMask++)
        {
            const Language right = smallLanguage(rightMask);
            const Language expectedOr = bitsetNegation(bitsetConjunction(
                bitsetNegation(left), bitsetNegation(right)));
            assert(bitsetDisjunction(left, right) == expectedOr);
            assert(bitsetConjunction(left, right) ==
                   bitsetConjunction(right, left));
            assert(left == smallLanguage(leftMask));
            assert(right == smallLanguage(rightMask));
        }
    }

    Language shortClassOne(1, std::vector<std::bitset<SIZE>>(2));
    shortClassOne[0][0].set(1);
    Language fullClassOne(1, std::vector<std::bitset<SIZE>>(2));
    fullClassOne[0][0].set(SIZE - 1);
    const Language capacityBoundary = bitsetConjunction(
        shortClassOne, fullClassOne);
    assert(capacityBoundary[0][0][SIZE - 1]);

    Language longTrue(1, std::vector<std::bitset<SIZE>>(2));
    longTrue[0][1].set(SIZE - 2);
    bool rejected = false;
    try
    {
        (void)bitsetConjunction(longTrue, longTrue);
    }
    catch (const std::length_error &)
    {
        rejected = true;
    }
    assert(rejected);

    Language longFalse(1, std::vector<std::bitset<SIZE>>(2));
    longFalse[0][0].set(SIZE - 2);
    rejected = false;
    try
    {
        (void)bitsetDisjunction(longFalse, longFalse);
    }
    catch (const std::length_error &)
    {
        rejected = true;
    }
    assert(rejected);

    rejected = false;
    try
    {
        (void)bitsetConjunction(Language{}, smallLanguage(0));
    }
    catch (const std::invalid_argument &)
    {
        rejected = true;
    }
    assert(rejected);

    rejected = false;
    try
    {
        Language malformed(1, std::vector<std::bitset<SIZE>>(1));
        (void)bitsetDisjunction(malformed, malformed);
    }
    catch (const std::invalid_argument &)
    {
        rejected = true;
    }
    assert(rejected);

    auto captureLengthError = [](auto &&operation)
    {
        try
        {
            return std::make_pair(false, operation());
        }
        catch (const std::length_error &)
        {
            return std::make_pair(true, Language{});
        }
    };

    std::mt19937 generator(0x5eed1234U);
    constexpr int randomizedTrials = 1500;
    for (int trial = 0; trial < randomizedTrials; trial++)
    {
        const int segmentCount = 1 + generator() % 3;
        Language left(
            segmentCount, std::vector<std::bitset<SIZE>>(2));
        Language right(
            segmentCount, std::vector<std::bitset<SIZE>>(2));
        const bool dense = trial % 2 != 0;

        auto populate = [&](Language &language)
        {
            for (auto &segment : language)
            {
                for (auto &bucket : segment)
                {
                    if (dense)
                    {
                        for (size_t index = 0; index < SIZE; index++)
                        {
                            if (generator() % 8 == 0)
                            {
                                bucket[index] = true;
                            }
                        }
                    }
                    else
                    {
                        const int draws = generator() % 5;
                        for (int draw = 0; draw < draws; draw++)
                        {
                            bucket[generator() % SIZE] = true;
                        }
                    }
                }
            }
        };
        populate(left);
        populate(right);

        if (trial % 6 == 2 || trial % 6 == 5)
        {
            for (auto &segment : right)
            {
                segment[0].reset();
                segment[1].reset();
                segment[trial % 6 == 2 ? 1 : 0][0] = true;
            }
        }

        const Language leftBefore = left;
        const Language rightBefore = right;

        const auto directOr = captureLengthError([&]
        {
            return bitsetDisjunction(left, right);
        });
        const auto deMorganOr = captureLengthError([&]
        {
            return bitsetNegation(bitsetConjunction(
                bitsetNegation(left), bitsetNegation(right)));
        });
        assert(directOr.first == deMorganOr.first);
        if (!directOr.first)
        {
            assert(directOr.second == deMorganOr.second);
        }

        const auto forwardAnd = captureLengthError([&]
        {
            return bitsetConjunction(left, right);
        });
        const auto reverseAnd = captureLengthError([&]
        {
            return bitsetConjunction(right, left);
        });
        assert(forwardAnd.first == reverseAnd.first);
        if (!forwardAnd.first)
        {
            assert(forwardAnd.second == reverseAnd.second);
        }

        const auto directPhi4 = captureLengthError([&]
        {
            return bitsetAlways(bitsetDisjunction(
                bitsetNegation(left), bitsetEventually(right)));
        });
        const auto legacyPhi4 = captureLengthError([&]
        {
            return bitsetAlways(bitsetNegation(bitsetConjunction(
                left, bitsetAlways(bitsetNegation(right)))));
        });
        assert(directPhi4.first == legacyPhi4.first);
        if (!directPhi4.first)
        {
            assert(directPhi4.second == legacyPhi4.second);
        }

        assert(left == leftBefore);
        assert(right == rightBefore);
    }
}

void testSharedPaperAbstraction()
{
    using Signals = std::vector<std::vector<std::pair<long long, double>>>;

    {
        const Signals signals{{{0, 0.0}, {5, 1.0}}};
        const auto regions = computeUncertaintyIntervals(signals, 2, 10);
        assert(regions[0][0].empty());
        assert((regions[0][1] == std::vector<long long>{3, 7}));
        const auto segmentation = computeCanonicalSegmentation(
            signals, regions, 10);
        assert((segmentation == std::vector<long long>{0, 3, 7, 10}));
        const auto gamma = computeValueExpressions(
            signals, regions, segmentation);
        assert((gamma[0][1] == std::set<std::string>{"0;1"}));
        assert((gamma[0][2] == std::set<std::string>{"1"}));

        const auto coarse = computeValueExpressionsCoarse(
            signals, regions, segmentation);
        assert((coarse[0][2] == std::set<std::string>{"1"}));
    }

    {
        const Signals signals{
            {{0, 0.0}, {1, 1.0}, {2, 2.0}},
            {{0, 0.0}, {3, 1.0}},
        };
        const auto regions = computeUncertaintyIntervals(signals, 2, 4);
        assert((regions[0][1] == std::vector<long long>{0, 3}));
        assert((regions[0][2] == std::vector<long long>{1, 4}));
        const auto segmentation = computeCanonicalSegmentation(
            signals, regions, 4);
        assert((segmentation == std::vector<long long>{0, 1, 3, 4}));
        const auto gamma = computeValueExpressions(
            signals, regions, segmentation);
        assert((gamma[0][1] == std::set<std::string>{
            "0;1", "0;1;2", "1", "1;2"}));
        const auto coarse = computeValueExpressionsCoarse(
            signals, regions, segmentation);
        assert((coarse[0][1] == std::set<std::string>{"0", "2"}));
    }

    {
        const Signals signals{{{0, 0.0}, {5, 1.0}, {6, 0.0}}};
        const auto regions = computeUncertaintyIntervals(signals, 10, 11);
        assert((regions[0][1] == std::vector<long long>{0, 10}));
        assert((regions[0][2] == std::vector<long long>{1, 11}));
        const auto segmentation = computeCanonicalSegmentation(
            signals, regions, 11);
        assert((segmentation == std::vector<long long>{0, 1, 10, 11}));
        const auto gamma = computeValueExpressions(
            signals, regions, segmentation);
        assert((gamma[0][0] == std::set<std::string>{"0", "0;1"}));
        assert((gamma[0][1] == std::set<std::string>{
            "0;1", "0;1;0", "1", "1;0"}));
        assert((gamma[0][2] == std::set<std::string>{"0", "1;0"}));
    }

    {
        const Signals signals{{{0, 4.0}, {1, 4.0}, {2, 5.0}}};
        const auto regions = computeUncertaintyIntervals(signals, 2, 5);
        assert(regions[0][1].empty());
        assert((regions[0][2] == std::vector<long long>{0, 4}));
    }

    {
        const Signals signals{{
            {0, 0.0}, {1, 0.0}, {8, 1.0}, {9, 0.0}}};
        const auto regions = computeUncertaintyIntervals(signals, 10, 10);
        assert(regions[0][0].empty());
        assert(regions[0][1].empty());
        assert((regions[0][2] == std::vector<long long>{0, 9}));
        assert((regions[0][3] == std::vector<long long>{1, 10}));
        const auto segmentation = computeCanonicalSegmentation(
            signals, regions, 10);
        assert((segmentation == std::vector<long long>{0, 1, 9, 10}));
        const auto gamma = computeValueExpressions(
            signals, regions, segmentation);
        assert((gamma[0][0] == std::set<std::string>{"0", "0;1"}));
        assert((gamma[0][1] == std::set<std::string>{
            "0;1", "0;1;0", "1", "1;0"}));
        assert((gamma[0][2] == std::set<std::string>{"0", "1;0"}));
    }

    {
        const Signals signals{{
            {0, 0.0}, {1, 1.0}, {2, 0.0}, {3, 1.0}}};
        const auto regions = computeUncertaintyIntervals(
            signals, std::numeric_limits<long long>::max(), 4);
        assert((regions[0][1] == std::vector<long long>{0, 2}));
        assert((regions[0][2] == std::vector<long long>{1, 3}));
        assert((regions[0][3] == std::vector<long long>{2, 4}));
    }

    {
        const Signals signals{
            {{0, 0.0}, {2, 1.0}, {3, 0.0}},
            {{0, 4.0}},
        };
        const std::set<int> exactSignals{0};
        const auto regions = computeUncertaintyIntervals(
            signals, 1, 5, exactSignals);
        const auto segmentation = computeCanonicalSegmentation(
            signals, regions, 5, exactSignals);
        assert((segmentation == std::vector<long long>{0, 2, 3, 5}));
        const auto gamma = computeRelativeValueExpressions(
            exactSignals, signals, regions, segmentation);
        assert((gamma[0][0] == std::set<std::string>{"0"}));
        assert((gamma[0][1] == std::set<std::string>{"1"}));
        assert((gamma[0][2] == std::set<std::string>{"0"}));
    }
}

void testSyntheticFormulaRegressions()
{
    auto signalFromWord = [](const std::string &word)
    {
        assert(word.size() == 4);
        std::vector<std::pair<long long, double>> signal{
            {0, static_cast<double>(word[0] - '0')}};
        for (std::size_t index = 1; index < word.size(); ++index)
        {
            if (word[index] != word[index - 1])
            {
                signal.push_back({
                    static_cast<long long>(index) * 1000,
                    static_cast<double>(word[index] - '0')});
            }
        }
        return signal;
    };

    auto evaluate = [&](const std::string &formula,
                        const std::string &pWord,
                        const std::string &qWord)
    {
        const std::vector<std::vector<std::pair<long long, double>>> signals{
            signalFromWord(pWord), signalFromWord(qWord)};
        const auto regions = computeUncertaintyIntervals(signals, 2000, 4000);
        const auto segmentation = computeCanonicalSegmentation(
            signals, regions, 4000);
        const auto gamma = computeValueExpressions(
            signals, regions, segmentation);
        const auto propositions = convertSignalsToAtomicPropositions(gamma, 0.0);
        const auto &p = propositions[0];
        const auto &q = propositions[1];

        std::vector<std::vector<std::bitset<SIZE>>> result;
        if (formula == "phi1")
        {
            result = bitsetAlways(bitsetConjunction(p, q));
        }
        else if (formula == "phi2")
        {
            result = bitsetAlways(bitsetDisjunction(p, q));
        }
        else if (formula == "phi3")
        {
            result = bitsetUnboundedUntil(p, q, segmentation, 0, true);
        }
        else if (formula == "phi4")
        {
            result = bitsetAlways(bitsetDisjunction(
                bitsetNegation(p), bitsetEventually(q)));
        }
        else
        {
            const long long upper = formula == "phi5" ? 1000 : 2000;
            result = bitsetAlways(bitsetNegation(bitsetConjunction(
                p,
                bitsetBoundedAlways(
                    bitsetNegation(q),
                    segmentation,
                    0,
                    upper,
                    true,
                    false))));
        }

        const bool falsePossible = result[0][0].any();
        const bool truePossible = result[0][1].any();
        if (falsePossible && truePossible)
        {
            return 2;
        }
        return truePossible ? 1 : 0;
    };

    evenMask = generateBitmask(0);
    oddMask = generateBitmask(1);
    allExceptFirstMask = evenMask | oddMask;
    allExceptFirstMask[0] = false;

    assert(evaluate("phi1", "1010", "1010") == 0);
    assert(evaluate("phi2", "0000", "1010") == 0);
    assert(evaluate("phi3", "0000", "1111") == 1);
    assert(evaluate("phi3", "1111", "0000") == 0);
    assert(evaluate("phi3", "1000", "1010") == 1);
    assert(evaluate("phi4", "1010", "0000") == 0);
    assert(evaluate("phi5", "1010", "0000") == 0);
    assert(evaluate("phi6", "1010", "0000") == 0);
}

void testWaterTankProperty()
{
    const std::vector<Variant> variants{
        Variant::Adm,
        Variant::AdmFine,
        Variant::AdmFineRelative,
        Variant::AdmCoarse,
        Variant::AdmCoarseRelative,
    };

    for (const Variant variant : variants)
    {
        assertVerdict(evaluateWaterTanks(
            {scalar({5.0}), scalar({5.0})}, 1, 2, variant), 0);
        assertVerdict(evaluateWaterTanks(
            {scalar({6.0}), scalar({5.0})}, 1, 2, variant), 1);
    }

    const ScalarSignals exactBoundary{
        scalar({3.3333336}),
        scalar({3.3333336}),
        scalar({3.3333328}),
    };
    assertVerdict(evaluateWaterTanks(exactBoundary, 1, 2), 0);

    const ScalarSignals barelySafe{
        scalar({5.0000004}),
        scalar({5.0000004}),
    };
    assertVerdict(evaluateWaterTanks(barelySafe, 1, 2), 1);

    const ScalarSignals wholeRegion{
        scalar({4.0, 4.0, 4.0}),
        scalar({8.0, 4.0, 8.0}),
    };
    assertVerdict(evaluateWaterTanks(wholeRegion, 3, 3), 0);

    const ScalarSignals orderDependent{
        scalar({6.0, 0.0}),
        scalar({6.0, 12.0}),
    };
    assertVerdict(evaluateWaterTanks(orderDependent, 1, 2), 2);

    const ScalarSignals relativeViolation{
        {{0, 0.0}, {10, 100.0}},
        {{0, 100.0}, {5, 0.0}},
    };
    assertVerdict(evaluateWaterTanks(
        relativeViolation, 1, 20, Variant::AdmFineRelative, 50.0), 0);
    assertVerdict(evaluateWaterTanks(
        relativeViolation, 1, 20, Variant::AdmCoarseRelative, 50.0), 0);

    const ScalarSignals relativeSafe{
        {{0, 8.0}, {8, 2.0}},
        {{0, 3.0}, {6, 9.0}},
    };
    assertVerdict(evaluateWaterTanks(
        relativeSafe, 1, 10, Variant::AdmFineRelative), 1);

    const ScalarSignals signedConstants{
        scalar({-5.0}), scalar({-7.0})};
    assertVerdict(evaluateWaterTanks(
        signedConstants, 1, 2, Variant::AdmFine, -1.0), 0);
    assertVerdict(evaluateWaterTanks(
        signedConstants, 1, 2, Variant::AdmCoarse, -1.0), 0);
}

void testCoarseMatchesFine()
{
    const std::vector<double> values{-3.0, 0.0, 4.0, 12.0};
    for (const double firstBefore : values)
    {
        for (const double firstAfter : values)
        {
            for (const double secondBefore : values)
            {
                for (const double secondAfter : values)
                {
                    const ScalarSignals signals{
                        scalar({firstBefore, firstAfter}),
                        scalar({secondBefore, secondAfter}),
                    };
                    for (const double threshold : {-1.0, 10.0})
                    {
                        for (const bool relative : {false, true})
                        {
                            const Variant fine = relative
                                ? Variant::AdmFineRelative
                                : Variant::AdmFine;
                            const Variant coarse = relative
                                ? Variant::AdmCoarseRelative
                                : Variant::AdmCoarse;
                            const Evaluation fineResult = evaluateWaterTanks(
                                signals, 1, 3, fine, threshold);
                            const Evaluation coarseResult = evaluateWaterTanks(
                                signals, 1, 3, coarse, threshold);
                            assert(fineResult.falsePossible() ==
                                   coarseResult.falsePossible());
                            assert(fineResult.truePossible() ==
                                   coarseResult.truePossible());
                        }
                    }
                }
            }
        }
    }
}

void testMutualSeparationProperty()
{
    const std::vector<Variant> variants{
        Variant::Adm,
        Variant::AdmFine,
        Variant::AdmFineRelative,
    };

    for (const Variant variant : variants)
    {
        assertVerdict(evaluateMutualSeparation(
            {points({{0.0, 0.0, 0.0}}),
             points({{0.0, 0.0, 0.0}})},
            1, 2, variant), 0);
        assertVerdict(evaluateMutualSeparation(
            {points({{0.0, 0.0, 0.0}}),
             points({{0.0001, 0.0, 0.0}})},
            1, 2, variant), 1);

        assertVerdict(evaluateMutualSeparation(
            {points({{10.0, 10.0, 10.0}}),
             points({{0.0, 0.0, 0.0}}),
             points({{0.0, 0.0, 0.0}})},
            1, 2, variant), 0);

        assertVerdict(evaluateMutualSeparation(
            {points({{10.0, 10.0, 10.0}}),
             points({{5.0, 5.0, 5.0}}),
             points({{0.0, 0.0, 0.0}}),
             points({{0.0, 0.0, 0.0}})},
            1, 2, variant), 0);
    }

    const PointSignals finalCollision{
        points({{0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}}),
        points({{2.0, 2.0, 2.0}, {1.0, 1.0, 1.0}}),
    };
    assertVerdict(evaluateMutualSeparation(
        finalCollision, 1, 3, Variant::Adm), 0);

    bool rejected = false;
    try
    {
        (void)evaluateMutualSeparation(
            {points({{0.0, 0.0, 0.0}}), points({{1.0, 0.0, 0.0}})},
            1, 2, Variant::AdmCoarse);
    }
    catch (const std::invalid_argument &)
    {
        rejected = true;
    }
    assert(rejected);
}

void testNumericAndVerdictRepresentation()
{
    const double value = 1.0 / 10.0;
    assert(std::stod(formatDouble(value)) == value);
    assert(formatDouble(-0.0) == "0");

    const std::set<std::string> stuttered = stutter2kStr({"1;2"}, 3);
    assert(!stuttered.empty());
    for (const std::string &word : stuttered)
    {
        assert(std::count(word.begin(), word.end(), ';') == 2);
    }

    bool rejected = false;
    try
    {
        Evaluation empty;
        (void)empty.verdict();
    }
    catch (const std::logic_error &)
    {
        rejected = true;
    }
    assert(rejected);
}

void testRightEndpointMarker()
{
    ScalarSignal scalar{{0, 1.0}, {1000, 2.0}};
    discardRightEndpointMarker(scalar, 1000);
    assert(scalar.size() == 1);
    assert(scalar.front().first == 0);
    assert(scalar.front().second == 1.0);

    PointSignal pointsWithMarker{
        {0, {1.0, 2.0, 3.0}},
        {1000, {4.0, 5.0, 6.0}},
    };
    discardRightEndpointMarker(pointsWithMarker, 1000);
    assert(pointsWithMarker.size() == 1);

    PointSignal rightOpen{
        {0, {1.0, 2.0, 3.0}},
        {950, {4.0, 5.0, 6.0}},
    };
    discardRightEndpointMarker(rightOpen, 1000);
    assert(rightOpen.size() == 2);

    ScalarSignal beyond{{0, 1.0}, {1001, 2.0}};
    discardRightEndpointMarker(beyond, 1000);
    bool rejected = false;
    try
    {
        validateApproximateSignal(beyond, 1000);
    }
    catch (const std::invalid_argument &)
    {
        rejected = true;
    }
    assert(rejected);
}

} // namespace

int main()
{
    testSharedPaperAbstraction();
    testSyntheticFormulaRegressions();
    testApproximateBitsetKernels();
    testWaterTankProperty();
    testCoarseMatchesFine();
    testMutualSeparationProperty();
    testNumericAndVerdictRepresentation();
    testRightEndpointMarker();
    std::cout << "offline approximate tests passed\n";
    return 0;
}
