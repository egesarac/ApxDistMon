#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>
#include <random>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

#include "monitoring_case_studies.hpp"

namespace
{

using caseStudyDetail::Arithmetic;
using caseStudyDetail::NumericDag;
using caseStudyDetail::NumericWord;
using caseStudyDetail::WordLanguage;

std::string encodeWord(const NumericWord &word)
{
    assert(!word.empty());
    std::string result = formatDouble(word.front());
    for (std::size_t index = 1; index < word.size(); ++index)
    {
        result += ";" + formatDouble(word[index]);
    }
    return result;
}

std::set<std::string> encodeLanguage(const WordLanguage &language)
{
    std::set<std::string> result;
    for (const NumericWord &word : language)
    {
        result.insert(encodeWord(word));
    }
    return result;
}

std::set<std::string> legacyProduct(
    const WordLanguage &left,
    const WordLanguage &right,
    Arithmetic arithmetic)
{
    const std::vector<std::set<std::string>> legacyLeft{
        encodeLanguage(left)};
    const std::vector<std::set<std::string>> legacyRight{
        encodeLanguage(right)};
    const auto result = arithmetic == Arithmetic::Sum
        ? asyncProdStrSum(legacyLeft, legacyRight)
        : asyncProdStrDiffSqr(legacyLeft, legacyRight);
    assert(result.size() == 1);
    return result.front();
}

void assertCanonical(const WordLanguage &language)
{
    for (const NumericWord &word : language)
    {
        assert(!word.empty());
        for (std::size_t index = 0; index < word.size(); ++index)
        {
            assert(std::isfinite(word[index]));
            if (word[index] == 0.0)
            {
                assert(!std::signbit(word[index]));
            }
            if (index != 0)
            {
                assert(word[index - 1] != word[index]);
            }
        }
    }
}

WordLanguage typedProduct(
    const WordLanguage &left,
    const WordLanguage &right,
    Arithmetic arithmetic)
{
    const WordLanguage result = caseStudyDetail::combineWordLanguages(
        left, right, arithmetic);
    assertCanonical(result);
    return result;
}

WordLanguage enumerateDag(const NumericDag &language)
{
    WordLanguage result;
    NumericWord prefix;
    std::function<void(caseStudyDetail::NodeId)> visit =
        [&](caseStudyDetail::NodeId node)
    {
        const auto &current = language.nodes[node];
        const bool appended = prefix.empty() || prefix.back() != current.value;
        if (appended)
        {
            prefix.push_back(current.value);
        }
        if (current.accepting)
        {
            result.push_back(prefix);
        }
        for (caseStudyDetail::NodeId child : current.children)
        {
            visit(child);
        }
        if (appended)
        {
            prefix.pop_back();
        }
    };

    for (caseStudyDetail::NodeId root : language.roots)
    {
        visit(root);
    }
    std::sort(result.begin(), result.end());
    result.erase(std::unique(result.begin(), result.end()), result.end());
    assertCanonical(result);
    return result;
}

NumericDag dagFromLanguage(const WordLanguage &language)
{
    NumericDag result;
    for (const NumericWord &word : language)
    {
        caseStudyDetail::NodeId previous = 0;
        for (std::size_t token = 0; token < word.size(); ++token)
        {
            const auto node = static_cast<caseStudyDetail::NodeId>(
                result.nodes.size());
            result.nodes.push_back({
                caseStudyDetail::normalizeValue(word[token]),
                {},
                token + 1 == word.size()});
            if (token == 0)
            {
                result.roots.push_back(node);
            }
            else
            {
                result.nodes[previous].children.push_back(node);
            }
            previous = node;
        }
    }
    return caseStudyDetail::determinizeDag(result);
}

void assertEquivalent(
    const WordLanguage &left,
    const WordLanguage &right,
    Arithmetic arithmetic)
{
    const WordLanguage leftBefore = left;
    const WordLanguage rightBefore = right;
    const WordLanguage typed = typedProduct(left, right, arithmetic);
    assert(encodeLanguage(typed) == legacyProduct(left, right, arithmetic));
    assert(left == leftBefore);
    assert(right == rightBefore);
}

void generateCanonicalWords(
    NumericWord &prefix,
    std::size_t targetLength,
    WordLanguage &words)
{
    static constexpr std::array<double, 3> alphabet{-1.0, 0.0, 1.0};
    if (prefix.size() == targetLength)
    {
        words.push_back(prefix);
        return;
    }

    for (double value : alphabet)
    {
        if (!prefix.empty() && prefix.back() == value)
        {
            continue;
        }
        prefix.push_back(value);
        generateCanonicalWords(prefix, targetLength, words);
        prefix.pop_back();
    }
}

WordLanguage allCanonicalWords()
{
    WordLanguage words;
    NumericWord prefix;
    for (std::size_t length = 1; length <= 3; ++length)
    {
        generateCanonicalWords(prefix, length, words);
    }
    std::sort(words.begin(), words.end());
    words.erase(std::unique(words.begin(), words.end()), words.end());
    assert(words.size() == 21);
    assertCanonical(words);
    return words;
}

void testExhaustiveCanonicalWords()
{
    const WordLanguage words = allCanonicalWords();
    for (const NumericWord &leftWord : words)
    {
        for (const NumericWord &rightWord : words)
        {
            const WordLanguage left{leftWord};
            const WordLanguage right{rightWord};
            assertEquivalent(left, right, Arithmetic::Sum);
            assertEquivalent(left, right, Arithmetic::SquaredDifference);
        }
    }

    // Also exercise union and deduplication across a complete language.
    assertEquivalent(words, words, Arithmetic::Sum);
    assertEquivalent(words, words, Arithmetic::SquaredDifference);
}

void testDagProducts()
{
    const WordLanguage words = allCanonicalWords();
    for (const NumericWord &leftWord : words)
    {
        const NumericDag left = dagFromLanguage({leftWord});
        assert(enumerateDag(left) == WordLanguage{leftWord});
        for (const NumericWord &rightWord : words)
        {
            const NumericDag right = dagFromLanguage({rightWord});
            for (Arithmetic arithmetic : {
                     Arithmetic::Sum,
                     Arithmetic::SquaredDifference})
            {
                const WordLanguage expected = typedProduct(
                    {leftWord}, {rightWord}, arithmetic);
                const NumericDag product = caseStudyDetail::productDags(
                    left, right, arithmetic);
                assert(enumerateDag(product) == expected);
            }
        }
    }
}

void assertExpected(
    const NumericWord &leftWord,
    const NumericWord &rightWord,
    Arithmetic arithmetic,
    const std::set<std::string> &expected)
{
    const WordLanguage left{leftWord};
    const WordLanguage right{rightWord};
    assertEquivalent(left, right, arithmetic);
    assert(encodeLanguage(typedProduct(left, right, arithmetic)) == expected);
}

void testDiagonalAndCancellationCases()
{
    assertExpected(
        {0.0, 1.0},
        {0.0, 1.0},
        Arithmetic::SquaredDifference,
        {"0", "0;1;0"});

    assertExpected(
        {0.0, 1.0},
        {0.0, -1.0},
        Arithmetic::Sum,
        {"0", "0;-1;0", "0;1;0"});

    assertExpected(
        {1.0},
        {-1.0},
        Arithmetic::Sum,
        {"0"});
}

void testZeroIdentity()
{
    const WordLanguage zero{{0.0}};
    const WordLanguage language{
        {-1.0, 0.0, 1.0},
        {0.0, 1.0},
        {1.0, 0.0, -1.0}};
    assertCanonical(zero);
    assertCanonical(language);

    const WordLanguage rightIdentity = caseStudyDetail::combineWordLanguages(
        language, zero, Arithmetic::Sum);
    const WordLanguage leftIdentity = caseStudyDetail::combineWordLanguages(
        zero, language, Arithmetic::Sum);
    assert(rightIdentity == language);
    assert(leftIdentity == language);
    assertEquivalent(language, zero, Arithmetic::Sum);
    assertEquivalent(zero, language, Arithmetic::Sum);
    assertEquivalent(zero, zero, Arithmetic::Sum);
}

void testFloatingNormalization()
{
    assert(std::signbit(-0.0));
    assert(!std::signbit(caseStudyDetail::normalizeValue(-0.0)));

    const WordLanguage cancellation = typedProduct(
        {{-1.0}}, {{1.0}}, Arithmetic::Sum);
    assert(cancellation == WordLanguage{{0.0}});
    assert(!std::signbit(cancellation.front().front()));
    assertEquivalent({{-1.0}}, {{1.0}}, Arithmetic::Sum);

    // Avoid the zero-language shortcut so that negative-zero arithmetic is
    // normalized by every grid cell before prefix interning.
    const WordLanguage signedZero = typedProduct(
        {{-0.0, 1.0}}, {{-0.0, -1.0}}, Arithmetic::Sum);
    assertCanonical(signedZero);
    assertEquivalent(
        {{-0.0, 1.0}}, {{-0.0, -1.0}}, Arithmetic::Sum);

    const WordLanguage fraction = typedProduct(
        {{0.1}}, {{0.2}}, Arithmetic::Sum);
    assert(fraction.size() == 1);
    assert(fraction.front().size() == 1);
    assert(formatDouble(fraction.front().front()) ==
           "0.30000000000000004");
    assertEquivalent({{0.1}}, {{0.2}}, Arithmetic::Sum);

    const double neighbor = std::nextafter(0.1, 1.0);
    assertEquivalent({{neighbor}}, {{-0.1}}, Arithmetic::Sum);

    // The mathematical square is nonzero, but binary64 underflows to zero.
    const double tiny = std::numeric_limits<double>::min();
    const WordLanguage underflow = typedProduct(
        {{tiny}}, {{0.0}}, Arithmetic::SquaredDifference);
    assert(underflow == WordLanguage{{0.0}});
    assert(!std::signbit(underflow.front().front()));
    assertEquivalent(
        {{tiny}}, {{0.0}}, Arithmetic::SquaredDifference);
}

template <typename Action>
bool throwsOverflow(Action action)
{
    try
    {
        action();
    }
    catch (const std::overflow_error &)
    {
        return true;
    }
    return false;
}

void assertOverflow(
    const WordLanguage &left,
    const WordLanguage &right,
    Arithmetic arithmetic)
{
    const bool typedOverflow = throwsOverflow([&]
    {
        (void)caseStudyDetail::combineWordLanguages(
            left, right, arithmetic);
    });
    const bool legacyOverflow = throwsOverflow([&]
    {
        (void)legacyProduct(left, right, arithmetic);
    });
    assert(typedOverflow);
    assert(legacyOverflow);
}

void testOverflow()
{
    const double maximum = std::numeric_limits<double>::max();
    assertOverflow({{maximum}}, {{maximum}}, Arithmetic::Sum);
    assertOverflow(
        {{maximum}}, {{-maximum}}, Arithmetic::SquaredDifference);
    assert(throwsOverflow([&]
    {
        (void)caseStudyDetail::normalizeValue(
            std::numeric_limits<double>::infinity());
    }));
}

void assertSegmentWordsMatch(
    const ScalarSignals &signals,
    long long epsilon,
    long long duration,
    const std::set<int> &exactSignals = {})
{
    const auto uncertainties = computeUncertaintyIntervals(
        signals, epsilon, duration, exactSignals);
    const auto segmentation = computeCanonicalSegmentation(
        signals, uncertainties, duration, exactSignals);
    const auto legacy = computeValueExpressions(
        signals, uncertainties, segmentation, exactSignals);
    const auto prepared = caseStudyDetail::prepareCaseStudy(
        signals, epsilon, duration, exactSignals);

    assert(prepared.segmentation == segmentation);
    assert(prepared.signals.size() == signals.size());
    assert(legacy.size() == signals.size());
    std::vector<caseStudyDetail::SweepCursor> cursors(signals.size());

    for (std::size_t segment = 0;
         segment + 1 < segmentation.size();
         ++segment)
    {
        const long long start = segmentation[segment];
        const long long end = segmentation[segment + 1];
        for (std::size_t signal = 0; signal < signals.size(); ++signal)
        {
            const caseStudyDetail::SegmentView view =
                caseStudyDetail::advanceSweep(
                    prepared.signals[signal], cursors[signal], start, end);
            const WordLanguage typed = caseStudyDetail::segmentWords(
                prepared.signals[signal], view, start, end);
            assertCanonical(typed);
            assert(encodeLanguage(typed) == legacy[signal][segment]);
            const NumericDag dag = caseStudyDetail::segmentDag(
                prepared.signals[signal], view, start, end);
            assert(enumerateDag(dag) == typed);
        }
    }
}

void testSyntheticSegmentWords()
{
    assertSegmentWordsMatch(
        {{{0, 0.0}, {1, 1.0}, {2, 2.0}},
         {{0, -0.0}, {3, -1.0}},
         {{0, 0.1}, {1, 0.1}, {4, 0.2}}},
        2,
        6);

    assertSegmentWordsMatch(
        {{{0, 0.0}, {2, 1.0}, {3, 0.0}},
         {{0, 4.0}, {1, 4.0}, {4, -1.0}}},
        2,
        6,
        {0});

    assertSegmentWordsMatch(
        {{{0, 0.0}, {5, 1.0}, {6, 0.0}}},
        10,
        11);
}

void testRandomSegmentWords()
{
    std::mt19937 generator(0x4e554d45U);
    static constexpr std::array<double, 6> values{
        -1.0, -0.0, 0.0, 0.1, 1.0, 2.0};

    for (int trial = 0; trial < 80; ++trial)
    {
        const long long duration = 6 + generator() % 5;
        const long long epsilon = 1 + generator() % duration;
        const std::size_t signalCount = 1 + generator() % 3;
        ScalarSignals signals(signalCount);
        std::set<int> exactSignals;

        for (std::size_t signal = 0; signal < signalCount; ++signal)
        {
            signals[signal].push_back({
                0, values[generator() % values.size()]});
            long long timestamp = 1 + generator() % 2;
            while (timestamp < duration && signals[signal].size() < 5)
            {
                signals[signal].push_back({
                    timestamp, values[generator() % values.size()]});
                timestamp += 1 + generator() % 3;
            }
            if (generator() % 4 == 0)
            {
                exactSignals.insert(static_cast<int>(signal));
            }
        }

        assertSegmentWordsMatch(
            signals, epsilon, duration, exactSignals);
    }
}

void testMultiInputBooleanSum()
{
    const WordLanguage words = allCanonicalWords();
    std::mt19937 generator(0x53554d44U);
    for (int trial = 0; trial < 300; ++trial)
    {
        const std::size_t operandCount = 2 + generator() % 3;
        std::vector<NumericDag> operands;
        operands.reserve(operandCount);
        WordLanguage expected{
            words[generator() % words.size()]};
        operands.push_back(dagFromLanguage(expected));
        for (std::size_t operand = 1; operand < operandCount; ++operand)
        {
            const WordLanguage next{
                words[generator() % words.size()]};
            operands.push_back(dagFromLanguage(next));
            expected = typedProduct(expected, next, Arithmetic::Sum);
        }

        const double threshold = static_cast<int>(generator() % 5) - 2.0;
        std::vector<std::bitset<SIZE>> actual(2);
        std::vector<std::bitset<SIZE>> reference(2);
        caseStudyDetail::includeBooleanSumProduct(
            actual, operands, threshold);
        for (const NumericWord &word : expected)
        {
            caseStudyDetail::includeBooleanWord(
                reference, word, threshold);
        }
        assert(actual == reference);
    }
}

caseStudyDetail::RawLanguage rawLanguageFromPaths(
    const std::vector<NumericWord> &paths,
    bool includeEmpty)
{
    caseStudyDetail::RawLanguage result{{caseStudyDetail::RawState{}}, 0};
    std::vector<std::size_t> tails;
    for (const NumericWord &path : paths)
    {
        assert(!path.empty());
        const std::size_t branch = result.states.size();
        result.states.push_back({});
        result.states[0].transitions.push_back({branch, 0.0, true});
        std::size_t current = branch;
        for (std::size_t token = 0; token < path.size(); ++token)
        {
            const std::size_t target = result.states.size();
            result.states.push_back({});
            result.states[current].transitions.push_back({
                target, path[token], false});
            current = target;
            if (token + 1 != path.size())
            {
                const std::size_t bridge = result.states.size();
                result.states.push_back({});
                result.states[current].transitions.push_back({
                    bridge, 0.0, true});
                current = bridge;
            }
        }
        tails.push_back(current);
    }

    result.finalState = result.states.size();
    result.states.push_back({});
    for (std::size_t tail : tails)
    {
        result.states[tail].transitions.push_back({
            result.finalState, 0.0, true});
    }
    if (includeEmpty)
    {
        result.states[0].transitions.push_back({
            result.finalState, 0.0, true});
    }
    return result;
}

NumericWord canonicalizeRawPath(const NumericWord &path)
{
    NumericWord result;
    for (double value : path)
    {
        value = caseStudyDetail::normalizeValue(value);
        if (result.empty() || result.back() != value)
        {
            result.push_back(value);
        }
    }
    return result;
}

void assertSparseSumMatches(
    const std::vector<std::vector<NumericWord>> &rawPaths,
    double threshold)
{
    std::vector<caseStudyDetail::RawLanguage> rawOperands;
    std::vector<WordLanguage> wordOperands;
    for (const auto &operandPaths : rawPaths)
    {
        rawOperands.push_back(rawLanguageFromPaths(operandPaths, true));
        WordLanguage words;
        for (const NumericWord &path : operandPaths)
        {
            words.push_back(canonicalizeRawPath(path));
        }
        std::sort(words.begin(), words.end());
        words.erase(std::unique(words.begin(), words.end()), words.end());
        wordOperands.push_back(std::move(words));
    }

    WordLanguage expected = wordOperands.front();
    for (std::size_t operand = 1; operand < wordOperands.size(); ++operand)
    {
        expected = typedProduct(
            expected, wordOperands[operand], Arithmetic::Sum);
    }

    std::vector<std::bitset<SIZE>> actual(2);
    std::vector<std::bitset<SIZE>> reference(2);
    caseStudyDetail::includeBooleanSparseSum(
        actual, rawOperands, threshold);
    for (const NumericWord &word : expected)
    {
        caseStudyDetail::includeBooleanWord(reference, word, threshold);
    }
    assert(actual == reference);
}

void testSparseMultiInputBooleanSum()
{
    const std::vector<std::vector<NumericWord>> fixed{
        {{-1.0, -1.0, 1.0}, {0.0, 0.0, 2.0}},
        {{1.0, 1.0, -1.0}, {0.0, 1.0, 1.0}},
        {{0.0, 0.0}, {2.0, 2.0, 0.0}},
        {{0.0, -2.0, -2.0}, {1.0, 1.0}}};
    for (double threshold : {-2.0, -1.0, 0.0, 1.0, 2.0})
    {
        assertSparseSumMatches(fixed, threshold);
    }

    std::mt19937 generator(0x53505253U);
    for (std::size_t operandCount = 2; operandCount <= 4; ++operandCount)
    {
        for (int trial = 0; trial < 40; ++trial)
        {
            std::vector<std::vector<NumericWord>> rawPaths(operandCount);
            for (auto &operandPaths : rawPaths)
            {
                for (int alternative = 0; alternative < 2; ++alternative)
                {
                    NumericWord path;
                    const std::size_t length = 1 + generator() % 2;
                    for (std::size_t token = 0; token < length; ++token)
                    {
                        const double value =
                            static_cast<int>(generator() % 5) - 2.0;
                        path.push_back(value);
                        if (generator() % 2 == 0)
                        {
                            path.push_back(value);
                        }
                    }
                    operandPaths.push_back(std::move(path));
                }
            }
            const double threshold = static_cast<int>(generator() % 5) - 2.0;
            assertSparseSumMatches(rawPaths, threshold);
        }
    }
}

namespace caseStudyLegacy
{

using NumericExpressions = std::vector<std::vector<std::set<std::string>>>;
using NumericLanguage = std::vector<std::set<std::string>>;

struct Abstraction
{
    std::vector<long long> segmentation;
    NumericExpressions valueExpressions;
};

Abstraction buildAbstraction(
    const ScalarSignals &signals,
    long long epsilon,
    long long duration,
    const std::set<int> &exactSignals)
{
    const auto uncertainties = computeUncertaintyIntervals(
        signals, epsilon, duration, exactSignals);
    std::vector<long long> segmentation = computeCanonicalSegmentation(
        signals, uncertainties, duration, exactSignals);
    NumericExpressions expressions = computeValueExpressions(
        signals, uncertainties, segmentation, exactSignals);
    return {std::move(segmentation), std::move(expressions)};
}

NumericLanguage sumLanguages(
    const NumericExpressions &expressions,
    bool fine)
{
    NumericLanguage result = expressions.front();
    for (std::size_t signal = 1; signal < expressions.size(); ++signal)
    {
        result = fine
            ? abstProdStrSum(result, expressions[signal])
            : asyncProdStrSum(result, expressions[signal]);
    }
    return result;
}

NumericLanguage coarseSum(const NumericExpressions &expressions)
{
    const std::size_t segmentCount = expressions.front().size();
    NumericLanguage result(segmentCount);
    for (std::size_t segment = 0; segment < segmentCount; ++segment)
    {
        double sumOfMinima = 0.0;
        double sumOfMaxima = 0.0;
        for (const NumericLanguage &signal : expressions)
        {
            double minimum = std::numeric_limits<double>::infinity();
            double maximum = -std::numeric_limits<double>::infinity();
            for (const std::string &expression : signal[segment])
            {
                for (double value : splitSet(expression, ";"))
                {
                    minimum = std::min(minimum, value);
                    maximum = std::max(maximum, value);
                }
            }
            sumOfMinima += minimum;
            sumOfMaxima += maximum;
        }
        result[segment].insert(formatDouble(sumOfMinima));
        result[segment].insert(formatDouble(sumOfMaxima));
    }
    return result;
}

BooleanLanguage invariantLanguage(
    const NumericLanguage &expressions,
    double threshold)
{
    const auto propositions = convertSignalsToAtomicPropositions(
        NumericExpressions{expressions}, threshold);
    return bitsetAlways(propositions[0]);
}

Evaluation evaluateWaterTanks(
    const ScalarSignals &signals,
    long long epsilon,
    long long duration,
    Variant variant,
    double threshold)
{
    const std::set<int> exactSignals = isRelative(variant)
        ? std::set<int>{0}
        : std::set<int>{};
    Abstraction abstraction = buildAbstraction(
        signals, epsilon, duration, exactSignals);
    NumericLanguage sum = isCoarse(variant)
        ? coarseSum(abstraction.valueExpressions)
        : sumLanguages(abstraction.valueExpressions, isFine(variant));

    evenMask = generateBitmask(0);
    oddMask = generateBitmask(1);
    return {
        std::move(abstraction.segmentation),
        invariantLanguage(sum, threshold)};
}

ScalarSignals flattenPointSignals(const PointSignals &signals)
{
    ScalarSignals flattened(signals.size() * 3);
    for (std::size_t agent = 0; agent < signals.size(); ++agent)
    {
        for (const auto &[timestamp, point] : signals[agent])
        {
            for (std::size_t axis = 0; axis < 3; ++axis)
            {
                flattened[3 * agent + axis].push_back({
                    timestamp, point[axis]});
            }
        }
    }
    return flattened;
}

NumericLanguage squaredDistance(
    const NumericExpressions &expressions,
    std::size_t left,
    std::size_t right,
    bool fine)
{
    auto add = [fine](
        const NumericLanguage &first,
        const NumericLanguage &second)
    {
        return fine
            ? abstProdStrSum(first, second)
            : asyncProdStrSum(first, second);
    };

    std::array<NumericLanguage, 3> squaredCoordinates;
    for (std::size_t axis = 0; axis < 3; ++axis)
    {
        const NumericLanguage &first = expressions[3 * left + axis];
        const NumericLanguage &second = expressions[3 * right + axis];
        squaredCoordinates[axis] = fine
            ? abstProdStrDiffSqr(first, second)
            : asyncProdStrDiffSqr(first, second);
    }
    return add(
        add(squaredCoordinates[0], squaredCoordinates[1]),
        squaredCoordinates[2]);
}

Evaluation evaluateMutualSeparation(
    const PointSignals &signals,
    long long epsilon,
    long long duration,
    Variant variant,
    double threshold)
{
    ScalarSignals flattened = flattenPointSignals(signals);
    const std::set<int> exactSignals = isRelative(variant)
        ? std::set<int>{0, 1, 2}
        : std::set<int>{};
    Abstraction abstraction = buildAbstraction(
        flattened, epsilon, duration, exactSignals);

    evenMask = generateBitmask(0);
    oddMask = generateBitmask(1);
    const double squaredThreshold = threshold * threshold;
    std::optional<BooleanLanguage> result;
    for (std::size_t left = 0; left < signals.size(); ++left)
    {
        for (std::size_t right = left + 1;
             right < signals.size();
             ++right)
        {
            BooleanLanguage pairInvariant = invariantLanguage(
                squaredDistance(
                    abstraction.valueExpressions,
                    left,
                    right,
                    isFine(variant)),
                squaredThreshold);
            result = result
                ? std::optional<BooleanLanguage>(
                    bitsetConjunction(*result, pairInvariant))
                : std::optional<BooleanLanguage>(
                    std::move(pairInvariant));
        }
    }
    return {std::move(abstraction.segmentation), std::move(*result)};
}

} // namespace caseStudyLegacy

void assertSameEvaluation(
    const Evaluation &actual,
    const Evaluation &reference)
{
    assert(actual.segmentation == reference.segmentation);
    assert(actual.language == reference.language);
    assert(actual.verdict() == reference.verdict());
}

void testWaterAdmUniformTruthValues()
{
    const auto assertMatchesLegacy = [](
        const ScalarSignals &signals,
        long long epsilon,
        long long duration,
        double threshold)
    {
        const Evaluation actual = evaluateWaterTanks(
            signals, epsilon, duration, Variant::Adm, threshold);
        const Evaluation reference = caseStudyLegacy::evaluateWaterTanks(
            signals, epsilon, duration, Variant::Adm, threshold);
        assertSameEvaluation(actual, reference);
    };

    assertMatchesLegacy(
        {{{0, 2.0}}, {{0, 3.0}}},
        1,
        2,
        4.0);
    assertMatchesLegacy(
        {{{0, 1.0}}, {{0, 1.0}}},
        1,
        2,
        2.0);

    // The lower bound equals the strict threshold, so both truth values are
    // possible under the strict comparison.
    assertMatchesLegacy(
        {{{0, 0.0}, {1, 1.0}}, {{0, 0.0}}},
        2,
        3,
        0.0);
}

void testWholeEvaluatorDifferential()
{
    std::mt19937 generator(0x4556414cU);
    static constexpr std::array<double, 5> values{
        -2.0, -0.0, 0.0, 1.0, 3.0};
    const std::array<Variant, 5> waterVariants{
        Variant::Adm,
        Variant::AdmFine,
        Variant::AdmFineRelative,
        Variant::AdmCoarse,
        Variant::AdmCoarseRelative};

    for (int trial = 0; trial < 18; ++trial)
    {
        const std::size_t agentCount = 2 + generator() % 3;
        ScalarSignals signals(agentCount);
        for (ScalarSignal &signal : signals)
        {
            signal.push_back({0, values[generator() % values.size()]});
            signal.push_back({2, values[generator() % values.size()]});
        }
        for (Variant variant : waterVariants)
        {
            const Evaluation actual = evaluateWaterTanks(
                signals, 1, 4, variant, 0.5);
            const Evaluation reference = caseStudyLegacy::evaluateWaterTanks(
                signals, 1, 4, variant, 0.5);
            assertSameEvaluation(actual, reference);
        }
    }

    for (int trial = 0; trial < 6; ++trial)
    {
        ScalarSignals signals(2);
        for (ScalarSignal &signal : signals)
        {
            signal.push_back({0, values[generator() % values.size()]});
            signal.push_back({2, values[generator() % values.size()]});
            signal.push_back({4, values[generator() % values.size()]});
        }
        for (Variant variant : waterVariants)
        {
            const Evaluation actual = evaluateWaterTanks(
                signals, 2, 7, variant, 0.5);
            const Evaluation reference = caseStudyLegacy::evaluateWaterTanks(
                signals, 2, 7, variant, 0.5);
            assertSameEvaluation(actual, reference);
        }
    }

    const std::array<Variant, 3> separationVariants{
        Variant::Adm,
        Variant::AdmFine,
        Variant::AdmFineRelative};
    for (int trial = 0; trial < 20; ++trial)
    {
        PointSignals signals(2 + generator() % 3);
        for (PointSignal &signal : signals)
        {
            Point first{};
            Point second{};
            for (std::size_t axis = 0; axis < 3; ++axis)
            {
                first[axis] = values[generator() % values.size()];
                second[axis] = values[generator() % values.size()];
            }
            signal.push_back({0, first});
            signal.push_back({2, second});
        }
        for (Variant variant : separationVariants)
        {
            const Evaluation actual = evaluateMutualSeparation(
                signals, 1, 4, variant, 0.5);
            const Evaluation reference =
                caseStudyLegacy::evaluateMutualSeparation(
                signals, 1, 4, variant, 0.5);
            assertSameEvaluation(actual, reference);
        }
    }
}

} // namespace

int main()
{
    testExhaustiveCanonicalWords();
    testDagProducts();
    testDiagonalAndCancellationCases();
    testZeroIdentity();
    testFloatingNormalization();
    testOverflow();
    testSyntheticSegmentWords();
    testRandomSegmentWords();
    testMultiInputBooleanSum();
    testSparseMultiInputBooleanSum();
    testWaterAdmUniformTruthValues();
    testWholeEvaluatorDifferential();
    std::cout << "case-study numeric tests passed\n";
    return 0;
}
