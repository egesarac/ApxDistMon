#include <bitset>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "monitoring.hpp"

namespace
{

using Signal = std::vector<std::pair<long long, double>>;
using Signals = std::vector<Signal>;
using SegmentedLanguage = std::vector<std::vector<std::bitset<SIZE>>>;

constexpr long long kMillisecondsPerSample = 1000;

Signal signalFromWord(const std::string &word)
{
    if (word.empty())
    {
        throw std::invalid_argument("Boolean words must not be empty");
    }
    if (word[0] != '0' && word[0] != '1')
    {
        throw std::invalid_argument("Boolean words may contain only 0 and 1");
    }

    Signal signal{{0, static_cast<double>(word[0] - '0')}};
    for (std::size_t index = 1; index < word.size(); ++index)
    {
        if (word[index] != '0' && word[index] != '1')
        {
            throw std::invalid_argument("Boolean words may contain only 0 and 1");
        }
        if (word[index] != word[index - 1])
        {
            signal.push_back({
                static_cast<long long>(index) * kMillisecondsPerSample,
                static_cast<double>(word[index] - '0')});
        }
    }
    return signal;
}

SegmentedLanguage evaluateFormula(
    const std::string &formula,
    const std::vector<std::vector<std::vector<std::bitset<SIZE>>>> &aps,
    const std::vector<long long> &segmentation)
{
    const auto &p = aps[0];
    const auto &q = aps[1];

    if (formula == "phi1")
    {
        return bitsetAlways(bitsetConjunction(p, q));
    }
    if (formula == "phi2")
    {
        return bitsetAlways(bitsetNegation(
            bitsetConjunction(bitsetNegation(p), bitsetNegation(q))));
    }
    if (formula == "phi3")
    {
        return bitsetUnboundedUntil(p, q, segmentation, 0, true);
    }
    if (formula == "phi4")
    {
        return bitsetAlways(bitsetNegation(
            bitsetConjunction(p, bitsetAlways(bitsetNegation(q)))));
    }

    long long upper = 0;
    if (formula == "phi5")
    {
        upper = 1000;
    }
    else if (formula == "phi6")
    {
        upper = 2000;
    }
    else
    {
        throw std::invalid_argument("unknown formula: " + formula);
    }

    return bitsetAlways(bitsetNegation(bitsetConjunction(
        p,
        bitsetBoundedAlways(
            bitsetNegation(q),
            segmentation,
            0,
            upper,
            true,
            false))));
}

struct Result
{
    int segments;
    bool falsePossible;
    bool truePossible;

    int verdict() const
    {
        if (falsePossible && truePossible)
        {
            return 2;
        }
        if (truePossible)
        {
            return 1;
        }
        if (falsePossible)
        {
            return 0;
        }
        throw std::logic_error("approximate monitor produced no verdict");
    }
};

Result evaluate(
    const std::string &formula,
    long long epsilonSeconds,
    const std::string &pWord,
    const std::string &qWord)
{
    if (pWord.size() != qWord.size())
    {
        throw std::invalid_argument("p and q words must have equal lengths");
    }
    if (epsilonSeconds <= 0)
    {
        throw std::invalid_argument("epsilon must be positive");
    }

    const long long duration =
        static_cast<long long>(pWord.size()) * kMillisecondsPerSample;
    const long long epsilon = epsilonSeconds * kMillisecondsPerSample;
    const Signals signals{signalFromWord(pWord), signalFromWord(qWord)};
    const auto regions = computeUncertaintyIntervals(
        signals, epsilon, duration);
    const auto segmentation = computeCanonicalSegmentation(
        signals, regions, duration);
    const auto gamma = computeValueExpressions(
        signals, regions, segmentation);
    const auto aps = convertSignalsToAtomicPropositions(gamma, 0.0);
    const auto root = evaluateFormula(formula, aps, segmentation);

    return {
        static_cast<int>(segmentation.size()) - 1,
        root[0][0].any(),
        root[0][1].any(),
    };
}

} // namespace

int main()
{
    evenMask = generateBitmask(0);
    oddMask = generateBitmask(1);
    allExceptFirstMask = evenMask | oddMask;
    allExceptFirstMask[0] = false;

    std::size_t query = 0;
    std::string formula;
    long long epsilon = 0;
    std::string pWord;
    std::string qWord;
    while (std::cin >> query >> formula >> epsilon >> pWord >> qWord)
    {
        try
        {
            const Result result = evaluate(formula, epsilon, pWord, qWord);
            std::cout << query << ' '
                      << result.segments << ' '
                      << static_cast<int>(result.falsePossible) << ' '
                      << static_cast<int>(result.truePossible) << ' '
                      << result.verdict() << '\n';
        }
        catch (const std::exception &error)
        {
            std::cerr << "query " << query << ": " << error.what() << '\n';
            return 1;
        }
    }
    return 0;
}
