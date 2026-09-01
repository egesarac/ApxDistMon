#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <random>
#include <tuple>
#include <map>
#include <set>
#include <bitset>
#include <algorithm>
#include <limits>
#include <cfloat>
#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <initializer_list>
#include <optional>
#include <stdexcept>
#include "monitoring.hpp"
using namespace std;

namespace fs = std::filesystem;

namespace
{

struct Options
{
    optional<int> duration;
    optional<int> epsilon;
    optional<int> chunkSize;
    optional<int> sample;
    int repetitions = 5;
    fs::path dataDirectory = "data/signals";
    fs::path outputDirectory = "results/online/incremental";
    bool quiet = false;
};

void printUsage(const char *program)
{
    cout
        << "Usage: " << program << " [options]\n"
        << "  --duration <16|32|64|128|256|512>  Restrict trace duration in seconds\n"
        << "  --epsilon <2|4|8>                  Restrict clock skew in seconds\n"
        << "  --chunk <4|8|16|32>                Restrict chunk size in seconds\n"
        << "  --sample <0..99>                   Restrict the random sample\n"
        << "  --repeat <count>                   Timed repetitions (default: 5)\n"
        << "  --data-dir <path>                  Signal dataset directory\n"
        << "  --output-dir <path>                Generated result directory\n"
        << "  --quiet                            Do not echo result rows\n"
        << "  --help                             Show this message\n";
}

int parseInteger(const string &value, const string &option)
{
    size_t parsed = 0;
    int result = 0;
    try
    {
        result = stoi(value, &parsed);
    }
    catch (const exception &)
    {
        throw runtime_error(option + " expects an integer: " + value);
    }
    if (parsed != value.size())
    {
        throw runtime_error(option + " expects an integer: " + value);
    }
    return result;
}

Options parseOptions(int argc, char **argv)
{
    Options options;
    for (int index = 1; index < argc; ++index)
    {
        const string argument = argv[index];
        auto valueFor = [&](const string &option)
        {
            if (index + 1 >= argc)
            {
                throw runtime_error("missing value for " + option);
            }
            return string(argv[++index]);
        };

        if (argument == "--duration")
        {
            options.duration = parseInteger(valueFor(argument), argument);
        }
        else if (argument == "--epsilon")
        {
            options.epsilon = parseInteger(valueFor(argument), argument);
        }
        else if (argument == "--chunk")
        {
            options.chunkSize = parseInteger(valueFor(argument), argument);
        }
        else if (argument == "--sample")
        {
            options.sample = parseInteger(valueFor(argument), argument);
        }
        else if (argument == "--repeat")
        {
            options.repetitions = parseInteger(valueFor(argument), argument);
        }
        else if (argument == "--data-dir")
        {
            options.dataDirectory = valueFor(argument);
        }
        else if (argument == "--output-dir")
        {
            options.outputDirectory = valueFor(argument);
        }
        else if (argument == "--quiet")
        {
            options.quiet = true;
        }
        else if (argument == "--help")
        {
            printUsage(argv[0]);
            exit(0);
        }
        else
        {
            throw runtime_error("unknown option: " + argument);
        }
    }

    if (options.repetitions < 1)
    {
        throw runtime_error("--repeat must be positive");
    }
    if (options.sample && (*options.sample < 0 || *options.sample > 99))
    {
        throw runtime_error("--sample must be between 0 and 99");
    }
    return options;
}

vector<long long> selectedMilliseconds(
    initializer_list<int> defaults,
    const optional<int> &selected,
    const string &option)
{
    if (selected && find(defaults.begin(), defaults.end(), *selected) == defaults.end())
    {
        throw runtime_error(
            "unsupported value for " + option + ": " + to_string(*selected));
    }

    vector<long long> values;
    if (selected)
    {
        values.push_back(static_cast<long long>(*selected) * 1000);
        return values;
    }
    for (const int value : defaults)
    {
        values.push_back(static_cast<long long>(value) * 1000);
    }
    return values;
}

void run(const Options &options)
{
    const vector<long long> D = selectedMilliseconds(
        {16, 32, 64, 128, 256, 512}, options.duration, "--duration");
    const vector<long long> EPS = selectedMilliseconds(
        {2, 4, 8}, options.epsilon, "--epsilon");
    const vector<long long> CHUNK = selectedMilliseconds(
        {4, 8, 16, 32}, options.chunkSize, "--chunk");

    const int nSignals = 2;
    const long long H_PHI = 0;

    fs::create_directories(options.outputDirectory);
    const fs::path outputPath = options.outputDirectory / "onlineIncr_ec.txt";
    ofstream results(outputPath);
    if (!results)
    {
        throw runtime_error("cannot open output file: " + outputPath.string());
    }

    const int firstSample = options.sample.value_or(0);
    const int lastSample = options.sample.value_or(99);

    for (auto d : D)
    {
        for (auto eps : EPS)
        {
            for (auto CHUNK_SZ : CHUNK)
            {
                if (d < eps)
                    continue;

                if (CHUNK_SZ <= eps || CHUNK_SZ > d)
                    continue;

                for (int c = firstSample; c <= lastSample; ++c)
                {
                    vector<vector<pair<long long, double>>> signalsReal(nSignals);
                    signalsReal[0] = getData(
                        (options.dataDirectory / (to_string(d / 1000) + "_" + to_string(c) + ".txt")).string(),
                        d / 1000);
                    signalsReal[1] = getData(
                        (options.dataDirectory / (to_string(d / 1000) + "_" + to_string(c + 100) + ".txt")).string(),
                        d / 1000);

                    evenMask = generateBitmask(0);
                    oddMask = generateBitmask(1);
                    allExceptFirstMask = evenMask | oddMask;
                    allExceptFirstMask[0] = 0;

                    auto runMonitor = [&]()
                    {
                        string sequence;
                        long long dPrevious = 0;
                        long long qPrevious = 0;

                        vector<vector<pair<long long, double>>> retainedSignals(nSignals);
                        vector<optional<long long>> predecessorLowerBounds(nSignals);
                        vector<int> nextSample(nSignals, 1);
                        for (int i = 0; i < nSignals; i++)
                        {
                            if (signalsReal[i].empty())
                            {
                                throw logic_error("an input signal has no initial value");
                            }
                            double value = signalsReal[i][0].second > 0 ? 1 : 0;
                            retainedSignals[i].push_back(make_pair(0, value));
                        }

                        vector<long long> retainedSegmentation{0};
                        vector<vector<bitset<SIZE>>> retainedP;
                        vector<vector<bitset<SIZE>>> retainedQ;
                        vector<vector<bitset<SIZE>>> retainedConjunction;
                        vector<vector<bitset<SIZE>>> retainedRoot;

                        while (dPrevious < d)
                        {
                            long long dCurrent = min(d, dPrevious + CHUNK_SZ);

                            int readEnd = dCurrent / 1000;
                            for (int i = 0; i < nSignals; i++)
                            {
                                while (nextSample[i] < readEnd)
                                {
                                    bool previous = signalsReal[i][nextSample[i] - 1].second > 0;
                                    bool current = signalsReal[i][nextSample[i]].second > 0;
                                    if (previous != current)
                                    {
                                        pair<long long, double> edge = make_pair(
                                            signalsReal[i][nextSample[i]].first,
                                            current ? 1 : 0);
                                        retainedSignals[i].push_back(edge);
                                    }
                                    nextSample[i]++;
                                }
                            }

                            long long qCurrent = max(0LL, dCurrent - eps);
                            if (qCurrent < qPrevious)
                            {
                                throw logic_error("the finalized watermark moved backwards");
                            }

                            const long long supportHorizon =
                                onlineUncertaintyHorizon(dCurrent, eps);
                            auto uncertainties = computeUncertaintyIntervals(
                                retainedSignals,
                                eps,
                                supportHorizon,
                                {},
                                predecessorLowerBounds);

                            if (qCurrent == qPrevious)
                            {
                                if (qCurrent == 0)
                                {
                                    sequence.push_back('2');
                                }
                                else
                                {
                                    sequence.push_back(bitsetVerdict(retainedRoot.back()));
                                }

                                long long inputCutoff = max(0LL, qCurrent - eps);
                                for (int i = 0; i < nSignals; i++)
                                {
                                    auto retainedInput = retainApproximateSignalFrom(
                                        retainedSignals[i],
                                        uncertainties[i],
                                        inputCutoff,
                                        predecessorLowerBounds[i]);
                                    retainedSignals[i] = move(retainedInput.signal);
                                    predecessorLowerBounds[i] =
                                        retainedInput.predecessorLowerBound;
                                }

                                dPrevious = dCurrent;
                                continue;
                            }

                            auto canonical = computeCanonicalSegmentation(
                                retainedSignals, uncertainties, supportHorizon);
                            auto newSegmentation = refineSegmentation(
                                canonical,
                                {},
                                qPrevious,
                                qCurrent);

                            if (newSegmentation.size() < 2 ||
                                newSegmentation.front() != qPrevious ||
                                newSegmentation.back() != qCurrent)
                            {
                                throw logic_error("invalid refined incremental segmentation");
                            }

                            auto valExprs = computeValueExpressions(
                                retainedSignals,
                                uncertainties,
                                newSegmentation);
                            auto aps = convertSignalsToAtomicPropositions(valExprs, 0.0);
                            auto newP = aps[0];
                            auto newQ = aps[1];

                            auto newConjunction = bitsetConjunction(newP, newQ);

                            bool root0 = true;
                            bool root1 = false;
                            if (!retainedRoot.empty())
                            {
                                tie(root0, root1) = bitsetLastBits(retainedRoot.back());
                            }
                            auto newRoot = bitsetEventuallyPast(
                                newConjunction,
                                0,
                                newConjunction.size(),
                                root0,
                                root1);

                            retainedSegmentation.insert(
                                retainedSegmentation.end(),
                                newSegmentation.begin() + 1,
                                newSegmentation.end());
                            retainedP.insert(retainedP.end(), newP.begin(), newP.end());
                            retainedQ.insert(retainedQ.end(), newQ.begin(), newQ.end());
                            retainedConjunction.insert(
                                retainedConjunction.end(),
                                newConjunction.begin(),
                                newConjunction.end());
                            retainedRoot.insert(retainedRoot.end(), newRoot.begin(), newRoot.end());

                            int segmentCount = retainedSegmentation.size() - 1;
                            if (retainedP.size() != segmentCount ||
                                retainedQ.size() != segmentCount ||
                                retainedConjunction.size() != segmentCount ||
                                retainedRoot.size() != segmentCount ||
                                retainedSegmentation.back() != qCurrent)
                            {
                                throw logic_error("incremental formula caches are misaligned");
                            }

                            sequence.push_back(bitsetVerdict(retainedRoot.back()));

                            long long semanticCutoff = qCurrent - H_PHI;
                            int keep = 0;
                            while (keep < retainedRoot.size() &&
                                   retainedSegmentation[keep + 1] < semanticCutoff)
                            {
                                keep++;
                            }

                            if (keep > 0)
                            {
                                retainedSegmentation.erase(
                                    retainedSegmentation.begin(),
                                    retainedSegmentation.begin() + keep);
                                retainedP.erase(retainedP.begin(), retainedP.begin() + keep);
                                retainedQ.erase(retainedQ.begin(), retainedQ.begin() + keep);
                                retainedConjunction.erase(
                                    retainedConjunction.begin(),
                                    retainedConjunction.begin() + keep);
                                retainedRoot.erase(retainedRoot.begin(), retainedRoot.begin() + keep);
                            }

                            long long inputCutoff = max(0LL, qCurrent - eps);
                            for (int i = 0; i < nSignals; i++)
                            {
                                auto retainedInput = retainApproximateSignalFrom(
                                    retainedSignals[i],
                                    uncertainties[i],
                                    inputCutoff,
                                    predecessorLowerBounds[i]);
                                retainedSignals[i] = move(retainedInput.signal);
                                predecessorLowerBounds[i] =
                                    retainedInput.predecessorLowerBound;
                            }

                            dPrevious = dCurrent;
                            qPrevious = qCurrent;
                        }

                        return sequence;
                    };

                    chrono::duration<double, milli> total(0);
                    string seq;
                    for (int rep = 0; rep < options.repetitions; rep++)
                    {
                        auto tBeg = chrono::high_resolution_clock::now();
                        string current = runMonitor();
                        auto tEnd = chrono::high_resolution_clock::now();
                        total += tEnd - tBeg;

                        if (rep == 0)
                        {
                            seq = current;
                        }
                        else if (current != seq)
                        {
                            throw logic_error("incremental repetitions produced different verdict sequences");
                        }
                    }

                    string w = to_string(d / 1000) + " " +
                               to_string(eps / 1000) + " " +
                               to_string(CHUNK_SZ / 1000) + " " +
                               to_string(c) + " " +
                               to_string((total.count() / options.repetitions) / 1000) + " " +
                               seq;
                    results << w << endl;
                    if (!options.quiet)
                    {
                        cout << w << endl;
                    }
                }
            }
        }
    }

    results.close();
}

} // namespace

int main(int argc, char **argv)
{
    try
    {
        run(parseOptions(argc, argv));
    }
    catch (const exception &error)
    {
        cerr << "error: " << error.what() << '\n';
        return 1;
    }
    return 0;
}
