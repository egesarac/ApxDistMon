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
    fs::path outputDirectory = "results/online/naive";
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

    fs::create_directories(options.outputDirectory);
    const fs::path outputPath = options.outputDirectory / "onlineNaive_ec.txt";
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
                        vector<long long> finalizationPoints;
                        long long dPrevious = 0;

                        while (dPrevious < d)
                        {
                            long long dCurrent = min(d, dPrevious + CHUNK_SZ);
                            long long qCurrent = max(0LL, dCurrent - eps);
                            finalizationPoints.push_back(qCurrent);

                            auto signals = convertSignalsToBool(
                                signalsReal,
                                dCurrent / 1000);

                            if (qCurrent == 0)
                            {
                                sequence.push_back('2');
                                dPrevious = dCurrent;
                                continue;
                            }

                            const long long supportHorizon =
                                onlineUncertaintyHorizon(dCurrent, eps);
                            auto uncertainties = computeUncertaintyIntervals(
                                signals, eps, supportHorizon);
                            auto canonical = computeCanonicalSegmentation(
                                signals, uncertainties, supportHorizon);
                            auto segmentation = refineSegmentation(
                                canonical,
                                finalizationPoints,
                                0,
                                qCurrent);

                            if (segmentation.size() < 2 ||
                                segmentation.front() != 0 ||
                                segmentation.back() != qCurrent)
                            {
                                throw logic_error("invalid refined naive segmentation");
                            }

                            auto valExprs = computeValueExpressions(signals, uncertainties, segmentation);
                            auto aps = convertSignalsToAtomicPropositions(valExprs, 0.0);

                            auto conjunction = bitsetConjunction(aps[0], aps[1]);
                            auto root = bitsetEventuallyPast(conjunction);

                            int last = root.size() - 1;
                            if (segmentation[last + 1] != qCurrent)
                            {
                                throw logic_error("the final naive segment does not end at the watermark");
                            }
                            sequence.push_back(bitsetVerdict(root[last]));

                            dPrevious = dCurrent;
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
                            throw logic_error("naive repetitions produced different verdict sequences");
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
