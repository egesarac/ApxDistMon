#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

#include "monitoring.hpp"

namespace fs = std::filesystem;

namespace
{

struct Formula
{
    std::string id;
    std::string description;
};

const std::vector<Formula> kFormulas{
    {"phi1", "always (p and q)"},
    {"phi2", "always (p or q)"},
    {"phi3", "p until q"},
    {"phi4", "always (p implies eventually q)"},
    {"phi5", "always (p implies eventually_[0,1) q)"},
    {"phi6", "always (p implies eventually_[0,2) q)"},
};

struct Options
{
    std::string formula = "all";
    std::optional<int> duration;
    std::optional<int> epsilon;
    std::optional<int> sample;
    int repetitions = 10;
    fs::path dataDirectory = "data/signals";
    fs::path outputDirectory = "results/offline/approximate/random";
    bool quiet = false;
    bool list = false;
};

void printUsage(const char *program)
{
    std::cout
        << "Usage: " << program << " [options]\n"
        << "  --formula <phi1..phi6|all>  Formula to evaluate (default: all)\n"
        << "  --duration <4|8|16|32>    Restrict trace duration in seconds\n"
        << "  --epsilon <1|2|4|8>        Restrict clock skew in seconds\n"
        << "  --sample <0..99>           Restrict the random sample\n"
        << "  --repeat <count>            Timed repetitions (default: 10)\n"
        << "  --data-dir <path>           Signal dataset directory\n"
        << "  --output-dir <path>         Generated result directory\n"
        << "  --quiet                     Do not echo result rows\n"
        << "  --list                      List formulas and exit\n"
        << "  --help                      Show this message\n";
}

int parseInteger(const std::string &value, const std::string &option)
{
    std::size_t parsed = 0;
    int result = 0;
    try
    {
        result = std::stoi(value, &parsed);
    }
    catch (const std::exception &)
    {
        throw std::runtime_error(option + " expects an integer: " + value);
    }

    if (parsed != value.size())
    {
        throw std::runtime_error(option + " expects an integer: " + value);
    }
    return result;
}

Options parseOptions(int argc, char **argv)
{
    Options options;

    for (int i = 1; i < argc; ++i)
    {
        const std::string argument = argv[i];
        auto valueFor = [&](const std::string &option) -> std::string
        {
            if (i + 1 >= argc)
            {
                throw std::runtime_error("missing value for " + option);
            }
            return argv[++i];
        };

        if (argument == "--formula")
        {
            options.formula = valueFor(argument);
        }
        else if (argument == "--duration")
        {
            options.duration = parseInteger(valueFor(argument), argument);
        }
        else if (argument == "--epsilon")
        {
            options.epsilon = parseInteger(valueFor(argument), argument);
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
        else if (argument == "--list")
        {
            options.list = true;
        }
        else if (argument == "--help")
        {
            printUsage(argv[0]);
            std::exit(0);
        }
        else
        {
            throw std::runtime_error("unknown option: " + argument);
        }
    }

    if (options.repetitions < 1)
    {
        throw std::runtime_error("--repeat must be positive");
    }
    if (options.sample && (*options.sample < 0 || *options.sample > 99))
    {
        throw std::runtime_error("--sample must be between 0 and 99");
    }

    return options;
}

std::vector<Formula> selectFormulas(const std::string &id)
{
    if (id == "all")
    {
        return kFormulas;
    }

    for (const Formula &formula : kFormulas)
    {
        if (formula.id == id)
        {
            return {formula};
        }
    }

    throw std::runtime_error("unknown formula: " + id);
}

std::vector<vector<bitset<SIZE>>> evaluateFormula(
    const std::string &formula,
    const vector<vector<vector<bitset<SIZE>>>> &atomicPropositions,
    const vector<long long> &segmentation)
{
    const auto &p = atomicPropositions[0];
    const auto &q = atomicPropositions[1];

    if (formula == "phi1")
    {
        return bitsetAlways(bitsetConjunction(p, q));
    }
    if (formula == "phi2")
    {
        return bitsetAlways(bitsetDisjunction(p, q));
    }
    if (formula == "phi3")
    {
        return bitsetUnboundedUntil(p, q, segmentation, 0, true);
    }
    if (formula == "phi4")
    {
        return bitsetAlways(
            bitsetDisjunction(bitsetNegation(p), bitsetEventually(q)));
    }
    if (formula == "phi5")
    {
        return bitsetAlways(
            bitsetNegation(
                bitsetConjunction(
                    p,
                    bitsetBoundedAlways(
                        bitsetNegation(q), segmentation, 0, 1000, true, false))));
    }
    if (formula == "phi6")
    {
        return bitsetAlways(
            bitsetNegation(
                bitsetConjunction(
                    p,
                    bitsetBoundedAlways(
                        bitsetNegation(q), segmentation, 0, 2000, true, false))));
    }

    throw std::runtime_error("unsupported formula: " + formula);
}

int verdict(const vector<vector<bitset<SIZE>>> &value)
{
    if (value[0][0].any() && value[0][1].any())
    {
        return 2;
    }
    if (value[0][1].any())
    {
        return 1;
    }
    return 0;
}

std::vector<int> selectValues(
    const std::vector<int> &defaults,
    const std::optional<int> &selected,
    const std::string &option)
{
    if (!selected)
    {
        return defaults;
    }
    if (std::find(defaults.begin(), defaults.end(), *selected) == defaults.end())
    {
        throw std::runtime_error("unsupported value for " + option + ": " + std::to_string(*selected));
    }
    return {*selected};
}

void runFormula(const Formula &formula, const Options &options)
{
    const std::vector<int> durations = selectValues({4, 8, 16, 32}, options.duration, "--duration");
    const std::vector<int> epsilons = selectValues({1, 2, 4, 8}, options.epsilon, "--epsilon");
    const int firstSample = options.sample.value_or(0);
    const int lastSample = options.sample.value_or(99);

    fs::create_directories(options.outputDirectory);
    const fs::path outputPath = options.outputDirectory / (formula.id + ".txt");
    std::ofstream results(outputPath);
    if (!results)
    {
        throw std::runtime_error("cannot open output file: " + outputPath.string());
    }

    for (const int duration : durations)
    {
        const long long durationMs = static_cast<long long>(duration) * 1000;

        for (const int epsilon : epsilons)
        {
            if (epsilon > duration)
            {
                continue;
            }
            const long long epsilonMs = static_cast<long long>(epsilon) * 1000;

            for (int sample = firstSample; sample <= lastSample; ++sample)
            {
                vector<vector<pair<long long, double>>> realSignals(2);
                realSignals[0] = getData(
                    (options.dataDirectory / (std::to_string(duration) + "_" + std::to_string(sample) + ".txt")).string(),
                    duration);
                realSignals[1] = getData(
                    (options.dataDirectory / (std::to_string(duration) + "_" + std::to_string(sample + 100) + ".txt")).string(),
                    duration);

                if (realSignals[0].size() != static_cast<std::size_t>(duration) ||
                    realSignals[1].size() != static_cast<std::size_t>(duration))
                {
                    throw std::runtime_error("missing or incomplete signal data for sample " + std::to_string(sample));
                }

                evenMask = generateBitmask(0);
                oddMask = generateBitmask(1);

                int segmentCount = 0;
                vector<vector<bitset<SIZE>>> result;
                const auto start = std::chrono::high_resolution_clock::now();

                for (int repetition = 0; repetition < options.repetitions; ++repetition)
                {
                    const auto signals = convertSignalsToBool(realSignals, duration);
                    const auto uncertainties = computeUncertaintyIntervals(signals, epsilonMs, durationMs);
                    const auto segmentation = computeCanonicalSegmentation(signals, uncertainties, durationMs);
                    segmentCount = static_cast<int>(segmentation.size()) - 1;
                    const auto valueExpressions = computeValueExpressions(signals, uncertainties, segmentation);
                    const auto atomicPropositions = convertSignalsToAtomicPropositions(valueExpressions, 0.0);
                    result = evaluateFormula(formula.id, atomicPropositions, segmentation);
                }

                const auto end = std::chrono::high_resolution_clock::now();
                const std::chrono::duration<double> elapsed = end - start;
                const double seconds = elapsed.count() / options.repetitions;

                const std::string row =
                    std::to_string(duration) + " " +
                    std::to_string(epsilon) + " 0 " +
                    std::to_string(sample) + " " +
                    std::to_string(segmentCount) + " " +
                    std::to_string(seconds) + " " +
                    std::to_string(verdict(result));

                results << row << '\n';
                if (!options.quiet)
                {
                    std::cout << formula.id << " " << row << '\n';
                }
            }
        }
    }
}

} // namespace

int main(int argc, char **argv)
{
    try
    {
        const Options options = parseOptions(argc, argv);
        if (options.list)
        {
            for (const Formula &formula : kFormulas)
            {
                std::cout << formula.id << ": " << formula.description << '\n';
            }
            return 0;
        }

        for (const Formula &formula : selectFormulas(options.formula))
        {
            runFormula(formula, options);
        }
    }
    catch (const std::exception &error)
    {
        std::cerr << "error: " << error.what() << '\n';
        return 1;
    }

    return 0;
}
