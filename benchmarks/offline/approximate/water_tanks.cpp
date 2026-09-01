#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "monitoring_case_studies.hpp"

namespace fs = std::filesystem;

namespace
{

constexpr int kDefaultRepetitions = 10;

struct Options
{
    std::optional<int> agents;
    std::optional<int> epsilonMultiplier;
    std::string variant = "adm";
    int repetitions = kDefaultRepetitions;
    fs::path dataDirectory = "data/case_studies/water_tanks";
    fs::path outputDirectory = "results/offline/approximate/case_studies";
    bool quiet = false;
};

void printUsage(const char *program)
{
    std::cout
        << "Usage: " << program << " [options]\n"
        << "  --agents <2|3|4>             Restrict the tank count\n"
        << "  --epsilon-multiplier <n>     Restrict epsilon to n * 50 ms\n"
        << "  --variant <name|all>         adm, adm-f, adm-fr, adm-c, adm-cr\n"
        << "                               (default: adm)\n"
        << "  --repeat <count>             Timed repetitions (default: 10)\n"
        << "  --data-dir <path>            Water-tank dataset directory\n"
        << "  --output-dir <path>          Generated result directory\n"
        << "  --quiet                      Do not echo result rows\n"
        << "  --help                       Show this message\n";
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
    for (int index = 1; index < argc; ++index)
    {
        const std::string argument = argv[index];
        auto valueFor = [&](const std::string &option)
        {
            if (index + 1 >= argc)
            {
                throw std::runtime_error("missing value for " + option);
            }
            return std::string(argv[++index]);
        };

        if (argument == "--agents")
        {
            options.agents = parseInteger(valueFor(argument), argument);
        }
        else if (argument == "--epsilon-multiplier")
        {
            options.epsilonMultiplier = parseInteger(valueFor(argument), argument);
        }
        else if (argument == "--variant")
        {
            options.variant = valueFor(argument);
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
            std::exit(0);
        }
        else
        {
            throw std::runtime_error("unknown option: " + argument);
        }
    }

    if (options.agents && (*options.agents < 2 || *options.agents > 4))
    {
        throw std::runtime_error("--agents must be 2, 3, or 4");
    }
    if (options.epsilonMultiplier && *options.epsilonMultiplier < 1)
    {
        throw std::runtime_error("--epsilon-multiplier must be positive");
    }
    if (options.repetitions < 1)
    {
        throw std::runtime_error("--repeat must be positive");
    }
    return options;
}

std::vector<Variant> selectVariants(const std::string &selection)
{
    if (selection == "all")
    {
        return {Variant::Adm,
                Variant::AdmFine,
                Variant::AdmFineRelative,
                Variant::AdmCoarse,
                Variant::AdmCoarseRelative};
    }
    return {parseVariant(selection)};
}

std::vector<int> selectedValues(
    const std::vector<int> &defaults,
    const std::optional<int> &selected)
{
    return selected ? std::vector<int>{*selected} : defaults;
}

fs::path outputPath(const Options &options, Variant variant)
{
    if (variant == Variant::Adm)
    {
        return options.outputDirectory / "water_tanks.txt";
    }
    return options.outputDirectory /
           ("water_tanks_" + variantFileSuffix(variant) + ".txt");
}

void runVariant(
    Variant variant,
    const Options &options,
    const ScalarSignals &allSignals)
{
    const std::vector<int> agentCounts = selectedValues({2, 3, 4}, options.agents);
    const std::vector<int> epsilonMultipliers = selectedValues(
        {1, 2, 4, 8}, options.epsilonMultiplier);

    const fs::path path = outputPath(options, variant);
    std::ofstream output(path);
    if (!output)
    {
        throw std::runtime_error("cannot open output file: " + path.string());
    }

    for (const int agents : agentCounts)
    {
        const ScalarSignals signals(allSignals.begin(), allSignals.begin() + agents);
        for (const int multiplier : epsilonMultipliers)
        {
            const long long epsilon =
                static_cast<long long>(multiplier) * kSamplingPeriodMs;
            Evaluation evaluation;
            const auto started = std::chrono::steady_clock::now();
            for (int repetition = 0; repetition < options.repetitions; ++repetition)
            {
                evaluation = evaluateWaterTanks(
                    signals, epsilon, kDurationMs, variant);
            }
            const auto finished = std::chrono::steady_clock::now();
            const std::chrono::duration<double> elapsed = finished - started;
            const double seconds = elapsed.count() / options.repetitions;

            const std::string row =
                "1 " + std::to_string(agents) + " " +
                std::to_string(multiplier) + " " +
                std::to_string(evaluation.segmentation.size() - 1) + " " +
                formatDouble(seconds) + " " +
                std::to_string(evaluation.falsePossible()) + " " +
                std::to_string(evaluation.truePossible()) + " " +
                std::to_string(evaluation.verdict());
            output << row << '\n';
            if (!options.quiet)
            {
                std::cout << variantName(variant) << " " << row << '\n';
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
        const std::vector<Variant> variants = selectVariants(options.variant);
        const std::vector<int> agentCounts = selectedValues({2, 3, 4}, options.agents);
        const int maximumAgents = *std::max_element(agentCounts.begin(), agentCounts.end());

        ScalarSignals signals;
        signals.reserve(maximumAgents);
        for (int agent = 0; agent < maximumAgents; ++agent)
        {
            ScalarSignal signal = loadScalarSignal(
                options.dataDirectory / ("s5_tank_" + std::to_string(agent)));
            discardRightEndpointMarker(signal, kDurationMs);
            signals.push_back(std::move(signal));
        }
        fs::create_directories(options.outputDirectory);
        for (const Variant variant : variants)
        {
            runVariant(variant, options, signals);
        }
    }
    catch (const std::exception &error)
    {
        std::cerr << "error: " << error.what() << '\n';
        return 1;
    }
    return 0;
}
