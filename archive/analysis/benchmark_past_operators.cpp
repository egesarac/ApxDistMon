#include <algorithm>
#include <bitset>
#include <chrono>
#include <cstdlib>
#include <functional>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "monitoring.hpp"

namespace
{

using Language = std::vector<std::vector<std::bitset<SIZE>>>;
using Clock = std::chrono::steady_clock;

volatile std::size_t observedBits = 0;

Language makeLanguage(int segments, const std::string &pattern, int phase)
{
    Language language(
        static_cast<std::size_t>(segments),
        std::vector<std::bitset<SIZE>>(2));
    for (int i = 0; i < segments; i++)
    {
        const int value = (i + phase) % 2;
        if (pattern == "constant")
        {
            language[i][value][0] = true;
        }
        else if (pattern == "ambiguous")
        {
            language[i][0][0] = true;
            language[i][1][0] = true;
        }
        else if (pattern == "rich")
        {
            language[i][0][0] = true;
            language[i][1][0] = true;
            language[i][value][1] = true;
            language[i][value ^ 1][2] = true;
        }
        else
        {
            throw std::invalid_argument(
                "pattern must be constant, ambiguous, or rich");
        }
    }
    return language;
}

std::size_t countBits(const Language &language)
{
    std::size_t total = 0;
    for (const auto &segment : language)
    {
        total += segment[0].count();
        total += segment[1].count();
    }
    return total;
}

double medianMilliseconds(const std::function<Language()> &operation)
{
    observedBits += countBits(operation());
    std::vector<double> samples;
    samples.reserve(7);
    for (int repetition = 0; repetition < 7; repetition++)
    {
        const auto begin = Clock::now();
        const Language result = operation();
        const auto end = Clock::now();
        observedBits += countBits(result);
        samples.push_back(std::chrono::duration<double, std::milli>(
            end - begin).count());
    }
    std::sort(samples.begin(), samples.end());
    return samples[samples.size() / 2];
}

void report(
    const std::string &name,
    const std::function<Language()> &optimized,
    const std::function<Language()> &legacy)
{
    const Language optimizedResult = optimized();
    const Language legacyResult = legacy();
    if (optimizedResult != legacyResult)
    {
        throw std::runtime_error(
            "optimized and legacy results differ for " + name);
    }

    const double optimizedMilliseconds = medianMilliseconds(optimized);
    const double legacyMilliseconds = medianMilliseconds(legacy);
    std::cout << std::left << std::setw(24) << name
              << std::right << std::fixed << std::setprecision(3)
              << " optimized=" << optimizedMilliseconds << " ms"
              << " legacy=" << legacyMilliseconds << " ms"
              << " speedup="
              << legacyMilliseconds / optimizedMilliseconds << "x\n";
}

} // namespace

int main(int argc, char **argv)
{
    const int segments = argc > 1 ? std::atoi(argv[1]) : 256;
    const std::string pattern = argc > 2 ? argv[2] : "rich";
    if (segments <= 0)
    {
        throw std::invalid_argument("segment count must be positive");
    }

    evenMask = generateBitmask(0);
    oddMask = generateBitmask(1);
    allExceptFirstMask = evenMask | oddMask;
    allExceptFirstMask[0] = false;

    const Language lhs = makeLanguage(segments, pattern, 0);
    const Language rhs = makeLanguage(segments, pattern, 1);
    std::vector<long long> segmentation(
        static_cast<std::size_t>(segments) + 1);
    for (int i = 0; i <= segments; i++)
    {
        segmentation[i] = i;
    }

    std::cout << "segments=" << segments
              << " pattern=" << pattern << '\n';
    report("untimed historically", [&]
    {
        return bitsetAlwaysPast(lhs);
    }, [&]
    {
        return bitsetAlwaysPastLegacy(lhs);
    });
    report("untimed once", [&]
    {
        return bitsetEventuallyPast(rhs);
    }, [&]
    {
        return bitsetEventuallyPastLegacy(rhs);
    });
    report("untimed since", [&]
    {
        return bitsetSinceNonStrict(lhs, rhs);
    }, [&]
    {
        return bitsetSinceNonStrictLegacy(lhs, rhs);
    });
    report("timed historically", [&]
    {
        return bitsetTimedUnary(
            lhs, segmentation, 0, 8, false,
            true, false, true, true, 0, -1);
    }, [&]
    {
        return bitsetTimedUnaryLegacy(
            lhs, segmentation, 0, 8, false,
            true, false, true, true, 0, -1);
    });
    report("timed once", [&]
    {
        return bitsetTimedUnary(
            rhs, segmentation, 0, 8, false,
            true, false, false, true, 0, -1);
    }, [&]
    {
        return bitsetTimedUnaryLegacy(
            rhs, segmentation, 0, 8, false,
            true, false, false, true, 0, -1);
    });
    report("timed since", [&]
    {
        return bitsetTimedSince(
            lhs, rhs, segmentation, 0, 8, false,
            true, false, 0, -1);
    }, [&]
    {
        return bitsetTimedSinceLegacy(
            lhs, rhs, segmentation, 0, 8, false,
            true, false, 0, -1);
    });
    report("timed historically [1,inf)", [&]
    {
        return bitsetTimedUnary(
            lhs, segmentation, 1, 0, true,
            true, false, true, true, 0, -1);
    }, [&]
    {
        return bitsetTimedUnaryLegacy(
            lhs, segmentation, 1, 0, true,
            true, false, true, true, 0, -1);
    });
    report("timed once [1,inf)", [&]
    {
        return bitsetTimedUnary(
            rhs, segmentation, 1, 0, true,
            true, false, false, true, 0, -1);
    }, [&]
    {
        return bitsetTimedUnaryLegacy(
            rhs, segmentation, 1, 0, true,
            true, false, false, true, 0, -1);
    });
    report("timed since [1,inf)", [&]
    {
        return bitsetTimedSince(
            lhs, rhs, segmentation, 1, 0, true,
            true, false, 0, -1);
    }, [&]
    {
        return bitsetTimedSinceLegacy(
            lhs, rhs, segmentation, 1, 0, true,
            true, false, 0, -1);
    });
    std::cerr << "observed-bits=" << observedBits << '\n';
}
