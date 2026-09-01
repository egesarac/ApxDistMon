#pragma once

#include <cassert>
#include <functional>
#include <queue>
#include <unordered_map>
#include <unordered_set>

#include "monitoring.hpp"

using ScalarSignal = vector<pair<long long, double>>;
using ScalarSignals = vector<ScalarSignal>;
using Point = array<double, 3>;
using PointSignal = vector<pair<long long, Point>>;
using PointSignals = vector<PointSignal>;
using BooleanLanguage = vector<vector<bitset<SIZE>>>;

constexpr long long kDurationMs = 1000;
constexpr long long kSamplingPeriodMs = 50;
constexpr double kWaterTankThreshold = 10.0;
constexpr double kSeparationThreshold = 0.0;

template <typename Value>
inline void discardRightEndpointMarker(
    vector<pair<long long, Value>> &signal,
    long long duration)
{
    if (!signal.empty() && signal.back().first == duration)
    {
        signal.pop_back();
    }
}

enum class Variant
{
    Adm,
    AdmFine,
    AdmFineRelative,
    AdmCoarse,
    AdmCoarseRelative,
};

inline string variantName(Variant variant)
{
    switch (variant)
    {
        case Variant::Adm:
            return "adm";
        case Variant::AdmFine:
            return "adm-f";
        case Variant::AdmFineRelative:
            return "adm-fr";
        case Variant::AdmCoarse:
            return "adm-c";
        case Variant::AdmCoarseRelative:
            return "adm-cr";
    }
    throw logic_error("unknown approximate case-study variant");
}

inline string variantFileSuffix(Variant variant)
{
    string suffix = variantName(variant);
    replace(suffix.begin(), suffix.end(), '-', '_');
    return suffix;
}

inline Variant parseVariant(const string &name)
{
    if (name == "adm")
    {
        return Variant::Adm;
    }
    if (name == "adm-f")
    {
        return Variant::AdmFine;
    }
    if (name == "adm-fr")
    {
        return Variant::AdmFineRelative;
    }
    if (name == "adm-c")
    {
        return Variant::AdmCoarse;
    }
    if (name == "adm-cr")
    {
        return Variant::AdmCoarseRelative;
    }
    throw runtime_error("unknown approximate variant: " + name);
}

inline bool isRelative(Variant variant)
{
    return variant == Variant::AdmFineRelative ||
           variant == Variant::AdmCoarseRelative;
}

inline bool isFine(Variant variant)
{
    return variant == Variant::AdmFine ||
           variant == Variant::AdmFineRelative;
}

inline bool isCoarse(Variant variant)
{
    return variant == Variant::AdmCoarse ||
           variant == Variant::AdmCoarseRelative;
}

struct Evaluation
{
    vector<long long> segmentation;
    BooleanLanguage language;

    bool falsePossible() const
    {
        validate();
        return language[0][0].any();
    }

    bool truePossible() const
    {
        validate();
        return language[0][1].any();
    }

    int verdict() const
    {
        validate();
        const bool canBeFalse = language[0][0].any();
        const bool canBeTrue = language[0][1].any();
        if (canBeFalse && canBeTrue)
        {
            return 2;
        }
        if (canBeTrue)
        {
            return 1;
        }
        if (canBeFalse)
        {
            return 0;
        }
        throw logic_error("approximate monitor produced an empty verdict");
    }

private:
    void validate() const
    {
        if (language.empty() || language[0].size() != 2)
        {
            throw logic_error("approximate monitor produced no boundary language");
        }
    }
};

inline void validateAgentCount(size_t count)
{
    if (count < 2 || count > 4)
    {
        throw invalid_argument("case studies require two, three, or four agents");
    }
}

inline long long timestampMilliseconds(
    double seconds,
    const filesystem::path &path,
    size_t lineNumber)
{
    if (!isfinite(seconds) || seconds < 0.0)
    {
        throw runtime_error(
            "invalid timestamp in " + path.string() + ":" +
            to_string(lineNumber));
    }

    const double milliseconds = seconds * 1000.0;
    const long long rounded = llround(milliseconds);
    if (fabs(milliseconds - static_cast<double>(rounded)) > 1e-6)
    {
        throw runtime_error(
            "timestamps must have whole-millisecond precision in " +
            path.string() + ":" + to_string(lineNumber));
    }
    return rounded;
}

inline ScalarSignal loadScalarSignal(const filesystem::path &path)
{
    ifstream input(path);
    if (!input)
    {
        throw runtime_error("cannot open signal data: " + path.string());
    }

    ScalarSignal signal;
    string line;
    size_t lineNumber = 0;
    while (getline(input, line))
    {
        ++lineNumber;
        if (line.find_first_not_of(" \t\r") == string::npos)
        {
            continue;
        }

        istringstream row(line);
        double seconds = 0.0;
        double value = 0.0;
        string extra;
        if (!(row >> seconds >> value) || (row >> extra))
        {
            throw runtime_error(
                "expected <time value> in " + path.string() + ":" +
                to_string(lineNumber));
        }
        signal.push_back({
            timestampMilliseconds(seconds, path, lineNumber), value});
    }
    if (signal.empty())
    {
        throw runtime_error("signal data is empty: " + path.string());
    }
    return signal;
}

inline PointSignal loadPointSignal(const filesystem::path &path)
{
    ifstream input(path);
    if (!input)
    {
        throw runtime_error("cannot open signal data: " + path.string());
    }

    PointSignal signal;
    string line;
    size_t lineNumber = 0;
    while (getline(input, line))
    {
        ++lineNumber;
        if (line.find_first_not_of(" \t\r") == string::npos)
        {
            continue;
        }

        istringstream row(line);
        double seconds = 0.0;
        Point point{};
        string extra;
        if (!(row >> seconds >> point[0] >> point[1] >> point[2]) ||
            (row >> extra))
        {
            throw runtime_error(
                "expected <time x y z> in " + path.string() + ":" +
                to_string(lineNumber));
        }
        signal.push_back({
            timestampMilliseconds(seconds, path, lineNumber), point});
    }
    if (signal.empty())
    {
        throw runtime_error("signal data is empty: " + path.string());
    }
    return signal;
}

namespace caseStudyDetail
{

struct PreparedEdge
{
    long long lower;
    long long upper;
    double before;
    double after;
};

struct PreparedSignal
{
    const ScalarSignal *samples;
    vector<PreparedEdge> edges;
    bool exact;
};

struct PreparedCaseStudy
{
    vector<long long> segmentation;
    vector<PreparedSignal> signals;
};

struct SweepCursor
{
    size_t sample = 0;
    size_t firstEdge = 0;
    size_t pastLastEdge = 0;
};

struct SegmentView
{
    double stableValue;
    size_t firstEdge;
    size_t pastLastEdge;
};

struct Bounds
{
    double minimum;
    double maximum;
};

using ValueSupport = vector<double>;

enum class Arithmetic
{
    Sum,
    SquaredDifference,
};

inline double normalizeValue(double value)
{
    if (!isfinite(value))
    {
        throw overflow_error("non-finite value in numeric value expression");
    }
    return value == 0.0 ? 0.0 : value;
}

inline double applyArithmetic(double left, double right, Arithmetic arithmetic)
{
    if (arithmetic == Arithmetic::Sum)
    {
        return normalizeValue(left + right);
    }
    const double difference = left - right;
    return normalizeValue(difference * difference);
}

inline PreparedCaseStudy prepareCaseStudy(
    const ScalarSignals &signals,
    long long epsilon,
    long long duration,
    const set<int> &exactSignals)
{
    for (const ScalarSignal &signal : signals)
    {
        if (signal.empty() || signal.front().first != 0)
        {
            throw invalid_argument(
                "case-study signals must start at timestamp zero");
        }
    }

    const auto uncertainties = computeUncertaintyIntervals(
        signals, epsilon, duration, exactSignals);
    PreparedCaseStudy prepared;
    prepared.segmentation = computeCanonicalSegmentation(
        signals, uncertainties, duration, exactSignals);
    prepared.signals.reserve(signals.size());

    for (size_t signalIndex = 0; signalIndex < signals.size(); ++signalIndex)
    {
        const bool exact = exactSignals.count(
            static_cast<int>(signalIndex)) != 0;
        PreparedSignal signal{&signals[signalIndex], {}, exact};
        if (!exact)
        {
            signal.edges.reserve(signals[signalIndex].size());
            for (size_t edge = 1; edge < signals[signalIndex].size(); ++edge)
            {
                if (!isActualEdge(signals[signalIndex], edge))
                {
                    continue;
                }
                signal.edges.push_back({
                    uncertainties[signalIndex][edge][0],
                    uncertainties[signalIndex][edge][1],
                    normalizeValue(signals[signalIndex][edge - 1].second),
                    normalizeValue(signals[signalIndex][edge].second)});
            }
        }
        prepared.signals.push_back(move(signal));
    }
    return prepared;
}

inline SegmentView advanceSweep(
    const PreparedSignal &signal,
    SweepCursor &cursor,
    long long start,
    long long end)
{
    assert(start < end);
    const ScalarSignal &samples = *signal.samples;
    while (cursor.sample + 1 < samples.size() &&
           samples[cursor.sample + 1].first <= start)
    {
        ++cursor.sample;
    }

    if (!signal.exact)
    {
        while (cursor.firstEdge < signal.edges.size() &&
               signal.edges[cursor.firstEdge].upper <= start)
        {
            ++cursor.firstEdge;
        }
        cursor.pastLastEdge = max(
            cursor.pastLastEdge, cursor.firstEdge);
        while (cursor.pastLastEdge < signal.edges.size() &&
               signal.edges[cursor.pastLastEdge].lower < end)
        {
            ++cursor.pastLastEdge;
        }
    }

    return {
        normalizeValue(samples[cursor.sample].second),
        cursor.firstEdge,
        cursor.pastLastEdge};
}

inline Bounds segmentBounds(
    const PreparedSignal &signal,
    const SegmentView &segment)
{
    if (signal.exact || segment.firstEdge == segment.pastLastEdge)
    {
        return {segment.stableValue, segment.stableValue};
    }

    const PreparedEdge &first = signal.edges[segment.firstEdge];
    Bounds result{
        min(first.before, first.after),
        max(first.before, first.after)};
    for (size_t edge = segment.firstEdge + 1;
         edge < segment.pastLastEdge;
         ++edge)
    {
        const PreparedEdge &current = signal.edges[edge];
        result.minimum = min(result.minimum, min(current.before, current.after));
        result.maximum = max(result.maximum, max(current.before, current.after));
    }
    return result;
}

inline Bounds sumSegmentBounds(
    const PreparedCaseStudy &prepared,
    const vector<SegmentView> &segments)
{
    assert(!prepared.signals.empty());
    assert(prepared.signals.size() == segments.size());
    Bounds total = segmentBounds(prepared.signals[0], segments[0]);
    for (size_t signal = 1; signal < prepared.signals.size(); ++signal)
    {
        const Bounds next = segmentBounds(
            prepared.signals[signal], segments[signal]);
        total.minimum = applyArithmetic(
            total.minimum, next.minimum, Arithmetic::Sum);
        total.maximum = applyArithmetic(
            total.maximum, next.maximum, Arithmetic::Sum);
    }
    return total;
}

inline ValueSupport segmentSupport(
    const PreparedSignal &signal,
    const SegmentView &segment)
{
    if (signal.exact || segment.firstEdge == segment.pastLastEdge)
    {
        return {segment.stableValue};
    }

    ValueSupport support;
    support.reserve(2 * (segment.pastLastEdge - segment.firstEdge));
    for (size_t edge = segment.firstEdge;
         edge < segment.pastLastEdge;
         ++edge)
    {
        support.push_back(signal.edges[edge].before);
        support.push_back(signal.edges[edge].after);
    }
    sort(support.begin(), support.end());
    support.erase(unique(support.begin(), support.end()), support.end());
    return support;
}

inline bool isZeroSupport(const ValueSupport &support)
{
    return support.size() == 1 && support[0] == 0.0;
}

inline ValueSupport combineSupports(
    ValueSupport left,
    ValueSupport right,
    Arithmetic arithmetic)
{
    if (arithmetic == Arithmetic::Sum)
    {
        if (isZeroSupport(right))
        {
            return left;
        }
        if (isZeroSupport(left))
        {
            return right;
        }
    }

    ValueSupport result;
    result.reserve(left.size() * right.size());
    for (double leftValue : left)
    {
        for (double rightValue : right)
        {
            result.push_back(applyArithmetic(
                leftValue, rightValue, arithmetic));
        }
    }
    sort(result.begin(), result.end());
    result.erase(unique(result.begin(), result.end()), result.end());
    return result;
}

inline void includeBooleanValue(
    vector<bitset<SIZE>> &segment,
    double value,
    double threshold)
{
    segment[value > threshold ? 1 : 0].set(0);
}

inline BooleanLanguage evaluateWaterPropositions(
    const PreparedCaseStudy &prepared,
    Variant variant,
    double threshold)
{
    const size_t segmentCount = prepared.segmentation.size() - 1;
    BooleanLanguage propositions(
        segmentCount, vector<bitset<SIZE>>(2));
    vector<SweepCursor> cursors(prepared.signals.size());
    vector<SegmentView> segments(prepared.signals.size());

    for (size_t segment = 0; segment < segmentCount; ++segment)
    {
        const long long start = prepared.segmentation[segment];
        const long long end = prepared.segmentation[segment + 1];
        for (size_t signal = 0; signal < prepared.signals.size(); ++signal)
        {
            segments[signal] = advanceSweep(
                prepared.signals[signal], cursors[signal], start, end);
        }

        if (isCoarse(variant))
        {
            const Bounds total = sumSegmentBounds(prepared, segments);
            includeBooleanValue(
                propositions[segment], total.minimum, threshold);
            includeBooleanValue(
                propositions[segment], total.maximum, threshold);
            continue;
        }

        ValueSupport total = segmentSupport(
            prepared.signals[0], segments[0]);
        for (size_t signal = 1; signal < prepared.signals.size(); ++signal)
        {
            total = combineSupports(
                move(total),
                segmentSupport(prepared.signals[signal], segments[signal]),
                Arithmetic::Sum);
        }
        for (double value : total)
        {
            includeBooleanValue(propositions[segment], value, threshold);
        }
    }
    return propositions;
}

inline ScalarSignals flattenPointSignalsTyped(const PointSignals &signals)
{
    ScalarSignals flattened(signals.size() * 3);
    for (size_t agent = 0; agent < signals.size(); ++agent)
    {
        for (size_t axis = 0; axis < 3; ++axis)
        {
            flattened[3 * agent + axis].reserve(signals[agent].size());
        }
        for (const auto &[timestamp, point] : signals[agent])
        {
            for (size_t axis = 0; axis < 3; ++axis)
            {
                flattened[3 * agent + axis].push_back({timestamp, point[axis]});
            }
        }
    }
    return flattened;
}

inline vector<pair<size_t, size_t>> agentPairs(size_t agentCount)
{
    vector<pair<size_t, size_t>> pairs;
    pairs.reserve(agentCount * (agentCount - 1) / 2);
    for (size_t left = 0; left < agentCount; ++left)
    {
        for (size_t right = left + 1; right < agentCount; ++right)
        {
            pairs.push_back({left, right});
        }
    }
    return pairs;
}

inline vector<BooleanLanguage> evaluateFineSeparationPropositions(
    const PreparedCaseStudy &prepared,
    size_t agentCount,
    double threshold)
{
    const size_t segmentCount = prepared.segmentation.size() - 1;
    const auto pairs = agentPairs(agentCount);
    vector<BooleanLanguage> propositions(
        pairs.size(),
        BooleanLanguage(segmentCount, vector<bitset<SIZE>>(2)));
    vector<SweepCursor> cursors(prepared.signals.size());
    vector<SegmentView> segments(prepared.signals.size());

    for (size_t segment = 0; segment < segmentCount; ++segment)
    {
        const long long start = prepared.segmentation[segment];
        const long long end = prepared.segmentation[segment + 1];
        vector<ValueSupport> supports;
        supports.reserve(prepared.signals.size());
        for (size_t signal = 0; signal < prepared.signals.size(); ++signal)
        {
            segments[signal] = advanceSweep(
                prepared.signals[signal], cursors[signal], start, end);
            supports.push_back(segmentSupport(
                prepared.signals[signal], segments[signal]));
        }

        for (size_t pair = 0; pair < pairs.size(); ++pair)
        {
            const auto [left, right] = pairs[pair];
            array<ValueSupport, 3> squaredCoordinates;
            for (size_t axis = 0; axis < 3; ++axis)
            {
                squaredCoordinates[axis] = combineSupports(
                    supports[3 * left + axis],
                    supports[3 * right + axis],
                    Arithmetic::SquaredDifference);
            }
            ValueSupport distance = combineSupports(
                move(squaredCoordinates[0]),
                move(squaredCoordinates[1]),
                Arithmetic::Sum);
            distance = combineSupports(
                move(distance),
                move(squaredCoordinates[2]),
                Arithmetic::Sum);
            for (double value : distance)
            {
                includeBooleanValue(
                    propositions[pair][segment], value, threshold);
            }
        }
    }
    return propositions;
}

using NumericWord = vector<double>;
using WordLanguage = vector<NumericWord>;
using NodeId = uint32_t;

struct PrefixKey
{
    NodeId parent;
    double value;

    bool operator==(const PrefixKey &other) const
    {
        return parent == other.parent && value == other.value;
    }
};

struct PrefixKeyHash
{
    size_t operator()(const PrefixKey &key) const
    {
        const size_t parentHash = hash<NodeId>{}(key.parent);
        const size_t valueHash = hash<double>{}(key.value);
        return parentHash ^ (valueHash + 0x9e3779b9U +
                             (parentHash << 6) + (parentHash >> 2));
    }
};

struct PrefixNode
{
    NodeId parent;
    double value;
    size_t length;
};

class PrefixArena
{
public:
    PrefixArena()
        : nodes{{0, 0.0, 0}}
    {
        interned.reserve(1024);
    }

    NodeId append(NodeId parent, double value)
    {
        value = normalizeValue(value);
        if (parent != 0 && nodes[parent].value == value)
        {
            return parent;
        }

        const PrefixKey key{parent, value};
        const auto found = interned.find(key);
        if (found != interned.end())
        {
            return found->second;
        }

        const NodeId id = static_cast<NodeId>(nodes.size());
        nodes.push_back({parent, value, nodes[parent].length + 1});
        interned.emplace(key, id);
        return id;
    }

    NumericWord materialize(NodeId id) const
    {
        NumericWord word(nodes[id].length);
        for (size_t position = word.size(); position > 0; --position)
        {
            word[position - 1] = nodes[id].value;
            id = nodes[id].parent;
        }
        return word;
    }

private:
    vector<PrefixNode> nodes;
    unordered_map<PrefixKey, NodeId, PrefixKeyHash> interned;
};

inline NumericWord appendWord(
    NumericWord prefix,
    const NumericWord &suffix)
{
    prefix.reserve(prefix.size() + suffix.size());
    for (double value : suffix)
    {
        if (prefix.empty() || prefix.back() != value)
        {
            prefix.push_back(value);
        }
    }
    return prefix;
}

inline vector<NumericWord> edgeContributions(
    const PreparedEdge &edge,
    long long start,
    long long end)
{
    const NumericWord whole{edge.before, edge.after};
    if (start <= edge.lower && edge.upper <= end)
    {
        return {whole};
    }
    if (start <= edge.lower && end < edge.upper)
    {
        return {{}, {edge.before}, whole};
    }
    if (edge.lower < start && edge.upper <= end)
    {
        return {{}, {edge.after}, whole};
    }
    return {{}, {edge.before}, {edge.after}, whole};
}

inline WordLanguage segmentWords(
    const PreparedSignal &signal,
    const SegmentView &segment,
    long long start,
    long long end)
{
    if (signal.exact || segment.firstEdge == segment.pastLastEdge)
    {
        return {{segment.stableValue}};
    }

    WordLanguage language(1);
    for (size_t edge = segment.firstEdge;
         edge < segment.pastLastEdge;
         ++edge)
    {
        const auto contributions = edgeContributions(
            signal.edges[edge], start, end);
        WordLanguage next;
        next.reserve(language.size() * contributions.size());
        for (const NumericWord &prefix : language)
        {
            for (const NumericWord &suffix : contributions)
            {
                next.push_back(appendWord(prefix, suffix));
            }
        }
        sort(next.begin(), next.end());
        next.erase(unique(next.begin(), next.end()), next.end());
        language = move(next);
    }
    language.erase(
        remove_if(
            language.begin(),
            language.end(),
            [](const NumericWord &word) { return word.empty(); }),
        language.end());
    return language;
}

inline bool isZeroLanguage(const WordLanguage &language)
{
    return language.size() == 1 && language[0].size() == 1 &&
           language[0][0] == 0.0;
}

inline void appendGridPredecessors(
    vector<NodeId> &cell,
    const vector<NodeId> &predecessors,
    double value,
    PrefixArena &arena)
{
    cell.reserve(cell.size() + predecessors.size());
    for (NodeId prefix : predecessors)
    {
        cell.push_back(arena.append(prefix, value));
    }
}

inline void productWordPair(
    const NumericWord &left,
    const NumericWord &right,
    Arithmetic arithmetic,
    PrefixArena &arena,
    vector<NodeId> &terminals)
{
    assert(!left.empty() && !right.empty());
    vector<vector<NodeId>> previous(right.size());
    vector<vector<NodeId>> current(right.size());

    for (size_t leftIndex = 0; leftIndex < left.size(); ++leftIndex)
    {
        for (size_t rightIndex = 0; rightIndex < right.size(); ++rightIndex)
        {
            vector<NodeId> &cell = current[rightIndex];
            cell.clear();
            const double value = applyArithmetic(
                left[leftIndex], right[rightIndex], arithmetic);

            if (leftIndex == 0 && rightIndex == 0)
            {
                cell.push_back(arena.append(0, value));
            }
            else
            {
                if (leftIndex > 0)
                {
                    appendGridPredecessors(
                        cell, previous[rightIndex], value, arena);
                }
                if (rightIndex > 0)
                {
                    appendGridPredecessors(
                        cell, current[rightIndex - 1], value, arena);
                }
                if (leftIndex > 0 && rightIndex > 0)
                {
                    appendGridPredecessors(
                        cell, previous[rightIndex - 1], value, arena);
                }
                sort(cell.begin(), cell.end());
                cell.erase(unique(cell.begin(), cell.end()), cell.end());
            }
        }
        previous.swap(current);
    }

    const vector<NodeId> &result = previous.back();
    terminals.insert(terminals.end(), result.begin(), result.end());
}

inline WordLanguage combineWordLanguages(
    const WordLanguage &left,
    const WordLanguage &right,
    Arithmetic arithmetic)
{
    if (arithmetic == Arithmetic::Sum)
    {
        if (isZeroLanguage(right))
        {
            return left;
        }
        if (isZeroLanguage(left))
        {
            return right;
        }
    }

    PrefixArena arena;
    vector<NodeId> terminals;
    for (const NumericWord &leftWord : left)
    {
        for (const NumericWord &rightWord : right)
        {
            productWordPair(
                leftWord, rightWord, arithmetic, arena, terminals);
        }
    }
    sort(terminals.begin(), terminals.end());
    terminals.erase(unique(terminals.begin(), terminals.end()), terminals.end());

    WordLanguage result;
    result.reserve(terminals.size());
    for (NodeId terminal : terminals)
    {
        result.push_back(arena.materialize(terminal));
    }
    sort(result.begin(), result.end());
    result.erase(unique(result.begin(), result.end()), result.end());
    return result;
}

inline void includeBooleanWord(
    vector<bitset<SIZE>> &segment,
    const NumericWord &word,
    double threshold)
{
    assert(!word.empty());
    const bool initial = word.front() > threshold;
    bool previous = initial;
    size_t transitions = 0;
    for (size_t token = 1; token < word.size(); ++token)
    {
        const bool current = word[token] > threshold;
        if (current != previous)
        {
            ++transitions;
            previous = current;
        }
    }
    if (transitions >= SIZE)
    {
        throw length_error("atomic-proposition word exceeds bitset capacity");
    }
    segment[initial ? 1 : 0].set(transitions);
}

struct NumericDagNode
{
    double value;
    vector<NodeId> children;
    bool accepting;
};

struct NumericDag
{
    vector<NumericDagNode> nodes;
    vector<NodeId> roots;
};

inline NumericDag minimizeDag(const NumericDag &language);
inline NumericDag determinizeDag(const NumericDag &language);

struct RawTransition
{
    size_t target;
    double value;
    bool epsilon;
};

struct RawState
{
    vector<RawTransition> transitions;
};

struct RawLanguage
{
    vector<RawState> states;
    size_t finalState;
};

inline RawLanguage segmentRawLanguage(
    const PreparedSignal &signal,
    const SegmentView &segment,
    long long start,
    long long end)
{
    if (signal.exact || segment.firstEdge == segment.pastLastEdge)
    {
        RawLanguage result{{RawState{}, RawState{}}, 1};
        result.states[0].transitions.push_back({
            1, segment.stableValue, false});
        return result;
    }

    RawLanguage result{{RawState{}}, 0};
    size_t current = 0;
    for (size_t edge = segment.firstEdge;
         edge < segment.pastLastEdge;
         ++edge)
    {
        const auto alternatives = edgeContributions(
            signal.edges[edge], start, end);
        vector<size_t> middle(alternatives.size(), 0);
        for (size_t alternative = 0;
             alternative < alternatives.size();
             ++alternative)
        {
            assert(alternatives[alternative].size() <= 2);
            if (alternatives[alternative].size() == 2)
            {
                middle[alternative] = result.states.size();
                result.states.push_back({});
            }
        }
        const size_t next = result.states.size();
        result.states.push_back({});

        for (size_t alternative = 0;
             alternative < alternatives.size();
             ++alternative)
        {
            const NumericWord &word = alternatives[alternative];
            if (word.empty())
            {
                result.states[current].transitions.push_back({
                    next, 0.0, true});
            }
            else if (word.size() == 1)
            {
                result.states[current].transitions.push_back({
                    next, word[0], false});
            }
            else
            {
                result.states[current].transitions.push_back({
                    middle[alternative], word[0], false});
                result.states[middle[alternative]].transitions.push_back({
                    next, word[1], false});
            }
        }
        current = next;
    }
    result.finalState = current;
    return result;
}

inline NumericDag segmentDag(
    const PreparedSignal &signal,
    const SegmentView &segment,
    long long start,
    long long end)
{
    const RawLanguage rawLanguage = segmentRawLanguage(
        signal, segment, start, end);
    const vector<RawState> &raw = rawLanguage.states;
    const size_t current = rawLanguage.finalState;

    NumericDag result;
    vector<size_t> tokenTargets;
    vector<vector<NodeId>> directTokens(raw.size());
    for (size_t state = 0; state < raw.size(); ++state)
    {
        for (const RawTransition &transition : raw[state].transitions)
        {
            if (transition.epsilon)
            {
                continue;
            }
            const NodeId token = static_cast<NodeId>(result.nodes.size());
            result.nodes.push_back({
                normalizeValue(transition.value), {}, false});
            tokenTargets.push_back(transition.target);
            directTokens[state].push_back(token);
        }
    }

    vector<unsigned char> canFinish(raw.size(), 0);
    vector<vector<NodeId>> nextTokens(raw.size());
    canFinish[current] = 1;
    for (size_t state = raw.size(); state-- > 0;)
    {
        nextTokens[state] = directTokens[state];
        for (const RawTransition &transition : raw[state].transitions)
        {
            if (!transition.epsilon)
            {
                continue;
            }
            canFinish[state] = canFinish[state] || canFinish[transition.target];
            nextTokens[state].insert(
                nextTokens[state].end(),
                nextTokens[transition.target].begin(),
                nextTokens[transition.target].end());
        }
        sort(nextTokens[state].begin(), nextTokens[state].end());
        nextTokens[state].erase(
            unique(nextTokens[state].begin(), nextTokens[state].end()),
            nextTokens[state].end());
    }

    for (NodeId token = 0; token < result.nodes.size(); ++token)
    {
        const size_t target = tokenTargets[token];
        result.nodes[token].children = nextTokens[target];
        result.nodes[token].accepting = canFinish[target] != 0;
    }
    result.roots = move(nextTokens[0]);
    return determinizeDag(result);
}

inline bool isZeroDag(const NumericDag &language)
{
    return language.nodes.size() == 1 && language.roots.size() == 1 &&
           language.roots[0] == 0 && language.nodes[0].value == 0.0 &&
           language.nodes[0].accepting &&
           language.nodes[0].children.empty();
}

struct DagSignature
{
    double value;
    vector<NodeId> children;
    bool accepting;

    bool operator==(const DagSignature &other) const
    {
        return value == other.value && children == other.children &&
               accepting == other.accepting;
    }
};

struct DagSignatureHash
{
    size_t operator()(const DagSignature &signature) const
    {
        size_t result = hash<double>{}(signature.value);
        result ^= static_cast<size_t>(signature.accepting) +
                  0x9e3779b9U + (result << 6) + (result >> 2);
        for (NodeId child : signature.children)
        {
            result ^= hash<NodeId>{}(child) + 0x9e3779b9U +
                      (result << 6) + (result >> 2);
        }
        return result;
    }
};

inline NumericDag minimizeDag(const NumericDag &language)
{
    NumericDag result;
    unordered_map<DagSignature, NodeId, DagSignatureHash> interned;
    vector<NodeId> canonical(
        language.nodes.size(), numeric_limits<NodeId>::max());
    interned.reserve(language.nodes.size());
    result.nodes.reserve(language.nodes.size());

    function<NodeId(NodeId)> reduce = [&](NodeId node) -> NodeId
    {
        if (canonical[node] != numeric_limits<NodeId>::max())
        {
            return canonical[node];
        }

        vector<NodeId> children;
        children.reserve(language.nodes[node].children.size());
        for (NodeId child : language.nodes[node].children)
        {
            children.push_back(reduce(child));
        }
        sort(children.begin(), children.end());
        children.erase(unique(children.begin(), children.end()), children.end());

        DagSignature signature{
            language.nodes[node].value,
            move(children),
            language.nodes[node].accepting};
        const auto found = interned.find(signature);
        if (found != interned.end())
        {
            canonical[node] = found->second;
            return found->second;
        }

        const NodeId id = static_cast<NodeId>(result.nodes.size());
        result.nodes.push_back({
            signature.value,
            signature.children,
            signature.accepting});
        interned.emplace(move(signature), id);
        canonical[node] = id;
        return id;
    };

    result.roots.reserve(language.roots.size());
    for (NodeId root : language.roots)
    {
        result.roots.push_back(reduce(root));
    }
    sort(result.roots.begin(), result.roots.end());
    result.roots.erase(
        unique(result.roots.begin(), result.roots.end()),
        result.roots.end());
    return result;
}

struct NodeSubsetHash
{
    size_t operator()(const vector<NodeId> &nodes) const
    {
        size_t result = 0;
        for (NodeId node : nodes)
        {
            result ^= hash<NodeId>{}(node) + 0x9e3779b9U +
                      (result << 6) + (result >> 2);
        }
        return result;
    }
};

inline NumericDag determinizeDag(const NumericDag &language)
{
    NumericDag result;
    unordered_map<vector<NodeId>, NodeId, NodeSubsetHash> subsets;
    vector<uint32_t> marks(language.nodes.size(), 0);
    uint32_t generation = 0;

    auto sameValueClosure = [&](vector<NodeId> seeds, double value)
    {
        ++generation;
        vector<NodeId> closure;
        closure.reserve(seeds.size());
        for (NodeId seed : seeds)
        {
            if (marks[seed] != generation)
            {
                marks[seed] = generation;
                closure.push_back(seed);
            }
        }
        for (size_t index = 0; index < closure.size(); ++index)
        {
            for (NodeId child : language.nodes[closure[index]].children)
            {
                if (language.nodes[child].value == value &&
                    marks[child] != generation)
                {
                    marks[child] = generation;
                    closure.push_back(child);
                }
            }
        }
        sort(closure.begin(), closure.end());
        return closure;
    };

    function<NodeId(vector<NodeId>)> build =
        [&](vector<NodeId> states) -> NodeId
    {
        const auto found = subsets.find(states);
        if (found != subsets.end())
        {
            return found->second;
        }

        const NodeId id = static_cast<NodeId>(result.nodes.size());
        subsets.emplace(states, id);
        const double value = language.nodes[states[0]].value;
        bool accepting = false;
        vector<pair<double, NodeId>> outgoing;
        for (NodeId state : states)
        {
            accepting = accepting || language.nodes[state].accepting;
            for (NodeId child : language.nodes[state].children)
            {
                if (language.nodes[child].value != value)
                {
                    outgoing.push_back({language.nodes[child].value, child});
                }
            }
        }
        result.nodes.push_back({value, {}, accepting});

        sort(outgoing.begin(), outgoing.end());
        vector<NodeId> children;
        for (size_t first = 0; first < outgoing.size();)
        {
            size_t last = first + 1;
            vector<NodeId> seeds{outgoing[first].second};
            while (last < outgoing.size() &&
                   outgoing[last].first == outgoing[first].first)
            {
                seeds.push_back(outgoing[last].second);
                ++last;
            }
            sort(seeds.begin(), seeds.end());
            seeds.erase(unique(seeds.begin(), seeds.end()), seeds.end());
            children.push_back(build(sameValueClosure(
                move(seeds), outgoing[first].first)));
            first = last;
        }
        sort(children.begin(), children.end());
        children.erase(unique(children.begin(), children.end()), children.end());
        result.nodes[id].children = move(children);
        return id;
    };

    vector<pair<double, NodeId>> roots;
    roots.reserve(language.roots.size());
    for (NodeId root : language.roots)
    {
        roots.push_back({language.nodes[root].value, root});
    }
    sort(roots.begin(), roots.end());
    for (size_t first = 0; first < roots.size();)
    {
        size_t last = first + 1;
        vector<NodeId> seeds{roots[first].second};
        while (last < roots.size() && roots[last].first == roots[first].first)
        {
            seeds.push_back(roots[last].second);
            ++last;
        }
        sort(seeds.begin(), seeds.end());
        seeds.erase(unique(seeds.begin(), seeds.end()), seeds.end());
        result.roots.push_back(build(sameValueClosure(
            move(seeds), roots[first].first)));
        first = last;
    }
    sort(result.roots.begin(), result.roots.end());
    result.roots.erase(
        unique(result.roots.begin(), result.roots.end()),
        result.roots.end());
    return minimizeDag(result);
}

inline uint64_t productKey(NodeId left, NodeId right)
{
    return (static_cast<uint64_t>(left) << 32) |
           static_cast<uint64_t>(right);
}

inline NumericDag productDags(
    const NumericDag &left,
    const NumericDag &right,
    Arithmetic arithmetic)
{
    NumericDag result;
    result.nodes.reserve(left.nodes.size() + right.nodes.size());
    unordered_map<uint64_t, NodeId> products;
    products.reserve(left.nodes.size() + right.nodes.size());

    function<NodeId(NodeId, NodeId)> build =
        [&](NodeId leftNode, NodeId rightNode) -> NodeId
    {
        const uint64_t key = productKey(leftNode, rightNode);
        const auto found = products.find(key);
        if (found != products.end())
        {
            return found->second;
        }

        const NodeId id = static_cast<NodeId>(result.nodes.size());
        products.emplace(key, id);
        result.nodes.push_back({
            applyArithmetic(
                left.nodes[leftNode].value,
                right.nodes[rightNode].value,
                arithmetic),
            {},
            left.nodes[leftNode].accepting &&
                right.nodes[rightNode].accepting});

        vector<NodeId> children;
        const auto &leftChildren = left.nodes[leftNode].children;
        const auto &rightChildren = right.nodes[rightNode].children;
        children.reserve(
            leftChildren.size() + rightChildren.size() +
            leftChildren.size() * rightChildren.size());
        for (NodeId child : leftChildren)
        {
            children.push_back(build(child, rightNode));
        }
        for (NodeId child : rightChildren)
        {
            children.push_back(build(leftNode, child));
        }
        for (NodeId leftChild : leftChildren)
        {
            for (NodeId rightChild : rightChildren)
            {
                children.push_back(build(leftChild, rightChild));
            }
        }
        sort(children.begin(), children.end());
        children.erase(unique(children.begin(), children.end()), children.end());
        result.nodes[id].children = move(children);
        return id;
    };

    result.roots.reserve(left.roots.size() * right.roots.size());
    for (NodeId leftRoot : left.roots)
    {
        for (NodeId rightRoot : right.roots)
        {
            result.roots.push_back(build(leftRoot, rightRoot));
        }
    }
    sort(result.roots.begin(), result.roots.end());
    result.roots.erase(
        unique(result.roots.begin(), result.roots.end()),
        result.roots.end());
    return determinizeDag(result);
}

struct BooleanDagState
{
    NodeId node;
    size_t transitions;
    bool initial;

    bool operator==(const BooleanDagState &other) const
    {
        return node == other.node && transitions == other.transitions &&
               initial == other.initial;
    }
};

struct BooleanDagStateHash
{
    size_t operator()(const BooleanDagState &state) const
    {
        size_t result = hash<NodeId>{}(state.node);
        result ^= hash<size_t>{}(state.transitions) + 0x9e3779b9U +
                  (result << 6) + (result >> 2);
        result ^= static_cast<size_t>(state.initial);
        return result;
    }
};

inline void includeBooleanDag(
    vector<bitset<SIZE>> &segment,
    const NumericDag &language,
    double threshold)
{
    vector<BooleanDagState> pending;
    unordered_set<BooleanDagState, BooleanDagStateHash> visited;
    pending.reserve(language.roots.size());
    visited.reserve(language.nodes.size());
    for (NodeId root : language.roots)
    {
        const bool truth = language.nodes[root].value > threshold;
        const BooleanDagState state{root, 0, truth};
        if (visited.insert(state).second)
        {
            pending.push_back(state);
        }
    }

    while (!pending.empty())
    {
        const BooleanDagState state = pending.back();
        pending.pop_back();
        const NumericDagNode &node = language.nodes[state.node];
        const bool currentTruth = node.value > threshold;
        if (node.accepting)
        {
            segment[state.initial ? 1 : 0].set(state.transitions);
        }

        for (NodeId child : node.children)
        {
            const bool truth = language.nodes[child].value > threshold;
            const size_t transitions = state.transitions +
                static_cast<size_t>(truth != currentTruth);
            if (transitions >= SIZE)
            {
                throw length_error(
                    "atomic-proposition word exceeds bitset capacity");
            }
            const BooleanDagState next{
                child, transitions, state.initial};
            if (visited.insert(next).second)
            {
                pending.push_back(next);
            }
        }
    }
}

struct MultiBooleanState
{
    array<NodeId, 4> nodes{};
    size_t transitions;
    bool initial;

    bool operator==(const MultiBooleanState &other) const
    {
        return nodes == other.nodes && transitions == other.transitions &&
               initial == other.initial;
    }
};

struct MultiBooleanStateHash
{
    size_t operator()(const MultiBooleanState &state) const
    {
        size_t result = hash<size_t>{}(state.transitions);
        for (NodeId node : state.nodes)
        {
            result ^= hash<NodeId>{}(node) + 0x9e3779b9U +
                      (result << 6) + (result >> 2);
        }
        result ^= static_cast<size_t>(state.initial);
        return result;
    }
};

inline double sumNodeValues(
    const vector<const NumericDag *> &operands,
    const array<NodeId, 4> &nodes)
{
    double value = operands[0]->nodes[nodes[0]].value;
    for (size_t operand = 1; operand < operands.size(); ++operand)
    {
        value = applyArithmetic(
            value,
            operands[operand]->nodes[nodes[operand]].value,
            Arithmetic::Sum);
    }
    return value;
}

inline void includeBooleanSumProduct(
    vector<bitset<SIZE>> &segment,
    const vector<NumericDag> &languages,
    double threshold)
{
    vector<const NumericDag *> operands;
    operands.reserve(languages.size());
    for (const NumericDag &language : languages)
    {
        if (!isZeroDag(language))
        {
            operands.push_back(&language);
        }
    }
    assert(operands.size() <= 4);
    if (operands.empty())
    {
        includeBooleanValue(segment, 0.0, threshold);
        return;
    }
    if (operands.size() == 1)
    {
        includeBooleanDag(segment, *operands[0], threshold);
        return;
    }

    vector<MultiBooleanState> pending;
    unordered_set<MultiBooleanState, MultiBooleanStateHash> visited;
    array<NodeId, 4> selected{};

    auto enqueueInitial = [&]()
    {
        const bool truth = sumNodeValues(operands, selected) > threshold;
        const MultiBooleanState state{selected, 0, truth};
        if (visited.insert(state).second)
        {
            pending.push_back(state);
        }
    };

    function<void(size_t)> selectRoots = [&](size_t operand)
    {
        if (operand == operands.size())
        {
            enqueueInitial();
            return;
        }
        for (NodeId root : operands[operand]->roots)
        {
            selected[operand] = root;
            selectRoots(operand + 1);
        }
    };
    selectRoots(0);

    while (!pending.empty())
    {
        const MultiBooleanState state = pending.back();
        pending.pop_back();

        bool accepting = true;
        for (size_t operand = 0; operand < operands.size(); ++operand)
        {
            accepting = accepting &&
                operands[operand]->nodes[state.nodes[operand]].accepting;
        }
        if (accepting)
        {
            segment[state.initial ? 1 : 0].set(state.transitions);
        }

        const bool currentTruth =
            sumNodeValues(operands, state.nodes) > threshold;
        selected = state.nodes;
        function<void(size_t, bool)> selectSuccessors =
            [&](size_t operand, bool advanced)
        {
            if (operand == operands.size())
            {
                if (!advanced)
                {
                    return;
                }
                const bool truth =
                    sumNodeValues(operands, selected) > threshold;
                const size_t transitions = state.transitions +
                    static_cast<size_t>(truth != currentTruth);
                if (transitions >= SIZE)
                {
                    throw length_error(
                        "atomic-proposition word exceeds bitset capacity");
                }
                const MultiBooleanState next{
                    selected, transitions, state.initial};
                if (visited.insert(next).second)
                {
                    pending.push_back(next);
                }
                return;
            }

            selected[operand] = state.nodes[operand];
            selectSuccessors(operand + 1, advanced);
            const auto &children =
                operands[operand]->nodes[state.nodes[operand]].children;
            for (NodeId child : children)
            {
                selected[operand] = child;
                selectSuccessors(operand + 1, true);
            }
        };
        selectSuccessors(0, false);
    }
}

struct RawCursor
{
    size_t state;
    double value;

    bool operator==(const RawCursor &other) const
    {
        return state == other.state && value == other.value;
    }
};

struct SparseCursorTuple
{
    array<RawCursor, 4> cursors{};

    bool operator==(const SparseCursorTuple &other) const
    {
        return cursors == other.cursors;
    }
};

struct SparseCursorTupleHash
{
    size_t operator()(const SparseCursorTuple &tuple) const
    {
        size_t result = 0;
        for (const RawCursor &cursor : tuple.cursors)
        {
            result ^= hash<size_t>{}(cursor.state) + 0x9e3779b9U +
                      (result << 6) + (result >> 2);
            result ^= hash<double>{}(cursor.value) + 0x9e3779b9U +
                      (result << 6) + (result >> 2);
        }
        return result;
    }
};

struct SparseCountState
{
    SparseCursorTuple tuple;
    // One reachable transition-count set for each initial truth value.
    array<bitset<SIZE>, 2> counts;
};

inline vector<RawCursor> firstRawCursors(const RawLanguage &language)
{
    vector<RawCursor> result;
    vector<size_t> pending{0};
    vector<unsigned char> visited(language.states.size(), 0);
    visited[0] = 1;
    while (!pending.empty())
    {
        const size_t state = pending.back();
        pending.pop_back();
        for (const RawTransition &transition :
             language.states[state].transitions)
        {
            if (transition.epsilon)
            {
                if (!visited[transition.target])
                {
                    visited[transition.target] = 1;
                    pending.push_back(transition.target);
                }
            }
            else
            {
                result.push_back({
                    transition.target,
                    normalizeValue(transition.value)});
            }
        }
    }
    sort(
        result.begin(),
        result.end(),
        [](const RawCursor &left, const RawCursor &right)
        {
            return left.state < right.state ||
                   (left.state == right.state && left.value < right.value);
        });
    result.erase(unique(result.begin(), result.end()), result.end());
    return result;
}

inline double sumRawValues(
    const array<RawCursor, 4> &cursors,
    size_t operandCount)
{
    double value = cursors[0].value;
    for (size_t operand = 1; operand < operandCount; ++operand)
    {
        value = applyArithmetic(
            value, cursors[operand].value, Arithmetic::Sum);
    }
    return value;
}

inline void includeBooleanSparseSum(
    vector<bitset<SIZE>> &segment,
    const vector<RawLanguage> &operands,
    double threshold)
{
    assert(!operands.empty() && operands.size() <= 4);
    vector<vector<RawCursor>> firstChoices;
    firstChoices.reserve(operands.size());
    for (const RawLanguage &operand : operands)
    {
        firstChoices.push_back(firstRawCursors(operand));
    }

    vector<SparseCountState> states;
    unordered_map<SparseCursorTuple, size_t, SparseCursorTupleHash> stateIds;
    priority_queue<
        pair<size_t, size_t>,
        vector<pair<size_t, size_t>>,
        greater<pair<size_t, size_t>>> pending;
    array<RawCursor, 4> selected{};

    auto rank = [&](const SparseCursorTuple &tuple)
    {
        size_t result = 0;
        for (size_t operand = 0; operand < operands.size(); ++operand)
        {
            result += tuple.cursors[operand].state;
        }
        return result;
    };

    auto enqueue = [&](const SparseCursorTuple &tuple,
                       const array<bitset<SIZE>, 2> &counts)
    {
        auto [position, inserted] = stateIds.emplace(tuple, states.size());
        if (inserted)
        {
            states.push_back({tuple, {}});
            pending.push({rank(tuple), position->second});
        }
        SparseCountState &state = states[position->second];
        state.counts[0] |= counts[0];
        state.counts[1] |= counts[1];
    };

    function<void(size_t)> selectFirst = [&](size_t operand)
    {
        if (operand == operands.size())
        {
            const bool truth =
                sumRawValues(selected, operands.size()) > threshold;
            array<bitset<SIZE>, 2> counts{};
            counts[truth ? 1 : 0].set(0);
            enqueue({selected}, counts);
            return;
        }
        for (const RawCursor &choice : firstChoices[operand])
        {
            selected[operand] = choice;
            selectFirst(operand + 1);
        }
    };
    selectFirst(0);

    while (!pending.empty())
    {
        const auto [currentRank, stateId] = pending.top();
        pending.pop();
        const SparseCursorTuple tuple = states[stateId].tuple;
        const array<bitset<SIZE>, 2> counts = states[stateId].counts;

        bool accepting = true;
        for (size_t operand = 0; operand < operands.size(); ++operand)
        {
            accepting = accepting &&
                tuple.cursors[operand].state == operands[operand].finalState;
        }
        if (accepting)
        {
            segment[0] |= counts[0];
            segment[1] |= counts[1];
        }

        const bool currentTruth =
            sumRawValues(tuple.cursors, operands.size()) > threshold;

        for (size_t operand = 0; operand < operands.size(); ++operand)
        {
            const auto &transitions = operands[operand]
                .states[tuple.cursors[operand].state].transitions;
            for (const RawTransition &transition : transitions)
            {
                if (!transition.epsilon)
                {
                    continue;
                }
                SparseCursorTuple next = tuple;
                next.cursors[operand].state = transition.target;
                assert(rank(next) > currentRank);
                enqueue(next, counts);
            }
        }

        selected = tuple.cursors;
        function<void(size_t, bool)> selectLabels =
            [&](size_t operand, bool advanced)
        {
            if (operand == operands.size())
            {
                if (!advanced)
                {
                    return;
                }
                const bool truth =
                    sumRawValues(selected, operands.size()) > threshold;
                array<bitset<SIZE>, 2> nextCounts = counts;
                if (truth != currentTruth)
                {
                    if (counts[0].test(SIZE - 1) ||
                        counts[1].test(SIZE - 1))
                    {
                        throw length_error(
                            "atomic-proposition word exceeds bitset capacity");
                    }
                    nextCounts[0] <<= 1;
                    nextCounts[1] <<= 1;
                }
                const SparseCursorTuple next{selected};
                assert(rank(next) > currentRank);
                enqueue(next, nextCounts);
                return;
            }

            selected[operand] = tuple.cursors[operand];
            selectLabels(operand + 1, advanced);
            const auto &transitions = operands[operand]
                .states[tuple.cursors[operand].state].transitions;
            for (const RawTransition &transition : transitions)
            {
                if (transition.epsilon)
                {
                    continue;
                }
                selected[operand] = {
                    transition.target,
                    normalizeValue(transition.value)};
                selectLabels(operand + 1, true);
            }
        };
        selectLabels(0, false);
    }
}

inline BooleanLanguage evaluateAdmWaterPropositions(
    const PreparedCaseStudy &prepared,
    double threshold)
{
    const size_t segmentCount = prepared.segmentation.size() - 1;
    BooleanLanguage propositions(
        segmentCount, vector<bitset<SIZE>>(2));
    vector<SweepCursor> cursors(prepared.signals.size());
    vector<SegmentView> segments(prepared.signals.size());

    for (size_t segment = 0; segment < segmentCount; ++segment)
    {
        const long long start = prepared.segmentation[segment];
        const long long end = prepared.segmentation[segment + 1];
        for (size_t signal = 0; signal < prepared.signals.size(); ++signal)
        {
            segments[signal] = advanceSweep(
                prepared.signals[signal], cursors[signal], start, end);
        }

        vector<NumericDag> sourceLanguages;
        sourceLanguages.reserve(prepared.signals.size());
        for (size_t signal = 0;
             signal < prepared.signals.size();
             ++signal)
        {
            sourceLanguages.push_back(segmentDag(
                prepared.signals[signal],
                segments[signal],
                start,
                end));
        }
        includeBooleanSumProduct(
            propositions[segment], sourceLanguages, threshold);
    }
    return propositions;
}

inline vector<BooleanLanguage> evaluateAdmSeparationPropositions(
    const PreparedCaseStudy &prepared,
    size_t agentCount,
    double threshold)
{
    const size_t segmentCount = prepared.segmentation.size() - 1;
    const auto pairs = agentPairs(agentCount);
    vector<BooleanLanguage> propositions(
        pairs.size(),
        BooleanLanguage(segmentCount, vector<bitset<SIZE>>(2)));
    vector<SweepCursor> cursors(prepared.signals.size());
    vector<SegmentView> segments(prepared.signals.size());

    for (size_t segment = 0; segment < segmentCount; ++segment)
    {
        const long long start = prepared.segmentation[segment];
        const long long end = prepared.segmentation[segment + 1];
        vector<NumericDag> sourceLanguages;
        sourceLanguages.reserve(prepared.signals.size());
        for (size_t signal = 0; signal < prepared.signals.size(); ++signal)
        {
            segments[signal] = advanceSweep(
                prepared.signals[signal], cursors[signal], start, end);
            sourceLanguages.push_back(segmentDag(
                prepared.signals[signal], segments[signal], start, end));
        }

        for (size_t pair = 0; pair < pairs.size(); ++pair)
        {
            const auto [left, right] = pairs[pair];
            vector<NumericDag> squaredCoordinates;
            squaredCoordinates.reserve(3);
            for (size_t axis = 0; axis < 3; ++axis)
            {
                squaredCoordinates.push_back(productDags(
                    sourceLanguages[3 * left + axis],
                    sourceLanguages[3 * right + axis],
                    Arithmetic::SquaredDifference));
            }
            includeBooleanSumProduct(
                propositions[pair][segment],
                squaredCoordinates,
                threshold);
        }
    }
    return propositions;
}

} // namespace caseStudyDetail

inline Evaluation evaluateWaterTanks(
    const ScalarSignals &signals,
    long long epsilon,
    long long duration = kDurationMs,
    Variant variant = Variant::Adm,
    double threshold = kWaterTankThreshold)
{
    validateAgentCount(signals.size());
    if (!isfinite(threshold))
    {
        throw invalid_argument("water-tank threshold must be finite");
    }
    const set<int> exactSignals = isRelative(variant)
        ? set<int>{0}
        : set<int>{};
    auto prepared = caseStudyDetail::prepareCaseStudy(
        signals, epsilon, duration, exactSignals);
    evenMask = generateBitmask(0);
    oddMask = generateBitmask(1);
    BooleanLanguage propositions = variant == Variant::Adm
        ? caseStudyDetail::evaluateAdmWaterPropositions(prepared, threshold)
        : caseStudyDetail::evaluateWaterPropositions(
            prepared, variant, threshold);
    BooleanLanguage result = bitsetAlways(propositions);
    return {move(prepared.segmentation), move(result)};
}

inline Evaluation evaluateMutualSeparation(
    const PointSignals &signals,
    long long epsilon,
    long long duration = kDurationMs,
    Variant variant = Variant::Adm,
    double threshold = kSeparationThreshold)
{
    validateAgentCount(signals.size());
    if (isCoarse(variant))
    {
        throw invalid_argument(
            "coarse variants are unsound for non-monotone squared separation");
    }
    if (!isfinite(threshold) || threshold < 0.0)
    {
        throw invalid_argument(
            "separation threshold must be finite and nonnegative");
    }
    ScalarSignals flattened = caseStudyDetail::flattenPointSignalsTyped(signals);
    const set<int> exactSignals = isRelative(variant)
        ? set<int>{0, 1, 2}
        : set<int>{};
    auto prepared = caseStudyDetail::prepareCaseStudy(
        flattened, epsilon, duration, exactSignals);
    evenMask = generateBitmask(0);
    oddMask = generateBitmask(1);

    auto propositions = variant == Variant::Adm
        ? caseStudyDetail::evaluateAdmSeparationPropositions(
            prepared, signals.size(), threshold * threshold)
        : caseStudyDetail::evaluateFineSeparationPropositions(
            prepared, signals.size(), threshold * threshold);
    BooleanLanguage result = bitsetAlways(propositions[0]);
    for (size_t pair = 1; pair < propositions.size(); ++pair)
    {
        result = bitsetConjunction(
            result, bitsetAlways(propositions[pair]));
    }
    return {move(prepared.segmentation), move(result)};
}
