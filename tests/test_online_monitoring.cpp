#include <cassert>
#include <iostream>
#include <set>
#include <random>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>
#include "monitoring.hpp"
using namespace std;

#ifndef RANDOM_TESTS
#define RANDOM_TESTS 2000
#endif
#ifndef PAST_CARRY_TESTS
#define PAST_CARRY_TESTS 500
#endif
#ifndef BOUNDED_RANGE_TESTS
#define BOUNDED_RANGE_TESTS 300
#endif
#ifndef BOUNDED_RESPONSE_TESTS
#define BOUNDED_RESPONSE_TESTS 100
#endif

using Language = vector<bitset<SIZE>>;
using SegmentedLanguage = vector<Language>;
using Signal = vector<pair<long long, double>>;

struct ResponseSnapshot
{
    long long q;
    vector<long long> segmentation;
    SegmentedLanguage p;
    SegmentedLanguage qLanguage;
    SegmentedLanguage once;
    SegmentedLanguage implication;
    SegmentedLanguage root;
    char verdict;
};

vector<vector<pair<long long, double>>> observedPrefix(
    const vector<vector<pair<long long, double>>> &signals,
    long long right)
{
    vector<vector<pair<long long, double>>> prefix(signals.size());
    for (int i = 0; i < signals.size(); i++)
    {
        for (const auto &edge : signals[i])
        {
            if (edge.first >= right)
            {
                break;
            }
            prefix[i].push_back(edge);
        }
        assert(!prefix[i].empty());
    }
    return prefix;
}

vector<ResponseSnapshot> runNaiveResponse(
    const vector<vector<pair<long long, double>>> &completeSignals,
    long long d,
    long long eps,
    long long chunk)
{
    vector<ResponseSnapshot> snapshots;
    vector<long long> finalizationPoints;
    long long dPrevious = 0;

    while (dPrevious < d)
    {
        long long dCurrent = min(d, dPrevious + chunk);
        long long qCurrent = max(0LL, dCurrent - eps);
        finalizationPoints.push_back(qCurrent);

        if (qCurrent == 0)
        {
            snapshots.push_back({qCurrent, {}, {}, {}, {}, {}, {}, '2'});
            dPrevious = dCurrent;
            continue;
        }

        auto signals = observedPrefix(completeSignals, dCurrent);
        const long long supportHorizon =
            onlineUncertaintyHorizon(dCurrent, eps);
        auto uncertainties = computeUncertaintyIntervals(
            signals, eps, supportHorizon);
        auto canonical = computeCanonicalSegmentation(
            signals, uncertainties, supportHorizon);
        auto segmentation = refineSegmentation(canonical, finalizationPoints, 0, qCurrent);
        auto valExprs = computeValueExpressions(signals, uncertainties, segmentation);
        auto aps = convertSignalsToAtomicPropositions(valExprs, 0.0);
        auto once = bitsetEventuallyPast(aps[1]);
        auto implication = bitsetDisjunction(bitsetNegation(aps[0]), once);
        auto root = bitsetAlwaysPast(implication);

        snapshots.push_back({qCurrent,
                             segmentation,
                             aps[0],
                             aps[1],
                             once,
                             implication,
                             root,
                             bitsetVerdict(root.back())});
        dPrevious = dCurrent;
    }

    return snapshots;
}

void assertSegmentLanguageEquals(const SegmentedLanguage &left,
                                 const SegmentedLanguage &right,
                                 const string &name)
{
    if (left != right)
    {
        cerr << "language mismatch in " << name << endl;
        cerr << "left segments: " << left.size() << ", right segments: " << right.size() << endl;
        for (int i = 0; i < min(left.size(), right.size()); i++)
        {
            if (left[i] != right[i])
            {
                cerr << "first mismatch at segment " << i << endl;
                cerr << "left: ";
                for (const auto &word : bitset2stringset(left[i])) cerr << word << " ";
                cerr << "\nright: ";
                for (const auto &word : bitset2stringset(right[i])) cerr << word << " ";
                cerr << endl;
                break;
            }
        }
        throw runtime_error("differential language mismatch");
    }
}

void compareNewRange(const ResponseSnapshot &naive,
                     long long qPrevious,
                     const vector<long long> &newSegmentation,
                     const SegmentedLanguage &newP,
                     const SegmentedLanguage &newQ,
                     const SegmentedLanguage &newOnce,
                     const SegmentedLanguage &newImplication,
                     const SegmentedLanguage &newRoot)
{
    auto it = find(naive.segmentation.begin(), naive.segmentation.end(), qPrevious);
    assert(it != naive.segmentation.end());
    int start = it - naive.segmentation.begin();

    vector<long long> expectedSegmentation(
        naive.segmentation.begin() + start,
        naive.segmentation.end());
    assert(expectedSegmentation == newSegmentation);

    auto suffix = [start](const SegmentedLanguage &language)
    {
        return SegmentedLanguage(language.begin() + start, language.end());
    };

    assertSegmentLanguageEquals(newP, suffix(naive.p), "new p");
    assertSegmentLanguageEquals(newQ, suffix(naive.qLanguage), "new q");
    assertSegmentLanguageEquals(newOnce, suffix(naive.once), "new once");
    assertSegmentLanguageEquals(newImplication, suffix(naive.implication), "new implication");
    assertSegmentLanguageEquals(newRoot, suffix(naive.root), "new root");
}

void compareRetained(const ResponseSnapshot &naive,
                     const vector<long long> &retainedSegmentation,
                     const SegmentedLanguage &retainedP,
                     const SegmentedLanguage &retainedQ,
                     const SegmentedLanguage &retainedOnce,
                     const SegmentedLanguage &retainedImplication,
                     const SegmentedLanguage &retainedRoot)
{
    assert(retainedSegmentation.size() == retainedRoot.size() + 1);

    for (int i = 0; i < retainedRoot.size(); i++)
    {
        auto it = find(naive.segmentation.begin(), naive.segmentation.end(), retainedSegmentation[i]);
        assert(it != naive.segmentation.end());
        int index = it - naive.segmentation.begin();
        assert(index + 1 < naive.segmentation.size());
        assert(naive.segmentation[index + 1] == retainedSegmentation[i + 1]);
        assert(retainedP[i] == naive.p[index]);
        assert(retainedQ[i] == naive.qLanguage[index]);
        assert(retainedOnce[i] == naive.once[index]);
        assert(retainedImplication[i] == naive.implication[index]);
        assert(retainedRoot[i] == naive.root[index]);
    }
}

void runIncrementalResponseAndCompare(
    const vector<vector<pair<long long, double>>> &completeSignals,
    long long d,
    long long eps,
    long long chunk,
    const vector<ResponseSnapshot> &naive)
{
    const int nSignals = completeSignals.size();
    long long dPrevious = 0;
    long long qPrevious = 0;
    int portion = 0;

    vector<vector<pair<long long, double>>> retainedSignals(nSignals);
    vector<optional<long long>> predecessorLowerBounds(nSignals);
    vector<int> nextEdge(nSignals, 1);
    for (int i = 0; i < nSignals; i++)
    {
        retainedSignals[i].push_back(completeSignals[i][0]);
    }

    vector<long long> retainedSegmentation{0};
    SegmentedLanguage retainedP;
    SegmentedLanguage retainedQ;
    SegmentedLanguage retainedOnce;
    SegmentedLanguage retainedImplication;
    SegmentedLanguage retainedRoot;
    set<tuple<string, long long, long long>> evaluated;

    while (dPrevious < d)
    {
        long long dCurrent = min(d, dPrevious + chunk);
        for (int i = 0; i < nSignals; i++)
        {
            while (nextEdge[i] < completeSignals[i].size() &&
                   completeSignals[i][nextEdge[i]].first < dCurrent)
            {
                const auto &edge = completeSignals[i][nextEdge[i]];
                retainedSignals[i].push_back(edge);
                nextEdge[i]++;
            }
        }

        long long qCurrent = max(0LL, dCurrent - eps);
        assert(qCurrent == naive[portion].q);

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
            char verdict = qCurrent == 0 ? '2' : bitsetVerdict(retainedRoot.back());
            assert(verdict == naive[portion].verdict);
            long long cutoff = max(0LL, qCurrent - eps);
            for (int i = 0; i < nSignals; i++)
            {
                auto retainedInput = retainApproximateSignalFrom(
                    retainedSignals[i],
                    uncertainties[i],
                    cutoff,
                    predecessorLowerBounds[i]);
                retainedSignals[i] = move(retainedInput.signal);
                predecessorLowerBounds[i] =
                    retainedInput.predecessorLowerBound;
            }
            dPrevious = dCurrent;
            portion++;
            continue;
        }

        auto canonical = computeCanonicalSegmentation(
            retainedSignals, uncertainties, supportHorizon);
        auto newSegmentation = refineSegmentation(canonical, {}, qPrevious, qCurrent);
        auto valExprs = computeValueExpressions(retainedSignals, uncertainties, newSegmentation);
        auto aps = convertSignalsToAtomicPropositions(valExprs, 0.0);
        auto newP = aps[0];
        auto newQ = aps[1];

        bool once0 = true;
        bool once1 = false;
        if (!retainedOnce.empty())
        {
            tie(once0, once1) = bitsetLastBits(retainedOnce.back());
        }
        auto newOnce = bitsetEventuallyPast(newQ, 0, newQ.size(), once0, once1);
        auto newImplication = bitsetDisjunction(bitsetNegation(newP), newOnce);

        bool root0 = false;
        bool root1 = true;
        if (!retainedRoot.empty())
        {
            tie(root0, root1) = bitsetLastBits(retainedRoot.back());
        }
        auto newRoot = bitsetAlwaysPast(
            newImplication,
            0,
            newImplication.size(),
            root0,
            root1);

        try
        {
            compareNewRange(
                naive[portion],
                qPrevious,
                newSegmentation,
                newP,
                newQ,
                newOnce,
                newImplication,
                newRoot);
        }
        catch (...)
        {
            cerr << "portion " << portion << ", qPrevious=" << qPrevious
                 << ", qCurrent=" << qCurrent << "\nnew segmentation: ";
            for (auto endpoint : newSegmentation) cerr << endpoint << " ";
            cerr << "\nnaive segmentation: ";
            for (auto endpoint : naive[portion].segmentation) cerr << endpoint << " ";
            cerr << "\nretained signals:\n";
            for (int component = 0; component < retainedSignals.size(); component++)
            {
                cerr << component << ": ";
                for (const auto &edge : retainedSignals[component])
                    cerr << "(" << edge.first << "," << edge.second << ") ";
                cerr << "\nuncertainties: ";
                for (const auto &interval : uncertainties[component])
                {
                    cerr << "[";
                    for (auto value : interval) cerr << value << ",";
                    cerr << "] ";
                }
                cerr << endl;
            }
            throw;
        }

        for (int i = 0; i < newRoot.size(); i++)
        {
            long long left = newSegmentation[i];
            long long right = newSegmentation[i + 1];
            for (const string &name : {"p", "q", "once", "implication", "root"})
            {
                assert(evaluated.insert(make_tuple(name, left, right)).second);
            }
        }

        retainedSegmentation.insert(
            retainedSegmentation.end(),
            newSegmentation.begin() + 1,
            newSegmentation.end());
        retainedP.insert(retainedP.end(), newP.begin(), newP.end());
        retainedQ.insert(retainedQ.end(), newQ.begin(), newQ.end());
        retainedOnce.insert(retainedOnce.end(), newOnce.begin(), newOnce.end());
        retainedImplication.insert(
            retainedImplication.end(),
            newImplication.begin(),
            newImplication.end());
        retainedRoot.insert(retainedRoot.end(), newRoot.begin(), newRoot.end());

        compareRetained(
            naive[portion],
            retainedSegmentation,
            retainedP,
            retainedQ,
            retainedOnce,
            retainedImplication,
            retainedRoot);
        assert(retainedRoot.back() == naive[portion].root.back());
        assert(bitsetVerdict(retainedRoot.back()) == naive[portion].verdict);

        long long semanticCutoff = qCurrent;
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
            retainedOnce.erase(retainedOnce.begin(), retainedOnce.begin() + keep);
            retainedImplication.erase(
                retainedImplication.begin(),
                retainedImplication.begin() + keep);
            retainedRoot.erase(retainedRoot.begin(), retainedRoot.begin() + keep);
        }

        compareRetained(
            naive[portion],
            retainedSegmentation,
            retainedP,
            retainedQ,
            retainedOnce,
            retainedImplication,
            retainedRoot);
        assert(retainedRoot.size() == 1);

        long long inputCutoff = max(0LL, qCurrent - eps);
        for (int component = 0; component < nSignals; component++)
        {
            auto retainedInput = retainApproximateSignalFrom(
                retainedSignals[component],
                uncertainties[component],
                inputCutoff,
                predecessorLowerBounds[component]);
            retainedSignals[component] = move(retainedInput.signal);
            predecessorLowerBounds[component] =
                retainedInput.predecessorLowerBound;
            assert(retainedSignals[component].front().first == inputCutoff);
            for (int i = 1; i < retainedSignals[component].size(); i++)
            {
                assert(retainedSignals[component][i].first > inputCutoff);
            }
        }

        dPrevious = dCurrent;
        qPrevious = qCurrent;
        portion++;
    }

    assert(portion == naive.size());
}

Language constantLanguage(bool value)
{
    Language language(2);
    language[value][0] = true;
    return language;
}

void testHelpers()
{
    Language language(2);
    language[0][0] = true;
    assert(bitsetLastBits(language) == make_pair(true, false));
    assert(bitsetVerdict(language) == '0');

    language.assign(2, bitset<SIZE>());
    language[0][1] = true;
    assert(bitsetLastBits(language) == make_pair(false, true));

    language.assign(2, bitset<SIZE>());
    language[1][1] = true;
    assert(bitsetLastBits(language) == make_pair(true, false));

    language.assign(2, bitset<SIZE>());
    language[1][0] = true;
    assert(bitsetLastBits(language) == make_pair(false, true));

    language[0][0] = true;
    assert(bitsetVerdict(language) == '2');

    vector<long long> canonical{0, 3, 7, 10};
    vector<long long> finalization{4, 7, 9, 12};
    assert((refineSegmentation(canonical, finalization, 0, 10) ==
            vector<long long>{0, 3, 4, 7, 9, 10}));
    assert((refineSegmentation(canonical, finalization, 4, 4) ==
            vector<long long>{4}));

    Signal signal{{0, 0}, {3, 1}, {5, 0}, {8, 1}};
    assert((retainSignalFrom(signal, 5) == Signal{{5, 0}, {8, 1}}));
    assert((retainSignalFrom(signal, 4) == Signal{{4, 1}, {5, 0}, {8, 1}}));
    assert((retainSignalFrom(signal, 0) == signal));
}

void testOnlineUncertaintyState()
{
    assert(onlineUncertaintyHorizon(48, 4) == 52);
    bool overflowed = false;
    try
    {
        onlineUncertaintyHorizon(numeric_limits<long long>::max(), 1);
    }
    catch (const overflow_error &)
    {
        overflowed = true;
    }
    assert(overflowed);

    Signal stable{{0, 0}, {40, 1}, {41, 0}, {42, 1},
                  {43, 0}, {45, 1}, {46, 0}};
    vector<vector<pair<long long, double>>> stableSignals{stable};
    const auto at48 = computeUncertaintyIntervals(
        stableSignals, 4, onlineUncertaintyHorizon(48, 4));
    const auto at51 = computeUncertaintyIntervals(
        stableSignals, 4, onlineUncertaintyHorizon(51, 4));
    assert(at48 == at51);

    Signal dense{{0, 0}, {1, 1}, {2, 0}, {3, 1}, {4, 0}, {5, 1}};
    vector<vector<pair<long long, double>>> denseSignals{dense};
    const auto fullRegions = computeUncertaintyIntervals(denseSignals, 2, 8);
    auto retained = retainApproximateSignalFrom(
        dense, fullRegions[0], 1);
    assert(retained.predecessorLowerBound == optional<long long>(0));
    vector<vector<pair<long long, double>>> retainedSignals{
        retained.signal};
    const vector<optional<long long>> predecessorLowerBounds{
        retained.predecessorLowerBound};
    const auto retainedRegions = computeUncertaintyIntervals(
        retainedSignals,
        2,
        8,
        {},
        predecessorLowerBounds);

    assert(retained.signal.size() + 1 == dense.size());
    for (size_t edge = 1; edge < retained.signal.size(); ++edge)
    {
        assert(retainedRegions[0][edge] == fullRegions[0][edge + 1]);
    }
}

void testSinceRangeCompilation()
{
    SegmentedLanguage lhs{constantLanguage(true), constantLanguage(false), constantLanguage(true)};
    SegmentedLanguage rhs{constantLanguage(false), constantLanguage(true), constantLanguage(false)};

    auto full = bitsetSinceNonStrict(lhs, rhs);
    auto first = bitsetSinceNonStrict(
        SegmentedLanguage{lhs[0], lhs[1]},
        SegmentedLanguage{rhs[0], rhs[1]});
    auto [since0, since1] = bitsetLastBits(first.back());
    auto [lhs0, lhs1] = bitsetLastBits(lhs[1]);
    bool carry0 = since0 || lhs0;
    bool carry1 = since1 && lhs1;
    auto last = bitsetSinceNonStrict(
        SegmentedLanguage{lhs[2]},
        SegmentedLanguage{rhs[2]},
        0,
        1,
        carry0,
        carry1);

    if (last[0] != full[2])
    {
        cerr << "since full[2]: ";
        for (const auto &word : bitset2stringset(full[2])) cerr << word << " ";
        cerr << "\nfirst.back: ";
        for (const auto &word : bitset2stringset(first.back())) cerr << word << " ";
        cerr << "\nlhs previous: ";
        for (const auto &word : bitset2stringset(lhs[1])) cerr << word << " ";
        cerr << "\ncarry " << carry0 << " " << carry1;
        cerr << "\nlast: ";
        for (const auto &word : bitset2stringset(last[0])) cerr << word << " ";
        cerr << endl;
    }
    assert(last[0] == full[2]);

    auto ranged = bitsetSinceNonStrict(lhs, rhs, 1, 3, true, false);
    assert(ranged.size() == lhs.size());
}

void testDifferentialResponse()
{
    vector<vector<pair<long long, double>>> signals(2);
    signals[0] = {{0, 0}, {1500, 1}, {4200, 0}, {7000, 1}, {9500, 0}};
    signals[1] = {{0, 0}, {500, 1}, {2800, 0}, {6200, 1}, {8400, 0}};

    for (const auto &[d, eps, chunk] : vector<tuple<long long, long long, long long>>{
             {12000, 2000, 3000},
             {12000, 4000, 3000},
             {12000, 2000, 5000},
             {12000, 1000, 2500}})
    {
        auto naive = runNaiveResponse(signals, d, eps, chunk);
        runIncrementalResponseAndCompare(signals, d, eps, chunk, naive);
    }
}

vector<ResponseSnapshot> runNaiveBoundedResponse(
    const vector<vector<pair<long long, double>>> &completeSignals,
    long long d,
    long long eps,
    long long chunk,
    long long upper)
{
    vector<ResponseSnapshot> snapshots;
    vector<long long> finalizationPoints;
    long long dPrevious = 0;

    while (dPrevious < d)
    {
        long long dCurrent = min(d, dPrevious + chunk);
        long long qCurrent = max(0LL, dCurrent - eps);
        finalizationPoints.push_back(qCurrent);

        if (qCurrent == 0)
        {
            snapshots.push_back({qCurrent, {}, {}, {}, {}, {}, {}, '2'});
            dPrevious = dCurrent;
            continue;
        }

        auto signals = observedPrefix(completeSignals, dCurrent);
        const long long supportHorizon =
            onlineUncertaintyHorizon(dCurrent, eps);
        auto uncertainties = computeUncertaintyIntervals(
            signals, eps, supportHorizon);
        auto canonical = computeCanonicalSegmentation(
            signals, uncertainties, supportHorizon);
        auto segmentation = refineSegmentation(canonical, finalizationPoints, 0, qCurrent);
        auto valExprs = computeValueExpressions(signals, uncertainties, segmentation);
        auto aps = convertSignalsToAtomicPropositions(valExprs, 0.0);
        auto once = bitsetBoundedEventuallyPast(
            aps[1],
            segmentation,
            0,
            upper,
            true,
            false);
        auto implication = bitsetDisjunction(bitsetNegation(aps[0]), once);
        auto root = bitsetAlwaysPast(implication);

        snapshots.push_back({qCurrent,
                             segmentation,
                             aps[0],
                             aps[1],
                             once,
                             implication,
                             root,
                             bitsetVerdict(root.back())});
        dPrevious = dCurrent;
    }

    return snapshots;
}

void runIncrementalBoundedResponseAndCompare(
    const vector<vector<pair<long long, double>>> &completeSignals,
    long long d,
    long long eps,
    long long chunk,
    long long upper,
    const vector<ResponseSnapshot> &naive)
{
    const int nSignals = completeSignals.size();
    long long dPrevious = 0;
    long long qPrevious = 0;
    int portion = 0;

    vector<vector<pair<long long, double>>> retainedSignals(nSignals);
    vector<optional<long long>> predecessorLowerBounds(nSignals);
    vector<int> nextEdge(nSignals, 1);
    for (int i = 0; i < nSignals; i++)
    {
        retainedSignals[i].push_back(completeSignals[i][0]);
    }

    vector<long long> retainedSegmentation{0};
    SegmentedLanguage retainedP;
    SegmentedLanguage retainedQ;
    SegmentedLanguage retainedOnce;
    SegmentedLanguage retainedImplication;
    SegmentedLanguage retainedRoot;
    set<tuple<string, long long, long long>> evaluated;

    while (dPrevious < d)
    {
        long long dCurrent = min(d, dPrevious + chunk);
        for (int i = 0; i < nSignals; i++)
        {
            while (nextEdge[i] < completeSignals[i].size() &&
                   completeSignals[i][nextEdge[i]].first < dCurrent)
            {
                const auto &edge = completeSignals[i][nextEdge[i]];
                retainedSignals[i].push_back(edge);
                nextEdge[i]++;
            }
        }

        long long qCurrent = max(0LL, dCurrent - eps);
        assert(qCurrent == naive[portion].q);

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
            char verdict = qCurrent == 0 ? '2' : bitsetVerdict(retainedRoot.back());
            assert(verdict == naive[portion].verdict);
            long long cutoff = max(0LL, qCurrent - eps);
            for (int i = 0; i < nSignals; i++)
            {
                auto retainedInput = retainApproximateSignalFrom(
                    retainedSignals[i],
                    uncertainties[i],
                    cutoff,
                    predecessorLowerBounds[i]);
                retainedSignals[i] = move(retainedInput.signal);
                predecessorLowerBounds[i] =
                    retainedInput.predecessorLowerBound;
            }
            dPrevious = dCurrent;
            portion++;
            continue;
        }

        auto canonical = computeCanonicalSegmentation(
            retainedSignals, uncertainties, supportHorizon);
        auto newSegmentation = refineSegmentation(canonical, {}, qPrevious, qCurrent);
        auto valExprs = computeValueExpressions(retainedSignals, uncertainties, newSegmentation);
        auto aps = convertSignalsToAtomicPropositions(valExprs, 0.0);
        auto newP = aps[0];
        auto newQ = aps[1];

        int historyCount = retainedQ.size();
        auto combinedSegmentation = retainedSegmentation;
        combinedSegmentation.insert(
            combinedSegmentation.end(),
            newSegmentation.begin() + 1,
            newSegmentation.end());
        auto combinedQ = retainedQ;
        combinedQ.insert(combinedQ.end(), newQ.begin(), newQ.end());
        auto combinedOnce = bitsetBoundedEventuallyPast(
            combinedQ,
            combinedSegmentation,
            0,
            upper,
            true,
            false,
            historyCount,
            combinedQ.size());
        SegmentedLanguage newOnce(
            combinedOnce.begin() + historyCount,
            combinedOnce.end());

        auto newImplication = bitsetDisjunction(bitsetNegation(newP), newOnce);
        bool root0 = false;
        bool root1 = true;
        if (!retainedRoot.empty())
        {
            tie(root0, root1) = bitsetLastBits(retainedRoot.back());
        }
        auto newRoot = bitsetAlwaysPast(
            newImplication,
            0,
            newImplication.size(),
            root0,
            root1);

        compareNewRange(
            naive[portion],
            qPrevious,
            newSegmentation,
            newP,
            newQ,
            newOnce,
            newImplication,
            newRoot);

        for (int i = 0; i < newRoot.size(); i++)
        {
            long long left = newSegmentation[i];
            long long right = newSegmentation[i + 1];
            for (const string &name : {"p", "q", "bounded-once", "implication", "root"})
            {
                assert(evaluated.insert(make_tuple(name, left, right)).second);
            }
        }

        retainedSegmentation.insert(
            retainedSegmentation.end(),
            newSegmentation.begin() + 1,
            newSegmentation.end());
        retainedP.insert(retainedP.end(), newP.begin(), newP.end());
        retainedQ.insert(retainedQ.end(), newQ.begin(), newQ.end());
        retainedOnce.insert(retainedOnce.end(), newOnce.begin(), newOnce.end());
        retainedImplication.insert(
            retainedImplication.end(),
            newImplication.begin(),
            newImplication.end());
        retainedRoot.insert(retainedRoot.end(), newRoot.begin(), newRoot.end());

        compareRetained(
            naive[portion],
            retainedSegmentation,
            retainedP,
            retainedQ,
            retainedOnce,
            retainedImplication,
            retainedRoot);
        assert(retainedRoot.back() == naive[portion].root.back());
        assert(bitsetVerdict(retainedRoot.back()) == naive[portion].verdict);

        long long semanticCutoff = qCurrent - upper;
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
            retainedOnce.erase(retainedOnce.begin(), retainedOnce.begin() + keep);
            retainedImplication.erase(
                retainedImplication.begin(),
                retainedImplication.begin() + keep);
            retainedRoot.erase(retainedRoot.begin(), retainedRoot.begin() + keep);
        }

        compareRetained(
            naive[portion],
            retainedSegmentation,
            retainedP,
            retainedQ,
            retainedOnce,
            retainedImplication,
            retainedRoot);

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
        portion++;
    }

    assert(portion == naive.size());
}

void testBoundedResponseDifferential()
{
    vector<vector<pair<long long, double>>> signals(2);
    signals[0] = {{0, 0}, {3, 1}, {8, 0}, {13, 1}, {19, 0}, {24, 1}};
    signals[1] = {{0, 1}, {2, 0}, {5, 1}, {11, 0}, {16, 1}, {21, 0}};

    for (long long upper : {3LL, 6LL, 10LL})
    {
        auto naive = runNaiveBoundedResponse(signals, 30, 4, 7, upper);
        runIncrementalBoundedResponseAndCompare(
            signals,
            30,
            4,
            7,
            upper,
            naive);
    }

    mt19937 generator(99173);
    for (int test = 0; test < BOUNDED_RESPONSE_TESTS; test++)
    {
        long long d = 20 + generator() % 31;
        long long eps = 1 + generator() % 8;
        long long chunk = 1 + generator() % 12;
        long long upper = 1 + generator() % 10;
        vector<vector<pair<long long, double>>> randomSignals(2);
        for (int component = 0; component < 2; component++)
        {
            int value = generator() % 2;
            randomSignals[component].push_back(make_pair(0, value));
            for (long long time = 1; time < d; time++)
            {
                if (generator() % 6 == 0)
                {
                    value = 1 - value;
                    randomSignals[component].push_back(make_pair(time, value));
                }
            }
        }

        auto naive = runNaiveBoundedResponse(
            randomSignals,
            d,
            eps,
            chunk,
            upper);
        runIncrementalBoundedResponseAndCompare(
            randomSignals,
            d,
            eps,
            chunk,
            upper,
            naive);
    }
}

Language randomLanguage(mt19937 &generator)
{
    Language language(2);
    int words = 1 + generator() % 5;
    for (int i = 0; i < words; i++)
    {
        int first = generator() % 2;
        int length = 1 + generator() % 5;
        language[first][length - 1] = true;
    }
    return language;
}

void testPastOperatorCarries()
{
    mt19937 generator(1701);

    for (int test = 0; test < PAST_CARRY_TESTS; test++)
    {
        int count = 2 + generator() % 6;
        int split = 1 + generator() % (count - 1);
        SegmentedLanguage lhs(count);
        SegmentedLanguage rhs(count);
        for (int i = 0; i < count; i++)
        {
            lhs[i] = randomLanguage(generator);
            rhs[i] = randomLanguage(generator);
        }

        auto fullOnce = bitsetEventuallyPast(rhs);
        auto firstOnce = bitsetEventuallyPast(
            SegmentedLanguage(rhs.begin(), rhs.begin() + split));
        auto [once0, once1] = bitsetLastBits(firstOnce.back());
        auto secondOnce = bitsetEventuallyPast(
            SegmentedLanguage(rhs.begin() + split, rhs.end()),
            0,
            count - split,
            once0,
            once1);
        for (int i = split; i < count; i++)
        {
            assert(secondOnce[i - split] == fullOnce[i]);
        }

        auto fullHistorically = bitsetAlwaysPast(rhs);
        auto firstHistorically = bitsetAlwaysPast(
            SegmentedLanguage(rhs.begin(), rhs.begin() + split));
        auto [historically0, historically1] = bitsetLastBits(firstHistorically.back());
        auto secondHistorically = bitsetAlwaysPast(
            SegmentedLanguage(rhs.begin() + split, rhs.end()),
            0,
            count - split,
            historically0,
            historically1);
        for (int i = split; i < count; i++)
        {
            assert(secondHistorically[i - split] == fullHistorically[i]);
        }

        auto fullSince = bitsetSinceNonStrict(lhs, rhs);
        auto firstSince = bitsetSinceNonStrict(
            SegmentedLanguage(lhs.begin(), lhs.begin() + split),
            SegmentedLanguage(rhs.begin(), rhs.begin() + split));
        auto [since0, since1] = bitsetLastBits(firstSince.back());
        auto [lhs0, lhs1] = bitsetLastBits(lhs[split - 1]);
        bool carry0 = since0 || lhs0;
        bool carry1 = since1 && lhs1;
        auto secondSince = bitsetSinceNonStrict(
            SegmentedLanguage(lhs.begin() + split, lhs.end()),
            SegmentedLanguage(rhs.begin() + split, rhs.end()),
            0,
            count - split,
            carry0,
            carry1);
        for (int i = split; i < count; i++)
        {
            assert(secondSince[i - split] == fullSince[i]);
        }

        // Strict Since is intentionally excluded from the maintained monitor.
#if 0
        auto fullSinceStrict = bitsetSinceStrict(lhs, rhs);
        auto firstSinceStrict = bitsetSinceStrict(
            SegmentedLanguage(lhs.begin(), lhs.begin() + split),
            SegmentedLanguage(rhs.begin(), rhs.begin() + split));
        auto [sinceStrict0, sinceStrict1] = bitsetLastBits(firstSinceStrict.back());
        auto secondSinceStrict = bitsetSinceStrict(
            SegmentedLanguage(lhs.begin() + split, lhs.end()),
            SegmentedLanguage(rhs.begin() + split, rhs.end()),
            0,
            count - split,
            sinceStrict0,
            sinceStrict1);
        for (int i = split; i < count; i++)
        {
            assert(secondSinceStrict[i - split] == fullSinceStrict[i]);
        }
#endif
    }
}

void testBoundedPastRange()
{
    mt19937 generator(8119);

    for (int test = 0; test < BOUNDED_RANGE_TESTS; test++)
    {
        int count = 2 + generator() % 6;
        vector<long long> segmentation{0};
        SegmentedLanguage child(count);
        for (int i = 0; i < count; i++)
        {
            segmentation.push_back(segmentation.back() + 1 + generator() % 4);
            child[i] = randomLanguage(generator);
        }

        long long upper = 1 + generator() % 8;
        int start = 1 + generator() % (count - 1);
        auto full = bitsetBoundedEventuallyPast(
            child,
            segmentation,
            0,
            upper,
            true,
            false);
        auto ranged = bitsetBoundedEventuallyPast(
            child,
            segmentation,
            0,
            upper,
            true,
            false,
            start,
            count);
        for (int i = start; i < count; i++)
        {
            assert(ranged[i] == full[i]);
        }

        long long cutoff = segmentation[start] - upper;
        int keep = 0;
        while (keep < start && segmentation[keep + 1] < cutoff)
        {
            keep++;
        }

        vector<long long> retainedSegmentation(
            segmentation.begin() + keep,
            segmentation.end());
        SegmentedLanguage retainedChild(
            child.begin() + keep,
            child.end());
        int retainedStart = start - keep;
        auto retainedResult = bitsetBoundedEventuallyPast(
            retainedChild,
            retainedSegmentation,
            0,
            upper,
            true,
            false,
            retainedStart,
            retainedChild.size());
        for (int i = start; i < count; i++)
        {
            assert(retainedResult[i - keep] == full[i]);
        }
    }
}

void testRandomizedDifferentialResponse()
{
    mt19937 generator(20260819);

    for (int test = 0; test < RANDOM_TESTS; test++)
    {
        long long d = 20 + generator() % 41;
        long long eps = 1 + generator() % 10;
        long long chunk = 1 + generator() % 15;

        vector<vector<pair<long long, double>>> signals(2);
        for (int component = 0; component < 2; component++)
        {
            int value = generator() % 2;
            signals[component].push_back(make_pair(0, value));

            vector<long long> candidates;
            for (long long time = 1; time < d; time++)
            {
                if (generator() % 5 == 0)
                {
                    candidates.push_back(time);
                }
            }

            for (const auto &time : candidates)
            {
                value = 1 - value;
                signals[component].push_back(make_pair(time, value));
            }
        }

        try
        {
            auto naive = runNaiveResponse(signals, d, eps, chunk);
            runIncrementalResponseAndCompare(signals, d, eps, chunk, naive);
        }
        catch (const exception &error)
        {
            cerr << "randomized test " << test
                 << " failed with d=" << d
                 << ", eps=" << eps
                 << ", chunk=" << chunk
                 << ": " << error.what() << endl;
            for (int component = 0; component < 2; component++)
            {
                cerr << "signal " << component << ": ";
                for (const auto &edge : signals[component])
                {
                    cerr << "(" << edge.first << "," << edge.second << ") ";
                }
                cerr << endl;
            }
            throw;
        }
    }
}

int main()
{
    evenMask = generateBitmask(0);
    oddMask = generateBitmask(1);
    allExceptFirstMask = evenMask | oddMask;
    allExceptFirstMask[0] = 0;

    testHelpers();
    testOnlineUncertaintyState();
    testSinceRangeCompilation();
    testPastOperatorCarries();
    testBoundedPastRange();
    testDifferentialResponse();
    testBoundedResponseDifferential();
    testRandomizedDifferentialResponse();

    cout << "online monitoring tests passed" << endl;
    return 0;
}
