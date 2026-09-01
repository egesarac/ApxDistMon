#include <cassert>
#include <functional>
#include <iostream>
#include <set>
#include <string>
#include <vector>

#include "monitoring.hpp"

using BoolLanguage = vector<bitset<SIZE>>;
using SegmentedLanguage = vector<vector<bitset<SIZE>>>;

void initializeMasks()
{
    evenMask = generateBitmask(0);
    oddMask = generateBitmask(1);
    allExceptFirstMask.set();
    allExceptFirstMask[0] = false;
}

string alternatingWord(const int &start, const int &length)
{
    string word;
    for (int i = 0; i < length; i++)
    {
        word += (char)('0' + (start ^ (i % 2)));
    }
    return word;
}

BoolLanguage makeLanguage(const set<string> &words)
{
    BoolLanguage out(2);
    for (const auto &word : words)
    {
        assert(!word.empty());
        for (std::size_t i = 1; i < word.size(); i++)
        {
            assert(word[i - 1] != word[i]);
        }
        out[word[0] - '0'][word.size() - 1] = true;
    }
    return out;
}

BoolLanguage makeLanguage(initializer_list<string> words)
{
    return makeLanguage(set<string>(words));
}

set<string> firstOracle(const set<string> &language)
{
    set<string> out;
    for (const auto &word : language)
    {
        out.insert(word.substr(0, 1));
    }
    return out;
}

set<string> prefixOracle(const set<string> &language)
{
    set<string> out;
    for (const auto &word : language)
    {
        for (std::size_t length = 1; length <= word.size(); length++)
        {
            out.insert(word.substr(0, length));
        }
    }
    return out;
}

set<string> suffixOracle(const set<string> &language)
{
    set<string> out;
    for (const auto &word : language)
    {
        for (std::size_t start = 0; start < word.size(); start++)
        {
            out.insert(word.substr(start));
        }
    }
    return out;
}

set<string> infixOracle(const set<string> &language)
{
    set<string> out;
    for (const auto &word : language)
    {
        for (std::size_t start = 0; start < word.size(); start++)
        {
            for (std::size_t length = 1; start + length <= word.size(); length++)
            {
                out.insert(word.substr(start, length));
            }
        }
    }
    return out;
}

string destutterWord(const string &word)
{
    string out;
    for (const auto &letter : word)
    {
        if (out.empty() || out.back() != letter)
        {
            out += letter;
        }
    }
    return out;
}

set<string> concatOracle(const set<string> &left, const set<string> &right)
{
    if (left.empty())
    {
        return right;
    }

    set<string> out;
    for (const auto &first : left)
    {
        for (const auto &second : right)
        {
            out.insert(destutterWord(first + second));
        }
    }
    return out;
}

template <typename Exception>
void expectException(const function<void()> &action)
{
    bool thrown = false;
    try
    {
        action();
    }
    catch (const Exception &)
    {
        thrown = true;
    }
    assert(thrown);
}

void testExistingLanguageHelpers()
{
    vector<string> words;
    for (int length = 1; length <= 6; length++)
    {
        for (int start = 0; start < 2; start++)
        {
            words.push_back(alternatingWord(start, length));
        }
    }

    for (unsigned int mask = 1; mask < (1U << words.size()); mask++)
    {
        set<string> language;
        for (std::size_t i = 0; i < words.size(); i++)
        {
            if ((mask & (1U << i)) != 0)
            {
                language.insert(words[i]);
            }
        }

        BoolLanguage encoded = makeLanguage(language);
        assert(bitset2stringset(segmentFirstBit(encoded)) == firstOracle(language));
        assert(bitset2stringset(segmentPrefix(encoded)) == prefixOracle(language));
        assert(bitset2stringset(segmentSuffix(encoded)) == suffixOracle(language));
        assert(bitset2stringset(segmentInfix(encoded)) == infixOracle(language));
    }

    vector<string> shortWords;
    for (int length = 1; length <= 4; length++)
    {
        for (int start = 0; start < 2; start++)
        {
            shortWords.push_back(alternatingWord(start, length));
        }
    }

    const unsigned int subsetCount = 1U << shortWords.size();
    for (unsigned int leftMask = 0; leftMask < subsetCount; leftMask++)
    {
        set<string> left;
        for (std::size_t i = 0; i < shortWords.size(); i++)
        {
            if ((leftMask & (1U << i)) != 0)
            {
                left.insert(shortWords[i]);
            }
        }

        for (unsigned int rightMask = 1; rightMask < subsetCount; rightMask++)
        {
            set<string> right;
            for (std::size_t i = 0; i < shortWords.size(); i++)
            {
                if ((rightMask & (1U << i)) != 0)
                {
                    right.insert(shortWords[i]);
                }
            }

            BoolLanguage actual = bitsetConcat(makeLanguage(left), makeLanguage(right));
            assert(bitset2stringset(actual) == concatOracle(left, right));
        }
    }
}

void testGeometry()
{
    TimedRange closed{0, 4, true, true};
    TimedRange halfOpen{1, 3, false, true};
    TimedRange intersection = timedIntersection(closed, halfOpen);

    assert(intersection.left == 1 && intersection.right == 3);
    assert(!intersection.leftClosed && intersection.rightClosed);
    assert(timedContains(intersection, 2));
    assert(!timedContains(intersection, 1));
    assert(timedContains(intersection, 3));

    TimedRange singleton{2, 2, true, true};
    assert(!timedEmpty(singleton));
    assert(timedContains(singleton, 2));

    TimedRange openSingleton{2, 2, true, false};
    assert(timedEmpty(openSingleton));

    TimedRange disjoint = timedIntersection(
        TimedRange{0, 1, true, true},
        TimedRange{2, 3, true, true});
    assert(timedEmpty(disjoint));
}
void testProfiles()
{
    vector<long long> oneSegment{0, 2};
    SegmentedLanguage value{makeLanguage({"010"})};

    TimedRange full{0, 4, true, false};
    assert(bitset2stringset(timedProfile(value, oneSegment, full)) == set<string>{"010"});

    TimedRange firstPoint{0, 0, true, true};
    assert(bitset2stringset(timedProfile(value, oneSegment, firstPoint)) == set<string>{"0"});

    TimedRange prefix{0, 2, true, true};
    assert(bitset2stringset(timedProfile(value, oneSegment, prefix)) == set<string>({"0", "01", "010"}));

    TimedRange suffix{2, 4, true, false};
    assert(bitset2stringset(timedProfile(value, oneSegment, suffix)) == set<string>({"0", "10", "010"}));

    TimedRange interior{1, 3, true, true};
    assert(bitset2stringset(timedProfile(value, oneSegment, interior)) == set<string>({"0", "1", "01", "10", "010"}));

    vector<long long> segmentation{0, 1, 3};
    SegmentedLanguage proposition{makeLanguage({"0"}), makeLanguage({"0", "01"})};
    TimedRange manuscriptHorizon{1, 3, true, true};
    assert(bitset2stringset(timedProfile(proposition, segmentation, manuscriptHorizon)) == set<string>({"0", "01"}));

    SegmentedLanguage boundaryValues{makeLanguage({"0"}), makeLanguage({"10"})};
    TimedRange closedBoundary{1, 2, true, true};
    TimedRange openBoundary{1, 2, true, false};
    assert(bitset2stringset(timedProfile(boundaryValues, segmentation, closedBoundary)) == set<string>{"01"});
    assert(bitset2stringset(timedProfile(boundaryValues, segmentation, openBoundary)) == set<string>{"0"});

    vector<long long> endSegmentation{0, 1};
    SegmentedLanguage endValue{makeLanguage({"0"})};
    TimedRange endsAtD{1, 2, true, true};
    assert(bitset2stringset(timedProfile(endValue, endSegmentation, timedIntersection(endsAtD, TimedRange{0, 2, true, false}))) == set<string>{"0"});
}

void compositions(const int &remaining, const int &parts, const int &position,
                  vector<int> &current, vector<vector<int>> &out)
{
    if (position + 1 == parts)
    {
        if (remaining >= 1)
        {
            current[position] = remaining;
            out.push_back(current);
        }
        return;
    }

    int laterParts = parts - position - 1;
    for (int length = 1; length <= remaining - laterParts; length++)
    {
        current[position] = length;
        compositions(remaining - length, parts, position + 1, current, out);
    }
}

vector<vector<int>> concreteSequences(const string &word, const int &atomCount)
{
    vector<vector<int>> out;
    if ((int)(word.size()) > atomCount)
    {
        return out;
    }

    vector<vector<int>> runLengths;
    vector<int> current(word.size());
    compositions(atomCount, (int)(word.size()), 0, current, runLengths);

    for (const auto &lengths : runLengths)
    {
        vector<int> sequence;
        for (std::size_t i = 0; i < lengths.size(); i++)
        {
            for (int j = 0; j < lengths[i]; j++)
            {
                sequence.push_back(word[i] - '0');
            }
        }
        out.push_back(sequence);
    }

    return out;
}

struct TestTimedCell
{
    bool singleton;
    bool inWindow;
};

vector<TestTimedCell> testTimedCells(const TimedRange &horizon,
                                    const TimedRange &window)
{
    vector<__int128> cuts{horizon.left, horizon.right, window.left, window.right};
    sort(cuts.begin(), cuts.end());
    cuts.erase(unique(cuts.begin(), cuts.end()), cuts.end());

    vector<TestTimedCell> cells;
    for (int i = 0; i < cuts.size(); i++)
    {
        if (timedContains(horizon, cuts[i]))
        {
            cells.push_back({true, timedContains(window, cuts[i])});
        }

        if (i + 1 < cuts.size())
        {
            TimedRange open{cuts[i], cuts[i + 1], false, false};
            if (!timedEmpty(timedIntersection(open, horizon)))
            {
                cells.push_back({false, !timedEmpty(timedIntersection(open, window))});
            }
        }
    }
    return cells;
}

array<bool, 2> slowPossibleBits(const string &lhsWord, const string &rhsWord,
                                const vector<TestTimedCell> &cells)
{
    int transitionCount = (int)(lhsWord.size() + rhsWord.size()) - 2;
    vector<bool> inWindow;
    vector<bool> pointAtom;

    for (const auto &cell : cells)
    {
        if (cell.singleton)
        {
            inWindow.push_back(cell.inWindow);
            pointAtom.push_back(true);
        }
        else
        {
            for (int event = 0; event < transitionCount; event++)
            {
                inWindow.push_back(cell.inWindow);
                pointAtom.push_back(false);
                inWindow.push_back(cell.inWindow);
                pointAtom.push_back(true);
            }
            inWindow.push_back(cell.inWindow);
            pointAtom.push_back(false);
        }
    }

    const int windowAtomCount = count(
        inWindow.begin(), inWindow.end(), true);
    if (windowAtomCount == 0)
    {
        return {true, false};
    }

    vector<vector<int>> lhsSignals = concreteSequences(lhsWord, (int)(inWindow.size()));
    vector<vector<int>> rhsSignals = concreteSequences(rhsWord, windowAtomCount);
    array<bool, 2> possible{false, false};

    for (const auto &lhs : lhsSignals)
    {
        for (const auto &rhs : rhsSignals)
        {
            bool lhsGood = true;
            bool witness = false;
            int rhsIndex = 0;
            for (std::size_t i = 0; i < inWindow.size(); i++)
            {
                if (inWindow[i])
                {
                    bool canWitnessHere = pointAtom[i] || lhs[i] == 1;
                    if (lhsGood && canWitnessHere && rhs[rhsIndex] == 1)
                    {
                        witness = true;
                    }
                    rhsIndex++;
                }
                lhsGood = lhsGood && lhs[i] == 1;
            }
            possible[witness ? 1 : 0] = true;
        }
    }

    return possible;
}

void assertPossibleBitsMatchOracle(
    const string &lhsWord,
    const string &rhsWord,
    const TimedRange &horizon,
    const TimedRange &window)
{
    const array<bool, 2> expected = slowPossibleBits(
        lhsWord, rhsWord, testTimedCells(horizon, window));
    const BoolLanguage lhs = makeLanguage({lhsWord});
    const BoolLanguage rhs = makeLanguage({rhsWord});
    assert(possibleUntilBitsLegacy(lhs, rhs, horizon, window) == expected);
    assert(possibleUntilBitsBitParallel(lhs, rhs, horizon, window) == expected);
}

void testPossibleBitSolver()
{
    vector<string> words;
    for (int length = 1; length <= 4; length++)
    {
        for (int start = 0; start < 2; start++)
        {
            words.push_back(alternatingWord(start, length));
        }
    }

    struct TimedCase
    {
        TimedRange horizon;
        TimedRange window;
    };

    vector<TimedCase> cases{
        {{0, 0, true, true}, {0, 0, true, true}},
        {{0, 0, true, true}, {1, 0, true, true}},
        {{0, 2, true, false}, {0, 2, true, false}},
        {{0, 2, true, false}, {0, 2, false, false}},
        {{0, 2, true, true}, {0, 2, true, true}},
        {{0, 2, true, true}, {0, 2, false, false}},
        {{0, 2, true, true}, {0, 2, false, true}},
        {{0, 2, true, true}, {0, 2, true, false}},
        {{0, 4, true, true}, {1, 3, true, true}},
        {{0, 4, true, true}, {1, 3, false, true}},
        {{0, 4, true, true}, {1, 3, true, false}},
        {{0, 4, true, true}, {1, 3, false, false}},
        {{0, 2, true, true}, {1, 0, true, true}}
    };

    for (const auto &testCase : cases)
    {
        vector<TestTimedCell> cells = testTimedCells(
            testCase.horizon, testCase.window);
        for (const auto &lhsWord : words)
        {
            for (const auto &rhsWord : words)
            {
                const array<bool, 2> expected = slowPossibleBits(
                    lhsWord, rhsWord, cells);
                const BoolLanguage lhs = makeLanguage({lhsWord});
                const BoolLanguage rhs = makeLanguage({rhsWord});
                assert(possibleUntilBitsLegacy(
                           lhs, rhs, testCase.horizon, testCase.window) ==
                       expected);
                assert(possibleUntilBitsBitParallel(
                           lhs, rhs, testCase.horizon, testCase.window) ==
                       expected);
            }
        }
    }

    vector<string> smallerWords;
    for (int length = 1; length <= 3; length++)
    {
        for (int start = 0; start < 2; start++)
        {
            smallerWords.push_back(alternatingWord(start, length));
        }
    }

    TimedRange subsetHorizon{0, 2, true, true};
    TimedRange subsetWindow{0, 2, false, false};
    vector<TestTimedCell> cells = testTimedCells(subsetHorizon, subsetWindow);
    vector<vector<array<bool, 2>>> pairResults(
        smallerWords.size(), vector<array<bool, 2>>(smallerWords.size()));
    for (std::size_t i = 0; i < smallerWords.size(); i++)
    {
        for (std::size_t j = 0; j < smallerWords.size(); j++)
        {
            pairResults[i][j] = slowPossibleBits(
                smallerWords[i], smallerWords[j], cells);
        }
    }

    const unsigned int subsetCount = 1U << smallerWords.size();
    for (unsigned int lhsMask = 1; lhsMask < subsetCount; lhsMask++)
    {
        set<string> lhsWords;
        for (std::size_t i = 0; i < smallerWords.size(); i++)
        {
            if ((lhsMask & (1U << i)) != 0)
            {
                lhsWords.insert(smallerWords[i]);
            }
        }

        for (unsigned int rhsMask = 1; rhsMask < subsetCount; rhsMask++)
        {
            set<string> rhsWords;
            array<bool, 2> expected{false, false};
            for (std::size_t i = 0; i < smallerWords.size(); i++)
            {
                if ((lhsMask & (1U << i)) == 0)
                {
                    continue;
                }
                for (std::size_t j = 0; j < smallerWords.size(); j++)
                {
                    if ((rhsMask & (1U << j)) != 0)
                    {
                        rhsWords.insert(smallerWords[j]);
                        expected[0] = expected[0] || pairResults[i][j][0];
                        expected[1] = expected[1] || pairResults[i][j][1];
                    }
                }
            }

            const BoolLanguage lhs = makeLanguage(lhsWords);
            const BoolLanguage rhs = makeLanguage(rhsWords);
            assert(possibleUntilBitsLegacy(
                       lhs, rhs, subsetHorizon, subsetWindow) == expected);
            assert(possibleUntilBitsBitParallel(
                       lhs, rhs, subsetHorizon, subsetWindow) == expected);
        }
    }

    TimedRange point{0, 0, true, true};
    assertPossibleBitsMatchOracle("0", "1", point, point);
    assertPossibleBitsMatchOracle("01", "1", subsetHorizon, subsetWindow);
    assertPossibleBitsMatchOracle("010", "1", subsetHorizon, subsetWindow);

    const TimedRange shiftedHorizon{0, 8, true, true};
    const TimedRange closedShiftedWindow{2, 8, true, true};
    const TimedRange openUpperShiftedWindow{2, 8, true, false};
    assert(slowPossibleBits(
               "1", "10",
               testTimedCells(shiftedHorizon, closedShiftedWindow)) ==
           (array<bool, 2>{false, true}));
    assertPossibleBitsMatchOracle(
        "1", "10", shiftedHorizon, closedShiftedWindow);
    assert(slowPossibleBits(
               "1", "01",
               testTimedCells(shiftedHorizon, openUpperShiftedWindow)) ==
           (array<bool, 2>{false, true}));
    assertPossibleBitsMatchOracle(
        "1", "01", shiftedHorizon, openUpperShiftedWindow);

    const TimedRange singletonWindow{4, 4, true, true};
    assert(slowPossibleBits(
               "1", "01",
               testTimedCells(shiftedHorizon, singletonWindow)) ==
           (array<bool, 2>{false, false}));
    assertPossibleBitsMatchOracle(
        "1", "01", shiftedHorizon, singletonWindow);

    const TimedRange emptyWindow{4, 4, true, false};
    assert(slowPossibleBits(
               "0101", "0101",
               testTimedCells(shiftedHorizon, emptyWindow)) ==
           (array<bool, 2>{true, false}));
    assertPossibleBitsMatchOracle(
        "0101", "0101", shiftedHorizon, emptyWindow);
}
void testConcatenationOverflow()
{
    BoolLanguage longWord(2);
    longWord[0][SIZE - 1] = true;
    expectException<overflow_error>([&]()
    {
        bitsetConcat(longWord, makeLanguage({"01"}));
    });

    BoolLanguage dense(2);
    dense[0].set();
    dense[1].set();
    SegmentedLanguage truth{makeLanguage({"1"})};
    SegmentedLanguage rhs{dense};
    expectException<overflow_error>([&]()
    {
        bitsetBoundedUntil(
            truth,
            rhs,
            vector<long long>{0, 1},
            0,
            0,
            true,
            true);
    });
}

void testPublicInterface()
{
    vector<long long> segmentation{0, 1};
    SegmentedLanguage lhs{makeLanguage({"1"})};
    SegmentedLanguage rhs{makeLanguage({"1"})};

    auto closed = bitsetBoundedUntil(
        lhs, rhs, segmentation, 0, 1, true, true);
    assert(bitset2stringset(closed[0]) == set<string>{"1"});
    assert(bitset2stringset(bitsetBoundedUntil(
        lhs, rhs, segmentation, 0, 1, true, false)[0]) == set<string>{"1"});
    assert(bitset2stringset(bitsetBoundedUntil(
        lhs, rhs, segmentation, 0, 1, false, true)[0]) == set<string>{"1"});
    assert(bitset2stringset(bitsetBoundedUntil(
        lhs, rhs, segmentation, 0, 1, false, false)[0]) == set<string>{"1"});
    assert(bitset2stringset(bitsetBoundedUntil(
        lhs, rhs, segmentation, 0, 0, true, true)[0]) == set<string>{"1"});
    assert(bitset2stringset(bitsetUnboundedUntil(
        lhs, rhs, segmentation, 0, true)[0]) == set<string>{"1"});
    assert(bitset2stringset(bitsetUnboundedUntil(
        lhs, rhs, segmentation, 0, false)[0]) == set<string>{"1"});

    auto direct = bitsetTimedUntil(
        lhs, rhs, segmentation, 0, 1, false, true, true, 0, -1);
    assert(direct == closed);
}
void testEndToEnd()
{
    vector<long long> segmentation{0, 1, 3};
    SegmentedLanguage truth{makeLanguage({"1"}), makeLanguage({"1"})};
    SegmentedLanguage proposition{makeLanguage({"0"}), makeLanguage({"0", "01"})};

    SegmentedLanguage manuscript = bitsetBoundedUntil(
        truth, proposition, segmentation, 0, 1, true, false);
    assert(bitset2stringset(manuscript[0]) == set<string>({"0", "01", "010", "0101"}));

    SegmentedLanguage eventually = bitsetBoundedEventually(
        proposition, segmentation, 0, 1, true, false);
    assert(eventually == manuscript);

    vector<long long> twoSegments{0, 1, 2};
    SegmentedLanguage falseSignal{makeLanguage({"0"}), makeLanguage({"0"})};
    SegmentedLanguage trueSignal{makeLanguage({"1"}), makeLanguage({"1"})};

    SegmentedLanguage includesNow = bitsetBoundedUntil(
        falseSignal, trueSignal, twoSegments, 0, 1, true, true, 0, 1);
    assert(bitset2stringset(includesNow[0]) == set<string>{"1"});
    assert(includesNow[1][0].none() && includesNow[1][1].none());

    SegmentedLanguage excludesNow = bitsetBoundedUntil(
        falseSignal, trueSignal, twoSegments, 0, 1, false, true, 0, 1);
    assert(bitset2stringset(excludesNow[0]) == set<string>{"0"});

    SegmentedLanguage boundaryRhs{makeLanguage({"0"}), makeLanguage({"10"})};
    SegmentedLanguage includesUpper = bitsetBoundedUntil(
        trueSignal, boundaryRhs, twoSegments, 0, 1, true, true, 0, 1);
    SegmentedLanguage excludesUpper = bitsetBoundedUntil(
        trueSignal, boundaryRhs, twoSegments, 0, 1, true, false, 0, 1);
    assert(firstOracle(bitset2stringset(includesUpper[0])) == set<string>{"1"});
    assert(bitset2stringset(excludesUpper[0]) == set<string>{"01"});

    SegmentedLanguage clippedEmpty = bitsetBoundedUntil(
        trueSignal, trueSignal, twoSegments, 3, 4, true, true);
    assert(bitset2stringset(clippedEmpty[0]) == set<string>{"0"});
    assert(bitset2stringset(clippedEmpty[1]) == set<string>{"0"});

    SegmentedLanguage punctual = bitsetBoundedUntil(
        falseSignal, trueSignal, twoSegments, 0, 0, true, true);
    assert(bitset2stringset(punctual[0]) == set<string>{"1"});
    assert(bitset2stringset(punctual[1]) == set<string>{"1"});

    SegmentedLanguage untimed = bitsetUntilNonStrict(truth, proposition);
    SegmentedLanguage unbounded = bitsetUnboundedUntil(
        truth, proposition, segmentation, 0, true);
    assert(unbounded == untimed);

    SegmentedLanguage unboundedSubrange = bitsetUnboundedUntil(
        truth, proposition, segmentation, 0, true, 1, 2);
    assert(unboundedSubrange[0][0].none() && unboundedSubrange[0][1].none());
    assert(unboundedSubrange[1] == untimed[1]);

    SegmentedLanguage closedUnbounded = bitsetUnboundedUntil(
        falseSignal, trueSignal, twoSegments, 0, true, 0, 1);
    assert(bitset2stringset(closedUnbounded[0]) == set<string>{"1"});

    SegmentedLanguage openUnbounded = bitsetUnboundedUntil(
        falseSignal, trueSignal, twoSegments, 0, false, 0, 1);
    assert(bitset2stringset(openUnbounded[0]) == set<string>{"0"});

    SegmentedLanguage delayedUnbounded = bitsetUnboundedUntil(
        trueSignal, trueSignal, twoSegments, 1, true, 0, 1);
    assert(bitset2stringset(delayedUnbounded[0]) == set<string>{"1"});
}

void testRefinedWindowProfileEndToEnd()
{
    const vector<long long> segmentation{0, 1, 3, 4, 5, 7, 8};
    const SegmentedLanguage truth(
        segmentation.size() - 1, makeLanguage({"1"}));
    const SegmentedLanguage proposition{
        makeLanguage({"0"}),
        makeLanguage({"0", "01"}),
        makeLanguage({"0", "01", "1"}),
        makeLanguage({"01", "010", "1", "10"}),
        makeLanguage({"0", "1", "10"}),
        makeLanguage({"0", "10"})};

    const SegmentedLanguage optimized = bitsetBoundedUntil(
        truth, proposition, segmentation, 0, 2, true, false, 2, 3);
    const SegmentedLanguage legacy = bitsetTimedUntilLegacy(
        truth, proposition, segmentation, 0, 2, false,
        true, false, 2, 3);

    assert(optimized == legacy);
    assert(bitset2stringset(optimized[2]) == set<string>{"1"});
    assert(optimized == bitsetBoundedEventually(
        proposition, segmentation, 0, 2, true, false, 2, 3));
    for (std::size_t segment = 0; segment < optimized.size(); segment++)
    {
        if (segment != 2)
        {
            assert(optimized[segment][0].none());
            assert(optimized[segment][1].none());
        }
    }
}

void testLhsHorizonOffsetLengthBound()
{
    const vector<long long> segmentation{0, 1, 2, 3};
    const SegmentedLanguage lhs{
        makeLanguage({"01"}),
        makeLanguage({"1"}),
        makeLanguage({"1"})};
    const BoolLanguage ambiguous = makeLanguage({"0", "1"});
    const SegmentedLanguage rhs{
        makeLanguage({"0"}), ambiguous, ambiguous};
    const set<string> expected{"0", "01", "010", "0101"};

    // For J = [1,2], offset zero belongs to M_J but not C_J. The change in
    // lhs on [0,1) must therefore contribute to the open placement's bound.
    const SegmentedLanguage legacy = bitsetTimedUntilLegacy(
        lhs, rhs, segmentation, 1, 2, false,
        true, true, 0, 1);
    const SegmentedLanguage optimized = bitsetBoundedUntil(
        lhs, rhs, segmentation, 1, 2, true, true, 0, 1);
    assert(optimized == legacy);
    assert(bitset2stringset(optimized[0]) == expected);
}

void testEventuallyMetamorphism()
{
    const vector<long long> segmentation{0, 1, 3, 4, 6};
    const SegmentedLanguage truth(
        segmentation.size() - 1, makeLanguage({"1"}));
    const SegmentedLanguage proposition{
        makeLanguage({"0", "01"}),
        makeLanguage({"0", "1", "10"}),
        makeLanguage({"01", "010", "1"}),
        makeLanguage({"0", "10"})};

    struct FiniteInterval
    {
        long long a;
        long long b;
        bool leftClosed;
        bool rightClosed;
    };

    const vector<FiniteInterval> finiteIntervals{
        {0, 2, true, true},
        {0, 2, true, false},
        {0, 2, false, true},
        {0, 2, false, false},
        {1, 3, true, true},
        {1, 3, true, false},
        {1, 3, false, true},
        {1, 3, false, false},
        {2, 2, true, true},
        {2, 2, true, false},
        {7, 9, true, true}};

    for (const auto &interval : finiteIntervals)
    {
        assert(bitsetBoundedUntil(
                   truth, proposition, segmentation,
                   interval.a, interval.b,
                   interval.leftClosed, interval.rightClosed) ==
               bitsetBoundedEventually(
                   proposition, segmentation,
                   interval.a, interval.b,
                   interval.leftClosed, interval.rightClosed));
        assert(bitsetBoundedUntil(
                   truth, proposition, segmentation,
                   interval.a, interval.b,
                   interval.leftClosed, interval.rightClosed, 1, 3) ==
               bitsetBoundedEventually(
                   proposition, segmentation,
                   interval.a, interval.b,
                   interval.leftClosed, interval.rightClosed, 1, 3));
    }

    struct UnboundedInterval
    {
        long long a;
        bool leftClosed;
    };

    const vector<UnboundedInterval> unboundedIntervals{
        {0, true},
        {0, false},
        {2, true},
        {2, false},
        {7, true}};

    for (const auto &interval : unboundedIntervals)
    {
        assert(bitsetUnboundedUntil(
                   truth, proposition, segmentation,
                   interval.a, interval.leftClosed) ==
               bitsetUnboundedEventually(
                   proposition, segmentation,
                   interval.a, interval.leftClosed));
        assert(bitsetUnboundedUntil(
                   truth, proposition, segmentation,
                   interval.a, interval.leftClosed, 1, 3) ==
               bitsetUnboundedEventually(
                   proposition, segmentation,
                   interval.a, interval.leftClosed, 1, 3));
    }
}

void testEmptyWindowSkipsProfiles()
{
    const int segments = SIZE + 1;
    vector<long long> segmentation(segments + 1);
    SegmentedLanguage alternating(segments, BoolLanguage(2));
    for (int segment = 0; segment <= segments; segment++)
    {
        segmentation[segment] = segment;
        if (segment < segments)
        {
            alternating[segment][segment & 1][0] = true;
        }
    }

    const TimedRange fullDomain{
        0, 2 * (__int128)(segments), true, false};
    expectException<overflow_error>([&]()
    {
        timedProfileLegacy(alternating, segmentation, fullDomain);
    });

    const long long a = segments + 1;
    const long long b = segments + 2;
    const SegmentedLanguage legacy = bitsetTimedUntilLegacy(
        alternating, alternating, segmentation,
        a, b, false, true, true, 0, 1);
    const SegmentedLanguage optimized = bitsetTimedUntil(
        alternating, alternating, segmentation,
        a, b, false, true, true, 0, 1);

    assert(optimized == legacy);
    assert(bitset2stringset(optimized[0]) == set<string>{"0"});
    for (int segment = 1; segment < segments; segment++)
    {
        assert(optimized[segment][0].none());
        assert(optimized[segment][1].none());
    }
}

void testUntimedUntilRewriteRegressions()
{
    static_assert(SIZE >= 5,
                  "untimed-until capacity tests require at least five bits");

    SegmentedLanguage boundaryLhs(1, BoolLanguage(2));
    SegmentedLanguage boundaryRhs(1, BoolLanguage(2));
    if (SIZE % 2 == 0)
    {
        // weak class 1 grows from length SIZE - 1 to exactly SIZE
        boundaryLhs[0][0][1] = true;
        boundaryRhs[0][0][SIZE - 2] = true;
    }
    else
    {
        // weak class 3 grows from length SIZE - 1 to exactly SIZE
        boundaryLhs[0][1][2] = true;
        boundaryRhs[0][1][SIZE - 2] = true;
    }
    const SegmentedLanguage boundaryLhsBefore = boundaryLhs;
    const SegmentedLanguage boundaryRhsBefore = boundaryRhs;
    SegmentedLanguage boundary = bitsetUntilNonStrict(
        boundaryLhs, boundaryRhs, false, true);
    if (SIZE % 2 == 0)
    {
        assert(boundary[0][0][SIZE - 1]);
    }
    else
    {
        assert(boundary[0][1][SIZE - 1]);
    }
    assert(boundaryLhs == boundaryLhsBefore);
    assert(boundaryRhs == boundaryRhsBefore);

    SegmentedLanguage overflowLhs(1, BoolLanguage(2));
    SegmentedLanguage overflowRhs(1, BoolLanguage(2));
    if (SIZE % 2 == 0)
    {
        // weak class 3 would require output index SIZE
        overflowLhs[0][1][2] = true;
        overflowRhs[0][1][SIZE - 1] = true;
    }
    else
    {
        // weak class 1 would require output index SIZE
        overflowLhs[0][0][1] = true;
        overflowRhs[0][0][SIZE - 1] = true;
    }
    const SegmentedLanguage overflowLhsBefore = overflowLhs;
    const SegmentedLanguage overflowRhsBefore = overflowRhs;
    expectException<overflow_error>([&]()
    {
        bitsetUntilNonStrict(overflowLhs, overflowRhs, false, true);
    });
    assert(overflowLhs == overflowLhsBefore);
    assert(overflowRhs == overflowRhsBefore);

    // The same extrema are representable in the strong branch. Capacity
    // validation must be local to the active carry branch.
    SegmentedLanguage strong = bitsetUntilNonStrict(
        overflowLhs, overflowRhs, true, false);
    if (SIZE % 2 == 0)
    {
        assert(strong[0][1][SIZE - 1]);
    }
    else
    {
        assert(strong[0][0][SIZE - 1]);
    }
    assert(overflowLhs == overflowLhsBefore);
    assert(overflowRhs == overflowRhsBefore);

    const SegmentedLanguage oneSegment{makeLanguage({"1"})};
    const SegmentedLanguage noSegments;
    expectException<invalid_argument>([&]()
    {
        bitsetUntilNonStrict(oneSegment, noSegments);
    });

    SegmentedLanguage malformed{BoolLanguage(1)};
    expectException<invalid_argument>([&]()
    {
        bitsetUntilNonStrict(malformed, oneSegment);
    });
    expectException<out_of_range>([&]()
    {
        bitsetUntilNonStrict(oneSegment, oneSegment, true, false, -1, 1);
    });
    expectException<out_of_range>([&]()
    {
        bitsetUntilNonStrict(oneSegment, oneSegment, true, false, 1, 0);
    });
    expectException<out_of_range>([&]()
    {
        bitsetUntilNonStrict(oneSegment, oneSegment, true, false, 0, 2);
    });

    const vector<long long> segmentation{0, 1, 2, 3};
    const SegmentedLanguage truth{
        makeLanguage({"1"}), makeLanguage({"1"}), makeLanguage({"1"})};
    const SegmentedLanguage delayedTruth{
        makeLanguage({"0"}), makeLanguage({"0"}), makeLanguage({"1"})};
    const SegmentedLanguage full = bitsetUnboundedUntil(
        truth, delayedTruth, segmentation, 0, true);
    const SegmentedLanguage partial = bitsetUnboundedUntil(
        truth, delayedTruth, segmentation, 0, true, 1, 2);
    const SegmentedLanguage clipped = bitsetUntilNonStrict(
        truth, delayedTruth, true, false, 1, 2);

    assert(partial[1] == full[1]);
    assert(partial[1] != clipped[1]);
    assert(partial[0][0].none() && partial[0][1].none());
    assert(partial[2][0].none() && partial[2][1].none());
}

void performanceSmokeTest()
{
    BoolLanguage dense(2);
    for (int i = 0; i < 80; i++)
    {
        dense[0][i] = true;
        dense[1][i] = true;
    }

    TimedRange horizon{0, 2, true, true};
    TimedRange window{0, 2, false, false};
    array<bool, 2> result = possibleUntilBits(dense, dense, horizon, window);
    assert(result[0] && result[1]);
}

int main()
{
    initializeMasks();

    testExistingLanguageHelpers();
    testGeometry();
    testProfiles();
    testPossibleBitSolver();
    testConcatenationOverflow();
    testPublicInterface();
    testEndToEnd();
    testRefinedWindowProfileEndToEnd();
    testLhsHorizonOffsetLengthBound();
    testEventuallyMetamorphism();
    testEmptyWindowSkipsProfiles();
    testUntimedUntilRewriteRegressions();
    performanceSmokeTest();

    cout << "All timed-until tests passed." << endl;
    return 0;
}
