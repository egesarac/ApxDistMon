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
#include "functions.cpp"
using namespace std;

int main()
{
    const vector<long long> D{16000, 32000, 64000, 128000, 256000, 512000};
    const vector<long long> EPS{2000, 4000, 8000};
    const vector<long long> CHUNK{4000, 8000, 16000, 32000};

    const int nSignals = 2;
    const int DEL = 1;
    const int REP = 5;

    ofstream results("onlineIncr_response01.txt");

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

                for (int c = 0; c < 100; ++c)
                {
                    vector<vector<pair<long long, double>>> signalsReal(nSignals);
                    signalsReal[0] = getData("dataNew/" + to_string(d / 1000) + "_" + to_string(c) + ".txt", d / 1000);
                    signalsReal[1] = getData("dataNew/" + to_string(d / 1000) + "_" + to_string(c + 100) + ".txt", d / 1000);

                    evenMask = generateBitmask(0);
                    oddMask = generateBitmask(1);
                    allExceptFirstMask = evenMask | oddMask;
                    allExceptFirstMask[0] = 0;

                    chrono::time_point<chrono::high_resolution_clock> tBeg, tEnd;

                    string seq;
                    bool lastBit0 = false;
                    bool lastBit1 = true;
                    bool lastBit0sub = true;
                    bool lastBit1sub = false;

                    int processed = 0;
                    int chunkIdx = 0;

                    tBeg = chrono::high_resolution_clock::now();
                    int lastSegmentEndPointFromPrevPortion = 0;
                    while (processed < d)
                    {
                        const int readStart = max<int>(0, lastSegmentEndPointFromPrevPortion - 2 * eps - 1000) / 1000;
                        const int realStart = lastSegmentEndPointFromPrevPortion;
                        const int readEnd = min<int>(d, processed + CHUNK_SZ) / 1000;
                        const int realEnd = min<int>(d, processed + CHUNK_SZ - eps);

                        auto signals = convertSignalsToBoolPartial(signalsReal, readStart, readEnd);

                        auto uncertainties = computeUncertaintyIntervals(signals, eps, DEL);
                        auto segmentation = computeCanonicalSegmentation(signals, uncertainties, d);
                        auto valExprs = computeValueExpressions(signals, uncertainties, segmentation);

                        auto itS = lower_bound(segmentation.begin(), segmentation.end(), realStart);
                        int processStart = int(itS - segmentation.begin());
                        if (segmentation[processStart] > realStart)
                        {
                            processStart--;
                        }

                        auto itE = lower_bound(segmentation.begin(), segmentation.end(), realEnd);
                        int processEnd = int(itE - segmentation.begin());
                        if (segmentation[processEnd] > realEnd)
                        {
                            processEnd--;
                        }
                        lastSegmentEndPointFromPrevPortion = segmentation[processEnd];

                        auto aps = convertSignalsToAtomicPropositions(valExprs, 0.0);
                        auto testsub = bitsetBoundedEventuallyPast(aps[1], segmentation, 0, 1000, true, false, processStart, processEnd);
                        auto test = bitsetAlwaysPast(bitsetDisjunction(bitsetNegation(aps[0]), testsub), processStart, processEnd, lastBit0, lastBit1);
                        int last = max(0, processEnd - 1);

                        if (processStart == processEnd)
                        {
                            test[last][0] = lastBit0;
                            test[last][1] = lastBit1;
                            testsub[last][0] = lastBit0sub;
                            testsub[last][1] = lastBit1sub;
                        }

                        lastBit0 = ((test[last][0] & evenMask).any() || (test[last][1] & oddMask).any());
                        lastBit1 = ((test[last][0] & oddMask).any() || (test[last][1] & evenMask).any());
                        lastBit0sub = ((testsub[last][0] & evenMask).any() || (testsub[last][1] & oddMask).any());
                        lastBit1sub = ((testsub[last][0] & oddMask).any() || (testsub[last][1] & evenMask).any());

                        char out;
                        if (lastBit0 && lastBit1)
                            out = '2';
                        else if (lastBit1)
                            out = '1';
                        else
                            out = '0';

                        seq.push_back(out);

                        processed += CHUNK_SZ;
                        ++chunkIdx;
                    }
                    tEnd = chrono::high_resolution_clock::now();

                    chrono::duration<double, milli> t = tEnd - tBeg;

                    string w = to_string(d / 1000) + " " + to_string(eps / 1000) + " " + to_string(CHUNK_SZ / 1000) + " " + to_string(c) + " " + to_string((t.count() / REP) / 1000) + " " + seq;
                    results << w << endl;
                    cout << w << endl;
                }
            }
        }
    }

    results.close();
    return 0;
}