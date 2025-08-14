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

    ofstream results("onlineNaive_since.txt");

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

                    int processed = 0;
                    int chunkIdx = 0;

                    tBeg = chrono::high_resolution_clock::now();
                    while (processed < d)
                    {
                        const int realStart = processed;

                        const int realEnd = min<int>(d, realStart + CHUNK_SZ - eps);
                        const int readEnd = min<int>(d, realStart + CHUNK_SZ) / 1000;

                        auto signals = convertSignalsToBool(signalsReal, readEnd);

                        auto uncertainties = computeUncertaintyIntervals(signals, eps, DEL);
                        auto segmentation = computeCanonicalSegmentation(signals, uncertainties, d);
                        auto valExprs = computeValueExpressions(signals, uncertainties, segmentation);

                        auto it_n = lower_bound(segmentation.begin(), segmentation.end(), realEnd);
                        int processEnd = int(it_n - segmentation.begin());

                        if (segmentation[processEnd] > realEnd)
                        {
                            processEnd--;
                        }

                        auto aps = convertSignalsToAtomicPropositions(valExprs, 0.0);
                        auto test = bitsetSinceStrict(aps[0], aps[1], 0, processEnd);
                        int last     = max(0, processEnd - 1);

                        if (processEnd == 0) {
                            test[last][0] = true;
                            test[last][1] = false;
                        }
    
                        bool lastBit0 = ((test[last][0] & evenMask).any() || (test[last][1] & oddMask).any());
                        bool lastBit1 = ((test[last][0] & oddMask).any() || (test[last][1] & evenMask).any());
    
                        char out;
                        if (lastBit0 && lastBit1)      out = '2';
                        else if (lastBit1)                   out = '1';
                        else                                        out = '0';

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