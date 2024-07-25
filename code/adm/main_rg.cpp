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
#include <set>
#include <bitset>
#include <algorithm>
#include <limits>
#include <cfloat>
#include <chrono>
#include "functions.cpp"
using namespace std;

#define REP 100

int main()
{
    /* set variables */
    vector<long long> D{4000, 8000, 16000, 32000};
    vector<long long> EPS{1000, 2000, 4000, 8000};
    int n = 2;
    int del = 1;

    /* for timed formulas */
    // long long a = 0;
    // long long b = 1000;
    // bool leftClosed = true;
    // bool rightClosed = false;

    ofstream results;
    string filename = "results.txt";
    results.open(filename);

    // for (long long b = 1000; b <= 16000; b += 1000)
    // {
    for (auto d : D)
    {
        // if (d <= b)
        // {
        //     continue;
        // }
        for (auto eps : EPS)
        {
            if (d < eps)
            {
                continue;
            }

            for (int c = 0; c < 100; c++)
            {
                /* read signal data */
                vector<vector<pair<long long, double>>> signalsReal(n);
                signalsReal[0] = getData("../../data/rg/" + to_string(d / 1000) + "_" + to_string(c) + ".txt", d / 1000);
                signalsReal[1] = getData("../../data/rg/" + to_string(d / 1000) + "_" + to_string(c + 100) + ".txt", d / 1000);
                
                /* prepare the bit masks */
                evenMask = generateBitmask(0);
                oddMask = generateBitmask(1);

                // bool invalidSignal = false;
                int numSegments;
                vector<vector<bitset<SIZE>>> test;
                chrono::time_point<chrono::system_clock> starttime;
                chrono::time_point<chrono::system_clock> endtime; 
                    
                starttime = chrono::system_clock::now();
                for (int rep = 0; rep < REP; rep++)
                {
                    /* read the signals as atomic propositions to improve performance */
                    vector<vector<pair<long long, double>>> signals = convertSignalsToBool(signalsReal);
                    
                    /* compute the uncertainty intervals and check if the signal is valid*/
                    vector<vector<vector<long long>>> uncertainties = computeUncertaintyIntervals(signals, eps, del);

                    /* compute the canonical segmentation */
                    vector<long long> segmentation = computeCanonicalSegmentation(signals, uncertainties, d);
                    numSegments = segmentation.size() - 1;

                    /* compute the value expressions */
                    vector<vector<set<string>>> valExprs = computeValueExpressions(signals, uncertainties, segmentation);

                    /* translate signals to atomic propositions */
                    vector<vector<vector<bitset<SIZE>>>> aps = convertSignalsToAtomicPropositions(valExprs, 0.0);

                    /* evaluate the formula */
                    // always(x1 and x2)
                    test = bitsetAlways(bitsetConjunction(aps[0], aps[1]));
                    
                    // always(x1 or x2)
                    // test = bitsetAlways(bitsetNegation(bitsetConjunction(bitsetNegation(aps[0]), bitsetNegation(aps[1]))));

                    // eventually(x1 and x2)
                    // test = bitsetEventually(bitsetConjunction(aps[0], aps[1]));

                    // eventually(x1 or x2)
                    // test = bitsetEventually(bitsetNegation(bitsetConjunction(bitsetNegation(aps[0]), bitsetNegation(aps[1]))));

                    // x1 until x2
                    // test = bitsetUntilStrict(aps[0], aps[1]);

                    // always(x1 implies eventually(x2))
                    // test = bitsetAlways(bitsetNegation(bitsetConjunction(aps[0], bitsetAlways(bitsetNegation(aps[1])))));
                       
                    // always(x1 implies eventually_[0,1)(x2))
                    // test = bitsetAlways(bitsetNegation(bitsetConjunction(aps[0], bitsetBoundedAlways(bitsetNegation(aps[1]), segmentation, 0, 1000, true, false))));
                        
                    // always(x1 implies eventually_[0,2)(x2))
                    // test = bitsetAlways(bitsetNegation(bitsetConjunction(aps[0], bitsetBoundedAlways(bitsetNegation(aps[1]), segmentation, 0, 2000, true, false))));
                        
                    // eventually_I(x1 or x2)
                    // test = bitsetBoundedEventually(bitsetNegation(bitsetConjunction(bitsetNegation(aps[0]), bitsetNegation(aps[1]))), segmentation, a, b, leftClosed, rightClosed);
                }

                endtime = chrono::system_clock::now();
                chrono::duration<double, milli> totalTime = endtime - starttime;

                string out = "";
                if (test[0][0].any() && test[0][1].any())
                {
                    out = "2";
                }
                else if (test[0][1].any())
                {
                    out = "1";
                }
                else
                {
                    out = "0";
                }

                string w = to_string(d / 1000) + " " + to_string(eps / 1000) + " " + to_string(del / 1000) + " " + to_string(c) + " ";
                // string w = to_string(b) + " " + to_string(d / 1000) + " " + to_string(eps / 1000) + " " + to_string(del / 1000) + " " + to_string(c) + " ";
                w = w + to_string(numSegments) + " " + to_string((totalTime.count() / REP) / 1000) + " " + out;
                results << w << endl;
                cout << w << endl;
            }
        }
    }
    // }

    results.close();

    return 0;
}
