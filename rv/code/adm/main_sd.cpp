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

#define REP 1

int main()
{
    /* set variables */
    vector<long long> N{2,3,4};
    vector<long long> EPS{2};
    int k = 50;
    int d = 1000;
    int del = 1;
    
   for (auto n : N) {
    for (auto eps : EPS) {
        eps = eps * k;

    /* read signal data */
    vector<vector<pair<long long, double>>> signalsReal(3 * n);
    for (int i = 0; i < 3 * n; i = i + 3)
    {
        signalsReal[i] = getData3d("../../data/sd/s1_uav_" + to_string(i / 3), 0);
        signalsReal[i + 1] = getData3d("../../data/sd/s1_uav_" + to_string(i / 3), 1);
        signalsReal[i + 2] = getData3d("../../data/sd/s1_uav_" + to_string(i / 3), 2);
    }
    
    /* prepare the bit masks */
    evenMask = generateBitmask(0);
    oddMask = generateBitmask(1);

    int numSegments;
    vector<vector<bitset<SIZE>>> test;
    chrono::time_point<chrono::system_clock> starttime;
    chrono::time_point<chrono::system_clock> endtime; 
    chrono::duration<double, milli> totalTime;
    
    starttime = chrono::system_clock::now();
    for (int rep = 0; rep < REP; rep++)
    {
        /* read the signals as atomic propositions to improve performance */
        //vector<vector<pair<long long, double>>> signals = convertSignalsToBool(signalsReal);

        /* compute the uncertainty intervals */
        vector<vector<vector<long long>>> uncertainties = computeUncertaintyIntervals(signalsReal, eps, del);

        /* compute the canonical segmentation */
        vector<long long> segmentation = computeCanonicalSegmentation(signalsReal, uncertainties, d);
        numSegments = segmentation.size() - 1;

        /* compute the value expressions */
        /* centralized monitor */
        vector<vector<set<string>>> valExprs = computeValueExpressions(signalsReal, uncertainties, segmentation);
        /* decentralized monitor (relative) */
        // set<int> leaders = {0, 1, 2};
        // vector<vector<set<string>>> valExprs = computeRelativeValueExpressions(leaders, signalsReal, uncertainties, segmentation);
       
        vector<vector<set<string>>> valExprsFinal;
        vector<set<string>> ve = valExprs[0];
        for (int i = 1; i < n; i++)
        {
            /* asynchronous product */
            vector<set<string>> vex = asyncProdStrDiffSqr(valExprs[0], valExprs[3 * i]); 
            vector<set<string>> vey = asyncProdStrDiffSqr(valExprs[1], valExprs[3 * i + 1]);
            vector<set<string>> vez = asyncProdStrDiffSqr(valExprs[2], valExprs[3 * i + 2]);
            vector<set<string>> vexy = asyncProdStrSum(vex, vey);
            vector<set<string>> vexyz = asyncProdStrSum(vexy, vez);

            /* fine */
            // vector<set<string>> vex = abstProdStrDiffSqr(valExprs[0], valExprs[3 * i]);
            // vector<set<string>> vey = abstProdStrDiffSqr(valExprs[1], valExprs[3 * i + 1]);
            // vector<set<string>> vez = abstProdStrDiffSqr(valExprs[2], valExprs[3 * i + 2]);
            // vector<set<string>> vexy = abstProdStrSum(vex, vey);
            // vector<set<string>> vexyz = abstProdStrSum(vexy, vez);

            valExprsFinal.push_back(vexyz);
        }

        /* translate signals to atomic propositions */
        vector<vector<vector<bitset<SIZE>>>> aps = convertSignalsToAtomicPropositions(valExprsFinal, 0.0);

        /* evaluate the formula */
        vector<vector<vector<bitset<SIZE>>>> testVec;
        for (int i = 0; i < n - 1; i++) 
        {
            testVec.push_back(bitsetAlways(aps[i]));
        }
        for (int i = 1; i < n - 1; i++) 
        {
            testVec[0] = bitsetConjunction(testVec[0], testVec[i]);
        }

        test = testVec[0];
    }
    
    endtime = chrono::system_clock::now();
    totalTime += endtime - starttime;

    cout << d << " " << eps << " " << del << " " << numSegments << " " << (totalTime.count() / REP) / 1000 << " " << test[0][0].any() << " " << test[0][1].any() << endl;
     }
   }
    return 0;
}