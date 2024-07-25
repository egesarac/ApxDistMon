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

#define REP 10

int main()
{
    /* set variables */
    vector<long long> N{2, 3, 4};
    vector<long long> EPS{1, 2, 4, 8};
    int k = 50;
    int d = 1000;
    int del = 1;
    
   for (auto n : N) {
    for (auto eps : EPS) {
        eps = eps * k;

    /* read signal data */
    vector<vector<pair<long long, double>>> signalsReal(n);
    for (int i = 0; i < n; i++)
    {
        signalsReal[i] = getData("../../data/wt/s5_tank_" + to_string(i), d / k);
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
        /* compute the uncertainty intervals */
        vector<vector<vector<long long>>> uncertainties = computeUncertaintyIntervals(signalsReal, eps, del);

        /* compute the canonical segmentation */
        vector<long long> segmentation = computeCanonicalSegmentation(signalsReal, uncertainties, d);
        numSegments = segmentation.size() - 1;

        /* compute the value expressions */
        /* centralized monitor */
        vector<vector<set<string>>> valExprs = computeValueExpressions(signalsReal, uncertainties, segmentation);
        /* decentralized monitor (relative) */
        // set<int> leaders = {0};
        // vector<vector<set<string>>> valExprs = computeRelativeValueExpressions(leaders, signalsReal, uncertainties, segmentation);
    
        vector<set<string>> ve = valExprs[0];
        for (int i = 1; i < n; i++)
        {   
            /* asynchronous product */
            ve = asyncProdStrSum(ve, valExprs[i]);

            /* fine */
            // ve = abstProdStrSum(ve, valExprs[i]);

            /* coarse */
            // ve = abstProdCoarseStrSum(ve, valExprs[i]);

        }
        valExprs = {ve};

        /* translate signals to atomic propositions */
        vector<vector<vector<bitset<SIZE>>>> aps = convertSignalsToAtomicPropositions(valExprs, 10.0);

        /* evaluate the formula */
        test = bitsetAlways(aps[0]);
    }
    
    endtime = chrono::system_clock::now();
    totalTime += endtime - starttime;

    cout << d << " " << eps << " " << del << " " << numSegments << " " << (totalTime.count() / REP) / 1000 << " " << test[0][0].any() << " " << test[0][1].any() << endl;
     }
   }
    return 0;
}