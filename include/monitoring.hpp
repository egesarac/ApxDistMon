#pragma once

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <filesystem>
#include <sstream>
#include <vector>
#include <string>
#include <random>
#include <tuple>
#include <map>
#include <set>
#include <array>
#include <cstdint>
// #include <bit>
#include <bitset>
#include <algorithm>
#include <limits>
#include <cfloat>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <optional>
#include <stdexcept>
#include <utility>
using namespace std;

#define SIZE 1000 // TODO: this should be as tight as possible

bitset<SIZE> evenMask;
bitset<SIZE> oddMask;
bitset<SIZE> allExceptFirstMask;

string formatDouble(double value)
{
    if (!isfinite(value))
    {
        throw overflow_error("non-finite value in numeric value expression");
    }
    if (value == 0.0)
    {
        return "0";
    }

    ostringstream output;
    output << setprecision(numeric_limits<double>::max_digits10) << value;
    return output.str();
}

template <std::size_t N>
int msb(const bitset<N> &bs)
{
    for (int i = bs.size() - 1; i >= 0; i--)
    {
        if (bs[i] == true)
        {
            return i;
        }
    }

    return -1;
}

template <std::size_t N>
int lsb(const bitset<N> &bs)
{
    for (int i = 0; i < bs.size(); i++)
    {
        if (bs[i] == true)
        {
            return i;
        }
    }

    return -1;
}

inline int mostSignificantBit(unsigned long long word)
{
    int result = -1;
    while (word != 0)
    {
        word >>= 1;
        result++;
    }
    return result;
}

template <std::size_t N>
int wordOrientedMsb(bitset<N> bits)
{
    static_assert(N > 0, "bitsets require positive capacity");
    static_assert(N <= static_cast<size_t>(numeric_limits<int>::max()),
                  "bitset capacity must fit in an int");
    constexpr size_t wordBits =
        numeric_limits<unsigned long long>::digits;
    const bitset<N> lowWordMask(
        numeric_limits<unsigned long long>::max());

    int result = -1;
    size_t offset = 0;
    while (bits.any())
    {
        const unsigned long long word =
            (bits & lowWordMask).to_ullong();
        const int local = mostSignificantBit(word);
        if (local >= 0)
        {
            result = static_cast<int>(offset) + local;
        }
        bits >>= wordBits;
        offset += wordBits;
    }
    return result;
}

struct AlternatingEndpointSummary
{
    array<int, 4> nonSingletonLengths{};
    bool singletonFalse = false;
    bool singletonTrue = false;
};

template <std::size_t N>
AlternatingEndpointSummary summarizeAlternatingEndpoints(
    const bitset<N> &startsFalse,
    const bitset<N> &startsTrue,
    bool negated = false)
{
    static_assert(N > 0, "alternating languages require positive capacity");
    static_assert(N <= static_cast<size_t>(numeric_limits<int>::max()),
                  "alternating language capacity must fit in an int");

    constexpr size_t wordBits =
        numeric_limits<unsigned long long>::digits;
    static const array<unsigned long long, 2> parityMasks = []
    {
        array<unsigned long long, 2> result{};
        for (size_t index = 0; index < wordBits; index++)
        {
            result[index % 2] |= 1ULL << index;
        }
        return result;
    }();
    const bitset<N> lowWordMask(
        numeric_limits<unsigned long long>::max());

    bitset<N> logicalFalse = negated ? startsTrue : startsFalse;
    bitset<N> logicalTrue = negated ? startsFalse : startsTrue;
    AlternatingEndpointSummary summary;
    summary.singletonFalse = logicalFalse[0];
    summary.singletonTrue = logicalTrue[0];

    size_t offset = 0;
    while (logicalFalse.any() || logicalTrue.any())
    {
        unsigned long long falseWord =
            (logicalFalse & lowWordMask).to_ullong();
        unsigned long long trueWord =
            (logicalTrue & lowWordMask).to_ullong();
        if (offset == 0)
        {
            falseWord &= ~1ULL;
            trueWord &= ~1ULL;
        }

        const array<unsigned long long, 4> classes{
            falseWord & parityMasks[0],
            falseWord & parityMasks[1],
            trueWord & parityMasks[1],
            trueWord & parityMasks[0]};
        for (size_t endpointClass = 0; endpointClass < 4;
             endpointClass++)
        {
            const int local = mostSignificantBit(classes[endpointClass]);
            if (local >= 0)
            {
                summary.nonSingletonLengths[endpointClass] =
                    static_cast<int>(offset) + local + 1;
            }
        }

        logicalFalse >>= wordBits;
        logicalTrue >>= wordBits;
        offset += wordBits;
    }
    return summary;
}

template <std::size_t N>
array<bitset<N>, 2> reverseAlternatingSegment(
    const bitset<N> &startsFalse,
    const bitset<N> &startsTrue)
{
    static_assert(N > 0, "alternating languages require positive capacity");

    static const array<bitset<N>, 2> parityMasks = []
    {
        array<bitset<N>, 2> result;
        for (size_t index = 0; index < N; index++)
        {
            result[index % 2].set(index);
        }
        return result;
    }();

    return {
        (startsFalse & parityMasks[0]) |
            (startsTrue & parityMasks[1]),
        (startsTrue & parityMasks[0]) |
            (startsFalse & parityMasks[1])};
}

template <std::size_t N>
array<bitset<N>, 2> reverseAlternatingSegment(
    const array<bitset<N>, 2> &language)
{
    return reverseAlternatingSegment(language[0], language[1]);
}

inline AlternatingEndpointSummary reverseAlternatingSummary(
    const AlternatingEndpointSummary &summary)
{
    AlternatingEndpointSummary reversed = summary;
    swap(reversed.nonSingletonLengths[1],
         reversed.nonSingletonLengths[2]);
    return reversed;
}

set<double> splitSet(string str, string delimiter)
{
    set<double> v;
    if (!str.empty()) {
        int start = 0;
        do {
            // Find the index of occurrence
            int idx = str.find(delimiter, start);
            if (idx == string::npos) {
                break;
            }

            // If found add the substring till that
            // occurrence in the vector
            int length = idx - start;
            v.insert(stod(str.substr(start, length)));
            start += (length + delimiter.size());
        } while (true);
        v.insert(stod(str.substr(start)));
    }

    return v;
}

vector<double> split(string str, string delimiter)
{
    vector<double> v;
    if (!str.empty()) {
        int start = 0;
        do {
            // Find the index of occurrence
            int idx = str.find(delimiter, start);
            if (idx == string::npos) {
                break;
            }

            // If found add the substring till that
            // occurrence in the vector
            int length = idx - start;
            v.push_back(stod(str.substr(start, length)));
            start += (length + delimiter.size());
        } while (true);
        v.push_back(stod(str.substr(start)));
    }

    return v;
}

set<string> stutterOnce(const set<string> &S)
{
    set<string> out;

    for (const auto &s : S)
    {
        for (int i = 0; i < s.length(); i++)
        {
            string sl = s.substr(0, i);
            string sr = s.substr(i, s.length());

            out.insert(sl + s[i] + sr);
        }
    }

    return out;
}

set<string> stutterOnceStr(const set<string> &S)
{
    set<string> out;

    for (const auto &s : S)
    {
        vector<int> ind;
        ind.push_back(-1);
        for (int i = 0; i < s.length(); i++)
        {
            if (s[i] == ';')
            {
                ind.push_back(i);
            }
        }
        ind.push_back(s.length());

        int prev = 0;
        for (int j = 1; j < ind.size(); j++) {
            string sl = s.substr(0, ind[j]);
            string sr = s.substr(ind[j - 1] + 1, s.length());
            out.insert(sl + ";" + sr);
        }
    }

    return out;
}

set<string> stutter2k(const set<string> &S, const int &k)
{
    int len = S.begin()->length();
    int iter = k - len;

    set<string> out = S;
    for (int i = 0; i < iter; i++)
    {
        set<string> temp = stutterOnce(out);
        out = temp;
    }

    return out;
}

set<string> stutter2kStr(const set<string> &S, const int &k)
{
    int len = count(S.begin()->begin(), S.begin()->end(), ';') + 1;
    int iter = k - len;

    set<string> out = S;
    for (int i = 0; i < iter; i++)
    {
        set<string> temp = stutterOnceStr(out);
        out = temp;
    }

    return out;
}

template <std::size_t N>
set<string> bitset2stringset(const vector<bitset<N>> &v)
{
    set<string> out;

    for (int i = 0; i < N; i++)
    {
        if (v[0][i] == true)
        {
            string s = "";

            for (int j = 0; j <= i; j++)
            {
                if (j % 2 == 0)
                {
                    s += "0";
                }
                else
                {
                    s += "1";
                }
            }

            out.insert(s);
        }

        if (v[1][i] == true)
        {
            string s = "";

            for (int j = 0; j <= i; j++)
            {
                if (j % 2 == 0)
                {
                    s += "1";
                }
                else
                {
                    s += "0";
                }
            }

            out.insert(s);
        }
    }

    return out;
}

template <std::size_t N>
vector<set<string>> bitset2stringset_withSegments(const vector<vector<bitset<N>>> &v)
{
    vector<set<string>> out(v.size());

    for (int i = 0; i < v.size(); i++)
    {
        out[i] = bitset2stringset(v[i]);
    }

    return out;
}

string destutterStr(const string &t1)
{
    vector<int> ind1;
    ind1.push_back(-1);
    for (int i = 0; i < t1.length(); i++)
    {
        if (t1[i] == ';')
        {
            ind1.push_back(i);
        }
    }
    ind1.push_back(t1.length());

    string s1p = t1.substr(ind1[0] + 1, ind1[1] - ind1[0] - 1);
    string x = s1p;
    string s1;
    for (int i = 1; i < ind1.size() - 1; i++) {
        s1 = t1.substr(ind1[i] + 1, ind1[i + 1] - ind1[i] - 1);

        if (s1 != s1p) {
            x += ";" + s1;
        }

        s1p = s1;
    }

    return x;
}

pair<string, string> destutterPairStr(const string &t1, const string &t2)
{
    vector<int> ind1;
    ind1.push_back(-1);
    for (int i = 0; i < t1.length(); i++)
    {
        if (t1[i] == ';')
        {
            ind1.push_back(i);
        }
    }
    ind1.push_back(t1.length());

    vector<int> ind2;
    ind2.push_back(-1);
    for (int i = 0; i < t2.length(); i++)
    {
        if (t2[i] == ';')
        {
            ind2.push_back(i);
        }
    }
    ind2.push_back(t2.length());

    string s1p = t1.substr(ind1[0] + 1, ind1[1] - ind1[0] - 1);
    string s2p = t2.substr(ind2[0] + 1, ind2[1] - ind2[0] - 1);
    string x = s1p;
    string y = s2p;
    string s1, s2;
    for (int i = 1; i < ind1.size() - 1; i++) {
        s1 = t1.substr(ind1[i] + 1, ind1[i + 1] - ind1[i] - 1);
        s2 = t2.substr(ind2[i] + 1, ind2[i + 1] - ind2[i] - 1);

        if (s1 != s1p || s2 != s2p) {
            x += ";" + s1;
            y += ";" + s2;
        }

        s1p = s1;
        s2p = s2;
    }

    return make_pair(x, y);
}

pair<string, string> destutterPair(const string &t1, const string &t2)
{
    string s1 = t1.substr(0, 1);
    string s2 = t2.substr(0, 1);
    int len = t1.length();

    for (int i = 1; i < len; i++)
    {
        if (t1[i - 1] != t1[i] || t2[i - 1] != t2[i])
        {
            s1 += t1[i];
            s2 += t2[i];
        }
    }

    return make_pair(s1, s2);
}

vector<set<string>> abstProdStrSum(const vector<set<string>> &v1, const vector<set<string>> &v2)
{
    vector<set<string>> out(v1.size());

    for (int i = 0; i < v1.size(); i++)
    {
        for (auto s1 : v1[i])
        {
            if (!s1.empty())
            {
                set<double> vals1 = splitSet(s1, ";");

                for (auto s2 : v2[i])
                {
                    if(!s2.empty())
                    {
                        set<double> vals2 = splitSet(s2, ";");

                        for (auto e1 : vals1)
                        {
                            for (auto e2 : vals2)
                            {
                                out[i].insert(formatDouble(e1 + e2));
                            }
                        }
                    }
                }
            }
        }
    }

    return out;
}

// vector<set<string>> abstProdCoarseStrSum(const vector<set<string>> &v1, const vector<set<string>> &v2)
// {
//     vector<set<string>> out(v1.size());

//     for (int i = 0; i < v1.size(); i++)
//     {
//         double m1 = DBL_MAX;
//         double m2 = DBL_MAX;
//         double M1 = DBL_MIN;
//         double M2 = DBL_MIN;

//         for (auto s1 : v1[i])
//         {
//             if (!s1.empty())
//             {
//                 vector<double> vals1 = split(s1, ";");
//                 auto [min1, max1] = minmax_element(vals1.begin(), vals1.end());

//                 if (*min1 < m1)
//                 {
//                     m1 = *min1;
//                 }

//                 if (*max1 > M1)
//                 {
//                     M1 = *max1;
//                 }
//             }
//         }

//         for (auto s2 : v2[i])
//         {
//             if(!s2.empty())
//             {
//                 vector<double> vals2 = split(s2, ";");
//                 auto [min2, max2] = minmax_element(vals2.begin(), vals2.end());

//                 if (*min2 < m2)
//                 {
//                     m2 = *min2;
//                 }

//                 if (*max2 > M2)
//                 {
//                     M2 = *max2;
//                 }
//             }
//         }

//         out[i].insert(formatDouble(M1 + M2));
//         out[i].insert(formatDouble(m1 + m2));
//     }

//     return out;
// }

vector<set<string>> abstProdStrDiffSqr(const vector<set<string>> &v1, const vector<set<string>> &v2)
{
    vector<set<string>> out(v1.size());

    for (int i = 0; i < v1.size(); i++)
    {
        for (auto s1 : v1[i])
        {
            if (!s1.empty())
            {
                set<double> vals1 = splitSet(s1, ";");

                for (auto s2 : v2[i])
                {
                    if(!s2.empty())
                    {
                        set<double> vals2 = splitSet(s2, ";");

                        for (auto e1 : vals1)
                        {
                            for (auto e2 : vals2)
                            {
                                out[i].insert(formatDouble((e1 - e2) * (e1 - e2)));
                            }
                        }
                    }
                }
            }
        }
    }

    return out;
}


vector<set<string>> asyncProdStrDiffSqr(const vector<set<string>> &v1, const vector<set<string>> &v2)
{
    vector<set<string>> out(v1.size());

    for (int i = 0; i < v1.size(); i++)
    {
        for (auto s1 : v1[i])
        {
            for (auto s2 : v2[i])
            {
                set<string> temp1, temp2;

                temp1.insert(s1);
                temp2.insert(s2);

                auto n1 = count(s1.begin(), s1.end(), ';') + 1;
                auto n2 = count(s2.begin(), s2.end(), ';') + 1;

                temp1 = stutter2kStr(temp1, n1 + n2 - 1);
                temp2 = stutter2kStr(temp2, n1 + n2 - 1);

                for (auto t1 : temp1)
                {
                    for (auto t2 : temp2)
                    {
                        vector<int> ind1;
                        ind1.push_back(-1);
                        for (int i = 0; i < t1.length(); i++)
                        {
                            if (t1[i] == ';')
                            {
                                ind1.push_back(i);
                            }
                        }
                        ind1.push_back(t1.length());

                        vector<int> ind2;
                        ind2.push_back(-1);
                        for (int i = 0; i < t2.length(); i++)
                        {
                            if (t2[i] == ';')
                            {
                                ind2.push_back(i);
                            }
                        }
                        ind2.push_back(t2.length());



                        string r = "";
                        for (int i = 0; i < ind1.size() - 1; i++) {
                            string r1 = t1.substr(ind1[i] + 1, ind1[i + 1] - ind1[i] - 1);
                            string r2 = t2.substr(ind2[i] + 1, ind2[i + 1] - ind2[i] - 1);

                            if (!r1.empty() && !r2.empty())
                            {
                                r += formatDouble((stod(r1) - stod(r2)) * (stod(r1) - stod(r2))) + ";";
                            }
                        }


                        out[i].insert(destutterStr(r.substr(0, r.length() - 1)));
                    }
                }
            }
        }
    }

    return out;
}

vector<set<string>> asyncProdStrDiff(const vector<set<string>> &v1, const vector<set<string>> &v2)
{
    vector<set<string>> out(v1.size());

    for (int i = 0; i < v1.size(); i++)
    {
        for (auto s1 : v1[i])
        {
            for (auto s2 : v2[i])
            {
                set<string> temp1, temp2;

                temp1.insert(s1);
                temp2.insert(s2);

                auto n1 = count(s1.begin(), s1.end(), ';') + 1;
                auto n2 = count(s2.begin(), s2.end(), ';') + 1;

                temp1 = stutter2kStr(temp1, n1 + n2 - 1);
                temp2 = stutter2kStr(temp2, n1 + n2 - 1);

                for (auto t1 : temp1)
                {
                    for (auto t2 : temp2)
                    {
                        vector<int> ind1;
                        ind1.push_back(-1);
                        for (int i = 0; i < t1.length(); i++)
                        {
                            if (t1[i] == ';')
                            {
                                ind1.push_back(i);
                            }
                        }
                        ind1.push_back(t1.length());

                        vector<int> ind2;
                        ind2.push_back(-1);
                        for (int i = 0; i < t2.length(); i++)
                        {
                            if (t2[i] == ';')
                            {
                                ind2.push_back(i);
                            }
                        }
                        ind2.push_back(t2.length());



                        string r = "";
                        for (int i = 0; i < ind1.size() - 1; i++) {
                            string r1 = t1.substr(ind1[i] + 1, ind1[i + 1] - ind1[i] - 1);
                            string r2 = t2.substr(ind2[i] + 1, ind2[i + 1] - ind2[i] - 1);

                            if (!r1.empty() && !r2.empty())
                            {
                                r += formatDouble(stod(r1) - stod(r2)) + ";";
                            }
                        }


                        out[i].insert(destutterStr(r.substr(0, r.length() - 1)));
                    }
                }
            }
        }
    }

    return out;
}

vector<set<string>> asyncProdStrSum(const vector<set<string>> &v1, const vector<set<string>> &v2)
{
    vector<set<string>> out(v1.size());

    for (int i = 0; i < v1.size(); i++)
    {
        for (auto s1 : v1[i])
        {
            for (auto s2 : v2[i])
            {
                set<string> temp1, temp2;

                temp1.insert(s1);
                temp2.insert(s2);

                auto n1 = count(s1.begin(), s1.end(), ';') + 1;
                auto n2 = count(s2.begin(), s2.end(), ';') + 1;

                temp1 = stutter2kStr(temp1, n1 + n2 - 1);
                temp2 = stutter2kStr(temp2, n1 + n2 - 1);

                for (auto t1 : temp1)
                {
                    for (auto t2 : temp2)
                    {
                        vector<int> ind1;
                        ind1.push_back(-1);
                        for (int i = 0; i < t1.length(); i++)
                        {
                            if (t1[i] == ';')
                            {
                                ind1.push_back(i);
                            }
                        }
                        ind1.push_back(t1.length());

                        vector<int> ind2;
                        ind2.push_back(-1);
                        for (int i = 0; i < t2.length(); i++)
                        {
                            if (t2[i] == ';')
                            {
                                ind2.push_back(i);
                            }
                        }
                        ind2.push_back(t2.length());



                        string r = "";
                        for (int i = 0; i < ind1.size() - 1; i++) {
                            string r1 = t1.substr(ind1[i] + 1, ind1[i + 1] - ind1[i] - 1);
                            string r2 = t2.substr(ind2[i] + 1, ind2[i + 1] - ind2[i] - 1);

                            if (!r1.empty() && !r2.empty())
                            {
                                r += formatDouble(stod(r1) + stod(r2)) + ";";
                            }
                        }


                        out[i].insert(destutterStr(r.substr(0, r.length() - 1)));
                    }
                }
            }
        }
    }

    return out;
}

template <std::size_t N>
vector<set<pair<string, string>>> asyncProd(const vector<vector<bitset<N>>> &v1, const vector<vector<bitset<N>>> &v2)
{
    vector<set<pair<string, string>>> out(v1.size());

    for (int i = 0; i < v1.size(); i++)
    {
        set<string> S1 = bitset2stringset(v1[i]);
        set<string> S2 = bitset2stringset(v2[i]);

        for (auto s1 : S1)
        {
            for (auto s2 : S2)
            {
                set<string> temp1, temp2;

                temp1.insert(s1);
                temp2.insert(s2);

                temp1 = stutter2k(temp1, s1.length() + s2.length() - 1);
                temp2 = stutter2k(temp2, s1.length() + s2.length() - 1);

                for (auto t1 : temp1)
                {
                    for (auto t2 : temp2)
                    {
                        out[i].insert(destutterPair(t1, t2));
                    }
                }
            }
        }
    }

    return out;
}

vector<set<string>> prodAlways(const vector<set<string>> &product)
{
    vector<set<string>> out(product.size());

    bool firstBit0 = false;
    bool firstBit1 = true;

    for (int i = product.size() - 1; i >= 0; i--)
    {
        if (firstBit0 == true)
        {
            out[i].insert("0");
        }

        if (firstBit1 == true)
        {
            bool f1 = false;
            bool f2 = false;
            bool f3 = false;

            for (const auto &s : product[i])
            {
                if (s == "1")
                {
                    out[i].insert("1");
                    f1 = true;
                }

                else if (s == "0" || s.substr(s.length() - 2, 2) == "10")
                {
                    out[i].insert("0");
                    f2 = true;
                }

                else
                {
                    out[i].insert("01");
                    f3 = true;
                }

                if (f1 && f2 && f3)
                {
                    break;
                }
            }
        }

        if (out[i].find("0") != out[i].end() || out[i].find("01") != out[i].end())
        {
            firstBit0 = true;
        }
        else
        {
            firstBit0 = false;
        }

        if (out[i].find("1") != out[i].end())
        {
            firstBit1 = true;
        }
        else
        {
            firstBit1 = false;
        }
    }

    return out;
}

vector<set<string>> prodEventuallyPast(const vector<set<string>> &product)
{
    vector<set<string>> out(product.size());

    bool firstBit0 = true;
    bool firstBit1 = false;

    for (int i = 0; i < product.size(); i++)
    {

        if (firstBit0 == true)
        {
            bool f1 = false;
            bool f2 = false;
            bool f3 = false;

            for (const auto &s : product[i])
            {
                if (s == "0")
                {
                    out[i].insert("0");
                    f1 = true;
                }

                else if (s == "1" || s.substr(0, 2) == "10")
                {
                    out[i].insert("1");
                    f2 = true;
                }

                else
                {
                    out[i].insert("01");
                    f3 = true;
                }

                if (f1 && f2 && f3)
                {
                    break;
                }
            }
        }

        if (firstBit1 == true)
        {
            out[i].insert("1");
        }

        if (out[i].find("0") != out[i].end())
        {
            firstBit0 = true;
        }
        else
        {
            firstBit0 = false;
        }

        if (out[i].find("1") != out[i].end() || out[i].find("01") != out[i].end())
        {
            firstBit1 = true;
        }
        else
        {
            firstBit1 = false;
        }
    }

    return out;
}

vector<set<string>> prodEventually(const vector<set<string>> &product)
{
    vector<set<string>> out(product.size());

    bool firstBit0 = true;
    bool firstBit1 = false;

    for (int i = product.size() - 1; i >= 0; i--)
    {

        if (firstBit0 == true)
        {
            bool f1 = false;
            bool f2 = false;
            bool f3 = false;

            for (const auto &s : product[i])
            {
                if (s == "0")
                {
                    out[i].insert("0");
                    f1 = true;
                }

                else if (s == "1" || s.substr(s.length() - 2, 2) == "01")
                {
                    out[i].insert("1");
                    f2 = true;
                }

                else
                {
                    out[i].insert("10");
                    f3 = true;
                }

                if (f1 && f2 && f3)
                {
                    break;
                }
            }
        }

        if (firstBit1 == true)
        {
            out[i].insert("1");
        }

        if (out[i].find("0") != out[i].end())
        {
            firstBit0 = true;
        }
        else
        {
            firstBit0 = false;
        }

        if (out[i].find("1") != out[i].end() || out[i].find("10") != out[i].end())
        {
            firstBit1 = true;
        }
        else
        {
            firstBit1 = false;
        }
    }

    return out;
}

vector<set<string>> prodSinceStrict(const vector<set<pair<string, string>>> &product, bool lastBit)
{
    vector<set<string>> out(product.size());
    vector<set<char>> lastBits(product.size() + 1);

    char lb = lastBit ? '1' : '0';
    lastBits[0].insert(lb);

    for (int i = 0; i < product.size(); i++)
    {
        set<string> temp;

        for (auto pp : product[i])
        {
            for (const auto b : lastBits[i])
            {
                auto p = pp;
                int len = p.second.length() + 1;
                string s = "";
                if (b == '0')
                {
                    p.first = "0" + p.first;
                    p.second = "0" + p.second;
                    s += "0";
                }
                else {
                    p.first = "1" + p.first;
                    p.second = "1" + p.second;
                    s += "1";
                }

                for (int j = 1; j <= len - 1; j++)
                {
                    if (p.first[j] == '0')
                    {
                        s += "0";
                    }
                    else if (p.second[j] == '1')
                    {
                        s += "1";
                    }
                    else {
                        if (p.second[j - 1] == '1')
                        {
                            s += "1";
                        }
                        else if (p.first[j - 1] == '1')
                        {
                            s += s[s.length() - 1];
                        }
                        else
                        {
                            s += "0";
                        }
                    }
                }

                s = s.substr(1, s.length());
                string ss = s.substr(0, 1);
                for (int j = 1; j < len - 1; j++)
                {
                    if (s[j - 1] != s[j])
                    {
                        ss += s[j];
                    }
                }

                // string sss(ss);
                // std::reverse(sss.begin(), sss.end());

                // temp.insert(sss);
                temp.insert(ss);
            }
        }

        for (const auto &t : temp)
        {
            lastBits[i+1].insert(t[t.length()-1]);
        }

        out[i] = temp;
    }

    return out;
}

vector<set<string>> prodUntilStrict(const vector<set<pair<string, string>>> &product, bool firstBit)
{
    vector<set<string>> out(product.size());
    vector<set<char>> firstBits(product.size() + 1);

    char fb = firstBit ? '1' : '0';
    firstBits[product.size()].insert(fb);

    for (int i = product.size() - 1; i >= 0; i--)
    {
        set<string> temp;

        for (auto pp : product[i])
        {
            for (const auto b : firstBits[i + 1])
            {
                auto p = pp;
                int len = p.second.length() + 1;
                string s = "";
                if (b == '0')
                {
                    p.first = p.first + "0";
                    p.second = p.second + "0";
                    s += "0";
                }
                else {
                    p.first = p.first + "1";
                    p.second = p.second + "1";
                    s += "1";
                }

                for (int j = len - 2; j >= 0; j--)
                {
                    if (p.first[j] == '0')
                    {
                        s += "0";
                    }
                    else if (p.second[j] == '1')
                    {
                        s += "1";
                    }
                    else {
                        if (p.second[j + 1] == '1')
                        {
                            s += "1";
                        }
                        else if (p.first[j + 1] == '1')
                        {
                            s += s[s.length() - 1];
                        }
                        else
                        {
                            s += "0";
                        }
                    }
                }

                s = s.substr(1, s.length());
                string ss = s.substr(0, 1);
                for (int j = 1; j < len - 1; j++)
                {
                    if (s[j - 1] != s[j])
                    {
                        ss += s[j];
                    }
                }

                string sss(ss);
                std::reverse(sss.begin(), sss.end());

                temp.insert(sss);
            }
        }

        for (const auto &t : temp)
        {
            firstBits[i].insert(t[0]);
        }

        out[i] = temp;
    }

    return out;
}

vector<set<string>> prodSinceStrict_reverse(const vector<set<pair<string, string>>> &product, bool lastBit)
{
    vector<set<pair<string, string>>> r(product.size());

    // reverse(r.begin(),r.end());

    int n = product.size() - 1;

    for (int t = n; t >= 0; t--) {
        for (auto e : product[t]) {
            pair<string, string> temp(string(e.first.rbegin(), e.first.rend()), string(e.second.rbegin(), e.second.rend()));
            r[n-t].insert(temp);
        }
    }

    // for (int i = 0; i < r.size(); i++) {
    //     for (auto &e : r[i]){
    //         e.first = string(e.first.rbegin(), e.first.rend());
    //         e.second = string(e.second.rbegin(), e.second.rend());
    //         // reverse(e.first.begin(), e.first.end());
    //         // reverse(e.second.begin(), e.second.end());
    //     }
    // }

    return prodUntilStrict(r, lastBit); // check
}

vector<set<string>> prodUntilNonStrict(const vector<set<pair<string, string>>> &product, bool firstBit)
{
    vector<set<string>> out(product.size());
    vector<set<char>> firstBits(product.size() + 1);

    char fb = firstBit ? '1' : '0';
    firstBits[product.size()].insert(fb);

    for (int i = product.size() - 1; i >= 0; i--)
    {
        set<string> temp;

        for (const auto &p : product[i])
        {
            int len = p.second.length();

            for (const auto b : firstBits[i + 1])
            {
                string s = "";
                s += b;

                for (int j = len - 1; j >= 0; j--)
                {
                    if (p.second[j] == '1' || (p.first[j] == '1' && s[len - 1 - j] == '1'))
                    {
                        s += "1";
                    }
                    else
                    {
                        s += "0";
                    }
                }

                s = s.substr(1, s.length());
                string ss = s.substr(0, 1);
                for (int j = 1; j < len; j++)
                {
                    if (s[j - 1] != s[j])
                    {
                        ss += s[j];
                    }
                }

                string sss(ss);
                std::reverse(sss.begin(), sss.end());

                temp.insert(sss);
            }
        }

        for (const auto &t : temp)
        {
            firstBits[i].insert(t[0]);
        }

        out[i] = temp;
    }

    return out;
}

vector<set<string>> prodConjunction(const vector<set<pair<string, string>>> &product)
{
    vector<set<string>> out(product.size());

    for (int i = 0; i < product.size(); i++)
    {
        set<string> temp;

        for (const auto &p : product[i])
        {
            int len = p.first.length();
            string s = "";

            for (int j = 0; j < len; j++)
            {
                if (p.first[j] == '1' && p.second[j] == '1')
                {
                    s += "1";
                }
                else
                {
                    s += "0";
                }
            }

            string ss = s.substr(0, 1);
            for (int j = 1; j < len; j++)
            {
                if (s[j - 1] != s[j])
                {
                    ss += s[j];
                }
            }

            temp.insert(ss);
        }

        out[i] = temp;
    }

    return out;
}

vector<set<string>> prodNegation(const vector<set<string>> &product)
{
    vector<set<string>> out(product.size());

    for (int i = 0; i < product.size(); i++)
    {
        for (const auto &p : product[i])
        {
            int len = p.length();
            string s = "";

            for (int j = 0; j < len; j++)
            {
                if (p[j] == '1')
                {
                    s += "0";
                }
                else
                {
                    s += "1";
                }
            }

            out[i].insert(s);
        }
    }

    return out;
}

template <std::size_t N>
bool isEqual(vector<set<string>> productResult, vector<vector<bitset<N>>> bitsetResult)
{
    bool flag = true;

    for (int i = 0; i < productResult.size(); i++)
    {
        set<string> S = bitset2stringset(bitsetResult[i]);

        if (S != productResult[i])
        {
            flag = false;
            // cout << "Error in segment " << i << endl;
        }
    }

    return flag;
}

bool isEqualProd(vector<set<string>> productResult, vector<set<string>> productResultNew)
{
    bool flag = true;

    for (int i = 0; i < productResult.size(); i++)
    {
        if (productResultNew[i] != productResult[i])
        {
            flag = false;
            // cout << "Error in segment " << i << endl;
        }
    }

    return flag;
}

template <std::size_t N>
vector<bitset<N>> segmentFirstBit(vector<bitset<N>> v1)
{
    bool hasZero = v1[0].any();
    bool hasOne = v1[1].any();

    v1[0].reset();
    v1[1].reset();
    v1[0][0] = hasZero;
    v1[1][0] = hasOne;

    return v1;
}

template <std::size_t N>
vector<bitset<N>> segmentLastBit(vector<bitset<N>> v1)
{
    bool hasZero = (v1[0] & evenMask).any() || (v1[1] & oddMask).any();
    bool hasOne = (v1[0] & oddMask).any() || (v1[1] & evenMask).any();

    v1[0].reset();
    v1[1].reset();
    v1[0][0] = hasZero;
    v1[1][0] = hasOne;

    return v1;
}

template <std::size_t N>
pair<bool, bool> bitsetLastBits(const vector<bitset<N>> &language)
{
    if (language.size() != 2)
    {
        throw invalid_argument("a Boolean language must have two bitset buckets");
    }

    bool hasZero = (language[0] & evenMask).any() ||
                   (language[1] & oddMask).any();
    bool hasOne = (language[0] & oddMask).any() ||
                  (language[1] & evenMask).any();

    return make_pair(hasZero, hasOne);
}

template <std::size_t N>
char bitsetVerdict(const vector<bitset<N>> &language)
{
    auto [hasZero, hasOne] = bitsetLastBits(language);

    if (!hasZero && !hasOne)
    {
        throw logic_error("cannot extract a verdict from an empty language");
    }
    if (hasZero && hasOne)
    {
        return '2';
    }
    if (hasOne)
    {
        return '1';
    }
    return '0';
}

vector<long long> refineSegmentation(const vector<long long> &canonical,
                                     const vector<long long> &finalizationPoints,
                                     long long left,
                                     long long right)
{
    if (left > right)
    {
        throw invalid_argument("segmentation interval has negative length");
    }

    vector<long long> refined;
    refined.push_back(left);

    for (const auto &endpoint : canonical)
    {
        if (left < endpoint && endpoint < right)
        {
            refined.push_back(endpoint);
        }
    }

    for (const auto &endpoint : finalizationPoints)
    {
        if (left < endpoint && endpoint < right)
        {
            refined.push_back(endpoint);
        }
    }

    if (left < right)
    {
        refined.push_back(right);
    }

    sort(refined.begin(), refined.end());
    refined.erase(unique(refined.begin(), refined.end()), refined.end());

    return refined;
}

vector<pair<long long, double>> retainSignalFrom(
    const vector<pair<long long, double>> &signal,
    long long cutoff)
{
    if (signal.empty())
    {
        throw invalid_argument("cannot prune an empty signal");
    }

    double value = signal.front().second;
    for (const auto &edge : signal)
    {
        if (edge.first > cutoff)
        {
            break;
        }
        value = edge.second;
    }

    vector<pair<long long, double>> retained;
    retained.push_back(make_pair(cutoff, value));
    for (const auto &edge : signal)
    {
        if (edge.first > cutoff)
        {
            retained.push_back(edge);
        }
    }

    return retained;
}

struct RetainedApproximateSignal
{
    vector<pair<long long, double>> signal;
    optional<long long> predecessorLowerBound;
};

// Pruning replaces the discarded prefix with its value at the cutoff. Keep
// the propagated lower bound of the last discarded actual edge separately so
// that recomputing the remaining regions preserves same-signal edge order.
RetainedApproximateSignal retainApproximateSignalFrom(
    const vector<pair<long long, double>> &signal,
    const vector<vector<long long>> &uncertainties,
    long long cutoff,
    optional<long long> predecessorLowerBound = nullopt)
{
    if (signal.size() != uncertainties.size())
    {
        throw invalid_argument(
            "a signal and its uncertainty regions are misaligned");
    }

    for (size_t edge = 1; edge < signal.size(); ++edge)
    {
        if (signal[edge].first > cutoff)
        {
            break;
        }
        if (signal[edge - 1].second == signal[edge].second)
        {
            continue;
        }
        if (uncertainties[edge].size() != 2)
        {
            throw invalid_argument("an actual edge has no uncertainty region");
        }
        predecessorLowerBound = uncertainties[edge][0];
    }

    return {
        retainSignalFrom(signal, cutoff),
        predecessorLowerBound};
}

// Online evaluation finalizes only through observedUntil - eps. Extending the
// region-construction horizon by eps keeps every observed edge's natural upper
// bound stable as the observed prefix grows.
long long onlineUncertaintyHorizon(long long observedUntil, long long eps)
{
    if (observedUntil <= 0)
    {
        throw invalid_argument("the observed horizon must be positive");
    }
    if (eps <= 0)
    {
        throw invalid_argument("clock skew must be positive");
    }
    if (observedUntil > numeric_limits<long long>::max() - eps)
    {
        throw overflow_error("online uncertainty horizon overflows");
    }
    return observedUntil + eps;
}


// Restrict a language to the same interval without its first point.
// The first value either persists to the right or occurs only at that point.
template <std::size_t N>
vector<bitset<N>> segmentWithoutFirstPoint(const vector<bitset<N>> &v1)
{
    vector<bitset<N>> out = v1;
    for (int first = 0; first < 2; first++)
    {
        for (int i = 1; i <= msb(v1[first]); i++)
        {
            if (v1[first][i])
            {
                out[1 - first][i - 1] = true;
            }
        }
    }
    return out;
}

template <std::size_t N>
vector<bitset<N>> segmentInfix(vector<bitset<N>> v1)
{
    int maxInd = max(msb(v1[0] & (evenMask | oddMask)), msb(v1[1] & (evenMask | oddMask)));

    for (int j = 0; j < maxInd; j++)
    {
        v1[0][j] = true;
        v1[1][j] = true;
    }

    return v1;
}

template <std::size_t N>
vector<bitset<N>> segmentSuffix(vector<bitset<N>> v1)
{
    int maxInd0 = max(msb(v1[0] & evenMask), msb(v1[1] & oddMask));
    int maxInd1 = max(msb(v1[0] & oddMask), msb(v1[1] & evenMask));

    for (int j = maxInd0; j >= 0; j--)
    {
        if (j % 2 == 0)
        {
            v1[0][j] = true;
        }
        else
        {
            v1[1][j] = true;
        }
    }

    for (int j = maxInd1; j >= 0; j--)
    {
        if (j % 2 == 0)
        {
            v1[1][j] = true;
        }
        else
        {
            v1[0][j] = true;
        }
    }

    return v1;
}

template <std::size_t N>
vector<bitset<N>> segmentPrefix(vector<bitset<N>> v1)
{
    int maxInd0 = msb(v1[0] & (evenMask | oddMask));
    int maxInd1 = msb(v1[1] & (evenMask | oddMask));

    for (int j = 0; j <= maxInd0; j++)
    {
        v1[0][j] = true;
    }

    for (int j = 0; j <= maxInd1; j++)
    {
        v1[1][j] = true;
    }

    return v1;
}

template <std::size_t N>
vector<vector<bitset<N>>> bitsetInfix(vector<vector<bitset<N>>> v1)
{
    for (int i = 0; i < v1.size(); i++)
    {
        int maxInd = max(msb(v1[i][0] & (evenMask | oddMask)), msb(v1[i][1] & (evenMask | oddMask)));

        for (int j = 0; j < maxInd; j++)
        {
            v1[i][0][j] = true;
            v1[i][1][j] = true;
        }
    }

    return v1;
}

template <std::size_t N>
vector<vector<bitset<N>>> bitsetSuffix(vector<vector<bitset<N>>> v1)
{
    for (int i = 0; i < v1.size(); i++)
    {
        int maxInd0 = max(msb(v1[i][0] & evenMask), msb(v1[i][1] & oddMask));
        int maxInd1 = max(msb(v1[i][0] & oddMask), msb(v1[i][1] & evenMask));

        for (int j = maxInd0; j >= 0; j--)
        {
            if (j % 2 == 0)
            {
                v1[i][0][j] = true;
            }
            else
            {
                v1[i][1][j] = true;
            }
        }

        for (int j = maxInd1; j >= 0; j--)
        {
            if (j % 2 == 0)
            {
                v1[i][1][j] = true;
            }
            else
            {
                v1[i][0][j] = true;
            }
        }
    }

    return v1;
}

template <std::size_t N>
vector<vector<bitset<N>>> bitsetPrefix(vector<vector<bitset<N>>> v1)
{
    for (int i = 0; i < v1.size(); i++)
    {
        int maxInd0 = msb(v1[i][0] & (evenMask | oddMask));
        int maxInd1 = msb(v1[i][1] & (evenMask | oddMask));

        for (int j = 0; j <= maxInd0; j++)
        {
            v1[i][0][j] = true;
        }

        for (int j = 0; j <= maxInd1; j++)
        {
            v1[i][1][j] = true;
        }
    }

    return v1;
}

template <std::size_t N>
vector<bitset<N>> bitsetConcatLegacy(const vector<bitset<N>> &v1, const vector<bitset<N>> &v2)
{
    if (v1[0].none() && v1[1].none())
    {
        return v2;
    }

    vector<bitset<N>> out(2);
    for (int first = 0; first < 2; first++)
    {
        int leftMax = msb(v1[first]);
        for (int i = 0; i <= leftMax; i++)
        {
            if (!v1[first][i])
            {
                continue;
            }

            int last = first ^ (i % 2);
            for (int second = 0; second < 2; second++)
            {
                int rightMax = msb(v2[second]);
                for (int j = 0; j <= rightMax; j++)
                {
                    if (v2[second][j])
                    {
                        int index = i + j + (last != second);
                        if (index >= N)
                        {
                            throw overflow_error("value expression exceeds the fixed bitset size");
                        }
                        out[first][index] = true;
                    }
                }
            }
        }
    }

    return out;
}

// Timed evaluation uses doubled coordinates. Integer times are even, while
// the representative of an open interval (x,y) is the exact midpoint x+y.
struct TimedRange
{
    __int128 left;
    __int128 right;
    bool leftClosed;
    bool rightClosed;
};

bool timedEmpty(const TimedRange &range)
{
    return range.left > range.right ||
           (range.left == range.right && !(range.leftClosed && range.rightClosed));
}

bool timedContains(const TimedRange &range, const __int128 &point)
{
    if (timedEmpty(range) || point < range.left || point > range.right)
    {
        return false;
    }
    if (point == range.left && !range.leftClosed)
    {
        return false;
    }
    if (point == range.right && !range.rightClosed)
    {
        return false;
    }
    return true;
}

TimedRange timedIntersection(const TimedRange &x, const TimedRange &y)
{
    TimedRange out;
    out.left = max(x.left, y.left);
    out.right = min(x.right, y.right);
    out.leftClosed = (x.left < out.left || x.leftClosed) &&
                     (y.left < out.left || y.leftClosed);
    out.rightClosed = (x.right > out.right || x.rightClosed) &&
                      (y.right > out.right || y.rightClosed);
    return out;
}

template <std::size_t N>
vector<bitset<N>> timedProfileLegacy(const vector<vector<bitset<N>>> &v,
                               const vector<long long> &segmentation,
                               const TimedRange &horizon)
{
    vector<bitset<N>> profile(2);

    for (int i = 0; i < v.size(); i++)
    {
        TimedRange segment{2 * (__int128)(segmentation[i]),
                           2 * (__int128)(segmentation[i + 1]), true, false};
        TimedRange part = timedIntersection(segment, horizon);

        if (timedEmpty(part))
        {
            continue;
        }

        vector<bitset<N>> restricted;
        if (part.left == segment.left && part.right == segment.right &&
            part.leftClosed == segment.leftClosed && part.rightClosed == segment.rightClosed)
        {
            restricted = v[i];
        }
        else if (part.left == part.right && part.left == segment.left)
        {
            restricted = segmentFirstBit(v[i]);
        }
        else if (part.left == segment.left)
        {
            if (part.right == segment.right)
            {
                restricted = v[i];
            }
            else
            {
                restricted = segmentPrefix(v[i]);
            }

            if (!part.leftClosed)
            {
                restricted = segmentWithoutFirstPoint(restricted);
            }
        }
        else if (part.right == segment.right)
        {
            restricted = segmentSuffix(v[i]);
        }
        else
        {
            restricted = segmentInfix(v[i]);
        }

        profile = bitsetConcatLegacy(profile, restricted);
    }

    return profile;
}

template <std::size_t N>
array<bool, 2> possibleUntilBitsLegacy(const vector<bitset<N>> &lhs,
                                 const vector<bitset<N>> &rhs,
                                 const TimedRange &horizon,
                                 const TimedRange &window)
{
    // lhs describes the whole horizon, while rhs describes only the witness
    // window. Keep rhs fixed whenever the horizon search is outside window.
    if (timedEmpty(window))
    {
        return {true, false};
    }

    struct Cell
    {
        bool point;
        bool inWindow;
    };

    vector<__int128> cuts{horizon.left, horizon.right, window.left, window.right};
    sort(cuts.begin(), cuts.end());
    cuts.erase(unique(cuts.begin(), cuts.end()), cuts.end());

    vector<Cell> cells;
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

    array<bool, 2> possible{false, false};

    // All alternating words with the same first bit are prefixes of one
    // alternating sequence. One search therefore covers every represented
    // length in the corresponding bitset bucket.
    for (int lhsFirst = 0; lhsFirst < 2; lhsFirst++)
    {
        int lhsMax = msb(lhs[lhsFirst]);
        if (lhsMax < 0)
        {
            continue;
        }

        for (int rhsFirst = 0; rhsFirst < 2; rhsFirst++)
        {
            int rhsMax = msb(rhs[rhsFirst]);
            if (rhsMax < 0)
            {
                continue;
            }

            int lhsCount = lhsMax + 1;
            int rhsCount = rhsMax + 1;
            vector<char> seen(cells.size() * lhsCount * rhsCount * 4, false);
            vector<int> todo;

            auto push = [&](int cell, int lhsIndex, int rhsIndex,
                            bool lhsGood, bool witness)
            {
                if (lhsIndex >= lhsCount || rhsIndex >= rhsCount)
                {
                    return;
                }

                int state = (((cell * lhsCount + lhsIndex) * rhsCount + rhsIndex) * 4) +
                            2 * lhsGood + witness;
                if (!seen[state])
                {
                    seen[state] = true;
                    todo.push_back(state);
                }
            };

            push(0, 0, 0, true, false);

            while (!todo.empty())
            {
                int state = todo.back();
                todo.pop_back();

                int flags = state % 4;
                state /= 4;
                int rhsIndex = state % rhsCount;
                state /= rhsCount;
                int lhsIndex = state % lhsCount;
                int cellIndex = state / lhsCount;

                bool lhsGood = (flags & 2) != 0;
                bool witness = (flags & 1) != 0;
                bool lhsValue = lhsFirst ^ (lhsIndex % 2);
                bool rhsValue = rhsFirst ^ (rhsIndex % 2);
                Cell cell = cells[cellIndex];

                // At a point s, check the witness before lhs(s), because the
                // until condition reads lhs on [t,s). Inside an open piece,
                // lhs must also hold on the part immediately before s.
                if (cell.point)
                {
                    witness = witness || (cell.inWindow && lhsGood && rhsValue);
                }
                else
                {
                    witness = witness || (cell.inWindow && lhsGood && lhsValue && rhsValue);
                }
                lhsGood = lhsGood && lhsValue;

                if (cellIndex + 1 == cells.size())
                {
                    if (lhs[lhsFirst][lhsIndex] && rhs[rhsFirst][rhsIndex])
                    {
                        possible[witness] = true;
                    }
                }
                else
                {
                    int rhsAdvanceLimit =
                        cell.inWindow && cells[cellIndex + 1].inWindow ? 1 : 0;
                    for (int lhsAdvance = 0; lhsAdvance <= 1; lhsAdvance++)
                    {
                        for (int rhsAdvance = 0;
                             rhsAdvance <= rhsAdvanceLimit; rhsAdvance++)
                        {
                            push(cellIndex + 1,
                                 lhsIndex + lhsAdvance,
                                 rhsIndex + rhsAdvance,
                                 lhsGood,
                                 witness);
                        }
                    }
                }

                if (!cell.point)
                {
                    // Further changes may occur inside the same open piece. Two
                    // advances of one word represent an isolated point value.
                    int rhsAdvanceLimit = cell.inWindow ? 2 : 0;
                    for (int lhsAdvance = 0; lhsAdvance <= 2; lhsAdvance++)
                    {
                        for (int rhsAdvance = 0;
                             rhsAdvance <= rhsAdvanceLimit; rhsAdvance++)
                        {
                            if (lhsAdvance == 0 && rhsAdvance == 0)
                            {
                                continue;
                            }

                            int lhsChoices = lhsAdvance == 1 ? 2 : 1;
                            int rhsChoices = rhsAdvance == 1 ? 2 : 1;

                            for (int x = 0; x < lhsChoices; x++)
                            {
                                bool lhsPoint = lhsValue ^
                                    (lhsAdvance == 2 || (lhsAdvance == 1 && x == 1));

                                for (int y = 0; y < rhsChoices; y++)
                                {
                                    bool rhsPoint = rhsValue ^
                                        (rhsAdvance == 2 || (rhsAdvance == 1 && y == 1));

                                    bool witnessAtPoint = witness ||
                                        (cell.inWindow && lhsGood && rhsPoint);
                                    bool lhsGoodAtPoint = lhsGood && lhsPoint;

                                    push(cellIndex,
                                         lhsIndex + lhsAdvance,
                                         rhsIndex + rhsAdvance,
                                         lhsGoodAtPoint,
                                         witnessAtPoint);
                                }
                            }
                        }
                    }
                }

                if (possible[0] && possible[1])
                {
                    return possible;
                }
            }
        }
    }

    return possible;
}

template <std::size_t N>
vector<vector<bitset<N>>> bitsetTimedUntilLegacy(const vector<vector<bitset<N>>> &v1,
                                           const vector<vector<bitset<N>>> &v2,
                                           const vector<long long> &segmentation,
                                           const long long &a,
                                           const long long &b,
                                           const bool &upperInfinite,
                                           const bool &leftClosed,
                                           const bool &rightClosed,
                                           int s,
                                           int e)
{
    if (e == -1)
    {
        e = v1.size();
    }

    vector<vector<bitset<N>>> vv(v1.size());
    for (auto &v : vv)
    {
        v.resize(2);
    }

    vector<long long> windowOffsets{a};
    if (!upperInfinite)
    {
        windowOffsets.push_back(b);
    }
    sort(windowOffsets.begin(), windowOffsets.end());
    windowOffsets.erase(
        unique(windowOffsets.begin(), windowOffsets.end()),
        windowOffsets.end());

    vector<long long> allOffsets = windowOffsets;
    allOffsets.push_back(0);
    sort(allOffsets.begin(), allOffsets.end());
    allOffsets.erase(
        unique(allOffsets.begin(), allOffsets.end()), allOffsets.end());

    long long d = segmentation.back();
    TimedRange domain{0, 2 * (__int128)(d), true, false};

    // E_psi(T + offset) from the manuscript.
    auto changes = [&](const vector<vector<bitset<N>>> &v,
                       long long left, long long right, long long offset)
    {
        TimedRange range{2 * ((__int128)(left) + offset),
                         2 * ((__int128)(right) + offset), false, false};
        long long result = 0;

        for (int i = 1; i + 1 < segmentation.size(); i++)
        {
            if (timedContains(range, 2 * (__int128)(segmentation[i])))
            {
                result++;
            }
        }

        for (int i = 0; i < v.size(); i++)
        {
            TimedRange segment{2 * (__int128)(segmentation[i]),
                               2 * (__int128)(segmentation[i + 1]), true, false};
            if (!timedEmpty(timedIntersection(range, segment)))
            {
                result += max(msb(v[i][0]), msb(v[i][1]));
            }
        }

        return result;
    };

    for (int i = s; i < e; i++)
    {
        vector<long long> critical;
        for (const auto &endpoint : segmentation)
        {
            for (const auto &offset : allOffsets)
            {
                __int128 point = (__int128)(endpoint) - offset;
                if (point >= segmentation[i] && point < segmentation[i + 1])
                {
                    critical.push_back((long long)(point));
                }
            }
        }
        sort(critical.begin(), critical.end());
        critical.erase(unique(critical.begin(), critical.end()), critical.end());

        // (x,x) denotes the singleton {x}; (x,y) with x<y denotes (x,y).
        vector<pair<long long, long long>> placements;
        long long cursor = segmentation[i];
        for (const auto &point : critical)
        {
            if (cursor < point)
            {
                placements.push_back({cursor, point});
            }
            placements.push_back({point, point});
            cursor = point;
        }
        if (cursor < segmentation[i + 1])
        {
            placements.push_back({cursor, segmentation[i + 1]});
        }

        vector<bitset<N>> result(2);
        for (const auto &placement : placements)
        {
            bool point = placement.first == placement.second;
            __int128 t = point ? 2 * (__int128)(placement.first)
                               : (__int128)(placement.first) + placement.second;

            TimedRange horizon{t,
                               upperInfinite ? domain.right : t + 2 * (__int128)(b),
                               true,
                               !upperInfinite};
            TimedRange window{t + 2 * (__int128)(a),
                              upperInfinite ? domain.right : t + 2 * (__int128)(b),
                              leftClosed,
                              upperInfinite ? false : rightClosed};
            horizon = timedIntersection(horizon, domain);
            window = timedIntersection(window, domain);

            array<bool, 2> bits{true, false};
            if (!timedEmpty(window))
            {
                vector<bitset<N>> lhsProfile =
                    timedProfileLegacy(v1, segmentation, horizon);
                vector<bitset<N>> rhsProfile =
                    timedProfileLegacy(v2, segmentation, window);
                bits = possibleUntilBitsLegacy(
                    lhsProfile, rhsProfile, horizon, window);
            }

            vector<bitset<N>> region(2);
            if (bits[0] && bits[1])
            {
                long long bound = 1;
                if (!point)
                {
                    long long sum = 0;
                    for (const auto &offset : allOffsets)
                    {
                        sum += changes(v1, placement.first, placement.second, offset);
                    }
                    for (const auto &offset : windowOffsets)
                    {
                        sum += changes(v2, placement.first, placement.second, offset);
                    }
                    bound += 2 * sum;
                }

                if (bound > N)
                {
                    throw overflow_error("timed until exceeds the fixed bitset size");
                }
                for (int length = 0; length < bound; length++)
                {
                    region[0][length] = true;
                    region[1][length] = true;
                }
            }
            else
            {
                if (bits[0])
                {
                    region[0][0] = true;
                }
                if (bits[1])
                {
                    region[1][0] = true;
                }
            }

            result = bitsetConcatLegacy(result, region);
        }

        vv[i] = result;
    }

    return vv;
}

// Optimized timed-operator kernels. This file is included from monitoring.hpp
// after the legacy differential-oracle implementations.

enum class TimedDirection
{
    Future,
    Past
};

inline TimedRange reverseTimedRange(const TimedRange &range)
{
    return {
        -range.right,
        -range.left,
        range.rightClosed,
        range.leftClosed};
}

template <std::size_t N>
vector<bitset<N>> bitsetConcat(
    const vector<bitset<N>> &left,
    const vector<bitset<N>> &right)
{
    static_assert(N > 0, "timed languages require positive bitset capacity");
    if (left[0].none() && left[1].none())
    {
        return right;
    }

    vector<bitset<N>> output(2);
    const array<int, 2> rightMaximum{msb(right[0]), msb(right[1])};
    constexpr int scalarSpanLimit = 8;

    for (int first = 0; first < 2; first++)
    {
        const int leftMaximum = msb(left[first]);
        for (int leftIndex = 0; leftIndex <= leftMaximum; leftIndex++)
        {
            if (!left[first][leftIndex])
            {
                continue;
            }

            const int last = first ^ (leftIndex % 2);
            for (int second = 0; second < 2; second++)
            {
                const int maximum = rightMaximum[second];
                if (maximum < 0)
                {
                    continue;
                }

                const int shift = leftIndex + (last != second);
                if (maximum + shift >= static_cast<int>(N))
                {
                    throw overflow_error(
                        "value expression exceeds the fixed bitset size");
                }

                // A word shift is substantially cheaper for a broad right
                // bucket. Preserve the scalar path only for short active
                // prefixes where scanning at most eight indices wins.
                if (maximum + 1 > scalarSpanLimit)
                {
                    output[first] |= right[second] << shift;
                }
                else
                {
                    for (int rightIndex = 0; rightIndex <= maximum;
                         rightIndex++)
                    {
                        if (right[second][rightIndex])
                        {
                            output[first][shift + rightIndex] = true;
                        }
                    }
                }
            }
        }
    }

    return output;
}

inline pair<size_t, size_t> timedIntersectingSegments(
    const vector<long long> &segmentation,
    const TimedRange &range)
{
    const size_t segmentCount = segmentation.empty()
        ? 0
        : segmentation.size() - 1;
    if (segmentCount == 0 || timedEmpty(range))
    {
        return {0, 0};
    }

    // First segment whose open right endpoint lies strictly after range.left.
    size_t low = 1;
    size_t high = segmentation.size();
    while (low < high)
    {
        const size_t middle = low + (high - low) / 2;
        if (2 * (__int128)(segmentation[middle]) <= range.left)
        {
            low = middle + 1;
        }
        else
        {
            high = middle;
        }
    }
    if (low == segmentation.size())
    {
        return {segmentCount, segmentCount};
    }
    const size_t first = low - 1;

    // A closed range endpoint includes a segment beginning at that point; an
    // open endpoint does not.
    low = 0;
    high = segmentCount;
    while (low < high)
    {
        const size_t middle = low + (high - low) / 2;
        const __int128 point = 2 * (__int128)(segmentation[middle]);
        const bool included = range.rightClosed
            ? point <= range.right
            : point < range.right;
        if (included)
        {
            low = middle + 1;
        }
        else
        {
            high = middle;
        }
    }
    const size_t last = max(first, low);
    return {first, min(last, segmentCount)};
}

template <std::size_t N>
vector<bitset<N>> timedRestrictedSegment(
    const vector<bitset<N>> &value,
    const TimedRange &segment,
    const TimedRange &part)
{
    if (part.left == segment.left && part.right == segment.right &&
        part.leftClosed == segment.leftClosed &&
        part.rightClosed == segment.rightClosed)
    {
        return value;
    }
    if (part.left == part.right && part.left == segment.left)
    {
        return segmentFirstBit(value);
    }
    if (part.left == segment.left)
    {
        vector<bitset<N>> restricted = part.right == segment.right
            ? value
            : segmentPrefix(value);
        if (!part.leftClosed)
        {
            restricted = segmentWithoutFirstPoint(restricted);
        }
        return restricted;
    }
    if (part.right == segment.right)
    {
        return segmentSuffix(value);
    }
    return segmentInfix(value);
}

template <std::size_t N>
vector<bitset<N>> timedProfile(
    const vector<vector<bitset<N>>> &values,
    const vector<long long> &segmentation,
    const TimedRange &horizon)
{
    vector<bitset<N>> profile(2);
    const auto [first, last] = timedIntersectingSegments(
        segmentation, horizon);
    for (size_t i = first; i < last; i++)
    {
        const TimedRange segment{
            2 * (__int128)(segmentation[i]),
            2 * (__int128)(segmentation[i + 1]), true, false};
        const TimedRange part = timedIntersection(segment, horizon);
        if (timedEmpty(part))
        {
            continue;
        }
        profile = bitsetConcat(
            profile,
            timedRestrictedSegment(values[i], segment, part));
    }
    return profile;
}

// Whole-range products are associative under the monitor invariant that every
// factor language is nonempty. bitsetTimedUntil checks that invariant and
// falls back to the legacy fold before constructing this index otherwise.
template <std::size_t N>
class TimedProfileIndex
{
public:
    TimedProfileIndex(
        const vector<vector<bitset<N>>> &values,
        const vector<long long> &segmentation)
        : segmentation_(segmentation)
    {
        while (base_ < values.size())
        {
            base_ *= 2;
        }
        products_.resize(2 * base_);
        for (size_t i = 0; i < values.size(); i++)
        {
            products_[base_ + i] = valueProduct(values[i]);
        }
        for (size_t i = base_ - 1; i > 0; i--)
        {
            products_[i] = combine(
                products_[2 * i], products_[2 * i + 1]);
        }
    }

    vector<bitset<N>> query(const TimedRange &horizon) const
    {
        vector<bitset<N>> profile(2);
        const auto [first, last] = timedIntersectingSegments(
            segmentation_, horizon);
        if (first == last)
        {
            return profile;
        }

        profile = bitsetConcat(
            profile, restricted(first, horizon));
        if (first + 1 < last - 1)
        {
            append(profile, wholeRange(first + 1, last - 1));
        }
        if (first + 1 < last)
        {
            profile = bitsetConcat(
                profile, restricted(last - 1, horizon));
        }
        return profile;
    }

private:
    enum class ProductState
    {
        identity,
        value,
        overflow
    };

    struct Product
    {
        ProductState state = ProductState::identity;
        vector<bitset<N>> language;
    };

    static Product valueProduct(const vector<bitset<N>> &language)
    {
        return {ProductState::value, language};
    }

    static Product combine(const Product &left, const Product &right)
    {
        if (left.state == ProductState::overflow ||
            right.state == ProductState::overflow)
        {
            return {ProductState::overflow, {}};
        }
        if (left.state == ProductState::identity)
        {
            return right;
        }
        if (right.state == ProductState::identity)
        {
            return left;
        }
        try
        {
            return {
                ProductState::value,
                bitsetConcat(left.language, right.language)};
        }
        catch (const overflow_error &)
        {
            // Eager cache construction must not make an unused wide range
            // fail. Propagate the overflow only if a query selects it.
            return {ProductState::overflow, {}};
        }
    }

    static void append(
        vector<bitset<N>> &profile,
        const Product &product)
    {
        if (product.state == ProductState::overflow)
        {
            throw overflow_error(
                "value expression exceeds the fixed bitset size");
        }
        if (product.state == ProductState::value)
        {
            profile = bitsetConcat(profile, product.language);
        }
    }

    Product wholeRange(size_t first, size_t last) const
    {
        Product left;
        Product right;
        first += base_;
        last += base_;
        while (first < last)
        {
            if (first & 1U)
            {
                left = combine(left, products_[first++]);
            }
            if (last & 1U)
            {
                right = combine(products_[--last], right);
            }
            first /= 2;
            last /= 2;
        }
        return combine(left, right);
    }

    vector<bitset<N>> restricted(
        size_t segmentIndex,
        const TimedRange &horizon) const
    {
        const TimedRange segment{
            2 * (__int128)(segmentation_[segmentIndex]),
            2 * (__int128)(segmentation_[segmentIndex + 1]),
            true,
            false};
        return timedRestrictedSegment(
            products_[base_ + segmentIndex].language,
            segment,
            timedIntersection(segment, horizon));
    }

    const vector<long long> &segmentation_;
    size_t base_ = 1;
    vector<Product> products_;
};

inline vector<vector<long long>> timedCriticalPoints(
    const vector<long long> &segmentation,
    const vector<long long> &offsets,
    bool includeSegmentStarts,
    TimedDirection direction)
{
    if (segmentation.size() < 2)
    {
        return {};
    }

    const size_t segmentCount = segmentation.size() - 1;
    vector<vector<long long>> buckets(segmentCount);
    if (includeSegmentStarts)
    {
        for (size_t i = 0; i < segmentCount; i++)
        {
            buckets[i].push_back(segmentation[i]);
        }
    }

    // For a fixed offset the shifted endpoints remain sorted, so a single
    // monotone segment pointer buckets all events in O(S). Across k distinct
    // offsets this is O(kS), rather than one binary search per event.
    for (const long long offset : offsets)
    {
        size_t segment = 0;
        for (const long long endpoint : segmentation)
        {
            const __int128 shifted = direction == TimedDirection::Future
                ? (__int128)(endpoint) - (__int128)(offset)
                : (__int128)(endpoint) + (__int128)(offset);
            if (shifted < segmentation.front())
            {
                continue;
            }
            if (shifted >= segmentation.back())
            {
                break;
            }
            while (segment + 1 < segmentation.size() &&
                   shifted >= segmentation[segment + 1])
            {
                segment++;
            }
            const long long point = static_cast<long long>(shifted);
            buckets[segment].push_back(point);
        }
    }

    for (auto &bucket : buckets)
    {
        sort(bucket.begin(), bucket.end());
        bucket.erase(unique(bucket.begin(), bucket.end()), bucket.end());
    }
    return buckets;
}

inline vector<vector<long long>> timedFutureCriticalPoints(
    const vector<long long> &segmentation,
    const vector<long long> &offsets,
    bool includeSegmentStarts)
{
    return timedCriticalPoints(
        segmentation, offsets, includeSegmentStarts,
        TimedDirection::Future);
}

template <std::size_t N>
class TimedChangeIndex
{
public:
    TimedChangeIndex(
        const vector<vector<bitset<N>>> &values,
        const vector<long long> &segmentation)
        : segmentation_(segmentation), prefix_(values.size() + 1, 0)
    {
        for (size_t i = 0; i < values.size(); i++)
        {
            prefix_[i + 1] = prefix_[i] +
                max(msb(values[i][0]), msb(values[i][1]));
        }
    }

    long long changes(
        long long left,
        long long right,
        long long offset,
        TimedDirection direction) const
    {
        const __int128 lower = direction == TimedDirection::Future
            ? (__int128)(left) + (__int128)(offset)
            : (__int128)(left) - (__int128)(offset);
        const __int128 upper = direction == TimedDirection::Future
            ? (__int128)(right) + (__int128)(offset)
            : (__int128)(right) - (__int128)(offset);
        if (lower >= upper)
        {
            return 0;
        }

        size_t firstEndpoint = 1;
        size_t high = segmentation_.size() - 1;
        while (firstEndpoint < high)
        {
            const size_t middle = firstEndpoint +
                (high - firstEndpoint) / 2;
            if ((__int128)(segmentation_[middle]) <= lower)
            {
                firstEndpoint = middle + 1;
            }
            else
            {
                high = middle;
            }
        }

        size_t lastEndpoint = 1;
        high = segmentation_.size() - 1;
        while (lastEndpoint < high)
        {
            const size_t middle = lastEndpoint +
                (high - lastEndpoint) / 2;
            if ((__int128)(segmentation_[middle]) < upper)
            {
                lastEndpoint = middle + 1;
            }
            else
            {
                high = middle;
            }
        }
        const long long endpointCount = lastEndpoint > firstEndpoint
            ? static_cast<long long>(lastEndpoint - firstEndpoint)
            : 0;

        const TimedRange range{
            2 * lower, 2 * upper, false, false};
        const auto [firstSegment, lastSegment] =
            timedIntersectingSegments(segmentation_, range);
        return endpointCount +
            prefix_[lastSegment] - prefix_[firstSegment];
    }

    long long changes(
        long long left,
        long long right,
        long long offset) const
    {
        return changes(
            left, right, offset, TimedDirection::Future);
    }

private:
    const vector<long long> &segmentation_;
    vector<long long> prefix_;
};

inline void validateTimedArguments(
    size_t firstSize,
    optional<size_t> secondSize,
    const vector<long long> &segmentation,
    long long a,
    long long b,
    bool upperInfinite,
    int s,
    int &e)
{
    if (secondSize && *secondSize != firstSize)
    {
        throw invalid_argument("timed operands must have equal segment counts");
    }
    if (segmentation.size() != firstSize + 1)
    {
        throw invalid_argument(
            "timed segmentation must contain one more endpoint than segments");
    }
    for (size_t i = 1; i < segmentation.size(); i++)
    {
        if (segmentation[i - 1] >= segmentation[i])
        {
            throw invalid_argument(
                "timed segmentation endpoints must be strictly increasing");
        }
    }
    if (a < 0 || (!upperInfinite && b < a))
    {
        throw invalid_argument("invalid timed interval bounds");
    }
    if (e == -1)
    {
        e = static_cast<int>(firstSize);
    }
    if (s < 0 || e < s || static_cast<size_t>(e) > firstSize)
    {
        throw out_of_range("invalid timed segment range");
    }
}

inline int untilTrailingZeroBits(std::uint64_t word)
{
    int count = 0;
    while ((word & 1) == 0)
    {
        word >>= 1;
        count++;
    }
    return count;
}

template <std::size_t N>
class UntilPackedGrid
{
public:
    UntilPackedGrid(int rows, int columns)
        : rows_(rows), columns_(columns), words_((columns + 63) / 64),
          data_(static_cast<std::size_t>(rows) * 4 * words_, 0)
    {
    }

    std::uint64_t *bits(int row, int flags)
    {
        return data_.data() +
               (static_cast<std::size_t>(row) * 4 + flags) * words_;
    }

    const std::uint64_t *bits(int row, int flags) const
    {
        return data_.data() +
               (static_cast<std::size_t>(row) * 4 + flags) * words_;
    }

    int rows() const
    {
        return rows_;
    }

    int columns() const
    {
        return columns_;
    }

    int words() const
    {
        return words_;
    }

    void clear()
    {
        std::fill(data_.begin(), data_.end(), 0);
    }

    void swap(UntilPackedGrid &other)
    {
        data_.swap(other.data_);
    }

private:
    int rows_;
    int columns_;
    int words_;
    std::vector<std::uint64_t> data_;
};

template <std::size_t N>
std::array<bool, 2> possibleUntilBitsBitParallel(
    const std::vector<std::bitset<N>> &lhs,
    const std::vector<std::bitset<N>> &rhs,
    const TimedRange &horizon,
    const TimedRange &window)
{
    if (timedEmpty(window))
    {
        return {true, false};
    }

    struct Cell
    {
        bool point;
        bool inWindow;
    };

    std::vector<__int128> cuts{horizon.left, horizon.right,
                               window.left, window.right};
    std::sort(cuts.begin(), cuts.end());
    cuts.erase(std::unique(cuts.begin(), cuts.end()), cuts.end());

    std::vector<Cell> cells;
    cells.reserve(7);
    for (std::size_t index = 0; index < cuts.size(); ++index)
    {
        if (timedContains(horizon, cuts[index]))
        {
            cells.push_back({true, timedContains(window, cuts[index])});
        }
        if (index + 1 < cuts.size())
        {
            TimedRange open{cuts[index], cuts[index + 1], false, false};
            if (!timedEmpty(timedIntersection(open, horizon)))
            {
                cells.push_back({false,
                    !timedEmpty(timedIntersection(open, window))});
            }
        }
    }

    constexpr std::uint64_t evenPositions = 0x5555555555555555ULL;
    constexpr std::uint64_t oddPositions = 0xAAAAAAAAAAAAAAAAULL;
    constexpr int maximumWords = static_cast<int>((N + 63) / 64);
    std::array<bool, 2> possible{false, false};

    for (int lhsFirst = 0; lhsFirst < 2; ++lhsFirst)
    {
        const int lhsMax = msb(lhs[lhsFirst]);
        if (lhsMax < 0)
        {
            continue;
        }

        for (int rhsFirst = 0; rhsFirst < 2; ++rhsFirst)
        {
            const int rhsMax = msb(rhs[rhsFirst]);
            if (rhsMax < 0)
            {
                continue;
            }

            const int lhsCount = lhsMax + 1;
            const int rhsCount = rhsMax + 1;
            UntilPackedGrid<N> current(lhsCount, rhsCount);
            UntilPackedGrid<N> next(lhsCount, rhsCount);
            const int words = current.words();
            const int tailBits = rhsCount & 63;
            const std::uint64_t tailMask = tailBits == 0
                ? ~std::uint64_t{0}
                : (std::uint64_t{1} << tailBits) - 1;
            const std::uint64_t rhsOneMask = rhsFirst
                ? evenPositions : oddPositions;
            const std::uint64_t rhsZeroMask = ~rhsOneMask;

            std::array<std::uint64_t, maximumWords> rhsAccepted{};
            for (int index = 0; index < rhsCount; ++index)
            {
                if (rhs[rhsFirst][index])
                {
                    rhsAccepted[index / 64] |=
                        std::uint64_t{1} << (index & 63);
                }
            }

            auto trim = [&](std::uint64_t *bits)
            {
                bits[words - 1] &= tailMask;
            };

            auto orShiftMasked = [&](std::uint64_t *target,
                                     const std::uint64_t *source,
                                     std::uint64_t mask,
                                     int shift)
            {
                if (shift == 0)
                {
                    for (int word = 0; word < words; ++word)
                    {
                        target[word] |= source[word] & mask;
                    }
                }
                else
                {
                    // Descending order also makes the operation safe when
                    // source and target alias.
                    for (int word = words - 1; word >= 0; --word)
                    {
                        const std::uint64_t value = source[word] & mask;
                        target[word] |= value << shift;
                        if (word + 1 < words)
                        {
                            target[word + 1] |= value >> (64 - shift);
                        }
                    }
                }
                trim(target);
            };

            auto orSuffixAfter = [&](std::uint64_t *target,
                                     const std::uint64_t *source1,
                                     const std::uint64_t *source2)
            {
                int firstWord = -1;
                std::uint64_t firstValue = 0;
                for (int word = 0; word < words; ++word)
                {
                    firstValue = source1[word] |
                                 (source2 == nullptr ? 0 : source2[word]);
                    if (firstValue != 0)
                    {
                        firstWord = word;
                        break;
                    }
                }
                if (firstWord < 0)
                {
                    return;
                }

                const int firstBit = untilTrailingZeroBits(firstValue);
                if (firstBit != 63)
                {
                    target[firstWord] |=
                        ~std::uint64_t{0} << (firstBit + 1);
                }
                for (int word = firstWord + 1; word < words; ++word)
                {
                    target[word] = ~std::uint64_t{0};
                }
                trim(target);
            };

            auto horizontalClosure = [&](int lhsIndex, bool lhsValue,
                                         bool inWindow)
            {
                std::uint64_t *state0 = current.bits(lhsIndex, 0);
                std::uint64_t *state1 = current.bits(lhsIndex, 1);
                std::uint64_t *state2 = current.bits(lhsIndex, 2);
                std::uint64_t *state3 = current.bits(lhsIndex, 3);

                if (!inWindow)
                {
                    return;
                }

                if (!lhsValue)
                {
                    orSuffixAfter(state0, state0, state2);
                    orSuffixAfter(state1, state1, state3);
                }
                else
                {
                    orSuffixAfter(state0, state0, nullptr);
                    orSuffixAfter(state1, state1, nullptr);
                    orSuffixAfter(state3, state3, state2);
                    orShiftMasked(state2, state2, rhsZeroMask, 1);
                }
            };

            auto evaluatedFlags = [](int flags, bool lhsValue,
                                     bool rhsValue, bool point,
                                     bool inWindow)
            {
                bool lhsGood = (flags & 2) != 0;
                bool witness = (flags & 1) != 0;
                if (point)
                {
                    witness = witness ||
                        (inWindow && lhsGood && rhsValue);
                }
                else
                {
                    witness = witness ||
                        (inWindow && lhsGood && lhsValue && rhsValue);
                }
                lhsGood = lhsGood && lhsValue;
                return 2 * static_cast<int>(lhsGood) + witness;
            };

            auto pointFlags = [](int flags, bool lhsPoint,
                                 bool rhsPoint, bool inWindow)
            {
                bool lhsGood = (flags & 2) != 0;
                bool witness = (flags & 1) != 0;
                witness = witness ||
                    (inWindow && lhsGood && rhsPoint);
                lhsGood = lhsGood && lhsPoint;
                return 2 * static_cast<int>(lhsGood) + witness;
            };

            auto hasIntersection = [&](const std::uint64_t *bits)
            {
                for (int word = 0; word < words; ++word)
                {
                    if ((bits[word] & rhsAccepted[word]) != 0)
                    {
                        return true;
                    }
                }
                return false;
            };

            current.bits(0, 2)[0] = 1;

            for (std::size_t cellIndex = 0;
                 cellIndex < cells.size(); ++cellIndex)
            {
                const Cell cell = cells[cellIndex];
                const bool lastCell = cellIndex + 1 == cells.size();

                for (int lhsIndex = 0; lhsIndex < lhsCount; ++lhsIndex)
                {
                    const bool lhsValue = lhsFirst ^ (lhsIndex & 1);
                    if (!cell.point)
                    {
                        horizontalClosure(lhsIndex, lhsValue, cell.inWindow);
                    }

                    std::array<std::array<std::uint64_t, maximumWords>, 4>
                        evaluated{};
                    for (int flags = 0; flags < 4; ++flags)
                    {
                        const std::uint64_t *source =
                            current.bits(lhsIndex, flags);
                        for (int rhsValue = 0; rhsValue < 2; ++rhsValue)
                        {
                            const int targetFlags = evaluatedFlags(
                                flags, lhsValue, rhsValue,
                                cell.point, cell.inWindow);
                            const std::uint64_t mask = rhsValue
                                ? rhsOneMask : rhsZeroMask;
                            for (int word = 0; word < words; ++word)
                            {
                                evaluated[targetFlags][word] |=
                                    source[word] & mask;
                            }
                        }
                    }

                    if (lastCell)
                    {
                        if (lhs[lhsFirst][lhsIndex])
                        {
                            possible[0] = possible[0] ||
                                hasIntersection(evaluated[0].data()) ||
                                hasIntersection(evaluated[2].data());
                            possible[1] = possible[1] ||
                                hasIntersection(evaluated[1].data()) ||
                                hasIntersection(evaluated[3].data());
                        }
                    }
                    else
                    {
                        const int rhsAdvanceLimit =
                            cell.inWindow && cells[cellIndex + 1].inWindow
                                ? 1 : 0;
                        for (int flags = 0; flags < 4; ++flags)
                        {
                            for (int lhsAdvance = 0;
                                 lhsAdvance <= 1; ++lhsAdvance)
                            {
                                if (lhsIndex + lhsAdvance >= lhsCount)
                                {
                                    continue;
                                }
                                for (int rhsAdvance = 0;
                                     rhsAdvance <= rhsAdvanceLimit;
                                     ++rhsAdvance)
                                {
                                    orShiftMasked(
                                        next.bits(lhsIndex + lhsAdvance, flags),
                                        evaluated[flags].data(),
                                        ~std::uint64_t{0}, rhsAdvance);
                                }
                            }
                        }
                    }

                    if (!cell.point)
                    {
                        for (int lhsAdvance = 1;
                             lhsAdvance <= 2; ++lhsAdvance)
                        {
                            if (lhsIndex + lhsAdvance >= lhsCount)
                            {
                                continue;
                            }
                            const int lhsChoices =
                                lhsAdvance == 1 ? 2 : 1;
                            const int rhsAdvanceLimit =
                                cell.inWindow ? 2 : 0;
                            for (int rhsAdvance = 0;
                                 rhsAdvance <= rhsAdvanceLimit; ++rhsAdvance)
                            {
                                const int rhsChoices =
                                    rhsAdvance == 1 ? 2 : 1;
                                for (int lhsChoice = 0;
                                     lhsChoice < lhsChoices; ++lhsChoice)
                                {
                                    const bool lhsPoint = lhsValue ^
                                        (lhsAdvance == 2 ||
                                         (lhsAdvance == 1 && lhsChoice == 1));
                                    for (int rhsChoice = 0;
                                         rhsChoice < rhsChoices; ++rhsChoice)
                                    {
                                        for (int flags = 0;
                                             flags < 4; ++flags)
                                        {
                                            for (int rhsValue = 0;
                                                 rhsValue < 2; ++rhsValue)
                                            {
                                                const bool rhsPoint = rhsValue ^
                                                    (rhsAdvance == 2 ||
                                                     (rhsAdvance == 1 &&
                                                      rhsChoice == 1));
                                                const int targetFlags =
                                                    pointFlags(
                                                        flags, lhsPoint,
                                                        rhsPoint,
                                                        cell.inWindow);
                                                orShiftMasked(
                                                    current.bits(
                                                        lhsIndex + lhsAdvance,
                                                        targetFlags),
                                                    evaluated[flags].data(),
                                                    rhsValue ? rhsOneMask
                                                             : rhsZeroMask,
                                                    rhsAdvance);
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }

                if (possible[0] && possible[1])
                {
                    return possible;
                }
                if (!lastCell)
                {
                    current.swap(next);
                    next.clear();
                }
            }
        }
    }

    return possible;
}

template <std::size_t N>
array<bool, 2> possibleUntilBits(
    const vector<bitset<N>> &lhs,
    const vector<bitset<N>> &rhs,
    const TimedRange &horizon,
    const TimedRange &window)
{
    return possibleUntilBitsBitParallel(lhs, rhs, horizon, window);
}

template <std::size_t N>
vector<vector<bitset<N>>> bitsetTimedUntil(
    const vector<vector<bitset<N>>> &lhs,
    const vector<vector<bitset<N>>> &rhs,
    const vector<long long> &segmentation,
    long long a,
    long long b,
    bool upperInfinite,
    bool leftClosed,
    bool rightClosed,
    int s,
    int e)
{
    validateTimedArguments(
        lhs.size(), rhs.size(), segmentation,
        a, b, upperInfinite, s, e);

    for (size_t i = 0; i < lhs.size(); i++)
    {
        if (lhs[i].size() != 2 || rhs[i].size() != 2)
        {
            throw invalid_argument(
                "timed languages require exactly two buckets per segment");
        }
        if ((lhs[i][0].none() && lhs[i][1].none()) ||
            (rhs[i][0].none() && rhs[i][1].none()))
        {
            return bitsetTimedUntilLegacy(
                lhs, rhs, segmentation, a, b, upperInfinite,
                leftClosed, rightClosed, s, e);
        }
    }

    vector<vector<bitset<N>>> output(
        lhs.size(), vector<bitset<N>>(2));
    vector<long long> windowOffsets{a};
    if (!upperInfinite)
    {
        windowOffsets.push_back(b);
    }
    sort(windowOffsets.begin(), windowOffsets.end());
    windowOffsets.erase(
        unique(windowOffsets.begin(), windowOffsets.end()),
        windowOffsets.end());

    vector<long long> allOffsets = windowOffsets;
    allOffsets.push_back(0);
    sort(allOffsets.begin(), allOffsets.end());
    allOffsets.erase(
        unique(allOffsets.begin(), allOffsets.end()), allOffsets.end());

    const vector<vector<long long>> criticalPoints =
        timedCriticalPoints(
            segmentation, allOffsets, false, TimedDirection::Future);
    const TimedProfileIndex<N> lhsProfiles(lhs, segmentation);
    const TimedProfileIndex<N> rhsProfiles(rhs, segmentation);
    const TimedChangeIndex<N> lhsChanges(lhs, segmentation);
    const TimedChangeIndex<N> rhsChanges(rhs, segmentation);
    const TimedRange domain{
        0, 2 * (__int128)(segmentation.back()), true, false};

    for (int i = s; i < e; i++)
    {
        vector<pair<long long, long long>> placements;
        placements.reserve(2 * criticalPoints[i].size() + 1);
        long long cursor = segmentation[i];
        for (const long long point : criticalPoints[i])
        {
            if (cursor < point)
            {
                placements.push_back({cursor, point});
            }
            placements.push_back({point, point});
            cursor = point;
        }
        if (cursor < segmentation[i + 1])
        {
            placements.push_back({cursor, segmentation[i + 1]});
        }

        vector<bitset<N>> result(2);
        for (const auto &placement : placements)
        {
            const bool point = placement.first == placement.second;
            const __int128 t = point
                ? 2 * (__int128)(placement.first)
                : (__int128)(placement.first) + placement.second;
            TimedRange horizon{
                t,
                upperInfinite
                    ? domain.right
                    : t + 2 * (__int128)(b),
                true,
                !upperInfinite};
            TimedRange window{
                t + 2 * (__int128)(a),
                upperInfinite
                    ? domain.right
                    : t + 2 * (__int128)(b),
                leftClosed,
                upperInfinite ? false : rightClosed};
            horizon = timedIntersection(horizon, domain);
            window = timedIntersection(window, domain);

            array<bool, 2> bits{true, false};
            if (!timedEmpty(window))
            {
                const vector<bitset<N>> lhsProfile =
                    lhsProfiles.query(horizon);
                const vector<bitset<N>> rhsProfile =
                    rhsProfiles.query(window);
                bits = possibleUntilBits(
                    lhsProfile, rhsProfile, horizon, window);
            }

            vector<bitset<N>> region(2);
            if (bits[0] && bits[1])
            {
                long long bound = 1;
                if (!point)
                {
                    long long sum = 0;
                    for (const long long offset : allOffsets)
                    {
                        sum += lhsChanges.changes(
                            placement.first, placement.second, offset,
                            TimedDirection::Future);
                    }
                    for (const long long offset : windowOffsets)
                    {
                        sum += rhsChanges.changes(
                            placement.first, placement.second, offset,
                            TimedDirection::Future);
                    }
                    bound += 2 * sum;
                }
                if (bound > static_cast<long long>(N))
                {
                    throw overflow_error(
                        "timed until exceeds the fixed bitset size");
                }
                for (int length = 0; length < bound; length++)
                {
                    region[0][length] = true;
                    region[1][length] = true;
                }
            }
            else
            {
                region[0][0] = bits[0];
                region[1][0] = bits[1];
            }
            result = bitsetConcat(result, region);
        }
        output[i] = result;
    }
    return output;
}

template <std::size_t N>
vector<vector<bitset<N>>> bitsetBoundedUntil(vector<vector<bitset<N>>> v1,
                                             vector<vector<bitset<N>>> v2,
                                             vector<long long> segmentation,
                                             const long long &a,
                                             const long long &b,
                                             const bool &leftClosed,
                                             const bool &rightClosed,
                                             int s = 0,
                                             int e = -1)
{
    return bitsetTimedUntil(v1, v2, segmentation, a, b, false,
                            leftClosed, rightClosed, s, e);
}

template <std::size_t N>
vector<vector<bitset<N>>> bitsetTimedSinceLegacy(
    const vector<vector<bitset<N>>> &v1,
    const vector<vector<bitset<N>>> &v2,
    const vector<long long> &segmentation,
    const long long &a,
    const long long &b,
    const bool &upperInfinite,
    const bool &leftClosed,
    const bool &rightClosed,
    int s,
    int e)
{
    if (e == -1)
    {
        e = v1.size();
    }

    vector<vector<bitset<N>>> vv(v1.size());
    for (auto &v : vv)
    {
        v.resize(2);
    }

    vector<long long> windowOffsets{a};
    if (!upperInfinite)
    {
        windowOffsets.push_back(b);
    }
    sort(windowOffsets.begin(), windowOffsets.end());
    windowOffsets.erase(
        unique(windowOffsets.begin(), windowOffsets.end()),
        windowOffsets.end());

    vector<long long> allOffsets = windowOffsets;
    allOffsets.push_back(0);
    sort(allOffsets.begin(), allOffsets.end());
    allOffsets.erase(
        unique(allOffsets.begin(), allOffsets.end()), allOffsets.end());

    long long d = segmentation.back();
    TimedRange domain{0, 2 * (__int128)(d), true, false};

    // E_psi(T - offset) from the manuscript.
    auto changes = [&](const vector<vector<bitset<N>>> &v,
                       long long left, long long right, long long offset)
    {
        TimedRange range{2 * ((__int128)(left) - offset),
                         2 * ((__int128)(right) - offset), false, false};
        long long result = 0;

        for (int i = 1; i + 1 < segmentation.size(); i++)
        {
            if (timedContains(range, 2 * (__int128)(segmentation[i])))
            {
                result++;
            }
        }

        for (int i = 0; i < v.size(); i++)
        {
            TimedRange segment{2 * (__int128)(segmentation[i]),
                               2 * (__int128)(segmentation[i + 1]), true, false};
            if (!timedEmpty(timedIntersection(range, segment)))
            {
                result += max(msb(v[i][0]), msb(v[i][1]));
            }
        }

        return result;
    };

    auto reverseLanguage = [](vector<bitset<N>> language)
    {
        bitset<N> zero = language[0];
        bitset<N> one = language[1];
        language[0] = (zero & evenMask) | (one & oddMask);
        language[1] = (one & evenMask) | (zero & oddMask);
        return language;
    };

    for (int i = s; i < e; i++)
    {
        vector<long long> critical;
        for (const auto &endpoint : segmentation)
        {
            for (const auto &offset : allOffsets)
            {
                __int128 point = (__int128)(endpoint) + offset;
                if (point >= segmentation[i] && point < segmentation[i + 1])
                {
                    critical.push_back((long long)(point));
                }
            }
        }
        sort(critical.begin(), critical.end());
        critical.erase(unique(critical.begin(), critical.end()), critical.end());

        // (x,x) denotes the singleton {x}; (x,y) with x<y denotes (x,y).
        vector<pair<long long, long long>> placements;
        long long cursor = segmentation[i];
        for (const auto &point : critical)
        {
            if (cursor < point)
            {
                placements.push_back({cursor, point});
            }
            placements.push_back({point, point});
            cursor = point;
        }
        if (cursor < segmentation[i + 1])
        {
            placements.push_back({cursor, segmentation[i + 1]});
        }

        vector<bitset<N>> result(2);
        for (const auto &placement : placements)
        {
            bool point = placement.first == placement.second;
            __int128 t = point ? 2 * (__int128)(placement.first)
                               : (__int128)(placement.first) + placement.second;

            TimedRange horizon{upperInfinite ? domain.left : t - 2 * (__int128)(b),
                               t,
                               true,
                               true};
            TimedRange window{upperInfinite ? domain.left : t - 2 * (__int128)(b),
                              t - 2 * (__int128)(a),
                              upperInfinite ? true : rightClosed,
                              leftClosed};
            horizon = timedIntersection(horizon, domain);
            window = timedIntersection(window, domain);

            array<bool, 2> bits{true, false};
            if (!timedEmpty(window))
            {
                vector<bitset<N>> lhsProfile = reverseLanguage(
                    timedProfileLegacy(v1, segmentation, horizon));
                vector<bitset<N>> rhsProfile = reverseLanguage(
                    timedProfileLegacy(v2, segmentation, window));
                TimedRange reverseHorizon{-horizon.right,
                                          -horizon.left,
                                          horizon.rightClosed,
                                          horizon.leftClosed};
                TimedRange reverseWindow{-window.right,
                                         -window.left,
                                         window.rightClosed,
                                         window.leftClosed};
                bits = possibleUntilBitsLegacy(
                    lhsProfile, rhsProfile, reverseHorizon, reverseWindow);
            }

            vector<bitset<N>> region(2);
            if (bits[0] && bits[1])
            {
                long long bound = 1;
                if (!point)
                {
                    long long sum = 0;
                    for (const auto &offset : allOffsets)
                    {
                        sum += changes(v1, placement.first, placement.second, offset);
                    }
                    for (const auto &offset : windowOffsets)
                    {
                        sum += changes(v2, placement.first, placement.second, offset);
                    }
                    bound += 2 * sum;
                }

                if (bound > N)
                {
                    throw overflow_error("timed since exceeds the fixed bitset size");
                }
                for (int length = 0; length < bound; length++)
                {
                    region[0][length] = true;
                    region[1][length] = true;
                }
            }
            else
            {
                if (bits[0])
                {
                    region[0][0] = true;
                }
                if (bits[1])
                {
                    region[1][0] = true;
                }
            }

            result = bitsetConcatLegacy(result, region);
        }

        vv[i] = result;
    }

    return vv;
}

// Optimized timed non-strict Since kernel. This file is included from
// monitoring.hpp after the frozen scalar Since implementation and the shared
// timed-future infrastructure.

template <std::size_t N>
vector<vector<bitset<N>>> bitsetTimedSince(
    const vector<vector<bitset<N>>> &lhs,
    const vector<vector<bitset<N>>> &rhs,
    const vector<long long> &segmentation,
    const long long &a,
    const long long &b,
    const bool &upperInfinite,
    const bool &leftClosed,
    const bool &rightClosed,
    int s,
    int e)
{
    static_assert(N > 0, "timed since requires positive bitset capacity");
    validateTimedArguments(
        lhs.size(), rhs.size(), segmentation,
        a, b, upperInfinite, s, e);

    for (size_t segment = 0; segment < lhs.size(); segment++)
    {
        if (lhs[segment].size() != 2 || rhs[segment].size() != 2)
        {
            throw invalid_argument(
                "timed languages require exactly two buckets per segment");
        }
        // Cached concatenation products require nonempty factors. Preserve
        // the established behavior for valid callers outside that invariant.
        if ((lhs[segment][0].none() && lhs[segment][1].none()) ||
            (rhs[segment][0].none() && rhs[segment][1].none()))
        {
            return bitsetTimedSinceLegacy(
                lhs, rhs, segmentation, a, b, upperInfinite,
                leftClosed, rightClosed, s, e);
        }
    }

    vector<vector<bitset<N>>> output(
        lhs.size(), vector<bitset<N>>(2));
    vector<long long> windowOffsets{a};
    if (!upperInfinite)
    {
        windowOffsets.push_back(b);
    }
    sort(windowOffsets.begin(), windowOffsets.end());
    windowOffsets.erase(
        unique(windowOffsets.begin(), windowOffsets.end()),
        windowOffsets.end());

    vector<long long> allOffsets = windowOffsets;
    allOffsets.push_back(0);
    sort(allOffsets.begin(), allOffsets.end());
    allOffsets.erase(
        unique(allOffsets.begin(), allOffsets.end()), allOffsets.end());

    const vector<vector<long long>> criticalPoints = timedCriticalPoints(
        segmentation, allOffsets, false, TimedDirection::Past);
    const TimedProfileIndex<N> lhsProfiles(lhs, segmentation);
    const TimedProfileIndex<N> rhsProfiles(rhs, segmentation);
    const TimedChangeIndex<N> lhsChanges(lhs, segmentation);
    const TimedChangeIndex<N> rhsChanges(rhs, segmentation);
    const TimedRange domain{
        0, 2 * (__int128)(segmentation.back()), true, false};

    const auto reverseProfile = [](const vector<bitset<N>> &profile)
    {
        const array<bitset<N>, 2> reversed =
            reverseAlternatingSegment(profile[0], profile[1]);
        return vector<bitset<N>>{reversed[0], reversed[1]};
    };

    for (int segment = s; segment < e; segment++)
    {
        vector<pair<long long, long long>> placements;
        placements.reserve(2 * criticalPoints[segment].size() + 1);
        long long cursor = segmentation[segment];
        for (const long long point : criticalPoints[segment])
        {
            if (cursor < point)
            {
                placements.push_back({cursor, point});
            }
            placements.push_back({point, point});
            cursor = point;
        }
        if (cursor < segmentation[segment + 1])
        {
            placements.push_back({cursor, segmentation[segment + 1]});
        }

        vector<bitset<N>> result(2);
        for (const auto &placement : placements)
        {
            const bool point = placement.first == placement.second;
            const __int128 t = point
                ? 2 * (__int128)(placement.first)
                : (__int128)(placement.first) + placement.second;

            TimedRange horizon{
                upperInfinite
                    ? domain.left
                    : t - 2 * (__int128)(b),
                t,
                true,
                true};
            TimedRange window{
                upperInfinite
                    ? domain.left
                    : t - 2 * (__int128)(b),
                t - 2 * (__int128)(a),
                upperInfinite ? true : rightClosed,
                leftClosed};
            horizon = timedIntersection(horizon, domain);
            window = timedIntersection(window, domain);

            array<bool, 2> bits{true, false};
            if (!timedEmpty(window))
            {
                const vector<bitset<N>> lhsProfile = reverseProfile(
                    lhsProfiles.query(horizon));
                const vector<bitset<N>> rhsProfile = reverseProfile(
                    rhsProfiles.query(window));
                const TimedRange reversedHorizon =
                    reverseTimedRange(horizon);
                const TimedRange reversedWindow =
                    reverseTimedRange(window);
                bits = possibleUntilBits(
                    lhsProfile, rhsProfile,
                    reversedHorizon, reversedWindow);
            }

            vector<bitset<N>> region(2);
            if (bits[0] && bits[1])
            {
                __int128 bound = 1;
                if (!point)
                {
                    __int128 changes = 0;
                    for (const long long offset : allOffsets)
                    {
                        changes += lhsChanges.changes(
                            placement.first, placement.second, offset,
                            TimedDirection::Past);
                    }
                    for (const long long offset : windowOffsets)
                    {
                        changes += rhsChanges.changes(
                            placement.first, placement.second, offset,
                            TimedDirection::Past);
                    }
                    bound += 2 * changes;
                }
                if (bound > static_cast<__int128>(N))
                {
                    throw overflow_error(
                        "timed since exceeds the fixed bitset size");
                }
                for (size_t length = 0;
                     length < static_cast<size_t>(bound); length++)
                {
                    region[0].set(length);
                    region[1].set(length);
                }
            }
            else
            {
                region[0][0] = bits[0];
                region[1][0] = bits[1];
            }
            result = bitsetConcat(result, region);
        }
        output[segment] = result;
    }
    return output;
}

template <std::size_t N>
vector<vector<bitset<N>>> bitsetBoundedSince(vector<vector<bitset<N>>> v1,
                                             vector<vector<bitset<N>>> v2,
                                             vector<long long> segmentation,
                                             const long long &a,
                                             const long long &b,
                                             const bool &leftClosed,
                                             const bool &rightClosed,
                                             int s = 0,
                                             int e = -1)
{
    return bitsetTimedSince(v1, v2, segmentation, a, b, false,
                            leftClosed, rightClosed, s, e);
}


template <std::size_t N>
vector<vector<bitset<N>>> getProfiles(vector<vector<bitset<N>>> v, vector<long long> segmentation, const long long &t1, const long long &t2, const long long &a, const long long &b, const bool &leftClosed, const bool &rightClosed, bool forward = true)
{
    vector<vector<bitset<N>>> profiles;

    // to handle the cases where the window falls out of the signal
    segmentation.push_back(INT32_MAX);

    if (forward){
        // this assumes the interval is trimmed to fit the temporal domain (equivalently, signal is assumed to extend beyond its domain)
        v.push_back(v[v.size() - 1]);
    }
    else{
        // this assumes the signal is false
        vector<bitset<N>> temp(2);
        temp[0][0] = true;
        v.push_back(temp);
    }

    vector<vector<long long>> breakpoints(2);

    // TODO: check and improve
    long long last = segmentation[segmentation.size() - 1];
    breakpoints[0].push_back(min(t1 + a, last));
    int low0 = upper_bound(segmentation.begin(), segmentation.end(), min(t1 + a, last)) - segmentation.begin();
    int high0 = lower_bound(segmentation.begin(), segmentation.end(), min(t2 + a, last)) - segmentation.begin() - 1;
    for (int i = low0; i <= high0; i++)
    {
        breakpoints[0].push_back(segmentation[i] - 1);
        breakpoints[0].push_back(segmentation[i]);
    }
    if (segmentation[high0] < min(t2 + a, last))
    {
        breakpoints[0].push_back(min(t2 + a, last) - 1);
        breakpoints[0].push_back(min(t2 + a, last));
    }

    breakpoints[1].push_back(min(t1 + b, last));
    int low1 = upper_bound(segmentation.begin(), segmentation.end(), min(t1 + b, last)) - segmentation.begin();
    int high1 = lower_bound(segmentation.begin(), segmentation.end(), min(t2 + b, last)) - segmentation.begin() - 1;
    for (int i = low1; i <= high1; i++)
    {
        breakpoints[1].push_back(segmentation[i] - 1);
        breakpoints[1].push_back(segmentation[i]);
    }
    if (segmentation[high1] < min(t2 + b, last))
    {
        breakpoints[1].push_back(min(t2 + b, last)- 1);
        breakpoints[1].push_back(min(t2 + b, last));
    }

    // to ignore the last interval which includes the point t2+b
    int l0 = breakpoints[0].size() - 1;
    int l1 = breakpoints[1].size() - 1;

    int i0 = 0;
    int i1 = 0;

    while (i0 < l0 && i1 <= l1)
    {
        // find the relation of the current window to the segments, determine which actions to carry for the profile (prefix, suffix, etc)

        long long left = breakpoints[0][i0];
        long long right = breakpoints[1][i1];

        int xind = upper_bound(segmentation.begin() + low0 - 1, segmentation.begin() + high1 + 1, left) - segmentation.begin() - 1;
        int yind = upper_bound(segmentation.begin() + low0 - 1, segmentation.begin() + high1 + 1, right) - segmentation.begin() - 1;

        vector<bitset<N>> pr(2);

        if (yind - xind == 0)
        { // the two ends of the window belong to the same segment
            if (left == segmentation[xind])
            { // the window starts together with the segment
                // concat: prefix
                pr = bitsetConcat(pr, segmentPrefix(v[xind]));
            }
            else
            { // the window starts after the beginning of the segment
                // concat: infix
                pr = bitsetConcat(pr, segmentInfix(v[xind]));
            }
        }
        else
        { // the two ends of the window fall in different segments
            if (left == segmentation[xind])
            {
                // concat: entire segment
                pr = bitsetConcat(pr, v[xind]);
            }
            else
            {
                // concat: suffix
                pr = bitsetConcat(pr, segmentSuffix(v[xind]));
            }

            // concat: every full segment between two end points
            for (int xx = xind + 1; xx < yind; xx++)
            {
                pr = bitsetConcat(pr, v[xx]);
            }

            // TODO: fix below
            if (right == segmentation[yind] && rightClosed)
            {
                // concat: firstBit
                pr = bitsetConcat(pr, segmentFirstBit(v[yind]));
            }
            else if (right != segmentation[yind])
            {
                // concat: prefix
                pr = bitsetConcat(pr, segmentPrefix(v[yind]));
            }
        }

        profiles.push_back(pr);

        // slide the window
        if (breakpoints[0][i0 + 1] - breakpoints[0][i0] < breakpoints[1][i1 + 1] - breakpoints[1][i1])
        {
            breakpoints[1][i1] += (breakpoints[0][i0 + 1] - breakpoints[0][i0]);
            i0++;
        }
        else if (breakpoints[0][i0 + 1] - breakpoints[0][i0] > breakpoints[1][i1 + 1] - breakpoints[1][i1])
        {
            breakpoints[0][i0] += (breakpoints[1][i1 + 1] - breakpoints[1][i1]);
            i1++;
        }
        else
        {
            breakpoints[1][i1] += (breakpoints[0][i0 + 1] - breakpoints[0][i0]);
            i0++;
            i1++;
        }
    }

    return profiles;
}


template <std::size_t N>
vector<vector<bitset<N>>> getBackwardProfiles(vector<vector<bitset<N>>> v, vector<long long> segmentation, const long long &t1, const long long &t2, const long long &a, const long long &b, const bool &leftClosed, const bool &rightClosed)
{
    std::vector<long long> segM;
    segM.reserve(segmentation.size());
    for (auto it = segmentation.rbegin(); it != segmentation.rend(); ++it)
        segM.push_back(-(*it));

    std::vector<std::vector<std::bitset<N>>> vM(v.rbegin(), v.rend());
    for (auto& bs : vM) {
        auto temp0 = bs[0];
        auto temp1 = bs[1];
        bs[0] = (temp0 & evenMask) | (temp1 & oddMask);
        bs[1] = (temp1 & evenMask) | (temp0 & oddMask);
    }

    // auto profM = getProfiles<N>(v, segM, -t2, -t1, a, b, leftClosed, rightClosed, false);
    auto profM = getProfiles<N>(vM, segM, -t2, -t1, a, b, leftClosed, rightClosed, true);
    std::reverse(profM.begin(), profM.end());
    for (auto& bs : profM) {
        auto temp0 = bs[0];
        auto temp1 = bs[1];
        bs[0] = (temp0 & evenMask) | (temp1 & oddMask);
        bs[1] = (temp1 & evenMask) | (temp0 & oddMask);
    }
    // for (auto& pr : profM)
    //     std::reverse(pr.begin(), pr.end());

    return profM;

    // vector<vector<bitset<N>>> profiles;
    // vector<vector<long long>> breakpoints(2);

    // /* TODO: CHECK UPPERBOUND LOWERBOUND USAGE AND INDICES */

    // long long first = segmentation[0];
    // breakpoints[0].push_back(max(t1 - a, first));
    // int low0 = upper_bound(segmentation.begin(), segmentation.end(), max(t1 - a, first)) - segmentation.begin();
    // int high0 = lower_bound(segmentation.begin(), segmentation.end(), max(t2 - a, first)) - segmentation.begin() - 1;
    // for (int i = high0; i >= low0; i--)
    // {
    //     breakpoints[0].push_back(segmentation[i] + 1);
    //     breakpoints[0].push_back(segmentation[i]);
    // }
    // if (segmentation[low0] < max(t2 - a, first))
    // {
    //     breakpoints[0].push_back(max(t2 - a, first) + 1);
    //     breakpoints[0].push_back(max(t2 - a, first));
    // }

    // breakpoints[1].push_back(max(t1 - b, first));
    // int low1 = upper_bound(segmentation.begin(), segmentation.end(), max(t1 - b, first)) - segmentation.begin();
    // int high1 = lower_bound(segmentation.begin(), segmentation.end(), max(t2 - b, first)) - segmentation.begin() - 1;
    // for (int i = high1; i >= low1; i--)
    // {
    //     breakpoints[1].push_back(segmentation[i] + 1);
    //     breakpoints[1].push_back(segmentation[i]);
    // }
    // if (segmentation[low1] < max(t2 - b, first))
    // {
    //     breakpoints[1].push_back(max(t2 - b, first) + 1);
    //     breakpoints[1].push_back(max(t2 - b, first));
    // }


    // /* TODO: CHECK AND UPDATE BELOW */
    // // to ignore the last interval which includes the point t2-b
    // int l0 = breakpoints[0].size() - 1;
    // int l1 = breakpoints[1].size() - 1;

    // int i0 = l0 - 1;
    // int i1 = l1 - 1;

    // while (i0 > 0 && i1 >= 0)
    // {
    //     // find the relation of the current window to the segments, determine which actions to carry for the profile (prefix, suffix, etc)

    //     long long left = breakpoints[0][i0];
    //     long long right = breakpoints[1][i1];

    //     int xind = upper_bound(segmentation.begin() + low0 - 1, segmentation.begin() + high1 + 1, left) - segmentation.begin() - 1;
    //     int yind = upper_bound(segmentation.begin() + low0 - 1, segmentation.begin() + high1 + 1, right) - segmentation.begin() - 1;

    //     vector<bitset<N>> pr(2);

    //     if (yind - xind == 0)
    //     { // the two ends of the window belong to the same segment
    //         if (left == segmentation[xind])
    //         { // the window starts together with the segment
    //             // concat: prefix
    //             pr = bitsetConcat(pr, segmentPrefix(v[xind]));
    //         }
    //         else
    //         { // the window starts after the beginning of the segment
    //             // concat: infix
    //             pr = bitsetConcat(pr, segmentInfix(v[xind]));
    //         }
    //     }
    //     else
    //     { // the two ends of the window fall in different segments
    //         if (left == segmentation[xind])
    //         {
    //             // concat: entire segment
    //             pr = bitsetConcat(pr, v[xind]);
    //         }
    //         else
    //         {
    //             // concat: suffix
    //             pr = bitsetConcat(pr, segmentSuffix(v[xind]));
    //         }

    //         // concat: every full segment between two end points
    //         for (int xx = xind + 1; xx < yind; xx++)
    //         {
    //             pr = bitsetConcat(pr, v[xx]);
    //         }

    //         if (right == segmentation[yind] && rightClosed)
    //         {
    //             // concat: firstBit
    //             pr = bitsetConcat(pr, segmentFirstBit(v[yind]));
    //         }
    //         else if (right != segmentation[yind])
    //         {
    //             // concat: prefix
    //             pr = bitsetConcat(pr, segmentPrefix(v[yind]));
    //         }
    //     }

    //     profiles.push_back(pr);

    //     // slide the window
    //     if (breakpoints[0][i0 + 1] - breakpoints[0][i0] < breakpoints[1][i1 + 1] - breakpoints[1][i1])
    //     {
    //         breakpoints[1][i1] += (breakpoints[0][i0 + 1] - breakpoints[0][i0]);
    //         i0++;
    //     }
    //     else if (breakpoints[0][i0 + 1] - breakpoints[0][i0] > breakpoints[1][i1 + 1] - breakpoints[1][i1])
    //     {
    //         breakpoints[0][i0] += (breakpoints[1][i1 + 1] - breakpoints[1][i1]);
    //         i1++;
    //     }
    //     else
    //     {
    //         breakpoints[1][i1] += (breakpoints[0][i0 + 1] - breakpoints[0][i0]);
    //         i0++;
    //         i1++;
    //     }
    // }

    // return profiles;
}

template <std::size_t N>
vector<vector<bitset<N>>> bitsetAlwaysPastLegacy(
    vector<vector<bitset<N>>> v1,
    int s = 0,
    int e = -1,
    bool b0 = false,
    bool b1 = true)
{
    if (e == -1) {
        e = v1.size();
    }

    vector<vector<bitset<N>>> vv(v1.size());

    for (auto &v : vv)
    {
        v.resize(2);
    }

    bool lastBit0 = b0;
    bool lastBit1 = b1;

    for (int i = s; i < e; i++)
    {
        if (lastBit1 == true)
        {
            vector<int> len1(4);
            len1[0] = msb(v1[i][0] & evenMask) + 1; // 0...0
            len1[1] = msb(v1[i][0] & oddMask) + 1;  // 0...1
            len1[2] = msb(v1[i][1] & oddMask) + 1;  // 1...0
            len1[3] = msb(v1[i][1] & evenMask) + 1; // 1...1

            vv[i][0][0] = (max(len1[0], len1[1]) > 0);
            vv[i][1][1] = (max(len1[2], len1[3] - 2) > 0);
            vv[i][1][0] = v1[i][1][0];
        }

        if (lastBit0 == true)
        {
            vv[i][0][0] = true;
        }

        if ((vv[i][0][0] == true) || (vv[i][1][1] == true))
        {
            lastBit0 = true;
        }
        else
        {
            lastBit0 = false;
        }

        if (vv[i][1][0] == true)
        {
            lastBit1 = true;
        }
        else
        {
            lastBit1 = false;
        }
    }

    return vv;
}

template <std::size_t N>
void validateUntimedPastUnaryArguments(
    const vector<vector<bitset<N>>> &language,
    int &e,
    int s,
    const char *operatorName)
{
    static_assert(N > 0,
                  "untimed past operators require positive bitset capacity");
    if (e == -1)
    {
        e = static_cast<int>(language.size());
    }
    if (s < 0 || e < s || static_cast<size_t>(e) > language.size())
    {
        throw out_of_range(string("invalid untimed ") + operatorName +
                           " segment range");
    }
    for (const auto &segment : language)
    {
        if (segment.size() != 2)
        {
            throw invalid_argument(string("each untimed ") + operatorName +
                                   " segment must have two start buckets");
        }
    }
}

template <std::size_t N>
vector<vector<bitset<N>>> bitsetAlwaysPast(
    const vector<vector<bitset<N>>> &v1,
    int s = 0,
    int e = -1,
    bool b0 = false,
    bool b1 = true)
{
    validateUntimedPastUnaryArguments(v1, e, s, "always-past");

    vector<vector<bitset<N>>> result(
        v1.size(), vector<bitset<N>>(2));
    bool lastBit0 = b0;
    bool lastBit1 = b1;

    for (int i = s; i < e; i++)
    {
        if (lastBit1)
        {
            result[i][0][0] = v1[i][0].any();
            result[i][1][0] = v1[i][1][0];
            if constexpr (N > 1)
            {
                result[i][1][1] = (v1[i][1] >> 1).any();
            }
        }

        if (lastBit0)
        {
            result[i][0][0] = true;
        }

        lastBit0 = result[i][0][0] ||
            (N > 1 && result[i][1][1]);
        lastBit1 = result[i][1][0];
    }

    return result;
}

template <std::size_t N>
vector<vector<bitset<N>>> bitsetAlways(const vector<vector<bitset<N>>> &v1)
{
    vector<vector<bitset<N>>> vv(v1.size());
    auto nonSingletonEvenMask = evenMask;
    nonSingletonEvenMask[0] = false;

    for (auto &v : vv)
    {
        v.resize(2);
    }

    bool firstBit0 = false;
    bool firstBit1 = true;

    for (size_t i = v1.size(); i-- > 0;)
    {
        if (firstBit1 == true)
        {
            vv[i][0][0] = (v1[i][0] & evenMask).any() ||
                          (v1[i][1] & oddMask).any();
            vv[i][0][1] = (v1[i][0] & oddMask).any() ||
                          (v1[i][1] & nonSingletonEvenMask).any();
            vv[i][1][0] = v1[i][1][0];
        }

        if (firstBit0 == true)
        {
            vv[i][0][0] = true;
        }

        if ((vv[i][0][0] == true) || (vv[i][0][1] == true))
        {
            firstBit0 = true;
        }
        else
        {
            firstBit0 = false;
        }

        if (vv[i][1][0] == true)
        {
            firstBit1 = true;
        }
        else
        {
            firstBit1 = false;
        }
    }

    return vv;
}

template <std::size_t N>
vector<vector<bitset<N>>> bitsetEventuallyPastLegacy(
    const vector<vector<bitset<N>>> &v1,
    int s = 0,
    int e = -1,
    bool b0 = true,
    bool b1 = false)
{
    if (e == -1) {
        e = v1.size();
    }

    vector<vector<bitset<N>>> vv(v1.size());

    for (auto &v : vv)
    {
        v.resize(2);
    }

    bool lastBit0 = b0;
    bool lastBit1 = b1;

    for (int i = s; i < e; i++)
    {
        if (lastBit0 == true)
        {
            vector<int> len1(4);
            len1[0] = msb(v1[i][0] & evenMask) + 1; // 0...0
            len1[1] = msb(v1[i][0] & oddMask) + 1;  // 0...1
            len1[2] = msb(v1[i][1] & oddMask) + 1;  // 1...0
            len1[3] = msb(v1[i][1] & evenMask) + 1; // 1...1

            vv[i][0][0] = v1[i][0][0];
            vv[i][1][0] = (max(len1[2], len1[3]) > 0);
            vv[i][0][1] = (max(len1[0] - 2, len1[1]) > 0);
        }

        if (lastBit1 == true)
        {
            vv[i][1][0] = true;
        }

        if (vv[i][0][0] == true)
        {
            lastBit0 = true;
        }
        else
        {
            lastBit0 = false;
        }

        if ((vv[i][1][0] == true) || (vv[i][0][1] == true))
        {
            lastBit1 = true;
        }
        else
        {
            lastBit1 = false;
        }
    }

    return vv;
}

template <std::size_t N>
vector<vector<bitset<N>>> bitsetEventuallyPast(
    const vector<vector<bitset<N>>> &v1,
    int s = 0,
    int e = -1,
    bool b0 = true,
    bool b1 = false)
{
    validateUntimedPastUnaryArguments(v1, e, s, "eventually-past");

    vector<vector<bitset<N>>> result(
        v1.size(), vector<bitset<N>>(2));
    bool lastBit0 = b0;
    bool lastBit1 = b1;

    for (int i = s; i < e; i++)
    {
        if (lastBit0)
        {
            result[i][0][0] = v1[i][0][0];
            result[i][1][0] = v1[i][1].any();
            if constexpr (N > 1)
            {
                result[i][0][1] = (v1[i][0] >> 1).any();
            }
        }

        if (lastBit1)
        {
            result[i][1][0] = true;
        }

        lastBit0 = result[i][0][0];
        lastBit1 = result[i][1][0] ||
            (N > 1 && result[i][0][1]);
    }

    return result;
}

template <std::size_t N>
vector<vector<bitset<N>>> bitsetEventually(const vector<vector<bitset<N>>> &v1)
{
    vector<vector<bitset<N>>> vv(v1.size());
    auto nonSingletonEvenMask = evenMask;
    nonSingletonEvenMask[0] = false;

    for (auto &v : vv)
    {
        v.resize(2);
    }

    bool firstBit0 = true;
    bool firstBit1 = false;

    for (size_t i = v1.size(); i-- > 0;)
    {
        if (firstBit0 == true)
        {
            vv[i][0][0] = v1[i][0][0];
            vv[i][1][0] = (v1[i][0] & oddMask).any() ||
                          (v1[i][1] & evenMask).any();
            vv[i][1][1] = (v1[i][0] & nonSingletonEvenMask).any() ||
                          (v1[i][1] & oddMask).any();
        }

        if (firstBit1 == true)
        {
            vv[i][1][0] = true;
        }

        if (vv[i][0][0] == true)
        {
            firstBit0 = true;
        }
        else
        {
            firstBit0 = false;
        }

        if ((vv[i][1][0] == true) || (vv[i][1][1] == true))
        {
            firstBit1 = true;
        }
        else
        {
            firstBit1 = false;
        }
    }

    return vv;
}

// Shared geometry for timed eventually and timed always in both directions.
template <std::size_t N>
vector<vector<bitset<N>>> bitsetTimedUnaryLegacy(const vector<vector<bitset<N>>> &v1,
                                           const vector<long long> &segmentation,
                                           const long long &a,
                                           const long long &b,
                                           const bool &upperInfinite,
                                           const bool &leftClosed,
                                           const bool &rightClosed,
                                           const bool &always,
                                           const bool &past,
                                           int s,
                                           int e)
{
    if (e == -1)
    {
        e = v1.size();
    }

    vector<vector<bitset<N>>> vv(v1.size());
    for (auto &v : vv)
    {
        v.resize(2);
    }

    // Only the moving boundaries of the quantified window matter here.
    vector<long long> offsets{a};
    if (!upperInfinite)
    {
        offsets.push_back(b);
    }
    sort(offsets.begin(), offsets.end());
    offsets.erase(unique(offsets.begin(), offsets.end()), offsets.end());

    long long d = segmentation.back();
    TimedRange domain{0, 2 * (__int128)(d), true, false};

    auto changes = [&](long long left, long long right, long long offset)
    {
        TimedRange range{2 * (past ? (__int128)(left) - offset
                                   : (__int128)(left) + offset),
                         2 * (past ? (__int128)(right) - offset
                                   : (__int128)(right) + offset), false, false};
        long long result = 0;

        for (int i = 1; i + 1 < segmentation.size(); i++)
        {
            if (timedContains(range, 2 * (__int128)(segmentation[i])))
            {
                result++;
            }
        }

        for (int i = 0; i < v1.size(); i++)
        {
            TimedRange segment{2 * (__int128)(segmentation[i]),
                               2 * (__int128)(segmentation[i + 1]), true, false};
            if (!timedEmpty(timedIntersection(range, segment)))
            {
                result += max(msb(v1[i][0]), msb(v1[i][1]));
            }
        }

        return result;
    };

    for (int i = s; i < e; i++)
    {
        // Keep the segment's left endpoint as a singleton so every remaining
        // placement can use the existing open-interval representation.
        vector<long long> critical{segmentation[i]};
        for (const auto &endpoint : segmentation)
        {
            for (const auto &offset : offsets)
            {
                __int128 point = past ? (__int128)(endpoint) + offset
                                       : (__int128)(endpoint) - offset;
                if (point >= segmentation[i] && point < segmentation[i + 1])
                {
                    critical.push_back((long long)(point));
                }
            }
        }
        sort(critical.begin(), critical.end());
        critical.erase(unique(critical.begin(), critical.end()), critical.end());

        // (x,x) denotes the singleton {x}; (x,y) with x<y denotes (x,y).
        vector<pair<long long, long long>> placements;
        long long cursor = segmentation[i];
        for (const auto &point : critical)
        {
            if (cursor < point)
            {
                placements.push_back({cursor, point});
            }
            placements.push_back({point, point});
            cursor = point;
        }
        if (cursor < segmentation[i + 1])
        {
            placements.push_back({cursor, segmentation[i + 1]});
        }

        vector<bitset<N>> result(2);
        for (const auto &placement : placements)
        {
            bool point = placement.first == placement.second;
            __int128 t = point ? 2 * (__int128)(placement.first)
                               : (__int128)(placement.first) + placement.second;

            TimedRange window;
            if (past)
            {
                window = upperInfinite
                    ? TimedRange{domain.left, t - 2 * (__int128)(a),
                                 true, leftClosed}
                    : TimedRange{t - 2 * (__int128)(b),
                                 t - 2 * (__int128)(a),
                                 rightClosed, leftClosed};
            }
            else
            {
                window = TimedRange{t + 2 * (__int128)(a),
                                    upperInfinite ? domain.right
                                                  : t + 2 * (__int128)(b),
                                    leftClosed,
                                    upperInfinite ? false : rightClosed};
            }
            window = timedIntersection(window, domain);

            array<bool, 2> bits{false, false};
            if (timedEmpty(window))
            {
                if (always)
                {
                    bits[1] = true;
                }
                else
                {
                    bits[0] = true;
                }
            }
            else
            {
                vector<bitset<N>> profile = timedProfileLegacy(v1, segmentation, window);
                if (window.left == window.right)
                {
                    profile = segmentFirstBit(profile);
                }

                vector<vector<bitset<N>>> oneWindow(1, profile);
                vector<vector<bitset<N>>> untimed;
                if (past)
                {
                    untimed = always
                        ? bitsetAlwaysPastLegacy(oneWindow)
                        : bitsetEventuallyPastLegacy(oneWindow);
                    vector<bitset<N>> last = segmentLastBit(untimed[0]);
                    bits[0] = last[0][0];
                    bits[1] = last[1][0];
                }
                else
                {
                    untimed = always
                        ? bitsetAlways(oneWindow)
                        : bitsetEventually(oneWindow);
                    bits[0] = untimed[0][0].any();
                    bits[1] = untimed[0][1].any();
                }
            }

            vector<bitset<N>> region(2);
            if (bits[0] && bits[1])
            {
                long long bound = 1;
                if (!point)
                {
                    long long sum = 0;
                    for (const auto &offset : offsets)
                    {
                        sum += changes(placement.first, placement.second, offset);
                    }
                    bound += 2 * sum;
                }

                if (bound > N)
                {
                    throw overflow_error(always
                        ? "timed always exceeds the fixed bitset size"
                        : "timed eventually exceeds the fixed bitset size");
                }
                for (int length = 0; length < bound; length++)
                {
                    region[0][length] = true;
                    region[1][length] = true;
                }
            }
            else
            {
                if (bits[0])
                {
                    region[0][0] = true;
                }
                if (bits[1])
                {
                    region[1][0] = true;
                }
            }

            result = bitsetConcatLegacy(result, region);
        }

        vv[i] = result;
    }

    return vv;
}

// Optimized timed Always/Eventually/Historically/Once kernel. This file is
// included from monitoring.hpp after bitsetTimedUnaryLegacy.

struct TimedTruthSummary
{
    bool allZero = true;
    bool allOne = true;
    bool hasZero = false;
    bool hasOne = false;
};

inline void combineTimedTruthSummary(
    TimedTruthSummary &left,
    const TimedTruthSummary &right)
{
    left.allZero = left.allZero && right.allZero;
    left.allOne = left.allOne && right.allOne;
    left.hasZero = left.hasZero || right.hasZero;
    left.hasOne = left.hasOne || right.hasOne;
}

template <std::size_t N>
TimedTruthSummary summarizeTimedTruth(
    const vector<bitset<N>> &language)
{
    const size_t falseSingleton = language[0][0] ? 1 : 0;
    const size_t trueSingleton = language[1][0] ? 1 : 0;
    return {
        language[0][0],
        language[1][0],
        language[0].any() || language[1].count() > trueSingleton,
        language[1].any() || language[0].count() > falseSingleton};
}

template <std::size_t N>
class TimedTruthSummaryIndex
{
public:
    TimedTruthSummaryIndex(
        const vector<vector<bitset<N>>> &values,
        const vector<long long> &segmentation)
        : values_(values), segmentation_(segmentation)
    {
        for (auto &prefix : prefixes_)
        {
            prefix.assign(values.size() + 1, 0);
        }
        for (size_t i = 0; i < values.size(); i++)
        {
            const TimedTruthSummary summary = summarizeTimedTruth(values[i]);
            prefixes_[0][i + 1] = prefixes_[0][i] + summary.allZero;
            prefixes_[1][i + 1] = prefixes_[1][i] + summary.allOne;
            prefixes_[2][i + 1] = prefixes_[2][i] + summary.hasZero;
            prefixes_[3][i + 1] = prefixes_[3][i] + summary.hasOne;
        }
    }

    TimedTruthSummary query(const TimedRange &window) const
    {
        const auto [first, last] = timedIntersectingSegments(
            segmentation_, window);
        TimedTruthSummary result;
        if (first == last)
        {
            // A geometrically nonempty window can still lie before the first
            // supplied segment when retained history starts after time zero.
            // Match the empty materialized profile: neither truth value is
            // represented. The default summary remains the combine identity.
            return {false, false, false, false};
        }

        if (first + 1 == last)
        {
            combineTimedTruthSummary(
                result, restrictedSummary(first, window));
            return result;
        }

        combineTimedTruthSummary(
            result, restrictedSummary(first, window));
        if (first + 1 < last - 1)
        {
            combineTimedTruthSummary(
                result, fullRange(first + 1, last - 1));
        }
        combineTimedTruthSummary(
            result, restrictedSummary(last - 1, window));
        return result;
    }

private:
    TimedTruthSummary fullRange(size_t first, size_t last) const
    {
        const int length = static_cast<int>(last - first);
        const auto count = [&](size_t property)
        {
            return prefixes_[property][last] -
                   prefixes_[property][first];
        };
        return {
            count(0) == length,
            count(1) == length,
            count(2) != 0,
            count(3) != 0};
    }

    TimedTruthSummary restrictedSummary(
        size_t segmentIndex,
        const TimedRange &window) const
    {
        const TimedRange segment{
            2 * (__int128)(segmentation_[segmentIndex]),
            2 * (__int128)(segmentation_[segmentIndex + 1]),
            true,
            false};
        const TimedRange part = timedIntersection(segment, window);
        return summarizeTimedTruth(timedRestrictedSegment(
            values_[segmentIndex], segment, part));
    }

    const vector<vector<bitset<N>>> &values_;
    const vector<long long> &segmentation_;
    array<vector<int>, 4> prefixes_;
};

template <std::size_t N>
vector<vector<bitset<N>>> bitsetTimedUnary(
    const vector<vector<bitset<N>>> &values,
    const vector<long long> &segmentation,
    const long long &a,
    const long long &b,
    const bool &upperInfinite,
    const bool &leftClosed,
    const bool &rightClosed,
    const bool &always,
    const bool &past,
    int s,
    int e)
{
    static_assert(N > 0, "timed languages require positive bitset capacity");

    validateTimedArguments(
        values.size(), nullopt, segmentation,
        a, b, upperInfinite, s, e);
    for (const auto &segment : values)
    {
        if (segment.size() != 2)
        {
            throw invalid_argument(
                "timed languages require exactly two buckets per segment");
        }
        // The profile algebra assumes the valid monitor invariant that every
        // factor language is nonempty. Retain exact legacy behavior for test
        // or external callers that violate that invariant.
        if (segment[0].none() && segment[1].none())
        {
            return bitsetTimedUnaryLegacy(
                values, segmentation, a, b, upperInfinite,
                leftClosed, rightClosed, always, past, s, e);
        }
    }

    vector<vector<bitset<N>>> output(
        values.size(), vector<bitset<N>>(2));
    vector<long long> offsets{a};
    if (!upperInfinite)
    {
        offsets.push_back(b);
    }
    sort(offsets.begin(), offsets.end());
    offsets.erase(unique(offsets.begin(), offsets.end()), offsets.end());

    const TimedDirection direction = past
        ? TimedDirection::Past
        : TimedDirection::Future;
    const vector<vector<long long>> criticalPoints =
        timedCriticalPoints(segmentation, offsets, true, direction);
    const TimedTruthSummaryIndex<N> summaries(values, segmentation);
    const TimedChangeIndex<N> changeIndex(values, segmentation);
    const TimedRange domain{
        0, 2 * (__int128)(segmentation.back()), true, false};

    for (int i = s; i < e; i++)
    {
        vector<pair<long long, long long>> placements;
        placements.reserve(2 * criticalPoints[i].size() + 1);
        long long cursor = segmentation[i];
        for (const long long point : criticalPoints[i])
        {
            if (cursor < point)
            {
                placements.push_back({cursor, point});
            }
            placements.push_back({point, point});
            cursor = point;
        }
        if (cursor < segmentation[i + 1])
        {
            placements.push_back({cursor, segmentation[i + 1]});
        }

        vector<bitset<N>> result(2);
        for (const auto &placement : placements)
        {
            const bool point = placement.first == placement.second;
            const __int128 t = point
                ? 2 * (__int128)(placement.first)
                : (__int128)(placement.first) + placement.second;
            TimedRange window;
            if (past)
            {
                window = upperInfinite
                    ? TimedRange{
                        domain.left,
                        t - 2 * (__int128)(a),
                        true,
                        leftClosed}
                    : TimedRange{
                        t - 2 * (__int128)(b),
                        t - 2 * (__int128)(a),
                        rightClosed,
                        leftClosed};
            }
            else
            {
                window = {
                    t + 2 * (__int128)(a),
                    upperInfinite
                        ? domain.right
                        : t + 2 * (__int128)(b),
                    leftClosed,
                    upperInfinite ? false : rightClosed};
            }
            window = timedIntersection(window, domain);

            array<bool, 2> bits{false, false};
            if (timedEmpty(window))
            {
                bits[always ? 1 : 0] = true;
            }
            else
            {
                TimedTruthSummary summary;
                if (window.left == window.right)
                {
                    summary = summarizeTimedTruth(segmentFirstBit(
                        timedProfile(values, segmentation, window)));
                }
                else
                {
                    summary = summaries.query(window);
                }
                bits = always
                    ? array<bool, 2>{summary.hasZero, summary.allOne}
                    : array<bool, 2>{summary.allZero, summary.hasOne};
            }

            vector<bitset<N>> region(2);
            if (bits[0] && bits[1])
            {
                __int128 bound = 1;
                if (!point)
                {
                    __int128 sum = 0;
                    for (const long long offset : offsets)
                    {
                        sum += changeIndex.changes(
                            placement.first, placement.second, offset,
                            direction);
                    }
                    bound += 2 * sum;
                }
                if (bound > static_cast<__int128>(N))
                {
                    throw overflow_error(always
                        ? "timed always exceeds the fixed bitset size"
                        : "timed eventually exceeds the fixed bitset size");
                }
                for (size_t length = 0;
                     length < static_cast<size_t>(bound);
                     length++)
                {
                    region[0][length] = true;
                    region[1][length] = true;
                }
            }
            else
            {
                region[0][0] = bits[0];
                region[1][0] = bits[1];
            }
            result = bitsetConcat(result, region);
        }
        output[i] = result;
    }
    return output;
}

template <std::size_t N>
vector<vector<bitset<N>>> bitsetBoundedEventually(vector<vector<bitset<N>>> v1,
                                                  vector<long long> segmentation,
                                                  const long long &a,
                                                  const long long &b,
                                                  const bool &leftClosed,
                                                  const bool &rightClosed,
                                                  int s = 0,
                                                  int e = -1)
{
    return bitsetTimedUnary(v1, segmentation, a, b, false,
                            leftClosed, rightClosed, false, false, s, e);
}

template <std::size_t N>
vector<vector<bitset<N>>> bitsetBoundedAlways(vector<vector<bitset<N>>> v1,
                                              vector<long long> segmentation,
                                              const long long &a,
                                              const long long &b,
                                              const bool &leftClosed,
                                              const bool &rightClosed,
                                              int s = 0,
                                              int e = -1)
{
    return bitsetTimedUnary(v1, segmentation, a, b, false,
                            leftClosed, rightClosed, true, false, s, e);
}

template <std::size_t N>
vector<vector<bitset<N>>> bitsetUnboundedEventually(
    const vector<vector<bitset<N>>> &v1,
    const vector<long long> &segmentation,
    long long a,
    bool leftClosed,
    int s = 0,
    int e = -1)
{
    if (!(a == 0 && leftClosed))
    {
        return bitsetTimedUnary(v1, segmentation, a, 0, true,
                                leftClosed, false, false, false, s, e);
    }

    if (e == -1)
    {
        e = static_cast<int>(v1.size());
    }
    if (segmentation.size() != v1.size() + 1)
    {
        throw invalid_argument(
            "unbounded eventually segmentation does not match language");
    }
    if (s < 0 || e < s || static_cast<size_t>(e) > v1.size())
    {
        throw out_of_range("invalid unbounded eventually segment range");
    }
    for (const auto &segment : v1)
    {
        if (segment.size() != 2)
        {
            throw invalid_argument(
                "an unbounded eventually segment must have two buckets");
        }
    }

    vector<vector<bitset<N>>> full = bitsetEventually(v1);
    if (s == 0 && static_cast<size_t>(e) == v1.size())
    {
        return full;
    }
    for (int i = 0; i < s; i++)
    {
        full[i][0].reset();
        full[i][1].reset();
    }
    for (size_t i = static_cast<size_t>(e); i < full.size(); i++)
    {
        full[i][0].reset();
        full[i][1].reset();
    }
    return full;
}

template <std::size_t N>
vector<vector<bitset<N>>> bitsetUnboundedAlways(
    const vector<vector<bitset<N>>> &v1,
    const vector<long long> &segmentation,
    long long a,
    bool leftClosed,
    int s = 0,
    int e = -1)
{
    if (!(a == 0 && leftClosed))
    {
        return bitsetTimedUnary(v1, segmentation, a, 0, true,
                                leftClosed, false, true, false, s, e);
    }

    if (e == -1)
    {
        e = static_cast<int>(v1.size());
    }
    if (segmentation.size() != v1.size() + 1)
    {
        throw invalid_argument(
            "unbounded always segmentation does not match language");
    }
    if (s < 0 || e < s || static_cast<size_t>(e) > v1.size())
    {
        throw out_of_range("invalid unbounded always segment range");
    }
    for (const auto &segment : v1)
    {
        if (segment.size() != 2)
        {
            throw invalid_argument(
                "an unbounded always segment must have two buckets");
        }
    }

    vector<vector<bitset<N>>> full = bitsetAlways(v1);
    if (s == 0 && static_cast<size_t>(e) == v1.size())
    {
        return full;
    }
    for (int i = 0; i < s; i++)
    {
        full[i][0].reset();
        full[i][1].reset();
    }
    for (size_t i = static_cast<size_t>(e); i < full.size(); i++)
    {
        full[i][0].reset();
        full[i][1].reset();
    }
    return full;
}

template <std::size_t N>
vector<vector<bitset<N>>> bitsetBoundedEventuallyPast(vector<vector<bitset<N>>> v1,
                                                      vector<long long> segmentation,
                                                      const long long &a,
                                                      const long long &b,
                                                      const bool &leftClosed,
                                                      const bool &rightClosed,
                                                      int s = 0,
                                                      int e = -1)
{
    return bitsetTimedUnary(v1, segmentation, a, b, false,
                            leftClosed, rightClosed, false, true, s, e);
}

template <std::size_t N>
vector<vector<bitset<N>>> bitsetBoundedAlwaysPast(vector<vector<bitset<N>>> v1,
                                                  vector<long long> segmentation,
                                                  const long long &a,
                                                  const long long &b,
                                                  const bool &leftClosed,
                                                  const bool &rightClosed,
                                                  int s = 0,
                                                  int e = -1)
{
    return bitsetTimedUnary(v1, segmentation, a, b, false,
                            leftClosed, rightClosed, true, true, s, e);
}

template <std::size_t N>
vector<vector<bitset<N>>> bitsetUnboundedEventuallyPast(
    const vector<vector<bitset<N>>> &v1,
    const vector<long long> &segmentation,
    long long a,
    bool leftClosed,
    int s = 0,
    int e = -1)
{
    if (!(a == 0 && leftClosed))
    {
        return bitsetTimedUnary(v1, segmentation, a, 0, true,
                                leftClosed, false, false, true, s, e);
    }

    if (segmentation.size() != v1.size() + 1)
    {
        throw invalid_argument(
            "unbounded past eventually segmentation does not match language");
    }
    if (s < 0)
    {
        throw out_of_range("invalid unbounded past eventually segment range");
    }

    vector<vector<bitset<N>>> result = bitsetEventuallyPast(v1, 0, e);
    if (e == -1)
    {
        e = static_cast<int>(v1.size());
    }
    if (e < s)
    {
        throw out_of_range("invalid unbounded past eventually segment range");
    }
    for (int i = 0; i < s; i++)
    {
        result[i][0].reset();
        result[i][1].reset();
    }
    return result;
}

template <std::size_t N>
vector<vector<bitset<N>>> bitsetUnboundedAlwaysPast(
    const vector<vector<bitset<N>>> &v1,
    const vector<long long> &segmentation,
    long long a,
    bool leftClosed,
    int s = 0,
    int e = -1)
{
    if (!(a == 0 && leftClosed))
    {
        return bitsetTimedUnary(v1, segmentation, a, 0, true,
                                leftClosed, false, true, true, s, e);
    }

    if (segmentation.size() != v1.size() + 1)
    {
        throw invalid_argument(
            "unbounded past always segmentation does not match language");
    }
    if (s < 0)
    {
        throw out_of_range("invalid unbounded past always segment range");
    }

    vector<vector<bitset<N>>> result = bitsetAlwaysPast(v1, 0, e);
    if (e == -1)
    {
        e = static_cast<int>(v1.size());
    }
    if (e < s)
    {
        throw out_of_range("invalid unbounded past always segment range");
    }
    for (int i = 0; i < s; i++)
    {
        result[i][0].reset();
        result[i][1].reset();
    }
    return result;
}

// Strict future Until is intentionally excluded from the maintained monitor.
#if 0
template <std::size_t N>
vector<vector<bitset<N>>> bitsetUntilStrict(vector<vector<bitset<N>>> v1, vector<vector<bitset<N>>> v2, bool b0 = true, bool b1 = false, int s = 0, int e = -1)
{
    if (e == -1) {
        e = v1.size();
    }

    vector<vector<bitset<N>>> vv(v1.size());

    for (auto &v : vv)
    {
        v.resize(2);
    }

    bool firstBit0 = b0;
    bool firstBit1 = b1;

    for (int i = e - 1; i >= s; i--)
    {
        // TODO: improve this
        vector<bool> save(4);
        save[0] = v1[i][0][0];
        save[1] = v1[i][1][0];
        save[2] = v2[i][0][0];
        save[3] = v2[i][1][0];

        if (firstBit0 == true)
        {
            // handling the corner cases of 0Ux=0, 1Ux=Ex, xU0=0, xU1=x
            if (v1[i][0][0] == true)
            {
                vv[i][0][0] = true;
                v1[i][0][0] = false;
            }

            if (v1[i][1][0] == true)
            {
                vector<int> temp(4);
                temp[0] = msb(v2[i][0] & evenMask) + 1; // 0...0
                temp[1] = msb(v2[i][0] & oddMask) + 1;  // 0...1
                temp[2] = msb(v2[i][1] & oddMask) + 1;  // 1...0
                temp[3] = msb(v2[i][1] & evenMask) + 1; // 1...1

                vv[i][0][0] = vv[i][0][0] || v2[i][0][0];
                vv[i][1][0] = vv[i][1][0] || (max(temp[1], temp[3]) > 0);
                vv[i][1][1] = vv[i][1][1] || (max(temp[0] - 2, temp[2]) > 0);

                v1[i][1][0] = false;
            }

            if (v2[i][0][0] == true)
            {
                vv[i][0][0] = true;
                v2[i][0][0] = false;
            }

            if (v2[i][1][0] == true)
            {
                vv[i][0] = vv[i][0] | v1[i][0];
                vv[i][1] = vv[i][1] | v1[i][1];
                v2[i][1][0] = false;
            }

            if ((v1[i][0].none() && v1[i][1].none()) || (v2[i][0].none() && v2[i][1].none()))
            {
                goto FIRSTBIT1_STRICT;
            }

            vector<int> len1(4);
            len1[0] = msb(v1[i][0] & evenMask) + 1; // 0...0
            len1[1] = msb(v1[i][0] & oddMask) + 1;  // 0...1
            len1[2] = msb(v1[i][1] & oddMask) + 1;  // 1...0
            len1[3] = msb(v1[i][1] & evenMask) + 1; // 1...1

            vector<int> len2(4);
            len2[0] = msb(v2[i][0] & evenMask) + 1; // 0...0
            len2[1] = msb(v2[i][0] & oddMask) + 1;  // 0...1
            len2[2] = msb(v2[i][1] & oddMask) + 1;  // 1...0
            len2[3] = msb(v2[i][1] & evenMask) + 1; // 1...1

            vector<int> len_strong(4);
            if (len2[0] > 0)
            {
                len_strong[0] = max({len1[0], len1[1] + 1, len1[2] - 1, len1[3]});
            }
            else if (len2[1] > 0 && len2[2] > 0)
            {
                if (len1[1] > 0)
                {
                    len_strong[0] = max({len1[0], len1[1] + 1, len1[2] - 1});
                }
                else
                {
                len_strong[0] = max({len1[0], len1[2] - 1});
                }
            }
            else if (len2[1] > 0 && len2[2] == 0)
            {
                len_strong[0] = max({len1[0], len1[2] - 1});
            }
            else if (len2[1] == 0 && len2[2] > 0)
            {
                if (len1[1] > 0)
                {
                    len_strong[0] = max({len1[0], len1[1] + 1});
                }
                else
                {
                    len_strong[0] = len1[0];
                }
            }
            else if (len2[3] > 0)
            {
                len_strong[0] = len1[0];
            }
            len_strong[1] = (len2[1] > 0) ? max(len1[1], len1[3] - 1) : len1[1];
            len_strong[2] = (max(len2[0], len2[2]) > 0) ? max(len1[2], len1[3] + 1) : len1[2];
            len_strong[3] = len1[3];
            // len_strong[0] = len2[0];
            // len_strong[1] = len2[1];
            // len_strong[2] = (max(len1[2], len1[3]) > 0) ? max(len2[2], len2[0] - 1) : len2[2];
            // len_strong[3] = (max(len1[2], len1[3]) > 0) ? max(len2[3], len2[1] - 1) : len2[3];

            int ctr;

            // 0...0
            ctr = len_strong[0] - 1;
            while (ctr > -1)
            {
                vv[i][0][ctr] = true;
                ctr -= 2;
            }

            // 0...1
            ctr = len_strong[1] - 1;
            if (len2[1] > 0 || len2[3] > 0)
            {
                while (ctr > 0)
                {
                    vv[i][0][ctr] = true;
                    ctr -= 2;
                }
            }

            // 1...0
            ctr = len_strong[2] - 1;
            if (ctr > 0)
            {
                if (len2[0] > 0 || len2[1] >= 4 || len2[2] > 0 || len2[3] > 0)
                {
                    while (ctr > 0)
                    {
                        vv[i][1][ctr] = true;
                        ctr -= 2;
                    }
                }
                else if (len2[1] == 2)
                {
                    while (ctr > 0)
                    {
                        vv[i][1][ctr] = v1[i][1][ctr];
                        ctr -= 2;
                    }
                }
            }

            // 1...1
            ctr = len_strong[3] - 1;
            if (ctr > 0)
            {
                if (len2[3] >= 3 || len2[1] >= 4)
                {
                    while (ctr > 1)
                    {
                        vv[i][1][ctr] = true;
                        ctr -= 2;
                    }
                }
                else if (len2[1] == 2)
                {
                    while (ctr > 1)
                    {
                        vv[i][1][ctr] = v1[i][1][ctr];
                        ctr -= 2;
                    }
                }
            }
        }

        FIRSTBIT1_STRICT:
        if (firstBit1 == true)
        {
            v1[i][0][0] = save[0];
            v1[i][1][0] = save[1];
            v2[i][0][0] = save[2];
            v2[i][1][0] = save[3];

            // handling the corner cases of 0Ux=0, 1Ux=1, xU0=Ax, xU1=x
            if (v1[i][0][0] == true)
            {
                vv[i][0][0] = true;
                v1[i][0][0] = false;
            }

            if (v1[i][1][0] == true)
            {
                vv[i][1][0] = true;
                v1[i][1][0] = false;
            }

            if (v2[i][0][0] == true)
            {
                v1[i][0][0] = save[0];
                v1[i][1][0] = save[1];

                vector<int> temp(4);
                temp[0] = msb(v1[i][0] & evenMask) + 1; // 0...0
                temp[1] = msb(v1[i][0] & oddMask) + 1;  // 0...1
                temp[2] = msb(v1[i][1] & oddMask) + 1;  // 1...0
                temp[3] = msb(v1[i][1] & evenMask) + 1; // 1...1

                vv[i][0][0] = vv[i][0][0] || (max(temp[0], temp[2]) > 0);
                vv[i][0][1] = vv[i][0][1] || (max(temp[1], temp[3] - 2) > 0);
                vv[i][1][0] = vv[i][1][0] || v1[i][1][0];

                v2[i][0][0] = false;
                v1[i][0][0] = false;
                v1[i][1][0] = false;
            }

            if (v2[i][1][0] == true)
            {
                vv[i][0] = vv[i][0] | v1[i][0];
                vv[i][1] = vv[i][1] | v1[i][1];
                v2[i][1][0] = false;
            }

            if ((v1[i][0].none() && v1[i][1].none()) || (v2[i][0].none() && v2[i][1].none()))
            {
                goto UPDATE_FIRSTBITS_STRICT;
            }

            vector<int> len1(4);
            len1[0] = msb(v1[i][0] & evenMask) + 1; // 0...0
            len1[1] = msb(v1[i][0] & oddMask) + 1;  // 0...1
            len1[2] = msb(v1[i][1] & oddMask) + 1;  // 1...0
            len1[3] = msb(v1[i][1] & evenMask) + 1; // 1...1

            vector<int> len2(4);
            len2[0] = msb(v2[i][0] & evenMask) + 1; // 0...0
            len2[1] = msb(v2[i][0] & oddMask) + 1;  // 0...1
            len2[2] = msb(v2[i][1] & oddMask) + 1;  // 1...0
            len2[3] = msb(v2[i][1] & evenMask) + 1; // 1...1

            vector<int> len_weak(4);
            len_weak[0] = (max(len2[0], len2[1]) > 0) ? max(len1[0], len1[2] - 1) : len1[0];
            len_weak[1] = (max(len2[0], len2[1]) > 0) ? max(len1[1], len1[3] - 1) : len1[1];
            len_weak[2] = len1[2];
            len_weak[3] = len1[3];

            int ctr;

            // 0...0
            ctr = len_weak[0] - 1;
            while (ctr >= 0)
            {
                vv[i][0][ctr] = true;
                ctr -= 2;
            }

            // 0...1
            ctr = len_weak[1] - 1;
            while (ctr > 0)
            {
                vv[i][0][ctr] = true;
                ctr -= 2;
            }

            // 1...0
            ctr = len_weak[2] - 1;
            if (ctr > -1)
            {
                if (len2[0] > 0 || len2[1] >= 4 || len2[2] > 0 || len2[3] > 0)
                {
                    while (ctr > 0)
                    {
                        vv[i][1][ctr] = true;
                        ctr -= 2;
                    }
                }
                else if (len2[1] == 2)
                {
                    while (ctr > 0)
                    {
                        vv[i][1][ctr] = v1[i][1][ctr];
                        ctr -= 2;
                    }
                }
            }

            // 1...1
            ctr = len_weak[3] - 1;
            if (ctr > -1)
            {
                if (len2[0] > 0 || len2[1] >= 4 || len2[2] > 0 || len2[3] > 0)
                {
                    while (ctr > 0)
                    {
                        vv[i][1][ctr] = true;
                        ctr -= 2;
                    }
                }
                else if (len2[1] == 2)
                {
                    while (ctr > 0)
                    {
                        vv[i][1][ctr] = v1[i][1][ctr];
                        ctr -= 2;
                    }
                }
            }
        }

        UPDATE_FIRSTBITS_STRICT:
        if (vv[i][0].any())
        {
            firstBit0 = true;
        }
        else
        {
            firstBit0 = false;
        }

        if (vv[i][1].any())
        {
            firstBit1 = true;
        }
        else
        {
            firstBit1 = false;
        }
    }

    return vv;
}
#endif

// Frozen pre-refactor scalar implementation used only as an independent
// differential oracle for the optimized Until and Since paths.
template <std::size_t N>
vector<vector<bitset<N>>> bitsetUntilNonStrictLegacy(vector<vector<bitset<N>>> v1,
                                               vector<vector<bitset<N>>> v2,
                                               bool b0 = true,
                                               bool b1 = false,
                                               int s = 0,
                                               int e = -1)
{
    if (e == -1)
    {
        e = v1.size();
    }

    vector<vector<bitset<N>>> vv(v1.size());

    for (auto &v : vv)
    {
        v.resize(2);
    }

    bool firstBit0 = b0;
    bool firstBit1 = b1;

    for (int i = e - 1; i >= s; i--)
    {
        // TODO: improve this
        vector<bool> save(4);
        save[0] = v1[i][0][0];
        save[1] = v1[i][1][0];
        save[2] = v2[i][0][0];
        save[3] = v2[i][1][0];

        if (firstBit0 == true)
        {
            // handling the corner cases of 0 U x, 1 U x, x U 0, x U 1
            if (v1[i][0][0] == true)
            {
                vv[i][0] = vv[i][0] | v2[i][0];
                vv[i][1] = vv[i][1] | v2[i][1];
                v1[i][0][0] = false;
            }

            if (v1[i][1][0] == true)
            {
                vector<int> temp(4);
                temp[0] = msb(v2[i][0] & evenMask) + 1; // 0...0
                temp[1] = msb(v2[i][0] & oddMask) + 1;  // 0...1
                temp[2] = msb(v2[i][1] & oddMask) + 1;  // 1...0
                temp[3] = msb(v2[i][1] & evenMask) + 1; // 1...1

                vv[i][0][0] = vv[i][0][0] || v2[i][0][0];
                vv[i][1][0] = vv[i][1][0] || (max(temp[1], temp[3]) > 0);
                vv[i][1][1] = vv[i][1][1] || (max(temp[0] - 2, temp[2]) > 0);

                v1[i][1][0] = false;
            }

            if (v2[i][0][0] == true)
            {
                vv[i][0][0] = true;
                v2[i][0][0] = false;
            }

            if (v2[i][1][0] == true)
            {
                vv[i][1][0] = true;
                v2[i][1][0] = false;
            }

            if ((v1[i][0].none() && v1[i][1].none()) || (v2[i][0].none() && v2[i][1].none()))
            {
                goto FIRSTBIT1_NONSTRICT;
            }

            vector<int> len1(4);
            len1[0] = msb(v1[i][0] & evenMask) + 1; // 0...0
            len1[1] = msb(v1[i][0] & oddMask) + 1;  // 0...1
            len1[2] = msb(v1[i][1] & oddMask) + 1;  // 1...0
            len1[3] = msb(v1[i][1] & evenMask) + 1; // 1...1

            vector<int> len2(4);
            len2[0] = msb(v2[i][0] & evenMask) + 1; // 0...0
            len2[1] = msb(v2[i][0] & oddMask) + 1;  // 0...1
            len2[2] = msb(v2[i][1] & oddMask) + 1;  // 1...0
            len2[3] = msb(v2[i][1] & evenMask) + 1; // 1...1

            vector<int> len_strong(4);
            len_strong[0] = len2[0];
            len_strong[1] = len2[1];
            len_strong[2] = (max(len1[2], len1[3]) > 0) ? max(len2[2], len2[0] - 1) : len2[2];
            len_strong[3] = (max(len1[2], len1[3]) > 0) ? max(len2[3], len2[1] - 1) : len2[3];

            int ctr;

            // 0...0
            ctr = len_strong[0] - 1;
            if (ctr > -1)
            {
                if (len1[2] == 2 && len1[0] == 0 && len1[1] == 0 && len1[3] == 0)
                {
                    while (ctr > 0)
                    {
                        vv[i][0][ctr] = v2[i][0][ctr];
                        ctr -= 2;
                    }
                }
                else
                {
                    while (ctr > 0)
                    {
                        vv[i][0][ctr] = true;
                        ctr -= 2;
                    }
                }
            }

            // 0...1
            ctr = len_strong[1] - 1;
            if (ctr > -1)
            {
                if (len1[2] == 2 && len1[0] == 0 && len1[1] == 0 && len1[3] == 0)
                {
                    while (ctr > 0)
                    {
                        vv[i][0][ctr] = v2[i][0][ctr];
                        ctr -= 2;
                    }
                }
                else
                {
                    while (ctr > 0)
                    {
                        vv[i][0][ctr] = true;
                        ctr -= 2;
                    }
                }
            }

            // 1...0
            ctr = len_strong[2] - 1;
            while (ctr > 0)
            {
                vv[i][1][ctr] = true;
                ctr -= 2;
            }

            // 1...1
            ctr = len_strong[3] - 1;
            while (ctr >= 0)
            {
                vv[i][1][ctr] = true;
                ctr -= 2;
            }
        }

        FIRSTBIT1_NONSTRICT:
        if (firstBit1 == true)
        {
            v1[i][0][0] = save[0];
            v1[i][1][0] = save[1];
            v2[i][0][0] = save[2];
            v2[i][1][0] = save[3];

            // handling the corner cases of 0Ux=x, 1Ux=1, xU0=Ax, xU1=1
            if (v2[i][0][0] == true)
            {
                vector<int> temp(4);
                temp[0] = msb(v1[i][0] & evenMask) + 1; // 0...0
                temp[1] = msb(v1[i][0] & oddMask) + 1;  // 0...1
                temp[2] = msb(v1[i][1] & oddMask) + 1;  // 1...0
                temp[3] = msb(v1[i][1] & evenMask) + 1; // 1...1

                vv[i][0][0] = vv[i][0][0] || (max(temp[0], temp[2]) > 0);
                vv[i][0][1] = vv[i][0][1] || (max(temp[1], temp[3] - 2) > 0);
                vv[i][1][0] = vv[i][1][0] || v1[i][1][0];

                v2[i][0][0] = false;
            }

            if (v1[i][0][0] == true)
            {
                v2[i][0][0] = save[2];

                vv[i][0] = vv[i][0] | v2[i][0];
                vv[i][1] = vv[i][1] | v2[i][1];
                v1[i][0][0] = false;

                v2[i][0][0] = false;
            }

            if (v1[i][1][0] == true)
            {
                vv[i][1][0] = true;
                v1[i][1][0] = false;
            }

            if (v2[i][1][0] == true)
            {
                vv[i][1][0] = true;
                v2[i][1][0] = false;
            }

            if ((v1[i][0].none() && v1[i][1].none()) || (v2[i][0].none() && v2[i][1].none()))
            {
                goto UPDATE_FIRSTBITS_NONSTRICT;
            }

            vector<int> len1(4);
            len1[0] = msb(v1[i][0] & evenMask) + 1; // 0...0
            len1[1] = msb(v1[i][0] & oddMask) + 1;  // 0...1
            len1[2] = msb(v1[i][1] & oddMask) + 1;  // 1...0
            len1[3] = msb(v1[i][1] & evenMask) + 1; // 1...1

            vector<int> len2(4);
            len2[0] = msb(v2[i][0] & evenMask) + 1; // 0...0
            len2[1] = msb(v2[i][0] & oddMask) + 1;  // 0...1
            len2[2] = msb(v2[i][1] & oddMask) + 1;  // 1...0
            len2[3] = msb(v2[i][1] & evenMask) + 1; // 1...1

            vector<int> len_weak(4);
            len_weak[0] = (max(len1[0], len1[2]) > 0) ? (len2[0]) : 0;
            len_weak[1] = ((max(len1[1], len1[3]) > 0) && len2[0] > 0) ? (max(len2[0] + 1, len2[1])) : len2[1];

            if (len1[0] > 0 || len1[2] > 0)
            {
                if (len1[2] > 0 && len2[0] > 0)
                {
                    len_weak[2] = max(len2[2], len2[0] - 1);
                }
                else
                {
                    len_weak[2] = len2[2];
                }
            }
            else
            {
                len_weak[2] = 0;
            }

            if (len1[3] > 0)
            {
                if (len2[2] > 0)
                {
                    len_weak[3] = max({len2[0], len2[1] - 1, len2[2] + 1, len2[3]});
                }
                else
                {
                    len_weak[3] = max({len2[0], len2[1] - 1, len2[3]});
                }
            }
            else if (len1[2] > 0 && len1[1] > 0)
            {
                if (len2[2] > 0)
                {
                    len_weak[3] = max({len2[1] - 1, len2[2] + 1, len2[3]});
                }
                else
                {
                    len_weak[3] = max(len2[1] - 1, len2[3]);
                }
            }
            else if (len1[2] > 0)
            {
                len_weak[3] = max(len2[1] - 1, len2[3]);
            }
            else if (len1[1] > 0)
            {
                if (len2[2] > 0)
                {
                    len_weak[3] = max(len2[2] + 1, len2[3]);
                }
                else
                {
                    len_weak[3] = len2[3];
                }
            }
            else
            {
                len_weak[3] = len2[3];
            }

            int ctr;

            // 0...0
            ctr = len_weak[0] - 1;
            if (ctr > -1)
            {
                if (len1[2] == 2 && len1[0] == 0)
                {
                    while (ctr > 0)
                    {
                        vv[i][0][ctr] = v2[i][0][ctr];
                        ctr -= 2;
                    }
                }
                else
                {
                    while (ctr > 0)
                    {
                        vv[i][0][ctr] = true;
                        ctr -= 2;
                    }
                }
            }

            // 0...1
            ctr = len_weak[1] - 1;
            if (ctr > -1)
            {
                if (len1[2] == 2 && len1[0] == 0 && len1[1] == 0 && len1[3] == 0)
                {
                    while (ctr > 0)
                    {
                        vv[i][0][ctr] = v2[i][0][ctr];
                        ctr -= 2;
                    }
                }
                else
                {
                    while (ctr > 0)
                    {
                        vv[i][0][ctr] = true;
                        ctr -= 2;
                    }
                }
            }

            // 1...0
            ctr = len_weak[2] - 1;
            if (ctr > -1)
            {
                if (len1[2] == 2 && len1[0] == 0 && len2[2] == 0)
                {
                    while (ctr > 0)
                    {
                        vv[i][1][ctr] = v2[i][1][ctr];
                        ctr -= 2;
                    }
                    int cc = len2[0] - 2;
                    while (cc > 0)
                    {
                        vv[i][1][cc] = true;
                        cc -= 2;
                    }
                }
                else
                {
                    while (ctr > 0)
                    {
                        vv[i][1][ctr] = true;
                        ctr -= 2;
                    }
                }
            }

            // 1...1
            ctr = len_weak[3] - 1;
            while (ctr >= 0)
            {
                vv[i][1][ctr] = true;
                ctr -= 2;
            }
        }

        UPDATE_FIRSTBITS_NONSTRICT:
        if (vv[i][0].any())
        {
            firstBit0 = true;
        }
        else
        {
            firstBit0 = false;
        }

        if (vv[i][1].any())
        {
            firstBit1 = true;
        }
        else
        {
            firstBit1 = false;
        }
    }

    return vv;
}

template <std::size_t N, typename LeftBuckets, typename RightBuckets>
array<bitset<N>, 2> bitsetUntilNonStrictSegment(
    const LeftBuckets &v1,
    const RightBuckets &v2,
    const AlternatingEndpointSummary &left,
    const AlternatingEndpointSummary &right,
    bool firstBit0,
    bool firstBit1)
{
    static_assert(N > 0, "untimed until requires positive bitset capacity");
    array<bitset<N>, 2> result;
    if (!firstBit0 && !firstBit1)
    {
        return result;
    }

    const array<int, 4> &len1 = left.nonSingletonLengths;
    const array<int, 4> &len2 = right.nonSingletonLengths;
    const bool leftHasNonSingleton =
        any_of(len1.begin(), len1.end(), [](int length)
        {
            return length > 0;
        });
    const bool rightHasNonSingleton =
        any_of(len2.begin(), len2.end(), [](int length)
        {
            return length > 0;
        });

    auto validateLengths = [](const array<long long, 4> &lengths)
    {
        for (const long long length : lengths)
        {
            if (length > static_cast<long long>(N))
            {
                throw overflow_error(
                    "untimed until exceeds the fixed bitset size");
            }
        }
    };

    if (firstBit0)
    {
        if (left.singletonFalse)
        {
            result[0] |= v2[0];
            result[1] |= v2[1];
        }
        if (left.singletonTrue)
        {
            result[0][0] = result[0][0] || right.singletonFalse;
            result[1][0] = result[1][0] ||
                right.singletonTrue || len2[1] > 0 || len2[3] > 0;
            if (len2[0] > 0 || len2[2] > 0)
            {
                result[1].set(1);
            }
        }
        if (right.singletonFalse)
        {
            result[0][0] = true;
        }
        if (right.singletonTrue)
        {
            result[1][0] = true;
        }

        if (leftHasNonSingleton && rightHasNonSingleton)
        {
            array<long long, 4> strongLengths;
            strongLengths[0] = len2[0];
            strongLengths[1] = len2[1];
            strongLengths[2] = max(len1[2], len1[3]) > 0
                ? max(len2[2], len2[0] - 1)
                : len2[2];
            strongLengths[3] = max(len1[2], len1[3]) > 0
                ? max(len2[3], len2[1] - 1)
                : len2[3];
            validateLengths(strongLengths);

            int index = static_cast<int>(strongLengths[0]) - 1;
            if (len1[2] == 2 && len1[0] == 0 &&
                len1[1] == 0 && len1[3] == 0)
            {
                while (index > 0)
                {
                    result[0][index] = v2[0][index];
                    index -= 2;
                }
            }
            else
            {
                while (index > 0)
                {
                    result[0][index] = true;
                    index -= 2;
                }
            }

            index = static_cast<int>(strongLengths[1]) - 1;
            if (len1[2] == 2 && len1[0] == 0 &&
                len1[1] == 0 && len1[3] == 0)
            {
                while (index > 0)
                {
                    result[0][index] = v2[0][index];
                    index -= 2;
                }
            }
            else
            {
                while (index > 0)
                {
                    result[0][index] = true;
                    index -= 2;
                }
            }

            index = static_cast<int>(strongLengths[2]) - 1;
            while (index > 0)
            {
                result[1][index] = true;
                index -= 2;
            }

            index = static_cast<int>(strongLengths[3]) - 1;
            while (index >= 0)
            {
                result[1][index] = true;
                index -= 2;
            }
        }
    }

    if (firstBit1)
    {
        if (right.singletonFalse)
        {
            result[0][0] = result[0][0] ||
                left.singletonFalse || len1[0] > 0 || len1[2] > 0;
            if (len1[1] > 0 || len1[3] > 0)
            {
                result[0].set(1);
            }
            result[1][0] = result[1][0] || left.singletonTrue;
        }
        if (left.singletonFalse)
        {
            result[0] |= v2[0];
            result[1] |= v2[1];
        }
        if (left.singletonTrue || right.singletonTrue)
        {
            result[1][0] = true;
        }

        if (leftHasNonSingleton && rightHasNonSingleton)
        {
            array<long long, 4> weakLengths{};
            weakLengths[0] = max(len1[0], len1[2]) > 0
                ? len2[0]
                : 0;
            weakLengths[1] =
                max(len1[1], len1[3]) > 0 && len2[0] > 0
                ? static_cast<long long>(max(
                    static_cast<long long>(len2[0]) + 1,
                    static_cast<long long>(len2[1])))
                : len2[1];

            if (len1[0] > 0 || len1[2] > 0)
            {
                weakLengths[2] = len1[2] > 0 && len2[0] > 0
                    ? max(len2[2], len2[0] - 1)
                    : len2[2];
            }

            if (len1[3] > 0)
            {
                weakLengths[3] = len2[2] > 0
                    ? max({
                        static_cast<long long>(len2[0]),
                        static_cast<long long>(len2[1]) - 1,
                        static_cast<long long>(len2[2]) + 1,
                        static_cast<long long>(len2[3])})
                    : max({
                        static_cast<long long>(len2[0]),
                        static_cast<long long>(len2[1]) - 1,
                        static_cast<long long>(len2[3])});
            }
            else if (len1[2] > 0 && len1[1] > 0)
            {
                weakLengths[3] = len2[2] > 0
                    ? max({
                        static_cast<long long>(len2[1]) - 1,
                        static_cast<long long>(len2[2]) + 1,
                        static_cast<long long>(len2[3])})
                    : max(
                        static_cast<long long>(len2[1]) - 1,
                        static_cast<long long>(len2[3]));
            }
            else if (len1[2] > 0)
            {
                weakLengths[3] = max(
                    static_cast<long long>(len2[1]) - 1,
                    static_cast<long long>(len2[3]));
            }
            else if (len1[1] > 0)
            {
                weakLengths[3] = len2[2] > 0
                    ? max(
                        static_cast<long long>(len2[2]) + 1,
                        static_cast<long long>(len2[3]))
                    : len2[3];
            }
            else
            {
                weakLengths[3] = len2[3];
            }
            validateLengths(weakLengths);

            int index = static_cast<int>(weakLengths[0]) - 1;
            if (len1[2] == 2 && len1[0] == 0)
            {
                while (index > 0)
                {
                    result[0][index] = v2[0][index];
                    index -= 2;
                }
            }
            else
            {
                while (index > 0)
                {
                    result[0][index] = true;
                    index -= 2;
                }
            }

            index = static_cast<int>(weakLengths[1]) - 1;
            if (len1[2] == 2 && len1[0] == 0 &&
                len1[1] == 0 && len1[3] == 0)
            {
                while (index > 0)
                {
                    result[0][index] = v2[0][index];
                    index -= 2;
                }
            }
            else
            {
                while (index > 0)
                {
                    result[0][index] = true;
                    index -= 2;
                }
            }

            index = static_cast<int>(weakLengths[2]) - 1;
            if (len1[2] == 2 && len1[0] == 0 && len2[2] == 0)
            {
                while (index > 0)
                {
                    result[1][index] = v2[1][index];
                    index -= 2;
                }
                int secondary = len2[0] - 2;
                while (secondary > 0)
                {
                    result[1][secondary] = true;
                    secondary -= 2;
                }
            }
            else
            {
                while (index > 0)
                {
                    result[1][index] = true;
                    index -= 2;
                }
            }

            index = static_cast<int>(weakLengths[3]) - 1;
            while (index >= 0)
            {
                result[1][index] = true;
                index -= 2;
            }
        }
    }

    return result;
}


template <std::size_t N>
vector<vector<bitset<N>>> bitsetUntilNonStrict(
    const vector<vector<bitset<N>>> &v1,
    const vector<vector<bitset<N>>> &v2,
    bool b0 = true,
    bool b1 = false,
    int s = 0,
    int e = -1)
{
    static_assert(N > 0, "untimed until requires positive bitset capacity");
    if (v1.size() != v2.size())
    {
        throw invalid_argument(
            "untimed until operands must have equal segment counts");
    }
    if (e == -1)
    {
        e = static_cast<int>(v1.size());
    }
    if (s < 0 || e < s || static_cast<size_t>(e) > v1.size())
    {
        throw out_of_range("invalid untimed until segment range");
    }
    for (size_t segment = 0; segment < v1.size(); segment++)
    {
        if (v1[segment].size() != 2 || v2[segment].size() != 2)
        {
            throw invalid_argument(
                "each untimed until segment must have two start buckets");
        }
    }

    vector<vector<bitset<N>>> result(
        v1.size(), vector<bitset<N>>(2));
    bool firstBit0 = b0;
    bool firstBit1 = b1;

    for (int i = e - 1; i >= s; i--)
    {
        if (!firstBit0 && !firstBit1)
        {
            continue;
        }

        const AlternatingEndpointSummary left =
            summarizeAlternatingEndpoints(v1[i][0], v1[i][1]);
        const AlternatingEndpointSummary right =
            summarizeAlternatingEndpoints(v2[i][0], v2[i][1]);
        const array<bitset<N>, 2> segmentResult =
            bitsetUntilNonStrictSegment<N>(
                v1[i], v2[i], left, right, firstBit0, firstBit1);
        result[i][0] = segmentResult[0];
        result[i][1] = segmentResult[1];
        firstBit0 = segmentResult[0].any();
        firstBit1 = segmentResult[1].any();
    }

    return result;
}



template <std::size_t N>
vector<vector<bitset<N>>> bitsetUnboundedUntil(
    const vector<vector<bitset<N>>> &v1,
    const vector<vector<bitset<N>>> &v2,
    const vector<long long> &segmentation,
    long long a,
    bool leftClosed,
    int s = 0,
    int e = -1)
{
    if (!(a == 0 && leftClosed))
    {
        return bitsetTimedUntil(v1, v2, segmentation, a, 0, true,
                                leftClosed, false, s, e);
    }

    if (e == -1)
    {
        e = static_cast<int>(v1.size());
    }
    if (s < 0 || e < s || static_cast<size_t>(e) > v1.size())
    {
        throw out_of_range("invalid unbounded until segment range");
    }

    vector<vector<bitset<N>>> full = bitsetUntilNonStrict(v1, v2);
    if (s == 0 && static_cast<size_t>(e) == v1.size())
    {
        return full;
    }
    for (int i = 0; i < s; i++)
    {
        full[i][0].reset();
        full[i][1].reset();
    }
    for (size_t i = static_cast<size_t>(e); i < full.size(); i++)
    {
        full[i][0].reset();
        full[i][1].reset();
    }
    return full;
}

// Strict Since is intentionally excluded from the maintained monitor.
#if 0
template <std::size_t N>
vector<vector<bitset<N>>> bitsetSinceStrict(vector<vector<bitset<N>>> v1, vector<vector<bitset<N>>> v2, int s = 0, int e = -1, bool b0 = true, bool b1 = false)
{
    if (e == -1) {
        e = v1.size();
    }

    reverse(v1.begin(), v1.end());
    for (auto &e : v1) {
        bitset<N> zero = e[0];
        bitset<N> one = e[1];
        e[0] = (zero & evenMask) | (one & oddMask);
        e[1] = (one & evenMask) | (zero & oddMask);
    }

    reverse(v2.begin(), v2.end());
    for (auto &e : v2) {
        bitset<N> zero = e[0];
        bitset<N> one = e[1];
        e[0] = (zero & evenMask) | (one & oddMask);
        e[1] = (one & evenMask) | (zero & oddMask);
    }

    int ss = v1.size() - s;
    int ee = v1.size() - e;

    vector<vector<bitset<N>>> vv = bitsetUntilStrict(v1, v2, b0, b1, ee, ss);
    reverse(vv.begin(), vv.end());
    for (auto &e : vv) {
        bitset<N> zero = e[0];
        bitset<N> one = e[1];
        e[0] = (zero & evenMask) | (one & oddMask);
        e[1] = (one & evenMask) | (zero & oddMask);
    }

    return vv;
}
#endif

// Process since from left to right so the carry between segments is exactly
// last(previous since) AND last(previous left operand).
template <std::size_t N>
vector<vector<bitset<N>>> bitsetSinceNonStrictLegacy(
    vector<vector<bitset<N>>> v1,
    vector<vector<bitset<N>>> v2,
    int s = 0,
    int e = -1,
    bool b0 = true,
    bool b1 = false)
{
    if (e == -1)
    {
        e = v1.size();
    }

    vector<vector<bitset<N>>> vv(v1.size(), vector<bitset<N>>(2));

    auto reverseLanguage = [](vector<bitset<N>> language)
    {
        bitset<N> zero = language[0];
        bitset<N> one = language[1];
        language[0] = (zero & evenMask) | (one & oddMask);
        language[1] = (one & evenMask) | (zero & oddMask);
        return language;
    };

    bool carry0 = b0;
    bool carry1 = b1;

    for (int i = s; i < e; i++)
    {
        vector<vector<bitset<N>>> lhs(1, reverseLanguage(v1[i]));
        vector<vector<bitset<N>>> rhs(1, reverseLanguage(v2[i]));
        auto result = bitsetUntilNonStrictLegacy(
            lhs, rhs, carry0, carry1, 0, 1);
        vv[i] = reverseLanguage(result[0]);

        auto [since0, since1] = bitsetLastBits(vv[i]);
        auto [lhs0, lhs1] = bitsetLastBits(v1[i]);
        carry0 = since0 || lhs0;
        carry1 = since1 && lhs1;
    }

    return vv;
}

template <std::size_t N>
vector<vector<bitset<N>>> bitsetSinceNonStrict(
    const vector<vector<bitset<N>>> &v1,
    const vector<vector<bitset<N>>> &v2,
    int s = 0,
    int e = -1,
    bool b0 = true,
    bool b1 = false)
{
    static_assert(N > 0, "untimed since requires positive bitset capacity");
    if (v1.size() != v2.size())
    {
        throw invalid_argument(
            "untimed since operands must have equal segment counts");
    }
    if (e == -1)
    {
        e = static_cast<int>(v1.size());
    }
    if (s < 0 || e < s || static_cast<size_t>(e) > v1.size())
    {
        throw out_of_range("invalid untimed since segment range");
    }
    for (size_t segment = 0; segment < v1.size(); segment++)
    {
        if (v1[segment].size() != 2 || v2[segment].size() != 2)
        {
            throw invalid_argument(
                "each untimed since segment must have two start buckets");
        }
    }

    vector<vector<bitset<N>>> result(
        v1.size(), vector<bitset<N>>(2));
    bool carry0 = b0;
    bool carry1 = b1;

    for (int i = s; i < e; i++)
    {
        const array<bitset<N>, 2> reversedLeft =
            reverseAlternatingSegment(v1[i][0], v1[i][1]);
        const array<bitset<N>, 2> reversedRight =
            reverseAlternatingSegment(v2[i][0], v2[i][1]);
        const AlternatingEndpointSummary left = reverseAlternatingSummary(
            summarizeAlternatingEndpoints(v1[i][0], v1[i][1]));
        const AlternatingEndpointSummary right = reverseAlternatingSummary(
            summarizeAlternatingEndpoints(v2[i][0], v2[i][1]));
        const array<bitset<N>, 2> reversedResult =
            bitsetUntilNonStrictSegment<N>(
                reversedLeft, reversedRight, left, right, carry0, carry1);

        const array<bitset<N>, 2> segmentResult =
            reverseAlternatingSegment(reversedResult);
        result[i][0] = segmentResult[0];
        result[i][1] = segmentResult[1];

        const bool sinceLastFalse = reversedResult[0].any();
        const bool sinceLastTrue = reversedResult[1].any();
        const bool leftLastFalse = reversedLeft[0].any();
        const bool leftLastTrue = reversedLeft[1].any();
        carry0 = sinceLastFalse || leftLastFalse;
        carry1 = sinceLastTrue && leftLastTrue;
    }

    return result;
}

template <std::size_t N>
vector<vector<bitset<N>>> bitsetUnboundedSince(
    const vector<vector<bitset<N>>> &v1,
    const vector<vector<bitset<N>>> &v2,
    const vector<long long> &segmentation,
    long long a,
    bool leftClosed,
    int s = 0,
    int e = -1)
{
    if (!(a == 0 && leftClosed))
    {
        return bitsetTimedSince(v1, v2, segmentation, a, 0, true,
                                leftClosed, false, s, e);
    }

    if (s < 0)
    {
        throw out_of_range("invalid unbounded since segment range");
    }

    vector<vector<bitset<N>>> result =
        bitsetSinceNonStrict(v1, v2, 0, e);
    if (e == -1)
    {
        e = static_cast<int>(v1.size());
    }
    if (e < s)
    {
        throw out_of_range("invalid unbounded since segment range");
    }
    for (int i = 0; i < s; i++)
    {
        result[i][0].reset();
        result[i][1].reset();
    }
    return result;
}


template <std::size_t N>
vector<vector<bitset<N>>> bitsetConjunctionWithPolarity(
    const vector<vector<bitset<N>>> &v1,
    const vector<vector<bitset<N>>> &v2,
    bool negated)
{
    if (v1.size() != v2.size())
    {
        throw invalid_argument(
            "bitset binary operands must have equal segment counts");
    }
    for (size_t segment = 0; segment < v1.size(); segment++)
    {
        if (v1[segment].size() != 2 || v2[segment].size() != 2)
        {
            throw invalid_argument(
                "each bitset language segment must have two start buckets");
        }
    }

    vector<vector<bitset<N>>> result(
        v1.size(), vector<bitset<N>>(2));
    constexpr long long missingLength =
        numeric_limits<long long>::min() / 4;
    auto logicalBucket = [negated](
        const vector<bitset<N>> &language,
        size_t value) -> const bitset<N> &
    {
        return language[value ^ static_cast<size_t>(negated)];
    };

    for (size_t segment = 0; segment < v1.size(); segment++)
    {
        const AlternatingEndpointSummary left =
            summarizeAlternatingEndpoints(
                v1[segment][0], v1[segment][1], negated);
        const AlternatingEndpointSummary right =
            summarizeAlternatingEndpoints(
                v2[segment][0], v2[segment][1], negated);
        auto outputBucket = [&](size_t value) -> bitset<N> &
        {
            return result[segment][
                value ^ static_cast<size_t>(negated)];
        };

        if (left.singletonFalse || right.singletonFalse)
        {
            outputBucket(0)[0] = true;
        }
        if (left.singletonTrue)
        {
            outputBucket(0) |= logicalBucket(v2[segment], 0);
            outputBucket(1) |= logicalBucket(v2[segment], 1);
        }
        if (right.singletonTrue)
        {
            outputBucket(0) |= logicalBucket(v1[segment], 0);
            outputBucket(1) |= logicalBucket(v1[segment], 1);
        }

        array<long long, 4> leftLengths;
        array<long long, 4> rightLengths;
        for (size_t endpointClass = 0; endpointClass < 4;
             endpointClass++)
        {
            leftLengths[endpointClass] =
                left.nonSingletonLengths[endpointClass] == 0
                ? missingLength
                : left.nonSingletonLengths[endpointClass];
            rightLengths[endpointClass] =
                right.nonSingletonLengths[endpointClass] == 0
                ? missingLength
                : right.nonSingletonLengths[endpointClass];
        }

        array<long long, 4> lengths;
        lengths[0] = max({
            leftLengths[0] + rightLengths[0] - 2,
            leftLengths[0] + rightLengths[1] - 1,
            leftLengths[0] + rightLengths[2] - 1,
            leftLengths[0] + rightLengths[3],
            leftLengths[1] + rightLengths[0] - 1,
            leftLengths[1] + rightLengths[2],
            leftLengths[2] + rightLengths[0] - 1,
            leftLengths[2] + rightLengths[1],
            leftLengths[3] + rightLengths[0]}) - 1;
        lengths[1] = max({
            leftLengths[1] + rightLengths[1] - 1,
            leftLengths[1] + rightLengths[3],
            leftLengths[3] + rightLengths[1]}) - 1;
        lengths[2] = max({
            leftLengths[2] + rightLengths[2] - 1,
            leftLengths[2] + rightLengths[3],
            leftLengths[3] + rightLengths[2]}) - 1;
        lengths[3] = leftLengths[3] + rightLengths[3] - 1;

        for (const long long length : lengths)
        {
            if (length > static_cast<long long>(N))
            {
                throw length_error(
                    negated
                    ? "bitset disjunction result exceeds bitset capacity"
                    : "bitset conjunction result exceeds bitset capacity");
            }
        }

        for (size_t endpointClass = 0; endpointClass < 4;
             endpointClass++)
        {
            const long long lowestIndex = endpointClass == 3 ? 1 : 0;
            for (long long index = lengths[endpointClass] - 1;
                 index >= lowestIndex;
                 index -= 2)
            {
                outputBucket(endpointClass / 2)[
                    static_cast<size_t>(index)] = true;
            }
        }
    }

    return result;
}

template <std::size_t N>
vector<vector<bitset<N>>> bitsetConjunction(
    const vector<vector<bitset<N>>> &v1,
    const vector<vector<bitset<N>>> &v2)
{
    return bitsetConjunctionWithPolarity(v1, v2, false);
}

template <std::size_t N>
vector<vector<bitset<N>>> bitsetNegation(vector<vector<bitset<N>>> v)
{
    for (auto &vv : v)
    {
        swap(vv[0], vv[1]);
    }

    return v;
}

template <std::size_t N>
vector<vector<bitset<N>>> bitsetDisjunction(
    const vector<vector<bitset<N>>> &v1,
    const vector<vector<bitset<N>>> &v2)
{
    return bitsetConjunctionWithPolarity(v1, v2, true);
}




// this only works for boolean signals
template <std::size_t N>
vector<vector<vector<bitset<N>>>> convertIntoBitset(vector<vector<set<string>>> &valExprs)
{
    vector<vector<vector<bitset<N>>>> aps(valExprs.size());

    for (int i = 0; i < valExprs.size(); i++)
    {
        aps[i].resize(valExprs[i].size());

        for (int j = 0; j < valExprs[i].size(); j++)
        {
            aps[i][j].resize(2);

            for (auto s : valExprs[i][j])
            {
                if (!s.empty())
                {
                    string temp = s;
                    temp.erase(remove(temp.begin(), temp.end(), ';'), temp.end());

                    int k = (temp[0] == '0') ? 0 : 1;
                    int l = temp.length() - 1;
                    aps[i][j][k].set(l, true);
                }
            }
        }
    }

    return aps;
}

set<string> concatWithDestutter(const set<string> &set1, const set<string> &set2)
{
    if (set1.empty())
    {
        return set2;
    }

    set<string> result;

    for (const auto &s1 : set1)
    {
        for (const auto &s2 : set2)
        {
            if (!s1.empty() && !s2.empty())
            {
                if (s1.rfind(';') == string::npos && s2.find(';') == string::npos)
                {
                    if (s1 != s2)
                    {
                        result.insert(s1 + ";" + s2);
                    }
                    else
                    {
                        result.insert(s1);
                    }
                }
                else if (s1.rfind(';') == string::npos && s2.find(';') != string::npos)
                {
                    if (s1 != s2.substr(0, s2.find(';')))
                    {
                        result.insert(s1 + ";" + s2);
                    }
                    else
                    {
                        result.insert(s2);
                    }
                }
                else if (s1.rfind(';') != string::npos && s2.find(';') == string::npos)
                {
                    if (s1.substr(s1.rfind(';') + 1, s1.length()) != s2)
                    {
                        result.insert(s1 + ";" + s2);
                    }
                    else
                    {
                        result.insert(s1);
                    }
                }
                else
                {
                    if (s1.substr(s1.rfind(';') + 1, s1.length()) != s2.substr(0, s2.find(';')))
                    {
                        result.insert(s1 + ";" + s2);
                    }
                    else
                    {
                        result.insert(s1 + ";" + s2.substr(s2.find(';') + 1, s2.length()));
                    }
                }
            }
            else
            {
                if (!s1.empty())
                {
                    result.insert(s1);
                }
                else
                {
                    result.insert(s2);
                }
            }
        }
    }

    return result;
}

pair<int, int> convertIntoBoolWithDestutter(string str, string delimiter, double thr)
{
    vector<bool> v;

    if (!str.empty())
    {
        int start = 0;
        do
        {
            // Find the index of occurrence
            int idx = str.find(delimiter, start);
            if (idx == string::npos)
            {
                break;
            }

            int length = idx - start;
            if (stod(str.substr(start, length)) > thr)
            {
                if (v.empty() || v[v.size() - 1] != true)
                {
                    v.push_back(true);
                }
            }
            else
            {
                if (v.empty() || v[v.size() - 1] != false)
                {
                    v.push_back(false);
                }
            }
            start += (length + delimiter.size());
        } while (true);

        if (stod(str.substr(start, str.length())) > thr)
        {
            if (v.empty() || v[v.size() - 1] != true)
            {
                v.push_back(true);
            }
        }
        else
        {
            if (v.empty() || v[v.size() - 1] != false)
            {
                v.push_back(false);
            }
        }
    }

    int x = -1;
    int y = -1;

    if (!v.empty())
    {
        if (v[0] == true)
        {
            x = 1;
        }
        else
        {
            x = 0;
        }

        y = v.size() - 1;
    }

    return make_pair(x, y);
}

vector<pair<long long, double>> getData(const string &fileName, const int &len)
{
    vector<pair<long long, double>> signalReal;
    ifstream sigdata(fileName);
    string line;
    int ctr = 0;
    while (getline(sigdata, line) && ctr < len)
    {
        stringstream linestream(line);
        double t, v;
        linestream >> t >> v;
        signalReal.push_back(make_pair((long long)(t * 1000), v));
        ctr++;
    }
    sigdata.close();

    return signalReal;
}

vector<pair<long long, double>> getData3d(const string &fileName, const int &i)
{
    vector<pair<long long, double>> signalReal;
    ifstream sigdata(fileName);
    string line;
    while (getline(sigdata, line))
    {
        stringstream linestream(line);
        vector<double> temp(3);
        double t;
        linestream >> t >> temp[0] >> temp[1] >> temp[2];
        signalReal.push_back(make_pair((long long)(t * 1000), temp[i]));
    }
    sigdata.close();

    return signalReal;
}

vector<vector<pair<long long, double>>> convertSignalsToBool(const vector<vector<pair<long long, double>>> &signalsReal, int size)
{
    int n = signalsReal.size();
    vector<vector<pair<long long, double>>> signals(n);

    for (int i = 0; i < n; i++)
    {
        if (signalsReal[i][0].second > 0)
        {
            signals[i].push_back(make_pair(0, 1));
        }
        else
        {
            signals[i].push_back(make_pair(0, 0));
        }

        for (int j = 1; j < size; j++)
        {
            if (signalsReal[i][j].second > 0 && signalsReal[i][j - 1].second <= 0)
            {
                signals[i].push_back(make_pair(signalsReal[i][j].first, 1));
            }
            else if (signalsReal[i][j].second <= 0 && signalsReal[i][j - 1].second > 0)
            {
                signals[i].push_back(make_pair(signalsReal[i][j].first, 0));
            }
        }
    }

    return signals;
}

vector<vector<pair<long long, double>>> convertSignalsToBoolPartial(const vector<vector<pair<long long, double>>> &signalsReal, int start, int end)
{
    int n = signalsReal.size();
    vector<vector<pair<long long, double>>> signals(n);

    for (int i = 0; i < n; i++)
    {
        // if (signalsReal[i][start].second != signalsReal[i][max(0, start - 1)].second) {
        //     if (signalsReal[i][start].second > 0)
        //     {
        //         signals[i].push_back(make_pair(signalsReal[i][start].first, 1));
        //     }
        //     else
        //     {
        //         signals[i].push_back(make_pair(signalsReal[i][start].first, 0));
        //     }
        // }

        // else {
            if (signalsReal[i][start].second > 0)
            {
                signals[i].push_back(make_pair(signalsReal[i][start].first, 1));
            }
            else
            {
                signals[i].push_back(make_pair(signalsReal[i][start].first, 0));
            }
        // }

        for (int j = start + 1; j < end; j++)
        {
            if (signalsReal[i][j].second > 0 && signalsReal[i][j - 1].second <= 0)
            {
                signals[i].push_back(make_pair(signalsReal[i][j].first, 1));
            }
            else if (signalsReal[i][j].second <= 0 && signalsReal[i][j - 1].second > 0)
            {
                signals[i].push_back(make_pair(signalsReal[i][j].first, 0));
            }
        }
    }

    return signals;
}

void validateApproximateSignal(
    const vector<pair<long long, double>> &signal,
    long long duration)
{
    if (signal.empty())
    {
        throw invalid_argument("signals must contain an initial value");
    }

    long long previousTimestamp = -1;
    for (const auto &[timestamp, value] : signal)
    {
        if (timestamp < 0 || timestamp >= duration)
        {
            throw invalid_argument("signal timestamps must lie in [0,d)");
        }
        if (timestamp <= previousTimestamp)
        {
            throw invalid_argument("signal timestamps must be strictly increasing");
        }
        if (!isfinite(value))
        {
            throw invalid_argument("signal values must be finite");
        }
        previousTimestamp = timestamp;
    }
}

bool isActualEdge(
    const vector<pair<long long, double>> &signal,
    size_t index)
{
    return index > 0 && signal[index - 1].second != signal[index].second;
}

double signalValueAtUnchecked(
    const vector<pair<long long, double>> &signal,
    long long timestamp)
{
    double value = signal.front().second;
    for (size_t index = 1; index < signal.size(); ++index)
    {
        if (signal[index].first > timestamp)
        {
            break;
        }
        value = signal[index].second;
    }
    return value;
}

double signalValueAt(
    const vector<pair<long long, double>> &signal,
    long long timestamp)
{
    if (signal.empty())
    {
        throw invalid_argument("cannot evaluate an empty signal");
    }
    return signalValueAtUnchecked(signal, timestamp);
}

// Start from the paper's uncertainty region, clipped to the supplied analysis
// horizon, for every actual edge.
// Edge positions use integral timestamps, so strict happened-before order for
// consecutive edges of one signal is represented by separating coincident
// bounds by one representable tick. This is an ordering tie-break, not an
// application-level minimum-delay parameter. Empty entries represent the
// initial value or repeated samples, neither of which is an edge.
vector<vector<vector<long long>>> computeUncertaintyIntervals(
    const vector<vector<pair<long long, double>>> &signals,
    long long eps,
    long long duration,
    const set<int> &exactSignals = {},
    const vector<optional<long long>> &predecessorLowerBounds = {})
{
    if (eps <= 0)
    {
        throw invalid_argument("clock skew must be positive");
    }
    if (duration <= 0)
    {
        throw invalid_argument("signal duration must be positive");
    }
    for (int exactSignal : exactSignals)
    {
        if (exactSignal < 0 || exactSignal >= static_cast<int>(signals.size()))
        {
            throw invalid_argument("exact signal index is out of range");
        }
    }
    if (!predecessorLowerBounds.empty() &&
        predecessorLowerBounds.size() != signals.size())
    {
        throw invalid_argument(
            "signals and predecessor lower bounds are misaligned");
    }

    vector<vector<vector<long long>>> uncertainties(signals.size());
    for (size_t signalIndex = 0; signalIndex < signals.size(); ++signalIndex)
    {
        validateApproximateSignal(signals[signalIndex], duration);
        uncertainties[signalIndex].resize(signals[signalIndex].size());

        if (exactSignals.count(static_cast<int>(signalIndex)) != 0)
        {
            continue;
        }

        vector<size_t> actualEdges;
        for (size_t edge = 1; edge < signals[signalIndex].size(); ++edge)
        {
            if (isActualEdge(signals[signalIndex], edge))
            {
                actualEdges.push_back(edge);
            }
        }

        vector<long long> lowerBounds(signals[signalIndex].size());
        vector<long long> upperBounds(signals[signalIndex].size());
        constexpr long long strictOrderTick = 1;

        bool havePrevious = !predecessorLowerBounds.empty() &&
                            predecessorLowerBounds[signalIndex].has_value();
        long long previousLower = havePrevious
                                      ? *predecessorLowerBounds[signalIndex]
                                      : 0;
        if (havePrevious &&
            (previousLower < 0 || previousLower >= duration))
        {
            throw invalid_argument(
                "predecessor lower bound lies outside the horizon");
        }
        for (const size_t edge : actualEdges)
        {
            const long long timestamp = signals[signalIndex][edge].first;
            long long lower = eps >= timestamp ? 0 : timestamp - eps;
            if (havePrevious)
            {
                if (previousLower >= timestamp)
                {
                    throw logic_error(
                        "predecessor lower bound excludes an observed edge");
                }
                lower = max(lower, previousLower + strictOrderTick);
            }
            if (lower >= timestamp)
            {
                throw logic_error(
                    "ordered lower bound excludes its observed edge");
            }
            lowerBounds[edge] = lower;
            previousLower = lower;
            havePrevious = true;
        }

        bool haveNext = false;
        long long nextUpper = 0;
        for (auto edge = actualEdges.rbegin(); edge != actualEdges.rend(); ++edge)
        {
            const long long timestamp = signals[signalIndex][*edge].first;
            const long long remaining = duration - timestamp;
            long long upper = eps >= remaining ? duration : timestamp + eps;
            if (haveNext)
            {
                upper = min(upper, nextUpper - strictOrderTick);
            }
            upperBounds[*edge] = upper;
            nextUpper = upper;
            haveNext = true;
        }

        for (const size_t edge : actualEdges)
        {
            uncertainties[signalIndex][edge] = {
                lowerBounds[edge], upperBounds[edge]};
        }
    }
    return uncertainties;
}

vector<long long> computeCanonicalSegmentation(
    const vector<vector<pair<long long, double>>> &signals,
    const vector<vector<vector<long long>>> &uncertainties,
    long long duration,
    const set<int> &exactSignals = {})
{
    // The signals and regions are the aligned output of uncertainty construction.
    set<long long> endpoints{0, duration};
    for (size_t signalIndex = 0; signalIndex < signals.size(); ++signalIndex)
    {
        const bool exact = exactSignals.count(static_cast<int>(signalIndex)) != 0;
        for (size_t edge = 1; edge < signals[signalIndex].size(); ++edge)
        {
            if (!isActualEdge(signals[signalIndex], edge))
            {
                continue;
            }
            if (exact)
            {
                endpoints.insert(signals[signalIndex][edge].first);
                continue;
            }
            endpoints.insert(uncertainties[signalIndex][edge][0]);
            endpoints.insert(uncertainties[signalIndex][edge][1]);
        }
    }

    vector<long long> segmentation(endpoints.begin(), endpoints.end());
    if (segmentation.size() < 2 || segmentation.front() != 0 ||
        segmentation.back() != duration)
    {
        throw logic_error("invalid canonical segmentation");
    }
    return segmentation;
}

set<string> computeEdgeContributionUnchecked(
    double before,
    double after,
    const vector<long long> &region,
    long long segmentStart,
    long long segmentEnd)
{
    const long long lower = region[0];
    const long long upper = region[1];
    const string first = formatDouble(before);
    const string last = formatDouble(after);
    const string whole = destutterStr(first + ";" + last);

    if (segmentStart <= lower && upper <= segmentEnd)
    {
        return {whole};
    }
    if (segmentStart <= lower && segmentEnd < upper)
    {
        return {"", first, whole};
    }
    if (lower < segmentStart && upper <= segmentEnd)
    {
        return {"", last, whole};
    }
    return {"", first, last, whole};
}

set<string> computeEdgeContribution(
    double before,
    double after,
    const vector<long long> &region,
    long long segmentStart,
    long long segmentEnd)
{
    if (region.size() != 2)
    {
        throw invalid_argument("an uncertainty region must have two endpoints");
    }
    return computeEdgeContributionUnchecked(
        before, after, region, segmentStart, segmentEnd);
}

vector<vector<set<string>>> computeValueExpressions(
    const vector<vector<pair<long long, double>>> &signals,
    const vector<vector<vector<long long>>> &uncertainties,
    const vector<long long> &segmentation,
    const set<int> &exactSignals = {})
{
    // Inputs retain the validated alignment established by the pipeline.
    const size_t numSegments = segmentation.size() - 1;
    vector<vector<set<string>>> valExprs(
        signals.size(), vector<set<string>>(numSegments));

    for (size_t signalIndex = 0; signalIndex < signals.size(); ++signalIndex)
    {
        const bool exact = exactSignals.count(static_cast<int>(signalIndex)) != 0;
        for (size_t segment = 0; segment < numSegments; ++segment)
        {
            const long long start = segmentation[segment];
            const long long end = segmentation[segment + 1];

            if (exact)
            {
                string expression = formatDouble(
                    signalValueAtUnchecked(signals[signalIndex], start));
                for (size_t edge = 1; edge < signals[signalIndex].size(); ++edge)
                {
                    const long long timestamp = signals[signalIndex][edge].first;
                    if (timestamp <= start)
                    {
                        continue;
                    }
                    if (timestamp >= end)
                    {
                        break;
                    }
                    if (isActualEdge(signals[signalIndex], edge))
                    {
                        expression += ";" + formatDouble(
                            signals[signalIndex][edge].second);
                    }
                }
                valExprs[signalIndex][segment].insert(destutterStr(expression));
                continue;
            }

            bool intersects = false;
            set<string> expressions;
            for (size_t edge = 1; edge < signals[signalIndex].size(); ++edge)
            {
                if (!isActualEdge(signals[signalIndex], edge))
                {
                    continue;
                }
                const vector<long long> &region = uncertainties[signalIndex][edge];
                if (region[0] >= end || region[1] <= start)
                {
                    continue;
                }

                intersects = true;
                expressions = concatWithDestutter(
                    expressions,
                    computeEdgeContributionUnchecked(
                        signals[signalIndex][edge - 1].second,
                        signals[signalIndex][edge].second,
                        region,
                        start,
                        end));
            }

            if (intersects)
            {
                expressions.erase("");
            }
            else
            {
                expressions.insert(formatDouble(
                    signalValueAtUnchecked(signals[signalIndex], start)));
            }
            valExprs[signalIndex][segment] = move(expressions);
        }
    }
    return valExprs;
}

vector<vector<set<string>>> computeRelativeValueExpressions(
    const set<int> &leaders,
    const vector<vector<pair<long long, double>>> &signals,
    const vector<vector<vector<long long>>> &uncertainties,
    const vector<long long> &segmentation)
{
    return computeValueExpressions(signals, uncertainties, segmentation, leaders);
}
// vector<vector<set<string>>> adjustValueExpressions(const vector<vector<set<string>>> &valExprs, const vector<long long> &segmentation, const int &realStart, const int &realEnd)
// {
//     vector<vector<set<string>>> valExprsAdjusted = valExprs;

//     auto it = std::lower_bound(segmentation.begin(), segmentation.end(), realStart);
//     auto indStart = it - segmentation.begin();
//     it = std::upper_bound(segmentation.begin(), segmentation.end(), realEnd);
//     it--;
//     auto indEnd = it - segmentation.begin();

//     for (int i = 0; i < valExprs.size(); i++) {
//         for (int j = 0; j < indStart; j++) {
//             for (auto e : valExprs[i][j]) {
//                 // insert suffixes of e to valExprsAdjusted[i][j]

//             }
//         }

//         for (int j = indEnd; j < segmentation.size(); j++) {
//             for (auto e : valExprs[i][j]) {
//                 // insert prefixes of e to valExprsAdjusted[i][j]

//             }
//         }
//     }


//     return valExprsAdjusted;
// }

bitset<SIZE> generateBitmask(const int &parity)
{
    bitset<SIZE> mask;

    for (int i = parity % 2; i < SIZE; i = i + 2)
    {
        mask[i] = true;
    }

    return mask;
}

vector<vector<vector<bitset<SIZE>>>> convertSignalsToAtomicPropositions(const vector<vector<set<string>>> &valExprs, const double &thr)
{
        int n = valExprs.size();
        int numSegments = valExprs[0].size();
        vector<vector<vector<bitset<SIZE>>>> aps(n);

        for (auto &v : aps)
        {
            v.resize(numSegments);
            for (auto &vv : v)
            {
                vv.resize(2);
            }
        }

        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < numSegments; j++)
            {
                for (auto expr : valExprs[i][j])
                {
                    if (expr != "")
                    {
                        pair<int, int> xy = convertIntoBoolWithDestutter(expr, ";", thr);
                        if (xy.first >= 0 && xy.second >= 0)
                        {
                            aps[i][j][xy.first][xy.second] = true;
                        }
                    }
                }
            }
        }

    return aps;
}


vector<vector<vector<bitset<SIZE>>>> adjustAtomicPropositions(const vector<vector<vector<bitset<SIZE>>>>  &aps, const vector<long long> &segmentation, const int &realStart, const int &realEnd)
{
    auto it = std::lower_bound(segmentation.begin(), segmentation.end(), realStart);
    int indStart = it - segmentation.begin();
    it = std::upper_bound(segmentation.begin(), segmentation.end(), realEnd);
    it--;
    int indEnd = it - segmentation.begin();

    // if (indStart >= indEnd) {
    //     return aps;
    // }

    vector<vector<vector<bitset<SIZE>>>> apsAdjusted = aps;
    int s = aps.size();

    for (int i = 0; i < s; i++) {
        int t = aps[i].size();
        // for (int j = 0; j < indStart; j++) {
        //     if (aps[i][min(t-1,j+1)][0].any()) {
        //         apsAdjusted[i][j][0] = aps[i][j][0] | ((aps[i][j][0] & oddMask) << 1);
        //     }
        //     if (aps[i][min(t-1,j+1)][1].any()) {
        //         apsAdjusted[i][j][1] = aps[i][j][1] | ((aps[i][j][1] & evenMask) << 1);
        //     }

        //     apsAdjusted[i][j] = segmentSuffix(apsAdjusted[i][j]);
        // }

        for (int j = 0; j < min(indStart, t - 1); j++) {
            apsAdjusted[i][j] = segmentSuffix(apsAdjusted[i][j]);
        }


        for (int j = t - 1; j >= max(indEnd, 0); j--) {
            apsAdjusted[i][j] = segmentPrefix(apsAdjusted[i][j]);
        }



        // for (int j = indEnd; j < segmentation.size() - 1; j++) {
        //     if (aps[i][max(0,j)][0].any()) {
        //         bitset<SIZE> asdf = aps[i][j-1][1] & evenMask;
        //         asdf = asdf << 1;
        //         asdf = asdf | aps[i][j-1][0];
        //         apsAdjusted[i][j-1][0] = asdf;
        //     }
        //     if (aps[i][max(0,j)][1].any()) {
        //         bitset<SIZE> asdf = aps[i][j-1][0] & oddMask;
        //         asdf = asdf << 1;
        //         asdf = asdf | aps[i][j-1][1];
        //         apsAdjusted[i][j-1][0] = asdf;
        //     }

        //     apsAdjusted[i][j] = segmentPrefix(apsAdjusted[i][j]);
        // }
    }


    return apsAdjusted;
}

vector<vector<set<string>>> computeValueExpressionsCoarse(
    const vector<vector<pair<long long, double>>> &signals,
    const vector<vector<vector<long long>>> &uncertainties,
    const vector<long long> &segmentation,
    const set<int> &exactSignals = {})
{
    // Inputs retain the validated alignment established by the pipeline.
    vector<vector<set<string>>> coarse(
        signals.size(), vector<set<string>>(segmentation.size() - 1));
    for (size_t signalIndex = 0; signalIndex < signals.size(); ++signalIndex)
    {
        const bool exact = exactSignals.count(static_cast<int>(signalIndex)) != 0;
        for (size_t segment = 0; segment + 1 < segmentation.size(); ++segment)
        {
            const long long start = segmentation[segment];
            const long long end = segmentation[segment + 1];
            double minimum = numeric_limits<double>::infinity();
            double maximum = -numeric_limits<double>::infinity();
            auto include = [&](double value)
            {
                minimum = min(minimum, value);
                maximum = max(maximum, value);
            };

            if (exact)
            {
                include(signalValueAtUnchecked(signals[signalIndex], start));
                for (size_t edge = 1; edge < signals[signalIndex].size(); ++edge)
                {
                    const long long timestamp = signals[signalIndex][edge].first;
                    if (timestamp <= start)
                    {
                        continue;
                    }
                    if (timestamp >= end)
                    {
                        break;
                    }
                    if (isActualEdge(signals[signalIndex], edge))
                    {
                        include(signals[signalIndex][edge].second);
                    }
                }
            }
            else
            {
                bool intersects = false;
                for (size_t edge = 1; edge < signals[signalIndex].size(); ++edge)
                {
                    if (!isActualEdge(signals[signalIndex], edge))
                    {
                        continue;
                    }
                    const vector<long long> &region = uncertainties[signalIndex][edge];
                    if (region[0] >= end || region[1] <= start)
                    {
                        continue;
                    }
                    intersects = true;
                    include(signals[signalIndex][edge - 1].second);
                    include(signals[signalIndex][edge].second);
                }
                if (!intersects)
                {
                    include(signalValueAtUnchecked(signals[signalIndex], start));
                }
            }

            if (!isfinite(minimum) || !isfinite(maximum))
            {
                throw logic_error("coarse abstraction contains no values");
            }
            coarse[signalIndex][segment].insert(formatDouble(minimum));
            coarse[signalIndex][segment].insert(formatDouble(maximum));
        }
    }
    return coarse;
}
