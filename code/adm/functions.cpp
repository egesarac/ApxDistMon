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
using namespace std;

#define SIZE 1000 // TODO: this should be as tight as possible

bitset<SIZE> evenMask;
bitset<SIZE> oddMask;

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
    int len = count(S.begin()->begin(),S.begin()->end(), ';');
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
                                out[i].insert(to_string(e1 + e2));
                            }
                        }
                    }
                }
            }
        }
    }

    return out;
}

vector<set<string>> abstProdCoarseStrSum(const vector<set<string>> &v1, const vector<set<string>> &v2)
{
    vector<set<string>> out(v1.size());

    for (int i = 0; i < v1.size(); i++)
    {
        double m1 = DBL_MAX;
        double m2 = DBL_MAX;
        double M1 = DBL_MIN;
        double M2 = DBL_MIN;

        for (auto s1 : v1[i])
        {
            if (!s1.empty())
            {
                vector<double> vals1 = split(s1, ";");
                auto [min1, max1] = minmax_element(vals1.begin(), vals1.end());

                if (*min1 < m1)
                {
                    m1 = *min1;
                }

                if (*max1 > M1) 
                {
                    M1 = *max1;
                }
            }
        }

        for (auto s2 : v2[i])
        {
            if(!s2.empty())
            {
                vector<double> vals2 = split(s2, ";");
                auto [min2, max2] = minmax_element(vals2.begin(), vals2.end());

                if (*min2 < m2)
                {
                    m2 = *min2;
                }

                if (*max2 > M2) 
                {
                    M2 = *max2;
                }
            }
        }

        out[i].insert(to_string(M1 + M2));
        out[i].insert(to_string(m1 + m2));
    }

    return out;
}

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
                                out[i].insert(to_string((e1 - e2) * (e1 - e2)));
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
                                r += to_string((stod(r1) - stod(r2)) * (stod(r1) - stod(r2))) + ";";
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
                                r += to_string(stod(r1) - stod(r2)) + ";";
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
                                r += to_string(stod(r1) + stod(r2)) + ";";
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

template <std::size_t N>
vector<bitset<N>> segmentFirstBit(vector<bitset<N>> v1)
{
    v1[0] = v1[0] & bitset<N>(1);
    v1[1] = v1[1] & bitset<N>(1);

    return v1;
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
vector<bitset<N>> bitsetConcat(const vector<bitset<N>> &v1, const vector<bitset<N>> &v2)
{
    if (v1[0].none() && v1[1].none())
    {
        return v2;
    }

    vector<bitset<N>> out(2);
    for (int i = 0; i < 2; i++)
    {
        for (int j = 0; j <= msb(v1[i]); j++)
        {
            if (v1[i][j] == true)
            {
                for (int x = 0; x < 2; x++)
                {
                    for (int y = 0; y <= msb(v2[x]); y++)
                    {
                        if (v2[x][y] == true)
                        {
                            if ((i == x && j % 2 == 1) || (i != x && j % 2 == 0))
                            {
                                out[i][j + y + 1] = true;
                            }
                            else
                            {
                                out[i][j + y] = true;
                            }
                        }
                    }
                }
            }
        }
    }

    return out;
}

template <std::size_t N>
vector<vector<bitset<N>>> getProfiles(vector<vector<bitset<N>>> v, vector<long long> segmentation, const long long &t1, const long long &t2, const long long &a, const long long &b, const bool &leftClosed, const bool &rightClosed)
{
    vector<vector<bitset<N>>> profiles;

    // // to handle the cases where the window falls out of the signal
    // segmentation.push_back(INT32_MAX);

    // // this assumes the signal is false
    // // vector<bitset<N>> temp(2);
    // // temp[0][0] = true;
    // // v.push_back(temp);

    // // this assumes the interval is trimmed to fit the temporal domain (equivalently, signal is assumed to extend beyond its domain)
    // v.push_back(v[v.size() - 1]);

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
vector<vector<bitset<N>>> bitsetAlways(const vector<vector<bitset<N>>> &v1)
{
    vector<vector<bitset<N>>> vv(v1.size());

    for (auto &v : vv)
    {
        v.resize(2);
    }

    bool firstBit0 = false;
    bool firstBit1 = true;

    for (int i = v1.size() - 1; i >= 0; i--)
    {
        if (firstBit1 == true)
        {
            vector<int> len1(4);
            len1[0] = msb(v1[i][0] & evenMask) + 1; // 0...0
            len1[1] = msb(v1[i][0] & oddMask) + 1;  // 0...1
            len1[2] = msb(v1[i][1] & oddMask) + 1;  // 1...0
            len1[3] = msb(v1[i][1] & evenMask) + 1; // 1...1

            vv[i][0][0] = (max(len1[0], len1[2]) > 0);
            vv[i][0][1] = (max(len1[1], len1[3] - 2) > 0);
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
vector<vector<bitset<N>>> bitsetEventually(const vector<vector<bitset<N>>> &v1)
{
    vector<vector<bitset<N>>> vv(v1.size());

    for (auto &v : vv)
    {
        v.resize(2);
    }

    bool firstBit0 = true;
    bool firstBit1 = false;

    for (int i = v1.size() - 1; i >= 0; i--)
    {
        if (firstBit0 == true)
        {
            vector<int> len1(4);
            len1[0] = msb(v1[i][0] & evenMask) + 1; // 0...0
            len1[1] = msb(v1[i][0] & oddMask) + 1;  // 0...1
            len1[2] = msb(v1[i][1] & oddMask) + 1;  // 1...0
            len1[3] = msb(v1[i][1] & evenMask) + 1; // 1...1

            vv[i][0][0] = v1[i][0][0];
            vv[i][1][0] = (max(len1[1], len1[3]) > 0);
            vv[i][1][1] = (max(len1[0] - 2, len1[2]) > 0);
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

template <std::size_t N>
vector<vector<bitset<N>>> bitsetBoundedAlways(vector<vector<bitset<N>>> v1, vector<long long> segmentation, const long long &a, const long long &b, const bool &leftClosed, const bool &rightClosed)
{
    vector<vector<bitset<N>>> vv(v1.size());

    for (auto &v : vv)
    {
        v.resize(2);
    }

    for (int i = 0; i < v1.size(); i++)
    {
        vector<vector<bitset<N>>> temp = getProfiles(v1, segmentation, segmentation[i], segmentation[i + 1], a, b, leftClosed, rightClosed);

        for (int k = 0; k < temp.size(); k++)
        {
            vector<vector<bitset<N>>> tempPerSegment;
            tempPerSegment.push_back(temp[k]);
            tempPerSegment = bitsetAlways(tempPerSegment);
            temp[k] = tempPerSegment[0];
        }

        for (int j = 1; j < temp.size(); j++)
        {
            temp[0] = bitsetConcat(temp[0], temp[j]);
        }

        vv[i] = temp[0];
    }

    return vv;
}

template <std::size_t N>
vector<vector<bitset<N>>> bitsetBoundedEventually(vector<vector<bitset<N>>> v1, vector<long long> segmentation, const long long &a, const long long &b, const bool &leftClosed, const bool &rightClosed)
{
    vector<vector<bitset<N>>> vv(v1.size());

    for (auto &v : vv)
    {
        v.resize(2);
    }

    for (int i = 0; i < v1.size(); i++)
    {
        vector<vector<bitset<N>>> temp = getProfiles(v1, segmentation, segmentation[i], segmentation[i + 1], a, b, leftClosed, rightClosed);

        for (int k = 0; k < temp.size(); k++)
        {
            vector<vector<bitset<N>>> tempPerSegment;
            tempPerSegment.push_back(temp[k]);
            tempPerSegment = bitsetEventually(tempPerSegment);
            temp[k] = tempPerSegment[0];
        }

        for (int j = 1; j < temp.size(); j++)
        {
            temp[0] = bitsetConcat(temp[0], temp[j]);
        }

        vv[i] = temp[0];
    }

    return vv;
}

template <std::size_t N>
vector<vector<bitset<N>>> bitsetUntilStrict(vector<vector<bitset<N>>> v1, vector<vector<bitset<N>>> v2)
{
    vector<vector<bitset<N>>> vv(v1.size());

    for (auto &v : vv)
    {
        v.resize(2);
    }

    bool firstBit0 = true;
    bool firstBit1 = false;

    for (int i = v1.size() - 1; i >= 0; i--)
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


template <std::size_t N>
vector<vector<bitset<N>>> bitsetUntilNonStrict(vector<vector<bitset<N>>> v1, vector<vector<bitset<N>>> v2)
{
    vector<vector<bitset<N>>> vv(v1.size());

    for (auto &v : vv)
    {
        v.resize(2);
    }

    bool firstBit0 = true;
    bool firstBit1 = false;

    for (int i = v1.size() - 1; i >= 0; i--)
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

template <std::size_t N>
vector<vector<bitset<N>>> bitsetConjunction(vector<vector<bitset<N>>> v1, vector<vector<bitset<N>>> v2)
{
    vector<vector<bitset<N>>> vv(v1.size());

    for (auto &v : vv)
    {
        v.resize(2);
    }

    for (int i = 0; i < v1.size(); i++)
    {
        // handling the corner cases of conjunction with 0 or 1
        if (v1[i][0][0] == true)
        {
            vv[i][0][0] = true;
            v1[i][0][0] = false;
        }
        if (v2[i][0][0] == true)
        {
            vv[i][0][0] = true;
            v2[i][0][0] = false;
        }
        if (v1[i][1][0] == true)
        {
            vv[i][0] = vv[i][0] | v2[i][0];
            vv[i][1] = vv[i][1] | v2[i][1];
            v1[i][1][0] = false;
        }
        if (v2[i][1][0] == true)
        {
            vv[i][0] = vv[i][0] | v1[i][0];
            vv[i][1] = vv[i][1] | v1[i][1];
            v2[i][1][0] = false;
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

        for (int i = 0; i < 4; i++)
        {
            if (len1[i] == 0)
            {
                len1[i] = INT32_MIN / 2 + 2;
            }

            if (len2[i] == 0)
            {
                len2[i] = INT32_MIN / 2 + 2;
            }
        }

        vector<int> len(4);
        len[0] = max({len1[0] + len2[0] - 2, len1[0] + len2[1] - 1, len1[0] + len2[2] - 1, len1[0] + len2[3], len1[1] + len2[0] - 1, len1[1] + len2[2], len1[2] + len2[0] - 1, len1[2] + len2[1], len1[3] + len2[0]}) - 1;
        len[1] = max({len1[1] + len2[1] - 1, len1[1] + len2[3], len1[3] + len2[1]}) - 1;
        len[2] = max({len1[2] + len2[2] - 1, len1[2] + len2[3], len1[3] + len2[2]}) - 1;
        len[3] = len1[3] + len2[3] - 1;

        int ctr;
        for (int j = 0; j < 3; j++)
        {
            ctr = len[j] - 1;

            while (ctr >= 0)
            {
                vv[i][j / 2][ctr] = true;
                ctr = ctr - 2;
            }
        }
        ctr = len[3] - 1;
        while (ctr > 0)
        {
            vv[i][1][ctr] = true;
            ctr = ctr - 2;
        }
    }

    return vv;
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

pair<int, int> convertIntoBoolWithDestutter(string str, string delimiter, float thr)
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
            if (stof(str.substr(start, length)) > thr)
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

        if (stof(str.substr(start, str.length())) > thr)
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

vector<vector<pair<long long, double>>> convertSignalsToBool(const vector<vector<pair<long long, double>>> &signalsReal)
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

            for (int j = 1; j < signalsReal[i].size(); j++)
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

vector<vector<vector<long long>>> computeUncertaintyIntervals(const vector<vector<pair<long long, double>>> &signals, const int &eps, const int &del)
{
    vector<vector<vector<long long>>> uncertainties(signals.size());

    for (int i = 0; i < signals.size(); i++)
        {
            int gi = signals[i].size();
            uncertainties[i].resize(gi);

            uncertainties[i][0].push_back(-del);
            for (int j = 1; j < gi; j++)
            {
                // if (signals[i][j].second != signals[i][j - 1].second)
                // {
                    long long tj = signals[i][j].first;
                    long long low = max(uncertainties[i][j - 1][0] + del, tj - eps);
                    uncertainties[i][j].push_back(low);                   
                // }

            }

            uncertainties[i][gi - 1].push_back(signals[i][gi - 1].first + eps);
            for (int j = gi - 2; j > 0; j--)
            {
                // if (signals[i][j].second != signals[i][j + 1].second)
                // {
                    long long tj = signals[i][j].first;
                    long long high = min(uncertainties[i][j + 1][1] - del, tj + eps);
                    uncertainties[i][j].push_back(high);
                // }
            }
        }
    
    return uncertainties;
}

vector<long long> computeCanonicalSegmentation(const vector<vector<pair<long long, double>>> &signals, const vector<vector<vector<long long>>> &uncertainties, const int &d)
{
    set<long long> segmentation_temp;

    segmentation_temp.insert(0);
    segmentation_temp.insert(max(signals[0][signals[0].size() - 1].first, (long long)(d)));

    for (const auto &s : uncertainties)
    {
        for (int j = 1; j < s.size(); j++)
        {
            for (auto b : s[j])
            {
                segmentation_temp.insert(b);
            }
        }
    }

    vector<long long> segmentation;
    for (const auto &b : segmentation_temp)
    {
        segmentation.push_back(b);
    }

    return segmentation;
}

vector<vector<set<string>>> computeRelativeValueExpressions(const set<int> &leaders, const vector<vector<pair<long long, double>>> &signals, const vector<vector<vector<long long>>> &uncertainties, const vector<long long> &segmentation)
{
    vector<vector<set<string>>> valExprs(signals.size());
    int numSegments = segmentation.size() - 1;

    for (auto leader : leaders) {
        valExprs[leader].resize(numSegments);

        for (int j = 0; j < numSegments; j++) 
        {
            set<string> temp;
            temp.insert(""); // ??

            string expr = "";

            for (auto p : signals[leader])
            {
                if (segmentation[j] <= p.first) {
                    if (p.first < segmentation[j + 1]) {
                        expr += to_string(p.second) + ";";
                    }
                    else {
                        break;
                    }
                }
            }

            if (expr == "") {
                expr = to_string(signals[leader][signals[leader].size() - 1].second) + ";";
            }

            temp.insert(destutterStr(expr.substr(0, expr.length() - 1)));
            valExprs[leader][j] = temp;
        }
    }


    for (int i = 0; i < signals.size(); i++)
        {
            if (leaders.find(i) == leaders.end()) {
                valExprs[i].resize(numSegments);

                for (int j = 0; j < numSegments; j++)
                {
                    //valExprs[i][j].insert("");

                    for (int k = 1; k < uncertainties[i].size(); k++)
                    {
                        if (segmentation[j] == uncertainties[i][k][0] && segmentation[j + 1] == uncertainties[i][k][1])
                        {
                            // entire expression
                            set<string> temp;
                            string expr;
                            if (signals[i][k - 1].second == signals[i][k].second) {
                                expr = to_string(signals[i][k - 1].second);
                            }
                            else
                            {
                                expr = to_string(signals[i][k - 1].second) + ";" + to_string(signals[i][k].second);
                            }
                            if (!(k < uncertainties[i].size() - 1 && uncertainties[i][k][0] == uncertainties[i][k + 1][0]))
                                temp.insert("");
                            temp.insert(expr);
                            valExprs[i][j] = concatWithDestutter(valExprs[i][j], temp);
                        }
                        else if (segmentation[j] == uncertainties[i][k][0] && segmentation[j + 1] < uncertainties[i][k][1])
                        {
                            // prefix
                            set<string> temp;
                            string expr;
                            if (signals[i][k - 1].second == signals[i][k].second) {
                                expr = to_string(signals[i][k - 1].second);
                            }
                            else
                            {
                                expr = to_string(signals[i][k - 1].second) + ";" + to_string(signals[i][k].second);
                            }
                            if (!(k < uncertainties[i].size() - 1 && uncertainties[i][k][0] == uncertainties[i][k + 1][0]))
                                temp.insert("");
                            temp.insert(expr);
                            temp.insert(expr.substr(0, expr.find(';')));
                            valExprs[i][j] = concatWithDestutter(valExprs[i][j], temp);
                        }
                        else if (segmentation[j] > uncertainties[i][k][0] && segmentation[j + 1] == uncertainties[i][k][1])
                        {
                            // suffix
                            set<string> temp;
                            string expr;
                            if (signals[i][k - 1].second == signals[i][k].second) {
                                expr = to_string(signals[i][k - 1].second);
                            }
                            else
                            {
                                expr = to_string(signals[i][k - 1].second) + ";" + to_string(signals[i][k].second);
                            }
                            if (!(k < uncertainties[i].size() - 1 && uncertainties[i][k][0] == uncertainties[i][k + 1][0]))
                                temp.insert("");
                            temp.insert(expr);
                            temp.insert(expr.substr(expr.find(';') + 1, expr.length()));
                            valExprs[i][j] = concatWithDestutter(valExprs[i][j], temp);
                        }
                        else if (segmentation[j] > uncertainties[i][k][0] && segmentation[j + 1] < uncertainties[i][k][1])
                        {
                            // infix
                            set<string> temp;
                            string expr;
                            if (signals[i][k - 1].second == signals[i][k].second) {
                                expr = to_string(signals[i][k - 1].second);
                            }
                            else
                            {
                                expr = to_string(signals[i][k - 1].second) + ";" + to_string(signals[i][k].second);
                            }
                            if (!(k < uncertainties[i].size() - 1 && uncertainties[i][k][0] == uncertainties[i][k + 1][0]))
                                temp.insert("");
                            temp.insert(expr);
                            temp.insert(expr.substr(0, expr.find(';')));
                            temp.insert(expr.substr(expr.find(';') + 1, expr.length()));
                            valExprs[i][j] = concatWithDestutter(valExprs[i][j], temp);
                        }
                    }

                    if (valExprs[i][j].empty())
                    {
                        int kk = 1;
                        while (kk < uncertainties[i].size() && uncertainties[i][kk][1] <= segmentation[j])
                        {
                            kk++;
                        }
                        valExprs[i][j].insert(to_string(signals[i][kk - 1].second));
                    }
                }

            }
        }

    return valExprs;
}


// TODO: check and improve
vector<vector<set<string>>> computeValueExpressions(const vector<vector<pair<long long, double>>> &signals, const vector<vector<vector<long long>>> &uncertainties, const vector<long long> &segmentation)
{
    vector<vector<set<string>>> valExprs(signals.size());
    int numSegments = segmentation.size() - 1;

    for (int i = 0; i < signals.size(); i++)
        {
            valExprs[i].resize(numSegments);

            for (int j = 0; j < numSegments; j++)
            {
                //valExprs[i][j].insert("");

                for (int k = 1; k < uncertainties[i].size(); k++)
                {
                    if (segmentation[j] == uncertainties[i][k][0] && segmentation[j + 1] == uncertainties[i][k][1])
                    {
                        // entire expression
                        set<string> temp;
                        string expr;
                        if (signals[i][k - 1].second == signals[i][k].second) {
                            expr = to_string(signals[i][k - 1].second);
                        }
                        else
                        {
                            expr = to_string(signals[i][k - 1].second) + ";" + to_string(signals[i][k].second);
                        }
                        if (!(k < uncertainties[i].size() - 1 && uncertainties[i][k][0] == uncertainties[i][k + 1][0]))
                            temp.insert("");
                        temp.insert(expr);
                        valExprs[i][j] = concatWithDestutter(valExprs[i][j], temp);
                    }
                    else if (segmentation[j] == uncertainties[i][k][0] && segmentation[j + 1] < uncertainties[i][k][1])
                    {
                        // prefix
                        set<string> temp;
                        string expr;
                        if (signals[i][k - 1].second == signals[i][k].second) {
                            expr = to_string(signals[i][k - 1].second);
                        }
                        else
                        {
                            expr = to_string(signals[i][k - 1].second) + ";" + to_string(signals[i][k].second);
                        }
                        if (!(k < uncertainties[i].size() - 1 && uncertainties[i][k][0] == uncertainties[i][k + 1][0]))
                            temp.insert("");
                        temp.insert(expr);
                        temp.insert(expr.substr(0, expr.find(';')));
                        valExprs[i][j] = concatWithDestutter(valExprs[i][j], temp);
                    }
                    else if (segmentation[j] > uncertainties[i][k][0] && segmentation[j + 1] == uncertainties[i][k][1])
                    {
                        // suffix
                        set<string> temp;
                        string expr;
                        if (signals[i][k - 1].second == signals[i][k].second) {
                            expr = to_string(signals[i][k - 1].second);
                        }
                        else
                        {
                            expr = to_string(signals[i][k - 1].second) + ";" + to_string(signals[i][k].second);
                        }
                        if (!(k < uncertainties[i].size() - 1 && uncertainties[i][k][0] == uncertainties[i][k + 1][0]))
                            temp.insert("");
                        temp.insert(expr);
                        temp.insert(expr.substr(expr.find(';') + 1, expr.length()));
                        valExprs[i][j] = concatWithDestutter(valExprs[i][j], temp);
                    }
                    else if (segmentation[j] > uncertainties[i][k][0] && segmentation[j + 1] < uncertainties[i][k][1])
                    {
                        // infix
                        set<string> temp;
                        string expr;
                        if (signals[i][k - 1].second == signals[i][k].second) {
                            expr = to_string(signals[i][k - 1].second);
                        }
                        else
                        {
                            expr = to_string(signals[i][k - 1].second) + ";" + to_string(signals[i][k].second);
                        }
                        if (!(k < uncertainties[i].size() - 1 && uncertainties[i][k][0] == uncertainties[i][k + 1][0]))
                            temp.insert("");
                        temp.insert(expr);
                        temp.insert(expr.substr(0, expr.find(';')));
                        temp.insert(expr.substr(expr.find(';') + 1, expr.length()));
                        valExprs[i][j] = concatWithDestutter(valExprs[i][j], temp);
                    }
                }

                if (valExprs[i][j].empty())
                {
                    int kk = 1;
                    while (kk < uncertainties[i].size() && uncertainties[i][kk][1] <= segmentation[j])
                    {
                        kk++;
                    }
                    valExprs[i][j].insert(to_string(signals[i][kk - 1].second));
                }
            }
        }

    return valExprs;
}

bitset<SIZE> generateBitmask(const int &parity)
{
    bitset<SIZE> mask;

    for (int i = parity % 2; i < SIZE; i = i + 2)
    {
        mask[i] = true;
    }

    return mask;
}

vector<vector<vector<bitset<SIZE>>>> convertSignalsToAtomicPropositions(const vector<vector<set<string>>> &valExprs, const float &thr)
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

vector<vector<set<string>>> computeValueExpressionsCoarse(const vector<vector<pair<long long, double>>> &signals, const vector<vector<vector<long long>>> &uncertainties, const vector<long long> &segmentation)
{
    vector<vector<set<string>>> valExprs(signals.size());
    int numSegments = segmentation.size() - 1;

    for (int i = 0; i < signals.size(); i++)
        {
            valExprs[i].resize(numSegments);

            for (int j = 0; j < numSegments; j++)
            {
                //valExprs[i][j].insert("");
                double m = DBL_MAX;
                double M = -DBL_MAX;

                for (int k = 1; k < uncertainties[i].size(); k++)
                {
                    if (segmentation[j] < uncertainties[i][k][1] && segmentation[j + 1] > uncertainties[i][k][0])
                    {
                        if (m > min(signals[i][k - 1].second, signals[i][k].second))
                        {
                            m = min(signals[i][k - 1].second, signals[i][k].second);
                        }

                        if (M < max(signals[i][k - 1].second, signals[i][k].second))
                        {
                            M = max(signals[i][k - 1].second, signals[i][k].second);
                        }
                    }
                }

                if (m == DBL_MAX) 
                {
                    m = signals[i][j].second;
                    M = m;
                }

                valExprs[i][j].insert(to_string(m));
                valExprs[i][j].insert(to_string(M));
            }
        }

    return valExprs;
}