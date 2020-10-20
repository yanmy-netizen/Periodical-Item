#include <bits/stdc++.h>
#include "BOBHash32.h"
#include "BOBHash64.h"

using namespace std;

#define inf 1e9
#define T 20
#define MAXL 30000

class Alg
{    
public:
    int w, h, nT, M;//nT次hash 共t*M个位置
    int l;//在一整段时间中取l个时间点
    int _K;//返回最大的_K个数
    BOBHash32* bobhash_sk[T];
    BOBHash32* bobhash_hg;
    double b, CenterL, CenterR, delta;
    float center[31];

    double sk[MAXL];//A部分，一个时间戳
    uint64_t hg_id[31][MAXL];
    int hg_count[31][MAXL];

    Alg(int nT, int M, int w, int h, double b, int l, int _K, double CenterL, double CenterR, double delta) : nT(nT), M(M), w(w), h(h), b(b), l(l), _K(_K), CenterL(CenterL), CenterR(CenterR), delta(delta) {
        srand(time(NULL));
        int x = 123;
        for (int k = 0; k < nT; ++k) bobhash_sk[k] = new BOBHash32(x + k);
        bobhash_hg = new BOBHash32(x);
        memset(sk, 0, sizeof(sk));
        memset(hg_id, 0, sizeof(hg_id));
        memset(hg_count, 0, sizeof(hg_count));
        if (l == 1) {
            center[0] = (double)(CenterR + CenterL) / 2; 
            return;
        }
        double de = (double)(CenterR - CenterL - 2 * delta) / (l-1);
        center[0] = CenterL + delta;
        for (int i = 1; i < l; ++i) {
            center[i] = center[i - 1] + de;
        }
    }

    void insert(const uint64_t& x, const double& y) {
        double last_time = y + 1; 
        for (int k = 0; k < nT; ++k) {
            unsigned long long int H = bobhash_sk[k] -> run((char *)&x, 8);
            unsigned long long int pos = H % M;
            
            last_time = min(last_time, sk[k * M + pos]);
            sk[k * M + pos] = y;
        }
        unsigned long long int H = bobhash_hg -> run((char *)&x, 8);
        unsigned long long int pos = H % w;
        for (int i = 0; i < l; i++) {
        
            if (y - last_time > center[i] + delta) continue;
            if (y - last_time < center[i] - delta) continue;

            int res = inf;
            int posi = 0;
            bool con = false;
            for (int k = 0; k < h; ++k) {
                if (!hg_id[i][pos * h + k]) {
                    hg_id[i][pos * h + k] = x;
                    hg_count[i][pos * h + k] = 1;
                    con = true;
                    break;
                } else if (hg_id[i][pos * h + k] == x) {
                    ++hg_count[i][pos * h + k];
                    con = true;
                    break;
                } else if (hg_count[i][pos * h + k] < res){
                    res = hg_count[i][pos * h + k];
                    posi = k;
                }
            }
            if (con) continue;
            double e = 1.0 / pow(b, res);
            if (double(rand()) / RAND_MAX < e) {
                --hg_count[i][pos * h + posi];
                if (!hg_count[i][pos * h + posi]) {
                    hg_id[i][pos * h + posi] = x;
                    hg_count[i][pos * h + posi] = 1;
                }
            }
        }
    }

    int query(const uint64_t& x) {
        int ret = 0;
        unsigned long long int H = bobhash_hg -> run((char *)&x, 8);
        unsigned long long int pos = H % w;

        for (int i = 0; i < l; ++i) {
            for (int k = 0; k < h; ++k) {
                if (hg_count[i][pos * h + k]) {
                    if (hg_id[i][pos * h + k] == x)
                        ret = max(ret, hg_count[i][pos * h + k]);
                } 
            }
        }
        return ret;
    }

    void getTopK(set <pair<int, uint64_t> > & s) {
        set <pair<int, uint64_t> > topK;
        map <uint64_t, int> num;
        s.clear();
        for (int i = 0; i < l; ++i) {
            int cur = 0;
            int sze = w * h;
            while (cur < sze) {
                uint64_t x = hg_id[i][cur];
                int y = hg_count[i][cur];
                if (num.find(x) != num.end()) {
                    if (num[x] < y) {
                        num[x] = y;
                    }
                } else {
                    num[x] = y;
                }
                ++cur;
            }
        }
        
        for (map<uint64_t, int>::iterator it = num.begin(); it != num.end(); ++it){
            if (s.size() < _K) {
                s.insert(make_pair(it->second, it->first));
            }
            else {
                pair<int, uint64_t> p = make_pair(it->second, it->first);
                set<pair<int, uint64_t> >::iterator iter = s.begin(); 
                if (p > *iter) {
                    s.erase(iter);
                    s.insert(p);
                }
            }
        }
    }
};
