#include <bits/stdc++.h>
#include "detector_level2_105.h"

using namespace std;

#define inf 1e9
#define T 20
#define CenterL 0.04
#define CenterR 0.2
#define delta 0.005
#define rangeL CenterL-delta
#define rangeR CenterR+delta

string datapath[60] = {"CAIDA2018/dataset/130000.dat",
                "CAIDA2018/dataset/130100.dat",
                "CAIDA2018/dataset/130200.dat",
                "CAIDA2018/dataset/130300.dat",
                "CAIDA2018/dataset/130400.dat",
                "CAIDA2018/dataset/130500.dat",
                "CAIDA2018/dataset/130600.dat",
                "CAIDA2018/dataset/130700.dat",
                "CAIDA2018/dataset/130800.dat",
                "CAIDA2018/dataset/130900.dat",
                "CAIDA2018/dataset/131000.dat",
                "CAIDA2018/dataset/131100.dat",
                "CAIDA2018/dataset/131200.dat",
                "CAIDA2018/dataset/131300.dat",
                "CAIDA2018/dataset/131400.dat",
                "CAIDA2018/dataset/131500.dat",
                "CAIDA2018/dataset/131600.dat",
                "CAIDA2018/dataset/131700.dat",
                "CAIDA2018/dataset/131800.dat",
                "CAIDA2018/dataset/131900.dat",
                "CAIDA2018/dataset/132000.dat",
                "CAIDA2018/dataset/132100.dat",
                "CAIDA2018/dataset/132200.dat",
                "CAIDA2018/dataset/132300.dat",
                "CAIDA2018/dataset/132400.dat",
                "CAIDA2018/dataset/132500.dat",
                "CAIDA2018/dataset/132600.dat",
                "CAIDA2018/dataset/132700.dat",
                "CAIDA2018/dataset/132800.dat",
                "CAIDA2018/dataset/132900.dat",
                "CAIDA2018/dataset/133000.dat",
                "CAIDA2018/dataset/133100.dat",
                "CAIDA2018/dataset/133200.dat",
                "CAIDA2018/dataset/133300.dat",
                "CAIDA2018/dataset/133400.dat",
                "CAIDA2018/dataset/133500.dat",
                "CAIDA2018/dataset/133600.dat",
                "CAIDA2018/dataset/133700.dat",
                "CAIDA2018/dataset/133800.dat",
                "CAIDA2018/dataset/133900.dat",
                "CAIDA2018/dataset/134000.dat",
                "CAIDA2018/dataset/134100.dat",
                "CAIDA2018/dataset/134200.dat",
                "CAIDA2018/dataset/134300.dat",
                "CAIDA2018/dataset/134400.dat",
                "CAIDA2018/dataset/134500.dat",
                "CAIDA2018/dataset/134600.dat",
                "CAIDA2018/dataset/134700.dat",
                "CAIDA2018/dataset/134800.dat",
                "CAIDA2018/dataset/134900.dat",
                "CAIDA2018/dataset/135000.dat",
                "CAIDA2018/dataset/135100.dat",
                "CAIDA2018/dataset/135200.dat",
                "CAIDA2018/dataset/135300.dat",
                "CAIDA2018/dataset/135400.dat",
                "CAIDA2018/dataset/135500.dat",
                "CAIDA2018/dataset/135600.dat", 
                "CAIDA2018/dataset/135700.dat",
                "CAIDA2018/dataset/135800.dat",
                "CAIDA2018/dataset/135900.dat"};

ifstream fin;

const int M = 2e7;
const int _K = 5000;

pair <uint64_t, double> Read()//新CAIDA
{   
    static bool isfirstRead = true;
    static int curFinID = 0;
    static double offset = 0;
    static double lastT = 0;

	double t; uint64_t s, _s;
    if (isfirstRead) {
        isfirstRead = false;
        fin.open(datapath[curFinID], std :: ios :: binary);
    }

	if (fin.eof()) {
        fin.close();
        fin.open(datapath[++curFinID], std :: ios :: binary);
        if (curFinID > 60) {
            fin.close();
            exit(0);
        }
    }
    fin.read((char*)&s, sizeof(uint64_t)); //srcip(4)+dstip(4) 
    fin.read((char*)&_s, 5);//srcport destport protcol
    fin.read((char*)&t, sizeof(double));

    
    t += offset;
    if(t < lastT) {
        offset += lastT - t;
        t += (lastT - t);
    }
    lastT = t;

	return make_pair(s, t * 2); 
}

pair <uint64_t, double> input[M + 7];

struct getGT{
    map <uint64_t, set <pair<double, bool> > > table;
    
    void init() {
        table.clear();
    }

    void insert(uint64_t id, double key) {
        if (key < rangeL + delta) return;
        if (key > rangeR - delta) return;
        table[id].insert(make_pair(key - delta, true));
        table[id].insert(make_pair(key + delta, false));
    }

    int query(const uint64_t& id) {
        if (table.find(id) == table.end()) return 0;
        int ret = 0;
        int cur = 0;
        for(auto it = table.find(id)->second.begin(); it != table.find(id)->second.end(); ++it){
            if (it->second) {
                ++cur;
                if (cur > ret) {
                    ret = cur;
                }
            } else {
                --cur;
            }
        }
        return ret;
    }
}intervalGT;

map <uint64_t, double> timeStamp;
map <uint64_t, int> intervalAnswer;

int total_ele, total_save;
pair <int, uint64_t> ele[M + 7], save[M + 7];
set <uint64_t> topKid, allid;
set<pair<int, uint64_t> >topKid_ours;

void GroundTruth() {
    for (int i = 0; i < M; ++i) {
        auto e = input[i]; 

        uint64_t id = e.first;
        double curTime = e.second;
        double lastTime = -1;
        
        if (timeStamp.find(id) != timeStamp.end()) {
            lastTime = timeStamp[id];
            intervalGT.insert(id, curTime - lastTime);
        }
        
        timeStamp[id] = curTime;
    }

    intervalAnswer.clear();
    total_ele = 0;
    for (int i = 0; i < M; ++i){
        auto e = input[i]; 
        uint64_t id = e.first;
        
        if (intervalAnswer.find(id) != intervalAnswer.end()) continue;

        allid.insert(id);

        int result = intervalGT.query(id);
        intervalAnswer[id] = result;
        ele[++total_ele] = make_pair(result, id);
    }

    sort(ele + 1, ele + total_ele + 1);
    for (int k = total_ele, i = 1; i <= _K; ++i, --k) {
        //cerr << '#' << '#' << "rank" << i <<' ' << ele[k].first << endl;
        topKid.insert(ele[k].second);
    }

}

int cnt;

double areK[2], aaeK[2];
double are[2], aae[2];

void OURS(int _L, int nK, int per_mem) {
    int memory = (_L + 1) * per_mem;
    Alg* Detector = new Alg(3, memory * 1024/((_L+1)*8*3), memory * 1024/((_L+1)*12*4), 4, 1.08, _L, _K, CenterL, CenterR, delta); 
    int topk_size = _K;

    for (int i = 0; i < M; ++i) {
        auto e = input[i]; 
        uint64_t id = e.first;
        double curTime = e.second;
        Detector -> insert(id, curTime);
    }

    Detector -> getTopK(topKid_ours);

    int both = 0, sze = 0;
    for (auto it = topKid_ours.begin(); it != topKid_ours.end(); ++it, ++sze) 
        if (topKid.find(it->second) != topKid.end()) ++both;
    
    printf("%d,%.6lf,%d,%d",memory, (double)both / _K, topk_size, _L);

    delete Detector;
}


int main() {
    freopen("level2_10599.csv", "w", stdout);
    srand(0);
    for (int i = 0; i < M + 1; ++i) input[i] = Read();
    
    GroundTruth();
    printf("CenterL:%.3lf,CenterR:%.3lf,delta:%.3lf\n",CenterL, CenterR, delta);
    puts("memory/KB,precision/recall,size,_L");

    for (int j = 100; j <= 500; j += 100) 
        for (int i = 1; i <= 16; i++)
            OURS(i, _K, j), putchar('\n');
	return 0;
}