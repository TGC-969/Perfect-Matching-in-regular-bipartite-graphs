#include<bits/stdc++.h>
using namespace std;
#define NMAX 100010
#define endl "\n"

/*
    n = |A| = |B|
    d is degree of each vertex

    s = 0, t = 2 * n + 1

    m = nd

    vector<int> p; stores the walk temporarily


*/



class RandomNumberGenerator
{
    // static mt19937 gen(std::random_device{}());
    static mt19937 gen;

    public:

    static int generateRand(int n){
        uniform_int_distribution<> distrib(0, n - 1);
        return distrib(gen);
    }
    static vector<int> generateRandPerm(int n){
        vector<int> p;
        for (int i = 1; i <= n; i++) p.push_back(i);
        shuffle(p.begin(), p.end(), gen);
        return p;
    }  
};
mt19937 RandomNumberGenerator::gen(std::random_device{}());

class Stats
{
    vector<int> values;
    public:
    void addVal(int val)
    {
        values.push_back(val);
    }
    long long getTotal()
    {
        long long res=0;
        for(auto val:values)res+=val;
        return res;
    }
    double getMean()
    {
        double res=getTotal();
        return res/values.size();
    }
    double getVariance()
    {
        double x=0;
        double xSquare=0;
        for(auto val:values)
        {
            x+=val;
            xSquare += (double)val * val;
        }
        x/=values.size();
        xSquare/=values.size();
        return xSquare - x * x;
    }
    double getSD()
    {
        return sqrt(getVariance());
    }
    int getMax()
    {
        int res=values[0];
        for(auto val:values)res= max(res,val);
        return res;
    }
    int getMin()
    {
        int res=values[0];
        for(auto val:values)res= min(res,val);
        return res;
    }
    void reset()
    {
        values.clear();
    }
    void printAll()
    {
        for(auto val : values){
            cout<<val<<" ";
        }
        cout<<endl;
    }
    void printfreqs()
    {
        map<int, int> freq;
        for (int v : values) freq[v]++;
        for (int i = 1; i <= 30; i++) cout << i << ' ';
        cout << endl;
        for (int i = 1; i <= 30; i++) cout << freq[i] << ' ';
        cout << endl;
    }
};


class RandomizedPerfectMatching
{   
    int n, d;
    // vector<int> adj[NMAX];
    vector<vector<int>> adj;

    int s, t;
    int k;
    vector<int> p;
    // int M[NMAX];
    vector<int> M;



    int truncatedWalkCnt = 0;
    int totalCnt=0;
    vector<int> came;
    vector<int> remainingVertices;
    void updateRemainingVertices(){
        if(remainingVertices.size() >= (k<<1))
        {
            vector<int> freeVertices;
            for(auto v:remainingVertices)if(!M[v])freeVertices.push_back(v);
            swap(freeVertices,remainingVertices);
        }
        // int sz = remainingVertices.size();
        // int idx = RandomNumberGenerator::generateRand(sz);
        // while(M[remainingVertices[idx]])
        // {
        //     idx = RandomNumberGenerator::generateRand(sz);
        // }
        // return remainingVertices[idx];
    }
    // int getUnmatched(int x){
    //     int c = 0;
    //     for (int i = 1; i <= n; i++){
    //         if (M[i] == 0) c++;
    //         if (c == x) return i;
    //     }
    //     return -1;
    // }


    int SAMPLE_OUT_EDGE(int u){
        // return t;
        int v = -1;
        int samplingCnt=0;
        if (u == s){
            // int idx = 1 + RandomNumberGenerator::generateRand(k);
            // v = getUnmatched(idx);
            // v = getUnmatched();
            updateRemainingVertices();
            int sz = remainingVertices.size();
            int idx = RandomNumberGenerator::generateRand(sz);

            int cnt=1;

            while(M[remainingVertices[idx]])
            {
                // cout << idx << endl;
                idx = RandomNumberGenerator::generateRand(sz);
                cnt++;
            }
            unmatchedSampling.addVal(cnt);
            samplingCnt += cnt;
            v = remainingVertices[idx];
        }
        else if (1 <= u && u <= n){
            v = M[u];
            int cnt=1;
            while (M[u] == v){
                int idx = RandomNumberGenerator::generateRand(d);
                // cout << "idx = " << idx << endl;
                v = adj[u][idx];
                cnt++;
            }
            nextVertexSampling.addVal(cnt);
            samplingCnt += cnt;
            if (M[v]) v = M[v];
        }
        else{
            v = t;
        }
        totalSampling.addVal(samplingCnt);
        totalCnt += samplingCnt;
        // cout << "G = " << u << ' ' << v << endl;
        return v;
    }


    bool TRUNCATED_WALK(int u, int b){
        truncatedWalkCnt++;
        if (u == t) return true;
        int v = SAMPLE_OUT_EDGE(u);
        // comment below later
        // if (v == -1) return true;
        b--;
        if (b < 0) return false;
        p.push_back(v);
        return TRUNCATED_WALK(v, b);
    }

    void removeloops(vector<int>& p){
        vector<int> up;
        for (int v : p){
            if (came[v]){
                while (!up.empty() && up.back() != v) up.pop_back();
            }
            else{
                came[v] = 1;
                up.push_back(v);
            }
        }
        p.clear();
        for (int v : up) p.push_back(v);
        for (int v : p) came[v]=0;
    }

    public:

    Stats unmatchedSampling;
    Stats nextVertexSampling;
    Stats totalSampling;
    Stats truncatedWalk;
    Stats totalTime;
    Stats totalTimeConstant;
    
    void runAlgorithm(){
        // step 1
        resetStats();
        totalCnt = 0;
        int j = 0;
        for (int i = s; i <= t; i++) M[i] = 0;
        k = n;
        for (j = 0; j < n; j++){
            // step 2
            truncatedWalkCnt=0;
            while (true){
                p.clear();
                int bj = 2 * (2 + (n / (n - j)));
                if (TRUNCATED_WALK(s, bj)) break;
            }

            truncatedWalk.addVal(truncatedWalkCnt);
            totalCnt += truncatedWalkCnt;
            // step 3
            // cout << "P = ";
            // for (int v : p) cout << v << ' ';
            // cout << endl;
            totalCnt += p.size();
            removeloops(p);
            for (int i = 0; i < p.size() - 2; i++){
                int mp = (p[i + 1] > n ? p[i + 1] : M[p[i + 1]]);
                M[p[i]] = mp;
                M[mp] = p[i];
            }
            // cout << "M = ";
            // for (int i = 1; i <= 2 * n; i++) cout << M[i] << ' ';
            // cout << endl;
            k--; // one more vertex gets matched
            //step 4
        }
        totalTime.addVal(totalCnt);
        totalTimeConstant.addVal(1 + (double)totalCnt / n / log2(n));
    }

    void runNonTruncatedAlgorithm()
    {
        // resetStats();
        initGraph(n,d,adj);
        // step 1
        totalCnt = 0;
        int j = 0;
        for (int i = s; i <= t; i++) M[i] = 0;
        k = n;
        for (j = 0; j < n; j++){
            // step 2
            truncatedWalkCnt=0;
            p.clear();
            int u = s;
            while( u != t)
            {
                u = SAMPLE_OUT_EDGE(u);
                p.push_back(u);
                truncatedWalkCnt++;
            }

            truncatedWalk.addVal(truncatedWalkCnt);
            totalCnt += truncatedWalkCnt;
            // step 3
            // cout << "P = ";
            // for (int v : p) cout << v << ' ';
            // cout << endl;
            totalCnt += p.size();
            removeloops(p);
            for (int i = 0; i < p.size() - 2; i++){
                int mp = (p[i + 1] > n ? p[i + 1] : M[p[i + 1]]);
                M[p[i]] = mp;
                M[mp] = p[i];
            }
            // cout << "M = ";
            // for (int i = 1; i <= 2 * n; i++) cout << M[i] << ' ';
            // cout << endl;
            k--; // one more vertex gets matched
            //step 4
        }
        totalTime.addVal(totalCnt);
    }

    /*
        n, d
        give adj of 1
        give adj of 2
        give adj of ..
        give adj of n

        all vertices in adj of i must be in n + 1 <= v <= 2 * n and distinct

    */
    // void getGraph(){
    //     cin >> n >> d;
    //     s = 0, t = 2 * n + 1;
    //     M.assign(2*n+2,0);
    //     adj.assign(n+1,vector<int>());
    //     for (int u = 1; u <= n; u++){
    //         // adj[u].clear();
    //         for (int j = 0; j < d; j++){
    //             int v; cin >> v;
    //             adj[u].push_back(v);
    //         }
    //     }
    // }

    void initGraph(int n,int d, vector<vector<int>> adj)
    {
        this->n = n;
        this->d = d;
        this->adj = adj;
        M.assign(2*n+2,0);
        came.assign(2*n+2,0);
        s = 0, t = 2 * n + 1;
        remainingVertices.clear();
        for(int i=1;i<=n;i++)remainingVertices.push_back(i);

        resetStats();
    }

    void resetStats()
    {
        truncatedWalk.reset();
        totalSampling.reset();
        unmatchedSampling.reset();
        nextVertexSampling.reset();
    }

    bool isPerfect(){
        set<int> V;
        for (int u = 1; u <= n; u++){
            V.insert(u);
            if (M[u] == s || M[u] == t) return false;
            V.insert(M[u]);
        }
        return V.size() == 2 * n;
    }

    vector<int> getMatching()
    {
        return M;
    }
    void printMatching()
    {
        for (int u = 1; u <= n; u++){
            cout << u << ' ' << M[u] << endl;
        }
    }
    void printRunStats()
    {
        cout << "---------------------------------------------------------------------------------" << endl;
        cout << "n: " << n << " d: " << d << endl;
        cout << "Total Run Time: " << totalCnt << endl;
        cout << "Constant in O(n log(n)): " << (double)totalCnt / double(n) / log2(n) << endl;
        cout << "Sampling Time Mean: " << totalSampling.getMean() << " Max: " << totalSampling.getMax() << " Sigma(SD): " << totalSampling.getSD() << endl;
        cout << "Truncated Walk Time Mean: " << truncatedWalk.getMean() << " Max: " << truncatedWalk.getMax() << " Constant in O(N): " << (double)truncatedWalk.getMean() / n << " Sigma(SD): " << truncatedWalk.getSD() << endl;
        cout << "---------------------------------------------------------------------------------" << endl;
    }

    void printAllRunStats()
    {
        cout << "---------------------------------------------------------------------------------" << endl;
        cout << "Printing All Runs Stats" << endl;
        cout << "Mean Run Time: " << totalTime.getMean() << endl;
        cout << "Signma(SD) Run Time: " << totalTime.getSD() << endl;
        cout << "Constant in O(n log(n)) Mean: " << totalTimeConstant.getMean() << endl;
        cout << "Constant Min: " << totalTimeConstant.getMin() << endl;
        cout << "Constant Max: " << totalTimeConstant.getMax() << endl;
        cout << "Constant Sigma: " << totalTimeConstant.getSD() << endl;
        cout << "---------------------------------------------------------------------------------" << endl;
    }
};

class InputGenerator
{

    /*
        n, d
        give adj of 1
        give adj of 2
        give adj of ..
        give adj of n

        all vertices in adj of i must be in n + 1 <= v <= 2 * n and distinct

    */
    public:
    int n,d;
    vector<vector<int>> adj;

    void getGraph(){
        cin >> n >> d;
        // s = 0, t = 2 * n + 1;
        // M.assign(2*n+2,0);
        adj.assign(n+1,vector<int>());
        for (int u = 1; u <= n; u++){
            // adj[u].clear();
            for (int j = 0; j < d; j++){
                int v; cin >> v;
                adj[u].push_back(v);
            }
        }
    }

    void getRandomInput(int n,int d)
    {   
        this->n = n;
        this->d = d;
        adj.assign(n+1, vector<int>());
        for (int i = 0; i < d; i++){
            vector<int> p = RandomNumberGenerator::generateRandPerm(n);
            // for (int j : p) cout << j << ' ';
            // cout << endl;
            for (int u = 1; u <= n; u++) adj[u].push_back(p[u - 1] + n);
        }
    }
};


int32_t main()
{
    RandomizedPerfectMatching randomizedPerfectMatching;
    InputGenerator inputGen;
    // inputGen.getGraph();
    // inputGen.getRandomInput(1000, 999);
    // // randomizedPerfectMatching.getGraph();
    // randomizedPerfectMatching.initGraph(inputGen.n,inputGen.d,inputGen.adj);
    // randomizedPerfectMatching.runAlgorithm();
    // // auto M = randomizedPerfectMatching.getMatching();

    // randomizedPerfectMatching.printMatching();
    // // for (int u = 1; u <= n; u++){
    // //     cout << u << ' ' << M[u] << endl;
    // // }

    // if (randomizedPerfectMatching.isPerfect()){
    //     cout << "M is perfect" << endl;
    // }
    // else cout << "M is not perfect" << endl;

    int n = 1000;
    int dCnt = 100;
    
    cout << "n = " << 1000 << endl;
    for (int d = 1; d <= n; d++){
        inputGen.getRandomInput(n, d);
        randomizedPerfectMatching.initGraph(inputGen.n,inputGen.d,inputGen.adj);
        
        randomizedPerfectMatching.runAlgorithm();
        cout << d << ' ';
    }
    cout << endl;
    randomizedPerfectMatching.totalTime.printAll();

    // cin >> n >> dCnt;
    // vector<int> l1, l2;
    // for(int i=1 ; i<1000 ; i++){
    //     int n=i+1;
    //     int d=max(i/2,1); 
    //     inputGen.getRandomInput(n, d);
    //     randomizedPerfectMatching.initGraph(inputGen.n,inputGen.d,inputGen.adj);

    //     randomizedPerfectMatching.runAlgorithm();
    //     l1.push_back(randomizedPerfectMatching.truncatedWalk.getMax());
    //     randomizedPerfectMatching.runNonTruncatedAlgorithm();
    //     l2.push_back(randomizedPerfectMatching.truncatedWalk.getMax());
    //     // cout<<n<<" ";
    // }
    // for (int i = 1; i < 1000; i++) cout << i << ' ';
    // cout << endl;
    // for (int v : l1) cout << v << ' ';
    // cout << endl;
    // for (int v : l2) cout << v << ' ';
    // cout << endl;

    // cout<<endl;
    // int cmax = randomizedPerfectMatching.totalTimeConstant.getMax();
    // int cmax = 2;
    // cout << cmax << endl;
    // for (int i = 2; i <= 1000; i++){
    //     cout << cmax * i * log2(i) << ' ';
    // }
    // cout << endl;
    // randomizedPerfectMatching.totalTime.printAll();


    // randomizedPerfectMatching.totalSampling.printfreqs();
    // cout << "Variance = " << randomizedPerfectMatching.totalSampling.getVariance() << endl;
    // cout << "Mean = " << randomizedPerfectMatching.totalSampling.getMean() << endl;
    // cout << "Max = " << randomizedPerfectMatching.totalSampling.getMax() << endl;

    // for(int d = 1; d <= n; d += n / dCnt)
    // {
    //     inputGen.getRandomInput(n, d);

    //     randomizedPerfectMatching.initGraph(inputGen.n,inputGen.d,inputGen.adj);

    //     cout << "Running Truncated Algorithm" << endl;

    //     randomizedPerfectMatching.runAlgorithm();

    //     if (randomizedPerfectMatching.isPerfect()){
    //         cout << "M is perfect" << endl;
    //     }
    //     else cout << "M is not perfect" << endl;

    //     randomizedPerfectMatching.printRunStats();

    //     cout << "Running Non Truncated Algorithm" << endl;

    //     randomizedPerfectMatching.runNonTruncatedAlgorithm();

    //     if (randomizedPerfectMatching.isPerfect()){
    //         cout << "M is perfect" << endl;
    //     }
    //     else cout << "M is not perfect" << endl;

    //     randomizedPerfectMatching.printRunStats();
    // }

    // randomizedPerfectMatching.printAllRunStats();
}