#include "iostream"
#include "vector"
#include "ctime"
#include "omp.h"

using namespace std;

extern vector<int> calcPrefixSum(vector<int> input, int num_threads);

int main() {
    int num_threads;
    cin >> num_threads;
    int n;
    cin >> n;
    vector<int> input;
    input.resize(n);
    for (int i = 0; i < n; i++) {
        cin >> input[i];
    }

    vector<int> prefixSum;
    double start_time = omp_get_wtime();
    prefixSum = calcPrefixSum(input, num_threads);
    double time_taken = omp_get_wtime() - start_time;

    // Printing stats and results
    cout<< time_taken << endl;
    cout<< prefixSum.size() << endl;
    n = prefixSum.size();
    for (int i = 0; i < n; i++){
        cout << prefixSum[i] << " " ;
    }
    cout << endl;

    return 0;
}