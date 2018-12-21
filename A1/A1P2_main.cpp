#include "iostream"
#include "vector"
#include "ctime"
#include "omp.h"

using namespace std;

extern vector< pair<int, int> > calcConvexHull(vector< vector<int> > image, int num_threads);

int main() {
    string file_headers;
    getline(cin, file_headers);
    getline(cin, file_headers);

    int num_threads;
    cin >> num_threads;
    int m, n;
    cin >> m >> n;
    vector< vector<int> > image;
    image.resize(n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            char a;
            cin >> a;
            image[i].push_back((int)a-'0');
        }
    }

    vector< pair<int, int> > convexHull;
    double start_time = omp_get_wtime();
    convexHull = calcConvexHull(image, num_threads);
    double time_taken = omp_get_wtime() - start_time;

    // Printing stats and results
    cout<< time_taken << endl;
    cout<< convexHull.size() << endl;

    for (int i = 0; i < (int)convexHull.size(); i++){
        cout << convexHull[i].first << " " << convexHull[i].second << endl;
    }

    return 0;
}