#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
using namespace std;

vector<int> calcPrefixSum ( vector<int> input, int num_threads)
{
	long long int n = input.size();
	long long int jump = 2, start = 0;
	omp_set_num_threads(num_threads);
	
	for (jump = 2; jump < n; jump *= 2)
	{
		#pragma omp parallel for schedule(static) firstprivate(jump, n)
			for(long long int j = start; j < n; j += jump)
			{
				input[j+jump/2] = input[j] + input[j+jump/2];
			}
		start = start + jump/2;
	}

	for (; jump > 1 ; jump /= 2)
	{
		#pragma omp parallel for schedule(static) firstprivate(jump, n)
			for(long long int j = start + jump/2; j < n; j += jump)
			{
				input[j] = input[j] + input[j-jump/2];
			}
		start = (start - jump/2)<0?start:(start - jump/2);
	}
	return input;
}