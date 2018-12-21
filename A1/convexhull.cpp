#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <vector>
#include <math.h>
#include <iostream>

using namespace std;


void findhull(int, int, vector<int> );

vector<pair<int, int> > hull;
vector<pair< pair<int, int> , int> > pts;
int size;
int width, height;
int NUM_THREADS;
int counter;
// Using a 0 indexed Cartesian coordinate system


vector< pair<int, int> > calcConvexHull(vector< vector<int> > image, int num_threads)
{
	width = image[0].size();
	height = image.size();
	NUM_THREADS = num_threads;
	counter = 1;
	size = 0;
	for (int j = 0; j < width; ++j)
	{
		for (int i = 0; i < height; ++i)
		{
			if(image[i][j]==1)
			{
				pts.push_back(make_pair(make_pair(i, j), 0));	// first argument is y, second is x
				size++;
			}
		}
	}
	hull.push_back(make_pair(pts[0].first.first, pts[0].first.second));
	hull.push_back(make_pair(pts[size-1].first.first, pts[size-1].first.second));
	pts[0].second = INT_MAX;
	pts[size-1].second = INT_MAX;
	

	int j1 = pts[0].first.first;
	int i1 = pts[0].first.second;
	int j2 = pts[size-1].first.first;
	int i2 = pts[size-1].first.second;
	
	
	vector<int> greater, less;
	
	for (int j = 1; j < size-1; ++j)
	{
		int j3 = pts[j].first.first;
		int i3 = pts[j].first.second;
		if( (i2 - i1) * (j3 - j1) - (j2 - j1) * (i3 - i1) > 0) // (Bx - Ax) * (Y - Ay) - (By - Ay) * (X - Ax)
		{
			greater.push_back(j);
		}
		else if ((i2 - i1) * (j3 - j1) - (j2 - j1) * (i3 - i1) < 0)
		{
			less.push_back(j);
		}
	}
	omp_set_nested(1);
	
	bool required = (greater.size() > 0) && (less.size() > 0);
	if((num_threads > 1) && required)
	{
		counter++;
		#pragma omp parallel for num_threads(2)
		for (int i = 0; i < 2; ++i)
		{
			if(i==0 && greater.size() > 0)
				findhull(0, size-1, greater);
			if(i==1 && less.size() > 0)
				findhull(size-1, 0, less);
		}
	}
	else
	{
		if(greater.size() > 0)
			findhull(0, size-1, greater);
		if(less.size() > 0)
			findhull(size-1, 0, less);
	}

	return hull;
}

void findhull(int a, int b, vector<int> set_pts)
{
	int j1 = pts[a].first.first;
	int i1 = pts[a].first.second;
	int j2 = pts[b].first.first;
	int i2 = pts[b].first.second;
	int maxpt = 0;
	float maxdist = 0;
	float s12 = (i1==i2)?INT_MAX:(float)(j2-j1)/(float)(i2-i1);
	float c = (float)j1 - s12*(float)(i1);
	
	int n = set_pts.size();
	
	for(int i=0; i < n; i++)
	{
		int j3 = pts[set_pts[i]].first.first;
		int i3 = pts[set_pts[i]].first.second;
		if(pts[set_pts[i]].second == INT_MAX) // i.e it is already included in the hull
			continue;
		float dist = (i1==i2)?abs(i3-i1):abs((float)(j3 - s12*(float)(i3) - c)/(float)(sqrt(s12*s12 + 1)));

		if(dist > maxdist)
		{
			maxdist = dist;
			maxpt = set_pts[i];
		}
	}
	if(pts[maxpt].second != INT_MAX) // if not already included in the hull
	{
		#pragma omp critical
		{
			hull.push_back(make_pair(pts[maxpt].first.first, pts[maxpt].first.second));
		}
		pts[maxpt].second = INT_MAX;
	}
	
	int j4 = pts[maxpt].first.first;
	int i4 = pts[maxpt].first.second;
	vector<int> one, two;
	for (int j = 0; j < n; ++j)
	{
		int j3 = pts[set_pts[j]].first.first;
		int i3 = pts[set_pts[j]].first.second;
		if( (i4 - i1) * (j3 - j1) - (j4 - j1) * (i3 - i1) > 0) // (Bx - Ax) * (Y - Ay) - (By - Ay) * (X - Ax)
		{
			one.push_back(set_pts[j]);
		}
		if ((i2 - i4) * (j3 - j4) - (j2 - j4) * (i3 - i4) > 0)
		{
			two.push_back(set_pts[j]);
		}
	}

	bool canAllocate = false;
	bool required = (one.size() > 0) && (two.size() > 0);
	#pragma omp critical
	{
		if(counter < NUM_THREADS)
		{
			counter++;
			canAllocate = true;
		}
	}
	if(canAllocate && required)
	{
		#pragma omp parallel for num_threads(2)
		for (int i = 0; i < 2; ++i)
		{
			if(i==0 && one.size()>0)
				findhull(a, maxpt, one);
			if(i==1 && two.size()>0)
				findhull(maxpt, b, two);
		}
	}
	else
	{
		if(one.size() > 0)
			findhull(a, maxpt, one);
		if(two.size() > 0)
			findhull(maxpt, b, two);
	}
	if(canAllocate)
	{
		counter--;
	}
}