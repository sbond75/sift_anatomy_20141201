/*
IPOL SIFT
Copyright (C) 2014, Ives Rey-Otero, CMLA ENS Cachan
<ives.rey-otero@cmla.ens-cachan.fr>

Version 20140911 (September 11th, 2014)

== Patent Warning and License =================================================

The SIFT method is patented 

    [3] "Method and apparatus for identifying scale invariant features
      in an image."
        David G. Lowe
        Patent number: 6711293
        Filing date: Mar 6, 2000
        Issue date: Mar 23, 2004
        Application number: 09/519,89
  
 These source codes are made available for the exclusive aim of serving as
 scientific tool to verify the soundness and completeness of the algorithm
 description. Compilation, execution and redistribution of this file may
 violate patents rights in certain countries. The situation being different
 for every country and changing over time, it is your responsibility to
 determine which patent rights restrictions apply to you before you compile,
 use, modify, or redistribute this file. A patent lawyer is qualified to make
 this determination. If and only if they don't conflict with any patent terms,
 you can benefit from the following license terms attached to this file.


This program is free software: you can use, modify and/or
redistribute it under the terms of the simplified BSD
License. You should have received a copy of this license along
this program. If not, see
<http://www.opensource.org/licenses/bsd-license.html>.

*/
/**
 * @file sift_matching.c
 * @brief data structures to store information relative to a pair of keypoints
 *
 * @li struct keypointPr     : Pair of keypoint data structure.
 * @li struct keypointPr_list : List of pairs.
 * @li print,save, read for lists of pairs.
 *
 * @author Ives Rey-Otero <ives.rey-otero@cmla.ens-cachan.fr>
 */

#include <algorithm>
#include <vector>
#include <cfloat>
#include <iostream>

extern "C" {

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "lib_keypoint.h"
#include "lib_matching.h"
#include "lib_util.h"

static void compute_keypoints_distance(float* dist,
                                       const struct sift_keypoints *k1,
                                       const struct sift_keypoints *k2)
{
    int n_hist = k1->list[0]->n_hist;
    int n_ori  = k1->list[0]->n_ori;
    int dim = n_hist*n_hist*n_ori;
    int n1 = k1->size;
    int n2 = k2->size;
    for(int i = 0; i < n1; i++){
        const float * v1 =  k1->list[i]->descr; // An array of descriptors
        for(int j = 0; j < n2; j++){
            const float * v2 =  k2->list[j]->descr; // An array of descriptors
            float d = euclidean_distance(v1, v2, dim /* length of v1 and v2 arrays */);
            dist[i*n2+j] = d;
        }
    }
}


static void find_the_two_nearest_keys(const float* dist, int n1, int n2,
                                     int* indexA, int* indexB,
                                     float* distA, float* distB)
{
    for(int i = 0; i < n1; i++){
        int iA, iB;
        float dA, dB;
        find_array_two_min(&dist[i*n2], n2, &dA, &dB, &iA, &iB);
        indexA[i] = iA;
        indexB[i] = iB;
        distA[i] = dA;
        distB[i] = dB;
    }
}

} // End extern "C"
template<typename F>
void forEachMatch(F f,
                 struct sift_keypoints *k1, // Keypoints from the first image
                 struct sift_keypoints *k2, // Keypoints from the second image
                 float thresh,
                 int flag,
                 float* distA,
                 float* distB,
                 int* indexA,
                 int* indexB
) {
    int n1 = k1->size;
    int n2 = k2->size;

    double average = 0;
    double stddev = 0;
    size_t processed = 0;
    double min = DBL_MAX, max = DBL_MIN;

    for(int i = 0; i < n1; i++){
        float val;
        val = (flag == 1 ? distA[i]/distB[i] : distA[i]);
        if (val < thresh){
            // The indices of the matching keypoints:
            int iA = indexA[i];
            int iB = indexB[i];
            // Points ptA and ptB are matching points:
            struct keypoint* k;
            struct keypoint* ptA = k2->list[iA]; // The second matched keypoint
            struct keypoint* ptB = k2->list[iB]; // Unknown meaning of this point
            
            f(k1->list[i], ptA, ptB);
            
            processed++;
        }
    }
}

extern "C" {
void matching(struct sift_keypoints *k1, // Keypoints from the first image
              struct sift_keypoints *k2, // Keypoints from the second image
              struct sift_keypoints *out_k1, // Matched keypoints in first image
              struct sift_keypoints *out_k2A, // Matched keypoints in second image
              struct sift_keypoints *out_k2B,
              float thresh,
              int flag)
{
    int n1 = k1->size;
    int n2 = k2->size;

    float* dist  = (float*)xmalloc(n1*n2*sizeof(float)); // Holds all possible distances checked across both arrays k1->list and k2->list
    float* distA = (float*)xmalloc(n1*sizeof(float));
    float* distB = (float*)xmalloc(n1*sizeof(float));
    int* indexA  = (int*)xmalloc(n1*sizeof(int));
    int* indexB  = (int*)xmalloc(n1*sizeof(int));

    compute_keypoints_distance(dist, k1, k2); // Find the distance between all keypoints in k1 to those in k2, for every possible pair that can be made across both arrays
    find_the_two_nearest_keys(dist, n1, n2, indexA, indexB, distA, distB);

//    std::vector<float> distances;
//    std::transform(k2->list->, std::back_inserter(Y), [](const std::vector<int>& value) {
//        return value.size();
//    });
//    cv::meanStdDev(
    // Compute average
    double average = 0;
    double min = DBL_MAX, max = DBL_MIN;
    size_t processed = 0;
#define DIST(ptA, pt1) sqrt(pow(ptA->x - pt1->x, 2) + pow(ptA->y - pt1->y, 2))
    forEachMatch([&](struct keypoint* pt1, struct keypoint* ptA, struct keypoint* ptB) {
        double dist = DIST(ptA, pt1);
        
        average += dist;
        
        // Maintain min and max
        if (dist < min) {
            min = dist;
        }
        if (dist > max) {
            max = dist;
        }
        
        processed++;
    }, k1, k2, thresh, flag, distA, distB, indexA, indexB);
    average /= processed;
    
    // Compute standard deviation
    double stddev = 0;
    processed = 0;
    forEachMatch([&](struct keypoint* pt1, struct keypoint* ptA, struct keypoint* ptB) {
        double dist = DIST(ptA, pt1);
        
        stddev += pow(dist - average, 2);
        
        processed++;
    }, k1, k2, thresh, flag, distA, distB, indexA, indexB);
    stddev /= processed;
    stddev = sqrt(stddev);
    double threshold = stddev * 3;
//    double threshold = stddev * 6;
    
    // Save matches
    forEachMatch([&](struct keypoint* pt1, struct keypoint* ptA, struct keypoint* ptB) {
        double dist = DIST(ptA, pt1);

        /* if (sqrt(pow(ptA->x - ptB->x, 2) + pow(ptA->y - ptB->y, 2)) > 60) { //100) { */
        /*     continue; */
        /* } */
        if (abs(dist - average) > threshold) {
            printf("Discarding %f, greater than %f\n", dist, threshold);
            return;
        }
        struct keypoint* k = sift_malloc_keypoint_from_model_and_copy(pt1); // The first matched keypoint
        sift_add_keypoint_to_list(k, out_k1); // The first matched keypoints
        k = sift_malloc_keypoint_from_model_and_copy(ptA);
        sift_add_keypoint_to_list(k, out_k2A); // The second matched keypoints
        k = sift_malloc_keypoint_from_model_and_copy(ptB);
        sift_add_keypoint_to_list(k, out_k2B); // Unknown meaning of this point
        
        processed++;
    }, k1, k2, thresh, flag, distA, distB, indexA, indexB);

    free(dist);
    free(indexA);
    free(indexB);
    free(distA);
    free(distB);
}

void print_pairs(const struct sift_keypoints *k1,
                 const struct sift_keypoints *k2)
{
    if (k1->size > 0){
        int n = k1->size;
        for(int i = 0; i < n ;i++){
            fprintf_one_keypoint(stdout, k1->list[i], 0, 0, 0);
            fprintf_one_keypoint(stdout, k2->list[i], 0, 0, 0);
            fprintf(stdout, "\n");
        }
    }
}

void save_pairs_extra(const char* name,
                      const struct sift_keypoints *k1,
                      const struct sift_keypoints *k2A,
                      const struct sift_keypoints *k2B)
{
    FILE* f = fopen(name,"w");
    
    if (k1->size > 0){

        int n_hist = k1->list[0]->n_hist;
        int n_ori = k1->list[0]->n_ori;
        int dim = n_hist*n_hist*n_ori;
        int n_bins  = k1->list[0]->n_bins;
        int n = k1->size;
        for(int i = 0; i < n; i++){
            fprintf_one_keypoint(f, k1->list[i], dim, n_bins, 2);
            fprintf_one_keypoint(f, k2A->list[i], dim, n_bins, 2);
            fprintf_one_keypoint(f, k2B->list[i], dim, n_bins, 2);
            fprintf(f, "\n");
        }
    }
    fclose(f);
}

}
