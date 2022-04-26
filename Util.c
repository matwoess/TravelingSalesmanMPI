#include <stdlib.h>
#include <stdbool.h>
#include <stdio.h>
#include "Util.h"

#define ERR_ALLOCATE_MEM -2
#define MAX_DISTANCE 100
static unsigned int SEED = 42;

int *get_random_edge_matrix(int nNodes, int populationPercentage) {
    int *edgeMatrix;
    allocate_int_array(&edgeMatrix, nNodes, nNodes);
    srand(SEED);
    for (int row = 0; row < nNodes; row++) {
        for (int col = 0; col < nNodes; col++) {
            if (rand_r(&SEED) % 100 < populationPercentage) {
                edgeMatrix[row * nNodes + col] = rand_r(&SEED) % MAX_DISTANCE;
            } else {
                edgeMatrix[row * nNodes + col] = 0;
            }
        }
    }
    return edgeMatrix;
}


void log_msg(bool verbose, char *s) {
    if (verbose) printf("%s\n", s);
}

void log_prune(bool verbose, int newDist, int currBest) {
    if (verbose)
        printf(
                "pruning path because %d is more than current best %d\n",
                newDist,
                currBest);
}

void log_curr_best_dist(bool verbose, int totalDist, int currBest) {
    if (verbose)
        printf(
                "totalDist: %d, currentBest %d\n",
                totalDist,
                currBest);
}

void logt_msg(bool verbose, int rank, char *s) {
    if (verbose) printf("T%d: %s\n", rank, s);
}

void logt_prune(bool verbose, int rank, int newDist, int currBest) {
    if (verbose)
        printf(
                "T%d: pruning path because %d is more than current best %d\n",
                rank,
                newDist,
                currBest);
}

void logt_curr_best_dist(bool verbose, int rank, int totalDist, int currBest) {
    if (verbose)
        printf(
                "T%d: totalDist: %d, currentBest %d\n",
                rank,
                totalDist,
                currBest);
}

void allocate_int_array(int **array, int rows, int columns) {
    int *a = (int *) malloc(rows * columns * sizeof(int));
    if (a == NULL) exit(ERR_ALLOCATE_MEM);
    *array = a;
}

void print_edge_matrix(int **edgeMatrix, int n) {
    printf("Edge Matrix:\n");
    for (int row = 0; row < n; row++) {
        printf("|");
        for (int col = 0; col < n; col++) {
            printf(" %1$2d", (*edgeMatrix)[row * n + col]);
        }
        printf("|\n");
    }
}

void print_path(bool verbose, int *path, int pathLength, int distance) {
    if (!verbose) return;
    for (int i = 0; i < pathLength; i++) {
        printf("%d", path[i]);
        if (i != pathLength - 1) printf("->");
    }
    printf(",%d\n", distance);
}

void printt_path(bool verbose, int rank, int *path, int pathLength, int distance) {
    if (!verbose) return;
    printf("T%d: ", rank);
    for (int i = 0; i < pathLength; i++) {
        printf("%d", path[i]);
        if (i != pathLength - 1) printf("->");
    }
    printf(",%d\n", distance);
}

void print_comm_buffer(int *commBuffer, int commBufferSize, int rank) {
    printf("T%d buffer: ", rank);
    for (int i = 0; i < commBufferSize; i++) {
        printf("%d ", commBuffer[i]);
    }
    printf("\n");
}