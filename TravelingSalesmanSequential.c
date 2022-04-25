#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <limits.h>
#include <mpi.h>
#include <stdlib.h>
#include <time.h>
#include "Util.h"

#define ERR_INVALID_ARGS -1

bool parse_args(int argc, char **argv);

void init_globals();

int get_path_length(int *p);

void set_path_length(int *p, int len);

int get_path_dist(int *p);

void set_path_dist(int *p, int dist);

int get_last_node(int *path);

int *init_path();

bool is_already_in_path(int *path, int node);

void add_path(int *path);

void remove_path(int *path);

void solve();

int add_node(int *path, int i);

void remove_node(int *path, int w);

void update_result(int *path);

bool prune = true;
bool verbose = false;

int N;
int *edgeMatrix;
int bestDistance;
int *bestPath;
int *paths;
int pathsInStack;


static const int EXAMPLE_EDGES[][4] = {
        {0, 1,  3,  8},
        {5, 0,  2,  6},
        {1, 18, 0,  10},
        {7, 4,  12, 0},
};
#define EXAMPLE_N_NODES 4

int main(int argc, char *argv[]) {
    if (!parse_args(argc, argv)) return ERR_INVALID_ARGS;
    print_edge_matrix(&edgeMatrix, N);
    init_globals();
    //MPI_Init(&argc, &argv);
    time_t t1 = time(NULL);
    solve();
    time_t t2 = time(NULL);
    if (bestDistance == INT_MAX) {
        printf("No solution possible for current graph!\n");
        return -1;
    }
    printf("\nBest path:\n");
    print_path(true, bestPath, N + 1, bestDistance);
    printf("\nAlgorithm took %.3fs\n", (double) (t2 - t1));
    return 0;
}

bool parse_args(int argc, char *argv[]) {
    bool applyExample = argc >= 2 && strcmp(argv[1], "example") == 0;
    if (applyExample) {
        edgeMatrix = (int *) EXAMPLE_EDGES;
        N = EXAMPLE_N_NODES;
    } else if (argc < 3) {
        printf("Not enough arguments!\n");
        printf("Usage: <program> <nNodes> <%%population> [-noprune] [-verbose]\n");
        printf("Alternative: <program> \"example\" [-noprune] [-verbose]\n");
        return false;
    } else {
        char *endptr;
        int nodesArg = (int) strtol(argv[1], &endptr, 10);
        int populationArg = (int) strtol(argv[2], &endptr, 10);
        N = nodesArg;
        edgeMatrix = get_random_edge_matrix(nodesArg, populationArg);
    }
    for (int a = 0; a < argc; a++) {
        if (strcmp(argv[a], "-noprune") == 0) {
            prune = false;
        } else if (strcmp(argv[a], "-verbose") == 0) {
            verbose = true;
        }
    }
    return true;
}


void init_globals() {
    allocate_int_array(&paths, N * (N - 1) / 2, N + 3);
    pathsInStack = 0;
    allocate_int_array(&bestPath, 1, N + 3);
    bestDistance = INT_MAX;
}

int get_path_length(int *p) {
    return p[N + 1];
}

void set_path_length(int *p, int len) {
    p[N + 1] = len;
}

int get_path_dist(int *p) {
    return p[N + 2];
}

void set_path_dist(int *p, int dist) {
    p[N + 2] = dist;
}

int get_last_node(int *path) {
    return path[get_path_length(path) - 1];
}

int *init_path() {
    int *initialPath;
    allocate_int_array(&initialPath, 1, N + 3);
    initialPath[0] = 0;
    set_path_length(initialPath, 1);
    set_path_dist(initialPath, 0);
    return initialPath;
}

bool is_already_in_path(int *path, int node) {
    for (int i = 0; i < get_path_length(path); i++)
        if (path[i] == node) return true;
    return false;
}

void add_path(int *path) {
    memcpy(&paths[pathsInStack * (N + 3)], path, (N + 3) * sizeof(int));
    pathsInStack++;
}

void remove_path(int *path) {
    pathsInStack--;
    memcpy(path, &paths[pathsInStack * (N + 3)], (N + 3) * sizeof(int));
}


void solve() {
    int *path = init_path();
    add_path(path);
    while (pathsInStack > 0) {
        remove_path(path);
        print_path(verbose, path, get_path_length(path), get_path_dist(path));
        if (get_path_length(path) == N) {
            update_result(path);
            continue;
        }
        for (int i = 0; i < N; i++) {
            int w = add_node(path, i);
            if (w < 0) continue;
            add_path(path);
            remove_node(path, w);
        }
    }
    free(path);
}

int add_node(int *path, int i) {
    int currPathNode = get_last_node(path);
    if (i == currPathNode) return -1;
    if (is_already_in_path(path, i)) return -1;
    int dist = edgeMatrix[currPathNode * N + i];
    if (dist == 0) return -1;
    int newDist = get_path_dist(path) + dist;
    if (prune && newDist > bestDistance) {
        log_prune(verbose, newDist, bestDistance);
        return -1;
    }
    int pathLength = get_path_length(path);
    path[pathLength] = i;
    set_path_length(path, pathLength + 1);
    set_path_dist(path, newDist);
    return dist;
}

void remove_node(int *path, int w) {
    set_path_dist(path, get_path_dist(path) - w);
    set_path_length(path, get_path_length(path) - 1);
}

void update_result(int *path) {
    int currPathNode = get_last_node(path);
    int distTo0 = edgeMatrix[currPathNode * N + 0];
    if (distTo0 == 0) {
        log_msg(verbose, "missing way back to 0!");
        return;
    }
    path[N] = 0;
    int totalDist = get_path_dist(path) + distTo0;
    set_path_dist(path, totalDist);
    set_path_length(path, N + 1);
    print_path(verbose, path, get_path_length(path), get_path_dist(path));
    log_curr_best_dist(verbose, totalDist, bestDistance);
    if (totalDist < bestDistance) {
        log_msg(verbose, "new best!");
        bestDistance = totalDist;
        memcpy(bestPath, path, (N + 3) * sizeof(int));
    }
}

