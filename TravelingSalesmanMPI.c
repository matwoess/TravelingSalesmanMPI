#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <limits.h>
#include <mpi.h>
#include <stdlib.h>
#include <time.h>
#include "Util.h"

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
void split_work(int *path);
void solve(int *path);
int add_node(int *path, int i);
void remove_node(int *path, int w);
void update_result(int *path);
int *get_path_from_manager();
void send_path_to_worker(int *path, int dest);
void listen_for_messages();
bool get_done_flag(int *buf);
void set_best_dist(int *buf, int distance);
void set_done_flag(int *buf, int flag);
int get_best_dist(int *buf);
bool all_threads_terminated();
void send_done_to_worker(int dest);
void send_result_to_manager();
void freeGlobals();
void send_beset_distance_to_worker(int dest);

bool prune = true;
bool verbose = false;

int N;
int *edgeMatrix;
int bestDistance;
int *bestPath;
int *paths;
int pathsInStack;

int nThreads, rank;
int *commBuffer;
int commBufferSize;
int doneFlag;
int workerThreadsTerminated;

static const int EXAMPLE_EDGES[][4] = {
        {0, 1,  3,  8},
        {5, 0,  2,  6},
        {1, 18, 0,  10},
        {7, 4,  12, 0},
};
#define EXAMPLE_N_NODES 4
#define ERR_INVALID_ARGS (-1)
#define OFFSET_PATH_LEN 1
#define OFFSET_PATH_DIST 2
#define MANAGER 0
#define OFFSET_BEST_DIST 3
#define OFFSET_DONE_FLAG 4

static const int TAG_REQUEST_PATH = 97;
static const int TAG_NEW_BEST = 98;
static const int TAG_BEST_DIST = 99;

int main(int argc, char *argv[]) {
    if (!parse_args(argc, argv)) return ERR_INVALID_ARGS;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nThreads);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) {
        print_edge_matrix(&edgeMatrix, N);
    }
    init_globals();
    clock_t t = clock();
    if (rank == 0) {
        int *path = init_path();
        add_path(path);
        free(path);
        while (true) {
            listen_for_messages();
            if (all_threads_terminated()) break;
        }
    } else {
        while (true) {
            int *path = get_path_from_manager();
            if (doneFlag) break;
            solve(path);
        }
    }
    t = clock() - t;
    double timeTaken = ((double) t) / CLOCKS_PER_SEC;
    if (rank == 0) {
        if (bestDistance == INT_MAX) {
            printf("No solution possible for current graph!\n");
            MPI_Finalize();
            return -1;
        }
        printf("\nBest path:\n");
        printt_path(true, rank, bestPath, N + 1, bestDistance);
        printf("\nAlgorithm took %.3fs\n", timeTaken);
    } else {
        logt_msg(true, rank, "thread exiting...");
    }
    freeGlobals();
    MPI_Finalize();
    return 0;
}

void init_globals() {
    allocate_int_array(&paths, N * (N - 1) / 2, N + 3);
    pathsInStack = 0;
    allocate_int_array(&bestPath, 1, N + 3);
    bestDistance = INT_MAX;
    commBufferSize = N + 5;
    allocate_int_array(&commBuffer, 1, commBufferSize);
    doneFlag = false;
    if (rank == 0) {
        workerThreadsTerminated = 0;
    }
}

bool all_threads_terminated() {
    return workerThreadsTerminated == nThreads - 1;
}

void freeGlobals() {
    free(paths);
    free(bestPath);
    free(commBuffer);
}

void listen_for_messages() {
    MPI_Status status;
    logt_msg(verbose, rank, "waiting for requests from workers...");
    MPI_Recv(commBuffer, commBufferSize, MPI_INT,
             MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
    if (status.MPI_TAG == TAG_REQUEST_PATH) {
        logt_msg(verbose, rank, "received request was for a new path...");
        if (doneFlag) {
            send_done_to_worker(status.MPI_SOURCE);
            workerThreadsTerminated++;
        } else {
            int *path = commBuffer;
            remove_path(path);
            if (pathsInStack == 0) {
                split_work(path);
                if (pathsInStack > 1) {
                    remove_path(path);
                }
            }
            send_path_to_worker(path, status.MPI_SOURCE);
            if (pathsInStack == 0) {
                doneFlag = true;
            }
        }
    } else if (status.MPI_TAG == TAG_NEW_BEST) {
        logt_msg(verbose, rank, "received request was for a new best score...");
        int newDist = get_best_dist(commBuffer);
        logt_curr_best_dist(verbose, rank, newDist, bestDistance);
        if (newDist < bestDistance) {
            logt_msg(verbose, rank, "new best. copying to bestPath...");
            memcpy(bestPath, commBuffer, (N + 3) * sizeof(int));
            bestDistance = newDist;
        }
        send_beset_distance_to_worker(status.MPI_SOURCE);
    } else {
        logt_msg(verbose, rank, "unknown tag received!");
    }
}

void send_done_to_worker(int dest) {
    logt_msg(verbose, rank, "informing worker that work is done.");
    set_best_dist(commBuffer, bestDistance);
    set_done_flag(commBuffer, doneFlag);
    MPI_Send(commBuffer, commBufferSize, MPI_INT,
             dest, TAG_REQUEST_PATH, MPI_COMM_WORLD);
}

void send_path_to_worker(int *path, int dest) {
    logt_msg(verbose, rank, "sending path to worker...");
    printt_path(verbose, rank, path, get_path_length(path), get_path_dist(path));
    set_best_dist(commBuffer, bestDistance);
    set_done_flag(commBuffer, doneFlag);
    MPI_Send(commBuffer, commBufferSize, MPI_INT,
             dest, TAG_REQUEST_PATH, MPI_COMM_WORLD);
}

void send_beset_distance_to_worker(int dest) {
    logt_msg(verbose, rank, "sending back best distance by now...");
    set_best_dist(commBuffer, bestDistance);
    set_done_flag(commBuffer, doneFlag);
    MPI_Send(commBuffer, commBufferSize, MPI_INT,
             dest, TAG_BEST_DIST, MPI_COMM_WORLD);
}

int *get_path_from_manager() {
    logt_msg(verbose, rank, "requesting path from manager...");
    MPI_Send(commBuffer, commBufferSize, MPI_INT,
             MANAGER, TAG_REQUEST_PATH, MPI_COMM_WORLD);
    MPI_Status status;
    logt_msg(verbose, rank, "waiting for path from manager...");
    MPI_Recv(commBuffer, commBufferSize, MPI_INT,
             MANAGER, TAG_REQUEST_PATH, MPI_COMM_WORLD, &status);
    doneFlag = get_done_flag(commBuffer);
    if (doneFlag) {
        logt_msg(verbose, rank, "received work is done. exiting...");
        return NULL;
    }
    logt_msg(verbose, rank, "received path from manager.");
    int *path = commBuffer;
    bestDistance = get_best_dist(path);
    return path;
}

void send_result_to_manager() {
    logt_msg(verbose, rank, "sending path to manager...");
    set_best_dist(commBuffer, bestDistance);
    MPI_Send(commBuffer, commBufferSize, MPI_INT,
             MANAGER, TAG_NEW_BEST, MPI_COMM_WORLD);
    MPI_Status status;
    logt_msg(verbose, rank, "waiting for best distance from manager...");
    MPI_Recv(commBuffer, commBufferSize, MPI_INT,
             MANAGER, TAG_BEST_DIST, MPI_COMM_WORLD, &status);
    bestDistance = get_best_dist(commBuffer);
}

int *init_path() {
    int *initialPath;
    allocate_int_array(&initialPath, 1, N + 3);
    initialPath[0] = 0;
    set_path_length(initialPath, 1);
    set_path_dist(initialPath, 0);
    return initialPath;
}


void add_path(int *path) {
    memcpy(&paths[pathsInStack * (N + 3)], path, (N + 3) * sizeof(int));
    pathsInStack++;
}

void remove_path(int *path) {
    pathsInStack--;
    memcpy(path, &paths[pathsInStack * (N + 3)], (N + 3) * sizeof(int));
}

void split_work(int *path) {
    logt_msg(verbose, rank, "expanding path: ");
    printt_path(verbose, rank, path, get_path_length(path), get_path_dist(path));
    for (int i = 0; i < N; i++) {
        int w = add_node(path, i);
        if (w < 0) continue;
        add_path(path);
        remove_node(path, w);
    }
}

void solve(int *path) {
    add_path(path);
    while (pathsInStack > 0) {
        remove_path(path);
        printt_path(verbose, rank, path, get_path_length(path), get_path_dist(path));
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
}

int add_node(int *path, int i) {
    int currPathNode = get_last_node(path);
    if (i == currPathNode) return -1;
    if (is_already_in_path(path, i)) return -1;
    int dist = edgeMatrix[currPathNode * N + i];
    if (dist == 0) return -1;
    int newDist = get_path_dist(path) + dist;
    if (prune && newDist > bestDistance) {
        logt_prune(verbose, rank, newDist, bestDistance);
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
        logt_msg(verbose, rank, "missing way back to 0!");
        return;
    }
    path[N] = 0;
    int totalDist = get_path_dist(path) + distTo0;
    set_path_dist(path, totalDist);
    set_path_length(path, N + 1);
    printt_path(verbose, rank, path, get_path_length(path), get_path_dist(path));
    logt_curr_best_dist(verbose, rank, totalDist, bestDistance);
    if (totalDist < bestDistance) {
        logt_msg(verbose, rank, "new best!");
        bestDistance = totalDist;
        memcpy(bestPath, path, (N + 3) * sizeof(int));
        send_result_to_manager();
    }
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

int get_path_length(int *p) {
    return p[N + OFFSET_PATH_LEN];
}

void set_path_length(int *p, int len) {
    p[N + OFFSET_PATH_LEN] = len;
}

int get_path_dist(int *p) {
    return p[N + OFFSET_PATH_DIST];
}

void set_path_dist(int *p, int dist) {
    p[N + OFFSET_PATH_DIST] = dist;
}

int get_last_node(int *path) {
    return path[get_path_length(path) - 1];
}

bool is_already_in_path(int *path, int node) {
    for (int i = 0; i < get_path_length(path); i++)
        if (path[i] == node) return true;
    return false;
}


int get_best_dist(int *buf) {
    return buf[N + OFFSET_BEST_DIST];
}

void set_best_dist(int *buf, int distance) {
    buf[N + OFFSET_BEST_DIST] = distance;
}

bool get_done_flag(int *buf) {
    return buf[N + OFFSET_DONE_FLAG];
}

void set_done_flag(int *buf, int flag) {
    buf[N + OFFSET_DONE_FLAG] = flag;
}