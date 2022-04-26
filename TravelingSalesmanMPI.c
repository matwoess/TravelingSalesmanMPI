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

void solve(int *path);

int add_node(int *path, int i);

void remove_node(int *path, int w);

void update_result(int *path);

int *get_path_from_manager();

void *send_path_to_worker(int *path, int dest);

void receive_result_from_workers();

void print_comm_buffer();

bool get_done_flag(int *buf);

void set_best_dist(int *buf, int distance);

void set_done_flag(int *buf, int flag);

int get_best_dist(int *buf);

bool prune = true;
bool verbose = false;

int N;
int *edgeMatrix;
int bestDistance;
int *bestPath;
int *paths;
int pathsInStack;
int doneFlag;

int size, rank;
int *commBuffer;
int commBufferSize;
int *threadsIdle;

static const int EXAMPLE_EDGES[][4] = {
        {0, 1,  3,  8},
        {5, 0,  2,  6},
        {1, 18, 0,  10},
        {7, 4,  12, 0},
};
#define EXAMPLE_N_NODES 4
#define MANAGER 0

static const int TAG_NEW_PATH = 99;
static const int TAG_NEW_BEST = 100;

int main(int argc, char *argv[]) {
    if (!parse_args(argc, argv)) return ERR_INVALID_ARGS;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) {
        print_edge_matrix(&edgeMatrix, N);
    }
    init_globals();
    time_t t1 = time(NULL);
    if (rank == 0) {
        int *path = init_path();
        pathsInStack = 1;
        send_path_to_worker(path, 1);
        do {
            receive_result_from_workers();
            if (pathsInStack == 0) {
                doneFlag = true;
                send_path_to_worker(path, 1);
                break;
            }
        } while (pathsInStack > 0);
    } else {
        while (true) {
            int *path = get_path_from_manager();
            if (doneFlag) break;
            solve(path);
        }
    }
    time_t t2 = time(NULL);
    if (rank == 0) {
        if (bestDistance == INT_MAX) {
            printf("No solution possible for current graph!\n");
            MPI_Finalize();
            return -1;
        }
        printf("\nBest path:\n");
        printt_path(true, rank, bestPath, N + 1, bestDistance);
        printf("\nAlgorithm took %.3fs\n", (double) (t2 - t1));
    }
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
        allocate_int_array(&threadsIdle, 1, N);
    }
}

void *send_path_to_worker(int *path, int dest) {
    logt_msg(verbose, rank, "sending path to worker...");
    memcpy(commBuffer, path, (N + 3) * sizeof(int));
    set_best_dist(commBuffer, bestDistance);
    set_done_flag(commBuffer, doneFlag);
    MPI_Send(commBuffer, commBufferSize, MPI_INT, dest, TAG_NEW_PATH, MPI_COMM_WORLD);
    logt_msg(verbose, rank, "sent path to worker.");
}

int *get_path_from_manager() {
    int *path;
    allocate_int_array(&path, 1, N + 3);
    MPI_Status status;
    logt_msg(verbose, rank, "waiting for path from manager...");
    MPI_Recv(commBuffer, commBufferSize, MPI_INT, MANAGER, TAG_NEW_PATH, MPI_COMM_WORLD, &status);
    logt_msg(verbose, rank, "received path from manager.");
    memcpy(path, commBuffer, (N + 3) * sizeof(int));
    printf("len of received path: %d\n", get_path_length(path));
    printf("dist of received path: %d\n", get_path_dist(path));
    bestDistance = get_best_dist(commBuffer);
    doneFlag = get_done_flag(commBuffer);
    return path;
}

void *send_result_to_manager(int *path) {
    logt_msg(verbose, rank, "sending path to manager...");
    memcpy(commBuffer, &path, (N + 3) * sizeof(int));
    set_best_dist(commBuffer, bestDistance);
    set_done_flag(commBuffer, doneFlag);
    MPI_Send(commBuffer, commBufferSize, MPI_INT, MANAGER, TAG_NEW_BEST, MPI_COMM_WORLD);
    logt_msg(verbose, rank, "sent path to manger.");
}

void receive_result_from_workers() {
    MPI_Status status;
    logt_msg(verbose, rank, "waiting for answer from workers...");
    MPI_Recv(commBuffer, commBufferSize, MPI_INT, MPI_ANY_SOURCE, TAG_NEW_BEST, MPI_COMM_WORLD, &status);
    logt_msg(verbose, rank, "received answer from worker.");
    int newDist = get_best_dist(commBuffer);
    if (newDist < bestDistance) {
        memcpy(bestPath, commBuffer, (N + 3) * sizeof(int));
        bestDistance = newDist;
    }
    if (get_done_flag(commBuffer) == true) {
        threadsIdle[status.MPI_SOURCE] = true;
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
    free(path);
    doneFlag = true;
    send_result_to_manager(bestPath);
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
        send_result_to_manager(bestPath);
    }
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

bool is_already_in_path(int *path, int node) {
    for (int i = 0; i < get_path_length(path); i++)
        if (path[i] == node) return true;
    return false;
}


int get_best_dist(int *buf) {
    return buf[N + 3];
}

void set_best_dist(int *buf, int distance) {
    buf[N + 3] = distance;
}

bool get_done_flag(int *buf) {
    return buf[N + 4];
}

void set_done_flag(int *buf, int flag) {
    buf[N + 4] = flag;
}


void print_comm_buffer() {
    printf("T%d buffer: ", rank);
    for (int i = 0; i < commBufferSize; i++) {
        printf("%d ", commBuffer[i]);
    }
    printf("\n");
}