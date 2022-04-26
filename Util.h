#ifndef TRAVELINGSALESMAN_UTIL_H
#define TRAVELINGSALESMAN_UTIL_H

void allocate_int_array(int **array, int rows, int columns);
void print_edge_matrix(int **edgeMatrix, int n);
int *get_random_edge_matrix(int nNodes, int populationPercentage);
void log_msg(bool verbose, char *s);
void log_prune(bool verbose, int newDist, int currBest);
void log_curr_best_dist(bool verbose, int totalDist, int currBest);
void logt_msg(bool verbose, int rank, char *s);
void logt_prune(bool verbose, int rank, int newDist, int currBest);
void logt_curr_best_dist(bool verbose, int rank, int totalDist, int currBest);
void print_path(bool verbose, int* path, int pathLength, int distance);
void printt_path(bool verbose, int rank, int* path, int pathLength, int distance);
void print_comm_buffer(int *commBuffer, int commBufferSize, int rank);

#endif //TRAVELINGSALESMAN_UTIL_H