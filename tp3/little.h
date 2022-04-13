#ifndef LITTLE_H
#define LITTLE_H

#include <stddef.h>

#define NBR_TOWNS 10

extern double* dist;
extern int starting_town[NBR_TOWNS];
extern int ending_town[NBR_TOWNS];
extern int* best_solution;
extern double best_eval;
extern float coord[NBR_TOWNS][2];

extern double* compute_distance(float coords[][2], size_t n_coords);
void print_matrix(double* d, size_t n_cities);
extern void print_solution(int* sol, double eval, size_t n_towns);
extern double evaluate(int* sol, size_t n_cities);

struct nn_t {
    double evaluation;
    int* solution;
};
extern struct nn_t build_nearest_neighbor_sub(double* dist, size_t n_towns);
extern double build_nearest_neighbor(double* dist, size_t n_towns);

extern void build_solution();
extern void little_algorithm(
    double* d0,
    int iteration,
    double eval_node_parent,
    size_t n_towns
);

#endif // LITTLE_H
