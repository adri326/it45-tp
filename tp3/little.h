#ifndef LITTLE_H
#define LITTLE_H

#include <stddef.h>
#include <stdbool.h>

#define MAX_CITIES 52

extern float berlin52_coords[][2];

struct solution_t {
    double evaluation;
    int* solution;
};
typedef struct solution_t solution_t;

extern double* compute_distance(float coords[][2], size_t n_coords);
void print_matrix(double* d, size_t n_cities);
extern void print_solution(solution_t solution, size_t n_towns);
extern void free_solution(solution_t solution);
extern double evaluate(const double* dist, int* sol, size_t n_cities);

extern solution_t build_nearest_neighbor_sub(const double* dist, size_t n_towns);
extern double build_nearest_neighbor(const double* dist, size_t n_towns);

extern bool is_valid_full_path(const int* path, size_t n_cities);
extern bool is_valid_sub_path(const int* path, size_t iteration, size_t n_cities);

extern void build_solution();
extern void little_algorithm_rec(
    const double* g_dist,
    double* current_dist,
    size_t iteration,
    double eval_node_parent,
    size_t n_towns,
    bool little_plus
);

extern solution_t little_algorithm(
    const double* dist,
    size_t n_towns
);

#endif // LITTLE_H
