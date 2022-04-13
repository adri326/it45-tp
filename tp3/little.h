#ifndef LITTLE_H
#define LITTLE_H

#define NBR_TOWNS 6

extern double dist[NBR_TOWNS][NBR_TOWNS];
extern int starting_town[NBR_TOWNS];
extern int ending_town[NBR_TOWNS];
extern int best_solution[NBR_TOWNS];
extern double best_eval;
extern float coord[NBR_TOWNS][2];

extern double* compute_distance(float coords[][2], size_t n_coords);
extern void print_matrix(double d[NBR_TOWNS][NBR_TOWNS]);
extern void print_solution(int* sol, double eval);
extern double evaluate(int* sol);

struct nn_t {
    double evaluation;
    int* solution;
};
extern struct nn_t build_nearest_neighbor_sub(double* dist, size_t n_towns);
extern double build_nearest_neighbor(double* dist, size_t n_towns);

extern void build_solution();
extern void little_algorithm(
    double d0[NBR_TOWNS][NBR_TOWNS],
    int iteration,
    double eval_node_parent
);

#endif // LITTLE_H
