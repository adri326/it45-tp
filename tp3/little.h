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
extern double evaluation_solution(int* sol);
extern double build_nearest_neighbour();
extern void build_solution();
extern void little_algorithm(double d0[NBR_TOWNS][NBR_TOWNS], int iteration, double eval_node_parent);

#endif // LITTLE_H
