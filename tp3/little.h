#ifndef LITTLE_H
#define LITTLE_H

#define NBR_TOWNS 6

static double dist[NBR_TOWNS][NBR_TOWNS];
static int starting_town[NBR_TOWNS];
static int ending_town[NBR_TOWNS];
static int best_solution[NBR_TOWNS];
static double best_eval;
static float coord[NBR_TOWNS][2];

extern void print_matrix(double d[NBR_TOWNS][NBR_TOWNS]);
extern void print_solution(int* sol, double eval);
extern double evaluation_solution(int* sol);
extern double build_nearest_neighbour();
extern void build_solution();
extern void little_algorithm(double d0[NBR_TOWNS][NBR_TOWNS], int iteration, double eval_node_parent);

#endif // LITTLE_H
