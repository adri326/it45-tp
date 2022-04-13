/**
 * Projec : gtsp (voyageur de commerce)
 *
 * Date   : 07/04/2014
 * Author : Olivier Grunder
 */

#include "little.h"
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Distance matrix */
double* dist = NULL;

/* Each edge has a starting and ending node */
// TODO: malloc this
int starting_town[NBR_TOWNS];
int ending_town[NBR_TOWNS];

/* no comment */
int* best_solution = NULL;
double best_eval = -1.0;

// TODO: replace `float[2]` by a bitcast-able struct
/**
 * Berlin52 :
 *  6 towns : Best solution (2315.15): 0 1 2 3 5 4
 * 10 towns : Best solution (2826.50): 0 1 6 2 7 8 9 3 5 4
 */
float coord[NBR_TOWNS][2] = {
    {565.0, 575.0},
    {25.0, 185.0},
    {345.0, 750.0},
    {945.0, 685.0},
    {845.0, 655.0},
    {880.0, 660.0},
    {25.0, 230.0},
    {525.0, 1000.0},
    {580.0, 1175.0},
    {650.0, 1130.0}
};

double dist2(float a[2], float b[2]) {
    float dx = a[0] - b[0];
    float dy = a[1] - b[1];
    return sqrt(dx * dx + dy * dy);
}

/// Computes the `dist` matrix, returning a 2d array of dimensions (n_coords,
/// n_coords)
double* compute_distance(float coords[][2], size_t n_coords) {
    double* res = malloc(sizeof(double) * n_coords * n_coords);
    for (size_t x = 0; x < n_coords; x++) {
        res[x * n_coords + x] = -1;
        for (size_t y = 0; y < x; y++) {
            res[y * n_coords + x] = dist2(coords[x], coords[y]);
            res[x * n_coords + y] = res[y * n_coords + x];
        }
    }

    return res;
}

/**
 * print a matrix
 */
void print_matrix(double* d, size_t n_cities) {
    for (size_t y = 0; y < n_cities; y++) {
        printf("[%zu]:", y);
        for (size_t x = 0; x < n_cities; x++) {
            printf("%6.1f ", d[y * n_cities + x]);
        }
        printf("\n");
    }
}

/**
 * print a solution
 */
void print_solution(int* sol, double eval, size_t n_cities) {
    printf("(%.2f): ", eval);
    for (size_t i = 0; i < n_cities; i++) {
        printf("%d ", sol[i]);
    }
    printf("\n");
}

/**
 * evaluation of a solution
 */
double evaluate(int* sol, size_t n_cities) {
    double eval = 0;
    int i;
    for (i = 0; i < n_cities - 1; i++) {
        eval += dist[sol[i + 1] * n_cities + sol[i]];
    }
    eval += dist[sol[0] * n_cities + sol[n_cities - 1]];

    return eval;
}

/** @fn build_nearest_neighbor_sub(dist, n_cities)
    Function called by `build_nearest_neighbor`.
**/
struct nn_t build_nearest_neighbor_sub(double* dist, size_t n_cities) {
    // Build the solution
    int* solution = malloc(sizeof(int) * n_cities);
    for (size_t i = 0; i < NBR_TOWNS; i++)
        solution[i] = -1;

    size_t current = 0;
    solution[0] = 0;

    double evaluation = 0.0;

    for (size_t n = 1; n < n_cities; n++) {
        int candidate = -1;
        double candidate_score = 0.0;
        for (size_t j = 1; j < n_cities; j++) {
            bool already_reached = false;
            if (j == current || dist[current * n_cities + j] < 0)
                continue; // Ignore invalid paths
            for (size_t i = 0; i < n; i++) {
                already_reached = already_reached || (solution[i] == (int)j);
            }
            if (already_reached)
                continue; // Ignore cities already explored

            if (dist[current * n_cities + j] < candidate_score ||
                candidate == -1) {
                candidate = j;
                candidate_score = dist[current * n_cities + j];
            }
        }

        if (candidate == -1) {
            fprintf(stderr, "No city reachable from city %zu!\n", current);
            if (n > 1)
                fprintf(stderr, "Path so far:\n");
            for (size_t i = 1; i < n; i++) {
                fprintf(stderr, "%d -> %d\n", solution[i - 1], solution[i]);
            }

            struct nn_t res = {.solution = NULL, .evaluation = -1.0};
            return res;
        }

        current = candidate;
        solution[n] = candidate;
        evaluation += candidate_score;
    }

    evaluation += dist[current * n_cities + 0];

    // for (size_t i = 1; i < n_cities; i++) {
    //     printf("%d -> %d\n", solution[i - 1], solution[i]);
    // }

    struct nn_t res = {.solution = solution, .evaluation = evaluation};
    return res;
}

/**
    @fn double build_nearest_neighbor(double* dist, size_t n_cities)
    Computes the path governed by the nearest neighbor heuristic; this is a
quick way to get an upper bound of the solution.
    @arg dist - the distance matrix; is expected to be allocated such that
`dist[y * n_cities + x]` is the distance from city `y` to city `x`.
    @arg n_cities - the number of towns; `dist` should have as size `n_cities *
n_cities`.
    @return The total distance travelled in the solution.
**/
double build_nearest_neighbor(double* dist, size_t n_cities) {
    struct nn_t res = build_nearest_neighbor_sub(dist, n_cities);

    free(res.solution);

    return res.evaluation;
}

/**
 *  Build final solution
 */
void build_solution(int* starting_town, int* ending_town, size_t n_cities) {
    int* solution = malloc(sizeof(int) * n_cities);
    int current = 0;

    for (size_t index = 0; index < NBR_TOWNS; index++) {
        solution[index] = current;

        // Test si le cycle est hamiltonien (n*O(n))
        for (size_t i = 0; i < index; i++) {
            if (solution[i] == current) {
                // fprintf(stderr, "Non-hamiltonian cycle!\n");
                // for (size_t n = 0; n < n_cities; n++) {
                //     fprintf(stderr, "%zu: (%d -> %d)\n", n, starting_town[n], ending_town[n]);
                // }
                free(solution);
                return;
            }
        }

        // Recherche de la ville suivante
        bool found = false;
        for (size_t i = 0; i < n_cities; i++) {
            if (starting_town[i] == current) {
                current = ending_town[i];
                found = true;
                break;
            }
        }

        if (!found) {
            fprintf(stderr, "No path starting at %d!\n", current);
            for (size_t n = 0; n < n_cities; n++) {
                fprintf(stderr, "%zu: (%d -> %d)\n", n, starting_town[n], ending_town[n]);
            }
            free(solution);
            return;
        }
    }

    double eval = evaluate(solution, n_cities);

    if (best_eval < 0 || eval < best_eval) {
        best_eval = eval;
        if (best_solution != NULL) free(best_solution);
        best_solution = solution;

        printf("New best solution: ");
        print_solution(solution, best_eval, n_cities);
    } else {
        free(solution);
    }
}

/**
 *  Little Algorithm
 */
void little_algorithm(
    double* dist,
    int iteration,
    double eval_node_parent,
    size_t n_cities
) {
    if (iteration == NBR_TOWNS) {
        build_solution(starting_town, ending_town, NBR_TOWNS);
        return;
    }

    // We do the modifications directly on dist;
    // the motive for this is that the parent call will have to have created a copy of dist anyways

    double eval = eval_node_parent; // TODO: remove and operate on `eval_node_parent` directly

    /**
     * substract the min of the rows and the min of the columns
     * and update the evaluation of the current node
     */

    // Zero-out columns
    for (size_t x = 0; x < n_cities; x++) {
        double min = INFINITY;
        for (size_t y = 0; y < n_cities; y++) {
            double current = dist[y * n_cities + x];
            if (current >= 0 && current < min) {
                min = current;
            }
        }
        if (isinf(min)) continue;

        eval += min;
        for (size_t y = 0; y < n_cities; y++) {
            // Should be optimizable with cmov or the f registers
            if (dist[y * n_cities + x] >= 0) dist[y * n_cities + x] -= min;
        }
    }

    // Zero-out rows
    for (size_t y = 0; y < n_cities; y++) {
        double min = INFINITY;
        for (size_t x = 0; x < n_cities; x++) {
            double current = dist[y * n_cities + x];
            if (current >= 0 && current < min) {
                min = current;
            }
        }
        if (isinf(min)) continue;

        eval += min;
        for (size_t x = 0; x < n_cities; x++) {
            // Should be optimizable with cmov or the f registers
            if (dist[y * n_cities + x] >= 0) dist[y * n_cities + x] -= min;
        }
    }

    /* Cut : stop the exploration of this node */
    if (best_eval >= 0 && eval >= best_eval) {
        #ifdef VERBOSE
        printf("[%d, %.2lf, %.2lf] Cut!\n", iteration, eval, best_eval);
        #endif
        return;
    }

    /**
     *  Compute the penalities to identify the zero with max penalty
     *  If no zero in the matrix, then return, solution infeasible
     */

    // row and column of the zero with the max penalty
    double best_penalty = -1.0;
    size_t x_zero = 0, y_zero = 0;

    for (size_t y = 0; y < n_cities; y++) {
        for (size_t x = 0; x < n_cities; x++) {
            if (x == y || dist[y * n_cities + x] != 0.0) continue;
            double penalty = 0.0;
            // Compute Σ_{n≠i} dist[n, y] and Σ_{n≠j} dist[x, n] at once
            for (size_t n = 0; n < n_cities; n++) {
                if (x != n && dist[y * n_cities + n] > 0) penalty += dist[y * n_cities + n];
                if (y != n && dist[n * n_cities + x] > 0) penalty += dist[n * n_cities + x];
            }

            if (penalty > best_penalty) {
                best_penalty = penalty;
                x_zero = x;
                y_zero = y;
            }
        }
    }


    // No zero, thus no feasible solution
    if (best_penalty == -1.0) {
        return;
    }

    /**
     *  Store the row and column of the zero with max penalty in
     *  starting_town and ending_town
     */
    starting_town[iteration] = x_zero;
    ending_town[iteration] = y_zero;

    // Do the modification on a copy of the distance matrix
    double* dist2 = malloc(sizeof(double) * n_cities * n_cities);
    memcpy(dist2, dist, sizeof(double) * n_cities * n_cities);

    /**
     *  Modify the matrix d2 according to the choice of the zero with the max
     * penalty
     */

    for (size_t n = 0; n < n_cities; n++) {
        dist2[n * n_cities + x_zero] = -1;
        dist2[y_zero * n_cities + n] = -1;
    }

    // Prevent n=2 cycles
    dist2[x_zero * n_cities + y_zero] = -1;

    #ifdef VERBOSE
    printf("[%d, %.2lf, %.2lf] Exploring branch %zu→%zu\n", iteration, eval, best_eval, x_zero, y_zero);
    #endif
    /* Explore left child node according to given choice */
    little_algorithm(dist2, iteration + 1, eval, n_cities);

    free(dist2);

    // The right branch is computed on dist, so we don't have to re-allocate more memory

    /**
     *  Modify the dist matrix to explore the other possibility : the non-choice
     *  of the zero with the max penalty
     */

    dist[y_zero * n_cities + x_zero] = -1;

    #ifdef VERBOSE
    printf("[%d, %.2lf, %.2lf] Exploring branch ¬%zu→%zu\n", iteration, eval, best_eval, x_zero, y_zero);
    #endif
    /* Explore right child node according to non-choice */
    little_algorithm(dist, iteration, eval, n_cities);
}

// We provide an alternative main() if we're in unit testing mode
#ifndef UNIT_TEST

/**
 *
 */
int main(int argc, char* argv[]) {

    best_eval = -1;

    /* Print problem informations */
    printf("Points coordinates:\n");
    int i;
    for (i = 0; i < NBR_TOWNS; i++) {
        printf("Node %d: x=%f, y=%f\n", i, coord[i][0], coord[i][1]);
    }
    printf("\n");

    /* Calcul de la matrice des distances */
    /**
     *  TO COMPLETE
     *  ...
     *  ...
     */

    dist = compute_distance(coord, NBR_TOWNS);
    double* dist2 = malloc(sizeof(double) * NBR_TOWNS * NBR_TOWNS);
    memcpy(dist2, dist, sizeof(double) * NBR_TOWNS * NBR_TOWNS);

    printf("Distance Matrix:\n");
    print_matrix(dist, NBR_TOWNS);
    printf("\n");

    double nearest_neighbour = build_nearest_neighbor((double*)dist, NBR_TOWNS);

    printf("%lf\n", nearest_neighbour);

    /** Little : uncomment when needed
     *
     *  int iteration = 0 ;
     *  double lowerbound = 0.0 ;
     *
     *  little_algorithm(dist, iteration, lowerbound) ;
     *
     *  printf("Best solution:") ;
     *  print_solution (best_solution, best_eval, NBR_TOWNS) ;
     */

    // TODO: allocate starting_town and ending_town; move little_algorithm to a subfunction?

    best_eval = nearest_neighbour + 0.001;
    little_algorithm(dist2, 0, 0.0, NBR_TOWNS);

    if (best_solution != NULL) {
        printf("Best solution:\n");
        print_solution(best_solution, best_eval, NBR_TOWNS);

        free(best_solution);
    }

    free(dist);
    free(dist2);

    // printf("Hit RETURN!\n");
    // getchar();

    return 0;
}

#endif // UNIT_TEST
