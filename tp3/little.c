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
// double* dist = NULL;

/* Each edge has a starting and ending node */
// TODO: malloc this
int* starting_town = NULL;
int* ending_town = NULL;

/* no comment */
int* best_solution = NULL;
double best_eval = -1.0;

// TODO: replace `float[2]` by a bitcast-able struct
/**
 * Berlin52 :
 *  6 towns : Best solution (2315.15): 0 1 2 3 5 4
 * 10 towns : Best solution (2826.50): 0 1 6 2 7 8 9 3 5 4
 */
float berlin52_coords[][2] = {
    {565.0, 575.0},
    {25.0, 185.0},
    {345.0, 750.0},
    {945.0, 685.0},
    {845.0, 655.0},
    {880.0, 660.0},
    {25.0, 230.0},
    {525.0, 1000.0},
    {580.0, 1175.0},
    {650.0, 1130.0},
    {1605.0, 620.0},
    {1220.0, 580.0},
    {1465.0, 200.0},
    {1530.0, 5.0},
    {845.0, 680.0},
    {725.0, 370.0},
    {145.0, 665.0},
    {415.0, 635.0},
    {510.0, 875.0},
    {560.0, 365.0},
    {300.0, 465.0},
    {520.0, 585.0},
    {480.0, 415.0},
    {835.0, 625.0},
    {975.0, 580.0},
    {1215.0, 245.0},
    {1320.0, 315.0},
    {1250.0, 400.0},
    {660.0, 180.0},
    {410.0, 250.0},
    {420.0, 555.0},
    {575.0, 665.0},
    {1150.0, 1160.0},
    {700.0, 580.0},
    {685.0, 595.0},
    {685.0, 610.0},
    {770.0, 610.0},
    {795.0, 645.0},
    {720.0, 635.0},
    {760.0, 650.0},
    {475.0, 960.0},
    {95.0, 260.0},
    {875.0, 920.0},
    {700.0, 500.0},
    {555.0, 815.0},
    {830.0, 485.0},
    {1170.0, 65.0},
    {830.0, 610.0},
    {605.0, 625.0},
    {595.0, 360.0},
    {1340.0, 725.0},
    {1740.0, 245.0}
};

double dist2(float a[2], float b[2]) {
    float dx = a[0] - b[0];
    float dy = a[1] - b[1];
    return sqrt(dx * dx + dy * dy);
}

/// Computes the `dist` matrix, returning a 2d array of dimensions
/// `(n_coords, n_coords)`
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
void print_solution(solution_t solution, size_t n_cities) {
    if (solution.solution != NULL) {
        printf("(%.2f): ", solution.evaluation);
        for (size_t i = 0; i < n_cities; i++) {
            printf("%02d ", solution.solution[i]);
        }
    } else {
        printf("(No solution)");
    }
    printf("\n");
}

/**
 * Free a solution
 **/
void free_solution(solution_t solution) {
    if (solution.solution != NULL) {
        free(solution.solution);
    }
}

/**
 * evaluation of a solution
 */
double evaluate(const double* dist, int* sol, size_t n_cities) {
    double eval = 0;
    for (size_t i = 0; i < n_cities - 1; i++) {
        eval += dist[sol[i + 1] * n_cities + sol[i]];
    }
    eval += dist[sol[0] * n_cities + sol[n_cities - 1]];

    return eval;
}

/** @fn build_nearest_neighbor_sub(dist, n_cities)
    Function called by `build_nearest_neighbor`.
**/
solution_t build_nearest_neighbor_sub(const double* dist, size_t n_cities) {
    // Build the solution
    int* solution = malloc(sizeof(int) * n_cities);
    for (size_t i = 0; i < n_cities; i++)
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

            solution_t res = {.solution = NULL, .evaluation = -1.0};
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

    solution_t res = {.solution = solution, .evaluation = evaluation};
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
double build_nearest_neighbor(const double* dist, size_t n_cities) {
    solution_t res = build_nearest_neighbor_sub(dist, n_cities);

    free(res.solution);

    return res.evaluation;
}

/**
 * Check whether the full path doesn't contain any sub-cycle; it is of time complexity Θ(n²)
 **/
bool is_valid_full_path(const int* path, size_t n_cities) {
    for (size_t n = 0; n < n_cities; n++) {
        for (size_t i = 0; i < n; i++) {
            if (path[i] == path[n]) return false;
        }
    }

    return true;
}

bool is_valid_sub_solution(int* starting_town, int* ending_town, size_t iteration, size_t n_cities) {
    // PERF TODO: don't malloc here
    int* indirect = malloc(sizeof(int) * n_cities);
    bool* visited = malloc(sizeof(bool) * n_cities);

    for (size_t index = 0; index < n_cities; index++) indirect[index] = -1;

    for (size_t index = 0; index < iteration; index++) {
        indirect[starting_town[index]] = ending_town[index];
    }

    for (size_t base = 0; base < n_cities; base++) {
        size_t ptr = base;
        while (indirect[ptr] >= 0) {
            // Node already checked
            //@variant sum(visited)
            if (visited[ptr]) break;

            visited[ptr] = true;

            // as per the loop condition,
            //@check indirect[ptr] >= 0
            ptr = (size_t)indirect[ptr];

            // Cycle detected!
            if (ptr == base) return false;
        }
    }

    free(indirect);
    free(visited);

    return true;
}

/**
 *  Build final solution
 */
void build_solution(const double* dist, int* starting_town, int* ending_town, size_t n_cities) {
    // PERF TODO: don't malloc here, use static instead?
    int* solution = malloc(sizeof(int) * n_cities);
    int current = 0;

    for (size_t index = 0; index < n_cities; index++) {
        solution[index] = current;

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

    if (!is_valid_full_path(solution, n_cities)) {
        free(solution);
        return;
    }

    double eval = evaluate(dist, solution, n_cities);

    if (best_eval < 0 || eval < best_eval) {
        best_eval = eval;
        if (best_solution != NULL) free(best_solution);
        best_solution = solution;

        #ifdef VERBOSE
        solution_t best_solution = {
            .solution = solution,
            .evaluation = best_eval
        };
        printf("New best solution: ");
        print_solution(best_solution, n_cities);
        #endif // VERBOSE
    } else {
        free(solution);
    }
}

solution_t little_algorithm(
    const double* dist,
    size_t n_cities,
    bool little_plus
) {
    double* dist2 = malloc(sizeof(double) * n_cities * n_cities);
    memcpy(dist2, dist, sizeof(double) * n_cities * n_cities);

    double nearest_neighbour = build_nearest_neighbor((double*)dist2, n_cities);

    // best_eval = nextafter(nearest_neighbour, INFINITY); // nearest_neighbour += ε
    best_eval = nearest_neighbour + 0.1;

    starting_town = malloc(sizeof(int) * n_cities);
    ending_town = malloc(sizeof(int) * n_cities);

    little_algorithm_rec(dist, dist2, 0, 0.0, n_cities, little_plus);

    free(dist2);
    free(starting_town);
    free(ending_town);

    // Monomorphize solution: if solution[1] > solution[n_cities - 1], swap the order around, check that the evaluation is the same and return that
    if (best_solution && n_cities > 2 && best_solution[1] > best_solution[n_cities - 1]) {
        int* mono_solution = malloc(sizeof(int) * n_cities);
        mono_solution[0] = best_solution[0];
        for (size_t n = 0; n < (n_cities - 1) / 2; n++) {
            mono_solution[1 + n] = best_solution[n_cities - n - 1];
            mono_solution[n_cities - n - 1] = best_solution[1 + n];
        }

        if (n_cities % 2 == 0) {
            mono_solution[(n_cities - 1) / 2 + 1] = best_solution[(n_cities - 1) / 2 + 1];
        }

        double mono_eval = evaluate(dist, mono_solution, n_cities);
        if (mono_eval <= best_eval) {
            free(best_solution);
            best_solution = mono_solution;
            best_eval = mono_eval;
        }
    }

    solution_t res = {
        .solution = best_solution,
        .evaluation = best_solution != NULL ? best_eval : INFINITY
    };

    best_solution = NULL;

    return res;
}

/**
 *  Little Algorithm
 */
void little_algorithm_rec(
    const double* dist,
    double* current_dist,
    size_t iteration,
    double eval_node_parent,
    size_t n_cities,
    bool little_plus
) {
    if (iteration == n_cities) {
        build_solution(dist, starting_town, ending_town, n_cities);
        return;
    } else if (little_plus && iteration > 0) {
        if (!is_valid_sub_solution(starting_town, ending_town, iteration, n_cities)) {
            return;
        }
    }

    // We do the modifications directly on current_dist;
    // the motive for this is that the parent call will have to have created a copy of current_dist anyways

    double eval = eval_node_parent; // TODO: remove and operate on `eval_node_parent` directly

    /**
     * substract the min of the rows and the min of the columns
     * and update the evaluation of the current node
     */

    // Zero-out columns
    for (size_t x = 0; x < n_cities; x++) {
        double min = INFINITY;
        for (size_t y = 0; y < n_cities; y++) {
            double current = current_dist[y * n_cities + x];
            if (current >= 0 && current < min) {
                min = current;
            }
        }
        if (isinf(min)) continue;

        eval += min;
        for (size_t y = 0; y < n_cities; y++) {
            // Should be optimizable with cmov or the f registers
            if (current_dist[y * n_cities + x] >= 0) current_dist[y * n_cities + x] -= min;
        }
    }

    // Zero-out rows
    for (size_t y = 0; y < n_cities; y++) {
        double min = INFINITY;
        for (size_t x = 0; x < n_cities; x++) {
            double current = current_dist[y * n_cities + x];
            if (current >= 0 && current < min) {
                min = current;
            }
        }
        if (isinf(min)) continue;

        eval += min;
        for (size_t x = 0; x < n_cities; x++) {
            // Should be optimizable with cmov or the f registers
            if (current_dist[y * n_cities + x] >= 0) current_dist[y * n_cities + x] -= min;
        }
    }

    /* Cut : stop the exploration of this node */
    if (best_eval >= 0 && eval >= best_eval) {
        #ifdef VERBOSE
        printf("[%zu, %.2lf, %.2lf] Cut!\n", iteration, eval, best_eval);
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
            if (x == y || current_dist[y * n_cities + x] != 0.0) continue;
            double penalty = 0.0;
            // Compute Σ_{n≠i} current_dist[n, y] and Σ_{n≠j} current_dist[x, n] at once
            for (size_t n = 0; n < n_cities; n++) {
                if (x != n && current_dist[y * n_cities + n] > 0) penalty += current_dist[y * n_cities + n];
                if (y != n && current_dist[n * n_cities + x] > 0) penalty += current_dist[n * n_cities + x];
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
    double* current_dist2 = malloc(sizeof(double) * n_cities * n_cities);
    memcpy(current_dist2, current_dist, sizeof(double) * n_cities * n_cities);

    /**
     *  Modify the matrix d2 according to the choice of the zero with the max
     * penalty
     */

    for (size_t n = 0; n < n_cities; n++) {
        current_dist2[n * n_cities + x_zero] = -1;
        current_dist2[y_zero * n_cities + n] = -1;
    }

    // Prevent n=2 cycles
    current_dist2[x_zero * n_cities + y_zero] = -1;

    #ifdef VERBOSE
    printf("[%zu, %.2lf, %.2lf] Exploring branch %zu→%zu\n", iteration, eval, best_eval, x_zero, y_zero);
    #endif
    /* Explore left child node according to given choice */
    little_algorithm_rec(dist, current_dist2, iteration + 1, eval, n_cities, little_plus);

    free(current_dist2);

    // The right branch is computed on dist, so we don't have to re-allocate more memory

    /**
     *  Modify the dist matrix to explore the other possibility : the non-choice
     *  of the zero with the max penalty
     */

    current_dist[y_zero * n_cities + x_zero] = -1;

    #ifdef VERBOSE
    printf("[%zu, %.2lf, %.2lf] Exploring branch ¬%zu→%zu\n", iteration, eval, best_eval, x_zero, y_zero);
    #endif
    /* Explore right child node according to non-choice */
    little_algorithm_rec(dist, current_dist, iteration, eval, n_cities, little_plus);
}
