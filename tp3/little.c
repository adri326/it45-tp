/**
 * Projec : gtsp (voyageur de commerce)
 *
 * Date   : 07/04/2014
 * Author : Olivier Grunder
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdbool.h>
#include "little.h"


/* Distance matrix */
double dist[NBR_TOWNS][NBR_TOWNS] ;

/* Each edge has a starting and ending node */
int starting_town[NBR_TOWNS] ;
int ending_town[NBR_TOWNS] ;

/* no comment */
int best_solution[NBR_TOWNS] ;
double best_eval=-1.0 ;


// TODO: replace `float[2]` by a bitcast-able struct
/**
 * Berlin52 :
 *  6 towns : Best solution (2315.15): 0 1 2 3 5 4
 * 10 towns : Best solution (2826.50): 0 1 6 2 7 8 9 3 5 4
 */
float coord[NBR_TOWNS][2]=
{
    {565.0,  575.0},
    { 25.0,  185.0},
    {345.0,  750.0},
    {945.0,  685.0},
    {845.0,  655.0},
    {880.0,  660.0}
} ;

double dist2(float a[2], float b[2]) {
    float dx = a[0] - b[0];
    float dy = a[1] - b[1];
    return sqrt(dx * dx + dy * dy);
}

/// Computes the `dist` matrix, returning a 2d array of dimensions (n_coords, n_coords)
double* compute_distance(float coords[][2], size_t n_coords) {
    double* res = malloc(sizeof(double) * n_coords * n_coords);
    for (size_t x = 0; x < n_coords; x++) {
        res[x*n_coords + x] = -1;
        for (size_t y = 0; y < x; y++) {
            res[y*n_coords + x] = dist2(coords[x], coords[y]);
            res[x*n_coords + y] = res[y*n_coords + x];
        }
    }

    return res;
}

/**
 * print a matrix
 */
void print_matrix(double d[NBR_TOWNS][NBR_TOWNS])
{
    int i, j ;
    for (i=0; i<NBR_TOWNS; i++)
    {
        printf ("%d:", i) ;
        for (j=0; j<NBR_TOWNS; j++)
        {
            printf ("%6.1f ", d[i][j]) ;
        }
        printf ("\n") ;
    }
}



/**
 * print a solution
 */
void print_solution(int* sol, double eval)
{
    int i ;
    printf ("(%.2f): ", eval) ;
    for (i=0; i<NBR_TOWNS; i++)
        printf ("%d ",sol[i]);
    printf("\n") ;
}




/**
 * evaluation of a solution
 */
double evaluate(int* sol)
{
    double eval=0 ;
    int i ;
    for (i=0; i<NBR_TOWNS-1; i++)
    {
        eval += dist[sol[i]][sol[i+1]] ;
    }
    eval += dist[sol[NBR_TOWNS-1]][sol[0]] ;

    return eval ;

}




/**
    @fn double build_nearest_neighbor(double* dist, size_t n_towns)
    Computes the path governed by the nearest neighbor heuristic; this is a quick way to get an upper bound of the solution.
    @arg dist - the distance matrix; is expected to be allocated such that `dist[y * n_towns + x]` is
        the distance from city `y` to city `x`.
    @arg n_towns - the number of towns; `dist` should have as size `n_towns * n_towns`.
    @return The distance travalled by the solution.
        In unit test mode, a struct containing the solution and the distance is returned; its `solution` field should be freed.
**/
#ifdef UNIT_TEST
    struct nn_t build_nearest_neighbor(double* dist, size_t n_towns) {
#else
    double build_nearest_neighbor(double* dist, size_t n_towns) {
#endif
    // Build the solution
    int* solution = malloc(sizeof(int) * n_towns);
    for (size_t i = 0; i < NBR_TOWNS; i++) solution[i] = -1;

    size_t current = 0;
    solution[0] = 0;

    double evaluation = 0.0;

    for (size_t n = 1; n < n_towns; n++) {
        int candidate = -1;
        double candidate_score = 0.0;
        for (size_t j = 1; j < n_towns; j++) {
            bool already_reached = false;
            if (j == current || dist[current * n_towns + j] < 0) continue; // Ignore invalid paths
            for (size_t i = 0; i < n; i++) {
                already_reached = already_reached || (solution[i] == (int)j);
            }
            if (already_reached) continue; // Ignore cities already explored

            if (dist[current * n_towns + j] < candidate_score || candidate == -1) {
                candidate = j;
                candidate_score = dist[current * n_towns + j];
            }
        }

        if (candidate == -1) {
            fprintf(stderr, "No city reachable from city %zu!\n", current);
            if (n > 1) fprintf(stderr, "Path so far:\n");
            for (size_t i = 1; i < n; i++) {
                fprintf(stderr, "%d -> %d\n", solution[i - 1], solution[i]);
            }

            #ifdef UNIT_TEST
                struct nn_t res = {
                    .solution = NULL,
                    .evaluation = -1.0
                };
                return res;
            #else
                return -1.0;
            #endif
        }

        current = candidate;
        solution[n] = candidate;
        evaluation += candidate_score;
    }

    evaluation += dist[current * n_towns + 0];

    // for (size_t i = 1; i < n_towns; i++) {
    //     printf("%d -> %d\n", solution[i - 1], solution[i]);
    // }

#ifdef UNIT_TEST
    struct nn_t res = {
        .solution = solution,
        .evaluation = evaluation
    };
    return res;
#else
    // Only return the evaluation of the solution
    free(solution);
    return evaluation;
#endif
}




/**
 *  Build final solution
 */
void build_solution()
{
    int i, solution[NBR_TOWNS] ;

    int indiceCour = 0;
    int villeCour = 0;

    while (indiceCour < NBR_TOWNS)
    {

        solution[indiceCour] = villeCour ;

        // Test si le cycle est hamiltonien
        for (i = 0; i < indiceCour; i++)
        {
            if (solution[i] == villeCour)
            {
                /* printf ("cycle non hamiltonien\n") ; */
                return;
            }
        }
        // Recherche de la ville suivante
        int trouve = 0;
        int i = 0;
        while ((!trouve) && (i < NBR_TOWNS))
        {
            if (starting_town[i] == villeCour)
            {
                trouve = 1;
                villeCour = ending_town[i];
            }
            i++;
        }
        indiceCour++;
    }

    double eval = evaluate(solution) ;

    if (best_eval<0 || eval < best_eval)
    {
        best_eval = eval ;
        for (i=0; i<NBR_TOWNS; i++)
            best_solution[i] = solution[i] ;
        printf ("New best solution: ") ;
        print_solution (solution, best_eval) ;
    }
    return;
}




/**
 *  Little Algorithm
 */
void little_algorithm(double d0[NBR_TOWNS][NBR_TOWNS], int iteration, double eval_node_parent)
{

    if (iteration == NBR_TOWNS)
    {
        build_solution ();
        return;
    }

    /* Do the modification on a copy of the distance matrix */
    double d[NBR_TOWNS][NBR_TOWNS] ;
    memcpy (d, d0, NBR_TOWNS*NBR_TOWNS*sizeof(double)) ;

    // int i, j ;

    double eval_node_child = eval_node_parent;

    /**
     * substract the min of the rows and the min of the columns
     * and update the evaluation of the current node
     *  TO COMPLETE
     *  ...
     *  ...
     */


    /* Cut : stop the exploration of this node */
    if (best_eval>=0 && eval_node_child >= best_eval)
        return;


    /**
     *  Compute the penalities to identify the zero with max penalty
     *  If no zero in the matrix, then return, solution infeasible
     *  TO COMPLETE
     *  ...
     *  ...
     */
    /* row and column of the zero with the max penalty */
    int izero=-1, jzero=-1 ;

    /**
     *  Store the row and column of the zero with max penalty in
     *  starting_town and ending_town
     *  TO COMPLETE
     *  ...
     *  ...
     */

    /* Do the modification on a copy of the distance matrix */
    double d2[NBR_TOWNS][NBR_TOWNS] ;
    memcpy (d2, d, NBR_TOWNS*NBR_TOWNS*sizeof(double)) ;

    /**
     *  Modify the matrix d2 according to the choice of the zero with the max penalty
     *  TO COMPLETE
     *  ...
     *  ...
     */

    /* Explore left child node according to given choice */
    little_algorithm(d2, iteration + 1, eval_node_child);

    /* Do the modification on a copy of the distance matrix */
    memcpy (d2, d, NBR_TOWNS*NBR_TOWNS*sizeof(double)) ;

    /**
     *  Modify the dist matrix to explore the other possibility : the non-choice
     *  of the zero with the max penalty
     *  TO COMPLETE
     *  ...
     *  ...
     */

    /* Explore right child node according to non-choice */
    little_algorithm(d2, iteration, eval_node_child);

}

// We provide an alternative main() if we're in unit testing mode
#ifndef UNIT_TEST

/**
 *
 */
int main (int argc, char* argv[])
{

    best_eval = -1 ;

    /* Print problem informations */
    printf ("Points coordinates:\n") ;
    int i ;
    for (i=0; i<NBR_TOWNS; i++)
    {
        printf ("Node %d: x=%f, y=%f\n", i, coord[i][0], coord[i][1]) ;
    }
    printf ("\n") ;


    /* Calcul de la matrice des distances */
    /**
     *  TO COMPLETE
     *  ...
     *  ...
     */

    double* dist_tmp = compute_distance(coord, NBR_TOWNS);
    memcpy(dist, dist_tmp, sizeof(double) * NBR_TOWNS * NBR_TOWNS);
    free(dist_tmp);

    printf("Distance Matrix:\n") ;
    print_matrix(dist) ;
    printf("\n") ;

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
     *  print_solution (best_solution, best_eval) ;
     */

    printf ("Hit RETURN!\n") ;
    getchar() ;

    return 0 ;
}

#endif // UNIT_TEST
