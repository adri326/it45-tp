#include <stdlib.h>
#include <stdio.h>
#include "little.h"

int main(int argc, char* argv[]) {
    const int arg = argc > 1 ? atoi(argv[1]) : 10;
    if (arg <= 0 || arg > MAX_CITIES) {
        fprintf(stderr, "Invalid number of cities given as parameter: expected a number between 1 and %d, got %d\n", MAX_CITIES, arg);
    }

    int better_approx = argc > 2 ? atoi(argv[2]) : 0;
    if (better_approx < 0) better_approx = 0;

    const size_t n_cities = (size_t)arg;

    /* Print problem informations */
    printf("Points coordinates:\n");
    for (size_t i = 0; i < n_cities; i++) {
        printf("Node %zu: x=%f, y=%f\n", i, berlin52_coords[i][0], berlin52_coords[i][1]);
    }
    printf("\n");

    // Calcul de la matrice des distances
    double* dist = compute_distance(berlin52_coords, n_cities);

    printf("Distance Matrix:\n");
    print_matrix(dist, n_cities);
    printf("\n");

    double nearest_neighbour = build_nearest_neighbor(dist, n_cities);
    printf("%lf\n", nearest_neighbour);

    solution_t solution = little_algorithm(dist, n_cities, true, better_approx);

    printf("Best solution:\n");
    print_solution(solution, n_cities);

    free_solution(solution);
    free(dist);

    return 0;
}
