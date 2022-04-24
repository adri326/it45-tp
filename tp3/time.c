#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "little.h"

int main(int argc, char* argv[]) {
    for (size_t n_cities = 3; n_cities <= MAX_CITIES; n_cities++) {
        double* dist = compute_distance(berlin52_coords, n_cities);
        clock_t start = clock();

        solution_t solution = little_algorithm(dist, n_cities, true, 0);

        clock_t diff = clock() - start;

        printf("%zu, %f\n", n_cities, (float)diff / (float)CLOCKS_PER_SEC);
        fflush(stdout); // for tee to write to stdout

        free(dist);
        free_solution(solution);
    }
}
