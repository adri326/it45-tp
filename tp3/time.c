#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "little.h"

int main(int argc, char* argv[]) {
    const int MAX_TIME = argc > 0 ? atoi(argv[1]) : 60;
    printf("\"n\", \"Time (s)\"\n");

    for (size_t n_cities = 3; n_cities <= MAX_CITIES; n_cities++) {
        double* dist = compute_distance(berlin52_coords, n_cities);

        int n = 0;
        double sum = 0.0;

        for (int i = 0; i < 10; i++) {
            clock_t start = clock();
            solution_t solution = little_algorithm(dist, n_cities, true, 0);
            clock_t diff = clock() - start;

            n++;
            sum += (double)diff / (double)CLOCKS_PER_SEC;
            free_solution(solution);

            if (sum > MAX_TIME) {
                break;
            }
        }

        printf("%zu, %f\n", n_cities, sum / (double)n);
        fflush(stdout); // for tee to write to stdout

        free(dist);
    }
}
