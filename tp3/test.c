#include "little.h"
#include <check.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

START_TEST(test_distance) {
    float coords[2][2] = {{565.0, 575.0}, {25.0, 185.0}};

    double* distances = compute_distance(coords, 2);

    ck_assert_double_eq_tol(distances[0 * 2 + 0], -1.0, 0.0001);
    ck_assert_double_eq_tol(distances[1 * 2 + 1], -1.0, 0.0001);
    ck_assert_double_eq_tol(distances[1 * 2 + 0], 666.10809, 0.0001);
    ck_assert_double_eq_tol(distances[0 * 2 + 1], 666.10809, 0.0001);

    free(distances);

    // Check permutation-independance

    float coords2[2][2] = {{25.0, 185.0}, {565.0, 575.0}};

    distances = compute_distance(coords2, 2);

    ck_assert_double_eq_tol(distances[0 * 2 + 0], -1.0, 0.0001);
    ck_assert_double_eq_tol(distances[1 * 2 + 1], -1.0, 0.0001);
    ck_assert_double_eq_tol(distances[1 * 2 + 0], 666.10809, 0.0001);
    ck_assert_double_eq_tol(distances[0 * 2 + 1], 666.10809, 0.0001);

    free(distances);
}
END_TEST

START_TEST(test_evaluation) {
    const size_t n_cities = 8;
    int* solution = malloc(sizeof(int) * n_cities);
    double expected = 0.0;
    double* dist = malloc(sizeof(double) * n_cities * n_cities);

    for (size_t x = 0; x < n_cities; x++) {
        for (size_t y = 0; y <= x; y++) {
            dist[y * n_cities + x] = dist[x * n_cities + y] = 1000.0 * rand() / (double)RAND_MAX;
        }
    }

    for (size_t n = 0; n < n_cities; n++) {
        solution[n] = n;
        if (n < n_cities - 1)
            expected += dist[(n + 1) * n_cities + n];
    }
    expected += dist[0 + solution[n_cities - 1]];

    ck_assert_double_eq_tol(evaluate(dist, solution, n_cities), expected, 0.0001);

    // Invert order
    for (size_t n = 1; n < n_cities; n++) {
        solution[n] = n_cities - n;
    }

    ck_assert_double_eq_tol(evaluate(dist, solution, n_cities), expected, 0.0001);

    free(solution);
    free(dist);
}
END_TEST

void check_solution(
    solution_t* solution,
    int* expected,
    double* dist,
    size_t n_cities
) {
    ck_assert_int_eq(solution->solution[0], 0);

    if (n_cities > 1 && solution->solution[1] == expected[1]) {
        // prograde traversal
        for (size_t n = 1; n < n_cities; n++) {
            if (solution->solution[n] != expected[n]) {
                fprintf(stderr, "Solution differs: expected is ");
                for (size_t i = 0; i < n_cities; i++) {
                    fprintf(stderr, "%02d ", expected[i]);
                }
                fprintf(stderr, "instead got ");
                for (size_t i = 0; i < n_cities; i++) {
                    fprintf(stderr, "%02d ", solution->solution[i]);
                }
                fprintf(stderr, "\n");
            }
            ck_assert_int_eq(solution->solution[n], expected[n]);
        }
    } else {
        // retrograde traversal
        for (size_t n = 1; n < n_cities; n++) {
            if (solution->solution[n] != expected[n_cities - n]) {
                fprintf(stderr, "Solution differs: expected is ");
                for (size_t i = 0; i < n_cities; i++) {
                    fprintf(stderr, "%02d ", expected[i]);
                }
                fprintf(stderr, "instead got ");
                for (size_t i = 0; i < n_cities; i++) {
                    fprintf(stderr, "%02d ", solution->solution[i]);
                }
                fprintf(stderr, "\n");
            }
            ck_assert_int_eq(solution->solution[n], expected[n_cities - n]);
        }
    }

    // Compute and check distance
    double sum = 0;
    for (size_t n = 1; n < n_cities; n++) {
        int prev = expected[n - 1];
        int next = expected[n];
        sum += dist[prev * n_cities + next];
    }
    sum += dist[expected[n_cities - 1] * n_cities + 0];

    ck_assert_double_eq_tol(
        solution->evaluation, sum, 0.000001 * sqrt(n_cities));
}

START_TEST(test_nearest_neighbor_simple) {
    float coords[4][2] = {{1.0, 0.0}, {1.0, 1.0}, {0.0, 0.0}, {0.0, 1.0}};
    double* dist = compute_distance(coords, 4);
    int expected[] = {0, 1, 3, 2};

    solution_t solution = build_nearest_neighbor_sub(dist, 4);

    check_solution(&solution, expected, dist, 4);

    free(solution.solution);
    free(dist);
}
END_TEST

/** @test test_nearest_neighbor_random
    Robust testing of the `build_nearest_neighbor_sub` function, using a set of
`n` cities placed around a circle. The order of the cities is shuffled between
two rounds, alongside the radius and phase of the regular polygon. Two solutions
are possible for each iteration, given the symmetries of the polygon.
`check_solution` takes care of this.
**/
START_TEST(test_nearest_neighbor_random) {
    srand(time(0));
    for (size_t n_cities = 5; n_cities < 20; n_cities++) {
        for (size_t sample = 0; sample < n_cities * n_cities / 4; sample++) {
            int* order = malloc(sizeof(int) * n_cities);
            // Fill with 0..n
            for (size_t n = 0; n < n_cities; n++)
                order[n] = (int)n;

            // Shuffle order
            for (size_t n = 1; n < n_cities - 1; n++) {
                size_t offset = rand() % (n_cities - 1 - n);
                int tmp = order[n];
                order[n] = order[n + offset];
                order[n + offset] = tmp;
            }

            float(*coords)[2] = malloc(sizeof(float) * 2 * n_cities);
            // Fill based on order
            float offset = ((float)rand() / (float)RAND_MAX) * M_PI * 2.0;
            float radius = ((float)rand() / (float)RAND_MAX) * 6.0 + 1.0;
            for (size_t n = 0; n < n_cities; n++) {
                coords[n][0] = radius *
                    cos(M_PI * 2.0 * (float)order[n] / (float)n_cities +
                        offset);
                coords[n][1] = radius *
                    sin(M_PI * 2.0 * (float)order[n] / (float)n_cities +
                        offset);
            }

            // Construct expected = order^{-1}
            int* expected = malloc(sizeof(int) * n_cities);
            for (size_t n = 0; n < n_cities; n++) {
                expected[order[n]] = n;
            }

            // Compute distance
            double* dist = compute_distance(coords, n_cities);

            // Build nearest neighbor solution
            solution_t solution = build_nearest_neighbor_sub(dist, n_cities);

            // Compare with expected solution
            check_solution(&solution, expected, dist, n_cities);

            free(solution.solution);
            free(dist);
            free(order);
            free(expected);
            free(coords);
        }
    }
}
END_TEST

// Tests the solution of little's algorithm against known results
START_TEST(test_little_berlin) {
    double* dist = compute_distance(berlin52_coords, 6);

    solution_t sol = little_algorithm(dist, 6, false, 0);
    int expected_6[6] = {0, 1, 2, 3, 5, 4};

    check_solution(&sol, expected_6, dist, 6);

    free_solution(sol);
    free(dist);

    dist = compute_distance(berlin52_coords, 10);

    sol = little_algorithm(dist, 10, false, 0);
    int expected_10[10] = {0, 1, 6, 2, 7, 8, 9, 3, 5, 4};

    check_solution(&sol, expected_10, dist, 10);

    free_solution(sol);
    free(dist);
}
END_TEST

// Verifies that little's algorithm gives a result that is always less than or equal to the nearest neighbor approximation
START_TEST(test_little_better) {
    int seed = time(0);
    srand(seed);
    printf("Seed: %d\n", seed);

    for (size_t i = 0; i < 1000; i++) {
        size_t n_cities = rand() % 9 + 4;

        double* dist = malloc(sizeof(double) * n_cities * n_cities);

        for (size_t x = 0; x < n_cities; x++) {
            dist[x * n_cities + x] = -1.0;

            for (size_t y = 0; y < x; y++) {
                dist[y * n_cities + x] = dist[x * n_cities + y] = 1000.0 * rand() / (double)RAND_MAX;
            }
        }

        // print_matrix(dist, n_cities);
        // printf("\n");

        solution_t nn_sol = build_nearest_neighbor_sub(dist, n_cities);
        solution_t little_sol = little_algorithm(dist, n_cities, false, 0);

        ck_assert_double_le_tol(little_sol.evaluation, nn_sol.evaluation, 0.0001);

        free_solution(nn_sol);
        free_solution(little_sol);
    }
}
END_TEST

// Verifies that is_valid_full_path correctly identifies solutions with cycles
START_TEST(test_is_valid_full_path) {
    int path_valid[] = {0, 1, 4, 3, 5, 2};
    ck_assert(is_valid_full_path(path_valid, 6));

    int path_invalid[] = {0, 1, 2, 1};
    ck_assert(!is_valid_full_path(path_invalid, 4));
}
END_TEST

// Verify that little+ gives the same result as little
START_TEST(test_little_plus) {
    int seed = time(0);
    srand(seed);
    printf("Seed: %d\n", seed);

    for (size_t i = 0; i < 1000; i++) {
        size_t n_cities = rand() % 9 + 4;

        double* dist = malloc(sizeof(double) * n_cities * n_cities);

        for (size_t x = 0; x < n_cities; x++) {
            dist[x * n_cities + x] = -1.0;

            for (size_t y = 0; y < x; y++) {
                dist[y * n_cities + x] = dist[x * n_cities + y] = 1000.0 * rand() / (double)RAND_MAX;
            }
        }

        // print_matrix(dist, n_cities);
        // printf("\n");

        solution_t nn_sol = build_nearest_neighbor_sub(dist, n_cities);
        solution_t little_sol = little_algorithm(dist, n_cities, false, 0);
        solution_t little_plus_sol = little_algorithm(dist, n_cities, true, 0);

        ck_assert_double_le_tol(little_sol.evaluation, nn_sol.evaluation, 0.0001);
        ck_assert_double_eq_tol(little_sol.evaluation, little_plus_sol.evaluation, 0.0001);

        check_solution(&little_plus_sol, little_sol.solution, dist, n_cities);

        free_solution(nn_sol);
        free_solution(little_sol);
        free_solution(little_plus_sol);
    }
}
END_TEST

Suite* little_suite() {
    Suite* s = suite_create("little");
    TCase* tc_core = tcase_create("core");

    tcase_add_test(tc_core, test_distance);
    tcase_add_test(tc_core, test_evaluation);

    tcase_add_test(tc_core, test_nearest_neighbor_simple);
    tcase_add_test(tc_core, test_nearest_neighbor_random);

    tcase_add_test(tc_core, test_little_berlin);
    tcase_add_test(tc_core, test_little_better);
    tcase_add_test(tc_core, test_is_valid_full_path);

    tcase_add_test(tc_core, test_little_plus);
    suite_add_tcase(s, tc_core);

    return s;
}

int main() {
    Suite* s = little_suite();
    SRunner* sr = srunner_create(s);

    srunner_run_all(sr, CK_VERBOSE);
    int number_failed = srunner_ntests_failed(sr);
    srunner_free(sr);

    return number_failed == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
