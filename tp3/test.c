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
    int solution[NBR_TOWNS];
    double expected = 0.0;
    for (size_t n = 0; n < NBR_TOWNS; n++) {
        solution[n] = n;
        if (n < NBR_TOWNS - 1)
            expected += dist[(n + 1) * NBR_TOWNS + n];
    }
    ck_assert_double_eq_tol(evaluate(solution, NBR_TOWNS), expected, 0.0001);

    // Invert order
    for (size_t n = 0; n < NBR_TOWNS; n++) {
        solution[n] = NBR_TOWNS - n - 1;
    }

    ck_assert_double_eq_tol(evaluate(solution, NBR_TOWNS), expected, 0.0001);
}
END_TEST

void check_nn_solution(
    struct nn_t* solution,
    int* expected,
    double* dist,
    size_t n_cities
) {
    ck_assert_int_eq(solution->solution[0], 0);

    if (n_cities > 1 && solution->solution[1] == expected[1]) {
        // prograde traversal
        for (size_t n = 1; n < n_cities; n++) {
            ck_assert_int_eq(solution->solution[n], expected[n]);
        }
    } else {
        // retrograde traversal
        for (size_t n = 1; n < n_cities; n++) {
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

    struct nn_t solution = build_nearest_neighbor_sub(dist, 4);

    check_nn_solution(&solution, expected, dist, 4);

    free(solution.solution);
    free(dist);
}
END_TEST

/** @test test_nearest_neighbor_random
    Robust testing of the `build_nearest_neighbor_sub` function, using a set of
`n` cities placed around a circle. The order of the cities is shuffled between
two rounds, alongside the radius and phase of the regular polygon. Two solutions
are possible for each iteration, given the symmetries of the polygon.
`check_nn_solution` takes care of this.
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
            struct nn_t solution = build_nearest_neighbor_sub(dist, n_cities);

            // Compare with expected solution
            check_nn_solution(&solution, expected, dist, n_cities);

            free(solution.solution);
            free(dist);
            free(order);
            free(expected);
        }
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
