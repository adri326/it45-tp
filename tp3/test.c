#include <stdlib.h>
#include <stdio.h>
#include <check.h>
#include "little.h"

START_TEST(test_distance) {

}
END_TEST

START_TEST(test_evaluation) {
    int solution[NBR_TOWNS];
    double expected = 0.0;
    for (size_t n = 0; n < NBR_TOWNS; n++) {
        solution[n] = n;
        if (n < NBR_TOWNS - 1) expected += dist[n][n + 1];
    }
    ck_assert_double_eq_tol(evaluation_solution(solution), expected, 0.0001);

    // Invert order
    for (size_t n = 0; n < NBR_TOWNS; n++) {
        solution[n] = NBR_TOWNS - n - 1;
    }

    ck_assert_double_eq_tol(evaluation_solution(solution), expected, 0.0001);
}
END_TEST

Suite* little_suite() {
    Suite* s = suite_create("little");
    TCase* tc_core = tcase_create("core");

    tcase_add_test(tc_core, test_evaluation);
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
