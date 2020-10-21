#include <gtest/gtest.h>
#include "../src/MIPSolver.h"
#include "utils.h"

#if 0

// 3 11
// 8 4   v1 w1
// 10 5  v2 w2
// 11 6  v3 w3
// max x1 * v1 + x2 * v2 + x3 * v3
// s.t. x1 * w1 + x2 * w2 + x3 * w3 <= 11
TEST(MIPSolver, Test1) {
    Eigen::VectorXd object;
    Eigen::MatrixXd constraint;
    Eigen::VectorXd maximum;

    object.resize(3);
    object << 8, 10, 11;

    constraint.resize(1, 3);
    constraint << 4, 4, 6;

    maximum.resize(1);
    maximum << 11;

    MIPSolver solver;
    solver.init(object, constraint, maximum, 3);

    auto [has_solution, maximize, variables] = solver.solve();

    EXPECT_TRUE(has_solution);
    match_result(maximize, 21, variables, {0, 1, 1});
}
#endif
