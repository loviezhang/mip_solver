#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include "../src/LPSolver.h"

using ::testing::ElementsAre;

// max  7x + 5y
// s.t. 2x + 3y <= 90
//      3x + 2y <= 120
// solution
//      x = 36
//      y = 6
//      maximum = 282
TEST(LPSolver, Test1) {
    Eigen::VectorXd object;
    Eigen::MatrixXd constraint;
    Eigen::VectorXd maximum;

    object.resize(2);
    object << 7, 5;

    constraint.resize(2, 2);
    constraint << 2, 3,
                  3, 2;

    maximum.resize(2);
    maximum << 90, 120;

    LPSolver solver;
    solver.init(object, constraint, maximum);


    auto [has_solution, maximize, variables] = solver.solve();

    ASSERT_TRUE(has_solution);
    ASSERT_EQ(maximize, 282);
    ASSERT_THAT(variables, ElementsAre(36, 6));
}

// max  7x + 5y
// s.t. 2x + 3y >= 90
//      3x + 2y >= 120
// solution
//      no solution
TEST(LPSolver, Test2) {
    Eigen::VectorXd object;
    Eigen::MatrixXd constraint;
    Eigen::VectorXd maximum;

    object.resize(2);
    object << 7, 5;

    constraint.resize(2, 2);
    constraint << -2, -3,
                  -3, -2;

    maximum.resize(2);
    maximum << -90, -120;

    LPSolver solver;
    solver.init(object, constraint, maximum);


    auto [has_solution, dummy1, dummy2] = solver.solve();

    ASSERT_FALSE(has_solution);
}

// min   x - 10y
// s.t. -x +  5y <= 25
//      6x +  5y <= 60
//       x +   y >= 2
// solution
//      x = 5
//      y = 6
//      minimize = -55
TEST(LPSolver, Test3) {
    Eigen::VectorXd object;
    Eigen::MatrixXd constraint;
    Eigen::VectorXd maximum;

    // 乘-1，转换为求最大值
    object.resize(2);
    object << -1, 10;

    constraint.resize(3, 2);
    constraint << -1,  5,
                   6,  5,
                  -1, -1;  // x + y >= 2 转为 -x - y <= -2

    maximum.resize(3);
    maximum << 25, 60, -2;

    LPSolver solver;
    solver.init(object, constraint, maximum);


    auto [has_solution, maximize, variables] = solver.solve();

    ASSERT_TRUE(has_solution);
    ASSERT_EQ(-1 * maximize, -55);
    ASSERT_THAT(variables, ElementsAre(5, 6));
}

// min  12x + 9y + 12z
// s.t.   x      +   z <= 40
//             y       <= 40
//                   z >= 20
//       5x + 4y - 10z  = 0
// solution
//      x = 8
//      y = 40
//      z = 20
//      minimize = 696
TEST(LPSolver, Test4) {
    Eigen::VectorXd object;
    Eigen::MatrixXd constraint;
    Eigen::VectorXd maximum;

    // 乘-1，转换为求最大值
    object.resize(3);
    object << -12, -9, -12;

    constraint.resize(5, 3);
    constraint <<  1,   0,   1,
                   0,   1,   0,
                   0,   0,  -1, // 变为 -z <= -20
                   5,   4, -10, // 变为5x + 4y - 10z <= 0 && -1 * (5x + 4y - 10z) >= 0
                  -5,  -4,  10;

    maximum.resize(5);
    maximum << 40, 40, -20, 0, 0;

    LPSolver solver;
    solver.init(object, constraint, maximum);


    auto [has_solution, maximize, variables] = solver.solve();

    ASSERT_TRUE(has_solution);
    ASSERT_EQ(-1 * maximize, 696);
    ASSERT_THAT(variables, ElementsAre(8, 40, 20));
}
