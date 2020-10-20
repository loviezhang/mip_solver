#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include "../src/MIPSolver.h"

using ::testing::ElementsAre;

// max  7x + 5y
// s.t. 2x + 3y <= 90
//      3x + 2y <= 120
// solution
//      x = 36
//      y = 6
//      maximum = 282
TEST(MIPSolver, Test1) {
    ASSERT_TRUE(true);
}
