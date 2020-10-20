#pragma once

#include <gtest/gtest.h>

#define ABS_ERROR 0.0001

inline void match_result(double maximize,
        double maximize_expect,
        std::vector<double> variables,
        std::vector<double> variables_expect) {
    EXPECT_NEAR(maximize, maximize_expect, ABS_ERROR);
    EXPECT_EQ(variables.size(), variables_expect.size());
    for (auto i = 0UL; i < variables.size(); i++) {
        EXPECT_NEAR(variables[i], variables_expect[i], ABS_ERROR);
    }
}
