/**
 * @file MIPSolver.h
 * @brief 求解混合整数规划问题，方案参考：https://www.coursera.org/learn/discrete-optimization
 * @author zhangwuxu
 * @version 1.0
 * @date 2020-10-18
 */

#pragma once

#include <vector>
#include <tuple>
#include <Eigen/Dense>

class MIPSolver {
public:
    MIPSolver();
    ~MIPSolver();

    /**
     * @brief 初始化求解器
     *
     * @param object
     *   目标函数。只支持求最大值，若要求的最小值，可以通过乘-1改为求最大值。
     * @param constraint
     *   约束条件，constraint * x <= maximum。
     *   对于最小值约束，参考目标函数的处理方法，改为最大值约束。
     *   对于等号约束，如x = 4，可以改为 x <= 4 && -1 * x <= -4。
     *   注意，这里还有个隐含的约束条件是x>=0。
     * @param maximum
     *   同constraint。
     * @param bin_var_num
     *   0-1变量的数量，建模时要保证0-1变量都在前面。如x1...x5是0-1变量，x6...x7是非0-1变量
     */
    void init(const Eigen::VectorXd& object,
            const Eigen::MatrixXd& constraint,
            const Eigen::VectorXd& maximum,
            int bin_var_num);

    /**
     * @brief 求解
     *
     * @return has_solution, maximize, variables
     */
    std::tuple<bool, float, std::vector<float>> solve();

protected:

protected:
};
