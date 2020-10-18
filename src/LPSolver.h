/**
 * @file LPSolver.h
 * @brief 线性规划问题求解，实现方案参考：https://brilliant.org/wiki/linear-programming/
 * @author zhangwuxu
 * @version 1.0
 * @date 2020-10-17
 */

#pragma once

#include <tuple>
#include <vector>
#include <Eigen/Dense>

/**
 * @brief 一个很大的数，需保证目标函数中所有的系数均小于这个数
 */
#define BIG_M 10000.0

class LPSolver {
public:
    LPSolver();
    ~LPSolver();

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
     */
    void init(const Eigen::VectorXd& object,
            const Eigen::MatrixXd& constraint,
            const Eigen::VectorXd& maximum);

    /**
     * @brief 求解
     *
     * @return has_solution, maximize, variables
     */
    std::tuple<bool, float, std::vector<float>> solve();

protected:
    /**
     * @brief 找到第一个可行解
     */
    void to_feasible_region();

    /**
     * @brief 是否是最优解
     *
     * @return 
     */
    bool optimal();

    /**
     * @brief 当前的解是否是有效解
     *
     * @return 
     */
    bool feasible_solution();

    /**
     * @brief 将(row, col)变为基础变量
     *
     * @param row
     * @param col
     */
    void pivot(int row, int col);

    /**
     * @brief 找到下一个基础变量
     *
     * @return row, col, has_solution
     */
    std::tuple<int, int, bool> select_basic_variable();

    /**
     * @brief 判断是否是基础变量
     *
     * @param col 所在的列
     *
     * @return is_basic, row
     */
    std::tuple<bool, int> is_basic_variable(int col);

    /**
     * @brief 找到目标函数中最小的系数
     *
     * @return coefficient, col
     */
    std::tuple<float, int> min_object_coeff();

protected:
    Eigen::MatrixXd equations_;
    int var_num_;
    int slack_var_num_;
    int artificial_var_num_;
};
