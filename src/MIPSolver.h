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
#include <memory>
#include <unordered_map>
#include <Eigen/Dense>

class Node : std::enable_shared_from_this<Node> {
public:
    Node(const Eigen::VectorXd& object,
            const Eigen::MatrixXd& constraint,
            const Eigen::VectorXd& maximum,
            const int& bin_var_num);

    ~Node();

    // 是否可行
    bool feasible() const {
        return feasible_;
    }

    // 是否是叶子节点
    bool is_leaf() const {
        return bin_var_num_ == 0;
    }

    // 是否已完成搜索
    bool has_explored() const {
        return left_explored_ && right_explored_;
    }

    // 预估值，对于叶子节点，预估值等于真实值
    double estimite() const {
        return estimite_ + value_;
    }

    const std::vector<double>& variables() const {
        return variables_;
    }

    int split_point() const {
        return split_point_;
    }

    bool on_left_side() const {
        return on_left_side_;
    }

    // 回溯到父节点，若当前已经是root，返回nullptr
    std::shared_ptr<Node> backtrack() const {
        return parent_;
    }

    // 继续探索，返回探索到的新节点，若当前已经是叶子节点，返回nullptr
    std::shared_ptr<Node> spawn();

protected:
    void init();

    // 找到分裂点
    // return go_left, col
    std::tuple<bool, int> calc_split_point() const;

protected:
    Eigen::VectorXd object_;
    Eigen::MatrixXd constraint_;
    Eigen::VectorXd maximum_;
    int bin_var_num_;

    double estimite_;
    double value_;
    std::vector<double> variables_;

    bool feasible_;

    std::shared_ptr<Node> parent_;

    int split_point_;
    bool left_explored_;
    bool right_explored_;

    // 是否是父节点的左子节点
    bool on_left_side_;
};


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
     *   对于0-1变量，不需要额外添加x <= 1或是x >= 0的约束
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
    std::tuple<bool, double, std::vector<double>> solve();

protected:
    void update_biggest_node(std::shared_ptr<Node>& cur);

    bool need_backtrack(std::shared_ptr<Node>& cur);

    std::tuple<bool, double, std::vector<double>> get_result() const;

protected:
    Eigen::VectorXd object_;
    Eigen::MatrixXd constraint_;
    Eigen::VectorXd maximum_;
    int bin_var_num_;

    std::shared_ptr<Node> biggest_;
};
