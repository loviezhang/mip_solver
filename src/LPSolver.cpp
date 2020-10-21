#include <iostream>
#include <cassert>
#include "LPSolver.h"
#include "utils.h"

LPSolver::LPSolver() {
}

LPSolver::~LPSolver() {
}

void LPSolver::init(const Eigen::VectorXd& object,
        const Eigen::MatrixXd& constraint,
        const Eigen::VectorXd& maximum) {
    if (object.rows() == 0) {
        equations_.resize(maximum.rows()+1, 1);
        equations_ << 0, maximum;
        return;
    }

    Eigen::MatrixXd full_constrant(constraint.rows(), constraint.cols()+1);
    full_constrant << constraint, maximum;
    std::cout << "object: " << object.transpose() << std::endl;
    std::cout << "constraint: " << std::endl << full_constrant << std::endl;

    /**
     * 加入松弛变量和人工变量后，表达式变为
     * (1) z - c * x = y (y初始为0)
     * (2) A * x = b
     * (3) x >= 0
     * 我们要找一组解，使z值最大
     *
     * 通过初等行变换，使得c >= 0
     * 因x >= 0，令所有非基础变量的值为0，此时 z = y 为最大值
     */
    Eigen::VectorXd c;
    Eigen::MatrixXd A;
    Eigen::VectorXd b;

    // 真正要求的变量
    var_num_ = constraint.cols();
    // 松弛变量
    slack_var_num_ = constraint.rows();
    // 人工变量
    artificial_var_num_ = 0;
    for (auto i = 0; i < maximum.rows(); i++) {
        if (maximum(i) < 0) {
            artificial_var_num_++;
        }
    }

    // 设置目标函数
    c.resize(var_num_ + slack_var_num_ + artificial_var_num_);
    c << -1 * object, Eigen::VectorXd::Zero(slack_var_num_),
       BIG_M * Eigen::VectorXd::Ones(artificial_var_num_);

    // 设置A, b
    A.resize(constraint.rows(), var_num_ + slack_var_num_ + artificial_var_num_);
    A << constraint,
      Eigen::MatrixXd::Identity(slack_var_num_, slack_var_num_),
      Eigen::MatrixXd::Zero(constraint.rows(), artificial_var_num_);
    b = maximum;

    // 设置人工变量
    auto j = 0;
    for (auto i = 0; i < b.rows(); i++) {
        if (b(i) >= 0) {
            continue;
        }
        A.row(i) *= -1;
        b.row(i) *= -1;
        A(i, constraint.cols() + slack_var_num_ + j) = 1;
        j++;
    }

    // 写成矩阵形式
    equations_ = Eigen::MatrixXd::Zero(A.rows() + 1, A.cols() + 1);
    equations_ << c.transpose(), 0,
                  A, b;

    std::cout << "equations: " << std::endl << equations_ << std::endl;
}

void LPSolver::set_variable(int idx, double value) {
}

void LPSolver::to_feasible_region() {
    if (artificial_var_num_ == 0) {
        return;
    }
    // 消除为负值的基础变量
    for (auto i = 1; i < equations_.rows(); i++) {
        for (auto j = var_num_ + slack_var_num_; j < equations_.cols() - 1; j++) {
            if (equations_(i, j) != 0) {
                pivot(i, j);
                break;
            }
        }
    }
}

std::tuple<bool, double, std::vector<double>> LPSolver::solve() {
    std::vector<double> vars;

    if (equations_.cols() == 1) {
        return std::make_tuple(equations_.minCoeff() >= 0, 0, vars);
    }

    to_feasible_region();
    // std::cout << "to_feasible_region: " << std::endl << equations_ << std::endl;

    int iterations = 0;
    while (!optimal()) {
        auto [row, col, has_solution] = next_basic_variable();
        if (!has_solution) {
            // std::cout << "LP no solution" << std::endl;
            return std::make_tuple(false, 0, vars);
        }
        pivot(row, col);

        // std::cout << "iteration " << iterations << ":" << std::endl << equations_ << std::endl;
        iterations++;
    }

    if (!feasible_solution()) {
        // std::cout << "LP no solution" << std::endl;
        return std::make_tuple(false, 0, vars);
    }

    auto maximize = equations_(0, equations_.cols() - 1);

    for (auto i = 0; i < var_num_; i++) {
        auto [is_basic, row] = is_basic_variable(i);
        if (is_basic) {
            auto var = equations_(row, equations_.cols()-1) / equations_(row, i);
            vars.push_back(var);
        } else {
            // 非基础变量全部制为0
            vars.push_back(0);
        }
    }

    std::cout << "============ linear programming solution =============" << std::endl;
    std::cout << "maximize: " << maximize << std::endl;
    std::cout << "variables: ";
    for (auto& var : vars) {
        std::cout << var << " ";
    }
    std::cout << std::endl;

    return std::make_tuple(true, maximize, vars);
}

std::tuple<double, int> LPSolver::min_object_coeff() {
    int col = 0;
    auto min = equations_.row(0).head(equations_.cols()-1).minCoeff(&col);
    return std::make_tuple(min, col);
}

std::tuple<bool, int> LPSolver::is_basic_variable(int col) {
    int row = 0;
    int non_zero_coeff_num = 0;
    for (auto i = 0; i < equations_.rows(); i++) {
        if (is_zero(equations_(i, col))) {
            continue;
        }
        row = i;
        non_zero_coeff_num++;
    }
    return std::make_tuple(non_zero_coeff_num == 1, row);
}

std::tuple<int, int, bool> LPSolver::next_basic_variable() {
    auto [dummy, col] = min_object_coeff();

    int row = 0;
    int b_col = equations_.cols() - 1;
    double min = 0;
    for (auto i = 1; i < equations_.rows(); i++) {
        if (equations_(i, col) <= 0) {
            continue;
        }
        double v = equations_(i, b_col) / equations_(i, col);
        if (row > 0 && min <= v) {
            continue;
        }
        row = i;
        min = v;
    }

    return std::make_tuple(row, col, row > 0);
}

bool LPSolver::optimal() {
    // 目标函数中所有系数都大于0，此时为最优解
    auto [min, dummy] = min_object_coeff();
    return min >= 0;
}

bool LPSolver::feasible_solution() {
    // 需要人工变量都不是基础变量，否则无解
    auto basic_vars_num = 0;
    for (auto i = 0; i < var_num_ + slack_var_num_; i++) {
        auto [is_basic, dummy] = is_basic_variable(i);
        if (is_basic) {
            basic_vars_num++;
        }
    }
    return basic_vars_num >= var_num_;
}

void LPSolver::pivot(int row, int col) {
    for (auto i = 0; i < equations_.rows(); i++) {
        if (i == row) {
            continue;
        }
        equations_.row(i) -= equations_.row(row) * equations_(i, col) / equations_(row, col);
    }
}
