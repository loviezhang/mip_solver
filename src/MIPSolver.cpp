#include <iostream>
#include <algorithm>
#include <list>
#include "MIPSolver.h"
#include "LPSolver.h"

Node::Node(const Eigen::VectorXd& object,
        const Eigen::MatrixXd& constraint,
        const Eigen::VectorXd& maximum,
        const int& bin_var_num)
    : object_(object)
      , constraint_(constraint)
      , maximum_(maximum)
      , bin_var_num_(bin_var_num)
      , estimite_(0)
      , value_(0)
      , feasible_(true)
      , split_point_(0)
      , left_explored_(false)
      , right_explored_(false)
      , on_left_side_(false) {
    init();
}

Node::~Node() {
}

void Node::init() {
    LPSolver solver;
    solver.init(object_, constraint_, maximum_);

    std::tie(feasible_, estimite_, variables_) = solver.solve();

    std::cout << "object: " << object_.transpose();
    std::cout << ", maximum: " << maximum_.transpose();
    std::cout << ", feasible: " << feasible_;
    std::cout << ", estimite: " << estimite_;
    std::cout << ", variables: ";
    for (auto& var : variables_) {
        std::cout << var << " ";
    }
    std::cout << std::endl;
}

std::tuple<bool, int> Node::calc_split_point() const {
    if (left_explored_) {
        return std::make_tuple(false, split_point_);
    }
    if (right_explored_) {
        return std::make_tuple(true, split_point_);
    }

    // 选取最接近0或1的变量
    double diff_min = 1;
    int col = 0;
    for (int i = 0; i < bin_var_num_; i++) {
        double diff = std::min(variables_[i], 1-variables_[i]);
        if (diff_min > diff) {
            diff_min = diff;
            col = i;
        }
    }

    bool go_left = true;
    if (variables_[col] < 1-variables_[col]) {
        go_left = false;
    }

    return std::make_tuple(go_left, col);
}

std::shared_ptr<Node> Node::spawn() {
    auto [go_left, split_point_] = calc_split_point();

    // 将对应位置的变量设置成已知数
    Eigen::VectorXd object(object_.rows()-1);
    Eigen::MatrixXd constraint(constraint_.rows(), constraint_.cols()-1);
    Eigen::VectorXd maximum(maximum_.rows());
    int bin_var_num = bin_var_num_ - 1;

    object << object_.head(split_point_), object_.tail(object_.rows()-split_point_-1);
    constraint << constraint_.block(0, 0, constraint_.rows(), split_point_),
                  constraint_.block(0, split_point_+1, constraint_.rows(), constraint_.cols()-split_point_-1);
    if (go_left) {
        maximum = maximum_ - constraint_.col(split_point_);
    } else {
        maximum = maximum_;
    }

    auto node = std::make_shared<Node>(object, constraint, maximum, bin_var_num);
    node->parent_ = shared_from_this();
    if (go_left) {
        left_explored_ = true;
        node->value_ = value_ + object_(split_point_);
    } else {
        right_explored_ = true;
        node->value_ = value_;
    }
    node->on_left_side_ = go_left;

    if (go_left) {
        std::cout << "go left, split point " << split_point_
            << " estimite " << node->estimite() << std::endl;
    } else {
        std::cout << "go right, split point " << split_point_
            << " estimite " << node->estimite() << std::endl;
    }

    return node;
}


MIPSolver::MIPSolver() {
}

MIPSolver::~MIPSolver() {
}

void MIPSolver::init(const Eigen::VectorXd& object,
        const Eigen::MatrixXd& constraint,
        const Eigen::VectorXd& maximum,
        int bin_var_num) {
    object_ = object;
    bin_var_num_ = bin_var_num;

    // 添加0-1变量的约束
    constraint_.resize(constraint.rows() + bin_var_num, constraint.cols());
    constraint_ << constraint,
                   Eigen::MatrixXd::Identity(bin_var_num, bin_var_num),
                   Eigen::MatrixXd::Zero(bin_var_num, constraint.cols() - bin_var_num);
    maximum_.resize(maximum.rows() + bin_var_num);
    maximum_ << maximum, Eigen::VectorXd::Ones(bin_var_num);
}

void MIPSolver::update_biggest_node(std::shared_ptr<Node>& cur) {
    if (!cur->is_leaf()) {
        return;
    }
    if (!cur->feasible()) {
        return;
    }
    // 对于叶子节点，estimite就是结果
    if (biggest_ == nullptr || biggest_->estimite() < cur->estimite()) {
        std::cout << "set biggest node" << std::endl;
        biggest_ = cur;
    }
}

bool MIPSolver::need_backtrack(std::shared_ptr<Node>& cur) {
    // 走到叶子节点
    if (cur->is_leaf()) {
        std::cout << "is_leaf" << std::endl;
        return true;
    }
    // 不可解
    if (!cur->feasible()) {
        std::cout << "infeasible" << std::endl;
        return true;
    }
    // 已搜索完成
    if (cur->has_explored()) {
        std::cout << "has explored" << std::endl;
        return true;
    }
    // 预估值小于biggest_node
    if (biggest_ != nullptr && cur->estimite() < biggest_->estimite()) {
        std::cout << "less than biggest" << std::endl;
        return true;
    }
    return false;
}

std::tuple<bool, double, std::vector<double>> MIPSolver::get_result() const {
    std::vector<double> vars;

    if (biggest_ == nullptr) {
        // 无解
        std::cout << "MIP no solution" << std::endl;
        return std::make_tuple(false, 0, vars);
    }

    std::list<std::shared_ptr<Node>> path;
    auto cur = biggest_;
    while (cur != nullptr) {
        path.push_front(cur);
        cur = cur->backtrack();
    }

    std::vector<int> var_index(bin_var_num_);
    for (auto i = 0; i < bin_var_num_; i++) {
        var_index[i] = i;
    }
    vars.resize(object_.rows(), -1.0);

    cur = path.front();
    path.pop_front();
    for (auto& child : path) {
        auto split_point = cur->split_point();
        for (auto i = 0; i < split_point || vars[i] != -1; i++) {
            if (vars[i] != -1) {
                split_point++;
            }
        }

        vars[split_point] = child->on_left_side() ? 1.0 : 0.0;

        cur = child;
    }

    std::copy(biggest_->variables().begin(), biggest_->variables().end(),
            vars.begin() + bin_var_num_);

    std::cout << "============ mixed integer programming solution =============" << std::endl;
    std::cout << "maximize: " << biggest_->estimite() << std::endl;
    std::cout << "variables: ";
    for (auto& var : vars) {
        std::cout << var << " ";
    }
    std::cout << std::endl;

    return std::make_tuple(true, biggest_->estimite(), vars);
}

std::tuple<bool, double, std::vector<double>> MIPSolver::solve() {
    auto cur = std::make_shared<Node>(object_, constraint_, maximum_, bin_var_num_);

    // 深度优先搜索
    while (cur != nullptr) {
        update_biggest_node(cur);

        // 探索到底部，回溯
        if (need_backtrack(cur)) {
            std::cout << "backtrack" << std::endl;
            cur = cur->backtrack();
            continue;
        }

        std::cout << "spawn" << std::endl;
        cur = cur->spawn();
    }
    std::cout << "search finish" << std::endl;

    return get_result();
}
