#include <iostream>
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
      , feasible_(false)
      , terminated_(false) {
    init();
}

Node::~Node() {
}

void Node::init() {
    LPSolver solver;
    solver.init(object_, constraint_, maximum_);

    std::vector<double> vars;
    std::tie(feasible_, estimite_, variables_) = solver.solve();
}

void Node::terminate() {
    terminated_ = true;
    left_ = nullptr;
    right_ = nullptr;
}

std::tuple<bool, int> Node::split_point() const {
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
    if (variables_[col] > 1-variables_[col]) {
        go_left = false;
    }

    return std::make_tuple(go_left, col);
}

std::shared_ptr<Node> Node::spawn() {
    auto [go_left, col] = split_point();

    // 将对应位置的变量设置成已知数
    Eigen::VectorXd object(object_.rows()-1);
    Eigen::MatrixXd constraint(constraint_.rows(), constraint_.cols()-1);
    Eigen::VectorXd maximum(maximum_.rows());
    int bin_var_num = bin_var_num_ - 1;

    object << object_.head(col), object_.tail(object_.rows()-col-1);
    constraint << constraint_.block(0, 0, constraint_.rows(), col),
                  constraint_.block(0, col+1, constraint_.rows(), constraint_.cols()-col-1);
    if (go_left) {
        maximum = maximum_ - constraint_.col(col);
    } else {
        maximum = maximum_;
    }

    auto node = std::make_shared<Node>(object, constraint, maximum, bin_var_num);
    node->parent_ = shared_from_this();
    if (go_left) {
        left_ = node;
        node->value_ = value_ + object_(col);
    } else {
        right_ = node;
        node->value_ = value_;
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
    if (biggest_ == nullptr || biggest_->estimite() < cur->estimite()) {
        biggest_ = cur;
    }
}

bool MIPSolver::need_backtrack(std::shared_ptr<Node>& cur) {
    // 走到叶子节点
    if (cur->is_leaf()) {
        return true;
    }
    // 不可解
    if (cur->feasible()) {
        return true;
    }
    // 已设置为终止搜索
    if (cur->terminated()) {
        return true;
    }
    // 预估值小于biggest_node
    if (biggest_ != nullptr && cur->estimite() < biggest_->estimite()) {
        return true;
    }
    return false;
}

std::tuple<bool, float, std::vector<float>> MIPSolver::solve() {
    auto root = std::make_shared<Node>(object_,
            constraint_, maximum_, bin_var_num_);
    auto cur = root;

    // 深度优先搜索
    while (cur != nullptr) {
        update_biggest_node(cur);

        // 探索到底部，回溯
        if (need_backtrack(cur)) {
            cur->terminate();
            cur = cur->backtrack();
            continue;
        }

        cur = cur->spawn();
    }

    std::vector<float> vars;
    return std::make_tuple(false, 0, vars);
}
