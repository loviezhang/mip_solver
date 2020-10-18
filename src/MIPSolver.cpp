#include <iostream>
#include "MIPSolver.h"
#include "LPSolver.h"

class Node {
public:
    Node(const Eigen::VectorXd& object,
            const Eigen::MatrixXd& constraint,
            const Eigen::VectorXd& maximum,
            const int& bin_var_num);

    ~Node();

    // 是否可行
    bool feasible() const;

    // 是否是叶子节点
    bool leaf() const;

    // 是否已终止继续向下搜索
    bool terminated() const;

    // 设置为终止状态
    void terminate();

protected:
    const Eigen::VectorXd& object_;
    const Eigen::MatrixXd& constraint_;
    const Eigen::VectorXd& maximum_;
    const int& bin_var_num_;

    std::shared_ptr<Node> parent_;
    std::shared_ptr<Node> left_;
    std::shared_ptr<Node> right_;

    double value_;
    double estimite_;

    std::unordered_map<int, double> knowns_;
};

Node::Node(const Eigen::VectorXd& object,
        const Eigen::MatrixXd& constraint,
        const Eigen::VectorXd& maximum,
        const int& bin_var_num)
    : object_(object)
      , constraint_(constraint)
      , maximum_(maximum_)
      , bin_var_num_(bin_var_num_)
      , value_(0)
      , estimite_(0) {
}

Node::~Node() {
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

std::tuple<bool, float, std::vector<float>> MIPSolver::solve() {
    std::shared_ptr<Node> root;
    std::shared_ptr<Node> biggest;

    root = std::make_shared<Node>(object_, constraint_, maximum_, bin_var_num_);

    std::vector<float> vars;
    return std::make_tuple(false, 0, vars);
}
