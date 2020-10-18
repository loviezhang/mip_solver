#include <iostream>
#include "MIPSolver.h"

MIPSolver::MIPSolver() {
}

MIPSolver::~MIPSolver() {
}

void MIPSolver::init(const Eigen::VectorXd& object,
        const Eigen::MatrixXd& constraint,
        const Eigen::VectorXd& maximum,
        int bin_var_num) {
}

std::tuple<bool, float, std::vector<float>> MIPSolver::solve() {
    std::vector<float> vars;
    return std::make_tuple(false, 0, vars);
}
