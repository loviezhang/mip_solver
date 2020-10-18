#include <iostream>
#include <vector>
#include <set>
#include <Eigen/Dense>

#include "lp_solver.h"

int main(int argc, char** argv) {
    Eigen::VectorXf object;
    Eigen::MatrixXf constraint;
    Eigen::VectorXf maximum;

    /*
    object.resize(2);
    object << 7, 5;

    constraint.resize(2, 2);
    constraint << 2, 3,
                  3, 2;

    maximum.resize(2);
    maximum << 90, 120;
    */

    object.resize(2);
    object << -1, 10;

    constraint.resize(3, 2);
    constraint << -1, 5,
                  6, 5,
                  -1, -1;

    maximum.resize(3);
    maximum << 25, 60, -2;

    LPSolver solver;
    solver.init(object, constraint, maximum);
    solver.solve();

    return 0;
}

#if 0
#include "mip_solver.h"

struct Item {
    float score;
    std::set<int> tags;
};

struct Rule {
    int tag;
    float cost;
    int limit;
    int start;
    int end;
};

void add_item(std::vector<Item>& items, float score, std::set<int>&& tags) {
    Item item;
    item.score = score;
    item.tags = std::move(tags);
    items.push_back(item);
}

void add_rule(std::vector<Rule>& rules, int tag, float cost, int limit, int start, int end) {
    Rule rule;
    rule.tag = tag;
    rule.cost = cost;
    rule.limit = limit;
    rule.start = start;
    rule.end = end;
    rules.push_back(rule);
}

MIPSolver init_solver(const std::vector<Item>& items,
        const std::vector<Rule>& rules,
        int output_size) {
    Eigen::VectorXf score;
    Eigen::VectorXd start;
    Eigen::VectorXd end;
    Eigen::VectorXd limit;
    Eigen::VectorXf cost;
    Eigen::MatrixXd match;

    score.resize(items.size());
    for (auto i = 0UL; i < items.size(); i++) {
        score(i) = items[i].score;
    }

    start.resize(rules.size());
    end.resize(rules.size());
    limit.resize(rules.size());
    cost.resize(rules.size());
    for (auto i = 0UL; i < rules.size(); i++) {
        start(i) = rules[i].start;
        end(i) = rules[i].end;
        limit(i) = rules[i].limit;
        cost(i) = rules[i].cost;
    }

    match.resize(items.size(), rules.size());
    for (auto i = 0UL; i < items.size(); i++) {
        for (auto j = 0UL; j < rules.size(); j++) {
            if (items[i].tags.find(rules[j].tag) != items[i].tags.end()) {
                match(i, j) = 1;
            }
        }
    }

    std::cout << "score: " << score.transpose() << std::endl;
    std::cout << "start: " << start.transpose() << std::endl;
    std::cout << "end: " << end.transpose() << std::endl;
    std::cout << "limit: " << limit.transpose() << std::endl;
    std::cout << "cost: " << cost.transpose() << std::endl;
    std::cout << "match: " << std::endl << match << std::endl;

    MIPSolver solver;
    solver.init(items.size(), rules.size(), output_size, score, start, end, limit, cost, match);
    return solver;
}

void test_mip_solver() {
    std::vector<Item> items;
    add_item(items, 0.9, {1, 2, 3});
    add_item(items, 0.8, {3, 4, 5});
    add_item(items, 0.7, {1, 4, 5});
    add_item(items, 0.6, {2, 3});

    std::vector<Rule> rules;
    add_rule(rules, 1, 100, 1, 0, 1);
    add_rule(rules, 2, 100, 1, 0, 1);
    add_rule(rules, 3, 100, 2, 0, 1);
    add_rule(rules, 4, 100, 1, 0, 1);
    add_rule(rules, 5, 100, 1, 0, 1);

    int output_size = 2;

    init_solver(items, rules, output_size);
}

int main(int argc, char** argv) {
    /*
    LPSolver solver;
    solver.init();
    solver.solve();
    */
    test_mip_solver();

    return 0;
}
#endif
