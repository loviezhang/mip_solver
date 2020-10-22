#include <gtest/gtest.h>
#include <vector>
#include <set>

#include "utils.h"
#include "../src/MIPSolver.h"
#include "RankRule.h"

#if 1

/*
struct Item {
    double score;
    std::set<int> tags;
};

struct Rule {
    int tag;
    double cost;
    int limit;
    int start;
    int end;
};

void add_item(std::vector<Item>& items, double score, std::set<int>&& tags) {
    Item item;
    item.score = score;
    item.tags = std::move(tags);
    items.push_back(item);
}

void add_rule(std::vector<Rule>& rules, int tag, double cost, int limit, int start, int end) {
    Rule rule;
    rule.tag = tag;
    rule.cost = cost;
    rule.limit = limit;
    rule.start = start;
    rule.end = end;
    rules.push_back(rule);
}

MIPSolver init(int item_size,
        int rule_size,
        int output_size,
        const Eigen::VectorXd& score,
        const Eigen::VectorXd& start,
        const Eigen::VectorXd& end,
        const Eigen::VectorXd& limit,
        const Eigen::VectorXd& cost,
        const Eigen::MatrixXd& match) {

    MIPSolver solver;

    auto bin_var_num = item_size * output_size;
    auto double_var_num = rule_size;
    auto constraint_num = 2
        + item_size
        + output_size
        + rule_size
        ;

    Eigen::VectorXd object(bin_var_num + double_var_num);
    Eigen::MatrixXd constraint(constraint_num, bin_var_num + double_var_num);
    Eigen::VectorXd maximum(constraint_num);

    // 目标函数
    for (auto i = 0; i < item_size; i++) {
        for (auto j = 0; j < output_size; j++) {
            object(i*output_size+j) = score(i);
        }
    }
    for (auto i = 0; i < rule_size; i++) {
        object(item_size*output_size+i) = -1.0 * cost(i);
    }

    auto idx = 0;

    // 输出数量 = output_size
    for (auto i = 0; i < item_size; i++) {
        for (auto j = 0; j < output_size; j++) {
            constraint(idx+0, i*output_size+j) = 1;
            constraint(idx+1, i*output_size+j) = -1;
        }
    }
    maximum(idx+0) = output_size;
    maximum(idx+1) = -1 * output_size;

    idx += 2;

    // 一个视频只能放一个位置
    for (auto i = 0; i < item_size; i++) {
        for (auto j = 0; j < output_size; j++) {
            constraint(idx+i, i*output_size+j) = 1;
        }
        maximum(idx+i) = 1;
    }

    idx += item_size;

    // 一个位置只能放一个视频
    for (auto j = 0; j < output_size; j++) {
        for (auto i = 0; i < item_size; i++) {
            constraint(idx+j, i*output_size+j) = 1;
        }
        maximum(idx+j) = 1;
    }
    idx += output_size;

    // 规则限制
    for (auto k = 0; k < rule_size; k++) {
        for (auto i = 0; i < item_size; i++) {
            for (auto t = start(k); t <= end(k); t++) {
                constraint(idx+k, i*output_size+t) = match(i, k);
            }
        }
        constraint(idx+k, item_size*output_size+k) = -1;
        maximum(idx+k) = limit(k);
    }

    idx += rule_size;

    solver.init(object, constraint, maximum, bin_var_num);
    return solver;
}

MIPSolver init(const std::vector<Item>& items,
        const std::vector<Rule>& rules,
        int output_size) {
    Eigen::VectorXd score;
    Eigen::VectorXd start;
    Eigen::VectorXd end;
    Eigen::VectorXd limit;
    Eigen::VectorXd cost;
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
            } else {
                match(i, j) = 0;
            }
        }
    }

    std::cout << "score: " << score.transpose() << std::endl;
    std::cout << "start: " << start.transpose() << std::endl;
    std::cout << "end: " << end.transpose() << std::endl;
    std::cout << "limit: " << limit.transpose() << std::endl;
    std::cout << "cost: " << cost.transpose() << std::endl;
    std::cout << "match: " << std::endl << match << std::endl;

    return init(items.size(), rules.size(), output_size,
            score, start, end, limit, cost, match);
}
*/

TEST(RankRule, Test1) {
    /*
    std::vector<Item> items;
    add_item(items, 0.9, {1});
    add_item(items, 0.8, {1, 2});
    add_item(items, 0.7, {2});

    std::vector<Rule> rules;
    add_rule(rules, 1, 100, 1, 0, 1);
    add_rule(rules, 2, 100, 1, 0, 1);

    int output_size = 2;

    auto solver = init(items, rules, output_size);
    auto [has_solution, maximize, variables] = solver.solve();
    EXPECT_TRUE(has_solution);
    */

    RankRule rank_rule;
    rank_rule.load("./data/202010120932330101151761700A1FEB66", "./data/rank_rule");
    EXPECT_TRUE(true);
}

#endif
