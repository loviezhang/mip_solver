#include <iostream>
#include <fstream>
#include <functional>
#include "RankRule.h"

std::vector<std::string> split(const std::string& s, char seperator) {
    std::vector<std::string> output;

    std::string::size_type prev_pos = 0, pos = 0;

    while((pos = s.find(seperator, pos)) != std::string::npos)
    {
        std::string substring( s.substr(prev_pos, pos-prev_pos) );

        output.push_back(substring);

        prev_pos = ++pos;
    }

    output.push_back(s.substr(prev_pos, pos-prev_pos)); // Last word

    return output;
}

bool RankRule::load(const std::string& item_file, const std::string& rule_file) {
    if (!load_item(item_file)) {
        return false;
    }
    if (!load_rule(rule_file)) {
        return false;
    }
    expand_rule();

    return true;
}

void RankRule::expand_rule() {
}

bool RankRule::load_item(const std::string& file) {
    std::ifstream ifs;
    ifs.open(file, std::fstream::in);
    std::string line;
    while (std::getline(ifs, line)) {
        auto words = split(line, ' ');

        Item item;
        item.gid = std::stoul(words[0]);
        item.score = std::stod(words[1]);
        auto tags = split(words[2], ',');
        for (auto& tag : tags) {
            item.tags.insert(std::stoul(tag));
        }
        auto cat_tags = split(words[3], ',');
        for (auto& s : cat_tags) {
            auto cat_tag = split(s, '|');
            item.cat_tags.push_back(std::make_tuple(std::stoul(cat_tag[0]), std::stoul(cat_tag[1])));
        }

        items_.push_back(std::move(item));
    }

    if (items_.empty()) {
        return true;
    }

    // 分数归一，修复分数太大的问题
    double max = items_[0].score;
    double min = items_[0].score;
    for (auto& item : items_) {
        if (max < item.score) {
            max = item.score;
        }
        if (min > item.score) {
            min = item.score;
        }
    }
    for (auto& item : items_) {
        item.score = (item.score - min / (max - min)) / 2.0 + 0.5;
    }

    return true;
}

bool RankRule::load_rule(const std::string& file) {
    std::ifstream ifs;
    ifs.open(file, std::fstream::in);

    std::string id;
    while (std::getline(ifs, id)) {
        std::string cat_or_tag, type, size, max_or_min;
        std::getline(ifs, cat_or_tag);
        std::getline(ifs, type);
        std::getline(ifs, size);
        std::getline(ifs, max_or_min);

        auto w1 = split(id, ' ');
        auto w2 = split(cat_or_tag, ' ');
        auto w3 = split(type, ' ');
        auto w4 = split(size, ' ');
        auto w5 = split(max_or_min, ' ');

        Rule rule;
        rule.id = std::stoi(w1[1]);
        if (w2[0] == "tags") {
            rule.tag = std::hash<std::string>{}(w2[1]);
        } else if (w2[0] == "category") {
            rule.cat = std::stoul(w2[1]);
        }
        rule.slide = w3[1] == "SLIDE";
        rule.size = std::stoi(w4[1]);
        if (w5[0] == "min") {
            rule.min = std::stoi(w5[1]);
        } else if (w5[0] == "max") {
            rule.max = std::stoi(w5[1]);
        }

        std::cout << "id:" << rule.id
            << ", cat:" << rule.cat
            << ", tag:" << rule.tag
            << ", min:" << rule.min
            << ", max:" << rule.max
            << ", size:" << rule.size
            << ", slide:" << rule.slide
            << std::endl;

        rules_.push_back(rule);
    }

    return true;
}
