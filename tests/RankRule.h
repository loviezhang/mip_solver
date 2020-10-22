#pragma once

#include <vector>
#include <tuple>
#include <set>

struct Item {
    uint64_t gid;
    double score;
    std::vector<std::tuple<uint64_t, uint64_t>> cat_tags;
    std::set<uint64_t> tags;

    Item() {
        gid = 0;
        score = 0;
    }
};

struct Rule {
    int id;
    uint64_t cat;
    uint64_t tag;
    int min;
    int max;
    int size;
    bool slide;

    Rule() {
        id = 0;
        cat = 0;
        tag = 0;
        min = 0;
        max = 0;
        size = 0;
        slide = false;
    }
};

class RankRule {
public:
    bool load(const std::string& item_file, const std::string& rule_file);

    const std::vector<Rule>& get_rule() const {
        return rules_;
    }

    const std::vector<Item>& get_item() const {
        return items_;
    }

private:
    bool load_item(const std::string& file);
    bool load_rule(const std::string& file);

    void expand_rule();

private:
    std::vector<Rule> rules_;
    std::vector<Item> items_;
};
