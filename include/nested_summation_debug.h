#pragma once
#include <vector>
#include <string>

template<typename indexing_t, class T> struct summation_debug_state;

template<typename indexing_t, class T>
struct summation_debug_item {
    bool is_zero = false;
    indexing_t current_index_value;
    T scaling;
    T inner_value;
    T new_item;
    summation_debug_state<indexing_t, T> items;
};

template<typename indexing_t, class T>
struct summation_debug_state {
    std::string class_name;
    std::string index_variable_name;
    T final_value;
    std::vector<summation_debug_item<indexing_t, T> > items;
};

// Track the status of the summation for debug purposes
template<typename indexing_t, class T>
class Summation_Debug_State {

public:
    void push_state(const char *class_name, const char *index_variable_name, indexing_t from, indexing_t to, indexing_t step)
    {
        current_variables.push_back({class_name, index_variable_name});
    }
    void add_item(indexing_t current_index_value, T scaling, T inner_value, T new_item) {}
    void add_zero_item(indexing_t current_index_value, T scaling) {}
    void pop_state(T total_sum)
    {
        current_variables.pop_back();
    }

private:
    std::vector< std::pair<std::string, std::string> > current_variables;
};
