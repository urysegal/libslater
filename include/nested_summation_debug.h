#pragma once
#include <vector>
#include <string>
#include <boost/core/demangle.hpp>

#include <iostream>
#include <typeinfo>
#include <assert.h>

namespace slater {

void print_tabs(int t);


template<typename indexing_t, class T>
struct summation_debug_state;

template<typename indexing_t, class T>
struct summation_debug_item {
    bool is_zero = false;
    bool is_loaded = false;
    indexing_t current_index_value;
    T scaling;
    T inner_value;
    T new_value; /// == scaling * inner_value
    summation_debug_state<indexing_t, T> *item = nullptr;

    void dump(const std::string &var_name, int tab_depth)
    {
        print_tabs(tab_depth);
        assert(is_loaded);
        if ( is_zero ) {
            std::cerr << var_name << " = " << current_index_value <<  ": 0.0 - scaling was " << scaling << std::endl;
        } else {
            std::cerr << var_name << " = " << current_index_value <<  ": " << new_value << " = ( " << scaling << " * " << inner_value << ")" << std::endl;
            item->dump(tab_depth+1);
        }
    }

    void cleanup()
    {
        if ( item ) {
            item->cleanup();
            delete item;
        }
    }

};

template<typename indexing_t, class T>
struct summation_debug_state {
    std::string class_name;
    std::string index_variable_name;
    indexing_t from;
    indexing_t to;
    indexing_t step;
    T final_value;
    summation_debug_state<indexing_t, T> *parent = nullptr;
    std::vector<summation_debug_item<indexing_t, T> *> items;

    void dump(int tab_depth)
    {
        print_tabs(tab_depth);
        std::cerr << class_name << " Result " << final_value
            << " index " << index_variable_name << " from " << from  << " to " << to
            << " step of " << step << std::endl;
        for ( auto const &it: items ) {
            it->dump(index_variable_name, tab_depth+1);
        }
    }

    void cleanup()
    {
        for ( auto it: items ) {
            it->cleanup();
            delete it;
        }
    }
};

// Track the status of the summation for debug purposes
template<typename indexing_t, class T>
class Summation_Debug_State {

public:
    void
    push_state(const char *class_name, const char *index_variable_name, indexing_t from, indexing_t to, indexing_t step)
    {
        if ( paused ) {
            return;
        }
        auto new_state = new summation_debug_state<indexing_t, T>();

        if ( not current_state ) {
            top_state = new_state;
        }
        new_state->index_variable_name = index_variable_name;
        new_state->class_name = boost::core::demangle(class_name);
        new_state->from = from;
        new_state->to = to;
        new_state->step = step;
        new_state->parent = current_state;
        current_state = new_state;
    }

    void add_item(indexing_t current_index_value, T scaling, T inner_value, T new_value)
    {
        if ( paused ) {
            return;
        }

        if ( not this->current_state->items.empty() ) {
            auto new_item = this->current_state->items.back();
            if (new_item->is_loaded == false) {
                new_item->current_index_value = current_index_value;
                new_item->scaling = scaling;
                new_item->inner_value = inner_value;
                new_item->new_value = new_value;
                new_item->is_loaded = true;
            }
        }
    }

    void add_zero_item(indexing_t current_index_value, T scaling)
    {
        if ( paused ) {
            return;
        }

        auto new_item = new summation_debug_item<indexing_t, T>();
        new_item->is_loaded = true;
        new_item->is_zero = true;
        new_item->current_index_value = current_index_value;
        new_item->scaling = scaling;
        current_state->items.push_back(new_item) ;
    }

    void pop_state(T total_sum)
    {
        if ( paused ) {
            return;
        }

        current_state->final_value = total_sum;
        if (current_state != top_state) {
            auto new_item = new summation_debug_item<indexing_t, T>();
            new_item->item = current_state;
            current_state->parent->items.push_back(new_item);
            current_state = current_state->parent;
        } else {
            dump();
            cleanup();
        }
    }

    void dump() {
        top_state->dump(0);
    }

    void cleanup()
    {
        top_state->cleanup();
        delete top_state;
        top_state = nullptr;
        current_state = nullptr;
    }

    void pause() { paused = true; }
    void resume() { paused = false; }


private:
    summation_debug_state<indexing_t, T> *current_state = nullptr;
    summation_debug_state<indexing_t, T> *top_state = nullptr;

    bool paused = false;
};

}