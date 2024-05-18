#include "domain.h"

void Domain::_bind_methods() {
    ClassDB::bind_method(D_METHOD("_m_verify_g", "state", "method", "state_var", "arg", "desired_val", "depth"), &Domain::_m_verify_g);
    ClassDB::bind_static_method("Domain", D_METHOD("_goals_not_achieved", "state", "multigoal"), &Domain::_goals_not_achieved);
    ClassDB::bind_method(D_METHOD("_m_verify_mg", "state", "method", "multigoal", "depth"), &Domain::_m_verify_mg);
    ClassDB::bind_method(D_METHOD("display"), &Domain::display);
}

Variant Domain::_m_verify_g(Dictionary state, String method, String state_var, String arg, Variant desired_val, int depth) {
    Dictionary state_dict = state[state_var]; 

    if (state_dict[arg] != desired_val) {
        if (verbose >= 3) {
            print_line(vformat("Depth %d: method %s didn't achieve\nGoal %s[%s] = %s", depth, method, state_var, arg, desired_val));
        }
        return false;
    }

    // if (!stn->is_consistent()) {
    //     if (verbose >= 3) {
    //         print_line(vformat("Depth %d: method %s resulted in inconsistent STN for %s", depth, method));
    //     }
    //     return false;
    // }

    if (verbose >= 3) {
        print_line(vformat("Depth %d: method %s achieved\nGoal %s[%s] = %s", depth, method, state_var, arg, desired_val));
    }
    return Array();
}

Dictionary Domain::_goals_not_achieved(Dictionary state, Ref<Multigoal> multigoal) {
    Dictionary unachieved;
    Array keys = multigoal->get_state().keys();
    for (int i = 0; i < keys.size(); ++i) {
        String n = keys[i];
        Dictionary sub_dict = multigoal->get_state()[n];
        Array sub_keys = sub_dict.keys();
        for (int j = 0; j < sub_keys.size(); ++j) {
            String arg = sub_keys[j];
            Variant val = sub_dict[arg];
            if (state[n].get_type() == Variant::DICTIONARY && Dictionary(state[n]).has(arg) && val != Dictionary(state[n])[arg]) {
                if (!unachieved.has(n)) {
                    unachieved[n] = Dictionary();
                }
                Dictionary temp = unachieved[n];
                temp[arg] = val;
                unachieved[n] = temp;
            }
        }
    }
    return unachieved;
}

Variant Domain::_m_verify_mg(Dictionary state, String method, Ref<Multigoal> multigoal, int depth) {
    Dictionary goal_dict = _goals_not_achieved(state, multigoal);
    if (!goal_dict.is_empty()) {
        if (verbose >= 3) {
            print_line(vformat("Depth %d: method %s didn't achieve %s", depth, method, multigoal));
        }
        return false;
    }

    // if (!stn->is_consistent()) {
    //     if (verbose >= 3) {
    //         print_line(vformat("Depth %d: method %s resulted in inconsistent STN for %s", depth, method));
    //     }
    //     return false;
    // }

    if (verbose >= 3) {
        print_line(vformat("Depth %d: method %s achieved %s", depth, method, multigoal));
    }
    return Array();
}


void Domain::display() {
    print_line(to_string());
}
