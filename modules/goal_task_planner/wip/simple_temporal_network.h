/**************************************************************************/
/*  simple_temporal_network.h                                             */
/**************************************************************************/
/*                         This file is part of:                          */
/*                             GODOT ENGINE                               */
/*                        https://godotengine.org                         */
/**************************************************************************/
/* Copyright (c) 2014-present Godot Engine contributors (see AUTHORS.md). */
/* Copyright (c) 2007-2014 Juan Linietsky, Ariel Manzur.                  */
/*                                                                        */
/* Permission is hereby granted, free of charge, to any person obtaining  */
/* a copy of this software and associated documentation files (the        */
/* "Software"), to deal in the Software without restriction, including    */
/* without limitation the rights to use, copy, modify, merge, publish,    */
/* distribute, sublicense, and/or sell copies of the Software, and to     */
/* permit persons to whom the Software is furnished to do so, subject to  */
/* the following conditions:                                              */
/*                                                                        */
/* The above copyright notice and this permission notice shall be         */
/* included in all copies or substantial portions of the Software.        */
/*                                                                        */
/* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,        */
/* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF     */
/* MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. */
/* IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY   */
/* CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,   */
/* TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE      */
/* SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.                 */
/**************************************************************************/

#ifndef SIMPLE_TEMPORAL_NETWORK_H
#define SIMPLE_TEMPORAL_NETWORK_H
// Copyright (c) 2023-present. This file is part of V-Sekai https://v-sekai.org/.
// K. S. Ernest (Fire) Lee & Contributors (see .all-contributorsrc).
// simple_temporal_network.h
// SPDX-License-Identifier: MIT

#include "core/io/resource.h"
#include "core/math/vector2i.h"
#include "temporal_constraint.h"

class SimpleTemporalNetwork : public Resource {
	GDCLASS(SimpleTemporalNetwork, Resource);

private:
	Array constraints;
	int num_nodes = 0;
	Array node_intervals;
	Dictionary node_indices;

	Dictionary node_index_cache;

public:
	Array getConstraints() const { return constraints; }
	int getNumNodes() const { return num_nodes; }
	Array getNodeIntervals() const { return node_intervals; }
	Dictionary getNodeIndices() const { return node_indices; }
	Dictionary getNodeIndexCache() const { return node_index_cache; }

	void setConstraints(const Array &value) { constraints = value; }
	void setNumNodes(int value) { num_nodes = value; }
	void setNodeIntervals(const Array &value) { node_intervals = value; }
	void setNodeIndices(const Dictionary &value) { node_indices = value; }
	void setNodeIndexCache(const Dictionary &value) { node_index_cache = value; }

	String _to_string();
	int get_node_index(int time_point);
	bool check_overlap(Ref<TemporalConstraint> new_constraint);
	bool add_temporal_constraint(Ref<TemporalConstraint> from_constraint, Ref<TemporalConstraint> to_constraint = nullptr, float min_gap = 0, float max_gap = 0);
	void update_constraints_list(Ref<TemporalConstraint> constraint, Ref<TemporalConstraint> node);
	void add_constraints_to_list(Ref<TemporalConstraint> from_constraint, Ref<TemporalConstraint> to_constraint);
	Ref<TemporalConstraint> process_constraint(Ref<TemporalConstraint> constraint);
	Ref<TemporalConstraint> get_temporal_constraint_by_name(String constraint_name);
	bool is_consistent();
	Array enumerate_decompositions(Ref<TemporalConstraint> vertex);
	bool is_leaf(Ref<TemporalConstraint> vertex);
	bool is_or(Ref<TemporalConstraint> vertex);
	Array get_children(Ref<TemporalConstraint> vertex);
	Array cartesian_product(Array arrays);
	void update_state(Dictionary state);

protected:
	static void _bind_methods();
};

#endif // SIMPLE_TEMPORAL_NETWORK_H
