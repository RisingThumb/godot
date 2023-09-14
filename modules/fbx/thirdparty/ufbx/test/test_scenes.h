#undef UFBXT_TEST_GROUP
#define UFBXT_TEST_GROUP "scenes"

UFBXT_FILE_TEST_FLAGS(maya_slime, UFBXT_FILE_TEST_FLAG_HEAVY_TO_FUZZ)
#if UFBXT_IMPL
{
	ufbx_node *node_high = ufbx_find_node(scene, "Slime_002:Slime_Body_high");
	ufbxt_assert(node_high);
	ufbxt_assert(!node_high->visible);
}
#endif

UFBXT_FILE_TEST(blender_293_barbarian)
#if UFBXT_IMPL
{
}
#endif

UFBXT_FILE_TEST_ALT(evaluate_alloc_fail, blender_293_barbarian)
#if UFBXT_IMPL
{
	for (size_t max_temp = 1; max_temp < 10000; max_temp++) {
		ufbx_evaluate_opts opts = { 0 };
		opts.temp_allocator.huge_threshold = 1;
		opts.temp_allocator.allocation_limit = max_temp;
		opts.evaluate_skinning = true;

		ufbxt_hintf("Temp limit: %zu", max_temp);

		ufbx_error error;
		ufbx_scene *eval_scene = ufbx_evaluate_scene(scene, NULL, 0.2, &opts, &error);
		if (eval_scene) {
			ufbxt_logf(".. Tested up to %zu temporary allocations", max_temp);
			ufbx_free_scene(eval_scene);
			break;
		}
		ufbxt_assert(error.type == UFBX_ERROR_ALLOCATION_LIMIT);
	}

	for (size_t max_result = 1; max_result < 10000; max_result++) {
		ufbx_evaluate_opts opts = { 0 };
		opts.result_allocator.huge_threshold = 1;
		opts.result_allocator.allocation_limit = max_result;
		opts.evaluate_skinning = true;

		ufbxt_hintf("Result limit: %zu", max_result);

		ufbx_error error;
		ufbx_scene *eval_scene = ufbx_evaluate_scene(scene, NULL, 0.2, &opts, &error);
		if (eval_scene) {
			ufbxt_logf(".. Tested up to %zu result allocations", max_result);
			ufbx_free_scene(eval_scene);
			break;
		}
		ufbxt_assert(error.type == UFBX_ERROR_ALLOCATION_LIMIT);
	}
}
#endif

UFBXT_FILE_TEST(maya_kenney_character)
#if UFBXT_IMPL
{
	ufbxt_check_frame(scene, err, false, "maya_kenney_character_4", NULL, 4.0/24.0);
	ufbxt_check_frame(scene, err, false, "maya_kenney_character_9", NULL, 9.0/24.0);

	{
		ufbx_node *node = ufbx_find_node(scene, "characterMedium");
		ufbxt_assert(node && node->mesh);
		ufbx_mesh *mesh = node->mesh;
		ufbxt_assert(mesh->skin_deformers.count == 1);
	}

	for (int frame = 3; frame <= 15; frame += 10) {
		double time = (double)frame / scene->settings.frames_per_second;

		ufbx_scene *state = ufbx_evaluate_scene(scene, NULL, time, NULL, NULL);
		ufbxt_assert(state);

		ufbx_node *node = ufbx_find_node(scene, "characterMedium");
		ufbxt_assert(node && node->mesh);
		ufbx_mesh *mesh = node->mesh;
		ufbxt_assert(mesh->skin_deformers.count == 1);

		for (size_t level = 0; level <= 2; level++) {
			ufbx_mesh *sub_mesh = mesh;
			if (level > 0) {
				ufbx_subdivide_opts opts = { 0 };
				opts.evaluate_source_vertices = true;
				opts.evaluate_skin_weights = true;
				sub_mesh = ufbx_subdivide_mesh(mesh, level, &opts, NULL);
				ufbxt_assert(sub_mesh);

				ufbxt_check_source_vertices(sub_mesh, mesh, err);
			} else {
				ufbx_retain_mesh(sub_mesh);
			}

			ufbx_skin_deformer *skin = mesh->skin_deformers.data[0];

			for (size_t vi = 0; vi < sub_mesh->num_vertices; vi++) {
				ufbx_vec3 skin_pos = sub_mesh->skinned_position.values.data[vi];

				if (level == 0) {
					ufbx_matrix mat = ufbx_get_skin_vertex_matrix(skin, vi, NULL);
					ufbx_vec3 local_pos = sub_mesh->vertices.data[vi];
					ufbx_vec3 world_pos = ufbx_transform_position(&mat, local_pos);
					ufbxt_assert_close_vec3(err, world_pos, skin_pos);

					ufbx_skin_vertex skin_vertex = skin->vertices.data[vi];
					ufbx_matrix sum = { 0 };
					for (size_t wi = 0; wi < skin_vertex.num_weights; wi++) {
						ufbx_skin_weight weight = skin->weights.data[skin_vertex.weight_begin + wi];
						ufbx_skin_cluster *cluster = skin->clusters.data[weight.cluster_index];
						for (size_t i = 0; i < 12; i++) {
							sum.v[i] += cluster->geometry_to_world.v[i] * weight.weight;
						}
					}

					ufbx_vec3 manual_pos = ufbx_transform_position(&sum, local_pos);
					ufbxt_assert_close_vec3(err, manual_pos, skin_pos);
				} else {
					ufbx_subdivision_result *sub_res = sub_mesh->subdivision_result;
					ufbxt_assert(sub_res);

					ufbx_subdivision_weight_range range = sub_res->skin_cluster_ranges.data[vi];
					ufbx_matrix sum = { 0 };
					for (size_t wi = 0; wi < range.num_weights; wi++) {
						ufbx_subdivision_weight weight = sub_res->skin_cluster_weights.data[range.weight_begin + wi];
						ufbx_skin_cluster *cluster = skin->clusters.data[weight.index];
						for (size_t i = 0; i < 12; i++) {
							sum.v[i] += cluster->geometry_to_world.v[i] * weight.weight;
						}
					}

					ufbx_vec3 local_pos = sub_mesh->vertices.data[vi];
					ufbx_vec3 manual_pos = ufbx_transform_position(&sum, local_pos);
					ufbxt_assert_close_vec3_threshold(err, manual_pos, skin_pos, 10.0f);
					manual_pos = manual_pos;
				}
			}

			ufbx_free_mesh(sub_mesh);
		}

		ufbx_free_scene(state);
	}
}
#endif
