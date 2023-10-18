#define The_4_2_PARTITION
auto ADHM::custom_BPS() -> void {
  auto path_it = get_iterator_by_id(new_paths_, 0);
  overwrite_path(path_it, cutoffPoint);
  evolve_path(path_it, 4 * kD4Cutoff);
  // save_data(0);

  path_it = get_iterator_by_id(new_paths_, 1);
  evolve_path(path_it, kCutoff);

  self_intersection_handler(1, true, -1, 0, false, false);

  path_it = get_iterator_by_id(new_paths_, 1);
  save_data(1);

  path_it = get_iterator_by_id(new_paths_, 2);
  evolve_path(path_it, 1.0 / 4 * kCutoff);

  path_it = get_iterator_by_id(new_paths_, 3);
  evolve_path(path_it, 1.0 / 4 * kCutoff);

  two_path_intersection_handler(2, 3, true, true, 0, 0, true, false);
  save_data(3);

  path_it = get_iterator_by_id(new_paths_, 4);
  auto endpoint = path_it->get_endpoint();
  invert_state(endpoint);
  add_new_path(endpoint);
  evolve_path(path_it, kCutoff);
  two_path_intersection_handler(1, 4, false, true, -1, 0, false, false);
  save_data(4);

  path_it = get_iterator_by_id(new_paths_, 5);
  evolve_path(path_it, kCutoff);

  path_it = get_iterator_by_id(new_paths_, 6);
  evolve_path(path_it, kCutoff);

  two_path_intersection_handler(5, 6, false, true, -1, 0, true, false);
  save_data(6);

  path_it = get_iterator_by_id(new_paths_, 6);
  endpoint = path_it->get_endpoint();

  path_it = get_iterator_by_id(new_paths_, 7);
  overwrite_path(path_it, endpoint);

  evolve_path(path_it, kCutoff);
  two_path_intersection_handler(2, 7, false, true, 0, 0, true, false);
  save_data(7);

  path_it = get_iterator_by_id(new_paths_, 8);
  evolve_path(path_it, kCutoff);
  two_path_intersection_handler(1, 8, false, true, -1, 0, false, false);
  save_data(8);

  path_it = get_iterator_by_id(new_paths_, 9);
  evolve_path(path_it, kCutoff);
  save_data(2);

  two_path_intersection_handler(5, 9, true, true, -1, 0, true, false);

  save_data(5);
  save_data(8);

  path_it = get_iterator_by_id(new_paths_, 10);
  evolve_path(path_it, kCutoff);
  two_path_intersection_handler(2, 10, false, true, 0, 0, true, false);
  save_data(9);
  save_data(10);

  path_it = get_iterator_by_id(new_paths_, 11);
  evolve_path(path_it, kCutoff);
  save_data(11);
}

#define THE_3_1_PARTITION

auto ADHM::custom_BPS() -> void {
  auto path_it = get_iterator_by_id(new_paths_, 0);
  overwrite_path(path_it, cutoffPoint);
  evolve_path(path_it, 4 * kD4Cutoff);

  path_it = get_iterator_by_id(new_paths_, 1);
  evolve_path(path_it, kCutoff);

  self_intersection_handler(1, true, -1, 0, false, false);

  path_it = get_iterator_by_id(new_paths_, 1);
  save_data(1);

  path_it = get_iterator_by_id(new_paths_, 2);
  evolve_path(path_it, 1.0 / 4 * kCutoff);

  path_it = get_iterator_by_id(new_paths_, 3);
  evolve_path(path_it, 1.0 / 4 * kCutoff);

  two_path_intersection_handler(2, 3, true, true, 0, 0, true, false);
  save_data(2);
  save_data(3);

  path_it = get_iterator_by_id(new_paths_, 4);
  auto endpoint = path_it->get_endpoint();
  invert_state(endpoint);
  add_new_path(endpoint);
  evolve_path(path_it, kCutoff);
  two_path_intersection_handler(1, 4, false, true, -1, 0, false, false);
  save_data(4);

  path_it = get_iterator_by_id(new_paths_, 5);
  evolve_path(path_it, kCutoff);

  path_it = get_iterator_by_id(new_paths_, 6);
  evolve_path(path_it, kCutoff);

  two_path_intersection_handler(6, 5, true, true, -1, 0, false, false);
  save_data(6);
  save_data(5);

  path_it = get_iterator_by_id(new_paths_, 7);
  endpoint = path_it->get_endpoint();
  endpoint.at(kIndexY1) -= 2 * pi * J;
  // invert_state(endpoint);

  path_it = get_iterator_by_id(new_paths_, 7);
  overwrite_path(path_it, endpoint);

  evolve_path(path_it, kCutoff);
  two_path_intersection_handler(2, 7, false, true, 0, 0, true, false);
  save_data(7);

  path_it = get_iterator_by_id(new_paths_, 8);
  evolve_path(path_it, kCutoff);
  // two_path_intersection_handler(1,8,false,true, -1, 0, false, false);
  save_data(8);
}

#define THE_4_1_PARTITIOM

auto ADHM::custom_BPS() -> void {
  auto path_it = get_iterator_by_id(new_paths_, 0);
  overwrite_path(path_it, cutoffPoint);
  evolve_path(path_it, 4 * kD4Cutoff);

  path_it = get_iterator_by_id(new_paths_, 1);
  evolve_path(path_it, kCutoff);

  self_intersection_handler(1, true, -1, 0, false, false);

  path_it = get_iterator_by_id(new_paths_, 1);
  save_data(1);

  path_it = get_iterator_by_id(new_paths_, 2);
  evolve_path(path_it, 1.0 / 4 * kCutoff);

  path_it = get_iterator_by_id(new_paths_, 3);
  evolve_path(path_it, 1.0 / 4 * kCutoff);

  two_path_intersection_handler(2, 3, true, true, 0, 0, true, false);
  save_data(2);
  save_data(3);

  path_it = get_iterator_by_id(new_paths_, 4);
  evolve_path(path_it, kCutoff);
  two_path_intersection_handler(1, 4, false, true, -1, 0, false, false);
  save_data(4);

  path_it = get_iterator_by_id(new_paths_, 5);
  evolve_path(path_it, kCutoff);
  two_path_intersection_handler(2, 5, false, true, 0, 0, true, false);
  save_data(5);

  path_it = get_iterator_by_id(new_paths_, 6);
  auto endpoint = path_it->get_endpoint();
  invert_state(endpoint);
  add_new_path(endpoint);
  evolve_path(path_it, kCutoff);
  two_path_intersection_handler(1, 6, false, true, -1, 0, false, false);
  save_data(6);

  path_it = get_iterator_by_id(new_paths_, 7);
  evolve_path(path_it, kCutoff);

  path_it = get_iterator_by_id(new_paths_, 8);
  evolve_path(path_it, kCutoff);

  two_path_intersection_handler(8, 7, true, true, -1, 0, false, false);
  save_data(8);
  save_data(7);

  path_it = get_iterator_by_id(new_paths_, 9);
  endpoint = path_it->get_endpoint();
  endpoint.at(kIndexY1) -= 2 * pi * J;
  // invert_state(endpoint);

  overwrite_path(path_it, endpoint);

  evolve_path(path_it, kCutoff);
  two_path_intersection_handler(2, 9, false, true, 0, 0, true, false);
  save_data(9);

  path_it = get_iterator_by_id(new_paths_, 10);
  evolve_path(path_it, kCutoff);
  // two_path_intersection_handler(1,8,false,true, -1, 0, false, false);
  save_data(10);
}

#define THE_DEFORMED_D0
auto ADHM::custom_BPS() -> void {
  auto path_it = get_iterator_by_id(new_paths_, 1);
  evolve_path(path_it, 0.1 * kCutoff);
  path_it->truncate(0, 1000);

  path_it = get_iterator_by_id(new_paths_, 2);
  // overwrite_path(path_it, cutoffPoint);
  evolve_path(path_it, 0.1 * kCutoff);
  auto endpoint = path_it->get_point(880);
  add_new_path(endpoint);
  path_it = get_iterator_by_id(new_paths_, 3);
  endpoint.at(kIndexY1) -= 2 * pi * J;
  overwrite_path(path_it, endpoint);
  evolve_path(path_it, 0.1 * kCutoff);
  // self_intersection_handler(1, true, -1, 1, false, false);
  two_path_intersection_handler(1, 3, true, true, -1, 0, false, false);
  path_it = get_iterator_by_id(new_paths_, 4);

  evolve_path(path_it, 0.1 * kCutoff);
  two_path_intersection_handler(2, 4, true, true, -1, 0, false, false);
  path_it = get_iterator_by_id(new_paths_, 2);
  save_data(1);
  save_data(2);
  save_data(3);
  save_data(4);
}

#define THE_4_2_2_PARTITION
auto ADHM::custom_BPS() -> void {
  spdlog::debug("Drawing state.");
  auto path_it = get_iterator_by_id(new_paths_, 0);
  overwrite_path(path_it, cutoffPoint);
  evolve_path(path_it, 4 * kD4Cutoff);

  path_it = get_iterator_by_id(new_paths_, 1);
  evolve_path(path_it, kCutoff);

  self_intersection_handler(1, true, -1, 0, false, false);

  path_it = get_iterator_by_id(new_paths_, 1);
  save_data(1);

  path_it = get_iterator_by_id(new_paths_, 2);
  evolve_path(path_it, 1.0 / 4 * kCutoff);

  path_it = get_iterator_by_id(new_paths_, 3);
  evolve_path(path_it, 1.0 / 4 * kCutoff);

  two_path_intersection_handler(2, 3, false, false, 0, 0, false, false);
  two_path_intersection_handler(2, 3, true, true, 0, 0, false, false);
  save_data(2);
  save_data(3);

  uint32_t current_index = 4;
  path_it = get_iterator_by_id(new_paths_, current_index);
  evolve_path(path_it, 1.0 / 4 * kCutoff);
  current_index++;
  path_it = get_iterator_by_id(new_paths_, current_index);
  auto endpoint = path_it->get_endpoint();
  endpoint.at(kIndexY1) -= 2 * pi * J;
  invert_state(endpoint);
  overwrite_path(path_it, endpoint);
  evolve_path(path_it, kCutoff);
  endpoint = path_it->get_point(17);
  endpoint.at(kIndexX) =
      0.5 * endpoint.at(kIndexX) + 0.5 * path_it->get_point(18).at(kIndexX);
  curve_->match_fiber(endpoint);
  invert_state(endpoint);
  endpoint.at(kIndexY1) += 2 * pi * J;
  add_new_path(endpoint);
  current_index++;
  path_it = get_iterator_by_id(new_paths_, current_index);
  evolve_path(path_it, kCutoff);
  two_path_intersection_handler(1, current_index, false, true, -1, 0, false,
                                false);
  save_data(current_index);
  current_index++;
  path_it = get_iterator_by_id(new_paths_, current_index);
  evolve_path(path_it, kCutoff);
  two_path_intersection_handler(5, current_index, true, true, 0, 0, false,
                                false);
  save_data(5);
  save_data(current_index);
  two_path_intersection_handler(1, 4, false, true, -2, 0, true, true);
  save_data(4);
  current_index += 2;
  path_it = get_iterator_by_id(new_paths_, current_index);
  evolve_path(path_it, kCutoff);
  two_path_intersection_handler(2, current_index, false, true, 0, 0, true,
                                false);
  save_data(current_index);
  current_index++;
  path_it = get_iterator_by_id(new_paths_, current_index);
  evolve_path(path_it, kCutoff);
  two_path_intersection_handler(1, current_index, false, true, -1, 0, true,
                                true);
  save_data(current_index);
  current_index++;
  path_it = get_iterator_by_id(new_paths_, current_index);
  evolve_path(path_it, kCutoff);
  two_path_intersection_handler(2, current_index, false, true, 0, 0, true,
                                false);
  save_data(current_index);
  current_index++;
  path_it = get_iterator_by_id(new_paths_, current_index);
  evolve_path(path_it, kCutoff);
  save_data(current_index);
  two_path_intersection_handler(0, 2, true, false, 0, 0, false, false);
  save_data(0);
}

auto ADHM::custom_BPS() -> void {
  spdlog::debug("Drawing state.");
  auto path_it = get_iterator_by_id(new_paths_, 0);
  overwrite_path(path_it, cutoffPoint);
  evolve_path(path_it, 4 * kD4Cutoff);

  path_it = get_iterator_by_id(new_paths_, 1);
  evolve_path(path_it, kCutoff);

  self_intersection_handler(1, true, -1, 0, false, false);

  path_it = get_iterator_by_id(new_paths_, 1);
  save_data(1);

  path_it = get_iterator_by_id(new_paths_, 2);
  evolve_path(path_it, 1.0 / 4 * kCutoff);

  path_it = get_iterator_by_id(new_paths_, 3);
  evolve_path(path_it, 1.0 / 4 * kCutoff);

  two_path_intersection_handler(2, 3, false, false, 0, 0, false, false);
  two_path_intersection_handler(2, 3, true, true, 0, 0, false, false);
  save_data(2);
  save_data(3);

  uint32_t current_index = 4;
  path_it = get_iterator_by_id(new_paths_, current_index);
  evolve_path(path_it, 1.0 / 4 * kCutoff);
  current_index++;
  path_it = get_iterator_by_id(new_paths_, current_index);
  auto endpoint = path_it->get_endpoint();
  endpoint.at(kIndexY1) -= 2 * pi * J;
  invert_state(endpoint);
  overwrite_path(path_it, endpoint);
  evolve_path(path_it, kCutoff);

  // Newton method by hand...
  endpoint = path_it->get_point(17);
  endpoint.at(kIndexX) =
      0.5 * endpoint.at(kIndexX) + 0.5 * path_it->get_point(18).at(kIndexX);

  curve_->match_fiber(endpoint);
  invert_state(endpoint);
  endpoint.at(kIndexY1) += 2 * pi * J;
  add_new_path(endpoint);
  current_index++;
  path_it = get_iterator_by_id(new_paths_, current_index);

  // New Baseline!
  endpoint = path_it->get_endpoint();
  invert_state(endpoint);

  evolve_path(path_it, kCutoff);
  two_path_intersection_handler(1, current_index, false, true, -1, 0, false,
                                false);
  save_data(current_index);
  current_index++;

  path_it = get_iterator_by_id(new_paths_, current_index);
  evolve_path(path_it, kCutoff);
  two_path_intersection_handler(5, current_index, true, true, 0, 0, false,
                                false);
  save_data(5);
  save_data(current_index);
  current_index++;

  auto baseline = current_index;
  path_it = get_iterator_by_id(new_paths_, baseline);
  overwrite_path(path_it, endpoint);
  evolve_path(path_it, kCutoff);
  save_data(baseline);
  two_path_intersection_handler(1, 4, false, true, -1, 0, false, false);
  save_data(4);
  current_index++;
  path_it = get_iterator_by_id(new_paths_, current_index);
  evolve_path(path_it, kCutoff);
  // save_data(current_index);
  two_path_intersection_handler(baseline, current_index, false, true, -2, 0,
                                true, false);
  save_data(current_index);
  current_index++;
  path_it = get_iterator_by_id(new_paths_, current_index);
  evolve_path(path_it, kCutoff);
  save_data(current_index);
  two_path_intersection_handler(2, current_index, false, true, 0, 0, true,
                                false);
  save_data(current_index);
  current_index++;
  path_it = get_iterator_by_id(new_paths_, current_index);
  spdlog::debug("Current index: {}", current_index);
  evolve_path(path_it, kCutoff);
  save_data(current_index);
  two_path_intersection_handler(1, current_index, false, true, -1, 0, false,
                                false);
  save_data(current_index);
  current_index++;
  path_it = get_iterator_by_id(new_paths_, current_index);
  evolve_path(path_it, kCutoff);
  two_path_intersection_handler(baseline, current_index, true, true, -1, 0,
                                true, false);
  save_data(current_index);
  save_data(baseline);
  current_index++;
  path_it = get_iterator_by_id(new_paths_, current_index);
  evolve_path(path_it, kCutoff);
  two_path_intersection_handler(2, current_index, false, true, 0, 0, true,
                                false);
  save_data(current_index);
  current_index++;
  path_it = get_iterator_by_id(new_paths_, current_index);
  evolve_path(path_it, kCutoff);
  save_data(current_index);
}

#define THE_2_PARTITION_GENERIC_STEP_1_2
spdlog::debug("Drawing state.");
auto path_it = get_iterator_by_id(new_paths_, 0);
overwrite_path(path_it, cutoffPoint);
evolve_path(path_it, 4 * kD4Cutoff);
// save_data(0);

path_it = get_iterator_by_id(new_paths_, 1);
evolve_path(path_it, kCutoff);

// self_intersection_handler(1, true, -1, 0, false, false);

path_it = get_iterator_by_id(new_paths_, 1);
// save_data(1);

const cplx x_coord = -1.4805;
auto fiber = curve_->get_fiber(x_coord);
std::transform(fiber.cbegin(), fiber.cend(), fiber.begin(),
               [](cplx y_coord) { return normalize(std::log(y_coord)); });
state_type state_1 = {x_coord, fiber.at(0), fiber.at(1)};
state_type state_2 = {x_coord, fiber.at(1), fiber.at(0)};
curve_->match_fiber(state_1);
curve_->match_fiber(state_2);

path_it = get_iterator_by_id(new_paths_, 2);
overwrite_path(path_it, state_1);
evolve_path(path_it, 1.0 / 4 * kCutoff);

// path_it = get_iterator_by_id(new_paths_, 3);
// evolve_path(path_it, 1.0/4 * kCutoff);

// two_path_intersection_handler(2, 3, false, false, 0, 0, false, false);
// two_path_intersection_handler(2,3,true, true, 0,0,false,false);
// save_data(3);
uint32_t current_index = 3;
add_new_path(state_2);
path_it = get_iterator_by_id(new_paths_, current_index);
evolve_path(path_it, 1.0 / 4 * kCutoff);
save_data(current_index);
two_path_intersection_handler(1, current_index, true, true, -1, 0, false,
                              false);
save_data(1);
save_data(current_index);
current_index++;
path_it = get_iterator_by_id(new_paths_, current_index);
evolve_path(path_it, 1.0 / 4 * kCutoff);
two_path_intersection_handler(2, current_index, true, true, 0, 0, true, false);
save_data(2);
save_data(current_index);
current_index++;
path_it = get_iterator_by_id(new_paths_, current_index);
evolve_path(path_it, 1.0 / 4 * kCutoff);
save_data(current_index);
// two_path_intersection_handler(current_index-1, current_index, true, false, 0,
// 0, true, false); save_data(current_index - 1);

#define FULL_N_EQ_2_NETWROK_WITH_CONES
auto ADHM::custom_BPS() -> void {
  spdlog::debug("Drawing state.");

  auto path_it = get_iterator_by_id(new_paths_, 0);
  overwrite_path(path_it, cutoffPoint);
  evolve_path(path_it, 4 * kD4Cutoff);
  // path_it->truncate(0,1480);
  // save_data(0);
  // save_data(0);

  path_it = get_iterator_by_id(new_paths_, 1);
  evolve_path(path_it, kCutoff);
  // path_it->truncate(0, 1540);
  self_intersection_handler(1, false, -1, 0, false, false);
  self_intersection_handler(1, true, -1, 0, true, true);
  save_data(1);

  path_it = get_iterator_by_id(new_paths_, 2);
  evolve_path(path_it, kCutoff);
  // self_intersection_handler(2, true, 1, 0, false, false);
  // save_data(2);

  // self_intersection_handler(1, true, -1, 0, false, false);

  path_it = get_iterator_by_id(new_paths_, 3);
  evolve_path(path_it, kCutoff);
  two_path_intersection_handler(2, 3, true, true, 0, 0, true, false);
  path_it = get_iterator_by_id(new_paths_, 4);
  evolve_path(path_it, kCutoff);

  path_it = get_iterator_by_id(new_paths_, 5);
  evolve_path(path_it, kCutoff);
  two_path_intersection_handler(1, 5, false, true, -1, 0, false, false);
  path_it = get_iterator_by_id(new_paths_, 6);
  evolve_path(path_it, kCutoff);
  two_path_intersection_handler(2, 6, false, true, 0, 0, true, false);
  path_it = get_iterator_by_id(new_paths_, 7);
  evolve_path(path_it, kCutoff);
  two_path_intersection_handler(2, 4, false, true, 0, 0, true, false);
  save_data(2);
  save_data(3);
  save_data(4);
  save_data(5);
  save_data(6);
  save_data(7);

  path_it = get_iterator_by_id(new_paths_, 1);
  add_new_path(path_it->get_endpoint());
  path_it = get_iterator_by_id(new_paths_, 9);
  evolve_path(path_it, 1);
  path_it->truncate(0, 30);
  save_data(9);

  path_it = get_iterator_by_id(new_paths_, 5);
  add_new_path(path_it->get_endpoint());
  path_it = get_iterator_by_id(new_paths_, 10);
  evolve_path(path_it, 1);
  path_it->truncate(0, 30);
  save_data(10);

  path_it = get_iterator_by_id(new_paths_, 5);
  auto startpoint = path_it->get_point(0);
  invert_state(startpoint);
  add_new_path(startpoint);
  path_it = get_iterator_by_id(new_paths_, 11);
  evolve_path(path_it, 1);
  path_it->truncate(0, 30);
  save_data(11);

  path_it = get_iterator_by_id(new_paths_, 7);
  startpoint = path_it->get_point(0);
  invert_state(startpoint);
  add_new_path(startpoint);
  path_it = get_iterator_by_id(new_paths_, 12);
  evolve_path(path_it, 1);
  path_it->truncate(0, 30);
  save_data(12);

  path_it = get_iterator_by_id(new_paths_, 2);
  add_new_path(path_it->get_endpoint());
  path_it = get_iterator_by_id(new_paths_, 13);
  evolve_path(path_it, 1);
  path_it->truncate(0, 30);
  save_data(13);
}

#define DEFORMED_2_1
auto ADHM::custom_BPS() -> void {
  spdlog::debug("Drawing state.");
  auto path_it = get_iterator_by_id(new_paths_, 0);
  overwrite_path(path_it, cutoffPoint);
  evolve_path(path_it, 4 * kD4Cutoff);
  path_it->truncate(0, 2270);
  save_data(0);

  path_it = get_iterator_by_id(new_paths_, 1);
  evolve_path(path_it, kCutoff);
  path_it->truncate(0, 1950);
  save_data(1);

  auto endpoint = path_it->get_endpoint();
  auto startpoint_3 = endpoint;
  // invert_state(startpoint_3);
  startpoint_3.at(kIndexY2) -= 2 * pi * J;
  auto startpoint_4 = endpoint;
  // invert_state(startpoint_4);
  startpoint_4.at(kIndexY2) = startpoint_4.at(kIndexY1) + 2 * pi * J;

  add_new_path(startpoint_3);
  add_new_path(startpoint_4);

  path_it = get_iterator_by_id(new_paths_, 3);
  evolve_path(path_it, kCutoff);
  path_it->truncate(0, 1000);
  save_data(3);

  // two_path_intersection_handler(2,3, true,true, 0,0,true,false);
  path_it = get_iterator_by_id(new_paths_, 4);
  evolve_path(path_it, kCutoff);
  path_it->truncate(0, 400);
  save_data(4);

  auto startpoint_5 = path_it->get_point(20);
  path_it = get_iterator_by_id(new_paths_, 1);
  auto fibers = path_it->get_endpoint();
  startpoint_5.at(kIndexY1) = fibers.at(kIndexY1);
  startpoint_5.at(kIndexY2) = fibers.at(kIndexY2);
  curve_->match_fiber(startpoint_5);
  invert_state(startpoint_5);
  auto startpoint_6 = startpoint_5;
  startpoint_6.at(kIndexY2) += 2 * pi * J;
  invert_state(startpoint_6);
  add_new_path(startpoint_5);
  add_new_path(startpoint_6);

  path_it = get_iterator_by_id(new_paths_, 5);
  evolve_path(path_it, kCutoff);
  path_it->truncate(0, 900);
  save_data(5);

  path_it = get_iterator_by_id(new_paths_, 6);
  evolve_path(path_it, kCutoff);
  path_it->truncate(0, 600);
  save_data(6);

  two_path_intersection_handler(4, 5, true, true, -1, 0, true, true);
  save_data(4);
  save_data(5);
  path_it = get_iterator_by_id(new_paths_, 5);
  endpoint = path_it->get_endpoint();
  endpoint.at(kIndexY2) += 4 * pi * J;
  // invert_state(endpoint);
  path_it = get_iterator_by_id(new_paths_, 7);
  overwrite_path(path_it, endpoint);
  evolve_path(path_it, kCutoff);
  path_it->truncate(0, 200);
  save_data(7);

  two_path_intersection_handler(6, 7, true, true, 1, 0, true, false);
  save_data(6);
  save_data(7);

  path_it = get_iterator_by_id(new_paths_, 8);
  endpoint = path_it->get_endpoint();
  endpoint.at(kIndexY2) = endpoint.at(kIndexY1) + 2 * pi * J;
  overwrite_path(path_it, endpoint);
  evolve_path(path_it, kCutoff);

  two_path_intersection_handler(3, 8, true, true, 0, 0, false, false);
  save_data(3);
  save_data(8);

  path_it = get_iterator_by_id(new_paths_, 9);
  evolve_path(path_it, kCutoff);
  save_data(9);
}