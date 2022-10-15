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
  evolve_path(path_it, 1.0/4 * kCutoff);
 
  path_it = get_iterator_by_id(new_paths_, 3);
  evolve_path(path_it, 1.0/4 * kCutoff);
  
  two_path_intersection_handler(2,3,true, true, 0,0,true,false);
  save_data(3);

  path_it = get_iterator_by_id(new_paths_, 4);
  auto endpoint = path_it->get_endpoint();
  invert_state(endpoint);
  add_new_path(endpoint);
  evolve_path(path_it, kCutoff);
  two_path_intersection_handler(1,4,false, true, -1, 0, false, false);
  save_data(4);
  
  path_it = get_iterator_by_id(new_paths_, 5);
  evolve_path(path_it, kCutoff);


  path_it = get_iterator_by_id(new_paths_, 6);
  evolve_path(path_it, kCutoff);

  two_path_intersection_handler(5, 6, false, true, -1, 0, true, false);
  save_data(6);

  path_it = get_iterator_by_id(new_paths_, 6);
  endpoint = path_it-> get_endpoint();

  path_it = get_iterator_by_id(new_paths_, 7);
  overwrite_path(path_it, endpoint);
    
  evolve_path(path_it, kCutoff);
  two_path_intersection_handler(2,7,false,true, 0, 0, true, false);
  save_data(7);
  

  path_it = get_iterator_by_id(new_paths_, 8);
  evolve_path(path_it, kCutoff);
  two_path_intersection_handler(1,8,false,true, -1, 0, false, false);
  save_data(8);

  path_it = get_iterator_by_id(new_paths_, 9);
  evolve_path(path_it, kCutoff);
  save_data(2);

  two_path_intersection_handler(5, 9, true, true, -1, 0, true, false);

  save_data(5);
  save_data(8);

  
  path_it =get_iterator_by_id(new_paths_, 10);
  evolve_path(path_it, kCutoff);
  two_path_intersection_handler(2, 10, false, true, 0,0, true, false);
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
  evolve_path(path_it, 1.0/4 * kCutoff);
 
  path_it = get_iterator_by_id(new_paths_, 3);
  evolve_path(path_it, 1.0/4 * kCutoff);
  
  two_path_intersection_handler(2,3,true, true, 0,0,true,false);
  save_data(2);
  save_data(3);

  path_it = get_iterator_by_id(new_paths_, 4);
  auto endpoint = path_it->get_endpoint();
  invert_state(endpoint);
  add_new_path(endpoint);
  evolve_path(path_it, kCutoff);
  two_path_intersection_handler(1,4,false, true, -1, 0, false, false);
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
  two_path_intersection_handler(2,7,false,true, 0, 0, true, false);
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
  evolve_path(path_it, 1.0/4 * kCutoff);
 
  path_it = get_iterator_by_id(new_paths_, 3);
  evolve_path(path_it, 1.0/4 * kCutoff);
  
  two_path_intersection_handler(2,3,true, true, 0,0,true,false);
  save_data(2);
  save_data(3);
  
  path_it = get_iterator_by_id(new_paths_, 4);
  evolve_path(path_it, kCutoff);
  two_path_intersection_handler(1,4,false, true, -1, 0, false, false);
  save_data(4);
  
  path_it = get_iterator_by_id(new_paths_, 5);
  evolve_path(path_it, kCutoff);
  two_path_intersection_handler(2,5,false, true, 0,0,true,false);
  save_data(5);

  path_it = get_iterator_by_id(new_paths_, 6);
  auto endpoint = path_it->get_endpoint();
  invert_state(endpoint);
  add_new_path(endpoint);
  evolve_path(path_it, kCutoff);
  two_path_intersection_handler(1,6,false, true, -1, 0, false, false);
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
  two_path_intersection_handler(2,9,false,true, 0, 0, true, false);
  save_data(9);
  

  path_it = get_iterator_by_id(new_paths_, 10);
  evolve_path(path_it, kCutoff);
  // two_path_intersection_handler(1,8,false,true, -1, 0, false, false);
  save_data(10);
}
  

