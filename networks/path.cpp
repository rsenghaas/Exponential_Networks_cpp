#include "path.hpp"

auto get_path_point(state_type v) -> path_point {
  path_point pp;
  auto coord_real = static_cast<int32_t>(
      (v.at(kIndexX).real() - kMapRangeReal.at(kIndexLowerBound)) /
      (kMapRangeReal.at(kIndexUpperBound) -
       kMapRangeReal.at(kIndexLowerBound)) *
      kMapResolutionReal);
  auto coord_imag = static_cast<int32_t>(
      (v.at(kIndexX).imag() - kMapRangeImag.at(kIndexLowerBound)) /
      (kMapRangeImag.at(kIndexUpperBound) -
       kMapRangeImag.at(kIndexLowerBound)) *
      kMapResolutionImag);
  pp.coordinate_real = coord_real;
  pp.coordinate_imag = coord_imag;
  return pp;
}

auto integrator_thread(std::shared_ptr<SW_curve> pCurve, const state_type& v0,
                       std::promise<std::vector<state_type>> v,
                       std::promise<std::vector<double>> masses, double theta,
                       double cutoff) -> void {
  spdlog::debug("Enter Integrator_thread");
  ODE_differential diff(std::move(pCurve), theta);
  ODE_integrator integrator(diff, v0);
  integrator.integrate_ode(cutoff);
  v.set_value(integrator.get_v());
  masses.set_value(integrator.get_masses());
  spdlog::debug("Exit Integrator_thread");
}

auto Path::get_endpoint() -> state_type { return v_.back(); }

auto Path::append_data(std::vector<state_type>& v, std::vector<double>& masses)
    -> void {
  spdlog::debug("Path {} was integrated with {} steps.", path_id_, v.size());
  double current_mass = masses_.back();
  std::transform(masses.begin(), masses.end(), masses.begin(),
                 [&](double m) { return m + current_mass; });
  masses_.insert(masses_.end(), ++masses.begin(), masses.end());
  v_.insert(v_.end(), ++v.begin(), v.end());
}

auto Path::update(future_type& future) -> void {
  try {
    std::vector<state_type> future_v = std::get<0>(future).get();
    std::vector<double> future_m = std::get<1>(future).get();
    append_data(future_v, future_m);
  } catch (const std::exception& e) {
    spdlog::debug("Future caught exception {} with validity v: {}, masses: {}",
                  e.what(), std::get<0>(future).valid(),
                  std::get<1>(future).valid());
  }
}

auto Path::integrate(std::shared_ptr<SW_curve> pCurve, double theta,
                     double cutoff) -> future_type {
  state_type v0 = get_endpoint();
  std::promise<std::vector<state_type>> prom_v;
  std::promise<std::vector<double>> prom_masses;
  std::future<std::vector<state_type>> future_v = prom_v.get_future();
  std::future<std::vector<double>> future_masses = prom_masses.get_future();

  std::thread{integrator_thread,      pCurve, v0,    std::move(prom_v),
              std::move(prom_masses), theta,  cutoff}
      .detach();
  future_type future;
  std::get<0>(future) = std::move(future_v);
  std::get<1>(future) = std::move(future_masses);
  return future;
}

auto draw_line(path_point pp1, path_point pp2) -> std::vector<path_point> {
  std::vector<path_point> line;
  path_point next_pp;
  next_pp.id = pp1.id;
  next_pp.pp_vec_index = pp2.pp_vec_index;
  next_pp.coordinate_real = pp1.coordinate_real;
  next_pp.coordinate_imag = pp1.coordinate_imag;
  next_pp.t = {pp1.t.at(1), pp2.t.at(0)};
  int32_t start_coord_real = pp1.coordinate_real;
  int32_t start_coord_imag = pp1.coordinate_imag;
  int32_t end_coord_real = pp2.coordinate_real;
  int32_t end_coord_imag = pp2.coordinate_imag;

  // spdlog::debug("Enter connecting Pixel ({},{}) and ({}, {}).",
  // pp1.coordinate_real, pp1.coordinate_imag, pp2.coordinate_real,
  // pp2.coordinate_imag);
  double ghost_coord_real{0};
  double ghost_coord_imag{0};

  auto distance = static_cast<int32_t>(
      kLineStepsPerUnit *
      (std::pow(std::pow(end_coord_real - start_coord_real, 2) +
                    std::pow(end_coord_imag - start_coord_imag, 2),
                1.0 / 2)));
  if (distance == 0) {
    return line;
  }
  for (int32_t counter = 0; counter <= distance; counter++) {
    ghost_coord_real = (1 - counter * 1.0 / distance) * (start_coord_real) +
                       counter * 1.0 / distance * (end_coord_real) +
                       1.0 / 2;  // Draw line from middle of the pixel.
    ghost_coord_imag = (1 - counter * 1.0 / distance) * (start_coord_imag) +
                       counter * 1.0 / distance * (end_coord_imag) +
                       1.0 / 2;  //

    if (static_cast<int32_t>(ghost_coord_real + kNumOffset) !=
            next_pp.coordinate_real ||
        static_cast<int32_t>(ghost_coord_imag + kNumOffset) !=
            next_pp.coordinate_imag) {
      next_pp.coordinate_real =
          static_cast<int32_t>(ghost_coord_real + kNumOffset);
      next_pp.coordinate_imag =
          static_cast<int32_t>(ghost_coord_imag + kNumOffset);
      line.push_back(next_pp);
    }
  }
  // spdlog::debug("Exit Connecting pixel");
  return line;
}

auto Path::compute_map_points(std::vector<path_point>& pp_vec,
                              std::vector<uint32_t>& index_vec) -> void {
  path_point start_pp = get_path_point(v_.front());
  start_pp.t = {0, 0};
  start_pp.id = path_id_;
  start_pp.pp_vec_index = 0;
  pp_vec.push_back(start_pp);
  index_vec.push_back(0);

  path_point new_pp;
  for (uint32_t i = 1; i < v_.size(); ++i) {
    new_pp = get_path_point(v_.at(i));
    new_pp.t = {i, i};
    new_pp.id = path_id_;
    new_pp.pp_vec_index = pp_vec.size();
    if (new_pp.coordinate_real != pp_vec.back().coordinate_real ||
        new_pp.coordinate_imag != pp_vec.back().coordinate_imag) {
      auto line = draw_line(pp_vec.back(), new_pp);
      pp_vec.insert(std::end(pp_vec), std::begin(line), std::end(line));
      index_vec.push_back(pp_vec.size() - 1);
    } else {
      pp_vec.back().t.at(kIndexEndTime)++;
    }
  }
  spdlog::debug("Path {} is represented by {} points on the map.", path_id_,
                pp_vec.size());
}

auto Path::print_data() -> void {
  for (uint32_t i = 0; i < v_.size(); i++) {
    std::cout << "Step: " << i
              << "x: " << complex_to_string(v_.at(i).at(kIndexX))
              << " y1: " << complex_to_string(v_.at(i).at(kIndexY1))
              << " y2: " << complex_to_string(v_.at(i).at(kIndexY2))
              << " with mass " << masses_.at(i) << std::endl;
  }
}

auto Path::save_data() -> void {
  std::fstream data_file;
  std::string filename =
      fmt::format("data/path_data/path_data_{}.csv", path_id_);
  data_file.open(filename, std::ios::out);
  if (!data_file) {
    spdlog::debug("{} could not be created.", filename);
  } else {
    for (const auto& zipped : boost::combine(v_, masses_)) {
      state_type line_v;
      double line_mass;
      boost::tie(line_v, line_mass) = zipped;
      std::string output_line =
          fmt::format("{},{},{},{}\n", complex_to_string(line_v.at(kIndexX)),
                      complex_to_string(line_v.at(kIndexY1)),
                      complex_to_string(line_v.at(kIndexY2)), line_mass);
      data_file << output_line;
    }
    spdlog::debug("Data saved to {}.", filename);
    data_file.close();
  }
}

auto Path::get_point(uint32_t t) -> state_type { return v_.at(t); }

auto Path::truncate(uint32_t t_start, uint32_t t_end) -> void {
  if (t_end <= t_start) {
    spdlog::debug("End time greater than start time. Skip truncation.");
    return;
  }
  spdlog::debug("Erasing");
  spdlog::debug("Id: {}, v size: {}", path_id_, v_.size());
  print_state_type(v_.at(t_start));
  v_.erase(v_.begin() + t_end, v_.end());
  spdlog::debug("{}", v_.size());
  masses_.erase(masses_.begin() + t_end, masses_.end());
  v_.erase(v_.begin(), v_.begin() + t_start);
  masses_.erase(masses_.begin(), masses_.begin() + t_start);
  double m0 = masses_.front();
  for (auto& mass : masses_) {
    mass -= m0;
  }
  print_state_type(v_.at(0));
  spdlog::debug("Done erasing.");
}

auto Path::add_single_point(state_type pt) -> void {
  state_type dv = pt + (-1 * v_.back());
  v_.push_back(pt);
  double dm = compute_dm(pt, dv);
  masses_.push_back(masses_.back() + dm);
}
