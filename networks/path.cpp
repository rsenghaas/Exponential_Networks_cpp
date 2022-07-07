#include "path.hpp"

auto integrator_thread(std::shared_ptr<SW_curve> pCurve, const state_type &v0, std::promise<std::vector<state_type>> v, std::promise<std::vector<double>> masses) -> void {
    spdlog::debug("Enter Integrator_thread");
    std::shared_ptr<SW_curve> thread_curve = pCurve;
    ODE_differential diff(thread_curve);
    ODE_integrator integrator(diff, v0);
    integrator.integrate_ode();
    v.set_value(integrator.get_v());
    masses.set_value(integrator.get_masses());
    spdlog::debug("Exit Integrator_thread");
}

auto Path::get_endpoint() -> state_type {
    return v_.back();
}

auto Path::append_data(std::vector<state_type>& v, std::vector<double>& masses) -> void {
    spdlog::debug("Path {} was integrated with {} steps.", id_, v.size());
    double current_mass = masses_.back();
    std::transform(masses.begin(), masses.end(), masses.begin(), [&](double m){return m + current_mass;});
    masses_.insert(masses_.end(), ++masses.begin(), masses.end());
    v_.insert(v_.end(), ++v.begin(), v.end());
}

auto Path::update(future_type& future) -> void {
    try {
        std::vector<state_type> future_v = std::get<0>(future).get();
        std::vector<double> future_m = std::get<1>(future).get();
        append_data(future_v, future_m);
    } catch(const std::exception& e) {
        spdlog::debug("Future caught exception {} with validity v: {}, masses: {}", e.what(), std::get<0>(future).valid(), std::get<1>(future).valid());
    }
}

auto Path::integrate(std::shared_ptr<SW_curve> pCurve) -> future_type {
    state_type v0 = get_endpoint();
    std::promise<std::vector<state_type>> prom_v_;
    std::promise<std::vector<double>> prom_masses_;
    std::future<std::vector<state_type>> future_v_ = prom_v_.get_future();
    std::future<std::vector<double>> future_masses_ = prom_masses_.get_future();

    std::thread{integrator_thread, pCurve, v0, std::move(prom_v_), std::move(prom_masses_)}.detach();
    future_type future;
    std::get<0>(future) = std::move(future_v_);
    std::get<1>(future) = std::move(future_masses_);
    return future;
}

auto Path::print_data() -> void {
    for (uint32_t i = 0; i < v_.size(); i++)
    {
        std::cout << "Step: " << i 
            << "x: " << complex_to_string(v_.at(i).at(kIndexX) )
            << " y1: " << complex_to_string(v_.at(i).at(kIndexY1)) 
            << " y2: " << complex_to_string(v_.at(i).at(kIndexY2)) 
            << " with mass " << masses_.at(i) << std::endl;
    }
}

auto Path::save_data() -> void {
    std::fstream data_file;
    std::string file_name = fmt::format("data/path_data/path_data_{}.csv", id_);
    data_file.open(file_name, std::ios::out);
    if (!data_file) {
        spdlog::debug("{} could not be created.", file_name);
    } else {
        for (const auto& zipped : boost::combine(v_, masses_)) {
            state_type line_v;
            double line_mass;
            boost::tie(line_v, line_mass) = zipped;
            std::string output_line = fmt::format("{},{},{},{}\n", 
                complex_to_string(line_v.at(kIndexX)), 
                complex_to_string(line_v.at(kIndexY1)), 
                complex_to_string(line_v.at(kIndexY2)), 
                line_mass);
            data_file << output_line;
        }
        spdlog::debug("Date saved to {}.", file_name);
        data_file.close();
    }
}
