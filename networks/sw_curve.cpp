#include "sw_curve.hpp"

auto SW_curve::compute_derivatives() -> void {
    dH_dx_ = H_.diff(x_);
    dH_dy_ = H_.diff(y_);
    d2H_dy2_ = H_.diff(y_, 2);
}

auto SW_curve::eval_H(const cplx &x, const cplx &y) -> cplx {
    GiNaC::ex x_ex = complex_to_ex(x);
    GiNaC::ex y_ex = complex_to_ex(y);
    return numeric_to_complex(GiNaC::ex_to<GiNaC::numeric>(H_.subs(GiNaC::lst{x_ == x_ex, y_ == y_ex}).evalf()));
}

auto SW_curve::eval_dH_dx(const cplx &x, const cplx &y) -> cplx {
    GiNaC::ex x_ex = complex_to_ex(x);
    GiNaC::ex y_ex = complex_to_ex(y);
    return numeric_to_complex(GiNaC::ex_to<GiNaC::numeric>(dH_dx_.subs(GiNaC::lst{x_ == x_ex, y_ == y_ex}).evalf()));
}

auto SW_curve::eval_dH_dy(const cplx &x, const cplx &y) -> cplx {
    GiNaC::ex x_ex = complex_to_ex(x);
    GiNaC::ex y_ex = complex_to_ex(y);
    return numeric_to_complex(GiNaC::ex_to<GiNaC::numeric>(dH_dy_.subs(GiNaC::lst{x_ == x_ex, y_ == y_ex}).evalf()));
}

auto SW_curve::eval_d2H_dy2(const cplx &x, const cplx &y) -> cplx {
    GiNaC::ex x_ex = complex_to_ex(x);
    GiNaC::ex y_ex = complex_to_ex(y);
    return numeric_to_complex(GiNaC::ex_to<GiNaC::numeric>(d2H_dy2_.subs(GiNaC::lst{x_ == x_ex, y_ == y_ex}).evalf()));
}

auto SW_curve::get_branch_points() -> std::vector<cplx> {
    std::vector<cplx> dummy{22.0/7};
    GiNaC::ex disc = discriminant(H_, y_);
    
    // GiNaC::ex f(GiNaC::pow(x_, 3) * 3 + x_ - 1);
    std::vector<cplx> branch_points = roots(disc, x_);
    return branch_points;
}


auto SW_curve::get_ramification_points() -> std::vector<std::array<cplx, 2>> {
    std::vector<cplx> branch_points = get_branch_points();
    std::vector<std::array<cplx, 2>> ramification_points;
    for (auto& b : branch_points)
    {
        std::array<cplx, 2> r;
        r.at(kIndexX) = b;
        r.at(kIndexY) = get_branched_sheet(b);
        ramification_points.push_back(r);
    }
    return ramification_points;
}

// At the moment this only does exponential networks
// It might to return to expressions and build the differential outside the class.
auto SW_curve::sw_differential(const state_type &v, state_type &dv) ->  void {
    dv.at(kIndexX) = v.at(kIndexX) / (v.at(kIndexY2) - v.at(kIndexY1));
    dv.at(kIndexY1) = -eval_dH_dx(v.at(kIndexX), std::exp(v.at(kIndexY1))) / eval_dH_dy(v.at(kIndexX), std::exp(v.at(kIndexY1))) * dv.at(kIndexX) / std::exp(v.at(kIndexY1));
    dv.at(kIndexY2) = -eval_dH_dx(v.at(kIndexX), std::exp(v.at(kIndexY2))) / eval_dH_dy(v.at(kIndexX), std::exp(v.at(kIndexY2))) * dv.at(kIndexX) / std::exp(v.at(kIndexY2));
}

auto SW_curve::get_fiber(const cplx &x) -> std::vector<cplx> {
    GiNaC::ex fiber_poly = H_.subs({x_ == complex_to_ex(x)});
    std::vector<cplx> fiber = roots(fiber_poly, y_);
    return fiber;
}

auto SW_curve::match_fiber(std::vector<state_type>::iterator v_it) -> void {
    cplx x = v_it->at(kIndexX);
    std::vector<cplx> fiber = get_fiber(x);
    for (auto it = std::next(v_it->begin()); it != v_it->end(); ++it) {
        cplx nearest_fiber = fiber.at(0);
        for (auto &f : fiber) {
            if (std::abs(std::exp(*it) - f) < std::abs(std::exp(*it) - nearest_fiber))
            {
                nearest_fiber = f;
            }
        }
        *it = std::log(nearest_fiber) + 2 * std::numbers::pi * static_cast<int>(((*it  + J*std::numbers::pi) / (2.0* J*std::numbers::pi)).real());
    }
    return;
}

auto SW_curve::get_branched_sheet(const cplx &x) -> cplx {
    std::vector<cplx> fiber = get_fiber(x);
    double diff = std::numeric_limits<double>::max();
    cplx sheet;
    for (auto it1 = std::rbegin(fiber); it1 != std::rend(fiber); it1++)
    {
        for (auto it2 = it1 + 1; it2 != std::rend(fiber); it2++){
            double new_diff = std::abs(*it1 - *it2);
            if (new_diff < diff)
            {
                diff = new_diff;
                sheet = (*it1 + *it2) / 2.0;
            }
        }
    }
    return sheet;
}
