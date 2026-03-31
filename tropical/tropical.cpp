/*
 * tropical_gen.cpp
 * Generates tropical diagram segment data and writes to tropical_segments.csv
 *
 * Compile:  g++ -O2 -std=c++17 -o tropical_gen tropical_gen.cpp
 * Run:      ./tropical_gen
 */

#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <tuple>
#include <algorithm>
#include <cmath>
#include <string>
#include <limits>
#include <map>

#include "tropical.hpp"

// ─── 2-D vector helpers ──────────────────────────────────────────────────────
struct Vec2 { double x, y; };
Vec2 operator+(Vec2 a, Vec2 b)         { return {a.x+b.x, a.y+b.y}; }
Vec2 operator*(double t, Vec2 v)       { return {t*v.x,   t*v.y};   }
Vec2 operator*(Vec2 v, double t)       { return t*v; }

// ─── Linear form  f(x,y) = ax + by + c ───────────────────────────────────────
struct LinearForm { double a, b, c; };

double eval(const LinearForm& f, double x, double y) {
    return f.a*x + f.b*y + f.c;
}

// ─── Parametric line  P + t*V ────────────────────────────────────────────────
struct Line {
    int i, j;       // indices of the two linear forms that are equal here
    Vec2 p, v;      // base point and direction
};

// ─── Solve 2x2 linear system  A*[t1;t2] = b  (returns false if singular) ────
bool solve2x2(double a00, double a01, double a10, double a11,
              double b0,  double b1,
              double& t1, double& t2)
{
    double det = a00*a11 - a01*a10;
    if (std::abs(det) < 1e-12) return false;
    t1 = ( a11*b0 - a01*b1) / det;
    t2 = (-a10*b0 + a00*b1) / det;
    return true;
}

// ─── Build linear forms from tropical monomials ──────────────────────────────
std::vector<LinearForm> build_linear_forms(double Q_num)
{
    std::vector<std::array<int,3>> monomials = monomials_5_1;

    std::vector<LinearForm> forms;
    for (auto& m : monomials) {
        double a = m[0];
        double b = m[1];
        double c = m[2] * Q_num;
        forms.push_back({a, b, c});
    }
    return forms;
}

// ─── Dominated-monomial elimination ──────────────────────────────────────────
std::vector<LinearForm> eliminate_dominated(
        const std::vector<LinearForm>& forms,
        double bound,
        double eps = 1e-6)
{
    int n = (int)forms.size();

    std::vector<std::pair<double,double>> pts;

    for (double sx : {-1.0, 0.0, 1.0})
        for (double sy : {-1.0, 0.0, 1.0})
            pts.push_back({sx*bound, sy*bound});

    for (int i = 0; i < n; ++i) {
        for (int j = i+1; j < n; ++j) {
            double da = forms[i].a - forms[j].a;
            double db = forms[i].b - forms[j].b;
            double dc = forms[i].c - forms[j].c;

            for (int k = j+1; k < n; ++k) {
                double ea = forms[i].a - forms[k].a;
                double eb = forms[i].b - forms[k].b;
                double ec = forms[i].c - forms[k].c;
                double det = da*eb - db*ea;
                if (std::abs(det) < 1e-12) continue;
                double xi = (-dc*eb + db*ec) / det;
                double yi = ( da*ec - (-dc)*ea) / det;
                if (std::abs(xi) > bound || std::abs(yi) > bound) continue;

                pts.push_back({xi,      yi     });
                pts.push_back({xi+eps,  yi     });
                pts.push_back({xi-eps,  yi     });
                pts.push_back({xi,      yi+eps });
                pts.push_back({xi,      yi-eps });
            }
        }
    }

    std::vector<bool> active(n, false);
    for (auto [x, y] : pts) {
        double max_val = -std::numeric_limits<double>::infinity();
        for (int k = 0; k < n; ++k)
            max_val = std::max(max_val, eval(forms[k], x, y));

        for (int k = 0; k < n; ++k)
            if (std::abs(eval(forms[k], x, y) - max_val) < 1e-9)
                active[k] = true;
    }

    std::vector<LinearForm> survivors;
    int pruned = 0;
    for (int k = 0; k < n; ++k) {
        if (active[k]) {
            survivors.push_back(forms[k]);
        } else {
            ++pruned;
            std::cout << "  [pruned] f" << k << ": "
                      << forms[k].a << "*x + "
                      << forms[k].b << "*y + "
                      << forms[k].c << "\n";
        }
    }
    std::cout << "Eliminated " << pruned << " dominated form(s); "
              << survivors.size() << " remain.\n";
    return survivors;
}

// ─── Segment merging ─────────────────────────────────────────────────────────
// Strategy:
//   Each maximal segment lives on a unique carrier line.  Two segments belong
//   to the same carrier iff they have the same (normalised) direction AND their
//   base points differ only along that direction (i.e. they are collinear).
//
//   For every carrier we collect 1-D intervals [t_start, t_end] obtained by
//   projecting each segment's endpoints onto the direction vector, then merge
//   all overlapping / touching / gap-smaller-than-eps intervals into one.
//   The result is written back as 2-D segments via  P_base + t * dir.
// ─────────────────────────────────────────────────────────────────────────────

struct RawSegment { Vec2 start, end; };

// Canonical key for a carrier line:
//   • direction normalised so that the first nonzero component is positive
//   • offset = component of any point on the line perpendicular to the direction
// Both quantities are rounded to a grid of width KEY_EPS so that floating-point
// segments on the same geometric line map to the same bucket.
static const double KEY_EPS = 1e-7;

static long long quantise(double v) {
    return static_cast<long long>(std::round(v / KEY_EPS));
}

struct CarrierKey {
    long long dx, dy;   // normalised direction (both integers after quantisation)
    long long offset;   // perpendicular offset

    bool operator<(const CarrierKey& o) const {
        if (dx != o.dx) return dx < o.dx;
        if (dy != o.dy) return dy < o.dy;
        return offset < o.offset;
    }
};

// Normalise a direction vector to a canonical representative:
// scale so that gcd(|dx_q|,|dy_q|) == 1 and the first nonzero component > 0.
static std::pair<long long,long long> normalise_dir(double vx, double vy)
{
    // Quantise with a coarser grid for direction (slopes are rational)
    const double D_EPS = 1e-9;
    long long qx = static_cast<long long>(std::round(vx / D_EPS));
    long long qy = static_cast<long long>(std::round(vy / D_EPS));

    // Reduce by gcd
    long long g = std::__gcd(std::abs(qx), std::abs(qy));
    if (g == 0) g = 1;
    qx /= g;  qy /= g;

    // Make canonical: first nonzero component positive
    if (qx < 0 || (qx == 0 && qy < 0)) { qx = -qx; qy = -qy; }

    return {qx, qy};
}

std::vector<RawSegment> merge_segments(const std::vector<RawSegment>& raw,
                                       double merge_gap = 1e-6)
{
    // ── 1. Bucket segments by carrier ────────────────────────────────────────
    // For each bucket we store:
    //   • the canonical base point (first segment's start)
    //   • the unit direction vector
    //   • a list of (t_lo, t_hi) intervals projected onto that direction
    struct BucketData {
        Vec2  base;
        Vec2  dir;          // unit vector (for back-projection)
        std::vector<std::pair<double,double>> intervals;
    };
    std::map<CarrierKey, BucketData> buckets;

    for (auto& seg : raw) {
        double vx = seg.end.x - seg.start.x;
        double vy = seg.end.y - seg.start.y;
        double len = std::hypot(vx, vy);
        if (len < 1e-15) continue;   // degenerate point — skip
        double ux = vx / len, uy = vy / len;

        auto [qx, qy] = normalise_dir(ux, uy);

        // Perpendicular offset: for direction (ux,uy) the perp is (-uy, ux).
        // offset = -uy * x0 + ux * y0  (any point on the line gives the same value)
        double perp = -uy * seg.start.x + ux * seg.start.y;

        CarrierKey key{ qx, qy, quantise(perp) };

        auto it = buckets.find(key);
        if (it == buckets.end()) {
            BucketData bd;
            bd.base = seg.start;
            bd.dir  = {ux, uy};
            // Project both endpoints onto direction
            double t0 = ux * seg.start.x + uy * seg.start.y;
            double t1 = ux * seg.end.x   + uy * seg.end.y;
            if (t0 > t1) std::swap(t0, t1);
            bd.intervals.push_back({t0, t1});
            buckets[key] = std::move(bd);
        } else {
            auto& bd = it->second;
            // Re-derive ux,uy from the stored unit vector for consistency
            double t0 = bd.dir.x * seg.start.x + bd.dir.y * seg.start.y;
            double t1 = bd.dir.x * seg.end.x   + bd.dir.y * seg.end.y;
            if (t0 > t1) std::swap(t0, t1);
            bd.intervals.push_back({t0, t1});
        }
    }

    // ── 2. Merge intervals within each bucket ────────────────────────────────
    // Back-projection: given parameter t along unit dir, the 2D point is
    //   P = base + (t - t_base) * dir
    // where t_base = dot(base, dir) is the projection of the stored base point.
    // This correctly reconstructs a point on the carrier regardless of where
    // base sits on the line.
    std::vector<RawSegment> merged;

    for (auto& [key, bd] : buckets) {
        auto& ivs = bd.intervals;
        std::sort(ivs.begin(), ivs.end());

        // Precompute projection of base onto dir so back-projection is cheap
        double t_base = bd.dir.x * bd.base.x + bd.dir.y * bd.base.y;

        auto make_point = [&](double t) -> Vec2 {
            double dt = t - t_base;
            return { bd.base.x + dt * bd.dir.x,
                     bd.base.y + dt * bd.dir.y };
        };

        // Sweep: merge intervals whose gap is <= merge_gap
        double cur_lo = ivs[0].first, cur_hi = ivs[0].second;
        for (int k = 1; k < (int)ivs.size(); ++k) {
            if (ivs[k].first <= cur_hi + merge_gap) {
                // Overlapping or touching — extend
                cur_hi = std::max(cur_hi, ivs[k].second);
            } else {
                // Gap — emit current interval and start a new one
                merged.push_back({ make_point(cur_lo), make_point(cur_hi) });
                cur_lo = ivs[k].first;
                cur_hi = ivs[k].second;
            }
        }
        // Emit last interval
        merged.push_back({ make_point(cur_lo), make_point(cur_hi) });
    }

    return merged;
}

// ─── Main ─────────────────────────────────────────────────────────────────────
int main()
{
    const double Q_NUM   = 1.0;
    const double T_BOUND = 10.0;
    const int    N_PTS   = 50;     // kept for compatibility; unused here

    // 1. Build linear forms
    auto forms_all = build_linear_forms(Q_NUM);
    int n_all = (int)forms_all.size();
    std::cout << "Linear forms before pruning (" << n_all << "):\n";
    for (int k = 0; k < n_all; ++k)
        std::cout << "  f" << k << ": " << forms_all[k].a << "*x + "
                  << forms_all[k].b << "*y + " << forms_all[k].c << "\n";

    // 1b. Eliminate dominated monomials
    auto forms = eliminate_dominated(forms_all, T_BOUND);
    int n = (int)forms.size();
    std::cout << "Active linear forms (" << n << "):\n";
    for (int k = 0; k < n; ++k)
        std::cout << "  f" << k << ": " << forms[k].a << "*x + "
                  << forms[k].b << "*y + " << forms[k].c << "\n";

    // 2. Build parametric lines for each pair (i,j) where f_i = f_j
    std::vector<Line> lines;
    for (int i = 0; i < n; ++i) {
        for (int j = i+1; j < n; ++j) {
            double da = forms[i].a - forms[j].a;
            double db = forms[i].b - forms[j].b;
            double dc = forms[i].c - forms[j].c;

            if (std::abs(da) < 1e-12 && std::abs(db) < 1e-12) continue;

            Line L;
            L.i = i; L.j = j;
            if (std::abs(db) > 1e-12) {
                double y0    = -dc / db;
                double slope = -da  / db;
                L.p = {0.0, y0};
                L.v = {1.0, slope};
            } else {
                double x0 = -dc / da;
                L.p = {x0,  0.0};
                L.v = {0.0, 1.0};
            }
            lines.push_back(L);
        }
    }
    std::cout << "Total parametric lines: " << lines.size() << "\n";

    // 3. Compute all pairwise intersections
    int Ln = (int)lines.size();
    std::vector<std::vector<double>> t_hits(Ln);

    for (int a = 0; a < Ln; ++a) {
        for (int b = a+1; b < Ln; ++b) {
            double rhs0 = lines[b].p.x - lines[a].p.x;
            double rhs1 = lines[b].p.y - lines[a].p.y;
            double t1, t2;
            if (solve2x2( lines[a].v.x, -lines[b].v.x,
                          lines[a].v.y, -lines[b].v.y,
                          rhs0, rhs1, t1, t2))
            {
                t_hits[a].push_back(t1);
                t_hits[b].push_back(t2);
            }
        }
    }

    // 4. Check each sub-segment for maximality and collect raw segments
    std::vector<RawSegment> raw_segments;

    for (int idx = 0; idx < Ln; ++idx) {
        auto& li = lines[idx];
        auto tv = t_hits[idx];
        std::sort(tv.begin(), tv.end());

        std::vector<double> ts;
        ts.push_back(-T_BOUND);
        for (double t : tv)
            if (t > -T_BOUND && t < T_BOUND) ts.push_back(t);
        ts.push_back( T_BOUND);

        for (int k = 0; k+1 < (int)ts.size(); ++k) {
            double t_mid = 0.5*(ts[k] + ts[k+1]);
            Vec2 pt = li.p + t_mid*li.v;
            double xm = pt.x, ym = pt.y;

            double max_val = -std::numeric_limits<double>::infinity();
            for (auto& f : forms) max_val = std::max(max_val, eval(f, xm, ym));

            bool fi_max = std::abs(eval(forms[li.i], xm, ym) - max_val) < 1e-9;
            bool fj_max = std::abs(eval(forms[li.j], xm, ym) - max_val) < 1e-9;

            if (fi_max || fj_max) {
                Vec2 s = li.p + ts[k  ]*li.v;
                Vec2 e = li.p + ts[k+1]*li.v;
                raw_segments.push_back({s, e});
            }
        }
    }
    std::cout << "Raw maximal segments (before merging): " << raw_segments.size() << "\n";

    // 5. Merge collinear, overlapping / touching segments  <─── NEW STEP
    auto merged_segments = merge_segments(raw_segments, /*merge_gap=*/1e-6);
    std::cout << "Merged segments (effective lines): " << merged_segments.size() << "\n";

    // 6. Write CSV:  x0,y0,x1,y1  (one merged segment per row)
    std::ofstream out("/home/spamdoodler/Uni/Exponential_Networks_cpp/tropical/tropical_segments.csv");
    out << "x0,y0,x1,y1\n";
    for (auto& seg : merged_segments)
        out << seg.start.x << "," << seg.start.y << ","
            << seg.end.x   << "," << seg.end.y   << "\n";
    out.close();
    std::cout << "Wrote tropical_segments.csv\n";

    return 0;
}
