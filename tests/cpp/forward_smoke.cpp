#include "tpsf/forward.hpp"
#include <cmath>
#include <iostream>
#include <vector>
int main() {
    constexpr int n = 10;
    constexpr double diameter = 0.5;
    constexpr double area = 0.19634954084936207;
    tpsf::Configuration cfg{
        n, 0,
        pipe_sim::TransientParameters{100., 1., 0.65, diameter, area,
                                      288.15, 0.9, 518.3, 1.1e-5, 1.0e-5},
        pipe_sim::DemandRamp{5.0e6, 10., 10., 0.},
        pipe_sim::NewtonOptions{1.0e-8, 1.0e-10, 30, false}
    };
    tpsf::InitialState initial{
        std::vector<double>(n + 1, 5.0e6),
        std::vector<double>(n, 10.0), 10.0
    };
    auto r = tpsf::run_forward(cfg, initial);
    if (r.info != 0 || r.completed_steps != 0 ||
        r.linepack.size() != 1 ||
        !std::isfinite(r.linepack[0]) || r.linepack[0] <= 0.0) {
        std::cerr << "Forward smoke test failed\n";
        return 1;
    }
    std::cout << "PASS initial linepack=" << r.linepack[0] << " kg\n";
}
