#pragma once

#include <cmath>
#include <iostream>
#include <string>

// burgers equation
// u0(x) = a + b sin(w x + phi) in [xl,xr]
// periodic boundary condition
class BurgersSinExact {
public:
    BurgersSinExact(double a, double b, double w, double phi, double ep)
        : m_a(a), m_b(b), m_w(w), m_phi(phi), m_ep(ep) {}

    double init_value(double x) const {
        return m_a + (m_b * sin((m_w * x) + m_phi));
    }

    double eval(double x, double t) const {
        double x2 = (m_w * x) + m_phi - (m_a * m_w * t);
        double t2 = m_b * m_w * t;
        return m_a + (m_b * eval_kernel(x2, t2, m_ep));
    }

    double eval_with_check(double x, double t) const {
        if (t >= get_tb())
            raise_error(x, t, "t >= tb: " + std::to_string(get_tb()));

        return eval(x, t);
    }

    double get_tb() const { return std::abs(1.0 / (m_b * m_w)); }

private:
    // u0(x) = sin(x)
    static double eval_kernel(double x, double t, double ep) {
        // check input
        if (t < 0) { raise_error(x, t, "t < 0: " + std::to_string(t)); }
        if (ep <= 0) { ep = 1e-6; }

        // keep x in [-pi,pi]
        while (x < -pi) { x += 2 * pi; }
        while (x >= pi) { x -= 2 * pi; }

        // initial value
        double k = 1.0 / (pi / 2 + t);
        double u = k * x;

        // Newton iteration
        // G(u) = u - sin(x - u * t) = 0
        const std::size_t iter_max = 1000000;
        std::size_t iter = 0;
        while (iter < iter_max) {
            double tmp = 1 + (cos(x - (u * t)) * t);
            double du = (u - sin(x - (u * t))) / tmp;
            u = u - du;

            iter++;
            if (std::abs(du) < ep) break;
        }

        return u;
    }

    static void raise_error(double x, double t, const std::string &msg) {
        std::string loc =
            " at (" + std::to_string(x) + "," + std::to_string(t) + ")";

        std::cerr << "BurgersExact: " + msg + loc;
        exit(1);
    }

    static constexpr double pi = 3.14159265358979323846;

    const double m_a;
    const double m_b;
    const double m_w;
    const double m_phi;
    const double m_ep;
};
