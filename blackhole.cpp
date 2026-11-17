#include <cmath>
#include <stdexcept>
#include "data.hpp"
using namespace std;

class schwarzschild_blackhole{
public:
    schwarzschild_blackhole(double mass_solar) { set_mass(mass_solar); }

    void set_mass(double mass_solar) {
        mass_kg_ = mass_solar * fixed::M_SUN;
        Rs = schwarzschild_radius(mass_kg_);
    }

    //코시슈바르츠 반지름 2(G*kg)/(c^2)
    double schwarzschild_radius(double mass_kg){
    return 2*(fixed::G * mass_kg) / pow(fixed::c, 2);
    }

    // 호킹 온도 T = hbar*c^3/(8*pi*G*M*kB)
    double hoking_temp (){
    return (fixed::HBAR*(pow(fixed::c, 3)) / (8*fixed::PI*fixed::G*mass_kg_*fixed::KB));
    }

    //r 반지름 좌표
    double metric_f(double r) const {
        if (r <= Rs){
            throw std::domain_error("r must be greater than r_s");}
            
        return 1.0 - (Rs / r);
        }

    // 사건의 지평선 반지름 (m)
    double event_horizon() const { return Rs; }

    // 시간 지연율: 먼 곳 시계가 1초 갈 때, r 에 가만히 있는 시계가 가는 시간 = sqrt(f(r))   (0 ~ 1)
    double time_dilation_rate(double r) const {
        return std::sqrt(metric_f(r));
    }

    double g_tt(double r) const { return metric_f(r) * pow(fixed::c, 2); }
    double g_rr(double r) const { return -1.0 / metric_f(r); }
    double g_thth(double r) const { return -r * r; }
    double g_phph(double r, double theta) const {
        const double s = std::sin(theta);
        return -(pow(r, 2))*(pow(s,2));
    }

    // ds^2 = g_tt dt^2 + g_rr dr^2 + g_thth dtheta^2 + g_phph dphi^2
    double interval(double r, double theta, double dt, double dr,
                    double dtheta, double dphi) const {
        return g_tt(r) * pow(dt, 2) + g_rr(r) * pow(dr, 2) +
            g_thth(r) * pow(dtheta, 2) + g_phph(r, theta) * pow(dphi, 2);
    }

private:
    double mass_kg_ = 0.0;
    double Rs = 0.0;

};