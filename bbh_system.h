#ifndef BBH_SYSTEM_H
#define BBH_SYSTEM_H

#include <vector>
#include "common.h"

struct Segment {
  double t;
  Bbh_coords coords;
  Segment(double t, const Bbh_coords& coords) : t(t), coords(coords) {}
  Segment() {}
};

class Bbh_system {
public:
  Bbh_system(const double q, const double R0, const double vx, const double vy);
  Bbh_coords coords;
  double f1 = 0, f2 = 0, ff1 = 0, ff2 = 0;
  double q = -999;
  std::array<double,4> interp_pos(const double t) const;
  std::array<double,4> interp_vel(const double t) const;
  std::tuple<double,double> energy_and_angmom();
  double pot(const Coords& star_coords, const double t) const;
  std::vector<Segment> segment;
};

#endif