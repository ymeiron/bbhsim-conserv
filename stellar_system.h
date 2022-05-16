#ifndef STELLAR_SYSTEM_H
#define STELLAR_SYSTEM_H

#include <algorithm>
#include <vector>
#include "common.h"
#include "orbit_control.h"

class Stellar_system {
public:
  template<class Generator>
  Stellar_system(size_t N, const Generator& generator, double total_mass = 1) :
    m(total_mass/N), coords(N), stat(N,0), E(N), L(N) {
      std::generate(begin(coords), end(coords), generator);
    }
  void minimize_rotation();
  void update_energy_and_angmom(const Orbit_control& control, const double t);

  double m;
  std::vector<Coords> coords;
  std::vector<int> stat;
  std::vector<double> E;
  std::vector<double> L;
  struct {
    double E_tot, L_tot;
    int diverged, crashed;
  } properties;
};

#endif