#ifndef ORBIT_CONTROL_H
#define ORBIT_CONTROL_H

#include "common.h"
#include "bbh_system.h"
#include "stellar_model.h"

class Orbit_control {
public:
  Orbit_control(const Bbh_system& bhs, const Stellar_model& stellar_model)
    : bhs(bhs), stellar_model(stellar_model) {}
  double energy(const Coords& V, const double t) const;
  static double angular_momentum(const Coords& V);
  static double kinetic_energy(const Coords& star_coords);
  int status(const Coords& coords, const double t) const;

/////////private:
  const Bbh_system& bhs;
  const Stellar_model& stellar_model;
};

#endif