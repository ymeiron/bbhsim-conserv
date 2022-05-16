#include <cmath>

#include "definitions.h"
#include "orbit_control.h"

double RMAX = 500;
double ORBIT_INFINITY = 510;

double Orbit_control::energy(const Coords& V, const double t) const
{
  double E = kinetic_energy(V); // make kinetic_energy a member
  E += bhs.pot(V, t);
  const double& x = V.x;
  const double& y = V.y;
  const double& z = V.z;
  double r = sqrt(x*x + y*y + z*z);
  E += stellar_model.pot(r);
  return E;
}

double Orbit_control::angular_momentum(const Coords& V)
{
  const double& x = V.x;
  const double& y = V.y;
  const double& vx = V.vx;
  const double& vy = V.vy;
  return x*vy - y*vx; 
}

double Orbit_control::kinetic_energy(const Coords& star_coords) {
  const double& vx = star_coords.vx;
  const double& vy = star_coords.vy;
  const double& vz = star_coords.vz;
  return 0.5*(vx*vx + vy*vy + vz*vz);
}

int Orbit_control::status(const Coords& coords, const double t) const
{
  auto BHCoordinates = bhs.interp_pos(t);

  double r_sqr  = coords.x*coords.x + coords.y*coords.y + coords.z*coords.z;
  if (r_sqr > ORBIT_INFINITY*ORBIT_INFINITY) return ORBIT_STAT_VERY_FAR;

  double x1 = (coords[0] - BHCoordinates[0]);
  double y1 = (coords[1] - BHCoordinates[1]);
  double x2 = (coords[0] - BHCoordinates[2]);
  double y2 = (coords[1] - BHCoordinates[3]);

  double r1_sqr = x1*x1 + y1*y1 + coords.z*coords.z;
  if (r1_sqr < CLOSENESS*CLOSENESS) return ORBIT_STAT_CLOSETOBH1;
  double r2_sqr = x2*x2 + y2*y2 + coords.z*coords.z;
  if (r2_sqr < CLOSENESS*CLOSENESS) return ORBIT_STAT_CLOSETOBH2;
  return ORBIT_STAT_STABLE;
}