#include "stellar_system.h"

void Stellar_system::minimize_rotation()
{
  int N = coords.size();
  for (int i = 0; i < N/2; i++) {
    for (int j = 0; j < 6; j++) {
      int sign = (j<3)?-1:+1;
      coords[N/2 + i][j] = sign*coords[i][j];
    }
  }
}

void Stellar_system::update_energy_and_angmom(const Orbit_control& control, const double t)
{
  properties.E_tot    = 0;
  properties.L_tot    = 0;
  properties.diverged = 0;
  properties.crashed  = 0;
  for (int i = 0; i < coords.size(); i++) {
    E[i] = control.energy(coords[i], t);
    L[i] = control.angular_momentum(coords[i]);
    properties.E_tot += E[i];
    properties.L_tot += L[i];
    if (stat[i] == ORBIT_STAT_DIVERGENT) properties.diverged++;
    else if ((stat[i] == ORBIT_STAT_FALLTOBH2) || (stat[i] == ORBIT_STAT_FALLTOBH2)) properties.crashed++;
  }
  properties.E_tot *= m;
  properties.L_tot *= m;
}