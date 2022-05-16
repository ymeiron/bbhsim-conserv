#include "stellar_model.h"
#include <cmath>

double Plummer_model::pot(const double r) const 
{
  return - mass / sqrt(r*r + a*a);
}

double Plummer_model::acc(const double r) const 
{
  return mass * r * pow(r*r + a*a, -1.5);
}

Coords Plummer_model::generate() const
{
  constexpr double RMAX = 10;
  double R;
  double X[8];
  do {
    X[1] = distribution(random_engine);
    R = 1/sqrt(pow(X[1], -2./3.) - 1);
  } while (R > RMAX);
  do {
    X[4] = distribution(random_engine);
    X[5] = distribution(random_engine);
  } while (0.1*X[5] >= X[4]*X[4]*pow(1-X[4]*X[4], 3.5));
  for (auto i : {2, 3, 6, 7}) X[i] = distribution(random_engine);

  double z = (1 - 2*X[2])*R;
  double x = sqrt(R*R - z*z) * cos(2*M_PI*X[3]);
  double y = sqrt(R*R - z*z) * sin(2*M_PI*X[3]);

  double Ve = M_SQRT2*pow(1 + R*R, -0.25);
  double V = Ve*X[4];

  double vz = (1 - 2*X[6])*V;
  double vx = sqrt(V*V - vz*vz) * cos(2*M_PI*X[7]);
  double vy = sqrt(V*V - vz*vz) * sin(2*M_PI*X[7]);

  return {x*a, y*a, z*a, vx/sqrt(a/mass), vy/sqrt(a/mass), vz/sqrt(a/mass)};
}
