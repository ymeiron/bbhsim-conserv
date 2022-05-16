#include <cmath>
#include <tuple>

#include "bbh_system.h"

Bbh_system::Bbh_system(const double q, const double R0, const double vx, const double vy)
{
  this->q = q;
  if (q < 1) coords = {0, 0, 0, 0, R0, 0, vx, vy};
  else coords = {R0, 0, vx, vy, -R0, 0, -vx, -vy};
  segment.clear();
  segment.emplace_back(0, coords);
}

std::tuple<double,double> Bbh_system::energy_and_angmom() { //improve!! redundant copies
  const auto& bh1 = segment.back().coords.bh1;
  const auto& bh2 = segment.back().coords.bh2;
  const auto [x1, y1, vx1, vy1] = std::array {bh1.x, bh1.y, bh1.vx, bh1.vy};
  const auto [x2, y2, vx2, vy2] = std::array {bh2.x, bh2.y, bh2.vx, bh2.vy};
  const double r = sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2));
  const double E = 0.5 * (vx1*vx1 + vy1*vy1) + 0.5 * q * (vx2*vx2 + vy2*vy2) - q / r;
  const double L = (x1*vy1 - y1*vx1) + q * (x2*vy2 - y2*vx2);
  return {E, L};
}

double Bbh_system::pot(const Coords& star_coords, const double t) const
{
  auto [x1, y1, x2, y2] = interp_pos(t);
  double x = star_coords.x;
  double y = star_coords.y;
  double z = star_coords.z;
  double r1 = sqrt((x-x1)*(x-x1) + (y-y1)*(y-y1) + z*z);
  double r2 = sqrt((x-x2)*(x-x2) + (y-y2)*(y-y2) + z*z);
  double E = - M/r1 - q*M/r2;
  return E;
}

std::array<double,4> Bbh_system::interp_pos(const double t) const
{
  std::array<double,4> pos1, pos2, vel1, vel2, output;

  double t_segment_begin = segment.front().t;
  double t_segment_end   = segment.back().t;
  if (t <= t_segment_begin) return segment.front().coords.extract_pos();
  if (t >= t_segment_end)   return segment.back().coords.extract_pos();

  int i;
  for (i = 0; i < segment.size(); i++) // that's kind of a search
    if (segment[i].t > t)
      break;

  pos1 = segment[i-1].coords.extract_pos();
  pos2 = segment[ i ].coords.extract_pos();
  vel1 = segment[i-1].coords.extract_vel();
  vel2 = segment[ i ].coords.extract_vel();

  double t1 = segment[i-1].t;
  double t2 = segment[i].t;
  double dt = t2 - t1;
  double dt2 = dt*dt;
  double dt3 = dt2*dt;

  for (i = 0; i < 4; i++) {
    double a0 = ( dt * (vel1[i]+vel2[i])      - 2. * (pos2[i]-pos1[i])) / dt3;
    double b0 = (-dt * (2. * vel1[i]+vel2[i]) + 3. * (pos2[i]-pos1[i])) / dt2;
    output[i] = pos1[i] + a0 * (t-t1) * (t-t1) * (t-t1) + b0 * (t-t1) * (t-t1) + vel1[i] * (t-t1);
  }
  return output;
}

std::array<double,4> Bbh_system::interp_vel(const double t) const
{
  std::array<double,4> vel1, vel2, output;

  double t_segment_begin = segment.front().t;
  double t_segment_end   = segment.back().t;
  if (t <= t_segment_begin) return segment.front().coords.extract_vel();
  if (t >= t_segment_end)   return segment.back().coords.extract_vel();

  int i;
  for (i = 0; i < segment.size(); i++) // that's kind of a search
    if (segment[i].t > t)
      break;

  vel1 = segment[i-1].coords.extract_vel();
  vel2 = segment[ i ].coords.extract_vel();

  double t1 = segment[i-1].t;
  double t2 = segment[i].t;
  double dt = t2 - t1;

  for (i = 0; i < 4; i++) {
    double slope = (vel2[i]-vel1[i]) / dt;
    output[i] = vel1[i] + slope * (t-t1);
  }
  return output;
}