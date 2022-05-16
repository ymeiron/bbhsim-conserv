#ifndef COMMON_H
#define COMMON_H

#include <array>

#include "definitions.h"

struct Coords {
  union {
    struct {
      double x, y, z, vx, vy, vz;
    };
    double data[6];
  };
  Coords() {}
  Coords(const double x, const double y, const double z, const double vx,
         const double vy, const double vz)
      : x(x), y(y), z(z), vx(vx), vy(vy), vz(vz) {}
  double &operator[](int i) { return data[i]; }
  const double operator[](int i) const { return data[i]; }
};

struct Coords2d {
  double x, y, vx, vy;
};

struct Bbh_coords {
  union {
    struct {
      Coords2d bh1, bh2;
    };
    double data[8];
  };
  double &operator[](int i) { return data[i]; }
  const double operator[](int i) const { return data[i]; }
  std::array<double,4> extract_pos() const
  {
    return {bh1.x, bh1.y, bh2.x, bh2.y};
  }
  std::array<double,4> extract_vel() const
  {
    return {bh1.vx, bh1.vy, bh2.vx, bh2.vy};
  }
};

#endif