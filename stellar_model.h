#ifndef STELLAR_MODEL_H
#define STELLAR_MODEL_H

#include <cmath>
#include <random>
#include "common.h"

class Stellar_model {
public:
  virtual double pot(const double r) const = 0;
  virtual double acc(const double r) const = 0;
  virtual Coords generate() const = 0;
};


class Plummer_model : public Stellar_model {
public:
  Plummer_model(const double a, const double mass, const unsigned seed)
    : a(a), mass(mass), distribution(0,1), random_engine(seed) {}
  double pot(const double r) const override;
  double acc(const double r) const override;
  Coords generate() const override;

private:
  double a, mass;
  mutable std::default_random_engine random_engine;
  mutable std::uniform_real_distribution<double> distribution;
};

#endif