#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>

class Integrator {
public:
  using Deriv_type = int (*) (double, const double[], double[], void*);
  Integrator(Deriv_type deriv, const size_t dim, const double eps_abs, double eps_rel, void* params)
  {
    s = gsl_odeiv2_step_alloc(step_type, dim);
    c = gsl_odeiv2_control_y_new(eps_abs, eps_rel);
    e = gsl_odeiv2_evolve_alloc(dim);
    sys = {deriv, nullptr, dim, params};
  }
  Integrator(Deriv_type deriv, const size_t dim, const double eps_abs=1.0E-6, double eps_rel=1.0E-6)
    : Integrator(deriv, dim, eps_abs, eps_rel, (void*)nullptr) {}
  ~Integrator()
  {
    gsl_odeiv2_evolve_free(e);
    gsl_odeiv2_control_free(c);
    gsl_odeiv2_step_free(s);
  }
  int step(double& t, double t_end, double& h, double y[])
  {
    return gsl_odeiv2_evolve_apply(e, c, s, &sys, &t, t_end, &h, y);
  }
  void reset()
  {
    gsl_odeiv2_step_reset(s);
    gsl_odeiv2_evolve_reset(e);
  }
private:
  const gsl_odeiv2_step_type* step_type = gsl_odeiv2_step_rkck;
  gsl_odeiv2_step* s;
  gsl_odeiv2_control* c;
  gsl_odeiv2_evolve* e;
  gsl_odeiv2_system sys;
};

#endif