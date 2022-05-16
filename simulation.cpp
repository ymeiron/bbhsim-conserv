#include <cmath>
#include <unistd.h>

#include "common.h"
#include "config.h"
#include "definitions.h"
#include "kepler.h"
#include "integrator.h"
#include "io.h"
#include "stellar_model.h"
#include "stellar_system.h"

/******************************************************************************/

int deriv_bhs(double t, const double V[], double dVdt[], void* params)
{
  auto& bhs = *(Bbh_system*)params;

  double r_rel;
  if (bhs.q < 1) {
    dVdt[0] = 0;
    dVdt[1] = 0;
    dVdt[2] = 0;
    dVdt[3] = 0;
    dVdt[4] = V[6];
    dVdt[5] = V[7];
    r_rel = sqrt(V[4] * V[4] + V[5] * V[5]);
  } else {
    dVdt[0] = V[2];
    dVdt[1] = V[3];
    dVdt[4] = V[6];
    dVdt[5] = V[7];
    r_rel = sqrt((V[0] - V[4]) * (V[0] - V[4]) + (V[1] - V[5]) * (V[1] - V[5]));
  }

  double r_rel_cube = r_rel * r_rel * r_rel;
  if (r_rel < TIDAL_RADIUS) {
    return GSL_EFAILED;
  }

  if (bhs.q < 1) {
    dVdt[6] = -V[4] / r_rel_cube;
    dVdt[7] = -V[5] / r_rel_cube;
  } else {
    dVdt[2] = -(V[0] - V[4]) / r_rel_cube;
    dVdt[3] = -(V[1] - V[5]) / r_rel_cube;
    // Newton's third law:
    dVdt[6] = -dVdt[2];
    dVdt[7] = -dVdt[3];
  }

  // The external force (friction like)
  double absv = sqrt(V[6] * V[6] + V[7] * V[7]);
  dVdt[6] += -bhs.f2 * V[6] / (absv * bhs.q);
  dVdt[7] += -bhs.f2 * V[7] / (absv * bhs.q);
  dVdt[6] += bhs.ff2 * V[7] / (absv * bhs.q);
  dVdt[7] += -bhs.ff2 * V[6] / (absv * bhs.q);

  if (bhs.q == 1) {
    dVdt[2] += -bhs.f1 * V[2] / absv;
    dVdt[3] += -bhs.f1 * V[3] / absv;
    dVdt[2] += bhs.ff1 * V[3] / absv;
    dVdt[3] += -bhs.ff1 * V[2] / absv;
  }

  return GSL_SUCCESS;
}

void evolve_bhs(Bbh_system& bhs, Integrator& integrator, double t0, double dt) {
  double t = t0;

  double h = 1;
  bhs.segment.clear();
  auto V = bhs.coords;
  while (t < t0+dt) {
    bhs.segment.emplace_back(t, V);
    int status = integrator.step(t, t0+dt, h, V.data);
    if (status != GSL_SUCCESS) {
      fprintf(stderr, "ODE solver failed... status=%d\n", status);
      exit(1);
    }
  }
  bhs.segment.emplace_back(t, V);
  bhs.coords = V;
}

int deriv_star(double t, const double V[], double dVdt[], void* params)
/*************************************
 * INPUT:
 * double t : time of the quary
 * double y[] : some generalized coordinate vector in phase space (vector 1x6)
 * OUTPUT:
 * double dydx[] : the RHS of the equations of motion for generalized
 *coordinates y[] (vector 1x6) default return value : 0
 *************************************/
{
  auto& [bhs, stellar_model] = *(std::tuple<Bbh_system&, Stellar_model&>*)params;

  auto BHCoordinates = bhs.interp_pos(t);
  double q = bhs.q;

  dVdt[0] = V[3];
  dVdt[1] = V[4];
  dVdt[2] = V[5];

  const double x1 = (V[0] - BHCoordinates[0]);
  const double y1 = (V[1] - BHCoordinates[1]);
  const double x2 = (V[0] - BHCoordinates[2]);
  const double y2 = (V[1] - BHCoordinates[3]);

  double r = sqrt(V[0]*V[0] + V[1]*V[1] + V[2]*V[2]);
  double r1 = sqrt(x1*x1 + y1*y1 + V[2]*V[2]);
  double r2 = sqrt(x2*x2 + y2*y2 + V[2]*V[2]);

  double r1cube = r1*r1*r1;
  double r2cube = r2*r2*r2;
  dVdt[3] = -M * x1 / r1cube;
  dVdt[4] = -M * y1 / r1cube;
  dVdt[5] = -M * V[2] / r1cube;
  dVdt[3] += -q * M * x2 / r2cube;
  dVdt[4] += -q * M * y2 / r2cube;
  dVdt[5] += -q * M * V[2] / r2cube;

  // ** Bulge Potential **
  double acc = stellar_model.acc(r);
  dVdt[3] += - acc * V[0]/r;
  dVdt[4] += - acc * V[1]/r;
  dVdt[5] += - acc * V[2]/r;

  return 0;
}

int analytical_solver(const Bbh_system& bhs, double V[6], int bh_num, double *t, double dt) {
  // this mostly does the transformation and calls AdvanceKepler
  // BH = 0 or 1
  double Vtmp[6];
  double hdid;
  int ReturnValue;
  int offset = bh_num * 2;
  double Mass;

  fprintf(stderr, "This is AnalyticalSolver: DeltaT=%f", dt);

  auto bhs_pos = bhs.interp_pos(*t); // obviously this check is done in derivs, we'll use *stat to transfer the info
  auto bhs_vel = bhs.interp_vel(*t);
  Vtmp[0] = V[0] - bhs_pos[offset];
  Vtmp[1] = V[1] - bhs_pos[offset+1];
  Vtmp[2] = V[2];
  Vtmp[3] = V[3] - bhs_vel[offset];
  Vtmp[4] = V[4] - bhs_vel[offset+1];
  Vtmp[5] = V[5];
  Mass = (bh_num == 0) ? 1 : bhs.q;
  ReturnValue = evolve_kepler(Mass, Vtmp, dt, CLOSENESS + 0.0001, &hdid);
  *t += hdid;

  bhs_pos = bhs.interp_pos(*t); // thats the new t
  bhs_vel  = bhs.interp_vel(*t);
  V[0] = Vtmp[0] + bhs_pos[offset];
  V[1] = Vtmp[1] + bhs_pos[offset+1];
  V[2] = Vtmp[2];
  V[3] = Vtmp[3] + bhs_vel[offset];
  V[4] = Vtmp[4] + bhs_vel[offset+1];
  V[5] = Vtmp[5];

  return ReturnValue;
}

void evolve_stars(Stellar_system& stars, Integrator& integrator, const Orbit_control& control, const double t0, const double dt) {
  int Star;
  double t = t0;
  double h;
#define AllowedSteps 350000000

// I think that the pragma should also include r1cube, r2cube as privates...
#pragma omp parallel for private(V, dVdt, t, h, Steps, r, r1, r2, E_current, L_current)
  for (Star = 0; Star < stars.coords.size(); Star++) {
    int Steps = 0;
    int &stat = stars.stat[Star];
    if (stat >= ORBIT_STAT_TERMINATE) continue;
    Coords& V = stars.coords[Star];

    t = t0;
    while ((t < t0 + dt) && (stat < ORBIT_STAT_TERMINATE)) {
      h = dt + t0 - t;
      stat = control.status(V, t);
      if (stat != 0) { // there are special cases: very close to 1, 2 or very
                       // far from the centre (then check energy). derivsStar
                       // won't return anything >= ORBIT_STAT_TERMINATE
        if (stat == ORBIT_STAT_VERY_FAR) {
          if (stars.E[Star] > 0) {
            stat = ORBIT_STAT_DIVERGENT;
            break;
          }
          /// NOTE This should be changed to something like StarsE[Star] >
          /// EscapeEnergy when there is a bulge. Note that when the potential
          /// diverges there is little point in doing so and the point is that
          /// the stars won't go back to the centre due to orbital change due to
          /// unsmooth potential.
          // otherwise we are very far but haven't diverged yet.
        } else if ((stat == ORBIT_STAT_CLOSETOBH1) ||
                   (stat == ORBIT_STAT_CLOSETOBH2)) {
          if (dt - t + t0 > 100)
            fprintf(stderr, "XXXXXXXXX %d\n", Star);
          if (stat == ORBIT_STAT_CLOSETOBH1)
            stat = analytical_solver(control.bhs, V.data, 0, &t, dt - t + t0);
          else if (stat == ORBIT_STAT_CLOSETOBH2)
            stat = analytical_solver(control.bhs, V.data, 1, &t, dt - t + t0);
          continue; // at any rate, don't do the rest of this iteration if you
                    // entered here.
        }
      }
      int integration_stat = integrator.step(t, t0+dt, h, V.data);
      integrator.reset();

      if (integration_stat != GSL_SUCCESS) {
        stat = ORBIT_INTEGRATOR_FAIL;
        break;
      }
      if (++Steps > AllowedSteps) {
        stat = ORBIT_STAT_EXCEED_MAXIT;
        break;
      }
    }
  }
}

std::array<double,6> segment_integrals(Bbh_system& bhs) {
  Bbh_coords Vector1, Vector2, Derivative1, Derivative2;
  double S1, S2, P1, P2, Q1, Q2;
  double a0, b0, v1, v2, dvdt1, dvdt2, p1, p2, dpdt1, dpdt2,
      q1, q2, dqdt1,
      dqdt2; // Beware! q1 & q2 have nothing to do with the mass ratio q.
  S1 = 0;
  S2 = 0;
  P1 = 0;
  P2 = 0;
  Q1 = 0;
  Q2 = 0;
  for (int row = 0; row < bhs.segment.size() - 1; row++) {
    double t1 = bhs.segment[row].t;
    double t2 = bhs.segment[row+1].t;
    double Delta_t = t2 - t1;
    Vector1 = bhs.segment[row].coords;
    Vector2 = bhs.segment[row+1].coords;
    deriv_bhs(t1, Vector1.data, Derivative1.data, &bhs);
    deriv_bhs(t2, Vector2.data, Derivative2.data, &bhs);
    int Offset = (bhs.q < 1) ? 4 : 0;
    for (;;) {
      v1 = sqrt(Vector1[2 + Offset] * Vector1[2 + Offset] +
                Vector1[3 + Offset] * Vector1[3 + Offset]);
      v2 = sqrt(Vector2[2 + Offset] * Vector2[2 + Offset] +
                Vector2[3 + Offset] * Vector2[3 + Offset]);
      dvdt1 = (Derivative1[0 + Offset] * Derivative1[2 + Offset] +
               Derivative1[1 + Offset] * Derivative1[3 + Offset]) /
              v1;
      dvdt2 = (Derivative2[0 + Offset] * Derivative2[2 + Offset] +
               Derivative2[1 + Offset] * Derivative2[3 + Offset]) /
              v2;
      a0 = (Delta_t * (dvdt1 + dvdt2) - 2.0 * (v2 - v1)) / pow(Delta_t, 3);
      b0 = (-Delta_t * (2.0 * dvdt1 + dvdt2) + 3.0 * (v2 - v1)) /
           pow(Delta_t, 2);
      if (!Offset)
        S1 += v1 * Delta_t + 0.25 * a0 * pow(Delta_t, 4) +
               (1.0 / 3.0) * b0 * pow(Delta_t, 3) +
               0.5 * dvdt1 * pow(Delta_t, 2);
      else
        S2 += v1 * Delta_t + 0.25 * a0 * pow(Delta_t, 4) +
               (1.0 / 3.0) * b0 * pow(Delta_t, 3) +
               0.5 * dvdt1 * pow(Delta_t, 2);

      p1 = (Vector1[0 + Offset] * Vector1[3 + Offset] -
            Vector1[1 + Offset] * Vector1[2 + Offset]) /
           v1;
      p2 = (Vector2[0 + Offset] * Vector2[3 + Offset] -
            Vector2[1 + Offset] * Vector2[2 + Offset]) /
           v2;
      dpdt1 = (Vector1[0 + Offset] * Derivative1[3 + Offset] -
               Vector1[1 + Offset] * Derivative1[2 + Offset]) /
                  v1 -
              (p1 / v1) * dvdt1;
      dpdt2 = (Vector2[0 + Offset] * Derivative2[3 + Offset] -
               Vector2[1 + Offset] * Derivative2[2 + Offset]) /
                  v2 -
              (p2 / v2) * dvdt2;
      a0 = (Delta_t * (dpdt1 + dpdt2) - 2.0 * (p2 - p1)) / pow(Delta_t, 3);
      b0 = (-Delta_t * (2.0 * dpdt1 + dpdt2) + 3.0 * (p2 - p1)) /
           pow(Delta_t, 2);
      if (!Offset)
        P1 += p1 * Delta_t + 0.25 * a0 * pow(Delta_t, 4) +
               (1.0 / 3.0) * b0 * pow(Delta_t, 3) +
               0.5 * dpdt1 * pow(Delta_t, 2);
      else
        P2 += p1 * Delta_t + 0.25 * a0 * pow(Delta_t, 4) +
               (1.0 / 3.0) * b0 * pow(Delta_t, 3) +
               0.5 * dpdt1 * pow(Delta_t, 2);

      q1 = (Vector1[0 + Offset] * Vector1[2 + Offset] +
            Vector1[1 + Offset] * Vector1[3 + Offset]) /
           v1;
      q2 = (Vector2[0 + Offset] * Vector2[2 + Offset] +
            Vector2[1 + Offset] * Vector2[3 + Offset]) /
           v2;
      dqdt1 = v1 +
              (Vector1[0 + Offset] * Derivative1[2 + Offset] +
               Vector1[1 + Offset] * Derivative1[3 + Offset]) /
                  v1 -
              (q1 / v1) * dvdt1;
      dqdt2 = v2 +
              (Vector2[0 + Offset] * Derivative2[2 + Offset] +
               Vector2[1 + Offset] * Derivative2[3 + Offset]) /
                  v2 -
              (q2 / v2) * dvdt2;
      a0 = (Delta_t * (dqdt1 + dqdt2) - 2.0 * (q2 - q1)) / pow(Delta_t, 3);
      b0 = (-Delta_t * (2.0 * dqdt1 + dqdt2) + 3.0 * (q2 - q1)) /
           pow(Delta_t, 2);
      if (!Offset) {
        Q1 += q1 * Delta_t + 0.25 * a0 * pow(Delta_t, 4) +
               (1.0 / 3.0) * b0 * pow(Delta_t, 3) +
               0.5 * dqdt1 * pow(Delta_t, 2);
        Offset = 4;
      } else {
        Q2 += q1 * Delta_t + 0.25 * a0 * pow(Delta_t, 4) +
               (1.0 / 3.0) * b0 * pow(Delta_t, 3) +
               0.5 * dqdt1 * pow(Delta_t, 2);
        break;
      }
    }
  }
  return {S1, S2, P1, P2, Q1, Q2};
}

#define MAX1 0.15
#define MAX2 0.15
//#define MAX1 5000
//#define MAX2 5000
bool calculate_friction(Bbh_system& bhs, double *F1, double *FF1, double *F2, double *FF2,
                       double DeltaE, double DeltaL) {
  double f1_tmp, ff1_tmp;
  auto [S1, S2, P1, P2, Q1, Q2] = segment_integrals(bhs);

  if (bhs.q < 1) {
    f1_tmp =
        DeltaE / S2; // I'm using f1_tmp but it's actually the force on BH2.
    ff1_tmp = (DeltaL * S2 - DeltaE * P2) / (S2 * Q2);
  } else {
    f1_tmp = DeltaE / (2 * S1);
    ff1_tmp = (DeltaL * S1 - DeltaE * P1) / (2 * S1 * Q1);
  }

  if ((fabs(ff1_tmp) > MAX1) || (fabs(f1_tmp) > MAX2)) // This segment is bad.
  {
    if (bhs.q < 1) {
      *F2 = DeltaL / P2;
      *FF2 = 0;
    } else {
      *F1 = bhs.f1;
      *F2 = bhs.f2;
      *FF1 = bhs.ff1;
      *FF2 = bhs.ff2;
    }
    //fprintf(stderr, "!!!!!!!!!!!!!!!!!BAD SEG!!!!!!!!!!!!!!!!!!\n");
    return false;
  } else {
    //fprintf(stderr, "!!!!!!!!!!!!!!!!!OK SEG!!!!!!!!!!!!!!!!!!\n");
  }
  if (bhs.q == 1)
    *F1 = f1_tmp;
  *F2 = f1_tmp;
  if (bhs.q == 1)
    *FF1 = ff1_tmp;
  *FF2 = ff1_tmp;
  return true;
}

std::string time_string()
{
  time_t rawtime;
  time(&rawtime);
  auto timeinfo = localtime(&rawtime);
  std::string output(asctime(timeinfo));
  output.resize(24);
  return output;
}

int main(int argc, const char *argv[]) {
  double E, E_prev, dE, Etot;
  double L, L_prev, dL, Ltot;
  double current_time = 0;
  int bad_segment_count = 0;
  int random_seed = 0;
  srand(random_seed);

  std::string simulation_name = "default";
  if (argc == 2) simulation_name = argv[1];

  std::string conf_filename = simulation_name + ".conf";

  Config config(conf_filename);
  if ((config.q <= 0) || (config.q > 1)) {
    fprintf(stderr, "Mass ratio has to be greater than zero and less than or equal to one.\n");
    exit(EXIT_FAILURE);
  }

  Bbh_system bhs(config.q, config.R0, config.vx, config.vy);

  constexpr double h0 = 1;
  Plummer_model stellar_model(config.a_plummer, config.m_plummer, 0);
  auto generator = [&stellar_model](){ return stellar_model.generate(); };
  Stellar_system stars(config.N_stars, generator, config.m_plummer);
  stars.minimize_rotation();

  Orbit_control control(bhs, stellar_model);

  std::string log_filename = simulation_name + ".log";
  FILE* log_file = fopen(log_filename.c_str(), "w");

  FILE* data_file;
  if (config.output_mode == "file") {
    std::string data_filename = simulation_name + ".csv";
    data_file = fopen(data_filename.c_str(), "w");
  } else if (config.output_mode == "screen") {
    data_file = stdout;
  } else {
    fprintf(stderr, "The \"output_mode\" setting can be \"file\" or \"screen\".\n");
    exit(EXIT_FAILURE);
  }


  fprintf(log_file, "==========================================\n");
  fprintf(log_file, "Starting BSMBH, %s.\n", time_string().c_str());

  char hostname[128];
  gethostname(hostname, 128);
  fprintf(log_file, "Machine: %s (Linux) with %ld CPU(s).\n", hostname, sysconf(_SC_NPROCESSORS_ONLN));
  fprintf(log_file, "t_end = %d\n", config.t_end);
  fprintf(log_file, "t = %f\n", current_time);

  fprintf(log_file, "Results are directly written to standard output.\n");




  /************************* Simulation Specific Stuff **************************/
  fprintf(log_file, "=== Simulation specific ===\n");

  fprintf(log_file, "R0 = %f\n", config.R0);
  fprintf(log_file, "vy = %f\n", config.vy);
  fprintf(log_file, "vx = %f\n", config.vx);
  fprintf(log_file, "N  = %d\n", config.N_stars);
  fprintf(log_file, "===========================\n");


  int number_of_stops = config.t_end / config.dt_fric;
  std::vector<double> friction_update_times(number_of_stops);
  for (int i = 0; i < number_of_stops; i++) friction_update_times[i] = i * config.dt_fric;

  Integrator integrator_bbh(deriv_bhs, 8, config.bh_eps_abs, config.bh_eps_rel, &bhs);
  std::tuple<Bbh_system&, Stellar_model&> bhs_and_stellar_model = {bhs, stellar_model};
  Integrator integrator_star(deriv_star, 6, config.star_eps_abs, config.star_eps_rel, &bhs_and_stellar_model);

  /******************************************************************************/
  
  stars.update_energy_and_angmom(control, 0);
  E = stars.properties.E_tot;
  L = stars.properties.L_tot;
  int diverged = stars.properties.diverged;
  int crashed = stars.properties.crashed;

  E_prev = E;
  L_prev = L;
  auto [Ebh, Lbh] = bhs.energy_and_angmom();
  Etot = E + Ebh;
  Ltot = L + Lbh;
  fprintf(log_file, "Initial values: (E=%f, L=%f)\n", E_prev + Ebh, L_prev + Lbh);
  fflush(log_file);

  fprintf(data_file, "%f,%f,%f,%f,%f,%f,%f,%f,%f,%.18e,%.18e,%.18e,%.18e,%.18e,%.18e,%.18e,%.18e,%d\n",
         friction_update_times[0], bhs.coords[0], bhs.coords[1], bhs.coords[2], bhs.coords[3], bhs.coords[4], bhs.coords[5],
         bhs.coords[6], bhs.coords[7], Ebh, E, Lbh, L, bhs.f1, bhs.ff1, bhs.f2, bhs.ff2, bad_segment_count);
  fflush(stdout);

  int progress = 0;
  for (int i = 0; i < number_of_stops - 1; i++) {
    double dt = friction_update_times[i + 1] - friction_update_times[i];

    evolve_bhs(bhs, integrator_bbh, friction_update_times[i], dt);
    evolve_stars(stars, integrator_star, control, friction_update_times[i], dt);

    // Here calculate the energy and angular momentum of the star
    stars.update_energy_and_angmom(control, friction_update_times[i+1]);
    E = stars.properties.E_tot;
    L = stars.properties.L_tot;
    diverged = stars.properties.diverged;
    crashed = stars.properties.crashed;

    std::tie(Ebh, Lbh) = bhs.energy_and_angmom();
    //dE = E - E_prev;
    //dL = L - L_prev;
    //dE = - Etot + Ebh + E;
    //dL = - Ltot + Lbh + L;
    dE = 0.999 * (E - E_prev) + 0.001 * (-Etot + Ebh + E);
    dL = 0.999 * (L - L_prev) + 0.001 * (-Ltot + Lbh + L);
    bool good_segment = calculate_friction(bhs, &bhs.f1, &bhs.ff1, &bhs.f2, &bhs.ff2, dE, dL);

    E_prev = E;
    L_prev = L;
    bad_segment_count += !good_segment;

    for (int j = 1; j < bhs.segment.size(); j++) {
      fprintf(data_file, "%f,%f,%f,%f,%f,%f,%f,%f,%f,%.18e,%.18e,%.18e,%.18e,%.18e,%.18e,%.18e,%.18e,%d\n",
             bhs.segment[j].t,
             bhs.segment[j].coords[0],
             bhs.segment[j].coords[1],
             bhs.segment[j].coords[2],
             bhs.segment[j].coords[3],
             bhs.segment[j].coords[4],
             bhs.segment[j].coords[5],
             bhs.segment[j].coords[6],
             bhs.segment[j].coords[7], Ebh, E, Lbh, L, bhs.f1, bhs.ff1, bhs.f2,
             bhs.ff2, bad_segment_count);
      fflush(data_file);
    }

    if ((int)(i * 10.0 / (number_of_stops - 1)) > progress) {
      progress = (int)(i * 10.0 / (number_of_stops - 1));
      fprintf(log_file, "%d%% at %s\n", progress * 10, time_string().c_str());
      fflush(log_file);
    }
  }
#ifdef HAS_HDF5
  char h5_filename[128];
  sprintf(h5_filename, "%s_%04d.h5\0", simulation_name.c_str(), number_of_stops-1);
  save_simulation(bhs, stars, h5_filename, simulation_name, config.description);
#endif

  fprintf(log_file, "Finished at %s\n", time_string().c_str());
  std::tie(Ebh, Lbh) = bhs.energy_and_angmom();
  fprintf(log_file, "Final values: (E=%f, L=%f)\n", E + Ebh, L + Lbh);
  fprintf(log_file, "Diverging particles: %d\n", diverged);
  fprintf(log_file, "Crashing particles: %d\n", crashed);
  fprintf(log_file, "Bad segments: %d\n", bad_segment_count);
  fprintf(log_file, "==========================================\n");
  fclose(log_file);
  fclose(data_file);

  return 0;
}
