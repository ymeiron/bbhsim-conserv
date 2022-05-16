#ifndef CONFIG_H
#define CONFIG_H

#include <string>

struct Config {
  Config(const std::string& filename);
  double q;
  double R0;
  double vx;
  double vy;
  double dt_fric;
  double t_end;
  int N_stars;
  double a_plummer;
  double m_plummer;
  std::string description;

  double bh_eps_abs, bh_eps_rel, star_eps_abs, star_eps_rel;
  std::string output_mode;
};

#endif