#include <iostream>
#include <libconfig.h++>

#include "config.h"

template<typename T>
void parse_setting(const libconfig::Setting& setting, T& output)
{
  output = (T)setting;
}

template<>
void parse_setting<double>(const libconfig::Setting& setting, double& output)
{
  double result;
  switch(setting.getType()) {
    case libconfig::Setting::Type::TypeFloat:
      output = (double)setting;
      break;
    case libconfig::Setting::Type::TypeInt:
      output = (double)(int)setting;
      break;
    case libconfig::Setting::Type::TypeInt64:
      output = (double)(long)setting;
      break;
    default:
      throw libconfig::SettingTypeException(setting);
  }
}

template<typename T>
void get_setting(libconfig::Config& cfg, const std::string& key, T& output, const T default_value)
{
  try {
    libconfig::Setting& setting = cfg.lookup(key);
  } catch (const libconfig::SettingNotFoundException &e) {
    output = default_value;
    return;
  }

  try {
    libconfig::Setting& setting = cfg.lookup(key);
    parse_setting(setting, output);
  } catch (const libconfig::SettingTypeException &e) {
    std::cerr << "The setting \"" << key << "\" has an incorrect type.\n";
    exit(EXIT_FAILURE);
  }
}

template<typename T>
void get_setting(libconfig::Config& cfg, const std::string& key, T& output)
{
  try {
    libconfig::Setting& setting = cfg.lookup(key);
  } catch (const libconfig::SettingNotFoundException &e) {
    std::cerr << "Mandatory setting \"" << key << "\" not found in the configuration file.\n";
    exit(EXIT_FAILURE);
  }

  try {
    libconfig::Setting& setting = cfg.lookup(key);
    parse_setting(setting, output);
  } catch (const libconfig::SettingTypeException &e) {
    std::cerr << "The setting \"" << key << "\" has an incorrect type.\n";
    exit(EXIT_FAILURE);
  }
}

void get_setting(libconfig::Config& cfg, const std::string& key, std::string& output, const char default_value[])
{
  get_setting(cfg, key, output, std::string(default_value));
}

Config::Config(const std::string& filename)
{
  libconfig::Config cfg;
  try {
    cfg.readFile(filename);
  } catch (const libconfig::FileIOException &e) {
    std::cerr << "Error reading configuration file.\n";
    exit(EXIT_FAILURE);
  } catch (libconfig::ParseException &e) {
    std::cerr << "Parse error at " << e.getFile() << ":" << e.getLine() << " - " << e.getError() << "\n";
    exit(EXIT_FAILURE);
  }

  get_setting(cfg, "q",         q        );
  get_setting(cfg, "R0",        R0       );
  get_setting(cfg, "vx",        vx       );
  get_setting(cfg, "vy",        vy       );
  get_setting(cfg, "dt_fric",   dt_fric  );
  get_setting(cfg, "t_end",     t_end    );
  get_setting(cfg, "N_stars",   N_stars  );
  get_setting(cfg, "a_plummer", a_plummer);
  get_setting(cfg, "m_plummer", m_plummer);

  get_setting(cfg, "bh_eps_abs",   bh_eps_abs,   1.0E-06);
  get_setting(cfg, "bh_eps_rel",   bh_eps_rel,   1.0E-06);
  get_setting(cfg, "star_eps_abs", star_eps_rel, 1.0E-06);
  get_setting(cfg, "star_eps_rel", star_eps_rel, 1.0E-06);

  get_setting(cfg, "description", description, "");
  get_setting(cfg, "output_mode", output_mode, "file");
}