#ifndef IO_H
#define IO_H

#include <string>
#include "bbh_system.h"
#include "stellar_system.h"

void save_simulation(const Bbh_system& bhs, const Stellar_system& stars, std::string filename, std::string name, std::string description);

#endif