/*******************************************************************************
* Copyright 2010-2018 Stefano Sinigardi                                        *
* The program is distributed under the terms of the GNU General Public License *
*******************************************************************************/

/**************************************************************************
    This file is part of PICcol.

    PICcol is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    PICcol is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with PICcol.  If not, see <http://www.gnu.org/licenses/>.
**************************************************************************/



#ifndef __EVOLUTION_H
#define __EVOLUTION_H

#include <vector>
#include <fstream>
#include "datatypes.h"


void evolveRK4_nogrid(Data, std::vector<Particle>&, std::vector<Field>&, std::ofstream&);
void evolveRK4_withgrid_onthefly(Data, std::vector<Particle>&, std::vector<Field>&, std::ofstream&);
void evolveRK4_withgrid(Data, std::vector<Particle>&, std::vector<Field>&, std::vector<Field>&, std::ofstream&);
void evolveLPF_nogrid(Data, std::vector<Particle>&, std::vector<Field>&, std::ofstream&);
void evolveLPF_withgrid_onthefly(Data, std::vector<Particle>&, std::vector<Field>&, std::ofstream&);
void evolveLPF_withgrid(Data, std::vector<Particle>&, std::vector<Field>&, std::vector<Field>&, std::ofstream&);
void evolveLPF4_nogrid_OLD(Data, std::vector<Particle>&, std::vector<Field>&, std::ofstream&);
void evolveLPF4_nogrid(Data, std::vector<Particle>&, std::ofstream&);
void evolveLPF4_withgrid_onthefly(Data, std::vector<Particle>&, std::vector<Field>&, std::ofstream&);
void evolveLPF4_withgrid(Data, std::vector<Particle>&, std::vector<Field>&, std::vector<Field>&, std::ofstream&);
void evolveLPF4_nogrid_TURCHETTI(Data, std::vector<Particle>&, std::vector<Field>&, std::ofstream&);

Field obtainFOP(Data, Particle, std::vector<Field>&);
Field interpolation2D_linear(Data, Particle, std::vector<Field>&);
Field interpolation2D_quadratic(Data, Particle, std::vector<Field>&);

void evolveRK4(Data, Particle, Field, std::vector<Particle>&, std::vector<Particle>&);
void evolveLPF(Data, Particle&, double, double);

void calcolaCampoInUnIstante(Data, std::vector<Particle>&, std::vector<Field>&, std::ofstream&);
void calcolaEvoluzioneCampoSuParticelle(Data, std::vector<Particle>&, std::vector<Field>&, std::ofstream&);

#endif //__EVOLUTION_HPP
