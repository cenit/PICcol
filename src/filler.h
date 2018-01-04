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

#ifndef __FILLER_H
#define __FILLER_H

#include "datatypes.h"
#include <vector>

void riempiPuntiGrigliaConCampoAnalitico1D(Data, std::vector<Field>&);
void riempiPuntiGrigliaConPotenzialeAnalitico1D(Data, std::vector<Field>&);
void calcolaCampoAnalitico1DSuParticelle(Data, std::vector<Particle>&, std::vector<Field>&);
void calcolaCampoAnalitico1DSuParticella(Data, Particle, Field&);
void calcolaPotenzialeAnalitico1DSuParticella(Data, Particle, Field&);
void creaVettoreParticelle(Data, std::vector<Particle>&);

#endif //__FILLER_HPP
