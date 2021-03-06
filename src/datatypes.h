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

#ifndef __DATATYPES_H
#define __DATATYPES_H

#include <fstream>
#include <string>

class Data {
  public:
  Data();
  int fillParametersFromFile(char*, std::ofstream&);
  void setNelectrons(int);
  void setNdim(int);
  void setThermalV(double);
  void setXaxis(double, double);
  void setYaxis(double, double);
  void setZaxis(double, double);
  void setNsteps(int);
  void setDeltaT(double);
  void setEvoType(int);
  void setParticleFillingMethod(int);

  void setK(double);
  void setA(double);

  void aumentaT(double);
  void resettaT();
  void setInitialT(double);

  int inputType;
  int n_electrons;
  int n_ions;
  int n_dim;

  double k, A;

  int evoType;
  int analyticType;
  int particleFillingMethod;

  int nGridPointsX;
  int nGridPointsY;
  int nGridPointsZ;
  double dimX;
  double dimY;
  double dimZ;
  double deltaX;
  double deltaY;
  double deltaZ;
  double thermal_vel;
  double dt, t;
  int nSteps;
  std::string filename;
};

class Field {
  public:
  Field();
  void setEx(double);
  void setEy(double);
  void setEz(double);
  void setBx(double);
  void setBy(double);
  void setBz(double);

  double ex;
  double ey;
  double ez;
  double e4; //e4 allocated to keep a double4 dimension
  double bx;
  double by;
  double bz;
  double b4; //b4 allocated to keep a double4/double8 dimension, speedup even if not used
};

class Particle {
  public:
  Particle();

  void setParticleX(double);
  void setParticleY(double);
  void setParticleZ(double);
  void setParticlePX(double);
  void setParticlePY(double);
  void setParticlePZ(double);
  void setParticleQ(int);
  void setParticleM(int);

  double x;
  double y;
  double z;
  double t; //t allocated to keep a double4 dimension
  double px;
  double py;
  double pz;
  double pt; //pt allocated to keep a double4/double8 dimension

  int q;
  int m;
  int cell;
};

#endif //__DATATYPES_H
