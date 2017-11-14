/*******************************************************************************
* Copyright 2010      Stefano Sinigardi                                        *
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


class Data
{
public:
  Data();
  void setInputType(char*);
  int getInputType();
  int fillParametersFromFile(char*, std::ofstream&);
  void fillData();
  int getNelectrons();
  void setNelectrons_UNSECURED(int);
  void setNelectrons();
  int getNdim();
  void setNdim_UNSECURED(int);
  void setNdim();
  double getThermalV();
  void setThermalV_UNSECURED(double);
  void setThermalV(double);

  void setXaxis_UNSECURED(double, double);
  void setXaxis_NOGRID_UNSECURED(double);
  void setXaxis(int);
  int getNgridPointsX();

  void setYaxis_UNSECURED(double, double);
  void setYaxis_NOGRID_UNSECURED(double);
  void setYaxis(int);
  int getNgridPointsY();

  void setZaxis_UNSECURED(double, double);
  void setZaxis_NOGRID_UNSECURED(double);
  void setZaxis(int);
  int getNgridPointsZ();

  int getNsteps();
  void setNsteps_UNSECURED(int);
  void setNsteps(int);

  double getDeltaT();
  void setDeltaT_UNSECURED(double);
  void setDeltaT(double);

  int getEvoType();
  void setEvoType_UNSECURED(int);
  void setEvoType();

  int getParticleFillingMethod();
  void setParticleFillingMethod_UNSECURED(int);
  void setParticleFillingMethod();

  double getDimX();
  double getDimY();
  double getDimZ();
  double getDeltaX();
  double getDeltaY();
  double getDeltaZ();

  void setK();
  void setK_UNSECURED(double);
  double getK();
  void setA();
  void setA_UNSECURED(double);
  double getA();

  //  void setAnalyticType();
  //  void setAnalyticType_UNSECURED(int );
  //  int getAnalyticType();

  void aumentaT(double);
  void resettaT();
  void impostaT(double);
  void impostaT_UNSECURED(double);
  double getT();

  ~Data();

private:
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

};

class Field
{
public:
  Field();
  void setEx(double);
  void setEy(double);
  void setEz(double);
  double getEx();
  double getEy();
  double getEz();
  void setBx(double);
  void setBy(double);
  void setBz(double);
  double getBx();
  double getBy();
  double getBz();
  ~Field();

private:
  double ex;
  double ey;
  double ez;
  double e4;  //e4 allocated to keep a double4 dimension
  double bx;
  double by;
  double bz;
  double b4;  //b4 allocated to keep a double4/double8 dimension, speedup even if not used
};

class Particle
{
public:
  Particle();
  void setParticleX(double);
  void setParticleX_fromUser(double, double);
  void setParticleY(double);
  void setParticleY_fromUser(double, double);
  void setParticleZ(double);
  void setParticleZ_fromUser(double, double);
  void setParticlePX(double);
  void setParticlePX_fromUser(double);
  void setParticlePY(double);
  void setParticlePY_fromUser(double);
  void setParticlePZ(double);
  void setParticlePZ_fromUser(double);
  void setParticleQ(int);
  void setParticleM(int);
  void setParticleCell(Data);

  double getParticleX();
  double getParticleY();
  double getParticleZ();
  double getParticlePX();
  double getParticlePY();
  double getParticlePZ();
  int getParticleQ();
  int getParticleM();
  int getParticleCell();
  ~Particle();
  Particle &operator*(const double);


private:
  double x;
  double y;
  double z;
  double t; //t allocated to keep a double4 dimension
  double px;
  double py;
  double pz;
  double pt;  //pt allocated to keep a double4/double8 dimension

  int q;
  int m;
  int cell;
};



#endif //__DATATYPES_HPP
