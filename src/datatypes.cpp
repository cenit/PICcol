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

#include "datatypes.h"
#include "jsoncons/json.hpp"
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <ios>
#include <iostream>
#include <limits>

using namespace std;

Data::Data()
{
  inputType = 0;
  n_electrons = 0;
  n_ions = 0;
  n_dim = 0;
  nSteps = 0;
  k = 0;
  A = 0;

  nGridPointsX = 0;
  nGridPointsY = 0;
  nGridPointsZ = 0;

  dimX = 0.;
  dimY = 0.;
  dimZ = 0.;
  deltaX = 0.;
  deltaY = 0.;
  deltaZ = 0.;

  thermal_vel = -1.;
  dt = 0.;
  t = 0.;
  evoType = 0;
  analyticType = 0;
  particleFillingMethod = -1;
  filename = "";
}

int Data::fillParametersFromFile(char* nomefile, ofstream& outputLOG)
{
  ifstream inputstreamParameters;
  jsoncons::json parameters;
  inputstreamParameters.open(nomefile);
  if (inputstreamParameters.is_open()) {
    outputLOG << "Lettura parametri da " << nomefile << " in corso..." << endl;
    inputstreamParameters.close();
    try {
      parameters = jsoncons::json::parse_file(nomefile);
    } catch (std::exception& e) {
      outputLOG << e.what() << std::endl;
      exit(-1);
    }
  }

  jsoncons::json empty_json;
  jsoncons::json json_lattice_elements = parameters.has_member("Magnetic_elements") ? parameters["Magnetic_elements"] : empty_json;

  // filename = parameters.has_member("phasespace_file") ?
  // parameters["phasespace_file"].as<string>() : "input.dat";

  setNelectrons(parameters.has_member("Nelectrons") ? parameters["Nelectrons"].as<int>() : 0);
  outputLOG << "Numero particelle da simulare: " << n_electrons << endl;

  setNdim(parameters.has_member("Ndim") ? parameters["Ndim"].as<int>() : 0);
  outputLOG << "Sistema " << n_dim << "D" << endl;

  setThermalV(parameters.has_member("ThermalV") ? parameters["ThermalV"].as<double>() : 0.0);
  outputLOG << "Velocita' termica: " << thermal_vel << endl;

  double DimX = parameters.has_member("DimX") ? parameters["DimX"].as<double>() : 0.0;
  double dx = parameters.has_member("dx") ? parameters["dx"].as<double>() : 0.0;
  setXaxis(DimX, dx);
  outputLOG << "Parametri asse x: dimX: " << dimX << ", dx: " << deltaX << endl;

  double DimY = parameters.has_member("DimY") ? parameters["DimY"].as<double>() : 0.0;
  double dy = parameters.has_member("dy") ? parameters["dy"].as<double>() : 0.0;
  setXaxis(DimY, dy);
  outputLOG << "Parametri asse y: dimY: " << dimY << ", dy: " << deltaY << endl;

  double DimZ = parameters.has_member("DimZ") ? parameters["DimZ"].as<double>() : 0.0;
  double dz = parameters.has_member("dz") ? parameters["dz"].as<double>() : 0.0;
  setXaxis(DimZ, dz);
  outputLOG << "Parametri asse z: dimZ: " << dimZ << ", dz: " << deltaZ << endl;

  setInitialT(parameters.has_member("InitialT") ? parameters["InitialT"].as<double>() : 0.0);
  outputLOG << "t iniziale: " << t << endl;

  setDeltaT(parameters.has_member("dt") ? parameters["dt"].as<double>() : 0.0);
  outputLOG << "dt: " << dt << endl;

  setNsteps(parameters.has_member("nSteps") ? parameters["nSteps"].as<int>() : 0);
  outputLOG << "nSteps: " << nSteps << endl;

  setParticleFillingMethod(parameters.has_member("ParticleFillingMethod")
          ? parameters["ParticleFillingMethod"].as<int>()
          : 0);
  if (particleFillingMethod == 1)
    outputLOG << "Generazione random particelle" << endl;
  else if (particleFillingMethod == 2)
    outputLOG << "Particelle generate in 0,0,0" << endl;
  else
    outputLOG << "Metodo non valido" << endl;

  setEvoType(parameters.has_member("EvoType") ? parameters["EvoType"].as<int>() : 0);
  switch (evoType) {
  case 1:
    outputLOG << "Distribuzione campi su particelle in un istante definito" << endl;
    break;
  case 2:
    outputLOG << "Evoluzione campo su una particella nel tempo" << endl;
    break;
  case 3:
    outputLOG << "Evoluzione RungeKutta4 no griglia (campo analitico calcolato sulla particella)" << endl;
    break;
  case 4:
    outputLOG << "Evoluzione RungeKutta4 griglia on-the-fly" << endl;
    break;
  case 5:
    outputLOG << "Evoluzione RungeKutta4 griglia persistente" << endl;
    break;
  case 6:
    outputLOG << "Evoluzione Leapfrog2 no griglia" << endl;
    break;
  case 7:
    outputLOG << "Evoluzione Leapfrog2 griglia on-the-fly" << endl;
    break;
  case 8:
    outputLOG << "Evoluzione Leapfrog2 griglia persistente" << endl;
    break;
  case 9:
    outputLOG << "Evoluzione Leapfrog4 no griglia" << endl;
    break;
  case 10:
    outputLOG << "Evoluzione Leapfrog4 griglia on-the-fly" << endl;
    break;
  case 11:
    outputLOG << "Evoluzione Leapfrog4 griglia persistente" << endl;
    break;
  case 12:
    outputLOG << "Evoluzione Leapfrog4 no griglia - algoritmo TURCHETTI" << endl;
  default:
    outputLOG << "Evoluzione non riconosciuta" << endl;
    break;
  }

  setA(parameters.has_member("A") ? parameters["A"].as<double>() : 0.0);
  setK(parameters.has_member("K") ? parameters["K"].as<double>() : 0.0);
  outputLOG << "E = -" << A << " " << k << " sin[" << k << "(x - ct)]" << endl;
  outputLOG << "A = " << A << " cos[" << k << "(x - ct)]" << endl;
  outputLOG << "Lettura parametri completata!" << endl;
  return 0;
}

void Data::setNelectrons(int n_e)
{
  n_electrons = n_e;
  n_ions = n_e;
}

void Data::setNdim(int n_dims) { n_dim = n_dims; }

void Data::setThermalV(double t_vel) { thermal_vel = t_vel; }

void Data::setXaxis(double dimensioneX, double dx)
{
  dimX = dimensioneX;
  deltaX = dx;
  nGridPointsX = (int)((dimX / deltaX) + 1);
  if (n_dim == 1) {
    dimZ = dimY = 1;
    deltaZ = deltaY = 1.0;
    nGridPointsZ = nGridPointsY = 1;
  }
}

void Data::setYaxis(double dimensioneY, double dy)
{
  dimY = dimensioneY;
  deltaY = dy;
  nGridPointsY = (int)((dimY / deltaY) + 2);
  if (n_dim == 2) {
    dimZ = 1;
    deltaZ = 1.0;
    nGridPointsZ = 1;
  }
}

void Data::setZaxis(double dimensioneZ, double dz)
{
  dimZ = dimensioneZ;
  deltaZ = dz;
  nGridPointsZ = (int)((dimZ / deltaZ) + 2);
}

void Data::setEvoType(int et) { evoType = et; }

void Data::setNsteps(int ns) { nSteps = ns; }

void Data::setDeltaT(double deltaT) { dt = deltaT; }

void Data::setParticleFillingMethod(int pfm) { particleFillingMethod = pfm; }

void Data::setK(double K) { k = K; }

void Data::setA(double a) { A = a; }

void Data::aumentaT(double deltaT) { t += deltaT; }

void Data::resettaT() { t = 0.; }

void Data::setInitialT(double tempo) { t = tempo; }

Field::Field() { ex = ey = ez = e4 = bx = by = bz = b4 = 0.; }

Particle::Particle()
{
  x = y = z = t = px = py = pz = pt = 0.;

  m = q = cell = 0;
}

void Field::setEx(double Ex) { ex = Ex; }

void Field::setEy(double Ey) { ey = Ey; }

void Field::setEz(double Ez) { ez = Ez; }

void Field::setBx(double Bx) { bx = Bx; }

void Field::setBy(double By) { by = By; }

void Field::setBz(double Bz) { bz = Bz; }

void Particle::setParticleX(double X) { x = X; }

void Particle::setParticleY(double Y) { y = Y; }

void Particle::setParticleZ(double Z) { z = Z; }

void Particle::setParticlePX(double PX) { px = PX; }

void Particle::setParticlePY(double PY) { py = PY; }

void Particle::setParticlePZ(double PZ) { pz = PZ; }

void Particle::setParticleQ(int Q) { q = Q; }

void Particle::setParticleM(int mass) { m = mass; }
