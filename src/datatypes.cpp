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



#define DEBUG

#include <iostream>
#include <cstdio>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <fstream>
#include <limits> // for declaration of 'numeric_limits'
#include <ios>    // for declaration of 'streamsize'

#include "datatypes.h"

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
}


Data::~Data() { }


int Data::fillParametersFromFile(char* nomefile, ofstream &outputLOG)
{
  ifstream inputstreamParameters;
  inputstreamParameters.open(nomefile);
  if (inputstreamParameters.is_open())
  {
    outputLOG << "Lettura parametri da " << nomefile << " in corso..." << endl;
    double genericdouble, anothergenericdouble;
    int genericInt;
    inputstreamParameters >> genericInt;
    setNelectrons_UNSECURED(genericInt);
    outputLOG << "Numero particelle da simulare: " << getNelectrons() << endl;
    inputstreamParameters >> genericInt;
    setNdim_UNSECURED(genericInt);
    outputLOG << "Sistema " << getNdim() << "D" << endl;

    inputstreamParameters >> genericdouble;
    setThermalV_UNSECURED(genericdouble);
    outputLOG << "Velocita' termica: " << getThermalV() << endl;

    inputstreamParameters >> genericdouble;
    inputstreamParameters >> anothergenericdouble;
    setXaxis_UNSECURED(genericdouble, anothergenericdouble);

    outputLOG << "Parametri asse x: dimX: " << getDimX() << ", dx: " << getDeltaX() << endl;

    if (getNdim() == 2 || getNdim() == 3)
    {
      inputstreamParameters >> genericdouble;
      inputstreamParameters >> anothergenericdouble;
      setYaxis_UNSECURED(genericdouble, anothergenericdouble);
      outputLOG << "Parametri asse y: dimY: " << getDimY() << ", dy: " << getDeltaY() << endl;
    }

    if (getNdim() == 3)
    {
      inputstreamParameters >> genericdouble;
      inputstreamParameters >> anothergenericdouble;
      setZaxis_UNSECURED(genericdouble, anothergenericdouble);
      outputLOG << "Parametri asse z: dimZ: " << getDimZ() << ", dz: " << getDeltaZ() << endl;
    }

    inputstreamParameters >> genericdouble;
    impostaT_UNSECURED(genericdouble);
    outputLOG << "t iniziale: " << getT() << endl;

    inputstreamParameters >> genericdouble;
    setDeltaT_UNSECURED(genericdouble);
    outputLOG << "dt: " << getDeltaT() << endl;

    inputstreamParameters >> genericInt;
    setNsteps_UNSECURED(genericInt);
    outputLOG << "nSteps: " << getNsteps() << endl;

    inputstreamParameters >> genericInt;
    setParticleFillingMethod_UNSECURED(genericInt);
    if (getParticleFillingMethod() == 1)
      outputLOG << "Generazione random particelle" << endl;
    else if (getParticleFillingMethod() == 2)
      outputLOG << "Particelle generate in 0,0,0" << endl;

    inputstreamParameters >> genericInt;
    setEvoType_UNSECURED(genericInt);
    if (getEvoType() == 1)
    {
      outputLOG << "Distribuzione campi su particelle in un istante definito" << endl;
      inputstreamParameters >> genericdouble;
      setA_UNSECURED(genericdouble);
      inputstreamParameters >> genericdouble;
      setK_UNSECURED(genericdouble);
      outputLOG << "E = -" << getA() << " " << getK() << " sin[" << getK() << "(x - ct)]" << endl;
      outputLOG << "A = " << getA() << " cos[" << getK() << "(x - ct)]" << endl;
      //      inputstreamParameters >> genericInt;
      //      setAnalyticType_UNSECURED(genericInt);
      //      if (getAnalyticType() == 1) outputLOG << "Distribuzione campi su particelle in un istante definito" << endl;
      //      else if (getAnalyticType() == 2) outputLOG << "Evoluzione campo su una particella nel tempo" << endl;
      //      else if (getAnalyticType() == 3) outputLOG << "Evoluzione particelle metodo RK4 con calcolo campo analitico sulla particella" << endl;
      //      else outputLOG << "Metodo non riconosciuto" << endl;
    }
    else if (getEvoType() == 2)
    {
      outputLOG << "Evoluzione campo su una particella nel tempo" << endl;
      inputstreamParameters >> genericdouble;
      setA_UNSECURED(genericdouble);
      inputstreamParameters >> genericdouble;
      setK_UNSECURED(genericdouble);
      outputLOG << "E = -" << getA() << " " << getK() << " sin[" << getK() << "(x - ct)]" << endl;
      outputLOG << "A = " << getA() << " cos[" << getK() << "(x - ct)]" << endl;
    }
    else if (getEvoType() == 3)
    {
      outputLOG << "Evoluzione RungeKutta4 no griglia (campo analitico calcolato sulla particella)" << endl;
      inputstreamParameters >> genericdouble;
      setA_UNSECURED(genericdouble);
      inputstreamParameters >> genericdouble;
      setK_UNSECURED(genericdouble);
      outputLOG << "E = -" << getA() << " " << getK() << " sin[" << getK() << "(x - ct)]" << endl;
      outputLOG << "A = " << getA() << " cos[" << getK() << "(x - ct)]" << endl;
    }
    else if (getEvoType() == 4)
    {
      outputLOG << "Evoluzione RungeKutta4 griglia on-the-fly" << endl;
      inputstreamParameters >> genericdouble;
      setA_UNSECURED(genericdouble);
      inputstreamParameters >> genericdouble;
      setK_UNSECURED(genericdouble);
      outputLOG << "E = -" << getA() << " " << getK() << " sin[" << getK() << "(x - ct)]" << endl;
      outputLOG << "A = " << getA() << " cos[" << getK() << "(x - ct)]" << endl;
    }
    else if (getEvoType() == 5)
    {
      outputLOG << "Evoluzione RungeKutta4 griglia persistente" << endl;
      inputstreamParameters >> genericdouble;
      setA_UNSECURED(genericdouble);
      inputstreamParameters >> genericdouble;
      setK_UNSECURED(genericdouble);
      outputLOG << "E = -" << getA() << " " << getK() << " sin[" << getK() << "(x - ct)]" << endl;
      outputLOG << "A = " << getA() << " cos[" << getK() << "(x - ct)]" << endl;
    }
    else if (getEvoType() == 6)
    {
      outputLOG << "Evoluzione Leapfrog2 no griglia" << endl;
      inputstreamParameters >> genericdouble;
      setA_UNSECURED(genericdouble);
      inputstreamParameters >> genericdouble;
      setK_UNSECURED(genericdouble);
      outputLOG << "E = -" << getA() << " " << getK() << " sin[" << getK() << "(x - ct)]" << endl;
      outputLOG << "A = " << getA() << " cos[" << getK() << "(x - ct)]" << endl;
    }
    else if (getEvoType() == 7)
    {
      outputLOG << "Evoluzione Leapfrog2 griglia on-the-fly" << endl;
      inputstreamParameters >> genericdouble;
      setA_UNSECURED(genericdouble);
      inputstreamParameters >> genericdouble;
      setK_UNSECURED(genericdouble);
      outputLOG << "E = -" << getA() << " " << getK() << " sin[" << getK() << "(x - ct)]" << endl;
      outputLOG << "A = " << getA() << " cos[" << getK() << "(x - ct)]" << endl;
    }
    else if (getEvoType() == 8)
    {
      outputLOG << "Evoluzione Leapfrog2 griglia persistente" << endl;
      inputstreamParameters >> genericdouble;
      setA_UNSECURED(genericdouble);
      inputstreamParameters >> genericdouble;
      setK_UNSECURED(genericdouble);
      outputLOG << "E = -" << getA() << " " << getK() << " sin[" << getK() << "(x - ct)]" << endl;
      outputLOG << "A = " << getA() << " cos[" << getK() << "(x - ct)]" << endl;
    }
    else if (getEvoType() == 9)
    {
      outputLOG << "Evoluzione Leapfrog4 no griglia" << endl;
      inputstreamParameters >> genericdouble;
      setA_UNSECURED(genericdouble);
      inputstreamParameters >> genericdouble;
      setK_UNSECURED(genericdouble);
      outputLOG << "E = -" << getA() << " " << getK() << " sin[" << getK() << "(x - ct)]" << endl;
      outputLOG << "A = " << getA() << " cos[" << getK() << "(x - ct)]" << endl;
    }
    else if (getEvoType() == 10)
    {
      outputLOG << "Evoluzione Leapfrog4 griglia on-the-fly" << endl;
      inputstreamParameters >> genericdouble;
      setA_UNSECURED(genericdouble);
      inputstreamParameters >> genericdouble;
      setK_UNSECURED(genericdouble);
      outputLOG << "E = -" << getA() << " " << getK() << " sin[" << getK() << "(x - ct)]" << endl;
      outputLOG << "A = " << getA() << " cos[" << getK() << "(x - ct)]" << endl;
    }
    else if (getEvoType() == 11)
    {
      outputLOG << "Evoluzione Leapfrog4 griglia persistente" << endl;
      inputstreamParameters >> genericdouble;
      setA_UNSECURED(genericdouble);
      inputstreamParameters >> genericdouble;
      setK_UNSECURED(genericdouble);
      outputLOG << "E = -" << getA() << " " << getK() << " sin[" << getK() << "(x - ct)]" << endl;
      outputLOG << "A = " << getA() << " cos[" << getK() << "(x - ct)]" << endl;
    }
    else if (getEvoType() == 12)
    {
      outputLOG << "Evoluzione Leapfrog4 no griglia - algoritmo TURCHETTI" << endl;
      inputstreamParameters >> genericdouble;
      setA_UNSECURED(genericdouble);
      inputstreamParameters >> genericdouble;
      setK_UNSECURED(genericdouble);
      outputLOG << "E = -" << getA() << " " << getK() << " sin[" << getK() << "(x - ct)]" << endl;
      outputLOG << "A = " << getA() << " cos[" << getK() << "(x - ct)]" << endl;
    }
    else outputLOG << "Evoluzione non riconosciuta" << endl;

    inputstreamParameters.close();
    outputLOG << "Lettura parametri completata!" << endl;
    return 1;
  }
  else
  {
    outputLOG << "Errore, impossibile leggere dal file " << nomefile << endl;
    return -254;
  }
}


void Data::setInputType(char* nomefile)
{
  cout << "Il file " << nomefile << " contiene parametri o particelle? Specificare:\n1 se e' un file di parametri\n2 se e' un file di particelle \n:";
  while (inputType != 1 && inputType != 2)
  {
    cin >> inputType;
    if (inputType != 1 && inputType != 2)
    {
      cin.clear();
      cin.ignore(numeric_limits<streamsize>::max(), '\n');
      cout << "Non valido!\n:";
    }
  }
}


int Data::getInputType()
{
  return inputType;
}


void Data::fillData()
{
  cout << "Da implementare" << endl;
}

int Data::getNelectrons()
{
  return n_electrons;
}


void Data::setNelectrons_UNSECURED(int n_e)
{
  n_electrons = n_e;
  n_ions = n_e;
}


void Data::setNelectrons()
{
  while (n_electrons < 1)
  {
    cout << "Inserire numero particelle da simulare: ";
    cin >> n_electrons;
    if (n_electrons < 1)
    {
      cin.clear();
      cin.ignore(numeric_limits<streamsize>::max(), '\n');
      cout << "Non valido!" << endl;
    }
  }
  n_ions = n_electrons;
}


int Data::getNdim()
{
  return n_dim;
}


void Data::setNdim_UNSECURED(int n_dims)
{
  n_dim = n_dims;
}


void Data::setNdim()
{
  while (n_dim != 1 && n_dim != 2 && n_dim != 3)
  {
    cout << "1: Simulazione 1D; 2: Simulazione 2D; 3: Simulazione 3D. Scelta: ";
    cin >> n_dim;
    if (n_dim != 1 && n_dim != 2 && n_dim != 3)
    {
      cin.clear();
      cin.ignore(numeric_limits<streamsize>::max(), '\n');
      cout << "Non valido!" << endl;
    }
  }
}


double Data::getThermalV()
{
  return thermal_vel;

}


void Data::setThermalV_UNSECURED(double t_vel)
{
  thermal_vel = t_vel;
}


void Data::setThermalV(double thermal_vel_min)
{
  while (thermal_vel < thermal_vel_min)
  {
    cout << "Inserisci la temperatura media delle particelle: ";
    cin >> thermal_vel;
    if (thermal_vel < thermal_vel_min)
    {
      cin.clear();
      cin.ignore(numeric_limits<streamsize>::max(), '\n');
      cout << "Non valido!" << endl;
    }
  }

}


void Data::setXaxis(int nGridPointsX_min)
{

  while (dimX <= 0.)
  {
    cout << "Consideriamo quindi l'asse x. Quanto e' lungo il sistema? ";
    cin >> dimX;
    if (dimX <= 0.)
    {
      cin.clear();
      cin.ignore(numeric_limits<streamsize>::max(), '\n');
      cout << "Non valido!" << endl;
    }
  }

  while (nGridPointsX < nGridPointsX_min)
  {
    cout << "Quanti punti di griglia? ";
    cin >> nGridPointsX;
    if (nGridPointsX < nGridPointsX_min)
    {
      cin.clear();
      cin.ignore(numeric_limits<streamsize>::max(), '\n');
      cout << "Non valido!" << endl;
    }
  }

  deltaX = (dimX / (nGridPointsX - 2));
  if (n_dim == 1)
  {
    dimZ = dimY = dimX;
    deltaZ = deltaY = deltaX;
    nGridPointsZ = nGridPointsY = nGridPointsX;
  }
}


void Data::setXaxis_UNSECURED(double dimensioneX, double dx)
{
  dimX = dimensioneX;
  deltaX = dx;
  nGridPointsX = (int)((dimX / deltaX) + 2);
  if (n_dim == 1)
  {
    dimZ = dimY = dimX;
    deltaZ = deltaY = deltaX;
    nGridPointsZ = nGridPointsY = nGridPointsX;
  }
}


void Data::setXaxis_NOGRID_UNSECURED(double dimensioneX)
{
  dimX = dimensioneX;
}


void Data::setYaxis(int nGridPointsY_min)
{

  while (dimY <= 0.)
  {
    cout << "Consideriamo quindi l'asse y. Quanto e' lungo il sistema? ";
    cin >> dimY;
    if (dimY <= 0.)
    {
      cin.clear();
      cin.ignore(numeric_limits<streamsize>::max(), '\n');
      cout << "Non valido!" << endl;
    }
  }

  while (nGridPointsY < nGridPointsY_min)
  {
    cout << "Quanti punti di griglia? ";
    cin >> nGridPointsY;
    if (nGridPointsY < nGridPointsY_min)
    {
      cin.clear();
      cin.ignore(numeric_limits<streamsize>::max(), '\n');
      cout << "Non valido!" << endl;
    }
  }
  deltaY = (dimY / (nGridPointsY - 2));
  if (n_dim == 2)
  {
    dimZ = dimY;
    deltaZ = deltaY;
    nGridPointsZ = nGridPointsY;
  }
}


void Data::setYaxis_UNSECURED(double dimensioneY, double dy)
{
  dimY = dimensioneY;
  deltaY = dy;
  nGridPointsY = (int)((dimY / deltaY) + 2);
  if (n_dim == 2)
  {
    dimZ = dimY;
    deltaZ = deltaY;
    nGridPointsZ = nGridPointsY;
  }
}


void Data::setYaxis_NOGRID_UNSECURED(double dimensioneY)
{
  dimY = dimensioneY;
}


void Data::setZaxis_UNSECURED(double dimensioneZ, double dz)
{
  dimZ = dimensioneZ;
  deltaZ = dz;
  nGridPointsZ = (int)((dimZ / deltaZ) + 2);
}


void Data::setZaxis_NOGRID_UNSECURED(double dimensioneZ)
{
  dimZ = dimensioneZ;
}


void Data::setZaxis(int nGridPointsZ_min)
{

  while (dimZ <= 0.)
  {
    cout << "Consideriamo quindi l'asse z. Quanto e' lungo il sistema? ";
    cin >> dimZ;
    if (dimZ <= 0.)
    {
      cin.clear();
      cin.ignore(numeric_limits<streamsize>::max(), '\n');
      cout << "Non valido!" << endl;
    }
  }

  while (nGridPointsZ < nGridPointsZ_min)
  {
    cout << "Quanti punti di griglia? ";
    cin >> nGridPointsZ;
    if (nGridPointsZ < nGridPointsZ_min)
    {
      cin.clear();
      cin.ignore(numeric_limits<streamsize>::max(), '\n');
      cout << "Non valido!" << endl;
    }
  }
  deltaZ = (dimZ / (nGridPointsZ - 2));
}


int Data::getNgridPointsX()
{
  return nGridPointsX;
}


int Data::getNgridPointsY()
{
  return nGridPointsY;
}


int Data::getNgridPointsZ()
{
  return nGridPointsZ;
}


int Data::getEvoType()
{
  return evoType;
}


void Data::setEvoType_UNSECURED(int et)
{
  evoType = et;
}


void Data::setEvoType()
{
  while (evoType != 1 && evoType != 2 && evoType != 3 && evoType != 4 && evoType != 5 && evoType != 6 && evoType != 7 && evoType != 8 && evoType != 9 && evoType != 10 && evoType != 11 && evoType != 12)
  {
    cout << "Ecco cosa puoi fare:" << endl;
    cout << "1-  Distribuzione campi su particelle in un istante definito" << endl;
    cout << "2-  Evoluzione campo su una particella nel tempo" << endl;
    cout << "3-  Evoluzione RungeKutta4 no griglia (campo analitico calcolato sulla particella)" << endl;
    cout << "4-  Evoluzione RungeKutta4 griglia on-the-fly" << endl;
    cout << "5-  Evoluzione RungeKutta4 griglia persistente" << endl;
    cout << "6-  Evoluzione Leapfrog2 no griglia" << endl;
    cout << "7-  Evoluzione Leapfrog2 griglia on-the-fly" << endl;
    cout << "8-  Evoluzione Leapfrog2 griglia persistente" << endl;
    cout << "9-  Evoluzione Leapfrog4 no griglia" << endl;
    cout << "10- Evoluzione Leapfrog4 griglia on-the-fly" << endl;
    cout << "11- Evoluzione Leapfrog4 griglia persistente" << endl;
    cout << "12- Evoluzione Leapfrog4 no griglia - algoritmo TURCHETTI" << endl;
    cout << "--------------------------------------------" << endl;
    cout << "Esprimi una scelta: ";
    cin >> evoType;
    if (evoType != 1 && evoType != 2 && evoType != 3 && evoType != 4 && evoType != 5 && evoType != 6 && evoType != 7 && evoType != 8 && evoType != 9 && evoType != 10 && evoType != 11 && evoType != 12)
    {
      cin.clear();
      cin.ignore(numeric_limits<streamsize>::max(), '\n');
      cout << "Non valido!" << endl;
      cout << "--------------------------------------------" << endl;
      cout << "============================================" << endl;

    }
  }
}


/*
void Data::setAnalyticType()
{
  while (analyticType != 1 && analyticType != 2 && analyticType != 3)
  {
    cout << "Cosa facciamo?" << endl;
    cout << "1- Distribuzione campi su particelle in un istante definito" << endl;
    cout << "2- Evoluzione campo su una particella nel tempo" << endl;
    cout << "3- Evoluzione particelle metodo RK4 con calcolo campo analitico sulla particella" << endl;
    cout << "Scegli: ";
    cin >> analyticType;
    if (analyticType != 1 && analyticType != 2 && analyticType != 3)
    {
      cin.clear();
      cin.ignore(numeric_limits<streamsize>::max(), '\n');
      cout << "Non valido!" << endl;
    }
  }
}


void Data::setAnalyticType_UNSECURED(int at)
{
  analyticType = at;
}


int Data::getAnalyticType()
{
  return analyticType;
}
*/


int Data::getNsteps()
{
  return nSteps;
}


void Data::setNsteps_UNSECURED(int ns)
{
  nSteps = ns;
}


void Data::setNsteps(int nSteps_min)
{
  while (nSteps < nSteps_min)
  {
    cout << "Inserisci il numero di step evolutivi che il sistema deve compiere: ";
    cin >> nSteps;
    if (nSteps < nSteps_min)
    {
      cin.clear();
      cin.ignore(numeric_limits<streamsize>::max(), '\n');
      cout << "Non valido!" << endl;
    }
  }
}


double Data::getDeltaT()
{
  return dt;
}


void Data::setDeltaT_UNSECURED(double deltaT)
{
  dt = deltaT;
}


void Data::setDeltaT(double dt_min)
{
  if (dt_min == 0.)
  {
    cout << "Inserisci la durata temporale dello step evolutivo (dt): ";
    cin >> dt;
    while (cin.fail() || dt <= 0.)
    {
      if (cin.fail())
      {
        cin.clear();
        cin.ignore(numeric_limits<streamsize>::max(), '\n');
      }
      cout << "Non valido!\n: ";
      cin >> dt;
    }
  }
  else
  {
    while (dt <= dt_min)
    {
      cout << "Inserisci la durata temporale dello step evolutivo (dt): ";
      cin >> dt;
      if (dt <= dt_min)
      {
        cin.clear();
        cin.ignore(numeric_limits<streamsize>::max(), '\n');
        cout << "Non valido!" << endl;
      }
    }
  }
}


int Data::getParticleFillingMethod()
{
  return particleFillingMethod;
}


void Data::setParticleFillingMethod_UNSECURED(int pfm)
{
  particleFillingMethod = pfm;
}


void Data::setParticleFillingMethod()
{
  while (particleFillingMethod != 0 && particleFillingMethod != 1 && particleFillingMethod != 2)
  {
    cout << "Particelle:" << endl;
    cout << "0: Inserimento manuale" << endl;
    cout << "1: Generazione random in posizione con momento definito dalla v termica" << endl;
    cout << "2: Tutte le particelle sono generate in 0,0,0 con momento definito dalla v termica" << endl;
    cout << "Scelta: ";
    cin >> particleFillingMethod;
    if (particleFillingMethod != 0 && particleFillingMethod != 1 && particleFillingMethod != 2)
    {
      cin.clear();
      cin.ignore(numeric_limits<streamsize>::max(), '\n');
      cout << "Non valido!" << endl;
    }
  }
}


double Data::getDimX()
{
  return dimX;
}


double Data::getDimY()
{
  return dimY;
}


double Data::getDimZ()
{
  return dimZ;
}


double Data::getDeltaX()
{
  return deltaX;
}


double Data::getDeltaY()
{
  return deltaY;
}


double Data::getDeltaZ()
{
  return deltaZ;
}


void Data::setK()
{
  if (A == 0) cout << "Ex = - A k sin[k(x - ct)]\n Inserisci k: ";
  else cout << " Inserisci k: ";
  cin >> k;
  while (cin.fail())
  {
    cin.clear();
    cin.ignore(numeric_limits<streamsize>::max(), '\n');
    cout << "Non valido!\n: " << endl;
    cin >> k;
  }
}


void Data::setK_UNSECURED(double K)
{
  k = K;
}


double Data::getK()
{
  return k;
}


void Data::setA()
{
  if (k == 0) cout << "Ex = - A k sin[k(x - ct)]\n Inserisci A: ";
  else cout << " Inserisci A: ";
  cin >> A;
  while (cin.fail())
  {
    cin.clear();
    cin.ignore(numeric_limits<streamsize>::max(), '\n');
    cout << "Non valido!\n: " << endl;
    cin >> A;
  }
}


void Data::setA_UNSECURED(double a)
{
  A = a;
}


double Data::getA()
{
  return A;
}


void Data::aumentaT(double deltaT)
{
  t += deltaT;
}


void Data::resettaT()
{
  t = 0.;
}


void Data::impostaT(double tempo)
{
  if (tempo == 0)
  {
    cout << "A che istante vuoi porre il tempo t? ";
    cin >> tempo;
    while (cin.fail())
    {
      cin.clear();
      cin.ignore(numeric_limits<streamsize>::max(), '\n');
      cin >> tempo;
      if (cin.fail())
      {
        cout << "Non valido!\n: " << endl;
      }
    }
  }
  t = tempo;
}


void Data::impostaT_UNSECURED(double tempo)
{
  t = tempo;
}


double Data::getT()
{
  return t;
}


Field::Field()
{
  ex = ey = ez = e4 = bx = by = bz = b4 = 0.;
}


void Field::setEx(double Ex)
{
  ex = Ex;
}


void Field::setEy(double Ey)
{
  ey = Ey;
}


void Field::setEz(double Ez)
{
  ez = Ez;
}


double Field::getEx()
{
  return ex;
}


double Field::getEy()
{
  return ey;
}


double Field::getEz()
{
  return ez;
}


void Field::setBx(double Bx)
{
  bx = Bx;
}


void Field::setBy(double By)
{
  by = By;
}


void Field::setBz(double Bz)
{
  bz = Bz;
}


double Field::getBx()
{
  return bx;
}


double Field::getBy()
{
  return by;
}


double Field::getBz()
{
  return bz;
}


Field::~Field() { }


Particle::Particle()
{
  x = y = z = t = px = py = pz = pt = 0.;

  m = q = cell = 0;
}


void Particle::setParticleX(double X)
{
  x = X;
}


void Particle::setParticleX_fromUser(double xmin, double xmax)
{
  cout << "Inserisci la coordinata x della particella: ";
  cin >> x;
  while (x < xmin || x > xmax)
  {
    cin.clear();
    cin.ignore(numeric_limits<streamsize>::max(), '\n');
    cout << "Non valido\n:";
    cin >> x;
  }
}


void Particle::setParticleY(double Y)
{
  y = Y;
}


void Particle::setParticleY_fromUser(double ymin, double ymax)
{
  cout << "Inserisci la coordinata y della particella: ";
  cin >> y;
  while (y < ymin || y > ymax)
  {
    cin.clear();
    cin.ignore(numeric_limits<streamsize>::max(), '\n');
    cout << "Non valido\n:";
    cin >> y;
  }
}


void Particle::setParticleZ(double Z)
{
  z = Z;
}


void Particle::setParticleZ_fromUser(double zmin, double zmax)
{
  cout << "Inserisci la coordinata z della particella: ";
  cin >> z;
  while (z < zmin || z > zmax)
  {
    cin.clear();
    cin.ignore(numeric_limits<streamsize>::max(), '\n');
    cout << "Non valido\n:";
    cin >> z;
  }
}


void Particle::setParticlePX(double PX)
{
  px = PX;
}


void Particle::setParticlePX_fromUser(double pmin)
{
  do
  {
    cin.clear();
    cin.ignore(numeric_limits<streamsize>::max(), '\n');

    cout << "Inserisci il momento lungo x della particella: ";
    cin >> px;
    if (cin.fail() || px < pmin)
    {
      cout << "Non valido oppure minore di " << pmin << " (impulso minimo)" << endl;
    }
  } while (cin.fail() || px < pmin);
}


void Particle::setParticlePY(double PY)
{
  py = PY;
}


void Particle::setParticlePY_fromUser(double pmin)
{
  do
  {
    cin.clear();
    cin.ignore(numeric_limits<streamsize>::max(), '\n');

    cout << "Inserisci il momento lungo y della particella: ";
    cin >> py;
    if (cin.fail() || px < pmin)
    {
      cout << "Non valido oppure minore di " << pmin << " (impulso minimo)" << endl;
    }
  } while (cin.fail() || px < pmin);
}


void Particle::setParticlePZ(double PZ)
{
  pz = PZ;
}


void Particle::setParticlePZ_fromUser(double pmin)
{
  do
  {
    cin.clear();
    cin.ignore(numeric_limits<streamsize>::max(), '\n');

    cout << "Inserisci il momento lungo z della particella: ";
    cin >> pz;
    if (cin.fail() || px < pmin)
    {
      cout << "Non valido oppure minore di " << pmin << " (impulso minimo)" << endl;
    }
  } while (cin.fail() || px < pmin);
}


void Particle::setParticleQ(int Q)
{
  q = Q;
}


void Particle::setParticleM(int mass)
{
  m = mass;
}


void Particle::setParticleCell(Data data)
{
  int posZ = 0;
  int posY = 0;
  int posX = 0;
  if (data.getDeltaZ()) posZ = (int)(z / data.getDeltaZ());
  if (data.getDeltaY()) posY = (int)(y / data.getDeltaY());
  if (data.getDeltaX()) posX = (int)(x / data.getDeltaX());

  cell = posZ * data.getNgridPointsX() * data.getNgridPointsY() + posY * data.getNgridPointsX() + posX + 1;
}


double Particle::getParticleX()
{
  return x;
}


double Particle::getParticleY()
{
  return y;
}


double Particle::getParticleZ()
{
  return z;
}


double Particle::getParticlePX()
{
  return px;
}


double Particle::getParticlePY()
{
  return py;
}


double Particle::getParticlePZ()
{
  return pz;
}


int Particle::getParticleQ()
{
  return q;
}


int Particle::getParticleM()
{
  return m;
}


int Particle::getParticleCell()
{
  return cell;
}

/*
Particle& Particle::operator*(const double value);
{

}
*/


Particle::~Particle() { }
