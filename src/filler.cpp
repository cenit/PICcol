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



#define C 1.00

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327
#endif

#include <iostream>
#include <cstdio>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <fstream>

#include "filler.h"

using namespace std;


void fillParticlesFromFile(char* nomefile, ofstream &outputLOG)     // NON E' ANCORA STATA IMPLEMENTATA!
{
  ifstream inputstreamParticles;
  inputstreamParticles.open(nomefile);
  if (inputstreamParticles.is_open())
  {
    outputLOG << "Lettura coordinate particelle da " << nomefile << " in corso..." << endl;
    outputLOG << "Funzione non implementata" << endl;
    cout << "Funzione non implementata" << endl;
    inputstreamParticles.close();
    outputLOG << "Lettura parametri completata!" << endl;
  }

  else
  {
    outputLOG << "Errore, impossibile leggere dal file " << nomefile << endl;
  }
}


void calcolaCampoAnalitico1DSuParticelle(Data data, vector<Particle> &particelle, vector<Field> &campoSuParticelle)
{
  Field *campoTemp;
  campoTemp = new Field[1];

  for (int i = 0; i < data.getNelectrons(); i++)
  {
    calcolaCampoAnalitico1DSuParticella(data, particelle[i], *campoTemp);
    campoSuParticelle.push_back(*campoTemp);
  }
}


void calcolaCampoAnalitico1DSuGriglia(Data data, vector<Particle> &posizioni, vector<Field> &campoSuPuntiGriglia)
{
  Field *campoTemp;
  campoTemp = new Field[1];

  int size = posizioni.size();

  for (int i = 0; i < size; i++)
  {
    calcolaCampoAnalitico1DSuParticella(data, posizioni[i], *campoTemp);
    campoSuPuntiGriglia.push_back(*campoTemp);
  }
}


void riempiPuntiGrigliaConCampoAnalitico1D(Data data, vector<Field> &campoSuPuntiGriglia)
{
  Field *campoTemp;
  campoTemp = new Field[1];

  Particle particellaTemp;
  //  for (int k = 0; k < data.getNgridPointsZ(); k++)
  //  {
  for (int j = 0; j < data.getNgridPointsY(); j++)
  {
    for (int i = 0; i < data.getNgridPointsX(); i++)
    {
      particellaTemp.setParticleX((i - 1)*data.getDeltaX());
      particellaTemp.setParticleY((j - 1)*data.getDeltaY());
      //        particellaTemp.setParticleZ((k-1)*data.getDeltaZ());
      particellaTemp.setParticleZ(0.);
      calcolaCampoAnalitico1DSuParticella(data, particellaTemp, *campoTemp);
      campoSuPuntiGriglia.push_back(*campoTemp);
    }
  }
  //  }
}


void riempiPuntiGrigliaConPotenzialeAnalitico1D(Data data, vector<Field> &campoSuPuntiGriglia)
{
  Field *campoTemp;
  campoTemp = new Field[1];

  Particle particellaTemp;
  //  for (int k = 0; k < data.getNgridPointsZ(); k++)
  //  {
  for (int j = 0; j < data.getNgridPointsY(); j++)
  {
    for (int i = 0; i < data.getNgridPointsX(); i++)
    {
      particellaTemp.setParticleX((i - 1)*data.getDeltaX());
      particellaTemp.setParticleY((j - 1)*data.getDeltaY());
      //        particellaTemp.setParticleZ((k-1)*data.getDeltaZ());
      particellaTemp.setParticleZ(0.);
      calcolaPotenzialeAnalitico1DSuParticella(data, particellaTemp, *campoTemp);
      campoSuPuntiGriglia.push_back(*campoTemp);
    }
  }
  //  }
}


void calcolaCampoAnalitico1DSuParticella(Data data, Particle particella, Field &campoSuParticella)
{
  //  /****************************************************************
    // PACCHETTO D'ONDA
  double f = 0., Amp = 0.;
  double tau = 5. * 2. * M_PI / data.getK();
  double x0 = -1.1 * tau;
  double arg = (particella.getParticleX() - C * data.getT() - x0) / tau;
  if (abs(arg) < 1)
  {
    Amp = pow(cos(arg*0.5*M_PI), 2);
    f = -data.getK() * sin(data.getK() * (particella.getParticleX() - C * data.getT())) * Amp
      - (0.5 * M_PI / tau) * cos(data.getK() * (particella.getParticleX() - C * data.getT())) * sin(M_PI * (particella.getParticleX() - C * data.getT() - x0) / tau);
  }
  f *= data.getA();

  campoSuParticella.setEx(0.);
  campoSuParticella.setEz(0.);
  campoSuParticella.setBx(0.);
  campoSuParticella.setBy(0.);

  campoSuParticella.setEy(f);
  campoSuParticella.setBz(f);
  //  ****************************************************************/

    /****************************************************************
    // PACCHETTO D'ONDA (TURCHETTI)
    double Amp, f=0., fp=0.;
    double x_foc = M_PI;
    double Ak = 2*M_PI;
    double ell = 5.;
    double beta = M_PI/(2*ell);  // non la capisco e quindi l'ho commentata e ridefinita pari a 1 sotto!
    double arg = beta * (particella.getParticleX() - x_foc - data.getT());
    if (abs(arg) < (0.5*M_PI))
    {
      Amp = pow(cos(arg),2);
      f = cos (Ak * (particella.getParticleX() - data.getT())) * Amp;
      fp = -Ak * sin (Ak * (particella.getParticleX() - data.getT())) * Amp;
      fp -= beta * cos (Ak * (particella.getParticleX() - data.getT())) * sin(2*arg);
    }
    f *= data.getA();
    fp *= data.getA();

    campoSuParticella.setEx(0.);
    campoSuParticella.setEz(0.);
    campoSuParticella.setBx(0.);
    campoSuParticella.setBy(0.);

    campoSuParticella.setEy(f);
    campoSuParticella.setBz(f);
    ****************************************************************/

    /****************************************************************
    // ONDA PIANA INFINITA
    campoSuParticella.setEx(0.);
    campoSuParticella.setEz(0.);
    campoSuParticella.setBx(0.);
    campoSuParticella.setBy(0.);

    campoSuParticella.setEy(- data.getK() * data.getA() * sin (data.getK() * (particella.getParticleX() - C * data.getT())));
    campoSuParticella.setBz(campoSuParticella.getEy());
    ****************************************************************/

    /****************************************************************
    // ONDA PIANA INFINITA CONFORME ALLE SCELTE DELLO SGATTO
    campoSuParticella.setEy(0.);
    campoSuParticella.setEz(0.);
    campoSuParticella.setBx(0.);
    campoSuParticella.setBz(0.);

    campoSuParticella.setEx(- data.getK() * data.getA() * sin (data.getK() * (particella.getParticleZ() - C * data.getT())));
    campoSuParticella.setBy(campoSuParticella.getEx());
    ****************************************************************/


    /****************************************************************
    // SOLO CAMPO MAGNETICO LUNGO Y
    campoSuParticella.setEy(0.);
    campoSuParticella.setBz(1.);
    ****************************************************************/
}


void calcolaPotenzialeAnalitico1DSuParticella(Data data, Particle particella, Field &potenzialeSuParticella)
{
  //  /****************************************************************
    // PACCHETTO D'ONDA
  double f = 0., Amp = 0.;
  double tau = 5. * 2. * M_PI / data.getK();
  double x0 = -1.1 * tau;
  double arg = (particella.getParticleX() - C * data.getT() - x0) / tau;
  if (abs(arg) < 1)
  {
    Amp = pow(cos(arg*0.5*M_PI), 2);
    f = cos(data.getK() * (particella.getParticleX() - C * data.getT())) * Amp;
  }
  f *= data.getA();

  potenzialeSuParticella.setEx(0.);
  potenzialeSuParticella.setEz(0.);
  potenzialeSuParticella.setBx(0.);
  potenzialeSuParticella.setBy(0.);
  potenzialeSuParticella.setBz(0.);

  potenzialeSuParticella.setEy(f);
  //  ****************************************************************/

    /****************************************************************
    // POTENZIALE PER PACCHETTO D'ONDA PIANA INFINITA (TURCHETTI)
    double Amp, f=0., fp=0.;
    double x_foc = M_PI;
    double Ak = 2*M_PI;
    double ell = 5.;
    double beta = M_PI/(2*ell);  // non la capisco e quindi l'ho commentata e ridefinita pari a 1 sotto!
    double arg = beta * (particella.getParticleX() - x_foc - data.getT());
    if (abs(arg) < (0.5*M_PI))
    {
      Amp = pow(cos(arg),2);
      f = cos (Ak * (particella.getParticleX() - data.getT())) * Amp;
      fp = -Ak * sin (Ak * (particella.getParticleX() - data.getT())) * Amp;
      fp -= beta * cos (Ak * (particella.getParticleX() - data.getT())) * sin(2*arg);
    }
    f *= data.getA();
    fp *= data.getA();

    potenzialeSuParticella.setEy(fp);
    potenzialeSuParticella.setEx(0.);
    potenzialeSuParticella.setEz(0.);
    potenzialeSuParticella.setBx(0.);   // nb: le 3 componenti del campo vettoriale le ho memorizzate in E; le tre componenti B di questo tipo di dato non sono utilizzate ma comunque poste a zero
    potenzialeSuParticella.setBy(0.);
    potenzialeSuParticella.setBz(0.);
    ****************************************************************/

    /****************************************************************
    // POTENZIALE PER ONDA PIANA INFINITA
    potenzialeSuParticella.setEy(data.getA() * cos (data.getK() * (particella.getParticleX() - C * data.getT())));
    potenzialeSuParticella.setEx(0.);
    potenzialeSuParticella.setEz(0.);
    potenzialeSuParticella.setBx(0.);   // nb: le 3 componenti del campo vettoriale le ho memorizzate in E; le tre componenti B di questo tipo di dato non sono utilizzate ma comunque poste a zero
    potenzialeSuParticella.setBy(0.);
    potenzialeSuParticella.setBz(0.);
    ****************************************************************/

    /****************************************************************
    // POTENZIALE PER ONDA PIANA INFINITA CONFORME ALLE SCELTE DELLO SGATTO
    potenzialeSuParticella.setEx(data.getA() * cos (data.getK() * (particella.getParticleZ() - C * data.getT())));
    potenzialeSuParticella.setEy(0.);
    potenzialeSuParticella.setEz(0.);
    potenzialeSuParticella.setBx(0.);   // nb: le 3 componenti del campo vettoriale le ho memorizzate in E; le tre componenti B di questo tipo di dato non sono utilizzate ma comunque poste a zero
    potenzialeSuParticella.setBy(0.);
    potenzialeSuParticella.setBz(0.);
    ****************************************************************/

}


void creaVettoreParticelle(Data data, vector<Particle> &particelle)
{
  Particle particellaTest;
  double normalizza = (1. / RAND_MAX);
  double randomX, randomY, randomZ;

  if (data.getParticleFillingMethod() == 0)
  {
    for (int i = 0; i < data.getNelectrons(); i++)
    {
      particellaTest.setParticleX_fromUser(0., data.getDimX());
      if (data.getNdim() == 2 || data.getNdim() == 3) particellaTest.setParticleY_fromUser(0., data.getDimY());
      else particellaTest.setParticleY(0.);
      if (data.getNdim() == 3) particellaTest.setParticleZ_fromUser(0., data.getDimZ());
      else particellaTest.setParticleZ(0.);

      if (data.getThermalV() == 0)
      {
        particellaTest.setParticlePX_fromUser(0.);
        if (data.getNdim() == 2 || data.getNdim() == 3) particellaTest.setParticlePY_fromUser(0.);
        else particellaTest.setParticlePY(0.);
        if (data.getNdim() == 3) particellaTest.setParticlePZ_fromUser(0.);
        else particellaTest.setParticlePZ(0.);
      }
      else
      {
        randomX = rand() * normalizza * data.getThermalV();
        particellaTest.setParticlePX(randomX);
        randomY = rand() * normalizza * data.getThermalV();
        particellaTest.setParticlePY(randomY);
        randomZ = rand() * normalizza * data.getThermalV();
        particellaTest.setParticlePZ(randomZ);
      }

      particellaTest.setParticleQ(1);
      particellaTest.setParticleM(1);
      particelle.push_back(particellaTest);
    }
  }
  else if (data.getParticleFillingMethod() == 1)
  {
    for (int i = 0; i < data.getNelectrons(); i++)
    {
      randomX = rand() * normalizza * data.getDimX();
      particellaTest.setParticleX(randomX);
      if (data.getNdim() == 2 || data.getNdim() == 3)
      {
        randomY = rand() * normalizza * data.getDimY();
        particellaTest.setParticleY(randomY);
      }
      else particellaTest.setParticleY(0.);
      if (data.getNdim() == 3)
      {
        randomZ = rand() * normalizza * data.getDimZ();
        particellaTest.setParticleZ(randomZ);
      }
      else particellaTest.setParticleZ(0.);

      randomX = rand() * normalizza * data.getThermalV();
      particellaTest.setParticlePX(randomX);
      if (data.getNdim() == 2 || data.getNdim() == 3)
      {
        randomY = rand() * normalizza * data.getThermalV();
        particellaTest.setParticlePY(randomY);
      }
      if (data.getNdim() == 3)
      {
        randomZ = rand() * normalizza * data.getThermalV();
        particellaTest.setParticlePZ(randomZ);
      }

      particellaTest.setParticleQ(1);
      particellaTest.setParticleM(1);
      particelle.push_back(particellaTest);
    }
  }
  else if (data.getParticleFillingMethod() == 2)
  {
    for (int i = 0; i < data.getNelectrons(); i++)
    {
      particellaTest.setParticleX(0.);
      particellaTest.setParticleY(0.);
      particellaTest.setParticleZ(0.);

      randomX = rand() * normalizza * data.getThermalV();
      particellaTest.setParticlePX(randomX);
      if (data.getNdim() == 2 || data.getNdim() == 3)
      {
        randomY = rand() * normalizza * data.getThermalV();
        particellaTest.setParticlePY(randomY);
      }
      if (data.getNdim() == 3)
      {
        randomZ = rand() * normalizza * data.getThermalV();
        particellaTest.setParticlePZ(randomZ);
      }


      particellaTest.setParticleQ(1);
      particellaTest.setParticleM(1);
      particelle.push_back(particellaTest);
    }
  }

}
