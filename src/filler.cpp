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

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>

#include "filler.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327
#endif

constexpr double SPEED_OF_LIGHT = 1.0;

using namespace std;

void calcolaCampoAnalitico1DSuParticelle(Data data, vector<Particle>& particelle, vector<Field>& campoSuParticelle)
{
    Field* campoTemp;
    campoTemp = new Field[1];

    for (int i = 0; i < data.n_electrons; i++) {
        calcolaCampoAnalitico1DSuParticella(data, particelle[i], *campoTemp);
        campoSuParticelle.push_back(*campoTemp);
    }
}

void calcolaCampoAnalitico1DSuGriglia(Data data, vector<Particle>& posizioni, vector<Field>& campoSuPuntiGriglia)
{
    Field* campoTemp;
    campoTemp = new Field[1];

    int size = posizioni.size();

    for (int i = 0; i < size; i++) {
        calcolaCampoAnalitico1DSuParticella(data, posizioni[i], *campoTemp);
        campoSuPuntiGriglia.push_back(*campoTemp);
    }
}

void riempiPuntiGrigliaConCampoAnalitico1D(Data data, vector<Field>& campoSuPuntiGriglia)
{
    Field* campoTemp;
    campoTemp = new Field[1];

    Particle particellaTemp;
    //  for (int k = 0; k < data.nGridPointsZ; k++)
    //  {
    for (int j = 0; j < data.nGridPointsY; j++) {
        for (int i = 0; i < data.nGridPointsX; i++) {
            particellaTemp.setParticleX((i - 1) * data.deltaX);
            particellaTemp.setParticleY((j - 1) * data.deltaY);
            //        particellaTemp.setParticleZ((k-1)*data.deltaZ);
            particellaTemp.setParticleZ(0.);
            calcolaCampoAnalitico1DSuParticella(data, particellaTemp, *campoTemp);
            campoSuPuntiGriglia.push_back(*campoTemp);
        }
    }
    //  }
}

void riempiPuntiGrigliaConPotenzialeAnalitico1D(Data data, vector<Field>& campoSuPuntiGriglia)
{
    Field* campoTemp;
    campoTemp = new Field[1];

    Particle particellaTemp;
    //  for (int k = 0; k < data.nGridPointsZ; k++)
    //  {
    for (int j = 0; j < data.nGridPointsY; j++) {
        for (int i = 0; i < data.nGridPointsX; i++) {
            particellaTemp.setParticleX((i - 1) * data.deltaX);
            particellaTemp.setParticleY((j - 1) * data.deltaY);
            //        particellaTemp.setParticleZ((k-1)*data.deltaZ);
            particellaTemp.setParticleZ(0.);
            calcolaPotenzialeAnalitico1DSuParticella(data, particellaTemp, *campoTemp);
            campoSuPuntiGriglia.push_back(*campoTemp);
        }
    }
    //  }
}

void calcolaCampoAnalitico1DSuParticella(Data data, Particle particella, Field& campoSuParticella)
{
    //  /****************************************************************
    // PACCHETTO D'ONDA
    double f = 0., Amp = 0.;
    double tau = 5. * 2. * M_PI / data.k;
    double x0 = -1.1 * tau;
    double arg = (particella.x - SPEED_OF_LIGHT * data.t - x0) / tau;
    if (abs(arg) < 1) {
        Amp = pow(cos(arg * 0.5 * M_PI), 2);
        f = -data.k * sin(data.k * (particella.x - SPEED_OF_LIGHT * data.t)) * Amp
            - (0.5 * M_PI / tau) * cos(data.k * (particella.x - SPEED_OF_LIGHT * data.t)) * sin(M_PI * (particella.x - SPEED_OF_LIGHT * data.t - x0) / tau);
    }
    f *= data.A;

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
    double arg = beta * (particella.x - x_foc - data.t);
    if (abs(arg) < (0.5*M_PI))
    {
      Amp = pow(cos(arg),2);
      f = cos (Ak * (particella.x - data.t)) * Amp;
      fp = -Ak * sin (Ak * (particella.x - data.t)) * Amp;
      fp -= beta * cos (Ak * (particella.x - data.t)) * sin(2*arg);
    }
    f *= data.A;
    fp *= data.A;

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

    campoSuParticella.setEy(- data.k * data.A * sin (data.k * (particella.x - SPEED_OF_LIGHT * data.t)));
    campoSuParticella.setBz(campoSuParticella.ey);
    ****************************************************************/

    /****************************************************************
    // ONDA PIANA INFINITA CONFORME ALLE SCELTE DELLO SGATTO
    campoSuParticella.setEy(0.);
    campoSuParticella.setEz(0.);
    campoSuParticella.setBx(0.);
    campoSuParticella.setBz(0.);

    campoSuParticella.setEx(- data.k * data.A * sin (data.k * (particella.z - SPEED_OF_LIGHT * data.t)));
    campoSuParticella.setBy(campoSuParticella.ex);
    ****************************************************************/

    /****************************************************************
    // SOLO CAMPO MAGNETICO LUNGO Y
    campoSuParticella.setEy(0.);
    campoSuParticella.setBz(1.);
    ****************************************************************/
}

void calcolaPotenzialeAnalitico1DSuParticella(Data data, Particle particella, Field& potenzialeSuParticella)
{
    // nb: le 3 componenti del campo vettoriale vengono memorizzate in E; le tre componenti B non sono utilizzate ma comunque poste a zero

    // POTENZIALE PER PACCHETTO D'ONDA
    double f = 0., Amp = 0.;
    double tau = 5. * 2. * M_PI / data.k;
    double x0 = -1.1 * tau;
    double arg = (particella.x - SPEED_OF_LIGHT * data.t - x0) / tau;
    if (abs(arg) < 1) {
        Amp = pow(cos(arg * 0.5 * M_PI), 2);
        f = cos(data.k * (particella.x - SPEED_OF_LIGHT * data.t)) * Amp;
    }
    f *= data.A;

    potenzialeSuParticella.ey = f;
    potenzialeSuParticella.ex = 0.;
    potenzialeSuParticella.ez = 0.;
    potenzialeSuParticella.bx = 0.;
    potenzialeSuParticella.by = 0.;
    potenzialeSuParticella.bz = 0.;

    /****************************************************************
    // POTENZIALE PER PACCHETTO D'ONDA (TURCHETTI)
    double Amp, f=0., fp=0.;
    double x_foc = M_PI;
    double Ak = 2*M_PI;
    double ell = 5.;
    double beta = M_PI/(2*ell);
    double arg = beta * (particella.x - x_foc - data.t);
    if (abs(arg) < (0.5*M_PI))
    {
      Amp = pow(cos(arg),2);
      f = cos (Ak * (particella.x - data.t)) * Amp;
      fp = -Ak * sin (Ak * (particella.x - data.t)) * Amp;
      fp -= beta * cos (Ak * (particella.x - data.t)) * sin(2*arg);
    }
    f *= data.A;
    fp *= data.A;
    potenzialeSuParticella.setEy(fp);
    potenzialeSuParticella.setEx(0.);
    potenzialeSuParticella.setEz(0.);
    potenzialeSuParticella.setBx(0.);
    potenzialeSuParticella.setBy(0.);
    potenzialeSuParticella.setBz(0.);
    ****************************************************************/

    /****************************************************************
    // POTENZIALE PER ONDA PIANA INFINITA
    potenzialeSuParticella.setEy(data.A * cos (data.k * (particella.x - SPEED_OF_LIGHT * data.t)));
    potenzialeSuParticella.setEx(0.);
    potenzialeSuParticella.setEz(0.);
    potenzialeSuParticella.setBx(0.);
    potenzialeSuParticella.setBy(0.);
    potenzialeSuParticella.setBz(0.);
    ****************************************************************/

    /****************************************************************
    // POTENZIALE PER ONDA PIANA INFINITA CONFORME ALLE SCELTE DELLO SGATTO
    potenzialeSuParticella.setEx(data.A * cos (data.k * (particella.z - SPEED_OF_LIGHT * data.t)));
    potenzialeSuParticella.setEy(0.);
    potenzialeSuParticella.setEz(0.);
    potenzialeSuParticella.setBx(0.);
    potenzialeSuParticella.setBy(0.);
    potenzialeSuParticella.setBz(0.);
    ****************************************************************/
}

void creaVettoreParticelle(Data data, vector<Particle>& particelle)
{
    Particle particellaTest;
    double normalizza = (1. / RAND_MAX);
    double randomX, randomY, randomZ;

    if (data.particleFillingMethod == 1) {
        for (int i = 0; i < data.n_electrons; i++) {
            randomX = rand() * normalizza * data.dimX;
            particellaTest.x = randomX;
            if (data.n_dim == 2 || data.n_dim == 3) {
                randomY = rand() * normalizza * data.dimY;
                particellaTest.setParticleY(randomY);
            } else
                particellaTest.setParticleY(0.);
            if (data.n_dim == 3) {
                randomZ = rand() * normalizza * data.dimZ;
                particellaTest.setParticleZ(randomZ);
            } else
                particellaTest.setParticleZ(0.);

            randomX = rand() * normalizza * data.thermal_vel;
            particellaTest.setParticlePX(randomX);
            if (data.n_dim == 2 || data.n_dim == 3) {
                randomY = rand() * normalizza * data.thermal_vel;
                particellaTest.setParticlePY(randomY);
            }
            if (data.n_dim == 3) {
                randomZ = rand() * normalizza * data.thermal_vel;
                particellaTest.setParticlePZ(randomZ);
            }

            particellaTest.setParticleQ(1);
            particellaTest.setParticleM(1);
            particelle.push_back(particellaTest);
        }
    } else if (data.particleFillingMethod == 2) {
        for (int i = 0; i < data.n_electrons; i++) {
            particellaTest.setParticleX(0.);
            particellaTest.setParticleY(0.);
            particellaTest.setParticleZ(0.);

            randomX = rand() * normalizza * data.thermal_vel;
            particellaTest.setParticlePX(randomX);
            if (data.n_dim == 2 || data.n_dim == 3) {
                randomY = rand() * normalizza * data.thermal_vel;
                particellaTest.setParticlePY(randomY);
            }
            if (data.n_dim == 3) {
                randomZ = rand() * normalizza * data.thermal_vel;
                particellaTest.setParticlePZ(randomZ);
            }

            particellaTest.setParticleQ(1);
            particellaTest.setParticleM(1);
            particelle.push_back(particellaTest);
        }
    }
}
