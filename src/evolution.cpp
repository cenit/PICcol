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
#include <iomanip>
#include <ios>
#include <iostream>
#include <limits>

#include "evolution.h"
#include "filler.h"

using namespace std;
constexpr double SPEED_OF_LIGHT = 1.0;

void evolveRK4_nogrid(Data data, vector<Particle>& particles, vector<Field>& fieldOnParticles, ofstream& outputstream)
{
  Particle x2, dx, dx2;
  double gamma;
  double zeta;
  double I1, I2, I3;

  Field* campoSuPunto;
  campoSuPunto = new Field[1];
  vector<Field> potenziale;

  for (int j = 0; j < data.nSteps; j++) {
    for (int i = 0; i < data.n_electrons; i++) {
      gamma = sqrt(1 + pow(particles[i].px, 2) + pow(particles[i].py, 2) + pow(particles[i].pz, 2));
      fieldOnParticles.clear();
      calcolaCampoAnalitico1DSuParticella(data, particles[i], *campoSuPunto);
      fieldOnParticles.push_back(*campoSuPunto);

      // PASSO #1
      // RIEMPIE dx = vx, vy, vz, Fx, Fy, Fz
      dx.setParticleX(particles[i].px / (gamma * particles[i].m));
      dx.setParticleY(particles[i].py / (gamma * particles[i].m));
      dx.setParticleZ(particles[i].pz / (gamma * particles[i].m));
      dx.setParticlePX(particles[i].q * fieldOnParticles[0].ex + particles[i].q * dx.y * fieldOnParticles[0].bz - particles[i].q * dx.z * fieldOnParticles[0].by);
      dx.setParticlePY(particles[i].q * fieldOnParticles[0].ey + particles[i].q * dx.z * fieldOnParticles[0].bx - particles[i].q * dx.x * fieldOnParticles[0].bz);
      dx.setParticlePZ(particles[i].q * fieldOnParticles[0].ez + particles[i].q * dx.x * fieldOnParticles[0].by - particles[i].q * dx.y * fieldOnParticles[0].bx);

      // RIEMPIE dx2 per il metodo RK4 con il dx del primo passo
      dx2 = dx;

      //RIEMPIE x2 con le posizioni avanzate di particles[i] dopo il primo passo
      x2.setParticleX(particles[i].x + dx.x * (data.dt * 0.5));
      x2.setParticleY(particles[i].y + dx.y * (data.dt * 0.5));
      x2.setParticleZ(particles[i].z + dx.z * (data.dt * 0.5));
      x2.setParticlePX(particles[i].px + dx.px * (data.dt * 0.5));
      x2.setParticlePY(particles[i].py + dx.py * (data.dt * 0.5));
      x2.setParticlePZ(particles[i].pz + dx.pz * (data.dt * 0.5));

      // PASSO #2
      // AVANZA TEMPO DI deltaT/2 e calcola campi nella posizione x2
      //fieldOnParticles.clear();
      data.aumentaT(data.dt * 0.5);
      calcolaCampoAnalitico1DSuParticella(data, x2, *campoSuPunto);
      fieldOnParticles.push_back(*campoSuPunto);

      gamma = sqrt(1 + pow(x2.px, 2) + pow(x2.py, 2) + pow(x2.pz, 2));

      // AGGIORNA dx con i nuovi valori
      dx.setParticleX(x2.px / (gamma * particles[i].m));
      dx.setParticleY(x2.py / (gamma * particles[i].m));
      dx.setParticleZ(x2.pz / (gamma * particles[i].m));
      dx.setParticlePX(particles[i].q * fieldOnParticles[1].ex + particles[i].q * dx.y * fieldOnParticles[1].bz - particles[i].q * dx.z * fieldOnParticles[1].by);
      dx.setParticlePY(particles[i].q * fieldOnParticles[1].ey + particles[i].q * dx.z * fieldOnParticles[1].bx - particles[i].q * dx.x * fieldOnParticles[1].bz);
      dx.setParticlePZ(particles[i].q * fieldOnParticles[1].ez + particles[i].q * dx.x * fieldOnParticles[1].by - particles[i].q * dx.y * fieldOnParticles[1].bx);

      // AGGIUNGE al dx2 i valori calcolati al secondo passo
      dx2.setParticleX(dx2.x + 2 * dx.x);
      dx2.setParticleY(dx2.y + 2 * dx.y);
      dx2.setParticleZ(dx2.z + 2 * dx.z);
      dx2.setParticlePX(dx2.px + 2 * dx.px);
      dx2.setParticlePY(dx2.py + 2 * dx.py);
      dx2.setParticlePZ(dx2.pz + 2 * dx.pz);

      // RIEMPIE x2 con le posizione avanzate di particles[i] dopo il secondo passo
      x2.setParticleX(particles[i].x + dx.x * (data.dt * 0.5));
      x2.setParticleY(particles[i].y + dx.y * (data.dt * 0.5));
      x2.setParticleZ(particles[i].z + dx.z * (data.dt * 0.5));
      x2.setParticlePX(particles[i].px + dx.px * (data.dt * 0.5));
      x2.setParticlePY(particles[i].py + dx.py * (data.dt * 0.5));
      x2.setParticlePZ(particles[i].pz + dx.pz * (data.dt * 0.5));

      // PASSO #3
      // CALCOLA campi nella posizione x2 appena aggiornata
      //fieldOnParticles.clear();
      calcolaCampoAnalitico1DSuParticella(data, x2, *campoSuPunto);
      fieldOnParticles.push_back(*campoSuPunto);

      gamma = sqrt(1 + pow(x2.px, 2) + pow(x2.py, 2) + pow(x2.pz, 2));

      //AGGIORNA dx con i nuovi valori
      dx.setParticleX(x2.px / (gamma * particles[i].m));
      dx.setParticleY(x2.py / (gamma * particles[i].m));
      dx.setParticleZ(x2.pz / (gamma * particles[i].m));
      dx.setParticlePX(particles[i].q * fieldOnParticles[2].ex + particles[i].q * dx.y * fieldOnParticles[2].bz - particles[i].q * dx.z * fieldOnParticles[2].by);
      dx.setParticlePY(particles[i].q * fieldOnParticles[2].ey + particles[i].q * dx.z * fieldOnParticles[2].bx - particles[i].q * dx.x * fieldOnParticles[2].bz);
      dx.setParticlePZ(particles[i].q * fieldOnParticles[2].ez + particles[i].q * dx.x * fieldOnParticles[2].by - particles[i].q * dx.y * fieldOnParticles[2].bx);

      //AGGIUNGE al dx2 i valori calcolati al terzo passo
      dx2.setParticleX(dx2.x + 2 * dx.x);
      dx2.setParticleY(dx2.y + 2 * dx.y);
      dx2.setParticleZ(dx2.z + 2 * dx.z);
      dx2.setParticlePX(dx2.px + 2 * dx.px);
      dx2.setParticlePY(dx2.py + 2 * dx.py);
      dx2.setParticlePZ(dx2.pz + 2 * dx.pz);

      //RIEMPIE x2 con le posizione avanzate di particles[i]
      x2.setParticleX(particles[i].x + dx.x * data.dt);
      x2.setParticleY(particles[i].y + dx.y * data.dt);
      x2.setParticleZ(particles[i].z + dx.z * data.dt);
      x2.setParticlePX(particles[i].px + dx.px * data.dt);
      x2.setParticlePY(particles[i].py + dx.py * data.dt);
      x2.setParticlePZ(particles[i].pz + dx.pz * data.dt);

      // PASSO #4
      //AVANZA TEMPO DI deltaT/2 e calcola campi nella posizione x2
      //fieldOnParticles.clear();
      data.aumentaT(data.dt * 0.5);
      calcolaCampoAnalitico1DSuParticella(data, x2, *campoSuPunto);
      fieldOnParticles.push_back(*campoSuPunto);
      gamma = sqrt(1 + pow(x2.px, 2) + pow(x2.py, 2) + pow(x2.pz, 2));

      //AGGIORNA dx con i nuovi valori
      dx.setParticleX(x2.px / (gamma * particles[i].m));
      dx.setParticleY(x2.py / (gamma * particles[i].m));
      dx.setParticleZ(x2.pz / (gamma * particles[i].m));
      dx.setParticlePX(particles[i].q * fieldOnParticles[3].ex + particles[i].q * dx.y * fieldOnParticles[3].bz - particles[i].q * dx.z * fieldOnParticles[3].by);
      dx.setParticlePY(particles[i].q * fieldOnParticles[3].ey + particles[i].q * dx.z * fieldOnParticles[3].bx - particles[i].q * dx.x * fieldOnParticles[3].bz);
      dx.setParticlePZ(particles[i].q * fieldOnParticles[3].ez + particles[i].q * dx.x * fieldOnParticles[3].by - particles[i].q * dx.y * fieldOnParticles[3].bx);

      //AGGIUNGE al dx2 i valori calcolati al terzo passo
      dx2.setParticleX(dx2.x + dx.x);
      dx2.setParticleY(dx2.y + dx.y);
      dx2.setParticleZ(dx2.z + dx.z);
      dx2.setParticlePX(dx2.px + dx.px);
      dx2.setParticlePY(dx2.py + dx.py);
      dx2.setParticlePZ(dx2.pz + dx.pz);

      particles[i].setParticleX(particles[i].x + (data.dt / 6.) * dx2.x);
      particles[i].setParticleY(particles[i].y + (data.dt / 6.) * dx2.y);
      particles[i].setParticleZ(particles[i].z + (data.dt / 6.) * dx2.z);
      particles[i].setParticlePX(particles[i].px + (data.dt / 6.) * dx2.px);
      particles[i].setParticlePY(particles[i].py + (data.dt / 6.) * dx2.py);
      particles[i].setParticlePZ(particles[i].pz + (data.dt / 6.) * dx2.pz);

      data.aumentaT(-data.dt); // faccio tornare il tempo indietro cosi` gli altri elettroni di questo step non hanno il tempo gia` avanzato
          // il tempo avanza definitivamente in ciascun ciclo 'for' esterno (sugli nsteps)

      // INTEGRALI PRIMI
      potenziale.clear();
      zeta = particles[i].x - SPEED_OF_LIGHT * data.t; // non usato nella mia struttura delle equazioni di campo
      gamma = sqrt(1 + pow(particles[i].px, 2) + pow(particles[i].py, 2) + pow(particles[i].pz, 2));
      calcolaPotenzialeAnalitico1DSuParticella(data, particles[i], *campoSuPunto);
      potenziale.push_back(*campoSuPunto);

      I1 = (particles[i].py / (particles[i].m * SPEED_OF_LIGHT)) + potenziale[0].ey;
      I2 = particles[i].pz / (particles[i].m * SPEED_OF_LIGHT);
      I3 = gamma - particles[i].px / (particles[i].m * SPEED_OF_LIGHT);

      potenziale.clear();

      outputstream
          << setw(12) << i + 1 << "\t"
          << setprecision(14) << setw(16) << particles[i].x << "\t"
          << setw(16) << particles[i].y << "\t"
          << setw(16) << particles[i].z << "\t"
          // << setw(16) << particles[i].getParticleCell() << "\t"
          << setw(16) << particles[i].px << "\t"
          << setw(16) << particles[i].py << "\t"
          << setw(16) << particles[i].pz << "\t"
          << setw(16) << zeta << "\t"
          << setw(16) << I1 << "\t"
          << setw(16) << I2 << "\t"
          << setw(16) << I3 << endl;
    }
    data.aumentaT(data.dt);
  }
}

void evolveRK4_withgrid_onthefly(Data data, vector<Particle>& particles, vector<Field>& fieldOnParticles, ofstream& outputstream)
{
  Particle x2, dx, dx2;
  double gamma;
  double zeta;
  double I1, I2, I3;

  Particle coordinatePrimoPuntoGriglia;
  Particle coordinateSecondoPuntoGriglia;
  Particle coordinateTerzoPuntoGriglia;
  Particle coordinateQuartoPuntoGriglia;
  Particle coordinateQuintoPuntoGriglia;
  Particle coordinateSestoPuntoGriglia;
  Particle coordinateSettimoPuntoGriglia;
  Particle coordinateOttavoPuntoGriglia;

  double gridXX, gridYY, gridZZ;
  int gridX, gridY, gridZ;

  Field FOP;

  Field* FOG;
  FOG = new Field[1];
  Field* campoSuPunto;
  campoSuPunto = new Field[1];

  vector<Field> fieldOnGrid;
  vector<Field> potenziale;

  for (int j = 0; j < data.nSteps; j++) {
    data.aumentaT(data.dt);
    for (int i = 0; i < data.n_electrons; i++) {

      fieldOnParticles.clear(); //SVUOTA il vettore dei campi sulle particelle ad ogni ciclo

      gridXX = particles[i].x / data.deltaX;
      if (gridXX < 0)
        gridX = ((int)gridXX) - 1; //il punto di griglia lungo l'asse x a sinistra della particella
      else
        gridX = (int)gridXX;

      gridYY = particles[i].y / data.deltaY;
      if (gridYY < 0)
        gridY = ((int)gridYY) - 1; //il punto di griglia lungo l'asse y a sinistra della particella
      else
        gridY = (int)gridYY;

      gridZZ = particles[i].z / data.deltaZ;
      if (gridZZ < 0)
        gridZ = ((int)gridZZ) - 1; //il punto di griglia lungo l'asse z a sinistra della particella
      else
        gridZ = (int)gridZZ;

      coordinatePrimoPuntoGriglia.setParticleX(gridX * data.deltaX);
      coordinatePrimoPuntoGriglia.setParticleY(gridY * data.deltaY);
      coordinatePrimoPuntoGriglia.setParticleZ(gridZ * data.deltaZ);
      coordinateSecondoPuntoGriglia.setParticleX(gridX * data.deltaX + data.deltaX);
      coordinateSecondoPuntoGriglia.setParticleY(gridY * data.deltaY);
      coordinateSecondoPuntoGriglia.setParticleZ(gridZ * data.deltaZ);
      coordinateTerzoPuntoGriglia.setParticleX(gridX * data.deltaX);
      coordinateTerzoPuntoGriglia.setParticleY(gridY * data.deltaY + data.deltaY);
      coordinateTerzoPuntoGriglia.setParticleZ(gridZ * data.deltaZ);
      coordinateQuartoPuntoGriglia.setParticleX(gridX * data.deltaX + data.deltaX);
      coordinateQuartoPuntoGriglia.setParticleY(gridY * data.deltaY + data.deltaY);
      coordinateQuartoPuntoGriglia.setParticleZ(gridZ * data.deltaZ);
      coordinateQuintoPuntoGriglia.setParticleX(gridX * data.deltaX);
      coordinateQuintoPuntoGriglia.setParticleY(gridY * data.deltaY);
      coordinateQuintoPuntoGriglia.setParticleZ(gridZ * data.deltaZ + data.deltaZ);
      coordinateSestoPuntoGriglia.setParticleX(gridX * data.deltaX + data.deltaX);
      coordinateSestoPuntoGriglia.setParticleY(gridY * data.deltaY);
      coordinateSestoPuntoGriglia.setParticleZ(gridZ * data.deltaZ + data.deltaZ);
      coordinateSettimoPuntoGriglia.setParticleX(gridX * data.deltaX);
      coordinateSettimoPuntoGriglia.setParticleY(gridY * data.deltaY + data.deltaY);
      coordinateSettimoPuntoGriglia.setParticleZ(gridZ * data.deltaZ + data.deltaZ);
      coordinateOttavoPuntoGriglia.setParticleX(gridX * data.deltaX + data.deltaX);
      coordinateOttavoPuntoGriglia.setParticleY(gridY * data.deltaY + data.deltaY);
      coordinateOttavoPuntoGriglia.setParticleZ(gridZ * data.deltaZ + data.deltaZ);

      calcolaCampoAnalitico1DSuParticella(data, coordinatePrimoPuntoGriglia, *FOG);
      fieldOnGrid.push_back(*FOG);
      calcolaCampoAnalitico1DSuParticella(data, coordinateSecondoPuntoGriglia, *FOG);
      fieldOnGrid.push_back(*FOG);
      calcolaCampoAnalitico1DSuParticella(data, coordinateTerzoPuntoGriglia, *FOG);
      fieldOnGrid.push_back(*FOG);
      calcolaCampoAnalitico1DSuParticella(data, coordinateQuartoPuntoGriglia, *FOG);
      fieldOnGrid.push_back(*FOG);
      calcolaCampoAnalitico1DSuParticella(data, coordinateQuintoPuntoGriglia, *FOG);
      fieldOnGrid.push_back(*FOG);
      calcolaCampoAnalitico1DSuParticella(data, coordinateSestoPuntoGriglia, *FOG);
      fieldOnGrid.push_back(*FOG);
      calcolaCampoAnalitico1DSuParticella(data, coordinateSettimoPuntoGriglia, *FOG);
      fieldOnGrid.push_back(*FOG);
      calcolaCampoAnalitico1DSuParticella(data, coordinateOttavoPuntoGriglia, *FOG);
      fieldOnGrid.push_back(*FOG);

      FOP = interpolation2D_linear(data, particles[i], fieldOnGrid);
      fieldOnParticles.push_back(FOP);

      fieldOnGrid.clear();

      gamma = sqrt(1 + pow(particles[i].px, 2) + pow(particles[i].py, 2) + pow(particles[i].pz, 2));

      // PASSO #1
      // RIEMPIE dx = vx, vy, vz, Fx, Fy, Fz
      dx.setParticleX(particles[i].px / (gamma * particles[i].m));
      dx.setParticleY(particles[i].py / (gamma * particles[i].m));
      dx.setParticleZ(particles[i].pz / (gamma * particles[i].m));
      dx.setParticlePX(particles[i].q * fieldOnParticles[0].ex + particles[i].q * dx.y * fieldOnParticles[0].bz - particles[i].q * dx.z * fieldOnParticles[0].by);
      dx.setParticlePY(particles[i].q * fieldOnParticles[0].ey + particles[i].q * dx.z * fieldOnParticles[0].bx - particles[i].q * dx.x * fieldOnParticles[0].bz);
      dx.setParticlePZ(particles[i].q * fieldOnParticles[0].ez + particles[i].q * dx.x * fieldOnParticles[0].by - particles[i].q * dx.y * fieldOnParticles[0].bx);

      // RIEMPIE dx2 per il metodo RK4 con il dx del primo passo
      dx2 = dx;

      // RIEMPIE x2 con le posizioni avanzate di particles[i] dopo il primo passo
      x2.setParticleX(particles[i].x + dx.x * (data.dt * 0.5));
      x2.setParticleY(particles[i].y + dx.y * (data.dt * 0.5));
      x2.setParticleZ(particles[i].z + dx.z * (data.dt * 0.5));
      x2.setParticlePX(particles[i].px + dx.px * (data.dt * 0.5));
      x2.setParticlePY(particles[i].py + dx.py * (data.dt * 0.5));
      x2.setParticlePZ(particles[i].pz + dx.pz * (data.dt * 0.5));

      // PASSO #2
      // AVANZA TEMPO DI deltaT/2 e calcola campi nella posizione x2
      gridXX = x2.x / data.deltaX;
      if (gridXX < 0)
        gridX = ((int)gridXX) - 1; //il punto di griglia lungo l'asse x a sinistra della particella
      else
        gridX = (int)gridXX;

      gridYY = x2.y / data.deltaY;
      if (gridYY < 0)
        gridY = ((int)gridYY) - 1; //il punto di griglia lungo l'asse y a sinistra della particella
      else
        gridY = (int)gridYY;

      gridZZ = x2.z / data.deltaZ;
      if (gridZZ < 0)
        gridZ = ((int)gridZZ) - 1; //il punto di griglia lungo l'asse z a sinistra della particella
      else
        gridZ = (int)gridZZ;

      coordinatePrimoPuntoGriglia.setParticleX(gridX * data.deltaX);
      coordinatePrimoPuntoGriglia.setParticleY(gridY * data.deltaY);
      coordinatePrimoPuntoGriglia.setParticleZ(gridZ * data.deltaZ);
      coordinateSecondoPuntoGriglia.setParticleX(gridX * data.deltaX + data.deltaX);
      coordinateSecondoPuntoGriglia.setParticleY(gridY * data.deltaY);
      coordinateSecondoPuntoGriglia.setParticleZ(gridZ * data.deltaZ);
      coordinateTerzoPuntoGriglia.setParticleX(gridX * data.deltaX);
      coordinateTerzoPuntoGriglia.setParticleY(gridY * data.deltaY + data.deltaY);
      coordinateTerzoPuntoGriglia.setParticleZ(gridZ * data.deltaZ);
      coordinateQuartoPuntoGriglia.setParticleX(gridX * data.deltaX + data.deltaX);
      coordinateQuartoPuntoGriglia.setParticleY(gridY * data.deltaY + data.deltaY);
      coordinateQuartoPuntoGriglia.setParticleZ(gridZ * data.deltaZ);
      coordinateQuintoPuntoGriglia.setParticleX(gridX * data.deltaX);
      coordinateQuintoPuntoGriglia.setParticleY(gridY * data.deltaY);
      coordinateQuintoPuntoGriglia.setParticleZ(gridZ * data.deltaZ + data.deltaZ);
      coordinateSestoPuntoGriglia.setParticleX(gridX * data.deltaX + data.deltaX);
      coordinateSestoPuntoGriglia.setParticleY(gridY * data.deltaY);
      coordinateSestoPuntoGriglia.setParticleZ(gridZ * data.deltaZ + data.deltaZ);
      coordinateSettimoPuntoGriglia.setParticleX(gridX * data.deltaX);
      coordinateSettimoPuntoGriglia.setParticleY(gridY * data.deltaY + data.deltaY);
      coordinateSettimoPuntoGriglia.setParticleZ(gridZ * data.deltaZ + data.deltaZ);
      coordinateOttavoPuntoGriglia.setParticleX(gridX * data.deltaX + data.deltaX);
      coordinateOttavoPuntoGriglia.setParticleY(gridY * data.deltaY + data.deltaY);
      coordinateOttavoPuntoGriglia.setParticleZ(gridZ * data.deltaZ + data.deltaZ);

      calcolaCampoAnalitico1DSuParticella(data, coordinatePrimoPuntoGriglia, *FOG);
      fieldOnGrid.push_back(*FOG);
      calcolaCampoAnalitico1DSuParticella(data, coordinateSecondoPuntoGriglia, *FOG);
      fieldOnGrid.push_back(*FOG);
      calcolaCampoAnalitico1DSuParticella(data, coordinateTerzoPuntoGriglia, *FOG);
      fieldOnGrid.push_back(*FOG);
      calcolaCampoAnalitico1DSuParticella(data, coordinateQuartoPuntoGriglia, *FOG);
      fieldOnGrid.push_back(*FOG);
      calcolaCampoAnalitico1DSuParticella(data, coordinateQuintoPuntoGriglia, *FOG);
      fieldOnGrid.push_back(*FOG);
      calcolaCampoAnalitico1DSuParticella(data, coordinateSestoPuntoGriglia, *FOG);
      fieldOnGrid.push_back(*FOG);
      calcolaCampoAnalitico1DSuParticella(data, coordinateSettimoPuntoGriglia, *FOG);
      fieldOnGrid.push_back(*FOG);
      calcolaCampoAnalitico1DSuParticella(data, coordinateOttavoPuntoGriglia, *FOG);
      fieldOnGrid.push_back(*FOG);

      data.aumentaT(data.dt * 0.5);
      FOP = interpolation2D_linear(data, x2, fieldOnGrid);
      fieldOnParticles.push_back(FOP);

      fieldOnGrid.clear();

      gamma = sqrt(1 + pow(x2.px, 2) + pow(x2.py, 2) + pow(x2.pz, 2));

      // AGGIORNA dx con i nuovi valori
      dx.setParticleX(x2.px / (gamma * particles[i].m));
      dx.setParticleY(x2.py / (gamma * particles[i].m));
      dx.setParticleZ(x2.pz / (gamma * particles[i].m));
      dx.setParticlePX(particles[i].q * fieldOnParticles[1].ex + particles[i].q * dx.y * fieldOnParticles[1].bz - particles[i].q * dx.z * fieldOnParticles[1].by);
      dx.setParticlePY(particles[i].q * fieldOnParticles[1].ey + particles[i].q * dx.z * fieldOnParticles[1].bx - particles[i].q * dx.x * fieldOnParticles[1].bz);
      dx.setParticlePZ(particles[i].q * fieldOnParticles[1].ez + particles[i].q * dx.x * fieldOnParticles[1].by - particles[i].q * dx.y * fieldOnParticles[1].bx);

      // AGGIUNGE al dx2 i valori calcolati al secondo passo
      dx2.setParticleX(dx2.x + 2 * dx.x);
      dx2.setParticleY(dx2.y + 2 * dx.y);
      dx2.setParticleZ(dx2.z + 2 * dx.z);
      dx2.setParticlePX(dx2.px + 2 * dx.px);
      dx2.setParticlePY(dx2.py + 2 * dx.py);
      dx2.setParticlePZ(dx2.pz + 2 * dx.pz);

      // RIEMPIE x2 con le posizione avanzate di particles[i] dopo il secondo passo
      x2.setParticleX(particles[i].x + dx.x * (data.dt * 0.5));
      x2.setParticleY(particles[i].y + dx.y * (data.dt * 0.5));
      x2.setParticleZ(particles[i].z + dx.z * (data.dt * 0.5));
      x2.setParticlePX(particles[i].px + dx.px * (data.dt * 0.5));
      x2.setParticlePY(particles[i].py + dx.py * (data.dt * 0.5));
      x2.setParticlePZ(particles[i].pz + dx.pz * (data.dt * 0.5));

      // PASSO #3
      // CALCOLA campi nella posizione x2 appena aggiornata
      gridXX = x2.x / data.deltaX;
      if (gridXX < 0)
        gridX = ((int)gridXX) - 1; //il punto di griglia lungo l'asse x a sinistra della particella
      else
        gridX = (int)gridXX;

      gridYY = x2.y / data.deltaY;
      if (gridYY < 0)
        gridY = ((int)gridYY) - 1; //il punto di griglia lungo l'asse y a sinistra della particella
      else
        gridY = (int)gridYY;

      gridZZ = x2.z / data.deltaZ;
      if (gridZZ < 0)
        gridZ = ((int)gridZZ) - 1; //il punto di griglia lungo l'asse z a sinistra della particella
      else
        gridZ = (int)gridZZ;

      coordinatePrimoPuntoGriglia.setParticleX(gridX * data.deltaX);
      coordinatePrimoPuntoGriglia.setParticleY(gridY * data.deltaY);
      coordinatePrimoPuntoGriglia.setParticleZ(gridZ * data.deltaZ);
      coordinateSecondoPuntoGriglia.setParticleX(gridX * data.deltaX + data.deltaX);
      coordinateSecondoPuntoGriglia.setParticleY(gridY * data.deltaY);
      coordinateSecondoPuntoGriglia.setParticleZ(gridZ * data.deltaZ);
      coordinateTerzoPuntoGriglia.setParticleX(gridX * data.deltaX);
      coordinateTerzoPuntoGriglia.setParticleY(gridY * data.deltaY + data.deltaY);
      coordinateTerzoPuntoGriglia.setParticleZ(gridZ * data.deltaZ);
      coordinateQuartoPuntoGriglia.setParticleX(gridX * data.deltaX + data.deltaX);
      coordinateQuartoPuntoGriglia.setParticleY(gridY * data.deltaY + data.deltaY);
      coordinateQuartoPuntoGriglia.setParticleZ(gridZ * data.deltaZ);
      coordinateQuintoPuntoGriglia.setParticleX(gridX * data.deltaX);
      coordinateQuintoPuntoGriglia.setParticleY(gridY * data.deltaY);
      coordinateQuintoPuntoGriglia.setParticleZ(gridZ * data.deltaZ + data.deltaZ);
      coordinateSestoPuntoGriglia.setParticleX(gridX * data.deltaX + data.deltaX);
      coordinateSestoPuntoGriglia.setParticleY(gridY * data.deltaY);
      coordinateSestoPuntoGriglia.setParticleZ(gridZ * data.deltaZ + data.deltaZ);
      coordinateSettimoPuntoGriglia.setParticleX(gridX * data.deltaX);
      coordinateSettimoPuntoGriglia.setParticleY(gridY * data.deltaY + data.deltaY);
      coordinateSettimoPuntoGriglia.setParticleZ(gridZ * data.deltaZ + data.deltaZ);
      coordinateOttavoPuntoGriglia.setParticleX(gridX * data.deltaX + data.deltaX);
      coordinateOttavoPuntoGriglia.setParticleY(gridY * data.deltaY + data.deltaY);
      coordinateOttavoPuntoGriglia.setParticleZ(gridZ * data.deltaZ + data.deltaZ);

      calcolaCampoAnalitico1DSuParticella(data, coordinatePrimoPuntoGriglia, *FOG);
      fieldOnGrid.push_back(*FOG);
      calcolaCampoAnalitico1DSuParticella(data, coordinateSecondoPuntoGriglia, *FOG);
      fieldOnGrid.push_back(*FOG);
      calcolaCampoAnalitico1DSuParticella(data, coordinateTerzoPuntoGriglia, *FOG);
      fieldOnGrid.push_back(*FOG);
      calcolaCampoAnalitico1DSuParticella(data, coordinateQuartoPuntoGriglia, *FOG);
      fieldOnGrid.push_back(*FOG);
      calcolaCampoAnalitico1DSuParticella(data, coordinateQuintoPuntoGriglia, *FOG);
      fieldOnGrid.push_back(*FOG);
      calcolaCampoAnalitico1DSuParticella(data, coordinateSestoPuntoGriglia, *FOG);
      fieldOnGrid.push_back(*FOG);
      calcolaCampoAnalitico1DSuParticella(data, coordinateSettimoPuntoGriglia, *FOG);
      fieldOnGrid.push_back(*FOG);
      calcolaCampoAnalitico1DSuParticella(data, coordinateOttavoPuntoGriglia, *FOG);
      fieldOnGrid.push_back(*FOG);

      FOP = interpolation2D_linear(data, x2, fieldOnGrid);
      fieldOnParticles.push_back(FOP);

      fieldOnGrid.clear();

      gamma = sqrt(1 + pow(x2.px, 2) + pow(x2.py, 2) + pow(x2.pz, 2));

      // AGGIORNA dx con i nuovi valori
      dx.setParticleX(x2.px / (gamma * particles[i].m));
      dx.setParticleY(x2.py / (gamma * particles[i].m));
      dx.setParticleZ(x2.pz / (gamma * particles[i].m));
      dx.setParticlePX(particles[i].q * fieldOnParticles[2].ex + particles[i].q * dx.y * fieldOnParticles[2].bz - particles[i].q * dx.z * fieldOnParticles[2].by);
      dx.setParticlePY(particles[i].q * fieldOnParticles[2].ey + particles[i].q * dx.z * fieldOnParticles[2].bx - particles[i].q * dx.x * fieldOnParticles[2].bz);
      dx.setParticlePZ(particles[i].q * fieldOnParticles[2].ez + particles[i].q * dx.x * fieldOnParticles[2].by - particles[i].q * dx.y * fieldOnParticles[2].bx);

      // AGGIUNGE al dx2 i valori calcolati al terzo passo
      dx2.setParticleX(dx2.x + 2 * dx.x);
      dx2.setParticleY(dx2.y + 2 * dx.y);
      dx2.setParticleZ(dx2.z + 2 * dx.z);
      dx2.setParticlePX(dx2.px + 2 * dx.px);
      dx2.setParticlePY(dx2.py + 2 * dx.py);
      dx2.setParticlePZ(dx2.pz + 2 * dx.pz);

      // RIEMPIE x2 con le posizione avanzate di particles[i]
      x2.setParticleX(particles[i].x + dx.x * data.dt);
      x2.setParticleY(particles[i].y + dx.y * data.dt);
      x2.setParticleZ(particles[i].z + dx.z * data.dt);
      x2.setParticlePX(particles[i].px + dx.px * data.dt);
      x2.setParticlePY(particles[i].py + dx.py * data.dt);
      x2.setParticlePZ(particles[i].pz + dx.pz * data.dt);

      // PASSO #4
      // AVANZA TEMPO DI deltaT/2 e calcola campi nella posizione x2
      gridXX = x2.x / data.deltaX;
      if (gridXX < 0)
        gridX = ((int)gridXX) - 1; //il punto di griglia lungo l'asse x a sinistra della particella
      else
        gridX = (int)gridXX;

      gridYY = x2.y / data.deltaY;
      if (gridYY < 0)
        gridY = ((int)gridYY) - 1; //il punto di griglia lungo l'asse y a sinistra della particella
      else
        gridY = (int)gridYY;

      gridZZ = x2.z / data.deltaZ;
      if (gridZZ < 0)
        gridZ = ((int)gridZZ) - 1; //il punto di griglia lungo l'asse z a sinistra della particella
      else
        gridZ = (int)gridZZ;

      coordinatePrimoPuntoGriglia.setParticleX(gridX * data.deltaX);
      coordinatePrimoPuntoGriglia.setParticleY(gridY * data.deltaY);
      coordinatePrimoPuntoGriglia.setParticleZ(gridZ * data.deltaZ);
      coordinateSecondoPuntoGriglia.setParticleX(gridX * data.deltaX + data.deltaX);
      coordinateSecondoPuntoGriglia.setParticleY(gridY * data.deltaY);
      coordinateSecondoPuntoGriglia.setParticleZ(gridZ * data.deltaZ);
      coordinateTerzoPuntoGriglia.setParticleX(gridX * data.deltaX);
      coordinateTerzoPuntoGriglia.setParticleY(gridY * data.deltaY + data.deltaY);
      coordinateTerzoPuntoGriglia.setParticleZ(gridZ * data.deltaZ);
      coordinateQuartoPuntoGriglia.setParticleX(gridX * data.deltaX + data.deltaX);
      coordinateQuartoPuntoGriglia.setParticleY(gridY * data.deltaY + data.deltaY);
      coordinateQuartoPuntoGriglia.setParticleZ(gridZ * data.deltaZ);
      coordinateQuintoPuntoGriglia.setParticleX(gridX * data.deltaX);
      coordinateQuintoPuntoGriglia.setParticleY(gridY * data.deltaY);
      coordinateQuintoPuntoGriglia.setParticleZ(gridZ * data.deltaZ + data.deltaZ);
      coordinateSestoPuntoGriglia.setParticleX(gridX * data.deltaX + data.deltaX);
      coordinateSestoPuntoGriglia.setParticleY(gridY * data.deltaY);
      coordinateSestoPuntoGriglia.setParticleZ(gridZ * data.deltaZ + data.deltaZ);
      coordinateSettimoPuntoGriglia.setParticleX(gridX * data.deltaX);
      coordinateSettimoPuntoGriglia.setParticleY(gridY * data.deltaY + data.deltaY);
      coordinateSettimoPuntoGriglia.setParticleZ(gridZ * data.deltaZ + data.deltaZ);
      coordinateOttavoPuntoGriglia.setParticleX(gridX * data.deltaX + data.deltaX);
      coordinateOttavoPuntoGriglia.setParticleY(gridY * data.deltaY + data.deltaY);
      coordinateOttavoPuntoGriglia.setParticleZ(gridZ * data.deltaZ + data.deltaZ);

      calcolaCampoAnalitico1DSuParticella(data, coordinatePrimoPuntoGriglia, *FOG);
      fieldOnGrid.push_back(*FOG);
      calcolaCampoAnalitico1DSuParticella(data, coordinateSecondoPuntoGriglia, *FOG);
      fieldOnGrid.push_back(*FOG);
      calcolaCampoAnalitico1DSuParticella(data, coordinateTerzoPuntoGriglia, *FOG);
      fieldOnGrid.push_back(*FOG);
      calcolaCampoAnalitico1DSuParticella(data, coordinateQuartoPuntoGriglia, *FOG);
      fieldOnGrid.push_back(*FOG);
      calcolaCampoAnalitico1DSuParticella(data, coordinateQuintoPuntoGriglia, *FOG);
      fieldOnGrid.push_back(*FOG);
      calcolaCampoAnalitico1DSuParticella(data, coordinateSestoPuntoGriglia, *FOG);
      fieldOnGrid.push_back(*FOG);
      calcolaCampoAnalitico1DSuParticella(data, coordinateSettimoPuntoGriglia, *FOG);
      fieldOnGrid.push_back(*FOG);
      calcolaCampoAnalitico1DSuParticella(data, coordinateOttavoPuntoGriglia, *FOG);
      fieldOnGrid.push_back(*FOG);

      data.aumentaT(data.dt * 0.5);
      FOP = interpolation2D_linear(data, x2, fieldOnGrid);
      fieldOnParticles.push_back(FOP);

      fieldOnGrid.clear();

      gamma = sqrt(1 + pow(x2.px, 2) + pow(x2.py, 2) + pow(x2.pz, 2));

      //AGGIORNA dx con i nuovi valori
      dx.setParticleX(x2.px / (gamma * particles[i].m));
      dx.setParticleY(x2.py / (gamma * particles[i].m));
      dx.setParticleZ(x2.pz / (gamma * particles[i].m));
      dx.setParticlePX(particles[i].q * fieldOnParticles[3].ex + particles[i].q * dx.y * fieldOnParticles[3].bz - particles[i].q * dx.z * fieldOnParticles[3].by);
      dx.setParticlePY(particles[i].q * fieldOnParticles[3].ey + particles[i].q * dx.z * fieldOnParticles[3].bx - particles[i].q * dx.x * fieldOnParticles[3].bz);
      dx.setParticlePZ(particles[i].q * fieldOnParticles[3].ez + particles[i].q * dx.x * fieldOnParticles[3].by - particles[i].q * dx.y * fieldOnParticles[3].bx);

      //AGGIUNGE al dx2 i valori calcolati al terzo passo
      dx2.setParticleX(dx2.x + dx.x);
      dx2.setParticleY(dx2.y + dx.y);
      dx2.setParticleZ(dx2.z + dx.z);
      dx2.setParticlePX(dx2.px + dx.px);
      dx2.setParticlePY(dx2.py + dx.py);
      dx2.setParticlePZ(dx2.pz + dx.pz);

      particles[i].setParticleX(particles[i].x + (data.dt / 6.) * dx2.x);
      particles[i].setParticleY(particles[i].y + (data.dt / 6.) * dx2.y);
      particles[i].setParticleZ(particles[i].z + (data.dt / 6.) * dx2.z);
      particles[i].setParticlePX(particles[i].px + (data.dt / 6.) * dx2.px);
      particles[i].setParticlePY(particles[i].py + (data.dt / 6.) * dx2.py);
      particles[i].setParticlePZ(particles[i].pz + (data.dt / 6.) * dx2.pz);

      data.aumentaT(-data.dt); // faccio tornare il tempo indietro cosi` gli altri elettroni di questo step non hanno il tempo gia` avanzato
          // il tempo avanza definitivamente in ciascun ciclo for esterno (sugli nsteps)

      //INTEGRALI PRIMI
      potenziale.clear();
      zeta = particles[i].x - SPEED_OF_LIGHT * data.t; // non usato nella mia struttura delle equazioni di campo
      gamma = sqrt(1 + pow(particles[i].px, 2) + pow(particles[i].py, 2) + pow(particles[i].pz, 2));

      gridXX = particles[i].x / data.deltaX;
      if (gridXX < 0)
        gridX = ((int)gridXX) - 1; //il punto di griglia lungo l'asse x a sinistra della particella
      else
        gridX = (int)gridXX;

      gridYY = particles[i].y / data.deltaY;
      if (gridYY < 0)
        gridY = ((int)gridYY) - 1; //il punto di griglia lungo l'asse y a sinistra della particella
      else
        gridY = (int)gridYY;

      gridZZ = particles[i].z / data.deltaZ;
      if (gridZZ < 0)
        gridZ = ((int)gridZZ) - 1; //il punto di griglia lungo l'asse z a sinistra della particella
      else
        gridZ = (int)gridZZ;

      coordinatePrimoPuntoGriglia.setParticleX(gridX * data.deltaX);
      coordinatePrimoPuntoGriglia.setParticleY(gridY * data.deltaY);
      coordinatePrimoPuntoGriglia.setParticleZ(gridZ * data.deltaZ);
      coordinateSecondoPuntoGriglia.setParticleX(gridX * data.deltaX + data.deltaX);
      coordinateSecondoPuntoGriglia.setParticleY(gridY * data.deltaY);
      coordinateSecondoPuntoGriglia.setParticleZ(gridZ * data.deltaZ);
      coordinateTerzoPuntoGriglia.setParticleX(gridX * data.deltaX);
      coordinateTerzoPuntoGriglia.setParticleY(gridY * data.deltaY + data.deltaY);
      coordinateTerzoPuntoGriglia.setParticleZ(gridZ * data.deltaZ);
      coordinateQuartoPuntoGriglia.setParticleX(gridX * data.deltaX + data.deltaX);
      coordinateQuartoPuntoGriglia.setParticleY(gridY * data.deltaY + data.deltaY);
      coordinateQuartoPuntoGriglia.setParticleZ(gridZ * data.deltaZ);
      coordinateQuintoPuntoGriglia.setParticleX(gridX * data.deltaX);
      coordinateQuintoPuntoGriglia.setParticleY(gridY * data.deltaY);
      coordinateQuintoPuntoGriglia.setParticleZ(gridZ * data.deltaZ + data.deltaZ);
      coordinateSestoPuntoGriglia.setParticleX(gridX * data.deltaX + data.deltaX);
      coordinateSestoPuntoGriglia.setParticleY(gridY * data.deltaY);
      coordinateSestoPuntoGriglia.setParticleZ(gridZ * data.deltaZ + data.deltaZ);
      coordinateSettimoPuntoGriglia.setParticleX(gridX * data.deltaX);
      coordinateSettimoPuntoGriglia.setParticleY(gridY * data.deltaY + data.deltaY);
      coordinateSettimoPuntoGriglia.setParticleZ(gridZ * data.deltaZ + data.deltaZ);
      coordinateOttavoPuntoGriglia.setParticleX(gridX * data.deltaX + data.deltaX);
      coordinateOttavoPuntoGriglia.setParticleY(gridY * data.deltaY + data.deltaY);
      coordinateOttavoPuntoGriglia.setParticleZ(gridZ * data.deltaZ + data.deltaZ);

      calcolaPotenzialeAnalitico1DSuParticella(data, coordinatePrimoPuntoGriglia, *FOG);
      fieldOnGrid.push_back(*FOG);
      calcolaPotenzialeAnalitico1DSuParticella(data, coordinateSecondoPuntoGriglia, *FOG);
      fieldOnGrid.push_back(*FOG);
      calcolaPotenzialeAnalitico1DSuParticella(data, coordinateTerzoPuntoGriglia, *FOG);
      fieldOnGrid.push_back(*FOG);
      calcolaPotenzialeAnalitico1DSuParticella(data, coordinateQuartoPuntoGriglia, *FOG);
      fieldOnGrid.push_back(*FOG);
      calcolaPotenzialeAnalitico1DSuParticella(data, coordinateQuintoPuntoGriglia, *FOG);
      fieldOnGrid.push_back(*FOG);
      calcolaPotenzialeAnalitico1DSuParticella(data, coordinateSestoPuntoGriglia, *FOG);
      fieldOnGrid.push_back(*FOG);
      calcolaPotenzialeAnalitico1DSuParticella(data, coordinateSettimoPuntoGriglia, *FOG);
      fieldOnGrid.push_back(*FOG);
      calcolaPotenzialeAnalitico1DSuParticella(data, coordinateOttavoPuntoGriglia, *FOG);
      fieldOnGrid.push_back(*FOG);

      FOP = interpolation2D_linear(data, particles[i], fieldOnGrid);
      potenziale.push_back(FOP);

      fieldOnGrid.clear();

      I1 = (particles[i].py / (particles[i].m * SPEED_OF_LIGHT)) + potenziale[0].ey;
      I2 = particles[i].pz / (particles[i].m * SPEED_OF_LIGHT);
      I3 = gamma - particles[i].px / (particles[i].m * SPEED_OF_LIGHT);

      potenziale.clear();

      outputstream
          << setw(12) << i + 1 << "\t"
          << setprecision(14) << setw(16) << particles[i].x << "\t"
          << setw(16) << particles[i].y << "\t"
          << setw(16) << particles[i].z << "\t"
          //        << setw(16) << particles[i].getParticleCell() << "\t"
          << setw(16) << particles[i].px << "\t"
          << setw(16) << particles[i].py << "\t"
          << setw(16) << particles[i].pz << "\t"
          << setw(16) << zeta << "\t"
          << setw(16) << I1 << "\t"
          << setw(16) << I2 << "\t"
          << setw(16) << I3 << endl;
    }
  }
}

void evolveRK4_withgrid(Data data, vector<Particle>& particles, vector<Field>& fieldOnParticles, vector<Field>& campoSuPuntiGriglia, ofstream& outputstream)
{

  Particle x2, dx, dx2;
  double gamma;
  double zeta;
  double I1, I2, I3;

  double gridXX, gridYY, gridZZ;
  int gridX, gridY, gridZ;
  Field FOP;
  vector<Field> fieldOnGrid;
  Field* FOG; //DA SISTEMARE, c'e` solo per avere un campo fittizio su una griglia on-the-fly per implementare l'algoritmo
  FOG = new Field[1];

  Field* campoSuPunto;
  campoSuPunto = new Field[1];
  vector<Field> potenziale;

  for (int j = 0; j < data.nSteps; j++) {
    data.aumentaT(data.dt);
    for (int i = 0; i < data.n_electrons; i++) {

      fieldOnParticles.clear(); //SVUOTA il vettore dei campi sulle particelle ad ogni ciclo

      gridXX = particles[i].x / data.deltaX;
      if (gridXX < 0)
        gridX = ((int)gridXX) - 1; //il punto di griglia lungo l'asse x a sinistra della particella
      else
        gridX = (int)gridXX;

      gridYY = particles[i].y / data.deltaY;
      if (gridYY < 0)
        gridY = ((int)gridYY) - 1; //il punto di griglia lungo l'asse y a sinistra della particella
      else
        gridY = (int)gridYY;

      gridZZ = particles[i].z / data.deltaZ;
      if (gridZZ < 0)
        gridZ = ((int)gridZZ) - 1; //il punto di griglia lungo l'asse z a sinistra della particella
      else
        gridZ = (int)gridZZ;

      if (gridX < 0 || gridY < 0 || gridZ < 0) {
        cout << "La particella e` uscita dai boundaries - simulazione interrotta" << endl;
        return;
      }

      campoSuPuntiGriglia.clear();
      riempiPuntiGrigliaConCampoAnalitico1D(data, campoSuPuntiGriglia);

      fieldOnGrid.push_back(campoSuPuntiGriglia.at(gridY * data.nGridPointsX + gridX)); // primo punto
      fieldOnGrid.push_back(campoSuPuntiGriglia.at(gridY * data.nGridPointsX + (gridX + 1))); // secondo punto
      fieldOnGrid.push_back(campoSuPuntiGriglia.at((gridY + 1) * data.nGridPointsX + gridX)); // terzo punto
      fieldOnGrid.push_back(campoSuPuntiGriglia.at((gridY + 1) * data.nGridPointsX + (gridX + 1))); // quarto punto

      FOP = interpolation2D_linear(data, particles[i], fieldOnGrid);
      fieldOnParticles.push_back(FOP);

      fieldOnGrid.clear();

      gamma = sqrt(1 + pow(particles[i].px, 2) + pow(particles[i].py, 2) + pow(particles[i].pz, 2));

      // PASSO #1
      //RIEMPIE dx = vx, vy, vz, Fx, Fy, Fz
      dx.setParticleX(particles[i].px / (gamma * particles[i].m));
      dx.setParticleY(particles[i].py / (gamma * particles[i].m));
      dx.setParticleZ(particles[i].pz / (gamma * particles[i].m));
      dx.setParticlePX(particles[i].q * fieldOnParticles[0].ex + particles[i].q * dx.y * fieldOnParticles[0].bz - particles[i].q * dx.z * fieldOnParticles[0].by);
      dx.setParticlePY(particles[i].q * fieldOnParticles[0].ey + particles[i].q * dx.z * fieldOnParticles[0].bx - particles[i].q * dx.x * fieldOnParticles[0].bz);
      dx.setParticlePZ(particles[i].q * fieldOnParticles[0].ez + particles[i].q * dx.x * fieldOnParticles[0].by - particles[i].q * dx.y * fieldOnParticles[0].bx);

      //RIEMPIE dx2 per il metodo RK4 con il dx del primo passo
      dx2 = dx;

      //RIEMPIE x2 con le posizioni avanzate di particles[i] dopo il primo passo
      x2.setParticleX(particles[i].x + dx.x * (data.dt * 0.5));
      x2.setParticleY(particles[i].y + dx.y * (data.dt * 0.5));
      x2.setParticleZ(particles[i].z + dx.z * (data.dt * 0.5));
      x2.setParticlePX(particles[i].px + dx.px * (data.dt * 0.5));
      x2.setParticlePY(particles[i].py + dx.py * (data.dt * 0.5));
      x2.setParticlePZ(particles[i].pz + dx.pz * (data.dt * 0.5));

      // PASSO #2
      //AVANZA TEMPO DI deltaT/2 e calcola campi nella posizione x2
      gridXX = x2.x / data.deltaX;
      if (gridXX < 0)
        gridX = ((int)gridXX) - 1; //il punto di griglia lungo l'asse x a sinistra della particella
      else
        gridX = (int)gridXX;

      gridYY = x2.y / data.deltaY;
      if (gridYY < 0)
        gridY = ((int)gridYY) - 1; //il punto di griglia lungo l'asse y a sinistra della particella
      else
        gridY = (int)gridYY;

      gridZZ = x2.z / data.deltaZ;
      if (gridZZ < 0)
        gridZ = ((int)gridZZ) - 1; //il punto di griglia lungo l'asse z a sinistra della particella
      else
        gridZ = (int)gridZZ;

      if (gridX < 0 || gridY < 0 || gridZ < 0) {
        cout << "La particella e` uscita dai boundaries - simulazione interrotta" << endl;
        return;
      }

      data.aumentaT(data.dt * 0.5);

      campoSuPuntiGriglia.clear();
      riempiPuntiGrigliaConCampoAnalitico1D(data, campoSuPuntiGriglia);

      fieldOnGrid.push_back(campoSuPuntiGriglia.at(gridY * data.nGridPointsX + gridX)); // primo punto
      fieldOnGrid.push_back(campoSuPuntiGriglia.at(gridY * data.nGridPointsX + (gridX + 1))); // secondo punto
      fieldOnGrid.push_back(campoSuPuntiGriglia.at((gridY + 1) * data.nGridPointsX + gridX)); // terzo punto
      fieldOnGrid.push_back(campoSuPuntiGriglia.at((gridY + 1) * data.nGridPointsX + (gridX + 1))); // quarto punto

      FOP = interpolation2D_linear(data, x2, fieldOnGrid);
      fieldOnParticles.push_back(FOP);

      fieldOnGrid.clear();

      gamma = sqrt(1 + pow(x2.px, 2) + pow(x2.py, 2) + pow(x2.pz, 2));

      //AGGIORNA dx con i nuovi valori
      dx.setParticleX(x2.px / (gamma * particles[i].m));
      dx.setParticleY(x2.py / (gamma * particles[i].m));
      dx.setParticleZ(x2.pz / (gamma * particles[i].m));
      dx.setParticlePX(particles[i].q * fieldOnParticles[1].ex + particles[i].q * dx.y * fieldOnParticles[1].bz - particles[i].q * dx.z * fieldOnParticles[1].by);
      dx.setParticlePY(particles[i].q * fieldOnParticles[1].ey + particles[i].q * dx.z * fieldOnParticles[1].bx - particles[i].q * dx.x * fieldOnParticles[1].bz);
      dx.setParticlePZ(particles[i].q * fieldOnParticles[1].ez + particles[i].q * dx.x * fieldOnParticles[1].by - particles[i].q * dx.y * fieldOnParticles[1].bx);

      //AGGIUNGE al dx2 i valori calcolati al secondo passo
      dx2.setParticleX(dx2.x + 2 * dx.x);
      dx2.setParticleY(dx2.y + 2 * dx.y);
      dx2.setParticleZ(dx2.z + 2 * dx.z);
      dx2.setParticlePX(dx2.px + 2 * dx.px);
      dx2.setParticlePY(dx2.py + 2 * dx.py);
      dx2.setParticlePZ(dx2.pz + 2 * dx.pz);

      //RIEMPIE x2 con le posizione avanzate di particles[i] dopo il secondo passo
      x2.setParticleX(particles[i].x + dx.x * (data.dt * 0.5));
      x2.setParticleY(particles[i].y + dx.y * (data.dt * 0.5));
      x2.setParticleZ(particles[i].z + dx.z * (data.dt * 0.5));
      x2.setParticlePX(particles[i].px + dx.px * (data.dt * 0.5));
      x2.setParticlePY(particles[i].py + dx.py * (data.dt * 0.5));
      x2.setParticlePZ(particles[i].pz + dx.pz * (data.dt * 0.5));

      // PASSO #3
      //CALCOLA campi nella posizione x2 appena aggiornata
      gridXX = x2.x / data.deltaX;
      if (gridXX < 0)
        gridX = ((int)gridXX) - 1; //il punto di griglia lungo l'asse x a sinistra della particella
      else
        gridX = (int)gridXX;

      gridYY = x2.y / data.deltaY;
      if (gridYY < 0)
        gridY = ((int)gridYY) - 1; //il punto di griglia lungo l'asse y a sinistra della particella
      else
        gridY = (int)gridYY;

      gridZZ = x2.z / data.deltaZ;
      if (gridZZ < 0)
        gridZ = ((int)gridZZ) - 1; //il punto di griglia lungo l'asse z a sinistra della particella
      else
        gridZ = (int)gridZZ;

      if (gridX < 0 || gridY < 0 || gridZ < 0) {
        cout << "La particella e` uscita dai boundaries - simulazione interrotta" << endl;
        return;
      }

      //      riempiPuntiGrigliaConCampoAnalitico1D(data, campoSuPuntiGriglia);   // non e` necessario ricalcolarlo perche' il tempo non e` variato

      fieldOnGrid.push_back(campoSuPuntiGriglia.at(gridY * data.nGridPointsX + gridX)); // primo punto
      fieldOnGrid.push_back(campoSuPuntiGriglia.at(gridY * data.nGridPointsX + (gridX + 1))); // secondo punto
      fieldOnGrid.push_back(campoSuPuntiGriglia.at((gridY + 1) * data.nGridPointsX + gridX)); // terzo punto
      fieldOnGrid.push_back(campoSuPuntiGriglia.at((gridY + 1) * data.nGridPointsX + (gridX + 1))); // quarto punto

      FOP = interpolation2D_linear(data, x2, fieldOnGrid);
      fieldOnParticles.push_back(FOP);

      fieldOnGrid.clear();

      gamma = sqrt(1 + pow(x2.px, 2) + pow(x2.py, 2) + pow(x2.pz, 2));

      //AGGIORNA dx con i nuovi valori
      dx.setParticleX(x2.px / (gamma * particles[i].m));
      dx.setParticleY(x2.py / (gamma * particles[i].m));
      dx.setParticleZ(x2.pz / (gamma * particles[i].m));
      dx.setParticlePX(particles[i].q * fieldOnParticles[2].ex + particles[i].q * dx.y * fieldOnParticles[2].bz - particles[i].q * dx.z * fieldOnParticles[2].by);
      dx.setParticlePY(particles[i].q * fieldOnParticles[2].ey + particles[i].q * dx.z * fieldOnParticles[2].bx - particles[i].q * dx.x * fieldOnParticles[2].bz);
      dx.setParticlePZ(particles[i].q * fieldOnParticles[2].ez + particles[i].q * dx.x * fieldOnParticles[2].by - particles[i].q * dx.y * fieldOnParticles[2].bx);

      //AGGIUNGE al dx2 i valori calcolati al terzo passo
      dx2.setParticleX(dx2.x + 2 * dx.x);
      dx2.setParticleY(dx2.y + 2 * dx.y);
      dx2.setParticleZ(dx2.z + 2 * dx.z);
      dx2.setParticlePX(dx2.px + 2 * dx.px);
      dx2.setParticlePY(dx2.py + 2 * dx.py);
      dx2.setParticlePZ(dx2.pz + 2 * dx.pz);

      //RIEMPIE x2 con le posizione avanzate di particles[i]
      x2.setParticleX(particles[i].x + dx.x * data.dt);
      x2.setParticleY(particles[i].y + dx.y * data.dt);
      x2.setParticleZ(particles[i].z + dx.z * data.dt);
      x2.setParticlePX(particles[i].px + dx.px * data.dt);
      x2.setParticlePY(particles[i].py + dx.py * data.dt);
      x2.setParticlePZ(particles[i].pz + dx.pz * data.dt);

      // PASSO #4
      //AVANZA TEMPO DI deltaT/2 e calcola campi nella posizione x2
      gridXX = x2.x / data.deltaX;
      if (gridXX < 0)
        gridX = ((int)gridXX) - 1; //il punto di griglia lungo l'asse x a sinistra della particella
      else
        gridX = (int)gridXX;

      gridYY = x2.y / data.deltaY;
      if (gridYY < 0)
        gridY = ((int)gridYY) - 1; //il punto di griglia lungo l'asse y a sinistra della particella
      else
        gridY = (int)gridYY;

      gridZZ = x2.z / data.deltaZ;
      if (gridZZ < 0)
        gridZ = ((int)gridZZ) - 1; //il punto di griglia lungo l'asse z a sinistra della particella
      else
        gridZ = (int)gridZZ;

      if (gridX < 0 || gridY < 0 || gridZ < 0) {
        cout << "La particella e` uscita dai boundaries - simulazione interrotta" << endl;
        return;
      }

      data.aumentaT(data.dt * 0.5);

      campoSuPuntiGriglia.clear();
      riempiPuntiGrigliaConCampoAnalitico1D(data, campoSuPuntiGriglia);

      fieldOnGrid.push_back(campoSuPuntiGriglia.at(gridY * data.nGridPointsX + gridX)); // primo punto
      fieldOnGrid.push_back(campoSuPuntiGriglia.at(gridY * data.nGridPointsX + (gridX + 1))); // secondo punto
      fieldOnGrid.push_back(campoSuPuntiGriglia.at((gridY + 1) * data.nGridPointsX + gridX)); // terzo punto
      fieldOnGrid.push_back(campoSuPuntiGriglia.at((gridY + 1) * data.nGridPointsX + (gridX + 1))); // quarto punto

      FOP = interpolation2D_linear(data, x2, fieldOnGrid);
      fieldOnParticles.push_back(FOP);

      fieldOnGrid.clear();

      gamma = sqrt(1 + pow(x2.px, 2) + pow(x2.py, 2) + pow(x2.pz, 2));

      //AGGIORNA dx con i nuovi valori
      dx.setParticleX(x2.px / (gamma * particles[i].m));
      dx.setParticleY(x2.py / (gamma * particles[i].m));
      dx.setParticleZ(x2.pz / (gamma * particles[i].m));
      dx.setParticlePX(particles[i].q * fieldOnParticles[3].ex + particles[i].q * dx.y * fieldOnParticles[3].bz - particles[i].q * dx.z * fieldOnParticles[3].by);
      dx.setParticlePY(particles[i].q * fieldOnParticles[3].ey + particles[i].q * dx.z * fieldOnParticles[3].bx - particles[i].q * dx.x * fieldOnParticles[3].bz);
      dx.setParticlePZ(particles[i].q * fieldOnParticles[3].ez + particles[i].q * dx.x * fieldOnParticles[3].by - particles[i].q * dx.y * fieldOnParticles[3].bx);

      //AGGIUNGE al dx2 i valori calcolati al terzo passo
      dx2.setParticleX(dx2.x + dx.x);
      dx2.setParticleY(dx2.y + dx.y);
      dx2.setParticleZ(dx2.z + dx.z);
      dx2.setParticlePX(dx2.px + dx.px);
      dx2.setParticlePY(dx2.py + dx.py);
      dx2.setParticlePZ(dx2.pz + dx.pz);

      particles[i].setParticleX(particles[i].x + (data.dt / 6.) * dx2.x);
      particles[i].setParticleY(particles[i].y + (data.dt / 6.) * dx2.y);
      particles[i].setParticleZ(particles[i].z + (data.dt / 6.) * dx2.z);
      particles[i].setParticlePX(particles[i].px + (data.dt / 6.) * dx2.px);
      particles[i].setParticlePY(particles[i].py + (data.dt / 6.) * dx2.py);
      particles[i].setParticlePZ(particles[i].pz + (data.dt / 6.) * dx2.pz);

      data.aumentaT(-data.dt); //faccio tornare il tempo indietro cosi` gli altri elettroni di questo step non hanno il tempo gia` avanzato
          //il tempo avanza definitivamente in ciascun ciclo for esterno (sugli nsteps)

      //INTEGRALI PRIMI da sistemare perche' calcolano il potenziale sulla particella e non sulla griglia e poi interpolando
      potenziale.clear();
      zeta = particles[i].x - SPEED_OF_LIGHT * data.t; // non usato nella mia struttura delle equazioni di campo
      gamma = sqrt(1 + pow(particles[i].px, 2) + pow(particles[i].py, 2) + pow(particles[i].pz, 2));

      gridXX = particles[i].x / data.deltaX;
      if (gridXX < 0)
        gridX = ((int)gridXX) - 1; //il punto di griglia lungo l'asse x a sinistra della particella
      else
        gridX = (int)gridXX;

      gridYY = particles[i].y / data.deltaY;
      if (gridYY < 0)
        gridY = ((int)gridYY) - 1; //il punto di griglia lungo l'asse y a sinistra della particella
      else
        gridY = (int)gridYY;

      gridZZ = particles[i].z / data.deltaZ;
      if (gridZZ < 0)
        gridZ = ((int)gridZZ) - 1; //il punto di griglia lungo l'asse z a sinistra della particella
      else
        gridZ = (int)gridZZ;

      if (gridX < 0 || gridY < 0 || gridZ < 0) {
        cout << "La particella e` uscita dai boundaries - simulazione interrotta" << endl;
        return;
      }

      campoSuPuntiGriglia.clear();
      riempiPuntiGrigliaConPotenzialeAnalitico1D(data, campoSuPuntiGriglia);

      fieldOnGrid.push_back(campoSuPuntiGriglia.at(gridY * data.nGridPointsX + gridX)); // primo punto
      fieldOnGrid.push_back(campoSuPuntiGriglia.at(gridY * data.nGridPointsX + (gridX + 1))); // secondo punto
      fieldOnGrid.push_back(campoSuPuntiGriglia.at((gridY + 1) * data.nGridPointsX + gridX)); // terzo punto
      fieldOnGrid.push_back(campoSuPuntiGriglia.at((gridY + 1) * data.nGridPointsX + (gridX + 1))); // quarto punto

      FOP = interpolation2D_linear(data, particles[i], fieldOnGrid);
      potenziale.push_back(FOP);

      fieldOnGrid.clear();

      I1 = (particles[i].py / (particles[i].m * SPEED_OF_LIGHT)) + potenziale[0].ey;
      I2 = particles[i].pz / (particles[i].m * SPEED_OF_LIGHT);
      I3 = gamma - particles[i].px / (particles[i].m * SPEED_OF_LIGHT);

      outputstream
          << setw(12) << i + 1 << "\t"
          << setprecision(14) << setw(16) << particles[i].x << "\t"
          << setw(16) << particles[i].y << "\t"
          << setw(16) << particles[i].z << "\t"
          //<< setw(16) << particles[i].getParticleCell() << "\t"
          << setw(16) << particles[i].px << "\t"
          << setw(16) << particles[i].py << "\t"
          << setw(16) << particles[i].pz << "\t"
          << setw(16) << zeta << "\t"
          << setw(16) << potenziale[0].ey << "\t"
          << setw(16) << I1 << "\t"
          << setw(16) << I2 << "\t"
          << setw(16) << I3 << endl;
    }
  }
}

Field obtainFOP(Data data, Particle particle, vector<Field>& FOG)
{
  size_t numberOfVertices = FOG.size();

  /*
  Se il numero di vertici (punti di griglia) dai quali bisogna interpolare il campo sulla particella sono:
  2 -> il primo e` alla sua sinistra, il secondo alla sua destra
  4 -> il primo e` il piu` vicino all'origine degli assi, il secondo si sposta lungo x, il terzo e` il primo spostato lungo y e l'ultimo e` il piu` lontano dagli assi
  8 -> con la stessa filosofia di prima, i primi 4 sono nello stesso ordine e sono quelli appartenenti al piano piu` vicino ad xy, gli altri 4 sono piu` lontani di un dz
  */

  Particle relativeCoord;
  double gridd;
  int grid;

  gridd = particle.x / data.deltaX;
  if (gridd < 0)
    grid = ((int)gridd) - 1;
  else
    grid = (int)gridd;
  relativeCoord.setParticleX(particle.x - data.deltaX * grid);

  gridd = particle.y / data.deltaY;
  if (gridd < 0)
    grid = ((int)gridd) - 1;
  else
    grid = (int)gridd;
  relativeCoord.setParticleY(particle.y - data.deltaY * grid);

  gridd = particle.z / data.deltaZ;
  if (gridd < 0)
    grid = ((int)gridd) - 1;
  else
    grid = (int)gridd;
  relativeCoord.setParticleZ(particle.z - data.deltaZ * grid);

  Field tempField;

  double x2_vsFirstGridPoint = (pow(relativeCoord.x, 2));
  double y2_vsFirstGridPoint = (pow(relativeCoord.y, 2));
  double z2_vsFirstGridPoint = (pow(relativeCoord.z, 2));
  double pos2_vsFirstGridPoint = x2_vsFirstGridPoint + y2_vsFirstGridPoint + z2_vsFirstGridPoint;

  double x2_vsSecondGridPoint = (pow(data.deltaX - relativeCoord.x, 2));
  double y2_vsSecondGridPoint = (pow(relativeCoord.y, 2));
  double z2_vsSecondGridPoint = (pow(relativeCoord.z, 2));
  double pos2_vsSecondGridPoint = x2_vsSecondGridPoint + y2_vsSecondGridPoint + z2_vsSecondGridPoint;

  double normal = pos2_vsFirstGridPoint + pos2_vsSecondGridPoint;

  double x2_vsThirdGridPoint;
  double y2_vsThirdGridPoint;
  double z2_vsThirdGridPoint;
  double pos2_vsThirdGridPoint;
  double x2_vsFourthGridPoint;
  double y2_vsFourthGridPoint;
  double z2_vsFourthGridPoint;
  double pos2_vsFourthGridPoint;
  double x2_vsFifthGridPoint;
  double y2_vsFifthGridPoint;
  double z2_vsFifthGridPoint;
  double pos2_vsFifthGridPoint;
  double x2_vsSixthGridPoint;
  double y2_vsSixthGridPoint;
  double z2_vsSixthGridPoint;
  double pos2_vsSixthGridPoint;
  double x2_vsSeventhGridPoint;
  double y2_vsSeventhGridPoint;
  double z2_vsSeventhGridPoint;
  double pos2_vsSeventhGridPoint;
  double x2_vsEighthGridPoint;
  double y2_vsEighthGridPoint;
  double z2_vsEighthGridPoint;
  double pos2_vsEighthGridPoint;

  if (numberOfVertices == 4 || numberOfVertices == 8) {
    x2_vsThirdGridPoint = (pow(relativeCoord.x, 2));
    y2_vsThirdGridPoint = (pow(data.deltaY - relativeCoord.z, 2));
    z2_vsThirdGridPoint = (pow(relativeCoord.z, 2));
    pos2_vsThirdGridPoint = x2_vsThirdGridPoint + y2_vsThirdGridPoint + z2_vsThirdGridPoint;

    x2_vsFourthGridPoint = (pow(data.deltaX - relativeCoord.x, 2));
    y2_vsFourthGridPoint = (pow(data.deltaY - relativeCoord.y, 2));
    z2_vsFourthGridPoint = (pow(relativeCoord.z, 2));
    pos2_vsFourthGridPoint = x2_vsFourthGridPoint + y2_vsFourthGridPoint + z2_vsFourthGridPoint;

    normal += pos2_vsThirdGridPoint + pos2_vsFourthGridPoint;
  }

  if (numberOfVertices == 8) {
    x2_vsFifthGridPoint = (pow(relativeCoord.x, 2));
    y2_vsFifthGridPoint = (pow(relativeCoord.y, 2));
    z2_vsFifthGridPoint = (pow(data.deltaZ - relativeCoord.z, 2));
    pos2_vsFifthGridPoint = x2_vsFifthGridPoint + y2_vsFifthGridPoint + z2_vsFifthGridPoint;

    x2_vsSixthGridPoint = (pow(data.deltaX - relativeCoord.x, 2));
    y2_vsSixthGridPoint = (pow(relativeCoord.y, 2));
    z2_vsSixthGridPoint = (pow(data.deltaZ - relativeCoord.z, 2));
    pos2_vsSixthGridPoint = x2_vsSixthGridPoint + y2_vsSixthGridPoint + z2_vsSixthGridPoint;

    x2_vsSeventhGridPoint = (pow(relativeCoord.x, 2));
    y2_vsSeventhGridPoint = (pow(data.deltaY - relativeCoord.z, 2));
    z2_vsSeventhGridPoint = (pow(data.deltaZ - relativeCoord.z, 2));
    pos2_vsSeventhGridPoint = x2_vsSeventhGridPoint + y2_vsSeventhGridPoint + z2_vsSeventhGridPoint;

    x2_vsEighthGridPoint = (pow(data.deltaX - relativeCoord.x, 2));
    y2_vsEighthGridPoint = (pow(data.deltaY - relativeCoord.y, 2));
    z2_vsEighthGridPoint = (pow(data.deltaZ - relativeCoord.z, 2));
    pos2_vsEighthGridPoint = x2_vsEighthGridPoint + y2_vsEighthGridPoint + z2_vsEighthGridPoint;

    normal += pos2_vsFifthGridPoint + pos2_vsSixthGridPoint + pos2_vsSeventhGridPoint + pos2_vsEighthGridPoint;
  }

  double normal_inv;
  normal_inv = 1.0 / normal;

  pos2_vsFirstGridPoint *= normal_inv;
  pos2_vsSecondGridPoint *= normal_inv;

  double pos_vsFirstGridPoint = sqrt(pos2_vsFirstGridPoint);
  double pos_vsSecondGridPoint = sqrt(pos2_vsSecondGridPoint);

  double fieldX_0 = FOG[0].ex * pos_vsFirstGridPoint;
  double fieldX_1 = FOG[1].ex * pos_vsSecondGridPoint;
  double fieldY_0 = FOG[0].ey * pos_vsFirstGridPoint;
  double fieldY_1 = FOG[1].ey * pos_vsSecondGridPoint;
  double fieldZ_0 = FOG[0].ez * pos_vsFirstGridPoint;
  double fieldZ_1 = FOG[1].ez * pos_vsSecondGridPoint;

  double fieldBX_0 = FOG[0].bx * pos_vsFirstGridPoint;
  double fieldBX_1 = FOG[1].bx * pos_vsSecondGridPoint;
  double fieldBY_0 = FOG[0].by * pos_vsFirstGridPoint;
  double fieldBY_1 = FOG[1].by * pos_vsSecondGridPoint;
  double fieldBZ_0 = FOG[0].bz * pos_vsFirstGridPoint;
  double fieldBZ_1 = FOG[1].bz * pos_vsSecondGridPoint;

  tempField.setEx(fieldX_0 + fieldX_1);
  tempField.setEy(fieldY_0 + fieldY_1);
  tempField.setEz(fieldZ_0 + fieldZ_1);
  tempField.setBx(fieldBX_0 + fieldBX_1);
  tempField.setBy(fieldBY_0 + fieldBY_1);
  tempField.setBz(fieldBZ_0 + fieldBZ_1);

  if (numberOfVertices == 4) {
    pos2_vsThirdGridPoint *= normal_inv;
    pos2_vsFourthGridPoint *= normal_inv;

    double pos_vsThirdGridPoint = sqrt(pos2_vsThirdGridPoint);
    double pos_vsFourthGridPoint = sqrt(pos2_vsFourthGridPoint);

    double fieldX_2 = FOG[2].ex * pos_vsThirdGridPoint;
    double fieldX_3 = FOG[3].ex * pos_vsFourthGridPoint;
    tempField.setEx(fieldX_0 + fieldX_1 + fieldX_2 + fieldX_3);
    double fieldBX_2 = FOG[2].bx * pos_vsThirdGridPoint;
    double fieldBX_3 = FOG[3].bx * pos_vsFourthGridPoint;
    tempField.setBx(fieldBX_0 + fieldBX_1 + fieldBX_2 + fieldBX_3);

    double fieldY_2 = FOG[2].ey * pos_vsThirdGridPoint;
    double fieldY_3 = FOG[3].ey * pos_vsFourthGridPoint;
    tempField.setEy(fieldY_0 + fieldY_1 + fieldY_2 + fieldY_3);
    double fieldBY_2 = FOG[2].by * pos_vsThirdGridPoint;
    double fieldBY_3 = FOG[3].by * pos_vsFourthGridPoint;
    tempField.setBy(fieldBY_0 + fieldBY_1 + fieldBY_2 + fieldBY_3);

    double fieldZ_2 = FOG[2].ez * pos_vsThirdGridPoint;
    double fieldZ_3 = FOG[3].ez * pos_vsFourthGridPoint;
    tempField.setEz(fieldZ_0 + fieldZ_1 + fieldZ_2 + fieldZ_3);
    double fieldBZ_2 = FOG[2].bz * pos_vsThirdGridPoint;
    double fieldBZ_3 = FOG[3].bz * pos_vsFourthGridPoint;
    tempField.setBz(fieldBZ_0 + fieldBZ_1 + fieldBZ_2 + fieldBZ_3);
  }

  if (numberOfVertices == 8) {
    pos2_vsThirdGridPoint *= normal_inv;
    pos2_vsFourthGridPoint *= normal_inv;
    pos2_vsFifthGridPoint *= normal_inv;
    pos2_vsSixthGridPoint *= normal_inv;
    pos2_vsSeventhGridPoint *= normal_inv;
    pos2_vsEighthGridPoint *= normal_inv;

    double pos_vsThirdGridPoint = sqrt(pos2_vsThirdGridPoint);
    double pos_vsFourthGridPoint = sqrt(pos2_vsFourthGridPoint);
    double pos_vsFifthGridPoint = sqrt(pos2_vsFifthGridPoint);
    double pos_vsSixthGridPoint = sqrt(pos2_vsSixthGridPoint);
    double pos_vsSeventhGridPoint = sqrt(pos2_vsSeventhGridPoint);
    double pos_vsEighthGridPoint = sqrt(pos2_vsEighthGridPoint);

    double fieldX_2 = FOG[2].ex * pos_vsThirdGridPoint;
    double fieldX_3 = FOG[3].ex * pos_vsFourthGridPoint;
    double fieldX_4 = FOG[4].ex * pos_vsFifthGridPoint;
    double fieldX_5 = FOG[5].ex * pos_vsSixthGridPoint;
    double fieldX_6 = FOG[6].ex * pos_vsSeventhGridPoint;
    double fieldX_7 = FOG[7].ex * pos_vsEighthGridPoint;
    tempField.setEx(fieldX_0 + fieldX_1 + fieldX_2 + fieldX_3 + fieldX_4 + fieldX_5 + fieldX_6 + fieldX_7);
    double fieldBX_2 = FOG[2].bx * pos_vsThirdGridPoint;
    double fieldBX_3 = FOG[3].bx * pos_vsFourthGridPoint;
    double fieldBX_4 = FOG[4].bx * pos_vsFifthGridPoint;
    double fieldBX_5 = FOG[5].bx * pos_vsSixthGridPoint;
    double fieldBX_6 = FOG[6].bx * pos_vsSeventhGridPoint;
    double fieldBX_7 = FOG[7].bx * pos_vsEighthGridPoint;
    tempField.setBx(fieldBX_0 + fieldBX_1 + fieldBX_2 + fieldBX_3 + fieldBX_4 + fieldBX_5 + fieldBX_6 + fieldBX_7);

    double fieldY_2 = FOG[2].ey * pos_vsThirdGridPoint;
    double fieldY_3 = FOG[3].ey * pos_vsFourthGridPoint;
    double fieldY_4 = FOG[4].ey * pos_vsFifthGridPoint;
    double fieldY_5 = FOG[5].ey * pos_vsSixthGridPoint;
    double fieldY_6 = FOG[6].ey * pos_vsSeventhGridPoint;
    double fieldY_7 = FOG[7].ey * pos_vsEighthGridPoint;
    tempField.setEy(fieldY_0 + fieldY_1 + fieldY_2 + fieldY_3 + fieldY_4 + fieldY_5 + fieldY_6 + fieldY_7);
    double fieldBY_2 = FOG[2].by * pos_vsThirdGridPoint;
    double fieldBY_3 = FOG[3].by * pos_vsFourthGridPoint;
    double fieldBY_4 = FOG[4].by * pos_vsFifthGridPoint;
    double fieldBY_5 = FOG[5].by * pos_vsSixthGridPoint;
    double fieldBY_6 = FOG[6].by * pos_vsSeventhGridPoint;
    double fieldBY_7 = FOG[7].by * pos_vsEighthGridPoint;
    tempField.setBy(fieldBY_0 + fieldBY_1 + fieldBY_2 + fieldBY_3 + fieldBY_4 + fieldBY_5 + fieldBY_6 + fieldBY_7);

    double fieldZ_2 = FOG[2].ez * pos_vsThirdGridPoint;
    double fieldZ_3 = FOG[3].ez * pos_vsFourthGridPoint;
    double fieldZ_4 = FOG[4].ez * pos_vsFifthGridPoint;
    double fieldZ_5 = FOG[5].ez * pos_vsSixthGridPoint;
    double fieldZ_6 = FOG[6].ez * pos_vsSeventhGridPoint;
    double fieldZ_7 = FOG[7].ez * pos_vsEighthGridPoint;
    tempField.setEz(fieldZ_0 + fieldZ_1 + fieldZ_2 + fieldZ_3 + fieldZ_4 + fieldZ_5 + fieldZ_6 + fieldZ_7);
    double fieldBZ_2 = FOG[2].bz * pos_vsThirdGridPoint;
    double fieldBZ_3 = FOG[3].bz * pos_vsFourthGridPoint;
    double fieldBZ_4 = FOG[4].bz * pos_vsFifthGridPoint;
    double fieldBZ_5 = FOG[5].bz * pos_vsSixthGridPoint;
    double fieldBZ_6 = FOG[6].bz * pos_vsSeventhGridPoint;
    double fieldBZ_7 = FOG[7].bz * pos_vsEighthGridPoint;
    tempField.setBz(fieldBZ_0 + fieldBZ_1 + fieldBZ_2 + fieldBZ_3 + fieldBZ_4 + fieldBZ_5 + fieldBZ_6 + fieldBZ_7);
  }

  /*
  PROTOTIPO
    - in input deve avere i dati riguardo la griglia e riguardo la posizione della particella, non un vettore ma solo 1 particella e i campi in quell'istante
    - da essi si calcola la distanza relativa tra i punti di griglia e la posizione della particella
    - in funzione di cio`, deve calcolare l'interpolazione dei campi sulla particella
  */
  return tempField;
}

Field interpolation2D_linear(Data data, Particle position, vector<Field>& valuesOnGrid)
{
  /*
  Se il numero di vertici (punti di griglia) dai quali bisogna interpolare il campo sulla particella sono diversi da 4, la funzione esce restituendo -1
  Se i vertici sono 4, vanno identificati in questo modo: il primo e` il piu` vicino all'origine degli assi, il secondo si sposta lungo x,
  il terzo e` il primo spostato lungo y e l'ultimo e` il piu` lontano dagli assi
  */

  Particle relativeCoord;

  double gridd;
  int grid = 0;
  double dx_inv = 1. / data.deltaX;
  double dy_inv = 1. / data.deltaY;

  gridd = position.x / data.deltaX;
  if (gridd < 0)
    grid = ((int)grid) - 1;
  else
    grid = (int)gridd;
  relativeCoord.setParticleX((position.x - data.deltaX * grid) * dx_inv);

  gridd = position.y / data.deltaY;
  if (gridd < 0)
    grid = ((int)grid) - 1;
  else
    grid = (int)gridd;
  relativeCoord.setParticleY((position.y - data.deltaY * grid) * dy_inv);

  relativeCoord.setParticleZ(0.);

  Field tempField;

  double weightOnFourthPoint = relativeCoord.x * relativeCoord.y;
  double weightOnFirstPoint = (1 - relativeCoord.x) * (1 - relativeCoord.y);
  double weightOnThirdPoint = relativeCoord.x * (1 - relativeCoord.y);
  double weightOnSecondPoint = (1 - relativeCoord.x) * relativeCoord.y;

  double fieldX_0 = valuesOnGrid[0].ex * weightOnFirstPoint;
  double fieldX_1 = valuesOnGrid[1].ex * weightOnSecondPoint;
  double fieldX_2 = valuesOnGrid[2].ex * weightOnThirdPoint;
  double fieldX_3 = valuesOnGrid[3].ex * weightOnFourthPoint;

  double fieldY_0 = valuesOnGrid[0].ey * weightOnFirstPoint;
  double fieldY_1 = valuesOnGrid[1].ey * weightOnSecondPoint;
  double fieldY_2 = valuesOnGrid[2].ey * weightOnThirdPoint;
  double fieldY_3 = valuesOnGrid[3].ey * weightOnFourthPoint;

  double fieldZ_0 = valuesOnGrid[0].ez * weightOnFirstPoint;
  double fieldZ_1 = valuesOnGrid[1].ez * weightOnSecondPoint;
  double fieldZ_2 = valuesOnGrid[2].ez * weightOnThirdPoint;
  double fieldZ_3 = valuesOnGrid[3].ez * weightOnFourthPoint;

  double fieldBX_0 = valuesOnGrid[0].bx * weightOnFirstPoint;
  double fieldBX_1 = valuesOnGrid[1].bx * weightOnSecondPoint;
  double fieldBX_2 = valuesOnGrid[2].bx * weightOnThirdPoint;
  double fieldBX_3 = valuesOnGrid[3].bx * weightOnFourthPoint;

  double fieldBY_0 = valuesOnGrid[0].by * weightOnFirstPoint;
  double fieldBY_1 = valuesOnGrid[1].by * weightOnSecondPoint;
  double fieldBY_2 = valuesOnGrid[2].by * weightOnThirdPoint;
  double fieldBY_3 = valuesOnGrid[3].by * weightOnFourthPoint;

  double fieldBZ_0 = valuesOnGrid[0].bz * weightOnFirstPoint;
  double fieldBZ_1 = valuesOnGrid[1].bz * weightOnSecondPoint;
  double fieldBZ_2 = valuesOnGrid[2].bz * weightOnThirdPoint;
  double fieldBZ_3 = valuesOnGrid[3].bz * weightOnFourthPoint;

  tempField.setEx(fieldX_0 + fieldX_1 + fieldX_2 + fieldX_3);
  tempField.setBx(fieldBX_0 + fieldBX_1 + fieldBX_2 + fieldBX_3);
  tempField.setEy(fieldY_0 + fieldY_1 + fieldY_2 + fieldY_3);
  tempField.setBy(fieldBY_0 + fieldBY_1 + fieldBY_2 + fieldBY_3);
  tempField.setEz(fieldZ_0 + fieldZ_1 + fieldZ_2 + fieldZ_3);
  tempField.setBz(fieldBZ_0 + fieldBZ_1 + fieldBZ_2 + fieldBZ_3);

  return tempField;
}

Field interpolation2D_quadratic(Data data, Particle position, vector<Field>& valuesOnGrid)
{
  /*
  1 2 3  |--->  x
  4 5 6  |
  7 8 9  V
         y
  */

  Particle relativeCoord;

  double gridd;
  int grid = 0;
  double dx_inv = 1 / data.deltaX;
  double dy_inv = 1 / data.deltaY;

  gridd = (position.x / data.deltaX) + 0.5;
  if (gridd < 0)
    grid = ((int)grid) - 1;
  else
    grid = (int)gridd;
  relativeCoord.setParticleX((position.x - data.deltaX * grid) * dx_inv);

  gridd = (position.y / data.deltaY) + 0.5;
  if (gridd < 0)
    grid = ((int)grid) - 1;
  else
    grid = (int)gridd;
  relativeCoord.setParticleY((position.y - data.deltaY * grid) * dy_inv);

  relativeCoord.setParticleZ(0.);

  Field tempField;

  double weightOnFirstPoint = 0.75 - pow(relativeCoord.x, 2.);
  double weightOnSecondPoint = 0.5 * (0.25 + pow(relativeCoord.x, 2.) + relativeCoord.x);
  double weightOnThirdPoint = 0.5 * (0.25 + pow(relativeCoord.x, 2.) - relativeCoord.x);
  double weightOnFourthPoint = 0.;
  double weightOnFifthPoint = 0.;
  double weightOnSixthPoint = 0.;
  double weightOnSeventhPoint = 0.;
  double weightOnEighthPoint = 0.;
  double weightOnNinethPoint = 0.;

  double fieldX_0 = valuesOnGrid[0].ex * weightOnFirstPoint;
  double fieldX_1 = valuesOnGrid[1].ex * weightOnSecondPoint;
  double fieldX_2 = valuesOnGrid[2].ex * weightOnThirdPoint;
  double fieldX_3 = valuesOnGrid[3].ex * weightOnFourthPoint;
  double fieldX_4 = valuesOnGrid[4].ex * weightOnFifthPoint;
  double fieldX_5 = valuesOnGrid[5].ex * weightOnSixthPoint;
  double fieldX_6 = valuesOnGrid[6].ex * weightOnSeventhPoint;
  double fieldX_7 = valuesOnGrid[7].ex * weightOnEighthPoint;
  double fieldX_8 = valuesOnGrid[8].ex * weightOnNinethPoint;

  double fieldY_0 = valuesOnGrid[0].ey * weightOnFirstPoint;
  double fieldY_1 = valuesOnGrid[1].ey * weightOnSecondPoint;
  double fieldY_2 = valuesOnGrid[2].ey * weightOnThirdPoint;
  double fieldY_3 = valuesOnGrid[3].ey * weightOnFourthPoint;
  double fieldY_4 = valuesOnGrid[4].ey * weightOnFifthPoint;
  double fieldY_5 = valuesOnGrid[5].ey * weightOnSixthPoint;
  double fieldY_6 = valuesOnGrid[6].ey * weightOnSeventhPoint;
  double fieldY_7 = valuesOnGrid[7].ey * weightOnEighthPoint;
  double fieldY_8 = valuesOnGrid[8].ey * weightOnNinethPoint;

  double fieldZ_0 = valuesOnGrid[0].ez * weightOnFirstPoint;
  double fieldZ_1 = valuesOnGrid[1].ez * weightOnSecondPoint;
  double fieldZ_2 = valuesOnGrid[2].ez * weightOnThirdPoint;
  double fieldZ_3 = valuesOnGrid[3].ez * weightOnFourthPoint;
  double fieldZ_4 = valuesOnGrid[4].ez * weightOnFifthPoint;
  double fieldZ_5 = valuesOnGrid[5].ez * weightOnSixthPoint;
  double fieldZ_6 = valuesOnGrid[6].ez * weightOnSeventhPoint;
  double fieldZ_7 = valuesOnGrid[7].ez * weightOnEighthPoint;
  double fieldZ_8 = valuesOnGrid[8].ez * weightOnNinethPoint;

  double fieldBX_0 = valuesOnGrid[0].bx * weightOnFirstPoint;
  double fieldBX_1 = valuesOnGrid[1].bx * weightOnSecondPoint;
  double fieldBX_2 = valuesOnGrid[2].bx * weightOnThirdPoint;
  double fieldBX_3 = valuesOnGrid[3].bx * weightOnFourthPoint;
  double fieldBX_4 = valuesOnGrid[4].bx * weightOnFifthPoint;
  double fieldBX_5 = valuesOnGrid[5].bx * weightOnSixthPoint;
  double fieldBX_6 = valuesOnGrid[6].bx * weightOnSeventhPoint;
  double fieldBX_7 = valuesOnGrid[7].bx * weightOnEighthPoint;
  double fieldBX_8 = valuesOnGrid[8].bx * weightOnNinethPoint;

  double fieldBY_0 = valuesOnGrid[0].by * weightOnFirstPoint;
  double fieldBY_1 = valuesOnGrid[1].by * weightOnSecondPoint;
  double fieldBY_2 = valuesOnGrid[2].by * weightOnThirdPoint;
  double fieldBY_3 = valuesOnGrid[3].by * weightOnFourthPoint;
  double fieldBY_4 = valuesOnGrid[4].by * weightOnFifthPoint;
  double fieldBY_5 = valuesOnGrid[5].by * weightOnSixthPoint;
  double fieldBY_6 = valuesOnGrid[6].by * weightOnSeventhPoint;
  double fieldBY_7 = valuesOnGrid[7].by * weightOnEighthPoint;
  double fieldBY_8 = valuesOnGrid[8].by * weightOnNinethPoint;

  double fieldBZ_0 = valuesOnGrid[0].bz * weightOnFirstPoint;
  double fieldBZ_1 = valuesOnGrid[1].bz * weightOnSecondPoint;
  double fieldBZ_2 = valuesOnGrid[2].bz * weightOnThirdPoint;
  double fieldBZ_3 = valuesOnGrid[3].bz * weightOnFourthPoint;
  double fieldBZ_4 = valuesOnGrid[4].bz * weightOnFifthPoint;
  double fieldBZ_5 = valuesOnGrid[5].bz * weightOnSixthPoint;
  double fieldBZ_6 = valuesOnGrid[6].bz * weightOnSeventhPoint;
  double fieldBZ_7 = valuesOnGrid[7].bz * weightOnEighthPoint;
  double fieldBZ_8 = valuesOnGrid[8].bz * weightOnNinethPoint;

  tempField.setEx(fieldX_0 + fieldX_1 + fieldX_2 + fieldX_3 + fieldX_4 + fieldX_5 + fieldX_6 + fieldX_7 + fieldX_8);
  tempField.setBx(fieldBX_0 + fieldBX_1 + fieldBX_2 + fieldBX_3 + fieldBX_4 + fieldBX_5 + fieldBX_6 + fieldBX_7 + fieldBX_8);
  tempField.setEy(fieldY_0 + fieldY_1 + fieldY_2 + fieldY_3 + fieldY_4 + fieldY_5 + fieldY_6 + fieldY_7 + fieldY_8);
  tempField.setBy(fieldBY_0 + fieldBY_1 + fieldBY_2 + fieldBY_3 + fieldBY_4 + fieldBY_5 + fieldBY_6 + fieldBY_7 + fieldBY_8);
  tempField.setEz(fieldZ_0 + fieldZ_1 + fieldZ_2 + fieldZ_3 + fieldZ_4 + fieldZ_5 + fieldZ_6 + fieldZ_7 + fieldZ_8);
  tempField.setBz(fieldBZ_0 + fieldBZ_1 + fieldBZ_2 + fieldBZ_3 + fieldBZ_4 + fieldBZ_5 + fieldBZ_6 + fieldBZ_7 + fieldBZ_8);

  return tempField;
}

void evolveLPF_nogrid(Data data, vector<Particle>& particles, vector<Field>& fieldOnParticles, ofstream& outputstream)
{
  Particle pstar, pHalfAdvanced;

  vector<Field> potenziale;
  Field b;
  Field* campoSuPunto;
  campoSuPunto = new Field[1];

  double gamma, gamma_inv, b2;
  double I1, I2, I3;
  double zeta;

  //ALGORITMO LEAPFROG PARTICELLE

  if (data.n_dim == 1) {
    for (int i = 0; i < data.nSteps; i++) {
      for (int j = 0; j < data.n_electrons; j++) {
        gamma = sqrt(1. + pow(particles[j].px, 2) + pow(particles[j].py, 2) + pow(particles[j].pz, 2));
        gamma_inv = 1. / gamma;

        particles[j].setParticleX(particles[j].x + 0.5 * data.dt * particles[j].px * gamma_inv);
        particles[j].setParticleY(particles[j].y + 0.5 * data.dt * particles[j].py * gamma_inv);
        particles[j].setParticleZ(particles[j].z + 0.5 * data.dt * particles[j].pz * gamma_inv);

        fieldOnParticles.clear(); //SVUOTA il vettore dei campi sulle particelle ad ogni ciclo

        //SICCOME NON HO ANCORA UN CAMPO SU GRIGLIA DEFINITO, TEMPORANEAMENTE CALCOLO IL CAMPO ANALITICO
        //DIRETTAMENTE SULLA PARTICELLA
        calcolaCampoAnalitico1DSuParticella(data, particles[j], *campoSuPunto);
        fieldOnParticles.push_back(*campoSuPunto);

        //pstar = p(n-1/2) + Q * E(n) * DeltaT/2
        pstar.setParticlePX(particles[j].px + 0.5 * data.dt * particles[j].q * fieldOnParticles[0].ex);
        pstar.setParticlePY(particles[j].py + 0.5 * data.dt * particles[j].q * fieldOnParticles[0].ey);
        pstar.setParticlePZ(particles[j].pz + 0.5 * data.dt * particles[j].q * fieldOnParticles[0].ez);
        gamma = sqrt(1. + pow(pstar.px, 2) + pow(pstar.py, 2) + pow(pstar.pz, 2));
        gamma_inv = 1. / gamma;

        //b = (DeltaT/2) * Q * B(n) / gamma
        b.setBx(0.5 * data.dt * fieldOnParticles[0].bx * particles[j].q * gamma_inv);
        b.setBy(0.5 * data.dt * fieldOnParticles[0].by * particles[j].q * gamma_inv);
        b.setBz(0.5 * data.dt * fieldOnParticles[0].bz * particles[j].q * gamma_inv);
        b2 = pow(b.bx, 2) + pow(b.by, 2) + pow(b.bz, 2);
        b2 += 1.;

        /****************************************************************************
        //p(n+1/2) = [pstar + pstar x b + (pstar x b) * b] / b2   SECONDO TURCHETTI  (ma i conti analitici diconon che e` sbagliato, anche perche' l'ultimo termine e` uno scalare!)
        pHalfAdvanced.setParticlePX( (pstar.px + pstar.py * b.bz - b.by * pstar.pz + (pstar.py * b.bz - b.by * pstar.pz) * b.bx ) / b2 );
        pHalfAdvanced.setParticlePY( (pstar.py + pstar.pz * b.bx - b.bz * pstar.px + (pstar.pz * b.bx - b.bz * pstar.px) * b.by ) / b2 );
        pHalfAdvanced.setParticlePZ( (pstar.pz + pstar.px * b.by - b.bx * pstar.py + (pstar.px * b.by - b.bx * pstar.py) * b.bz ) / b2 );
        ****************************************************************************/

        /****************************************************************************
        //p(n+1/2) = [pstar + pstar x b + b x (pstar x b)] / b2   SECONDO LONDRILLO
        pHalfAdvanced.setParticlePX( (pstar.px + pstar.py * b.bz - b.by * pstar.pz + b.bx * pstar.px * b.bx ) / b2 );
        pHalfAdvanced.setParticlePY( (pstar.py + pstar.pz * b.bx - b.bz * pstar.px + b.by * pstar.py * b.by ) / b2 );
        pHalfAdvanced.setParticlePZ( (pstar.pz + pstar.px * b.by - b.bx * pstar.py + b.bz * pstar.pz * b.bz ) / b2 );
        ****************************************************************************/

        //      /****************************************************************************
        //p(n+1/2) = [pstar + pstar x b] / b2   COME IMPLEMENTATO DA TURCHETTI NEL CODICE
        pHalfAdvanced.setParticlePX((pstar.px + pstar.py * b.bz) / b2);
        pHalfAdvanced.setParticlePY((pstar.py - pstar.px * b.bz) / b2);
        pHalfAdvanced.setParticlePZ(0.); //per come ho strutturato i pacchetti d'onda su z non dovrebbe succedere nulla e quindi lo fisso io a zero.
        //      ****************************************************************************/

        particles[j].setParticlePX(2. * pHalfAdvanced.px - particles[j].px);
        particles[j].setParticlePY(2. * pHalfAdvanced.py - particles[j].py);
        particles[j].setParticlePZ(2. * pHalfAdvanced.pz - particles[j].pz);

        gamma = sqrt(1. + pow(particles[j].px, 2) + pow(particles[j].py, 2) + pow(particles[j].pz, 2));
        gamma_inv = 1. / gamma;

        particles[j].setParticleX(particles[j].x + 0.5 * data.dt * particles[j].px * gamma_inv);
        particles[j].setParticleY(particles[j].y + 0.5 * data.dt * particles[j].py * gamma_inv);
        particles[j].setParticleZ(particles[j].z + 0.5 * data.dt * particles[j].pz * gamma_inv);

        //INTEGRALI PRIMI (da sistemare, li ho copiati dall'evoluzione sui campi analitici)
        potenziale.clear();
        zeta = particles[j].x - SPEED_OF_LIGHT * data.t; // non usato nella mia struttura delle equazioni di campo
        calcolaPotenzialeAnalitico1DSuParticella(data, particles[j], *campoSuPunto);
        potenziale.push_back(*campoSuPunto);

        I1 = (particles[j].py / (particles[j].m * SPEED_OF_LIGHT)) + potenziale[0].ey;
        I2 = particles[j].pz / (particles[j].m * SPEED_OF_LIGHT);
        I3 = gamma - particles[j].px / (particles[j].m * SPEED_OF_LIGHT);

        outputstream
            << setw(12) << i + 1 << "\t"
            << setprecision(14) << setw(16) << particles[j].x << "\t"
            << setw(16) << particles[j].y << "\t"
            << setw(16) << particles[j].z << "\t"
            //        << setw(16) << particles[j].getParticleCell() << "\t"
            << setw(16) << particles[j].px << "\t"
            << setw(16) << particles[j].py << "\t"
            << setw(16) << particles[j].pz << "\t"
            << setw(16) << fieldOnParticles[j].ey << "\t"
            << setw(16) << zeta << "\t"
            << setw(16) << I1 << "\t"
            << setw(16) << I2 << "\t"
            << setw(16) << I3 << endl;
      }
      data.aumentaT(data.dt);
    }
  }
}

void evolveLPF_withgrid_onthefly(Data data, vector<Particle>& particles, vector<Field>& fieldOnParticles, ofstream& outputstream)
{
  vector<Field> potenziale;

  Particle pstar, pHalfAdvanced;
  Particle coordinatePrimoPuntoGriglia;
  Particle coordinateSecondoPuntoGriglia;
  Particle coordinateTerzoPuntoGriglia;
  Particle coordinateQuartoPuntoGriglia;
  Particle coordinateQuintoPuntoGriglia;
  Particle coordinateSestoPuntoGriglia;
  Particle coordinateSettimoPuntoGriglia;
  Particle coordinateOttavoPuntoGriglia;

  int posX = 0, posY = 0, posZ = 0;
  Field b;
  Field FOP;
  vector<Field> fieldOnGrid;

  Field* FOG; //DA SISTEMARE, c'è solo per avere un campo fittizio su una griglia on-the-fly per implementare l'algoritmo
  FOG = new Field[1];

  double gridXX, gridYY, gridZZ;
  int gridX, gridY, gridZ;
  double gamma, gamma_inv, b2;
  double I1, I2, I3;
  double zeta;

  Field* campoSuPunto;
  campoSuPunto = new Field[1];

  //ALGORITMO LEAPFROG PARTICELLE

  if (data.n_dim == 1) {
    for (int i = 0; i < data.nSteps; i++) {

      data.aumentaT(data.dt);

      for (int j = 0; j < data.n_electrons; j++) {
        gamma = sqrt(1. + pow(particles[j].px, 2) + pow(particles[j].py, 2) + pow(particles[j].pz, 2));
        gamma_inv = 1. / gamma;

        particles[j].setParticleX(particles[j].x + 0.5 * data.dt * particles[j].px * gamma_inv);
        particles[j].setParticleY(particles[j].y + 0.5 * data.dt * particles[j].py * gamma_inv);
        particles[j].setParticleZ(particles[j].z + 0.5 * data.dt * particles[j].pz * gamma_inv);

        fieldOnParticles.clear();

        //QUI MI CALCOLO I PUNTI DI RIFERIMENTO PER LA PARTICELLA IN ESAME nella posizione attuale
        //forse c'e` un modo piu` furbo per farlo sfruttando il # di cella in cui si trova che e` noto
        gridXX = particles[j].x / data.deltaX;
        if (gridXX < 0)
          gridX = ((int)gridXX) - 1; //il punto di griglia lungo l'asse x a sinistra della particella
        else
          gridX = (int)gridXX;

        gridYY = particles[j].y / data.deltaY;
        if (gridYY < 0)
          gridY = ((int)gridYY) - 1; //il punto di griglia lungo l'asse y a sinistra della particella
        else
          gridY = (int)gridYY;

        gridZZ = particles[j].z / data.deltaZ;
        if (gridZZ < 0)
          gridZ = ((int)gridZZ) - 1; //il punto di griglia lungo l'asse z a sinistra della particella
        else
          gridZ = (int)gridZZ;

        coordinatePrimoPuntoGriglia.setParticleX(gridX * data.deltaX);
        coordinatePrimoPuntoGriglia.setParticleY(gridY * data.deltaY);
        coordinatePrimoPuntoGriglia.setParticleZ(gridZ * data.deltaZ);
        coordinateSecondoPuntoGriglia.setParticleX(gridX * data.deltaX + data.deltaX);
        coordinateSecondoPuntoGriglia.setParticleY(gridY * data.deltaY);
        coordinateSecondoPuntoGriglia.setParticleZ(gridZ * data.deltaZ);
        coordinateTerzoPuntoGriglia.setParticleX(gridX * data.deltaX);
        coordinateTerzoPuntoGriglia.setParticleY(gridY * data.deltaY + data.deltaY);
        coordinateTerzoPuntoGriglia.setParticleZ(gridZ * data.deltaZ);
        coordinateQuartoPuntoGriglia.setParticleX(gridX * data.deltaX + data.deltaX);
        coordinateQuartoPuntoGriglia.setParticleY(gridY * data.deltaY + data.deltaY);
        coordinateQuartoPuntoGriglia.setParticleZ(gridZ * data.deltaZ);
        coordinateQuintoPuntoGriglia.setParticleX(gridX * data.deltaX);
        coordinateQuintoPuntoGriglia.setParticleY(gridY * data.deltaY);
        coordinateQuintoPuntoGriglia.setParticleZ(gridZ * data.deltaZ + data.deltaZ);
        coordinateSestoPuntoGriglia.setParticleX(gridX * data.deltaX + data.deltaX);
        coordinateSestoPuntoGriglia.setParticleY(gridY * data.deltaY);
        coordinateSestoPuntoGriglia.setParticleZ(gridZ * data.deltaZ + data.deltaZ);
        coordinateSettimoPuntoGriglia.setParticleX(gridX * data.deltaX);
        coordinateSettimoPuntoGriglia.setParticleY(gridY * data.deltaY + data.deltaY);
        coordinateSettimoPuntoGriglia.setParticleZ(gridZ * data.deltaZ + data.deltaZ);
        coordinateOttavoPuntoGriglia.setParticleX(gridX * data.deltaX + data.deltaX);
        coordinateOttavoPuntoGriglia.setParticleY(gridY * data.deltaY + data.deltaY);
        coordinateOttavoPuntoGriglia.setParticleZ(gridZ * data.deltaZ + data.deltaZ);

        //SICCOME NON HO ANCORA UN CAMPO SU GRIGLIA DEFINITO, TEMPORANEAMENTE CALCOLO IL CAMPO ANALITICO
        //SUI PUNTI DI RIFERIMENTO DELLA GRIGLIA APPENA OTTENUTI
        //saranno questi i miei valori di En e Bn.
        calcolaCampoAnalitico1DSuParticella(data, coordinatePrimoPuntoGriglia, *FOG);
        fieldOnGrid.push_back(*FOG);
        calcolaCampoAnalitico1DSuParticella(data, coordinateSecondoPuntoGriglia, *FOG);
        fieldOnGrid.push_back(*FOG);
        calcolaCampoAnalitico1DSuParticella(data, coordinateTerzoPuntoGriglia, *FOG);
        fieldOnGrid.push_back(*FOG);
        calcolaCampoAnalitico1DSuParticella(data, coordinateQuartoPuntoGriglia, *FOG);
        fieldOnGrid.push_back(*FOG);
        calcolaCampoAnalitico1DSuParticella(data, coordinateQuintoPuntoGriglia, *FOG);
        fieldOnGrid.push_back(*FOG);
        calcolaCampoAnalitico1DSuParticella(data, coordinateSestoPuntoGriglia, *FOG);
        fieldOnGrid.push_back(*FOG);
        calcolaCampoAnalitico1DSuParticella(data, coordinateSettimoPuntoGriglia, *FOG);
        fieldOnGrid.push_back(*FOG);
        calcolaCampoAnalitico1DSuParticella(data, coordinateOttavoPuntoGriglia, *FOG);
        fieldOnGrid.push_back(*FOG);

        FOP = interpolation2D_linear(data, particles[j], fieldOnGrid);
        fieldOnParticles.push_back(FOP);

        fieldOnGrid.clear();

        //pstar = p(n-1/2) + Q * E(n) * DeltaT/2
        pstar.setParticlePX(particles[j].px + 0.5 * data.dt * particles[j].q * fieldOnParticles[0].ex);
        pstar.setParticlePY(particles[j].py + 0.5 * data.dt * particles[j].q * fieldOnParticles[0].ey);
        pstar.setParticlePZ(particles[j].pz + 0.5 * data.dt * particles[j].q * fieldOnParticles[0].ez);
        gamma = sqrt(1. + pow(pstar.px, 2) + pow(pstar.py, 2) + pow(pstar.pz, 2));
        gamma_inv = 1. / gamma;

        //b = (DeltaT/2) * Q * B(n) / gamma
        b.setBx(0.5 * data.dt * fieldOnParticles[0].bx * particles[j].q * gamma_inv);
        b.setBy(0.5 * data.dt * fieldOnParticles[0].by * particles[j].q * gamma_inv);
        b.setBz(0.5 * data.dt * fieldOnParticles[0].bz * particles[j].q * gamma_inv);
        b2 = pow(b.bx, 2) + pow(b.by, 2) + pow(b.bz, 2);
        b2 += 1.;

        /****************************************************************************
        //p(n+1/2) = [pstar + pstar x b + (pstar x b) * b] / b2   SECONDO TURCHETTI  (ma i conti analitici diconon che e` sbagliato, anche perche' l'ultimo termine e` uno scalare!)
        pHalfAdvanced.setParticlePX( (pstar.px + pstar.py * b.bz - b.by * pstar.pz + (pstar.py * b.bz - b.by * pstar.pz) * b.bx ) / b2 );
        pHalfAdvanced.setParticlePY( (pstar.py + pstar.pz * b.bx - b.bz * pstar.px + (pstar.pz * b.bx - b.bz * pstar.px) * b.by ) / b2 );
        pHalfAdvanced.setParticlePZ( (pstar.pz + pstar.px * b.by - b.bx * pstar.py + (pstar.px * b.by - b.bx * pstar.py) * b.bz ) / b2 );
        ****************************************************************************/

        /****************************************************************************
        //p(n+1/2) = [pstar + pstar x b + b x (pstar x b)] / b2   SECONDO LONDRILLO
        pHalfAdvanced.setParticlePX( (pstar.px + pstar.py * b.bz - b.by * pstar.pz + b.bx * pstar.px * b.bx ) / b2 );
        pHalfAdvanced.setParticlePY( (pstar.py + pstar.pz * b.bx - b.bz * pstar.px + b.by * pstar.py * b.by ) / b2 );
        pHalfAdvanced.setParticlePZ( (pstar.pz + pstar.px * b.by - b.bx * pstar.py + b.bz * pstar.pz * b.bz ) / b2 );
        ****************************************************************************/

        // /****************************************************************************
        //p(n+1/2) = [pstar + pstar x b] / b2   COME IMPLEMENTATO DA TURCHETTI NEL CODICE
        pHalfAdvanced.setParticlePX((pstar.px + pstar.py * b.bz) / b2);
        pHalfAdvanced.setParticlePY((pstar.py - pstar.px * b.bz) / b2);
        pHalfAdvanced.setParticlePZ(0.); //per come ho strutturato i pacchetti d'onda su z non dovrebbe succedere nulla e quindi lo fisso io a zero.
        // ****************************************************************************/

        particles[j].setParticlePX(2. * pHalfAdvanced.px - particles[j].px);
        particles[j].setParticlePY(2. * pHalfAdvanced.py - particles[j].py);
        particles[j].setParticlePZ(2. * pHalfAdvanced.pz - particles[j].pz);

        gamma = sqrt(1. + pow(particles[j].px, 2) + pow(particles[j].py, 2) + pow(particles[j].pz, 2));
        gamma_inv = 1. / gamma;

        particles[j].setParticleX(particles[j].x + 0.5 * data.dt * particles[j].px * gamma_inv);
        particles[j].setParticleY(particles[j].y + 0.5 * data.dt * particles[j].py * gamma_inv);
        particles[j].setParticleZ(particles[j].z + 0.5 * data.dt * particles[j].pz * gamma_inv);

        //INTEGRALI PRIMI (da sistemare, li ho copiati dall'evoluzione sui campi analitici)
        potenziale.clear();
        zeta = particles[j].x - SPEED_OF_LIGHT * data.t; // non usato nella mia struttura delle equazioni di campo

        calcolaPotenzialeAnalitico1DSuParticella(data, particles[j], *campoSuPunto);
        potenziale.push_back(*campoSuPunto);

        I1 = (particles[j].py / (particles[j].m * SPEED_OF_LIGHT)) + potenziale[0].ey;
        I2 = particles[j].pz / (particles[j].m * SPEED_OF_LIGHT);
        I3 = gamma - particles[j].px / (particles[j].m * SPEED_OF_LIGHT);

        outputstream
            << setw(12) << i + 1 << "\t"
            << setprecision(14) << setw(16) << particles[j].x << "\t"
            << setw(16) << particles[j].y << "\t"
            << setw(16) << particles[j].z << "\t"
            // << setw(16) << particles[j].getParticleCell() << "\t"
            << setw(16) << particles[j].px << "\t"
            << setw(16) << particles[j].py << "\t"
            << setw(16) << particles[j].pz << "\t"
            << setw(16) << fieldOnParticles[0].ey << "\t"
            << setw(16) << potenziale[0].ey << "\t"
            << setw(16) << zeta << "\t"
            << setw(16) << I1 << "\t"
            << setw(16) << I2 << "\t"
            << setw(16) << I3 << endl;
      }
    }
  }
}

void evolveLPF_withgrid(Data data, vector<Particle>& particles, vector<Field>& fieldOnParticles, vector<Field>& campoSuPuntiGriglia, ofstream& outputstream)
{
  vector<Field> potenzialeSuPuntiGriglia;

  vector<Field> potenziale;
  vector<Field> fieldOnGrid;

  Particle pstar, pHalfAdvanced;

  int posX = 0, posY = 0, posZ = 0;
  Field b;
  Field FOP;

  double gridXX, gridYY, gridZZ;
  int gridX, gridY, gridZ;
  double gamma, gamma_inv, b2;
  double I1, I2, I3;
  double zeta;

  Field* campoSuPunto;
  campoSuPunto = new Field[1];

  //ALGORITMO LEAPFROG PARTICELLE

  if (data.n_dim == 1) {
    for (int i = 0; i < data.nSteps; i++) {
      data.aumentaT(data.dt);

      campoSuPuntiGriglia.clear();
      riempiPuntiGrigliaConCampoAnalitico1D(data, campoSuPuntiGriglia);

      potenzialeSuPuntiGriglia.clear();
      riempiPuntiGrigliaConPotenzialeAnalitico1D(data, potenzialeSuPuntiGriglia);

      for (int j = 0; j < data.n_electrons; j++) {
        gamma = sqrt(1. + pow(particles[j].px, 2) + pow(particles[j].py, 2) + pow(particles[j].pz, 2));
        gamma_inv = 1. / gamma;

        particles[j].setParticleX(particles[j].x + 0.5 * data.dt * particles[j].px * gamma_inv);
        particles[j].setParticleY(particles[j].y + 0.5 * data.dt * particles[j].py * gamma_inv);
        particles[j].setParticleZ(particles[j].z + 0.5 * data.dt * particles[j].pz * gamma_inv);

        // QUI MI CALCOLO I PUNTI DI RIFERIMENTO PER LA PARTICELLA IN ESAME nella posizione attuale
        // forse c'e` un modo piu` furbo per farlo sfruttando il # di cella in cui si trova che e` noto

        fieldOnParticles.clear(); //SVUOTA il vettore dei campi sulle particelle ad ogni ciclo

        gridXX = particles[j].x / data.deltaX;
        if (gridXX < 0)
          gridX = ((int)gridXX) - 1; //il punto di griglia lungo l'asse x a sinistra della particella
        else
          gridX = (int)gridXX;

        gridYY = particles[j].y / data.deltaY;
        if (gridYY < 0)
          gridY = ((int)gridYY) - 1; //il punto di griglia lungo l'asse y a sinistra della particella
        else
          gridY = (int)gridYY;

        gridZZ = particles[j].z / data.deltaZ;
        if (gridZZ < 0)
          gridZ = ((int)gridZZ) - 1; //il punto di griglia lungo l'asse z a sinistra della particella
        else
          gridZ = (int)gridZZ;

        if (gridX < 0 || gridY < 0 || gridZ < 0) {
          cout << "La particella e` uscita dai boundaries - simulazione interrotta" << endl;
          return;
        }

        fieldOnGrid.push_back(campoSuPuntiGriglia.at(gridY * data.nGridPointsX + gridX)); // primo punto
        fieldOnGrid.push_back(campoSuPuntiGriglia.at(gridY * data.nGridPointsX + (gridX + 1))); // secondo punto
        fieldOnGrid.push_back(campoSuPuntiGriglia.at((gridY + 1) * data.nGridPointsX + gridX)); // terzo punto
        fieldOnGrid.push_back(campoSuPuntiGriglia.at((gridY + 1) * data.nGridPointsX + (gridX + 1))); // quarto punto

        FOP = interpolation2D_linear(data, particles[j], fieldOnGrid);
        fieldOnParticles.push_back(FOP);

        fieldOnGrid.clear();

        //pstar = p(n-1/2) + Q * E(n) * DeltaT/2
        pstar.setParticlePX(particles[j].px + 0.5 * data.dt * particles[j].q * fieldOnParticles[0].ex);
        pstar.setParticlePY(particles[j].py + 0.5 * data.dt * particles[j].q * fieldOnParticles[0].ey);
        pstar.setParticlePZ(particles[j].pz + 0.5 * data.dt * particles[j].q * fieldOnParticles[0].ez);
        gamma = sqrt(1. + pow(pstar.px, 2) + pow(pstar.py, 2) + pow(pstar.pz, 2));
        gamma_inv = 1. / gamma;

        //b = (DeltaT/2) * Q * B(n) / gamma
        b.setBx(0.5 * data.dt * fieldOnParticles[0].bx * particles[j].q * gamma_inv);
        b.setBy(0.5 * data.dt * fieldOnParticles[0].by * particles[j].q * gamma_inv);
        b.setBz(0.5 * data.dt * fieldOnParticles[0].bz * particles[j].q * gamma_inv);
        b2 = pow(b.bx, 2) + pow(b.by, 2) + pow(b.bz, 2);
        b2 += 1.;

        /****************************************************************************
        //p(n+1/2) = [pstar + pstar x b + (pstar x b) * b] / b2   SECONDO TURCHETTI  (ma i conti analitici diconon che e` sbagliato, anche perche' l'ultimo termine e` uno scalare!)
        pHalfAdvanced.setParticlePX( (pstar.px + pstar.py * b.bz - b.by * pstar.pz + (pstar.py * b.bz - b.by * pstar.pz) * b.bx ) / b2 );
        pHalfAdvanced.setParticlePY( (pstar.py + pstar.pz * b.bx - b.bz * pstar.px + (pstar.pz * b.bx - b.bz * pstar.px) * b.by ) / b2 );
        pHalfAdvanced.setParticlePZ( (pstar.pz + pstar.px * b.by - b.bx * pstar.py + (pstar.px * b.by - b.bx * pstar.py) * b.bz ) / b2 );
        ****************************************************************************/

        /****************************************************************************
        //p(n+1/2) = [pstar + pstar x b + b x (pstar x b)] / b2   SECONDO LONDRILLO
        pHalfAdvanced.setParticlePX( (pstar.px + pstar.py * b.bz - b.by * pstar.pz + b.bx * pstar.px * b.bx ) / b2 );
        pHalfAdvanced.setParticlePY( (pstar.py + pstar.pz * b.bx - b.bz * pstar.px + b.by * pstar.py * b.by ) / b2 );
        pHalfAdvanced.setParticlePZ( (pstar.pz + pstar.px * b.by - b.bx * pstar.py + b.bz * pstar.pz * b.bz ) / b2 );
        ****************************************************************************/

        //      /****************************************************************************
        //p(n+1/2) = [pstar + pstar x b] / b2   COME IMPLEMENTATO DA TURCHETTI NEL CODICE
        pHalfAdvanced.setParticlePX((pstar.px + pstar.py * b.bz) / b2);
        pHalfAdvanced.setParticlePY((pstar.py - pstar.px * b.bz) / b2);
        pHalfAdvanced.setParticlePZ(0.); //per come ho strutturato i pacchetti d'onda su z non dovrebbe succedere nulla e quindi lo fisso io a zero.
        //      ****************************************************************************/

        particles[j].setParticlePX(2. * pHalfAdvanced.px - particles[j].px);
        particles[j].setParticlePY(2. * pHalfAdvanced.py - particles[j].py);
        particles[j].setParticlePZ(2. * pHalfAdvanced.pz - particles[j].pz);

        gamma = sqrt(1. + pow(particles[j].px, 2) + pow(particles[j].py, 2) + pow(particles[j].pz, 2));
        gamma_inv = 1. / gamma;

        particles[j].setParticleX(particles[j].x + 0.5 * data.dt * particles[j].px * gamma_inv);
        particles[j].setParticleY(particles[j].y + 0.5 * data.dt * particles[j].py * gamma_inv);
        particles[j].setParticleZ(particles[j].z + 0.5 * data.dt * particles[j].pz * gamma_inv);

        //INTEGRALI PRIMI (da sistemare, li ho copiati dall'evoluzione sui campi analitici)
        potenziale.clear();
        zeta = particles[j].x - SPEED_OF_LIGHT * data.t; // non usato nella mia struttura delle equazioni di campo

        gridXX = particles[j].x / data.deltaX;
        if (gridXX < 0)
          gridX = ((int)gridXX) - 1; //il punto di griglia lungo l'asse x a sinistra della particella
        else
          gridX = (int)gridXX;

        gridYY = particles[j].y / data.deltaY;
        if (gridYY < 0)
          gridY = ((int)gridYY) - 1; //il punto di griglia lungo l'asse y a sinistra della particella
        else
          gridY = (int)gridYY;

        gridZZ = particles[j].z / data.deltaZ;
        if (gridZZ < 0)
          gridZ = ((int)gridZZ) - 1; //il punto di griglia lungo l'asse z a sinistra della particella
        else
          gridZ = (int)gridZZ;

        if (gridX < 0 || gridY < 0 || gridZ < 0) {
          cout << "La particella e` uscita dai boundaries - simulazione interrotta" << endl;
          return;
        }

        fieldOnGrid.clear();

        fieldOnGrid.push_back(potenzialeSuPuntiGriglia.at(gridY * data.nGridPointsX + gridX)); // primo punto
        fieldOnGrid.push_back(potenzialeSuPuntiGriglia.at(gridY * data.nGridPointsX + (gridX + 1))); // secondo punto
        fieldOnGrid.push_back(potenzialeSuPuntiGriglia.at((gridY + 1) * data.nGridPointsX + gridX)); // terzo punto
        fieldOnGrid.push_back(potenzialeSuPuntiGriglia.at((gridY + 1) * data.nGridPointsX + (gridX + 1))); // quarto punto

        FOP = interpolation2D_linear(data, particles[j], fieldOnGrid);
        potenziale.push_back(FOP);

        I1 = (particles[j].py / (particles[j].m * SPEED_OF_LIGHT)) + potenziale[0].ey;
        I2 = particles[j].pz / (particles[j].m * SPEED_OF_LIGHT);
        I3 = gamma - particles[j].px / (particles[j].m * SPEED_OF_LIGHT);

        outputstream
            << setw(12) << i + 1 << "\t"
            << setprecision(14) << setw(16) << particles[j].x << "\t"
            << setw(16) << particles[j].y << "\t"
            << setw(16) << particles[j].z << "\t"
            //        << setw(16) << particles[j].getParticleCell() << "\t"
            << setw(16) << particles[j].px << "\t"
            << setw(16) << particles[j].py << "\t"
            << setw(16) << particles[j].pz << "\t"
            << setw(16) << potenziale[0].ey << "\t"
            << setw(16) << zeta << "\t"
            << setw(16) << I1 << "\t"
            << setw(16) << I2 << "\t"
            << setw(16) << I3 << endl;

        /*
        //UPDATE CAMPI ATTORNO ALLA NUOVA POSIZIONE DELLA PARTICELLA  (al momento non li evolve, il set li pone pari al get precedente)
        gridX = (int) (particles[j].x/data.deltaX);   //il punto di griglia lungo l'asse x piu` vicino all'origine
        gridY = (int) (particles[j].y/data.deltaY);   //il punto di griglia lungo l'asse y piu` vicino all'origine
        gridZ = (int) (particles[j].z/data.deltaZ);   //il punto di griglia lungo l'asse z piu` vicino all'origine

        campoSuPuntiGriglia.at(gridZ * data.nGridPointsX * data.nGridPointsY + gridY * data.nGridPointsX + gridX).setEx( campoSuPuntiGriglia.at(gridZ * data.nGridPointsX * data.nGridPointsY + gridY * data.nGridPointsX + gridX).ex );       // primo punto
        campoSuPuntiGriglia.at(gridZ * data.nGridPointsX * data.nGridPointsY + gridY * data.nGridPointsX + gridX).setEy( campoSuPuntiGriglia.at(gridZ * data.nGridPointsX * data.nGridPointsY + gridY * data.nGridPointsX + gridX).ey );
        campoSuPuntiGriglia.at(gridZ * data.nGridPointsX * data.nGridPointsY + gridY * data.nGridPointsX + gridX).setEz( campoSuPuntiGriglia.at(gridZ * data.nGridPointsX * data.nGridPointsY + gridY * data.nGridPointsX + gridX).ez );
        campoSuPuntiGriglia.at(gridZ * data.nGridPointsX * data.nGridPointsY + gridY * data.nGridPointsX + gridX).setBx( campoSuPuntiGriglia.at(gridZ * data.nGridPointsX * data.nGridPointsY + gridY * data.nGridPointsX + gridX).bx );
        campoSuPuntiGriglia.at(gridZ * data.nGridPointsX * data.nGridPointsY + gridY * data.nGridPointsX + gridX).setBy( campoSuPuntiGriglia.at(gridZ * data.nGridPointsX * data.nGridPointsY + gridY * data.nGridPointsX + gridX).by );
        campoSuPuntiGriglia.at(gridZ * data.nGridPointsX * data.nGridPointsY + gridY * data.nGridPointsX + gridX).setBz( campoSuPuntiGriglia.at(gridZ * data.nGridPointsX * data.nGridPointsY + gridY * data.nGridPointsX + gridX).bz );

        campoSuPuntiGriglia.at(gridZ * data.nGridPointsX * data.nGridPointsY + gridY * data.nGridPointsX + (gridX+1)).setEx( campoSuPuntiGriglia.at(gridZ * data.nGridPointsX * data.nGridPointsY + gridY * data.nGridPointsX + (gridX+1)).ex );   // secondo punto
        campoSuPuntiGriglia.at(gridZ * data.nGridPointsX * data.nGridPointsY + gridY * data.nGridPointsX + (gridX+1)).setEy( campoSuPuntiGriglia.at(gridZ * data.nGridPointsX * data.nGridPointsY + gridY * data.nGridPointsX + (gridX+1)).ey );
        campoSuPuntiGriglia.at(gridZ * data.nGridPointsX * data.nGridPointsY + gridY * data.nGridPointsX + (gridX+1)).setEz( campoSuPuntiGriglia.at(gridZ * data.nGridPointsX * data.nGridPointsY + gridY * data.nGridPointsX + (gridX+1)).ez );
        campoSuPuntiGriglia.at(gridZ * data.nGridPointsX * data.nGridPointsY + gridY * data.nGridPointsX + (gridX+1)).setBx( campoSuPuntiGriglia.at(gridZ * data.nGridPointsX * data.nGridPointsY + gridY * data.nGridPointsX + (gridX+1)).bx );
        campoSuPuntiGriglia.at(gridZ * data.nGridPointsX * data.nGridPointsY + gridY * data.nGridPointsX + (gridX+1)).setBy( campoSuPuntiGriglia.at(gridZ * data.nGridPointsX * data.nGridPointsY + gridY * data.nGridPointsX + (gridX+1)).by );
        campoSuPuntiGriglia.at(gridZ * data.nGridPointsX * data.nGridPointsY + gridY * data.nGridPointsX + (gridX+1)).setBz( campoSuPuntiGriglia.at(gridZ * data.nGridPointsX * data.nGridPointsY + gridY * data.nGridPointsX + (gridX+1)).bz );

        campoSuPuntiGriglia.at(gridZ * data.nGridPointsX * data.nGridPointsY + (gridY+1) * data.nGridPointsX + gridX).setEx( campoSuPuntiGriglia.at(gridZ * data.nGridPointsX * data.nGridPointsY + (gridY+1) * data.nGridPointsX + gridX).ex );   // terzo punto
        campoSuPuntiGriglia.at(gridZ * data.nGridPointsX * data.nGridPointsY + (gridY+1) * data.nGridPointsX + gridX).setEy( campoSuPuntiGriglia.at(gridZ * data.nGridPointsX * data.nGridPointsY + (gridY+1) * data.nGridPointsX + gridX).ey );
        campoSuPuntiGriglia.at(gridZ * data.nGridPointsX * data.nGridPointsY + (gridY+1) * data.nGridPointsX + gridX).setEz( campoSuPuntiGriglia.at(gridZ * data.nGridPointsX * data.nGridPointsY + (gridY+1) * data.nGridPointsX + gridX).ez );
        campoSuPuntiGriglia.at(gridZ * data.nGridPointsX * data.nGridPointsY + (gridY+1) * data.nGridPointsX + gridX).setBx( campoSuPuntiGriglia.at(gridZ * data.nGridPointsX * data.nGridPointsY + (gridY+1) * data.nGridPointsX + gridX).bx );
        campoSuPuntiGriglia.at(gridZ * data.nGridPointsX * data.nGridPointsY + (gridY+1) * data.nGridPointsX + gridX).setBy( campoSuPuntiGriglia.at(gridZ * data.nGridPointsX * data.nGridPointsY + (gridY+1) * data.nGridPointsX + gridX).by );
        campoSuPuntiGriglia.at(gridZ * data.nGridPointsX * data.nGridPointsY + (gridY+1) * data.nGridPointsX + gridX).setBz( campoSuPuntiGriglia.at(gridZ * data.nGridPointsX * data.nGridPointsY + (gridY+1) * data.nGridPointsX + gridX).bz );

        campoSuPuntiGriglia.at(gridZ * data.nGridPointsX * data.nGridPointsY + (gridY+1) * data.nGridPointsX + (gridX+1)).setEx( campoSuPuntiGriglia.at(gridZ * data.nGridPointsX * data.nGridPointsY + (gridY+1) * data.nGridPointsX + (gridX+1)).ex );   // quarto punto
        campoSuPuntiGriglia.at(gridZ * data.nGridPointsX * data.nGridPointsY + (gridY+1) * data.nGridPointsX + (gridX+1)).setEy( campoSuPuntiGriglia.at(gridZ * data.nGridPointsX * data.nGridPointsY + (gridY+1) * data.nGridPointsX + (gridX+1)).ey );
        campoSuPuntiGriglia.at(gridZ * data.nGridPointsX * data.nGridPointsY + (gridY+1) * data.nGridPointsX + (gridX+1)).setEz( campoSuPuntiGriglia.at(gridZ * data.nGridPointsX * data.nGridPointsY + (gridY+1) * data.nGridPointsX + (gridX+1)).ez );
        campoSuPuntiGriglia.at(gridZ * data.nGridPointsX * data.nGridPointsY + (gridY+1) * data.nGridPointsX + (gridX+1)).setBx( campoSuPuntiGriglia.at(gridZ * data.nGridPointsX * data.nGridPointsY + (gridY+1) * data.nGridPointsX + (gridX+1)).bx );
        campoSuPuntiGriglia.at(gridZ * data.nGridPointsX * data.nGridPointsY + (gridY+1) * data.nGridPointsX + (gridX+1)).setBy( campoSuPuntiGriglia.at(gridZ * data.nGridPointsX * data.nGridPointsY + (gridY+1) * data.nGridPointsX + (gridX+1)).by );
        campoSuPuntiGriglia.at(gridZ * data.nGridPointsX * data.nGridPointsY + (gridY+1) * data.nGridPointsX + (gridX+1)).setBz( campoSuPuntiGriglia.at(gridZ * data.nGridPointsX * data.nGridPointsY + (gridY+1) * data.nGridPointsX + (gridX+1)).bz );

        campoSuPuntiGriglia.at((gridZ+1) * data.nGridPointsX * data.nGridPointsY + gridY * data.nGridPointsX + gridX).setEx( campoSuPuntiGriglia.at((gridZ+1) * data.nGridPointsX * data.nGridPointsY + gridY * data.nGridPointsX + gridX).ex );   // quinto punto
        campoSuPuntiGriglia.at((gridZ+1) * data.nGridPointsX * data.nGridPointsY + gridY * data.nGridPointsX + gridX).setEy( campoSuPuntiGriglia.at((gridZ+1) * data.nGridPointsX * data.nGridPointsY + gridY * data.nGridPointsX + gridX).ey );
        campoSuPuntiGriglia.at((gridZ+1) * data.nGridPointsX * data.nGridPointsY + gridY * data.nGridPointsX + gridX).setEz( campoSuPuntiGriglia.at((gridZ+1) * data.nGridPointsX * data.nGridPointsY + gridY * data.nGridPointsX + gridX).ez );
        campoSuPuntiGriglia.at((gridZ+1) * data.nGridPointsX * data.nGridPointsY + gridY * data.nGridPointsX + gridX).setBx( campoSuPuntiGriglia.at((gridZ+1) * data.nGridPointsX * data.nGridPointsY + gridY * data.nGridPointsX + gridX).bx );
        campoSuPuntiGriglia.at((gridZ+1) * data.nGridPointsX * data.nGridPointsY + gridY * data.nGridPointsX + gridX).setBy( campoSuPuntiGriglia.at((gridZ+1) * data.nGridPointsX * data.nGridPointsY + gridY * data.nGridPointsX + gridX).by );
        campoSuPuntiGriglia.at((gridZ+1) * data.nGridPointsX * data.nGridPointsY + gridY * data.nGridPointsX + gridX).setBz( campoSuPuntiGriglia.at((gridZ+1) * data.nGridPointsX * data.nGridPointsY + gridY * data.nGridPointsX + gridX).bz );

        campoSuPuntiGriglia.at((gridZ+1) * data.nGridPointsX * data.nGridPointsY + gridY * data.nGridPointsX + (gridX+1)).setEx( campoSuPuntiGriglia.at((gridZ+1) * data.nGridPointsX * data.nGridPointsY + gridY * data.nGridPointsX + (gridX+1)).ex );   // sesto punto
        campoSuPuntiGriglia.at((gridZ+1) * data.nGridPointsX * data.nGridPointsY + gridY * data.nGridPointsX + (gridX+1)).setEy( campoSuPuntiGriglia.at((gridZ+1) * data.nGridPointsX * data.nGridPointsY + gridY * data.nGridPointsX + (gridX+1)).ey );
        campoSuPuntiGriglia.at((gridZ+1) * data.nGridPointsX * data.nGridPointsY + gridY * data.nGridPointsX + (gridX+1)).setEz( campoSuPuntiGriglia.at((gridZ+1) * data.nGridPointsX * data.nGridPointsY + gridY * data.nGridPointsX + (gridX+1)).ez );
        campoSuPuntiGriglia.at((gridZ+1) * data.nGridPointsX * data.nGridPointsY + gridY * data.nGridPointsX + (gridX+1)).setBx( campoSuPuntiGriglia.at((gridZ+1) * data.nGridPointsX * data.nGridPointsY + gridY * data.nGridPointsX + (gridX+1)).bx );
        campoSuPuntiGriglia.at((gridZ+1) * data.nGridPointsX * data.nGridPointsY + gridY * data.nGridPointsX + (gridX+1)).setBy( campoSuPuntiGriglia.at((gridZ+1) * data.nGridPointsX * data.nGridPointsY + gridY * data.nGridPointsX + (gridX+1)).by );
        campoSuPuntiGriglia.at((gridZ+1) * data.nGridPointsX * data.nGridPointsY + gridY * data.nGridPointsX + (gridX+1)).setBz( campoSuPuntiGriglia.at((gridZ+1) * data.nGridPointsX * data.nGridPointsY + gridY * data.nGridPointsX + (gridX+1)).bz );

        campoSuPuntiGriglia.at((gridZ+1) * data.nGridPointsX * data.nGridPointsY + (gridY+1) * data.nGridPointsX + gridX).setEx( campoSuPuntiGriglia.at((gridZ+1) * data.nGridPointsX * data.nGridPointsY + (gridY+1) * data.nGridPointsX + gridX).ex );   // settimo punto
        campoSuPuntiGriglia.at((gridZ+1) * data.nGridPointsX * data.nGridPointsY + (gridY+1) * data.nGridPointsX + gridX).setEy( campoSuPuntiGriglia.at((gridZ+1) * data.nGridPointsX * data.nGridPointsY + (gridY+1) * data.nGridPointsX + gridX).ey );
        campoSuPuntiGriglia.at((gridZ+1) * data.nGridPointsX * data.nGridPointsY + (gridY+1) * data.nGridPointsX + gridX).setEz( campoSuPuntiGriglia.at((gridZ+1) * data.nGridPointsX * data.nGridPointsY + (gridY+1) * data.nGridPointsX + gridX).ez );
        campoSuPuntiGriglia.at((gridZ+1) * data.nGridPointsX * data.nGridPointsY + (gridY+1) * data.nGridPointsX + gridX).setBx( campoSuPuntiGriglia.at((gridZ+1) * data.nGridPointsX * data.nGridPointsY + (gridY+1) * data.nGridPointsX + gridX).bx );
        campoSuPuntiGriglia.at((gridZ+1) * data.nGridPointsX * data.nGridPointsY + (gridY+1) * data.nGridPointsX + gridX).setBy( campoSuPuntiGriglia.at((gridZ+1) * data.nGridPointsX * data.nGridPointsY + (gridY+1) * data.nGridPointsX + gridX).by );
        campoSuPuntiGriglia.at((gridZ+1) * data.nGridPointsX * data.nGridPointsY + (gridY+1) * data.nGridPointsX + gridX).setBz( campoSuPuntiGriglia.at((gridZ+1) * data.nGridPointsX * data.nGridPointsY + (gridY+1) * data.nGridPointsX + gridX).bz );

        campoSuPuntiGriglia.at((gridZ+1) * data.nGridPointsX * data.nGridPointsY + (gridY+1) * data.nGridPointsX + (gridX+1)).setEx( campoSuPuntiGriglia.at((gridZ+1) * data.nGridPointsX * data.nGridPointsY + (gridY+1) * data.nGridPointsX + (gridX+1)).ex );   // ottavo punto
        campoSuPuntiGriglia.at((gridZ+1) * data.nGridPointsX * data.nGridPointsY + (gridY+1) * data.nGridPointsX + (gridX+1)).setEy( campoSuPuntiGriglia.at((gridZ+1) * data.nGridPointsX * data.nGridPointsY + (gridY+1) * data.nGridPointsX + (gridX+1)).ex );
        campoSuPuntiGriglia.at((gridZ+1) * data.nGridPointsX * data.nGridPointsY + (gridY+1) * data.nGridPointsX + (gridX+1)).setEz( campoSuPuntiGriglia.at((gridZ+1) * data.nGridPointsX * data.nGridPointsY + (gridY+1) * data.nGridPointsX + (gridX+1)).ex );
        campoSuPuntiGriglia.at((gridZ+1) * data.nGridPointsX * data.nGridPointsY + (gridY+1) * data.nGridPointsX + (gridX+1)).setBx( campoSuPuntiGriglia.at((gridZ+1) * data.nGridPointsX * data.nGridPointsY + (gridY+1) * data.nGridPointsX + (gridX+1)).ex );
        campoSuPuntiGriglia.at((gridZ+1) * data.nGridPointsX * data.nGridPointsY + (gridY+1) * data.nGridPointsX + (gridX+1)).setBy( campoSuPuntiGriglia.at((gridZ+1) * data.nGridPointsX * data.nGridPointsY + (gridY+1) * data.nGridPointsX + (gridX+1)).ex );
        campoSuPuntiGriglia.at((gridZ+1) * data.nGridPointsX * data.nGridPointsY + (gridY+1) * data.nGridPointsX + (gridX+1)).setBz( campoSuPuntiGriglia.at((gridZ+1) * data.nGridPointsX * data.nGridPointsY + (gridY+1) * data.nGridPointsX + (gridX+1)).ex );
        */
      }
    }
  }
}

void calcolaCampoInUnIstante(Data data, vector<Particle>& particelle, vector<Field>& campiSuParticelle, ofstream& outputDATA)
{
  /*
  In questa funzione si genera un vettore contenente n particelle, n da definire prima di chiamare la funzione
  invocando data.setNelectrons() se si vuole avere qualcosa indietro, in posizione random. Su ciascuna viene calcolato
  il valore del campo analitico che deve essere gia` stato definito prima di chiamare la funzione invocando
  data.setK() e data.setA()  ( E= - k A sin[k(z-ct)] ) in un dato istante gia` definito tramite data.impostaT(0)
  */

  calcolaCampoAnalitico1DSuParticelle(data, particelle, campiSuParticelle);
  for (int i = 0; i < data.n_electrons; i++) {
    outputDATA
        << setw(12) << i + 1 << "\t"
        << setprecision(14) << setw(16) << particelle[i].x << "\t"
        << setw(16) << particelle[i].y << "\t"
        << setw(16) << particelle[i].z << "\t"
        << setw(16) << data.t << "\t"
        << setw(16) << campiSuParticelle[i].ex << "\t"
        << setw(16) << campiSuParticelle[i].ey << "\t"
        << setw(16) << campiSuParticelle[i].ez << "\t"
        << setw(16) << campiSuParticelle[i].bx << "\t"
        << setw(16) << campiSuParticelle[i].by << "\t"
        << setw(16) << campiSuParticelle[i].bz << endl;
  }
}

void calcolaEvoluzioneCampoSuParticelle(Data data, vector<Particle>& particelle, vector<Field>& campiSuParticelle, ofstream& outputDATA)
{
  /*
  In questa funzione si genera un vettore contenente n particelle, n da definire prima di chiamare la funzione
  invocando data.setNelectrons() se si vuole avere qualcosa indietro, in posizione random se sono 3+, manuale se sono di meno.
  Su ciascuna viene calcolato il valore del campo analitico, che deve essere gia` stato definito prima di chiamare la funzione invocando
  data.setK() e data.setA()  ( E = -k A sin[k(z-ct)] ) in un dato istante gia` definito tramite data.setT().
  L'operazione viene ripetuta per nsteps distanti ciascuno un deltaT, impostati rispettivamente con data.setNsteps(int nStepsMin)
  e data.setDeltaT(double deltaTmin)
  */

  data.setDeltaT(0.);
  data.setNsteps(1);

  Field* campoTemp;
  campoTemp = new Field[1];
  campiSuParticelle.clear();

  for (int i = 0; i < data.nSteps; i++) {
    for (int j = 0; j < data.n_electrons; j++) {
      calcolaCampoAnalitico1DSuParticella(data, particelle[j], *campoTemp);
      campiSuParticelle.push_back(*campoTemp);
      outputDATA
          << setw(12) << j + 1 << "\t"
          << setprecision(14) << setw(16) << particelle[j].x << "\t"
          << setw(16) << particelle[j].y << "\t"
          << setw(16) << particelle[j].z << "\t"
          << setw(16) << data.t << "\t"
          << setw(16) << campiSuParticelle[j].ex << "\t"
          << setw(16) << campiSuParticelle[j].ey << "\t"
          << setw(16) << campiSuParticelle[j].ez << "\t"
          << setw(16) << campiSuParticelle[j].bx << "\t"
          << setw(16) << campiSuParticelle[j].by << "\t"
          << setw(16) << campiSuParticelle[j].bz << endl;
    }
    data.aumentaT(data.dt);
    campiSuParticelle.clear();
  }
}

void evolveLPF4_nogrid_TURCHETTI(Data data, vector<Particle>& particles, vector<Field>& fieldOnParticles, ofstream& outputstream)
{
  Particle evolvedParticle;
  Particle pstar, pHalfAdvanced;

  Field b;
  Field* FOG;
  FOG = new Field[1];
  vector<Field> potenziale;
  Field* campoSuPunto;
  campoSuPunto = new Field[1];

  double gamma, gammastar, coeff;
  double I1, I2, I3;
  double zeta;

  double alpha = 1. / (2. - pow(2., (1. / 3.)));
  double dtx1 = alpha * data.dt * 0.5;
  double dtp1 = alpha * data.dt;
  double dtx2 = (1. - alpha) * 0.5 * data.dt;
  double dtp2 = (1. - 2. * alpha) * data.dt;
  double dtx3 = (1. - alpha) * 0.5 * data.dt;
  double dtp3 = alpha * data.dt;

  for (int i = 0; i < data.nSteps; i++) {
    for (int j = 0; j < data.n_electrons; j++) {
      fieldOnParticles.clear();

      // PASSO #1
      gamma = sqrt(1. + pow(particles[j].px, 2) + pow(particles[j].py, 2));
      particles[j].setParticleX(particles[j].x + dtx1 * particles[j].px / (gamma));
      particles[j].setParticleY(particles[j].y + dtx1 * particles[j].py / (gamma));

      calcolaCampoAnalitico1DSuParticella(data, particles[j], *FOG);
      fieldOnParticles.push_back(*FOG);

      pstar.setParticlePX(particles[j].px);
      pstar.setParticlePY(particles[j].py + dtp1 * 0.5 * fieldOnParticles[0].ey);

      gammastar = sqrt(1. + pow(pstar.px, 2) + pow(pstar.py, 2));
      b.setEy(0.5 * dtp1 * fieldOnParticles[0].ey / gammastar);
      coeff = 1. / (1. + pow(b.ey, 2));

      pHalfAdvanced.setParticlePX(coeff * (pstar.px + pstar.py * b.bz));
      pHalfAdvanced.setParticlePY(coeff * (pstar.py - pstar.px * b.bz));

      particles[j].setParticlePX(2. * pHalfAdvanced.px - particles[j].px);
      particles[j].setParticlePY(2. * pHalfAdvanced.py - particles[j].py);

      data.aumentaT(dtp1);

      // PASSO #2
      gamma = sqrt(1. + pow(particles[j].px, 2) + pow(particles[j].py, 2));
      particles[j].setParticleX(particles[j].x + dtx2 * particles[j].px / (gamma));
      particles[j].setParticleY(particles[j].y + dtx2 * particles[j].py / (gamma));

      calcolaCampoAnalitico1DSuParticella(data, particles[j], *FOG);
      fieldOnParticles.push_back(*FOG);

      pstar.setParticlePX(particles[j].px);
      pstar.setParticlePY(particles[j].py + dtp2 * 0.5 * fieldOnParticles[1].ey);

      gammastar = sqrt(1. + pow(pstar.px, 2) + pow(pstar.py, 2));
      b.setEy(0.5 * dtp2 * fieldOnParticles[1].ey / gammastar);
      coeff = 1. / (1. + pow(b.ey, 2));

      pHalfAdvanced.setParticlePX(coeff * (pstar.px + pstar.py * b.ey)); //?????
      pHalfAdvanced.setParticlePY(coeff * (pstar.py - pstar.px * b.ey)); //?????

      particles[j].setParticlePX(2. * pHalfAdvanced.px - particles[j].px);
      particles[j].setParticlePY(2. * pHalfAdvanced.py - particles[j].py);

      data.aumentaT(dtp2);

      // PASSO #3
      gamma = sqrt(1. + pow(particles[j].px, 2) + pow(particles[j].py, 2));
      particles[j].setParticleX(particles[j].x + dtx3 * particles[j].px / (gamma));
      particles[j].setParticleY(particles[j].y + dtx3 * particles[j].py / (gamma));

      calcolaCampoAnalitico1DSuParticella(data, particles[j], *FOG);
      fieldOnParticles.push_back(*FOG);

      pstar.setParticlePX(particles[j].px);
      pstar.setParticlePY(particles[j].py + dtp3 * 0.5 * fieldOnParticles[2].ey);

      gammastar = sqrt(1. + pow(pstar.px, 2) + pow(pstar.py, 2));
      b.setEy(0.5 * dtp3 * fieldOnParticles[2].ey / gammastar);
      coeff = 1. / (1. + pow(b.ey, 2));

      pHalfAdvanced.setParticlePX(coeff * (pstar.px + pstar.py * b.ey)); //?????
      pHalfAdvanced.setParticlePY(coeff * (pstar.py - pstar.px * b.ey)); //?????

      particles[j].setParticlePX(2. * pHalfAdvanced.px - particles[j].px);
      particles[j].setParticlePY(2. * pHalfAdvanced.py - particles[j].py);

      // PASSO #4 - SOLO EVOLUZIONE POSIZIONI
      gamma = sqrt(1. + pow(particles[j].px, 2) + pow(particles[j].py, 2));

      particles[j].setParticleX(particles[j].x + dtx1 * particles[j].px / gamma);
      particles[j].setParticleY(particles[j].y + dtx1 * particles[j].py / gamma);

      data.aumentaT(dtp3);

      //INTEGRALI PRIMI (da sistemare, li ho copiati dall'evoluzione sui campi analitici)
      potenziale.clear();
      zeta = particles[j].x - SPEED_OF_LIGHT * data.t; // non usato nella mia struttura delle equazioni di campo
      gamma = sqrt(1 + pow(particles[j].px, 2) + pow(particles[j].py, 2) + pow(particles[j].pz, 2));
      calcolaPotenzialeAnalitico1DSuParticella(data, particles[j], *campoSuPunto);
      potenziale.push_back(*campoSuPunto);

      I1 = (particles[j].py / (particles[j].m * SPEED_OF_LIGHT)) + potenziale[0].ey;
      I2 = particles[j].pz / (particles[j].m * SPEED_OF_LIGHT);
      I3 = gamma - particles[j].px / (particles[j].m * SPEED_OF_LIGHT);

      potenziale.clear();

      outputstream
          << setw(12) << i + 1 << "\t"
          << setprecision(14) << setw(16) << particles[j].x << "\t"
          << setw(16) << particles[j].y << "\t"
          << setw(16) << particles[j].z << "\t"
          //        << setw(16) << particles[j].getParticleCell() << "\t"
          << setw(16) << particles[j].px << "\t"
          << setw(16) << particles[j].py << "\t"
          << setw(16) << particles[j].pz << "\t"
          << setw(16) << fieldOnParticles[0].ey << "\t"
          << setw(16) << fieldOnParticles[1].ey << "\t"
          << setw(16) << fieldOnParticles[2].ey << "\t"
          << setw(16) << zeta << "\t"
          << setw(16) << I1 << "\t"
          << setw(16) << I2 << "\t"
          << setw(16) << I3 << endl;
    }
  }
}

void evolveLPF4_nogrid(Data data, vector<Particle>& particles, ofstream& outputstream)
{
  Field* campoSuPunto;
  campoSuPunto = new Field[1];

  double gamma, gamma_inv;
  double I1, I2, I3;
  double zeta;

  double alpha = 1. / (2. - pow(2., (1. / 3.)));
  double dtx1 = alpha * data.dt * 0.5;
  double dtp1 = alpha * data.dt;
  double dtx2 = (1. - alpha) * 0.5 * data.dt;
  double dtp2 = (1. - 2. * alpha) * data.dt;
  double dtx3 = (1. - alpha) * 0.5 * data.dt;
  double dtp3 = alpha * data.dt;

  //Il leapfrog di ordine 4 si ottiene componendo 3 schemi leapfrog invertiti e usando un peso alpha
  //   r(n+1/2) = r(n-1/2) + DeltaT(pesatoXciascunPassaggio) * p(n)/gamma(n)
  //   p(n+1) = p(n) + DeltaT(pesatoXciascunPassaggio) * F(n+1/2)
  //dove          F(n+1/2) = E(n+1/2) + p(n+1/2)/gamma(n+1/2) * B(n+1/2)
  //nel quale     p(n+1/2) = ( p(n) + p(n+1) ) / 2

  for (int i = 0; i < data.nSteps; i++) {
    for (int j = 0; j < data.n_electrons; j++) {
      // PASSO #1
      evolveLPF(data, particles[j], dtp1, dtx1);
      // PASSO #2
      evolveLPF(data, particles[j], dtp2, dtx2);
      // PASSO #3
      evolveLPF(data, particles[j], dtp3, dtx3);

      // PASSO #4 - evoluzione solo posizioni
      gamma = sqrt(1. + pow(particles[j].px, 2) + pow(particles[j].py, 2));
      gamma_inv = 1. / gamma;
      particles[j].setParticleX(particles[j].x + dtx1 * particles[j].px * gamma_inv);
      particles[j].setParticleY(particles[j].y + dtx1 * particles[j].py * gamma_inv);
      particles[j].setParticleZ(particles[j].z + dtx1 * particles[j].pz * gamma_inv);

      data.aumentaT(dtx1);

      zeta = particles[j].x - SPEED_OF_LIGHT * data.t; // non usato nella mia struttura delle equazioni di campo
      calcolaPotenzialeAnalitico1DSuParticella(data, particles[j], *campoSuPunto);

      I1 = (particles[j].py / (particles[j].m * SPEED_OF_LIGHT)) + (*campoSuPunto).ey;
      I2 = particles[j].pz / (particles[j].m * SPEED_OF_LIGHT);
      I3 = gamma - particles[j].px / (particles[j].m * SPEED_OF_LIGHT);

      outputstream
          << setw(12) << i + 1 << "\t"
          << setprecision(14)
          << setw(16) << particles[j].x << "\t"
          << setw(16) << particles[j].y << "\t"
          << setw(16) << particles[j].z << "\t"
          << setw(16) << particles[j].px << "\t"
          << setw(16) << particles[j].py << "\t"
          << setw(16) << particles[j].pz << "\t"
          << setw(16) << zeta << "\t"
          << setw(16) << I1 << "\t"
          << setw(16) << I2 << "\t"
          << setw(16) << I3 << endl;

      data.aumentaT(-data.dt); // cosi` un altro elettrone dello stesso step non si trova i tempi avanzati
    }
    data.aumentaT(data.dt); // e` questo il vero avanzamento di uno step temporale del tempo del sistema
  }
}

void evolveLPF4_withgrid_onthefly(Data data, vector<Particle>& particles, vector<Field>& fieldOnParticles, ofstream& outputstream)
{
  vector<Field> potenziale;

  Particle pstar, pHalfAdvanced;
  Particle coordinatePrimoPuntoGriglia;
  Particle coordinateSecondoPuntoGriglia;
  Particle coordinateTerzoPuntoGriglia;
  Particle coordinateQuartoPuntoGriglia;
  Particle coordinateQuintoPuntoGriglia;
  Particle coordinateSestoPuntoGriglia;
  Particle coordinateSettimoPuntoGriglia;
  Particle coordinateOttavoPuntoGriglia;

  int posX = 0, posY = 0, posZ = 0;
  Field b;
  Field FOP;
  vector<Field> fieldOnGrid;

  Field* FOG; // DA SISTEMARE, c'è solo per avere un campo fittizio su una griglia on-the-fly per implementare l'algoritmo
  FOG = new Field[1];

  double gridXX, gridYY, gridZZ;
  int gridX, gridY, gridZ;
  double gamma, gamma_inv, b2;
  double I1, I2, I3;
  double zeta;

  double alpha = 1. / (2. - pow(2., (1. / 3.)));
  double dtx1 = 0.5 * alpha * data.dt;
  double dtp1 = alpha * data.dt;
  double dtx2 = (1. - alpha) * 0.5 * data.dt;
  double dtp2 = (1. - 2 * alpha) * data.dt;
  double dtx3 = (1. - alpha) * 0.5 * data.dt;
  double dtp3 = alpha * data.dt;

  Field* campoSuPunto;
  campoSuPunto = new Field[1];

  //ALGORITMO LEAPFROG PARTICELLE

  for (int i = 0; i < data.nSteps; i++) {
    for (int j = 0; j < data.n_electrons; j++) {
      gamma = sqrt(1 + pow(particles[j].px, 2) + pow(particles[j].py, 2) + pow(particles[j].pz, 2));
      gamma_inv = 1. / gamma;

      particles[j].setParticleX(particles[j].x + dtx1 * particles[j].px * gamma_inv);
      particles[j].setParticleY(particles[j].y + dtx1 * particles[j].py * gamma_inv);
      particles[j].setParticleZ(particles[j].z + dtx1 * particles[j].pz * gamma_inv);

      fieldOnParticles.clear();

      // STEP #1

      gridXX = particles[j].x / data.deltaX;
      if (gridXX < 0)
        gridX = ((int)gridXX) - 1; //il punto di griglia lungo l'asse x a sinistra della particella
      else
        gridX = (int)gridXX;

      gridYY = particles[j].y / data.deltaY;
      if (gridYY < 0)
        gridY = ((int)gridYY) - 1; //il punto di griglia lungo l'asse y a sinistra della particella
      else
        gridY = (int)gridYY;

      gridZZ = particles[j].z / data.deltaZ;
      if (gridZZ < 0)
        gridZ = ((int)gridZZ) - 1; //il punto di griglia lungo l'asse z a sinistra della particella
      else
        gridZ = (int)gridZZ;

      coordinatePrimoPuntoGriglia.setParticleX(gridX * data.deltaX);
      coordinatePrimoPuntoGriglia.setParticleY(gridY * data.deltaY);
      coordinatePrimoPuntoGriglia.setParticleZ(gridZ * data.deltaZ);
      coordinateSecondoPuntoGriglia.setParticleX(gridX * data.deltaX + data.deltaX);
      coordinateSecondoPuntoGriglia.setParticleY(gridY * data.deltaY);
      coordinateSecondoPuntoGriglia.setParticleZ(gridZ * data.deltaZ);
      coordinateTerzoPuntoGriglia.setParticleX(gridX * data.deltaX);
      coordinateTerzoPuntoGriglia.setParticleY(gridY * data.deltaY + data.deltaY);
      coordinateTerzoPuntoGriglia.setParticleZ(gridZ * data.deltaZ);
      coordinateQuartoPuntoGriglia.setParticleX(gridX * data.deltaX + data.deltaX);
      coordinateQuartoPuntoGriglia.setParticleY(gridY * data.deltaY + data.deltaY);
      coordinateQuartoPuntoGriglia.setParticleZ(gridZ * data.deltaZ);
      /*
      coordinateQuintoPuntoGriglia.setParticleX(gridX*data.deltaX);
      coordinateQuintoPuntoGriglia.setParticleY(gridY*data.deltaY);
      coordinateQuintoPuntoGriglia.setParticleZ(gridZ*data.deltaZ+data.deltaZ);
      coordinateSestoPuntoGriglia.setParticleX(gridX*data.deltaX+data.deltaX);
      coordinateSestoPuntoGriglia.setParticleY(gridY*data.deltaY);
      coordinateSestoPuntoGriglia.setParticleZ(gridZ*data.deltaZ+data.deltaZ);
      coordinateSettimoPuntoGriglia.setParticleX(gridX*data.deltaX);
      coordinateSettimoPuntoGriglia.setParticleY(gridY*data.deltaY+data.deltaY);
      coordinateSettimoPuntoGriglia.setParticleZ(gridZ*data.deltaZ+data.deltaZ);
      coordinateOttavoPuntoGriglia.setParticleX(gridX*data.deltaX+data.deltaX);
      coordinateOttavoPuntoGriglia.setParticleY(gridY*data.deltaY+data.deltaY);
      coordinateOttavoPuntoGriglia.setParticleZ(gridZ*data.deltaZ+data.deltaZ);
      */

      data.aumentaT(dtx1);

      calcolaCampoAnalitico1DSuParticella(data, coordinatePrimoPuntoGriglia, *FOG);
      fieldOnGrid.push_back(*FOG);
      calcolaCampoAnalitico1DSuParticella(data, coordinateSecondoPuntoGriglia, *FOG);
      fieldOnGrid.push_back(*FOG);
      calcolaCampoAnalitico1DSuParticella(data, coordinateTerzoPuntoGriglia, *FOG);
      fieldOnGrid.push_back(*FOG);
      calcolaCampoAnalitico1DSuParticella(data, coordinateQuartoPuntoGriglia, *FOG);
      fieldOnGrid.push_back(*FOG);
      /*
      calcolaCampoAnalitico1DSuParticella(data, coordinateQuintoPuntoGriglia, *FOG);
      fieldOnGrid.push_back(*FOG);
      calcolaCampoAnalitico1DSuParticella(data, coordinateSestoPuntoGriglia, *FOG);
      fieldOnGrid.push_back(*FOG);
      calcolaCampoAnalitico1DSuParticella(data, coordinateSettimoPuntoGriglia, *FOG);
      fieldOnGrid.push_back(*FOG);
      calcolaCampoAnalitico1DSuParticella(data, coordinateOttavoPuntoGriglia, *FOG);
      fieldOnGrid.push_back(*FOG);
      */

      FOP = interpolation2D_linear(data, particles[j], fieldOnGrid);
      fieldOnParticles.push_back(FOP);

      fieldOnGrid.clear();

      //Il leapfrog di ordine 4 si ottiene componendo 3 schemi leapfrog invertiti e usando un peso alpha
      //   r(n+1/2) = r(n-1/2) + DeltaT(pesatoXciascunPassaggio) * p(n)/gamma(n)
      //   p(n+1) = p(n) + DeltaT(pesatoXciascunPassaggio) * F(n+1/2)
      //dove          F(n+1/2) = E(n+1/2) + p(n+1/2)/gamma(n+1/2) * B(n+1/2)
      //nel quale     p(n+1/2) = ( p(n) + p(n+1) ) / 2

      //pstar = p(n) + Q * E(n+1/2) * DeltaT/2
      pstar.setParticlePX(particles[j].px + 0.5 * dtp1 * particles[j].q * fieldOnParticles[0].ex);
      pstar.setParticlePY(particles[j].py + 0.5 * dtp1 * particles[j].q * fieldOnParticles[0].ey);
      pstar.setParticlePZ(particles[j].pz + 0.5 * dtp1 * particles[j].q * fieldOnParticles[0].ez);
      gamma = sqrt(1 + pow(pstar.px, 2) + pow(pstar.py, 2) + pow(pstar.pz, 2));
      gamma_inv = 1. / gamma;

      //b = (DeltaT/2) * Q * B(n+1/2) / gamma(aggiornato al pstar)
      b.setBx(0.5 * dtp1 * fieldOnParticles[0].bx * particles[j].q * gamma_inv);
      b.setBy(0.5 * dtp1 * fieldOnParticles[0].by * particles[j].q * gamma_inv);
      b.setBz(0.5 * dtp1 * fieldOnParticles[0].bz * particles[j].q * gamma_inv);
      b2 = 1. + pow(b.bx, 2) + pow(b.by, 2) + pow(b.bz, 2);

      /****************************************************************************
      //p(n+1/2) = [pstar + pstar x b + (pstar x b) * b] / b2   SECONDO TURCHETTI  (ma i conti analitici diconon che e` sbagliato, anche perche' l'ultimo termine e` uno scalare!)
      pHalfAdvanced.setParticlePX( (pstar.px + pstar.py * b.bz - b.by * pstar.pz + (pstar.py * b.bz - b.by * pstar.pz) * b.bx ) / b2 );
      pHalfAdvanced.setParticlePY( (pstar.py + pstar.pz * b.bx - b.bz * pstar.px + (pstar.pz * b.bx - b.bz * pstar.px) * b.by ) / b2 );
      pHalfAdvanced.setParticlePZ( (pstar.pz + pstar.px * b.by - b.bx * pstar.py + (pstar.px * b.by - b.bx * pstar.py) * b.bz ) / b2 );
      ****************************************************************************/

      /****************************************************************************
      //p(n+1/2) = [pstar + pstar x b + b x (pstar x b)] / b2   SECONDO LONDRILLO
      pHalfAdvanced.setParticlePX( (pstar.px + pstar.py * b.bz - b.by * pstar.pz + b.bx * pstar.px * b.bx ) / b2 );
      pHalfAdvanced.setParticlePY( (pstar.py + pstar.pz * b.bx - b.bz * pstar.px + b.by * pstar.py * b.by ) / b2 );
      pHalfAdvanced.setParticlePZ( (pstar.pz + pstar.px * b.by - b.bx * pstar.py + b.bz * pstar.pz * b.bz ) / b2 );
      ****************************************************************************/

      // /****************************************************************************
      //p(n+1/2) = [pstar + pstar x b] / b2   COME IMPLEMENTATO DA TURCHETTI NEL CODICE
      pHalfAdvanced.setParticlePX((pstar.px + pstar.py * b.bz) / b2);
      pHalfAdvanced.setParticlePY((pstar.py - pstar.px * b.bz) / b2);
      pHalfAdvanced.setParticlePZ(0.); //per come ho strutturato i pacchetti d'onda su z non dovrebbe succedere nulla e quindi lo fisso io a zero.
      // ****************************************************************************/

      particles[j].setParticlePX(2 * pHalfAdvanced.px - particles[j].px);
      particles[j].setParticlePY(2 * pHalfAdvanced.py - particles[j].py);
      particles[j].setParticlePZ(2 * pHalfAdvanced.pz - particles[j].pz);

      // STEP #2
      gamma = sqrt(1 + pow(particles[j].px, 2) + pow(particles[j].py, 2) + pow(particles[j].pz, 2));
      gamma_inv = 1. / gamma;

      particles[j].setParticleX(particles[j].x + dtx2 * particles[j].px * gamma_inv);
      particles[j].setParticleY(particles[j].y + dtx2 * particles[j].py * gamma_inv);
      particles[j].setParticleZ(particles[j].z + dtx2 * particles[j].pz * gamma_inv);

      gridXX = particles[j].x / data.deltaX;
      if (gridXX < 0)
        gridX = ((int)gridXX) - 1; //il punto di griglia lungo l'asse x a sinistra della particella
      else
        gridX = (int)gridXX;

      gridYY = particles[j].y / data.deltaY;
      if (gridYY < 0)
        gridY = ((int)gridYY) - 1; //il punto di griglia lungo l'asse y a sinistra della particella
      else
        gridY = (int)gridYY;

      gridZZ = particles[j].z / data.deltaZ;
      if (gridZZ < 0)
        gridZ = ((int)gridZZ) - 1; //il punto di griglia lungo l'asse z a sinistra della particella
      else
        gridZ = (int)gridZZ;

      coordinatePrimoPuntoGriglia.setParticleX(gridX * data.deltaX);
      coordinatePrimoPuntoGriglia.setParticleY(gridY * data.deltaY);
      coordinatePrimoPuntoGriglia.setParticleZ(gridZ * data.deltaZ);
      coordinateSecondoPuntoGriglia.setParticleX(gridX * data.deltaX + data.deltaX);
      coordinateSecondoPuntoGriglia.setParticleY(gridY * data.deltaY);
      coordinateSecondoPuntoGriglia.setParticleZ(gridZ * data.deltaZ);
      coordinateTerzoPuntoGriglia.setParticleX(gridX * data.deltaX);
      coordinateTerzoPuntoGriglia.setParticleY(gridY * data.deltaY + data.deltaY);
      coordinateTerzoPuntoGriglia.setParticleZ(gridZ * data.deltaZ);
      coordinateQuartoPuntoGriglia.setParticleX(gridX * data.deltaX + data.deltaX);
      coordinateQuartoPuntoGriglia.setParticleY(gridY * data.deltaY + data.deltaY);
      coordinateQuartoPuntoGriglia.setParticleZ(gridZ * data.deltaZ);

      data.aumentaT(dtx2);

      calcolaCampoAnalitico1DSuParticella(data, coordinatePrimoPuntoGriglia, *FOG);
      fieldOnGrid.push_back(*FOG);
      calcolaCampoAnalitico1DSuParticella(data, coordinateSecondoPuntoGriglia, *FOG);
      fieldOnGrid.push_back(*FOG);
      calcolaCampoAnalitico1DSuParticella(data, coordinateTerzoPuntoGriglia, *FOG);
      fieldOnGrid.push_back(*FOG);
      calcolaCampoAnalitico1DSuParticella(data, coordinateQuartoPuntoGriglia, *FOG);
      fieldOnGrid.push_back(*FOG);

      FOP = interpolation2D_linear(data, particles[j], fieldOnGrid);
      fieldOnParticles.push_back(FOP);

      fieldOnGrid.clear();

      //pstar = p(n) + Q * E(n+1/2) * DeltaT/2
      pstar.setParticlePX(particles[j].px + 0.5 * dtp2 * particles[j].q * fieldOnParticles[1].ex);
      pstar.setParticlePY(particles[j].py + 0.5 * dtp2 * particles[j].q * fieldOnParticles[1].ey);
      pstar.setParticlePZ(particles[j].pz + 0.5 * dtp2 * particles[j].q * fieldOnParticles[1].ez);
      gamma = sqrt(1 + pow(pstar.px, 2) + pow(pstar.py, 2) + pow(pstar.pz, 2));
      gamma_inv = 1. / gamma;

      //b = (DeltaT/2) * Q * B(n+1/2) / gamma(aggiornato al pstar)
      b.setBx(0.5 * dtp2 * fieldOnParticles[1].bx * particles[j].q * gamma_inv);
      b.setBy(0.5 * dtp2 * fieldOnParticles[1].by * particles[j].q * gamma_inv);
      b.setBz(0.5 * dtp2 * fieldOnParticles[1].bz * particles[j].q * gamma_inv);
      b2 = 1. + pow(b.bx, 2) + pow(b.by, 2) + pow(b.bz, 2);

      /****************************************************************************
      //p(n+1/2) = [pstar + pstar x b + (pstar x b) * b] / b2   SECONDO TURCHETTI  (ma i conti analitici diconon che e` sbagliato, anche perche' l'ultimo termine e` uno scalare!)
      pHalfAdvanced.setParticlePX( (pstar.px + pstar.py * b.bz - b.by * pstar.pz + (pstar.py * b.bz - b.by * pstar.pz) * b.bx ) / b2 );
      pHalfAdvanced.setParticlePY( (pstar.py + pstar.pz * b.bx - b.bz * pstar.px + (pstar.pz * b.bx - b.bz * pstar.px) * b.by ) / b2 );
      pHalfAdvanced.setParticlePZ( (pstar.pz + pstar.px * b.by - b.bx * pstar.py + (pstar.px * b.by - b.bx * pstar.py) * b.bz ) / b2 );
      ****************************************************************************/

      /****************************************************************************
      //p(n+1/2) = [pstar + pstar x b + b x (pstar x b)] / b2   SECONDO LONDRILLO
      pHalfAdvanced.setParticlePX( (pstar.px + pstar.py * b.bz - b.by * pstar.pz + b.bx * pstar.px * b.bx ) / b2 );
      pHalfAdvanced.setParticlePY( (pstar.py + pstar.pz * b.bx - b.bz * pstar.px + b.by * pstar.py * b.by ) / b2 );
      pHalfAdvanced.setParticlePZ( (pstar.pz + pstar.px * b.by - b.bx * pstar.py + b.bz * pstar.pz * b.bz ) / b2 );
      ****************************************************************************/

      // /****************************************************************************
      //p(n+1/2) = [pstar + pstar x b] / b2   COME IMPLEMENTATO DA TURCHETTI NEL CODICE
      pHalfAdvanced.setParticlePX((pstar.px + pstar.py * b.bz) / b2);
      pHalfAdvanced.setParticlePY((pstar.py - pstar.px * b.bz) / b2);
      pHalfAdvanced.setParticlePZ(0.); //per come ho strutturato i pacchetti d'onda su z non dovrebbe succedere nulla e quindi lo fisso io a zero.
      // ****************************************************************************/

      particles[j].setParticlePX(2 * pHalfAdvanced.px - particles[j].px);
      particles[j].setParticlePY(2 * pHalfAdvanced.py - particles[j].py);
      particles[j].setParticlePZ(2 * pHalfAdvanced.pz - particles[j].pz);

      // STEP #3
      gamma = sqrt(1 + pow(particles[j].px, 2) + pow(particles[j].py, 2) + pow(particles[j].pz, 2));
      gamma_inv = 1. / gamma;

      particles[j].setParticleX(particles[j].x + dtx3 * particles[j].px * gamma_inv);
      particles[j].setParticleY(particles[j].y + dtx3 * particles[j].py * gamma_inv);
      particles[j].setParticleZ(particles[j].z + dtx3 * particles[j].pz * gamma_inv);

      gridXX = particles[j].x / data.deltaX;
      if (gridXX < 0)
        gridX = ((int)gridXX) - 1; //il punto di griglia lungo l'asse x a sinistra della particella
      else
        gridX = (int)gridXX;

      gridYY = particles[j].y / data.deltaY;
      if (gridYY < 0)
        gridY = ((int)gridYY) - 1; //il punto di griglia lungo l'asse y a sinistra della particella
      else
        gridY = (int)gridYY;

      gridZZ = particles[j].z / data.deltaZ;
      if (gridZZ < 0)
        gridZ = ((int)gridZZ) - 1; //il punto di griglia lungo l'asse z a sinistra della particella
      else
        gridZ = (int)gridZZ;

      coordinatePrimoPuntoGriglia.setParticleX(gridX * data.deltaX);
      coordinatePrimoPuntoGriglia.setParticleY(gridY * data.deltaY);
      coordinatePrimoPuntoGriglia.setParticleZ(gridZ * data.deltaZ);
      coordinateSecondoPuntoGriglia.setParticleX(gridX * data.deltaX + data.deltaX);
      coordinateSecondoPuntoGriglia.setParticleY(gridY * data.deltaY);
      coordinateSecondoPuntoGriglia.setParticleZ(gridZ * data.deltaZ);
      coordinateTerzoPuntoGriglia.setParticleX(gridX * data.deltaX);
      coordinateTerzoPuntoGriglia.setParticleY(gridY * data.deltaY + data.deltaY);
      coordinateTerzoPuntoGriglia.setParticleZ(gridZ * data.deltaZ);
      coordinateQuartoPuntoGriglia.setParticleX(gridX * data.deltaX + data.deltaX);
      coordinateQuartoPuntoGriglia.setParticleY(gridY * data.deltaY + data.deltaY);
      coordinateQuartoPuntoGriglia.setParticleZ(gridZ * data.deltaZ);

      data.aumentaT(dtx3);

      calcolaCampoAnalitico1DSuParticella(data, coordinatePrimoPuntoGriglia, *FOG);
      fieldOnGrid.push_back(*FOG);
      calcolaCampoAnalitico1DSuParticella(data, coordinateSecondoPuntoGriglia, *FOG);
      fieldOnGrid.push_back(*FOG);
      calcolaCampoAnalitico1DSuParticella(data, coordinateTerzoPuntoGriglia, *FOG);
      fieldOnGrid.push_back(*FOG);
      calcolaCampoAnalitico1DSuParticella(data, coordinateQuartoPuntoGriglia, *FOG);
      fieldOnGrid.push_back(*FOG);

      FOP = interpolation2D_linear(data, particles[j], fieldOnGrid);
      fieldOnParticles.push_back(FOP);

      fieldOnGrid.clear();

      //pstar = p(n) + Q * E(n+1/2) * DeltaT/2
      pstar.setParticlePX(particles[j].px + 0.5 * dtp3 * particles[j].q * fieldOnParticles[2].ex);
      pstar.setParticlePY(particles[j].py + 0.5 * dtp3 * particles[j].q * fieldOnParticles[2].ey);
      pstar.setParticlePZ(particles[j].pz + 0.5 * dtp3 * particles[j].q * fieldOnParticles[2].ez);
      gamma = sqrt(1 + pow(pstar.px, 2) + pow(pstar.py, 2) + pow(pstar.pz, 2));
      gamma_inv = 1. / gamma;

      //b = (DeltaT/2) * Q * B(n+1/2) / gamma(aggiornato al pstar)
      b.setBx(0.5 * dtp3 * fieldOnParticles[2].bx * particles[j].q * gamma_inv);
      b.setBy(0.5 * dtp3 * fieldOnParticles[2].by * particles[j].q * gamma_inv);
      b.setBz(0.5 * dtp3 * fieldOnParticles[2].bz * particles[j].q * gamma_inv);
      b2 = 1. + pow(b.bx, 2) + pow(b.by, 2) + pow(b.bz, 2);

      /****************************************************************************
      //p(n+1/2) = [pstar + pstar x b + (pstar x b) * b] / b2   SECONDO TURCHETTI  (ma i conti analitici diconon che e` sbagliato, anche perche' l'ultimo termine e` uno scalare!)
      pHalfAdvanced.setParticlePX( (pstar.px + pstar.py * b.bz - b.by * pstar.pz + (pstar.py * b.bz - b.by * pstar.pz) * b.bx ) / b2 );
      pHalfAdvanced.setParticlePY( (pstar.py + pstar.pz * b.bx - b.bz * pstar.px + (pstar.pz * b.bx - b.bz * pstar.px) * b.by ) / b2 );
      pHalfAdvanced.setParticlePZ( (pstar.pz + pstar.px * b.by - b.bx * pstar.py + (pstar.px * b.by - b.bx * pstar.py) * b.bz ) / b2 );
      ****************************************************************************/

      /****************************************************************************
      //p(n+1/2) = [pstar + pstar x b + b x (pstar x b)] / b2   SECONDO LONDRILLO
      pHalfAdvanced.setParticlePX( (pstar.px + pstar.py * b.bz - b.by * pstar.pz + b.bx * pstar.px * b.bx ) / b2 );
      pHalfAdvanced.setParticlePY( (pstar.py + pstar.pz * b.bx - b.bz * pstar.px + b.by * pstar.py * b.by ) / b2 );
      pHalfAdvanced.setParticlePZ( (pstar.pz + pstar.px * b.by - b.bx * pstar.py + b.bz * pstar.pz * b.bz ) / b2 );
      ****************************************************************************/

      // /****************************************************************************
      //p(n+1/2) = [pstar + pstar x b] / b2   COME IMPLEMENTATO DA TURCHETTI NEL CODICE
      pHalfAdvanced.setParticlePX((pstar.px + pstar.py * b.bz) / b2);
      pHalfAdvanced.setParticlePY((pstar.py - pstar.px * b.bz) / b2);
      pHalfAdvanced.setParticlePZ(0.); //per come ho strutturato i pacchetti d'onda su z non dovrebbe succedere nulla e quindi lo fisso io a zero.
      // ****************************************************************************/

      particles[j].setParticlePX(2 * pHalfAdvanced.px - particles[j].px);
      particles[j].setParticlePY(2 * pHalfAdvanced.py - particles[j].py);
      particles[j].setParticlePZ(2 * pHalfAdvanced.pz - particles[j].pz);

      // PASSO #4 - evoluzione solo posizioni
      gamma = sqrt(1. + pow(particles[j].px, 2) + pow(particles[j].py, 2));
      gamma_inv = 1. / gamma;

      particles[j].setParticleX(particles[j].x + dtx1 * particles[j].px * gamma_inv);
      particles[j].setParticleY(particles[j].y + dtx1 * particles[j].py * gamma_inv);
      particles[j].setParticleZ(particles[j].z + dtx1 * particles[j].pz * gamma_inv);

      data.aumentaT(dtx1);

      //INTEGRALI PRIMI (da sistemare, li ho copiati dall'evoluzione sui campi analitici)
      potenziale.clear();
      zeta = particles[j].x - SPEED_OF_LIGHT * data.t; // non usato nella mia struttura delle equazioni di campo
      gamma = sqrt(1 + pow(particles[j].px, 2) + pow(particles[j].py, 2) + pow(particles[j].pz, 2));
      calcolaPotenzialeAnalitico1DSuParticella(data, particles[j], *campoSuPunto);
      potenziale.push_back(*campoSuPunto);

      I1 = (particles[j].py / (particles[j].m * SPEED_OF_LIGHT)) + potenziale[0].ey;
      I2 = particles[j].pz / (particles[j].m * SPEED_OF_LIGHT);
      I3 = gamma - particles[j].px / (particles[j].m * SPEED_OF_LIGHT);

      potenziale.clear();

      outputstream
          << setw(12) << i + 1 << "\t"
          << setprecision(14) << setw(16) << particles[j].x << "\t"
          << setw(16) << particles[j].y << "\t"
          << setw(16) << particles[j].z << "\t"
          // << setw(16) << particles[j].getParticleCell() << "\t"
          << setw(16) << particles[j].px << "\t"
          << setw(16) << particles[j].py << "\t"
          << setw(16) << particles[j].pz << "\t"
          << setw(16) << zeta << "\t"
          << setw(16) << I1 << "\t"
          << setw(16) << I2 << "\t"
          << setw(16) << I3 << endl;
    }
  }
}

void evolveLPF4_withgrid(Data data, vector<Particle>& particles, vector<Field>& fieldOnParticles, vector<Field>& campoSuPuntiGriglia, ofstream& outputstream)
{

  //AL MOMENTO NON E' CORRETTAMENTE IMPLEMENTATA, FA LE STESSE COSE DI evolveLPF_withgrid (_nemmeno_ quello di ordine 4, quindi!)
  //OLTRETUTTO non e` nemmeno la versione piu` aggiornata di evolveLPF_withgrid, quindi e` pure completamente buggato!

  vector<Field> potenzialeSuPuntiGriglia;

  vector<Field> potenziale;
  vector<Field> fieldOnGrid;

  Particle pstar, pHalfAdvanced;

  int posX = 0, posY = 0, posZ = 0;
  Field b;
  Field FOP;

  double gridXX, gridYY, gridZZ;
  int gridX, gridY, gridZ;
  double gamma, gamma_inv, b2;
  double I1, I2, I3;
  double zeta;

  Field* campoSuPunto;
  campoSuPunto = new Field[1];

  //ALGORITMO LEAPFROG PARTICELLE

  if (data.n_dim == 1) {
    for (int i = 0; i < data.nSteps; i++) {
      data.aumentaT(data.dt);

      campoSuPuntiGriglia.clear();
      riempiPuntiGrigliaConCampoAnalitico1D(data, campoSuPuntiGriglia);

      potenzialeSuPuntiGriglia.clear();
      riempiPuntiGrigliaConPotenzialeAnalitico1D(data, potenzialeSuPuntiGriglia);

      for (int j = 0; j < data.n_electrons; j++) {
        gamma = sqrt(1. + pow(particles[j].px, 2) + pow(particles[j].py, 2));
        gamma_inv = 1. / gamma;

        particles[j].setParticleX(particles[j].x + 0.5 * data.dt * particles[j].px * gamma_inv);
        particles[j].setParticleY(particles[j].y + 0.5 * data.dt * particles[j].py * gamma_inv);
        particles[j].setParticleZ(particles[j].z + 0.5 * data.dt * particles[j].pz * gamma_inv);

        //QUI MI CALCOLO I PUNTI DI RIFERIMENTO PER LA PARTICELLA IN ESAME nella posizione attuale
        //forse c'e` un modo piu` furbo per farlo sfruttando il # di cella in cui si trova che e` noto

        fieldOnParticles.clear(); //SVUOTA il vettore dei campi sulle particelle ad ogni ciclo

        gridXX = particles[j].x / data.deltaX;
        if (gridXX < 0)
          gridX = ((int)gridXX) - 1; //il punto di griglia lungo l'asse x a sinistra della particella
        else
          gridX = (int)gridXX;

        gridYY = particles[j].y / data.deltaY;
        if (gridYY < 0)
          gridY = ((int)gridYY) - 1; //il punto di griglia lungo l'asse y a sinistra della particella
        else
          gridY = (int)gridYY;

        gridZZ = particles[j].z / data.deltaZ;
        if (gridZZ < 0)
          gridZ = ((int)gridZZ) - 1; //il punto di griglia lungo l'asse z a sinistra della particella
        else
          gridZ = (int)gridZZ;

        if (gridX < 0 || gridY < 0 || gridZ < 0) {
          cout << "La particella e` uscita dai boundaries - simulazione interrotta" << endl;
          return;
        }

        fieldOnGrid.push_back(campoSuPuntiGriglia.at(gridY * data.nGridPointsX + gridX)); // primo punto
        fieldOnGrid.push_back(campoSuPuntiGriglia.at(gridY * data.nGridPointsX + (gridX + 1))); // secondo punto
        fieldOnGrid.push_back(campoSuPuntiGriglia.at((gridY + 1) * data.nGridPointsX + gridX)); // terzo punto
        fieldOnGrid.push_back(campoSuPuntiGriglia.at((gridY + 1) * data.nGridPointsX + (gridX + 1))); // quarto punto

        FOP = interpolation2D_linear(data, particles[j], fieldOnGrid);
        fieldOnParticles.push_back(FOP);

        fieldOnGrid.clear();

        //pstar = p(n-1/2) + Q * E(n) * DeltaT/2
        pstar.setParticlePX(particles[j].px + 0.5 * data.dt * particles[j].q * fieldOnParticles[0].ex);
        pstar.setParticlePY(particles[j].py + 0.5 * data.dt * particles[j].q * fieldOnParticles[0].ey);
        pstar.setParticlePZ(particles[j].pz + 0.5 * data.dt * particles[j].q * fieldOnParticles[0].ez);
        gamma = sqrt(1 + pow(pstar.px, 2) + pow(pstar.py, 2) + pow(pstar.pz, 2));
        gamma_inv = 1. / gamma;

        //b = (DeltaT/2) * Q * B(n) / gamma
        b.setBx(0.5 * data.dt * fieldOnParticles[0].bx * particles[j].q * gamma_inv);
        b.setBy(0.5 * data.dt * fieldOnParticles[0].by * particles[j].q * gamma_inv);
        b.setBz(0.5 * data.dt * fieldOnParticles[0].bz * particles[j].q * gamma_inv);
        b2 = pow(b.bx, 2) + pow(b.by, 2) + pow(b.bz, 2);
        b2 += 1.;

        /****************************************************************************
        //p(n+1/2) = [pstar + pstar x b + (pstar x b) * b] / b2   SECONDO TURCHETTI  (ma i conti analitici diconon che e` sbagliato, anche perche' l'ultimo termine e` uno scalare!)
        pHalfAdvanced.setParticlePX( (pstar.px + pstar.py * b.bz - b.by * pstar.pz + (pstar.py * b.bz - b.by * pstar.pz) * b.bx ) / b2 );
        pHalfAdvanced.setParticlePY( (pstar.py + pstar.pz * b.bx - b.bz * pstar.px + (pstar.pz * b.bx - b.bz * pstar.px) * b.by ) / b2 );
        pHalfAdvanced.setParticlePZ( (pstar.pz + pstar.px * b.by - b.bx * pstar.py + (pstar.px * b.by - b.bx * pstar.py) * b.bz ) / b2 );
        ****************************************************************************/

        /****************************************************************************
        //p(n+1/2) = [pstar + pstar x b + b x (pstar  x b)] / b2   SECONDO LONDRILLO
        pHalfAdvanced.setParticlePX( (pstar.px + pstar.py * b.bz - b.by * pstar.pz + b.bx * pstar.px * b.bx ) / b2 );
        pHalfAdvanced.setParticlePY( (pstar.py + pstar.pz * b.bx - b.bz * pstar.px + b.by * pstar.py * b.by ) / b2 );
        pHalfAdvanced.setParticlePZ( (pstar.pz + pstar.px * b.by - b.bx * pstar.py + b.bz * pstar.pz * b.bz ) / b2 );
        ****************************************************************************/

        //      /****************************************************************************
        //p(n+1/2) = [pstar + pstar x b] / b2   COME IMPLEMENTATO DA TURCHETTI NEL CODICE
        pHalfAdvanced.setParticlePX((pstar.px + pstar.py * b.bz) / b2);
        pHalfAdvanced.setParticlePY((pstar.py - pstar.px * b.bz) / b2);
        pHalfAdvanced.setParticlePZ(0.); //per come ho strutturato i pacchetti d'onda su z non dovrebbe succedere nulla e quindi lo fisso io a zero.
        //      ****************************************************************************/

        particles[j].setParticlePX(2 * pHalfAdvanced.px - particles[j].px);
        particles[j].setParticlePY(2 * pHalfAdvanced.py - particles[j].py);
        particles[j].setParticlePZ(2 * pHalfAdvanced.pz - particles[j].pz);

        gamma = sqrt(1. + pow(particles[j].px, 2) + pow(particles[j].py, 2));
        gamma_inv = 1. / gamma;

        particles[j].setParticleX(particles[j].x + 0.5 * data.dt * particles[j].px * gamma_inv);
        particles[j].setParticleY(particles[j].y + 0.5 * data.dt * particles[j].py * gamma_inv);
        particles[j].setParticleZ(particles[j].z + 0.5 * data.dt * particles[j].pz * gamma_inv);

        //INTEGRALI PRIMI (da sistemare, li ho copiati dall'evoluzione sui campi analitici)
        potenziale.clear();
        zeta = particles[j].x - SPEED_OF_LIGHT * data.t; // non usato nella mia struttura delle equazioni di campo
        gamma = sqrt(1 + pow(particles[j].px, 2) + pow(particles[j].py, 2) + pow(particles[j].pz, 2));

        gridXX = particles[j].x / data.deltaX;
        if (gridXX < 0)
          gridX = ((int)gridXX) - 1; //il punto di griglia lungo l'asse x a sinistra della particella
        else
          gridX = (int)gridXX;

        gridYY = particles[j].y / data.deltaY;
        if (gridYY < 0)
          gridY = ((int)gridYY) - 1; //il punto di griglia lungo l'asse y a sinistra della particella
        else
          gridY = (int)gridYY;

        gridZZ = particles[j].z / data.deltaZ;
        if (gridZZ < 0)
          gridZ = ((int)gridZZ) - 1; //il punto di griglia lungo l'asse z a sinistra della particella
        else
          gridZ = (int)gridZZ;

        if (gridX < 0 || gridY < 0 || gridZ < 0) {
          cout << "La particella e` uscita dai boundaries - simulazione interrotta" << endl;
          return;
        }

        fieldOnGrid.push_back(potenzialeSuPuntiGriglia.at(gridY * data.nGridPointsX + gridX)); // primo punto
        fieldOnGrid.push_back(potenzialeSuPuntiGriglia.at(gridY * data.nGridPointsX + (gridX + 1))); // secondo punto
        fieldOnGrid.push_back(potenzialeSuPuntiGriglia.at((gridY + 1) * data.nGridPointsX + gridX)); // terzo punto
        fieldOnGrid.push_back(potenzialeSuPuntiGriglia.at((gridY + 1) * data.nGridPointsX + (gridX + 1))); // quarto punto

        FOP = interpolation2D_linear(data, particles[j], fieldOnGrid);
        potenziale.push_back(FOP);

        fieldOnGrid.clear();

        I1 = (particles[j].py / (particles[j].m * SPEED_OF_LIGHT)) + potenziale[0].ey;
        I2 = particles[j].pz / (particles[j].m * SPEED_OF_LIGHT);
        I3 = gamma - particles[j].px / (particles[j].m * SPEED_OF_LIGHT);

        outputstream
            << setw(12) << i + 1 << "\t"
            << setprecision(14) << setw(16) << particles[j].x << "\t"
            << setw(16) << particles[j].y << "\t"
            << setw(16) << particles[j].z << "\t"
            // << setw(16) << particles[j].getParticleCell() << "\t"
            << setw(16) << particles[j].px << "\t"
            << setw(16) << particles[j].py << "\t"
            << setw(16) << particles[j].pz << "\t"
            << setw(16) << zeta << "\t"
            << setw(16) << potenziale[0].ey << "\t"
            << setw(16) << I1 << "\t"
            << setw(16) << I2 << "\t"
            << setw(16) << I3 << endl;

        /*
        //UPDATE CAMPI ATTORNO ALLA NUOVA POSIZIONE DELLA PARTICELLA  (al momento non li evolve, il set li pone pari al get precedente)
        gridX = (int) (particles[j].x/data.deltaX);   //il punto di griglia lungo l'asse x piu` vicino all'origine
        gridY = (int) (particles[j].y/data.deltaY);   //il punto di griglia lungo l'asse y piu` vicino all'origine
        gridZ = (int) (particles[j].z/data.deltaZ);   //il punto di griglia lungo l'asse z piu` vicino all'origine

        campoSuPuntiGriglia.at(gridZ * data.nGridPointsX * data.nGridPointsY + gridY * data.nGridPointsX + gridX).setEx( campoSuPuntiGriglia.at(gridZ * data.nGridPointsX * data.nGridPointsY + gridY * data.nGridPointsX + gridX).ex );       // primo punto
        campoSuPuntiGriglia.at(gridZ * data.nGridPointsX * data.nGridPointsY + gridY * data.nGridPointsX + gridX).setEy( campoSuPuntiGriglia.at(gridZ * data.nGridPointsX * data.nGridPointsY + gridY * data.nGridPointsX + gridX).ey );
        campoSuPuntiGriglia.at(gridZ * data.nGridPointsX * data.nGridPointsY + gridY * data.nGridPointsX + gridX).setEz( campoSuPuntiGriglia.at(gridZ * data.nGridPointsX * data.nGridPointsY + gridY * data.nGridPointsX + gridX).ez );
        campoSuPuntiGriglia.at(gridZ * data.nGridPointsX * data.nGridPointsY + gridY * data.nGridPointsX + gridX).setBx( campoSuPuntiGriglia.at(gridZ * data.nGridPointsX * data.nGridPointsY + gridY * data.nGridPointsX + gridX).bx );
        campoSuPuntiGriglia.at(gridZ * data.nGridPointsX * data.nGridPointsY + gridY * data.nGridPointsX + gridX).setBy( campoSuPuntiGriglia.at(gridZ * data.nGridPointsX * data.nGridPointsY + gridY * data.nGridPointsX + gridX).by );
        campoSuPuntiGriglia.at(gridZ * data.nGridPointsX * data.nGridPointsY + gridY * data.nGridPointsX + gridX).setBz( campoSuPuntiGriglia.at(gridZ * data.nGridPointsX * data.nGridPointsY + gridY * data.nGridPointsX + gridX).bz );

        campoSuPuntiGriglia.at(gridZ * data.nGridPointsX * data.nGridPointsY + gridY * data.nGridPointsX + (gridX+1)).setEx( campoSuPuntiGriglia.at(gridZ * data.nGridPointsX * data.nGridPointsY + gridY * data.nGridPointsX + (gridX+1)).ex );   // secondo punto
        campoSuPuntiGriglia.at(gridZ * data.nGridPointsX * data.nGridPointsY + gridY * data.nGridPointsX + (gridX+1)).setEy( campoSuPuntiGriglia.at(gridZ * data.nGridPointsX * data.nGridPointsY + gridY * data.nGridPointsX + (gridX+1)).ey );
        campoSuPuntiGriglia.at(gridZ * data.nGridPointsX * data.nGridPointsY + gridY * data.nGridPointsX + (gridX+1)).setEz( campoSuPuntiGriglia.at(gridZ * data.nGridPointsX * data.nGridPointsY + gridY * data.nGridPointsX + (gridX+1)).ez );
        campoSuPuntiGriglia.at(gridZ * data.nGridPointsX * data.nGridPointsY + gridY * data.nGridPointsX + (gridX+1)).setBx( campoSuPuntiGriglia.at(gridZ * data.nGridPointsX * data.nGridPointsY + gridY * data.nGridPointsX + (gridX+1)).bx );
        campoSuPuntiGriglia.at(gridZ * data.nGridPointsX * data.nGridPointsY + gridY * data.nGridPointsX + (gridX+1)).setBy( campoSuPuntiGriglia.at(gridZ * data.nGridPointsX * data.nGridPointsY + gridY * data.nGridPointsX + (gridX+1)).by );
        campoSuPuntiGriglia.at(gridZ * data.nGridPointsX * data.nGridPointsY + gridY * data.nGridPointsX + (gridX+1)).setBz( campoSuPuntiGriglia.at(gridZ * data.nGridPointsX * data.nGridPointsY + gridY * data.nGridPointsX + (gridX+1)).bz );

        campoSuPuntiGriglia.at(gridZ * data.nGridPointsX * data.nGridPointsY + (gridY+1) * data.nGridPointsX + gridX).setEx( campoSuPuntiGriglia.at(gridZ * data.nGridPointsX * data.nGridPointsY + (gridY+1) * data.nGridPointsX + gridX).ex );   // terzo punto
        campoSuPuntiGriglia.at(gridZ * data.nGridPointsX * data.nGridPointsY + (gridY+1) * data.nGridPointsX + gridX).setEy( campoSuPuntiGriglia.at(gridZ * data.nGridPointsX * data.nGridPointsY + (gridY+1) * data.nGridPointsX + gridX).ey );
        campoSuPuntiGriglia.at(gridZ * data.nGridPointsX * data.nGridPointsY + (gridY+1) * data.nGridPointsX + gridX).setEz( campoSuPuntiGriglia.at(gridZ * data.nGridPointsX * data.nGridPointsY + (gridY+1) * data.nGridPointsX + gridX).ez );
        campoSuPuntiGriglia.at(gridZ * data.nGridPointsX * data.nGridPointsY + (gridY+1) * data.nGridPointsX + gridX).setBx( campoSuPuntiGriglia.at(gridZ * data.nGridPointsX * data.nGridPointsY + (gridY+1) * data.nGridPointsX + gridX).bx );
        campoSuPuntiGriglia.at(gridZ * data.nGridPointsX * data.nGridPointsY + (gridY+1) * data.nGridPointsX + gridX).setBy( campoSuPuntiGriglia.at(gridZ * data.nGridPointsX * data.nGridPointsY + (gridY+1) * data.nGridPointsX + gridX).by );
        campoSuPuntiGriglia.at(gridZ * data.nGridPointsX * data.nGridPointsY + (gridY+1) * data.nGridPointsX + gridX).setBz( campoSuPuntiGriglia.at(gridZ * data.nGridPointsX * data.nGridPointsY + (gridY+1) * data.nGridPointsX + gridX).bz );

        campoSuPuntiGriglia.at(gridZ * data.nGridPointsX * data.nGridPointsY + (gridY+1) * data.nGridPointsX + (gridX+1)).setEx( campoSuPuntiGriglia.at(gridZ * data.nGridPointsX * data.nGridPointsY + (gridY+1) * data.nGridPointsX + (gridX+1)).ex );   // quarto punto
        campoSuPuntiGriglia.at(gridZ * data.nGridPointsX * data.nGridPointsY + (gridY+1) * data.nGridPointsX + (gridX+1)).setEy( campoSuPuntiGriglia.at(gridZ * data.nGridPointsX * data.nGridPointsY + (gridY+1) * data.nGridPointsX + (gridX+1)).ey );
        campoSuPuntiGriglia.at(gridZ * data.nGridPointsX * data.nGridPointsY + (gridY+1) * data.nGridPointsX + (gridX+1)).setEz( campoSuPuntiGriglia.at(gridZ * data.nGridPointsX * data.nGridPointsY + (gridY+1) * data.nGridPointsX + (gridX+1)).ez );
        campoSuPuntiGriglia.at(gridZ * data.nGridPointsX * data.nGridPointsY + (gridY+1) * data.nGridPointsX + (gridX+1)).setBx( campoSuPuntiGriglia.at(gridZ * data.nGridPointsX * data.nGridPointsY + (gridY+1) * data.nGridPointsX + (gridX+1)).bx );
        campoSuPuntiGriglia.at(gridZ * data.nGridPointsX * data.nGridPointsY + (gridY+1) * data.nGridPointsX + (gridX+1)).setBy( campoSuPuntiGriglia.at(gridZ * data.nGridPointsX * data.nGridPointsY + (gridY+1) * data.nGridPointsX + (gridX+1)).by );
        campoSuPuntiGriglia.at(gridZ * data.nGridPointsX * data.nGridPointsY + (gridY+1) * data.nGridPointsX + (gridX+1)).setBz( campoSuPuntiGriglia.at(gridZ * data.nGridPointsX * data.nGridPointsY + (gridY+1) * data.nGridPointsX + (gridX+1)).bz );

        campoSuPuntiGriglia.at((gridZ+1) * data.nGridPointsX * data.nGridPointsY + gridY * data.nGridPointsX + gridX).setEx( campoSuPuntiGriglia.at((gridZ+1) * data.nGridPointsX * data.nGridPointsY + gridY * data.nGridPointsX + gridX).ex );   // quinto punto
        campoSuPuntiGriglia.at((gridZ+1) * data.nGridPointsX * data.nGridPointsY + gridY * data.nGridPointsX + gridX).setEy( campoSuPuntiGriglia.at((gridZ+1) * data.nGridPointsX * data.nGridPointsY + gridY * data.nGridPointsX + gridX).ey );
        campoSuPuntiGriglia.at((gridZ+1) * data.nGridPointsX * data.nGridPointsY + gridY * data.nGridPointsX + gridX).setEz( campoSuPuntiGriglia.at((gridZ+1) * data.nGridPointsX * data.nGridPointsY + gridY * data.nGridPointsX + gridX).ez );
        campoSuPuntiGriglia.at((gridZ+1) * data.nGridPointsX * data.nGridPointsY + gridY * data.nGridPointsX + gridX).setBx( campoSuPuntiGriglia.at((gridZ+1) * data.nGridPointsX * data.nGridPointsY + gridY * data.nGridPointsX + gridX).bx );
        campoSuPuntiGriglia.at((gridZ+1) * data.nGridPointsX * data.nGridPointsY + gridY * data.nGridPointsX + gridX).setBy( campoSuPuntiGriglia.at((gridZ+1) * data.nGridPointsX * data.nGridPointsY + gridY * data.nGridPointsX + gridX).by );
        campoSuPuntiGriglia.at((gridZ+1) * data.nGridPointsX * data.nGridPointsY + gridY * data.nGridPointsX + gridX).setBz( campoSuPuntiGriglia.at((gridZ+1) * data.nGridPointsX * data.nGridPointsY + gridY * data.nGridPointsX + gridX).bz );

        campoSuPuntiGriglia.at((gridZ+1) * data.nGridPointsX * data.nGridPointsY + gridY * data.nGridPointsX + (gridX+1)).setEx( campoSuPuntiGriglia.at((gridZ+1) * data.nGridPointsX * data.nGridPointsY + gridY * data.nGridPointsX + (gridX+1)).ex );   // sesto punto
        campoSuPuntiGriglia.at((gridZ+1) * data.nGridPointsX * data.nGridPointsY + gridY * data.nGridPointsX + (gridX+1)).setEy( campoSuPuntiGriglia.at((gridZ+1) * data.nGridPointsX * data.nGridPointsY + gridY * data.nGridPointsX + (gridX+1)).ey );
        campoSuPuntiGriglia.at((gridZ+1) * data.nGridPointsX * data.nGridPointsY + gridY * data.nGridPointsX + (gridX+1)).setEz( campoSuPuntiGriglia.at((gridZ+1) * data.nGridPointsX * data.nGridPointsY + gridY * data.nGridPointsX + (gridX+1)).ez );
        campoSuPuntiGriglia.at((gridZ+1) * data.nGridPointsX * data.nGridPointsY + gridY * data.nGridPointsX + (gridX+1)).setBx( campoSuPuntiGriglia.at((gridZ+1) * data.nGridPointsX * data.nGridPointsY + gridY * data.nGridPointsX + (gridX+1)).bx );
        campoSuPuntiGriglia.at((gridZ+1) * data.nGridPointsX * data.nGridPointsY + gridY * data.nGridPointsX + (gridX+1)).setBy( campoSuPuntiGriglia.at((gridZ+1) * data.nGridPointsX * data.nGridPointsY + gridY * data.nGridPointsX + (gridX+1)).by );
        campoSuPuntiGriglia.at((gridZ+1) * data.nGridPointsX * data.nGridPointsY + gridY * data.nGridPointsX + (gridX+1)).setBz( campoSuPuntiGriglia.at((gridZ+1) * data.nGridPointsX * data.nGridPointsY + gridY * data.nGridPointsX + (gridX+1)).bz );

        campoSuPuntiGriglia.at((gridZ+1) * data.nGridPointsX * data.nGridPointsY + (gridY+1) * data.nGridPointsX + gridX).setEx( campoSuPuntiGriglia.at((gridZ+1) * data.nGridPointsX * data.nGridPointsY + (gridY+1) * data.nGridPointsX + gridX).ex );   // settimo punto
        campoSuPuntiGriglia.at((gridZ+1) * data.nGridPointsX * data.nGridPointsY + (gridY+1) * data.nGridPointsX + gridX).setEy( campoSuPuntiGriglia.at((gridZ+1) * data.nGridPointsX * data.nGridPointsY + (gridY+1) * data.nGridPointsX + gridX).ey );
        campoSuPuntiGriglia.at((gridZ+1) * data.nGridPointsX * data.nGridPointsY + (gridY+1) * data.nGridPointsX + gridX).setEz( campoSuPuntiGriglia.at((gridZ+1) * data.nGridPointsX * data.nGridPointsY + (gridY+1) * data.nGridPointsX + gridX).ez );
        campoSuPuntiGriglia.at((gridZ+1) * data.nGridPointsX * data.nGridPointsY + (gridY+1) * data.nGridPointsX + gridX).setBx( campoSuPuntiGriglia.at((gridZ+1) * data.nGridPointsX * data.nGridPointsY + (gridY+1) * data.nGridPointsX + gridX).bx );
        campoSuPuntiGriglia.at((gridZ+1) * data.nGridPointsX * data.nGridPointsY + (gridY+1) * data.nGridPointsX + gridX).setBy( campoSuPuntiGriglia.at((gridZ+1) * data.nGridPointsX * data.nGridPointsY + (gridY+1) * data.nGridPointsX + gridX).by );
        campoSuPuntiGriglia.at((gridZ+1) * data.nGridPointsX * data.nGridPointsY + (gridY+1) * data.nGridPointsX + gridX).setBz( campoSuPuntiGriglia.at((gridZ+1) * data.nGridPointsX * data.nGridPointsY + (gridY+1) * data.nGridPointsX + gridX).bz );

        campoSuPuntiGriglia.at((gridZ+1) * data.nGridPointsX * data.nGridPointsY + (gridY+1) * data.nGridPointsX + (gridX+1)).setEx( campoSuPuntiGriglia.at((gridZ+1) * data.nGridPointsX * data.nGridPointsY + (gridY+1) * data.nGridPointsX + (gridX+1)).ex );   // ottavo punto
        campoSuPuntiGriglia.at((gridZ+1) * data.nGridPointsX * data.nGridPointsY + (gridY+1) * data.nGridPointsX + (gridX+1)).setEy( campoSuPuntiGriglia.at((gridZ+1) * data.nGridPointsX * data.nGridPointsY + (gridY+1) * data.nGridPointsX + (gridX+1)).ex );
        campoSuPuntiGriglia.at((gridZ+1) * data.nGridPointsX * data.nGridPointsY + (gridY+1) * data.nGridPointsX + (gridX+1)).setEz( campoSuPuntiGriglia.at((gridZ+1) * data.nGridPointsX * data.nGridPointsY + (gridY+1) * data.nGridPointsX + (gridX+1)).ex );
        campoSuPuntiGriglia.at((gridZ+1) * data.nGridPointsX * data.nGridPointsY + (gridY+1) * data.nGridPointsX + (gridX+1)).setBx( campoSuPuntiGriglia.at((gridZ+1) * data.nGridPointsX * data.nGridPointsY + (gridY+1) * data.nGridPointsX + (gridX+1)).ex );
        campoSuPuntiGriglia.at((gridZ+1) * data.nGridPointsX * data.nGridPointsY + (gridY+1) * data.nGridPointsX + (gridX+1)).setBy( campoSuPuntiGriglia.at((gridZ+1) * data.nGridPointsX * data.nGridPointsY + (gridY+1) * data.nGridPointsX + (gridX+1)).ex );
        campoSuPuntiGriglia.at((gridZ+1) * data.nGridPointsX * data.nGridPointsY + (gridY+1) * data.nGridPointsX + (gridX+1)).setBz( campoSuPuntiGriglia.at((gridZ+1) * data.nGridPointsX * data.nGridPointsY + (gridY+1) * data.nGridPointsX + (gridX+1)).ex );
        */
      }
    }
  }
}

void evolveRK4(Data data, Particle particle, Field fieldOnParticle, vector<Particle>& dx2, vector<Particle>& x2)
{
  Particle dx, dx2a, x2a;
  double gamma = sqrt(1 + pow(particle.px, 2) + pow(particle.py, 2) + pow(particle.pz, 2));

  dx.setParticleX(particle.px / (gamma * particle.m));
  dx.setParticleY(particle.py / (gamma * particle.m));
  dx.setParticleZ(particle.pz / (gamma * particle.m));
  dx.setParticlePX(particle.q * fieldOnParticle.ex + particle.q * dx.y * fieldOnParticle.bz - particle.q * dx.z * fieldOnParticle.by);
  dx.setParticlePY(particle.q * fieldOnParticle.ey + particle.q * dx.z * fieldOnParticle.bx - particle.q * dx.x * fieldOnParticle.bz);
  dx.setParticlePZ(particle.q * fieldOnParticle.ez + particle.q * dx.x * fieldOnParticle.by - particle.q * dx.y * fieldOnParticle.bx);

  //RIEMPIE dx2 per il metodo RK4 con il dx del primo passo
  dx2a = dx;
  dx2.push_back(dx2a);

  //RIEMPIE x2 con le posizioni avanzate di particles[i] dopo il primo passo
  x2a.setParticleX(particle.x + dx.x * (data.dt * 0.5));
  x2a.setParticleY(particle.y + dx.y * (data.dt * 0.5));
  x2a.setParticleZ(particle.z + dx.z * (data.dt * 0.5));
  x2a.setParticlePX(particle.px + dx.px * (data.dt * 0.5));
  x2a.setParticlePY(particle.py + dx.py * (data.dt * 0.5));
  x2a.setParticlePZ(particle.pz + dx.pz * (data.dt * 0.5));

  x2.push_back(x2a);
}

void evolveLPF(Data data, Particle& particle, double dtp, double dtx)
{
  double gamma, gamma_inv, b2;
  Particle pstar, pHalfAdvanced;

  Field b;
  Field* campoSuPunto;
  campoSuPunto = new Field[1];

  gamma = sqrt(1. + pow(particle.px, 2) + pow(particle.py, 2) + pow(particle.pz, 2));
  gamma_inv = 1.0 / gamma;

  particle.setParticleX(particle.x + dtx * particle.px * gamma_inv);
  particle.setParticleY(particle.y + dtx * particle.py * gamma_inv);
  particle.setParticleZ(particle.z + dtx * particle.pz * gamma_inv);

  data.aumentaT(dtx);

  calcolaCampoAnalitico1DSuParticella(data, particle, *campoSuPunto);

  //pstar = p(n) + Q * E(n+1/2) * DeltaT/2
  pstar.setParticlePX(particle.px + 0.5 * dtp * particle.q * (*campoSuPunto).ex);
  pstar.setParticlePY(particle.py + 0.5 * dtp * particle.q * (*campoSuPunto).ey);
  pstar.setParticlePZ(particle.pz + 0.5 * dtp * particle.q * (*campoSuPunto).ez);
  gamma = sqrt(1. + pow(pstar.px, 2) + pow(pstar.py, 2) + pow(pstar.pz, 2));
  gamma_inv = 1.0 / gamma;

  //b = (DeltaT/2) * Q * B(n+1/2) / gamma(aggiornato al pstar)
  b.setBx(0.5 * dtp * (*campoSuPunto).bx * particle.q * gamma_inv);
  b.setBy(0.5 * dtp * (*campoSuPunto).by * particle.q * gamma_inv);
  b.setBz(0.5 * dtp * (*campoSuPunto).bz * particle.q * gamma_inv);
  b2 = pow(b.bx, 2) + pow(b.by, 2) + pow(b.bz, 2);
  b2 += 1.0;

  /****************************************************************************
  //p(n+1/2) = [pstar + pstar x b + (pstar x b) * b] / b2   SECONDO TURCHETTI  (ma i conti analitici diconon che e` sbagliato, anche perche' l'ultimo termine e` uno scalare!)
  pHalfAdvanced.setParticlePX( (pstar.px + pstar.py * b.bz - b.by * pstar.pz + (pstar.py * b.bz - b.by * pstar.pz) * b.bx ) / b2 );
  pHalfAdvanced.setParticlePY( (pstar.py + pstar.pz * b.bx - b.bz * pstar.px + (pstar.pz * b.bx - b.bz * pstar.px) * b.by ) / b2 );
  pHalfAdvanced.setParticlePZ( (pstar.pz + pstar.px * b.by - b.bx * pstar.py + (pstar.px * b.by - b.bx * pstar.py) * b.bz ) / b2 );
  ****************************************************************************/

  //  /****************************************************************************
  //p(n+1/2) = [pstar + pstar x b + b x (pstar x b)] / b2   SECONDO LONDRILLO
  pHalfAdvanced.setParticlePX((pstar.px + pstar.py * b.bz - b.by * pstar.pz + b.bx * pstar.px * b.bx) / b2);
  pHalfAdvanced.setParticlePY((pstar.py + pstar.pz * b.bx - b.bz * pstar.px + b.by * pstar.py * b.by) / b2);
  pHalfAdvanced.setParticlePZ((pstar.pz + pstar.px * b.by - b.bx * pstar.py + b.bz * pstar.pz * b.bz) / b2);
  //  ****************************************************************************/

  /****************************************************************************
    //p(n+1/2) = [pstar + pstar x b] / b2   COME IMPLEMENTATO DA TURCHETTI NEL CODICE
    pHalfAdvanced.setParticlePX( (pstar.px + pstar.py * b.bz) / b2 );
    pHalfAdvanced.setParticlePY( (pstar.py - pstar.px * b.bz) / b2 );
    pHalfAdvanced.setParticlePZ( 0. );  //per come ho strutturato i pacchetti d'onda su z non dovrebbe succedere nulla e quindi lo fisso io a zero.
    ****************************************************************************/

  particle.setParticlePX(2.0 * pHalfAdvanced.px - particle.px);
  particle.setParticlePY(2.0 * pHalfAdvanced.py - particle.py);
  particle.setParticlePZ(2.0 * pHalfAdvanced.pz - particle.pz);
}
