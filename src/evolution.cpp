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



#define OUTPUTALL

#include <iostream>
#include <iomanip>
#include <cstdio>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <fstream>
#include <limits> // for declaration of 'numeric_limits'
#include <ios>    // for declaration of 'streamsize'

#include "evolution.h"
#include "filler.h"

using namespace std;

const double C = 1.00;

void evolveRK4_nogrid(Data data, vector<Particle> &particles, vector<Field> &fieldOnParticles, ofstream &outputstream)
{
  Particle x2, dx, dx2;
  double gamma;
  double zeta;
  double I1, I2, I3;

  Field *campoSuPunto;
  campoSuPunto = new Field[1];
  vector<Field> potenziale;

  for (int j = 0; j < data.getNsteps(); j++)
  {
    for (int i = 0; i < data.getNelectrons(); i++)
    {
      gamma = sqrt(1 + pow(particles[i].getParticlePX(), 2) + pow(particles[i].getParticlePY(), 2) + pow(particles[i].getParticlePZ(), 2));
      fieldOnParticles.clear();
      calcolaCampoAnalitico1DSuParticella(data, particles[i], *campoSuPunto);
      fieldOnParticles.push_back(*campoSuPunto);


      // PASSO #1
      // RIEMPIE dx = vx, vy, vz, Fx, Fy, Fz
      dx.setParticleX(particles[i].getParticlePX() / (gamma * particles[i].getParticleM()));
      dx.setParticleY(particles[i].getParticlePY() / (gamma * particles[i].getParticleM()));
      dx.setParticleZ(particles[i].getParticlePZ() / (gamma * particles[i].getParticleM()));
      dx.setParticlePX(particles[i].getParticleQ() * fieldOnParticles[0].getEx() + particles[i].getParticleQ() * dx.getParticleY() * fieldOnParticles[0].getBz() - particles[i].getParticleQ() * dx.getParticleZ() * fieldOnParticles[0].getBy());  //CONTROLLARE I SEGNI!!!!
      dx.setParticlePY(particles[i].getParticleQ() * fieldOnParticles[0].getEy() + particles[i].getParticleQ() * dx.getParticleZ() * fieldOnParticles[0].getBx() - particles[i].getParticleQ() * dx.getParticleX() * fieldOnParticles[0].getBz());
      dx.setParticlePZ(particles[i].getParticleQ() * fieldOnParticles[0].getEz() + particles[i].getParticleQ() * dx.getParticleX() * fieldOnParticles[0].getBy() - particles[i].getParticleQ() * dx.getParticleY() * fieldOnParticles[0].getBx());

      // RIEMPIE dx2 per il metodo RK4 con il dx del primo passo
      dx2 = dx;

      //RIEMPIE x2 con le posizioni avanzate di particles[i] dopo il primo passo
      x2.setParticleX(particles[i].getParticleX() + dx.getParticleX() * (data.getDeltaT()*0.5));      // E LE FORZE???
      x2.setParticleY(particles[i].getParticleY() + dx.getParticleY() * (data.getDeltaT()*0.5));
      x2.setParticleZ(particles[i].getParticleZ() + dx.getParticleZ() * (data.getDeltaT()*0.5));
      x2.setParticlePX(particles[i].getParticlePX() + dx.getParticlePX() *  (data.getDeltaT()*0.5));
      x2.setParticlePY(particles[i].getParticlePY() + dx.getParticlePY() *  (data.getDeltaT()*0.5));
      x2.setParticlePZ(particles[i].getParticlePZ() + dx.getParticlePZ() *  (data.getDeltaT()*0.5));


      // PASSO #2
      // AVANZA TEMPO DI deltaT/2 e calcola campi nella posizione x2
//      fieldOnParticles.clear();
      data.aumentaT(data.getDeltaT() * 0.5);
      calcolaCampoAnalitico1DSuParticella(data, x2, *campoSuPunto);
      fieldOnParticles.push_back(*campoSuPunto);

      gamma = sqrt(1 + pow(x2.getParticlePX(), 2) + pow(x2.getParticlePY(), 2) + pow(x2.getParticlePZ(), 2));

      // AGGIORNA dx con i nuovi valori
      dx.setParticleX(x2.getParticlePX() / (gamma * particles[i].getParticleM()));
      dx.setParticleY(x2.getParticlePY() / (gamma * particles[i].getParticleM()));
      dx.setParticleZ(x2.getParticlePZ() / (gamma * particles[i].getParticleM()));
      dx.setParticlePX(particles[i].getParticleQ() * fieldOnParticles[1].getEx() + particles[i].getParticleQ() * dx.getParticleY() * fieldOnParticles[1].getBz() - particles[i].getParticleQ() * dx.getParticleZ() * fieldOnParticles[1].getBy());  //CONTROLLARE I SEGNI!!!!
      dx.setParticlePY(particles[i].getParticleQ() * fieldOnParticles[1].getEy() + particles[i].getParticleQ() * dx.getParticleZ() * fieldOnParticles[1].getBx() - particles[i].getParticleQ() * dx.getParticleX() * fieldOnParticles[1].getBz());
      dx.setParticlePZ(particles[i].getParticleQ() * fieldOnParticles[1].getEz() + particles[i].getParticleQ() * dx.getParticleX() * fieldOnParticles[1].getBy() - particles[i].getParticleQ() * dx.getParticleY() * fieldOnParticles[1].getBx());


      // AGGIUNGE al dx2 i valori calcolati al secondo passo
      dx2.setParticleX(dx2.getParticleX() + 2 * dx.getParticleX());
      dx2.setParticleY(dx2.getParticleY() + 2 * dx.getParticleY());
      dx2.setParticleZ(dx2.getParticleZ() + 2 * dx.getParticleZ());
      dx2.setParticlePX(dx2.getParticlePX() + 2 * dx.getParticlePX());
      dx2.setParticlePY(dx2.getParticlePY() + 2 * dx.getParticlePY());
      dx2.setParticlePZ(dx2.getParticlePZ() + 2 * dx.getParticlePZ());

      // RIEMPIE x2 con le posizione avanzate di particles[i] dopo il secondo passo  // DA SISTEMARE!!!
      x2.setParticleX(particles[i].getParticleX() + dx.getParticleX() * (data.getDeltaT()*0.5));
      x2.setParticleY(particles[i].getParticleY() + dx.getParticleY() * (data.getDeltaT()*0.5));
      x2.setParticleZ(particles[i].getParticleZ() + dx.getParticleZ() * (data.getDeltaT()*0.5));
      x2.setParticlePX(particles[i].getParticlePX() + dx.getParticlePX() *  (data.getDeltaT()*0.5));
      x2.setParticlePY(particles[i].getParticlePY() + dx.getParticlePY() *  (data.getDeltaT()*0.5));
      x2.setParticlePZ(particles[i].getParticlePZ() + dx.getParticlePZ() *  (data.getDeltaT()*0.5));



      // PASSO #3
      // CALCOLA campi nella posizione x2 appena aggiornata
      //fieldOnParticles.clear();
      calcolaCampoAnalitico1DSuParticella(data, x2, *campoSuPunto);
      fieldOnParticles.push_back(*campoSuPunto);

      gamma = sqrt(1 + pow(x2.getParticlePX(), 2) + pow(x2.getParticlePY(), 2) + pow(x2.getParticlePZ(), 2));

      //AGGIORNA dx con i nuovi valori  DA SISTEMARE!!!!!
      dx.setParticleX(x2.getParticlePX() / (gamma * particles[i].getParticleM()));
      dx.setParticleY(x2.getParticlePY() / (gamma * particles[i].getParticleM()));
      dx.setParticleZ(x2.getParticlePZ() / (gamma * particles[i].getParticleM()));
      dx.setParticlePX(particles[i].getParticleQ() * fieldOnParticles[2].getEx() + particles[i].getParticleQ() * dx.getParticleY() * fieldOnParticles[2].getBz() - particles[i].getParticleQ() * dx.getParticleZ() * fieldOnParticles[2].getBy());  //CONTROLLARE I SEGNI!!!!
      dx.setParticlePY(particles[i].getParticleQ() * fieldOnParticles[2].getEy() + particles[i].getParticleQ() * dx.getParticleZ() * fieldOnParticles[2].getBx() - particles[i].getParticleQ() * dx.getParticleX() * fieldOnParticles[2].getBz());
      dx.setParticlePZ(particles[i].getParticleQ() * fieldOnParticles[2].getEz() + particles[i].getParticleQ() * dx.getParticleX() * fieldOnParticles[2].getBy() - particles[i].getParticleQ() * dx.getParticleY() * fieldOnParticles[2].getBx());

      //AGGIUNGE al dx2 i valori calcolati al terzo passo
      dx2.setParticleX(dx2.getParticleX() + 2 * dx.getParticleX());
      dx2.setParticleY(dx2.getParticleY() + 2 * dx.getParticleY());
      dx2.setParticleZ(dx2.getParticleZ() + 2 * dx.getParticleZ());
      dx2.setParticlePX(dx2.getParticlePX() + 2 * dx.getParticlePX());
      dx2.setParticlePY(dx2.getParticlePY() + 2 * dx.getParticlePY());
      dx2.setParticlePZ(dx2.getParticlePZ() + 2 * dx.getParticlePZ());

      //RIEMPIE x2 con le posizione avanzate di particles[i]      // DA SISTEMARE!!!
      x2.setParticleX(particles[i].getParticleX() + dx.getParticleX() * data.getDeltaT());
      x2.setParticleY(particles[i].getParticleY() + dx.getParticleY() * data.getDeltaT());
      x2.setParticleZ(particles[i].getParticleZ() + dx.getParticleZ() * data.getDeltaT());
      x2.setParticlePX(particles[i].getParticlePX() + dx.getParticlePX() *  data.getDeltaT());
      x2.setParticlePY(particles[i].getParticlePY() + dx.getParticlePY() *  data.getDeltaT());
      x2.setParticlePZ(particles[i].getParticlePZ() + dx.getParticlePZ() *  data.getDeltaT());



      // PASSO #4
      //AVANZA TEMPO DI deltaT/2 e calcola campi nella posizione x2
      //fieldOnParticles.clear();
      data.aumentaT(data.getDeltaT() * 0.5);
      calcolaCampoAnalitico1DSuParticella(data, x2, *campoSuPunto);
      fieldOnParticles.push_back(*campoSuPunto);
      gamma = sqrt(1 + pow(x2.getParticlePX(), 2) + pow(x2.getParticlePY(), 2) + pow(x2.getParticlePZ(), 2));

      //AGGIORNA dx con i nuovi valori  DA SISTEMARE!!!!!
      dx.setParticleX(x2.getParticlePX() / (gamma * particles[i].getParticleM()));
      dx.setParticleY(x2.getParticlePY() / (gamma * particles[i].getParticleM()));
      dx.setParticleZ(x2.getParticlePZ() / (gamma * particles[i].getParticleM()));
      dx.setParticlePX(particles[i].getParticleQ() * fieldOnParticles[3].getEx() + particles[i].getParticleQ() * dx.getParticleY() * fieldOnParticles[3].getBz() - particles[i].getParticleQ() * dx.getParticleZ() * fieldOnParticles[3].getBy());  //CONTROLLARE I SEGNI!!!!
      dx.setParticlePY(particles[i].getParticleQ() * fieldOnParticles[3].getEy() + particles[i].getParticleQ() * dx.getParticleZ() * fieldOnParticles[3].getBx() - particles[i].getParticleQ() * dx.getParticleX() * fieldOnParticles[3].getBz());
      dx.setParticlePZ(particles[i].getParticleQ() * fieldOnParticles[3].getEz() + particles[i].getParticleQ() * dx.getParticleX() * fieldOnParticles[3].getBy() - particles[i].getParticleQ() * dx.getParticleY() * fieldOnParticles[3].getBx());


      //AGGIUNGE al dx2 i valori calcolati al terzo passo
      dx2.setParticleX(dx2.getParticleX() + dx.getParticleX());
      dx2.setParticleY(dx2.getParticleY() + dx.getParticleY());
      dx2.setParticleZ(dx2.getParticleZ() + dx.getParticleZ());
      dx2.setParticlePX(dx2.getParticlePX() + dx.getParticlePX());
      dx2.setParticlePY(dx2.getParticlePY() + dx.getParticlePY());
      dx2.setParticlePZ(dx2.getParticlePZ() + dx.getParticlePZ());


      particles[i].setParticleX(particles[i].getParticleX() + (data.getDeltaT() / 6.) * dx2.getParticleX());
      particles[i].setParticleY(particles[i].getParticleY() + (data.getDeltaT() / 6.) * dx2.getParticleY());
      particles[i].setParticleZ(particles[i].getParticleZ() + (data.getDeltaT() / 6.) * dx2.getParticleZ());
      particles[i].setParticlePX(particles[i].getParticlePX() + (data.getDeltaT() / 6.) * dx2.getParticlePX());
      particles[i].setParticlePY(particles[i].getParticlePY() + (data.getDeltaT() / 6.) * dx2.getParticlePY());
      particles[i].setParticlePZ(particles[i].getParticlePZ() + (data.getDeltaT() / 6.) * dx2.getParticlePZ());

      particles[i].setParticleCell(data);

      data.aumentaT(-data.getDeltaT());       // faccio tornare il tempo indietro cosi` gli altri elettroni di questo step non hanno il tempo gia` avanzato
                                              // il tempo avanza definitivamente in ciascun ciclo 'for' esterno (sugli nsteps)

      // INTEGRALI PRIMI
      potenziale.clear();
      zeta = particles[i].getParticleX() - C * data.getT();             // non usato nella mia struttura delle equazioni di campo
      gamma = sqrt(1 + pow(particles[i].getParticlePX(), 2) + pow(particles[i].getParticlePY(), 2) + pow(particles[i].getParticlePZ(), 2));
      calcolaPotenzialeAnalitico1DSuParticella(data, particles[i], *campoSuPunto);
      potenziale.push_back(*campoSuPunto);

      I1 = (particles[i].getParticlePY() / (particles[i].getParticleM()*C)) + potenziale[0].getEy();
      I2 = particles[i].getParticlePZ() / (particles[i].getParticleM()*C);
      I3 = gamma - particles[i].getParticlePX() / (particles[i].getParticleM()*C);

      potenziale.clear();

      outputstream
        << setw(12) << i + 1 << "\t"
        << setprecision(14) << setw(16) << particles[i].getParticleX() << "\t"
        << setw(16) << particles[i].getParticleY() << "\t"
        << setw(16) << particles[i].getParticleZ() << "\t"
        //        << setw(16) << particles[i].getParticleCell() << "\t"
        << setw(16) << particles[i].getParticlePX() << "\t"
        << setw(16) << particles[i].getParticlePY() << "\t"
        << setw(16) << particles[i].getParticlePZ() << "\t"
        << setw(16) << zeta << "\t"
        << setw(16) << I1 << "\t"
        << setw(16) << I2 << "\t"
        << setw(16) << I3 << endl;
    }
    data.aumentaT(data.getDeltaT());
  }
}


void evolveRK4_withgrid_onthefly(Data data, vector<Particle> &particles, vector<Field> &fieldOnParticles, ofstream &outputstream)
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

  Field *FOG;
  FOG = new Field[1];
  Field *campoSuPunto;
  campoSuPunto = new Field[1];

  vector<Field> fieldOnGrid;
  vector<Field> potenziale;

  for (int j = 0; j < data.getNsteps(); j++)
  {
    data.aumentaT(data.getDeltaT());
    for (int i = 0; i < data.getNelectrons(); i++)
    {

      fieldOnParticles.clear();       //SVUOTA il vettore dei campi sulle particelle ad ogni ciclo

      gridXX = particles[i].getParticleX() / data.getDeltaX();
      if (gridXX < 0) gridX = ((int)gridXX) - 1;           //il punto di griglia lungo l'asse x a sinistra della particella
      else gridX = (int)gridXX;                      //

      gridYY = particles[i].getParticleY() / data.getDeltaY();
      if (gridYY < 0) gridY = ((int)gridYY) - 1;           //il punto di griglia lungo l'asse y a sinistra della particella
      else gridY = (int)gridYY;                      //

      gridZZ = particles[i].getParticleZ() / data.getDeltaZ();
      if (gridZZ < 0) gridZ = ((int)gridZZ) - 1;           //il punto di griglia lungo l'asse z a sinistra della particella
      else gridZ = (int)gridZZ;                      //

      coordinatePrimoPuntoGriglia.setParticleX(gridX*data.getDeltaX());
      coordinatePrimoPuntoGriglia.setParticleY(gridY*data.getDeltaY());
      coordinatePrimoPuntoGriglia.setParticleZ(gridZ*data.getDeltaZ());
      coordinateSecondoPuntoGriglia.setParticleX(gridX*data.getDeltaX() + data.getDeltaX());
      coordinateSecondoPuntoGriglia.setParticleY(gridY*data.getDeltaY());
      coordinateSecondoPuntoGriglia.setParticleZ(gridZ*data.getDeltaZ());
      coordinateTerzoPuntoGriglia.setParticleX(gridX*data.getDeltaX());
      coordinateTerzoPuntoGriglia.setParticleY(gridY*data.getDeltaY() + data.getDeltaY());
      coordinateTerzoPuntoGriglia.setParticleZ(gridZ*data.getDeltaZ());
      coordinateQuartoPuntoGriglia.setParticleX(gridX*data.getDeltaX() + data.getDeltaX());
      coordinateQuartoPuntoGriglia.setParticleY(gridY*data.getDeltaY() + data.getDeltaY());
      coordinateQuartoPuntoGriglia.setParticleZ(gridZ*data.getDeltaZ());
      coordinateQuintoPuntoGriglia.setParticleX(gridX*data.getDeltaX());
      coordinateQuintoPuntoGriglia.setParticleY(gridY*data.getDeltaY());
      coordinateQuintoPuntoGriglia.setParticleZ(gridZ*data.getDeltaZ() + data.getDeltaZ());
      coordinateSestoPuntoGriglia.setParticleX(gridX*data.getDeltaX() + data.getDeltaX());
      coordinateSestoPuntoGriglia.setParticleY(gridY*data.getDeltaY());
      coordinateSestoPuntoGriglia.setParticleZ(gridZ*data.getDeltaZ() + data.getDeltaZ());
      coordinateSettimoPuntoGriglia.setParticleX(gridX*data.getDeltaX());
      coordinateSettimoPuntoGriglia.setParticleY(gridY*data.getDeltaY() + data.getDeltaY());
      coordinateSettimoPuntoGriglia.setParticleZ(gridZ*data.getDeltaZ() + data.getDeltaZ());
      coordinateOttavoPuntoGriglia.setParticleX(gridX*data.getDeltaX() + data.getDeltaX());
      coordinateOttavoPuntoGriglia.setParticleY(gridY*data.getDeltaY() + data.getDeltaY());
      coordinateOttavoPuntoGriglia.setParticleZ(gridZ*data.getDeltaZ() + data.getDeltaZ());

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


      gamma = sqrt(1 + pow(particles[i].getParticlePX(), 2) + pow(particles[i].getParticlePY(), 2) + pow(particles[i].getParticlePZ(), 2));


      // PASSO #1
      // RIEMPIE dx = vx, vy, vz, Fx, Fy, Fz
      dx.setParticleX(particles[i].getParticlePX() / (gamma * particles[i].getParticleM()));
      dx.setParticleY(particles[i].getParticlePY() / (gamma * particles[i].getParticleM()));
      dx.setParticleZ(particles[i].getParticlePZ() / (gamma * particles[i].getParticleM()));
      dx.setParticlePX(particles[i].getParticleQ() * fieldOnParticles[0].getEx() + particles[i].getParticleQ() * dx.getParticleY() * fieldOnParticles[0].getBz() - particles[i].getParticleQ() * dx.getParticleZ() * fieldOnParticles[0].getBy());  //CONTROLLARE I SEGNI!!!!
      dx.setParticlePY(particles[i].getParticleQ() * fieldOnParticles[0].getEy() + particles[i].getParticleQ() * dx.getParticleZ() * fieldOnParticles[0].getBx() - particles[i].getParticleQ() * dx.getParticleX() * fieldOnParticles[0].getBz());
      dx.setParticlePZ(particles[i].getParticleQ() * fieldOnParticles[0].getEz() + particles[i].getParticleQ() * dx.getParticleX() * fieldOnParticles[0].getBy() - particles[i].getParticleQ() * dx.getParticleY() * fieldOnParticles[0].getBx());

      // RIEMPIE dx2 per il metodo RK4 con il dx del primo passo
      dx2 = dx;

      // RIEMPIE x2 con le posizioni avanzate di particles[i] dopo il primo passo
      x2.setParticleX(particles[i].getParticleX() + dx.getParticleX() * (data.getDeltaT()*0.5));      // E LE FORZE???
      x2.setParticleY(particles[i].getParticleY() + dx.getParticleY() * (data.getDeltaT()*0.5));
      x2.setParticleZ(particles[i].getParticleZ() + dx.getParticleZ() * (data.getDeltaT()*0.5));
      x2.setParticlePX(particles[i].getParticlePX() + dx.getParticlePX() *  (data.getDeltaT()*0.5));
      x2.setParticlePY(particles[i].getParticlePY() + dx.getParticlePY() *  (data.getDeltaT()*0.5));
      x2.setParticlePZ(particles[i].getParticlePZ() + dx.getParticlePZ() *  (data.getDeltaT()*0.5));


      // PASSO #2
      // AVANZA TEMPO DI deltaT/2 e calcola campi nella posizione x2
      gridXX = x2.getParticleX() / data.getDeltaX();
      if (gridXX < 0) gridX = ((int)gridXX) - 1;           //il punto di griglia lungo l'asse x a sinistra della particella
      else gridX = (int)gridXX;

      gridYY = x2.getParticleY() / data.getDeltaY();
      if (gridYY < 0) gridY = ((int)gridYY) - 1;           //il punto di griglia lungo l'asse y a sinistra della particella
      else gridY = (int)gridYY;

      gridZZ = x2.getParticleZ() / data.getDeltaZ();
      if (gridZZ < 0) gridZ = ((int)gridZZ) - 1;           //il punto di griglia lungo l'asse z a sinistra della particella
      else gridZ = (int)gridZZ;

      coordinatePrimoPuntoGriglia.setParticleX(gridX*data.getDeltaX());
      coordinatePrimoPuntoGriglia.setParticleY(gridY*data.getDeltaY());
      coordinatePrimoPuntoGriglia.setParticleZ(gridZ*data.getDeltaZ());
      coordinateSecondoPuntoGriglia.setParticleX(gridX*data.getDeltaX() + data.getDeltaX());
      coordinateSecondoPuntoGriglia.setParticleY(gridY*data.getDeltaY());
      coordinateSecondoPuntoGriglia.setParticleZ(gridZ*data.getDeltaZ());
      coordinateTerzoPuntoGriglia.setParticleX(gridX*data.getDeltaX());
      coordinateTerzoPuntoGriglia.setParticleY(gridY*data.getDeltaY() + data.getDeltaY());
      coordinateTerzoPuntoGriglia.setParticleZ(gridZ*data.getDeltaZ());
      coordinateQuartoPuntoGriglia.setParticleX(gridX*data.getDeltaX() + data.getDeltaX());
      coordinateQuartoPuntoGriglia.setParticleY(gridY*data.getDeltaY() + data.getDeltaY());
      coordinateQuartoPuntoGriglia.setParticleZ(gridZ*data.getDeltaZ());
      coordinateQuintoPuntoGriglia.setParticleX(gridX*data.getDeltaX());
      coordinateQuintoPuntoGriglia.setParticleY(gridY*data.getDeltaY());
      coordinateQuintoPuntoGriglia.setParticleZ(gridZ*data.getDeltaZ() + data.getDeltaZ());
      coordinateSestoPuntoGriglia.setParticleX(gridX*data.getDeltaX() + data.getDeltaX());
      coordinateSestoPuntoGriglia.setParticleY(gridY*data.getDeltaY());
      coordinateSestoPuntoGriglia.setParticleZ(gridZ*data.getDeltaZ() + data.getDeltaZ());
      coordinateSettimoPuntoGriglia.setParticleX(gridX*data.getDeltaX());
      coordinateSettimoPuntoGriglia.setParticleY(gridY*data.getDeltaY() + data.getDeltaY());
      coordinateSettimoPuntoGriglia.setParticleZ(gridZ*data.getDeltaZ() + data.getDeltaZ());
      coordinateOttavoPuntoGriglia.setParticleX(gridX*data.getDeltaX() + data.getDeltaX());
      coordinateOttavoPuntoGriglia.setParticleY(gridY*data.getDeltaY() + data.getDeltaY());
      coordinateOttavoPuntoGriglia.setParticleZ(gridZ*data.getDeltaZ() + data.getDeltaZ());

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

      data.aumentaT(data.getDeltaT() * 0.5);
      FOP = interpolation2D_linear(data, x2, fieldOnGrid);
      fieldOnParticles.push_back(FOP);

      fieldOnGrid.clear();

      gamma = sqrt(1 + pow(x2.getParticlePX(), 2) + pow(x2.getParticlePY(), 2) + pow(x2.getParticlePZ(), 2));

      // AGGIORNA dx con i nuovi valori
      dx.setParticleX(x2.getParticlePX() / (gamma * particles[i].getParticleM()));
      dx.setParticleY(x2.getParticlePY() / (gamma * particles[i].getParticleM()));
      dx.setParticleZ(x2.getParticlePZ() / (gamma * particles[i].getParticleM()));
      dx.setParticlePX(particles[i].getParticleQ() * fieldOnParticles[1].getEx() + particles[i].getParticleQ() * dx.getParticleY() * fieldOnParticles[1].getBz() - particles[i].getParticleQ() * dx.getParticleZ() * fieldOnParticles[1].getBy());  //CONTROLLARE I SEGNI!!!!
      dx.setParticlePY(particles[i].getParticleQ() * fieldOnParticles[1].getEy() + particles[i].getParticleQ() * dx.getParticleZ() * fieldOnParticles[1].getBx() - particles[i].getParticleQ() * dx.getParticleX() * fieldOnParticles[1].getBz());
      dx.setParticlePZ(particles[i].getParticleQ() * fieldOnParticles[1].getEz() + particles[i].getParticleQ() * dx.getParticleX() * fieldOnParticles[1].getBy() - particles[i].getParticleQ() * dx.getParticleY() * fieldOnParticles[1].getBx());


      // AGGIUNGE al dx2 i valori calcolati al secondo passo
      dx2.setParticleX(dx2.getParticleX() + 2 * dx.getParticleX());
      dx2.setParticleY(dx2.getParticleY() + 2 * dx.getParticleY());
      dx2.setParticleZ(dx2.getParticleZ() + 2 * dx.getParticleZ());
      dx2.setParticlePX(dx2.getParticlePX() + 2 * dx.getParticlePX());
      dx2.setParticlePY(dx2.getParticlePY() + 2 * dx.getParticlePY());
      dx2.setParticlePZ(dx2.getParticlePZ() + 2 * dx.getParticlePZ());

      // RIEMPIE x2 con le posizione avanzate di particles[i] dopo il secondo passo  // DA SISTEMARE!!!
      x2.setParticleX(particles[i].getParticleX() + dx.getParticleX() * (data.getDeltaT()*0.5));
      x2.setParticleY(particles[i].getParticleY() + dx.getParticleY() * (data.getDeltaT()*0.5));
      x2.setParticleZ(particles[i].getParticleZ() + dx.getParticleZ() * (data.getDeltaT()*0.5));
      x2.setParticlePX(particles[i].getParticlePX() + dx.getParticlePX() *  (data.getDeltaT()*0.5));
      x2.setParticlePY(particles[i].getParticlePY() + dx.getParticlePY() *  (data.getDeltaT()*0.5));
      x2.setParticlePZ(particles[i].getParticlePZ() + dx.getParticlePZ() *  (data.getDeltaT()*0.5));



      // PASSO #3
      // CALCOLA campi nella posizione x2 appena aggiornata
      gridXX = x2.getParticleX() / data.getDeltaX();
      if (gridXX < 0) gridX = ((int)gridXX) - 1;           //il punto di griglia lungo l'asse x a sinistra della particella
      else gridX = (int)gridXX;                      //

      gridYY = x2.getParticleY() / data.getDeltaY();
      if (gridYY < 0) gridY = ((int)gridYY) - 1;           //il punto di griglia lungo l'asse y a sinistra della particella
      else gridY = (int)gridYY;                      //

      gridZZ = x2.getParticleZ() / data.getDeltaZ();
      if (gridZZ < 0) gridZ = ((int)gridZZ) - 1;           //il punto di griglia lungo l'asse z a sinistra della particella
      else gridZ = (int)gridZZ;                      //

      coordinatePrimoPuntoGriglia.setParticleX(gridX*data.getDeltaX());
      coordinatePrimoPuntoGriglia.setParticleY(gridY*data.getDeltaY());
      coordinatePrimoPuntoGriglia.setParticleZ(gridZ*data.getDeltaZ());
      coordinateSecondoPuntoGriglia.setParticleX(gridX*data.getDeltaX() + data.getDeltaX());
      coordinateSecondoPuntoGriglia.setParticleY(gridY*data.getDeltaY());
      coordinateSecondoPuntoGriglia.setParticleZ(gridZ*data.getDeltaZ());
      coordinateTerzoPuntoGriglia.setParticleX(gridX*data.getDeltaX());
      coordinateTerzoPuntoGriglia.setParticleY(gridY*data.getDeltaY() + data.getDeltaY());
      coordinateTerzoPuntoGriglia.setParticleZ(gridZ*data.getDeltaZ());
      coordinateQuartoPuntoGriglia.setParticleX(gridX*data.getDeltaX() + data.getDeltaX());
      coordinateQuartoPuntoGriglia.setParticleY(gridY*data.getDeltaY() + data.getDeltaY());
      coordinateQuartoPuntoGriglia.setParticleZ(gridZ*data.getDeltaZ());
      coordinateQuintoPuntoGriglia.setParticleX(gridX*data.getDeltaX());
      coordinateQuintoPuntoGriglia.setParticleY(gridY*data.getDeltaY());
      coordinateQuintoPuntoGriglia.setParticleZ(gridZ*data.getDeltaZ() + data.getDeltaZ());
      coordinateSestoPuntoGriglia.setParticleX(gridX*data.getDeltaX() + data.getDeltaX());
      coordinateSestoPuntoGriglia.setParticleY(gridY*data.getDeltaY());
      coordinateSestoPuntoGriglia.setParticleZ(gridZ*data.getDeltaZ() + data.getDeltaZ());
      coordinateSettimoPuntoGriglia.setParticleX(gridX*data.getDeltaX());
      coordinateSettimoPuntoGriglia.setParticleY(gridY*data.getDeltaY() + data.getDeltaY());
      coordinateSettimoPuntoGriglia.setParticleZ(gridZ*data.getDeltaZ() + data.getDeltaZ());
      coordinateOttavoPuntoGriglia.setParticleX(gridX*data.getDeltaX() + data.getDeltaX());
      coordinateOttavoPuntoGriglia.setParticleY(gridY*data.getDeltaY() + data.getDeltaY());
      coordinateOttavoPuntoGriglia.setParticleZ(gridZ*data.getDeltaZ() + data.getDeltaZ());

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


      gamma = sqrt(1 + pow(x2.getParticlePX(), 2) + pow(x2.getParticlePY(), 2) + pow(x2.getParticlePZ(), 2));

      // AGGIORNA dx con i nuovi valori  DA SISTEMARE!!!!!
      dx.setParticleX(x2.getParticlePX() / (gamma * particles[i].getParticleM()));
      dx.setParticleY(x2.getParticlePY() / (gamma * particles[i].getParticleM()));
      dx.setParticleZ(x2.getParticlePZ() / (gamma * particles[i].getParticleM()));
      dx.setParticlePX(particles[i].getParticleQ() * fieldOnParticles[2].getEx() + particles[i].getParticleQ() * dx.getParticleY() * fieldOnParticles[2].getBz() - particles[i].getParticleQ() * dx.getParticleZ() * fieldOnParticles[2].getBy());  //CONTROLLARE I SEGNI!!!!
      dx.setParticlePY(particles[i].getParticleQ() * fieldOnParticles[2].getEy() + particles[i].getParticleQ() * dx.getParticleZ() * fieldOnParticles[2].getBx() - particles[i].getParticleQ() * dx.getParticleX() * fieldOnParticles[2].getBz());
      dx.setParticlePZ(particles[i].getParticleQ() * fieldOnParticles[2].getEz() + particles[i].getParticleQ() * dx.getParticleX() * fieldOnParticles[2].getBy() - particles[i].getParticleQ() * dx.getParticleY() * fieldOnParticles[2].getBx());

      // AGGIUNGE al dx2 i valori calcolati al terzo passo
      dx2.setParticleX(dx2.getParticleX() + 2 * dx.getParticleX());
      dx2.setParticleY(dx2.getParticleY() + 2 * dx.getParticleY());
      dx2.setParticleZ(dx2.getParticleZ() + 2 * dx.getParticleZ());
      dx2.setParticlePX(dx2.getParticlePX() + 2 * dx.getParticlePX());
      dx2.setParticlePY(dx2.getParticlePY() + 2 * dx.getParticlePY());
      dx2.setParticlePZ(dx2.getParticlePZ() + 2 * dx.getParticlePZ());

      // RIEMPIE x2 con le posizione avanzate di particles[i]      // DA SISTEMARE!!!
      x2.setParticleX(particles[i].getParticleX() + dx.getParticleX() * data.getDeltaT());
      x2.setParticleY(particles[i].getParticleY() + dx.getParticleY() * data.getDeltaT());
      x2.setParticleZ(particles[i].getParticleZ() + dx.getParticleZ() * data.getDeltaT());
      x2.setParticlePX(particles[i].getParticlePX() + dx.getParticlePX() *  data.getDeltaT());
      x2.setParticlePY(particles[i].getParticlePY() + dx.getParticlePY() *  data.getDeltaT());
      x2.setParticlePZ(particles[i].getParticlePZ() + dx.getParticlePZ() *  data.getDeltaT());



      // PASSO #4
      // AVANZA TEMPO DI deltaT/2 e calcola campi nella posizione x2
      gridXX = x2.getParticleX() / data.getDeltaX();
      if (gridXX < 0) gridX = ((int)gridXX) - 1;           //il punto di griglia lungo l'asse x a sinistra della particella
      else gridX = (int)gridXX;

      gridYY = x2.getParticleY() / data.getDeltaY();
      if (gridYY < 0) gridY = ((int)gridYY) - 1;           //il punto di griglia lungo l'asse y a sinistra della particella
      else gridY = (int)gridYY;

      gridZZ = x2.getParticleZ() / data.getDeltaZ();
      if (gridZZ < 0) gridZ = ((int)gridZZ) - 1;           //il punto di griglia lungo l'asse z a sinistra della particella
      else gridZ = (int)gridZZ;

      coordinatePrimoPuntoGriglia.setParticleX(gridX*data.getDeltaX());
      coordinatePrimoPuntoGriglia.setParticleY(gridY*data.getDeltaY());
      coordinatePrimoPuntoGriglia.setParticleZ(gridZ*data.getDeltaZ());
      coordinateSecondoPuntoGriglia.setParticleX(gridX*data.getDeltaX() + data.getDeltaX());
      coordinateSecondoPuntoGriglia.setParticleY(gridY*data.getDeltaY());
      coordinateSecondoPuntoGriglia.setParticleZ(gridZ*data.getDeltaZ());
      coordinateTerzoPuntoGriglia.setParticleX(gridX*data.getDeltaX());
      coordinateTerzoPuntoGriglia.setParticleY(gridY*data.getDeltaY() + data.getDeltaY());
      coordinateTerzoPuntoGriglia.setParticleZ(gridZ*data.getDeltaZ());
      coordinateQuartoPuntoGriglia.setParticleX(gridX*data.getDeltaX() + data.getDeltaX());
      coordinateQuartoPuntoGriglia.setParticleY(gridY*data.getDeltaY() + data.getDeltaY());
      coordinateQuartoPuntoGriglia.setParticleZ(gridZ*data.getDeltaZ());
      coordinateQuintoPuntoGriglia.setParticleX(gridX*data.getDeltaX());
      coordinateQuintoPuntoGriglia.setParticleY(gridY*data.getDeltaY());
      coordinateQuintoPuntoGriglia.setParticleZ(gridZ*data.getDeltaZ() + data.getDeltaZ());
      coordinateSestoPuntoGriglia.setParticleX(gridX*data.getDeltaX() + data.getDeltaX());
      coordinateSestoPuntoGriglia.setParticleY(gridY*data.getDeltaY());
      coordinateSestoPuntoGriglia.setParticleZ(gridZ*data.getDeltaZ() + data.getDeltaZ());
      coordinateSettimoPuntoGriglia.setParticleX(gridX*data.getDeltaX());
      coordinateSettimoPuntoGriglia.setParticleY(gridY*data.getDeltaY() + data.getDeltaY());
      coordinateSettimoPuntoGriglia.setParticleZ(gridZ*data.getDeltaZ() + data.getDeltaZ());
      coordinateOttavoPuntoGriglia.setParticleX(gridX*data.getDeltaX() + data.getDeltaX());
      coordinateOttavoPuntoGriglia.setParticleY(gridY*data.getDeltaY() + data.getDeltaY());
      coordinateOttavoPuntoGriglia.setParticleZ(gridZ*data.getDeltaZ() + data.getDeltaZ());

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

      data.aumentaT(data.getDeltaT() * 0.5);
      FOP = interpolation2D_linear(data, x2, fieldOnGrid);
      fieldOnParticles.push_back(FOP);

      fieldOnGrid.clear();

      gamma = sqrt(1 + pow(x2.getParticlePX(), 2) + pow(x2.getParticlePY(), 2) + pow(x2.getParticlePZ(), 2));

      //AGGIORNA dx con i nuovi valori  DA SISTEMARE!!!!!
      dx.setParticleX(x2.getParticlePX() / (gamma * particles[i].getParticleM()));
      dx.setParticleY(x2.getParticlePY() / (gamma * particles[i].getParticleM()));
      dx.setParticleZ(x2.getParticlePZ() / (gamma * particles[i].getParticleM()));
      dx.setParticlePX(particles[i].getParticleQ() * fieldOnParticles[3].getEx() + particles[i].getParticleQ() * dx.getParticleY() * fieldOnParticles[3].getBz() - particles[i].getParticleQ() * dx.getParticleZ() * fieldOnParticles[3].getBy());  //CONTROLLARE I SEGNI!!!!
      dx.setParticlePY(particles[i].getParticleQ() * fieldOnParticles[3].getEy() + particles[i].getParticleQ() * dx.getParticleZ() * fieldOnParticles[3].getBx() - particles[i].getParticleQ() * dx.getParticleX() * fieldOnParticles[3].getBz());
      dx.setParticlePZ(particles[i].getParticleQ() * fieldOnParticles[3].getEz() + particles[i].getParticleQ() * dx.getParticleX() * fieldOnParticles[3].getBy() - particles[i].getParticleQ() * dx.getParticleY() * fieldOnParticles[3].getBx());


      //AGGIUNGE al dx2 i valori calcolati al terzo passo
      dx2.setParticleX(dx2.getParticleX() + dx.getParticleX());
      dx2.setParticleY(dx2.getParticleY() + dx.getParticleY());
      dx2.setParticleZ(dx2.getParticleZ() + dx.getParticleZ());
      dx2.setParticlePX(dx2.getParticlePX() + dx.getParticlePX());
      dx2.setParticlePY(dx2.getParticlePY() + dx.getParticlePY());
      dx2.setParticlePZ(dx2.getParticlePZ() + dx.getParticlePZ());


      particles[i].setParticleX(particles[i].getParticleX() + (data.getDeltaT() / 6.) * dx2.getParticleX());
      particles[i].setParticleY(particles[i].getParticleY() + (data.getDeltaT() / 6.) * dx2.getParticleY());
      particles[i].setParticleZ(particles[i].getParticleZ() + (data.getDeltaT() / 6.) * dx2.getParticleZ());
      particles[i].setParticlePX(particles[i].getParticlePX() + (data.getDeltaT() / 6.) * dx2.getParticlePX());
      particles[i].setParticlePY(particles[i].getParticlePY() + (data.getDeltaT() / 6.) * dx2.getParticlePY());
      particles[i].setParticlePZ(particles[i].getParticlePZ() + (data.getDeltaT() / 6.) * dx2.getParticlePZ());

      particles[i].setParticleCell(data);

      data.aumentaT(-data.getDeltaT());       // faccio tornare il tempo indietro cosi` gli altri elettroni di questo step non hanno il tempo gia` avanzato
                                              // il tempo avanza definitivamente in ciascun ciclo for esterno (sugli nsteps)

      //INTEGRALI PRIMI
      potenziale.clear();
      zeta = particles[i].getParticleX() - C * data.getT();       // non usato nella mia struttura delle equazioni di campo
      gamma = sqrt(1 + pow(particles[i].getParticlePX(), 2) + pow(particles[i].getParticlePY(), 2) + pow(particles[i].getParticlePZ(), 2));

      gridXX = particles[i].getParticleX() / data.getDeltaX();
      if (gridXX < 0) gridX = ((int)gridXX) - 1;           //il punto di griglia lungo l'asse x a sinistra della particella
      else gridX = (int)gridXX;                      //

      gridYY = particles[i].getParticleY() / data.getDeltaY();
      if (gridYY < 0) gridY = ((int)gridYY) - 1;           //il punto di griglia lungo l'asse y a sinistra della particella
      else gridY = (int)gridYY;                      //

      gridZZ = particles[i].getParticleZ() / data.getDeltaZ();
      if (gridZZ < 0) gridZ = ((int)gridZZ) - 1;           //il punto di griglia lungo l'asse z a sinistra della particella
      else gridZ = (int)gridZZ;                      //

      coordinatePrimoPuntoGriglia.setParticleX(gridX*data.getDeltaX());
      coordinatePrimoPuntoGriglia.setParticleY(gridY*data.getDeltaY());
      coordinatePrimoPuntoGriglia.setParticleZ(gridZ*data.getDeltaZ());
      coordinateSecondoPuntoGriglia.setParticleX(gridX*data.getDeltaX() + data.getDeltaX());
      coordinateSecondoPuntoGriglia.setParticleY(gridY*data.getDeltaY());
      coordinateSecondoPuntoGriglia.setParticleZ(gridZ*data.getDeltaZ());
      coordinateTerzoPuntoGriglia.setParticleX(gridX*data.getDeltaX());
      coordinateTerzoPuntoGriglia.setParticleY(gridY*data.getDeltaY() + data.getDeltaY());
      coordinateTerzoPuntoGriglia.setParticleZ(gridZ*data.getDeltaZ());
      coordinateQuartoPuntoGriglia.setParticleX(gridX*data.getDeltaX() + data.getDeltaX());
      coordinateQuartoPuntoGriglia.setParticleY(gridY*data.getDeltaY() + data.getDeltaY());
      coordinateQuartoPuntoGriglia.setParticleZ(gridZ*data.getDeltaZ());
      coordinateQuintoPuntoGriglia.setParticleX(gridX*data.getDeltaX());
      coordinateQuintoPuntoGriglia.setParticleY(gridY*data.getDeltaY());
      coordinateQuintoPuntoGriglia.setParticleZ(gridZ*data.getDeltaZ() + data.getDeltaZ());
      coordinateSestoPuntoGriglia.setParticleX(gridX*data.getDeltaX() + data.getDeltaX());
      coordinateSestoPuntoGriglia.setParticleY(gridY*data.getDeltaY());
      coordinateSestoPuntoGriglia.setParticleZ(gridZ*data.getDeltaZ() + data.getDeltaZ());
      coordinateSettimoPuntoGriglia.setParticleX(gridX*data.getDeltaX());
      coordinateSettimoPuntoGriglia.setParticleY(gridY*data.getDeltaY() + data.getDeltaY());
      coordinateSettimoPuntoGriglia.setParticleZ(gridZ*data.getDeltaZ() + data.getDeltaZ());
      coordinateOttavoPuntoGriglia.setParticleX(gridX*data.getDeltaX() + data.getDeltaX());
      coordinateOttavoPuntoGriglia.setParticleY(gridY*data.getDeltaY() + data.getDeltaY());
      coordinateOttavoPuntoGriglia.setParticleZ(gridZ*data.getDeltaZ() + data.getDeltaZ());

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

      I1 = (particles[i].getParticlePY() / (particles[i].getParticleM()*C)) + potenziale[0].getEy();
      I2 = particles[i].getParticlePZ() / (particles[i].getParticleM()*C);
      I3 = gamma - particles[i].getParticlePX() / (particles[i].getParticleM()*C);

      potenziale.clear();

      outputstream
        << setw(12) << i + 1 << "\t"
        << setprecision(14) << setw(16) << particles[i].getParticleX() << "\t"
        << setw(16) << particles[i].getParticleY() << "\t"
        << setw(16) << particles[i].getParticleZ() << "\t"
        //        << setw(16) << particles[i].getParticleCell() << "\t"
        << setw(16) << particles[i].getParticlePX() << "\t"
        << setw(16) << particles[i].getParticlePY() << "\t"
        << setw(16) << particles[i].getParticlePZ() << "\t"
        << setw(16) << zeta << "\t"
        << setw(16) << I1 << "\t"
        << setw(16) << I2 << "\t"
        << setw(16) << I3 << endl;
    }
  }
}


void evolveRK4_withgrid(Data data, vector<Particle> &particles, vector<Field> &fieldOnParticles, vector<Field> &campoSuPuntiGriglia, ofstream &outputstream)
{

  Particle x2, dx, dx2;
  double gamma;
  double zeta;
  double I1, I2, I3;

  double gridXX, gridYY, gridZZ;
  int gridX, gridY, gridZ;
  Field FOP;
  vector<Field> fieldOnGrid;
  Field *FOG;                             //DA SISTEMARE, c'e` solo per avere un campo fittizio su una griglia on-the-fly per implementare l'algoritmo
  FOG = new Field[1];

  Field *campoSuPunto;
  campoSuPunto = new Field[1];
  vector<Field> potenziale;

  for (int j = 0; j < data.getNsteps(); j++)
  {
    data.aumentaT(data.getDeltaT());
    for (int i = 0; i < data.getNelectrons(); i++)
    {

      fieldOnParticles.clear();       //SVUOTA il vettore dei campi sulle particelle ad ogni ciclo

      gridXX = particles[i].getParticleX() / data.getDeltaX();
      if (gridXX < 0) gridX = ((int)gridXX) - 1;           //il punto di griglia lungo l'asse x a sinistra della particella
      else gridX = (int)gridXX;                      //

      gridYY = particles[i].getParticleY() / data.getDeltaY();
      if (gridYY < 0) gridY = ((int)gridYY) - 1;           //il punto di griglia lungo l'asse y a sinistra della particella
      else gridY = (int)gridYY;                      //

      gridZZ = particles[i].getParticleZ() / data.getDeltaZ();
      if (gridZZ < 0) gridZ = ((int)gridZZ) - 1;           //il punto di griglia lungo l'asse z a sinistra della particella
      else gridZ = (int)gridZZ;                      //

      if (gridX < 0 || gridY < 0 || gridZ < 0)
      {
        cout << "La particella e` uscita dai boundaries - simulazione interrotta" << endl;
        return;
      }

      campoSuPuntiGriglia.clear();
      riempiPuntiGrigliaConCampoAnalitico1D(data, campoSuPuntiGriglia);

      fieldOnGrid.push_back(campoSuPuntiGriglia.at(gridY * data.getNgridPointsX() + gridX));                // primo punto
      fieldOnGrid.push_back(campoSuPuntiGriglia.at(gridY * data.getNgridPointsX() + (gridX + 1)));          // secondo punto
      fieldOnGrid.push_back(campoSuPuntiGriglia.at((gridY + 1) * data.getNgridPointsX() + gridX));          // terzo punto
      fieldOnGrid.push_back(campoSuPuntiGriglia.at((gridY + 1) * data.getNgridPointsX() + (gridX + 1)));    // quarto punto

      FOP = interpolation2D_linear(data, particles[i], fieldOnGrid);
      fieldOnParticles.push_back(FOP);

      fieldOnGrid.clear();


      gamma = sqrt(1 + pow(particles[i].getParticlePX(), 2) + pow(particles[i].getParticlePY(), 2) + pow(particles[i].getParticlePZ(), 2));


      // PASSO #1
      //RIEMPIE dx = vx, vy, vz, Fx, Fy, Fz
      dx.setParticleX(particles[i].getParticlePX() / (gamma * particles[i].getParticleM()));
      dx.setParticleY(particles[i].getParticlePY() / (gamma * particles[i].getParticleM()));
      dx.setParticleZ(particles[i].getParticlePZ() / (gamma * particles[i].getParticleM()));
      dx.setParticlePX(particles[i].getParticleQ() * fieldOnParticles[0].getEx() + particles[i].getParticleQ() * dx.getParticleY() * fieldOnParticles[0].getBz() - particles[i].getParticleQ() * dx.getParticleZ() * fieldOnParticles[0].getBy());  //CONTROLLARE I SEGNI!!!!
      dx.setParticlePY(particles[i].getParticleQ() * fieldOnParticles[0].getEy() + particles[i].getParticleQ() * dx.getParticleZ() * fieldOnParticles[0].getBx() - particles[i].getParticleQ() * dx.getParticleX() * fieldOnParticles[0].getBz());
      dx.setParticlePZ(particles[i].getParticleQ() * fieldOnParticles[0].getEz() + particles[i].getParticleQ() * dx.getParticleX() * fieldOnParticles[0].getBy() - particles[i].getParticleQ() * dx.getParticleY() * fieldOnParticles[0].getBx());

      //RIEMPIE dx2 per il metodo RK4 con il dx del primo passo
      dx2 = dx;

      //RIEMPIE x2 con le posizioni avanzate di particles[i] dopo il primo passo
      x2.setParticleX(particles[i].getParticleX() + dx.getParticleX() * (data.getDeltaT()*0.5));      // E LE FORZE???
      x2.setParticleY(particles[i].getParticleY() + dx.getParticleY() * (data.getDeltaT()*0.5));
      x2.setParticleZ(particles[i].getParticleZ() + dx.getParticleZ() * (data.getDeltaT()*0.5));
      x2.setParticlePX(particles[i].getParticlePX() + dx.getParticlePX() *  (data.getDeltaT()*0.5));
      x2.setParticlePY(particles[i].getParticlePY() + dx.getParticlePY() *  (data.getDeltaT()*0.5));
      x2.setParticlePZ(particles[i].getParticlePZ() + dx.getParticlePZ() *  (data.getDeltaT()*0.5));


      // PASSO #2
      //AVANZA TEMPO DI deltaT/2 e calcola campi nella posizione x2
      gridXX = x2.getParticleX() / data.getDeltaX();
      if (gridXX < 0) gridX = ((int)gridXX) - 1;           //il punto di griglia lungo l'asse x a sinistra della particella
      else gridX = (int)gridXX;                      //

      gridYY = x2.getParticleY() / data.getDeltaY();
      if (gridYY < 0) gridY = ((int)gridYY) - 1;           //il punto di griglia lungo l'asse y a sinistra della particella
      else gridY = (int)gridYY;                      //

      gridZZ = x2.getParticleZ() / data.getDeltaZ();
      if (gridZZ < 0) gridZ = ((int)gridZZ) - 1;           //il punto di griglia lungo l'asse z a sinistra della particella
      else gridZ = (int)gridZZ;                      //

      if (gridX < 0 || gridY < 0 || gridZ < 0)
      {
        cout << "La particella e` uscita dai boundaries - simulazione interrotta" << endl;
        return;
      }


      data.aumentaT(data.getDeltaT() * 0.5);

      campoSuPuntiGriglia.clear();
      riempiPuntiGrigliaConCampoAnalitico1D(data, campoSuPuntiGriglia);

      fieldOnGrid.push_back(campoSuPuntiGriglia.at(gridY * data.getNgridPointsX() + gridX));        // primo punto
      fieldOnGrid.push_back(campoSuPuntiGriglia.at(gridY * data.getNgridPointsX() + (gridX + 1)));      // secondo punto
      fieldOnGrid.push_back(campoSuPuntiGriglia.at((gridY + 1) * data.getNgridPointsX() + gridX));      // terzo punto
      fieldOnGrid.push_back(campoSuPuntiGriglia.at((gridY + 1) * data.getNgridPointsX() + (gridX + 1)));    // quarto punto

      FOP = interpolation2D_linear(data, x2, fieldOnGrid);
      fieldOnParticles.push_back(FOP);

      fieldOnGrid.clear();

      gamma = sqrt(1 + pow(x2.getParticlePX(), 2) + pow(x2.getParticlePY(), 2) + pow(x2.getParticlePZ(), 2));

      //AGGIORNA dx con i nuovi valori
      dx.setParticleX(x2.getParticlePX() / (gamma * particles[i].getParticleM()));
      dx.setParticleY(x2.getParticlePY() / (gamma * particles[i].getParticleM()));
      dx.setParticleZ(x2.getParticlePZ() / (gamma * particles[i].getParticleM()));
      dx.setParticlePX(particles[i].getParticleQ() * fieldOnParticles[1].getEx() + particles[i].getParticleQ() * dx.getParticleY() * fieldOnParticles[1].getBz() - particles[i].getParticleQ() * dx.getParticleZ() * fieldOnParticles[1].getBy());  //CONTROLLARE I SEGNI!!!!
      dx.setParticlePY(particles[i].getParticleQ() * fieldOnParticles[1].getEy() + particles[i].getParticleQ() * dx.getParticleZ() * fieldOnParticles[1].getBx() - particles[i].getParticleQ() * dx.getParticleX() * fieldOnParticles[1].getBz());
      dx.setParticlePZ(particles[i].getParticleQ() * fieldOnParticles[1].getEz() + particles[i].getParticleQ() * dx.getParticleX() * fieldOnParticles[1].getBy() - particles[i].getParticleQ() * dx.getParticleY() * fieldOnParticles[1].getBx());


      //AGGIUNGE al dx2 i valori calcolati al secondo passo
      dx2.setParticleX(dx2.getParticleX() + 2 * dx.getParticleX());
      dx2.setParticleY(dx2.getParticleY() + 2 * dx.getParticleY());
      dx2.setParticleZ(dx2.getParticleZ() + 2 * dx.getParticleZ());
      dx2.setParticlePX(dx2.getParticlePX() + 2 * dx.getParticlePX());
      dx2.setParticlePY(dx2.getParticlePY() + 2 * dx.getParticlePY());
      dx2.setParticlePZ(dx2.getParticlePZ() + 2 * dx.getParticlePZ());

      //RIEMPIE x2 con le posizione avanzate di particles[i] dopo il secondo passo  // DA SISTEMARE!!!
      x2.setParticleX(particles[i].getParticleX() + dx.getParticleX() * (data.getDeltaT()*0.5));
      x2.setParticleY(particles[i].getParticleY() + dx.getParticleY() * (data.getDeltaT()*0.5));
      x2.setParticleZ(particles[i].getParticleZ() + dx.getParticleZ() * (data.getDeltaT()*0.5));
      x2.setParticlePX(particles[i].getParticlePX() + dx.getParticlePX() *  (data.getDeltaT()*0.5));
      x2.setParticlePY(particles[i].getParticlePY() + dx.getParticlePY() *  (data.getDeltaT()*0.5));
      x2.setParticlePZ(particles[i].getParticlePZ() + dx.getParticlePZ() *  (data.getDeltaT()*0.5));



      // PASSO #3
      //CALCOLA campi nella posizione x2 appena aggiornata
      gridXX = x2.getParticleX() / data.getDeltaX();
      if (gridXX < 0) gridX = ((int)gridXX) - 1;           //il punto di griglia lungo l'asse x a sinistra della particella
      else gridX = (int)gridXX;                      //

      gridYY = x2.getParticleY() / data.getDeltaY();
      if (gridYY < 0) gridY = ((int)gridYY) - 1;           //il punto di griglia lungo l'asse y a sinistra della particella
      else gridY = (int)gridYY;                      //

      gridZZ = x2.getParticleZ() / data.getDeltaZ();
      if (gridZZ < 0) gridZ = ((int)gridZZ) - 1;           //il punto di griglia lungo l'asse z a sinistra della particella
      else gridZ = (int)gridZZ;                      //

      if (gridX < 0 || gridY < 0 || gridZ < 0)
      {
        cout << "La particella e` uscita dai boundaries - simulazione interrotta" << endl;
        return;
      }


      //      riempiPuntiGrigliaConCampoAnalitico1D(data, campoSuPuntiGriglia);   // non e` necessario ricalcolarlo perche' il tempo non e` variato

      fieldOnGrid.push_back(campoSuPuntiGriglia.at(gridY * data.getNgridPointsX() + gridX));        // primo punto
      fieldOnGrid.push_back(campoSuPuntiGriglia.at(gridY * data.getNgridPointsX() + (gridX + 1)));      // secondo punto
      fieldOnGrid.push_back(campoSuPuntiGriglia.at((gridY + 1) * data.getNgridPointsX() + gridX));      // terzo punto
      fieldOnGrid.push_back(campoSuPuntiGriglia.at((gridY + 1) * data.getNgridPointsX() + (gridX + 1)));    // quarto punto

      FOP = interpolation2D_linear(data, x2, fieldOnGrid);
      fieldOnParticles.push_back(FOP);

      fieldOnGrid.clear();


      gamma = sqrt(1 + pow(x2.getParticlePX(), 2) + pow(x2.getParticlePY(), 2) + pow(x2.getParticlePZ(), 2));

      //AGGIORNA dx con i nuovi valori  DA SISTEMARE!!!!!
      dx.setParticleX(x2.getParticlePX() / (gamma * particles[i].getParticleM()));
      dx.setParticleY(x2.getParticlePY() / (gamma * particles[i].getParticleM()));
      dx.setParticleZ(x2.getParticlePZ() / (gamma * particles[i].getParticleM()));
      dx.setParticlePX(particles[i].getParticleQ() * fieldOnParticles[2].getEx() + particles[i].getParticleQ() * dx.getParticleY() * fieldOnParticles[2].getBz() - particles[i].getParticleQ() * dx.getParticleZ() * fieldOnParticles[2].getBy());  //CONTROLLARE I SEGNI!!!!
      dx.setParticlePY(particles[i].getParticleQ() * fieldOnParticles[2].getEy() + particles[i].getParticleQ() * dx.getParticleZ() * fieldOnParticles[2].getBx() - particles[i].getParticleQ() * dx.getParticleX() * fieldOnParticles[2].getBz());
      dx.setParticlePZ(particles[i].getParticleQ() * fieldOnParticles[2].getEz() + particles[i].getParticleQ() * dx.getParticleX() * fieldOnParticles[2].getBy() - particles[i].getParticleQ() * dx.getParticleY() * fieldOnParticles[2].getBx());

      //AGGIUNGE al dx2 i valori calcolati al terzo passo
      dx2.setParticleX(dx2.getParticleX() + 2 * dx.getParticleX());
      dx2.setParticleY(dx2.getParticleY() + 2 * dx.getParticleY());
      dx2.setParticleZ(dx2.getParticleZ() + 2 * dx.getParticleZ());
      dx2.setParticlePX(dx2.getParticlePX() + 2 * dx.getParticlePX());
      dx2.setParticlePY(dx2.getParticlePY() + 2 * dx.getParticlePY());
      dx2.setParticlePZ(dx2.getParticlePZ() + 2 * dx.getParticlePZ());

      //RIEMPIE x2 con le posizione avanzate di particles[i]      // DA SISTEMARE!!!
      x2.setParticleX(particles[i].getParticleX() + dx.getParticleX() * data.getDeltaT());
      x2.setParticleY(particles[i].getParticleY() + dx.getParticleY() * data.getDeltaT());
      x2.setParticleZ(particles[i].getParticleZ() + dx.getParticleZ() * data.getDeltaT());
      x2.setParticlePX(particles[i].getParticlePX() + dx.getParticlePX() *  data.getDeltaT());
      x2.setParticlePY(particles[i].getParticlePY() + dx.getParticlePY() *  data.getDeltaT());
      x2.setParticlePZ(particles[i].getParticlePZ() + dx.getParticlePZ() *  data.getDeltaT());



      // PASSO #4
      //AVANZA TEMPO DI deltaT/2 e calcola campi nella posizione x2
      gridXX = x2.getParticleX() / data.getDeltaX();
      if (gridXX < 0) gridX = ((int)gridXX) - 1;           //il punto di griglia lungo l'asse x a sinistra della particella
      else gridX = (int)gridXX;

      gridYY = x2.getParticleY() / data.getDeltaY();
      if (gridYY < 0) gridY = ((int)gridYY) - 1;           //il punto di griglia lungo l'asse y a sinistra della particella
      else gridY = (int)gridYY;

      gridZZ = x2.getParticleZ() / data.getDeltaZ();
      if (gridZZ < 0) gridZ = ((int)gridZZ) - 1;           //il punto di griglia lungo l'asse z a sinistra della particella
      else gridZ = (int)gridZZ;

      if (gridX < 0 || gridY < 0 || gridZ < 0)
      {
        cout << "La particella e` uscita dai boundaries - simulazione interrotta" << endl;
        return;
      }

      data.aumentaT(data.getDeltaT() * 0.5);

      campoSuPuntiGriglia.clear();
      riempiPuntiGrigliaConCampoAnalitico1D(data, campoSuPuntiGriglia);

      fieldOnGrid.push_back(campoSuPuntiGriglia.at(gridY * data.getNgridPointsX() + gridX));        // primo punto
      fieldOnGrid.push_back(campoSuPuntiGriglia.at(gridY * data.getNgridPointsX() + (gridX + 1)));      // secondo punto
      fieldOnGrid.push_back(campoSuPuntiGriglia.at((gridY + 1) * data.getNgridPointsX() + gridX));      // terzo punto
      fieldOnGrid.push_back(campoSuPuntiGriglia.at((gridY + 1) * data.getNgridPointsX() + (gridX + 1)));    // quarto punto

      FOP = interpolation2D_linear(data, x2, fieldOnGrid);
      fieldOnParticles.push_back(FOP);

      fieldOnGrid.clear();

      gamma = sqrt(1 + pow(x2.getParticlePX(), 2) + pow(x2.getParticlePY(), 2) + pow(x2.getParticlePZ(), 2));

      //AGGIORNA dx con i nuovi valori  DA SISTEMARE!!!!!
      dx.setParticleX(x2.getParticlePX() / (gamma * particles[i].getParticleM()));
      dx.setParticleY(x2.getParticlePY() / (gamma * particles[i].getParticleM()));
      dx.setParticleZ(x2.getParticlePZ() / (gamma * particles[i].getParticleM()));
      dx.setParticlePX(particles[i].getParticleQ() * fieldOnParticles[3].getEx() + particles[i].getParticleQ() * dx.getParticleY() * fieldOnParticles[3].getBz() - particles[i].getParticleQ() * dx.getParticleZ() * fieldOnParticles[3].getBy());  //CONTROLLARE I SEGNI!!!!
      dx.setParticlePY(particles[i].getParticleQ() * fieldOnParticles[3].getEy() + particles[i].getParticleQ() * dx.getParticleZ() * fieldOnParticles[3].getBx() - particles[i].getParticleQ() * dx.getParticleX() * fieldOnParticles[3].getBz());
      dx.setParticlePZ(particles[i].getParticleQ() * fieldOnParticles[3].getEz() + particles[i].getParticleQ() * dx.getParticleX() * fieldOnParticles[3].getBy() - particles[i].getParticleQ() * dx.getParticleY() * fieldOnParticles[3].getBx());


      //AGGIUNGE al dx2 i valori calcolati al terzo passo
      dx2.setParticleX(dx2.getParticleX() + dx.getParticleX());
      dx2.setParticleY(dx2.getParticleY() + dx.getParticleY());
      dx2.setParticleZ(dx2.getParticleZ() + dx.getParticleZ());
      dx2.setParticlePX(dx2.getParticlePX() + dx.getParticlePX());
      dx2.setParticlePY(dx2.getParticlePY() + dx.getParticlePY());
      dx2.setParticlePZ(dx2.getParticlePZ() + dx.getParticlePZ());


      particles[i].setParticleX(particles[i].getParticleX() + (data.getDeltaT() / 6.) * dx2.getParticleX());
      particles[i].setParticleY(particles[i].getParticleY() + (data.getDeltaT() / 6.) * dx2.getParticleY());
      particles[i].setParticleZ(particles[i].getParticleZ() + (data.getDeltaT() / 6.) * dx2.getParticleZ());
      particles[i].setParticlePX(particles[i].getParticlePX() + (data.getDeltaT() / 6.) * dx2.getParticlePX());
      particles[i].setParticlePY(particles[i].getParticlePY() + (data.getDeltaT() / 6.) * dx2.getParticlePY());
      particles[i].setParticlePZ(particles[i].getParticlePZ() + (data.getDeltaT() / 6.) * dx2.getParticlePZ());

      particles[i].setParticleCell(data);

      data.aumentaT(-data.getDeltaT());       //faccio tornare il tempo indietro cosi` gli altri elettroni di questo step non hanno il tempo gia` avanzato
                              //il tempo avanza definitivamente in ciascun ciclo for esterno (sugli nsteps)

      //INTEGRALI PRIMI da sistemare perche' calcolano il potenziale sulla particella e non sulla griglia e poi interpolando
      potenziale.clear();
      zeta = particles[i].getParticleX() - C * data.getT();             // non usato nella mia struttura delle equazioni di campo
      gamma = sqrt(1 + pow(particles[i].getParticlePX(), 2) + pow(particles[i].getParticlePY(), 2) + pow(particles[i].getParticlePZ(), 2));

      gridXX = particles[i].getParticleX() / data.getDeltaX();
      if (gridXX < 0) gridX = ((int)gridXX) - 1;           //il punto di griglia lungo l'asse x a sinistra della particella
      else gridX = (int)gridXX;                      //

      gridYY = particles[i].getParticleY() / data.getDeltaY();
      if (gridYY < 0) gridY = ((int)gridYY) - 1;           //il punto di griglia lungo l'asse y a sinistra della particella
      else gridY = (int)gridYY;                      //

      gridZZ = particles[i].getParticleZ() / data.getDeltaZ();
      if (gridZZ < 0) gridZ = ((int)gridZZ) - 1;           //il punto di griglia lungo l'asse z a sinistra della particella
      else gridZ = (int)gridZZ;                      //

      if (gridX < 0 || gridY < 0 || gridZ < 0)
      {
        cout << "La particella e` uscita dai boundaries - simulazione interrotta" << endl;
        return;
      }

      campoSuPuntiGriglia.clear();
      riempiPuntiGrigliaConPotenzialeAnalitico1D(data, campoSuPuntiGriglia);

      fieldOnGrid.push_back(campoSuPuntiGriglia.at(gridY * data.getNgridPointsX() + gridX));                // primo punto
      fieldOnGrid.push_back(campoSuPuntiGriglia.at(gridY * data.getNgridPointsX() + (gridX + 1)));          // secondo punto
      fieldOnGrid.push_back(campoSuPuntiGriglia.at((gridY + 1) * data.getNgridPointsX() + gridX));          // terzo punto
      fieldOnGrid.push_back(campoSuPuntiGriglia.at((gridY + 1) * data.getNgridPointsX() + (gridX + 1)));    // quarto punto

      FOP = interpolation2D_linear(data, particles[i], fieldOnGrid);
      potenziale.push_back(FOP);

      fieldOnGrid.clear();

      I1 = (particles[i].getParticlePY() / (particles[i].getParticleM()*C)) + potenziale[0].getEy();
      I2 = particles[i].getParticlePZ() / (particles[i].getParticleM()*C);
      I3 = gamma - particles[i].getParticlePX() / (particles[i].getParticleM()*C);


      outputstream
        << setw(12) << i + 1 << "\t"
        << setprecision(14) << setw(16) << particles[i].getParticleX() << "\t"
        << setw(16) << particles[i].getParticleY() << "\t"
        << setw(16) << particles[i].getParticleZ() << "\t"
        //<< setw(16) << particles[i].getParticleCell() << "\t"
        << setw(16) << particles[i].getParticlePX() << "\t"
        << setw(16) << particles[i].getParticlePY() << "\t"
        << setw(16) << particles[i].getParticlePZ() << "\t"
        << setw(16) << zeta << "\t"
        << setw(16) << potenziale[0].getEy() << "\t"
        << setw(16) << I1 << "\t"
        << setw(16) << I2 << "\t"
        << setw(16) << I3 << endl;
    }
  }
}



Field obtainFOP(Data data, Particle particle, vector<Field> &FOG)     // Difettosa, era solo un primissimo tentativo di interpolazione
{
  int numberOfVertices = FOG.size();

  /*
  Se il numero di vertici (punti di griglia) dai quali bisogna interpolare il campo sulla particella sono:
  2 -> il primo e` alla sua sinistra, il secondo alla sua destra
  4 -> il primo e` il piu` vicino all'origine degli assi, il secondo si sposta lungo x, il terzo e` il primo spostato lungo y e l'ultimo e` il piu` lontano dagli assi
  8 -> con la stessa filosofia di prima, i primi 4 sono nello stesso ordine e sono quelli appartenenti al piano piu` vicino ad xy, gli altri 4 sono piu` lontani di un dz
  */

  Particle relativeCoord;
  double gridd;
  int grid;

  gridd = particle.getParticleX() / data.getDeltaX();
  if (gridd < 0) grid = ((int)gridd) - 1;
  else grid = (int)gridd;
  relativeCoord.setParticleX(particle.getParticleX() - data.getDeltaX()*grid);

  gridd = particle.getParticleY() / data.getDeltaY();
  if (gridd < 0) grid = ((int)gridd) - 1;
  else grid = (int)gridd;
  relativeCoord.setParticleY(particle.getParticleY() - data.getDeltaY()*grid);

  gridd = particle.getParticleZ() / data.getDeltaZ();
  if (gridd < 0) grid = ((int)gridd) - 1;
  else grid = (int)gridd;
  relativeCoord.setParticleZ(particle.getParticleZ() - data.getDeltaZ()*grid);

  Field tempField;

  double x2_vsFirstGridPoint = (pow(relativeCoord.getParticleX(), 2));
  double y2_vsFirstGridPoint = (pow(relativeCoord.getParticleY(), 2));
  double z2_vsFirstGridPoint = (pow(relativeCoord.getParticleZ(), 2));
  double pos2_vsFirstGridPoint = x2_vsFirstGridPoint + y2_vsFirstGridPoint + z2_vsFirstGridPoint;

  double x2_vsSecondGridPoint = (pow(data.getDeltaX() - relativeCoord.getParticleX(), 2));
  double y2_vsSecondGridPoint = (pow(relativeCoord.getParticleY(), 2));
  double z2_vsSecondGridPoint = (pow(relativeCoord.getParticleZ(), 2));
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

  //  if (numberOfVertices == 4 || numberOfVertices == 8)
  //  {
  x2_vsThirdGridPoint = (pow(relativeCoord.getParticleX(), 2));
  y2_vsThirdGridPoint = (pow(data.getDeltaY() - relativeCoord.getParticleZ(), 2));
  z2_vsThirdGridPoint = (pow(relativeCoord.getParticleZ(), 2));
  pos2_vsThirdGridPoint = x2_vsThirdGridPoint + y2_vsThirdGridPoint + z2_vsThirdGridPoint;

  x2_vsFourthGridPoint = (pow(data.getDeltaX() - relativeCoord.getParticleX(), 2));
  y2_vsFourthGridPoint = (pow(data.getDeltaY() - relativeCoord.getParticleY(), 2));
  z2_vsFourthGridPoint = (pow(relativeCoord.getParticleZ(), 2));
  pos2_vsFourthGridPoint = x2_vsFourthGridPoint + y2_vsFourthGridPoint + z2_vsFourthGridPoint;

  normal += pos2_vsThirdGridPoint + pos2_vsFourthGridPoint;
  //  }

  //  if (numberOfVertices == 8)
  //  {
  x2_vsFifthGridPoint = (pow(relativeCoord.getParticleX(), 2));
  y2_vsFifthGridPoint = (pow(relativeCoord.getParticleY(), 2));
  z2_vsFifthGridPoint = (pow(data.getDeltaZ() - relativeCoord.getParticleZ(), 2));
  pos2_vsFifthGridPoint = x2_vsFifthGridPoint + y2_vsFifthGridPoint + z2_vsFifthGridPoint;

  x2_vsSixthGridPoint = (pow(data.getDeltaX() - relativeCoord.getParticleX(), 2));
  y2_vsSixthGridPoint = (pow(relativeCoord.getParticleY(), 2));
  z2_vsSixthGridPoint = (pow(data.getDeltaZ() - relativeCoord.getParticleZ(), 2));
  pos2_vsSixthGridPoint = x2_vsSixthGridPoint + y2_vsSixthGridPoint + z2_vsSixthGridPoint;

  x2_vsSeventhGridPoint = (pow(relativeCoord.getParticleX(), 2));
  y2_vsSeventhGridPoint = (pow(data.getDeltaY() - relativeCoord.getParticleZ(), 2));
  z2_vsSeventhGridPoint = (pow(data.getDeltaZ() - relativeCoord.getParticleZ(), 2));
  pos2_vsSeventhGridPoint = x2_vsSeventhGridPoint + y2_vsSeventhGridPoint + z2_vsSeventhGridPoint;

  x2_vsEighthGridPoint = (pow(data.getDeltaX() - relativeCoord.getParticleX(), 2));
  y2_vsEighthGridPoint = (pow(data.getDeltaY() - relativeCoord.getParticleY(), 2));
  z2_vsEighthGridPoint = (pow(data.getDeltaZ() - relativeCoord.getParticleZ(), 2));
  pos2_vsEighthGridPoint = x2_vsEighthGridPoint + y2_vsEighthGridPoint + z2_vsEighthGridPoint;

  normal += pos2_vsFifthGridPoint + pos2_vsSixthGridPoint + pos2_vsSeventhGridPoint + pos2_vsEighthGridPoint;
  //  }


  double normal_inv;
  normal_inv = 1.0 / normal;

  pos2_vsFirstGridPoint *= normal_inv;
  pos2_vsSecondGridPoint *= normal_inv;

  double pos_vsFirstGridPoint = sqrt(pos2_vsFirstGridPoint);
  double pos_vsSecondGridPoint = sqrt(pos2_vsSecondGridPoint);

  double fieldX_0 = FOG[0].getEx() * pos_vsFirstGridPoint;
  double fieldX_1 = FOG[1].getEx() * pos_vsSecondGridPoint;
  double fieldY_0 = FOG[0].getEy() * pos_vsFirstGridPoint;
  double fieldY_1 = FOG[1].getEy() * pos_vsSecondGridPoint;
  double fieldZ_0 = FOG[0].getEz() * pos_vsFirstGridPoint;
  double fieldZ_1 = FOG[1].getEz() * pos_vsSecondGridPoint;

  double fieldBX_0 = FOG[0].getBx() * pos_vsFirstGridPoint;
  double fieldBX_1 = FOG[1].getBx() * pos_vsSecondGridPoint;
  double fieldBY_0 = FOG[0].getBy() * pos_vsFirstGridPoint;
  double fieldBY_1 = FOG[1].getBy() * pos_vsSecondGridPoint;
  double fieldBZ_0 = FOG[0].getBz() * pos_vsFirstGridPoint;
  double fieldBZ_1 = FOG[1].getBz() * pos_vsSecondGridPoint;

  //  if (numberOfVertices == 2)
  //  {
  tempField.setEx(fieldX_0 + fieldX_1);
  tempField.setEy(fieldY_0 + fieldY_1);
  tempField.setEz(fieldZ_0 + fieldZ_1);
  tempField.setBx(fieldBX_0 + fieldBX_1);
  tempField.setBy(fieldBY_0 + fieldBY_1);
  tempField.setBz(fieldBZ_0 + fieldBZ_1);
  //  }
  /*
  //  if (numberOfVertices == 4)
  //  {
      pos2_vsThirdGridPoint *= normal_inv;
      pos2_vsFourthGridPoint *= normal_inv;

      double pos_vsThirdGridPoint =   sqrt(pos2_vsThirdGridPoint);
      double pos_vsFourthGridPoint =   sqrt(pos2_vsFourthGridPoint);

      double fieldX_2 = FOG[2].getEx() * pos_vsThirdGridPoint;
      double fieldX_3 = FOG[3].getEx() * pos_vsFourthGridPoint;
      tempField.setEx(fieldX_0 + fieldX_1 + fieldX_2 + fieldX_3);
      double fieldBX_2 = FOG[2].getBx() * pos_vsThirdGridPoint;
      double fieldBX_3 = FOG[3].getBx() * pos_vsFourthGridPoint;
      tempField.setBx(fieldBX_0 + fieldBX_1 + fieldBX_2 + fieldBX_3);

      double fieldY_2 = FOG[2].getEy() * pos_vsThirdGridPoint;
      double fieldY_3 = FOG[3].getEy() * pos_vsFourthGridPoint;
      tempField.setEy(fieldY_0 + fieldY_1 + fieldY_2 + fieldY_3);
      double fieldBY_2 = FOG[2].getBy() * pos_vsThirdGridPoint;
      double fieldBY_3 = FOG[3].getBy() * pos_vsFourthGridPoint;
      tempField.setBy(fieldBY_0 + fieldBY_1 + fieldBY_2 + fieldBY_3);

      double fieldZ_2 = FOG[2].getEz() * pos_vsThirdGridPoint;
      double fieldZ_3 = FOG[3].getEz() * pos_vsFourthGridPoint;
      tempField.setEz(fieldZ_0 + fieldZ_1 + fieldZ_2 + fieldZ_3);
      double fieldBZ_2 = FOG[2].getBz() * pos_vsThirdGridPoint;
      double fieldBZ_3 = FOG[3].getBz() * pos_vsFourthGridPoint;
      tempField.setBz(fieldBZ_0 + fieldBZ_1 + fieldBZ_2 + fieldBZ_3);
    }
  */
  //  if (numberOfVertices == 8)
  //  {
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

  double fieldX_2 = FOG[2].getEx() * pos_vsThirdGridPoint;
  double fieldX_3 = FOG[3].getEx() * pos_vsFourthGridPoint;
  double fieldX_4 = FOG[4].getEx() * pos_vsFifthGridPoint;
  double fieldX_5 = FOG[5].getEx() * pos_vsSixthGridPoint;
  double fieldX_6 = FOG[6].getEx() * pos_vsSeventhGridPoint;
  double fieldX_7 = FOG[7].getEx() * pos_vsEighthGridPoint;
  tempField.setEx(fieldX_0 + fieldX_1 + fieldX_2 + fieldX_3 + fieldX_4 + fieldX_5 + fieldX_6 + fieldX_7);
  double fieldBX_2 = FOG[2].getBx() * pos_vsThirdGridPoint;
  double fieldBX_3 = FOG[3].getBx() * pos_vsFourthGridPoint;
  double fieldBX_4 = FOG[4].getBx() * pos_vsFifthGridPoint;
  double fieldBX_5 = FOG[5].getBx() * pos_vsSixthGridPoint;
  double fieldBX_6 = FOG[6].getBx() * pos_vsSeventhGridPoint;
  double fieldBX_7 = FOG[7].getBx() * pos_vsEighthGridPoint;
  tempField.setBx(fieldBX_0 + fieldBX_1 + fieldBX_2 + fieldBX_3 + fieldBX_4 + fieldBX_5 + fieldBX_6 + fieldBX_7);

  double fieldY_2 = FOG[2].getEy() * pos_vsThirdGridPoint;
  double fieldY_3 = FOG[3].getEy() * pos_vsFourthGridPoint;
  double fieldY_4 = FOG[4].getEy() * pos_vsFifthGridPoint;
  double fieldY_5 = FOG[5].getEy() * pos_vsSixthGridPoint;
  double fieldY_6 = FOG[6].getEy() * pos_vsSeventhGridPoint;
  double fieldY_7 = FOG[7].getEy() * pos_vsEighthGridPoint;
  tempField.setEy(fieldY_0 + fieldY_1 + fieldY_2 + fieldY_3 + fieldY_4 + fieldY_5 + fieldY_6 + fieldY_7);
  double fieldBY_2 = FOG[2].getBy() * pos_vsThirdGridPoint;
  double fieldBY_3 = FOG[3].getBy() * pos_vsFourthGridPoint;
  double fieldBY_4 = FOG[4].getBy() * pos_vsFifthGridPoint;
  double fieldBY_5 = FOG[5].getBy() * pos_vsSixthGridPoint;
  double fieldBY_6 = FOG[6].getBy() * pos_vsSeventhGridPoint;
  double fieldBY_7 = FOG[7].getBy() * pos_vsEighthGridPoint;
  tempField.setBy(fieldBY_0 + fieldBY_1 + fieldBY_2 + fieldBY_3 + fieldBY_4 + fieldBY_5 + fieldBY_6 + fieldBY_7);

  double fieldZ_2 = FOG[2].getEz() * pos_vsThirdGridPoint;
  double fieldZ_3 = FOG[3].getEz() * pos_vsFourthGridPoint;
  double fieldZ_4 = FOG[4].getEz() * pos_vsFifthGridPoint;
  double fieldZ_5 = FOG[5].getEz() * pos_vsSixthGridPoint;
  double fieldZ_6 = FOG[6].getEz() * pos_vsSeventhGridPoint;
  double fieldZ_7 = FOG[7].getEz() * pos_vsEighthGridPoint;
  tempField.setEz(fieldZ_0 + fieldZ_1 + fieldZ_2 + fieldZ_3 + fieldZ_4 + fieldZ_5 + fieldZ_6 + fieldZ_7);
  double fieldBZ_2 = FOG[2].getBz() * pos_vsThirdGridPoint;
  double fieldBZ_3 = FOG[3].getBz() * pos_vsFourthGridPoint;
  double fieldBZ_4 = FOG[4].getBz() * pos_vsFifthGridPoint;
  double fieldBZ_5 = FOG[5].getBz() * pos_vsSixthGridPoint;
  double fieldBZ_6 = FOG[6].getBz() * pos_vsSeventhGridPoint;
  double fieldBZ_7 = FOG[7].getBz() * pos_vsEighthGridPoint;
  tempField.setBz(fieldBZ_0 + fieldBZ_1 + fieldBZ_2 + fieldBZ_3 + fieldBZ_4 + fieldBZ_5 + fieldBZ_6 + fieldBZ_7);
  //  }

  /*
  PROTOTIPO
    - in input deve avere i dati riguardo la griglia e riguardo la posizione della particella, non un vettore ma solo 1 particella e i campi in quell'istante
    - da essi si calcola la distanza relativa tra i punti di griglia e la posizione della particella
    - in funzione di cio`, deve calcolare l'interpolazione dei campi sulla particella
  */
  return tempField;
}


Field interpolation2D_linear(Data data, Particle position, vector<Field> &valuesOnGrid)
{
  /*
  Se il numero di vertici (punti di griglia) dai quali bisogna interpolare il campo sulla particella sono diversi da 4, la funzione esce restituendo -1
  Se i vertici sono 4, vanno identificati in questo modo: il primo e` il piu` vicino all'origine degli assi, il secondo si sposta lungo x,
  il terzo e` il primo spostato lungo y e l'ultimo e` il piu` lontano dagli assi
  */

  Particle relativeCoord;

  double gridd;
  int grid = 0;
  double dx_inv = 1. / data.getDeltaX();
  double dy_inv = 1. / data.getDeltaY();

  gridd = position.getParticleX() / data.getDeltaX();
  if (gridd < 0) grid = ((int)grid) - 1;
  else grid = (int)gridd;
  relativeCoord.setParticleX((position.getParticleX() - data.getDeltaX()*grid) * dx_inv);

  gridd = position.getParticleY() / data.getDeltaY();
  if (gridd < 0) grid = ((int)grid) - 1;
  else grid = (int)gridd;
  relativeCoord.setParticleY((position.getParticleY() - data.getDeltaY()*grid) * dy_inv);

  relativeCoord.setParticleZ(0.);

  Field tempField;

  double weightOnFourthPoint = relativeCoord.getParticleX() * relativeCoord.getParticleY();
  double weightOnFirstPoint = (1 - relativeCoord.getParticleX()) * (1 - relativeCoord.getParticleY());
  double weightOnThirdPoint = relativeCoord.getParticleX() * (1 - relativeCoord.getParticleY());
  double weightOnSecondPoint = (1 - relativeCoord.getParticleX()) * relativeCoord.getParticleY();

  double fieldX_0 = valuesOnGrid[0].getEx() * weightOnFirstPoint;
  double fieldX_1 = valuesOnGrid[1].getEx() * weightOnSecondPoint;
  double fieldX_2 = valuesOnGrid[2].getEx() * weightOnThirdPoint;
  double fieldX_3 = valuesOnGrid[3].getEx() * weightOnFourthPoint;

  double fieldY_0 = valuesOnGrid[0].getEy() * weightOnFirstPoint;
  double fieldY_1 = valuesOnGrid[1].getEy() * weightOnSecondPoint;
  double fieldY_2 = valuesOnGrid[2].getEy() * weightOnThirdPoint;
  double fieldY_3 = valuesOnGrid[3].getEy() * weightOnFourthPoint;

  double fieldZ_0 = valuesOnGrid[0].getEz() * weightOnFirstPoint;
  double fieldZ_1 = valuesOnGrid[1].getEz() * weightOnSecondPoint;
  double fieldZ_2 = valuesOnGrid[2].getEz() * weightOnThirdPoint;
  double fieldZ_3 = valuesOnGrid[3].getEz() * weightOnFourthPoint;

  double fieldBX_0 = valuesOnGrid[0].getBx() * weightOnFirstPoint;
  double fieldBX_1 = valuesOnGrid[1].getBx() * weightOnSecondPoint;
  double fieldBX_2 = valuesOnGrid[2].getBx() * weightOnThirdPoint;
  double fieldBX_3 = valuesOnGrid[3].getBx() * weightOnFourthPoint;

  double fieldBY_0 = valuesOnGrid[0].getBy() * weightOnFirstPoint;
  double fieldBY_1 = valuesOnGrid[1].getBy() * weightOnSecondPoint;
  double fieldBY_2 = valuesOnGrid[2].getBy() * weightOnThirdPoint;
  double fieldBY_3 = valuesOnGrid[3].getBy() * weightOnFourthPoint;

  double fieldBZ_0 = valuesOnGrid[0].getBz() * weightOnFirstPoint;
  double fieldBZ_1 = valuesOnGrid[1].getBz() * weightOnSecondPoint;
  double fieldBZ_2 = valuesOnGrid[2].getBz() * weightOnThirdPoint;
  double fieldBZ_3 = valuesOnGrid[3].getBz() * weightOnFourthPoint;

  tempField.setEx(fieldX_0 + fieldX_1 + fieldX_2 + fieldX_3);
  tempField.setBx(fieldBX_0 + fieldBX_1 + fieldBX_2 + fieldBX_3);
  tempField.setEy(fieldY_0 + fieldY_1 + fieldY_2 + fieldY_3);
  tempField.setBy(fieldBY_0 + fieldBY_1 + fieldBY_2 + fieldBY_3);
  tempField.setEz(fieldZ_0 + fieldZ_1 + fieldZ_2 + fieldZ_3);
  tempField.setBz(fieldBZ_0 + fieldBZ_1 + fieldBZ_2 + fieldBZ_3);

  return tempField;
}


Field interpolation2D_quadratic(Data data, Particle position, vector<Field> &valuesOnGrid)
{
  /*
  1 2 3    |--->  x
  4 5 6  |
  7 8 9  V
         y
  */

  Particle relativeCoord;

  double gridd;
  int grid = 0;
  double dx_inv = 1 / data.getDeltaX();
  double dy_inv = 1 / data.getDeltaY();

  gridd = (position.getParticleX() / data.getDeltaX()) + 0.5;
  if (gridd < 0) grid = ((int)grid) - 1;
  else grid = (int)gridd;
  relativeCoord.setParticleX((position.getParticleX() - data.getDeltaX()*grid) * dx_inv);

  gridd = (position.getParticleY() / data.getDeltaY()) + 0.5;
  if (gridd < 0) grid = ((int)grid) - 1;
  else grid = (int)gridd;
  relativeCoord.setParticleY((position.getParticleY() - data.getDeltaY()*grid) * dy_inv);

  relativeCoord.setParticleZ(0.);


  Field tempField;

  // DA SISTEMARE!!!
  double weightOnFirstPoint = 0.75 - pow(relativeCoord.getParticleX(), 2.);
  double weightOnSecondPoint = 0.5 * (0.25 + pow(relativeCoord.getParticleX(), 2.) + relativeCoord.getParticleX());
  double weightOnThirdPoint = 0.5 * (0.25 + pow(relativeCoord.getParticleX(), 2.) - relativeCoord.getParticleX());
  double weightOnFourthPoint = 0.;
  double weightOnFifthPoint = 0.;
  double weightOnSixthPoint = 0.;
  double weightOnSeventhPoint = 0.;
  double weightOnEighthPoint = 0.;
  double weightOnNinethPoint = 0.;

  double fieldX_0 = valuesOnGrid[0].getEx() * weightOnFirstPoint;
  double fieldX_1 = valuesOnGrid[1].getEx() * weightOnSecondPoint;
  double fieldX_2 = valuesOnGrid[2].getEx() * weightOnThirdPoint;
  double fieldX_3 = valuesOnGrid[3].getEx() * weightOnFourthPoint;
  double fieldX_4 = valuesOnGrid[4].getEx() * weightOnFifthPoint;
  double fieldX_5 = valuesOnGrid[5].getEx() * weightOnSixthPoint;
  double fieldX_6 = valuesOnGrid[6].getEx() * weightOnSeventhPoint;
  double fieldX_7 = valuesOnGrid[7].getEx() * weightOnEighthPoint;
  double fieldX_8 = valuesOnGrid[8].getEx() * weightOnNinethPoint;

  double fieldY_0 = valuesOnGrid[0].getEy() * weightOnFirstPoint;
  double fieldY_1 = valuesOnGrid[1].getEy() * weightOnSecondPoint;
  double fieldY_2 = valuesOnGrid[2].getEy() * weightOnThirdPoint;
  double fieldY_3 = valuesOnGrid[3].getEy() * weightOnFourthPoint;
  double fieldY_4 = valuesOnGrid[4].getEy() * weightOnFifthPoint;
  double fieldY_5 = valuesOnGrid[5].getEy() * weightOnSixthPoint;
  double fieldY_6 = valuesOnGrid[6].getEy() * weightOnSeventhPoint;
  double fieldY_7 = valuesOnGrid[7].getEy() * weightOnEighthPoint;
  double fieldY_8 = valuesOnGrid[8].getEy() * weightOnNinethPoint;

  double fieldZ_0 = valuesOnGrid[0].getEz() * weightOnFirstPoint;
  double fieldZ_1 = valuesOnGrid[1].getEz() * weightOnSecondPoint;
  double fieldZ_2 = valuesOnGrid[2].getEz() * weightOnThirdPoint;
  double fieldZ_3 = valuesOnGrid[3].getEz() * weightOnFourthPoint;
  double fieldZ_4 = valuesOnGrid[4].getEz() * weightOnFifthPoint;
  double fieldZ_5 = valuesOnGrid[5].getEz() * weightOnSixthPoint;
  double fieldZ_6 = valuesOnGrid[6].getEz() * weightOnSeventhPoint;
  double fieldZ_7 = valuesOnGrid[7].getEz() * weightOnEighthPoint;
  double fieldZ_8 = valuesOnGrid[8].getEz() * weightOnNinethPoint;

  double fieldBX_0 = valuesOnGrid[0].getBx() * weightOnFirstPoint;
  double fieldBX_1 = valuesOnGrid[1].getBx() * weightOnSecondPoint;
  double fieldBX_2 = valuesOnGrid[2].getBx() * weightOnThirdPoint;
  double fieldBX_3 = valuesOnGrid[3].getBx() * weightOnFourthPoint;
  double fieldBX_4 = valuesOnGrid[4].getBx() * weightOnFifthPoint;
  double fieldBX_5 = valuesOnGrid[5].getBx() * weightOnSixthPoint;
  double fieldBX_6 = valuesOnGrid[6].getBx() * weightOnSeventhPoint;
  double fieldBX_7 = valuesOnGrid[7].getBx() * weightOnEighthPoint;
  double fieldBX_8 = valuesOnGrid[8].getBx() * weightOnNinethPoint;

  double fieldBY_0 = valuesOnGrid[0].getBy() * weightOnFirstPoint;
  double fieldBY_1 = valuesOnGrid[1].getBy() * weightOnSecondPoint;
  double fieldBY_2 = valuesOnGrid[2].getBy() * weightOnThirdPoint;
  double fieldBY_3 = valuesOnGrid[3].getBy() * weightOnFourthPoint;
  double fieldBY_4 = valuesOnGrid[4].getBy() * weightOnFifthPoint;
  double fieldBY_5 = valuesOnGrid[5].getBy() * weightOnSixthPoint;
  double fieldBY_6 = valuesOnGrid[6].getBy() * weightOnSeventhPoint;
  double fieldBY_7 = valuesOnGrid[7].getBy() * weightOnEighthPoint;
  double fieldBY_8 = valuesOnGrid[8].getBy() * weightOnNinethPoint;

  double fieldBZ_0 = valuesOnGrid[0].getBz() * weightOnFirstPoint;
  double fieldBZ_1 = valuesOnGrid[1].getBz() * weightOnSecondPoint;
  double fieldBZ_2 = valuesOnGrid[2].getBz() * weightOnThirdPoint;
  double fieldBZ_3 = valuesOnGrid[3].getBz() * weightOnFourthPoint;
  double fieldBZ_4 = valuesOnGrid[4].getBz() * weightOnFifthPoint;
  double fieldBZ_5 = valuesOnGrid[5].getBz() * weightOnSixthPoint;
  double fieldBZ_6 = valuesOnGrid[6].getBz() * weightOnSeventhPoint;
  double fieldBZ_7 = valuesOnGrid[7].getBz() * weightOnEighthPoint;
  double fieldBZ_8 = valuesOnGrid[8].getBz() * weightOnNinethPoint;

  tempField.setEx(fieldX_0 + fieldX_1 + fieldX_2 + fieldX_3 + fieldX_4 + fieldX_5 + fieldX_6 + fieldX_7 + fieldX_8);
  tempField.setBx(fieldBX_0 + fieldBX_1 + fieldBX_2 + fieldBX_3 + fieldBX_4 + fieldBX_5 + fieldBX_6 + fieldBX_7 + fieldBX_8);
  tempField.setEy(fieldY_0 + fieldY_1 + fieldY_2 + fieldY_3 + fieldY_4 + fieldY_5 + fieldY_6 + fieldY_7 + fieldY_8);
  tempField.setBy(fieldBY_0 + fieldBY_1 + fieldBY_2 + fieldBY_3 + fieldBY_4 + fieldBY_5 + fieldBY_6 + fieldBY_7 + fieldBY_8);
  tempField.setEz(fieldZ_0 + fieldZ_1 + fieldZ_2 + fieldZ_3 + fieldZ_4 + fieldZ_5 + fieldZ_6 + fieldZ_7 + fieldZ_8);
  tempField.setBz(fieldBZ_0 + fieldBZ_1 + fieldBZ_2 + fieldBZ_3 + fieldBZ_4 + fieldBZ_5 + fieldBZ_6 + fieldBZ_7 + fieldBZ_8);

  return tempField;
}



void evolveLPF_nogrid(Data data, vector<Particle> &particles, vector<Field> &fieldOnParticles, ofstream &outputstream)
{
  Particle pstar, pHalfAdvanced;

  vector<Field> potenziale;
  Field b;
  Field *campoSuPunto;
  campoSuPunto = new Field[1];

  double gamma, gamma_inv, b2;
  double I1, I2, I3;
  double zeta;


  //ALGORITMO LEAPFROG PARTICELLE

  if (data.getNdim() == 1)
  {
    for (int i = 0; i < data.getNsteps(); i++)
    {
      for (int j = 0; j < data.getNelectrons(); j++)
      {
        gamma = sqrt(1. + pow(particles[j].getParticlePX(), 2) + pow(particles[j].getParticlePY(), 2) + pow(particles[j].getParticlePZ(), 2));
        gamma_inv = 1. / gamma;

        particles[j].setParticleX(particles[j].getParticleX() + 0.5 * data.getDeltaT() * particles[j].getParticlePX() * gamma_inv);
        particles[j].setParticleY(particles[j].getParticleY() + 0.5 * data.getDeltaT() * particles[j].getParticlePY() * gamma_inv);
        particles[j].setParticleZ(particles[j].getParticleZ() + 0.5 * data.getDeltaT() * particles[j].getParticlePZ() * gamma_inv);

        fieldOnParticles.clear();                         //SVUOTA il vettore dei campi sulle particelle ad ogni ciclo

        //SICCOME NON HO ANCORA UN CAMPO SU GRIGLIA DEFINITO, TEMPORANEAMENTE CALCOLO IL CAMPO ANALITICO
        //DIRETTAMENTE SULLA PARTICELLA
        calcolaCampoAnalitico1DSuParticella(data, particles[j], *campoSuPunto);
        fieldOnParticles.push_back(*campoSuPunto);


        //pstar = p(n-1/2) + Q * E(n) * DeltaT/2
        pstar.setParticlePX(particles[j].getParticlePX() + 0.5*data.getDeltaT() * particles[j].getParticleQ() * fieldOnParticles[0].getEx());
        pstar.setParticlePY(particles[j].getParticlePY() + 0.5*data.getDeltaT() * particles[j].getParticleQ() * fieldOnParticles[0].getEy());
        pstar.setParticlePZ(particles[j].getParticlePZ() + 0.5*data.getDeltaT() * particles[j].getParticleQ() * fieldOnParticles[0].getEz());
        gamma = sqrt(1. + pow(pstar.getParticlePX(), 2) + pow(pstar.getParticlePY(), 2) + pow(pstar.getParticlePZ(), 2));
        gamma_inv = 1. / gamma;

        //b = (DeltaT/2) * Q * B(n) / gamma
        b.setBx(0.5*data.getDeltaT() * fieldOnParticles[0].getBx() * particles[j].getParticleQ() * gamma_inv);
        b.setBy(0.5*data.getDeltaT() * fieldOnParticles[0].getBy() * particles[j].getParticleQ() * gamma_inv);
        b.setBz(0.5*data.getDeltaT() * fieldOnParticles[0].getBz() * particles[j].getParticleQ() * gamma_inv);
        b2 = pow(b.getBx(), 2) + pow(b.getBy(), 2) + pow(b.getBz(), 2);
        b2 += 1.;

        /****************************************************************************
        //p(n+1/2) = [pstar + pstar x b + (pstar x b) * b] / b2   SECONDO TURCHETTI  (ma i conti analitici diconon che e` sbagliato, anche perche' l'ultimo termine e` uno scalare!)
        pHalfAdvanced.setParticlePX( (pstar.getParticlePX() + pstar.getParticlePY() * b.getBz() - b.getBy() * pstar.getParticlePZ() + (pstar.getParticlePY() * b.getBz() - b.getBy() * pstar.getParticlePZ()) * b.getBx() ) / b2 );
        pHalfAdvanced.setParticlePY( (pstar.getParticlePY() + pstar.getParticlePZ() * b.getBx() - b.getBz() * pstar.getParticlePX() + (pstar.getParticlePZ() * b.getBx() - b.getBz() * pstar.getParticlePX()) * b.getBy() ) / b2 );
        pHalfAdvanced.setParticlePZ( (pstar.getParticlePZ() + pstar.getParticlePX() * b.getBy() - b.getBx() * pstar.getParticlePY() + (pstar.getParticlePX() * b.getBy() - b.getBx() * pstar.getParticlePY()) * b.getBz() ) / b2 );
        ****************************************************************************/

        /****************************************************************************
        //p(n+1/2) = [pstar + pstar x b + b x (pstar x b)] / b2   SECONDO LONDRILLO
        pHalfAdvanced.setParticlePX( (pstar.getParticlePX() + pstar.getParticlePY() * b.getBz() - b.getBy() * pstar.getParticlePZ() + b.getBx() * pstar.getParticlePX() * b.getBx() ) / b2 );
        pHalfAdvanced.setParticlePY( (pstar.getParticlePY() + pstar.getParticlePZ() * b.getBx() - b.getBz() * pstar.getParticlePX() + b.getBy() * pstar.getParticlePY() * b.getBy() ) / b2 );
        pHalfAdvanced.setParticlePZ( (pstar.getParticlePZ() + pstar.getParticlePX() * b.getBy() - b.getBx() * pstar.getParticlePY() + b.getBz() * pstar.getParticlePZ() * b.getBz() ) / b2 );
        ****************************************************************************/

        //      /****************************************************************************
              //p(n+1/2) = [pstar + pstar x b] / b2   COME IMPLEMENTATO DA TURCHETTI NEL CODICE
        pHalfAdvanced.setParticlePX((pstar.getParticlePX() + pstar.getParticlePY() * b.getBz()) / b2);
        pHalfAdvanced.setParticlePY((pstar.getParticlePY() - pstar.getParticlePX() * b.getBz()) / b2);
        pHalfAdvanced.setParticlePZ(0.);  //per come ho strutturato i pacchetti d'onda su z non dovrebbe succedere nulla e quindi lo fisso io a zero.
  //      ****************************************************************************/

        particles[j].setParticlePX(2.* pHalfAdvanced.getParticlePX() - particles[j].getParticlePX());
        particles[j].setParticlePY(2.* pHalfAdvanced.getParticlePY() - particles[j].getParticlePY());
        particles[j].setParticlePZ(2.* pHalfAdvanced.getParticlePZ() - particles[j].getParticlePZ());

        gamma = sqrt(1. + pow(particles[j].getParticlePX(), 2) + pow(particles[j].getParticlePY(), 2) + pow(particles[j].getParticlePZ(), 2));
        gamma_inv = 1. / gamma;

        particles[j].setParticleX(particles[j].getParticleX() + 0.5 * data.getDeltaT() * particles[j].getParticlePX() * gamma_inv);
        particles[j].setParticleY(particles[j].getParticleY() + 0.5 * data.getDeltaT() * particles[j].getParticlePY() * gamma_inv);
        particles[j].setParticleZ(particles[j].getParticleZ() + 0.5 * data.getDeltaT() * particles[j].getParticlePZ() * gamma_inv);
        particles[j].setParticleCell(data);


        //INTEGRALI PRIMI (da sistemare, li ho copiati dall'evoluzione sui campi analitici)
        potenziale.clear();
        zeta = particles[j].getParticleX() - C * data.getT();             // non usato nella mia struttura delle equazioni di campo
        calcolaPotenzialeAnalitico1DSuParticella(data, particles[j], *campoSuPunto);
        potenziale.push_back(*campoSuPunto);

        I1 = (particles[j].getParticlePY() / (particles[j].getParticleM()*C)) + potenziale[0].getEy();
        I2 = particles[j].getParticlePZ() / (particles[j].getParticleM()*C);
        I3 = gamma - particles[j].getParticlePX() / (particles[j].getParticleM()*C);


        outputstream
          << setw(12) << i + 1 << "\t"
          << setprecision(14) << setw(16) << particles[j].getParticleX() << "\t"
          << setw(16) << particles[j].getParticleY() << "\t"
          << setw(16) << particles[j].getParticleZ() << "\t"
          //        << setw(16) << particles[j].getParticleCell() << "\t"
          << setw(16) << particles[j].getParticlePX() << "\t"
          << setw(16) << particles[j].getParticlePY() << "\t"
          << setw(16) << particles[j].getParticlePZ() << "\t"
          << setw(16) << fieldOnParticles[j].getEy() << "\t"
          << setw(16) << zeta << "\t"
          << setw(16) << I1 << "\t"
          << setw(16) << I2 << "\t"
          << setw(16) << I3 << endl;
      }
      data.aumentaT(data.getDeltaT());
    }
  }
}


void evolveLPF_withgrid_onthefly(Data data, vector<Particle> &particles, vector<Field> &fieldOnParticles, ofstream &outputstream)
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

  Field *FOG;                             //DA SISTEMARE, c'è solo per avere un campo fittizio su una griglia on-the-fly per implementare l'algoritmo
  FOG = new Field[1];                         //-------

  double gridXX, gridYY, gridZZ;
  int gridX, gridY, gridZ;
  double gamma, gamma_inv, b2;
  double I1, I2, I3;
  double zeta;

  Field *campoSuPunto;
  campoSuPunto = new Field[1];

  //ALGORITMO LEAPFROG PARTICELLE

  if (data.getNdim() == 1)
  {
    for (int i = 0; i < data.getNsteps(); i++)
    {

      data.aumentaT(data.getDeltaT());

      for (int j = 0; j < data.getNelectrons(); j++)
      {
        gamma = sqrt(1. + pow(particles[j].getParticlePX(), 2) + pow(particles[j].getParticlePY(), 2) + pow(particles[j].getParticlePZ(), 2));
        gamma_inv = 1. / gamma;

        particles[j].setParticleX(particles[j].getParticleX() + 0.5 * data.getDeltaT() * particles[j].getParticlePX() * gamma_inv);
        particles[j].setParticleY(particles[j].getParticleY() + 0.5 * data.getDeltaT() * particles[j].getParticlePY() * gamma_inv);
        particles[j].setParticleZ(particles[j].getParticleZ() + 0.5 * data.getDeltaT() * particles[j].getParticlePZ() * gamma_inv);


        fieldOnParticles.clear();

        //QUI MI CALCOLO I PUNTI DI RIFERIMENTO PER LA PARTICELLA IN ESAME nella posizione attuale
        //forse c'e` un modo piu` furbo per farlo sfruttando il # di cella in cui si trova che e` noto
        gridXX = particles[j].getParticleX() / data.getDeltaX();
        if (gridXX < 0) gridX = ((int)gridXX) - 1;           //il punto di griglia lungo l'asse x a sinistra della particella
        else gridX = (int)gridXX;

        gridYY = particles[j].getParticleY() / data.getDeltaY();
        if (gridYY < 0) gridY = ((int)gridYY) - 1;           //il punto di griglia lungo l'asse y a sinistra della particella
        else gridY = (int)gridYY;

        gridZZ = particles[j].getParticleZ() / data.getDeltaZ();
        if (gridZZ < 0) gridZ = ((int)gridZZ) - 1;           //il punto di griglia lungo l'asse z a sinistra della particella
        else gridZ = (int)gridZZ;

        coordinatePrimoPuntoGriglia.setParticleX(gridX*data.getDeltaX());
        coordinatePrimoPuntoGriglia.setParticleY(gridY*data.getDeltaY());
        coordinatePrimoPuntoGriglia.setParticleZ(gridZ*data.getDeltaZ());
        coordinateSecondoPuntoGriglia.setParticleX(gridX*data.getDeltaX() + data.getDeltaX());
        coordinateSecondoPuntoGriglia.setParticleY(gridY*data.getDeltaY());
        coordinateSecondoPuntoGriglia.setParticleZ(gridZ*data.getDeltaZ());
        coordinateTerzoPuntoGriglia.setParticleX(gridX*data.getDeltaX());
        coordinateTerzoPuntoGriglia.setParticleY(gridY*data.getDeltaY() + data.getDeltaY());
        coordinateTerzoPuntoGriglia.setParticleZ(gridZ*data.getDeltaZ());
        coordinateQuartoPuntoGriglia.setParticleX(gridX*data.getDeltaX() + data.getDeltaX());
        coordinateQuartoPuntoGriglia.setParticleY(gridY*data.getDeltaY() + data.getDeltaY());
        coordinateQuartoPuntoGriglia.setParticleZ(gridZ*data.getDeltaZ());
        coordinateQuintoPuntoGriglia.setParticleX(gridX*data.getDeltaX());
        coordinateQuintoPuntoGriglia.setParticleY(gridY*data.getDeltaY());
        coordinateQuintoPuntoGriglia.setParticleZ(gridZ*data.getDeltaZ() + data.getDeltaZ());
        coordinateSestoPuntoGriglia.setParticleX(gridX*data.getDeltaX() + data.getDeltaX());
        coordinateSestoPuntoGriglia.setParticleY(gridY*data.getDeltaY());
        coordinateSestoPuntoGriglia.setParticleZ(gridZ*data.getDeltaZ() + data.getDeltaZ());
        coordinateSettimoPuntoGriglia.setParticleX(gridX*data.getDeltaX());
        coordinateSettimoPuntoGriglia.setParticleY(gridY*data.getDeltaY() + data.getDeltaY());
        coordinateSettimoPuntoGriglia.setParticleZ(gridZ*data.getDeltaZ() + data.getDeltaZ());
        coordinateOttavoPuntoGriglia.setParticleX(gridX*data.getDeltaX() + data.getDeltaX());
        coordinateOttavoPuntoGriglia.setParticleY(gridY*data.getDeltaY() + data.getDeltaY());
        coordinateOttavoPuntoGriglia.setParticleZ(gridZ*data.getDeltaZ() + data.getDeltaZ());


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
        pstar.setParticlePX(particles[j].getParticlePX() + 0.5 * data.getDeltaT() * particles[j].getParticleQ() * fieldOnParticles[0].getEx());
        pstar.setParticlePY(particles[j].getParticlePY() + 0.5 * data.getDeltaT() * particles[j].getParticleQ() * fieldOnParticles[0].getEy());
        pstar.setParticlePZ(particles[j].getParticlePZ() + 0.5 * data.getDeltaT() * particles[j].getParticleQ() * fieldOnParticles[0].getEz());
        gamma = sqrt(1. + pow(pstar.getParticlePX(), 2) + pow(pstar.getParticlePY(), 2) + pow(pstar.getParticlePZ(), 2));
        gamma_inv = 1. / gamma;

        //b = (DeltaT/2) * Q * B(n) / gamma
        b.setBx(0.5 * data.getDeltaT() * fieldOnParticles[0].getBx() * particles[j].getParticleQ() * gamma_inv);
        b.setBy(0.5 * data.getDeltaT() * fieldOnParticles[0].getBy() * particles[j].getParticleQ() * gamma_inv);
        b.setBz(0.5 * data.getDeltaT() * fieldOnParticles[0].getBz() * particles[j].getParticleQ() * gamma_inv);
        b2 = pow(b.getBx(), 2) + pow(b.getBy(), 2) + pow(b.getBz(), 2);
        b2 += 1.;

        /****************************************************************************
        //p(n+1/2) = [pstar + pstar x b + (pstar x b) * b] / b2   SECONDO TURCHETTI  (ma i conti analitici diconon che e` sbagliato, anche perche' l'ultimo termine e` uno scalare!)
        pHalfAdvanced.setParticlePX( (pstar.getParticlePX() + pstar.getParticlePY() * b.getBz() - b.getBy() * pstar.getParticlePZ() + (pstar.getParticlePY() * b.getBz() - b.getBy() * pstar.getParticlePZ()) * b.getBx() ) / b2 );
        pHalfAdvanced.setParticlePY( (pstar.getParticlePY() + pstar.getParticlePZ() * b.getBx() - b.getBz() * pstar.getParticlePX() + (pstar.getParticlePZ() * b.getBx() - b.getBz() * pstar.getParticlePX()) * b.getBy() ) / b2 );
        pHalfAdvanced.setParticlePZ( (pstar.getParticlePZ() + pstar.getParticlePX() * b.getBy() - b.getBx() * pstar.getParticlePY() + (pstar.getParticlePX() * b.getBy() - b.getBx() * pstar.getParticlePY()) * b.getBz() ) / b2 );
        ****************************************************************************/

        /****************************************************************************
        //p(n+1/2) = [pstar + pstar x b + b x (pstar x b)] / b2   SECONDO LONDRILLO
        pHalfAdvanced.setParticlePX( (pstar.getParticlePX() + pstar.getParticlePY() * b.getBz() - b.getBy() * pstar.getParticlePZ() + b.getBx() * pstar.getParticlePX() * b.getBx() ) / b2 );
        pHalfAdvanced.setParticlePY( (pstar.getParticlePY() + pstar.getParticlePZ() * b.getBx() - b.getBz() * pstar.getParticlePX() + b.getBy() * pstar.getParticlePY() * b.getBy() ) / b2 );
        pHalfAdvanced.setParticlePZ( (pstar.getParticlePZ() + pstar.getParticlePX() * b.getBy() - b.getBx() * pstar.getParticlePY() + b.getBz() * pstar.getParticlePZ() * b.getBz() ) / b2 );
        ****************************************************************************/

        // /****************************************************************************
        //p(n+1/2) = [pstar + pstar x b] / b2   COME IMPLEMENTATO DA TURCHETTI NEL CODICE
        pHalfAdvanced.setParticlePX((pstar.getParticlePX() + pstar.getParticlePY() * b.getBz()) / b2);
        pHalfAdvanced.setParticlePY((pstar.getParticlePY() - pstar.getParticlePX() * b.getBz()) / b2);
        pHalfAdvanced.setParticlePZ(0.);  //per come ho strutturato i pacchetti d'onda su z non dovrebbe succedere nulla e quindi lo fisso io a zero.
        // ****************************************************************************/

        particles[j].setParticlePX(2. * pHalfAdvanced.getParticlePX() - particles[j].getParticlePX());
        particles[j].setParticlePY(2. * pHalfAdvanced.getParticlePY() - particles[j].getParticlePY());
        particles[j].setParticlePZ(2. * pHalfAdvanced.getParticlePZ() - particles[j].getParticlePZ());


        gamma = sqrt(1. + pow(particles[j].getParticlePX(), 2) + pow(particles[j].getParticlePY(), 2) + pow(particles[j].getParticlePZ(), 2));
        gamma_inv = 1. / gamma;

        particles[j].setParticleX(particles[j].getParticleX() + 0.5 * data.getDeltaT() * particles[j].getParticlePX() * gamma_inv);
        particles[j].setParticleY(particles[j].getParticleY() + 0.5 * data.getDeltaT() * particles[j].getParticlePY() * gamma_inv);
        particles[j].setParticleZ(particles[j].getParticleZ() + 0.5 * data.getDeltaT() * particles[j].getParticlePZ() * gamma_inv);
        particles[j].setParticleCell(data);


        //INTEGRALI PRIMI (da sistemare, li ho copiati dall'evoluzione sui campi analitici)
        potenziale.clear();
        zeta = particles[j].getParticleX() - C * data.getT();             // non usato nella mia struttura delle equazioni di campo

        calcolaPotenzialeAnalitico1DSuParticella(data, particles[j], *campoSuPunto);
        potenziale.push_back(*campoSuPunto);

        I1 = (particles[j].getParticlePY() / (particles[j].getParticleM()*C)) + potenziale[0].getEy();
        I2 = particles[j].getParticlePZ() / (particles[j].getParticleM()*C);
        I3 = gamma - particles[j].getParticlePX() / (particles[j].getParticleM()*C);

        outputstream
          << setw(12) << i + 1 << "\t"
          << setprecision(14) << setw(16) << particles[j].getParticleX() << "\t"
          << setw(16) << particles[j].getParticleY() << "\t"
          << setw(16) << particles[j].getParticleZ() << "\t"
          // << setw(16) << particles[j].getParticleCell() << "\t"
          << setw(16) << particles[j].getParticlePX() << "\t"
          << setw(16) << particles[j].getParticlePY() << "\t"
          << setw(16) << particles[j].getParticlePZ() << "\t"
          << setw(16) << fieldOnParticles[0].getEy() << "\t"
          << setw(16) << potenziale[0].getEy() << "\t"
          << setw(16) << zeta << "\t"
          << setw(16) << I1 << "\t"
          << setw(16) << I2 << "\t"
          << setw(16) << I3 << endl;
      }
    }
  }
}


void evolveLPF_withgrid(Data data, vector<Particle> &particles, vector<Field> &fieldOnParticles, vector<Field> &campoSuPuntiGriglia, ofstream &outputstream)
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

  Field *campoSuPunto;
  campoSuPunto = new Field[1];

  //ALGORITMO LEAPFROG PARTICELLE

  if (data.getNdim() == 1)
  {
    for (int i = 0; i < data.getNsteps(); i++)
    {
      data.aumentaT(data.getDeltaT());

      campoSuPuntiGriglia.clear();
      riempiPuntiGrigliaConCampoAnalitico1D(data, campoSuPuntiGriglia);

      potenzialeSuPuntiGriglia.clear();
      riempiPuntiGrigliaConPotenzialeAnalitico1D(data, potenzialeSuPuntiGriglia);

      for (int j = 0; j < data.getNelectrons(); j++)
      {
        gamma = sqrt(1. + pow(particles[j].getParticlePX(), 2) + pow(particles[j].getParticlePY(), 2) + pow(particles[j].getParticlePZ(), 2));
        gamma_inv = 1. / gamma;

        particles[j].setParticleX(particles[j].getParticleX() + 0.5 * data.getDeltaT() * particles[j].getParticlePX() * gamma_inv);
        particles[j].setParticleY(particles[j].getParticleY() + 0.5 * data.getDeltaT() * particles[j].getParticlePY() * gamma_inv);
        particles[j].setParticleZ(particles[j].getParticleZ() + 0.5 * data.getDeltaT() * particles[j].getParticlePZ() * gamma_inv);


        // QUI MI CALCOLO I PUNTI DI RIFERIMENTO PER LA PARTICELLA IN ESAME nella posizione attuale
        // forse c'e` un modo piu` furbo per farlo sfruttando il # di cella in cui si trova che e` noto

        fieldOnParticles.clear();       //SVUOTA il vettore dei campi sulle particelle ad ogni ciclo


        gridXX = particles[j].getParticleX() / data.getDeltaX();
        if (gridXX < 0) gridX = ((int)gridXX) - 1;           //il punto di griglia lungo l'asse x a sinistra della particella
        else gridX = (int)gridXX;                      //

        gridYY = particles[j].getParticleY() / data.getDeltaY();
        if (gridYY < 0) gridY = ((int)gridYY) - 1;           //il punto di griglia lungo l'asse y a sinistra della particella
        else gridY = (int)gridYY;                      //

        gridZZ = particles[j].getParticleZ() / data.getDeltaZ();
        if (gridZZ < 0) gridZ = ((int)gridZZ) - 1;           //il punto di griglia lungo l'asse z a sinistra della particella
        else gridZ = (int)gridZZ;                      //

        if (gridX < 0 || gridY < 0 || gridZ < 0)
        {
          cout << "La particella e` uscita dai boundaries - simulazione interrotta" << endl;
          return;
        }


        fieldOnGrid.push_back(campoSuPuntiGriglia.at(gridY * data.getNgridPointsX() + gridX));        // primo punto
        fieldOnGrid.push_back(campoSuPuntiGriglia.at(gridY * data.getNgridPointsX() + (gridX + 1)));      // secondo punto
        fieldOnGrid.push_back(campoSuPuntiGriglia.at((gridY + 1) * data.getNgridPointsX() + gridX));      // terzo punto
        fieldOnGrid.push_back(campoSuPuntiGriglia.at((gridY + 1) * data.getNgridPointsX() + (gridX + 1)));    // quarto punto

        FOP = interpolation2D_linear(data, particles[j], fieldOnGrid);
        fieldOnParticles.push_back(FOP);

        fieldOnGrid.clear();



        //pstar = p(n-1/2) + Q * E(n) * DeltaT/2
        pstar.setParticlePX(particles[j].getParticlePX() + 0.5 * data.getDeltaT() * particles[j].getParticleQ() * fieldOnParticles[0].getEx());
        pstar.setParticlePY(particles[j].getParticlePY() + 0.5 * data.getDeltaT() * particles[j].getParticleQ() * fieldOnParticles[0].getEy());
        pstar.setParticlePZ(particles[j].getParticlePZ() + 0.5 * data.getDeltaT() * particles[j].getParticleQ() * fieldOnParticles[0].getEz());
        gamma = sqrt(1. + pow(pstar.getParticlePX(), 2) + pow(pstar.getParticlePY(), 2) + pow(pstar.getParticlePZ(), 2));
        gamma_inv = 1. / gamma;

        //b = (DeltaT/2) * Q * B(n) / gamma
        b.setBx(0.5 * data.getDeltaT() * fieldOnParticles[0].getBx() * particles[j].getParticleQ() * gamma_inv);
        b.setBy(0.5 * data.getDeltaT() * fieldOnParticles[0].getBy() * particles[j].getParticleQ() * gamma_inv);
        b.setBz(0.5 * data.getDeltaT() * fieldOnParticles[0].getBz() * particles[j].getParticleQ() * gamma_inv);
        b2 = pow(b.getBx(), 2) + pow(b.getBy(), 2) + pow(b.getBz(), 2);
        b2 += 1.;

        /****************************************************************************
        //p(n+1/2) = [pstar + pstar x b + (pstar x b) * b] / b2   SECONDO TURCHETTI  (ma i conti analitici diconon che e` sbagliato, anche perche' l'ultimo termine e` uno scalare!)
        pHalfAdvanced.setParticlePX( (pstar.getParticlePX() + pstar.getParticlePY() * b.getBz() - b.getBy() * pstar.getParticlePZ() + (pstar.getParticlePY() * b.getBz() - b.getBy() * pstar.getParticlePZ()) * b.getBx() ) / b2 );
        pHalfAdvanced.setParticlePY( (pstar.getParticlePY() + pstar.getParticlePZ() * b.getBx() - b.getBz() * pstar.getParticlePX() + (pstar.getParticlePZ() * b.getBx() - b.getBz() * pstar.getParticlePX()) * b.getBy() ) / b2 );
        pHalfAdvanced.setParticlePZ( (pstar.getParticlePZ() + pstar.getParticlePX() * b.getBy() - b.getBx() * pstar.getParticlePY() + (pstar.getParticlePX() * b.getBy() - b.getBx() * pstar.getParticlePY()) * b.getBz() ) / b2 );
        ****************************************************************************/

        /****************************************************************************
        //p(n+1/2) = [pstar + pstar x b + b x (pstar x b)] / b2   SECONDO LONDRILLO
        pHalfAdvanced.setParticlePX( (pstar.getParticlePX() + pstar.getParticlePY() * b.getBz() - b.getBy() * pstar.getParticlePZ() + b.getBx() * pstar.getParticlePX() * b.getBx() ) / b2 );
        pHalfAdvanced.setParticlePY( (pstar.getParticlePY() + pstar.getParticlePZ() * b.getBx() - b.getBz() * pstar.getParticlePX() + b.getBy() * pstar.getParticlePY() * b.getBy() ) / b2 );
        pHalfAdvanced.setParticlePZ( (pstar.getParticlePZ() + pstar.getParticlePX() * b.getBy() - b.getBx() * pstar.getParticlePY() + b.getBz() * pstar.getParticlePZ() * b.getBz() ) / b2 );
        ****************************************************************************/

        //      /****************************************************************************
              //p(n+1/2) = [pstar + pstar x b] / b2   COME IMPLEMENTATO DA TURCHETTI NEL CODICE
        pHalfAdvanced.setParticlePX((pstar.getParticlePX() + pstar.getParticlePY() * b.getBz()) / b2);
        pHalfAdvanced.setParticlePY((pstar.getParticlePY() - pstar.getParticlePX() * b.getBz()) / b2);
        pHalfAdvanced.setParticlePZ(0.);  //per come ho strutturato i pacchetti d'onda su z non dovrebbe succedere nulla e quindi lo fisso io a zero.
  //      ****************************************************************************/

        particles[j].setParticlePX(2. * pHalfAdvanced.getParticlePX() - particles[j].getParticlePX());
        particles[j].setParticlePY(2. * pHalfAdvanced.getParticlePY() - particles[j].getParticlePY());
        particles[j].setParticlePZ(2. * pHalfAdvanced.getParticlePZ() - particles[j].getParticlePZ());

        gamma = sqrt(1. + pow(particles[j].getParticlePX(), 2) + pow(particles[j].getParticlePY(), 2) + pow(particles[j].getParticlePZ(), 2));
        gamma_inv = 1. / gamma;

        particles[j].setParticleX(particles[j].getParticleX() + 0.5 * data.getDeltaT() * particles[j].getParticlePX() * gamma_inv);
        particles[j].setParticleY(particles[j].getParticleY() + 0.5 * data.getDeltaT() * particles[j].getParticlePY() * gamma_inv);
        particles[j].setParticleZ(particles[j].getParticleZ() + 0.5 * data.getDeltaT() * particles[j].getParticlePZ() * gamma_inv);
        particles[j].setParticleCell(data);


        //INTEGRALI PRIMI (da sistemare, li ho copiati dall'evoluzione sui campi analitici)
        potenziale.clear();
        zeta = particles[j].getParticleX() - C * data.getT();             // non usato nella mia struttura delle equazioni di campo

        gridXX = particles[j].getParticleX() / data.getDeltaX();
        if (gridXX < 0) gridX = ((int)gridXX) - 1;           //il punto di griglia lungo l'asse x a sinistra della particella
        else gridX = (int)gridXX;                      //

        gridYY = particles[j].getParticleY() / data.getDeltaY();
        if (gridYY < 0) gridY = ((int)gridYY) - 1;           //il punto di griglia lungo l'asse y a sinistra della particella
        else gridY = (int)gridYY;                      //

        gridZZ = particles[j].getParticleZ() / data.getDeltaZ();
        if (gridZZ < 0) gridZ = ((int)gridZZ) - 1;           //il punto di griglia lungo l'asse z a sinistra della particella
        else gridZ = (int)gridZZ;                      //

        if (gridX < 0 || gridY < 0 || gridZ < 0)
        {
          cout << "La particella e` uscita dai boundaries - simulazione interrotta" << endl;
          return;
        }

        fieldOnGrid.clear();

        fieldOnGrid.push_back(potenzialeSuPuntiGriglia.at(gridY * data.getNgridPointsX() + gridX));       // primo punto
        fieldOnGrid.push_back(potenzialeSuPuntiGriglia.at(gridY * data.getNgridPointsX() + (gridX + 1)));     // secondo punto
        fieldOnGrid.push_back(potenzialeSuPuntiGriglia.at((gridY + 1) * data.getNgridPointsX() + gridX));     // terzo punto
        fieldOnGrid.push_back(potenzialeSuPuntiGriglia.at((gridY + 1) * data.getNgridPointsX() + (gridX + 1)));   // quarto punto

        FOP = interpolation2D_linear(data, particles[j], fieldOnGrid);
        potenziale.push_back(FOP);

        I1 = (particles[j].getParticlePY() / (particles[j].getParticleM()*C)) + potenziale[0].getEy();
        I2 = particles[j].getParticlePZ() / (particles[j].getParticleM()*C);
        I3 = gamma - particles[j].getParticlePX() / (particles[j].getParticleM()*C);



        outputstream
          << setw(12) << i + 1 << "\t"
          << setprecision(14) << setw(16) << particles[j].getParticleX() << "\t"
          << setw(16) << particles[j].getParticleY() << "\t"
          << setw(16) << particles[j].getParticleZ() << "\t"
          //        << setw(16) << particles[j].getParticleCell() << "\t"
          << setw(16) << particles[j].getParticlePX() << "\t"
          << setw(16) << particles[j].getParticlePY() << "\t"
          << setw(16) << particles[j].getParticlePZ() << "\t"
          << setw(16) << potenziale[0].getEy() << "\t"
          << setw(16) << zeta << "\t"
          << setw(16) << I1 << "\t"
          << setw(16) << I2 << "\t"
          << setw(16) << I3 << endl;


        /*
        //UPDATE CAMPI ATTORNO ALLA NUOVA POSIZIONE DELLA PARTICELLA  (al momento non li evolve, il set li pone pari al get precedente)
        gridX = (int) (particles[j].getParticleX()/data.getDeltaX());   //il punto di griglia lungo l'asse x piu` vicino all'origine
        gridY = (int) (particles[j].getParticleY()/data.getDeltaY());   //il punto di griglia lungo l'asse y piu` vicino all'origine
        gridZ = (int) (particles[j].getParticleZ()/data.getDeltaZ());   //il punto di griglia lungo l'asse z piu` vicino all'origine

        campoSuPuntiGriglia.at(gridZ * data.getNgridPointsX() * data.getNgridPointsY() + gridY * data.getNgridPointsX() + gridX).setEx( campoSuPuntiGriglia.at(gridZ * data.getNgridPointsX() * data.getNgridPointsY() + gridY * data.getNgridPointsX() + gridX).getEx() );       // primo punto
        campoSuPuntiGriglia.at(gridZ * data.getNgridPointsX() * data.getNgridPointsY() + gridY * data.getNgridPointsX() + gridX).setEy( campoSuPuntiGriglia.at(gridZ * data.getNgridPointsX() * data.getNgridPointsY() + gridY * data.getNgridPointsX() + gridX).getEy() );
        campoSuPuntiGriglia.at(gridZ * data.getNgridPointsX() * data.getNgridPointsY() + gridY * data.getNgridPointsX() + gridX).setEz( campoSuPuntiGriglia.at(gridZ * data.getNgridPointsX() * data.getNgridPointsY() + gridY * data.getNgridPointsX() + gridX).getEz() );
        campoSuPuntiGriglia.at(gridZ * data.getNgridPointsX() * data.getNgridPointsY() + gridY * data.getNgridPointsX() + gridX).setBx( campoSuPuntiGriglia.at(gridZ * data.getNgridPointsX() * data.getNgridPointsY() + gridY * data.getNgridPointsX() + gridX).getBx() );
        campoSuPuntiGriglia.at(gridZ * data.getNgridPointsX() * data.getNgridPointsY() + gridY * data.getNgridPointsX() + gridX).setBy( campoSuPuntiGriglia.at(gridZ * data.getNgridPointsX() * data.getNgridPointsY() + gridY * data.getNgridPointsX() + gridX).getBy() );
        campoSuPuntiGriglia.at(gridZ * data.getNgridPointsX() * data.getNgridPointsY() + gridY * data.getNgridPointsX() + gridX).setBz( campoSuPuntiGriglia.at(gridZ * data.getNgridPointsX() * data.getNgridPointsY() + gridY * data.getNgridPointsX() + gridX).getBz() );

        campoSuPuntiGriglia.at(gridZ * data.getNgridPointsX() * data.getNgridPointsY() + gridY * data.getNgridPointsX() + (gridX+1)).setEx( campoSuPuntiGriglia.at(gridZ * data.getNgridPointsX() * data.getNgridPointsY() + gridY * data.getNgridPointsX() + (gridX+1)).getEx() );   // secondo punto
        campoSuPuntiGriglia.at(gridZ * data.getNgridPointsX() * data.getNgridPointsY() + gridY * data.getNgridPointsX() + (gridX+1)).setEy( campoSuPuntiGriglia.at(gridZ * data.getNgridPointsX() * data.getNgridPointsY() + gridY * data.getNgridPointsX() + (gridX+1)).getEy() );
        campoSuPuntiGriglia.at(gridZ * data.getNgridPointsX() * data.getNgridPointsY() + gridY * data.getNgridPointsX() + (gridX+1)).setEz( campoSuPuntiGriglia.at(gridZ * data.getNgridPointsX() * data.getNgridPointsY() + gridY * data.getNgridPointsX() + (gridX+1)).getEz() );
        campoSuPuntiGriglia.at(gridZ * data.getNgridPointsX() * data.getNgridPointsY() + gridY * data.getNgridPointsX() + (gridX+1)).setBx( campoSuPuntiGriglia.at(gridZ * data.getNgridPointsX() * data.getNgridPointsY() + gridY * data.getNgridPointsX() + (gridX+1)).getBx() );
        campoSuPuntiGriglia.at(gridZ * data.getNgridPointsX() * data.getNgridPointsY() + gridY * data.getNgridPointsX() + (gridX+1)).setBy( campoSuPuntiGriglia.at(gridZ * data.getNgridPointsX() * data.getNgridPointsY() + gridY * data.getNgridPointsX() + (gridX+1)).getBy() );
        campoSuPuntiGriglia.at(gridZ * data.getNgridPointsX() * data.getNgridPointsY() + gridY * data.getNgridPointsX() + (gridX+1)).setBz( campoSuPuntiGriglia.at(gridZ * data.getNgridPointsX() * data.getNgridPointsY() + gridY * data.getNgridPointsX() + (gridX+1)).getBz() );

        campoSuPuntiGriglia.at(gridZ * data.getNgridPointsX() * data.getNgridPointsY() + (gridY+1) * data.getNgridPointsX() + gridX).setEx( campoSuPuntiGriglia.at(gridZ * data.getNgridPointsX() * data.getNgridPointsY() + (gridY+1) * data.getNgridPointsX() + gridX).getEx() );   // terzo punto
        campoSuPuntiGriglia.at(gridZ * data.getNgridPointsX() * data.getNgridPointsY() + (gridY+1) * data.getNgridPointsX() + gridX).setEy( campoSuPuntiGriglia.at(gridZ * data.getNgridPointsX() * data.getNgridPointsY() + (gridY+1) * data.getNgridPointsX() + gridX).getEy() );
        campoSuPuntiGriglia.at(gridZ * data.getNgridPointsX() * data.getNgridPointsY() + (gridY+1) * data.getNgridPointsX() + gridX).setEz( campoSuPuntiGriglia.at(gridZ * data.getNgridPointsX() * data.getNgridPointsY() + (gridY+1) * data.getNgridPointsX() + gridX).getEz() );
        campoSuPuntiGriglia.at(gridZ * data.getNgridPointsX() * data.getNgridPointsY() + (gridY+1) * data.getNgridPointsX() + gridX).setBx( campoSuPuntiGriglia.at(gridZ * data.getNgridPointsX() * data.getNgridPointsY() + (gridY+1) * data.getNgridPointsX() + gridX).getBx() );
        campoSuPuntiGriglia.at(gridZ * data.getNgridPointsX() * data.getNgridPointsY() + (gridY+1) * data.getNgridPointsX() + gridX).setBy( campoSuPuntiGriglia.at(gridZ * data.getNgridPointsX() * data.getNgridPointsY() + (gridY+1) * data.getNgridPointsX() + gridX).getBy() );
        campoSuPuntiGriglia.at(gridZ * data.getNgridPointsX() * data.getNgridPointsY() + (gridY+1) * data.getNgridPointsX() + gridX).setBz( campoSuPuntiGriglia.at(gridZ * data.getNgridPointsX() * data.getNgridPointsY() + (gridY+1) * data.getNgridPointsX() + gridX).getBz() );

        campoSuPuntiGriglia.at(gridZ * data.getNgridPointsX() * data.getNgridPointsY() + (gridY+1) * data.getNgridPointsX() + (gridX+1)).setEx( campoSuPuntiGriglia.at(gridZ * data.getNgridPointsX() * data.getNgridPointsY() + (gridY+1) * data.getNgridPointsX() + (gridX+1)).getEx() );   // quarto punto
        campoSuPuntiGriglia.at(gridZ * data.getNgridPointsX() * data.getNgridPointsY() + (gridY+1) * data.getNgridPointsX() + (gridX+1)).setEy( campoSuPuntiGriglia.at(gridZ * data.getNgridPointsX() * data.getNgridPointsY() + (gridY+1) * data.getNgridPointsX() + (gridX+1)).getEy() );
        campoSuPuntiGriglia.at(gridZ * data.getNgridPointsX() * data.getNgridPointsY() + (gridY+1) * data.getNgridPointsX() + (gridX+1)).setEz( campoSuPuntiGriglia.at(gridZ * data.getNgridPointsX() * data.getNgridPointsY() + (gridY+1) * data.getNgridPointsX() + (gridX+1)).getEz() );
        campoSuPuntiGriglia.at(gridZ * data.getNgridPointsX() * data.getNgridPointsY() + (gridY+1) * data.getNgridPointsX() + (gridX+1)).setBx( campoSuPuntiGriglia.at(gridZ * data.getNgridPointsX() * data.getNgridPointsY() + (gridY+1) * data.getNgridPointsX() + (gridX+1)).getBx() );
        campoSuPuntiGriglia.at(gridZ * data.getNgridPointsX() * data.getNgridPointsY() + (gridY+1) * data.getNgridPointsX() + (gridX+1)).setBy( campoSuPuntiGriglia.at(gridZ * data.getNgridPointsX() * data.getNgridPointsY() + (gridY+1) * data.getNgridPointsX() + (gridX+1)).getBy() );
        campoSuPuntiGriglia.at(gridZ * data.getNgridPointsX() * data.getNgridPointsY() + (gridY+1) * data.getNgridPointsX() + (gridX+1)).setBz( campoSuPuntiGriglia.at(gridZ * data.getNgridPointsX() * data.getNgridPointsY() + (gridY+1) * data.getNgridPointsX() + (gridX+1)).getBz() );

        campoSuPuntiGriglia.at((gridZ+1) * data.getNgridPointsX() * data.getNgridPointsY() + gridY * data.getNgridPointsX() + gridX).setEx( campoSuPuntiGriglia.at((gridZ+1) * data.getNgridPointsX() * data.getNgridPointsY() + gridY * data.getNgridPointsX() + gridX).getEx() );   // quinto punto
        campoSuPuntiGriglia.at((gridZ+1) * data.getNgridPointsX() * data.getNgridPointsY() + gridY * data.getNgridPointsX() + gridX).setEy( campoSuPuntiGriglia.at((gridZ+1) * data.getNgridPointsX() * data.getNgridPointsY() + gridY * data.getNgridPointsX() + gridX).getEy() );
        campoSuPuntiGriglia.at((gridZ+1) * data.getNgridPointsX() * data.getNgridPointsY() + gridY * data.getNgridPointsX() + gridX).setEz( campoSuPuntiGriglia.at((gridZ+1) * data.getNgridPointsX() * data.getNgridPointsY() + gridY * data.getNgridPointsX() + gridX).getEz() );
        campoSuPuntiGriglia.at((gridZ+1) * data.getNgridPointsX() * data.getNgridPointsY() + gridY * data.getNgridPointsX() + gridX).setBx( campoSuPuntiGriglia.at((gridZ+1) * data.getNgridPointsX() * data.getNgridPointsY() + gridY * data.getNgridPointsX() + gridX).getBx() );
        campoSuPuntiGriglia.at((gridZ+1) * data.getNgridPointsX() * data.getNgridPointsY() + gridY * data.getNgridPointsX() + gridX).setBy( campoSuPuntiGriglia.at((gridZ+1) * data.getNgridPointsX() * data.getNgridPointsY() + gridY * data.getNgridPointsX() + gridX).getBy() );
        campoSuPuntiGriglia.at((gridZ+1) * data.getNgridPointsX() * data.getNgridPointsY() + gridY * data.getNgridPointsX() + gridX).setBz( campoSuPuntiGriglia.at((gridZ+1) * data.getNgridPointsX() * data.getNgridPointsY() + gridY * data.getNgridPointsX() + gridX).getBz() );

        campoSuPuntiGriglia.at((gridZ+1) * data.getNgridPointsX() * data.getNgridPointsY() + gridY * data.getNgridPointsX() + (gridX+1)).setEx( campoSuPuntiGriglia.at((gridZ+1) * data.getNgridPointsX() * data.getNgridPointsY() + gridY * data.getNgridPointsX() + (gridX+1)).getEx() );   // sesto punto
        campoSuPuntiGriglia.at((gridZ+1) * data.getNgridPointsX() * data.getNgridPointsY() + gridY * data.getNgridPointsX() + (gridX+1)).setEy( campoSuPuntiGriglia.at((gridZ+1) * data.getNgridPointsX() * data.getNgridPointsY() + gridY * data.getNgridPointsX() + (gridX+1)).getEy() );
        campoSuPuntiGriglia.at((gridZ+1) * data.getNgridPointsX() * data.getNgridPointsY() + gridY * data.getNgridPointsX() + (gridX+1)).setEz( campoSuPuntiGriglia.at((gridZ+1) * data.getNgridPointsX() * data.getNgridPointsY() + gridY * data.getNgridPointsX() + (gridX+1)).getEz() );
        campoSuPuntiGriglia.at((gridZ+1) * data.getNgridPointsX() * data.getNgridPointsY() + gridY * data.getNgridPointsX() + (gridX+1)).setBx( campoSuPuntiGriglia.at((gridZ+1) * data.getNgridPointsX() * data.getNgridPointsY() + gridY * data.getNgridPointsX() + (gridX+1)).getBx() );
        campoSuPuntiGriglia.at((gridZ+1) * data.getNgridPointsX() * data.getNgridPointsY() + gridY * data.getNgridPointsX() + (gridX+1)).setBy( campoSuPuntiGriglia.at((gridZ+1) * data.getNgridPointsX() * data.getNgridPointsY() + gridY * data.getNgridPointsX() + (gridX+1)).getBy() );
        campoSuPuntiGriglia.at((gridZ+1) * data.getNgridPointsX() * data.getNgridPointsY() + gridY * data.getNgridPointsX() + (gridX+1)).setBz( campoSuPuntiGriglia.at((gridZ+1) * data.getNgridPointsX() * data.getNgridPointsY() + gridY * data.getNgridPointsX() + (gridX+1)).getBz() );

        campoSuPuntiGriglia.at((gridZ+1) * data.getNgridPointsX() * data.getNgridPointsY() + (gridY+1) * data.getNgridPointsX() + gridX).setEx( campoSuPuntiGriglia.at((gridZ+1) * data.getNgridPointsX() * data.getNgridPointsY() + (gridY+1) * data.getNgridPointsX() + gridX).getEx() );   // settimo punto
        campoSuPuntiGriglia.at((gridZ+1) * data.getNgridPointsX() * data.getNgridPointsY() + (gridY+1) * data.getNgridPointsX() + gridX).setEy( campoSuPuntiGriglia.at((gridZ+1) * data.getNgridPointsX() * data.getNgridPointsY() + (gridY+1) * data.getNgridPointsX() + gridX).getEy() );
        campoSuPuntiGriglia.at((gridZ+1) * data.getNgridPointsX() * data.getNgridPointsY() + (gridY+1) * data.getNgridPointsX() + gridX).setEz( campoSuPuntiGriglia.at((gridZ+1) * data.getNgridPointsX() * data.getNgridPointsY() + (gridY+1) * data.getNgridPointsX() + gridX).getEz() );
        campoSuPuntiGriglia.at((gridZ+1) * data.getNgridPointsX() * data.getNgridPointsY() + (gridY+1) * data.getNgridPointsX() + gridX).setBx( campoSuPuntiGriglia.at((gridZ+1) * data.getNgridPointsX() * data.getNgridPointsY() + (gridY+1) * data.getNgridPointsX() + gridX).getBx() );
        campoSuPuntiGriglia.at((gridZ+1) * data.getNgridPointsX() * data.getNgridPointsY() + (gridY+1) * data.getNgridPointsX() + gridX).setBy( campoSuPuntiGriglia.at((gridZ+1) * data.getNgridPointsX() * data.getNgridPointsY() + (gridY+1) * data.getNgridPointsX() + gridX).getBy() );
        campoSuPuntiGriglia.at((gridZ+1) * data.getNgridPointsX() * data.getNgridPointsY() + (gridY+1) * data.getNgridPointsX() + gridX).setBz( campoSuPuntiGriglia.at((gridZ+1) * data.getNgridPointsX() * data.getNgridPointsY() + (gridY+1) * data.getNgridPointsX() + gridX).getBz() );

        campoSuPuntiGriglia.at((gridZ+1) * data.getNgridPointsX() * data.getNgridPointsY() + (gridY+1) * data.getNgridPointsX() + (gridX+1)).setEx( campoSuPuntiGriglia.at((gridZ+1) * data.getNgridPointsX() * data.getNgridPointsY() + (gridY+1) * data.getNgridPointsX() + (gridX+1)).getEx() );   // ottavo punto
        campoSuPuntiGriglia.at((gridZ+1) * data.getNgridPointsX() * data.getNgridPointsY() + (gridY+1) * data.getNgridPointsX() + (gridX+1)).setEy( campoSuPuntiGriglia.at((gridZ+1) * data.getNgridPointsX() * data.getNgridPointsY() + (gridY+1) * data.getNgridPointsX() + (gridX+1)).getEx() );
        campoSuPuntiGriglia.at((gridZ+1) * data.getNgridPointsX() * data.getNgridPointsY() + (gridY+1) * data.getNgridPointsX() + (gridX+1)).setEz( campoSuPuntiGriglia.at((gridZ+1) * data.getNgridPointsX() * data.getNgridPointsY() + (gridY+1) * data.getNgridPointsX() + (gridX+1)).getEx() );
        campoSuPuntiGriglia.at((gridZ+1) * data.getNgridPointsX() * data.getNgridPointsY() + (gridY+1) * data.getNgridPointsX() + (gridX+1)).setBx( campoSuPuntiGriglia.at((gridZ+1) * data.getNgridPointsX() * data.getNgridPointsY() + (gridY+1) * data.getNgridPointsX() + (gridX+1)).getEx() );
        campoSuPuntiGriglia.at((gridZ+1) * data.getNgridPointsX() * data.getNgridPointsY() + (gridY+1) * data.getNgridPointsX() + (gridX+1)).setBy( campoSuPuntiGriglia.at((gridZ+1) * data.getNgridPointsX() * data.getNgridPointsY() + (gridY+1) * data.getNgridPointsX() + (gridX+1)).getEx() );
        campoSuPuntiGriglia.at((gridZ+1) * data.getNgridPointsX() * data.getNgridPointsY() + (gridY+1) * data.getNgridPointsX() + (gridX+1)).setBz( campoSuPuntiGriglia.at((gridZ+1) * data.getNgridPointsX() * data.getNgridPointsY() + (gridY+1) * data.getNgridPointsX() + (gridX+1)).getEx() );
        */
      }
    }
  }
}


void calcolaCampoInUnIstante(Data data, vector<Particle> &particelle, vector<Field> &campiSuParticelle, ofstream &outputDATA)
{
  /*
  In questa funzione si genera un vettore contenente n particelle, n da definire prima di chiamare la funzione
  invocando data.setNelectrons() se si vuole avere qualcosa indietro, in posizione random. Su ciascuna viene calcolato
  il valore del campo analitico che deve essere gia` stato definito prima di chiamare la funzione invocando
  data.setK() e data.setA()  ( E= - k A sin[k(z-ct)] ) in un dato istante gia` definito tramite data.impostaT(0)
  */

  calcolaCampoAnalitico1DSuParticelle(data, particelle, campiSuParticelle);
  for (int i = 0; i < data.getNelectrons(); i++)
  {
    outputDATA
      << setw(12) << i + 1 << "\t"
      << setprecision(14) << setw(16) << particelle[i].getParticleX() << "\t"
      << setw(16) << particelle[i].getParticleY() << "\t"
      << setw(16) << particelle[i].getParticleZ() << "\t"
      << setw(16) << data.getT() << "\t"
      << setw(16) << campiSuParticelle[i].getEx() << "\t"
      << setw(16) << campiSuParticelle[i].getEy() << "\t"
      << setw(16) << campiSuParticelle[i].getEz() << "\t"
      << setw(16) << campiSuParticelle[i].getBx() << "\t"
      << setw(16) << campiSuParticelle[i].getBy() << "\t"
      << setw(16) << campiSuParticelle[i].getBz() << endl;
  }
}


void calcolaEvoluzioneCampoSuParticelle(Data data, vector<Particle> &particelle, vector<Field> &campiSuParticelle, ofstream &outputDATA)
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


  Field *campoTemp;
  campoTemp = new Field[1];
  campiSuParticelle.clear();

  for (int i = 0; i < data.getNsteps(); i++)
  {
    for (int j = 0; j < data.getNelectrons(); j++)
    {
      calcolaCampoAnalitico1DSuParticella(data, particelle[j], *campoTemp);
      campiSuParticelle.push_back(*campoTemp);
      outputDATA
        << setw(12) << j + 1 << "\t"
        << setprecision(14) << setw(16) << particelle[j].getParticleX() << "\t"
        << setw(16) << particelle[j].getParticleY() << "\t"
        << setw(16) << particelle[j].getParticleZ() << "\t"
        << setw(16) << data.getT() << "\t"
        << setw(16) << campiSuParticelle[j].getEx() << "\t"
        << setw(16) << campiSuParticelle[j].getEy() << "\t"
        << setw(16) << campiSuParticelle[j].getEz() << "\t"
        << setw(16) << campiSuParticelle[j].getBx() << "\t"
        << setw(16) << campiSuParticelle[j].getBy() << "\t"
        << setw(16) << campiSuParticelle[j].getBz() << endl;
    }
    data.aumentaT(data.getDeltaT());
    campiSuParticelle.clear();
  }
}


void evolveLPF4_nogrid_TURCHETTI(Data data, vector<Particle> &particles, vector<Field> &fieldOnParticles, ofstream &outputstream)
{
  Particle evolvedParticle;
  Particle pstar, pHalfAdvanced;

  Field b;
  Field *FOG;
  FOG = new Field[1];
  vector<Field> potenziale;
  Field *campoSuPunto;
  campoSuPunto = new Field[1];

  double gamma, gammastar, coeff;
  double I1, I2, I3;
  double zeta;

  double alpha = 1. / (2. - pow(2., (1. / 3.)));
  double dtx1 = alpha*data.getDeltaT()*0.5;
  double dtp1 = alpha*data.getDeltaT();
  double dtx2 = (1. - alpha)*0.5*data.getDeltaT();
  double dtp2 = (1. - 2.*alpha)*data.getDeltaT();
  double dtx3 = (1. - alpha)*0.5*data.getDeltaT();
  double dtp3 = alpha*data.getDeltaT();

  for (int i = 0; i < data.getNsteps(); i++)
  {
    for (int j = 0; j < data.getNelectrons(); j++)
    {
      fieldOnParticles.clear();

      // PASSO #1
      gamma = sqrt(1. + pow(particles[j].getParticlePX(), 2) + pow(particles[j].getParticlePY(), 2));
      particles[j].setParticleX(particles[j].getParticleX() + dtx1 * particles[j].getParticlePX() / (gamma));
      particles[j].setParticleY(particles[j].getParticleY() + dtx1 * particles[j].getParticlePY() / (gamma));

      calcolaCampoAnalitico1DSuParticella(data, particles[j], *FOG);
      fieldOnParticles.push_back(*FOG);

      pstar.setParticlePX(particles[j].getParticlePX());
      pstar.setParticlePY(particles[j].getParticlePY() + dtp1 * 0.5 * fieldOnParticles[0].getEy());

      gammastar = sqrt(1. + pow(pstar.getParticlePX(), 2) + pow(pstar.getParticlePY(), 2));
      b.setEy(0.5 * dtp1 * fieldOnParticles[0].getEy() / gammastar);
      coeff = 1. / (1. + pow(b.getEy(), 2));

      pHalfAdvanced.setParticlePX(coeff * (pstar.getParticlePX() + pstar.getParticlePY() * b.getBz()));
      pHalfAdvanced.setParticlePY(coeff * (pstar.getParticlePY() - pstar.getParticlePX() * b.getBz()));

      particles[j].setParticlePX(2.* pHalfAdvanced.getParticlePX() - particles[j].getParticlePX());
      particles[j].setParticlePY(2.* pHalfAdvanced.getParticlePY() - particles[j].getParticlePY());

      data.aumentaT(dtp1);



      // PASSO #2
      gamma = sqrt(1. + pow(particles[j].getParticlePX(), 2) + pow(particles[j].getParticlePY(), 2));
      particles[j].setParticleX(particles[j].getParticleX() + dtx2 * particles[j].getParticlePX() / (gamma));
      particles[j].setParticleY(particles[j].getParticleY() + dtx2 * particles[j].getParticlePY() / (gamma));

      calcolaCampoAnalitico1DSuParticella(data, particles[j], *FOG);
      fieldOnParticles.push_back(*FOG);

      pstar.setParticlePX(particles[j].getParticlePX());
      pstar.setParticlePY(particles[j].getParticlePY() + dtp2 * 0.5 * fieldOnParticles[1].getEy());

      gammastar = sqrt(1. + pow(pstar.getParticlePX(), 2) + pow(pstar.getParticlePY(), 2));
      b.setEy(0.5 * dtp2 * fieldOnParticles[1].getEy() / gammastar);
      coeff = 1. / (1. + pow(b.getEy(), 2));

      pHalfAdvanced.setParticlePX(coeff * (pstar.getParticlePX() + pstar.getParticlePY() * b.getEy()));   //?????
      pHalfAdvanced.setParticlePY(coeff * (pstar.getParticlePY() - pstar.getParticlePX() * b.getEy()));   //?????

      particles[j].setParticlePX(2.* pHalfAdvanced.getParticlePX() - particles[j].getParticlePX());
      particles[j].setParticlePY(2.* pHalfAdvanced.getParticlePY() - particles[j].getParticlePY());

      data.aumentaT(dtp2);


      // PASSO #3
      gamma = sqrt(1. + pow(particles[j].getParticlePX(), 2) + pow(particles[j].getParticlePY(), 2));
      particles[j].setParticleX(particles[j].getParticleX() + dtx3 * particles[j].getParticlePX() / (gamma));
      particles[j].setParticleY(particles[j].getParticleY() + dtx3 * particles[j].getParticlePY() / (gamma));

      calcolaCampoAnalitico1DSuParticella(data, particles[j], *FOG);
      fieldOnParticles.push_back(*FOG);

      pstar.setParticlePX(particles[j].getParticlePX());
      pstar.setParticlePY(particles[j].getParticlePY() + dtp3 * 0.5 * fieldOnParticles[2].getEy());

      gammastar = sqrt(1. + pow(pstar.getParticlePX(), 2) + pow(pstar.getParticlePY(), 2));
      b.setEy(0.5 * dtp3 * fieldOnParticles[2].getEy() / gammastar);
      coeff = 1. / (1. + pow(b.getEy(), 2));

      pHalfAdvanced.setParticlePX(coeff * (pstar.getParticlePX() + pstar.getParticlePY() * b.getEy()));   //?????
      pHalfAdvanced.setParticlePY(coeff * (pstar.getParticlePY() - pstar.getParticlePX() * b.getEy()));   //?????

      particles[j].setParticlePX(2.* pHalfAdvanced.getParticlePX() - particles[j].getParticlePX());
      particles[j].setParticlePY(2.* pHalfAdvanced.getParticlePY() - particles[j].getParticlePY());


      // PASSO #4 - SOLO EVOLUZIONE POSIZIONI
      gamma = sqrt(1. + pow(particles[j].getParticlePX(), 2) + pow(particles[j].getParticlePY(), 2));

      particles[j].setParticleX(particles[j].getParticleX() + dtx1 * particles[j].getParticlePX() / gamma);
      particles[j].setParticleY(particles[j].getParticleY() + dtx1 * particles[j].getParticlePY() / gamma);

      data.aumentaT(dtp3);



      particles[j].setParticleCell(data);


      //INTEGRALI PRIMI (da sistemare, li ho copiati dall'evoluzione sui campi analitici)
      potenziale.clear();
      zeta = particles[j].getParticleX() - C * data.getT();             // non usato nella mia struttura delle equazioni di campo
      gamma = sqrt(1 + pow(particles[j].getParticlePX(), 2) + pow(particles[j].getParticlePY(), 2) + pow(particles[j].getParticlePZ(), 2));
      calcolaPotenzialeAnalitico1DSuParticella(data, particles[j], *campoSuPunto);
      potenziale.push_back(*campoSuPunto);

      I1 = (particles[j].getParticlePY() / (particles[j].getParticleM()*C)) + potenziale[0].getEy();
      I2 = particles[j].getParticlePZ() / (particles[j].getParticleM()*C);
      I3 = gamma - particles[j].getParticlePX() / (particles[j].getParticleM()*C);

      potenziale.clear();


      outputstream
        << setw(12) << i + 1 << "\t"
        << setprecision(14) << setw(16) << particles[j].getParticleX() << "\t"
        << setw(16) << particles[j].getParticleY() << "\t"
        << setw(16) << particles[j].getParticleZ() << "\t"
        //        << setw(16) << particles[j].getParticleCell() << "\t"
        << setw(16) << particles[j].getParticlePX() << "\t"
        << setw(16) << particles[j].getParticlePY() << "\t"
        << setw(16) << particles[j].getParticlePZ() << "\t"
        << setw(16) << fieldOnParticles[0].getEy() << "\t"
        << setw(16) << fieldOnParticles[1].getEy() << "\t"
        << setw(16) << fieldOnParticles[2].getEy() << "\t"
        << setw(16) << zeta << "\t"
        << setw(16) << I1 << "\t"
        << setw(16) << I2 << "\t"
        << setw(16) << I3 << endl;
    }
  }
}


void evolveLPF4_nogrid(Data data, vector<Particle> &particles, ofstream &outputstream)
{
  Field *campoSuPunto;
  campoSuPunto = new Field[1];

  double gamma, gamma_inv;
  double I1, I2, I3;
  double zeta;

  double alpha = 1. / (2. - pow(2., (1. / 3.)));
  double dtx1 = alpha * data.getDeltaT() * 0.5;
  double dtp1 = alpha * data.getDeltaT();
  double dtx2 = (1. - alpha) * 0.5 * data.getDeltaT();
  double dtp2 = (1. - 2. * alpha) * data.getDeltaT();
  double dtx3 = (1. - alpha) * 0.5 * data.getDeltaT();
  double dtp3 = alpha * data.getDeltaT();

  //Il leapfrog di ordine 4 si ottiene componendo 3 schemi leapfrog invertiti e usando un peso alpha
  //   r(n+1/2) = r(n-1/2) + DeltaT(pesatoXciascunPassaggio) * p(n)/gamma(n)
  //   p(n+1) = p(n) + DeltaT(pesatoXciascunPassaggio) * F(n+1/2)
  //dove          F(n+1/2) = E(n+1/2) + p(n+1/2)/gamma(n+1/2) * B(n+1/2)
  //nel quale     p(n+1/2) = ( p(n) + p(n+1) ) / 2


  for (int i = 0; i < data.getNsteps(); i++)
  {
    for (int j = 0; j < data.getNelectrons(); j++)
    {
      // PASSO #1
      evolveLPF(data, particles[j], dtp1, dtx1);
      // PASSO #2
      evolveLPF(data, particles[j], dtp2, dtx2);
      // PASSO #3
      evolveLPF(data, particles[j], dtp3, dtx3);

      // PASSO #4 - evoluzione solo posizioni
      gamma = sqrt(1. + pow(particles[j].getParticlePX(), 2) + pow(particles[j].getParticlePY(), 2));
      gamma_inv = 1. / gamma;
      particles[j].setParticleX(particles[j].getParticleX() + dtx1 * particles[j].getParticlePX() * gamma_inv);
      particles[j].setParticleY(particles[j].getParticleY() + dtx1 * particles[j].getParticlePY() * gamma_inv);
      particles[j].setParticleZ(particles[j].getParticleZ() + dtx1 * particles[j].getParticlePZ() * gamma_inv);
      particles[j].setParticleCell(data);

      data.aumentaT(dtx1);

      zeta = particles[j].getParticleX() - C * data.getT();             // non usato nella mia struttura delle equazioni di campo
      calcolaPotenzialeAnalitico1DSuParticella(data, particles[j], *campoSuPunto);

      I1 = (particles[j].getParticlePY() / (particles[j].getParticleM()*C)) + (*campoSuPunto).getEy();
      I2 = particles[j].getParticlePZ() / (particles[j].getParticleM()*C);
      I3 = gamma - particles[j].getParticlePX() / (particles[j].getParticleM()*C);


      outputstream
        << setw(12) << i + 1 << "\t"
        << setprecision(14)
        << setw(16) << particles[j].getParticleX() << "\t"
        << setw(16) << particles[j].getParticleY() << "\t"
        << setw(16) << particles[j].getParticleZ() << "\t"
        << setw(16) << particles[j].getParticlePX() << "\t"
        << setw(16) << particles[j].getParticlePY() << "\t"
        << setw(16) << particles[j].getParticlePZ() << "\t"
        << setw(16) << zeta << "\t"
        << setw(16) << I1 << "\t"
        << setw(16) << I2 << "\t"
        << setw(16) << I3 << endl;

      data.aumentaT(-data.getDeltaT());   // cosi` un altro elettrone dello stesso step non si trova i tempi avanzati
    }
    data.aumentaT(data.getDeltaT());      // e` questo il vero avanzamento di uno step temporale del tempo del sistema
  }
}



void evolveLPF4_withgrid_onthefly(Data data, vector<Particle> &particles, vector<Field> &fieldOnParticles, ofstream &outputstream)
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

  Field *FOG;                             // DA SISTEMARE, c'è solo per avere un campo fittizio su una griglia on-the-fly per implementare l'algoritmo
  FOG = new Field[1];

  double gridXX, gridYY, gridZZ;
  int gridX, gridY, gridZ;
  double gamma, gamma_inv, b2;
  double I1, I2, I3;
  double zeta;

  double alpha = 1. / (2. - pow(2., (1. / 3.)));
  double dtx1 = 0.5*alpha*data.getDeltaT();
  double dtp1 = alpha*data.getDeltaT();
  double dtx2 = (1. - alpha)*0.5*data.getDeltaT();
  double dtp2 = (1. - 2 * alpha)*data.getDeltaT();
  double dtx3 = (1. - alpha)*0.5*data.getDeltaT();
  double dtp3 = alpha*data.getDeltaT();

  Field *campoSuPunto;
  campoSuPunto = new Field[1];

  //ALGORITMO LEAPFROG PARTICELLE

  for (int i = 0; i < data.getNsteps(); i++)
  {
    for (int j = 0; j < data.getNelectrons(); j++)
    {
      gamma = sqrt(1 + pow(particles[j].getParticlePX(), 2) + pow(particles[j].getParticlePY(), 2) + pow(particles[j].getParticlePZ(), 2));
      gamma_inv = 1. / gamma;

      particles[j].setParticleX(particles[j].getParticleX() + dtx1 * particles[j].getParticlePX() * gamma_inv);
      particles[j].setParticleY(particles[j].getParticleY() + dtx1 * particles[j].getParticlePY() * gamma_inv);
      particles[j].setParticleZ(particles[j].getParticleZ() + dtx1 * particles[j].getParticlePZ() * gamma_inv);


      fieldOnParticles.clear();

      // STEP #1

      gridXX = particles[j].getParticleX() / data.getDeltaX();
      if (gridXX < 0) gridX = ((int)gridXX) - 1;           //il punto di griglia lungo l'asse x a sinistra della particella
      else gridX = (int)gridXX;                      //

      gridYY = particles[j].getParticleY() / data.getDeltaY();
      if (gridYY < 0) gridY = ((int)gridYY) - 1;           //il punto di griglia lungo l'asse y a sinistra della particella
      else gridY = (int)gridYY;                      //

      gridZZ = particles[j].getParticleZ() / data.getDeltaZ();
      if (gridZZ < 0) gridZ = ((int)gridZZ) - 1;           //il punto di griglia lungo l'asse z a sinistra della particella
      else gridZ = (int)gridZZ;                      //

      coordinatePrimoPuntoGriglia.setParticleX(gridX*data.getDeltaX());
      coordinatePrimoPuntoGriglia.setParticleY(gridY*data.getDeltaY());
      coordinatePrimoPuntoGriglia.setParticleZ(gridZ*data.getDeltaZ());
      coordinateSecondoPuntoGriglia.setParticleX(gridX*data.getDeltaX() + data.getDeltaX());
      coordinateSecondoPuntoGriglia.setParticleY(gridY*data.getDeltaY());
      coordinateSecondoPuntoGriglia.setParticleZ(gridZ*data.getDeltaZ());
      coordinateTerzoPuntoGriglia.setParticleX(gridX*data.getDeltaX());
      coordinateTerzoPuntoGriglia.setParticleY(gridY*data.getDeltaY() + data.getDeltaY());
      coordinateTerzoPuntoGriglia.setParticleZ(gridZ*data.getDeltaZ());
      coordinateQuartoPuntoGriglia.setParticleX(gridX*data.getDeltaX() + data.getDeltaX());
      coordinateQuartoPuntoGriglia.setParticleY(gridY*data.getDeltaY() + data.getDeltaY());
      coordinateQuartoPuntoGriglia.setParticleZ(gridZ*data.getDeltaZ());
      /*
      coordinateQuintoPuntoGriglia.setParticleX(gridX*data.getDeltaX());
      coordinateQuintoPuntoGriglia.setParticleY(gridY*data.getDeltaY());
      coordinateQuintoPuntoGriglia.setParticleZ(gridZ*data.getDeltaZ()+data.getDeltaZ());
      coordinateSestoPuntoGriglia.setParticleX(gridX*data.getDeltaX()+data.getDeltaX());
      coordinateSestoPuntoGriglia.setParticleY(gridY*data.getDeltaY());
      coordinateSestoPuntoGriglia.setParticleZ(gridZ*data.getDeltaZ()+data.getDeltaZ());
      coordinateSettimoPuntoGriglia.setParticleX(gridX*data.getDeltaX());
      coordinateSettimoPuntoGriglia.setParticleY(gridY*data.getDeltaY()+data.getDeltaY());
      coordinateSettimoPuntoGriglia.setParticleZ(gridZ*data.getDeltaZ()+data.getDeltaZ());
      coordinateOttavoPuntoGriglia.setParticleX(gridX*data.getDeltaX()+data.getDeltaX());
      coordinateOttavoPuntoGriglia.setParticleY(gridY*data.getDeltaY()+data.getDeltaY());
      coordinateOttavoPuntoGriglia.setParticleZ(gridZ*data.getDeltaZ()+data.getDeltaZ());
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
      pstar.setParticlePX(particles[j].getParticlePX() + 0.5*dtp1 * particles[j].getParticleQ() * fieldOnParticles[0].getEx());
      pstar.setParticlePY(particles[j].getParticlePY() + 0.5*dtp1 * particles[j].getParticleQ() * fieldOnParticles[0].getEy());
      pstar.setParticlePZ(particles[j].getParticlePZ() + 0.5*dtp1 * particles[j].getParticleQ() * fieldOnParticles[0].getEz());
      gamma = sqrt(1 + pow(pstar.getParticlePX(), 2) + pow(pstar.getParticlePY(), 2) + pow(pstar.getParticlePZ(), 2));
      gamma_inv = 1. / gamma;

      //b = (DeltaT/2) * Q * B(n+1/2) / gamma(aggiornato al pstar)
      b.setBx(0.5*dtp1 * fieldOnParticles[0].getBx() * particles[j].getParticleQ() * gamma_inv);
      b.setBy(0.5*dtp1 * fieldOnParticles[0].getBy() * particles[j].getParticleQ() * gamma_inv);
      b.setBz(0.5*dtp1 * fieldOnParticles[0].getBz() * particles[j].getParticleQ() * gamma_inv);
      b2 = 1. + pow(b.getBx(), 2) + pow(b.getBy(), 2) + pow(b.getBz(), 2);

      /****************************************************************************
      //p(n+1/2) = [pstar + pstar x b + (pstar x b) * b] / b2   SECONDO TURCHETTI  (ma i conti analitici diconon che e` sbagliato, anche perche' l'ultimo termine e` uno scalare!)
      pHalfAdvanced.setParticlePX( (pstar.getParticlePX() + pstar.getParticlePY() * b.getBz() - b.getBy() * pstar.getParticlePZ() + (pstar.getParticlePY() * b.getBz() - b.getBy() * pstar.getParticlePZ()) * b.getBx() ) / b2 );
      pHalfAdvanced.setParticlePY( (pstar.getParticlePY() + pstar.getParticlePZ() * b.getBx() - b.getBz() * pstar.getParticlePX() + (pstar.getParticlePZ() * b.getBx() - b.getBz() * pstar.getParticlePX()) * b.getBy() ) / b2 );
      pHalfAdvanced.setParticlePZ( (pstar.getParticlePZ() + pstar.getParticlePX() * b.getBy() - b.getBx() * pstar.getParticlePY() + (pstar.getParticlePX() * b.getBy() - b.getBx() * pstar.getParticlePY()) * b.getBz() ) / b2 );
      ****************************************************************************/

      /****************************************************************************
      //p(n+1/2) = [pstar + pstar x b + b x (pstar x b)] / b2   SECONDO LONDRILLO
      pHalfAdvanced.setParticlePX( (pstar.getParticlePX() + pstar.getParticlePY() * b.getBz() - b.getBy() * pstar.getParticlePZ() + b.getBx() * pstar.getParticlePX() * b.getBx() ) / b2 );
      pHalfAdvanced.setParticlePY( (pstar.getParticlePY() + pstar.getParticlePZ() * b.getBx() - b.getBz() * pstar.getParticlePX() + b.getBy() * pstar.getParticlePY() * b.getBy() ) / b2 );
      pHalfAdvanced.setParticlePZ( (pstar.getParticlePZ() + pstar.getParticlePX() * b.getBy() - b.getBx() * pstar.getParticlePY() + b.getBz() * pstar.getParticlePZ() * b.getBz() ) / b2 );
      ****************************************************************************/

      // /****************************************************************************
            //p(n+1/2) = [pstar + pstar x b] / b2   COME IMPLEMENTATO DA TURCHETTI NEL CODICE
      pHalfAdvanced.setParticlePX((pstar.getParticlePX() + pstar.getParticlePY() * b.getBz()) / b2);
      pHalfAdvanced.setParticlePY((pstar.getParticlePY() - pstar.getParticlePX() * b.getBz()) / b2);
      pHalfAdvanced.setParticlePZ(0.);  //per come ho strutturato i pacchetti d'onda su z non dovrebbe succedere nulla e quindi lo fisso io a zero.
      // ****************************************************************************/

      particles[j].setParticlePX(2 * pHalfAdvanced.getParticlePX() - particles[j].getParticlePX());
      particles[j].setParticlePY(2 * pHalfAdvanced.getParticlePY() - particles[j].getParticlePY());
      particles[j].setParticlePZ(2 * pHalfAdvanced.getParticlePZ() - particles[j].getParticlePZ());


      // STEP #2
      gamma = sqrt(1 + pow(particles[j].getParticlePX(), 2) + pow(particles[j].getParticlePY(), 2) + pow(particles[j].getParticlePZ(), 2));
      gamma_inv = 1. / gamma;

      particles[j].setParticleX(particles[j].getParticleX() + dtx2 * particles[j].getParticlePX() * gamma_inv);
      particles[j].setParticleY(particles[j].getParticleY() + dtx2 * particles[j].getParticlePY() * gamma_inv);
      particles[j].setParticleZ(particles[j].getParticleZ() + dtx2 * particles[j].getParticlePZ() * gamma_inv);


      gridXX = particles[j].getParticleX() / data.getDeltaX();
      if (gridXX < 0) gridX = ((int)gridXX) - 1;           //il punto di griglia lungo l'asse x a sinistra della particella
      else gridX = (int)gridXX;

      gridYY = particles[j].getParticleY() / data.getDeltaY();
      if (gridYY < 0) gridY = ((int)gridYY) - 1;           //il punto di griglia lungo l'asse y a sinistra della particella
      else gridY = (int)gridYY;

      gridZZ = particles[j].getParticleZ() / data.getDeltaZ();
      if (gridZZ < 0) gridZ = ((int)gridZZ) - 1;           //il punto di griglia lungo l'asse z a sinistra della particella
      else gridZ = (int)gridZZ;

      coordinatePrimoPuntoGriglia.setParticleX(gridX*data.getDeltaX());
      coordinatePrimoPuntoGriglia.setParticleY(gridY*data.getDeltaY());
      coordinatePrimoPuntoGriglia.setParticleZ(gridZ*data.getDeltaZ());
      coordinateSecondoPuntoGriglia.setParticleX(gridX*data.getDeltaX() + data.getDeltaX());
      coordinateSecondoPuntoGriglia.setParticleY(gridY*data.getDeltaY());
      coordinateSecondoPuntoGriglia.setParticleZ(gridZ*data.getDeltaZ());
      coordinateTerzoPuntoGriglia.setParticleX(gridX*data.getDeltaX());
      coordinateTerzoPuntoGriglia.setParticleY(gridY*data.getDeltaY() + data.getDeltaY());
      coordinateTerzoPuntoGriglia.setParticleZ(gridZ*data.getDeltaZ());
      coordinateQuartoPuntoGriglia.setParticleX(gridX*data.getDeltaX() + data.getDeltaX());
      coordinateQuartoPuntoGriglia.setParticleY(gridY*data.getDeltaY() + data.getDeltaY());
      coordinateQuartoPuntoGriglia.setParticleZ(gridZ*data.getDeltaZ());

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
      pstar.setParticlePX(particles[j].getParticlePX() + 0.5*dtp2 * particles[j].getParticleQ() * fieldOnParticles[1].getEx());
      pstar.setParticlePY(particles[j].getParticlePY() + 0.5*dtp2 * particles[j].getParticleQ() * fieldOnParticles[1].getEy());
      pstar.setParticlePZ(particles[j].getParticlePZ() + 0.5*dtp2 * particles[j].getParticleQ() * fieldOnParticles[1].getEz());
      gamma = sqrt(1 + pow(pstar.getParticlePX(), 2) + pow(pstar.getParticlePY(), 2) + pow(pstar.getParticlePZ(), 2));
      gamma_inv = 1. / gamma;

      //b = (DeltaT/2) * Q * B(n+1/2) / gamma(aggiornato al pstar)
      b.setBx(0.5*dtp2 * fieldOnParticles[1].getBx() * particles[j].getParticleQ() * gamma_inv);
      b.setBy(0.5*dtp2 * fieldOnParticles[1].getBy() * particles[j].getParticleQ() * gamma_inv);
      b.setBz(0.5*dtp2 * fieldOnParticles[1].getBz() * particles[j].getParticleQ() * gamma_inv);
      b2 = 1. + pow(b.getBx(), 2) + pow(b.getBy(), 2) + pow(b.getBz(), 2);

      /****************************************************************************
      //p(n+1/2) = [pstar + pstar x b + (pstar x b) * b] / b2   SECONDO TURCHETTI  (ma i conti analitici diconon che e` sbagliato, anche perche' l'ultimo termine e` uno scalare!)
      pHalfAdvanced.setParticlePX( (pstar.getParticlePX() + pstar.getParticlePY() * b.getBz() - b.getBy() * pstar.getParticlePZ() + (pstar.getParticlePY() * b.getBz() - b.getBy() * pstar.getParticlePZ()) * b.getBx() ) / b2 );
      pHalfAdvanced.setParticlePY( (pstar.getParticlePY() + pstar.getParticlePZ() * b.getBx() - b.getBz() * pstar.getParticlePX() + (pstar.getParticlePZ() * b.getBx() - b.getBz() * pstar.getParticlePX()) * b.getBy() ) / b2 );
      pHalfAdvanced.setParticlePZ( (pstar.getParticlePZ() + pstar.getParticlePX() * b.getBy() - b.getBx() * pstar.getParticlePY() + (pstar.getParticlePX() * b.getBy() - b.getBx() * pstar.getParticlePY()) * b.getBz() ) / b2 );
      ****************************************************************************/

      /****************************************************************************
      //p(n+1/2) = [pstar + pstar x b + b x (pstar x b)] / b2   SECONDO LONDRILLO
      pHalfAdvanced.setParticlePX( (pstar.getParticlePX() + pstar.getParticlePY() * b.getBz() - b.getBy() * pstar.getParticlePZ() + b.getBx() * pstar.getParticlePX() * b.getBx() ) / b2 );
      pHalfAdvanced.setParticlePY( (pstar.getParticlePY() + pstar.getParticlePZ() * b.getBx() - b.getBz() * pstar.getParticlePX() + b.getBy() * pstar.getParticlePY() * b.getBy() ) / b2 );
      pHalfAdvanced.setParticlePZ( (pstar.getParticlePZ() + pstar.getParticlePX() * b.getBy() - b.getBx() * pstar.getParticlePY() + b.getBz() * pstar.getParticlePZ() * b.getBz() ) / b2 );
      ****************************************************************************/

      // /****************************************************************************
      //p(n+1/2) = [pstar + pstar x b] / b2   COME IMPLEMENTATO DA TURCHETTI NEL CODICE
      pHalfAdvanced.setParticlePX((pstar.getParticlePX() + pstar.getParticlePY() * b.getBz()) / b2);
      pHalfAdvanced.setParticlePY((pstar.getParticlePY() - pstar.getParticlePX() * b.getBz()) / b2);
      pHalfAdvanced.setParticlePZ(0.);  //per come ho strutturato i pacchetti d'onda su z non dovrebbe succedere nulla e quindi lo fisso io a zero.
      // ****************************************************************************/


      particles[j].setParticlePX(2 * pHalfAdvanced.getParticlePX() - particles[j].getParticlePX());
      particles[j].setParticlePY(2 * pHalfAdvanced.getParticlePY() - particles[j].getParticlePY());
      particles[j].setParticlePZ(2 * pHalfAdvanced.getParticlePZ() - particles[j].getParticlePZ());


      // STEP #3
      gamma = sqrt(1 + pow(particles[j].getParticlePX(), 2) + pow(particles[j].getParticlePY(), 2) + pow(particles[j].getParticlePZ(), 2));
      gamma_inv = 1. / gamma;

      particles[j].setParticleX(particles[j].getParticleX() + dtx3 * particles[j].getParticlePX() * gamma_inv);
      particles[j].setParticleY(particles[j].getParticleY() + dtx3 * particles[j].getParticlePY() * gamma_inv);
      particles[j].setParticleZ(particles[j].getParticleZ() + dtx3 * particles[j].getParticlePZ() * gamma_inv);

      gridXX = particles[j].getParticleX() / data.getDeltaX();
      if (gridXX < 0) gridX = ((int)gridXX) - 1;           //il punto di griglia lungo l'asse x a sinistra della particella
      else gridX = (int)gridXX;

      gridYY = particles[j].getParticleY() / data.getDeltaY();
      if (gridYY < 0) gridY = ((int)gridYY) - 1;           //il punto di griglia lungo l'asse y a sinistra della particella
      else gridY = (int)gridYY;

      gridZZ = particles[j].getParticleZ() / data.getDeltaZ();
      if (gridZZ < 0) gridZ = ((int)gridZZ) - 1;           //il punto di griglia lungo l'asse z a sinistra della particella
      else gridZ = (int)gridZZ;

      coordinatePrimoPuntoGriglia.setParticleX(gridX*data.getDeltaX());
      coordinatePrimoPuntoGriglia.setParticleY(gridY*data.getDeltaY());
      coordinatePrimoPuntoGriglia.setParticleZ(gridZ*data.getDeltaZ());
      coordinateSecondoPuntoGriglia.setParticleX(gridX*data.getDeltaX() + data.getDeltaX());
      coordinateSecondoPuntoGriglia.setParticleY(gridY*data.getDeltaY());
      coordinateSecondoPuntoGriglia.setParticleZ(gridZ*data.getDeltaZ());
      coordinateTerzoPuntoGriglia.setParticleX(gridX*data.getDeltaX());
      coordinateTerzoPuntoGriglia.setParticleY(gridY*data.getDeltaY() + data.getDeltaY());
      coordinateTerzoPuntoGriglia.setParticleZ(gridZ*data.getDeltaZ());
      coordinateQuartoPuntoGriglia.setParticleX(gridX*data.getDeltaX() + data.getDeltaX());
      coordinateQuartoPuntoGriglia.setParticleY(gridY*data.getDeltaY() + data.getDeltaY());
      coordinateQuartoPuntoGriglia.setParticleZ(gridZ*data.getDeltaZ());

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
      pstar.setParticlePX(particles[j].getParticlePX() + 0.5*dtp3 * particles[j].getParticleQ() * fieldOnParticles[2].getEx());
      pstar.setParticlePY(particles[j].getParticlePY() + 0.5*dtp3 * particles[j].getParticleQ() * fieldOnParticles[2].getEy());
      pstar.setParticlePZ(particles[j].getParticlePZ() + 0.5*dtp3 * particles[j].getParticleQ() * fieldOnParticles[2].getEz());
      gamma = sqrt(1 + pow(pstar.getParticlePX(), 2) + pow(pstar.getParticlePY(), 2) + pow(pstar.getParticlePZ(), 2));
      gamma_inv = 1. / gamma;

      //b = (DeltaT/2) * Q * B(n+1/2) / gamma(aggiornato al pstar)
      b.setBx(0.5*dtp3 * fieldOnParticles[2].getBx() * particles[j].getParticleQ() * gamma_inv);
      b.setBy(0.5*dtp3 * fieldOnParticles[2].getBy() * particles[j].getParticleQ() * gamma_inv);
      b.setBz(0.5*dtp3 * fieldOnParticles[2].getBz() * particles[j].getParticleQ() * gamma_inv);
      b2 = 1. + pow(b.getBx(), 2) + pow(b.getBy(), 2) + pow(b.getBz(), 2);


      /****************************************************************************
      //p(n+1/2) = [pstar + pstar x b + (pstar x b) * b] / b2   SECONDO TURCHETTI  (ma i conti analitici diconon che e` sbagliato, anche perche' l'ultimo termine e` uno scalare!)
      pHalfAdvanced.setParticlePX( (pstar.getParticlePX() + pstar.getParticlePY() * b.getBz() - b.getBy() * pstar.getParticlePZ() + (pstar.getParticlePY() * b.getBz() - b.getBy() * pstar.getParticlePZ()) * b.getBx() ) / b2 );
      pHalfAdvanced.setParticlePY( (pstar.getParticlePY() + pstar.getParticlePZ() * b.getBx() - b.getBz() * pstar.getParticlePX() + (pstar.getParticlePZ() * b.getBx() - b.getBz() * pstar.getParticlePX()) * b.getBy() ) / b2 );
      pHalfAdvanced.setParticlePZ( (pstar.getParticlePZ() + pstar.getParticlePX() * b.getBy() - b.getBx() * pstar.getParticlePY() + (pstar.getParticlePX() * b.getBy() - b.getBx() * pstar.getParticlePY()) * b.getBz() ) / b2 );
      ****************************************************************************/

      /****************************************************************************
      //p(n+1/2) = [pstar + pstar x b + b x (pstar x b)] / b2   SECONDO LONDRILLO
      pHalfAdvanced.setParticlePX( (pstar.getParticlePX() + pstar.getParticlePY() * b.getBz() - b.getBy() * pstar.getParticlePZ() + b.getBx() * pstar.getParticlePX() * b.getBx() ) / b2 );
      pHalfAdvanced.setParticlePY( (pstar.getParticlePY() + pstar.getParticlePZ() * b.getBx() - b.getBz() * pstar.getParticlePX() + b.getBy() * pstar.getParticlePY() * b.getBy() ) / b2 );
      pHalfAdvanced.setParticlePZ( (pstar.getParticlePZ() + pstar.getParticlePX() * b.getBy() - b.getBx() * pstar.getParticlePY() + b.getBz() * pstar.getParticlePZ() * b.getBz() ) / b2 );
      ****************************************************************************/

      // /****************************************************************************
      //p(n+1/2) = [pstar + pstar x b] / b2   COME IMPLEMENTATO DA TURCHETTI NEL CODICE
      pHalfAdvanced.setParticlePX((pstar.getParticlePX() + pstar.getParticlePY() * b.getBz()) / b2);
      pHalfAdvanced.setParticlePY((pstar.getParticlePY() - pstar.getParticlePX() * b.getBz()) / b2);
      pHalfAdvanced.setParticlePZ(0.);  //per come ho strutturato i pacchetti d'onda su z non dovrebbe succedere nulla e quindi lo fisso io a zero.
      // ****************************************************************************/

      particles[j].setParticlePX(2 * pHalfAdvanced.getParticlePX() - particles[j].getParticlePX());
      particles[j].setParticlePY(2 * pHalfAdvanced.getParticlePY() - particles[j].getParticlePY());
      particles[j].setParticlePZ(2 * pHalfAdvanced.getParticlePZ() - particles[j].getParticlePZ());



      // PASSO #4 - evoluzione solo posizioni
      gamma = sqrt(1. + pow(particles[j].getParticlePX(), 2) + pow(particles[j].getParticlePY(), 2));
      gamma_inv = 1. / gamma;

      particles[j].setParticleX(particles[j].getParticleX() + dtx1 * particles[j].getParticlePX() * gamma_inv);
      particles[j].setParticleY(particles[j].getParticleY() + dtx1 * particles[j].getParticlePY() * gamma_inv);
      particles[j].setParticleZ(particles[j].getParticleZ() + dtx1 * particles[j].getParticlePZ() * gamma_inv);
      particles[j].setParticleCell(data);

      data.aumentaT(dtx1);

      //INTEGRALI PRIMI (da sistemare, li ho copiati dall'evoluzione sui campi analitici)
      potenziale.clear();
      zeta = particles[j].getParticleX() - C * data.getT();             // non usato nella mia struttura delle equazioni di campo
      gamma = sqrt(1 + pow(particles[j].getParticlePX(), 2) + pow(particles[j].getParticlePY(), 2) + pow(particles[j].getParticlePZ(), 2));
      calcolaPotenzialeAnalitico1DSuParticella(data, particles[j], *campoSuPunto);
      potenziale.push_back(*campoSuPunto);

      I1 = (particles[j].getParticlePY() / (particles[j].getParticleM()*C)) + potenziale[0].getEy();
      I2 = particles[j].getParticlePZ() / (particles[j].getParticleM()*C);
      I3 = gamma - particles[j].getParticlePX() / (particles[j].getParticleM()*C);

      potenziale.clear();


      outputstream
        << setw(12) << i + 1 << "\t"
        << setprecision(14) << setw(16) << particles[j].getParticleX() << "\t"
        << setw(16) << particles[j].getParticleY() << "\t"
        << setw(16) << particles[j].getParticleZ() << "\t"
        // << setw(16) << particles[j].getParticleCell() << "\t"
        << setw(16) << particles[j].getParticlePX() << "\t"
        << setw(16) << particles[j].getParticlePY() << "\t"
        << setw(16) << particles[j].getParticlePZ() << "\t"
        << setw(16) << zeta << "\t"
        << setw(16) << I1 << "\t"
        << setw(16) << I2 << "\t"
        << setw(16) << I3 << endl;
    }
  }
}


void evolveLPF4_withgrid(Data data, vector<Particle> &particles, vector<Field> &fieldOnParticles, vector<Field> &campoSuPuntiGriglia, ofstream &outputstream)
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

  Field *campoSuPunto;
  campoSuPunto = new Field[1];

  //ALGORITMO LEAPFROG PARTICELLE

  if (data.getNdim() == 1)
  {
    for (int i = 0; i < data.getNsteps(); i++)
    {
      data.aumentaT(data.getDeltaT());

      campoSuPuntiGriglia.clear();
      riempiPuntiGrigliaConCampoAnalitico1D(data, campoSuPuntiGriglia);

      potenzialeSuPuntiGriglia.clear();
      riempiPuntiGrigliaConPotenzialeAnalitico1D(data, potenzialeSuPuntiGriglia);

      for (int j = 0; j < data.getNelectrons(); j++)
      {
        gamma = sqrt(1. + pow(particles[j].getParticlePX(), 2) + pow(particles[j].getParticlePY(), 2));
        gamma_inv = 1. / gamma;

        particles[j].setParticleX(particles[j].getParticleX() + 0.5 * data.getDeltaT() * particles[j].getParticlePX() * gamma_inv);
        particles[j].setParticleY(particles[j].getParticleY() + 0.5 * data.getDeltaT() * particles[j].getParticlePY() * gamma_inv);
        particles[j].setParticleZ(particles[j].getParticleZ() + 0.5 * data.getDeltaT() * particles[j].getParticlePZ() * gamma_inv);


        //QUI MI CALCOLO I PUNTI DI RIFERIMENTO PER LA PARTICELLA IN ESAME nella posizione attuale
        //forse c'e` un modo piu` furbo per farlo sfruttando il # di cella in cui si trova che e` noto

        fieldOnParticles.clear();       //SVUOTA il vettore dei campi sulle particelle ad ogni ciclo


        gridXX = particles[j].getParticleX() / data.getDeltaX();
        if (gridXX < 0) gridX = ((int)gridXX) - 1;           //il punto di griglia lungo l'asse x a sinistra della particella
        else gridX = (int)gridXX;

        gridYY = particles[j].getParticleY() / data.getDeltaY();
        if (gridYY < 0) gridY = ((int)gridYY) - 1;           //il punto di griglia lungo l'asse y a sinistra della particella
        else gridY = (int)gridYY;

        gridZZ = particles[j].getParticleZ() / data.getDeltaZ();
        if (gridZZ < 0) gridZ = ((int)gridZZ) - 1;           //il punto di griglia lungo l'asse z a sinistra della particella
        else gridZ = (int)gridZZ;

        if (gridX < 0 || gridY < 0 || gridZ < 0)
        {
          cout << "La particella e` uscita dai boundaries - simulazione interrotta" << endl;
          return;
        }


        fieldOnGrid.push_back(campoSuPuntiGriglia.at(gridY * data.getNgridPointsX() + gridX));        // primo punto
        fieldOnGrid.push_back(campoSuPuntiGriglia.at(gridY * data.getNgridPointsX() + (gridX + 1)));      // secondo punto
        fieldOnGrid.push_back(campoSuPuntiGriglia.at((gridY + 1) * data.getNgridPointsX() + gridX));      // terzo punto
        fieldOnGrid.push_back(campoSuPuntiGriglia.at((gridY + 1) * data.getNgridPointsX() + (gridX + 1)));    // quarto punto

        FOP = interpolation2D_linear(data, particles[j], fieldOnGrid);
        fieldOnParticles.push_back(FOP);

        fieldOnGrid.clear();



        //pstar = p(n-1/2) + Q * E(n) * DeltaT/2
        pstar.setParticlePX(particles[j].getParticlePX() + 0.5*data.getDeltaT() * particles[j].getParticleQ() * fieldOnParticles[0].getEx());
        pstar.setParticlePY(particles[j].getParticlePY() + 0.5*data.getDeltaT() * particles[j].getParticleQ() * fieldOnParticles[0].getEy());
        pstar.setParticlePZ(particles[j].getParticlePZ() + 0.5*data.getDeltaT() * particles[j].getParticleQ() * fieldOnParticles[0].getEz());
        gamma = sqrt(1 + pow(pstar.getParticlePX(), 2) + pow(pstar.getParticlePY(), 2) + pow(pstar.getParticlePZ(), 2));
        gamma_inv = 1. / gamma;

        //b = (DeltaT/2) * Q * B(n) / gamma
        b.setBx(0.5*data.getDeltaT() * fieldOnParticles[0].getBx() * particles[j].getParticleQ() * gamma_inv);
        b.setBy(0.5*data.getDeltaT() * fieldOnParticles[0].getBy() * particles[j].getParticleQ() * gamma_inv);
        b.setBz(0.5*data.getDeltaT() * fieldOnParticles[0].getBz() * particles[j].getParticleQ() * gamma_inv);
        b2 = pow(b.getBx(), 2) + pow(b.getBy(), 2) + pow(b.getBz(), 2);
        b2 += 1.;

        /****************************************************************************
        //p(n+1/2) = [pstar + pstar x b + (pstar x b) * b] / b2   SECONDO TURCHETTI  (ma i conti analitici diconon che e` sbagliato, anche perche' l'ultimo termine e` uno scalare!)
        pHalfAdvanced.setParticlePX( (pstar.getParticlePX() + pstar.getParticlePY() * b.getBz() - b.getBy() * pstar.getParticlePZ() + (pstar.getParticlePY() * b.getBz() - b.getBy() * pstar.getParticlePZ()) * b.getBx() ) / b2 );
        pHalfAdvanced.setParticlePY( (pstar.getParticlePY() + pstar.getParticlePZ() * b.getBx() - b.getBz() * pstar.getParticlePX() + (pstar.getParticlePZ() * b.getBx() - b.getBz() * pstar.getParticlePX()) * b.getBy() ) / b2 );
        pHalfAdvanced.setParticlePZ( (pstar.getParticlePZ() + pstar.getParticlePX() * b.getBy() - b.getBx() * pstar.getParticlePY() + (pstar.getParticlePX() * b.getBy() - b.getBx() * pstar.getParticlePY()) * b.getBz() ) / b2 );
        ****************************************************************************/

        /****************************************************************************
        //p(n+1/2) = [pstar + pstar x b + b x (pstar  x b)] / b2   SECONDO LONDRILLO
        pHalfAdvanced.setParticlePX( (pstar.getParticlePX() + pstar.getParticlePY() * b.getBz() - b.getBy() * pstar.getParticlePZ() + b.getBx() * pstar.getParticlePX() * b.getBx() ) / b2 );
        pHalfAdvanced.setParticlePY( (pstar.getParticlePY() + pstar.getParticlePZ() * b.getBx() - b.getBz() * pstar.getParticlePX() + b.getBy() * pstar.getParticlePY() * b.getBy() ) / b2 );
        pHalfAdvanced.setParticlePZ( (pstar.getParticlePZ() + pstar.getParticlePX() * b.getBy() - b.getBx() * pstar.getParticlePY() + b.getBz() * pstar.getParticlePZ() * b.getBz() ) / b2 );
        ****************************************************************************/

        //      /****************************************************************************
              //p(n+1/2) = [pstar + pstar x b] / b2   COME IMPLEMENTATO DA TURCHETTI NEL CODICE
        pHalfAdvanced.setParticlePX((pstar.getParticlePX() + pstar.getParticlePY() * b.getBz()) / b2);
        pHalfAdvanced.setParticlePY((pstar.getParticlePY() - pstar.getParticlePX() * b.getBz()) / b2);
        pHalfAdvanced.setParticlePZ(0.);  //per come ho strutturato i pacchetti d'onda su z non dovrebbe succedere nulla e quindi lo fisso io a zero.
  //      ****************************************************************************/


        particles[j].setParticlePX(2 * pHalfAdvanced.getParticlePX() - particles[j].getParticlePX());
        particles[j].setParticlePY(2 * pHalfAdvanced.getParticlePY() - particles[j].getParticlePY());
        particles[j].setParticlePZ(2 * pHalfAdvanced.getParticlePZ() - particles[j].getParticlePZ());

        gamma = sqrt(1. + pow(particles[j].getParticlePX(), 2) + pow(particles[j].getParticlePY(), 2));
        gamma_inv = 1. / gamma;

        particles[j].setParticleX(particles[j].getParticleX() + 0.5 * data.getDeltaT() * particles[j].getParticlePX() * gamma_inv);
        particles[j].setParticleY(particles[j].getParticleY() + 0.5 * data.getDeltaT() * particles[j].getParticlePY() * gamma_inv);
        particles[j].setParticleZ(particles[j].getParticleZ() + 0.5 * data.getDeltaT() * particles[j].getParticlePZ() * gamma_inv);


        //INTEGRALI PRIMI (da sistemare, li ho copiati dall'evoluzione sui campi analitici)
        potenziale.clear();
        zeta = particles[j].getParticleX() - C * data.getT();             // non usato nella mia struttura delle equazioni di campo
        gamma = sqrt(1 + pow(particles[j].getParticlePX(), 2) + pow(particles[j].getParticlePY(), 2) + pow(particles[j].getParticlePZ(), 2));

        gridXX = particles[j].getParticleX() / data.getDeltaX();
        if (gridXX < 0) gridX = ((int)gridXX) - 1;           //il punto di griglia lungo l'asse x a sinistra della particella
        else gridX = (int)gridXX;

        gridYY = particles[j].getParticleY() / data.getDeltaY();
        if (gridYY < 0) gridY = ((int)gridYY) - 1;           //il punto di griglia lungo l'asse y a sinistra della particella
        else gridY = (int)gridYY;

        gridZZ = particles[j].getParticleZ() / data.getDeltaZ();
        if (gridZZ < 0) gridZ = ((int)gridZZ) - 1;           //il punto di griglia lungo l'asse z a sinistra della particella
        else gridZ = (int)gridZZ;

        if (gridX < 0 || gridY < 0 || gridZ < 0)
        {
          cout << "La particella e` uscita dai boundaries - simulazione interrotta" << endl;
          return;
        }


        fieldOnGrid.push_back(potenzialeSuPuntiGriglia.at(gridY * data.getNgridPointsX() + gridX));       // primo punto
        fieldOnGrid.push_back(potenzialeSuPuntiGriglia.at(gridY * data.getNgridPointsX() + (gridX + 1)));     // secondo punto
        fieldOnGrid.push_back(potenzialeSuPuntiGriglia.at((gridY + 1) * data.getNgridPointsX() + gridX));     // terzo punto
        fieldOnGrid.push_back(potenzialeSuPuntiGriglia.at((gridY + 1) * data.getNgridPointsX() + (gridX + 1)));   // quarto punto

        FOP = interpolation2D_linear(data, particles[j], fieldOnGrid);
        potenziale.push_back(FOP);

        fieldOnGrid.clear();

        I1 = (particles[j].getParticlePY() / (particles[j].getParticleM()*C)) + potenziale[0].getEy();
        I2 = particles[j].getParticlePZ() / (particles[j].getParticleM()*C);
        I3 = gamma - particles[j].getParticlePX() / (particles[j].getParticleM()*C);



        outputstream
          << setw(12) << i + 1 << "\t"
          << setprecision(14) << setw(16) << particles[j].getParticleX() << "\t"
          << setw(16) << particles[j].getParticleY() << "\t"
          << setw(16) << particles[j].getParticleZ() << "\t"
          // << setw(16) << particles[j].getParticleCell() << "\t"
          << setw(16) << particles[j].getParticlePX() << "\t"
          << setw(16) << particles[j].getParticlePY() << "\t"
          << setw(16) << particles[j].getParticlePZ() << "\t"
          << setw(16) << zeta << "\t"
          << setw(16) << potenziale[0].getEy() << "\t"
          << setw(16) << I1 << "\t"
          << setw(16) << I2 << "\t"
          << setw(16) << I3 << endl;


        /*
        //UPDATE CAMPI ATTORNO ALLA NUOVA POSIZIONE DELLA PARTICELLA  (al momento non li evolve, il set li pone pari al get precedente)
        gridX = (int) (particles[j].getParticleX()/data.getDeltaX());   //il punto di griglia lungo l'asse x piu` vicino all'origine
        gridY = (int) (particles[j].getParticleY()/data.getDeltaY());   //il punto di griglia lungo l'asse y piu` vicino all'origine
        gridZ = (int) (particles[j].getParticleZ()/data.getDeltaZ());   //il punto di griglia lungo l'asse z piu` vicino all'origine

        campoSuPuntiGriglia.at(gridZ * data.getNgridPointsX() * data.getNgridPointsY() + gridY * data.getNgridPointsX() + gridX).setEx( campoSuPuntiGriglia.at(gridZ * data.getNgridPointsX() * data.getNgridPointsY() + gridY * data.getNgridPointsX() + gridX).getEx() );       // primo punto
        campoSuPuntiGriglia.at(gridZ * data.getNgridPointsX() * data.getNgridPointsY() + gridY * data.getNgridPointsX() + gridX).setEy( campoSuPuntiGriglia.at(gridZ * data.getNgridPointsX() * data.getNgridPointsY() + gridY * data.getNgridPointsX() + gridX).getEy() );
        campoSuPuntiGriglia.at(gridZ * data.getNgridPointsX() * data.getNgridPointsY() + gridY * data.getNgridPointsX() + gridX).setEz( campoSuPuntiGriglia.at(gridZ * data.getNgridPointsX() * data.getNgridPointsY() + gridY * data.getNgridPointsX() + gridX).getEz() );
        campoSuPuntiGriglia.at(gridZ * data.getNgridPointsX() * data.getNgridPointsY() + gridY * data.getNgridPointsX() + gridX).setBx( campoSuPuntiGriglia.at(gridZ * data.getNgridPointsX() * data.getNgridPointsY() + gridY * data.getNgridPointsX() + gridX).getBx() );
        campoSuPuntiGriglia.at(gridZ * data.getNgridPointsX() * data.getNgridPointsY() + gridY * data.getNgridPointsX() + gridX).setBy( campoSuPuntiGriglia.at(gridZ * data.getNgridPointsX() * data.getNgridPointsY() + gridY * data.getNgridPointsX() + gridX).getBy() );
        campoSuPuntiGriglia.at(gridZ * data.getNgridPointsX() * data.getNgridPointsY() + gridY * data.getNgridPointsX() + gridX).setBz( campoSuPuntiGriglia.at(gridZ * data.getNgridPointsX() * data.getNgridPointsY() + gridY * data.getNgridPointsX() + gridX).getBz() );

        campoSuPuntiGriglia.at(gridZ * data.getNgridPointsX() * data.getNgridPointsY() + gridY * data.getNgridPointsX() + (gridX+1)).setEx( campoSuPuntiGriglia.at(gridZ * data.getNgridPointsX() * data.getNgridPointsY() + gridY * data.getNgridPointsX() + (gridX+1)).getEx() );   // secondo punto
        campoSuPuntiGriglia.at(gridZ * data.getNgridPointsX() * data.getNgridPointsY() + gridY * data.getNgridPointsX() + (gridX+1)).setEy( campoSuPuntiGriglia.at(gridZ * data.getNgridPointsX() * data.getNgridPointsY() + gridY * data.getNgridPointsX() + (gridX+1)).getEy() );
        campoSuPuntiGriglia.at(gridZ * data.getNgridPointsX() * data.getNgridPointsY() + gridY * data.getNgridPointsX() + (gridX+1)).setEz( campoSuPuntiGriglia.at(gridZ * data.getNgridPointsX() * data.getNgridPointsY() + gridY * data.getNgridPointsX() + (gridX+1)).getEz() );
        campoSuPuntiGriglia.at(gridZ * data.getNgridPointsX() * data.getNgridPointsY() + gridY * data.getNgridPointsX() + (gridX+1)).setBx( campoSuPuntiGriglia.at(gridZ * data.getNgridPointsX() * data.getNgridPointsY() + gridY * data.getNgridPointsX() + (gridX+1)).getBx() );
        campoSuPuntiGriglia.at(gridZ * data.getNgridPointsX() * data.getNgridPointsY() + gridY * data.getNgridPointsX() + (gridX+1)).setBy( campoSuPuntiGriglia.at(gridZ * data.getNgridPointsX() * data.getNgridPointsY() + gridY * data.getNgridPointsX() + (gridX+1)).getBy() );
        campoSuPuntiGriglia.at(gridZ * data.getNgridPointsX() * data.getNgridPointsY() + gridY * data.getNgridPointsX() + (gridX+1)).setBz( campoSuPuntiGriglia.at(gridZ * data.getNgridPointsX() * data.getNgridPointsY() + gridY * data.getNgridPointsX() + (gridX+1)).getBz() );

        campoSuPuntiGriglia.at(gridZ * data.getNgridPointsX() * data.getNgridPointsY() + (gridY+1) * data.getNgridPointsX() + gridX).setEx( campoSuPuntiGriglia.at(gridZ * data.getNgridPointsX() * data.getNgridPointsY() + (gridY+1) * data.getNgridPointsX() + gridX).getEx() );   // terzo punto
        campoSuPuntiGriglia.at(gridZ * data.getNgridPointsX() * data.getNgridPointsY() + (gridY+1) * data.getNgridPointsX() + gridX).setEy( campoSuPuntiGriglia.at(gridZ * data.getNgridPointsX() * data.getNgridPointsY() + (gridY+1) * data.getNgridPointsX() + gridX).getEy() );
        campoSuPuntiGriglia.at(gridZ * data.getNgridPointsX() * data.getNgridPointsY() + (gridY+1) * data.getNgridPointsX() + gridX).setEz( campoSuPuntiGriglia.at(gridZ * data.getNgridPointsX() * data.getNgridPointsY() + (gridY+1) * data.getNgridPointsX() + gridX).getEz() );
        campoSuPuntiGriglia.at(gridZ * data.getNgridPointsX() * data.getNgridPointsY() + (gridY+1) * data.getNgridPointsX() + gridX).setBx( campoSuPuntiGriglia.at(gridZ * data.getNgridPointsX() * data.getNgridPointsY() + (gridY+1) * data.getNgridPointsX() + gridX).getBx() );
        campoSuPuntiGriglia.at(gridZ * data.getNgridPointsX() * data.getNgridPointsY() + (gridY+1) * data.getNgridPointsX() + gridX).setBy( campoSuPuntiGriglia.at(gridZ * data.getNgridPointsX() * data.getNgridPointsY() + (gridY+1) * data.getNgridPointsX() + gridX).getBy() );
        campoSuPuntiGriglia.at(gridZ * data.getNgridPointsX() * data.getNgridPointsY() + (gridY+1) * data.getNgridPointsX() + gridX).setBz( campoSuPuntiGriglia.at(gridZ * data.getNgridPointsX() * data.getNgridPointsY() + (gridY+1) * data.getNgridPointsX() + gridX).getBz() );

        campoSuPuntiGriglia.at(gridZ * data.getNgridPointsX() * data.getNgridPointsY() + (gridY+1) * data.getNgridPointsX() + (gridX+1)).setEx( campoSuPuntiGriglia.at(gridZ * data.getNgridPointsX() * data.getNgridPointsY() + (gridY+1) * data.getNgridPointsX() + (gridX+1)).getEx() );   // quarto punto
        campoSuPuntiGriglia.at(gridZ * data.getNgridPointsX() * data.getNgridPointsY() + (gridY+1) * data.getNgridPointsX() + (gridX+1)).setEy( campoSuPuntiGriglia.at(gridZ * data.getNgridPointsX() * data.getNgridPointsY() + (gridY+1) * data.getNgridPointsX() + (gridX+1)).getEy() );
        campoSuPuntiGriglia.at(gridZ * data.getNgridPointsX() * data.getNgridPointsY() + (gridY+1) * data.getNgridPointsX() + (gridX+1)).setEz( campoSuPuntiGriglia.at(gridZ * data.getNgridPointsX() * data.getNgridPointsY() + (gridY+1) * data.getNgridPointsX() + (gridX+1)).getEz() );
        campoSuPuntiGriglia.at(gridZ * data.getNgridPointsX() * data.getNgridPointsY() + (gridY+1) * data.getNgridPointsX() + (gridX+1)).setBx( campoSuPuntiGriglia.at(gridZ * data.getNgridPointsX() * data.getNgridPointsY() + (gridY+1) * data.getNgridPointsX() + (gridX+1)).getBx() );
        campoSuPuntiGriglia.at(gridZ * data.getNgridPointsX() * data.getNgridPointsY() + (gridY+1) * data.getNgridPointsX() + (gridX+1)).setBy( campoSuPuntiGriglia.at(gridZ * data.getNgridPointsX() * data.getNgridPointsY() + (gridY+1) * data.getNgridPointsX() + (gridX+1)).getBy() );
        campoSuPuntiGriglia.at(gridZ * data.getNgridPointsX() * data.getNgridPointsY() + (gridY+1) * data.getNgridPointsX() + (gridX+1)).setBz( campoSuPuntiGriglia.at(gridZ * data.getNgridPointsX() * data.getNgridPointsY() + (gridY+1) * data.getNgridPointsX() + (gridX+1)).getBz() );

        campoSuPuntiGriglia.at((gridZ+1) * data.getNgridPointsX() * data.getNgridPointsY() + gridY * data.getNgridPointsX() + gridX).setEx( campoSuPuntiGriglia.at((gridZ+1) * data.getNgridPointsX() * data.getNgridPointsY() + gridY * data.getNgridPointsX() + gridX).getEx() );   // quinto punto
        campoSuPuntiGriglia.at((gridZ+1) * data.getNgridPointsX() * data.getNgridPointsY() + gridY * data.getNgridPointsX() + gridX).setEy( campoSuPuntiGriglia.at((gridZ+1) * data.getNgridPointsX() * data.getNgridPointsY() + gridY * data.getNgridPointsX() + gridX).getEy() );
        campoSuPuntiGriglia.at((gridZ+1) * data.getNgridPointsX() * data.getNgridPointsY() + gridY * data.getNgridPointsX() + gridX).setEz( campoSuPuntiGriglia.at((gridZ+1) * data.getNgridPointsX() * data.getNgridPointsY() + gridY * data.getNgridPointsX() + gridX).getEz() );
        campoSuPuntiGriglia.at((gridZ+1) * data.getNgridPointsX() * data.getNgridPointsY() + gridY * data.getNgridPointsX() + gridX).setBx( campoSuPuntiGriglia.at((gridZ+1) * data.getNgridPointsX() * data.getNgridPointsY() + gridY * data.getNgridPointsX() + gridX).getBx() );
        campoSuPuntiGriglia.at((gridZ+1) * data.getNgridPointsX() * data.getNgridPointsY() + gridY * data.getNgridPointsX() + gridX).setBy( campoSuPuntiGriglia.at((gridZ+1) * data.getNgridPointsX() * data.getNgridPointsY() + gridY * data.getNgridPointsX() + gridX).getBy() );
        campoSuPuntiGriglia.at((gridZ+1) * data.getNgridPointsX() * data.getNgridPointsY() + gridY * data.getNgridPointsX() + gridX).setBz( campoSuPuntiGriglia.at((gridZ+1) * data.getNgridPointsX() * data.getNgridPointsY() + gridY * data.getNgridPointsX() + gridX).getBz() );

        campoSuPuntiGriglia.at((gridZ+1) * data.getNgridPointsX() * data.getNgridPointsY() + gridY * data.getNgridPointsX() + (gridX+1)).setEx( campoSuPuntiGriglia.at((gridZ+1) * data.getNgridPointsX() * data.getNgridPointsY() + gridY * data.getNgridPointsX() + (gridX+1)).getEx() );   // sesto punto
        campoSuPuntiGriglia.at((gridZ+1) * data.getNgridPointsX() * data.getNgridPointsY() + gridY * data.getNgridPointsX() + (gridX+1)).setEy( campoSuPuntiGriglia.at((gridZ+1) * data.getNgridPointsX() * data.getNgridPointsY() + gridY * data.getNgridPointsX() + (gridX+1)).getEy() );
        campoSuPuntiGriglia.at((gridZ+1) * data.getNgridPointsX() * data.getNgridPointsY() + gridY * data.getNgridPointsX() + (gridX+1)).setEz( campoSuPuntiGriglia.at((gridZ+1) * data.getNgridPointsX() * data.getNgridPointsY() + gridY * data.getNgridPointsX() + (gridX+1)).getEz() );
        campoSuPuntiGriglia.at((gridZ+1) * data.getNgridPointsX() * data.getNgridPointsY() + gridY * data.getNgridPointsX() + (gridX+1)).setBx( campoSuPuntiGriglia.at((gridZ+1) * data.getNgridPointsX() * data.getNgridPointsY() + gridY * data.getNgridPointsX() + (gridX+1)).getBx() );
        campoSuPuntiGriglia.at((gridZ+1) * data.getNgridPointsX() * data.getNgridPointsY() + gridY * data.getNgridPointsX() + (gridX+1)).setBy( campoSuPuntiGriglia.at((gridZ+1) * data.getNgridPointsX() * data.getNgridPointsY() + gridY * data.getNgridPointsX() + (gridX+1)).getBy() );
        campoSuPuntiGriglia.at((gridZ+1) * data.getNgridPointsX() * data.getNgridPointsY() + gridY * data.getNgridPointsX() + (gridX+1)).setBz( campoSuPuntiGriglia.at((gridZ+1) * data.getNgridPointsX() * data.getNgridPointsY() + gridY * data.getNgridPointsX() + (gridX+1)).getBz() );

        campoSuPuntiGriglia.at((gridZ+1) * data.getNgridPointsX() * data.getNgridPointsY() + (gridY+1) * data.getNgridPointsX() + gridX).setEx( campoSuPuntiGriglia.at((gridZ+1) * data.getNgridPointsX() * data.getNgridPointsY() + (gridY+1) * data.getNgridPointsX() + gridX).getEx() );   // settimo punto
        campoSuPuntiGriglia.at((gridZ+1) * data.getNgridPointsX() * data.getNgridPointsY() + (gridY+1) * data.getNgridPointsX() + gridX).setEy( campoSuPuntiGriglia.at((gridZ+1) * data.getNgridPointsX() * data.getNgridPointsY() + (gridY+1) * data.getNgridPointsX() + gridX).getEy() );
        campoSuPuntiGriglia.at((gridZ+1) * data.getNgridPointsX() * data.getNgridPointsY() + (gridY+1) * data.getNgridPointsX() + gridX).setEz( campoSuPuntiGriglia.at((gridZ+1) * data.getNgridPointsX() * data.getNgridPointsY() + (gridY+1) * data.getNgridPointsX() + gridX).getEz() );
        campoSuPuntiGriglia.at((gridZ+1) * data.getNgridPointsX() * data.getNgridPointsY() + (gridY+1) * data.getNgridPointsX() + gridX).setBx( campoSuPuntiGriglia.at((gridZ+1) * data.getNgridPointsX() * data.getNgridPointsY() + (gridY+1) * data.getNgridPointsX() + gridX).getBx() );
        campoSuPuntiGriglia.at((gridZ+1) * data.getNgridPointsX() * data.getNgridPointsY() + (gridY+1) * data.getNgridPointsX() + gridX).setBy( campoSuPuntiGriglia.at((gridZ+1) * data.getNgridPointsX() * data.getNgridPointsY() + (gridY+1) * data.getNgridPointsX() + gridX).getBy() );
        campoSuPuntiGriglia.at((gridZ+1) * data.getNgridPointsX() * data.getNgridPointsY() + (gridY+1) * data.getNgridPointsX() + gridX).setBz( campoSuPuntiGriglia.at((gridZ+1) * data.getNgridPointsX() * data.getNgridPointsY() + (gridY+1) * data.getNgridPointsX() + gridX).getBz() );

        campoSuPuntiGriglia.at((gridZ+1) * data.getNgridPointsX() * data.getNgridPointsY() + (gridY+1) * data.getNgridPointsX() + (gridX+1)).setEx( campoSuPuntiGriglia.at((gridZ+1) * data.getNgridPointsX() * data.getNgridPointsY() + (gridY+1) * data.getNgridPointsX() + (gridX+1)).getEx() );   // ottavo punto
        campoSuPuntiGriglia.at((gridZ+1) * data.getNgridPointsX() * data.getNgridPointsY() + (gridY+1) * data.getNgridPointsX() + (gridX+1)).setEy( campoSuPuntiGriglia.at((gridZ+1) * data.getNgridPointsX() * data.getNgridPointsY() + (gridY+1) * data.getNgridPointsX() + (gridX+1)).getEx() );
        campoSuPuntiGriglia.at((gridZ+1) * data.getNgridPointsX() * data.getNgridPointsY() + (gridY+1) * data.getNgridPointsX() + (gridX+1)).setEz( campoSuPuntiGriglia.at((gridZ+1) * data.getNgridPointsX() * data.getNgridPointsY() + (gridY+1) * data.getNgridPointsX() + (gridX+1)).getEx() );
        campoSuPuntiGriglia.at((gridZ+1) * data.getNgridPointsX() * data.getNgridPointsY() + (gridY+1) * data.getNgridPointsX() + (gridX+1)).setBx( campoSuPuntiGriglia.at((gridZ+1) * data.getNgridPointsX() * data.getNgridPointsY() + (gridY+1) * data.getNgridPointsX() + (gridX+1)).getEx() );
        campoSuPuntiGriglia.at((gridZ+1) * data.getNgridPointsX() * data.getNgridPointsY() + (gridY+1) * data.getNgridPointsX() + (gridX+1)).setBy( campoSuPuntiGriglia.at((gridZ+1) * data.getNgridPointsX() * data.getNgridPointsY() + (gridY+1) * data.getNgridPointsX() + (gridX+1)).getEx() );
        campoSuPuntiGriglia.at((gridZ+1) * data.getNgridPointsX() * data.getNgridPointsY() + (gridY+1) * data.getNgridPointsX() + (gridX+1)).setBz( campoSuPuntiGriglia.at((gridZ+1) * data.getNgridPointsX() * data.getNgridPointsY() + (gridY+1) * data.getNgridPointsX() + (gridX+1)).getEx() );
        */
      }
    }
  }
}


void evolveRK4(Data data, Particle particle, Field fieldOnParticle, vector<Particle> &dx2, vector<Particle> &x2)
{
  Particle dx, dx2a, x2a;
  double gamma = sqrt(1 + pow(particle.getParticlePX(), 2) + pow(particle.getParticlePY(), 2) + pow(particle.getParticlePZ(), 2));

  dx.setParticleX(particle.getParticlePX() / (gamma * particle.getParticleM()));
  dx.setParticleY(particle.getParticlePY() / (gamma * particle.getParticleM()));
  dx.setParticleZ(particle.getParticlePZ() / (gamma * particle.getParticleM()));
  dx.setParticlePX(particle.getParticleQ() * fieldOnParticle.getEx() + particle.getParticleQ() * dx.getParticleY() * fieldOnParticle.getBz() - particle.getParticleQ() * dx.getParticleZ() * fieldOnParticle.getBy());  //CONTROLLARE I SEGNI!!!!
  dx.setParticlePY(particle.getParticleQ() * fieldOnParticle.getEy() + particle.getParticleQ() * dx.getParticleZ() * fieldOnParticle.getBx() - particle.getParticleQ() * dx.getParticleX() * fieldOnParticle.getBz());
  dx.setParticlePZ(particle.getParticleQ() * fieldOnParticle.getEz() + particle.getParticleQ() * dx.getParticleX() * fieldOnParticle.getBy() - particle.getParticleQ() * dx.getParticleY() * fieldOnParticle.getBx());

  //RIEMPIE dx2 per il metodo RK4 con il dx del primo passo
  dx2a = dx;
  dx2.push_back(dx2a);


  //RIEMPIE x2 con le posizioni avanzate di particles[i] dopo il primo passo
  x2a.setParticleX(particle.getParticleX() + dx.getParticleX() * (data.getDeltaT()*0.5));
  x2a.setParticleY(particle.getParticleY() + dx.getParticleY() * (data.getDeltaT()*0.5));
  x2a.setParticleZ(particle.getParticleZ() + dx.getParticleZ() * (data.getDeltaT()*0.5));
  x2a.setParticlePX(particle.getParticlePX() + dx.getParticlePX() *  (data.getDeltaT()*0.5));
  x2a.setParticlePY(particle.getParticlePY() + dx.getParticlePY() *  (data.getDeltaT()*0.5));
  x2a.setParticlePZ(particle.getParticlePZ() + dx.getParticlePZ() *  (data.getDeltaT()*0.5));

  x2.push_back(x2a);
}


void evolveLPF(Data data, Particle &particle, double dtp, double dtx)
{
  double gamma, gamma_inv, b2;
  Particle pstar, pHalfAdvanced;

  Field b;
  Field *campoSuPunto;
  campoSuPunto = new Field[1];

  gamma = sqrt(1. + pow(particle.getParticlePX(), 2) + pow(particle.getParticlePY(), 2) + pow(particle.getParticlePZ(), 2));
  gamma_inv = 1.0 / gamma;

  particle.setParticleX(particle.getParticleX() + dtx * particle.getParticlePX() * gamma_inv);
  particle.setParticleY(particle.getParticleY() + dtx * particle.getParticlePY() * gamma_inv);
  particle.setParticleZ(particle.getParticleZ() + dtx * particle.getParticlePZ() * gamma_inv);

  data.aumentaT(dtx);

  calcolaCampoAnalitico1DSuParticella(data, particle, *campoSuPunto);

  //pstar = p(n) + Q * E(n+1/2) * DeltaT/2
  pstar.setParticlePX(particle.getParticlePX() + 0.5 * dtp * particle.getParticleQ() * (*campoSuPunto).getEx());
  pstar.setParticlePY(particle.getParticlePY() + 0.5 * dtp * particle.getParticleQ() * (*campoSuPunto).getEy());
  pstar.setParticlePZ(particle.getParticlePZ() + 0.5 * dtp * particle.getParticleQ() * (*campoSuPunto).getEz());
  gamma = sqrt(1. + pow(pstar.getParticlePX(), 2) + pow(pstar.getParticlePY(), 2) + pow(pstar.getParticlePZ(), 2));
  gamma_inv = 1.0 / gamma;

  //b = (DeltaT/2) * Q * B(n+1/2) / gamma(aggiornato al pstar)
  b.setBx(0.5 * dtp * (*campoSuPunto).getBx() * particle.getParticleQ() * gamma_inv);
  b.setBy(0.5 * dtp * (*campoSuPunto).getBy() * particle.getParticleQ() * gamma_inv);
  b.setBz(0.5 * dtp * (*campoSuPunto).getBz() * particle.getParticleQ() * gamma_inv);
  b2 = pow(b.getBx(), 2) + pow(b.getBy(), 2) + pow(b.getBz(), 2);
  b2 += 1.0;

  /****************************************************************************
  //p(n+1/2) = [pstar + pstar x b + (pstar x b) * b] / b2   SECONDO TURCHETTI  (ma i conti analitici diconon che e` sbagliato, anche perche' l'ultimo termine e` uno scalare!)
  pHalfAdvanced.setParticlePX( (pstar.getParticlePX() + pstar.getParticlePY() * b.getBz() - b.getBy() * pstar.getParticlePZ() + (pstar.getParticlePY() * b.getBz() - b.getBy() * pstar.getParticlePZ()) * b.getBx() ) / b2 );
  pHalfAdvanced.setParticlePY( (pstar.getParticlePY() + pstar.getParticlePZ() * b.getBx() - b.getBz() * pstar.getParticlePX() + (pstar.getParticlePZ() * b.getBx() - b.getBz() * pstar.getParticlePX()) * b.getBy() ) / b2 );
  pHalfAdvanced.setParticlePZ( (pstar.getParticlePZ() + pstar.getParticlePX() * b.getBy() - b.getBx() * pstar.getParticlePY() + (pstar.getParticlePX() * b.getBy() - b.getBx() * pstar.getParticlePY()) * b.getBz() ) / b2 );
  ****************************************************************************/

  //  /****************************************************************************
    //p(n+1/2) = [pstar + pstar x b + b x (pstar x b)] / b2   SECONDO LONDRILLO
  pHalfAdvanced.setParticlePX((pstar.getParticlePX() + pstar.getParticlePY() * b.getBz() - b.getBy() * pstar.getParticlePZ() + b.getBx() * pstar.getParticlePX() * b.getBx()) / b2);
  pHalfAdvanced.setParticlePY((pstar.getParticlePY() + pstar.getParticlePZ() * b.getBx() - b.getBz() * pstar.getParticlePX() + b.getBy() * pstar.getParticlePY() * b.getBy()) / b2);
  pHalfAdvanced.setParticlePZ((pstar.getParticlePZ() + pstar.getParticlePX() * b.getBy() - b.getBx() * pstar.getParticlePY() + b.getBz() * pstar.getParticlePZ() * b.getBz()) / b2);
  //  ****************************************************************************/

    /****************************************************************************
    //p(n+1/2) = [pstar + pstar x b] / b2   COME IMPLEMENTATO DA TURCHETTI NEL CODICE
    pHalfAdvanced.setParticlePX( (pstar.getParticlePX() + pstar.getParticlePY() * b.getBz()) / b2 );
    pHalfAdvanced.setParticlePY( (pstar.getParticlePY() - pstar.getParticlePX() * b.getBz()) / b2 );
    pHalfAdvanced.setParticlePZ( 0. );  //per come ho strutturato i pacchetti d'onda su z non dovrebbe succedere nulla e quindi lo fisso io a zero.
    ****************************************************************************/

  particle.setParticlePX(2.0 * pHalfAdvanced.getParticlePX() - particle.getParticlePX());
  particle.setParticlePY(2.0 * pHalfAdvanced.getParticlePY() - particle.getParticlePY());
  particle.setParticlePZ(2.0 * pHalfAdvanced.getParticlePZ() - particle.getParticlePZ());
}
