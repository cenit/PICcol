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



#include <iostream>
#include <iomanip>
#include <cstdio>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <fstream>
#include <vector>
#include <limits>
#include <ios>

#include "datatypes.h"      // include il file con la definizione dei tipi degli oggetti: Particle, Field, Data
#include "evolution.h"      // include i metodi di evoluzione del sistema
#include "filler.h"         // include i metodi di riempimento dei vettori di particelle e dei campi, iniziali e durante l'elaborazione

using namespace std;


int main(int argc, char*argv[])
{
  /************************************************************
  * INIZIALIZZAZIONE generatore random tramite seed temporale *
  ************************************************************/
  time_t now;
  time(&now);
  srand((unsigned int)now);

  /************************************************************
  * INIZIALIZZAZIONE file di output                           *
  * - outputLOG conterra` i parametri dell'esecuzione corrente *
  *   ed apre in modalita` append il file log.txt              *
  * - outputDATA conterra` invece tutti i dati delle           *
  *   particelle istante per istante e va a sovrascrivere ad  *
  *   ogni esecuzione del programma il file output.dat        *
  ************************************************************/
  ofstream outputLOG("log.txt", fstream::app);
  ofstream outputDATA("output.dat");
  outputLOG << "------------------------------------------------------" << endl;
  auto t = std::time(nullptr);
  auto tm = *std::localtime(&t);
  outputLOG << "Data e ora: " << std::put_time(&tm, "%d-%m-%Y %H-%M-%S") << endl;

  /************************************************************
  * INIZIALIZZAZIONE fondamenti del codice:                   *
  * - "data", di tipo Data, contiene tutti i parametri che    *
  *   possono servire ad una funzione di questo codice per    *
  *   sapere la configurazione del sistema in esame           *
  * - un vettore di particelle che raggruppa le stesse        *
  * - un vettore di campi calcolati nelle posizioni delle     *
  *   particelle                                              *
  * - un vettore di campi calcolati nei nodi della griglia    *
  *   definita dall'utente e sovrimposta al sistema (celle)   *
  ************************************************************/
  Data data;
  vector<Particle> particelle;
  vector<Field> campiSuParticelle;
  vector<Field> campiSuGriglia;

  /************************************************************
  * L'eseguibile generato deve essere in grado di ricevere    *
  * dall'utente in fase di runtime tutti i parametri di       *
  * lavoro o le coordinate delle particelle da utilizzare,    *
  * oppure anche essere in grado di leggere dei file di input *
  * che possono essere di vario tipo:                         *
  * - parametri del sistema, da utilizzarli per generare lo   *
  *   stesso (procedura completamente automatizzata)          *
  * - dati riguardo i valori dei campi e le coordinate delle  *
  *   particelle da utilizzare durante i calcoli; se non      *
  *   legge anche i parametri deve naturalmente calcolarsi    *
  *   quelli ovvi e chidere quelli mancanti                   *
  * - una combinazione dei precedenti (con controlli di       *
  *   coerenza                                                *
  ************************************************************/
  if (argc == 2)
  {
    data.setInputType(argv[1]);
    if (data.getInputType() == 1)       // data.getInputType == 1  --> l'utente ha definito che il file contiene parametri
    {
      data.fillParametersFromFile(argv[1], outputLOG);
    }
    else if (data.getInputType() == 2)      // data.getInputType == 2  --> 2: lettura particelle gia` generate (no campi)
    {
      fillParticlesFromFile(argv[2], outputLOG);
      return 3;
    }
    else return -249;
  }
  else if (argc == 3)
  {
    data.setInputType(argv[1]);
    data.fillParametersFromFile(argv[1], outputLOG);    // se vengono specificati due file, il primo e` dei parametri
    fillParticlesFromFile(argv[2], outputLOG);          // e il secondo delle coordinate delle particelle
    return 4;
  }
  else if (argc > 3)
  {
    cout << "Troppi parametri: eseguire senza parametri oppure specificare il nome di un file contenente i parametri o le particelle" << endl;
    return -255;
  }
  else
  {
    outputLOG << "Parametri impostati dall'utente in fase di runtime" << endl;
    data.setEvoType();
    if (data.getEvoType() == 1)
    {
      cout << "Verranno ora chiesti dati sull'asse x; la griglia costruita non viene utilizzata, \nserve solo per inserire inutilmente le particelle in una cella" << endl;
      data.setNdim();
      data.setXaxis(3);
      if (data.getNdim() == 2 || data.getNdim() == 3) data.setYaxis(3);
      if (data.getNdim() == 3) data.setZaxis(3);
      data.setDeltaT(0.);
      data.setNsteps(1);
      data.setA();
      data.setK();
      data.impostaT(0);
      data.setNelectrons();
      data.setThermalV(0.);
      //      data.setAnalyticType();
      data.setParticleFillingMethod();
    }


    else if (data.getEvoType() == 2)
    {
      //IMPOSTAZIONI NON CORRETTE, RICHIESTE SOLO PER LA MANCANZA DELL'IMPLEMENTAZIONE DEI CAMPI SU GRIGLIA!
      data.setNdim();
      data.setXaxis(3);
      if (data.getNdim() == 2 || data.getNdim() == 3) data.setYaxis(3);
      if (data.getNdim() == 3) data.setZaxis(3);
      data.setDeltaT(0.);
      data.setNsteps(1);
      data.setA();
      data.setK();
      data.impostaT(0);
      data.setNelectrons();
      data.setThermalV(0.);
      data.setParticleFillingMethod();
      //------------------------------------- VERIFICARE QUALI ANDRANNO TENUTE QUANDO CI SARA' UN CAMPO SU GRIGLIA!
    }


    else if (data.getEvoType() == 3)
    {
      //SONO RICHIESTE LE STESSE COSE DEL CASO 2 NEL QUALE NON C'ERA UN VERO CAMPO SU GRIGLIA; FORSE E' SBAGLIATA ANCHE QUESTA VERSIONE!
      data.setNdim();
      data.setXaxis(3);
      if (data.getNdim() == 2 || data.getNdim() == 3) data.setYaxis(3);
      if (data.getNdim() == 3) data.setZaxis(3);
      data.setDeltaT(0.);
      data.setNsteps(1);
      data.setA();
      data.setK();
      data.impostaT(0);
      data.setNelectrons();
      data.setThermalV(0.);
      data.setParticleFillingMethod();
    }


    else if (data.getEvoType() == 4)
    {
      //SONO RICHIESTE LE STESSE COSE DEL CASO 2 NEL QUALE NON C'ERA UN VERO CAMPO SU GRIGLIA; FORSE E' SBAGLIATA ANCHE QUESTA VERSIONE!
      data.setNdim();
      data.setXaxis(3);
      if (data.getNdim() == 2 || data.getNdim() == 3) data.setYaxis(3);
      if (data.getNdim() == 3) data.setZaxis(3);
      data.setDeltaT(0.);
      data.setNsteps(1);
      data.setA();
      data.setK();
      data.impostaT(0);
      data.setNelectrons();
      data.setThermalV(0.);
      data.setParticleFillingMethod();
    }


    else if (data.getEvoType() == 5)
    {
      //SONO RICHIESTE LE STESSE COSE DEL CASO 2 NEL QUALE NON C'ERA UN VERO CAMPO SU GRIGLIA; FORSE E' SBAGLIATA ANCHE QUESTA VERSIONE!
      data.setNdim();
      data.setXaxis(3);
      if (data.getNdim() == 2 || data.getNdim() == 3) data.setYaxis(3);
      if (data.getNdim() == 3) data.setZaxis(3);
      data.setDeltaT(0.);
      data.setNsteps(1);
      data.setA();
      data.setK();
      data.impostaT(0);
      data.setNelectrons();
      data.setThermalV(0.);
      data.setParticleFillingMethod();
    }


    else if (data.getEvoType() == 6)
    {
      //SONO RICHIESTE LE STESSE COSE DEL CASO 2 NEL QUALE NON C'ERA UN VERO CAMPO SU GRIGLIA; FORSE E' SBAGLIATA ANCHE QUESTA VERSIONE!
      data.setNdim();
      data.setXaxis(3);
      if (data.getNdim() == 2 || data.getNdim() == 3) data.setYaxis(3);
      if (data.getNdim() == 3) data.setZaxis(3);
      data.setDeltaT(0.);
      data.setNsteps(1);
      data.setA();
      data.setK();
      data.impostaT(0);
      data.setNelectrons();
      data.setThermalV(0.);
      data.setParticleFillingMethod();
    }


    else if (data.getEvoType() == 7)
    {
      //SONO RICHIESTE LE STESSE COSE DEL CASO 2 NEL QUALE NON C'ERA UN VERO CAMPO SU GRIGLIA; FORSE E' SBAGLIATA ANCHE QUESTA VERSIONE!
      data.setNdim();
      data.setXaxis(3);
      if (data.getNdim() == 2 || data.getNdim() == 3) data.setYaxis(3);
      if (data.getNdim() == 3) data.setZaxis(3);
      data.setDeltaT(0.);
      data.setNsteps(1);
      data.setA();
      data.setK();
      data.impostaT(0);
      data.setNelectrons();
      data.setThermalV(0.);
      data.setParticleFillingMethod();
    }


    else if (data.getEvoType() == 8)
    {
      //SONO RICHIESTE LE STESSE COSE DEL CASO 2 NEL QUALE NON C'ERA UN VERO CAMPO SU GRIGLIA; FORSE E' SBAGLIATA ANCHE QUESTA VERSIONE!
      data.setNdim();
      data.setXaxis(3);
      if (data.getNdim() == 2 || data.getNdim() == 3) data.setYaxis(3);
      if (data.getNdim() == 3) data.setZaxis(3);
      data.setDeltaT(0.);
      data.setNsteps(1);
      data.setA();
      data.setK();
      data.impostaT(0);
      data.setNelectrons();
      data.setThermalV(0.);
      data.setParticleFillingMethod();
    }


    else if (data.getEvoType() == 9)
    {
      //SONO RICHIESTE LE STESSE COSE DEL CASO 2 NEL QUALE NON C'ERA UN VERO CAMPO SU GRIGLIA; FORSE E' SBAGLIATA ANCHE QUESTA VERSIONE!
      data.setNdim();
      data.setXaxis(3);
      if (data.getNdim() == 2 || data.getNdim() == 3) data.setYaxis(3);
      if (data.getNdim() == 3) data.setZaxis(3);
      data.setDeltaT(0.);
      data.setNsteps(1);
      data.setA();
      data.setK();
      data.impostaT(0);
      data.setNelectrons();
      data.setThermalV(0.);
      data.setParticleFillingMethod();
    }


    else if (data.getEvoType() == 10)
    {
      //SONO RICHIESTE LE STESSE COSE DEL CASO 2 NEL QUALE NON C'ERA UN VERO CAMPO SU GRIGLIA; FORSE E' SBAGLIATA ANCHE QUESTA VERSIONE!
      data.setNdim();
      data.setXaxis(3);
      if (data.getNdim() == 2 || data.getNdim() == 3) data.setYaxis(3);
      if (data.getNdim() == 3) data.setZaxis(3);
      data.setDeltaT(0.);
      data.setNsteps(1);
      data.setA();
      data.setK();
      data.impostaT(0);
      data.setNelectrons();
      data.setThermalV(0.);
      data.setParticleFillingMethod();
    }


    else if (data.getEvoType() == 11)
    {
      //SONO RICHIESTE LE STESSE COSE DEL CASO 2 NEL QUALE NON C'ERA UN VERO CAMPO SU GRIGLIA; FORSE E' SBAGLIATA ANCHE QUESTA VERSIONE!
      data.setNdim();
      data.setXaxis(3);
      if (data.getNdim() == 2 || data.getNdim() == 3) data.setYaxis(3);
      if (data.getNdim() == 3) data.setZaxis(3);
      data.setDeltaT(0.);
      data.setNsteps(1);
      data.setA();
      data.setK();
      data.impostaT(0);
      data.setNelectrons();
      data.setThermalV(0.);
      data.setParticleFillingMethod();
    }


    else if (data.getEvoType() == 12)
    {
      //SONO RICHIESTE LE STESSE COSE DEL CASO 2 NEL QUALE NON C'ERA UN VERO CAMPO SU GRIGLIA; FORSE E' SBAGLIATA ANCHE QUESTA VERSIONE!
      data.setNdim();
      data.setXaxis(3);
      if (data.getNdim() == 2 || data.getNdim() == 3) data.setYaxis(3);
      if (data.getNdim() == 3) data.setZaxis(3);
      data.setDeltaT(0.);
      data.setNsteps(1);
      data.setA();
      data.setK();
      data.impostaT(0);
      data.setNelectrons();
      data.setThermalV(0.);
      data.setParticleFillingMethod();
    }


    else
    {
      cout << "Evoluzione non riconosciuta" << endl;
      outputLOG << "Evoluzione non riconosciuta" << endl;
      return -10;
    }


    outputLOG << "Numero particelle da simulare: " << data.getNelectrons() << endl;
    outputLOG << "Sistema " << data.getNdim() << "D" << endl;
    outputLOG << "Velocita' termica: " << data.getThermalV() << endl;
    outputLOG << "Parametri asse x: dimX: " << data.getDimX() << ", num punti griglia: " << data.getNgridPointsX() << endl;
    if (data.getNdim() == 2 || data.getNdim() == 3) outputLOG << "Parametri asse y: dimY: " << data.getDimY() << ", num punti griglia: " << data.getNgridPointsY() << endl;
    if (data.getNdim() == 3) outputLOG << "Parametri asse z: dimZ: " << data.getDimZ() << ", num punti griglia: " << data.getNgridPointsZ() << endl;
    outputLOG << "t iniziale: " << data.getT() << endl;
    outputLOG << "dt: " << data.getDeltaT() << endl;
    outputLOG << "nSteps: " << data.getNsteps() << endl;
    if (data.getParticleFillingMethod() == 1) outputLOG << "Generazione random particelle" << endl;
    else if (data.getParticleFillingMethod() == 2) outputLOG << "Particelle generate in 0,0,0" << endl;

    if (data.getEvoType() == 1)  outputLOG << "Distribuzione campi su particelle in un istante definito" << endl;
    else if (data.getEvoType() == 2)  outputLOG << "Evoluzione campo su una particella nel tempo" << endl;
    else if (data.getEvoType() == 3)  outputLOG << "Evoluzione RungeKutta4_no grid" << endl;
    else if (data.getEvoType() == 4)  outputLOG << "Evoluzione RungeKutta4_grid on the fly" << endl;
    else if (data.getEvoType() == 5)  outputLOG << "Evoluzione RungeKutta4_grid persistent" << endl;
    else if (data.getEvoType() == 6)  outputLOG << "Evoluzione Leapfrog2_no grid" << endl;
    else if (data.getEvoType() == 7)  outputLOG << "Evoluzione Leapfrog2_grid on the fly" << endl;
    else if (data.getEvoType() == 8)  outputLOG << "Evoluzione Leapfrog2_grid persistent" << endl;
    else if (data.getEvoType() == 9)  outputLOG << "Evoluzione Leapfrog4_no grid" << endl;
    else if (data.getEvoType() == 10) outputLOG << "Evoluzione Leapfrog4_grid on the fly" << endl;
    else if (data.getEvoType() == 11) outputLOG << "Evoluzione Leapfrog4_grid persistent" << endl;
    else outputLOG << "Metodo non riconosciuto" << endl;
  }

  /*************************************************************
  * Ora il programma dovrebbe aver saputo "tutto" dall'utente  *
  * e pertanto genera il vettore iniziale delle particelle     *
  * seguendo i parametri impostati e solo nel caso non fosse   *
  * stato gia` definito e letto un file contenente questi dati *
  *************************************************************/
  if (argc != 3) creaVettoreParticelle(data, particelle);


  /*************************************************************
  * In base alle scelte che l'utente ha fatto, il programma:   *
  * - calcola l'evoluzione con un metodo runge-kutta del       *
  *   quarto ordine e un campo analitico predefinito           *
  *   (modificabili solo i parametri) calcolato direttamente   *
  *   nelle posizioni delle particelle                         *
  * - calcola l'evoluzione con un metodo leapfrog e un campo   *
  *   analitico predefinito (modificabili solo i parametri)    *
  *   calcolato nei punti di griglia circostanti la            *
  *   particella ed interpolato su di essa                     *
  * - calcola l'evoluzione con un metodo leapfrog e un campo   *
  *   predeterminato sui punti di griglia, che a sua volta si  *
  *   evolve secondo precise regole (ad ora, i punti di        *
  *   griglia sono riempiti inizialmente con i valori del      *
  *   solito campo analitico in quelle posizioni ed evoluti    *
  *   ponendoli uguali al passo precedente, senza correnti!)   *
  * - altri metodi ancora non descritti in questo commento     *
  *************************************************************/
  if (data.getEvoType() == 1)  calcolaCampoInUnIstante(data, particelle, campiSuParticelle, outputDATA);
  else if (data.getEvoType() == 2)  calcolaEvoluzioneCampoSuParticelle(data, particelle, campiSuParticelle, outputDATA);
  else if (data.getEvoType() == 3)  evolveRK4_nogrid(data, particelle, campiSuParticelle, outputDATA);
  else if (data.getEvoType() == 4)  evolveRK4_withgrid_onthefly(data, particelle, campiSuParticelle, outputDATA);
  else if (data.getEvoType() == 5)  evolveRK4_withgrid(data, particelle, campiSuParticelle, campiSuGriglia, outputDATA);
  else if (data.getEvoType() == 6)  evolveLPF_nogrid(data, particelle, campiSuParticelle, outputDATA);
  else if (data.getEvoType() == 7)  evolveLPF_withgrid_onthefly(data, particelle, campiSuParticelle, outputDATA);
  else if (data.getEvoType() == 8)  evolveLPF_withgrid(data, particelle, campiSuParticelle, campiSuGriglia, outputDATA);
  else if (data.getEvoType() == 9)  evolveLPF4_nogrid(data, particelle, outputDATA);
  else if (data.getEvoType() == 10) evolveLPF4_withgrid_onthefly(data, particelle, campiSuParticelle, outputDATA);
  else if (data.getEvoType() == 11) evolveLPF4_withgrid(data, particelle, campiSuParticelle, campiSuGriglia, outputDATA);
  else if (data.getEvoType() == 12) evolveLPF4_nogrid_TURCHETTI(data, particelle, campiSuParticelle, outputDATA);

  cout << "Esecuzione completata!" << endl;


  /*************************************************************
  * Chiusura dei file di log (log.txt) e dei dati (output.dat) *
  *************************************************************/
  outputLOG.close();
  outputDATA.close();
  return 0;
}
