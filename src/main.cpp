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
#include <vector>

#include "datatypes.h" // include il file con la definizione dei tipi degli oggetti: Particle, Field, Data
#include "evolution.h" // include i metodi di evoluzione del sistema
#include "filler.h" // include i metodi di riempimento dei vettori di particelle e dei campi, iniziali e durante l'elaborazione

using namespace std;

int main(int argc, char* argv[])
{
    /************************************************************
  * INIZIALIZZAZIONE generatore random tramite seed temporale *
  ************************************************************/
    time_t now;
    time(&now);
    srand((unsigned int)now);

    /*************************************************************
  * INIZIALIZZAZIONE file di output                            *
  * - outputLOG conterra` i parametri dell'esecuzione corrente *
  *   ed apre in modalita` append il file log.txt              *
  * - outputDATA conterra` invece tutti i dati delle           *
  *   particelle istante per istante e va a sovrascrivere ad   *
  *   ogni esecuzione del programma il file output.dat         *
  *************************************************************/
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
    vector<Field> campiSuGriglia;

    /************************************************************
  * L'eseguibile legge l'input da un file json contenente     *
  * tutti i parametri necessari per l'esecuzione.             *
  * La sua assenza sulla riga di comando e` bloccante         *
  ************************************************************/
    if (argc < 2) {
        cout << "E` necessario specificare un file json contenente la configurazione della simulazione!" << endl;
        exit(-1);
    }

    data.fillParametersFromFile(argv[1], outputLOG);

    /*************************************************************
  * Ora il programma dovrebbe aver saputo "tutto" dall'utente  *
  * e pertanto genera il vettore iniziale delle particelle     *
  * seguendo i parametri impostati e solo nel caso non fosse   *
  * stato gia` definito e letto un file contenente questi dati *
  *************************************************************/

    creaVettoreParticelle(data, particelle);

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

    cout << "NB: la griglia costruita non viene utilizzata,\nserve solo per inserire inutilmente le particelle in una cella" << endl;

    if (data.getEvoType() == 1)
        calcolaCampoInUnIstante(data, particelle, campiSuParticelle, outputDATA);
    else if (data.getEvoType() == 2)
        calcolaEvoluzioneCampoSuParticelle(data, particelle, campiSuParticelle, outputDATA);
    else if (data.getEvoType() == 3)
        evolveRK4_nogrid(data, particelle, campiSuParticelle, outputDATA);
    else if (data.getEvoType() == 4)
        evolveRK4_withgrid_onthefly(data, particelle, campiSuParticelle, outputDATA);
    else if (data.getEvoType() == 5)
        evolveRK4_withgrid(data, particelle, campiSuParticelle, campiSuGriglia, outputDATA);
    else if (data.getEvoType() == 6)
        evolveLPF_nogrid(data, particelle, campiSuParticelle, outputDATA);
    else if (data.getEvoType() == 7)
        evolveLPF_withgrid_onthefly(data, particelle, campiSuParticelle, outputDATA);
    else if (data.getEvoType() == 8)
        evolveLPF_withgrid(data, particelle, campiSuParticelle, campiSuGriglia, outputDATA);
    else if (data.getEvoType() == 9)
        evolveLPF4_nogrid(data, particelle, outputDATA);
    else if (data.getEvoType() == 10)
        evolveLPF4_withgrid_onthefly(data, particelle, campiSuParticelle, outputDATA);
    else if (data.getEvoType() == 11)
        evolveLPF4_withgrid(data, particelle, campiSuParticelle, campiSuGriglia, outputDATA);
    else if (data.getEvoType() == 12)
        evolveLPF4_nogrid_TURCHETTI(data, particelle, campiSuParticelle, outputDATA);

    cout << "Esecuzione completata!" << endl;

    /*************************************************************
  * Chiusura dei file di log (log.txt) e dei dati (output.dat) *
  *************************************************************/
    outputLOG.close();
    outputDATA.close();
    return 0;
}
