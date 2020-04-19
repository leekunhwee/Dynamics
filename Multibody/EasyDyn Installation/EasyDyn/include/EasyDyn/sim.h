/*

Copyright (C) 2003 Olivier VERLINDEN
    Service de Mecanique rationnelle, Dynamique et Vibrations
    Faculte Polytechnique de Mons
    31, Bd Dolez, 7000 MONS (Belgium)
    Olivier.Verlinden@fpms.ac.be

This file is part of EasyDyn

EasyDyn is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the
Free Software Foundation; either version 2, or (at your option) any
later version.

EasyDyn is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with EasyDyn; see the file COPYING.  If not, write to the Free
Software Foundation, 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

*/

#include <fstream>
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <string.h>

using namespace std;

//---------------------------------------------------------------------------

/* Global variables */

#ifdef EASYDYNSIMMAIN
double t;
double *q=0,*qd=0,*qdd=0,*f=0,*u=0,*x,*xd,*y=0;
char *application=0;
int nbrdof=0;
int nbrflexdof=0;
int nbrinput=0;
int nbroutput=0;
int DEBUG=0;
double AbsTol=1E-6, RelTol=1E-6, hminmin=1E-8;
#else
extern double t;
extern double *q,*qd,*qdd,*f,*u,*x,*xd,*y;
extern char *application;
extern int nbrdof, nbrflexdof, nbrinput, nbroutput;
extern int DEBUG;
extern double AbsTol, RelTol,hminmin;
#endif

//---------------------------------------------------------------------------

// Procedures that the user must provide

void SaveData(ostream &OutFile);
void WriteDataHeader(ostream &OutFile);
void ComputeResidual();

//---------------------------------------------------------------------------

// Procedures that we provide

#ifdef _WIN32
int gettimeofday(struct timeval *tv, void *tz);
#endif

void InitEasyDynsim();
void EndEasyDynsim();
void StaticEquilibrium(int *doflocked=0);
int NewmarkOneStep(double h, double &errqd, int *doflocked=0, int hmodified=1);
void NewmarkPrediction(double h, double &errqd);
int NewmarkCorrection(double h, double &errqd, int *doflocked=0, int hmodified=1);
void NewmarkInterval(double tfinal, double &h, double hmax, int *doflocked=0);
void NewmarkIntervalold(double tfinal, double &h, double hmax, int *doflocked=0);
void NewmarkIntegration(double tfinal, double hsave, double hmax, 
			int *doflocked=0, int old=0);
void NewmarkIntegration2(double tfinal, double hsave, double hmax, 
			int *doflocked=0, int old=0);
void SetIntegrationPrecision(double NewAbsTol, double NewRelTol);
void SaveStateVariables(ostream &OutFile);
void WriteStateVariablesHeader(ostream &OutFile);
void SaveLinearizedSystem(int *doflocked=0);
void ComputePoles(int *doflocked=0, double freqmin=0,double freqmax=1E30);

//---------------------------------------------------------------------------
