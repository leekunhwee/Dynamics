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

#include <math.h>
#include <string.h>

// Simulation of a mass-spring system (eigenfrequency=1 Hz)

#define EASYDYNSIMMAIN // to declare the global variables
#include <EasyDyn/sim.h>

double L=0.5, R=2, Kb=0.02, J=0.0002, Km=0.02, Kf=0, Kp, Ki, Kd;

//--------------------------------------------------------------------

void ComputeResidual()
{
double w=50, erreur;
u[0]=-Kp*qd[2]-Ki*q[2]-Kd*qdd[2];
f[0]=L*qdd[0]+R*qd[0]+Kb*qd[1]-u[0];
f[1]=J*qdd[1]+Kf*qd[1]-Km*qd[0];
erreur=qd[1]-w;
f[2]=qdd[2]+1000*qd[2]-1000*erreur;
}

//--------------------------------------------------------------------

void WriteDataHeader(ostream &OutFile)
{
WriteStateVariablesHeader(OutFile);
OutFile << endl;
}

//--------------------------------------------------------------------

void SaveData(ostream &OutFile)
{
SaveStateVariables(OutFile);
OutFile << endl;
}

//--------------------------------------------------------------------

int main()
{
  // Initialisation and memory allocation
  nbrdof=3; nbrinput=1;
  application=new char[20]; strcpy(application,"dcmotorPID");
  InitEasyDynsim();
  // Initial configuration
  q[0]=0;
  // Let's go !
  cout << "Simulation of a DC motor with a PID control on velocity" << endl;
  cout << "Typical values of the controller" << endl;
  cout << "P   only: Kp=0.2  Ki=0    Kd=0" << endl;
  cout << "PI  only: Kp=0.18 Ki=0.1  Kd=0" << endl;
  cout << "PID     : Kp=0.24 Ki=0.17 Kd=0.03" << endl;
  cout << "Kp="; cin >> Kp;
  cout << "Ki="; cin >> Ki;
  cout << "Kd="; cin >> Kd;
  NewmarkIntegration(5,0.01,0.05);
  EndEasyDynsim();
}

//--------------------------------------------------------------------
