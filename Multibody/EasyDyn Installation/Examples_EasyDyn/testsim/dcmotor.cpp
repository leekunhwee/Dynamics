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

// Simulation of DC motor

#include <math.h>
#include <string.h>

#define EASYDYNSIMMAIN // to declare the global variables
#include <EasyDyn/sim.h>

double L=0.5, R=2, Kb=0.02, J=0.0002, Km=0.02, Kf=0; //0.002;

//--------------------------------------------------------------------

void ComputeResidual()
{
f[0]=L*qdd[0]+R*qd[0]+Kb*qd[1]-u[0];
f[1]=J*qdd[1]+Kf*qd[1]-Km*qd[0];
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
  nbrdof=2; nbrinput=1;
  application=new char[20]; strcpy(application,"dcmotor");
  InitEasyDynsim();
  // Initial configuration
  q[0]=0;
  // Export input-output system
  SaveLinearizedSystem();
  // Let's go !
  u[0]=1; // simulation of the step response
  NewmarkIntegration(10,0.1,0.05);
  EndEasyDynsim();
}

//--------------------------------------------------------------------
