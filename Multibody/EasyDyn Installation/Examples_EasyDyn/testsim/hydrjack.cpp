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

// Simulation of a hydraulic jack pushing a mass-spring system

#define EASYDYNSIMMAIN // to declare the global variables
#include <math.h>
#include <EasyDyn/sim.h>

double m=10, k=1E5, Kf=2E8, rho=860, pE1=1E6, pE2=1E6, S1=1E-3, S2=5E-4,
  V01=1E-4, V02=3E-4, Cd=0.611, Se1=1E-5, Se2=1E-5;

//--------------------------------------------------------------------

void ComputeResidual()
{
double p1,p2,p1d,p2d;
p1=1E5*qd[1]; p1d=1E5*qdd[1]; // In qd, the pressure is in bar 
p2=1E5*qd[2]; p2d=1E5*qdd[2];
pE1=1E6+9E6*(100*t);
if (t>0.01) pE1=1E7;
f[0]=m*qdd[0]+k*q[0]-p1*S1+p2*S2;
f[1]=S1*qd[0]+(V01+S1*q[0])*p1d/Kf
    -Se1*Cd*sqrt(fabs(pE1-p1)/rho)*(pE1-p1)/(1+fabs(pE1-p1));
f[2]=-S2*qd[0]+(V02-S2*q[0])*p2d/Kf
    -Se2*Cd*sqrt(fabs(pE2-p2)/rho)*(pE2-p2)/(1+fabs(pE2-p2));
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
  // Initialization and memory allocation
  nbrdof=3; 
  application=new char[20]; strcpy(application, "hydrjack");
  InitEasyDynsim();
  // Initial configuration
  qd[1]=pE2/1E5;
  qd[2]=pE2/1E5;
  // Let's go !
  DEBUG=1;
  NewmarkIntegration(1,0.001,0.0005);
  EndEasyDynsim();
}

//--------------------------------------------------------------------
