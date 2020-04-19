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
	
// Simulation of a double pendulum

#define EASYDYNSIMMAIN // to declare the global variables
#include <EasyDyn/sim.h>

double m1=1.1,m2=0.9,l1=1.2,l2=1.1,IG1zz=m1*l1*l1/12,IG2zz=m2*l2*l2/12,g=9.81;

//--------------------------------------------------------------------

void ComputeResidual()
{
f[0]=(m1*l1*l1/4+m2*l1*l1+IG1zz+m2*l1*l2*0.5*cos(q[1]))*qdd[0]
      +m2*l1*l2*0.5*cos(q[1])*qdd[1]
      -m2*l1*l2*0.5*sin(q[1])*(qd[0]+qd[1])*(qd[0]+qd[1])
      +(0.5*m1+m2)*g*l1*sin(q[0]);
f[1]=(m2*l2*l2/4+IG2zz+m2*l1*l2*0.5*cos(q[1]))*qdd[0]
     +(m2*l2*l2/4+IG2zz)*qdd[1]
     +m2*l1*l2*0.5*sin(q[1])*qd[0]*qd[0]
     +m2*g*l2*0.5*sin(q[0]+q[1]);
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
  nbrdof=2; 
  application=new char[20]; strcpy(application, "dp1");
  InitEasyDynsim();
  // Initial configuration
  q[1]=1;
  // Let's go !
  NewmarkIntegration(5,0.01,0.005);
  EndEasyDynsim();
}

//--------------------------------------------------------------------
