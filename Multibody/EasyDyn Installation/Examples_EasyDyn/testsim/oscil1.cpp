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

double w0=2*3.14159265;

//--------------------------------------------------------------------

void ComputeResidual()
{
f[0]=qdd[0]+w0*w0*q[0];
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
  nbrdof=1;
  application=new char[20]; strcpy(application, "oscil1");
  InitEasyDynsim();
  // Initial configuration
  q[0]=1;
  // Let's go !
  NewmarkIntegration(5,0.01,0.005);
  EndEasyDynsim();
}

//--------------------------------------------------------------------
