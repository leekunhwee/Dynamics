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

// Simulation of a double pendulum with relative coordinates

#define EASYDYNMBSMAIN // to declare the global variables
#include <EasyDyn/mbs.h>

// Definition of global application variables
double l1=1.2, l2=1.1;

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

void SetInertiaData()

{

// Inertia for body 0
body[0].mass=1.1;
body[0].PhiG.put(1,1,body[0].mass*l1*l1/12);
// Inertia for body 1
body[1].mass=0.9;
body[1].PhiG.put(1,1,body[1].mass*l2*l2/12);

}

//--------------------------------------------------------------------

void ComputeMotion()

{

// Kinematics for body 0
body[0].T0G=Trotz(q[0])*Tdisp(0,-0.5*l1,0);
body[0].omega.put(0,0,qd[0]);
body[0].omegad.put(0,0,qdd[0]);
vec OG1=body[0].T0G.e;
body[0].vG=(body[0].omega^OG1);
body[0].aG=(body[0].omegad^OG1)+(body[0].omega^(body[0].omega^OG1));
// Kinematics for body 1
body[1].TrefG=Tdisp(0,-0.5*l1,0)*Trotz(q[1])*Tdisp(0,-0.5*l2,0);
body[1].omegarel=vcoord(0,0,qd[1]);
body[1].omegadrel=vcoord(0,0,qdd[1]);
vec AG2=body[1].TrefG.R*vcoord(0,-0.5*l2,0);
body[1].vGrel=(body[1].omegarel^AG2);
body[1].aGrel=(body[1].omegadrel^AG2)
          +(body[1].omegarel^(body[1].omegarel^AG2));
ComposeMotion(1,0);
}

//--------------------------------------------------------------------

void AddAppliedEfforts()

{
// Contribution of external applied forces
vec gravity(0,-9.81,0);
AddGravityForces(gravity);
}

//--------------------------------------------------------------------

void ComputeResidual()

{
ComputeResidualmbs();
}

//--------------------------------------------------------------------

int main()

{
  // Initialisation and memory allocation
  nbrdof=2;
  nbrbody=2;
  application=new char [20]; strcpy(application, "dp4");
  InitEasyDynmbs();
  // Initial configuration
  q[1]=1; 
  // Let's go !
  NewmarkIntegration(5,0.01,0.005);
  // The clean way to finish !
  EndEasyDynmbs();
}

//--------------------------------------------------------------------
