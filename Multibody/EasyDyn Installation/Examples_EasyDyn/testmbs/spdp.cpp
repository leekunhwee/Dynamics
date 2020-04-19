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

#include <math.h>
#include <string.h>

*/

// Simulation of a mass-spring-damper system

#define EASYDYNMBSMAIN // to declare the global variables
#include <EasyDyn/mbs.h>

// Definition of global application variables
vec ut; // direction of motion
// The result should be the same whatever ut

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

// Inertia for body 0 (=ground in this example)
body[0].mass=1;
body[0].PhiG.put(1,1,1); // no impact on the result anyway
// Inertia for body 1
body[1].mass=1;
body[1].PhiG.put(1,1,1); // no impact on the result anyway

}

//--------------------------------------------------------------------

void ComputeMotion()

{

// Kinematics for body 0
// =ground

// Kinematics for body 1
body[1].T0G=Tdisp(q[0]*ut);
body[1].vG=qd[0]*ut;
body[1].aG=qdd[0]*ut;

}

//--------------------------------------------------------------------

void AddAppliedEfforts()

{
// Contribution of external applied forces
vec gravity=39.478417*ut;
AddGravityForces(gravity);
AddSpringForce(39.478417,2,1,vcoord(0,0,0),0,2*ut);
// Uncomment the following if you want to add a damper
//AddDamperForce(1,0,2*ut,1,vcoord(0,0,0));
}

//--------------------------------------------------------------------

void ComputeResidual()

{
ComputeResidualmbs();
}

//--------------------------------------------------------------------

int main()

{
  // Initialization and memory allocation
  nbrdof=1;
  nbrbody=2;
  application=new char[20]; strcpy(application, "spdp");
  InitEasyDynmbs();
  // Initial configuration
  // Everything =0, so nothing to do :-)
  // Let's go !
  ut.put(1.1,2.2,3.3); ut.unite();
  NewmarkIntegration(5,0.01,0.005);
  // the clean way to finish
  EndEasyDynmbs();
}

//--------------------------------------------------------------------
