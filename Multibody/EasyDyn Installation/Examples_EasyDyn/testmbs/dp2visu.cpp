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

// Simulation of a double pendulum with relative coordinates

#define EASYDYNMBSMAIN // to declare the global variables
#include <EasyDyn/mbs.h>
#include <EasyDyn/visu.h>

#include <math.h>
#include <string.h>
#include <fstream>
using namespace std;

scene thescene;
ofstream VanFile;

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
thescene.WriteCoord(VanFile);
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
body[1].T0G=Trotz(q[0])*Tdisp(0,-l1,0)*Trotz(q[1])*Tdisp(0,-0.5*l2,0);
body[1].omega=body[0].omega+vcoord(0,0,qd[1]);
body[1].omegad=body[0].omegad+vcoord(0,0,qdd[1]);
vec AG2=body[1].T0G.R*vcoord(0,-0.5*l2,0);
body[1].vG=2*body[0].vG+(body[1].omega^AG2);
body[1].aG=2*body[0].aG+(body[1].omegad^AG2)
          +(body[1].omega^(body[1].omega^AG2));

}

//--------------------------------------------------------------------

void AddAppliedEfforts()

{
// Contribution of applied external forces
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
  // Initialization and memory allocation
  nbrdof=2;
  nbrbody=2;
  application=new char[40];
  strcpy(application,"dp2visu");
  InitEasyDynmbs();
  // Initial configuration
  q[1]=1; 

  // Let's create the shapes
  shape *s1,*s2;
  vec e1(-0.1,-0.6,-0.1), e2(0.1,0.6,0.1);
  s1=new box(&(body[0].T0G),e1,e2,1,2);
  vec e3(-0.1,-0.55,-0.1),e4(0.1,0.55,0.1);
  s2=new box(&(body[1].T0G),e3,e4,1,2);

  // Let's add the shapes to the scene
  thescene.AddShape(s1);
  thescene.AddShape(s2);

  // uncomment the following line if you want a moving observer
  // thescene.SetVisuFrame(&Tref);

  // Let's save the structure of the scene
  ComputeMotion(); // to get relevant position matrices
  thescene.CreateVolFile("dp2visu.vol");

  // Let's open the animation file
  VanFile.open("dp2visu.van");

  // Let's go !
  NewmarkIntegration(5,0.01,0.005);
  // The clean way to finish !
  EndEasyDynmbs();
  VanFile.close();
}

//--------------------------------------------------------------------
