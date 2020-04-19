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

// Simulation of a tyre falling on the ground
// The wheel has an initial rotation velocity so that
// it should finally roll at a constant speed
// an initial slip or camber angle can be given to test the
// lateral behaviour too

#define EASYDYNMBSMAIN // to declare the global variables
#include <EasyDyn/mbs.h>
#include <EasyDyn/visu.h>

#include <fstream>
using namespace std;

scene thescene;
ofstream VanFile;

// Definition of global application variables
double R1=0.5, R2=0.1, camberini=0.3, slipini=0.1;
structtyre tyre1;

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
body[0].PhiG.put(1,1,body[0].mass*R1*R1/2);

}

//--------------------------------------------------------------------

void ComputeMotion()

{
// Kinematics for body 0
body[0].T0G=Tdisp(q[0],q[1],q[2])*Trotz(slipini)*Trotx(camberini)*Troty(q[3]);
body[0].omega=qd[3]*body[0].T0G.R.uy();
body[0].omegad=qdd[3]*body[0].T0G.R.uy();
body[0].vG.put(qd[0],qd[1],qd[2]);
body[0].aG.put(qdd[0],qdd[1],qdd[2]);
}

//--------------------------------------------------------------------

void AddAppliedEfforts()

{
// Contribution of applied external forces
vec gravity(0,0,-9.81);
AddGravityForces(gravity);
AddTyreEfforts(0,vcoord(0,1,0),tyre1);
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
  nbrdof=4;
  nbrbody=1;
  application=new char[20]; strcpy(application, "tyre2");
  InitEasyDynmbs();
  // Initial configuration
  tyre1.r1=R1;
  tyre1.r2=R2;
  tyre1.Kz=10000;
  tyre1.Cz=10;
  tyre1.Fznom=2500;
  tyre1.Clongnom=20000;
  tyre1.nlong=0.8;
  tyre1.Clatnom=20000;
  tyre1.nlat=0.8;
  tyre1.Ccambernom=2000;
  tyre1.ncamber=0.8;
  tyre1.fClbs=1.0;
  tyre1.fClbd=0.8;
	
  q[2]=0.5; 
  qd[3]=1;

  // Let's create the shapes
  shape *s1,*s2;
  vec e1(0,-R2,0), e2(0,R2,0);
  s1=new tyre(&(body[0].T0G),e1,e2,R1-2*R2,R1,12,7,0);
  vec e3(-1,-1,0),e4(1,1,-0.05);
  mth T0;
  s2=new box(&T0,e3,e4,6,6);

  // Let's add the shapes to the scene
  thescene.AddShape(s1);
  thescene.AddShape(s2);

  // uncomment the following line if you want a moving observer
  // thescene.SetVisuFrame(&Tref);

  // Let's save the structure of the scene
  ComputeMotion(); // to get relevant position matrices
  thescene.CreateVolFile("tyre2.vol");

  // Let's open the animation file
  VanFile.open("tyre2.van");

  // Let's go !
  NewmarkIntegration(2,0.01,0.005);
  // The clean way to finish !
  EndEasyDynmbs();
  VanFile.close();
}

//--------------------------------------------------------------------
