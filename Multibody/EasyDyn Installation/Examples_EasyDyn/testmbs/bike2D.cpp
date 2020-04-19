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

// Simulation of a bike (in 2D)
// The bike is launched with some initial velocity

#define EASYDYNMBSMAIN // to declare the global variables
#include <EasyDyn/mbs.h>
#include <EasyDyn/visu.h>

#include <fstream>
using namespace std;

scene thescene;
ofstream VanFile;

// Definition of global application variables
double R1=0.35, R2=0.05,L=1;
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

// Inertia for body 0 (frame)
body[0].mass=80;
body[0].PhiG.put(10,body[0].mass*L*L/12,10);
// Inertia for body 1 (front wheel)
body[1].mass=1.1;
body[1].PhiG.put(1,body[1].mass*R1*R1/2,1);
// Inertia for body 2 (rear wheel
body[2].mass=1.1;
body[2].PhiG.put(1,body[2].mass*R1*R1/2,1);

}

//--------------------------------------------------------------------

void ComputeMotion()

{
// Kinematics for body 0 (frame)
body[0].T0G=Tdisp(q[0],0,q[1])*Troty(q[2]);
body[0].omega.put(0,qd[2],0);
body[0].omegad.put(0,qdd[2],0);
body[0].vG.put(qd[0],0,qd[1]);
body[0].aG.put(qdd[0],0,qdd[1]);
// Kinematics for body 1 (front wheel)
body[1].T0G=body[0].T0G*Tdisp(0.5*L,0,0)*Troty(q[3]);
body[1].omega.put(0,qd[2]+qd[3],0);
body[1].omegad.put(0,qdd[2]+qdd[3],0);
body[1].vG=body[0].vG-qd[2]*0.5*L*body[0].T0G.R.uz();
body[1].aG=body[0].aG-qdd[2]*0.5*L*body[0].T0G.R.uz()
                      -qd[2]*qd[2]*0.5*L*body[0].T0G.R.ux();
// Kinematics for body 2 (rear wheel)
body[2].T0G=body[0].T0G*Tdisp(-0.5*L,0,0)*Troty(q[4]);
body[2].omega.put(0,qd[2]+qd[4],0);
body[2].omegad.put(0,qdd[2]+qdd[4],0);
body[2].vG=body[0].vG+qd[2]*0.5*L*body[0].T0G.R.uz();
body[2].aG=body[0].aG+qdd[2]*0.5*L*body[0].T0G.R.uz()
                      +qd[2]*qd[2]*0.5*L*body[0].T0G.R.ux();
}

//--------------------------------------------------------------------

void AddAppliedEfforts()

{
// Contribution of applied external forces
vec gravity(0,0,-9.81);
AddGravityForces(gravity);
AddTyreEfforts(1,vcoord(0,1,0),tyre1);
AddTyreEfforts(2,vcoord(0,1,0),tyre1);
// To accelerate the bike
body[2].MG+=vcoord(0,1,0);
body[0].MG-=vcoord(0,1,0); // The reaction
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
  nbrdof=5;
  nbrbody=3;
  application=new char[20]; strcpy(application, "bike2D");
  InitEasyDynmbs();
  // Initial configuration
  tyre1.r1=R1;
  tyre1.r2=R2;
  tyre1.Kz=100000;
  tyre1.Cz=30;
  tyre1.Fznom=500;
  tyre1.Clongnom=4000;
  tyre1.nlong=0.8;
  tyre1.Clatnom=4000;
  tyre1.nlat=0.8;
  tyre1.Ccambernom=400;
  tyre1.ncamber=0.8;
  tyre1.fClbs=1.0;
  tyre1.fClbd=0.8;
	
  q[1]=R1+0.01;
  qd[3]=1;
  qd[4]=1;
  qd[0]=R1*qd[3]; 

  // Let's create the shapes
  shape *s1,*s2,*s3;
  mth T0;
  vec e3(-1,-1,0),e4(15,1,-0.05);
  s1=new box(&T0,e3,e4,6,6);
  vec e1(0,-R2,0), e2(0,R2,0);
  s2=new tyre(&(body[1].T0G),e1,e2,R1-2*R2,R1,12,7,0);
  s3=new tyre(&(body[2].T0G),e1,e2,R1-2*R2,R1,12,7,0);

  // Let's add the shapes to the scene
  thescene.AddShape(s1);
  thescene.AddShape(s2);
  thescene.AddShape(s3);

  // uncomment the following line if you want a moving observer
  thescene.SetVisuFrame(&body[0].T0G,1);

  // Let's save the structure of the scene
  ComputeMotion(); // to get relevant position matrices
  thescene.CreateVolFile("bike2D.vol");

  // Let's open the animation file
  VanFile.open("bike2D.van");

  // Let's go !
  NewmarkIntegration(20,0.1,0.05);
  // The clean way to finish !
  EndEasyDynmbs();
  VanFile.close();
}

//--------------------------------------------------------------------
