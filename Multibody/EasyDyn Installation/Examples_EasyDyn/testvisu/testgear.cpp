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

// This file shows how to use the visu library to create
// animation files (.vol and .van)

#include <EasyDyn/vec.h>
#include <EasyDyn/visu.h>

#include <fstream>
using namespace std;

main()

{

// Here are the variables of the gear
int nd1=40, nd2=30;
double module=0.0254, R1=0.5*nd1*module, R2=0.5*nd2*module;
	
// Let's create the homogeneous transformations matrices the shapes 
// will be attached to
mth T1,T2=Tdisp(R1+R2,0,0),T3=Tdisp(-R1-R2,0,0),T4,Tref;
	
// Let's declare the shape pointers and the scene
shape *s1,*s2,*s3,*s4;
scene thescene;
	
// Let's create the shapes
vec e1(0,0,-4*module), e2(0,0,4*module);
s1=new gear(&T1,e1,e2,module,nd1,0,0.5*R1,1,2);
s2=new gear(&T2,e1,e2,module,nd2,3.1416/nd2,0.5*R2,1,4);
s3=new gear(&T3,e1,e2,module,nd2,3.1416/nd2,0.5*R2,1,4);
s4=new gear(&T4,e1,e2,module,nd1+2*nd2,3.1416/(nd1+2*nd2),1.2*(R1+2*R2),1,6);
	
// Let's add the shapes to the scene
thescene.AddShape(s1);
thescene.AddShape(s2);
thescene.AddShape(s3);
thescene.AddShape(s4);

// uncomment the following line if you want a moving observer
// thescene.SetVisuFrame(&Tref);

// Let's save the structure of the scene
thescene.CreateVolFile("testgear.vol");

// Let's create an animation by rotation of the pieces
ofstream VanFile("testgear.van");
int i;
for (i=0;i<10;i++)
     {
     T1=Trotz(i*0.1*6.2832/nd1);
     T2=Tdisp(R1+R2,0,0)*Trotz(-(i*0.1/nd1)*6.2832*R1/R2);
     T3=Tdisp(-R1-R2,0,0)*Trotz(-(i*0.1/nd1)*6.2832*R1/R2);
     T4=Trotz(-(i*0.1/nd1)*6.2832*R1/(R1+2*R2));
     Tref=Tdisp(0.5,0.5,0)*Troty(i*0.01*6.2832);
     thescene.WriteCoord(VanFile);
     }
VanFile.close();

cout << "Files testgear.vol and testgear.van created !" << endl;
	
}
