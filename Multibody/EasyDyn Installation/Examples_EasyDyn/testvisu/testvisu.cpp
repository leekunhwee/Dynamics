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

#include <math.h>
#include <fstream>

#include <EasyDyn/vec.h>
#include <EasyDyn/visu.h>

using namespace std;

main()

{
	
// Let's create the homogeneous transformations matrices the shapes 
// will be attached to
  mth T1,T2=Tdisp(5,0,0),T3=Tdisp(0,0,-5),T4=Tdisp(0,0,-5),
      T5=Tdisp(0,5,0), T6=Tdisp(0,5,-2.5),Tref=Tdisp(0.5,0.5,0);

// Let's declare the shpe pointers and the scene
shape *s1,*s2,*s3,*s4,*s5,*s6,*s7;
scene thescene;

// Let's create the shapes
vec e1(-1,-1,-1), e2(1,1,1);
s1=new box(&T1,e1,e2,1,2);
s2=new line(&T1,2*e1,2*e2,1);
vec e3(0,1,0),e4(0,-1,0);
s3=new frustum(&T2,e3,e4,0.7,1.4,12,3,4);
vec e5(1,0,0),e6(-1,0,2);
s4=new triangle(&T3,e3,e4,e5,5,6);
s5=new triangle2(&T3,e3,&T3,e4,&T4,e6,5,6);
vec e7(0,0,3), e8(0,0,1), e9(0,0,0);
s6=new spring(&T5,e7,&T6,e8,0.4,0.5,10,12,1);
s7=new sphere(&T6,e9,1.0,7,12,14,14);

// Let's add the shapes to the scene
thescene.AddShape(s1);
thescene.AddShape(s2);
thescene.AddShape(s3);
thescene.AddShape(s4);
thescene.AddShape(s5);
thescene.AddShape(s6);
thescene.AddShape(s7);

// uncomment the following line if you want a moving observer
// thescene.SetVisuFrame(&Tref);

// Let's save the structure of the scene
thescene.CreateVolFile("testvisu.vol");

// Let's create an animation by rotation of the pieces
ofstream VanFile("testvisu.van");
int i;
for (i=0;i<100;i++)
     {
     T1=Trotz(i*0.02*6.2832);
     T2=Tdisp(5,0,0)*Troty(i*0.03*6.2832);
     T3=Tdisp(0,0,-5)*Trotx(i*0.02*6.2832);
     T4=Tdisp(0,sin(0.02*i*3.14159265),-5);
     T6=Tdisp(0,5,-2.25+sin(0.02*i*3.14159265));
     Tref=Tdisp(0.5,0.5,0)*Trotz(i*0.01*6.2832);
     thescene.WriteCoord(VanFile);
     }
VanFile.close();

cout << "Files testvisu.vol and testvisu.van created !" << endl;
	
}
