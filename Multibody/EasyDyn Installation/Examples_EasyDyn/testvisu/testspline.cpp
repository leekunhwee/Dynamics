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

// Let's create the homogeneous transformation matrices the shapes 
// will be attached to
mth T1;
	
// Let's declare the shape pointers and the scene
shape *s1;
scene thescene;

// Let's create the shapes
vec e[2], t[2];
e[0].put(1,0,0); e[1].put(0,1,0);
t[0].put(0,1.567,0); t[1].put(-1.567,0,0);
s1=new grspline(&T1,e,t,1,10,3);

// Let's add the shapes to the scene
thescene.AddShape(s1);

// Let's save the structure of the scene
thescene.CreateVolFile("testspline.vol");

cout << "File testspline.vol created !" << endl;
	
}
