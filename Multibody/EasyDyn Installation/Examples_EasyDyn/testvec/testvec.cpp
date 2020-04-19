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

#include <stdio.h>
#include <math.h>
#include <EasyDyn/vec.h>

using namespace std;

int main()

{

/*************************************************************************/

cout<<"*******************************************************************\n";
cout<<"*           Validation of vector operators and methods            *\n";
cout<<"*******************************************************************\n";

// Test of assignments
vec a(1,2,3),b,c; 
b.put(0.3,0.1,0.2); 

// Test of output
cout << "Vector a= " << a << "\n";
cout << "Vector b= " << b << "\n";

// Test of method length()
cout << "Length of a= " << a.length() << "\n";
cout << "Length of a= " << sqrt(a*a) << "\n";

// Test of method unite()
cout << "a normalized= " << a/a.length() << "\n";
 c=a; c.unite();
cout << "a normalized= " << c << "\n";

// Test of left and right multiplication by a double
cout << "2*a=" << 2*a << "\n";
cout << "a*2=" << a*2 << "\n";

// Test of division by a double
cout << "a/2=" << a/2 << "\n";

// Test of operators +,-,*,^,+= and -= 
cout << "a+b=" << a+b << "\n";
cout << "a+b=" << vcoord(1,2,3)+vcoord(0.3,0.1,0.2) << "\n";
cout << "a-b=" << a-b << "\n";
cout << "a-b=" << vcoord(1,2,3)-vcoord(0.3,0.1,0.2) << "\n";
cout << "a*b=" << a*b << "\n";
cout << "a*b=" << vcoord(1,2,3)*vcoord(0.3,0.1,0.2) << "\n";
cout << "a^b=" << (a^b) << "\n";
cout << "a^b=" << (vcoord(1,2,3)^vcoord(0.3,0.1,0.2)) << "\n";
cout << "a*(a^b)=" << a*(a^b) << "\n";
cout << "b*(a^b)=" << b*(a^b) << "\n";
c=a; c+=a;
cout << "a+=a" << c << "\n";
c=a; c-=c;
cout << "a-=a" << c << "\n";

// Test of long expressions
cout << "The next two vectors should be the same\n";
cout << "a^(b^(b^a))=" << (a^(b^(b^a))) << "\n";
cout << "(a*b)*(a^b)=" << (a*b)*(a^b) << "\n";

/*************************************************************************/

cout<<"*******************************************************************\n";
cout<<"*       Test of rotation tensor utilities and multiplication      *\n";
cout<<"*******************************************************************\n";
 
trot R1,R2,R3;
R1=Rrotz(0.1)*Rroty(-0.2)*Rrotx(0.3);
R2=Rrotx(-0.3)*Rroty(0.2)*Rrotz(-0.1);
cout << "Here is a rotation tensor\n" << R2;
cout << "We should get the same rotation tensor\n" << R1.inv();
cout << "We should get the identity matrix\n" << R1*R2;
R3=R1;
R3*=R2;
cout << "We should get the identity matrix again\n" << R3;

// Test of multiplication of vectors by rotation tensors
c=R1*a;
cout << "The next two vectors should be the same\n";
cout << a << "\n" << R2*c << "\n";
cout << "The next two values should be the same\n";
cout << a.length() << "  " << c.length() << "\n";

// Test of the decomposition routines
double th1,th2,th3;
vec axe1,axe2,axe3;
R3.unite();
cout<<"We should get 6 times 0.1 -0.2 0.3 followed by twice the same vector\n";
getrotzyx(R3,R1,th1,th2,th3,axe1,axe2,axe3);
cout << th1 << " " << th2  << " " << th3 << endl;
cout << axe1 << endl << R3.uz() << endl;
cout << axe3 << endl << R1.ux() << endl;
R3=R1;
R1=R3*Rroty(0.1)*Rrotx(-0.2)*Rrotz(0.3);
getrotyxz(R3,R1,th1,th2,th3,axe1,axe2,axe3);
cout << th1 << " " << th2  << " " << th3 << endl;
cout << axe1 << endl << R3.uy() << endl;
cout << axe3 << endl << R1.uz() << endl;
R3=R1;
R1=R3*Rrotx(0.1)*Rrotz(-0.2)*Rroty(0.3);
getrotxzy(R3,R1,th1,th2,th3,axe1,axe2,axe3);
cout << th1 << " " << th2  << " " << th3 << endl;
cout << axe1 << endl << R3.ux() << endl;
cout << axe3 << endl << R1.uy() << endl;
R3=R1;
R1=R3*Rrotx(0.1)*Rroty(-0.2)*Rrotz(0.3);
getrotxyz(R3,R1,th1,th2,th3,axe1,axe2,axe3);
cout << th1 << " " << th2  << " " << th3 << endl;
cout << axe1 << endl << R3.ux() << endl;
cout << axe3 << endl << R1.uz() << endl;
R3=R1;
R1=R3*Rrotz(0.1)*Rrotx(-0.2)*Rroty(0.3);
getrotzxy(R3,R1,th1,th2,th3,axe1,axe2,axe3);
cout << th1 << " " << th2  << " " << th3 << endl;
cout << axe1 << endl << R3.uz() << endl;
cout << axe3 << endl << R1.uy() << endl;
R3=R1;
R1=R3*Rroty(0.1)*Rrotz(-0.2)*Rrotx(0.3);
getrotyzx(R3,R1,th1,th2,th3,axe1,axe2,axe3);
cout << th1 << " " << th2  << " " << th3 << endl;
cout << axe1 << endl << R3.uy() << endl;
cout << axe3 << endl << R1.ux() << endl;

/*************************************************************************/

cout<<"*******************************************************************\n";
cout<<"*   Validation of transformation matrices methods and operators   *\n"; 
cout<<"*******************************************************************\n";

mth T1,T2,T3;
T1=Trotz(0.1)*Tdisp(a)*Troty(-0.2)*Tdisp(b)*Trotx(0.3)*Tdisp(-0.2,-0.4,-0.7);
T2=Tdisp(0.2,0.4,0.7)*Trotx(-0.3)*Tdisp(-b)*Troty(0.2)*Tdisp(-a)*Trotz(-0.1);
cout << "Here is an homogeneous transformation matrix\n" << T2;
cout << "We should get the same matrix\n" << T1.inv();
cout << "We should get the identity matrix\n" << T1*T2;
T3=T1;
T3*=T2;
cout << "We should get the identity matrix again\n" << T3;

// Test of multiplication of vectors by transformation matrices
c=T1*a;
cout << "The next two vectors should be the same\n";
cout << a << "\n" << T2*c << "\n";

// Test of the decomposition routines
T3.unite();
cout<<"We should get 6 times 0.1 -0.2 0.3 followed by twice the same vector\n";
getrotzyx(T3.R,T1.R,th1,th2,th3,axe1,axe2,axe3);
cout << th1 << " " << th2  << " " << th3 << endl;
cout << axe1 << endl << T3.R.uz() << endl;
cout << axe3 << endl << T1.R.ux() << endl;
T3=T1;
T1=T3*Troty(0.1)*Trotx(-0.2)*Trotz(0.3);
getrotyxz(T3.R,T1.R,th1,th2,th3,axe1,axe2,axe3);
cout << th1 << " " << th2  << " " << th3 << endl;
cout << axe1 << endl << T3.R.uy() << endl;
cout << axe3 << endl << T1.R.uz() << endl;
T3=T1;
T1=T3*Trotx(0.1)*Trotz(-0.2)*Troty(0.3);
getrotxzy(T3.R,T1.R,th1,th2,th3,axe1,axe2,axe3);
cout << th1 << " " << th2  << " " << th3 << endl;
cout << axe1 << endl << T3.R.ux() << endl;
cout << axe3 << endl << T1.R.uy() << endl;
T3=T1;
T1=T3*Trotx(0.1)*Troty(-0.2)*Trotz(0.3);
getrotxyz(T3.R,T1.R,th1,th2,th3,axe1,axe2,axe3);
cout << th1 << " " << th2  << " " << th3 << endl;
cout << axe1 << endl << T3.R.ux() << endl;
cout << axe3 << endl << T1.R.uz() << endl;
T3=T1;
T1=T3*Trotz(0.1)*Trotx(-0.2)*Troty(0.3);
getrotzxy(T3.R,T1.R,th1,th2,th3,axe1,axe2,axe3);
cout << th1 << " " << th2  << " " << th3 << endl;
cout << axe1 << endl << T3.R.uz() << endl;
cout << axe3 << endl << T1.R.uy() << endl;
T3=T1;
T1=T3*Troty(0.1)*Trotz(-0.2)*Trotx(0.3);
getrotyzx(T3.R,T1.R,th1,th2,th3,axe1,axe2,axe3);
cout << th1 << " " << th2  << " " << th3 << endl;
cout << axe1 << endl << T3.R.uy() << endl;
cout << axe3 << endl << T1.R.ux() << endl;

// Test of routine decomposevec
double a1,a2,a3;
if (decomposevec(a,axe1,axe2,axe3,a1,a2,a3))
  {
    cout << "Vecteurs de la base coplanaires !!" << endl;
    exit(1);
  }
cout << "The next two vectors should be the same\n";
cout << a << endl << a1*axe1+a2*axe2+a3*axe3 << endl;

cout<<"*******************************************************************\n";
cout<<"*              Validation of inertia tensor routines              *\n";
cout<<"*******************************************************************\n";

// Test of initialisation

tiner Phi1(1,2,3),Phi2;
Phi2.put(1,2,3);

cout << "We should get '1 2 3 0 0 0' twice\n";
cout << Phi1 << "\n";
cout << Phi2 << "\n";

Phi1.put(1,2,3,0.1,0.2,0.3);
cout << "We should get '1 2 3 0.1 0.2 0.3'\n";
cout << Phi1 << "\n";
	
cout << "We should get '0.2 3 8.2'\n";
cout << "Phi1*a=" << Phi1*a << "\n";
	
}
