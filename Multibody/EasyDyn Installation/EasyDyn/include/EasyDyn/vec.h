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

/****************************************************************************/
/*                                                                          */
/*                             EasyDyn: vec.h                               */
/*                                                                          */
/*         Definition of various classes related to vector calculus         */
/*                                                                          */
/*      clas mth : homogeneous transformation matrix                        */
/*      clas vec : vector 3D                                                */
/*                                                                          */
/*      the ususal vector operators are also redefined :                    */
/*         matrix product, product of a vector by a matrix,                 */
/*         dot product, cross produit, addition of vectors                  */
/*                                                                          */
/****************************************************************************/

#ifndef EASYDYNVEC_H
#define EASYDYNVEC_H

#include <iostream>
#include <iomanip>
#include <math.h>
#include <stdlib.h>

using namespace std;

/****************************************************************************/

/****************************************************************************/
/*                                                                          */
/*                Definition of classes of the library                      */
/*                                                                          */
/****************************************************************************/

// Definition of class vec (3D vector)
class vec
  {
  public :
    double x,y,z;
    vec(double xinit=0.0, double yinit=0.0, double zinit=0.0);
    void put(double xcoord, double ycoord , double zcoord);
    void zero();
    void unite();
    double length();
    vec operator+(vec v2);
    vec operator+=(vec v2);
    vec operator-(vec v2);
    vec operator-=(vec v2);
    vec operator^(vec v2);
    double  operator*(vec v2);
    vec operator*(double d);
    vec operator*=(double d);
    vec operator/(double d);
    vec operator/=(double d);
    friend vec operator*(double d, vec v);
    friend vec operator-(vec v);
    friend ostream &operator<<(ostream &stream, vec v);
    friend istream &operator>>(istream &stream, vec &v);
  };

// Definition of class trot (rotation tensor)
class trot
  {
  public :
    double r11,r12,r13,r21,r22,r23,r31,r32,r33;
    trot();
    void unite();
    void rotx(double th);
    void roty(double th);
    void rotz(double th);
    void rotn(double nx, double ny, double nz, double th);
    vec ux();
    vec uy();
    vec uz();
    void setux(vec v);
    void setuy(vec v);
    void setuz(vec v);
    trot inv();
    trot operator*(trot tr2);
    trot operator*=(trot tr2);
    vec operator*(vec v);
    friend ostream &operator<<(ostream &stream, trot tr);
    friend istream &operator>>(istream &stream, trot &tr);
  };

// Definition of class tgen (general tensor)
class tgen
  {
  public :
    double Txx,Txy,Txz,Tyx,Tyy,Tyz,Tzx,Tzy,Tzz;
    tgen();
    tgen(double Jxx,double Jxy,double Jxz,
           double Jyx, double Jyy, double Jyz,
           double Jzx, double Jzy, double Jzz);
    void put(double Jxx,double Jxy,double Jxz,
        double Jyx, double Jyy, double Jyz,
        double Jzx, double Jzy, double Jzz);
    tgen trans();
    double trace();
    tgen rotate(trot R);
    tgen DiagRotate(trot R);
    tgen BeamTRrotate(trot R);
    friend ostream &operator<<(ostream &stream, tgen ti);
    friend istream &operator>>(istream &stream, tgen &ti);
    vec operator*(vec v);
    tgen operator+(tgen t);
    tgen operator*(double d);
  };

// Definition of class tiner (inertia tensor, symmetric)
class tiner
  {
  public :
    double Ixx,Iyy,Izz,Ixy,Ixz,Iyz;
    tiner();
    tiner(double Jxx,double Jyy,double Jzz,
          double Jxy=0, double Jxz=0, double Jyz=0);
    void put(double Jxx, double Jyy, double Jzz,
             double Jxy=0, double Jxz=0, double Jyz=0);
    double trace();
    tiner rotate(trot R);
    friend ostream &operator<<(ostream &stream, tiner ti);
    friend istream &operator>>(istream &stream, tiner &ti);
    vec operator*(vec v);
  };

// Definition of class tdiag (diagonal tensor)
class tdiag
  {
  public :
    double Txx,Tyy,Tzz;
    tdiag();
    tdiag(double Jxx,double Jyy,double Jzz);
    void put(double Jxx, double Jyy, double Jzz);
    double trace();
    tgen rotate(trot R);
    friend ostream &operator<<(ostream &stream, tdiag ti);
    friend istream &operator>>(istream &stream, tdiag &ti);
    vec operator*(vec v);
  };

// Definition of class mth (homogeneous transformation matrix)
class mth
  {
  public :
    trot R; vec e;
    mth();
    void unite();
    mth inv();
    mth operator*(mth m2);
    mth operator*=(mth m2);
    vec operator*(vec v);
    friend ostream &operator<<(ostream &stream, mth m);
    friend istream &operator>>(istream &stream, mth &m);
  };

/****************************************************************************/
/*                                                                          */
/*                Definition of some utility functions                      */
/*                                                                          */
/****************************************************************************/

vec vcoord(double x, double y, double z);

trot Rrotx(double th);
trot Rroty(double th);
trot Rrotz(double th);
trot Rrotn(vec v, double th);
trot Rrotn(double nx, double ny, double nz, double th);

mth Trotx(double th);
mth Troty(double th);
mth Trotz(double th);
mth Trotn(vec v, double th);
mth Trotn(double nx,double ny, double nz, double th);


mth Tdisp(vec v);
mth Tdisp(double x, double y, double z);

void getrotxyz(trot m1, trot m2, double & th1, double & th2, double & th3,
                            vec &axe1, vec &axe2, vec &axe3);
void getrotyzx(trot m1, trot m2, double & th1, double & th2, double & th3,
                            vec &axe1, vec &axe2, vec &axe3);
void getrotzxy(trot m1, trot m2, double & th1, double & th2, double & th3,
                            vec &axe1, vec &axe2, vec &axe3);
void getrotyxz(trot m1, trot m2, double & th1, double & th2, double & th3,
                            vec &axe1, vec &axe2, vec &axe3);
void getrotxzy(trot m1, trot m2, double & th1, double & th2, double & th3,
                            vec &axe1, vec &axe2, vec &axe3);
void getrotzyx(trot m1, trot m2, double & th1, double & th2, double & th3,
                            vec &axe1, vec &axe2, vec &axe3);
int decomposevec(vec a,vec axe1, vec axe2, vec axe3, 
		 double &a1,double &a2,double &a3);

#endif
