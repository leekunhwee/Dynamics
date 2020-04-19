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
#include <EasyDyn/vec.h>

/****************************************************************************/
/*                        Methods relative to vec                           */
/****************************************************************************/


/****************************************************************************/

vec::vec(double xinit, double yinit, double zinit)
  {
  x=xinit; y=yinit; z=zinit;
  }

/****************************************************************************/

void vec::put(double xcoord, double ycoord, double zcoord)
  {
  x=xcoord; y=ycoord; z=zcoord;
  }

/****************************************************************************/

void vec::zero()
  {
  x=0; y=0; z=0;
  }

/****************************************************************************/

void vec::unite()
  {
  double norm=sqrt(x*x+y*y+z*z);
  if (norm>1E-20) { x/=norm; y/=norm; z/=norm; }
  }

/****************************************************************************/

double vec::length()
  {
  return sqrt(x*x+y*y+z*z);
  }

/****************************************************************************/

vec vec::operator+(vec v2)
  {
  vec vtemp;
  vtemp.x=x+v2.x;
  vtemp.y=y+v2.y;
  vtemp.z=z+v2.z;
  return vtemp;
  }

/****************************************************************************/

vec vec::operator+=(vec v2)
  {
  x+=v2.x;
  y+=v2.y;
  z+=v2.z;
  return *this;
  }

/****************************************************************************/

vec vec::operator-(vec v2)
  {
  vec vtemp;
  vtemp.x=x-v2.x;
  vtemp.y=y-v2.y;
  vtemp.z=z-v2.z;
  return vtemp;
  }

/****************************************************************************/

vec vec::operator-=(vec v2)
  {
  x-=v2.x;
  y-=v2.y;
  z-=v2.z;
  return *this;
  }

/****************************************************************************/

vec vec::operator^(vec v2)
  {
  vec vtemp;
  vtemp.x=y*v2.z-z*v2.y;
  vtemp.y=z*v2.x-x*v2.z;
  vtemp.z=x*v2.y-y*v2.x;
  return vtemp;
  }

/****************************************************************************/

double vec::operator*(vec v2)
  {
  return (x*v2.x+y*v2.y+z*v2.z);
  }

/****************************************************************************/

vec vec::operator*(double d)
  {
  vec vtemp;
  vtemp.x=x*d;
  vtemp.y=y*d;
  vtemp.z=z*d;
  return vtemp;
  }

/****************************************************************************/

vec vec::operator*=(double d)
  {
  x*=d;
  y*=d;
  z*=d;
  return *this;
  }

/****************************************************************************/

vec vec::operator/(double d)
  {
  vec vtemp;
  vtemp.x=x/d;
  vtemp.y=y/d;
  vtemp.z=z/d;
  return vtemp;
  }

/****************************************************************************/

vec vec::operator/=(double d)
  {
  x/=d;
  y/=d;
  z/=d;
  return *this;
  }

/****************************************************************************/

vec operator*(double d, vec v)
  {
  vec vtemp;
  vtemp.x=v.x*d;
  vtemp.y=v.y*d;
  vtemp.z=v.z*d;
  return vtemp;
  }

/****************************************************************************/

vec operator-(vec v)
  {
  vec vtemp;
  vtemp.x=-v.x;
  vtemp.y=-v.y;
  vtemp.z=-v.z;
  return vtemp;
  }

/****************************************************************************/

ostream &operator<<(ostream &stream, vec v)
  {
  stream << setw(15) << v.x << " " << setw(15) << v.y
	 << setw(15) << v.z;
  return stream;
  }

/****************************************************************************/

istream &operator>>(istream &stream, vec &v)
  {
  stream >> v.x >> v.y >> v.z;
  return stream;
  }

/****************************************************************************/
/*                       Methods relative to trot                           */
/****************************************************************************/

trot::trot()
  {
  r11=1.0; r12=0.0; r13=0.0;
  r21=0.0; r22=1.0; r23=0.0;
  r31=0.0; r32=0.0; r33=1.0;
  }

/****************************************************************************/

void trot::unite()
  {
  r11=1.0; r12=0.0; r13=0.0;
  r21=0.0; r22=1.0; r23=0.0;
  r31=0.0; r32=0.0; r33=1.0;
  }

/****************************************************************************/

void trot::rotx(double th)
  {
  r11=1.0; r12=0.0;     r13=0.0;     
  r21=0.0; r22=cos(th); r23=-sin(th);
  r31=0.0; r32=sin(th); r33=cos(th); 
  }

/****************************************************************************/

void trot::roty(double th)
  {
  r11=cos(th);  r12=0.0; r13=sin(th);
  r21=0.0;      r22=1.0; r23=0.0;    
  r31=-sin(th); r32=0.0; r33=cos(th);
  }

/****************************************************************************/

void trot::rotz(double th)
  {
  r11=cos(th); r12=-sin(th); r13=0.0;
  r21=sin(th); r22=cos(th);  r23=0.0;
  r31=0.0;     r32=0.0;      r33=1.0;
  }

/****************************************************************************/

void trot::rotn(double nx, double ny, double nz, double th)
  {
  double C=cos(th), V=1-cos(th), S=sin(th);
  r11=nx*nx*V+C; r12=nx*ny*V-nz*S; r13=nx*nz*V+ny*S;
  r21=nx*ny*V+nz*S; r22=ny*ny*V+C; r23=ny*nz*V-nx*S;
  r31=nx*nz*V-ny*S; r32=ny*nz*V+nx*S; r33=nz*nz*V+C;
  }

/****************************************************************************/

vec trot::ux()
  {
  vec vtemp;
  vtemp.x=r11; vtemp.y=r21; vtemp.z=r31;
  return vtemp;
  }

/****************************************************************************/

vec trot::uy()
  {
  vec vtemp;
  vtemp.x=r12; vtemp.y=r22; vtemp.z=r32;
  return vtemp;
  }

/****************************************************************************/

vec trot::uz()
  {
  vec vtemp;
  vtemp.x=r13; vtemp.y=r23; vtemp.z=r33;
  return vtemp;
  }

/****************************************************************************/

void trot::setux(vec v)
  {
  r11=v.x; r21=v.y; r31=v.z; 
  }

/****************************************************************************/

void trot::setuy(vec v)
  {
  r12=v.x; r22=v.y; r32=v.z; 
  }

/****************************************************************************/

void trot::setuz(vec v)
  {
  r13=v.x; r23=v.y; r33=v.z; 
  }

/****************************************************************************/

trot trot::inv()
  {
  trot trtemp;
  trtemp.r11=r11; trtemp.r12=r21; trtemp.r13=r31;
  trtemp.r21=r12; trtemp.r22=r22; trtemp.r23=r32;
  trtemp.r31=r13; trtemp.r32=r23; trtemp.r33=r33;
  return trtemp;
  }

/****************************************************************************/

trot trot::operator*(trot tr2)
  {
  trot trtemp;
  trtemp.r11=r11*tr2.r11+r12*tr2.r21+r13*tr2.r31;
  trtemp.r12=r11*tr2.r12+r12*tr2.r22+r13*tr2.r32;
  trtemp.r13=r11*tr2.r13+r12*tr2.r23+r13*tr2.r33;
  trtemp.r21=r21*tr2.r11+r22*tr2.r21+r23*tr2.r31;
  trtemp.r22=r21*tr2.r12+r22*tr2.r22+r23*tr2.r32;
  trtemp.r23=r21*tr2.r13+r22*tr2.r23+r23*tr2.r33;
  trtemp.r31=r31*tr2.r11+r32*tr2.r21+r33*tr2.r31;
  trtemp.r32=r31*tr2.r12+r32*tr2.r22+r33*tr2.r32;
  trtemp.r33=r31*tr2.r13+r32*tr2.r23+r33*tr2.r33;
  return trtemp;
  }

/****************************************************************************/

trot trot::operator*=(trot tr2)
  {
  double r1,r2,r3;
  r1=r11*tr2.r11+r12*tr2.r21+r13*tr2.r31;
  r2=r11*tr2.r12+r12*tr2.r22+r13*tr2.r32;
  r3=r11*tr2.r13+r12*tr2.r23+r13*tr2.r33;
  r11=r1; r12=r2; r13=r3;
  r1=r21*tr2.r11+r22*tr2.r21+r23*tr2.r31;
  r2=r21*tr2.r12+r22*tr2.r22+r23*tr2.r32;
  r3=r21*tr2.r13+r22*tr2.r23+r23*tr2.r33;
  r21=r1; r22=r2; r23=r3;
  r1=r31*tr2.r11+r32*tr2.r21+r33*tr2.r31;
  r2=r31*tr2.r12+r32*tr2.r22+r33*tr2.r32;
  r3=r31*tr2.r13+r32*tr2.r23+r33*tr2.r33;
  r31=r1; r32=r2; r33=r3;
  return *this;
  }

/****************************************************************************/

vec trot::operator*(vec v)
  {
  vec vtemp;
  vtemp.x=r11*v.x+r12*v.y+r13*v.z;
  vtemp.y=r21*v.x+r22*v.y+r23*v.z;
  vtemp.z=r31*v.x+r32*v.y+r33*v.z;
  return vtemp;
  }

/****************************************************************************/

ostream &operator<<(ostream &stream, trot tr)
  {
  stream << tr.r11 << " " << tr.r12 << " " << tr.r13 << "\n";
  stream << tr.r21 << " " << tr.r22 << " " << tr.r23 << "\n";
  stream << tr.r31 << " " << tr.r32 << " " << tr.r33 << "\n";
  return stream;
  }

/****************************************************************************/

istream &operator>>(istream &stream, trot &tr)
  {
  stream >> tr.r11 >> tr.r12 >> tr.r13;
  stream >> tr.r21 >> tr.r22 >> tr.r23;
  stream >> tr.r31 >> tr.r32 >> tr.r33;
  return stream;
  }

/****************************************************************************/

/****************************************************************************/
/*                      Methods relative to tgen                            */
/****************************************************************************/

tgen::tgen()
  {
    Txx=0.0; Txy=0.0; Txz=0.0;
    Tyx=0.0; Tyy=0.0; Tyz=0.0;
    Tzx=0.0; Tzy=0.0; Tzz=0.0;
  }

/****************************************************************************/

tgen::tgen(double Jxx,double Jxy,double Jxz,
           double Jyx, double Jyy, double Jyz,
           double Jzx, double Jzy, double Jzz)
  {
  Txx=Jxx; Txy=Jxy; Txz=Jxz;
  Tyx=Jyx; Tyy=Jyy; Tyz=Jyz; 
  Tzx=Jzx; Tzy=Jzy; Tzz=Jzz; 
  }

/****************************************************************************/

void tgen::put(double Jxx,double Jxy,double Jxz,
           double Jyx, double Jyy, double Jyz,
           double Jzx, double Jzy, double Jzz)
  {
  Txx=Jxx; Txy=Jxy; Txz=Jxz;
  Tyx=Jyx; Tyy=Jyy; Tyz=Jyz; 
  Tzx=Jzx; Tzy=Jzy; Tzz=Jzz; 
  }

/****************************************************************************/

vec tgen::operator*(vec v)
  {
  vec vtemp;
  vtemp.x = Txx*v.x + Txy*v.y + Txz*v.z;
  vtemp.y = Tyx*v.x + Tyy*v.y + Tyz*v.z;
  vtemp.z = Tzx*v.x + Tzy*v.y + Tzz*v.z;
  return vtemp;
  }

/****************************************************************************/

tgen tgen::operator+(tgen t)
  {
  tgen tsrtemp;
  tsrtemp.Txx=Txx+t.Txx; tsrtemp.Txy=Txy+t.Txy; tsrtemp.Txz=Txz+t.Txz;
  tsrtemp.Tyx=Tyx+t.Tyx; tsrtemp.Tyy=Tyy+t.Tyy; tsrtemp.Tyz=Tyz+t.Tyz;
  tsrtemp.Tzx=Tzx+t.Tzx; tsrtemp.Tzy=Tzy+t.Tzy; tsrtemp.Tzz=Tzz+t.Tzz;
  return tsrtemp;
  }

/****************************************************************************/

tgen tgen::operator*(double d)
  {
  tgen tsrtemp;
  tsrtemp.Txx=Txx*d; tsrtemp.Txy=Txy*d; tsrtemp.Txz=Txz*d;
  tsrtemp.Tyx=Tyx*d; tsrtemp.Tyy=Tyy*d; tsrtemp.Tyz=Tyz*d;
  tsrtemp.Tzx=Tzx*d; tsrtemp.Tzy=Tzy*d; tsrtemp.Tzz=Tzz*d;
  return tsrtemp;
  }
  
/****************************************************************************/

tgen tgen::trans()
  {
  tgen tsrtemp;
  tsrtemp.Txx=Txx; tsrtemp.Txy=Tyx; tsrtemp.Txz=Tzx;
  tsrtemp.Tyx=Txy; tsrtemp.Tyy=Tyy; tsrtemp.Tyz=Tzy;
  tsrtemp.Tzx=Txz; tsrtemp.Tzy=Tyz; tsrtemp.Tzz=Tzz;
  return tsrtemp;
  }

/****************************************************************************/

double tgen::trace()
  {
    return (Txx+Tyy+Tzz);
  }

/****************************************************************************/

tgen tgen::rotate(trot R)
  {
  tgen tsrtemp;
  tsrtemp.Txx= R.r11*(R.r11*Txx + R.r12*Tyx + R.r13*Tzx)
             + R.r12*(R.r11*Txy + R.r12*Tyy + R.r13*Tzy)
             + R.r13*(R.r11*Txz + R.r12*Tyz + R.r13*Tzz);
  tsrtemp.Txy= R.r21*(R.r11*Txx + R.r12*Tyx + R.r13*Tzx)
             + R.r22*(R.r11*Txy + R.r12*Tyy + R.r13*Tzy)
             + R.r23*(R.r11*Txz + R.r12*Tyz + R.r13*Tzz);
  tsrtemp.Txz= R.r31*(R.r11*Txx + R.r12*Tyx + R.r13*Tzx)
             + R.r32*(R.r11*Txy + R.r12*Tyy + R.r13*Tzy)
             + R.r33*(R.r11*Txz + R.r12*Tyz + R.r13*Tzz);
  tsrtemp.Tyx= R.r11*(R.r21*Txx + R.r22*Tyx + R.r23*Tzx)
             + R.r12*(R.r21*Txy + R.r22*Tyy + R.r23*Tzy)
             + R.r13*(R.r21*Txz + R.r22*Tyz + R.r23*Tzz);
  tsrtemp.Tyy= R.r21*(R.r21*Txx + R.r22*Tyx + R.r23*Tzx)
             + R.r22*(R.r21*Txy + R.r22*Tyy + R.r23*Tzy)
             + R.r23*(R.r21*Txz + R.r22*Tyz + R.r23*Tzz);
  tsrtemp.Tyz= R.r31*(R.r21*Txx + R.r22*Tyx + R.r23*Tzx)
             + R.r32*(R.r21*Txy + R.r22*Tyy + R.r23*Tzy)
             + R.r33*(R.r21*Txz + R.r22*Tyz + R.r23*Tzz);
  tsrtemp.Tzx= R.r11*(R.r31*Txx + R.r32*Tyx + R.r33*Tzx)
             + R.r12*(R.r31*Txy + R.r32*Tyy + R.r33*Tzy)
             + R.r13*(R.r31*Txz + R.r32*Tyz + R.r33*Tzz);
  tsrtemp.Tzy= R.r21*(R.r31*Txx + R.r32*Tyx + R.r33*Tzx)
             + R.r22*(R.r31*Txy + R.r32*Tyy + R.r33*Tzy)
             + R.r23*(R.r31*Txz + R.r32*Tyz + R.r33*Tzz);
  tsrtemp.Tzz= R.r31*(R.r31*Txx + R.r32*Tyx + R.r33*Tzx)
             + R.r32*(R.r31*Txy + R.r32*Tyy + R.r33*Tzy)
             + R.r33*(R.r31*Txz + R.r32*Tyz + R.r33*Tzz);
  return tsrtemp;
  }
  
/****************************************************************************/

tgen tgen::DiagRotate(trot R)
  {
  tgen tsrtemp;
  tsrtemp.Txx= R.r11*(R.r11*Txx)
             + R.r12*(R.r12*Tyy)
             + R.r13*(R.r13*Tzz);
  tsrtemp.Txy= R.r21*(R.r11*Txx)
             + R.r22*(R.r12*Tyy)
             + R.r23*(R.r13*Tzz);
  tsrtemp.Txz= R.r31*(R.r11*Txx)
             + R.r32*(R.r12*Tyy)
             + R.r33*(R.r13*Tzz);
  tsrtemp.Tyx= R.r11*(R.r21*Txx)
             + R.r12*(R.r22*Tyy)
             + R.r13*(R.r23*Tzz);
  tsrtemp.Tyy= R.r21*(R.r21*Txx)
             + R.r22*(R.r22*Tyy)
             + R.r23*(R.r23*Tzz);
  tsrtemp.Tyz= R.r31*(R.r21*Txx)
             + R.r32*(R.r22*Tyy)
             + R.r33*(R.r23*Tzz);
  tsrtemp.Tzx= R.r11*(R.r31*Txx)
             + R.r12*(R.r32*Tyy)
             + R.r13*(R.r33*Tzz);
  tsrtemp.Tzy= R.r21*(R.r31*Txx)
             + R.r22*(R.r32*Tyy)
             + R.r23*(R.r33*Tzz);
  tsrtemp.Tzz= R.r31*(R.r31*Txx)
             + R.r32*(R.r32*Tyy)
             + R.r33*(R.r33*Tzz);
  return tsrtemp;
  }

/****************************************************************************/

tgen tgen::BeamTRrotate(trot R)
  {
  tgen tsrtemp;
  tsrtemp.Txx= R.r12*(R.r13*Tzy)
             + R.r13*(R.r12*Tyz);
  tsrtemp.Txy= R.r22*(R.r13*Tzy)
             + R.r23*(R.r12*Tyz);
  tsrtemp.Txz= R.r32*(R.r13*Tzy)
             + R.r33*(R.r12*Tyz);
  tsrtemp.Tyx= R.r12*(R.r23*Tzy)
             + R.r13*(R.r22*Tyz);
  tsrtemp.Tyy= R.r22*(R.r23*Tzy)
             + R.r23*(R.r22*Tyz);
  tsrtemp.Tyz= R.r32*(R.r23*Tzy)
             + R.r33*(R.r22*Tyz);
  tsrtemp.Tzx= R.r12*(R.r33*Tzy)
             + R.r13*(R.r32*Tyz);
  tsrtemp.Tzy= R.r22*(R.r33*Tzy)
             + R.r23*(R.r32*Tyz);
  tsrtemp.Tzz= R.r32*(R.r33*Tzy)
             + R.r33*(R.r32*Tyz);
  return tsrtemp;
  }

/****************************************************************************/

ostream &operator<<(ostream &stream, tgen ti)
  {
  stream << ti.Txx << " " << ti.Txy << " " << ti.Txz << endl;
  stream << ti.Tyx << " " << ti.Tyy << " " << ti.Tyz << endl;
  stream << ti.Tzx << " " << ti.Tzy << " " << ti.Tzz;
  return stream;
  }

/****************************************************************************/

istream &operator>>(istream &stream, tgen &ti)
  {
  stream >> ti.Txx >> ti.Txy >> ti.Txz
	 >> ti.Tyx >> ti.Tyy >> ti.Tyz
	 >> ti.Tzx >> ti.Tzy >> ti.Tzz;
  return stream;
  }

/****************************************************************************/

/****************************************************************************/
/*                      Methods relative to tiner                           */
/****************************************************************************/

tiner::tiner()
  {
  Ixx=0.0; Iyy=0.0; Izz=0.0; Ixy=0.0; Ixz=0.0; Iyz=0.0;
  }

/****************************************************************************/

tiner::tiner(double Jxx,double Jyy, double Jzz,
             double Jxy, double Jxz, double Jyz)
  {
  Ixx=Jxx; Iyy=Jyy; Izz=Jzz; 
  Ixy=Jxy; Ixz=Jxz; Iyz=Jyz;
  }

/****************************************************************************/

void tiner::put(double Jxx,double Jyy,double Jzz,
                double Jxy, double Jxz, double Jyz)
  {
  Ixx=Jxx; Iyy=Jyy; Izz=Jzz; 
  Ixy=Jxy; Ixz=Jxz; Iyz=Jyz;
  }

/****************************************************************************/

double tiner::trace()
  {
    return (Ixx+Iyy+Izz);
  }

/****************************************************************************/

tiner tiner::rotate(trot R)
  {
  tiner tinertemp;
  tinertemp.Ixx= R.r11*(Ixx*R.r11 - Ixy*R.r12 - Ixz*R.r13)
               + R.r12*(-Ixy*R.r11 + Iyy*R.r12 - Iyz*R.r13)
               + R.r13*(-Ixz*R.r11 - Iyz*R.r12 + Izz*R.r13);
  tinertemp.Iyy= R.r21*(Ixx*R.r21 - Ixy*R.r22 - Ixz*R.r23)
               + R.r22*(-Ixy*R.r21 + Iyy*R.r22 - Iyz*R.r23)
               + R.r23*(-Ixz*R.r21 - Iyz*R.r22 + Izz*R.r23);
  tinertemp.Izz= R.r31*(Ixx*R.r31 - Ixy*R.r32 - Ixz*R.r33)
               + R.r32*(-Ixy*R.r31 + Iyy*R.r32 - Iyz*R.r33)
               + R.r33*(-Ixz*R.r31 - Iyz*R.r32 + Izz*R.r33);
  tinertemp.Ixy= -R.r21*(Ixx*R.r11 - Ixy*R.r12 - Ixz*R.r13)
                - R.r22*(-Ixy*R.r11 + Iyy*R.r12 - Iyz*R.r13)
                - R.r23*(-Ixz*R.r11 - Iyz*R.r12 + Izz*R.r13);
  tinertemp.Ixz= -R.r31*(Ixx*R.r11 - Ixy*R.r12 - Ixz*R.r13)
                - R.r32*(-Ixy*R.r11 + Iyy*R.r12 - Iyz*R.r13)
                - R.r33*(-Ixz*R.r11 - Iyz*R.r12 + Izz*R.r13);
  tinertemp.Iyz= -R.r31*(Ixx*R.r21 - Ixy*R.r22 - Ixz*R.r23)
                - R.r32*(-Ixy*R.r21 + Iyy*R.r22 - Iyz*R.r23)
                - R.r33*(-Ixz*R.r21 - Iyz*R.r22 + Izz*R.r23);
  return tinertemp;
  }

/****************************************************************************/

vec tiner::operator*(vec v)
  {
  vec vtemp;
  vtemp.x=Ixx*v.x-Ixy*v.y-Ixz*v.z;
  vtemp.y=-Ixy*v.x+Iyy*v.y-Iyz*v.z;
  vtemp.z=-Ixz*v.x-Iyz*v.y+Izz*v.z;
  return vtemp;
  }

/****************************************************************************/

ostream &operator<<(ostream &stream, tiner ti)
  {
  stream << ti.Ixx << " " << ti.Iyy << " " << ti.Izz << " ";
  stream << ti.Ixy << " " << ti.Ixz << " " << ti.Iyz;
  return stream;
  }

/****************************************************************************/

istream &operator>>(istream &stream, tiner &ti)
  {
  stream >> ti.Ixx >> ti.Iyy >> ti.Izz >> ti.Ixy >> ti.Ixz >> ti.Iyz;
  return stream;
  }

/****************************************************************************/

/****************************************************************************/
/*                      Methods relative to tdiag                           */
/****************************************************************************/

tdiag::tdiag()
  {
  Txx=0.0; Tyy=0.0; Tzz=0.0;
  }

/****************************************************************************/

tdiag::tdiag(double Jxx,double Jyy, double Jzz)
  {
  Txx=Jxx; Tyy=Jyy; Tzz=Jzz; 
  }

/****************************************************************************/

void tdiag::put(double Jxx,double Jyy,double Jzz)
  {
  Txx=Jxx; Tyy=Jyy; Tzz=Jzz; 
  }

/****************************************************************************/

double tdiag::trace()
  {
    return (Txx+Tyy+Tzz);
  }

/****************************************************************************/

tgen tdiag::rotate(trot R)
  {
  tgen tsrtemp;
  tsrtemp.Txx= R.r11*R.r11*Txx + R.r12*R.r12*Tyy + R.r13*R.r13*Tzz;
  tsrtemp.Txy= R.r11*R.r21*Txx + R.r12*R.r22*Tyy + R.r13*R.r23*Tzz;
  tsrtemp.Txz= R.r11*R.r31*Txx + R.r12*R.r32*Tyy + R.r13*R.r33*Tzz;
  tsrtemp.Tyx= tsrtemp.Txy;
  tsrtemp.Tyy= R.r21*R.r21*Txx + R.r22*R.r22*Tyy + R.r23*R.r23*Tzz;
  tsrtemp.Tyz= R.r21*R.r31*Txx + R.r22*R.r32*Tyy + R.r23*R.r33*Tzz;
  tsrtemp.Tzx= tsrtemp.Txz;
  tsrtemp.Tzy= tsrtemp.Tyz;
  tsrtemp.Tzz= R.r31*R.r31*Txx + R.r32*R.r32*Tyy + R.r33*R.r33*Tzz;
  return tsrtemp;
  }

/****************************************************************************/

vec tdiag::operator*(vec v)
  {
  vec vtemp;
  vtemp.x=Txx*v.x;
  vtemp.y=Tyy*v.y;
  vtemp.z=Tzz*v.z;
  return vtemp;
  }

/****************************************************************************/

ostream &operator<<(ostream &stream, tdiag ti)
  {
  stream << ti.Txx << " " << ti.Tyy << " " << ti.Tzz;
  return stream;
  }

/****************************************************************************/

istream &operator>>(istream &stream, tdiag &ti)
  {
  stream >> ti.Txx >> ti.Tyy >> ti.Tzz;
  return stream;
  }

/****************************************************************************/

/****************************************************************************/
/*                       Methods relative to mth                            */
/****************************************************************************/

mth::mth()
  {
  }

/****************************************************************************/

void mth::unite()
  {
  R.unite();
  e.zero();
  }

/****************************************************************************/

mth mth::inv()
  {
  mth mtemp;
  mtemp.R.r11=R.r11; mtemp.R.r12=R.r21; mtemp.R.r13=R.r31;
  mtemp.R.r21=R.r12; mtemp.R.r22=R.r22; mtemp.R.r23=R.r32;
  mtemp.R.r31=R.r13; mtemp.R.r32=R.r23; mtemp.R.r33=R.r33;
  mtemp.e.x=-R.r11*e.x-R.r21*e.y-R.r31*e.z;
  mtemp.e.y=-R.r12*e.x-R.r22*e.y-R.r32*e.z;
  mtemp.e.z=-R.r13*e.x-R.r23*e.y-R.r33*e.z;
  return mtemp;
  }

/****************************************************************************/

mth mth::operator*(mth m2)
  {
  mth mtemp;
  mtemp.R.r11=R.r11*m2.R.r11+R.r12*m2.R.r21+R.r13*m2.R.r31;
  mtemp.R.r12=R.r11*m2.R.r12+R.r12*m2.R.r22+R.r13*m2.R.r32;
  mtemp.R.r13=R.r11*m2.R.r13+R.r12*m2.R.r23+R.r13*m2.R.r33;
  mtemp.e.x=R.r11*m2.e.x+R.r12*m2.e.y+R.r13*m2.e.z+e.x;
  mtemp.R.r21=R.r21*m2.R.r11+R.r22*m2.R.r21+R.r23*m2.R.r31;
  mtemp.R.r22=R.r21*m2.R.r12+R.r22*m2.R.r22+R.r23*m2.R.r32;
  mtemp.R.r23=R.r21*m2.R.r13+R.r22*m2.R.r23+R.r23*m2.R.r33;
  mtemp.e.y=R.r21*m2.e.x+R.r22*m2.e.y+R.r23*m2.e.z+e.y;
  mtemp.R.r31=R.r31*m2.R.r11+R.r32*m2.R.r21+R.r33*m2.R.r31;
  mtemp.R.r32=R.r31*m2.R.r12+R.r32*m2.R.r22+R.r33*m2.R.r32;
  mtemp.R.r33=R.r31*m2.R.r13+R.r32*m2.R.r23+R.r33*m2.R.r33;
  mtemp.e.z=R.r31*m2.e.x+R.r32*m2.e.y+R.r33*m2.e.z+e.z;
  return mtemp;
  }

/****************************************************************************/

mth mth::operator*=(mth m2)
  {
  double r1,r2,r3;
  r1=R.r11*m2.R.r11+R.r12*m2.R.r21+R.r13*m2.R.r31;
  r2=R.r11*m2.R.r12+R.r12*m2.R.r22+R.r13*m2.R.r32;
  r3=R.r11*m2.R.r13+R.r12*m2.R.r23+R.r13*m2.R.r33;
  e.x=R.r11*m2.e.x+R.r12*m2.e.y+R.r13*m2.e.z+e.x;
  R.r11=r1; R.r12=r2; R.r13=r3;
  r1=R.r21*m2.R.r11+R.r22*m2.R.r21+R.r23*m2.R.r31;
  r2=R.r21*m2.R.r12+R.r22*m2.R.r22+R.r23*m2.R.r32;
  r3=R.r21*m2.R.r13+R.r22*m2.R.r23+R.r23*m2.R.r33;
  e.y=R.r21*m2.e.x+R.r22*m2.e.y+R.r23*m2.e.z+e.y;
  R.r21=r1; R.r22=r2; R.r23=r3;
  r1=R.r31*m2.R.r11+R.r32*m2.R.r21+R.r33*m2.R.r31;
  r2=R.r31*m2.R.r12+R.r32*m2.R.r22+R.r33*m2.R.r32;
  r3=R.r31*m2.R.r13+R.r32*m2.R.r23+R.r33*m2.R.r33;
  e.z=R.r31*m2.e.x+R.r32*m2.e.y+R.r33*m2.e.z+e.z;
  R.r31=r1; R.r32=r2; R.r33=r3;
  return *this;
  }

/****************************************************************************/

vec mth::operator*(vec v)
  {
  vec vtemp;
  vtemp.x=R.r11*v.x+R.r12*v.y+R.r13*v.z+e.x;
  vtemp.y=R.r21*v.x+R.r22*v.y+R.r23*v.z+e.y;
  vtemp.z=R.r31*v.x+R.r32*v.y+R.r33*v.z+e.z;
  return vtemp;
  }

/****************************************************************************/

ostream &operator<<(ostream &stream, mth m)
  {
  stream << m.R.r11 <<" "<< m.R.r12 <<" "<< m.R.r13 <<" "<< m.e.x << "\n";
  stream << m.R.r21 <<" "<< m.R.r22 <<" "<< m.R.r23 <<" "<< m.e.y << "\n";
  stream << m.R.r31 <<" "<< m.R.r32 <<" "<< m.R.r33 <<" "<< m.e.z << "\n";
  return stream;
  }

/****************************************************************************/

istream &operator>>(istream &stream, mth &m)
  {
  stream >> m.R.r11 >> m.R.r12 >> m.R.r13 >> m.e.x;
  stream >> m.R.r21 >> m.R.r22 >> m.R.r23 >> m.e.y;
  stream >> m.R.r31 >> m.R.r32 >> m.R.r33 >> m.e.z;
  return stream;
  }

/****************************************************************************/



/****************************************************************************/
/*                             Utlility routines                            */
/****************************************************************************/

/****************************************************************************/

vec vcoord(double x, double y, double z)
  {
  vec vtemp;
  vtemp.x=x;
  vtemp.y=y;
  vtemp.z=z;
  return vtemp;
  }

/****************************************************************************/

trot Rrotx(double th)
  {
  trot trtemp;
  trtemp.rotx(th);
  return trtemp;
  }

/****************************************************************************/

trot Rroty(double th)
  {
  trot trtemp;
  trtemp.roty(th);
  return trtemp;
  }

/****************************************************************************/

trot Rrotz(double th)
  {
  trot trtemp;
  trtemp.rotz(th);
  return trtemp;
  }

/****************************************************************************/

trot Rrotn(vec v, double th)
  {
  trot trtemp;
  trtemp.rotn(v.x, v.y, v.z, th);
  return trtemp;
  }

/****************************************************************************/

trot Rrotn(double nx, double ny, double nz, double th)
  {
  trot trtemp;
  trtemp.rotn(nx, ny, nz, th);
  return trtemp;
  }

/****************************************************************************/

mth Trotx(double th)
  {
  mth mtemp;
  mtemp.R.rotx(th);
  return mtemp;
  }

/****************************************************************************/

mth Troty(double th)
  {
  mth mtemp;
  mtemp.R.roty(th);
  return mtemp;
  }

/****************************************************************************/

mth Trotz(double th)
  {
  mth mtemp;
  mtemp.R.rotz(th);
  return mtemp;
  }

/****************************************************************************/

mth Trotn(vec v, double th)
  {
  mth mtemp;
  mtemp.R.rotn(v.x, v.y, v.z, th);
  return mtemp;
  }

/****************************************************************************/

mth Trotn(double nx, double ny, double nz, double th)
  {
  mth mtemp;
  mtemp.R.rotn(nx, ny, nz, th);
  return mtemp;
  }

/****************************************************************************/

mth Tdisp(vec v)
  {
  mth mtemp;
  mtemp.e=v;
  return mtemp;
  }

/****************************************************************************/

mth Tdisp(double x, double y, double z)
  {
  mth mtemp;
  mtemp.e.x=x;
  mtemp.e.y=y;
  mtemp.e.z=z;
  return mtemp;
  }

/****************************************************************************/

/****************************************************************************/
/*                        Decomposition routines                            */
/****************************************************************************/

/****************************************************************************/

void getrotxyz(trot tr1, trot tr2, double &th1, double &th2, double &th3,
                            vec &axe1, vec &axe2, vec &axe3)
// Determination des angles menant de m1 a m2 par des rotations successives
// autour d'axes locaux x, y et z
  {
  // Copie des axes extremes
  axe1=tr1.ux();
  axe3=tr2.uz();
  // Determination de l'axe intermediaire
  axe2=axe3^axe1;  // ou l'oppose mais il faut bien choisir !
  if (axe2.length()<1e-6) axe2=tr2.uy(); // singularite !!
  else axe2/=axe2.length();
  vec zp=axe1^axe2; // axe z du 1er repere intermediaire
  // Determination des angles de rotation
  th1=atan2(axe2*tr1.uz(),axe2*tr1.uy());
  th2=atan2(axe3*axe1,axe3*zp);
  th3=atan2(axe2*tr2.ux(),axe2*tr2.uy());
  }

/****************************************************************************/

void getrotyzx(trot tr1, trot tr2, double &th1, double &th2, double &th3,
                            vec &axe1, vec &axe2, vec &axe3)
// Determination des angles menant de m1 a m2 par des rotations successives
// autour d'axes locaux y, z et x
  {
  // Copie des axes extremes
  axe1=tr1.uy();
  axe3=tr2.ux();
  // Determination de l'axe intermediaire
  axe2=axe3^axe1;  // ou l'oppose mais il faut bien choisir !
  if (axe2.length()<1e-6) axe2=tr2.uz(); // singularite !!
  else axe2/=axe2.length();
  vec xp=axe1^axe2; // axe x du 1er repere intermediaire
  // Determination des angles de rotation
  th1=atan2(axe2*tr1.ux(),axe2*tr1.uz());
  th2=atan2(axe3*axe1,axe3*xp);
  th3=atan2(axe2*tr2.uy(),axe2*tr2.uz());
  }

/****************************************************************************/

void getrotzxy(trot tr1, trot tr2, double &th1, double &th2, double &th3,
                            vec &axe1, vec &axe2, vec &axe3)
// Determination des angles menant de m1 a m2 par des rotations successives
// autour d'axes locaux z, x et y
  {
  // Copie des axes extremes
  axe1=tr1.uz();
  axe3=tr2.uy();
  // Determination de l'axe intermediaire
  axe2=axe3^axe1;  // ou l'oppose mais il faut bien choisir !
  if (axe2.length()<1e-6) axe2=tr2.ux(); // singularite !!
  else axe2/=axe2.length();
  vec yp=axe1^axe2; // axe y du 1er repere intermediaire
  // Determination des angles de rotation
  th1=atan2(axe2*tr1.uy(),axe2*tr1.ux());
  th2=atan2(axe3*axe1,axe3*yp);
  th3=atan2(axe2*tr2.uz(),axe2*tr2.ux());
  }

/****************************************************************************/

void getrotyxz(trot tr1, trot tr2, double &th1, double &th2, double &th3,
                            vec &axe1, vec &axe2, vec &axe3)
// Determination des angles menant de m1 a m2 par des rotations successives
// autour d'axes locaux y, x et z
  {
  // Copie des axes extremes
  axe1=tr1.uy();
  axe3=tr2.uz();
  // Determination de l'axe intermediaire
  axe2=axe1^axe3;  // ou l'oppose mais il faut bien choisir !
  if (axe2.length()<1e-6) axe2=tr2.ux(); // singularite !!
  else axe2/=axe2.length();
  vec zp=axe2^axe1; // axe z du 1er repere intermediaire
  // Determination des angles de rotation
  th1=atan2(-(axe2*tr1.uz()),axe2*tr1.ux());
  th2=atan2(-(axe3*axe1),axe3*zp);
  th3=atan2(-(axe2*tr2.uy()),axe2*tr2.ux());
  }

/****************************************************************************/

void getrotxzy(trot tr1, trot tr2, double &th1, double &th2, double &th3,
                            vec &axe1, vec &axe2, vec &axe3)
// Determination des angles menant de m1 a m2 par des rotations successives
// autour d'axes locaux x, z et y
  {
  // Copie des axes extremes
  axe1=tr1.ux();
  axe3=tr2.uy();
  // Determination de l'axe intermediaire
  axe2=axe1^axe3;  // ou l'oppose mais il faut bien choisir !
  if (axe2.length()<1e-6) axe2=tr2.uz(); // singularite !!
  else axe2/=axe2.length();
  vec yp=axe2^axe1; // axe y du 1er repere intermediaire
  // Determination des angles de rotation
  th1=atan2(-(axe2*tr1.uy()),axe2*tr1.uz());
  th2=atan2(-(axe3*axe1),axe3*yp);
  th3=atan2(-(axe2*tr2.ux()),axe2*tr2.uz());
  }

/****************************************************************************/

void getrotzyx(trot tr1, trot tr2, double &th1, double &th2, double &th3,
                            vec &axe1, vec &axe2, vec &axe3)
// Determination des angles menant de m1 a m2 par des rotations successives
// autour d'axes locaux z, y et x
  {
  // Copie des axes extremes
  axe1=tr1.uz();
  axe3=tr2.ux();
  // Determination de l'axe intermediaire
  axe2=axe1^axe3;  // ou l'oppose mais il faut bien choisir !
  if (axe2.length()<1e-6) axe2=tr2.uy(); // singularite !!
  else axe2/=axe2.length();
  vec xp=axe2^axe1; // axe x du 1er repere intermediaire
  // Determination des angles de rotation
  th1=atan2(-(axe2*tr1.ux()),axe2*tr1.uy());
  th2=atan2(-(axe3*axe1),axe3*xp);
  th3=atan2(-(axe2*tr2.uz()),axe2*tr2.uy());
  }

/****************************************************************************/

int decomposevec(vec a,vec axe1, vec axe2, vec axe3, 
		     double &a1,double &a2,double &a3)
  // Calcul des composantes d'un vecteur dans une base non cartesienne
  {
  // Calcul du volume construit sur les 3 vecteurs de la base
  double V=axe1*(axe2^axe3);
  if (fabs(V)<1e-6) return 1; // les vecteurs sont coplanaires !
  // calcul des composantes par la theorie des vecteurs reciproques
  a1=a*(axe2^axe3)/V;
  a2=a*(axe3^axe1)/V;
  a3=a*(axe1^axe2)/V;
  return 0;
  }

/****************************************************************************/
// Thats'all Folks !
