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
#include <stdlib.h>

#include <EasyDyn/toeol.h>
#include <EasyDyn/spline.h>

// fonctions utilitaires pour les poutres flexibles, le point sur sol
// et le point ou segment sur courbe correspondant aux fonctions de forme
// transversales des poutres flexibles et aux splines cubiques

double N22(double xi)     { return xi*xi*(3-2*xi); }
double N26(double xi)     { return xi*xi*(xi-1);   }
double N1(double u)       { return u*u*(3-2*u);    }
double N2(double u)       { return u*u*(u-1);      }
double dN1_du(double u)   { return 6*u*(1-u);      }
double dN2_du(double u)   { return u*(3*u-2);      }
double d2N1_du2(double u) { return 6-12*u;         }
double d2N2_du2(double u) { return 6*u-2;          }

/****************************************************************************/

void spline::set(int nbrtr2, double *Xnode,double *Ynode,
			     double *Xtnode, double *Ytnode)
  {
  nbrtr=nbrtr2;
  Xli=new double[nbrtr+1];
  Yli=new double[nbrtr+1];
  Xtli=new double[nbrtr+1];
  Ytli=new double[nbrtr+1];
  Qli=new double[nbrtr+1];
  if ((!Xli) || (!Yli) || (!Xtli) || (!Ytli) || (!Qli))
    {
    cout << "Allocation error in spline" << "\n";
    exit(1);
    }
  int itr;
  for (itr=0;itr<=nbrtr;itr++)
    {
    Xli[itr]=Xnode[itr];   Yli[itr]=Ynode[itr];
    Xtli[itr]=Xtnode[itr]; Ytli[itr]=Ytnode[itr];
    }
  InitQli();
  }

/****************************************************************************/

void spline::Get(int &nbrtr2, double *&Xnode,double *&Ynode,
			      double *&Xtnode, double *&Ytnode)
  {
  nbrtr2=nbrtr;
  Xnode=Xli;
  Ynode=Yli;
  Xtnode=Xtli;
  Ytnode=Ytli;
  }

/****************************************************************************/

void spline::ReadSpline(istream &stream)
  {
  cout << "Lecture de nbrtr" << "\n";
  stream >> nbrtr >> toeol;
  // Attribution de memoire pour les differents vecs
  Xli=new double[nbrtr+1];
  Yli=new double[nbrtr+1];
  Xtli=new double[nbrtr+1];
  Ytli=new double[nbrtr+1];
  Qli=new double[nbrtr+1];
  int itr;
  cout << "Lecture des coordonnees du spline" << "\n";
  for (itr=0; itr<=nbrtr; itr++)
    stream >> Xli[itr] >> Yli[itr] >> Xtli[itr] >> Ytli[itr] >> toeol;
  cout << "Initialisation de Qli" << "\n";
  InitQli();
  }

/****************************************************************************/

void spline::WriteSpline(ostream &stream)
  {
  // Ecriture nombre de troncons constituant la courbe
  stream << nbrtr << " =nbrtr" << "\n";
  int itr;
  for (itr=0; itr<=nbrtr; itr++)
    stream << Xli[itr] << " " << Yli[itr] << " " << Xtli[itr]
	   << " " << Ytli[itr] << " =X Y XT YT au noeud " << itr << "\n";
  }

/****************************************************************************/

// fonction d'initialisation du vec Qli et de deltaqth pour un spline

void spline::InitQli()
  {
  int itr,iu;
  double u,dX_du,dY_du,ds_du,d2X_du2,d2Y_du2,d2s_du2;
  Qli[0]=0.0;
  deltaqth=1E12;
  for (itr=1; itr<=nbrtr; itr++)
    {
    // Calcul de la contribution du premier point du segment
    dX_du=Xtli[itr-1];
    dY_du=Ytli[itr-1];
    ds_du=sqrt(dX_du*dX_du+dY_du*dY_du);
    d2X_du2=-6*Xli[itr-1]-4*Xtli[itr-1]+6*Xli[itr]-2*Xtli[itr];
    d2Y_du2=-6*Yli[itr-1]-4*Ytli[itr-1]+6*Yli[itr]-2*Ytli[itr];
    d2s_du2=(dX_du*d2X_du2+dY_du*d2Y_du2)/ds_du;
    Qli[itr]=Qli[itr-1]+0.005*ds_du+0.0001*d2s_du2/12;
    // Calcul de la contribution du dernier point du segment
    dX_du=Xtli[itr];
    dY_du=Ytli[itr];
    ds_du=sqrt(dX_du*dX_du+dY_du*dY_du);
    d2X_du2=6*Xli[itr-1]+2*Xtli[itr-1]-6*Xli[itr]+4*Xtli[itr];
    d2Y_du2=6*Yli[itr-1]+2*Ytli[itr-1]-6*Yli[itr]+4*Ytli[itr];
    d2s_du2=(dX_du*d2X_du2+dY_du*d2Y_du2)/ds_du;
    Qli[itr]+=0.005*ds_du-0.0001*d2s_du2/12;
    // Calcul de la contribution de 99 points dans le segment
    for (iu=1; iu<100; iu++)
      {
      u=iu*0.01;
      dX_du=-Xli[itr-1]*dN1_du(1-u)+Xli[itr]*dN1_du(u)
	    +Xtli[itr-1]*dN2_du(1-u)+Xtli[itr]*dN2_du(u);
      dY_du=-Yli[itr-1]*dN1_du(1-u)+Yli[itr]*dN1_du(u)
	    +Ytli[itr-1]*dN2_du(1-u)+Ytli[itr]*dN2_du(u);
      Qli[itr]+=0.01*sqrt(dX_du*dX_du+dY_du*dY_du);
      }
    if ((Qli[itr]-Qli[itr-1])<(20*deltaqth)) deltaqth=(Qli[itr]-Qli[itr-1])/20;
    }
  }

/****************************************************************************/

// fonction de calcul des caracteristiques en un point tangent d'un spline

void spline::GetPttg(double qs,
		     double &xpt, double &dx_ds,  double &d2x_ds2,
		     double &ypt, double &dy_ds,  double &d2y_ds2,
		     double &th,  double &dth_ds, double &d2th_ds2)
  {
  double qA=qs-deltaqth, qB=qs+deltaqth;
  if ((qs<0) || (qs>Qli[nbrtr]))
    {
    cout << "Spline out of bounds" << "\n";
    exit(1);
    }
  // Recherche des numeros des troncons ou se trouvent les points A,P et B
  int itr,numtr,numtrA=0,numtrB=nbrtr+1;
  for (itr=1; itr<=nbrtr; itr++)
    {
    if ((qs>=Qli[itr-1]) && (qs<=Qli[itr])) numtr=itr;
    if ((qA>=Qli[itr-1])   && (qA<=Qli[itr]))   numtrA=itr;
    if ((qB>=Qli[itr-1])   && (qB<=Qli[itr]))   numtrB=itr;
    }
  // Recherche des coordonnees exactes des points A,P et B
  double Xl0,Yl0,Xtl0,Ytl0,Xl1,Yl1,Xtl1,Ytl1,
	 dq_du, dx_du, dy_du, d2x_du2, d2y_du2,
	 dqA_du, dxA_du, dyA_du, d2xA_du2, d2yA_du2,
	 xA, yA, dxA_ds, dyA_ds, d2xA_ds2, d2yA_ds2,
	 dqB_du, dxB_du, dyB_du, d2xB_du2, d2yB_du2,
	 xB, yB, dxB_ds, dyB_ds, d2xB_ds2, d2yB_ds2,
	 dx1_du, dy1_du, ds1_du, d2x1_du2, d2y1_du2, d2s1_du2,
	 dx2_du, dy2_du, ds2_du, d2x2_du2, d2y2_du2, d2s2_du2,
	 u, uA, uB, usto, qtr, deltaq, xi;
  // Calcul des caracteristiques au point A s'il est en dehors du spline
  if (numtrA==0)
    {
    double utAx=Xtli[0]/sqrt(Xtli[0]*Xtli[0]+Ytli[0]*Ytli[0]),
	   utAy=Ytli[0]/sqrt(Xtli[0]*Xtli[0]+Ytli[0]*Ytli[0]);
    xA=Xli[0]+qA*utAx;
    yA=Yli[0]+qA*utAy;
    dxA_ds=utAx;
    dyA_ds=utAy;
    d2xA_ds2=0;
    d2yA_ds2=0;
    }
  // Recherche de la coordonnee uA correspondant a l'abscisse curviligne qA
  // s'il n'est pas dans le meme troncon que le point P
  // l'integration se fait ici en arriere
  if ((numtrA>0) && (numtrA<numtr))
    {
    Xl0=Xli[numtrA-1];   Yl0=Yli[numtrA-1];
    Xtl0=Xtli[numtrA-1]; Ytl0=Ytli[numtrA-1];
    Xl1=Xli[numtrA];     Yl1=Yli[numtrA];
    Xtl1=Xtli[numtrA];   Ytl1=Ytli[numtrA];
    uA=1.0; qtr=Qli[numtrA];
    dx2_du=Xtl1; dy2_du=Ytl1;
    ds2_du=sqrt(dx2_du*dx2_du+dy2_du*dy2_du);
    d2x2_du2=6*Xl0+2*Xtl0-6*Xl1+4*Xtl1;
    d2y2_du2=6*Yl0+2*Ytl0-6*Yl1+4*Ytl1;
    d2s2_du2=(dx2_du*d2x2_du2+dy2_du*d2y2_du2)/ds2_du;
    // initialisation
    deltaq=1;
    ds1_du=ds2_du;
    d2s1_du2=d2s2_du2;
    while (qtr>qA)
      {
      ds2_du=ds1_du;
      d2s2_du2=d2s1_du2;
      uA=uA-0.01;
      dx1_du=-Xl0*dN1_du(1-uA)+Xl1*dN1_du(uA)+Xtl0*dN2_du(1-uA)+Xtl1*dN2_du(uA);
      dy1_du=-Yl0*dN1_du(1-uA)+Yl1*dN1_du(uA)+Ytl0*dN2_du(1-uA)+Ytl1*dN2_du(uA);
      ds1_du=sqrt(dx1_du*dx1_du+dy1_du*dy1_du);
      d2x1_du2=Xl0*d2N1_du2(1-uA)+Xl1*d2N1_du2(uA)-Xtl0*d2N2_du2(1-uA)+Xtl1*d2N2_du2(uA);
      d2y1_du2=Yl0*d2N1_du2(1-uA)+Yl1*d2N1_du2(uA)-Ytl0*d2N2_du2(1-uA)+Ytl1*d2N2_du2(uA);
      d2s1_du2=(dx1_du*d2x1_du2+dy1_du*d2y1_du2)/ds1_du;
      deltaq=0.005*(ds1_du+ds2_du)+0.0001*(d2s1_du2-d2s2_du2)/12;
      qtr-=deltaq;
      }
    // uA=uA+0.01*(qA-qtr)/deltaq;
    xi=(qA-qtr)/deltaq;
    uA=uA*N1(1-xi)+(uA+0.01)*N1(xi)+deltaq*(-N2(1-xi)/ds1_du+N2(xi)/ds2_du);
    // Calcul de la position et des derivees au point uA
    xA=Xl0*N1(1-uA)+Xl1*N1(uA)-Xtl0*N2(1-uA)+Xtl1*N2(uA);
    yA=Yl0*N1(1-uA)+Yl1*N1(uA)-Ytl0*N2(1-uA)+Ytl1*N2(uA);
    dxA_du=-Xl0*dN1_du(1-uA)+Xl1*dN1_du(uA)+Xtl0*dN2_du(1-uA)+Xtl1*dN2_du(uA);
    dyA_du=-Yl0*dN1_du(1-uA)+Yl1*dN1_du(uA)+Ytl0*dN2_du(1-uA)+Ytl1*dN2_du(uA);
    d2xA_du2=Xl0*d2N1_du2(1-uA)+Xl1*d2N1_du2(uA)
	    -Xtl0*d2N2_du2(1-uA)+Xtl1*d2N2_du2(uA);
    d2yA_du2=Yl0*d2N1_du2(1-uA)+Yl1*d2N1_du2(uA)
	   -Ytl0*d2N2_du2(1-uA)+Ytl1*d2N2_du2(uA);
    dqA_du=sqrt(dxA_du*dxA_du+dyA_du*dyA_du);
    dxA_ds=dxA_du/dqA_du;
    dyA_ds=dyA_du/dqA_du;
    d2xA_ds2=(d2xA_du2*dyA_du*dyA_du-dxA_du*dyA_du*d2yA_du2)
	     /(dqA_du*dqA_du*dqA_du*dqA_du);
    d2yA_ds2=(d2yA_du2*dxA_du*dxA_du-dxA_du*dyA_du*d2xA_du2)
	     /(dqA_du*dqA_du*dqA_du*dqA_du);
    }
  // Recherche de la coordonnee u correspondant a l'abscisse curviligne qs
  Xl0=Xli[numtr-1];   Yl0=Yli[numtr-1];
  Xtl0=Xtli[numtr-1]; Ytl0=Ytli[numtr-1];
  Xl1=Xli[numtr];     Yl1=Yli[numtr];
  Xtl1=Xtli[numtr];   Ytl1=Ytli[numtr];
  u=0.0; qtr=Qli[numtr-1];
  dx1_du=Xtl0;
  dy1_du=Ytl0;
  ds1_du=sqrt(dx1_du*dx1_du+dy1_du*dy1_du);
  d2x1_du2=-6*Xl0-4*Xtl0+6*Xl1-2*Xtl1;
  d2y1_du2=-6*Yl0-4*Ytl0+6*Yl1-2*Ytl1;
  d2s1_du2=(dx1_du*d2x1_du2+dy1_du*d2y1_du2)/ds1_du;
  // initialisation
  deltaq=1; ds2_du=ds1_du; d2s2_du2=d2s1_du2;
  while (qtr<qs)
    {
    ds1_du=ds2_du;
    d2s1_du2=d2s2_du2;
    u=u+0.01;
    dx2_du=-Xl0*dN1_du(1-u)+Xl1*dN1_du(u)+Xtl0*dN2_du(1-u)+Xtl1*dN2_du(u);
    dy2_du=-Yl0*dN1_du(1-u)+Yl1*dN1_du(u)+Ytl0*dN2_du(1-u)+Ytl1*dN2_du(u);
    ds2_du=sqrt(dx2_du*dx2_du+dy2_du*dy2_du);
    d2x2_du2=Xl0*d2N1_du2(1-u)+Xl1*d2N1_du2(u)-Xtl0*d2N2_du2(1-u)+Xtl1*d2N2_du2(u);
    d2y2_du2=Yl0*d2N1_du2(1-u)+Yl1*d2N1_du2(u)-Ytl0*d2N2_du2(1-u)+Ytl1*d2N2_du2(u);
    d2s2_du2=(dx2_du*d2x2_du2+dy2_du*d2y2_du2)/ds2_du;
    deltaq=0.005*(ds1_du+ds2_du)+0.0001*(d2s1_du2-d2s2_du2)/12;
    qtr+=deltaq;
    if ((qA<=qtr) && (qA>(qtr-deltaq)))
      {
      // uA=u-0.01*(qtr-qA)/deltaq;
      xi=(qA-qtr+deltaq)/deltaq;
      uA=(u-0.01)*N1(1-xi)+u*N1(xi)+deltaq*(-N2(1-xi)/ds1_du+N2(xi)/ds2_du);
      // Calcul de la position et des derivees au point uA
      xA=Xl0*N1(1-uA)+Xl1*N1(uA)-Xtl0*N2(1-uA)+Xtl1*N2(uA);
      yA=Yl0*N1(1-uA)+Yl1*N1(uA)-Ytl0*N2(1-uA)+Ytl1*N2(uA);
      dxA_du=-Xl0*dN1_du(1-uA)+Xl1*dN1_du(uA)+Xtl0*dN2_du(1-uA)+Xtl1*dN2_du(uA);
      dyA_du=-Yl0*dN1_du(1-uA)+Yl1*dN1_du(uA)+Ytl0*dN2_du(1-uA)+Ytl1*dN2_du(uA);
      d2xA_du2=Xl0*d2N1_du2(1-uA)+Xl1*d2N1_du2(uA)
	      -Xtl0*d2N2_du2(1-uA)+Xtl1*d2N2_du2(uA);
      d2yA_du2=Yl0*d2N1_du2(1-uA)+Yl1*d2N1_du2(uA)
	      -Ytl0*d2N2_du2(1-uA)+Ytl1*d2N2_du2(uA);
      dqA_du=sqrt(dxA_du*dxA_du+dyA_du*dyA_du);
      dxA_ds=dxA_du/dqA_du;
      dyA_ds=dyA_du/dqA_du;
      d2xA_ds2=(d2xA_du2*dyA_du*dyA_du-dxA_du*dyA_du*d2yA_du2)
	       /(dqA_du*dqA_du*dqA_du*dqA_du);
      d2yA_ds2=(d2yA_du2*dxA_du*dxA_du-dxA_du*dyA_du*d2xA_du2)
	       /(dqA_du*dqA_du*dqA_du*dqA_du);
      }
    }
  usto=u; // stockage de u pour continuer l'integration si B est dans le
	  // meme troncon que P
  // calcul precis de u
  /*double u_I=usto-0.01, u_II=usto, s1=qtr-deltaq, s_I=s1, s_II=qtr,
	 dtmp, xi, xi2;
  while (fabs(qtr-qs)>1E-6)
    {
    u=u_I+(u_II-u_I)*(qs-s_I)/(s_II-s_I);
    xi=(u+0.01-usto)/0.01; xi2=xi*xi;
    dtmp=0.01*(ds1_du*(xi-xi*xi2+xi2*xi2/2)+ds2_du*(xi*xi2-xi2*xi2/2))
	+1E-4*(d2s1_du2*(xi2/2-2*xi*xi2/3+xi2*xi2/4)
	      +d2s2_du2*(xi2*xi2/4-xi*xi2/3));
    qtr=s1+dtmp;
    if (qtr>qs)  { s_II=qtr; u_II=u; }
    if (qtr<=qs) { s_I=qtr;  u_I=u; }
    } */
  xi=(qs-qtr+deltaq)/deltaq;
  u=(u-0.01)*N1(1-xi)+u*N1(xi)+deltaq*(-N2(1-xi)/ds1_du+N2(xi)/ds2_du);
  // Calcul de la position et des derivees au point u
  xpt=Xl0*N1(1-u)+Xl1*N1(u)-Xtl0*N2(1-u)+Xtl1*N2(u);
  ypt=Yl0*N1(1-u)+Yl1*N1(u)-Ytl0*N2(1-u)+Ytl1*N2(u);
  dx_du=-Xl0*dN1_du(1-u)+Xl1*dN1_du(u)+Xtl0*dN2_du(1-u)+Xtl1*dN2_du(u);
  dy_du=-Yl0*dN1_du(1-u)+Yl1*dN1_du(u)+Ytl0*dN2_du(1-u)+Ytl1*dN2_du(u);
  d2x_du2=Xl0*d2N1_du2(1-u)+Xl1*d2N1_du2(u)
	-Xtl0*d2N2_du2(1-u)+Xtl1*d2N2_du2(u);
  d2y_du2=Yl0*d2N1_du2(1-u)+Yl1*d2N1_du2(u)
	    -Ytl0*d2N2_du2(1-u)+Ytl1*d2N2_du2(u);
  dq_du=sqrt(dx_du*dx_du+dy_du*dy_du);
  dx_ds=dx_du/dq_du;
  dy_ds=dy_du/dq_du;
  d2x_ds2=(d2x_du2*dy_du*dy_du-dx_du*dy_du*d2y_du2)/(dq_du*dq_du*dq_du*dq_du);
  d2y_ds2=(d2y_du2*dx_du*dx_du-dx_du*dy_du*d2x_du2)/(dq_du*dq_du*dq_du*dq_du);
  // Recherche du point B s'il est dans le meme troncon que P
  uB=usto;
  if (numtrB==numtr) while (qtr<qB)
    {
    ds1_du=ds2_du;
    d2s1_du2=d2s2_du2;
    uB=uB+0.01;
    dx2_du=-Xl0*dN1_du(1-uB)+Xl1*dN1_du(uB)+Xtl0*dN2_du(1-uB)+Xtl1*dN2_du(uB);
    dy2_du=-Yl0*dN1_du(1-uB)+Yl1*dN1_du(uB)+Ytl0*dN2_du(1-uB)+Ytl1*dN2_du(uB);
    ds2_du=sqrt(dx2_du*dx2_du+dy2_du*dy2_du);
    d2x2_du2=Xl0*d2N1_du2(1-uB)+Xl1*d2N1_du2(uB)-Xtl0*d2N2_du2(1-uB)+Xtl1*d2N2_du2(uB);
    d2y2_du2=Yl0*d2N1_du2(1-uB)+Yl1*d2N1_du2(uB)-Ytl0*d2N2_du2(1-uB)+Ytl1*d2N2_du2(uB);
    d2s2_du2=(dx2_du*d2x2_du2+dy2_du*d2y2_du2)/ds2_du;
    deltaq=0.005*(ds1_du+ds2_du)+0.0001*(d2s1_du2-d2s2_du2)/12;
    qtr+=deltaq;
    }
  // uB=u-0.01*(qtr-qB)/deltaq;
  xi=(qB-qtr+deltaq)/deltaq;
  uB=(uB-0.01)*N1(1-xi)+uB*N1(xi)+deltaq*(-N2(1-xi)/ds1_du+N2(xi)/ds2_du);
  // Recherche de uB s'il n'est pas dans le meme troncon que P
  if ((numtrB<=nbrtr) && (numtrB>numtr))
    {
    Xl0=Xli[numtrB-1];   Yl0=Yli[numtrB-1];
    Xtl0=Xtli[numtrB-1]; Ytl0=Ytli[numtrB-1];
    Xl1=Xli[numtrB];     Yl1=Yli[numtrB];
    Xtl1=Xtli[numtrB];   Ytl1=Ytli[numtrB];
    uB=0.0; qtr=Qli[numtrB-1];
    dx1_du=Xtl0;
    dy1_du=Ytl0;
    ds1_du=sqrt(dx1_du*dx1_du+dy1_du*dy1_du);
    d2x1_du2=-6*Xl0-4*Xtl0+6*Xl1-2*Xtl1;
    d2y1_du2=-6*Yl0-4*Ytl0+6*Yl1-2*Ytl1;
    d2s1_du2=(dx1_du*d2x1_du2+dy1_du*d2y1_du2)/ds1_du;
    // initialisation
    deltaq=1; ds2_du=ds1_du; d2s2_du2=d2s1_du2;
    while (qtr<qB)
      {
      ds1_du=ds2_du;
      d2s1_du2=d2s2_du2;
      uB=uB+0.01;
      dx2_du=-Xl0*dN1_du(1-uB)+Xl1*dN1_du(uB)+Xtl0*dN2_du(1-uB)+Xtl1*dN2_du(uB);
      dy2_du=-Yl0*dN1_du(1-uB)+Yl1*dN1_du(uB)+Ytl0*dN2_du(1-uB)+Ytl1*dN2_du(uB);
      ds2_du=sqrt(dx2_du*dx2_du+dy2_du*dy2_du);
      d2x2_du2=Xl0*d2N1_du2(1-uB)+Xl1*d2N1_du2(uB)-Xtl0*d2N2_du2(1-uB)+Xtl1*d2N2_du2(uB);
      d2y2_du2=Yl0*d2N1_du2(1-uB)+Yl1*d2N1_du2(uB)-Ytl0*d2N2_du2(1-uB)+Ytl1*d2N2_du2(uB);
      d2s2_du2=(dx2_du*d2x2_du2+dy2_du*d2y2_du2)/ds2_du;
      deltaq=0.005*(ds1_du+ds2_du)+0.0001*(d2s1_du2-d2s2_du2)/12;
      qtr+=deltaq;
      }
    // uB=uB-0.01*(qtr-qB)/deltaq;
    xi=(qB-qtr+deltaq)/deltaq;
    uB=(uB-0.01)*N1(1-xi)+uB*N1(xi)+deltaq*(-N2(1-xi)/ds1_du+N2(xi)/ds2_du);
    }
  // Calcul de la position et des derivees au point uB s'il est dans
  // le spline (numtrB<=nbrtr)
  if (numtrB<=nbrtr)
    {
    xB=Xl0*N1(1-uB)+Xl1*N1(uB)-Xtl0*N2(1-uB)+Xtl1*N2(uB);
    yB=Yl0*N1(1-uB)+Yl1*N1(uB)-Ytl0*N2(1-uB)+Ytl1*N2(uB);
    dxB_du=-Xl0*dN1_du(1-uB)+Xl1*dN1_du(uB)+Xtl0*dN2_du(1-uB)+Xtl1*dN2_du(uB);
    dyB_du=-Yl0*dN1_du(1-uB)+Yl1*dN1_du(uB)+Ytl0*dN2_du(1-uB)+Ytl1*dN2_du(uB);
    d2xB_du2=Xl0*d2N1_du2(1-uB)+Xl1*d2N1_du2(uB)
	    -Xtl0*d2N2_du2(1-uB)+Xtl1*d2N2_du2(uB);
    d2yB_du2=Yl0*d2N1_du2(1-uB)+Yl1*d2N1_du2(uB)
	    -Ytl0*d2N2_du2(1-uB)+Ytl1*d2N2_du2(uB);
    dqB_du=sqrt(dxB_du*dxB_du+dyB_du*dyB_du);
    dxB_ds=dxB_du/dqB_du;
    dyB_ds=dyB_du/dqB_du;
    d2xB_ds2=(d2xB_du2*dyB_du*dyB_du-dxB_du*dyB_du*d2yB_du2)
	     /(dqB_du*dqB_du*dqB_du*dqB_du);
    d2yB_ds2=(d2yB_du2*dxB_du*dxB_du-dxB_du*dyB_du*d2xB_du2)
	     /(dqB_du*dqB_du*dqB_du*dqB_du);
    }
  // Calcul des caracteristiques au point B s'il est en dehors du spline
  if (numtrB==(nbrtr+1))
    {
    double utBx=Xtli[nbrtr]/sqrt(Xtli[nbrtr]*Xtli[nbrtr]+Ytli[nbrtr]*Ytli[nbrtr]),
	   utBy=Ytli[nbrtr]/sqrt(Xtli[nbrtr]*Xtli[nbrtr]+Ytli[nbrtr]*Ytli[nbrtr]);
    xB=Xli[nbrtr]+(qB-Qli[nbrtr])*utBx;
    yB=Yli[nbrtr]+(qB-Qli[nbrtr])*utBy;
    dxB_ds=utBx;
    dyB_ds=utBy;
    d2xB_ds2=0;
    d2yB_ds2=0;
    }
  // calcul du reste de tES et de t0S
  double dAB=sqrt((xB-xA)*(xB-xA)+(yB-yA)*(yB-yA)),
	 uABx=(xB-xA)/dAB,
	 uABy=(yB-yA)/dAB;
  // Calcul de th,dth_ds et d2th_ds2
  th=atan2(uABy,uABx);
  dth_ds=(uABx*(dyB_ds-dyA_ds)-uABy*(dxB_ds-dxA_ds))/dAB;
  d2th_ds2=(uABx*(d2yB_ds2-d2yA_ds2)-uABy*(d2xB_ds2-d2xA_ds2)
	    -2*(uABx*(dxB_ds-dxA_ds)+uABy*(dyB_ds-dyA_ds))*dth_ds)/dAB;
  }

/****************************************************************************/
