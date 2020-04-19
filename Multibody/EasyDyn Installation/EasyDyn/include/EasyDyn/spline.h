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

#ifndef EASYDYNSPLINE_H
#define EASYDYNSPLINE_H

#include<iostream>
using namespace std;

/****************************************************************************/
/*                         Module MR_SPLIN.H                                */
/*             Fonctions et classe relatives au spline cubique              */
/****************************************************************************/

double N22(double xi);
double N26(double xi);
double N1(double u);
double N2(double u);
double dN1_du(double u);
double dN2_du(double u);
double d2N1_du2(double u);
double d2N2_du2(double u);

class spline
  {
  int nbrtr;       // nombre de tron‡ons du spline;
  double deltaqth, // increment d'abscisse curviligne des deux points
		   // utilises pour determiner la tangente
	 *Xli,     // abscisses des noeuds du spline
	 *Yli,     // ordonnees des noeuds du spline
	 *Xtli,    // abscisses des vecs tangents aux noeuds
	 *Ytli,    // ordonnees des vecs tangents aux noeuds
	 *Qli;     // vec des abscisses curvilignes aux noeuds
  public:
  void set(int nbrtr, double *Xnode,double *Ynode,
		      double *Xtnode, double *Ytnode);
  void Get(int &nbrtr, double *&Xnode,double *&Ynode,
		       double *&Xtnode, double *&Ytnode);
  void ReadSpline(istream &stream);
  void WriteSpline(ostream &stream);
  void InitQli();
  void GetPttg(double qs,
	       double &xpt, double &dx_ds,  double &d2x_ds2,
	       double &ypt, double &dy_ds,  double &d2y_ds2,
	       double &th,  double &dth_ds, double &d2th_ds2);

  void GetPt(double qs,
	     double &xpt, double &dx_ds,  double &d2x_ds2,
	     double &ypt, double &dy_ds,  double &d2y_ds2);
  };

#endif
