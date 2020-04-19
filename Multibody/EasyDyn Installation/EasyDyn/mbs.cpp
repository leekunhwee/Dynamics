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
#include <string.h>
#include <stdio.h>
#include <vector>


#include <EasyDyn/mbs.h>
#include <gsl/gsl_math.h>
//#include <gsl/gsl_errno.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

gsl_matrix *MM=0;
gsl_vector_view fmbsview,qddmbsview;
gsl_permutation *permmbs;

//---------------------------------------------------------------------------

void EasyDynmbsError(const char *texte)
{
 cout << "EasdyDynmbs ERROR: " << texte << "\n";
 cout << "Hit ENTER to leave"; char dummy[3];
 char* useless=fgets(dummy,2,stdin);
 exit(1);
}

//--------------------------------------------------------------------

void InitEasyDynmbs()

{
if (nbrdof==0)
  EasyDynmbsError("nbrdof must be initialized before calling InitMbsSim");
if (application==0)
  EasyDynmbsError("application must be initialized before calling InitMbsSim");
// Memory allocation for variables deriving from sim (q,qd,qdd,f,..)
InitEasyDynsim();
// Memory allocation for dependent variables
if (nbrdep>0)
    {
    p=new double[nbrdep];
    pd=new double[nbrdep];
    pdd=new double[nbrdep];
    }
//
if ((nbrbody==0) && (nbrflexbody==0))
  EasyDynmbsError("InitEasyDynmbs:nbrbody or nbrflexbody must be different from 0");

 int ibody, iflexbody, iframe, jframe;

// Memory allocation for rigid bodies and bodies variables
if (DEBUG) cout << "Initializing bodies" << endl;
if (nbrbody>0)
    {
    body=new structbody[nbrbody];
    for (ibody=0; ibody<nbrbody; ibody++)
        {
        body[ibody].vGpartial=new vec[nbrdof];
        body[ibody].omegapartial=new vec[nbrdof];
        body[ibody].vGrelpartial=new vec[nbrdof];
        body[ibody].omegarelpartial=new vec[nbrdof];
        body[ibody].underdof=new int[nbrdof];
	}
    }

// Memory allocation for flexible bodies
if (DEBUG) cout << "Initializing flexible bodies" << endl;
if (nbrflexbody>0)
    {
    flexbody=new structflexbody[nbrflexbody];
    for (iflexbody=0; iflexbody<nbrflexbody; iflexbody++)
        {
        flexbody[iflexbody].nbrframe=nbrframesperflexbody[iflexbody];
        flexbody[iflexbody].length=0.0;
        flexbody[iflexbody].rho=0.0;
        flexbody[iflexbody].E=0.0;
        flexbody[iflexbody].nu=0.0;
        flexbody[iflexbody].section=0.0;
        flexbody[iflexbody].Iyy=0.0;
        flexbody[iflexbody].Izz=0.0;

        // allocation for mass and stiffness matrices
        for (iframe=0; iframe<flexbody[iflexbody].nbrframe; iframe++)
          {
          flexbody[iflexbody].MTT=new tgen*[flexbody[iflexbody].nbrframe];
          flexbody[iflexbody].MTR=new tgen*[flexbody[iflexbody].nbrframe];
          flexbody[iflexbody].MRT=new tgen*[flexbody[iflexbody].nbrframe];
          flexbody[iflexbody].MRR=new tgen*[flexbody[iflexbody].nbrframe];
          flexbody[iflexbody].KTT=new tgen*[flexbody[iflexbody].nbrframe];
          flexbody[iflexbody].KTR=new tgen*[flexbody[iflexbody].nbrframe];
          flexbody[iflexbody].KRT=new tgen*[flexbody[iflexbody].nbrframe];
          flexbody[iflexbody].KRR=new tgen*[flexbody[iflexbody].nbrframe];
          flexbody[iflexbody].CTT=new tgen*[flexbody[iflexbody].nbrframe];
          flexbody[iflexbody].CTR=new tgen*[flexbody[iflexbody].nbrframe];
          flexbody[iflexbody].CRT=new tgen*[flexbody[iflexbody].nbrframe];
          flexbody[iflexbody].CRR=new tgen*[flexbody[iflexbody].nbrframe];
          for (jframe=0; jframe<flexbody[iflexbody].nbrframe; jframe++)
             {
             flexbody[iflexbody].MTT[jframe]=new tgen[flexbody[iflexbody].nbrframe];
             flexbody[iflexbody].MTR[jframe]=new tgen[flexbody[iflexbody].nbrframe];
             flexbody[iflexbody].MRT[jframe]=new tgen[flexbody[iflexbody].nbrframe];
             flexbody[iflexbody].MRR[jframe]=new tgen[flexbody[iflexbody].nbrframe];
             flexbody[iflexbody].KTT[jframe]=new tgen[flexbody[iflexbody].nbrframe];
             flexbody[iflexbody].KTR[jframe]=new tgen[flexbody[iflexbody].nbrframe];
             flexbody[iflexbody].KRT[jframe]=new tgen[flexbody[iflexbody].nbrframe];
             flexbody[iflexbody].KRR[jframe]=new tgen[flexbody[iflexbody].nbrframe];
             flexbody[iflexbody].CTT[jframe]=new tgen[flexbody[iflexbody].nbrframe];
             flexbody[iflexbody].CTR[jframe]=new tgen[flexbody[iflexbody].nbrframe];
             flexbody[iflexbody].CRT[jframe]=new tgen[flexbody[iflexbody].nbrframe];
             flexbody[iflexbody].CRR[jframe]=new tgen[flexbody[iflexbody].nbrframe];
	     }
	     }
        flexbody[iflexbody].frame=new structframe[flexbody[iflexbody].nbrframe];
        for (iframe=0; iframe<flexbody[iflexbody].nbrframe; iframe++)
          {
          flexbody[iflexbody].frame[iframe].vFpartial=new vec[nbrdof];
          flexbody[iflexbody].frame[iframe].vFrelpartial=new vec[nbrdof];
          flexbody[iflexbody].frame[iframe].omegapartial=new vec[nbrdof];
          flexbody[iflexbody].frame[iframe].omegarelpartial=new vec[nbrdof];
          flexbody[iflexbody].frame[iframe].underdof=new int[nbrdof];
          }
        flexbody[iflexbody].omegacoropartial=new vec[nbrdof];
        flexbody[iflexbody].omegacororelpartial=new vec[nbrdof];
        flexbody[iflexbody].corounderdof=new int[nbrdof];
	}
    }

// Inertia data loading for rigid and flexible bodies through SetInertiaData
if (DEBUG) cout << "Call to SetInertiaData" << endl;
SetInertiaData();

//Rayleigh Damping through the Damping Matrix C = AlphaDamp * M + BetaDamp * K
if(AlphaDamp > 0.0 || BetaDamp > 0.0)
{
    for (iflexbody=0; iflexbody<nbrflexbody; iflexbody++)
    {
        for(int i=0; i< nbrframesperflexbody[iflexbody];i++) //tensor line frame
        {
            for(int j=0; j< nbrframesperflexbody[iflexbody];j++) //tensor column frame
            {
              flexbody[iflexbody].CTT[i][j]=flexbody[iflexbody].MTT[i][j]*AlphaDamp + flexbody[iflexbody].KTT[i][j]*BetaDamp;
              flexbody[iflexbody].CTR[i][j]=flexbody[iflexbody].MTR[i][j]*AlphaDamp + flexbody[iflexbody].KTR[i][j]*BetaDamp;
              flexbody[iflexbody].CRT[i][j]=flexbody[iflexbody].MRT[i][j]*AlphaDamp + flexbody[iflexbody].KRT[i][j]*BetaDamp;
              flexbody[iflexbody].CRR[i][j]=flexbody[iflexbody].MRR[i][j]*AlphaDamp + flexbody[iflexbody].KRR[i][j]*BetaDamp;
            }
        }
    }
}

// Verification of inertia data of rigid bodies
if (DEBUG) cout << "Verifying inertia data of rigid bodies" << endl;
if (nbrbody>0) for (ibody=0; ibody<nbrbody; ibody++)
    {
    if (body[ibody].mass<=0) EasyDynmbsError("Mass null or negative");
    if (body[ibody].PhiG.Ixx<=0) EasyDynmbsError("Ixx null or negative");
    if (body[ibody].PhiG.Iyy<=0) EasyDynmbsError("Iyy null or negative");
    if (body[ibody].PhiG.Izz<=0) EasyDynmbsError("Izz null or negative");
    if (body[ibody].PhiG.Ixx>(body[ibody].PhiG.Iyy+body[ibody].PhiG.Izz))
        EasyDynmbsError("Ixx > Iyy+Izz");
    if (body[ibody].PhiG.Iyy>(body[ibody].PhiG.Ixx+body[ibody].PhiG.Izz))
        EasyDynmbsError("Iyy > Ixx+Izz");
    if (body[ibody].PhiG.Izz>(body[ibody].PhiG.Ixx+body[ibody].PhiG.Iyy))
        EasyDynmbsError("Izz > Ixx+Iyy");
    }

fmbsview=gsl_vector_view_array(f,nbrdof);
qddmbsview=gsl_vector_view_array(qdd,nbrdof);
permmbs = gsl_permutation_alloc(nbrdof);
MM=gsl_matrix_calloc(nbrdof,nbrdof);
}

//---------------------------------------------------------------------------

void EndEasyDynmbs()
{
EndEasyDynsim();
delete body;
gsl_permutation_free(permmbs);
gsl_matrix_free(MM);
}

//--------------------------------------------------------------------

void CreateBeamMassMatrix(int iflexbody, double LL, double AA, double rho, double Iyy, double Izz)
{
     //create the Mass Matrix for the specified flexbody
     //iflexbody : number of the flexbody
     //LL length beam; AA section beam; rho density beam; Iyy Izz moment of inertia

    if(nbrframesperflexbody[iflexbody] > 2) //If several frames per elements, you must divide the Length of beam
    {
        LL = LL/ (double (nbrframesperflexbody[iflexbody] -1));
    }

    double mass = rho * AA * LL;

    flexbody[iflexbody].MTT[0][0].put(mass/3.0,0,0,0,13.0*mass/35.0,0,0,0,13.0*mass/35.0);
    flexbody[iflexbody].MTT[0][1].put(mass/6.0,0,0,0,9.0*mass/70.0,0,0,0,9.0*mass/70.0);
    flexbody[iflexbody].MTT[1][0].put(mass/6.0,0,0,0,9.0*mass/70.0,0,0,0,9.0*mass/70.0);
    flexbody[iflexbody].MTT[1][1].put(mass/3.0,0,0,0,13.0*mass/35.0,0,0,0,13.0*mass/35.0);
    flexbody[iflexbody].MTR[0][0].put(0,0,0,0,0,11.0*mass*LL/210.0,0,-11.0*mass*LL/210.0,0);
    flexbody[iflexbody].MTR[0][1].put(0,0,0,0,0,-13.0*mass*LL/420.0,0,13.0*mass*LL/420.0,0);
    flexbody[iflexbody].MTR[1][0].put(0,0,0,0,0,13.0*mass*LL/420.0,0,-13.0*mass*LL/420.0,0);
    flexbody[iflexbody].MTR[1][1].put(0,0,0,0,0,-11.0*mass*LL/210.0,0,11.0*mass*LL/210.0,0);
    flexbody[iflexbody].MRT[0][0].put(0,0,0,0,0,-11.0*mass*LL/210.0,0,11.0*mass*LL/210.0,0);
    flexbody[iflexbody].MRT[0][1].put(0,0,0,0,0,-13.0*mass*LL/420.0,0,13.0*mass*LL/420.0,0);
    flexbody[iflexbody].MRT[1][0].put(0,0,0,0,0,13.0*mass*LL/420.0,0,-13.0*mass*LL/420.0,0);
    flexbody[iflexbody].MRT[1][1].put(0,0,0,0,0,11.0*mass*LL/210.0,0,-11.0*mass*LL/210.0,0);
    flexbody[iflexbody].MRR[0][0].put(rho*LL*(Iyy+Izz)/3.0,0,0,0,mass*LL*LL/105.0,0,0,0,mass*LL*LL/105.0);
    flexbody[iflexbody].MRR[0][1].put(rho*LL*(Iyy+Izz)/6.0,0,0,0,-mass*LL*LL/140.0,0,0,0,-mass*LL*LL/140.0);
    flexbody[iflexbody].MRR[1][0].put(rho*LL*(Iyy+Izz)/6.0,0,0,0,-mass*LL*LL/140.0,0,0,0,-mass*LL*LL/140.0);
    flexbody[iflexbody].MRR[1][1].put(rho*LL*(Iyy+Izz)/3.0,0,0,0,mass*LL*LL/105.0,0,0,0,mass*LL*LL/105.0);
}

//--------------------------------------------------------------------

gsl_matrix *mbeam3d(int iflexbody, int initframe, int endframe, double AA, double rho, double Ixx, double x3, double y3, double z3)
{
     // Create the mass matrix of a 3D beam with its rotation (same as in EasyFEM)
     //iflexbody: number of the flexbody of which the elementary mass matrix is being computed
     //initframe: number of the starting frame of the flexbody iflexbody
     //endframe: number of the ending frame of the flexbody iflexbody
     //AA: section m2
     //rho: density kg/m3
     //Ix: moment of inertia x
     //x3: fictive point to uniquely orientate the beam in space

     ComputeMotion();
     double LL = 0.0; //length of the beam

     //Coordinates of first frame
     double x1 =flexbody[iflexbody].frame[initframe].T0F.e.x;
     double y1 =flexbody[iflexbody].frame[initframe].T0F.e.y;
     double z1 =flexbody[iflexbody].frame[initframe].T0F.e.z;

     //Coordinates of second frame
     double x2 =flexbody[iflexbody].frame[endframe].T0F.e.x;
     double y2 =flexbody[iflexbody].frame[endframe].T0F.e.y;
     double z2 =flexbody[iflexbody].frame[endframe].T0F.e.z;

     LL = sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)+(z2-z1)*(z2-z1)); //current length of the beam

     gsl_vector *d2 = gsl_vector_alloc(3);
     gsl_vector *d3 = gsl_vector_alloc(3);
     gsl_vector *ex = gsl_vector_alloc(3);
     gsl_vector *ey = gsl_vector_alloc(3);
     gsl_vector *ez = gsl_vector_alloc(3);
     gsl_matrix *T_Rot = gsl_matrix_alloc(12,12);
     gsl_matrix *MassMat = gsl_matrix_alloc(12,12);
     gsl_matrix *Mass_Rot = gsl_matrix_alloc(12,12);
     gsl_matrix *Mass_Matrix = gsl_matrix_alloc(12,12);

     gsl_vector_set(d2,0,x2-x1);
     gsl_vector_set(d2,1,y2-y1);
     gsl_vector_set(d2,2,z2-z1);

     gsl_vector_set(ex,0,x2-x1);
     gsl_vector_set(ex,1,y2-y1);
     gsl_vector_set(ex,2,z2-z1);

     gsl_vector_scale (ex, 1.0/LL);

     gsl_vector_set(d3,0,x3-x1);
     gsl_vector_set(d3,1,y3-y1);
     gsl_vector_set(d3,2,z3-z1);

     //cross product d3^d2
     gsl_vector_set(ey,0,(gsl_vector_get(d3,1)*gsl_vector_get(d2,2)) - (gsl_vector_get(d2,1)*gsl_vector_get(d3,2)) );
     gsl_vector_set(ey,1,(gsl_vector_get(d2,0)*gsl_vector_get(d3,2)) - (gsl_vector_get(d3,0)*gsl_vector_get(d2,2)) );
     gsl_vector_set(ey,2,(gsl_vector_get(d3,0)*gsl_vector_get(d2,1)) - (gsl_vector_get(d2,0)*gsl_vector_get(d3,1)) );

     gsl_vector_scale(ey,1.0/(LL*gsl_blas_dnrm2(d3)));

     //cross product ex^ey
     gsl_vector_set(ez,0,(gsl_vector_get(ex,1)*gsl_vector_get(ey,2)) - (gsl_vector_get(ey,1)*gsl_vector_get(ex,2)) );
     gsl_vector_set(ez,1,(gsl_vector_get(ey,0)*gsl_vector_get(ex,2)) - (gsl_vector_get(ex,0)*gsl_vector_get(ey,2)) );
     gsl_vector_set(ez,2,(gsl_vector_get(ex,0)*gsl_vector_get(ey,1)) - (gsl_vector_get(ey,0)*gsl_vector_get(ex,1)) );

     //Transformation matrix composed of rotation matrices (Block Diagonal Matrix)
     gsl_matrix_set_zero(T_Rot);
     gsl_matrix_set(T_Rot,0,0,gsl_vector_get(ex,0));
     gsl_matrix_set(T_Rot,0,1,gsl_vector_get(ex,1));
     gsl_matrix_set(T_Rot,0,2,gsl_vector_get(ex,2));
     gsl_matrix_set(T_Rot,1,0,gsl_vector_get(ey,0));
     gsl_matrix_set(T_Rot,1,1,gsl_vector_get(ey,1));
     gsl_matrix_set(T_Rot,1,2,gsl_vector_get(ey,2));
     gsl_matrix_set(T_Rot,2,0,gsl_vector_get(ez,0));
     gsl_matrix_set(T_Rot,2,1,gsl_vector_get(ez,1));
     gsl_matrix_set(T_Rot,2,2,gsl_vector_get(ez,2));

     gsl_matrix_set(T_Rot,3,3,gsl_vector_get(ex,0));
     gsl_matrix_set(T_Rot,3,4,gsl_vector_get(ex,1));
     gsl_matrix_set(T_Rot,3,5,gsl_vector_get(ex,2));
     gsl_matrix_set(T_Rot,4,3,gsl_vector_get(ey,0));
     gsl_matrix_set(T_Rot,4,4,gsl_vector_get(ey,1));
     gsl_matrix_set(T_Rot,4,5,gsl_vector_get(ey,2));
     gsl_matrix_set(T_Rot,5,3,gsl_vector_get(ez,0));
     gsl_matrix_set(T_Rot,5,4,gsl_vector_get(ez,1));
     gsl_matrix_set(T_Rot,5,5,gsl_vector_get(ez,2));

     gsl_matrix_set(T_Rot,6,6,gsl_vector_get(ex,0));
     gsl_matrix_set(T_Rot,6,7,gsl_vector_get(ex,1));
     gsl_matrix_set(T_Rot,6,8,gsl_vector_get(ex,2));
     gsl_matrix_set(T_Rot,7,6,gsl_vector_get(ey,0));
     gsl_matrix_set(T_Rot,7,7,gsl_vector_get(ey,1));
     gsl_matrix_set(T_Rot,7,8,gsl_vector_get(ey,2));
     gsl_matrix_set(T_Rot,8,6,gsl_vector_get(ez,0));
     gsl_matrix_set(T_Rot,8,7,gsl_vector_get(ez,1));
     gsl_matrix_set(T_Rot,8,8,gsl_vector_get(ez,2));

     gsl_matrix_set(T_Rot,9,9,gsl_vector_get(ex,0));
     gsl_matrix_set(T_Rot,9,10,gsl_vector_get(ex,1));
     gsl_matrix_set(T_Rot,9,11,gsl_vector_get(ex,2));
     gsl_matrix_set(T_Rot,10,9,gsl_vector_get(ey,0));
     gsl_matrix_set(T_Rot,10,10,gsl_vector_get(ey,1));
     gsl_matrix_set(T_Rot,10,11,gsl_vector_get(ey,2));
     gsl_matrix_set(T_Rot,11,9,gsl_vector_get(ez,0));
     gsl_matrix_set(T_Rot,11,10,gsl_vector_get(ez,1));
     gsl_matrix_set(T_Rot,11,11,gsl_vector_get(ez,2));

     //Set the Mass Matrix
     double mass = rho * AA * LL;
     double r2 = Ixx/AA;
     double m1 = mass/3.0;
     double m2 = mass/6.0;
     double m3 = 13.0*mass/35.0;
     double m4 = 11.0*mass*LL/210.0;
     double m5 = 9.0*mass/70.0;
     double m6 = -13.0*mass*LL/420.0;
     double m7 = mass*r2/3.0;
     double m8 = mass*r2/6.0;
     double m9 = mass*LL*LL/105.0;
     double m10 = -mass*LL*LL/420.0;
     double m11 = -mass*LL*LL/140.0;

     gsl_matrix_set_zero(MassMat);
     gsl_matrix_set(MassMat,0,0,m1);
     gsl_matrix_set(MassMat,6,0,m2);
     gsl_matrix_set(MassMat,1,1,m3);
     gsl_matrix_set(MassMat,5,1,m4);
     gsl_matrix_set(MassMat,7,1,m5);
     gsl_matrix_set(MassMat,11,1,m6);
     gsl_matrix_set(MassMat,2,2,m3);
     gsl_matrix_set(MassMat,4,2,-m4);
     gsl_matrix_set(MassMat,8,2,m5);
     gsl_matrix_set(MassMat,10,2,-m6);
     gsl_matrix_set(MassMat,3,3,m7);
     gsl_matrix_set(MassMat,9,3,m8);
     gsl_matrix_set(MassMat,2,4,-m4);
     gsl_matrix_set(MassMat,4,4,m9);
     gsl_matrix_set(MassMat,8,4,m6);
     gsl_matrix_set(MassMat,10,4,m10);
     gsl_matrix_set(MassMat,1,5,m4);
     gsl_matrix_set(MassMat,5,5,m9);
     gsl_matrix_set(MassMat,7,5,-m6);
     gsl_matrix_set(MassMat,11,5,m11);
     gsl_matrix_set(MassMat,0,6,m2);
     gsl_matrix_set(MassMat,6,6,m1);
     gsl_matrix_set(MassMat,1,7,m5);
     gsl_matrix_set(MassMat,5,7,-m6);
     gsl_matrix_set(MassMat,7,7,m3);
     gsl_matrix_set(MassMat,11,7,-m4);
     gsl_matrix_set(MassMat,2,8,m5);
     gsl_matrix_set(MassMat,4,8,m6);
     gsl_matrix_set(MassMat,8,8,m3);
     gsl_matrix_set(MassMat,10,8,m4);
     gsl_matrix_set(MassMat,3,9,m8);
     gsl_matrix_set(MassMat,9,9,m7);
     gsl_matrix_set(MassMat,2,10,-m6);
     gsl_matrix_set(MassMat,4,10,m10);
     gsl_matrix_set(MassMat,8,10,m4);
     gsl_matrix_set(MassMat,10,10,m9);
     gsl_matrix_set(MassMat,1,11,m6);
     gsl_matrix_set(MassMat,5,11,m11);
     gsl_matrix_set(MassMat,7,11,-m4);
     gsl_matrix_set(MassMat,11,11,m9);

     gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, MassMat, T_Rot,0.0, Mass_Rot);
     gsl_blas_dgemm (CblasTrans, CblasNoTrans,1.0, T_Rot, Mass_Rot, 0.0, Mass_Matrix);

     gsl_matrix_free(MassMat);
     gsl_matrix_free(Mass_Rot);
     gsl_vector_free(d2);
     gsl_vector_free(d3);
     gsl_vector_free(ex);
     gsl_vector_free(ey);
     gsl_vector_free(ez);

     return Mass_Matrix;
}

//--------------------------------------------------------------------

gsl_matrix *kbeam3d(int iflexbody, int initframe, int endframe, double E, double nu, double AA, double Ixx, double Iyy, double Izz, double x3, double y3, double z3)
{
     // Create the stiffness matrix of a 3D beam with its rotation (same as in EasyFEM)
     //iflexbody: number of the flexbody of which the elementary mass matrix is being computed
     //initframe: number of the starting frame of the flexbody iflexbody
     //endframe: number of the ending frame of the flexbody iflexbody
     //E: Young's Modulus
     //nu: Poisson's Modulus
     //AA: section m2
     //rho: density kg/m3
     //Ix: moment of inertia x
     //x3: fictive point to uniquely orientate the beam in space

     ComputeMotion();
     double LL = 0.0; //length of the beam
     double G = E/(2.0*(1.0+nu));

     //Coordinates of first frame
     double x1 =flexbody[iflexbody].frame[initframe].T0F.e.x;
     double y1 =flexbody[iflexbody].frame[initframe].T0F.e.y;
     double z1 =flexbody[iflexbody].frame[initframe].T0F.e.z;

     //Coordinates of second frame
     double x2 =flexbody[iflexbody].frame[endframe].T0F.e.x;
     double y2 =flexbody[iflexbody].frame[endframe].T0F.e.y;
     double z2 =flexbody[iflexbody].frame[endframe].T0F.e.z;

     LL = sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)+(z2-z1)*(z2-z1)); //current length of the beam

     gsl_vector *d2 = gsl_vector_alloc(3);
     gsl_vector *d3 = gsl_vector_alloc(3);
     gsl_vector *ex = gsl_vector_alloc(3);
     gsl_vector *ey = gsl_vector_alloc(3);
     gsl_vector *ez = gsl_vector_alloc(3);
     gsl_matrix *T_Rot = gsl_matrix_alloc(12,12);
     gsl_matrix *StiffMat = gsl_matrix_alloc(12,12);
     gsl_matrix *Stiff_Rot = gsl_matrix_alloc(12,12);
     gsl_matrix *Stiff_Matrix = gsl_matrix_alloc(12,12);

     gsl_vector_set(d2,0,x2-x1);
     gsl_vector_set(d2,1,y2-y1);
     gsl_vector_set(d2,2,z2-z1);

     gsl_vector_set(ex,0,x2-x1);
     gsl_vector_set(ex,1,y2-y1);
     gsl_vector_set(ex,2,z2-z1);

     gsl_vector_scale (ex, 1.0/LL);

     gsl_vector_set(d3,0,x3-x1);
     gsl_vector_set(d3,1,y3-y1);
     gsl_vector_set(d3,2,z3-z1);

     //cross product d3^d2
     gsl_vector_set(ey,0,(gsl_vector_get(d3,1)*gsl_vector_get(d2,2)) - (gsl_vector_get(d2,1)*gsl_vector_get(d3,2)) );
     gsl_vector_set(ey,1,(gsl_vector_get(d2,0)*gsl_vector_get(d3,2)) - (gsl_vector_get(d3,0)*gsl_vector_get(d2,2)) );
     gsl_vector_set(ey,2,(gsl_vector_get(d3,0)*gsl_vector_get(d2,1)) - (gsl_vector_get(d2,0)*gsl_vector_get(d3,1)) );

     gsl_vector_scale(ey,1.0/(LL*gsl_blas_dnrm2(d3)));

     //cross product ex^ey
     gsl_vector_set(ez,0,(gsl_vector_get(ex,1)*gsl_vector_get(ey,2)) - (gsl_vector_get(ey,1)*gsl_vector_get(ex,2)) );
     gsl_vector_set(ez,1,(gsl_vector_get(ey,0)*gsl_vector_get(ex,2)) - (gsl_vector_get(ex,0)*gsl_vector_get(ey,2)) );
     gsl_vector_set(ez,2,(gsl_vector_get(ex,0)*gsl_vector_get(ey,1)) - (gsl_vector_get(ey,0)*gsl_vector_get(ex,1)) );

     //Transformation matrix composed of rotation matrices (Block Diagonal Matrix)
     gsl_matrix_set_zero(T_Rot);
     gsl_matrix_set(T_Rot,0,0,gsl_vector_get(ex,0));
     gsl_matrix_set(T_Rot,0,1,gsl_vector_get(ex,1));
     gsl_matrix_set(T_Rot,0,2,gsl_vector_get(ex,2));
     gsl_matrix_set(T_Rot,1,0,gsl_vector_get(ey,0));
     gsl_matrix_set(T_Rot,1,1,gsl_vector_get(ey,1));
     gsl_matrix_set(T_Rot,1,2,gsl_vector_get(ey,2));
     gsl_matrix_set(T_Rot,2,0,gsl_vector_get(ez,0));
     gsl_matrix_set(T_Rot,2,1,gsl_vector_get(ez,1));
     gsl_matrix_set(T_Rot,2,2,gsl_vector_get(ez,2));

     gsl_matrix_set(T_Rot,3,3,gsl_vector_get(ex,0));
     gsl_matrix_set(T_Rot,3,4,gsl_vector_get(ex,1));
     gsl_matrix_set(T_Rot,3,5,gsl_vector_get(ex,2));
     gsl_matrix_set(T_Rot,4,3,gsl_vector_get(ey,0));
     gsl_matrix_set(T_Rot,4,4,gsl_vector_get(ey,1));
     gsl_matrix_set(T_Rot,4,5,gsl_vector_get(ey,2));
     gsl_matrix_set(T_Rot,5,3,gsl_vector_get(ez,0));
     gsl_matrix_set(T_Rot,5,4,gsl_vector_get(ez,1));
     gsl_matrix_set(T_Rot,5,5,gsl_vector_get(ez,2));

     gsl_matrix_set(T_Rot,6,6,gsl_vector_get(ex,0));
     gsl_matrix_set(T_Rot,6,7,gsl_vector_get(ex,1));
     gsl_matrix_set(T_Rot,6,8,gsl_vector_get(ex,2));
     gsl_matrix_set(T_Rot,7,6,gsl_vector_get(ey,0));
     gsl_matrix_set(T_Rot,7,7,gsl_vector_get(ey,1));
     gsl_matrix_set(T_Rot,7,8,gsl_vector_get(ey,2));
     gsl_matrix_set(T_Rot,8,6,gsl_vector_get(ez,0));
     gsl_matrix_set(T_Rot,8,7,gsl_vector_get(ez,1));
     gsl_matrix_set(T_Rot,8,8,gsl_vector_get(ez,2));

     gsl_matrix_set(T_Rot,9,9,gsl_vector_get(ex,0));
     gsl_matrix_set(T_Rot,9,10,gsl_vector_get(ex,1));
     gsl_matrix_set(T_Rot,9,11,gsl_vector_get(ex,2));
     gsl_matrix_set(T_Rot,10,9,gsl_vector_get(ey,0));
     gsl_matrix_set(T_Rot,10,10,gsl_vector_get(ey,1));
     gsl_matrix_set(T_Rot,10,11,gsl_vector_get(ey,2));
     gsl_matrix_set(T_Rot,11,9,gsl_vector_get(ez,0));
     gsl_matrix_set(T_Rot,11,10,gsl_vector_get(ez,1));
     gsl_matrix_set(T_Rot,11,11,gsl_vector_get(ez,2));

     //Set the Stiffness Matrix
     double k1 = E*AA/LL;
     double k2 = (12.0*E*Izz)/(LL*LL*LL);
     double k3 = (6.0*E*Izz)/(LL*LL);
     double k4 = (12.0*E*Iyy)/(LL*LL*LL);
     double k5 = (-6.0*E*Iyy)/(LL*LL);
     double k6 = (G*Ixx)/LL;
     double k7 = (4.0*E*Iyy)/LL;
     double k8 = (6.0*E*Iyy)/(LL*LL);
     double k9 = (2.0*E*Iyy)/LL;
     double k10 = (4.0*E*Izz)/(LL);
     double k11 = (-6*E*Izz)/(LL*LL);
     double k12 = (2*E*Izz)/LL;

     gsl_matrix_set_zero(StiffMat);
     gsl_matrix_set(StiffMat,0,0,k1);
     gsl_matrix_set(StiffMat,6,0,-k1);
     gsl_matrix_set(StiffMat,1,1,k2);
     gsl_matrix_set(StiffMat,5,1,k3);
     gsl_matrix_set(StiffMat,7,1,-k2);
     gsl_matrix_set(StiffMat,11,1,k3);
     gsl_matrix_set(StiffMat,2,2,k4);
     gsl_matrix_set(StiffMat,4,2,k5);
     gsl_matrix_set(StiffMat,8,2,-k4);
     gsl_matrix_set(StiffMat,10,2,k5);
     gsl_matrix_set(StiffMat,3,3,k6);
     gsl_matrix_set(StiffMat,9,3,-k6);
     gsl_matrix_set(StiffMat,2,4,k5);
     gsl_matrix_set(StiffMat,4,4,k7);
     gsl_matrix_set(StiffMat,8,4,k8);
     gsl_matrix_set(StiffMat,10,4,k9);
     gsl_matrix_set(StiffMat,1,5,k3);
     gsl_matrix_set(StiffMat,5,5,k10);
     gsl_matrix_set(StiffMat,7,5,k11);
     gsl_matrix_set(StiffMat,11,5,k12);
     gsl_matrix_set(StiffMat,0,6,-k1);
     gsl_matrix_set(StiffMat,6,6,k1);
     gsl_matrix_set(StiffMat,1,7,-k2);
     gsl_matrix_set(StiffMat,5,7,k11);
     gsl_matrix_set(StiffMat,7,7,k2);
     gsl_matrix_set(StiffMat,11,7,-k3);
     gsl_matrix_set(StiffMat,2,8,-k4);
     gsl_matrix_set(StiffMat,4,8,k8);
     gsl_matrix_set(StiffMat,8,8,k4);
     gsl_matrix_set(StiffMat,10,8,-k5);
     gsl_matrix_set(StiffMat,3,9,-k6);
     gsl_matrix_set(StiffMat,9,9,k6);
     gsl_matrix_set(StiffMat,2,10,k5);
     gsl_matrix_set(StiffMat,4,10,k9);
     gsl_matrix_set(StiffMat,8,10,-k5);
     gsl_matrix_set(StiffMat,10,10,k7);
     gsl_matrix_set(StiffMat,1,11,k3);
     gsl_matrix_set(StiffMat,5,11,k12);
     gsl_matrix_set(StiffMat,7,11,-k3);
     gsl_matrix_set(StiffMat,11,11,k10);

     gsl_blas_dgemm (CblasNoTrans, CblasNoTrans,1.0, StiffMat, T_Rot,0.0, Stiff_Rot);
     gsl_blas_dgemm (CblasTrans, CblasNoTrans,1.0, T_Rot, Stiff_Rot, 0.0, Stiff_Matrix);

     gsl_matrix_free(StiffMat);
     gsl_matrix_free(Stiff_Rot);
     gsl_vector_free(d2);
     gsl_vector_free(d3);
     gsl_vector_free(ex);
     gsl_vector_free(ey);
     gsl_vector_free(ez);

     return Stiff_Matrix;
}

//--------------------------------------------------------------------

void connect(int iflexbody, gsl_matrix *MassElement, gsl_matrix *StiffElement, int q0, int q1, int q2, int q3, int q4, int q5, int q6, int q7, int q8, int q9, int q10, int q11)
{
    // Assemble the mbeam3d and kneam3d matrices
   static int flexbodynumber = 0;

   static gsl_matrix *GlobalMassMatrix = gsl_matrix_alloc(nbrframesperflexbody[iflexbody]*6,nbrframesperflexbody[iflexbody]*6);
   static gsl_matrix *GlobalStiffMatrix = gsl_matrix_alloc(nbrframesperflexbody[iflexbody]*6,nbrframesperflexbody[iflexbody]*6);

   if (iflexbody != flexbodynumber)  //Reinitialize if it is another flexbody
   {
      gsl_matrix_set_zero(GlobalMassMatrix);
      gsl_matrix_set_zero(GlobalStiffMatrix);
   }
   flexbodynumber=iflexbody;

   gsl_vector *connect = gsl_vector_alloc(12);  //(nbrDOF,nbrBeams)
   gsl_vector_set(connect,0,q0);
   gsl_vector_set(connect,1,q1);
   gsl_vector_set(connect,2,q2);
   gsl_vector_set(connect,3,q3);
   gsl_vector_set(connect,4,q4);
   gsl_vector_set(connect,5,q5);
   gsl_vector_set(connect,6,q6);
   gsl_vector_set(connect,7,q7);
   gsl_vector_set(connect,8,q8);
   gsl_vector_set(connect,9,q9);
   gsl_vector_set(connect,10,q10);
   gsl_vector_set(connect,11,q11);

   int ii, jj;
   for(int i=0 ; i<12 ; i++)
   {
        for(int j=0 ; j<12 ; j++)
        {
              ii=(int) gsl_vector_get(connect,i);
              jj=(int) gsl_vector_get(connect,j);

              if ((ii >= 0) && (jj >= 0))
              {
                gsl_matrix_set(GlobalMassMatrix,ii,jj,gsl_matrix_get(GlobalMassMatrix,ii,jj)+gsl_matrix_get(MassElement,i,j));
              }
        }
   }
   for(int i=0 ; i<12 ; i++)
   {
        for(int j=0 ; j<12 ; j++)
        {
              ii=(int) gsl_vector_get(connect,i);
              jj=(int) gsl_vector_get(connect,j);

              if ((ii >= 0) && (jj >= 0))
              {
                gsl_matrix_set(GlobalStiffMatrix,ii,jj,gsl_matrix_get(GlobalStiffMatrix,ii,jj)+gsl_matrix_get(StiffElement,i,j));
              }
        }
   }

   for(int i=0; i<nbrframesperflexbody[iflexbody];i++)  //Transfer the global  mass matrix into tensors
   {
        for(int j=0; j<nbrframesperflexbody[iflexbody];j++)
        {
           //Mass Matrix
           flexbody[iflexbody].MTT[i][j].put(gsl_matrix_get(GlobalMassMatrix,0+(i*6),0+(j*6)),gsl_matrix_get(GlobalMassMatrix,0+(i*6),1+(j*6)),gsl_matrix_get(GlobalMassMatrix,0+(i*6),2+(j*6)),gsl_matrix_get(GlobalMassMatrix,1+(i*6),0+(j*6)),gsl_matrix_get(GlobalMassMatrix,1+(i*6),1+(j*6)),gsl_matrix_get(GlobalMassMatrix,1+(i*6),2+(j*6)),gsl_matrix_get(GlobalMassMatrix,2+(i*6),0+(j*6)),gsl_matrix_get(GlobalMassMatrix,2+(i*6),1+(j*6)),gsl_matrix_get(GlobalMassMatrix,2+(i*6),2+(j*6)));
           flexbody[iflexbody].MTR[i][j].put(gsl_matrix_get(GlobalMassMatrix,0+(i*6),3+(j*6)),gsl_matrix_get(GlobalMassMatrix,0+(i*6),4+(j*6)),gsl_matrix_get(GlobalMassMatrix,0+(i*6),5+(j*6)),gsl_matrix_get(GlobalMassMatrix,1+(i*6),3+(j*6)),gsl_matrix_get(GlobalMassMatrix,1+(i*6),4+(j*6)),gsl_matrix_get(GlobalMassMatrix,1+(i*6),5+(j*6)),gsl_matrix_get(GlobalMassMatrix,2+(i*6),3+(j*6)),gsl_matrix_get(GlobalMassMatrix,2+(i*6),4+(j*6)),gsl_matrix_get(GlobalMassMatrix,2+(i*6),5+(j*6)));
           flexbody[iflexbody].MRT[i][j].put(gsl_matrix_get(GlobalMassMatrix,3+(i*6),0+(j*6)),gsl_matrix_get(GlobalMassMatrix,3+(i*6),1+(j*6)),gsl_matrix_get(GlobalMassMatrix,3+(i*6),2+(j*6)),gsl_matrix_get(GlobalMassMatrix,4+(i*6),0+(j*6)),gsl_matrix_get(GlobalMassMatrix,4+(i*6),1+(j*6)),gsl_matrix_get(GlobalMassMatrix,4+(i*6),2+(j*6)),gsl_matrix_get(GlobalMassMatrix,5+(i*6),0+(j*6)),gsl_matrix_get(GlobalMassMatrix,5+(i*6),1+(j*6)),gsl_matrix_get(GlobalMassMatrix,5+(i*6),2+(j*6)));
           flexbody[iflexbody].MRR[i][j].put(gsl_matrix_get(GlobalMassMatrix,3+(i*6),3+(j*6)),gsl_matrix_get(GlobalMassMatrix,3+(i*6),4+(j*6)),gsl_matrix_get(GlobalMassMatrix,3+(i*6),5+(j*6)),gsl_matrix_get(GlobalMassMatrix,4+(i*6),3+(j*6)),gsl_matrix_get(GlobalMassMatrix,4+(i*6),4+(j*6)),gsl_matrix_get(GlobalMassMatrix,4+(i*6),5+(j*6)),gsl_matrix_get(GlobalMassMatrix,5+(i*6),3+(j*6)),gsl_matrix_get(GlobalMassMatrix,5+(i*6),4+(j*6)),gsl_matrix_get(GlobalMassMatrix,5+(i*6),5+(j*6)));

        }
   }

   for(int i=0; i<nbrframesperflexbody[iflexbody];i++)  //Transfer the global stiffness matrix into tensors
   {
        for(int j=0; j<nbrframesperflexbody[iflexbody];j++)
        {
           //Stiffness Matrix
           flexbody[iflexbody].KTT[i][j].put(gsl_matrix_get(GlobalStiffMatrix,0+(i*6),0+(j*6)),gsl_matrix_get(GlobalStiffMatrix,0+(i*6),1+(j*6)),gsl_matrix_get(GlobalStiffMatrix,0+(i*6),2+(j*6)),gsl_matrix_get(GlobalStiffMatrix,1+(i*6),0+(j*6)),gsl_matrix_get(GlobalStiffMatrix,1+(i*6),1+(j*6)),gsl_matrix_get(GlobalStiffMatrix,1+(i*6),2+(j*6)),gsl_matrix_get(GlobalStiffMatrix,2+(i*6),0+(j*6)),gsl_matrix_get(GlobalStiffMatrix,2+(i*6),1+(j*6)),gsl_matrix_get(GlobalStiffMatrix,2+(i*6),2+(j*6)));
           flexbody[iflexbody].KTR[i][j].put(gsl_matrix_get(GlobalStiffMatrix,0+(i*6),3+(j*6)),gsl_matrix_get(GlobalStiffMatrix,0+(i*6),4+(j*6)),gsl_matrix_get(GlobalStiffMatrix,0+(i*6),5+(j*6)),gsl_matrix_get(GlobalStiffMatrix,1+(i*6),3+(j*6)),gsl_matrix_get(GlobalStiffMatrix,1+(i*6),4+(j*6)),gsl_matrix_get(GlobalStiffMatrix,1+(i*6),5+(j*6)),gsl_matrix_get(GlobalStiffMatrix,2+(i*6),3+(j*6)),gsl_matrix_get(GlobalStiffMatrix,2+(i*6),4+(j*6)),gsl_matrix_get(GlobalStiffMatrix,2+(i*6),5+(j*6)));
           flexbody[iflexbody].KRT[i][j].put(gsl_matrix_get(GlobalStiffMatrix,3+(i*6),0+(j*6)),gsl_matrix_get(GlobalStiffMatrix,3+(i*6),1+(j*6)),gsl_matrix_get(GlobalStiffMatrix,3+(i*6),2+(j*6)),gsl_matrix_get(GlobalStiffMatrix,4+(i*6),0+(j*6)),gsl_matrix_get(GlobalStiffMatrix,4+(i*6),1+(j*6)),gsl_matrix_get(GlobalStiffMatrix,4+(i*6),2+(j*6)),gsl_matrix_get(GlobalStiffMatrix,5+(i*6),0+(j*6)),gsl_matrix_get(GlobalStiffMatrix,5+(i*6),1+(j*6)),gsl_matrix_get(GlobalStiffMatrix,5+(i*6),2+(j*6)));
           flexbody[iflexbody].KRR[i][j].put(gsl_matrix_get(GlobalStiffMatrix,3+(i*6),3+(j*6)),gsl_matrix_get(GlobalStiffMatrix,3+(i*6),4+(j*6)),gsl_matrix_get(GlobalStiffMatrix,3+(i*6),5+(j*6)),gsl_matrix_get(GlobalStiffMatrix,4+(i*6),3+(j*6)),gsl_matrix_get(GlobalStiffMatrix,4+(i*6),4+(j*6)),gsl_matrix_get(GlobalStiffMatrix,4+(i*6),5+(j*6)),gsl_matrix_get(GlobalStiffMatrix,5+(i*6),3+(j*6)),gsl_matrix_get(GlobalStiffMatrix,5+(i*6),4+(j*6)),gsl_matrix_get(GlobalStiffMatrix,5+(i*6),5+(j*6)));

        }
   }
}

//--------------------------------------------------------------------

void CreateBeamStiffnessMatrix(int iflexbody, double LL, double AA, double EE, double Iyy, double Izz, double nu)
{
     //create the Stiffness Matrix for the specified flexbody
     //iflexbody : number of the flexbody
     //LL length beam; AA section beam; EE Young's modulus beam; Iyy Izz moment of inertia; Ip polar inertia; nu Poisson ratio

    if(nbrframesperflexbody[iflexbody] > 2) //If several frames per elements, you must divide the Length of beam
    {
        LL = LL/ (double(nbrframesperflexbody[iflexbody] -1));
    }

     double G = EE / (2.0*(1.0+nu)) ;//Shear modulus
     double Ip = Iyy + Izz; //Polar inertia

    flexbody[iflexbody].KTT[0][0].put((EE*AA)/LL,0,0,0,(12.0*EE*Izz)/(LL*LL*LL),0,0,0,(12.0*EE*Iyy)/(LL*LL*LL));
    flexbody[iflexbody].KTT[0][1].put((-EE*AA)/LL,0,0,0,(-12.0*EE*Izz)/(LL*LL*LL),0,0,0,(-12.0*EE*Iyy)/(LL*LL*LL));
    flexbody[iflexbody].KTT[1][0].put((-EE*AA)/LL,0,0,0,(-12.0*EE*Izz)/(LL*LL*LL),0,0,0,(-12.0*EE*Iyy)/(LL*LL*LL));
    flexbody[iflexbody].KTT[1][1].put((EE*AA)/LL,0,0,0,(12.0*EE*Izz)/(LL*LL*LL),0,0,0,(12.0*EE*Iyy)/(LL*LL*LL));
    flexbody[iflexbody].KTR[0][0].put(0,0,0,0,0,(6.0*EE*Izz)/(LL*LL),0,(-6.0*EE*Iyy)/(LL*LL),0);
    flexbody[iflexbody].KTR[0][1].put(0,0,0,0,0,(6.0*EE*Izz)/(LL*LL),0,(-6.0*EE*Iyy)/(LL*LL),0);
    flexbody[iflexbody].KTR[1][0].put(0,0,0,0,0,(-6.0*EE*Izz)/(LL*LL),0,(6.0*EE*Iyy)/(LL*LL),0);
    flexbody[iflexbody].KTR[1][1].put(0,0,0,0,0,(-6.0*EE*Izz)/(LL*LL),0,(6.0*EE*Iyy)/(LL*LL),0);
    flexbody[iflexbody].KRT[0][0].put(0,0,0,0,0,(-6.0*EE*Iyy)/(LL*LL),0,(6.0*EE*Izz)/(LL*LL),0);
    flexbody[iflexbody].KRT[0][1].put(0,0,0,0,0,(6.0*EE*Iyy)/(LL*LL),0,(-6.0*EE*Izz)/(LL*LL),0);
    flexbody[iflexbody].KRT[1][0].put(0,0,0,0,0,(-6.0*EE*Iyy)/(LL*LL),0,(6.0*EE*Izz)/(LL*LL),0);
    flexbody[iflexbody].KRT[1][1].put(0,0,0,0,0,(6.0*EE*Iyy)/(LL*LL),0,(-6.0*EE*Izz)/(LL*LL),0);
    flexbody[iflexbody].KRR[0][0].put((G*Ip)/(LL),0,0,0,(4.0*EE*Iyy)/(LL),0,0,0,(4.0*EE*Izz)/(LL));
    flexbody[iflexbody].KRR[0][1].put((-G*Ip)/(LL),0,0,0,(2.0*EE*Iyy)/(LL),0,0,0,(2.0*EE*Izz)/(LL));
    flexbody[iflexbody].KRR[1][0].put((-G*Ip)/(LL),0,0,0,(2.0*EE*Iyy)/(LL),0,0,0,(2.0*EE*Izz)/(LL));
    flexbody[iflexbody].KRR[1][1].put((G*Ip)/(LL),0,0,0,(4.0*EE*Iyy)/(LL),0,0,0,(4.0*EE*Izz)/(LL));
}

//--------------------------------------------------------------------

void AssemblingBeamMK()
{
//Procedure to assemble mass and stiffness matrices for a beam element along the X axis
   if(nbrflexbody > 0)
   {
         for (int iflexbody=0; iflexbody<nbrflexbody; iflexbody++)
         {
             if((nbrframesperflexbody[iflexbody]) > 2)
             {
                for(int i=1; i<(nbrframesperflexbody[iflexbody]);i++)
                {
                  flexbody[iflexbody].MTT[i][i]=flexbody[iflexbody].MTT[1][1];
                  flexbody[iflexbody].MTR[i][i]=flexbody[iflexbody].MTR[1][1];
                  flexbody[iflexbody].MRT[i][i]=flexbody[iflexbody].MRT[1][1];
                  flexbody[iflexbody].MRR[i][i]=flexbody[iflexbody].MRR[1][1];
                  flexbody[iflexbody].KTT[i][i]=flexbody[iflexbody].KTT[1][1];
                  flexbody[iflexbody].KTR[i][i]=flexbody[iflexbody].KTR[1][1];
                  flexbody[iflexbody].KRT[i][i]=flexbody[iflexbody].KRT[1][1];
                  flexbody[iflexbody].KRR[i][i]=flexbody[iflexbody].KRR[1][1];
                }

                for(int i=0; i< nbrframesperflexbody[iflexbody];i++) //tensor line frame
                {
                  for(int j=0; j< nbrframesperflexbody[iflexbody];j++) //tensor column frame
                  {
                      if(i==j && i>0 && i<(nbrframesperflexbody[iflexbody]-1))
                      {
                          flexbody[iflexbody].MTT[i][j]=flexbody[iflexbody].MTT[i][j]+flexbody[iflexbody].MTT[0][0];
                          flexbody[iflexbody].MTR[i][j]=flexbody[iflexbody].MTR[i][j]+flexbody[iflexbody].MTR[0][0];
                          flexbody[iflexbody].MRT[i][j]=flexbody[iflexbody].MRT[i][j]+flexbody[iflexbody].MRT[0][0];
                          flexbody[iflexbody].MRR[i][j]=flexbody[iflexbody].MRR[i][j]+flexbody[iflexbody].MRR[0][0];
                          flexbody[iflexbody].KTT[i][j]=flexbody[iflexbody].KTT[i][j]+flexbody[iflexbody].KTT[0][0];
                          flexbody[iflexbody].KTR[i][j]=flexbody[iflexbody].KTR[i][j]+flexbody[iflexbody].KTR[0][0];
                          flexbody[iflexbody].KRT[i][j]=flexbody[iflexbody].KRT[i][j]+flexbody[iflexbody].KRT[0][0];
                          flexbody[iflexbody].KRR[i][j]=flexbody[iflexbody].KRR[i][j]+flexbody[iflexbody].KRR[0][0];
                      }

                      if(i+1==j && j>1)
                      {
                          flexbody[iflexbody].MTT[i][j]=flexbody[iflexbody].MTT[0][1];
                          flexbody[iflexbody].MTR[i][j]=flexbody[iflexbody].MTR[0][1];
                          flexbody[iflexbody].MRT[i][j]=flexbody[iflexbody].MRT[0][1];
                          flexbody[iflexbody].MRR[i][j]=flexbody[iflexbody].MRR[0][1];
                          flexbody[iflexbody].KTT[i][j]=flexbody[iflexbody].KTT[0][1];
                          flexbody[iflexbody].KTR[i][j]=flexbody[iflexbody].KTR[0][1];
                          flexbody[iflexbody].KRT[i][j]=flexbody[iflexbody].KRT[0][1];
                          flexbody[iflexbody].KRR[i][j]=flexbody[iflexbody].KRR[0][1];
                      }

                      if(j+1==i && i>1)
                      {
                          flexbody[iflexbody].MTT[i][j]=flexbody[iflexbody].MTT[1][0];
                          flexbody[iflexbody].MTR[i][j]=flexbody[iflexbody].MTR[1][0];
                          flexbody[iflexbody].MRT[i][j]=flexbody[iflexbody].MRT[1][0];
                          flexbody[iflexbody].MRR[i][j]=flexbody[iflexbody].MRR[1][0];
                          flexbody[iflexbody].KTT[i][j]=flexbody[iflexbody].KTT[1][0];
                          flexbody[iflexbody].KTR[i][j]=flexbody[iflexbody].KTR[1][0];
                          flexbody[iflexbody].KRT[i][j]=flexbody[iflexbody].KRT[1][0];
                          flexbody[iflexbody].KRR[i][j]=flexbody[iflexbody].KRR[1][0];
                      }
                  }
                }
             }
         }
   }

}

//--------------------------------------------------------------------

void ReadMassMatrix(int iflexbody, char name[100])
{
    //Read Mass Matrix from Extern File (.mm file)
    int matrix_size = nbrframesperflexbody[iflexbody] * 6; //(nframe * 6 X nframe * 6)
    double mass[matrix_size][matrix_size];
    #define MASS_MATRIX_EXT "mm"
    cout<<"\nRead Mass Matrix from an External File"<<endl;
    cout<<"--------------------------------"<<endl;
    cout<<"Name of File .("<<MASS_MATRIX_EXT<<") for flexbody "<<iflexbody<<": "<<endl;
    cout<<"->";
    char filename_current[100];
    if(name != NULL)
    {
      filename=name;
    }
    else
    {
      cin >> filename_current;  //entrer le fichier de données
      filename=filename_current;
    }

    char nameextension[104];
    strcpy(nameextension,filename);
    strcat(nameextension,".");
    strcat(nameextension,MASS_MATRIX_EXT);

    ifstream data_file(nameextension);
    if (data_file.is_open()) //Test to open file
    {
         vector<double> numbers;  // Count number of values in file
         double n;
         while(data_file >> n)  numbers.push_back(n);

         if(numbers.size() > matrix_size*matrix_size) //Test to know if the input matrix size is correct
         {
             cout<<"Too Many Arguments\n " <<endl;
             cout<<"Number of values in the file is: "<<numbers.size()<<endl;
             cout<<"-> "<<(numbers.size() - (matrix_size*matrix_size))<<" argument(s) in excess"<<endl;
             data_file.close();
             char dummy[3];
             char* useless=fgets(dummy,2,stdin);
             EasyDynmbsError("Exit Program");
         }
         if(numbers.size() < matrix_size*matrix_size) //Test to know if the input matrix size is correct
         {
             cout<<"Too few arguments\n " <<endl;
             cout<<"Number of values in the file is: "<<numbers.size()<<endl;
             cout<<"-> "<<((matrix_size*matrix_size) - numbers.size())<<" argument(s) missing"<<endl;
             data_file.close();
             char dummy[3];
             char* useless=fgets(dummy,2,stdin);
             EasyDynmbsError("Exit Program");
         }

         size_t index = 0;
         for(int i=0; i <matrix_size; i++) //Reading Mass Matrix
         {
            for(int j=0; j <matrix_size; j++)
            {
              mass[i][j] = numbers[index];
              index++;
            }
         }

         data_file.close();

         for(int i=0; i<(nbrframesperflexbody[iflexbody]);i++)  //Transfer the mass matrix into tensors
         {
            for(int j=0; j<(nbrframesperflexbody[iflexbody]);j++)
            {
               //Mass Matrix
               flexbody[iflexbody].MTT[i][j].put(mass[0+(i*6)][0+(j*6)],mass[0+(i*6)][1+(j*6)],mass[0+(i*6)][2+(j*6)],mass[1+(i*6)][0+(j*6)],mass[1+(i*6)][1+(j*6)],mass[1+(i*6)][2+(j*6)],mass[2+(i*6)][0+(j*6)],mass[2+(i*6)][1+(j*6)],mass[2+(i*6)][2+(j*6)]);
               flexbody[iflexbody].MTR[i][j].put(mass[0+(i*6)][3+(j*6)],mass[0+(i*6)][4+(j*6)],mass[0+(i*6)][5+(j*6)],mass[1+(i*6)][3+(j*6)],mass[1+(i*6)][4+(j*6)],mass[1+(i*6)][5+(j*6)],mass[2+(i*6)][3+(j*6)],mass[2+(i*6)][4+(j*6)],mass[2+(i*6)][5+(j*6)]);
               flexbody[iflexbody].MRT[i][j].put(mass[3+(i*6)][0+(j*6)],mass[3+(i*6)][1+(j*6)],mass[3+(i*6)][2+(j*6)],mass[4+(i*6)][0+(j*6)],mass[4+(i*6)][1+(j*6)],mass[4+(i*6)][2+(j*6)],mass[5+(i*6)][0+(j*6)],mass[5+(i*6)][1+(j*6)],mass[5+(i*6)][2+(j*6)]);
               flexbody[iflexbody].MRR[i][j].put(mass[3+(i*6)][3+(j*6)],mass[3+(i*6)][4+(j*6)],mass[3+(i*6)][5+(j*6)],mass[4+(i*6)][3+(j*6)],mass[4+(i*6)][4+(j*6)],mass[4+(i*6)][5+(j*6)],mass[5+(i*6)][3+(j*6)],mass[5+(i*6)][4+(j*6)],mass[5+(i*6)][5+(j*6)]);
            }
         }


    }
    else // Unable to open file
    {
         cout << "\nCannot Open File "<< nameextension <<endl;
         cout << "\nBad entry or no such file ... \n" <<endl;
         data_file.close();
         char dummy[3];
         char* useless=fgets(dummy,2,stdin);
         EasyDynmbsError("Exit Program");
    }
}

//--------------------------------------------------------------------

void ReadStiffnessMatrix(int iflexbody, char name[100])
{
    //Read Stiffness Matrix from Extern File (.kk file)
    int matrix_size = nbrframesperflexbody[iflexbody] * 6; //(nframe * 6 X nframe * 6)
    double stiffness[matrix_size][matrix_size];
    #define STIFF_MATRIX_EXT "kk"
    cout<<"\nRead Stiffness Matrix from an External File"<<endl;
    cout<<"--------------------------------"<<endl;
    cout<<"Name of File .("<<STIFF_MATRIX_EXT<<") for flexbody "<<iflexbody<<": "<<endl;
    cout<<"->";
    char filename_current[100];
    if(name != NULL)
    {
      filename=name;
    }
    else
    {
      cin >> filename_current;  //entrer le fichier de données
      filename=filename_current;
    }

    char nameextension[104];
    strcpy(nameextension,filename);
    strcat(nameextension,".");
    strcat(nameextension,STIFF_MATRIX_EXT);

    ifstream data_file(nameextension);
    if (data_file.is_open()) //Test to open file
    {
         vector<double> numbers;  // Count number of values in file
         double n;
         while(data_file >> n)  numbers.push_back(n);

         if(numbers.size() > matrix_size*matrix_size) //Test to know if the input matrix size is correct
         {
             cout<<"Too Many Arguments\n " <<endl;
             cout<<"Number of values in the file is: "<<numbers.size()<<endl;
             cout<<"-> "<<(numbers.size() - (matrix_size*matrix_size))<<" argument(s) in excess"<<endl;
             data_file.close();
             char dummy[3];
             char* useless=fgets(dummy,2,stdin);
             EasyDynmbsError("Exit Program");
         }
         if(numbers.size() < matrix_size*matrix_size) //Test to know if the input matrix size is correct
         {
             cout<<"Too few arguments\n " <<endl;
             cout<<"Number of values in the file is: "<<numbers.size()<<endl;
             cout<<"-> "<<((matrix_size*matrix_size) - numbers.size())<<" argument(s) missing"<<endl;
             data_file.close();
             char dummy[3];
             char* useless=fgets(dummy,2,stdin);
             EasyDynmbsError("Exit Program");
         }

         size_t index = 0;
         for(int i=0; i <matrix_size; i++) //Reading Stiffness Matrix
         {
            for(int j=0; j <matrix_size; j++)
            {
              stiffness[i][j] = numbers[index];
              index++;
            }
         }

         data_file.close();

         for(int i=0; i<(nbrframesperflexbody[iflexbody]);i++)  //Transfer the stiffness matrix into tensors
         {
            for(int j=0; j<(nbrframesperflexbody[iflexbody]);j++)
            {
               //Stiffness Matrix
               flexbody[iflexbody].KTT[i][j].put(stiffness[0+(i*6)][0+(j*6)],stiffness[0+(i*6)][1+(j*6)],stiffness[0+(i*6)][2+(j*6)],stiffness[1+(i*6)][0+(j*6)],stiffness[1+(i*6)][1+(j*6)],stiffness[1+(i*6)][2+(j*6)],stiffness[2+(i*6)][0+(j*6)],stiffness[2+(i*6)][1+(j*6)],stiffness[2+(i*6)][2+(j*6)]);
               flexbody[iflexbody].KTR[i][j].put(stiffness[0+(i*6)][3+(j*6)],stiffness[0+(i*6)][4+(j*6)],stiffness[0+(i*6)][5+(j*6)],stiffness[1+(i*6)][3+(j*6)],stiffness[1+(i*6)][4+(j*6)],stiffness[1+(i*6)][5+(j*6)],stiffness[2+(i*6)][3+(j*6)],stiffness[2+(i*6)][4+(j*6)],stiffness[2+(i*6)][5+(j*6)]);
               flexbody[iflexbody].KRT[i][j].put(stiffness[3+(i*6)][0+(j*6)],stiffness[3+(i*6)][1+(j*6)],stiffness[3+(i*6)][2+(j*6)],stiffness[4+(i*6)][0+(j*6)],stiffness[4+(i*6)][1+(j*6)],stiffness[4+(i*6)][2+(j*6)],stiffness[5+(i*6)][0+(j*6)],stiffness[5+(i*6)][1+(j*6)],stiffness[5+(i*6)][2+(j*6)]);
               flexbody[iflexbody].KRR[i][j].put(stiffness[3+(i*6)][3+(j*6)],stiffness[3+(i*6)][4+(j*6)],stiffness[3+(i*6)][5+(j*6)],stiffness[4+(i*6)][3+(j*6)],stiffness[4+(i*6)][4+(j*6)],stiffness[4+(i*6)][5+(j*6)],stiffness[5+(i*6)][3+(j*6)],stiffness[5+(i*6)][4+(j*6)],stiffness[5+(i*6)][5+(j*6)]);
            }
         }


    }
    else // Unable to open file
    {
         cout << "\nCannot Open File "<< nameextension <<endl;
         cout << "\nBad entry or no such file ... \n" <<endl;
         data_file.close();
         char dummy[3];
         char* useless=fgets(dummy,2,stdin);
         EasyDynmbsError("Exit Program");
    }
}

//--------------------------------------------------------------------

void DisplayMassMatrix(int iflexbody)
{
     //Dipslay the mass matrix of flexbody: iflexbody

     cout<<"Display Mass Matrix of flexbody ["<< iflexbody<<"]" <<endl;
     for(int i=0; i< nbrframesperflexbody[iflexbody];i++) //tensor line frame
     {
        for(int j=0; j< nbrframesperflexbody[iflexbody];j++) //tensor column frame
        {
         cout<<" "<<endl;
         cout<<"MTT ["<<i<<"]["<< j <<"]: "<<endl;
         cout<<flexbody[iflexbody].MTT[i][j]<< endl;
         cout<<" "<<endl;
         cout<<"MTR ["<<i<<"]["<< j <<"]: "<<endl;
         cout<<flexbody[iflexbody].MTR[i][j]<< endl;
         cout<<" "<<endl;
         cout<<"MRT ["<<i<<"]["<< j <<"]: "<<endl;
         cout<<flexbody[iflexbody].MRT[i][j]<< endl;
         cout<<" "<<endl;
         cout<<"MRR ["<<i<<"]["<< j <<"]: "<<endl;
         cout<<flexbody[iflexbody].MRR[i][j]<< endl;
        }
     }
}

//--------------------------------------------------------------------

void DisplayStiffnessMatrix(int iflexbody)
{
     //Dipslay the stiffness matrix of flexbody: iflexbody

    cout<<"Display Stiffness Matrix of flexbody ["<< iflexbody<<"]" <<endl;
    for(int i=0; i< nbrframesperflexbody[iflexbody];i++) //tensor line frame
    {
        for(int j=0; j< nbrframesperflexbody[iflexbody];j++) //tensor column frame
        {
         cout<<" "<<endl;
         cout<<"KTT ["<<i<<"]["<< j <<"]: "<<endl;
         cout<<flexbody[iflexbody].KTT[i][j]<< endl;
         cout<<" "<<endl;
         cout<<"KTR ["<<i<<"]["<< j <<"]: "<<endl;
         cout<<flexbody[iflexbody].KTR[i][j]<< endl;
         cout<<" "<<endl;
         cout<<"KRT ["<<i<<"]["<< j <<"]: "<<endl;
         cout<<flexbody[iflexbody].KRT[i][j]<< endl;
         cout<<" "<<endl;
         cout<<"KRR ["<<i<<"]["<< j <<"]: "<<endl;
         cout<<flexbody[iflexbody].KRR[i][j]<< endl;
        }
    }
}

//--------------------------------------------------------------------

void DisplayDampingMatrix(int iflexbody)
{
     //Dipslay the damping matrix of flexbody: iflexbody

    cout<<"Display Damping Matrix of flexbody ["<< iflexbody<<"]" <<endl;
    for(int i=0; i< nbrframesperflexbody[iflexbody];i++) //tensor line frame
    {
        for(int j=0; j< nbrframesperflexbody[iflexbody];j++) //tensor column frame
        {
         cout<<" "<<endl;
         cout<<"CTT ["<<i<<"]["<< j <<"]: "<<endl;
         cout<<flexbody[iflexbody].CTT[i][j]<< endl;
         cout<<" "<<endl;
         cout<<"CTR ["<<i<<"]["<< j <<"]: "<<endl;
         cout<<flexbody[iflexbody].CTR[i][j]<< endl;
         cout<<" "<<endl;
         cout<<"CRT ["<<i<<"]["<< j <<"]: "<<endl;
         cout<<flexbody[iflexbody].CRT[i][j]<< endl;
         cout<<" "<<endl;
         cout<<"CRR ["<<i<<"]["<< j <<"]: "<<endl;
         cout<<flexbody[iflexbody].CRR[i][j]<< endl;
        }
    }
}

//--------------------------------------------------------------------

void ComputeInertiaEfforts()

{
int ibody, iflexbody;
vec LG;
double angle_x, angle_y, angle_z; //for Getrotxyz
vec axe_x, axe_y, axe_z;
vec aF_omegacoro_vF, omegad_omegacoro_omega, omegacoro_omega; //for intermediate computations

  // Computation of reactions of inertia: rigid
  for (ibody=0; ibody<nbrbody; ibody++)
  {
  body[ibody].R=-body[ibody].mass*body[ibody].aG;
  // Computation of PHI*wd first in the local frame
  LG=body[ibody].PhiG*(body[ibody].T0G.R.inv()*body[ibody].omegad);
  body[ibody].MG=-(body[ibody].T0G.R*LG);
  LG=body[ibody].T0G.R*(body[ibody].PhiG*(body[ibody].T0G.R.inv()*body[ibody].omega));
  body[ibody].MG-=(body[ibody].omega^LG);
  }
  // Computation of reactions of inertia: flexible (corotational approach)
  if (nbrflexbody>0)
  {
     for (iflexbody=0; iflexbody<nbrflexbody; iflexbody++)
     {
        flexbody[iflexbody].MG=0;
        for(int i=0; i< nbrframesperflexbody[iflexbody];i++) //tensor line frame
        {
           flexbody[iflexbody].frame[i].R=0;
           flexbody[iflexbody].frame[i].MG=0;
           for(int j=0; j< nbrframesperflexbody[iflexbody];j++) //tensor column frame
           {
                {//Contributions of the Mass Matrix
                //Preliminary computations
                aF_omegacoro_vF = ((flexbody[iflexbody].frame[j].aF - flexgravity) - ((flexbody[iflexbody].omegacoro)^(flexbody[iflexbody].frame[j].vF)) );
                omegad_omegacoro_omega = ((flexbody[iflexbody].frame[j].omegad) - ((flexbody[iflexbody].omegacoro)^(flexbody[iflexbody].frame[j].omega)) );
                omegacoro_omega = (flexbody[iflexbody].omegacoro-flexbody[iflexbody].frame[i].omega);

                  //if(flexbody[iflexbody].length > 0.0) //Beam Case: speed up the computation with BeamTRrotate and DiagRotate
                  /*{
                    //Reactions and torques
                    flexbody[iflexbody].frame[i].R-=   ((flexbody[iflexbody].MTT[i][j].DiagRotate(flexbody[iflexbody].T0coro.R)*(aF_omegacoro_vF        )) + ((flexbody[iflexbody].omegacoro)^(flexbody[iflexbody].MTT[i][j].DiagRotate(flexbody[iflexbody].T0coro.R)*(flexbody[iflexbody].frame[j].vF))) ); //contribution Mtitj
                    flexbody[iflexbody].frame[i].R-=   ((flexbody[iflexbody].MTR[i][j].BeamTRrotate(flexbody[iflexbody].T0coro.R)*(omegad_omegacoro_omega )) + ((flexbody[iflexbody].omegacoro)^(flexbody[iflexbody].MTR[i][j].BeamTRrotate(flexbody[iflexbody].T0coro.R)*(flexbody[iflexbody].frame[j].omega)))); //contribution Mtirj
                    flexbody[iflexbody].frame[i].MG-=  ((flexbody[iflexbody].MRR[i][j].DiagRotate(flexbody[iflexbody].T0coro.R)*(omegad_omegacoro_omega )) + ((omegacoro_omega)^(flexbody[iflexbody].MRR[i][j].DiagRotate(flexbody[iflexbody].T0coro.R)*(flexbody[iflexbody].frame[j].omega)))); //contribution Mrirj
                    flexbody[iflexbody].frame[i].MG-=  ((flexbody[iflexbody].MRT[i][j].BeamTRrotate(flexbody[iflexbody].T0coro.R)*(aF_omegacoro_vF        )) + ((flexbody[iflexbody].MRT[i][j].BeamTRrotate(flexbody[iflexbody].T0coro.R) * (flexbody[iflexbody].frame[j].vF))^(flexbody[iflexbody].frame[i].omega-flexbody[iflexbody].omegacoro))); //contribution Mritj
                    //Torque on * (corotational frame)
                    flexbody[iflexbody].MG-=  (flexbody[iflexbody].frame[i].vF    ^(flexbody[iflexbody].MTT[i][j].DiagRotate(flexbody[iflexbody].T0coro.R)*(flexbody[iflexbody].frame[j].vF))); //contribtion Mtitj
                    flexbody[iflexbody].MG-=  (flexbody[iflexbody].frame[i].omega ^(flexbody[iflexbody].MRR[i][j].DiagRotate(flexbody[iflexbody].T0coro.R)*(flexbody[iflexbody].frame[j].omega))); //contribution Mrirj
                    flexbody[iflexbody].MG-=  ((flexbody[iflexbody].frame[i].omega^(flexbody[iflexbody].MRT[i][j].BeamTRrotate(flexbody[iflexbody].T0coro.R)*(flexbody[iflexbody].frame[j].vF)))    + ((flexbody[iflexbody].frame[i].vF)   ^(flexbody[iflexbody].MTR[i][j].BeamTRrotate(flexbody[iflexbody].T0coro.R)*(flexbody[iflexbody].frame[j].omega)))); //contribution Mtirj et Mritj
                  }*/
                  //else
                  //{
                    //Reactions and torques
                    flexbody[iflexbody].frame[i].R-=   ((flexbody[iflexbody].MTT[i][j].rotate(flexbody[iflexbody].T0coro.R)*(aF_omegacoro_vF        )) + ((flexbody[iflexbody].omegacoro)^(flexbody[iflexbody].MTT[i][j].rotate(flexbody[iflexbody].T0coro.R)*(flexbody[iflexbody].frame[j].vF))) ); //contribution Mtitj
                    flexbody[iflexbody].frame[i].R-=   ((flexbody[iflexbody].MTR[i][j].rotate(flexbody[iflexbody].T0coro.R)*(omegad_omegacoro_omega )) + ((flexbody[iflexbody].omegacoro)^(flexbody[iflexbody].MTR[i][j].rotate(flexbody[iflexbody].T0coro.R)*(flexbody[iflexbody].frame[j].omega)))); //contribution Mtirj
                    flexbody[iflexbody].frame[i].MG-=  ((flexbody[iflexbody].MRR[i][j].rotate(flexbody[iflexbody].T0coro.R)*(omegad_omegacoro_omega )) + ((omegacoro_omega)^(flexbody[iflexbody].MRR[i][j].rotate(flexbody[iflexbody].T0coro.R)*(flexbody[iflexbody].frame[j].omega)))); //contribution Mrirj
                    flexbody[iflexbody].frame[i].MG-=  ((flexbody[iflexbody].MRT[i][j].rotate(flexbody[iflexbody].T0coro.R)*(aF_omegacoro_vF        )) + ((flexbody[iflexbody].MRT[i][j].rotate(flexbody[iflexbody].T0coro.R) * (flexbody[iflexbody].frame[j].vF))^(flexbody[iflexbody].frame[i].omega-flexbody[iflexbody].omegacoro))); //contribution Mritj
                    //Torque on * (corotational frame)
                    flexbody[iflexbody].MG-=  (flexbody[iflexbody].frame[i].vF    ^(flexbody[iflexbody].MTT[i][j].rotate(flexbody[iflexbody].T0coro.R)*(flexbody[iflexbody].frame[j].vF))); //contribtion Mtitj
                    flexbody[iflexbody].MG-=  (flexbody[iflexbody].frame[i].omega ^(flexbody[iflexbody].MRR[i][j].rotate(flexbody[iflexbody].T0coro.R)*(flexbody[iflexbody].frame[j].omega))); //contribution Mrirj
                    flexbody[iflexbody].MG-=  ((flexbody[iflexbody].frame[i].omega^(flexbody[iflexbody].MRT[i][j].rotate(flexbody[iflexbody].T0coro.R)*(flexbody[iflexbody].frame[j].vF)))    + ((flexbody[iflexbody].frame[i].vF)   ^(flexbody[iflexbody].MTR[i][j].rotate(flexbody[iflexbody].T0coro.R)*(flexbody[iflexbody].frame[j].omega)))); //contribution Mtirj et Mritj
                  //}
                }
            }

           //Strains
           flexbody[iflexbody].frame[i].TcoroF0_res=flexbody[iflexbody].T0coro*flexbody[iflexbody].frame[i].TcoroF0;
           flexbody[iflexbody].frame[i].deform = flexbody[iflexbody].frame[i].T0F.e - flexbody[iflexbody].frame[i].TcoroF0_res.e ;
           getrotxyz(  flexbody[iflexbody].frame[i].TcoroF0_res.R, flexbody[iflexbody].frame[i].T0F.R, angle_x, angle_y, angle_z, axe_x, axe_y, axe_z);
           flexbody[iflexbody].frame[i].theta.put(angle_x,angle_y,angle_z);

           //Strain Velocities
           flexbody[iflexbody].frame[i].Vdeform = flexbody[iflexbody].frame[i].vF    - (flexbody[iflexbody].omegacoro ^ (flexbody[iflexbody].frame[i].T0F.e - flexbody[iflexbody].T0coro.e));
           flexbody[iflexbody].frame[i].Vtheta  = flexbody[iflexbody].frame[i].omega - (flexbody[iflexbody].omegacoro);

        }

        //Contributions of the Stiffness Matrix
        for(int i=0; i< nbrframesperflexbody[iflexbody];i++) //tensor line frame
        {
           for(int j=0; j< nbrframesperflexbody[iflexbody];j++) //tensor column frame
           {
                /*if(flexbody[iflexbody].length >  0.0) //Beam Case: speed up the computation with BeamTRrotate and DiagRotate
                {
                   flexbody[iflexbody].frame[i].R-=  ((flexbody[iflexbody].KTT[i][j].DiagRotate(flexbody[iflexbody].T0coro.R) * flexbody[iflexbody].frame[j].deform) + (flexbody[iflexbody].T0coro.R * (flexbody[iflexbody].KTR[i][j] * flexbody[iflexbody].frame[j].theta)) );
                   flexbody[iflexbody].frame[i].MG-= ((flexbody[iflexbody].KRT[i][j].BeamTRrotate(flexbody[iflexbody].T0coro.R) * flexbody[iflexbody].frame[j].deform) + (flexbody[iflexbody].T0coro.R * (flexbody[iflexbody].KRR[i][j] * flexbody[iflexbody].frame[j].theta)) );
                }
                else
                {*/
                   flexbody[iflexbody].frame[i].R-=  ((flexbody[iflexbody].KTT[i][j].rotate(flexbody[iflexbody].T0coro.R) * flexbody[iflexbody].frame[j].deform) + (flexbody[iflexbody].T0coro.R * (flexbody[iflexbody].KTR[i][j] * flexbody[iflexbody].frame[j].theta)) );
                   flexbody[iflexbody].frame[i].MG-= ((flexbody[iflexbody].KRT[i][j].rotate(flexbody[iflexbody].T0coro.R) * flexbody[iflexbody].frame[j].deform) + (flexbody[iflexbody].T0coro.R * (flexbody[iflexbody].KRR[i][j] * flexbody[iflexbody].frame[j].theta)) );
                //}

           }
        }

        //Rayleigh Damping through the Damping Matrix C = AlphaDamp * M + BetaDamp * K
        if(AlphaDamp > 0.0 || BetaDamp > 0.0)
        {
            //Contributions of the Damping Matrix
            for(int i=0; i< nbrframesperflexbody[iflexbody];i++) //tensor line frame
            {
               for(int j=0; j< nbrframesperflexbody[iflexbody];j++) //tensor column frame
               {
                   flexbody[iflexbody].frame[i].R-=  ((flexbody[iflexbody].CTT[i][j].rotate(flexbody[iflexbody].T0coro.R) * flexbody[iflexbody].frame[j].Vdeform) + (flexbody[iflexbody].T0coro.R * (flexbody[iflexbody].CTR[i][j] * flexbody[iflexbody].frame[j].Vtheta)) );
                   flexbody[iflexbody].frame[i].MG-= ((flexbody[iflexbody].CRT[i][j].rotate(flexbody[iflexbody].T0coro.R) * flexbody[iflexbody].frame[j].Vdeform) + (flexbody[iflexbody].T0coro.R * (flexbody[iflexbody].CRR[i][j] * flexbody[iflexbody].frame[j].Vtheta)) );

               }
            }
        }

     }

  }
}

//--------------------------------------------------------------------

void ComputePartialVelocitiesDefault()

{

ComputeMotion();

int ibody,idof,iflexbody,iframe;

// Store initial state
for (ibody=0; ibody<nbrbody; ibody++)
    {
    body[ibody].vGsto=body[ibody].vG;
    body[ibody].omegasto=body[ibody].omega;
    }
for (iflexbody=0; iflexbody<nbrflexbody; iflexbody++)
    {
    flexbody[iflexbody].omegacorosto=flexbody[iflexbody].omegacoro;
    for (iframe=0; iframe<flexbody[iflexbody].nbrframe; iframe++)
       {
       flexbody[iflexbody].frame[iframe].vFsto=flexbody[iflexbody].frame[iframe].vF;
       flexbody[iflexbody].frame[iframe].omegasto=flexbody[iflexbody].frame[iframe].omega;
       }
    }

// perform numerical derivation
for (idof=0; idof<nbrdof; idof++)
  {
  qd[idof]+=1.0;
  ComputeMotion();
  qd[idof]-=1.0;

    for (ibody=0; ibody<nbrbody; ibody++)
    {
    body[ibody].vGpartial[idof]=body[ibody].vG-body[ibody].vGsto;
    body[ibody].omegapartial[idof]=body[ibody].omega-body[ibody].omegasto;
    }
    for (iflexbody=0; iflexbody<nbrflexbody; iflexbody++)
    {
    flexbody[iflexbody].omegacoropartial[idof]=flexbody[iflexbody].omegacoro
                                              -flexbody[iflexbody].omegacorosto;
       for (iframe=0; iframe<flexbody[iflexbody].nbrframe; iframe++)
       {
       flexbody[iflexbody].frame[iframe].vFpartial[idof]=flexbody[iflexbody].frame[iframe].vF
	                                                -flexbody[iflexbody].frame[iframe].vFsto;;
       flexbody[iflexbody].frame[iframe].omegapartial[idof]=flexbody[iflexbody].frame[iframe].omega
	                                                -flexbody[iflexbody].frame[iframe].omegasto;
       }
    }
  }

// Restore intial state
ComputeMotion();

/*for (ibody=0; ibody<nbrbody; ibody++)
    {
    body[ibody].vG=body[ibody].vGsto;
    body[ibody].omega=body[ibody].omegasto;
    }
for (iflexbody=0; iflexbody<nbrflexbody; iflexbody++)
    {
    flexbody[iflexbody].omegacoro=flexbody[iflexbody].omegacorosto;
    for (iframe=0; iframe<flexbody[iflexbody].nbrframe; iframe++)
       {
       flexbody[iflexbody].frame[iframe].vF=flexbody[iflexbody].frame[iframe].vFsto;
       flexbody[iflexbody].frame[iframe].omega=flexbody[iflexbody].frame[iframe].omegasto;
       }
       }*/

}

//--------------------------------------------------------------------

void ComputeResidualmbs()

{
// Let's control the user assiduity
if ((body==0) && (flexbody==0))
  EasyDynmbsError("Probably InitMbsSim was not called (body or flexbody not allocated)");

ComputeMotion();
ComputePartialVelocities();
ComputeInertiaEfforts();
AddAppliedEfforts();

int ibody;

    for (int idof=0; idof<nbrdof; idof++)
    {
        f[idof]=0;
	// rigid bodies
        for (ibody=0; ibody<nbrbody; ibody++)
         f[idof]-=body[ibody].R*body[ibody].vGpartial[idof]
    	    +body[ibody].MG*body[ibody].omegapartial[idof];
	// and now flexible bodies
        if(nbrflexbody>0) //flexible bodies
        {
           for (int iflexbody=0; iflexbody<nbrflexbody; iflexbody++)
           {
              f[idof]-=flexbody[iflexbody].MG*flexbody[iflexbody].omegacoropartial[idof]; //corotational frame
              for (int i=0; i<nbrframesperflexbody[iflexbody]; i++)
              {
                  f[idof]-=flexbody[iflexbody].frame[i].R*flexbody[iflexbody].frame[i].vFpartial[idof]      //other frames
                  	   +flexbody[iflexbody].frame[i].MG*flexbody[iflexbody].frame[i].omegapartial[idof];
              }
           }
	    }
    }

}

//--------------------------------------------------------------------

void ComputeMassMatrix(int MotionToCompute=0)

{
int ibody,idof,jdof;
vec LG;
double jacterm;

if (MotionToCompute) ComputeMotion();


// Computation of mass matrix from partial velocities
for (idof=0; idof<nbrdof; idof++) for (jdof=0; jdof<=idof; jdof++)
  {
  jacterm=0;
  for (ibody=0; ibody<nbrbody; ibody++)
    {
    jacterm+=body[ibody].mass*(body[ibody].vGpartial[idof]
				      *body[ibody].vGpartial[jdof]);
    LG=body[ibody].T0G.R*(body[ibody].PhiG*
            (body[ibody].T0G.R.inv()*body[ibody].omegapartial[jdof]));
    jacterm+=body[ibody].omegapartial[idof]*LG;
    }
  gsl_matrix_set(MM,idof,jdof,jacterm);
  if (jdof<idof) gsl_matrix_set(MM,jdof,idof,jacterm);
  }
}

//--------------------------------------------------------------------

void ComputeAccelerations()

{
  int ibody,idof,jdof,s;

for (idof=0; idof<nbrdof; idof++) qdd[idof]=0;
ComputeResidualmbs();
ComputeMassMatrix();
gsl_linalg_LU_decomp(MM,permmbs,&s);
gsl_linalg_LU_solve(MM,permmbs,&fmbsview.vector,&qddmbsview.vector);
for (idof=0; idof<nbrdof; idof++)
   {
   qdd[idof]=-qdd[idof];
   xd[idof]=qd[idof];
   }
}

//--------------------------------------------------------------------

void ComposeMotionOld(int ibody, int ibodyref)

{
  if ((ibody<0) || (ibody>nbrbody))
   EasyDynmbsError("ComposeMotion, incorrect number for ibody");
  if ((ibodyref<0) || (ibodyref>nbrbody))
   EasyDynmbsError("ComposeMotion, incorrect number for ibodyref");
  if (ibodyref>ibody)
   EasyDynmbsError("ComposeMotion, ibody must be larger than ibodyref");

// Let's first project all relative terms in the absolute coordinate system
vec erel=body[ibodyref].T0G.R*body[ibody].TrefG.e,
    vrel=body[ibodyref].T0G.R*body[ibody].vGrel,
    arel=body[ibodyref].T0G.R*body[ibody].aGrel,
    wrel=body[ibodyref].T0G.R*body[ibody].omegarel,
    wdrel=body[ibodyref].T0G.R*body[ibody].omegadrel;

// Let's go for the composition
body[ibody].T0G=body[ibodyref].T0G*body[ibody].TrefG;
body[ibody].vG=body[ibodyref].vG+(body[ibodyref].omega^erel)+vrel;
body[ibody].omega=body[ibodyref].omega+wrel;
body[ibody].aG=body[ibodyref].aG+(body[ibodyref].omegad^erel)
             +(body[ibodyref].omega^(body[ibodyref].omega^erel))
             +2*(body[ibodyref].omega^vrel)+arel;
body[ibody].omegad=body[ibodyref].omegad+(body[ibodyref].omega^wrel)
                  +wdrel;
}

//--------------------------------------------------------------------

void ComposeMotion(int ibody, int ibodyref, int iframe, int iframeref)
{
     //par défaut: iframe = -1 et iframeref=-1 (dans mbs.h)

     //if ((ibody<0) || (ibody>nbrbody))
      //EasyDynmbsError("ComposeMotion, incorrect number for ibody");
     //if ((ibodyref<0) || (ibodyref>nbrbody))
      //EasyDynmbsError("ComposeMotion, incorrect number for ibodyref");
     //if (ibodyref>ibody)
      //EasyDynmbsError("ComposeMotion, ibody must be larger than ibodyref");

     //Case 1: ibody: rigid, ibodyref: rigid, iframe=-1, iframref=-1
     if(iframe < 0 && iframeref < 0)
     {
          // Let's first project all relative terms in the absolute coordinate system
          vec erel=body[ibodyref].T0G.R*body[ibody].TrefG.e,
                vrel=body[ibodyref].T0G.R*body[ibody].vGrel,
                arel=body[ibodyref].T0G.R*body[ibody].aGrel,
                wrel=body[ibodyref].T0G.R*body[ibody].omegarel,
                wdrel=body[ibodyref].T0G.R*body[ibody].omegadrel;

          // Let's go for the composition
          body[ibody].T0G=body[ibodyref].T0G*body[ibody].TrefG;
          body[ibody].vG=body[ibodyref].vG+(body[ibodyref].omega^erel)+vrel;
          body[ibody].omega=body[ibodyref].omega+wrel;
          body[ibody].aG=body[ibodyref].aG+(body[ibodyref].omegad^erel)
                         +(body[ibodyref].omega^(body[ibodyref].omega^erel))
                         +2*(body[ibodyref].omega^vrel)+arel;
          body[ibody].omegad=body[ibodyref].omegad+(body[ibodyref].omega^wrel)
                              +wdrel;
     }

     //Case 2: ibody: rigid, ibodyref: flexible, iframe=-1, iframref>=0
     if(iframe < 0 && iframeref >= 0)
     {
           if(nbrflexbody>0) //for struct flexbody
           {
              // Let's first project all relative terms in the absolute coordinate system
              vec erel=flexbody[ibodyref].frame[iframeref].T0F.R*body[ibody].TrefG.e,
                    vrel=flexbody[ibodyref].frame[iframeref].T0F.R*body[ibody].vGrel,
                    arel=flexbody[ibodyref].frame[iframeref].T0F.R*body[ibody].aGrel,
                    wrel=flexbody[ibodyref].frame[iframeref].T0F.R*body[ibody].omegarel,
                    wdrel=flexbody[ibodyref].frame[iframeref].T0F.R*body[ibody].omegadrel;

              // Let's go for the composition
              body[ibody].T0G=flexbody[ibodyref].frame[iframeref].T0F*body[ibody].TrefG;
              body[ibody].vG=flexbody[ibodyref].frame[iframeref].vF+(flexbody[ibodyref].frame[iframeref].omega^erel)+vrel;
              body[ibody].omega=flexbody[ibodyref].frame[iframeref].omega+wrel;
              body[ibody].aG=flexbody[ibodyref].frame[iframeref].aF+(flexbody[ibodyref].frame[iframeref].omegad^erel)
                             +(flexbody[ibodyref].frame[iframeref].omega^(flexbody[ibodyref].frame[iframeref].omega^erel))
                             +2*(flexbody[ibodyref].frame[iframeref].omega^vrel)+arel;
              body[ibody].omegad=flexbody[ibodyref].frame[iframeref].omegad+(flexbody[ibodyref].frame[iframeref].omega^wrel)
                                  +wdrel;
           }
     }

     //Case 3: ibody: flexible, ibodyref: rigid, iframe>=0, iframref=-1
     if(iframe >= 0 && iframeref < 0)
     {
           if(nbrflexbody>0) //for struct flexbody
           {
              // Let's first project all relative terms in the absolute coordinate system
              vec erel=body[ibodyref].T0G.R*flexbody[ibody].frame[iframe].TrefF.e,
                    vrel=body[ibodyref].T0G.R*flexbody[ibody].frame[iframe].vFrel,
                    arel=body[ibodyref].T0G.R*flexbody[ibody].frame[iframe].aFrel,
                    wrel=body[ibodyref].T0G.R*flexbody[ibody].frame[iframe].omegarel,
                    wdrel=body[ibodyref].T0G.R*flexbody[ibody].frame[iframe].omegadrel;

              // Let's go for the composition
              flexbody[ibody].frame[iframe].T0F=body[ibodyref].T0G*flexbody[ibody].frame[iframe].TrefF;
              flexbody[ibody].frame[iframe].vF=body[ibodyref].vG+(body[ibodyref].omega^erel)+vrel;
              flexbody[ibody].frame[iframe].omega=body[ibodyref].omega+wrel;
              flexbody[ibody].frame[iframe].aF=body[ibodyref].aG+(body[ibodyref].omegad^erel)
                             +(body[ibodyref].omega^(body[ibodyref].omega^erel))
                             +2*(body[ibodyref].omega^vrel)+arel;
              flexbody[ibody].frame[iframe].omegad=body[ibodyref].omegad+(body[ibodyref].omega^wrel)
                                  +wdrel;

              //Corotational Frame
              // Let's first project all relative terms in the absolute coordinate system
              vec wrelcoro=body[ibodyref].T0G.R*flexbody[ibody].omegacororel;

              // Let's go for the composition
              flexbody[ibody].T0coro=body[ibodyref].T0G*flexbody[ibody].T0cororef;
              flexbody[ibody].omegacoro=body[ibodyref].omega+wrelcoro;
           }
     }

     //Case 4: ibody: flexible, ibodyref: flexible, iframe>=0, iframref>=0  (default corotational frame wrt to the highest number of the frame of bodyref)
     if(iframe >= 0 && iframeref >= 0)
     {
           // 2 sous-cas: si c'est par rapport à une frame ou par rapport à un flexbody: ça va changer le calcul du repère corotationnel

           if(nbrflexbody>0) //for struct flexbody
           {
              // iframeref wrt iframe for the same flexbody: pas de calcul pour le repère corotationnel
              // Let's first project all relative terms in the absolute coordinate system
              vec erel=flexbody[ibodyref].frame[iframeref].T0F.R*flexbody[ibody].frame[iframe].TrefF.e,
                    vrel=flexbody[ibodyref].frame[iframeref].T0F.R*flexbody[ibody].frame[iframe].vFrel,
                    arel=flexbody[ibodyref].frame[iframeref].T0F.R*flexbody[ibody].frame[iframe].aFrel,
                    wrel=flexbody[ibodyref].frame[iframeref].T0F.R*flexbody[ibody].frame[iframe].omegarel,
                    wdrel=flexbody[ibodyref].frame[iframeref].T0F.R*flexbody[ibody].frame[iframe].omegadrel;

              // Let's go for the composition
              flexbody[ibody].frame[iframe].T0F=flexbody[ibodyref].frame[iframeref].T0F*flexbody[ibody].frame[iframe].TrefF;
              flexbody[ibody].frame[iframe].vF=flexbody[ibodyref].frame[iframeref].vF+(flexbody[ibodyref].frame[iframeref].omega^erel)+vrel;
              flexbody[ibody].frame[iframe].omega=flexbody[ibodyref].frame[iframeref].omega+wrel;
              flexbody[ibody].frame[iframe].aF=flexbody[ibodyref].frame[iframeref].aF+(flexbody[ibodyref].frame[iframeref].omegad^erel)
                             +(flexbody[ibodyref].frame[iframeref].omega^(flexbody[ibodyref].frame[iframeref].omega^erel))
                             +2*(flexbody[ibodyref].frame[iframeref].omega^vrel)+arel;
              flexbody[ibody].frame[iframe].omegad=flexbody[ibodyref].frame[iframeref].omegad+(flexbody[ibodyref].frame[iframeref].omega^wrel)
                                  +wdrel;

              if(ibody > ibodyref) // iflexbody i>j wrt jflexbodyref j: calcul pour le repère corotationnel par rapport à la frame
                                   // de plus haut indice du flexbody précédent
              {
                  // Let's first project all relative terms in the absolute coordinate system
                  int higherframe=(nbrframesperflexbody[ibodyref])-1; //higher frame of the current flexbody
                  vec wrel=flexbody[ibodyref].frame[higherframe].T0F.R*flexbody[ibody].omegacororel;

                  // Let's go for the composition
                  flexbody[ibody].T0coro=flexbody[ibodyref].frame[higherframe].T0F*flexbody[ibody].T0cororef;
                  flexbody[ibody].omegacoro=flexbody[ibodyref].frame[higherframe].omega+wrel;
              }
           }
     }
}

//--------------------------------------------------------------------

void ComposePartialVelocitiesOld(int ibody, int ibodyref)

{
  if ((ibody<0) || (ibody>nbrbody))
   EasyDynmbsError("ComposePartialVelocities, incorrect number for ibody");
  if ((ibodyref<0) || (ibodyref>nbrbody))
   EasyDynmbsError("ComposePartialVelocities, incorrect number for ibodyref");
  if (ibodyref>ibody)
   EasyDynmbsError("ComposePartialVelocities, ibody must be larger than ibodyref");

// Let's first project all relative terms in the absolute coordinate system
vec erel=body[ibodyref].T0G.R*body[ibody].TrefG.e;
vec vrel, wrel;
int idof;
for (idof=0; idof<nbrdof; idof++)
  {
  vrel=body[ibodyref].T0G.R*body[ibody].vGrelpartial[idof];
  wrel=body[ibodyref].T0G.R*body[ibody].omegarelpartial[idof];
  body[ibody].vGpartial[idof]=body[ibodyref].vGpartial[idof]
                             +(body[ibodyref].omegapartial[idof]^erel)
                             +vrel;
  body[ibody].omegapartial[idof]=body[ibodyref].omegapartial[idof]+wrel;
  }
}

//--------------------------------------------------------------------

void ComposePartialVelocities(int ibody, int ibodyref, int iframe, int iframeref)
{
     //par défaut: iframe = -1 et iframeref=-1 (dans mbs.h)

     //if ((ibody<0) || (ibody>nbrbody))
      //EasyDynmbsError("ComposePartialVelocities, incorrect number for ibody");
     //if ((ibodyref<0) || (ibodyref>nbrbody))
      //EasyDynmbsError("ComposePartialVelocities, incorrect number for ibodyref");
     //if (ibodyref>ibody)
      //EasyDynmbsError("ComposePartialVelocities, ibody must be larger than ibodyref");

     //Case 1: ibody: rigid, ibodyref: rigid, iframe=-1, iframref=-1
     if(iframe < 0 && iframeref < 0)
     {
          // Let's first project all relative terms in the absolute coordinate system
          vec erel=body[ibodyref].T0G.R*body[ibody].TrefG.e;
          vec vrel, wrel;
          int idof;
          for (idof=0; idof<nbrdof; idof++)
          {
               vrel=body[ibodyref].T0G.R*body[ibody].vGrelpartial[idof];
               wrel=body[ibodyref].T0G.R*body[ibody].omegarelpartial[idof];
               body[ibody].vGpartial[idof]=body[ibodyref].vGpartial[idof]
                                          +(body[ibodyref].omegapartial[idof]^erel)
                                          +vrel;
               body[ibody].omegapartial[idof]=body[ibodyref].omegapartial[idof]+wrel;
          }
     }

     //Case 2: ibody: rigid, ibodyref: flexible, iframe=-1, iframeref>=0
     if(iframe < 0 && iframeref >= 0)
     {
           if(nbrflexbody>0) //for struct flexbody
           {
              // Let's first project all relative terms in the absolute coordinate system
              vec erel=flexbody[ibodyref].frame[iframeref].T0F.R*body[ibody].TrefG.e;
              vec vrel, wrel;
              int idof;
              for (idof=0; idof<nbrdof; idof++)
              {
                   vrel=flexbody[ibodyref].frame[iframeref].T0F.R*body[ibody].vGrelpartial[idof];
                   wrel=flexbody[ibodyref].frame[iframeref].T0F.R*body[ibody].omegarelpartial[idof];
                   body[ibody].vGpartial[idof]=flexbody[ibodyref].frame[iframeref].vFpartial[idof]
                                              +(flexbody[ibodyref].frame[iframeref].omegapartial[idof]^erel)
                                              +vrel;
                   body[ibody].omegapartial[idof]=flexbody[ibodyref].frame[iframeref].omegapartial[idof]+wrel;
              }
           }
     }

     //Case 3: ibody: flexible, ibodyref: rigid, iframe>=0, iframref=-1
     if(iframe >= 0 && iframeref < 0)
     {
          if(nbrflexbody>0) //for struct flexbody
          {
              // Let's first project all relative terms in the absolute coordinate system
              vec erel=body[ibodyref].T0G.R*flexbody[ibody].frame[iframe].TrefF.e;
              vec vrel, wrel;
              int idof;
              for (idof=0; idof<nbrdof; idof++)
              {
                   vrel=body[ibodyref].T0G.R*flexbody[ibody].frame[iframe].vFrelpartial[idof];
                   wrel=body[ibodyref].T0G.R*flexbody[ibody].frame[iframe].omegarelpartial[idof];
                   flexbody[ibody].frame[iframe].vFpartial[idof]=body[ibodyref].vGpartial[idof]
                                              +(body[ibodyref].omegapartial[idof]^erel)
                                              +vrel;
                   flexbody[ibody].frame[iframe].omegapartial[idof]=body[ibodyref].omegapartial[idof]+wrel;
              }
          }
     }

     //Case 4: ibody: flexible, ibodyref: flexible, iframe>=0, iframref>=0  (default corotational frame wrt to the highest number of the frame of bodyref)
     if(iframe >= 0 && iframeref >= 0)
     {
           // 2 sous-cas: si c'est par rapport à une frame ou par rapport à un flexbody: ça va changer le calcul du repère corotationnel

           if(nbrflexbody>0) //for struct flexbody
           {
              // iframeref wrt iframe for the same flexbody: pas de calcul pour le repère corotationnel
              // Let's first project all relative terms in the absolute coordinate system
              vec erel=flexbody[ibodyref].frame[iframeref].T0F.R*flexbody[ibody].frame[iframe].TrefF.e;
              vec vrel, wrel;
              int idof;
              for (idof=0; idof<nbrdof; idof++)
              {
                   vrel=flexbody[ibodyref].frame[iframeref].T0F.R*flexbody[ibody].frame[iframe].vFrelpartial[idof];
                   wrel=flexbody[ibodyref].frame[iframeref].T0F.R*flexbody[ibody].frame[iframe].omegarelpartial[idof];
                   flexbody[ibody].frame[iframe].vFpartial[idof]=flexbody[ibodyref].frame[iframeref].vFpartial[idof]
                                              +(flexbody[ibodyref].frame[iframeref].omegapartial[idof]^erel)
                                              +vrel;
                   flexbody[ibody].frame[iframe].omegapartial[idof]=flexbody[ibodyref].frame[iframeref].omegapartial[idof]+wrel;
              }

              if(ibody > ibodyref) // iflexbody i>j wrt jflexbodyref j: calcul pour le repère corotationnel par rapport à la frame
                                   // de plus haut indice du flexbody précédent
              {
                  // Let's first project all relative terms in the absolute coordinate system
                  int higherframe=(nbrframesperflexbody[ibodyref])-1; //higher frame of the current flexbody
                  vec wrel;
                  int idof;
                  for (idof=0; idof<nbrdof; idof++)
                  {
                       wrel=flexbody[ibodyref].frame[higherframe].T0F.R*flexbody[ibody].omegacororelpartial[idof];
                       flexbody[ibody].omegacoropartial[idof]=flexbody[ibodyref].frame[higherframe].omegapartial[idof]+wrel;
                  }
              }
           }
     }
}

//--------------------------------------------------------------------

vec VelPointBody(int ibody, vec r)
{
  vec r2,vpb;
  r2=body[ibody].T0G.R*r; // projection in the global frame
  vpb=body[ibody].vG+(body[ibody].omega^r2);
  return vpb;
}

//--------------------------------------------------------------------

vec AccPointBody(int ibody, vec r)
{
  vec r2,apb;
  r2=body[ibody].T0G.R*r; // projection in the global frame
  apb=body[ibody].aG+(body[ibody].omegad^r2)
    +(body[ibody].omega^(body[ibody].omega^r2));
  return apb;
}

//--------------------------------------------------------------------

void AddGravityForces(vec grav)

{
int ibody;
  for (ibody=0; ibody<nbrbody; ibody++)
  {
    body[ibody].R+=body[ibody].mass*grav;
  }
}

//--------------------------------------------------------------------

void AddSpringForce(double K, double L0,
                    int ibodyA, vec rA, int ibodyB, vec rB)

{
vec uAB=body[ibodyB].T0G*rB-body[ibodyA].T0G*rA;
double L=uAB.length();
uAB.unite();
vec force=K*(L-L0)*uAB;
body[ibodyA].R+=force;
body[ibodyA].MG+=((body[ibodyA].T0G.R*rA)^force);
body[ibodyB].R-=force;
body[ibodyB].MG-=((body[ibodyB].T0G.R*rB)^force);
}

//--------------------------------------------------------------------


void AddSpringForceQB(double K, double L0,
                    int ibodyA, vec rA, int ibodyB, vec rB, vec *force)

{
vec uAB=body[ibodyB].T0G*rB-body[ibodyA].T0G*rA;
double L=uAB.length();
uAB.unite();
(*force)=K*(L-L0)*uAB;
body[ibodyA].R+=(*force);
body[ibodyA].MG+=((body[ibodyA].T0G.R*rA)^(*force));
body[ibodyB].R-=(*force);
body[ibodyB].MG-=((body[ibodyB].T0G.R*rB)^(*force));
}


//--------------------------------------------------------------------

void AddDamperForce(double C,
                    int ibodyA, vec rA, int ibodyB, vec rB)

{
vec uAB=body[ibodyB].T0G*rB-body[ibodyA].T0G*rA;
uAB.unite();
vec vrel=body[ibodyB].vG+(body[ibodyB].omega^(body[ibodyB].T0G.R*rB))
         -body[ibodyA].vG-(body[ibodyA].omega^(body[ibodyA].T0G.R*rA));
vec force=C*(vrel*uAB)*uAB;
body[ibodyA].R+=force;
body[ibodyA].MG+=((body[ibodyA].T0G.R*rA)^force);
body[ibodyB].R-=force;
body[ibodyB].MG-=((body[ibodyB].T0G.R*rB)^force);
}


//--------------------------------------------------------------------


void AddDamperForceQB(double C,
                    int ibodyA, vec rA, int ibodyB, vec rB, vec *force)

{
vec uAB=body[ibodyB].T0G*rB-body[ibodyA].T0G*rA;
uAB.unite();
vec vrel=body[ibodyB].vG+(body[ibodyB].omega^(body[ibodyB].T0G.R*rB))
         -body[ibodyA].vG-(body[ibodyA].omega^(body[ibodyA].T0G.R*rA));
(*force)=C*(vrel*uAB)*uAB;
body[ibodyA].R+=(*force);
body[ibodyA].MG+=((body[ibodyA].T0G.R*rA)^(*force));
body[ibodyB].R-=(*force);
body[ibodyB].MG-=((body[ibodyB].T0G.R*rB)^(*force));
}

//--------------------------------------------------------------------

void AddTyreEfforts(int ibody, vec axe, structtyre tyre)
{
// Determination of the conventional contact coordinate system
// z perpendicular to the ground
// y projection of the rotation axis on the ground
mth tSAE;
vec zSAE(0,0,1), uaxe=body[ibody].T0G.R*axe;
uaxe.unite();
vec xSAE=uaxe^zSAE;
xSAE.unite();
vec ySAE=zSAE^xSAE;
tSAE.R.setux(xSAE);
tSAE.R.setuy(ySAE);
tSAE.R.setuz(zSAE);
tSAE.e=body[ibody].T0G.e+(tyre.r1-tyre.r2)*(uaxe^xSAE)-tyre.r2*zSAE;

// Determination of penetration
double pen=-tSAE.e.z;

if (pen<=0) return;  // all efforts are of course equal to zero

// The SAE coordinate system is adapted so as to take into account the
// effective rolling radius
tSAE.e.z=tSAE.e.z*0.66667;
// Determination of contact point velocity
vec vP=body[ibody].vG+(body[ibody].omega^(tSAE.e-body[ibody].T0G.e));
// Determination of relevant coordinates
double vPxSAE=vP*xSAE, vPySAE=vP*ySAE, vPzSAE=vP*zSAE,
       vCxSAE=body[ibody].vG*xSAE;
// Determination of reference velocity
double vref;
if (vCxSAE*vPxSAE>0) vref=fabs(vCxSAE);  /* braking */
    else vref=fabs(vCxSAE-vPxSAE); /* traction */
if (vref<1E-3) vref=1E-3;
// Determination of slips
double slong=vPxSAE/vref, slat=vPySAE/vref, scamber=uaxe*zSAE;
// Determination of normal force
double Fz=pen*tyre.Kz-sqrt(pen*tyre.Kz/tyre.Fznom)*tyre.Cz*vPzSAE;

if (Fz<0) return;

double lc=sqrt(8*tyre.r1*pen);
// Determination of tyre efforts according to the model developed at
// the University of Arizona by Gim Gwanghun (cf Ph.D and doc ADAMS)

double Clong,Clat,Ccamber,Stot,B1,B2,B3,fClb,Sn,ln;
Stot=sqrt(slong*slong+slat*slat);
if (Stot<1E-6) Stot=1E-6;
// Determination of friction coefficient according to total slip
double Stot2=(Stot<1 ? Stot : 1);
fClb=tyre.fClbs*(1-Stot2)+tyre.fClbd*Stot2;
// Determination of stiffness coefficients in terms of vertical load
Clong=tyre.Clongnom*pow((Fz/tyre.Fznom),tyre.nlong);
Clat=tyre.Clatnom*pow((Fz/tyre.Fznom),tyre.nlat);
Ccamber=tyre.Ccambernom*pow((Fz/tyre.Fznom),tyre.ncamber);
// Is there skid ?
int skid=((Clong*Clong*slong*slong
   +(Clat*slat+3*Ccamber*scamber)*(Clat*slat+3*Ccamber*scamber))
	       > (9*fClb*fClb*Fz*Fz));
// Determination of efforts
double Fx=0,Fy=0,My=0,Mz=0;
if (skid) // in case of skid
      {
	// cout << "skid\n";
      Fx=-fClb*slong*Fz/Stot;
      Fy=-fClb*slat*Fz/Stot;
      My=0.0;
      Mz=-0.6*fClb*slong*fClb*slat*Fz*Fz*lc/(Stot*Stot*Clat);
      }
// Computation of non dimensional length of slip area
else // no skid
      {
      B1=9*(fClb*fClb*Fz*Fz-Ccamber*Ccamber*scamber*scamber);
      B2=3*Clat*slat*Ccamber*scamber;
      B3=-(Clat*Clat*slat*slat+Clong*Clong*slong*slong);
      Sn=(B2+sqrt(B2*B2-B1*B3))/B1; ln=1-Sn;
      // Calcul des efforts dans le cas lineaire des faibles glissements
        /*  if (Sn<1E-6)
        {
        cout << "\nPetits glissements Sn=" << Sn << endl;
        Eff.Fx=-Clong*slong;
        Eff.Fy=-Clat*slat;
        Eff.My=0.0;
        Eff.Mz=-2*Clong*slong*slat*lc/3+Clat*slat*lc/6;
        return Eff;
        }*/
      // Computation of contact forces in general
      Fx=-Clong*slong*ln*ln-fClb*slong*Fz*(1+ln*ln*(2*ln-3))/Stot;
      Fy=-Clat*slat*ln*ln
	     -fClb*slat*Fz*(1+ln*ln*(2*ln-3))/Stot
	     -Ccamber*scamber*ln*ln*(3-2*ln);
      My=0.0;
      Mz=-2*Clong*slong*slat*lc*ln*ln*ln/3
         -0.6*fClb*slong*fClb*slat*Fz*Fz*
	  lc*(1+ln*ln*ln*(-10+15*ln-6*ln*ln))/(Stot*Stot*Clat);
      // Differentiation of case 3 from cases 1 and 2 (cf. doc. ADAMS)
      if((Clat*Clat*slat*slat+Ccamber*scamber*Clat*slat)>0)
        {
        Mz+=(Clat*slat*(-0.5+2*ln/3)
	     +1.5*(fClb*slat*Fz/Stot-Ccamber*scamber)*Sn*Sn)*lc*ln*ln;
        }
      else
        {
        Mz+=Clat*slat*lc*ln/6;
        }
      }
vec Ftot=Fx*xSAE+Fy*ySAE+Fz*zSAE;
body[ibody].R+=Ftot;
// The SAE coordinate system is pulled at the level of the ground for
// the application of the forces
tSAE.e.z=0;
body[ibody].MG+=My*ySAE+Mz*zSAE+((tSAE.e-body[ibody].T0G.e)^Ftot);
}

//--------------------------------------------------------------------

void AddContactPlanePointForce(double Kn, double pk,
                               double Cdamp, double pd,
                               double ffrict, double vglim,
                               int ibodyplane, vec rplane, vec nplane,
		               int ibodypoint, vec rpoint)

{
vec eplane, vplane, epoint, vpoint, nglob;

eplane=body[ibodyplane].T0G*rplane;
epoint=body[ibodypoint].T0G*rpoint;
nglob=body[ibodyplane].T0G.R*nplane;
nglob.unite();
double pen=-(epoint-eplane)*nglob;
if (pen>0)
 {
 vec eforce=epoint+0.5*pen*nglob;
 vplane=body[ibodyplane].vG
      +(body[ibodyplane].omega^(eforce-body[ibodyplane].T0G.e));
 vpoint=body[ibodypoint].vG
      +(body[ibodypoint].omega^(eforce-body[ibodypoint].T0G.e));
 vec vr=(vpoint-vplane);
 double pend=-vr*nglob;
 vec vt=vr+pend*nglob;
 double Fn=Kn*pow(pen,pk)+Cdamp*pow(pen,pd)*pend;
 if (Fn<0) Fn=0;
 vec force=Fn*nglob;
 if (vt.length()>vglim) force-=ffrict*Fn*vt/vt.length();
 else force-=ffrict*Fn*vt/vglim;
 body[ibodypoint].R+=force;
 body[ibodypoint].MG+=((eforce-body[ibodypoint].T0G.e)^force);
 body[ibodyplane].R-=force;
 body[ibodyplane].MG-=((eforce-body[ibodyplane].T0G.e)^force);
 }
}

//--------------------------------------------------------------------

void AddContactPlaneSphereForce(double Kn, double pk,
                               double Cdamp, double pd,
                               double ffrict, double vglim,
                               int ibodyplane, vec rplane, vec nplane,
			       int ibodysphere, vec rcentre, double radius)

{
vec eplane, vplane, ecentre, vpoint, nglob;

eplane=body[ibodyplane].T0G*rplane;
ecentre=body[ibodysphere].T0G*rcentre;
nglob=body[ibodyplane].T0G.R*nplane;
nglob.unite();
double pen=-(ecentre-eplane)*nglob+radius;
if (pen>0)
 {
 vec eforce=ecentre+(0.5*pen-radius)*nglob;
 vplane=body[ibodyplane].vG
      +(body[ibodyplane].omega^(eforce-body[ibodyplane].T0G.e));
 vpoint=body[ibodysphere].vG
      +(body[ibodysphere].omega^(eforce-body[ibodysphere].T0G.e));
 vec vr=(vpoint-vplane);
 double pend=-vr*nglob;
 vec vt=vr+pend*nglob;
 double Fn=Kn*pow(pen,pk)+Cdamp*pow(pen,pd)*pend;
 if (Fn<0) Fn=0;
 vec force=Fn*nglob;
 if (vt.length()>vglim) force-=ffrict*Fn*vt/vt.length();
 else force-=ffrict*Fn*vt/vglim;
 body[ibodysphere].R+=force;
 body[ibodysphere].MG+=((eforce-body[ibodysphere].T0G.e)^force);
 body[ibodyplane].R-=force;
 body[ibodyplane].MG-=((eforce-body[ibodyplane].T0G.e)^force);
 }
}

//--------------------------------------------------------------------

void AddContactPlanePointForceQB(double Kn, double pk,
                               double Cdamp, double pd,
                               double ffrict, double vglim,
                               int ibodyplane, vec rplane, vec nplane,
		               int ibodypoint, vec rpoint, double *Fn, vec *force)

{
vec eplane, vplane, epoint, vpoint, nglob;

//cout << "rPlan:" << rplane << " nplan:" << nplane << endl;

eplane=body[ibodyplane].T0G*rplane;
epoint=body[ibodypoint].T0G*rpoint;
nglob=(body[ibodyplane].T0G.R)*nplane;
nglob.unite();
double pen=-(epoint-eplane)*nglob;
if (pen>0)
 {
 vec eforce=epoint+0.5*pen*nglob;
 vplane=body[ibodyplane].vG
      +(body[ibodyplane].omega^(eforce-body[ibodyplane].T0G.e));
 vpoint=body[ibodypoint].vG
      +(body[ibodypoint].omega^(eforce-body[ibodypoint].T0G.e));
 vec vr=(vpoint-vplane);
 double pend=-vr*nglob;
 vec vt=vr+pend*nglob;
 /*double*/ (*Fn)=Kn*pow(pen,pk)+Cdamp*pow(pen,pd)*pend;
 if ((*Fn)<0) (*Fn)=0;
 /*vec*/ (*force)=(*Fn)*nglob;
 if (vt.length()>vglim) (*force)-=ffrict*(*Fn)*vt/vt.length();
 else (*force)-=ffrict*(*Fn)*vt/vglim;
 body[ibodypoint].R+=(*force);
 body[ibodypoint].MG+=((eforce-body[ibodypoint].T0G.e)^(*force));
 body[ibodyplane].R-=(*force);
 body[ibodyplane].MG-=((eforce-body[ibodyplane].T0G.e)^(*force));
 }
}

//--------------------------------------------------------------------

void AddBushingEfforts(structbushing stbs, int ibodyA, mth TrelA,
                                           int ibodyB, mth TrelB)
{
mth T0A=body[ibodyA].T0G*TrelA, T0B=body[ibodyB].T0G*TrelB;
vec emoy=0.5*(T0A.e+T0B.e);
// Calcul des positions et vitesses relatives
double thx, thy, thz;
vec axex, axey, axez;
switch(stbs.AxesOrder)
  {
  case ED_AXES_XYZ: {getrotxyz(T0A.R,T0B.R,thx,thy,thz,axex,axey,axez); break;}
  case ED_AXES_YZX: {getrotyzx(T0A.R,T0B.R,thy,thz,thx,axey,axez,axex); break;}
  case ED_AXES_ZXY: {getrotzxy(T0A.R,T0B.R,thz,thx,thy,axez,axex,axey); break;}
  case ED_AXES_XZY: {getrotxzy(T0A.R,T0B.R,thx,thz,thy,axex,axez,axey); break;}
  case ED_AXES_YXZ: {getrotyxz(T0A.R,T0B.R,thy,thx,thz,axey,axex,axez); break;}
  case ED_AXES_ZYX: {getrotzyx(T0A.R,T0B.R,thz,thy,thx,axez,axey,axex); break;}
  default: {cout << "Bushing: AxesOrder must be between 1 and 6\n"; exit(1); }
  }
vec delta=T0B.e-T0A.e;
double deltax, deltay, deltaz;
decomposevec(delta,axex,axey,axez,deltax,deltay,deltaz);
vec vitrel=VelPointBody(ibodyB,TrelB.e)-VelPointBody(ibodyA,TrelA.e);
double vitrelx, vitrely, vitrelz;
decomposevec(vitrel,axex,axey,axez,vitrelx,vitrely,vitrelz);
vec omegarel=body[ibodyB].omega-body[ibodyA].omega;
double omegarelx, omegarely, omegarelz;
decomposevec(omegarel,axex,axey,axez,omegarelx,omegarely,omegarelz);
vec force=axex*(stbs.Kx*deltax+stbs.Cx*vitrelx)
         +axey*(stbs.Ky*deltay+stbs.Cy*vitrely)
         +axez*(stbs.Kz*deltaz+stbs.Cz*vitrelz);
vec moment=axex*(stbs.KTx*thx+stbs.CTx*omegarelx)
          +axey*(stbs.KTy*thy+stbs.CTy*omegarely)
          +axez*(stbs.KTz*thz+stbs.CTz*omegarelz);
  // report des efforts sur les elements
body[ibodyA].R+=force+T0A.R*stbs.PreForceOnA;
body[ibodyB].R-=force+T0A.R*stbs.PreForceOnA;
body[ibodyA].MG+=moment+((emoy-body[ibodyA].T0G.e)^force);
body[ibodyB].MG-=moment+((emoy-body[ibodyB].T0G.e)^force);
}

//--------------------------------------------------------------------

void CreateVmoFile(scene &sc)

{
// Storage of mth in reference configuration
ComputeMotion();
mth *T0Gref;
T0Gref=new mth[nbrbody];
int ibody;
for (ibody=0; ibody<nbrbody; ibody++) T0Gref[ibody]=body[ibody].T0G;

// Resaving the reference configuration to be sure
char *FileName;
FileName=new char[strlen(application)+5];
strcpy(FileName,application);
strcat(FileName,".vol");
sc.CreateVolFile(FileName);

// Opening the mod file
strcpy(FileName,application);
strcat(FileName,".mod");
ifstream ModFile(FileName);

// Opening the vmo file
strcpy(FileName,application);
strcat(FileName,".vmo");
ofstream VmoFile(FileName);

double alpha,omega,freq,amor,massmod,*qref,*qreal,*qimag,coeff,angle,anglemax;
int iddl,itmp;
trot Rrel;

// Storing the reference position
qref=new double[nbrdof];
for (iddl=0; iddl<nbrdof; iddl++) qref[iddl]=q[iddl];
qreal=new double[nbrdof];
qimag=new double[nbrdof];
while (!ModFile.rdstate())
  {
  ModFile >> itmp >> alpha >> omega >> freq >> amor >> massmod;
  coeff=0;
  // Load the mode and determine the amplitude
  for (iddl=0; iddl<nbrdof; iddl++) if (!ModFile.rdstate())
      {
      ModFile >> qreal[iddl] >> qimag[iddl];
      coeff+=qreal[iddl]*qreal[iddl]+qimag[iddl]*qimag[iddl];
      }
  if (!ModFile.rdstate()) // data treated only if no problem after reading
    {
    // Normalize the amplitude of the mode
    for (iddl=0; iddl<nbrdof; iddl++)
      {
      qreal[iddl]*=0.01/sqrt(coeff);
      qimag[iddl]*=0.01/sqrt(coeff);
      }
    // Determine the maximum angle for real part of mode
    for (iddl=0; iddl<nbrdof; iddl++) q[iddl]=qref[iddl]+qreal[iddl];
    ComputeMotion();
    anglemax=0;
    for (ibody=0; ibody<nbrbody; ibody++)
      {
      Rrel=T0Gref[ibody].R.inv()*body[ibody].T0G.R;
      angle=2*acos(0.5*sqrt(1+Rrel.r11+Rrel.r22+Rrel.r33));
      if (angle>anglemax) anglemax=angle;
      }
    // and go on for imaginary part of mode
    for (iddl=0; iddl<nbrdof; iddl++) q[iddl]=qref[iddl]+qimag[iddl];
    ComputeMotion();
    for (ibody=0; ibody<nbrbody; ibody++)
      {
      Rrel=T0Gref[ibody].R.inv()*body[ibody].T0G.R;
      angle=2*acos(0.5*sqrt(1+Rrel.r11+Rrel.r22+Rrel.r33));
      if (angle>anglemax) anglemax=angle;
      }
    // Adapt the amplitude of the mode
    for (iddl=0; iddl<nbrdof; iddl++)
      {
      qreal[iddl]*=0.3/anglemax;
      qimag[iddl]*=0.3/anglemax;
      }
    // Store mode in VmoFile
    VmoFile << alpha << " " << omega << " " << freq << " " << amor << endl;
    for (iddl=0; iddl<nbrdof; iddl++) q[iddl]=qref[iddl]+qreal[iddl];
    ComputeMotion();
    sc.WriteCoord(VmoFile);
    for (iddl=0; iddl<nbrdof; iddl++) q[iddl]=qref[iddl]+qimag[iddl];
    ComputeMotion();
    sc.WriteCoord(VmoFile);
    }
  }
ModFile.close();
VmoFile.close();

// restore initial position
for (iddl=0; iddl<nbrdof; iddl++) q[iddl]=qref[iddl];
ComputeMotion();

delete T0Gref;
delete qref;
delete qreal;
delete qimag;
delete FileName;

}

//--------------------------------------------------------------------

double GetKineticEnergy()

{
double Energy=0;
for (int ibody=0; ibody<nbrbody; ibody++)
  {
  Energy+=0.5*body[ibody].mass*(body[ibody].vG*body[ibody].vG);
  Energy+=0.5*body[ibody].omega*(body[ibody].PhiG*body[ibody].omega);
  }
return Energy;
}

//--------------------------------------------------------------------

double GetFlexibleKineticEnergy()
// Kinetic Energy of all flexible bodies
{
double Energy=0;
    if (nbrflexbody>0)
    {
       for (int iflexbody=0; iflexbody<nbrflexbody; iflexbody++)
       {
           for(int i=0; i< nbrframesperflexbody[iflexbody];i++) //tensor line frame
           {
               for(int j=0; j< nbrframesperflexbody[iflexbody];j++) //tensor column frame
               {
                      Energy+= (flexbody[iflexbody].frame[i].vF * (flexbody[iflexbody].T0coro.R * (flexbody[iflexbody].MTT[i][j] * (flexbody[iflexbody].T0coro.R.inv() * flexbody[iflexbody].frame[j].vF)))) * 0.5;
                      Energy+= (flexbody[iflexbody].frame[i].vF * (flexbody[iflexbody].T0coro.R * (flexbody[iflexbody].MTR[i][j] * (flexbody[iflexbody].T0coro.R.inv() * flexbody[iflexbody].frame[j].omega)))) * 0.5;
                      Energy+= (flexbody[iflexbody].frame[i].omega * (flexbody[iflexbody].T0coro.R * (flexbody[iflexbody].MRT[i][j] * (flexbody[iflexbody].T0coro.R.inv() * flexbody[iflexbody].frame[j].vF)))) * 0.5;
                      Energy+= (flexbody[iflexbody].frame[i].omega * (flexbody[iflexbody].T0coro.R * (flexbody[iflexbody].MRR[i][j] * (flexbody[iflexbody].T0coro.R.inv() * flexbody[iflexbody].frame[j].omega)))) * 0.5;
               }
           }
       }
    }
return Energy;
}

//--------------------------------------------------------------------

double GetGravitationalEnergy(vec grav)

{
double Energy=0;
for (int ibody=0; ibody<nbrbody; ibody++)
  {
  Energy-=body[ibody].mass*(body[ibody].T0G.e*grav);
  }
return Energy;
}

//--------------------------------------------------------------------
