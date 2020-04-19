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

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

#include <EasyDyn/vec.h>
#include <EasyDyn/visu.h>

struct structbody
{
double mass;
tiner PhiG;
mth T0G,  // homogeneous transformation matrix giving the position
            //of the body with respect to the global refernce frame
    TrefG; // idem for the relative position
vec vG,     // velocity of center of gravity
    aG,     // acceleration of center of gravity
    omega,  // rotation velocity of body
    omegad, // rotational acceleration of body
    // Relative motion terms, used in ComposeKinematics() and
    // assumed to be expressed in the coordinate system used as reference
    // for the relative motion
    vGrel,  // relative velocity of center of gravity
    aGrel,  // relative acceleration of center of gravity
    omegarel,  // relative rotation velocity of body
    omegadrel, // relative rotational acceleration of body
    // vectors used for forces
    R, MG,  // resultant force and moment at center of gravity
    vGsto, omegasto, // used for numerical derivation
    *vGpartial, *omegapartial,*vGrelpartial, *omegarelpartial;
int *underdof;
};

struct structframe
{
mth T0F, // homogeneous transformation matrices giving the situation
         // of the frame with respect to the global reference frame
  TrefF, // idem for the relative situation
  TcoroF0, // mth giving the undeformed situation of frame wrt to corot frame
           // which is necessary to compute elastic forces
  TcoroF0_res; //resultant mth after rigid displacement from the corotational frame
vec vF,     // velocities of frame
    aF,     // acceleration of frame
    omega,  // rotation velocity of frame
    omegad, // rotational acceleration of frame
    deform, theta, //strains and angular deformation to compute elastic forces
    Vdeform, Vtheta, //strains and angular deformation VELOCITY to compute elastic Damping
    // Relative motion terms, used in ComposeKinematics() and
    // assumed to be expressed in the coordinate system used as reference
    // for the relative motion
    vFrel,  // relative velocity of frame
    aFrel,  // relative acceleration of frame
    omegarel,  // relative rotation velocity of frame
    omegadrel, // relative rotational acceleration of frame
    // vectors used for forces
    R, MG,  // resultant force and moment on frame
    vFsto, omegasto, // used only for numerical derivation
    *vFpartial, *omegapartial,*vFrelpartial, *omegarelpartial;
    int *underdof;
};

struct structflexbody
{
int nbrframe;

double length; // length of the beam
double rho;   // material density of the beam
double E;     // Young's modulus of the beam
double nu;    // Poisson's ratio of the beam
double section; // section of the beam
double Iyy, Izz; // moments of inertia of the beam

tgen **MTT, **MTR, **MRT, **MRR;
tgen **KTT, **KTR, **KRT, **KRR;
tgen **CTT, **CTR, **CRT, **CRR;
structframe *frame;
mth T0coro, // homogeneous transformation matrix giving the situation
           // of the corot frame with respect to the global reference frame
    T0cororef; // relative situation of the corotational frame
vec omegacoro, // velocity of corotational frame
    omegacororel, // relative rotation velocity of corotational frame
    *omegacoropartial, // partial velocities of the corotational frame
    *omegacororelpartial, // relative partial velocities of the corotational frame
    omegacorosto, // used for numerical derivation
    MG; // resultant moment on the corotational frame
int *corounderdof;
};

struct structtyre
{
double r1; // outer radius
double r2; // equivalent torus radius
double Kz; // vertical stiffness
double Cz; // vertical damping
double Fznom; // nominal vertical force on tyre
double Clongnom; // nominal longitudinal stiffness
double nlong; // Clong=Clongnom*(Fz/Fznom)^nlong
double Clatnom; // nominal cornering stiffness
double nlat; // Clat=Clatnom*(Fz/Fznom)^nlat
double Ccambernom; // nominal camber stiffness
double ncamber; // Ccamber=Ccambernom*(Fz/Fznom)^ncamber
double fClbs, fClbd; // Coulomb friction coefficients (static and dynamic)
};

struct structbushing
{
double Kx, Ky, Kz;    // Translational stiffnesses
double KTx, KTy, KTz; // Rotational stiffnesses
double Cx, Cy, Cz;    // Translational stiffnesses
double CTx, CTy, CTz; // Rotational damping coefficients
vec PreForceOnA;
int AxesOrder;
};

const int ED_AXES_XYZ=1, ED_AXES_YZX=2, ED_AXES_ZXY=3,
          ED_AXES_XZY=4, ED_AXES_YXZ=5, ED_AXES_ZYX=6;

#ifdef EASYDYNMBSMAIN
#define EASYDYNSIMMAIN
structbody *body=0;
structflexbody *flexbody=0;
double *p,*pd,*pdd;
int nbrbody=0,nbrdep=0,nbrflexbody=0,*nbrframesperflexbody;
int mbsAdvanced=0;
vec flexgravity=vcoord(0.0,0.0,0.0);
char filename[100];
double AlphaDamp, BetaDamp;
#else
extern structbody *body;
extern structflexbody *flexbody;
extern double *p,*pd,*pdd;
extern int nbrbody,nbrdep,nbrflexbody,*nbrframesperflexbody;
extern int mbsAdvanced;
extern vec flexgravity;
extern char* filename;
extern double AlphaDamp, BetaDamp;
#endif

#include <EasyDyn/sim.h>

//--------------------------------------------------------------------
// Procedures that must be provided by the user

void SetInertiaData();
void ComputeMotion();
void ComputePartialVelocities();
void AddAppliedEfforts();
void SaveData(ostream &OutFile); // eventually empty
void WriteDataHeader(ostream &OutFile); // eventually empty
void AssemblingBeamMK(); // assembling matrix element for multiple frames beam representation (along the X axis)

//--------------------------------------------------------------------
// Procedures that we provide

void InitEasyDynmbs();
void EndEasyDynmbs();
void ComputeInertiaEfforts();
void ComputePartialVelocitiesDefault();
#ifdef EASYDYNMBSMAIN
#ifndef EASYDYNMBSADVANCED
void ComputePartialVelocities() { ComputePartialVelocitiesDefault(); }
#endif
#endif
void CreateBeamMassMatrix(int iflexbody, double LL, double AA, double rho, double Iyy, double Izz); // Create the mass matrix of a 3D beam along the X direction
void CreateBeamStiffnessMatrix(int iflexbody, double LL, double AA, double EE, double Iyy, double Izz, double nu); // Create the stiffness matrix of a 3D beam along the X direction
gsl_matrix *mbeam3d(int iflexbody, int initframe, int endframe, double AA, double rho, double Ixx, double x3, double y3, double z3); // Create the mass matrix of a 3D beam with its rotation
gsl_matrix *kbeam3d(int iflexbody, int initframe, int endframe, double E, double nu, double AA, double Ixx, double Iyy, double Izz, double x3, double y3, double z3); // Create the stiffness matrix of a 3D beam with its rotation
void connect(int iflexbody, gsl_matrix *MassElement, gsl_matrix *StiffElement, int q0, int q1, int q2, int q3, int q4, int q5, int q6, int q7, int q8, int q9, int q10, int q11); // Assemble the mbeam3d and kneam3d matrices
void ReadMassMatrix(int iflexbody, char *name=NULL); // Read a mass matrix from a .mm file
void ReadStiffnessMatrix(int iflexbody, char *name=NULL); // Read a stiffness matrix from a .kk file
void DisplayMassMatrix(int iflexbody); // Display the mass matrix of iflexbody
void DisplayStiffnessMatrix(int iflexbody); // Display the stiffness matrix of iflexbody
void DisplayDampingMatrix(int iflexbody); // Display the damping matrix of iflexbody

void ComputeResidualmbs();
void ComputeMassMatrix(int MotionToCompute);
void ComputeAccelerations();
void ComposeMotionOld(int ibody, int ibodyref); //Limited to rigid bodies
void ComposeMotion(int ibody, int ibodyref, int iframe=-1, int iframeref=-1); //Extended to flexible bodies

void ComposePartialVelocitiesOld(int ibody, int ibodyref); //Limited to rigid bodies
void ComposePartialVelocities(int ibody, int ibodyref, int iframe=-1, int iframeref=-1); //Extended to flexible bodies

vec VelPointBody(int ibody, vec r);
vec AccPointBody(int ibody, vec r);
void AddGravityForces(vec grav);
void AddSpringForce(double K, double L0,
                    int ibodyA, vec rA, int ibodyB, vec rB);
void AddSpringForceQB(double K, double L0,
                    int ibodyA, vec rA, int ibodyB, vec rB, vec *force);
void AddDamperForce(double C,
                    int ibodyA, vec rA, int ibodyB, vec rB);
void AddDamperForceQB(double C,
                    int ibodyA, vec rA, int ibodyB, vec rB, vec *force);
void AddContactPlanePointForce(double Kn, double pk,
                               double Cdamp, double pd,
                               double ffrict, double vglim,
                               int ibodyplane, vec rplane, vec nplane,
		               int ibodypoint, vec rpoint);
void AddContactPlanePointForceQB(double Kn, double pk,
                               double Cdamp, double pd,
                               double ffrict, double vglim,
                               int ibodyplane, vec rplane, vec nplane,
		               int ibodypoint, vec rpoint,
			       double *Fn, vec *force);
void AddContactPlaneSphereForce(double Kn, double pk, double Cdamp, double pd,
                               double ffrict, double vglim,
                               int ibodyplane, vec rplane, vec nplane,
			       int ibodysphere, vec rcentre, double radius);
void AddTyreEfforts(int ibody, vec axe, structtyre tyre);
void AddBushingEfforts(structbushing stbs,
                       int ibodyA, mth TrelA,
                       int ibodyB, mth TrelB);
void CreateVmoFile(scene &sc);
double GetKineticEnergy();
double GetFlexibleKineticEnergy();  //Kinetic Energy of all flexible bodies
double GetGravitationalEnergy(vec grav);

//--------------------------------------------------------------------
