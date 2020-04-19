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

#include <string.h>
#include <math.h>
#include <fstream>
#include <iostream>
#include <iomanip>
//#include <time.h>

#ifdef __unix__
#include <sys/time.h> // pour gettimeofday
#endif

#ifdef _WIN32
  #include <time.h>
  #include <windows.h>
  #if defined(_MSC_VER) || defined(_MSC_EXTENSIONS)
    #define DELTA_EPOCH_IN_MICROSECS  11644473600000000Ui64
  #else
    #define DELTA_EPOCH_IN_MICROSECS  11644473600000000ULL
  #endif
#endif

using namespace std;

#include <gsl/gsl_math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_blas.h>

#include <EasyDyn/sim.h>

gsl_vector_view fview,qview,qdview,qddview,uview,yview;
gsl_vector *qsto=0,*qdsto=0,*qddsto=0,*qddsto2=0,*qddcor=0;
gsl_permutation *perm;
gsl_matrix *jac=0;
ofstream DbgFile;

//---------------------------------------------------------------------------

// gettimeofday replacement for windows

#ifdef _WIN32

int gettimeofday(struct timeval *tv, void *tz)
{
  FILETIME ft;
  unsigned __int64 tmpres = 0;
  static int tzflag = 0;

  if (NULL != tv)
  {
    GetSystemTimeAsFileTime(&ft);

    tmpres |= ft.dwHighDateTime;
    tmpres <<= 32;
    tmpres |= ft.dwLowDateTime;

    tmpres /= 10;  /*convert into microseconds*/
    /*converting file time to unix epoch*/
    tmpres -= DELTA_EPOCH_IN_MICROSECS;
    tv->tv_sec = (long)(tmpres / 1000000UL);
    tv->tv_usec = (long)(tmpres % 1000000UL);
  }

  return 0;
}

#endif

//---------------------------------------------------------------------------

void EasyDynsimError(const char *texte)
{
 cout << "EasyDynsim ERROR: " << texte << "\n";
 exit(1);
}

//---------------------------------------------------------------------------

void InitEasyDynsim()

{
if (nbrdof==0)
  EasyDynsimError("nbrdof must be initialized before calling InitEasyDynSim");
if (application==0)
  EasyDynsimError("application must be initialized before calling InitEasyDynSim");
f=new double[nbrdof];
x=new double[2*nbrdof];
q=x;
qd=x+nbrdof;
xd=new double[2*nbrdof];
qdd=xd+nbrdof;
if (nbrinput>0) u=new double[nbrinput];
if (nbroutput>0) y=new double[nbroutput];
fview=gsl_vector_view_array(f,nbrdof);
qview=gsl_vector_view_array(q,nbrdof);
qdview=gsl_vector_view_array(qd,nbrdof);
qddview=gsl_vector_view_array(qdd,nbrdof);
if (nbrinput>0) uview=gsl_vector_view_array(u,nbrinput);
if (nbroutput>0) yview=gsl_vector_view_array(y,nbroutput);
qsto=gsl_vector_alloc(nbrdof);
qdsto=gsl_vector_alloc(nbrdof);
qddsto=gsl_vector_alloc(nbrdof);
qddsto2=gsl_vector_alloc(nbrdof);
qddcor=gsl_vector_alloc(nbrdof);
perm = gsl_permutation_alloc(nbrdof);
jac=gsl_matrix_calloc(nbrdof,nbrdof);

 int idof,iinput,ioutput;
for (idof=0; idof<nbrdof; idof++)
   {
     q[idof]=0.0; qd[idof]=0.0; qdd[idof]=0.0; xd[idof]=0.0;
   }
for (iinput=0; iinput<nbrinput; iinput++) u[iinput]=0;
for (ioutput=0; ioutput<nbroutput; ioutput++) y[ioutput]=0;
}

//---------------------------------------------------------------------------

void EndEasyDynsim()

{
delete f; delete x; delete xd; if (nbrinput>0) delete u;
gsl_vector_free(qsto);
gsl_vector_free(qdsto);
gsl_vector_free(qddsto);
gsl_vector_free(qddsto2);
gsl_vector_free(qddcor);
gsl_permutation_free(perm);
gsl_matrix_free(jac);
}

//---------------------------------------------------------------------------

void WriteStateVariablesHeader(ostream &OutFile)

{
int idof,iinput;
OutFile << " time ";
for (idof=0;idof<nbrdof; idof++)
  OutFile << " q" << idof << " qd" << idof << " qdd" << idof;
for (iinput=0;iinput<nbrinput; iinput++)
  OutFile << " u" << iinput;
OutFile << " ";
}

//--------------------------------------------------------------------

void SaveStateVariables(ostream &OutFile)

{
int idof,iinput;
OutFile << " " << t;
for (idof=0;idof<nbrdof; idof++)
  OutFile << " " << q[idof] << " " << qd[idof] << " " << qdd[idof];
for (iinput=0; iinput<nbrinput; iinput++)
  OutFile << " " << u[iinput];
OutFile << " ";
}

//--------------------------------------------------------------------

void ComputeResidual2(int *doflocked)

{
int idof;
ComputeResidual();
if (doflocked) for (idof=0; idof<nbrdof; idof++)
   if (doflocked[idof]) f[idof]=0.0;
}

//--------------------------------------------------------------------

void Calcule_Matrice_J(double alphaM, double alphaC, double alphaK, double tol,
                       int *doflocked)

{
  int iddl,jddl;
  double *fsto;
  fsto=new double[nbrdof];
  ComputeResidual2(doflocked);
  for (iddl=0; iddl<nbrdof; iddl++) fsto[iddl]=f[iddl];
  for (jddl=0; jddl<nbrdof; jddl++)
    {
    q[jddl]+=alphaK*tol;
    qd[jddl]+=alphaC*tol;
    qdd[jddl]+=alphaM*tol;
    if (doflocked) { if (!doflocked[jddl]) ComputeResidual2(doflocked); }
    else ComputeResidual2(doflocked);
     for (iddl=0; iddl<nbrdof; iddl++)
     {
      gsl_matrix_set(jac,iddl,jddl,(f[iddl]-fsto[iddl])/tol);
     }
    q[jddl]-=alphaK*tol;
    qd[jddl]-=alphaC*tol;
    qdd[jddl]-=alphaM*tol;
    }
  for (iddl=0; iddl<nbrdof; iddl++) f[iddl]=fsto[iddl];
  delete fsto;
  if (doflocked) for (iddl=0; iddl<nbrdof; iddl++) if (doflocked[iddl])
    {
    for (jddl=0; jddl<nbrdof; jddl++)
      {
      gsl_matrix_set(jac,iddl,jddl,0);
      gsl_matrix_set(jac,jddl,iddl,0);
      }
    gsl_matrix_set(jac,iddl,iddl,1);
    }
}

//---------------------------------------------------------------------------

void StaticEquilibrium(int *doflocked)

{
  int nstep,s;
  double err;
  nstep=0;
  cout << "Equilibrium iterations: ";
  do
    {
    nstep++;
    cout << "+";
    if ((nstep%10)==1)
      {
	Calcule_Matrice_J(100,0,1,1E-5,doflocked);
	gsl_linalg_LU_decomp(jac,perm,&s);
      }
    ComputeResidual2(doflocked);
    gsl_linalg_LU_solve(jac,perm,&fview.vector,qddcor);
    gsl_vector_sub(&qview.vector,qddcor);
    err=gsl_blas_dnrm2(qddcor);
    if (DEBUG) cout << "Err=" << err << "\n";
    }
  while ((nstep<1000) && (err>(1E-8*sqrt(1.0*nbrdof))));
  if(nstep==1000) EasyDynsimError("No convergence after 1000 iterations");
  cout << " Bingo" << "\n";
}

//---------------------------------------------------------------------------

int NewmarkOneStep(double h, double &errq, int *doflocked, int hmodified)
  {
  // Let's control the user assiduity
  if ((f==0) || (q==0) || (qd==0) || (qdd==0) || (qsto==0) || (qdsto==0)
      || (qddsto==0) || (qddsto2==0) || (qddcor==0) || (jac==0))
    EasyDynsimError("NewmarkOneStep: Newmark variables not allocated");
  int iddl,nstep,istepjac,s;
  double err,Beta=0.25,Gamma=0.5;
  // Copie de reference pour les accelerations pour calculer l'erreur
  //for(int i=0 ; i<nbrdof ; i++) cout<< "qddview "<<i<<" "<<gsl_vector_get(&qddview.vector,i)<<endl;
  gsl_vector_memcpy(qddsto2,&qddview.vector);
  // Prediction
  t+=h;
  for (iddl=0; iddl<nbrdof; iddl++)
      {
      q[iddl]+=h*qd[iddl]+0.5*h*h*qdd[iddl];
      qd[iddl]+=h*qdd[iddl];
      }

  if (DEBUG) cout << "Newmark iterations: ";
  // Correction
  istepjac=0;
  if (hmodified) istepjac=1;
  nstep=0;
  do
    {
    nstep++;
    if (DEBUG) cout << "+";
    if ((nstep%4)==istepjac)
      {
	Calcule_Matrice_J(1,Gamma*h,Beta*h*h,1E-2,doflocked);
        gsl_linalg_LU_decomp(jac,perm,&s);
      }

    ComputeResidual2(doflocked);

    //gsl_vector_fprintf(stdout,&residu.vector,"%g");

    /*for(int i=0; i<nbrdof; i++)
    {
         for(int j=0; j<nbrdof; j++)
         {
                 cout<<gsl_matrix_get(jac,i,j)<<" " ;
         }
         cout<< " " <<endl;
    }
    system("pause");*/

    gsl_linalg_LU_solve(jac,perm,&fview.vector,qddcor);
    //err=gsl_blas_dnrm2(qddcor)/sqrt(1.0*nbrdof);
    err=gsl_blas_dnrm2(qddcor)/(sqrt(1.0*nbrdof)*(1+gsl_blas_dnrm2(&qddview.vector)));
    gsl_vector_sub(&qddview.vector,qddcor);
    gsl_vector_scale(qddcor,Gamma*h);
    gsl_vector_sub(&qdview.vector,qddcor);
    gsl_vector_scale(qddcor,Beta*h/Gamma);
    gsl_vector_sub(&qview.vector,qddcor);
    //gsl_vector_fprintf(stdout,&qdcopy.vector,"%g");
    //err=0;
    //cout << "Err=" << err << "\n";
    }
  while ((nstep<100) && (err>(sqrt(AbsTol*RelTol)/100.0)));
  if(nstep==100) return 1;
  if (DEBUG) cout << " Bingo" << "\n";
  // Calcul de l'erreur a partir de la variation du vecteur acceleration
  gsl_vector_sub(qddsto2,&qddview.vector);
  errq=0;
  for (iddl=0; iddl<nbrdof; iddl++)
      {
      double errtmp, absqiddl;
      errtmp=h*h*fabs(gsl_vector_get(qddsto2,iddl))/12.0;
      absqiddl=fabs(q[iddl]);
      if (fabs(qd[iddl]*h)>absqiddl) absqiddl=h*fabs(qd[iddl]);
      if (fabs(qdd[iddl]*h*h)>absqiddl) absqiddl=h*h*fabs(qdd[iddl]);
      if ((absqiddl*RelTol)>AbsTol) errtmp=errtmp*AbsTol/(absqiddl*RelTol);
      errq+=errtmp*errtmp;
      }
  errq=sqrt(errq/nbrdof);

  if (!gsl_finite(errq)) return 2;
  return 0;
}
//---------------------------------------------------------------------------

void NewmarkPrediction(double h, double &errq)
{
  // Let's control the user assiduity
  if ((f==0) || (q==0) || (qd==0) || (qdd==0) || (qsto==0) || (qdsto==0)
      || (qddsto==0) || (qddsto2==0) || (qddcor==0) || (jac==0))
    EasyDynsimError("NewmarkOneStep: Newmark variables not allocated");
  // Copie de reference pour les accelerations pour calculer l'erreur
  //for(int i=0 ; i<nbrdof ; i++) cout<< "qddview "<<i<<" "<<gsl_vector_get(&qddview.vector,i)<<endl;
  gsl_vector_memcpy(qddsto2,&qddview.vector);
  // Prediction
  t+=h;
  for (int iddl=0; iddl<nbrdof; iddl++)
      {
      q[iddl]+=h*qd[iddl]+0.5*h*h*qdd[iddl];
      qd[iddl]+=h*qdd[iddl];
      }

  if (DEBUG) cout << "Newmark iterations: ";
}

//---------------------------------------------------------------------------

int NewmarkCorrection(double h, double &errq, int *doflocked, int hmodified)
{
   // Let's control the user assiduity
  if ((f==0) || (q==0) || (qd==0) || (qdd==0) || (qsto==0) || (qdsto==0)
      || (qddsto==0) || (qddsto2==0) || (qddcor==0) || (jac==0))
    EasyDynsimError("NewmarkOneStep: Newmark variables not allocated");
  int iddl,nstep,istepjac,s;
  double err,Beta=0.25,Gamma=0.5;
  // Correction
  istepjac=0;
  if (hmodified) istepjac=1;
  nstep=0;
  do
    {
    nstep++;
    if (DEBUG) cout << "+";
    if ((nstep%4)==istepjac)
      {
	Calcule_Matrice_J(1,Gamma*h,Beta*h*h,1E-2,doflocked);
        gsl_linalg_LU_decomp(jac,perm,&s);
      }

    ComputeResidual2(doflocked);

    //gsl_vector_fprintf(stdout,&residu.vector,"%g");

    gsl_linalg_LU_solve(jac,perm,&fview.vector,qddcor);
    //err=gsl_blas_dnrm2(qddcor)/sqrt(1.0*nbrdof);
    err=gsl_blas_dnrm2(qddcor)/(sqrt(1.0*nbrdof)*(1+gsl_blas_dnrm2(&qddview.vector)));
    gsl_vector_sub(&qddview.vector,qddcor);
    gsl_vector_scale(qddcor,Gamma*h);
    gsl_vector_sub(&qdview.vector,qddcor);
    gsl_vector_scale(qddcor,Beta*h/Gamma);
    gsl_vector_sub(&qview.vector,qddcor);
    //gsl_vector_fprintf(stdout,&qdcopy.vector,"%g");
    //err=0;
    //cout << "Err=" << err << "\n";
    }
  while ((nstep<100) && (err>(sqrt(AbsTol*RelTol)/100.0)));
  if(nstep==100)
           return 1;

  if (DEBUG) cout << " Bingo" << "\n";
  // Calcul de l'erreur a partir de la variation du vecteur acceleration
  gsl_vector_sub(qddsto2,&qddview.vector);
  errq=0;
  for (iddl=0; iddl<nbrdof; iddl++)
      {
      double errtmp, absqiddl;
      errtmp=h*h*fabs(gsl_vector_get(qddsto2,iddl))/12.0;
      absqiddl=fabs(q[iddl]);
      if (fabs(qd[iddl]*h)>absqiddl) absqiddl=h*fabs(qd[iddl]);
      if (fabs(qdd[iddl]*h*h)>absqiddl) absqiddl=h*h*fabs(qdd[iddl]);
      if ((absqiddl*RelTol)>AbsTol) errtmp=errtmp*AbsTol/(absqiddl*RelTol);
      errq+=errtmp*errtmp;
      }
  errq=sqrt(errq/nbrdof);

  if (!gsl_finite(errq)) return 2;
  return 0;
}
//---------------------------------------------------------------------------

void NewmarkInterval(double tfinal, double &h, double hmax, int *doflocked)
// Procede a l'integration jusque tfinal

{
  // Let's control the user assiduity
  if ((f==0) || (q==0) || (qd==0) || (qdd==0) || (qsto==0) || (qdsto==0)
      || (qddsto==0) || (qddsto2==0) || (qddcor==0) || (jac==0))
    EasyDynsimError("NewmarkInterval: Probably InitMbsSim was not called");
  int iddl, istep, code, hchanged=1;
  double errq,timesto;
  // Time step verification
  if (h>hmax) h=hmax;
  // Copie de reference pour les accelerations
  while (t<tfinal)
    {
    if ((t+1.4*h)>=tfinal) { h=tfinal-t; hchanged=1; }
    timesto=t;
    gsl_vector_memcpy(qsto,&qview.vector);
    gsl_vector_memcpy(qdsto,&qdview.vector);
    gsl_vector_memcpy(qddsto,&qddview.vector);
    code=NewmarkOneStep(h,errq,doflocked,hchanged);
    hchanged=0;
    if (DEBUG) cout << "time:" << t << " h=" << h
                    << " Erreur sur les positions: " << errq << "\n";
    if ((code) || (errq>AbsTol))  // la solution n'est pas convaincante
      {
        if (DEBUG)
          {
	    char chartmp;
          if (code==1) cout << "\nNo convergence - ";
          if (code==2) cout << "\nNumerical trouble - ";
          if (gsl_finite(errq)) if (errq>AbsTol)
			  cout << "\nLack of accuracy - ";
          //cin >> chartmp;
          }
        if ((code==1) || (code==2) || (!gsl_finite(errq))) h*=0.25;
        // The step is divided by a factor between 2 and 5
        // else h*=sqrt((0.21*AbsTol+0.04*errq)/errq);
        else
          {
            // the step is divided by a power of 2
	    int n2=ceil(log(sqrt(errq/AbsTol))/log(2.0));
            for (int in2=0; in2<n2; in2++) h=h/2.0;
	  }
        hchanged=1;
        if (DEBUG) cout << "Reducing time step to " << h << "\n";
        t=timesto;
        gsl_vector_memcpy(&qview.vector,qsto);
        gsl_vector_memcpy(&qdview.vector,qdsto);
        gsl_vector_memcpy(&qddview.vector,qddsto);
        if (h<hminmin)
          EasyDynsimError("\nNo convergence with the smallest time step - STOP");
      }
    else
      {
      // Let's save the information
      if (DEBUG) DbgFile << t << " " << h << " " << errq << "\n";
      // The time step is increased if the error is small
      if ((errq<0.25*AbsTol) && (h<hmax))
        {
        if (DEBUG) cout << "\nMmh.. Good ... Increasing time step\n";
        // the step is multiplied by a factor between 2 and 5
        // h*=sqrt(AbsTol/(2.1*errq+0.04*AbsTol));
        // the step is divided by a power of 2
	int n2=floor(log(sqrt(AbsTol/errq))/log(2.0));
        for (int in2=0; in2<n2; in2++) h=h*2.0;
        if (h>hmax) h=hmax;
        hchanged=1;
        }
      }
    }
}

//---------------------------------------------------------------------------

void NewmarkIntervalold(double tfinal, double &h, double hmax, int *doflocked)
// NewmarkInterval with time step management of EasyDyn-1.2.3
// Procede a l'integration jusque tfinal

{
  // Let's control the user assiduity
  if ((f==0) || (q==0) || (qd==0) || (qdd==0) || (qsto==0) || (qdsto==0)
      || (qddsto==0) || (qddsto2==0) || (qddcor==0) || (jac==0))
    EasyDynsimError("NewmarkInterval: Probably InitMbsSim was not called");
  int iddl, istep, code, hchanged=1;
  double errq,timesto;
  // Time step verification
  if (h>hmax) h=hmax;
  // Copie de reference pour les accelerations
  while (t<tfinal)
    {
    if ((t+1.4*h)>=tfinal) { h=tfinal-t; hchanged=1; }
    timesto=t;
    gsl_vector_memcpy(qsto,&qview.vector);
    gsl_vector_memcpy(qdsto,&qdview.vector);
    gsl_vector_memcpy(qddsto,&qddview.vector);
    code=NewmarkOneStep(h,errq,doflocked,hchanged);
    errq=errq/h; // So as to get back errqd
    hchanged=0;
    if (DEBUG) cout << "time:" << t << " h=" << h
                    << " Erreur sur les positions: " << errq << "\n";
    if ((code) || (errq>AbsTol))  // la solution n'est pas convaincante
      {
        if (DEBUG)
          {
	    char chartmp;
          if (code==1) cout << "\nNo convergence - ";
          if (code==2) cout << "\nNumerical trouble - ";
          if (gsl_finite(errq)) if (errq>AbsTol)
			  cout << "\nLack of accuracy - ";
          //cin >> chartmp;
          }
        if ((code==1) || (code==2) || (!gsl_finite(errq))) h*=0.25;
        // The step is divided by a factor between 2 and 5
        // else h*=sqrt((0.21*AbsTol+0.04*errq)/errq);
        else
          {
            // the step is divided by a power of 2
	    int n2=ceil(log(errq/AbsTol)/log(2.0));
            for (int in2=0; in2<n2; in2++) h=h/2.0;
	  }
        hchanged=1;
        if (DEBUG) cout << "Reducing time step to " << h << "\n";
        t=timesto;
        gsl_vector_memcpy(&qview.vector,qsto);
        gsl_vector_memcpy(&qdview.vector,qdsto);
        gsl_vector_memcpy(&qddview.vector,qddsto);
        if (h<hminmin)
          EasyDynsimError("\nNo convergence with the smallest time step - STOP");
      }
    else
      {
      // Let's save the information
      if (DEBUG) DbgFile << t << " " << h << " " << errq << "\n";
      // The time step is increased if the error is small
      if ((errq<(0.5*AbsTol)) && (h<hmax))
        {
        if (DEBUG) cout << "\nMmh.. Good ... Increasing time step\n";
        // the step is multiplied by a factor between 2 and 5
        // h*=sqrt(AbsTol/(2.1*errq+0.04*AbsTol));
        // the step is divided by a power of 2
	int n2=floor(log(AbsTol/errq)/log(2.0));
        for (int in2=0; in2<n2; in2++) h=h*2.0;
        if (h>hmax) h=hmax;
        hchanged=1;
        }
      }
    }
}

//---------------------------------------------------------------------------

void NewmarkIntegration(double tfinal, double hsave, double hmax,
			int *doflocked, int old)
// Procede a l'integration complete sur des intervalles successifs
// correspondant au pas de sauvegarde

{
  timeval tim;
  gettimeofday(&tim, NULL);
  double t1=tim.tv_sec+(tim.tv_usec/1000000.0);
  // Let's control the user assiduity
  if ((f==0) || (q==0) || (qd==0) || (qdd==0) || (qsto==0) || (qdsto==0)
      || (qddsto==0) || (qddsto2==0) || (qddcor==0) || (jac==0))
    EasyDynsimError("NewmarkIntegration: Probably InitMbsSim was not called");
  // Ouverture du fichier de resultats et ecriture de la ligne d'entete
  char *nom;
  nom=new char[strlen(application)+5];
  strcpy(nom,application);
  strcat(nom,".res");
  cout << "Opening file " << nom << "\n" << "\n";
  ofstream ResFile(nom);
  if (!ResFile) EasyDynsimError("Cannot open file");
  WriteDataHeader(ResFile);
  // Ouverture du fichier DEBUG
  if (DEBUG)
      {
      strcpy(nom,application);
      strcat(nom,".dbg");
      cout << "Opening file " << nom << "\n" << "\n";
      DbgFile.open(nom);
      if (!DbgFile) EasyDynsimError("Cannot open file");
      }
  // Calcul des accelerations initiales et sauvegarde
  t=0;
  double h=0,dummy;
  NewmarkOneStep(0,dummy,doflocked);
  cout << "Sauvegarde\n";
  SaveData(ResFile);
  // Integration proprement dite
  int ipas=0;
  h=hmax;
  while (t<tfinal)
    {
      ipas++;
      cout << "Process up to time  " << hsave*ipas << " s" ;
      if (!old) NewmarkInterval(hsave*ipas,h,hmax,doflocked);
      if (old) NewmarkIntervalold(hsave*ipas,h,hmax,doflocked);
      cout << " - Finished - Saving ";
      SaveData(ResFile);
      cout << " - Done" << "\n";
    }
  // Fermeture des fichiers
  ResFile.close();
  if (DEBUG) DbgFile.close();
  gettimeofday(&tim, NULL);
  double t2=tim.tv_sec+(tim.tv_usec/1000000.0);
  cout << "Time needed for integration: " << t2-t1 << " seconds" << endl;
}

//---------------------------------------------------------------------------

//Original NewmarkIntegration() routine without any "cout<< "
void NewmarkIntegration2(double tfinal, double hsave, double hmax,
			int *doflocked, int old)
// Procede a l'integration complete sur des intervalles successifs
// correspondant au pas de sauvegarde

{
  timeval tim;
  gettimeofday(&tim, NULL);
  double t1=tim.tv_sec+(tim.tv_usec/1000000.0);
  // Let's control the user assiduity
  if ((f==0) || (q==0) || (qd==0) || (qdd==0) || (qsto==0) || (qdsto==0)
      || (qddsto==0) || (qddsto2==0) || (qddcor==0) || (jac==0))
    EasyDynsimError("NewmarkIntegration: Probably InitMbsSim was not called");
  // Ouverture du fichier de resultats et ecriture de la ligne d'entete
  char *nom;
  nom=new char[strlen(application)+5];
  strcpy(nom,application);
  strcat(nom,".res");
  //cout << "Opening file " << nom << "\n" << "\n";
  ofstream ResFile(nom);
  if (!ResFile) EasyDynsimError("Cannot open file");
  WriteDataHeader(ResFile);
  // Ouverture du fichier DEBUG
  if (DEBUG)
      {
      strcpy(nom,application);
      strcat(nom,".dbg");
      cout << "Opening file " << nom << "\n" << "\n";
      DbgFile.open(nom);
      if (!DbgFile) EasyDynsimError("Cannot open file");
      }
  // Calcul des accelerations initiales et sauvegarde
  t=0;
  double h=0,dummy;
  NewmarkOneStep(0,dummy,doflocked);
  //cout << "Sauvegarde\n";
  SaveData(ResFile);
  // Integration proprement dite
  int ipas=0;
  h=hmax;
  while (t<tfinal)
    {
      ipas++;
      //cout << "Process up to time  " << hsave*ipas << " s" ;
      if (!old) NewmarkInterval(hsave*ipas,h,hmax,doflocked);
      if (old) NewmarkIntervalold(hsave*ipas,h,hmax,doflocked);
      //cout << " - Finished - Saving ";
      SaveData(ResFile);
      //cout << " - Done" << "\n";
    }
  // Fermeture des fichiers
  ResFile.close();
  if (DEBUG) DbgFile.close();
  gettimeofday(&tim, NULL);
  double t2=tim.tv_sec+(tim.tv_usec/1000000.0);
  //cout << "Time needed for integration: " << t2-t1 << " seconds" << endl;
}

//---------------------------------------------------------------------------

void SetIntegrationPrecision(double NewAbsTol, double NewRelTol)
{
  AbsTol=NewAbsTol;
  RelTol=NewRelTol;
}

//---------------------------------------------------------------------------

void SaveLinearizedSystem(int *doflocked)

// Saves mass, damping and stiffness matrices and eventually
// inputs contribution

{
  int idof,jdof;
  char *nom;
  nom=new char[strlen(application)+5];

  cout << "Constructing and saving stiffness matrix ...";
  Calcule_Matrice_J(0,0,1,1E-8,doflocked);
  strcpy(nom,application);
  strcat(nom,".kk");
  ofstream stifffile(nom);
  if (!stifffile) EasyDynsimError("Cannot create stiffness matrix file");
  for (idof=0; idof<nbrdof; idof++)
     {
     for (jdof=0;jdof<nbrdof; jdof++)
	 stifffile << " " << gsl_matrix_get(jac,idof,jdof);
     stifffile << endl;
     }
  stifffile.close();
  cout << " Done\n";

  cout << "Constructing and saving damping matrix ...";
  Calcule_Matrice_J(0,1,0,1E-3,doflocked);
  strcpy(nom,application);
  strcat(nom,".cc");
  ofstream dampfile(nom);
  if (!dampfile) EasyDynsimError("Cannot create damping matrix file");
  for (idof=0; idof<nbrdof; idof++)
     {
     for (jdof=0;jdof<nbrdof; jdof++)
	 dampfile << " " << gsl_matrix_get(jac,idof,jdof);
     dampfile << endl;
     }
  dampfile.close();
  cout << " Done\n";

  cout << "Constructing and saving mass matrix ...";
  Calcule_Matrice_J(1,0,0,1,doflocked);
  strcpy(nom,application);
  strcat(nom,".mm");
  ofstream massfile(nom);
  if (!massfile) EasyDynsimError("Cannot create mass matrix file");
  for (idof=0; idof<nbrdof; idof++)
     {
     for (jdof=0;jdof<nbrdof; jdof++)
	 massfile << " " << gsl_matrix_get(jac,idof,jdof);
     massfile << endl;
     }
  massfile.close();
  cout << " Done\n";

  if (nbrinput>0)
    {
    cout << "Constructing and saving input contribution matrix ...";
    double **FF, tol=1e-3;
    int iinput;
    FF=new double*[nbrdof];
    for (idof=0; idof<nbrdof;idof++) FF[idof]=new double[nbrinput];
    for (iinput=0; iinput<nbrinput; iinput++)
      {
      u[iinput]+=tol;
      ComputeResidual();
      for (idof=0; idof<nbrdof;idof++) FF[idof][iinput]=0.5*f[idof]/tol;
      u[iinput]-=2*tol;
      ComputeResidual();
      for (idof=0; idof<nbrdof;idof++) FF[idof][iinput]-=0.5*f[idof]/tol;
      u[iinput]+=tol;
      }
    strcpy(nom,application);
    strcat(nom,".ff");
    ofstream fffile(nom);
    if (!fffile) EasyDynsimError("Cannot create influence matrix file");
    for (idof=0;idof<nbrdof; idof++)
      {
      for (iinput=0; iinput<nbrinput; iinput++)
        fffile << " " << -FF[idof][iinput];
      fffile << endl;
      }
    fffile.close();
    cout << " Done\n";
    for (idof=0;idof<nbrdof; idof++) delete FF[idof];
    delete FF;
    }
}

//---------------------------------------------------------------------------

void ComputePoles(int *doflocked, double freqmin, double freqmax)
// Procede au calcul des poles pour la configuration actuelle

{
  int NN,INFO,*Liste;
  gsl_matrix *AA,*BB;
  gsl_matrix_complex *VEC;
  gsl_vector_complex *ALPHA;
  gsl_vector *BETA;
  gsl_eigen_genv_workspace *WORK;
  NN=2*nbrdof;
  gsl_matrix_complex *PSI;
  gsl_matrix_complex *MM_PSI;
  gsl_matrix_complex *MM;
  gsl_matrix_complex *MODMASS;
  gsl_complex mass_element;
  gsl_complex unit1 = gsl_complex_rect(1.0,0.0); // needed in method zgemm of gsl to compute the modal mass
  gsl_complex zero1 = gsl_complex_rect(0.0,0.0); // needed in method zgemm of gsl to compute the modal mass
  // Memory allocation
  AA=gsl_matrix_alloc(NN,NN);
  BB=gsl_matrix_alloc(NN,NN);
  VEC=gsl_matrix_complex_alloc(NN,NN);
  ALPHA=gsl_vector_complex_alloc(NN);
  BETA=gsl_vector_alloc(NN);
  WORK=gsl_eigen_genv_alloc(NN);
  Liste=new int[NN];
  PSI=gsl_matrix_complex_alloc(nbrdof,1);
  MM_PSI=gsl_matrix_complex_alloc(nbrdof,1);
  MM=gsl_matrix_complex_alloc(nbrdof,nbrdof);
  MODMASS=gsl_matrix_complex_alloc(1,1);
  mass_element=GSL_COMPLEX_ZERO;

  // Building eigenvalue problem matrices
  int iddl,jddl,ipole,jpole;
  // Zeroing matrices
  gsl_matrix_set_zero(AA);
  gsl_matrix_set_zero(BB);
  // Adding stiffness matrix contribution
  cout << "Constructing Stiffness matrix ...";
  Calcule_Matrice_J(0,0,1,1E-5,doflocked);
  cout << " Done\n";
  for (iddl=0; iddl<nbrdof; iddl++) for (jddl=0;jddl<nbrdof; jddl++)
     gsl_matrix_set(AA,iddl,jddl,-gsl_matrix_get(jac,iddl,jddl));

  // Adding damping matrix contribution
  cout << "Constructing Damping matrix ...";
  Calcule_Matrice_J(0,1,0,1E-3,doflocked);
  cout << " Done\n";
  for (iddl=0; iddl<nbrdof; iddl++) for (jddl=0;jddl<nbrdof; jddl++)
     gsl_matrix_set(BB,iddl,jddl,gsl_matrix_get(jac,iddl,jddl));

  // Adding mass matrix contribution
  cout << "Constructing Mass matrix ...";
  Calcule_Matrice_J(1,0,0,1,doflocked);
  cout << " Done\n";
  for (iddl=0; iddl<nbrdof; iddl++) for (jddl=0;jddl<nbrdof; jddl++)
     {
     gsl_matrix_set(AA,nbrdof+iddl,nbrdof+jddl,gsl_matrix_get(jac,iddl,jddl));
     gsl_matrix_set(BB,nbrdof+iddl,jddl,gsl_matrix_get(jac,iddl,jddl));
     gsl_matrix_set(BB,iddl,nbrdof+jddl,gsl_matrix_get(jac,iddl,jddl));
     mass_element=gsl_complex_rect(gsl_matrix_get(jac,iddl,jddl),0.0); // set the mass matrix into its complex form (imaginary part =0)
     gsl_matrix_complex_set(MM,iddl,jddl,mass_element);
     }

  // Calling GSL routine for nonsymmetric problem
  INFO=gsl_eigen_genv(AA,BB,ALPHA,BETA,VEC,WORK);
  gsl_eigen_genv_free(WORK); // immediately freeing work memory
  if (INFO!=0) // did everything happen properly ?
   {
   cout << "Code returned by gsl_eigen_genv: " << INFO << endl;
   EasyDynsimError("gsl_eigen_genv failed");
   }

  // Reporting eigen values to ALPHA
  double beta; gsl_complex alpha;
  for (ipole=0; ipole<NN; ipole++)
    {
    beta=gsl_vector_get(BETA,ipole);
    alpha=gsl_vector_complex_get(ALPHA,ipole);
    if (fabs(beta)<1E-15)
      {
      beta=1.0;
      alpha=gsl_complex_rect(1E30,1E30);
      }
    alpha=gsl_complex_div_real(alpha,beta);
    gsl_vector_set(BETA,ipole,beta);
    gsl_vector_complex_set(ALPHA,ipole,alpha);
    }

  // Building ordered list of poles
  double sigmai,sigmaj, omegai, omegaj, modalmass, modalmass2;
  for (ipole=0; ipole<NN; ipole++) Liste[ipole]=ipole;
  for (ipole=0; ipole<NN-1; ipole++) for (jpole=ipole; jpole<NN; jpole++)
    {
    alpha=gsl_vector_complex_get(ALPHA,Liste[ipole]);
    omegai=GSL_IMAG(alpha);
    sigmai=GSL_REAL(alpha);
    alpha=gsl_vector_complex_get(ALPHA,Liste[jpole]);
    omegaj=GSL_IMAG(alpha);
    sigmaj=GSL_REAL(alpha);
    if (omegaj<omegai)
           {
           int itmp=Liste[ipole];
           Liste[ipole]=Liste[jpole];
           Liste[jpole]=itmp;
           }
    if (fabs(omegaj-omegai)<1E-15)
      if (sigmaj<sigmai)
           {
           int itmp=Liste[ipole];
           Liste[ipole]=Liste[jpole];
           Liste[jpole]=itmp;
           }
    }

  // printing eigenvalues if in debug mode
  if (DEBUG) for (ipole=0; ipole<NN; ipole++)
    {
    alpha=gsl_vector_complex_get(ALPHA,Liste[ipole]);
    omegai=GSL_IMAG(alpha);
    sigmai=GSL_REAL(alpha);
    if (omegai>-1E-15)
        cout << "Pole " << Liste[ipole] << "=" << sigmai
             << "+-i" << omegai << "\n";
    }
  // Building complete results file
  int numpole, FREQMINMAX=0;
  double polenorm,freq,amor,coeff;
  char *nom;
  if ((freqmin!=0) || (freqmax!=1E30)) FREQMINMAX=1;

  nom=new char[strlen(application)+5];
  strcpy(nom,application);
  strcat(nom,".lst"); // file with eigenvalues
  ofstream lstfile(nom);
  if (!lstfile)
    {
    cout<<"\n Fichier "<<nom<<" impossible a ouvrir";
    exit(1);
    }
  lstfile << "  POLE       ALPHA        OMEGA   "
        << "FREQ(Hz)  DAMP_RATIO  MASS_MOD\n";

  ofstream ls2file;
  if (FREQMINMAX)
    {
    strcpy(nom,application);
    strcat(nom,".ls2");
    ls2file.open(nom);
    if (!ls2file)
      {
      cout<<"\n Cannot open " << nom;
      exit(1);
      }
    }
  strcpy(nom,application);
  strcat(nom,".mod");
  ofstream modfile(nom);
  if (!modfile)
    {
      cout<<"\n Cannot open " << nom << endl;
    exit(1);
    }

  for(ipole=0;ipole<NN;ipole++)
    {
    numpole=Liste[ipole];
    alpha=gsl_vector_complex_get(ALPHA,numpole);
    omegai=GSL_IMAG(alpha);
    sigmai=GSL_REAL(alpha);
    polenorm=hypot(sigmai,omegai);
    freq=omegai/(2*3.14159265);
    if (polenorm<1E-15) polenorm=1;
    amor=-sigmai/polenorm;
    if (omegai>-1E-15)  // conjugate poles are saved only once
      {
	// saving modal properties in lst file (but modal mass)
      lstfile << setw(5) << ipole+1
              << " " << setw(12) << setprecision(5) << sigmai
              << " " << setw(12) << setprecision(5) << omegai
              << " " << setw(10) << setprecision(4) << freq
              << " " << setw(11) << setprecision(3) << amor;
      if ((freq>=freqmin) && (freq<=freqmax)) // only selected frequency
	// range is saved in mod and ls2 files
        {
        // Constructing mode
        for (iddl=0; iddl<nbrdof; iddl++)
        {
	      alpha=gsl_matrix_complex_get(VEC,nbrdof+iddl,numpole);
	      gsl_matrix_complex_set(PSI,iddl,0,alpha);
        }

        // Computing Modal Mass for each pole (Psi_i^T * M * Psi_i = m_i)
        gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, unit1, MM, PSI, zero1, MM_PSI);
        gsl_blas_zgemm(CblasTrans, CblasNoTrans, unit1, PSI, MM_PSI, zero1, MODMASS);
        modalmass=GSL_REAL(gsl_matrix_complex_get(MODMASS,0,0));
        // Saving modal mass in lst file
        lstfile << " " << setw(9) << setprecision(3) << fabs(modalmass) << "\n";
	// Saving modal properties in mod and ls2 files
        if (FREQMINMAX) ls2file << setw(5) << ipole+1
                               << " " << setw(12) << setprecision(5) << sigmai
                               << " " << setw(12) << setprecision(5) << omegai
                               << " " << setw(10) << setprecision(4) << freq
                               << " " << setw(11) << setprecision(3) << amor
                               <<  " " << setw(9) << setprecision(3) << fabs(modalmass)
			                   << "\n";
        modfile << setw(5) << ipole+1
                << " " << setw(12) << setprecision(5) << sigmai
                << " " << setw(12) << setprecision(5) << omegai
                << " " << setw(10) << setprecision(4) << freq
                << " " << setw(11) << setprecision(3) << amor
                << " " << setw(9) << setprecision(3) << fabs(modalmass) << "\n";
         // Saving mode in mod file
         for (iddl=0; iddl<nbrdof; iddl++)
	     {
             alpha=gsl_matrix_complex_get(VEC,nbrdof+iddl,numpole);
             modfile << GSL_REAL(alpha) << " " << GSL_IMAG(alpha) << endl;
         }
        }
        else
        {            // Constructing mode
            for (iddl=0; iddl<nbrdof; iddl++)
	        {
	           alpha=gsl_matrix_complex_get(VEC,nbrdof+iddl,numpole);
	           gsl_matrix_complex_set(PSI,iddl,0,alpha);
            }
            // Computing Modal Mass for each pole (Psi_i^T * M * Psi_i = m_i)
            gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, unit1, MM, PSI, zero1, MM_PSI);
            gsl_blas_zgemm(CblasTrans, CblasNoTrans, unit1, PSI, MM_PSI, zero1, MODMASS);
            modalmass=GSL_REAL(gsl_matrix_complex_get(MODMASS,0,0));
            // Saving modal mass in lst file
            lstfile << " " << setw(9) << setprecision(3) << fabs(modalmass) << "\n";
        }
      }
    }
  lstfile.close();
  if (FREQMINMAX) ls2file.close();
  modfile.close();

  // Cleaning memeory
  gsl_matrix_free(AA);
  gsl_matrix_free(BB);
  gsl_matrix_complex_free(VEC);
  gsl_vector_complex_free(ALPHA);
  gsl_vector_free(BETA);
  gsl_matrix_complex_free(PSI);
  gsl_matrix_complex_free(MM_PSI);
  gsl_matrix_complex_free(MM);
  gsl_matrix_complex_free(MODMASS);

  delete Liste;
}

//---------------------------------------------------------------------------
