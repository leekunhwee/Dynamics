/*

Copyright (C) 2003 Olivier VERLINDEN
    Service de Mecanique rationnelle, Dynamique et Vibrations
    Faculte Polytechnique de Mons
    31, Bd Dolez, 7000 MONS (Belgium)
    Olivier.Verlinden@fpms.ac.be

This file is part of EasyAnim

EasyAnim is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the
Free Software Foundation; either version 2, or (at your option) any
later version.

EasyAnim is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with EasyAnim; see the file COPYING.  If not, write to the Free
Software Foundation, 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

*/

#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include <iostream>
#include "EasyAnimcnv.h"
#include <GL/glu.h>
#include <GL/gl.h>

#define _FATAL   1
#define _WARNING 0

def_repere* s1=NULL;
def_volvan* s2=NULL;
def_volvmo* s3=NULL;
def_uff* s4=NULL;

//*********************************************************************

istream &toeol(istream &stream)
  {
  char c;
  do { stream.get(c); } while ((c!='\r') && (c!='\n') && (!stream.rdstate()));
  return stream;
  }

/********************** Gestion des messages d'erreur *****************/

void ErreurMessage(int i,const char *message)
{
 if (i) std::cout<<"\n !!! ERREUR FATALE !!!\n";
	else std::cout<<"\n !!! ATTENTION !!!\n";
	std::cout<<message<<"\n";
 if (i) exit(1);
}

/*---------------------------------------------------------------------*/

/* impression des messages d'erreur */
/* si ErrType < 100 : arret du programme */
void ErreurMessage(int ErrType,const char *mot,int num_lig)
{

 if (ErrType < 100) printf("\n !!! ERREUR FATALE !!!\n");
switch (ErrType)
 {
  case 1 : printf("\t fichier %s impossible a ouvrir. ",mot);
	   break;
  case 2 : printf("\t nombre de points %d incorrect dans le fichier %d.",num_lig,mot);
	   break;
  case 3 : printf("\t numero de point %d incorrect dans le fichier %s.",num_lig,mot);
	   break;
  case 4 : printf("\t nombre d'arretes %d incorrect dans le fichier %d.",num_lig,mot);
	   break;
  case 5 : printf("\t numero d'arrete %d incorrect dans le fichier %s.",num_lig,mot);
	   break;
  case 6 : printf("\t nombre de faces %d incorrect dans le fichier %d.",num_lig,mot);
	   break;
  case 7 : printf("\t numero de faces %d incorrect dans le fichier %s.",num_lig,mot);
	   break;
  case 8 : printf("\t nombre de sommets dans la face %d incorrect dans le fichier %s.",num_lig,mot);
	   break;
  case 9 : printf("\t sommet non defini dans l'arrete %d du fichier %s.",num_lig,mot);
	   break;

  case 10 : printf("\t sommet non defini dans la face %d du fichier %s.",num_lig,mot);
	   break;

  case 101 : printf("\n ****** Lecture de %s avec succes ****** \n",mot);
	   break;
  case 102 : printf(" ****** Nombre de lignes : %4d      ****** \n",num_lig);
	     break;
  case 103 : printf("\n ****** Lecture de %s             ****** \n",mot);
	     break;
  default : printf("\t il y a une erreur a la ligne %d.",num_lig);
	    break;
 }
 if (ErrType < 100) exit(0);
}

/*---------------------------------------------------------------------*/
/*                    Methods relative to def_volvan                   */
/*---------------------------------------------------------------------*/

/*---------------------------------------------------------------------*/

def_volvan::def_volvan(const char *VolFileName, int &errcode)
{
 int ll=strlen(VolFileName);
 strcpy(FileName,VolFileName);
 FileName[ll-4]=0;
 errcode=0;
 /* lecture du fichier *.VOL */
 char titre[180];
 int i,j,k,l,num;
 float x,y,z;
 ifstream VolFile(VolFileName);
 if (!VolFile) { errcode=6; return; }

 //Reading the title
 VolFile.getline(titre,179);

 // Reading the nodes
 VolFile >> NbrNodes >> toeol;
 if (NbrNodes<=0) { errcode=4; return; }
 node=new NodeStruct[NbrNodes];
 if (node==NULL) { errcode=5; return; }
 xmin=1E10; xmax=-1E10; ymin=1E10; ymax=-1E10; zmin=1E10; zmax=-1E10;
 for (i=0;i<NbrNodes;i++)
   {
   VolFile >> node[i].name >> node[i].x >> node[i].y >> node[i].z >> toeol;
   if(node[i].x<xmin) xmin=node[i].x;
   if(node[i].x>xmax) xmax=node[i].x;
   if(node[i].y<ymin) ymin=node[i].y;
   if(node[i].y>ymax) ymax=node[i].y;
   if(node[i].z<zmin) zmin=node[i].z;
   if(node[i].z>zmax) zmax=node[i].z;
   }/*endfor i*/
  // Memory allocation for the spheres representing the nodes
  NodeSphere=gluNewQuadric();
  gluQuadricNormals(NodeSphere,GLU_FLAT);
  gluQuadricTexture(NodeSphere,GL_FALSE);
  gluQuadricOrientation(NodeSphere,GLU_OUTSIDE);
  gluQuadricDrawStyle(NodeSphere,GLU_FILL);

 // Reading the edges
 VolFile >> NbrEdges >> toeol;
 if (NbrEdges<0) { errcode=4; return; }
 if (NbrEdges) 
    {
    edge=new EdgeStruct[NbrEdges];
    if (edge==NULL) { errcode=5; return; }
    }
 for (i=0;i<NbrEdges;i++)
   {
   char Node1Name[40],Node2Name[40];
   VolFile >> edge[i].name >> Node1Name >> Node2Name >> edge[i].color;
   edge[i].node1=GetNumNodeOfName(Node1Name);
   edge[i].node2=GetNumNodeOfName(Node2Name);
   if (edge[i].node1==-1) { errcode=4; return; }
   if (edge[i].node2==-1) { errcode=4; return; }
   }/*endfor i*/

 // Reading the sides
 VolFile >> NbrSides >> toeol;
 if (NbrSides<0) { errcode=4; return; }
 if (NbrSides) 
    {
    side=new SideStruct[NbrSides];
    if (side==NULL) { errcode=5; return; }
    }
 for (i=0;i<NbrSides;i++)
    {
    VolFile >> side[i].name >> side[i].NbrNodes;
    if (side[i].NbrNodes<=0) { errcode=4; return; }
    side[i].NumNode=new int[side[i].NbrNodes];
    if (side[i].NumNode==NULL) { errcode=5; return; }
    for (k=0;k<side[i].NbrNodes;k++)
       {
       char NodeName[40];
       VolFile >> NodeName;
       side[i].NumNode[k]=GetNumNodeOfName(NodeName);
       if (side[i].NumNode[k]==-1) { errcode=4; return; }
       }
    VolFile >> side[i].color >> toeol;
    } /*endfor i*/
 VolFile.close();

 // Lecture du fichier van s'il existe
 char nomf[200]="";
 strcpy(nomf,FileName);
 strcat(nomf,".van");
 fpos.open(nomf);
 if (!fpos) errcode=2;
 NbrImages=0;
 CurImage=0;
 if (fpos)
   {
   while (!fpos.rdstate())
     {
     for (i=0; i<NbrNodes;i++) fpos >> x >> y >> z;
     if (!fpos.rdstate()) NbrImages++;
     }
   fpos.close();
   fpos.clear();
   fpos.open(nomf);
   }
 
 // Lecture du fichier cfg s'il existe
 strcpy(nomf,FileName);
 strcat(nomf,".cfg");
 if (Lecture_Config(nomf)) { errcode+=1; return; }

 }

/*---------------------------------------------------------------------*/

// destructeur
def_volvan::~def_volvan()
{
if (fpos) fpos.close();
}	

/*---------------------------------------------------------------------*/

void def_volvan::Next()
// Lecture des coordonnées
  {
  if (NbrImages) 
    {
    int i;
    if (CurImage<NbrImages) 
      {
      for (i=0;i<NbrNodes;i++) fpos>>node[i].x>>node[i].y>>node[i].z;
      CurImage++;
      }
    else
      {
      fpos.close();
      fpos.clear();
      char VanFileName[200]="";
      strcpy(VanFileName,FileName);
      strcat(VanFileName,".van");
      fpos.open(VanFileName);
      if (fpos)
	{
        for (i=0;i<NbrNodes;i++) fpos>>node[i].x>>node[i].y>>node[i].z;
        CurImage=1;
        }
      else std::cout << "Cannot reopen file " << VanFileName << "\n";
      }
    }
  }

/*---------------------------------------------------------------------*/

void def_volvan::LoadImage(int NumImage)
// Lecture des coordonnées
  {
  if (NbrImages) 
    {
    if (NumImage>NbrImages) NumImage=NbrImages;
    if (NumImage<1) NumImage=1;
    int i,j;
    if (NumImage>CurImage) while (CurImage<NumImage)
      {
      for (i=0;i<NbrNodes;i++) fpos>>node[i].x>>node[i].y>>node[i].z;
      CurImage++;
      }
    else
      {
      fpos.close();
      fpos.clear();
      char VanFileName[200]="";
      strcpy(VanFileName,FileName);
      strcat(VanFileName,".van");
      fpos.open(VanFileName);
      CurImage=0;
      if (fpos) while (CurImage<NumImage)
	{
        for (i=0;i<NbrNodes;i++) fpos>>node[i].x>>node[i].y>>node[i].z;
        CurImage++;
        }
      else std::cout << "Cannot reopen file " << VanFileName << "\n";
      }
    }
  }

/*---------------------------------------------------------------------*/

/*---------------------------------------------------------------------*/
/*                    Methods relative to def_volvmo                   */
/*---------------------------------------------------------------------*/

/*---------------------------------------------------------------------*/

def_volvmo::def_volvmo(const char *VolFileName, int &errcode)
{
 int ll=strlen(VolFileName);
 strcpy(FileName,VolFileName);
 FileName[ll-4]=0;
 errcode=0;
 /* lecture du fichier *.VOL */
 char titre[180];
 int i,j,k,l,num;
 ifstream VolFile(VolFileName);
 if (!VolFile) { errcode=6; return; }

 //Reading the title
 VolFile.getline(titre,179);

 // Reading the nodes
 VolFile >> NbrNodes >> toeol;
 if (NbrNodes<=0) { errcode=4; return; }
 node=new NodeStruct[NbrNodes];
 if (node==NULL) { errcode=5; return; }
 PosRef=new VecStruct[NbrNodes];
 if (PosRef==NULL) { errcode=5; return; }
 DepReal=new VecStruct[NbrNodes];
 if (DepReal==NULL) { errcode=5; return; }
 DepImag=new VecStruct[NbrNodes];
 if (DepImag==NULL) { errcode=5; return; }
 xmin=1E10; xmax=-1E10; ymin=1E10; ymax=-1E10; zmin=1E10; zmax=-1E10;
 for (i=0;i<NbrNodes;i++)
   {
   VolFile >> node[i].name >> node[i].x >> node[i].y >> node[i].z >> toeol;
   if(node[i].x<xmin) xmin=node[i].x;
   if(node[i].x>xmax) xmax=node[i].x;
   if(node[i].y<ymin) ymin=node[i].y;
   if(node[i].y>ymax) ymax=node[i].y;
   if(node[i].z<zmin) zmin=node[i].z;
   if(node[i].z>zmax) zmax=node[i].z;
   PosRef[i].x=node[i].x;
   PosRef[i].y=node[i].y;
   PosRef[i].z=node[i].z;
   DepReal[i].x=0.0; DepImag[i].x=0.0;
   DepReal[i].y=0.0; DepImag[i].y=0.0;
   DepReal[i].z=0.0; DepImag[i].z=0.0;
   }/*endfor i*/

  // Memory allocation for the spheres representing the nodes
  NodeSphere=gluNewQuadric();
  gluQuadricNormals(NodeSphere,GLU_FLAT);
  gluQuadricTexture(NodeSphere,GL_FALSE);
  gluQuadricOrientation(NodeSphere,GLU_OUTSIDE);
  gluQuadricDrawStyle(NodeSphere,GLU_FILL);

 // Reading the edges
 VolFile >> NbrEdges >> toeol;
 if (NbrEdges<0) { errcode=4; return; }
 if (NbrEdges) 
    {
    edge=new EdgeStruct[NbrEdges];
    if (edge==NULL) { errcode=5; return; }
    }
 for (i=0;i<NbrEdges;i++)
   {
   char Node1Name[40],Node2Name[40];
   VolFile >> edge[i].name >> Node1Name >> Node2Name >> edge[i].color;
   edge[i].node1=GetNumNodeOfName(Node1Name);
   edge[i].node2=GetNumNodeOfName(Node2Name);
   if (edge[i].node1==-1) { errcode=4; return; }
   if (edge[i].node2==-1) { errcode=4; return; }
   }/*endfor i*/

 // Reading the sides
 VolFile >> NbrSides >> toeol;
 if (NbrSides<0) { errcode=4; return; }
 if (NbrSides) 
    {
    side=new SideStruct[NbrSides];
    if (side==NULL) { errcode=5; return; }
    }
 for (i=0;i<NbrSides;i++)
    {
    VolFile >> side[i].name >> side[i].NbrNodes;
    if (side[i].NbrNodes<=0) { errcode=4; return; }
    side[i].NumNode=new int[side[i].NbrNodes];
    if (side[i].NumNode==NULL) { errcode=5; return; }
    for (k=0;k<side[i].NbrNodes;k++)
       {
       char NodeName[40];
       VolFile >> NodeName;
       side[i].NumNode[k]=GetNumNodeOfName(NodeName);
       if (side[i].NumNode[k]==-1) { errcode=4; return; }
       }
    VolFile >> side[i].color >> toeol;
    } /*endfor i*/
 VolFile.close();

 // Reading vmo file if existing
 if (LoadMode(1)) errcode=2;

 // reading cfg if existing
 char nomf[200];
 strcpy(nomf,FileName);
 strcat(nomf,".cfg");
 if (Lecture_Config(nomf)) { errcode+=1; return; }

 }

// destructeur
def_volvmo::~def_volvmo()
{
}	

/*---------------------------------------------------------------------*/

// destructeur
int def_volvmo::LoadMode(int num)
{
 char nomf[200]="";
 int imode, ipt,nbrmode;
 strcpy(nomf,FileName);
 strcat(nomf,".vmo");
 ifstream vmofile(nomf);
 // vmofile >> nbrmode;
 if (num<1) return 1;
 for (imode=0; imode<num; imode++)
   {
   // Lecture des infos sur le pole correspondant
   vmofile >> alpha >> beta >> freq >> amor;
   // Initialisation de DepMax pour eviter des divisions par zero
   DepMax=1E-9;
   // Lecture de la partie relle du mode 
   for (ipt=0; ipt<NbrNodes; ipt++)
     {
     vmofile >> DepReal[ipt].x >> DepReal[ipt].y >> DepReal[ipt].z;
     if (fabs(DepReal[ipt].x-PosRef[ipt].x)>DepMax) 
       DepMax=fabs(DepReal[ipt].x-PosRef[ipt].x);
     if (fabs(DepReal[ipt].y-PosRef[ipt].y)>DepMax) 
       DepMax=fabs(DepReal[ipt].y-PosRef[ipt].y);
     if (fabs(DepReal[ipt].z-PosRef[ipt].z)>DepMax) 
       DepMax=fabs(DepReal[ipt].z-PosRef[ipt].z);
     }
   // Lecture de la partie imaginaire du mode 
   for (ipt=0; ipt<NbrNodes; ipt++)
     {
     vmofile >> DepImag[ipt].x >> DepImag[ipt].y >> DepImag[ipt].z;
     if (fabs(DepImag[ipt].x-PosRef[ipt].x)>DepMax) 
       DepMax=fabs(DepImag[ipt].x-PosRef[ipt].x);
     if (fabs(DepImag[ipt].y-PosRef[ipt].y)>DepMax) 
       DepMax=fabs(DepImag[ipt].y-PosRef[ipt].y);
     if (fabs(DepImag[ipt].z-PosRef[ipt].z)>DepMax) 
       DepMax=fabs(DepImag[ipt].z-PosRef[ipt].z);
     }
   }
 if (vmofile.rdstate()) 
   {
   // Puisque ca n'a pas marche, on remet tout a zero
     alpha=0.0; beta=0.0; freq=0.0; amor=0.0;
   for (ipt=0; ipt<NbrNodes; ipt++)
     { DepReal[ipt].x=PosRef[ipt].x; DepImag[ipt].x=PosRef[ipt].x; 
       DepReal[ipt].y=PosRef[ipt].y; DepImag[ipt].y=PosRef[ipt].y; 
       DepReal[ipt].z=PosRef[ipt].z; DepImag[ipt].z=PosRef[ipt].z; 
       node[ipt].x=PosRef[ipt].x; node[ipt].y=PosRef[ipt].y;
       node[ipt].z=PosRef[ipt].z; }
   NbrImages=0;
   return 2;
   }
 vmofile.close();
 float lcaract;
 lcaract=sqrt((xmax-xmin)*(xmax-xmin)+(ymax-ymin)*(ymax-ymin)
               +(zmax-zmin)*(zmax-zmin));
 for (ipt=0;ipt<NbrNodes;ipt++)
      {
      node[ipt].x=PosRef[ipt].x
                +(DepReal[ipt].x-PosRef[ipt].x)*0.05*lcaract/DepMax;
      node[ipt].y=PosRef[ipt].y
                +(DepReal[ipt].y-PosRef[ipt].y)*0.05*lcaract/DepMax;
      node[ipt].z=PosRef[ipt].z
                +(DepReal[ipt].z-PosRef[ipt].z)*0.05*lcaract/DepMax;
      }
 CurImage=0;
 NbrImages=50;
 return 0;
}	

/*---------------------------------------------------------------------*/

void def_volvmo::Next()
// Construction de l'image modale
  {
  if (NbrImages)
    {
    int ipt;
    float lcaract, pi=3.14151695;
    lcaract=sqrt((xmax-xmin)*(xmax-xmin)+(ymax-ymin)*(ymax-ymin)
                +(zmax-zmin)*(zmax-zmin));
    CurImage=(CurImage++)%NbrImages;
    for (ipt=0;ipt<NbrNodes;ipt++)
      {
      node[ipt].x=PosRef[ipt].x
                 +(DepReal[ipt].x-PosRef[ipt].x)*0.05*lcaract
                       *cos(2.0*CurImage*pi/NbrImages)/DepMax
                 +(DepImag[ipt].x-PosRef[ipt].x)*0.05*lcaract
                       *sin(2.0*CurImage*pi/NbrImages)/DepMax;
      node[ipt].y=PosRef[ipt].y
                 +(DepReal[ipt].y-PosRef[ipt].y)*0.05*lcaract
                       *cos(2.0*CurImage*pi/NbrImages)/DepMax
                 +(DepImag[ipt].y-PosRef[ipt].y)*0.05*lcaract
                       *sin(2.0*CurImage*pi/NbrImages)/DepMax;
      node[ipt].z=PosRef[ipt].z
                 +(DepReal[ipt].z-PosRef[ipt].z)*0.05*lcaract
                       *cos(2.0*CurImage*pi/NbrImages)/DepMax
                 +(DepImag[ipt].z-PosRef[ipt].z)*0.05*lcaract
                       *sin(2.0*CurImage*pi/NbrImages)/DepMax;
      }
    }
  }

/*---------------------------------------------------------------------*/

void def_volvmo::LoadImage(int NumImage)
// Building a requested image
  {
  if (NbrImages) 
    {
    int ipt;
    float lcaract, pi=3.14151695;
    lcaract=sqrt((xmax-xmin)*(xmax-xmin)+(ymax-ymin)*(ymax-ymin)
                 +(zmax-zmin)*(zmax-zmin));
    CurImage=NumImage%NbrImages;
    for (ipt=0;ipt<NbrNodes;ipt++)
      {
      node[ipt].x=PosRef[ipt].x
                 +(DepReal[ipt].x-PosRef[ipt].x)*0.05*lcaract
                       *cos(2.0*CurImage*pi/NbrImages)/DepMax
                 +(DepImag[ipt].x-PosRef[ipt].x)*0.05*lcaract
                       *sin(2.0*CurImage*pi/NbrImages)/DepMax;
      node[ipt].y=PosRef[ipt].y
                 +(DepReal[ipt].y-PosRef[ipt].y)*0.05*lcaract
                       *cos(2.0*CurImage*pi/NbrImages)/DepMax
                 +(DepImag[ipt].y-PosRef[ipt].y)*0.05*lcaract
                       *sin(2.0*CurImage*pi/NbrImages)/DepMax;
      node[ipt].z=PosRef[ipt].z
                 +(DepReal[ipt].z-PosRef[ipt].z)*0.05*lcaract
                       *cos(2.0*CurImage*pi/NbrImages)/DepMax
                 +(DepImag[ipt].z-PosRef[ipt].z)*0.05*lcaract
                       *sin(2.0*CurImage*pi/NbrImages)/DepMax;
      }
    }
  }

/*---------------------------------------------------------------------*/

/*---------------------------------------------------------------------*/
/*                      Methods relative to def_uff                    */
/*---------------------------------------------------------------------*/

/*---------------------------------------------------------------------*/

int GetWords2(const char *Ligne, int &nbrmotslus, char **Mot)

/* Fonction permettant d'extraire d'une ligne tous les mots separes par un
   ou plusieurs espaces et tabulations */

  {

  char c;
  int ich,lncrtword=0; // longueur du mot courant;

  nbrmotslus=0; // nombre de mots
  for (ich=0; ich<strlen(Ligne); ich++) // boucle principale
    {
      c=Ligne[ich];
      if (isspace(c) && (lncrtword>0))  // cloture du mot si premier
        {                               // espace apres un mot
        Mot[nbrmotslus-1][lncrtword]='\0';
        lncrtword=0;
        }
      // construction du mot si caractere adequat
      if (!isspace(c))
        {
        if (lncrtword==0) nbrmotslus++; // nouveau mot si prem. car. apres espaces
	//if (nbrmotslus>DIMMOT) RunPrepError(nli,"Too many words in the command");
        lncrtword++;
        //if (lncrtword>=LONGMOT) RunPrepError(nli,"Word too long");
        Mot[nbrmotslus-1][lncrtword-1]=c;
        }
     }

} // fin de GetWords

/***************************************************************************/

int GetWordsInLine(istream &stream, int &nbrmotslus, char **Mot)

/* Fonction permettant d'extraire d'une ligne tous les mots separes par un
   ou plusieurs espaces et tabulations */

  {
  char c;
  int lncrtword=0; // longueur du mot courant;
  nbrmotslus=0; // nombre de mots
  if (stream.rdstate()) return 1;
  do
    {
    // Reading spaces
    while ((stream.peek()==' ') || (stream.peek()=='\t')) c=stream.get();
    // Reading word
    nbrmotslus++; lncrtword=0;
    while ((stream.peek()!=' ') && (stream.peek()!='\t') &&
           (stream.peek()!='\n') && (stream.peek()!='\r') &&
           (stream.peek()!='\f') && (stream.peek()!=EOF))
      {
      Mot[nbrmotslus-1][lncrtword]=stream.get();
      lncrtword++;
      }
    if (lncrtword==0) nbrmotslus--;
    else Mot[nbrmotslus-1][lncrtword]='\0';
    }
  while ((stream.peek()!='\n') && (stream.peek()!='\f') 
         &&(stream.peek()!='\r') && !stream.rdstate()); 
  // Reading end of line characters
  while ((stream.peek()=='\n') || (stream.peek()=='\f') 
         || (stream.peek()=='\r')) c=stream.get(); 
  // Looking for end of file before returning
  if (stream.peek()==EOF) c=stream.get(); 
  return 0;
  } // fin de GetWords

/***************************************************************************/

int FindNextMinusOne(istream &stream, char **Mot)
{
int NbrWords;
do
  {
  if (GetWordsInLine(stream,NbrWords,Mot)) return 1;
  }
  while ((NbrWords!=1) || (strcmp(Mot[0],"-1")!=0));
return 0;
}

/***************************************************************************/

int FindNextUff(istream &stream, int type, char **Mot)
{
int NbrWords;
 do
   {
   if (FindNextMinusOne(stream,Mot)) return 1;
   if (GetWordsInLine(stream,NbrWords,Mot)) return 1;
   if ((NbrWords==1) && (atoi(Mot[0])==type))
     return 0;
   else if (FindNextMinusOne(stream,Mot)) return 1;
   }
   while (!stream.rdstate());
 return 1;
}

/***************************************************************************/

def_uff::def_uff(const char *UffFileName1,const char *UffFileName2, 
                       int &errcode)
{
 char Ligne[181],**Mot;
 Mot=new char*[20];
 int im;
 for (im=0;im<20;im++) Mot[im]=new char[30];
 int NbrWords;
 // First reading for counting the nodes
 ifstream UffFile1(UffFileName1);
 if (!UffFile1) { errcode=6; return; }
 if (FindNextUff(UffFile1,15,Mot)) { errcode=8; return; }
 NbrNodes=0;
 do
  {
  if (GetWordsInLine(UffFile1,NbrWords,Mot)) { errcode=7; NbrNodes=0; return; }
  NbrNodes++;
  }
 while ((NbrWords!=1) || (strcmp(Mot[0],"-1")!=0));
 UffFile1.close();
 NbrNodes--;
 cout << "Number of nodes:" << NbrNodes <<  endl;
 if (NbrNodes==0) { errcode=7; return; }  

 // Second reading for counting the edges
 UffFile1.clear();
 UffFile1.open(UffFileName1);
 if (!UffFile1) { errcode=6; return; }
 if (FindNextUff(UffFile1,82,Mot)) { errcode=4; return; }
 // Read number of items
 GetWordsInLine(UffFile1,NbrWords,Mot);
 int Nentry=atoi(Mot[1]),in,inode1,inode2;
 // Drop one comment line
 if (GetWordsInLine(UffFile1,NbrWords,Mot)) { errcode=3; NbrEdges=0; return; }
 NbrEdges=0;
 UffFile1 >> inode2;
 for (in=1; in<Nentry; in++)
       {
       inode1=inode2;
       UffFile1 >> inode2;
       if ((inode1!=0) && (inode2!=0)) NbrEdges++;
       }
 cout << "NbrEdges=" << NbrEdges << endl;
 if (FindNextMinusOne(UffFile1,Mot)) { errcode=3; NbrEdges=0; return; } 
 UffFile1.close();

 // Allocating memory
 node=new NodeStruct[NbrNodes];
 if (node==NULL) { errcode=5; return; }
 PosRef=new VecStruct[NbrNodes];
 if (PosRef==NULL) { errcode=5; return; }
 DepReal=new VecStruct[NbrNodes];
 if (DepReal==NULL) { errcode=5; return; }
 DepImag=new VecStruct[NbrNodes];
 if (DepImag==NULL) { errcode=5; return; }

 if (NbrEdges) edge=new EdgeStruct[NbrEdges];

 // Memory allocation for the spheres representing the nodes
 NodeSphere=gluNewQuadric();
 gluQuadricNormals(NodeSphere,GLU_FLAT);
 gluQuadricTexture(NodeSphere,GL_FALSE);
 gluQuadricOrientation(NodeSphere,GLU_OUTSIDE);
 gluQuadricDrawStyle(NodeSphere,GLU_FILL);

 // Second reading to assign the nodes
 cout << "Third reading for the nodes" << endl;
 UffFile1.clear();
 UffFile1.open(UffFileName1);
 if (!UffFile1) { errcode=6; return; }
 cout << "Third reading" << endl;
 if (FindNextUff(UffFile1,15,Mot)) 
     { errcode=8; NbrNodes=0; NbrEdges=0; return; }
 xmin=1E10; xmax=-1E10; ymin=1E10; ymax=-1E10; zmin=1E10; zmax=-1E10;
 for (in=0; in<NbrNodes; in++)
       {
       if (GetWordsInLine(UffFile1,NbrWords,Mot))
           { errcode=7; NbrNodes=0; NbrEdges=0; return; }
       strcpy(node[in].name,Mot[0]);
       node[in].x=atof(Mot[4]);
       node[in].y=atof(Mot[5]);
       node[in].z=atof(Mot[6]);
       if(node[in].x<xmin) xmin=node[in].x;
       if(node[in].x>xmax) xmax=node[in].x;
       if(node[in].y<ymin) ymin=node[in].y;
       if(node[in].y>ymax) ymax=node[in].y;
       if(node[in].z<zmin) zmin=node[in].z;
       if(node[in].z>zmax) zmax=node[in].z;
       PosRef[in].x=node[in].x;
       PosRef[in].y=node[in].y;
       PosRef[in].z=node[in].z;
       DepReal[in].x=0.0; DepImag[in].x=0.0;
       DepReal[in].y=0.0; DepImag[in].y=0.0;
       DepReal[in].z=0.0; DepImag[in].z=0.0;
       }
 if (FindNextMinusOne(UffFile1,Mot)) 
    { errcode=7; NbrNodes=0; NbrEdges=0; return; } 
 UffFile1.close();

 // Third reading to assign the edges
 cout << "Fourth reading for the edges" << endl;
 UffFile1.clear();
 UffFile1.open(UffFileName1);
 if (!UffFile1) { errcode=6; return; }
 cout << "Fourth reading" << endl;
 if (FindNextUff(UffFile1,82,Mot))
    { errcode=4; NbrNodes=0; NbrEdges=0; return; } 
 // Read number of items
 if (GetWordsInLine(UffFile1,NbrWords,Mot))
    { errcode=3; NbrNodes=0; NbrEdges=0; return; } 
 Nentry=atoi(Mot[1]);
 // Drop one comment line
 if (GetWordsInLine(UffFile1,NbrWords,Mot))
    { errcode=3; NbrNodes=0; NbrEdges=0; return; } 
 NbrEdges=0;
 char nomnode1[20],nomnode2[20];
 UffFile1 >> nomnode2;
 for (in=1; in<Nentry; in++)
       {
       strcpy(nomnode1,nomnode2);
       UffFile1 >> nomnode2;
       if ((strcmp(nomnode1,"0")!=0) && (strcmp(nomnode2,"0")!=0)) 
           { 
           NbrEdges++; 
           edge[NbrEdges-1].node1=GetNumNodeOfName(nomnode1);
           edge[NbrEdges-1].node2=GetNumNodeOfName(nomnode2);
           edge[NbrEdges-1].color=1;
           if (edge[NbrEdges-1].node1==-1) { errcode=4; return; }
           if (edge[NbrEdges-1].node2==-1) { errcode=4; return; }
           cout << nomnode1 << " " << nomnode2 << endl;
           cout << edge[NbrEdges-1].node1 << " " 
                << edge[NbrEdges-1].node2 << endl;
	   }
       }
 cout << "NbrEdges=" << NbrEdges << endl;
 if (FindNextMinusOne(UffFile1,Mot))
    { errcode=3; NbrNodes=0; NbrEdges=0; return; } 
 UffFile1.close();

 cout << "Number of nodes:" << NbrNodes <<  endl;
 cout << "Number of edges:" << NbrEdges <<  endl;
 strcpy(FileName,UffFileName2);

 // Reading modes if existing
 errcode=LoadMode(1);
 
}

/*---------------------------------------------------------------------*/

// destructeur
def_uff::~def_uff()
{
}	

/*---------------------------------------------------------------------*/

// destructeur
int def_uff::LoadMode(int num)
{
 int imode, ipt,nbrmode;
 // Resetting everything
 for (ipt=0; ipt<NbrNodes; ipt++)
     { DepReal[ipt].x=PosRef[ipt].x; DepImag[ipt].x=PosRef[ipt].x; 
       DepReal[ipt].y=PosRef[ipt].y; DepImag[ipt].y=PosRef[ipt].y; 
       DepReal[ipt].z=PosRef[ipt].z; DepImag[ipt].z=PosRef[ipt].z; 
       node[ipt].x=PosRef[ipt].x; node[ipt].y=PosRef[ipt].y;
       node[ipt].z=PosRef[ipt].z; }
 NbrImages=0;

 // Controlling mode number
 if (num<1) return 1;

 ifstream UffFile(FileName);
 if (!UffFile) { return 1; }
 char **Mot;
 Mot=new char*[20];
 int im;
 for (im=0;im<20;im++) Mot[im]=new char[30];

 int ComplexData, WithRot,NbrWords;
 // Skipping the first blocks of type 55
 for (imode=1; imode<num; imode++)
   {
   if (FindNextUff(UffFile,55,Mot)) { NbrImages=0; return 1; }
   if (FindNextMinusOne(UffFile,Mot)) { NbrImages=0; return 1; }
   }
 // up to the right one 
 if (FindNextUff(UffFile,55,Mot)) { NbrImages=0; return 1; }
 // Skipping the 5 comment lines
 for (im=0; im<5; im++)  if (GetWordsInLine(UffFile,NbrWords,Mot))
   { NbrImages=0; return 2; } 
 // Reading data parameters
 if (GetWordsInLine(UffFile,NbrWords,Mot))
   { NbrImages=0; return 2; }  
 if (NbrWords!=6)
   { NbrImages=0; cout << "Record 6: incorrect number of values\n"; return 2; }  
 if (strcmp(Mot[0],"1")!=0)
   { NbrImages=0; cout << "Model not structural\n"; return 2; }  
 if (strcmp(Mot[1],"2")!=0)
   { NbrImages=0; cout << "Analysis not modal\n"; return 2; }  
 if ((strcmp(Mot[2],"2")!=0) && (strcmp(Mot[2],"3")!=0))
   { NbrImages=0; cout << "Bad data type\n"; return 2; }  
 if (!strcmp(Mot[2],"2")) WithRot=0;
 if (!strcmp(Mot[2],"3")) WithRot=1;
 if ((strcmp(Mot[4],"2")!=0) && (strcmp(Mot[4],"5")!=0))
   { NbrImages=0; cout << "Bad data type\n"; return 2; }  
 if (!strcmp(Mot[4],"2")) ComplexData=0;
 if (!strcmp(Mot[4],"5")) ComplexData=1;
 if (WithRot &&  (strcmp(Mot[5],"6")!=0))
   { NbrImages=0; cout << "NDV incorrect\n"; return 2; }  
 if ((WithRot==0) &&  (strcmp(Mot[5],"3")!=0))
   { NbrImages=0; cout << "NDV incorrect\n"; return 2; }  
 // Skipping record 7
 if (GetWordsInLine(UffFile,NbrWords,Mot))
   { NbrImages=0; return 2; }  
 // Reading record 8
 if (GetWordsInLine(UffFile,NbrWords,Mot))
   { NbrImages=0; return 2; }  
 if (NbrWords!=4)
   { NbrImages=0; cout << "Record 8: incorrect number of values\n"; return 2; }
 freq=atof(Mot[0]);
 amor=atof(Mot[2]);
 cout << "freq=" << freq << " amor=" << amor << endl;
 // Initialisation of DepMax to avoid divisions by zero
 DepMax=1E-9;
 do
  {
  if (GetWordsInLine(UffFile,NbrWords,Mot)) { NbrImages=-1; }
  if (strcmp(Mot[0],"-1")!=0)
   {
   ipt=GetNumNodeOfName(Mot[0]);
   cout << "Node " << Mot[0] << endl;
   if (ipt<0) { NbrImages=-1; }
   if(ipt!=-1)
     {
     if (GetWordsInLine(UffFile,NbrWords,Mot)) { NbrImages=-1; }
     cout << "NbrWords=" << NbrWords << " WithRot=" << WithRot 
          << " Cplx=" << ComplexData << endl;
     if (NbrWords!=(3*(1+ComplexData)*(1+WithRot)))
      { NbrImages=-1; cout<<"Record 10: incorrect number of data\n"; }
     if (ComplexData)
       {
       DepReal[ipt].x=PosRef[ipt].x+atof(Mot[0]);
       DepImag[ipt].x=PosRef[ipt].x+atof(Mot[1]);
       DepReal[ipt].y=PosRef[ipt].y+atof(Mot[2]); 
       DepImag[ipt].y=PosRef[ipt].y+atof(Mot[3]);
       DepReal[ipt].z=PosRef[ipt].z+atof(Mot[4]); 
       DepImag[ipt].z=PosRef[ipt].z+atof(Mot[6]);
       }
     else
       {
       DepReal[ipt].x=PosRef[ipt].x+atof(Mot[0]);
       DepReal[ipt].y=PosRef[ipt].y+atof(Mot[1]);
       DepReal[ipt].z=PosRef[ipt].z+atof(Mot[2]);
       }
     if (fabs(DepReal[ipt].x-PosRef[ipt].x)>DepMax) 
       DepMax=fabs(DepReal[ipt].x-PosRef[ipt].x);
     if (fabs(DepReal[ipt].y-PosRef[ipt].y)>DepMax) 
       DepMax=fabs(DepReal[ipt].y-PosRef[ipt].y);
     if (fabs(DepReal[ipt].z-PosRef[ipt].z)>DepMax) 
       DepMax=fabs(DepReal[ipt].z-PosRef[ipt].z);
     if (fabs(DepImag[ipt].x-PosRef[ipt].x)>DepMax) 
       DepMax=fabs(DepImag[ipt].x-PosRef[ipt].x);
     if (fabs(DepImag[ipt].y-PosRef[ipt].y)>DepMax) 
       DepMax=fabs(DepImag[ipt].y-PosRef[ipt].y);
     if (fabs(DepImag[ipt].z-PosRef[ipt].z)>DepMax) 
       DepMax=fabs(DepImag[ipt].z-PosRef[ipt].z);
     }
   }
  }
 while ((NbrWords!=1) || (strcmp(Mot[0],"-1")!=0));

 UffFile.close();
 if (NbrImages==-1) 
   {
   // Puisque ca n'a pas marche, on remet tout a zero
     alpha=0.0; beta=0.0; freq=0.0; amor=0.0;
   for (ipt=0; ipt<NbrNodes; ipt++)
     { DepReal[ipt].x=PosRef[ipt].x; DepImag[ipt].x=PosRef[ipt].x; 
       DepReal[ipt].y=PosRef[ipt].y; DepImag[ipt].y=PosRef[ipt].y; 
       DepReal[ipt].z=PosRef[ipt].z; DepImag[ipt].z=PosRef[ipt].z; 
       node[ipt].x=PosRef[ipt].x; node[ipt].y=PosRef[ipt].y;
       node[ipt].z=PosRef[ipt].z; }
   NbrImages=0;
   return 2;
   }

 float lcaract;
 lcaract=sqrt((xmax-xmin)*(xmax-xmin)+(ymax-ymin)*(ymax-ymin)
               +(zmax-zmin)*(zmax-zmin));
 for (ipt=0;ipt<NbrNodes;ipt++)
      {
      node[ipt].x=PosRef[ipt].x
                +(DepReal[ipt].x-PosRef[ipt].x)*0.05*lcaract/DepMax;
      node[ipt].y=PosRef[ipt].y
                +(DepReal[ipt].y-PosRef[ipt].y)*0.05*lcaract/DepMax;
      node[ipt].z=PosRef[ipt].z
                +(DepReal[ipt].z-PosRef[ipt].z)*0.05*lcaract/DepMax;
      }
 CurImage=0;
 NbrImages=50;
 return 0;
}	

/*---------------------------------------------------------------------*/

void def_uff::Next()
// Construction de l'image modale
  {
  if (NbrImages)
    {
    int ipt;
    float lcaract, pi=3.14151695;
    lcaract=sqrt((xmax-xmin)*(xmax-xmin)+(ymax-ymin)*(ymax-ymin)
                +(zmax-zmin)*(zmax-zmin));
    CurImage=(CurImage++)%NbrImages;
    for (ipt=0;ipt<NbrNodes;ipt++)
      {
      node[ipt].x=PosRef[ipt].x
                 +(DepReal[ipt].x-PosRef[ipt].x)*0.05*lcaract
                       *cos(2.0*CurImage*pi/NbrImages)/DepMax
                 +(DepImag[ipt].x-PosRef[ipt].x)*0.05*lcaract
                       *sin(2.0*CurImage*pi/NbrImages)/DepMax;
      node[ipt].y=PosRef[ipt].y
                 +(DepReal[ipt].y-PosRef[ipt].y)*0.05*lcaract
                       *cos(2.0*CurImage*pi/NbrImages)/DepMax
                 +(DepImag[ipt].y-PosRef[ipt].y)*0.05*lcaract
                       *sin(2.0*CurImage*pi/NbrImages)/DepMax;
      node[ipt].z=PosRef[ipt].z
                 +(DepReal[ipt].z-PosRef[ipt].z)*0.05*lcaract
                       *cos(2.0*CurImage*pi/NbrImages)/DepMax
                 +(DepImag[ipt].z-PosRef[ipt].z)*0.05*lcaract
                       *sin(2.0*CurImage*pi/NbrImages)/DepMax;
      }
    }
  }

/*---------------------------------------------------------------------*/

void def_uff::LoadImage(int NumImage)
// Building a requested image
  {
  if (NbrImages) 
    {
    int ipt;
    float lcaract, pi=3.14151695;
    lcaract=sqrt((xmax-xmin)*(xmax-xmin)+(ymax-ymin)*(ymax-ymin)
                 +(zmax-zmin)*(zmax-zmin));
    CurImage=NumImage%NbrImages;
    for (ipt=0;ipt<NbrNodes;ipt++)
      {
      node[ipt].x=PosRef[ipt].x
                 +(DepReal[ipt].x-PosRef[ipt].x)*0.05*lcaract
                       *cos(2.0*CurImage*pi/NbrImages)/DepMax
                 +(DepImag[ipt].x-PosRef[ipt].x)*0.05*lcaract
                       *sin(2.0*CurImage*pi/NbrImages)/DepMax;
      node[ipt].y=PosRef[ipt].y
                 +(DepReal[ipt].y-PosRef[ipt].y)*0.05*lcaract
                       *cos(2.0*CurImage*pi/NbrImages)/DepMax
                 +(DepImag[ipt].y-PosRef[ipt].y)*0.05*lcaract
                       *sin(2.0*CurImage*pi/NbrImages)/DepMax;
      node[ipt].z=PosRef[ipt].z
                 +(DepReal[ipt].z-PosRef[ipt].z)*0.05*lcaract
                       *cos(2.0*CurImage*pi/NbrImages)/DepMax
                 +(DepImag[ipt].z-PosRef[ipt].z)*0.05*lcaract
                       *sin(2.0*CurImage*pi/NbrImages)/DepMax;
      }
    }
  }

/*---------------------------------------------------------------------*/
/*                    Methods relative to def_repere                   */
/*---------------------------------------------------------------------*/

/*---------------------------------------------------------------------*/

//dessine un repere a partir du fichier repere.vol
def_repere::def_repere()
  {
  NbrNodes=40;
  NbrEdges=31;
  xmin=-0.3; xmax=1.0;
  ymin=-0.3; ymax=1.0;
  zmin=-0.3; zmax=1.0;
  node = new NodeStruct[NbrNodes];
  if (node==NULL) ErreurMessage(_FATAL,"Memory allocation problem");
  edge= new EdgeStruct[NbrEdges];
  if (edge==NULL) ErreurMessage(_FATAL,"Memory allocation problem");
	
  node[0].x=0.0f; node[0].y= 0.0f; node[0].z= 0.0f;
  node[1].x=1.0f; node[1].y= 0.0f; node[1].z= 0.0f;
  node[2].x=0.8f; node[2].y= 0.1f; node[2].z= 0.0f;
  node[3].x=0.8f; node[3].y= 0.0f; node[3].z= 0.1f;
  node[4].x=0.8f; node[4].y=-0.1f; node[4].z= 0.0f;
  node[5].x=0.8f; node[5].y= 0.0f; node[5].z=-0.1f;
  node[6].x=1.0f; node[6].y=-0.3f; node[6].z= 0.0f;
  node[7].x=0.8f; node[7].y=-0.2f; node[7].z= 0.0f;
  node[8].x=0.8f; node[8].y=-0.3f; node[8].z= 0.0f;
  node[9].x=1.0f; node[9].y=-0.2f; node[9].z= 0.0f;
  node[10].x= 1.0f; node[10].y= 0.0f; node[10].z=-0.2f;
  node[11].x= 0.8f; node[11].y= 0.0f; node[11].z=-0.3f;
  node[12].x= 0.8f; node[12].y= 0.0f; node[12].z=-0.2f; 
  node[13].x= 1.0f; node[13].y= 0.0f; node[13].z=-0.3f; 
  node[14].x= 0.0f; node[14].y= 1.0f; node[14].z= 0.0f;
  node[15].x= 0.1f; node[15].y= 0.8f; node[15].z= 0.0f;
  node[16].x= 0.0f; node[16].y= 0.8f; node[16].z= 0.1f;
  node[17].x=-0.1f; node[17].y= 0.8f; node[17].z= 0.0f;
  node[18].x= 0.0f; node[18].y= 0.8f; node[18].z=-0.1f;
  node[19].x=-0.25f;node[19].y= 0.8f; node[19].z= 0.0f;
  node[20].x=-0.2f; node[20].y= 1.0f; node[20].z= 0.0f;
  node[21].x=-0.3f; node[21].y= 1.0f; node[21].z= 0.0f;
  node[22].x=-0.25f;node[22].y= 0.9f; node[22].z= 0.0f;
  node[23].x= 0.0f; node[23].y= 1.0f; node[23].z= -0.2f;
  node[24].x= 0.0f; node[24].y= 1.0f; node[24].z= -0.3f;
  node[25].x= 0.0f; node[25].y= 0.8f; node[25].z= -0.25f;
  node[26].x= 0.0f; node[26].y= 0.9f; node[26].z= -0.25f;
  node[27].x= 0.0f; node[27].y= 0.0f; node[27].z= 1.0f;
  node[28].x= 0.0f; node[28].y= -0.1f;node[28].z= 0.8f;
  node[29].x= -0.1f;node[29].y= 0.0f; node[29].z= 0.8f; 
  node[30].x= 0.0f; node[30].y= 0.1f; node[30].z= 0.8f;
  node[31].x= 0.1f; node[31].y= 0.0f; node[31].z= 0.8f; 
  node[32].x= 0.0f; node[32].y= -0.2f;node[32].z= 0.8f;
  node[33].x= 0.0f; node[33].y= -0.3f;node[33].z= 0.8f;
  node[34].x= 0.0f; node[34].y= -0.2f;node[34].z= 1.0f;
  node[35].x= 0.0f; node[35].y= -0.3f;node[35].z= 1.0f;
  node[36].x= -0.2f;node[36].y= 0.0f; node[36].z= 0.8f;
  node[37].x= -0.3f;node[37].y= 0.0f; node[37].z= 0.8f;
  node[38].x= -0.2f;node[38].y= 0.0f; node[38].z= 1.0f;
  node[39].x= -0.3f;node[39].y= 0.0f; node[39].z= 1.0f;

  edge[0].node1=0; edge[0].node2=1; edge[0].color=1;
  edge[1].node1=1; edge[1].node2=2; edge[1].color=1;
  edge[2].node1=1; edge[2].node2=3; edge[2].color=1;
  edge[3].node1=1; edge[3].node2=4; edge[3].color=1;
  edge[4].node1=1; edge[4].node2=5; edge[4].color=1;
  edge[5].node1=6; edge[5].node2=7; edge[5].color=1;
  edge[6].node1=8; edge[6].node2=9; edge[6].color=1;
  edge[7].node1=10;edge[7].node2=11; edge[7].color=1;
  edge[8].node1=12;edge[8].node2=13; edge[8].color=1;
  edge[9].node1=0; edge[9].node2=14; edge[9].color=4;
  edge[10].node1=14; edge[10].node2=15; edge[10].color=4;
  edge[11].node1=14; edge[11].node2=16; edge[11].color=4;
  edge[12].node1=14; edge[12].node2=17; edge[12].color=4; 
  edge[13].node1=14; edge[13].node2=18; edge[13].color=4; 
  edge[14].node1=19; edge[14].node2=22; edge[14].color=4; 
  edge[15].node1=22; edge[15].node2=20; edge[15].color=4; 
  edge[16].node1=22; edge[16].node2=21; edge[16].color=4; 
  edge[17].node1=25; edge[17].node2=26; edge[17].color=4; 
  edge[18].node1=26; edge[18].node2=23; edge[18].color=4; 
  edge[19].node1=26; edge[19].node2=24; edge[19].color=4; 
  edge[20].node1=0;  edge[20].node2=27; edge[20].color=5;
  edge[21].node1=27; edge[21].node2=28; edge[21].color=5;
  edge[22].node1=27; edge[22].node2=29; edge[22].color=5;
  edge[23].node1=27; edge[23].node2=30; edge[23].color=5;
  edge[24].node1=27; edge[24].node2=31; edge[24].color=5;
  edge[25].node1=32; edge[25].node2=33; edge[25].color=5;
  edge[26].node1=33; edge[26].node2=34; edge[26].color=5;
  edge[27].node1=34; edge[27].node2=35; edge[27].color=5;
  edge[28].node1=36; edge[28].node2=37; edge[28].color=5;
  edge[29].node1=37; edge[29].node2=38; edge[29].color=5;
  edge[30].node1=38; edge[30].node2=39; edge[30].color=5;
  }

/*---------------------------------------------------------------------*/

// routine de dessin
void def_repere::Trace()
  {
  //creation de la liste GL
  int i,j,k;
  GLfloat  mat_ambient[] = { 0.0f, 0.0f, 0.0f, 1.0f};
  GLfloat  mat_diffuse[] = { 0.0f, 0.0f, 0.0f, 1.0f};
  GLfloat  mat_gray[] = { 0.2f, 0.2f, 0.2f, 1.0f};
  GLfloat  mat_black[] = { 0.0f, 0.0f, 0.0f, 0.0f};

  if (!hidden)
    {
    glBegin(GL_LINES);
    for (i=0;i<NbrEdges;i++) 
      {
      mat_ambient[0]=COULEUR[edge[i].color][0];
      mat_ambient[1]=COULEUR[edge[i].color][1];
      mat_ambient[2]=COULEUR[edge[i].color][2];
      glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, mat_ambient);
      glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, mat_black);
      glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, mat_black);
      glVertex3f(node[edge[i].node1].x,node[edge[i].node1].y,
                 node[edge[i].node1].z);/*origine*/
      glVertex3f(node[edge[i].node2].x,node[edge[i].node2].y,
                 node[edge[i].node2].z);/*extremite*/
      }
    glEnd();

    glFlush();
    }
  }

/*---------------------------------------------------------------------*/

// routine de dessin
void def_repere::LoadImage(int NumImage)
  {
  }

/*---------------------------------------------------------------------*/

// routine pour l'affichage ou le masquage du repère
void repere() 
  {
  if (s1->Hidden()) s1->Show();
  else s1->Hide();
  }

/*---------------------------------------------------------------------*/

// routines de dimensionnement du repère
void repere_plus()
  {	
  s1->Homothetie(1.25f);	
  }

/*---------------------------------------------------------------------*/

void repere_moins() 
  {
  s1->Homothetie(0.8f);
  }

/*---------------------------------------------------------------------*/
					       
//===================>>> myOGLCanvasPane::myOGLCanvasPane <<<====================
// constructeur du canvas
  myOGLCanvasPane::myOGLCanvasPane(unsigned int vGLmode, PaneType pt)
    : vBaseGLCanvasPane(vGLmode,pt)
  {
    initDone = 0;
  }

//===================>>> myOGLCanvasPane::~myOGLCanvasPane <<<====================
// destructeur du canvas
  myOGLCanvasPane::~myOGLCanvasPane()
  {
  }

//======================>>> myOGLCanvasPane::TimerAnimate <<<========================
// instruction exécutées par le timer
  void myOGLCanvasPane::TimerAnimate(void)
  {
  // **** Called by CmdWindow AuxTimer for animation.
  vglMakeCurrent();  // Typically done here
  anime();
  vglFlush();  // After you draw, typically flush
  }

//======================>>> myOGLCanvasPane::graphicsInit <<<========================
// Initialmisation des graphiques OpenGL
// Exécuté une fois lorsque le canvas est affiché
  void myOGLCanvasPane::graphicsInit(void)
  {
    vBaseGLCanvasPane::graphicsInit();	// Always call the superclass first!

	vglMakeCurrent();  // Typically done here
	GLSetupRC();
        s1=new def_repere;
        s2=NULL;
        Add_Scene(s1);
	GLResize(GetWidth(),GetHeight());
	GLRenderScene();
	vglFlush();  // After you draw, typically flush
	initDone = 1;
  }

//======================>>> myOGLCanvasPane::MouseDown <<<======================
  void myOGLCanvasPane::MouseDown(int X, int Y, int button)
  {
    vBaseGLCanvasPane::MouseDown(X,Y,button);
  }

//========================>>> myOGLCanvasPane::MouseUp <<<======================
  void myOGLCanvasPane::MouseUp(int X, int Y, int button)
  {
    vBaseGLCanvasPane::MouseUp(X,Y,button);
  }

//======================>>> myOGLCanvasPane::MouseMove <<<======================
  void myOGLCanvasPane::MouseMove(int x, int y, int button)
  {
    vBaseGLCanvasPane::MouseMove(x,y,button);
  }

//=========================>>> myOGLCanvasPane::Redraw <<<======================
// Routiine exécutée lorsque le canvas doir être redessiné
  void myOGLCanvasPane::Redraw(int x, int y, int w, int h)
  {
    static int inRedraw = 0;

    if (inRedraw || !initDone)  // Don't draw until initialized
        return;

    inRedraw = 1;  // Don't allow recursive redraws.
    vglMakeCurrent();  // Typically done here
    GLRenderScene();
    inRedraw = 0;  // Out of Redraw

  }


//======================>>> myOGLCanvasPane::Resize <<<======================
// routine exécutée lorsque le canvas est redimensionné
  void myOGLCanvasPane::Resize(int w, int h)
  {
    vBaseGLCanvasPane::Resize(w,h);
	GLResize(w,h);
	GLRenderScene();
//	}
  }

