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

// EasyAnimGL.cpp: the core routines for viewing with GL

#include <stdio.h>
#include <stdlib.h>
#ifdef WINDOWS
  #include <windows.h>
  #include <wingdi.h>
#endif
#include <GL/gl.h>
#include <GL/glu.h>
#include <math.h>
#include <fstream>
#include <iostream>

#include "EasyAnimGL.h"

// Declarations and initializations
int i;

float PI =3.141592654f;
float PI_2 =1.570796327f;
float PI_10 = 0.3141592654f;

// Variables related to the view

GLfloat _EYE_X=2886.75f,_EYE_Y=2886.75f,_EYE_Z=2886.75f;
GLfloat _CENTER_X=0.0f,_CENTER_Y=0.0f,_CENTER_Z=0.0f;
GLfloat _THETA=PI*0.25f,_PHI=0.6154797;
GLfloat _GAMMA=0.0f;
GLfloat _UP_X=-0.408f,_UP_Y=-0.408f,_UP_Z=0.577f;
GLfloat _RIGHT_X=-0.707f,_RIGHT_Y=0.707f,_RIGHT_Z=0.0f;
GLfloat _DIST=5000.0f;
static GLsizei _GLOB_W=0;
static GLsizei _GLOB_H=0;
static GLfloat _LARGEUR=10;

GLfloat COULEUR[16][3];
GLfloat MATERIAU[16][3][3];

int SHOWNODES=0;
int SHOWEDGES=1;
int SHOWSIDES=1;
int _ORTHO=1;
int _COULEUR_FOND = 0;
static int _N_SKIP_IMAGE=0;

// Variables controlling aspect of edges and nodes
GLdouble NodeRadius=0.1;
GLdouble EdgeRadius=0.05;

// Variables controlling animation and type of view
int _ANIME = 0;
int _MULTI;

// Variables related to the scenes
def_scene *SCENE[10];
int n_scene=0;

// Routine to normalize a vector
void ReduceToUnit(float vector[3])
	{
	float length;
	
	// Calculate the length of the vector		
	length = (float)sqrt((vector[0]*vector[0]) + 
						(vector[1]*vector[1]) +
						(vector[2]*vector[2]));

	// Keep the program from blowing up by providing an exceptable
	// value for vectors that may calculated too close to zero.
	if(length == 0.0f)
		length = 1.0f;

	// Dividing each element by the length will result in a
	// unit normal vector.
	vector[0] /= length;
	vector[1] /= length;
	vector[2] /= length;
	}

/*---------------------------------------------------------------------*/

// routine de calcul du vecteur normal d'une surface
void calcNormal(float v[3][3], float out[3])
	{
	float v1[3],v2[3];
	static const int x = 0;
	static const int y = 1;
	static const int z = 2;

	// Calcul deux vecteurs à partir des trois points
	v1[x] = v[0][x] - v[1][x];
	v1[y] = v[0][y] - v[1][y];
	v1[z] = v[0][z] - v[1][z];

	v2[x] = v[1][x] - v[2][x];
	v2[y] = v[1][y] - v[2][y];
	v2[z] = v[1][z] - v[2][z];

	// Take the cross product of the two vectors to get
	// the normal vector which will be stored in out
	out[x] = v1[y]*v2[z] - v1[z]*v2[y];
	out[y] = v1[z]*v2[x] - v1[x]*v2[z];
	out[z] = v1[x]*v2[y] - v1[y]*v2[x];

	// normalise le vecteur
	ReduceToUnit(out);
	}

/*---------------------------------------------------------------------*/

/*---------------------------------------------------------------------*/
/*     Premiere partie: procedures generales relatives aux scenes      */
/*---------------------------------------------------------------------*/

// destructeur
def_scene::def_scene()
  {
  NbrImages=0;
  CurImage=0;
  hidden=0;
  NbrNodes=0;
  NbrEdges=0;
  NbrSides=0;
  node=NULL;
  edge=NULL;
  side=NULL;
  xmin=0.0;
  xmax=0.0;
  ymin=0.0;
  ymax=0.0;
  zmin=0.0;
  zmax=0.0;
  }	

/*---------------------------------------------------------------------*/

// destructeur
def_scene::~def_scene()
  {
  int i;
  if (NbrNodes) delete node;
  if (NbrEdges) delete edge;
  if (NbrSides)
    {
    for (i=0;i<NbrSides;i++) delete side[i].NumNode;
    delete side;
    }
  NbrNodes = 0;
  NbrEdges = 0;
  NbrSides = 0;
  }	

/*---------------------------------------------------------------------*/

// Lecture des fichiers
int def_scene::Set(int npoint, float **point,int narr, int **arrete, 
                    int nface, int **face)
  {
  int i,j,k,l;
  float long_max=0.0f;
  if (npoint<=0) return 1;
  NbrNodes=npoint;
  node=new NodeStruct[NbrNodes];
  if (node==NULL) return 2;
  xmin=point[0][0]; xmax=point[0][0];
  ymin=point[0][1]; ymax=point[0][1];
  zmin=point[0][2]; zmax=point[0][2];
  for (i=0;i<NbrNodes;i++)
	{
	node[i].x=point[i][0];
	node[i].y=point[i][1];
	node[i].z=point[i][2];
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

  // Managing the edges
  if (narr<=0) return 3;
  NbrEdges=narr;
  edge=new EdgeStruct[NbrEdges];
  if (edge==NULL) return 5;
  for (i=0;i<NbrEdges;i++)
      {
      edge[i].node1=arrete[i][0];edge[i].node2=arrete[i][1];
      edge[i].color=arrete[i][2];
      if ((edge[i].node1<0)||(edge[i].node1>NbrNodes)) return 6;
      if ((edge[i].node2<0)||(edge[i].node2>NbrNodes)) return 7;
      }
  if (nface<0) return 8;
  NbrSides=nface;
  if (NbrSides!=0)
      {
      side = new SideStruct[NbrSides];
      if (side==NULL) return 9;
      for (i=0;i<NbrSides;i++)
	  {
	  side[i].NbrNodes=face[i][0];/*nombre de sommets*/
          side[i].NumNode = new int [side[i].NbrNodes];
          if (side[i].NumNode==NULL) return 10;
          for (k=0;k<=side[i].NbrNodes;k++)
            {
            side[i].NumNode[k]=face[i][k+1];/*numero des sommets*/
            if ((side[i].NumNode[k]<0)||(side[i].NumNode[k]>NbrNodes)) 
               return 11;
            }
          side[i].color=face[i][face[i][0]+1];/*numero de la couleur*/
          }/*endfor i*/
      }/*endif*/
  }

/*---------------------------------------------------------------------*/

void def_scene::GetDim(float &x1,float &x2,float &y1,float &y2,
                       float &z1,float &z2)
  {
  x1=xmin; x2=xmax; y1=ymin; y2=ymax; z1=zmin; z2=zmax;
  }

/*---------------------------------------------------------------------*/

int def_scene::Hidden()	
  {
  return hidden;
  }

/*---------------------------------------------------------------------*/

void def_scene::Hide()	
  {
  hidden=1;
  }  

/*---------------------------------------------------------------------*/

void def_scene::Show()	
  {
  hidden=0;
  }  

/*---------------------------------------------------------------------*/

void def_scene::Homothetie(float MULT)	
// Homothétie de la scène
  {
  int i;
  //modification de la taille de l'image
  if (MULT>0.0f)
    {
    for (i=0;i<NbrNodes;i++)
	{
	node[i].x*=MULT;
	node[i].y*=MULT;
	node[i].z*=MULT;
	}
    xmin*=MULT; xmax*=MULT;
    ymin*=MULT; ymax*=MULT;
    zmin*=MULT; zmax*=MULT;
    }
	
}

/*---------------------------------------------------------------------*/

int def_scene::GetNumNodeOfName(const char *NodeName)
  {
  int i,num=-1;
  for (i=0;i<NbrNodes;i++) if (!strcmp(NodeName,node[i].name)) num=i;
  return num;
  }

/*---------------------------------------------------------------------*/

//coordonnees d'un point de l'image
/*void def_scene::Get_Point(int num,float Pt[3])	
  {
  Pt[0]=0.0;
  Pt[1]=0.0;
  Pt[2]=0.0;
  if ((num>=0)&&(num<NbrNodes))
      {
      Pt[0]=POINT[num][0];
      Pt[1]=POINT[num][1];
      Pt[2]=POINT[num][2];
      }
      }*/

/*---------------------------------------------------------------------*/

int def_scene::GetNbrImages()	
  {
    return NbrImages;
  }

/*---------------------------------------------------------------------*/

int def_scene::GetCurImage()	
  {
    return CurImage;
  }

/*---------------------------------------------------------------------*/

void def_scene::Next()	
  {
  }

/*---------------------------------------------------------------------*/

void def_scene::LoadImage(int NumImage)	
  {
  }

/*---------------------------------------------------------------------*/

// routine de dessin
void def_scene::Trace()
  {
  //creation de la liste GL
  int i,j,k;
  GLfloat  mat_ambient[] = { 0.0f, 0.0f, 0.0f, 1.0f};
  GLfloat  mat_diffuse[] = { 0.0f, 0.0f, 0.0f, 1.0f};
  GLfloat  mat_gray[] = { 0.2f, 0.2f, 0.2f, 1.0f};
  GLfloat  mat_lightgray[] = { 0.4f, 0.4f, 0.4f, 0.0f};
  GLfloat  mat_black[] = { 0.0f, 0.0f, 0.0f, 0.0f};

  if (!hidden)
    {
    // Display the spheres at the nodes
    if (SHOWNODES)
	{
        glPushMatrix();
	mat_ambient[0]=0.2f; mat_ambient[1]=0.2f; mat_ambient[2]=0.6f;
        glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, mat_ambient);
	mat_diffuse[0]=0.2f; mat_diffuse[1]=0.2f; mat_diffuse[2]=0.6f;
        glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, mat_diffuse);
        glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, mat_lightgray);
        for (i=0; i<NbrNodes; i++)
          {
          glTranslatef(node[i].x,node[i].y,node[i].z);
          gluSphere(NodeSphere,NodeRadius,10,10);
          glTranslatef(-node[i].x,-node[i].y,-node[i].z);
          }
        glPopMatrix();
	}
      //glBegin(GL_LINES);
    if (SHOWEDGES)
        {
        float angle, dx, dy, dz, le;
        for (i=0;i<NbrEdges;i++) 
          {
	  ///mat_ambient[0]=COULEUR[edge[i].color][0];
	  ///mat_ambient[1]=COULEUR[edge[i].color][1];
	  ///mat_ambient[2]=COULEUR[edge[i].color][2];
          mat_ambient[0]=0.6*COULEUR[edge[i].color][0];
          mat_ambient[1]=0.6*COULEUR[edge[i].color][1];
          mat_ambient[2]=0.6*COULEUR[edge[i].color][2];
          mat_diffuse[0]=0.4*COULEUR[edge[i].color][0];
          mat_diffuse[1]=0.4*COULEUR[edge[i].color][1];
          mat_diffuse[2]=0.4*COULEUR[edge[i].color][2];
          glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, mat_ambient);
          //glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, mat_black);
          glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, mat_diffuse);
          glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, mat_black);
          //glVertex3f(node[edge[i].node1].x,node[edge[i].node1].y,
          //          node[edge[i].node1].z);/*origine*/
          //glVertex3f(node[edge[i].node2].x,node[edge[i].node2].y,
          //           node[edge[i].node2].z);/*extremite*/
          dx=node[edge[i].node2].x-node[edge[i].node1].x;
          dy=node[edge[i].node2].y-node[edge[i].node1].y;
          dz=node[edge[i].node2].z-node[edge[i].node1].z;
          le=sqrt(dx*dx+dy*dy+dz*dz);
          if (le>1E-6)
	    {
            glPushMatrix();
            glTranslatef(node[edge[i].node1].x,node[edge[i].node1].y,
                         node[edge[i].node1].z);
            if (fabs(dz/le)>0.99999) 
                {
		angle=90*(1-dz/le);
		glRotatef(angle,1.0,0.0,0.0);
                }
            else
	        {
		angle=57.29578*acos(dz/le);
                if (dy<0) glRotatef(angle,-dy,dx,0);
                else glRotatef(-angle,dy,-dx,0);
                }
            gluCylinder(NodeSphere, EdgeRadius, EdgeRadius,le,8,1);
            glPopMatrix();
            }
	  }
	}

    //glEnd();
    if (SHOWSIDES)
      {
      float normal[3];
      for (i=0;i<NbrSides;i++) 
        {
        mat_ambient[0]=0.6*COULEUR[side[i].color][0];
        mat_ambient[1]=0.6*COULEUR[side[i].color][1];
        mat_ambient[2]=0.6*COULEUR[side[i].color][2];
        mat_diffuse[0]=0.4*COULEUR[side[i].color][0];
        mat_diffuse[1]=0.4*COULEUR[side[i].color][1];
        mat_diffuse[2]=0.4*COULEUR[side[i].color][2];
        glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, mat_ambient);
        glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, mat_diffuse);
        glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, mat_gray);
        glBegin(GL_POLYGON);
        float v[3][3]=
        {{node[side[i].NumNode[0]].x,node[side[i].NumNode[0]].y,
          node[side[i].NumNode[0]].z},
         {node[side[i].NumNode[1]].x,node[side[i].NumNode[1]].y,
          node[side[i].NumNode[1]].z},
         {node[side[i].NumNode[2]].x,node[side[i].NumNode[2]].y,
	  node[side[i].NumNode[2]].z}};
        calcNormal(v,normal);
        for (k=0;k<side[i].NbrNodes;k++)
	  {
          glNormal3fv(normal);
	  glVertex3f(node[side[i].NumNode[k]].x,node[side[i].NumNode[k]].y,
                     node[side[i].NumNode[k]].z);/*sommet*/
	  }
        glEnd();
        }
      }
    glFlush();
    }
  }

/*---------------------------------------------------------------------*/

void Add_Scene(def_scene *sc)
  {
  n_scene++;
  SCENE[n_scene-1]=sc;
  }

/*---------------------------------------------------------------------*/

void Remove_Scene()
  {
  if (n_scene)
    {
    n_scene--;
    //SCENE[n_scene]=0;
    }
  }

/*---------------------------------------------------------------------*/
/* Deuxieme partie: procedures relatives a la visualisation des scenes */
/*---------------------------------------------------------------------*/

//initialisation des couleurs;
void Init_Couleur()
{
	//Noir			0
	COULEUR[0][0]=COULEUR[0][1]=COULEUR[0][2]=0.0f;
	//Bleu			1
	COULEUR[1][0]=COULEUR[1][1]=0.0f;COULEUR[1][2]=1.0f;
	//Vert			2
	COULEUR[2][0]=0.0f;COULEUR[2][1]=1.0f;COULEUR[2][2]=0.0f;
	//Cyan			3
	COULEUR[3][0]=0.0f;COULEUR[3][1]=1.0f;COULEUR[3][2]=1.0f;
	//Rouge			4
	COULEUR[4][0]=1.0f;COULEUR[4][1]=0.0f;COULEUR[4][2]=0.0f;
	//magenta		5
	COULEUR[5][0]=1.0f;COULEUR[5][1]=0.0f;COULEUR[5][2]=1.0f;
	//Brun			6 
	COULEUR[6][0]=0.50f;COULEUR[6][1]=0.2f;COULEUR[6][2]=0.0f;
	//Gris clair	7
	COULEUR[7][0]=0.8f;COULEUR[7][1]=0.8f;COULEUR[7][2]=0.8f;
	//Gris foncé	8
	COULEUR[8][0]=0.25f;COULEUR[8][1]=0.25f;COULEUR[8][2]=0.25f;
	//Bleu clair	9
	COULEUR[9][0]=0.0f;COULEUR[9][1]=0.5f;COULEUR[9][2]=1.0f;
	//vert clair	10
	COULEUR[10][0]=0.0f;COULEUR[10][1]=1.0f;COULEUR[10][2]=0.5f;
	//Cyan clair	11
	COULEUR[11][0]=0.5f;COULEUR[11][1]=1.0f;COULEUR[11][2]=1.0f;
	//rouge clair	12
	 COULEUR[12][0]=1.0f;COULEUR[12][1]=0.5f;COULEUR[12][2]=0.5f;
	//Magenta clair	13
	COULEUR[13][0]=1.0f;COULEUR[13][1]=0.5f;COULEUR[13][2]=1.0f;
	//Jaune		14
	COULEUR[14][0]=1.0f;COULEUR[14][1]=1.0f;COULEUR[14][2]=0.0f;
	//Blanc		15
	COULEUR[15][0]=1.0f;COULEUR[15][1]=1.0f;COULEUR[15][2]=1.0f;
}

/*---------------------------------------------------------------------*/

// changement de la couleur du canvas
void Change_Couleur_Fond(int c)
{
	if ((c<0)||(c>15)) 
	{
		std::cout<<"couleur de fond illicite";
		return;
	}
	glClearColor(COULEUR[c][0],COULEUR[c][1],COULEUR[c][2],0.0f);
	_COULEUR_FOND=c;
}										   

/*---------------------------------------------------------------------*/

// Routine de calcul des nouvelles coordonnées
void rotation_phi_theta()
  {
  float cst,csf,snt,snf,csr,snr;
  cst=(float)cos((double)_THETA);snt=(float)sin((double)_THETA);
  csf=(float)cos((double)_PHI);snf=(float)sin((double)_PHI);
  _EYE_X=_CENTER_X+cst*csf*_DIST;
  _EYE_Y=_CENTER_Y+csf*snt*_DIST;
  _EYE_Z=_CENTER_Z+snf*_DIST;
  csr=(float)cos((double)_GAMMA);snr=(float)sin((double)_GAMMA);
  _UP_X=snr*snt-csr*cst*snf;
  _UP_Y=-snr*cst-csr*snf*snt;
  _UP_Z=csr*csf;
  _RIGHT_X=-snt*csr-cst*snf*snr;
  _RIGHT_Y=cst*csr-snt*snf*snr;
  _RIGHT_Z=csf*snr;
  }

/*---------------------------------------------------------------------*/

// Routine contenant les initialisations OpenGL
void GLSetupRC() 
  {
  Init_Couleur();			//initialisation des couleurs
  rotation_phi_theta();
  glEnable(GL_DEPTH_TEST);//permet l'utilisation du Z buffer
  GLint depth;
  /*std::ofstream ftmp("del.del");
  glGetIntegerv(GL_DEPTH_BITS,&depth);
  ftmp << "depth=" << depth << "\n";
  glGetIntegerv(GL_RED_BITS,&depth);
  ftmp << "red_bits=" << depth << "\n";
  glGetIntegerv(GL_GREEN_BITS,&depth);
  ftmp << "green_bits=" << depth << "\n";
  glGetIntegerv(GL_BLUE_BITS,&depth);
  ftmp << "blue_bits=" << depth << "\n";
  glGetIntegerv(GL_ALPHA_BITS,&depth);
  ftmp << "alpha_bits=" << depth << "\n";
  ftmp.close();*/


  // Activation des valeurs pour les lumieres
  // Par defaut, toutes les composantes de la lumiere sont blanches
  glEnable(GL_LIGHTING);
  GLfloat  ambientLight[] = { 1.0f, 1.0f, 1.0f, 1.0f };
  glLightfv(GL_LIGHT0,GL_AMBIENT,ambientLight);
  GLfloat  diffuseLight[] = { 1.0f, 1.0f, 1.0f, 1.0f };
  glLightfv(GL_LIGHT0,GL_DIFFUSE,diffuseLight);
  GLfloat  specular[] = { 1.0f, 1.0f, 1.0f, 0.6f};
  glLightfv(GL_LIGHT0,GL_SPECULAR,specular);
  GLfloat  lightPos[] = { 0.0f, 150.0f, 150.0f, 0.0f };
  glLightfv(GL_LIGHT0,GL_POSITION,lightPos);
  glEnable(GL_LIGHT0);

  // Pour que les deux faces soient traitees sur le meme pied
  // glMightModeli(LIGHT_MODEL_TWO_SIDE, GL_TRUE);

  // Defintion de la brillance (les couleurs sont gerees durant le trace)
  glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 50.0);

  glClearColor (0.0f, 0.0f, 0.0f, 0.0f);
  _COULEUR_FOND=0;
  glShadeModel (GL_FLAT);
  } 

/*---------------------------------------------------------------------*/

// Routine de dessin
void GLRenderScene()
  {
  // declarations
  int i;
  //efface l'ecran et le buffer Z
  glClear (GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
  glLoadIdentity ();             // Nettoie la matrice courante
  // transformation de vue
  gluLookAt(_EYE_X,_EYE_Y,_EYE_Z,_CENTER_X,_CENTER_Y,_CENTER_Z,
            _UP_X,_UP_Y,_UP_Z); 

  // Dessin de la scène
  for (i=0;i<n_scene;i++)
    {
    if (_MULTI!=0) // pour les vues multiples
	{
	// Vue globale 
	glViewport(0,0,_GLOB_W/2,_GLOB_H/2); //changement de la zone de dessin
	glLoadIdentity(); //Nettoyage de la matrice courante
        // transformation de vue
        gluLookAt(_EYE_X,_EYE_Y,_EYE_Z,_CENTER_X,_CENTER_Y,_CENTER_Z,
                  _UP_X,_UP_Y,_UP_Z);
	SCENE[i]->Trace(); // tracé de la scène
        // Vue selon X
	glViewport(0,_GLOB_H/2,_GLOB_W/2,_GLOB_H/2);
	glLoadIdentity();
        gluLookAt(_CENTER_X+_DIST,_CENTER_Y,_CENTER_Z,_CENTER_X,_CENTER_Y,
                  _CENTER_Z,0.0,0.0,1.0);
	SCENE[i]->Trace();
	// Vue selon Y
	glViewport(_GLOB_W/2,_GLOB_H/2,_GLOB_W/2,_GLOB_H/2);
	glLoadIdentity();
        gluLookAt(_CENTER_X,_CENTER_Y+_DIST,_CENTER_Z,_CENTER_X,_CENTER_Y,
                  _CENTER_Z,0.0,0.0,1.0);
	SCENE[i]->Trace();
	// Vue selon Z
	glViewport(_GLOB_W/2,0,_GLOB_W/2,_GLOB_H/2);
	glLoadIdentity();
        gluLookAt(_CENTER_X,_CENTER_Y,_CENTER_Z+_DIST,_CENTER_X,_CENTER_Y,
                  _CENTER_Z,0.0,1.0,0.0);
	SCENE[i]->Trace();
	// retour a la vue de départ
	glLoadIdentity();
        gluLookAt(_EYE_X,_EYE_Y,_EYE_Z,_CENTER_X,_CENTER_Y,_CENTER_Z,_UP_X,
                  _UP_Y,_UP_Z);
	}
    else // une seule vue 
	{
        // zone de vue correspondant à tout l'écran
	glViewport(0, 0, _GLOB_W, _GLOB_H); 
	SCENE[i]->Trace();
        }
     //    glPopMatrix();//on va chercher la matrice
}

glFlush(); // oblige le tracé de la scène et du repère avant de continuer    

}

/*---------------------------------------------------------------------*/

// Routine utilisée lors du redimensionnement du canvas
void GLResize(GLsizei w, GLsizei h)
  {
  _GLOB_W=w;
  _GLOB_H=h;
  ortho_perspective();
  }

/*---------------------------------------------------------------------*/

// Routine passant de la projection parallèle à la projection avec perspective et inversément
void ortho_perspective() 
  {
  glMatrixMode (GL_PROJECTION);
  glLoadIdentity ();
  if (_ORTHO)
	{
	glOrtho(-_LARGEUR/2,_LARGEUR/2,-_LARGEUR*_GLOB_H/(_GLOB_W*2),
                _LARGEUR*_GLOB_H/(_GLOB_W*2),0.5*_DIST,1.5*_DIST);
	}
    else
	{
	glFrustum(-_LARGEUR/2,_LARGEUR/2,-_LARGEUR*_GLOB_H/(_GLOB_W*2),
                          _LARGEUR*_GLOB_H/(_GLOB_W*2),0.5*_DIST,1.5*_DIST);
	}
  glMatrixMode (GL_MODELVIEW);
  }

/*---------------------------------------------------------------------*/

// Routines de zoom

void zoom_plus() 
  {
  _LARGEUR*=1.1;
  ortho_perspective();
  }

/*---------------------------------------------------------------------*/

void zoom_moins() 
  {
  _LARGEUR/=1.1;
  ortho_perspective();
  }

/*---------------------------------------------------------------------*/

//Routines de changement de vue

void vue_selon_x_plus() 
  {
  _THETA=0.0f;_PHI=0.0f;_GAMMA=0; // affectation des nouvelles valeurs
  rotation_phi_theta(); // calcul des nouvelles coordonnées
  }

/*---------------------------------------------------------------------*/

void vue_selon_x_moins() 
  {
  _THETA=PI;_PHI=0.0;_GAMMA=0;
  rotation_phi_theta();
  }

/*---------------------------------------------------------------------*/

void vue_selon_y_plus() 
  {
  _THETA=PI_2;_PHI=0.0f;_GAMMA=0;
  rotation_phi_theta();
  }

/*---------------------------------------------------------------------*/

void vue_selon_y_moins() 
  {
  _THETA=-PI_2;_PHI=0.0f;_GAMMA=0;
  rotation_phi_theta();
  }

/*---------------------------------------------------------------------*/

void vue_selon_z_plus() 
  {
  _THETA=0.0f;_PHI=PI_2;_GAMMA=-PI_2;
  rotation_phi_theta();
  }

/*---------------------------------------------------------------------*/

void vue_selon_z_moins() 
  {
  _THETA=0.0f;_PHI=-PI_2;_GAMMA=-PI_2;
  rotation_phi_theta();
  }

/*---------------------------------------------------------------------*/

// routines pour les déplacements de la caméra
void distance_plus() 
  {
  _DIST*=1.1;
  rotation_phi_theta();
  }

/*---------------------------------------------------------------------*/

void distance_moins() 
  {
  if (_DIST>1.0) _DIST/=1.1;
  rotation_phi_theta();
  }

/*---------------------------------------------------------------------*/

void distance_plus_plus()
  {
  _DIST*=2;
  rotation_phi_theta();
  }

/*---------------------------------------------------------------------*/

void distance_moins_moins()
  {
  _DIST/=2;
  rotation_phi_theta();
  }

/*---------------------------------------------------------------------*/

// routines pour un changement de vue plus fin
void phi_plus() 
  {
  _PHI+=PI_10;
  rotation_phi_theta();
  }

/*---------------------------------------------------------------------*/

void phi_moins() 
  {
  _PHI-=PI_10;
  rotation_phi_theta();
  }

/*---------------------------------------------------------------------*/

void theta_plus() 
  {
  _THETA+=PI_10;
  rotation_phi_theta();
  }

/*---------------------------------------------------------------------*/

void theta_moins() 
  {
  _THETA-=PI_10;
  rotation_phi_theta();
  }

/*---------------------------------------------------------------------*/

void gamma_plus() 
  {
  _GAMMA+=PI_10;
  rotation_phi_theta();
  }

/*---------------------------------------------------------------------*/

void gamma_moins() 
  {
  _GAMMA-=PI_10;
  rotation_phi_theta();
  }

/*---------------------------------------------------------------------*/

// routine pour definir le rayon des spheres representant les noeuds
void SetNodeRadius(double NewRadius)
  {
  NodeRadius=NewRadius;
  }

/*---------------------------------------------------------------------*/

// routine pour multiplier le rayon des spheres representant les noeuds
void MultNodeRadius(double Factor)
  {
  NodeRadius*=Factor;
  }

/*---------------------------------------------------------------------*/

// routine pour definir le rayon des cylindres representant les lignes
void SetEdgeRadius(double NewRadius)
  {
  EdgeRadius=NewRadius;
  }

/*---------------------------------------------------------------------*/

// routine pour multiplier le rayon des cylindres representant les lignes
void MultEdgeRadius(double Factor)
  {
  EdgeRadius*=Factor;
  }

/*---------------------------------------------------------------------*/

// routine pour l'affichage ou le masquage des noeuds
void SetShowNodesOn()
  {
  SHOWNODES=1;
  }

/*---------------------------------------------------------------------*/

// routine pour l'affichage ou le masquage des noeuds
void SetShowNodesOff()
  {
  SHOWNODES=0;
  }

/*---------------------------------------------------------------------*/

// routine pour l'affichage ou le masquage des noeuds
void SwitchShowNodes()
  {
  if (SHOWNODES==0) SHOWNODES=1;
  else SHOWNODES=0;
  GLRenderScene();
  }

/*---------------------------------------------------------------------*/

// routine renvoyant le drapeau d'affichage des noeuds
int GetShowNodes()
  {
  return SHOWNODES;
  }

/*---------------------------------------------------------------------*/

// routine pour l'affichage ou le masquage des aretes
void SetShowEdgesOn()
  {
  SHOWEDGES=1;
  }

/*---------------------------------------------------------------------*/

// routine pour l'affichage ou le masquage des aretes
void SetShowEdgesOff()
  {
  SHOWEDGES=0;
  }

/*---------------------------------------------------------------------*/

// routine pour l'affichage ou le masquage des aretes
void SwitchShowEdges()
  {
  if (SHOWEDGES==0) SHOWEDGES=1;
  else SHOWEDGES=0;
  GLRenderScene();
  }

/*---------------------------------------------------------------------*/

// routine renvoyant le drapeau d'affichage des aretes
int GetShowEdges()
  {
  return SHOWEDGES;
  }

/*---------------------------------------------------------------------*/

// routine pour l'affichage ou le masquage des faces
void SetShowSidesOn()
  {
  SHOWSIDES=1;
  }

/*---------------------------------------------------------------------*/

// routine pour l'affichage ou le masquage des faces
void SetShowSidesOff()
  {
  SHOWSIDES=0;
  }

/*---------------------------------------------------------------------*/

// routine pour l'affichage ou le masquage des faces
void SwitchShowSides()
  {
  if (SHOWSIDES==0) SHOWSIDES=1;
  else SHOWSIDES=0;
  GLRenderScene();
  }

/*---------------------------------------------------------------------*/

// routine renvoyant le drapeau d'affichage des faces
int GetShowSides()
  {
  return SHOWSIDES;
  }

/*---------------------------------------------------------------------*/

// routines de déplacement de la scène par rapport aux axes
void depl_x_plus()
  {
  _CENTER_X-=0.1f*(float)sqrt(_LARGEUR);
  GLRenderScene();
  }

/*---------------------------------------------------------------------*/

void depl_x_moins()
  {
  _CENTER_X+=0.1f*(float)sqrt(_LARGEUR);
  GLRenderScene();
  }

/*---------------------------------------------------------------------*/

void depl_y_plus()
  {_CENTER_Y-=0.1f*(float)sqrt(_LARGEUR);
  GLRenderScene();
  }

/*---------------------------------------------------------------------*/

void depl_y_moins()
  {
  _CENTER_Y+=0.1f*(float)sqrt(_LARGEUR);
  GLRenderScene();
  }

/*---------------------------------------------------------------------*/

void depl_z_plus()
  {
  _CENTER_Z-=0.1f*(float)sqrt(_LARGEUR);
  GLRenderScene();
  }

/*---------------------------------------------------------------------*/

void depl_z_moins()
  {
  _CENTER_Z+=0.1f*(float)sqrt(_LARGEUR);
  GLRenderScene();
  }

/*---------------------------------------------------------------------*/

// routines de déplacement de la scène par rapport à l'écran 
void depl_haut()
  {
  _EYE_X-=_UP_X*(_LARGEUR)*0.02;
  _EYE_Y-=_UP_Y*(_LARGEUR)*0.02;
  _EYE_Z-=_UP_Z*(_LARGEUR)*0.02;
  _CENTER_X-=_UP_X*(_LARGEUR)*0.02;
  _CENTER_Y-=_UP_Y*(_LARGEUR)*0.02;
  _CENTER_Z-=_UP_Z*(_LARGEUR)*0.02;
  GLRenderScene();
  }

/*---------------------------------------------------------------------*/

void depl_bas()
  {
  _EYE_X+=_UP_X*(_LARGEUR)*0.02;
  _EYE_Y+=_UP_Y*(_LARGEUR)*0.02;
  _EYE_Z+=_UP_Z*(_LARGEUR)*0.02;
  _CENTER_X+=_UP_X*(_LARGEUR)*0.02;
  _CENTER_Y+=_UP_Y*(_LARGEUR)*0.02;
  _CENTER_Z+=_UP_Z*(_LARGEUR)*0.02;
  GLRenderScene();
  }

/*---------------------------------------------------------------------*/

void depl_gauche()
  {
  _EYE_X+=_RIGHT_X*(_LARGEUR)*0.02;
  _EYE_Y+=_RIGHT_Y*(_LARGEUR)*0.02;
  _EYE_Z+=_RIGHT_Z*(_LARGEUR)*0.02;
  _CENTER_X+=_RIGHT_X*(_LARGEUR)*0.02;
  _CENTER_Y+=_RIGHT_Y*(_LARGEUR)*0.02;
  _CENTER_Z+=_RIGHT_Z*(_LARGEUR)*0.02;
  GLRenderScene();
  }

/*---------------------------------------------------------------------*/

void depl_droite()
  {
  _EYE_X-=_RIGHT_X*(_LARGEUR)*0.02;
  _EYE_Y-=_RIGHT_Y*(_LARGEUR)*0.02;
  _EYE_Z-=_RIGHT_Z*(_LARGEUR)*0.02;
  _CENTER_X-=_RIGHT_X*(_LARGEUR)*0.02;
  _CENTER_Y-=_RIGHT_Y*(_LARGEUR)*0.02;
  _CENTER_Z-=_RIGHT_Z*(_LARGEUR)*0.02;
  GLRenderScene();
  }

/*---------------------------------------------------------------------*/

// routine pour l'animation
void anime()
 { 
 int j,i;
 if (_ANIME)
     {
     for (j=1;j<n_scene;j++) for (i=-1;i<_N_SKIP_IMAGE;i++) SCENE[j]->Next();
     GLRenderScene();
     }
  }

/*---------------------------------------------------------------------*/

void start_anime()
  {
  _ANIME=1;
  }

/*---------------------------------------------------------------------*/

void stop_anime()
  {
  _ANIME=0;
  }

/*---------------------------------------------------------------------*/

// routine de démarrage et d'arrêt de l'animation
void start_stop_anime()
  {
  if(_ANIME) _ANIME=0;
  else _ANIME=1;
  }

/*---------------------------------------------------------------------*/

void SetOneView()
  {
  _MULTI=0;
  }

/*---------------------------------------------------------------------*/

void SetMultiView()
  {
  _MULTI=1;
  }

/*---------------------------------------------------------------------*/

// routine de démarrage et d'arrêt de l'animation
void SwitchView()
  {
  if(_MULTI) _MULTI=0;
  else _MULTI=1;
  }

/*---------------------------------------------------------------------*/

// routine pour passer à l'image suivante
void image_suivante()
  {
  int j;
  _ANIME=0;
  for (j=1;j<n_scene;j++) SCENE[j]->Next();
  GLRenderScene();
  }

/*---------------------------------------------------------------------*/

// routine pour passer à l'image suivante
void image_precedente()
  {
  int j;
  _ANIME=0;
  for (j=1;j<n_scene;j++) if (SCENE[j]->GetCurImage()>1) 
      SCENE[j]->LoadImage(SCENE[j]->GetCurImage()-1);
  GLRenderScene();
  }

/*---------------------------------------------------------------------*/

// routine utilisée pour revenir en arrière dans l'animation

// routines changeant la vitesse de l'animation
void skip_image_plus()
  {
  _N_SKIP_IMAGE+=2; // diminue le nombre d'images affichées
  }

/*---------------------------------------------------------------------*/

void skip_image_moins()
  {
  _N_SKIP_IMAGE-=2; // augmente le nombre d'images affichées
  if(_N_SKIP_IMAGE<0) _N_SKIP_IMAGE=0;
  }

/*---------------------------------------------------------------------*/

// routine permettant de revenir à la vue de départ
void vue_depart()
  {
  _THETA= 0.785398f; _PHI= 0.6154797f; _GAMMA= 0.0f;
  rotation_phi_theta();
  }

/*---------------------------------------------------------------------*/

void vue_centre()
  {
  if (n_scene)
    {
    int j;
    float xmin,xmax,ymin,ymax,zmin,zmax;
    float xmin2,xmax2,ymin2,ymax2,zmin2,zmax2;
    SCENE[0]->GetDim(xmin,xmax,ymin,ymax,zmin,zmax);
    for (j=1; j<n_scene; j++)
      {
      SCENE[j]->GetDim(xmin2,xmax2,ymin2,ymax2,zmin2,zmax2);
      if (xmin2<xmin) xmin=xmin2;
      if (xmax2>xmax) xmax=xmax2;
      if (ymin2<ymin) ymin=ymin2;
      if (ymax2>ymax) ymax=ymax2;
      if (zmin2<zmin) zmin=zmin2;
      if (zmax2>zmax) zmax=zmax2;
      }
    _LARGEUR=sqrt((xmax-xmin)*(xmax-xmin)+(ymax-ymin)*(ymax-ymin)
                  +(zmax-zmin)*(zmax-zmin));
    ortho_perspective();
    _CENTER_X=0.5*(xmax+xmin);
    _CENTER_Y=0.5*(ymax+ymin);
    _CENTER_Z=0.5*(zmax+zmin);
    rotation_phi_theta();
    GLRenderScene();
    }
  }

/*---------------------------------------------------------------------*/

// routine de sauvegarde de configuration
int Sauve_Config(const char *nom)
  {
  if (nom==NULL) return 1;
  std::ofstream fout(nom);
  if(!fout) return 2;
  fout<<"_ORTHO " << _ORTHO<< " SHOWNODES " << SHOWNODES
      <<" SHOWEDGES " << SHOWEDGES << " SHOWSIDES " << SHOWSIDES << "\n";
  fout << "NodeRadius " << NodeRadius << " EdgeRadius " << EdgeRadius << "\n";
  fout<<"_CENTER_X "<<_CENTER_X<<" _CENTER_Y "<<_CENTER_Y
      <<" _CENTER_Z "<<_CENTER_Z<<"\n";
  fout<<"_THETA "<<_THETA<<" _PHI "<<_PHI<<" _GAMMA "<<_GAMMA<<"\n";
  fout<<"_DIST "<<_DIST<<"\n";
  fout<<" _LARGEUR "<<_LARGEUR<<"\n";
  fout<<"_N_SKIP_IMAGE "<<_N_SKIP_IMAGE<<"\n";
  fout<<"_COULEUR_FOND "<<_COULEUR_FOND<<"\n";
  fout.close();
  return 0;
  }

/*---------------------------------------------------------------------*/

// routine de lecture d'une configuration
int Lecture_Config(const char *nom_fichier)
  {
  char mot[1000];
  int _FILAIRE, OK;
  std::ifstream fin(nom_fichier);
  if(!fin) return 1;
  _FILAIRE=2;
  while (!fin.rdstate())
      {
      fin >> mot;
      if (!fin.rdstate())
        {
	OK=0;
        if (!strcmp(mot,"_ORTHO")) { fin >> _ORTHO; OK=1; }
        if (!strcmp(mot,"_FILAIRE")) { fin >> _FILAIRE; OK=1; }
        if (!strcmp(mot,"SHOWNODES")) { fin >> SHOWNODES; OK=1; }
        if (!strcmp(mot,"SHOWEDGES")) { fin >> SHOWEDGES; OK=1; }
        if (!strcmp(mot,"SHOWSIDES")) { fin >> SHOWSIDES; OK=1; }
        if (!strcmp(mot,"NodeRadius")) { fin >> NodeRadius; OK=1; }
        if (!strcmp(mot,"EdgeRadius")) { fin >> EdgeRadius; OK=1; }
        if (!strcmp(mot,"_CENTER_X")) { fin >> _CENTER_X; OK=1; }
        if (!strcmp(mot,"_CENTER_Y")) { fin >> _CENTER_Y; OK=1; }
        if (!strcmp(mot,"_CENTER_Z")) { fin >> _CENTER_Z; OK=1; }
        if (!strcmp(mot,"_THETA")) { fin >> _THETA; OK=1; }
        if (!strcmp(mot,"_PHI")) { fin >> _PHI; OK=1; }
        if (!strcmp(mot,"_GAMMA")) { fin >> _GAMMA; OK=1; }
        if (!strcmp(mot,"_DIST")) { fin >> _DIST; OK=1; }
        if (!strcmp(mot,"_LARGEUR")) { fin >> _LARGEUR; OK=1; }
        if (!strcmp(mot,"_N_SKIP_IMAGE")) { fin >> _N_SKIP_IMAGE; OK=1; }
        if (!strcmp(mot,"_COULEUR_FOND")) { fin >> _COULEUR_FOND; OK=1; }
        if (fin.rdstate()) OK=0;
        }
      }
  fin.close();
  if (!OK) return 2;
  if (_FILAIRE==0) 
      { 
      SHOWNODES=0; SHOWEDGES=1; SHOWSIDES=1; 
      NodeRadius=0.005*_LARGEUR; EdgeRadius=0.003*_LARGEUR;
      }
  if (_FILAIRE==1) 
      { 
      SHOWNODES=0; SHOWEDGES=1; SHOWSIDES=0; 
      NodeRadius=0.005*_LARGEUR; EdgeRadius=0.003*_LARGEUR;
      }
  rotation_phi_theta();
  ortho_perspective();
  Change_Couleur_Fond(_COULEUR_FOND);
  GLRenderScene();
  return 0;
  }

/*---------------------------------------------------------------------*/
