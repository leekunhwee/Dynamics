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

// Definition of the class relative to the scenes used for wiewing

#ifndef EasyAnimGLH
#define EasyAnimGLH

struct VecStruct
  {
    double x,y,z;
};

struct NodeStruct
  {
    double x,y,z;
    int color;
    char name[41];
};

struct EdgeStruct
  {
    int node1, node2;
    int color;
    char name[41];
};

struct SideStruct
  {
    int NbrNodes, *NumNode;
    int color;
    char name[41];
};

class def_scene
  {
  protected :
  int NbrNodes,NbrEdges,NbrSides,hidden,NbrImages,CurImage;
  NodeStruct *node; EdgeStruct *edge; SideStruct *side;
  float xmin,xmax,ymin,ymax,zmin,zmax;
  GLUquadricObj *NodeSphere;
  public:
  def_scene();
  ~def_scene();
  int Set(int npoint, float **point,int narr, int **arrete, 
           int nface, int **face);
  void GetDim(float &x1,float &x2,float &y1,float &y2,float &z1,float &z2);
  int Hidden();
  void Hide();
  void Show();
  void Homothetie(float MULT);	
  //void Get_Point(int num,float Pt[3]);
  int GetNumNodeOfName(const char *NodeName);
  int GetNumEdgeOfName(const char *EdgeName);
  int GetNumSideOfName(const char *SideName);
  int GetNbrImages();
  int GetCurImage();
  virtual void Next();
  virtual void LoadImage(int NumImage);
  virtual void Trace();
  void Trace2();
  };

extern GLfloat COULEUR[16][3];
extern GLfloat MATERIAU[16][3][3];

// Fonctions d'initialisation
void Init_Couleur();
void GLSetupRC(); 

// Fonctions permettant d'ajouter et de retirer des scenes
void Add_Scene(def_scene *sc);
void Remove_Scene();

// fonction de visualisation
void GLRenderScene();

// fonction de gestion de la fenetre
void GLResize(GLsizei w, GLsizei h);
void Change_Couleur_Fond(int c);
void SetOneView();
void SetMultiView();
void SwitchView();

//Fonctions pour gerer la dimension des noeuds et des lignes
void SetNodeRadius(double NewRadius);
void MultNodeRadius(double Factor);
void SetEdgeRadius(double NewRadius);
void MultEdgeRadius(double Factor);

// Fonctions pour activer ou desactiver la vision fil de fer
void SetShowNodesOn();
void SetShowNodesOff();
void SwitchShowNodes();
int  GetShowNodes();
void SetShowEdgesOn();
void SetShowEdgesOff();
void SwitchShowEdges();
int  GetShowEdges();
void SetShowSidesOn();
void SetShowSidesOff();
void SwitchShowSides();
int  GetShowSides();

// fonctions de gestion de la camera
void ortho_perspective();
void vue_depart();
void vue_centre();
void zoom_plus();
void zoom_moins();
void vue_selon_x_plus();
void vue_selon_x_moins();
void vue_selon_y_plus();
void vue_selon_y_moins();
void vue_selon_z_plus();
void vue_selon_z_moins();
void distance_plus();
void distance_moins();
void distance_plus_plus();
void distance_moins_moins();
void phi_plus();
void phi_moins();
void theta_plus();
void theta_moins();
void gamma_plus();
void gamma_moins();
void depl_x_plus();
void depl_x_moins();
void depl_y_plus();
void depl_y_moins();
void depl_z_plus();
void depl_z_moins();
void depl_haut();
void depl_bas();
void depl_gauche();
void depl_droite();
int Sauve_Config(const char *nom);
int Lecture_Config(const char *nom_fichier);

// fonctions de gestion de l'animation
void anime();
void start_anime();
void stop_anime();
void start_stop_anime();
void image_suivante();
void image_precedente();
void skip_image_plus();
void skip_image_moins();

#endif
