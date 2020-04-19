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

#include <EasyDyn/vec.h>

#include <iostream>
#include <stdlib.h>
#include <string.h>

using namespace std;

#ifndef EASYDYNVISU_H
#define EASYDYNVISU_H

class shape
  {
  protected:
  //mth **mthref;
  int nbrmth,*nbrnodes, nbredges, nbrsurf;
  int *edgecolor, *surfcolor;
  //vec **nodecoord;
  int *edgenode1, *edgenode2;
  int *nbrnodesurf, **numnodesurf;
  int security;
  public:
  mth **mthref;
  vec **nodecoord;
  shape(){};
  shape(mth **mt_per_pt,int n_pt,vec **pt,int n_edges, int *enode1,int *enode2, int *eclr,
	int n_surf,int *n_nodesurf,int **num_nodesurf,int *sclr,int protection=1);
  shape(mth **mt_per_pt,int n_pt,vec **pt,int protection=1);
  shape(mth **mt_per_pt, istream &stream);
  shape(int n_mth,mth **mt,int *n_pt_per_mth,vec **pt,int n_edges, int *enode1,
	int *enode2, int *ecolor,int n_surf,int *n_nodesurf,int **numnodesurf,int *sclor,int protection=1);
  shape(int n_mth,mth **mt,int *n_pt_per_mth,vec **pt,int protection=1);
  shape(int n_mth, mth **mt, istream &stream);
  ~shape();
  shape *nextshape;
  int GetNbrNodes();
  int GetNbrEdges();
  int GetNbrSurf();
  virtual void WriteNodes(ostream &stream, int initnode);
  virtual void WriteNodes(mth *mthvisu,ostream &stream, int initnode, int norot=0);
  void WriteEdges(ostream &stream, int initnode, int initedge);
  void WriteSurf(ostream &stream, int initnode, int initsurf);
  virtual void WriteCoord(ostream &stream);
  virtual void WriteCoord(mth *mthvisu, ostream &stream, int norot=0);
  virtual const char* GetShapeType() { return "SHAPE"; };
  virtual void WriteShape(ostream &stream);
  };

class line: public shape
  {
  private:
  void build();
  protected:
  vec point1, point2;
  int ecolor;
  public:
  line(mth *mt, vec pt1, vec pt2, int eclr);
  line(mth *mt, istream &stream);
  ~line();
  const char* GetShapeType() { return "LINE"; };
  virtual void WriteShape(ostream &stream);
  };

class triangle: public shape
  {
  private:
  void build();
  protected:
  vec corner1, corner2, corner3;
  int ecolor, scolor;
  public:
  triangle(mth *mt, vec c1, vec c2, vec c3,
	   int eclr, int sclr);
  triangle(mth *mt, istream &stream);
  ~triangle();
  const char* GetShapeType() { return "TRIANGLE"; };
  virtual void WriteShape(ostream &stream);
  };

class triangle2: public shape
  {
  private:
  void build();
  protected:
  vec corner1, corner2, corner3;
  int ecolor, scolor;
  public:
  triangle2(mth *mt1, vec c1, mth* mt2, vec c2, mth* mt3, vec c3, 
            int eclr, int sclr);
  triangle2(mth *mt, istream &stream);
  ~triangle2();
  const char* GetShapeType() { return "TRIANGLE2"; };
  virtual void WriteShape(ostream &stream);
  };

class box: public shape
  {
  private:
  void build();
  protected:
  vec corner1, corner2;
  int ecolor, scolor;
  public:
  box(mth *mt, vec c1, vec c2, int eclr, int sclr);
  box(mth *mt, istream &stream);
  ~box();
  const char* GetShapeType() { return "BOX"; };
  virtual void WriteShape(ostream &stream);
  };

class grspline: public shape
  {
  private:
  void build();
  protected:
  int nbrtr,nbrseg;
  vec *pt, *tgt;
  int ecolor, scolor;
  public:
  grspline(mth *mt, vec *p, vec *t, int nt, int ns, int eclr);
  grspline(mth *mt, istream &stream);
  ~grspline();
  const char* GetShapeType() { return "GRSPLINE"; };
  virtual void WriteShape(ostream &stream);
  };

class frustum: public shape
  {
  private:
  void build();
  protected:
  vec axispt1, axispt2;
  double radius1, radius2;
  int nbrsect, ecolor, scolor;
  public:
  frustum(mth *mt, vec ap1, vec ap2,
	 double r1, double r2, int nbs,
	 int eclr, int sclr);
  frustum(mth *mt, istream &stream);
  ~frustum();
  const char* GetShapeType() { return "FRUSTUM"; };
  virtual void WriteShape(ostream &stream);
  };

class tyre: public shape
  {
  private:
  void build();
  protected:
  vec axispt1, axispt2;
  double radext, radint;
  int nbrsect, ecolor, scolor;
  public:
  tyre(mth *mt, vec ap1, vec ap2,
	 double ri, double re, int nbs,
	 int eclr, int sclr);
  tyre(mth *mt, istream &stream);
  ~tyre();
  const char* GetShapeType() { return "TYRE"; };
  virtual void WriteShape(ostream &stream);
  };

class gear: public shape
  {
  private:
  void build();
  protected:
  vec axispt1, axispt2;
  double module,holeradius, beta0;
  int nbrdents, ecolor, scolor;
  public:
  gear(mth *mt, vec ap1, vec ap2,
       double mod, int nd, double b0, double r,
       int eclr, int sclr);
  gear(mth *mt, istream &stream);
  ~gear();
  const char* GetShapeType() { return "GEAR"; };
  virtual void WriteShape(ostream &stream);
  };

class sphere: public shape
  {
  private:
  void build();
  protected:
  vec center;
  double radius;
  int nbrslices, nbrsect, ecolor, scolor;
  public:
  sphere(mth *mt, vec ce, double r, int nsl, int nbs, int eclr, int sclr);
  sphere(mth *mt, istream &stream);
  ~sphere();
  const char* GetShapeType() { return "SPHERE"; };
  virtual void WriteShape(ostream &stream);
  };

class habfile: public shape
  {
  private:
  void build();
  protected:
  char *filename;
  public:
  habfile(mth *mt, const char *filename);
  habfile(mth *mt, istream &stream);
  ~habfile();
  const char* GetShapeType() { return "HABFILE"; };
  virtual void WriteShape(ostream &stream);
  };

class var_surface: public shape
  {
  private:
  void build();
  protected:
  int nbrpt,scolor;
  public:
  var_surface(mth **mt_per_pt,int n_pt,vec **pt,int sclr,int protection);
  var_surface(mth **mt_per_pt, istream &stream);
  ~var_surface();
  const char* GetShapeType() { return "var_SURFACE"; };
  virtual void WriteShape(ostream &stream);
  };

class var_path: public shape
  {
  private:
  void build();
  protected:
  int nbrpt,ecolor;
  public:
  var_path(mth **mt_per_pt,int n_pt,vec **pt,int eclr,int protection);
  var_path(mth **mt_per_pt, istream &stream);
  ~var_path();
  const char* GetShapeType() { return "var_PATH"; };
  virtual void WriteShape(ostream &stream);
  };

class spring: public shape
  {
  private:
  void buildnodes();
  void build();
  protected:
  vec end1, end2;
  int nbrspirals, nbrsect, ecolor;
  double radius, lrod;
  public:
  spring(mth* mt1,vec en1, mth* mt2, vec en2, double r, double l,
         int nsp,int ns,int eclr);
  spring(mth *mt, istream &stream);
  ~spring();
  void WriteNodes(ostream &stream, int initnode);
  void WriteNodes(mth *mthvisu,ostream &stream, int initnode, int norot=0);
  void WriteCoord(ostream &stream);
  void WriteCoord(mth *mthvisu, ostream &stream, int norot=0);
  const char* GetShapeType() { return "SPRING"; };
  virtual void WriteShape(ostream &stream);
  };

class scene
  {
  protected:
  int nbrshapes,nbrnodes,nbredges,nbrsurf;
  shape *firstshape,*curshape;
  mth *mthvisu;
  int norot;
  public:
  scene();
  void AddShape(shape *newshape);
  void CreateVolFile(const char *filename);
  void WriteCoord(ostream &stream);
  void SetVisuFrame(mth *mthtmp, int nrt=0);
  };

#endif 
