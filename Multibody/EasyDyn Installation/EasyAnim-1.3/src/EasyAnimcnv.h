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

#ifndef myCNV_H
#define myCNV_H

#include <fstream>
#include <v/vbglcnv.h>
#include "EasyAnimGL.h"

using namespace std;

class myOGLCanvasPane : public vBaseGLCanvasPane
  {
  public:		//---------------------------------------- public
    //  myOGLCanvasPane(unsigned int vGLmode = (vGL_RGB|vGL_DoubleBuffer|vGL_Depth),
  myOGLCanvasPane(unsigned int vGLmode = (vGL_RGB|vGL_DoubleBuffer|vGL_Depth),
                  PaneType pt = P_Canvas);
  virtual ~myOGLCanvasPane();

  virtual void graphicsInit(void);
  void TimerAnimate(void);	// for AuxTimer animation
	
  // Events
  virtual void MouseDown(int x, int y, int button);
  virtual void MouseUp(int x, int y, int button);
  virtual void MouseMove(int x, int y, int button);

  virtual void Redraw(int x, int y, int width, int height);
  virtual void Resize(int newW, int newH);

  protected:	//--------------------------------------- protected

  private:		//--------------------------------------- private
  int initDone;
  };

class def_volvan: public def_scene
  {
  ifstream fpos;
  char FileName[200];
  public :
  def_volvan(const char *nom_fichier_vol, int &errcode);
  ~def_volvan();
  void Next();
  void LoadImage(int NumImage);
  };

class def_volvmo: public def_scene
  {
  char FileName[200];
  VecStruct *PosRef,*DepReal,*DepImag;
  double DepMax;
  public :
  float alpha, beta, freq, amor;
  def_volvmo(const char *VolFileName, int &errcode);
  ~def_volvmo();
  void Next();
  int LoadMode(int num);
  void LoadImage(int NumImage);
  };

class def_uff: public def_scene
  {
  char FileName[200];
  VecStruct *PosRef,*DepReal,*DepImag;
  double DepMax;
  public :
  float alpha, beta, freq, amor;
  def_uff(const char *UffFileName1,const char *UffFileName2,
	     int &errcode);
  ~def_uff();
  void Next();
  int LoadMode(int num);
  void LoadImage(int NumImage);
  };

class def_repere: public def_scene
  {
  public :
  void Trace();
  def_repere();
  void LoadImage(int NumImage);
  };

extern def_repere* s1;
extern def_volvan* s2;
extern def_volvmo* s3;
extern def_uff* s4;

void repere();
void repere_plus();
void repere_moins();

#endif


