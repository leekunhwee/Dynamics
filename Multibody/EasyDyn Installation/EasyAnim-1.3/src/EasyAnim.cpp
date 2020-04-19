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

#include <stdio.h>
#include <fstream>
#include "EasyAnim.h"              // Header file

myCmdWindow* _myCmdWin;
//=========================>>> myApp::myApp <<<==========================
  myApp::myApp(char* name, int sdi) : vApp(name, sdi)
  {
	  // Constructor
  }

//=====================>>> myApp::NewAppWin <<<==========================
  vWindow* myApp::NewAppWin(vWindow* win, char* name,
    int w, int h, vAppWinInfo* winInfo)
  {
	vWindow* _myCmdWin = win;
	vAppWinInfo* awinfo = winInfo;
    char *appname = name;

    if (!*name)
      {
	 appname = "EasyAnim 1.3 - Faculte Polytechnique de Mons (www.mecara.fpms.ac.be)";		// Default name
      }
	
    UserDebug1(Build,"myApp::NewAppWin(%s)\n",appname);

    // Create the first window using provided CmdWindow

	if (!_myCmdWin)
      {
	_myCmdWin = new myCmdWindow(appname, w, h);
      }
	
	if (!awinfo)
	awinfo = new vAppWinInfo(appname);

    return vApp::NewAppWin(_myCmdWin, appname, w, h, awinfo);
  }

//============================>>> myApp::Exit <<<===========================
  void myApp::Exit(void)
  {
    // This is called to close all windows.

    UserDebug(Build,"myApp::Exit()\n");

    vApp::Exit();		// Default behavior
  }

//======================>>> myApp::CloseAppWin <<<===========================
  int myApp::CloseAppWin(vWindow* win)
  {
    // This will be called BEFORE a window has been unregistered or
    // closed.  Default behavior: unregister and close the window.

    UserDebug(Build,"myApp::CloseAppWin()\n");

    return vApp::CloseAppWin(win);
  }

//=====================>>> myApp::AppCommand <<<==============================
  void myApp::AppCommand(vWindow* win, ItemVal id, ItemVal val, CmdType cType)
  {
    // Commands not processed by the window will be passed here

    UserDebug1(Build,"myApp::AppCmd(ID: %d)\n",id);
    vApp::AppCommand(win, id, val, cType);
  }

//=========================>>> myApp::KeyIn <<<==============================
  void myApp::KeyIn(vWindow* win, vKey key, unsigned int shift)
  {
    // Key strokes not processed by the window will be passed here

    vApp::KeyIn(win, key, shift);
  }

//###########################################################################

  static myApp my_App(" ",1);	// The instance of the app

//============================>>> AppMain <<<==============================
  int AppMain(int argc, char** argv)
  {
    // Use AppMain to create the main window

    myCmdWindow* cw = (myCmdWindow*) theApp->NewAppWin(0,"",700,550);

  if (argc == 2)
  {
  	char *fichier_tmp = argv[1];
	_myCmdWin->Lecture(fichier_tmp);
  }

   return 0;           // 0 means OK 
  }
