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

#ifndef EasyAnimMDLG_H
#define EasyAnimMDLG_H

#include <v/vmodald.h>

    class myCmdWindow;

    class myModalDialog : public vModalDialog
      {
      public:		//---------------------------------------- public
	myModalDialog(vBaseWindow* bw, char* title = "Change Color");
	virtual ~myModalDialog();		// Destructor
	virtual void DialogCommand(ItemVal,ItemVal,CmdType); // action selected
	virtual int myAction(char* msg);

      protected:	//--------------------------------------- protected

      private:		//--------------------------------------- private

	myCmdWindow* _myCmdWin;
      };
#endif

