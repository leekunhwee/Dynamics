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

#include "EasyAnimmdlg.h"
#include "EasyAnimcmdw.h"
#include <v/vnotice.h>

int retval;
//@V@:BeginIDs
    enum {
	cbxid = 0,
    lblMainMsg = 1000,
    // add your id's here
  };
//@V@:EndIds

// Valeurs de la combo liste
char* comboList[] = 
	{
		"Black",
		"Blue",
		"Green",
		"Cyan",
		"Red",
		"Magenta",
		"Brown",
		"Light Gray",
		"Dark Gray",
		"Light Blue",
		"Light Green",
		"Light Cyan",
		"Light Red",
		"Light Magenta",
		"Yellow",
		"White",
		0
	};

//@V@:BeginDialogCmd DefaultCmds
    static DialogCmd DefaultCmds[] =
      {
	{C_Label, lblMainMsg, 0,"Change background color",NoList,CA_None,isSens,NoFrame, 0, 0},
	{C_Spinner, cbxid, cbxid,"",(void*)comboList,CA_Text,isSens,NoFrame,lblMainMsg,0,100},
	{C_Button, M_Cancel, 0, " Cancel ",NoList,CA_None,isSens,NoFrame,0, lblMainMsg},
	{C_Button, M_OK, 0, " OK ",NoList,CA_DefaultButton,isSens,NoFrame,M_Cancel,lblMainMsg},
	{C_EndOfList,0,0,0,0,CA_None,0,0,0}
    };
//@V@:EndDialogCmd


//======================>>> myModalDialog::myModalDialog <<<==================
  myModalDialog::myModalDialog(vBaseWindow* bw, char* title) :
    vModalDialog(bw, title)
  {
    UserDebug(Constructor,"myModalDialog::myModalDialog()\n")
    _myCmdWin = (myCmdWindow*) bw;
    AddDialogCmds(DefaultCmds);		// add the predefined commands
  }

//===================>>> myModalDialog::~myModalDialog <<<====================
  myModalDialog::~myModalDialog()
  {
    UserDebug(Destructor,"myModalDialog::~myModalDialog() destructor\n")
  }

//====================>>> myModalDialog::myAction <<<====================
  int myModalDialog::myAction(char* msg)
  {
    ItemVal ans,rval;
	ans = ShowModalDialog(msg,rval);
	retval = GetValue(cbxid);
    if (ans == M_Cancel)
	return 0;


	// *** Add code to process dialog values here

    return ans == M_OK;
  }

//====================>>> myModalDialog::DialogCommand <<<====================
  void myModalDialog::DialogCommand(ItemVal id, ItemVal retval, CmdType ctype)
  {
    UserDebug2(CmdEvents,"myModalDialog::DialogCommand(id:%d, val:%d)\n",id, retval)

      vModalDialog::DialogCommand(id,retval,ctype);
  }

