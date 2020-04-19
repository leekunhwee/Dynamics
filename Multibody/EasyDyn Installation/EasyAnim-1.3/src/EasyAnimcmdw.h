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

#ifndef EasyAnimCMDW_H
#define EasyAnimCMDW_H

#include <v/vcmdwin.h>	// So we can use vCmdWindow
#include <v/vmenu.h>	// For the menu pane
#include <v/vutil.h>	// For V Utilities
#include <v/vcmdpane.h> // command pane
#include <v/vstatusp.h> // status pane
#include <v/vtimer.h>	// Timer
#include <v/vicon.h>   // Icons
#include "zoompl.vbm"
#include "zoomm.vbm"
#include "xplus.vbm"
#include "xmoins.vbm"
#include "yplus.vbm"
#include "ymoins.vbm"
#include "zplus.vbm"
#include "zmoins.vbm"
#include "st_stop.vbm"
// #include "stop.vbm"
#include "suiv.vbm"
#include "prec.vbm"
// #include "fast.vbm"
// #include "slow.vbm"
// #include "save.vbm"
// #include "open.vbm"
// #include "print.vbm"
// #include "filaire.vbm"
#include "repere.vbm"
#include "reperep.vbm"
#include "reperem.vbm"
#include "noderadp.vbm"
#include "node.vbm"
#include "noderadm.vbm"
#include "line.vbm"
#include "lineradp.vbm"
#include "lineradm.vbm"
#include "side.vbm"
#include "phipl.vbm"
#include "phim.vbm"
#include "thetap.vbm"
#include "thetam.vbm"
#include "gammap.vbm"
#include "gammam.vbm"
//#include "distpl.vbm"
//#include "distm.vbm"
//#include "distpp.vbm"
//#include "distmm.vbm"
#include "dephaut.vbm"
#include "depbas.vbm"
#include "depg.vbm"
#include "depd.vbm"
#include "vuedep.vbm"
#include "multi.vbm"
#include "vcentr.vbm"

#ifdef vDEBUG
#include <v/vdebug.h>
#endif

#include "EasyAnimmdlg.h"     // my modal dialog
#include "EasyAnimcnv.h"      // myOGLCanvasPane

static vIcon ssicon(&st_stop_bits[0],st_stop_height,st_stop_width,st_stop_depth);
static vIcon precicon(&prec_bits[0],prec_height,prec_width,prec_depth);
static vIcon suivicon(&suiv_bits[0],suiv_height,suiv_width,suiv_depth);
//static vIcon fasticon(&fast_bits[0],fast_height,fast_width,fast_depth);
//static vIcon slowicon(&slow_bits[0],slow_height,slow_width,slow_depth);

static vIcon zpicon(&zoompl_bits[0],zoompl_height,zoompl_width,zoompl_depth);
static vIcon zmicon(&zoomm_bits[0],zoomm_height,zoomm_width,zoomm_depth);
//static vIcon filicon(&filaire_bits[0],filaire_height,filaire_width,filaire_depth);

static vIcon repicon(&repere_bits[0],repere_height,repere_width,repere_depth);
static vIcon reppicon(&reperep_bits[0],reperep_height,reperep_width,reperep_depth);
static vIcon repmicon(&reperem_bits[0],reperem_height,reperem_width,reperem_depth);
static vIcon nodeicon(&node_bits[0],node_height,node_width,node_depth);
static vIcon noderadpicon(&noderadp_bits[0],noderadp_height,noderadp_width,noderadp_depth);
static vIcon noderadmicon(&noderadm_bits[0],noderadm_height,noderadm_width,noderadm_depth);
static vIcon lineicon(&line_bits[0],line_height,line_width,line_depth);
static vIcon lineradpicon(&lineradp_bits[0],lineradp_height,lineradp_width,lineradp_depth);
static vIcon lineradmicon(&lineradm_bits[0],lineradm_height,lineradm_width,lineradm_depth);
static vIcon sideicon(&side_bits[0],side_height,side_width,side_depth);

static vIcon vuedepicon(&vuedep_bits[0],vuedep_height,vuedep_width,vuedep_depth);
static vIcon vuexpicon(&xplus_bits[0],xplus_height,xplus_width,xplus_depth);
static vIcon vuexmicon(&xmoins_bits[0],xmoins_height,xmoins_width,xmoins_depth);
static vIcon vueypicon(&yplus_bits[0],yplus_height,yplus_width,yplus_depth);
static vIcon vueymicon(&ymoins_bits[0],ymoins_height,ymoins_width,ymoins_depth);
static vIcon vuezpicon(&zplus_bits[0],zplus_height,zplus_width,zplus_depth);
static vIcon vuezmicon(&zmoins_bits[0],zmoins_height,zmoins_width,zmoins_depth);

static vIcon phiplusicon(&phipl_bits[0],phipl_height,phipl_width,phipl_depth);
static vIcon phimoinsicon(&phim_bits[0],phim_height,phim_width,phim_depth);
static vIcon thetaplusicon(&thetap_bits[0],thetap_height,thetap_width,thetap_depth);
static vIcon thetamoinsicon(&thetam_bits[0],thetam_height,thetam_width,thetam_depth);
static vIcon gammaplusicon(&gammap_bits[0],gammap_height,gammap_width,gammap_depth);
static vIcon gammamoinsicon(&gammam_bits[0],gammam_height,gammam_width,gammam_depth);

//static vIcon distplicon(&distpl_bits[0],distpl_height,distpl_width,distpl_depth);
//static vIcon distppicon(&distpp_bits[0],distpp_height,distpp_width,distpp_depth);
//static vIcon distmicon(&distm_bits[0],distm_height,distm_width,distm_depth);
//static vIcon distmmicon(&distmm_bits[0],distmm_height,distmm_width,distmm_depth);


static vIcon vuecentricon(&vcentr_bits[0],vcentr_height,vcentr_width,vcentr_depth);
static vIcon dephauticon(&dephaut_bits[0],dephaut_height,dephaut_width,dephaut_depth);
static vIcon depbasicon(&depbas_bits[0],depbas_height,depbas_width,depbas_depth);
static vIcon depgaucheicon(&depg_bits[0],depg_height,depg_width,depg_depth);
static vIcon depdroiteicon(&depd_bits[0],depd_height,depd_width,depd_depth);
static vIcon multiicon(&multi_bits[0],multi_height,multi_width,multi_depth);

	class myCmdWindow;
	
    class myAuxTimer : public vTimer
      {
      public:		//---------------------------------------- public
	myAuxTimer(myCmdWindow* cw) { cmdw = cw; }
	~myAuxTimer() {}
	virtual void TimerTick();
      private:		//--------------------------------------- private
	myCmdWindow* cmdw;
      };
	



    class myCmdWindow : public vCmdWindow
      {
	friend int AppMain(int, char**);	// allow AppMain access

      public:		//---------------------------------------- public
	myCmdWindow(char*, int, int);
	virtual ~myCmdWindow();
	virtual void WindowCommand(ItemVal id, ItemVal val, CmdType cType);
	void Lecture(char* fichier_config);
	myAuxTimer* GetmyAuxTimer() {return _auxTimer;}// Aux Timer
	  protected:	//--------------------------------------- protected

      private:		//--------------------------------------- private


	// Standard elements
	vMenuPane* myMenu;		// For the menu bar
	myOGLCanvasPane* myCanvas;		// For the canvas
        vCommandPane* DisplacePane;        // For the Command Pane
	vCommandPane* RotatePane;		// For the Command Pane
	vCommandPane* AnimPane;		// For the Command Pane
	vStatusPane* InfoPane;		// For the Command Pane
	myAuxTimer* _auxTimer;	// Aux Timer
 
	// Dialogs associated with CmdWindow
	myModalDialog* myMDlg;

      };
#endif
