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

#include <v/vnotice.h>	// for vNoticeDialog
#include <v/vkeys.h>	// to map keys
#include <v/vfilesel.h> // for file select
#include <v/vreply.h> // for file select

#include "EasyAnimcmdw.h"     // our header
#include "EasyAnimcnv.h"
#include <fstream>
#include <math.h>

// Déclarations et initialisations
extern int retval;

int _INIT_VOL=0;
int LOOP=0;

//	Start defines for the main window with 100

//@V@:BeginIDs
    enum {
	m_FirstCmd = 100, // Dummy Command
    blanc,
        m_open2,
        m_open3,
	m_vue_dep,
	zoom_in,
	zoom_out,
	m_start_stop,
	m_suiv,
	m_prec,
	m_fast,
	m_slow,
	m_shownodes,
	m_showedges,
	m_showsides,
	m_repere,
	m_reperep,
	m_reperem,
	m_noderadp,
	m_noderadm,
	m_lineradp,
	m_lineradm,
	m_vuexp,
	m_vuexm,
	m_vueyp,
	m_vueym,
	m_vuezp,
	m_vuezm,
	m_vuecentre,
	m_phiplus,
	m_phimoins,
	m_thetaplus,
	m_thetamoins,
	m_gammaplus,
	m_gammamoins,
	m_distpl,
	m_distpp,
	m_distm,
	m_distmm,
	m_dephaut,
	m_depbas,
	m_depgauche,
	m_depdroite,
	m_multi,
	m_couleur,
	m_loadmode,
	cmdAuxTimer,	// AuxTimer
	m_lblspeed,     // the label for the animation speed
	m_sldspeed,     // the slider for the animation speed
	m_lblprogress,  // the label for the animation progress bar
	m_sldprogress,  // the slider for the animation progress bar
        m_loop,
        m_info,
	blkLast		// Last item
      };
//@V@:EndIDs

//@V@:BeginPulldownMenu FileMenu
// Définition du Menu "File"
    static vMenu FileMenu[] =
      {
	{"&Open animation (vol-van)", M_Open, isSens, notChk, noKeyLbl, noKey, noSub},
	{"&Open modal analysis (vol-vmo)", m_open2, isSens, notChk, noKeyLbl, noKey, noSub},
	{"&Open modal analysis (universal)", m_open3, isSens, notChk, noKeyLbl, noKey, noSub},
	{"&Close...", M_CloseFile, isSens, notChk, noKeyLbl, noKey, noSub},
	{"-", M_Line, notSens, notChk, noKeyLbl, noKey, noSub},
	{"Save &Config", M_Save, isSens, notChk, noKeyLbl, noKey, noSub},
	{"-", M_Line, notSens, notChk, noKeyLbl, noKey, noSub},
	{"&Print", M_Print, notSens, notChk, noKeyLbl, noKey, noSub},
	{"-", M_Line, notSens, notChk, noKeyLbl, noKey, noSub},
	{"E&xit", M_Exit, isSens, notChk, noKeyLbl, noKey, noSub},
	{NULL} // dernière instruction des menus
      };
//@V@:EndPulldownMenu

//@V@:BeginPulldownMenu OptionMenu
// Définition du menu "Options"
	static vMenu OptionMenu[] =
	  {
	{"&Background Color...", m_couleur, isSens, notChk, noKeyLbl, noKey, noSub},
	{"&Load mode...", m_loadmode, isSens, notChk, noKeyLbl, noKey, noSub},
	//{"Hide/show &Nodes", m_shownodes, isSens, notChk, noKeyLbl, noKey, noSub},
        //{"Hide/show &Edges", m_showedges, isSens, notChk, noKeyLbl, noKey, noSub},
        //{"Hide/show &Sides", m_showsides, isSens, notChk, noKeyLbl, noKey, noSub},
	//{"Hide/show &Mark", m_repere, isSens, notChk, noKeyLbl, noKey, noSub},
	{NULL}
	  };
//@V@:EndPulldownMenu

//@V@:BeginMenu StandardMenu
// Définition de la barre de menu
    static vMenu StandardMenu[] =
      {
	{"&File", M_File, isSens, notUsed, notUsed, noKey, &FileMenu[0]},
	{"&Options", M_Options, isSens, notUsed, notUsed, noKey, &OptionMenu[0]},
	{NULL}
      };
//@V@:EndMenu


//@V@:BeginCmdPane ToolBar
// Définition de la barre d'outil
    static CommandObject DisplaceBar[] =
      {
	{C_IconButton,zoom_in,zoom_in,"",(void*)&zpicon,CA_None,isSens,NoFrame,0,0,0,"Zoom +"},
	{C_IconButton,zoom_out,zoom_out,"",(void*)&zmicon,CA_None,isSens,NoFrame,0,0,0,"Zoom -"},
	{C_Blank,blanc,0," ",NoList,CA_None,isSens,NoFrame,0,0,5},  
		{C_Blank,blanc,0," ",NoList,CA_None,isSens,NoFrame,0,0,5},
		{C_IconButton,m_vuecentre,m_vuecentre,"",(void*)&vuecentricon,CA_None,isSens,NoFrame,0,0,0,"Centered View"},
		{C_IconButton,m_dephaut,m_dephaut,"",(void*)&dephauticon,CA_None,isSens,NoFrame,0,0,0,"Up Move"},
		{C_IconButton,m_depbas,m_depbas,"",(void*)&depbasicon,CA_None,isSens,NoFrame,0,0,0,"Down move"},
		{C_IconButton,m_depgauche,m_depgauche,"",(void*)&depgaucheicon,CA_None,isSens,NoFrame,0,0,0,"Left Move"},
		{C_IconButton,m_depdroite,m_depdroite,"",(void*)&depdroiteicon,CA_None,isSens,NoFrame,0,0,0,"Right Move"},
	{C_ToggleIconButton,m_multi,0,"",(void*)&multiicon,CA_None,isSens,NoFrame,0,0,0,"MDI"},
	{C_Blank,blanc,0," ",NoList,CA_None,isSens,NoFrame,0,0,5},
	{C_ToggleIconButton,m_repere,m_repere,"",(void*)&repicon,CA_None,isSens,NoFrame,0,0,0,"Show marker"},
	{C_IconButton,m_reperep,m_reperep,"",(void*)&reppicon,CA_None,isSens,NoFrame,0,0,0,"Increase size of marker"},
	{C_IconButton,m_reperem,m_reperem,"",(void*)&repmicon,CA_None,isSens,NoFrame,0,0,0,"Decrease size of marker"},
	{C_Blank,blanc,0," ",NoList,CA_None,isSens,NoFrame,0,0,5},
	{C_ToggleIconButton,m_showedges,m_showedges,"",(void*)&lineicon,CA_None,isSens,NoFrame,0,0,0,"Show edges"},
	{C_IconButton,m_lineradp,m_lineradp,"",(void*)&lineradpicon,CA_None,isSens,NoFrame,0,0,0,"Increase radius of edges"},
	{C_IconButton,m_lineradm,m_lineradm,"",(void*)&lineradmicon,CA_None,isSens,NoFrame,0,0,0,"Decrease radius of edges"},
	{C_Blank,blanc,0," ",NoList,CA_None,isSens,NoFrame,0,0,5},
	{C_ToggleIconButton,m_shownodes,m_shownodes,"",(void*)&nodeicon,CA_None,isSens,NoFrame,0,0,0,"Show nodes"},
	{C_IconButton,m_noderadp,m_noderadp,"",(void*)&noderadpicon,CA_None,isSens,NoFrame,0,0,0,"Increase size of nodes"},
	{C_IconButton,m_noderadm,m_noderadm,"",(void*)&noderadmicon,CA_None,isSens,NoFrame,0,0,0,"Decrease size of nodes"},
	{C_Blank,blanc,0," ",NoList,CA_None,isSens,NoFrame,0,0,5},
	{C_ToggleIconButton,m_showsides,m_showsides,"",(void*)&sideicon,CA_None,isSens,NoFrame,0,0,0,"Show sides"},
	//{C_IconButton,m_distmm,m_distmm,"",(void*)&distmmicon,CA_None,isSens,NoFrame,0,0,0,"Decrease the distance of the camera"},
	//{C_IconButton,m_distm,m_distm,"",(void*)&distmicon,CA_None,isSens,NoFrame,0,0,0,"Decrease the distance of the camera"},
	//{C_IconButton,m_distpl,m_distpl,"",(void*)&distplicon,CA_None,isSens,NoFrame,0,0,0,"Increase the distance of the camera"},
	//{C_IconButton,m_distpp,m_distpp,"",(void*)&distppicon,CA_None,isSens,NoFrame,0,0,0,"Increase the distance of the camera"},
	{C_EndOfList,0,0,0,0,CA_None,0,0,0} // dernière instruction des listes
      };

//@V@:EndCmdPane
//@V@:BeginCmdPane CommandBar
	static CommandObject RotateBar[] =
	{ 
		{C_IconButton,m_vue_dep,m_vue_dep,"",(void*)&vuedepicon,CA_None,isSens,NoFrame,0,0,0,"Start View"},
		{C_Blank,blanc,0," ",NoList,CA_None,isSens,NoFrame,0,0,5},
		{C_IconButton,m_vuexp,m_vuexp,"",(void*)&vuexpicon,CA_None,isSens,NoFrame,0,0,0,"X + View"},
		{C_IconButton,m_vuexm,m_vuexm,"",(void*)&vuexmicon,CA_None,isSens,NoFrame,0,0,0,"X - View"},
		{C_IconButton,m_vueyp,m_vueyp,"",(void*)&vueypicon,CA_None,isSens,NoFrame,0,0,0,"Y + View"},
		{C_IconButton,m_vueym,m_vueym,"",(void*)&vueymicon,CA_None,isSens,NoFrame,0,0,0,"Y - View"},
		{C_IconButton,m_vuezp,m_vuezp,"",(void*)&vuezpicon,CA_None,isSens,NoFrame,0,0,0,"Z + View"},
		{C_IconButton,m_vuezm,m_vuezm,"",(void*)&vuezmicon,CA_None,isSens,NoFrame,0,0,0,"Z - View"},
		{C_Blank,blanc,0," ",NoList,CA_None,isSens,NoFrame,0,0,5},
		{C_IconButton,m_phiplus,m_phiplus,"",(void*)&phiplusicon,CA_None,isSens,NoFrame,0,0,0,"Phi + Rotation"},
		{C_IconButton,m_phimoins,m_phimoins,"",(void*)&phimoinsicon,CA_None,isSens,NoFrame,0,0,0,"Phi - Rotation"},
		{C_IconButton,m_thetaplus,m_thetaplus,"",(void*)&thetaplusicon,CA_None,isSens,NoFrame,0,0,0,"Theta + Rotation"},
		{C_IconButton,m_thetamoins,m_thetamoins,"",(void*)&thetamoinsicon,CA_None,isSens,NoFrame,0,0,0,"Theta - Rotation"},
		{C_IconButton,m_gammaplus,m_gammaplus,"",(void*)&gammaplusicon,CA_None,isSens,NoFrame,0,0,0,"Gamma + Rotation"},
		{C_IconButton,m_gammamoins,m_gammamoins,"",(void*)&gammamoinsicon,CA_None,isSens,NoFrame,0,0,0,"Gamma - Rotation"},
		{C_EndOfList,0,0,0,0,CA_None,0,0,0}
	};
//@V@:EndCmdPane

// Définition de la barre d'outil
    static CommandObject AnimBar[] =
      {
	//	{C_IconButton,m_slow,m_slow,"",(void*)&slowicon,CA_None,isSens,NoFrame,0,0,0,"Slow animation"},
	//	{C_IconButton,m_fast,m_fast,"",(void*)&fasticon,CA_None,isSens,NoFrame,0,0,0,"Fast animation"},
	{C_Blank,blanc,0," ",NoList,CA_None,isSens,NoFrame,0,0,5},   
	{C_Label,m_lblspeed,0,"Speed",NoList,CA_None,isSens,NoFrame,0,0,0,""},
	{C_Slider,m_sldspeed,0,"",NoList,CA_Horizontal,isSens,NoFrame,0,0,0,""},
	{C_Blank,blanc,0," ",NoList,CA_None,isSens,NoFrame,0,0,5},   
	{C_Label,m_lblprogress,0,"Progress",NoList,CA_None,isSens,NoFrame,0,0,0,""},
       	{C_IconButton,m_prec,m_prec,"",(void*)&precicon,CA_None,isSens,NoFrame,0,0,0,"Previous Frame"}, 
	{C_Slider,m_sldprogress,0,"",NoList,CA_Horizontal,isSens,NoFrame,0,0,0,""},
	{C_IconButton,m_suiv,m_suiv,"",(void*)&suivicon,CA_None,isSens,NoFrame,0,0,0,"Next Frame"},
	{C_Blank,blanc,0," ",NoList,CA_None,isSens,NoFrame,0,0,5},   
	{C_ToggleIconButton,m_start_stop,0,"",(void*)&ssicon,CA_None,isSens,NoFrame,0,0,0,"Start/Stop"},
	{C_Blank,blanc,0," ",NoList,CA_None,isSens,NoFrame,0,0,5},   
	{C_ToggleButton,m_loop,0,"LOOP",NoList,CA_None,isSens,NoFrame,0,0},   
	{C_EndOfList,0,0,0,0,CA_None,0,0,0} // derniere instruction des listes
      };

//@V@:BeginCmdPane InfoBar
static vStatus InfoBar[] =
 { 
   {"No application yet",m_info,CA_None,isSens,500},
   {0,0,0,0,0}
 };

//====================>>> myAuxTimer::TimerTick <<<====================
  void myAuxTimer::TimerTick()
  {
	cmdw->WindowCommand(cmdAuxTimer, cmdAuxTimer, C_Button); // update clock
  }


//====================>>> myCmdWindow::myCmdWindow <<<====================
// Définition du contenu de la fenêtre de commandes
  myCmdWindow::myCmdWindow(char* name, int width, int height) :
    vCmdWindow(name, width, height)
  {
    UserDebug1(Constructor,"myCmdWindow::myCmdWindow(%s) Constructor\n",name)
    
	// The Menu Bar
    myMenu = new vMenuPane(StandardMenu);
    AddPane(myMenu);

	// The Canvas
    myCanvas = new myOGLCanvasPane;
    AddPane(myCanvas);

	// The Timer
	_auxTimer = new myAuxTimer(this);	// create aux timer
    _auxTimer->TimerSet(1);		// 1/2 second intervals

    // The Command Pane
    DisplacePane = new vCommandPane(DisplaceBar);
    AddPane(DisplacePane);
     
    // The second Command Pane
    RotatePane = new vCommandPane(RotateBar);
    AddPane(RotatePane);
     
    // The second Command Pane
    AnimPane = new vCommandPane(AnimBar);
    AddPane(AnimPane);
     
    // The status Pane
    InfoPane = new vStatusPane(InfoBar);
    AddPane(InfoPane);
	
    // Associated dialogs
    myMDlg = new myModalDialog(this,name);

    // Show Window
    ShowWindow();
    WindowCommand(cmdAuxTimer,cmdAuxTimer,C_Button);	// update clock
    
    SetValue(m_sldspeed,100,Value);
    SetValue(m_sldprogress,0,Value);
    SetValue(m_multi,0,Value);
    SetValue(m_repere,1,Value);
    SetValue(m_shownodes,0,Value);
    SetValue(m_showsides,1,Value);
    SetValue(m_showedges,1,Value);
    }

//====================>>> myCmdWindow::~myCmdWindow <<<====================
// Destructeur de la fenêtre de commandes
	myCmdWindow::~myCmdWindow()
  {
    UserDebug(Destructor,"myCmdWindow::~myCmdWindow() destructor\n")

    // Now put a delete for each new in the constructor.

    delete myMenu;
    delete myCanvas;
    delete DisplacePane;
    delete RotatePane;
    delete AnimPane;
    delete myMDlg;
    _auxTimer->TimerStop();	// end it
    delete _auxTimer;	// free it


  }

//===================>>> myCmdWindow::Lecture <<<======================
  void myCmdWindow::Lecture(char* fichier_config)
  {
	  Lecture_Config(fichier_config);
  }

  
//===================>>> myCmdWindow::OpenFile <<<======================
// ouverture des fichier ".vol" et ".van"

  void myCmdWindow::WindowCommand(ItemVal id, ItemVal val, CmdType cType)
  {
  // Default: route menu and toolbar commands here

  static char nomfichiervol[200]="";
  vNoticeDialog note(this);
  char Utiltmp[40];
	
  UserDebug1(CmdEvents,"myCmdWindow:WindowCommand(%d)\n",id)

  // instructions exécutées lors de l'activation des menu, boutons, ...
  switch (id)
    {
    case M_Open:// This demos vFileSelect dialog
      {
      char name[200] = "";        // start out with null name
      vFileSelect fsel(this);     // an instance of vFileSelect
      int fI = 0;		      // Filter index
      static char* filter[]={"*.vol",0};    // Filter for file select
      int ans;  
      // Show the file select dialog
      ans = fsel.FileSelect("Open file",name,99,filter,fI);
      if (ans && *name)   // User picked a file name
	  {
          if (s2!=NULL) { Remove_Scene(); delete s2; s2=NULL; }
          if (s3!=NULL) { Remove_Scene(); delete s3; s3=NULL; }
          if (s4!=NULL) { Remove_Scene(); delete s4; s4=NULL; }
          strcpy(nomfichiervol,name);
          int ec;
          s2=new def_volvan(nomfichiervol,ec);
          if (s2==NULL) note.Notice("Impossible d'allouer s2");
          if ((s2!=NULL) && (ec<4))
            {
            Add_Scene(s2);
            float xmin,xmax,ymin,ymax,zmin,zmax,diag;
            // Remise a l'echelle du repere
            s1->GetDim(xmin,xmax,ymin,ymax,zmin,zmax);
            s1->Homothetie(1/xmax);
            s2->GetDim(xmin,xmax,ymin,ymax,zmin,zmax);
            diag=sqrt((xmax-xmin)*(xmax-xmin)+(ymax-ymin)*(ymax-ymin)
                    +(zmax-zmin)*(zmax-zmin));
            s1->Homothetie(0.2*diag);
            if ((ec==1)||(ec==3)) 
                {
                vue_centre(); 
                SetNodeRadius(0.005*diag); 
                SetEdgeRadius(0.003*diag);
		}
            GLRenderScene();
            SetValue(m_sldprogress,0,Value);
            SetValue(m_shownodes,GetShowNodes(),Value);
            SetValue(m_showedges,GetShowEdges(),Value);
            SetValue(m_showsides,GetShowSides(),Value);
            char Message[80];
            strcpy(Message,"Anim - Number of images=");
            sprintf(Utiltmp,"%4d",s2->GetNbrImages());
            strcat(Message,Utiltmp);
            SetString(m_info,Message);
            }
          else { delete s2; s2=NULL; }
          if ((ec==1)||(ec==3)) 
            note.Notice ("Configuration (cfg) file not found or not conform");
          if ((ec==2)||(ec==3)) 
              note.Notice ("Animation (van) file not found or corrupted");
          if ((ec==4)||(ec==6)) 
              note.Notice ("Image (vol) file not found or corrupted");
          if (ec==5) 
              note.Notice ("Memory allocation problems");
          }
      break;
      }

    case m_open2:// Opening a modal file
      {
      char name[200] = "";        // start out with null name
      vFileSelect fsel(this);     // an instance of vFileSelect
      int fI = 0;		      // Filter index
      static char* filter[]={"*.vol",0};    // Filter for file select
      int ans;  
      // Show the file select dialog
      ans = fsel.FileSelect("Open file",name,99,filter,fI);
      if (ans && *name)   // User picked a file name
	  {
          if (s2!=NULL) { Remove_Scene(); delete s2; s2=NULL; }
          if (s3!=NULL) { Remove_Scene(); delete s3; s3=NULL; }
          if (s4!=NULL) { Remove_Scene(); delete s4; s4=NULL; }
          strcpy(nomfichiervol,name);
          int ec;
          s3=new def_volvmo(nomfichiervol,ec);
          if (s3==NULL) note.Notice("Impossible d'allouer s3");
          if ((s3!=NULL) && (ec<4))
            {
            Add_Scene(s3);
            float xmin,xmax,ymin,ymax,zmin,zmax,diag;
            // Remise a l'echelle du repere
            s1->GetDim(xmin,xmax,ymin,ymax,zmin,zmax);
            s1->Homothetie(1/xmax);
            s3->GetDim(xmin,xmax,ymin,ymax,zmin,zmax);
            diag=sqrt((xmax-xmin)*(xmax-xmin)+(ymax-ymin)*(ymax-ymin)
                    +(zmax-zmin)*(zmax-zmin));
            s1->Homothetie(0.2*diag);
            if ((ec==1)||(ec==3)) 
                {
                vue_centre(); 
                SetNodeRadius(0.005*diag); 
                SetEdgeRadius(0.003*diag);
		}
            GLRenderScene();
            SetValue(m_sldprogress,0,Value);
            SetValue(m_shownodes,GetShowNodes(),Value);
            SetValue(m_showedges,GetShowEdges(),Value);
            SetValue(m_showsides,GetShowSides(),Value);
            if (ec<2)
              {
              char Message[80];
              strcpy(Message,"Modal - Freq(Hz)=");
              sprintf(Utiltmp,"%8.4f",s3->freq);
              strcat(Message,Utiltmp);
              strcat(Message," Damp(%)=");
              sprintf(Utiltmp,"%8.4f",100*s3->amor);
              strcat(Message,Utiltmp);
              SetString(m_info,Message);
              }
            else SetString(m_info,"Modal - No mode");
            }
          else { delete s3; s3=NULL; }
          if ((ec==1)||(ec==3)) 
            note.Notice ("Configuration (cfg) file not found or not conform");
          if ((ec==2)||(ec==3)) 
              note.Notice ("Mode (vmo) file not found or corrupted");
          if ((ec==4)||(ec==6)) 
              note.Notice ("Image (vol) file not found or corrupted");
          if (ec==5) 
              note.Notice ("Memory allocation problems");
          }
      break;
      }

    case m_open3:// Opening a modal file
      {
      char name[200] = "";        // start out with null name
      vFileSelect fsel(this);     // an instance of vFileSelect
      int fI = 0;		      // Filter index
      static char* filter[]={"*.unv;*.uff","*",0};    // Filter for file select
      int ans;  
      // Show the file select dialog
      ans = fsel.FileSelect("Open geometry file",name,99,filter,fI);
      if (ans && *name)   // User picked a file name
	  {
          if (s2!=NULL) { Remove_Scene(); delete s2; s2=NULL; }
          if (s3!=NULL) { Remove_Scene(); delete s3; s3=NULL; }
          if (s4!=NULL) { Remove_Scene(); delete s4; s4=NULL; }
          strcpy(nomfichiervol,name);
          note.Notice ("Now choose the uff modal file");
          char name2[200] = "";
          ans = fsel.FileSelect("Open modal file",name,99,filter,fI);
          if (ans && *name) { strcpy(name2,name); }
          int ec;
          s4=new def_uff(nomfichiervol,name2,ec);
          if (s4==NULL) note.Notice("Impossible d'allouer s4");
          if ((s4!=NULL) && (ec<2))
            {
            Add_Scene(s4);
            float xmin,xmax,ymin,ymax,zmin,zmax,diag;
            // Remise a l'echelle du repere
            s1->GetDim(xmin,xmax,ymin,ymax,zmin,zmax);
            s1->Homothetie(1/xmax);
            s4->GetDim(xmin,xmax,ymin,ymax,zmin,zmax);
            diag=sqrt((xmax-xmin)*(xmax-xmin)+(ymax-ymin)*(ymax-ymin)
                    +(zmax-zmin)*(zmax-zmin));
            s1->Homothetie(0.2*diag);
            vue_centre(); 
            SetNodeRadius(0.005*diag); 
            SetEdgeRadius(0.003*diag);
            GLRenderScene();
            SetValue(m_sldprogress,0,Value);
            SetValue(m_shownodes,GetShowNodes(),Value);
            SetValue(m_showedges,GetShowEdges(),Value);
            SetValue(m_showsides,GetShowSides(),Value);
            if (ec!=1)
              {
              char Message[80];
              strcpy(Message,"Modal - Freq(Hz)=");
              sprintf(Utiltmp,"%8.4f",s4->freq);
              strcat(Message,Utiltmp);
              strcat(Message," Damp(%)=");
              sprintf(Utiltmp,"%8.4f",100*s4->amor);
              strcat(Message,Utiltmp);
              SetString(m_info,Message);
              }
            else SetString(m_info,"Modal - No mode");
            }
          else { delete s4; s4=NULL; }
          if (ec==1) 
            note.Notice ("Cannot open modal file");
          if (ec==2) 
              note.Notice ("Error in modal file");
          if (ec==3) 
              note.Notice ("Uff block 82 corrupted");
          if (ec==4) 
              note.Notice ("Uff block 82 not found");
          if (ec==7) 
              note.Notice ("Uff block 15 corrupted");
          if (ec==8) 
              note.Notice ("Uff block 15 not found");
          if (ec==5) 
              note.Notice ("Memory allocation problems");
          if (ec==6) 
            note.Notice ("Cannot open geometry file");
          }
      break;
      }


	//@V@:Case M_Save
      case M_Save:
	  {
	  char nom_fichier_cfg[200]="";
	  if(strcmp(nomfichiervol,""))
	     {
             int ll=strlen(nomfichiervol);
	     strcpy(nom_fichier_cfg,nomfichiervol);
             nom_fichier_cfg[ll-3]=0;
	     strcat(nom_fichier_cfg,"cfg");				
             if(Sauve_Config(nom_fichier_cfg)) 
               note.Notice("Was not able to save configuration file");
	     }
	  else note.Notice("You must have opened an animation");
	  break;
	  }	//@V@:EndCase

	//@V@:Case M_CloseFile
	case M_CloseFile:
	  {
	  Change_Couleur_Fond(0);
          if (s2!=NULL) { Remove_Scene(); delete s2; s2=NULL; }
          if (s3!=NULL) { Remove_Scene(); delete s3; s3=NULL; }
          if (s4!=NULL) { Remove_Scene(); delete s4; s4=NULL; }
          strcpy(nomfichiervol,"");
          float xmin,xmax,ymin,ymax,zmin,zmax;
          s1->GetDim(xmin,xmax,ymin,ymax,zmin,zmax);
          s1->Homothetie(1/xmax);
	  vue_depart();
          vue_centre();
          break;
	  }	//@V@:EndCase


	//@V@:Case M_Exit
	case M_Exit:
	  {
            if (s1!=NULL) { delete s1; s1=NULL; }
            if (s2!=NULL) { delete s2; s2=NULL; }
            if (s3!=NULL) { delete s3; s3=NULL; }
            if (s4!=NULL) { delete s4; s4=NULL; }
	    theApp->Exit();
	    break;
	  }	//@V@:EndCase

        case m_loadmode:// Chargement d'un autre mode
          {
	  vReplyDialog rp(this);
          if (rp.Reply("Enter mode number",Utiltmp,9)==M_OK)
	    {
	    int modenum=atoi(Utiltmp),lstatus;
            if (s3!=NULL) lstatus=s3->LoadMode(modenum);
            if (s4!=NULL) lstatus=s4->LoadMode(modenum);
            if (lstatus==1) note.Notice("Invalid mode number");
            if (lstatus==2) note.Notice("Problem in vmo file");
            if (!lstatus)
              {
              GLRenderScene();
              char Message[80];
              strcpy(Message,"Modal - Freq(Hz)=");
              if (s3!=NULL) sprintf(Utiltmp,"%8.4f",s3->freq);
              if (s4!=NULL) sprintf(Utiltmp,"%8.4f",s4->freq);
              strcat(Message,Utiltmp);
              strcat(Message," Damp(%)=");
              if (s3!=NULL) sprintf(Utiltmp,"%8.4f",100*s3->amor);
              if (s4!=NULL) sprintf(Utiltmp,"%8.4f",100*s4->amor);
              strcat(Message,Utiltmp);
              SetString(m_info,Message);
              //sprintf(Utiltmp,"%8.4f",s2->alpha);
              //sprintf(Utiltmp,"%8.4f",s2->beta);
              }
            else SetString(m_info,"Modal - No mode");
            }
          break;
          }

	//@V@:Case m_couleur
	case m_couleur:
	  {
		  int retour;
		 retour = myMDlg->myAction("Sample Modal Dialog");
		 if (retour) Change_Couleur_Fond(retval);
		 break;
	  } //@V@:EndCase

	//@V@:Case Vue_depart 
	case m_vue_dep:
		{
		vue_depart();
		if (!GetValue(m_start_stop)) GLRenderScene();
		break;
		}

    //@V@:Case zoom in
		case zoom_in:
		{
		zoom_moins();
		if (!GetValue(m_start_stop)) GLRenderScene();
		break;
		}

	//@V@:Case zoom out
		case zoom_out:
		{
		zoom_plus();
		if (!GetValue(m_start_stop)) GLRenderScene();
		break;
		}

	//@V@:Case shownodes
		case m_shownodes:
		{
		if (GetValue(m_shownodes)) SetShowNodesOn(); 
                else SetShowNodesOff();
		if (!GetValue(m_start_stop)) GLRenderScene();
		break;
		}

	//@V@:Case showedges
		case m_showedges:
		{
		if (GetValue(m_showedges)) SetShowEdgesOn(); 
                else SetShowEdgesOff();
		if (!GetValue(m_start_stop)) GLRenderScene();
		break;
		}

	//@V@:Case showsides
		case m_showsides:
		{
		if (GetValue(m_showsides)) SetShowSidesOn(); 
                else SetShowSidesOff();
		if (!GetValue(m_start_stop)) GLRenderScene();
		break;
		}

	//@V@:Case m_suiv
		case m_suiv:
		{
		image_suivante();
		SetValue(m_start_stop,0,Value);	
		break;
		}

	//@V@:Case m_prec
		case m_prec:
		{
		stop_anime();
                image_precedente();
	        SetValue(m_start_stop,0,Value);
		break;
		}

	//@V@:Case m_fast
		case m_fast:
		{
			skip_image_plus();
			break;
		}

	//@V@:Case m_slow
		case m_slow:
		{
			skip_image_moins();
			break;
		}

	//@V@:Case start_stop
		case m_start_stop:
		{
		if (GetValue(m_start_stop)) start_anime();
                else stop_anime();
		break;
		}

	//@V@:Case sldprogress
		case m_sldprogress:
		{
		stop_anime();
		if (s2!=NULL) if (s2->GetNbrImages()>2)
                    s2->LoadImage((int) val*(s2->GetNbrImages()-1)/100+1);
		if (s3!=NULL) if (s3->GetNbrImages()>2)
                    s3->LoadImage((int) val*(s3->GetNbrImages()-1)/100+1);
		if (s4!=NULL) if (s4->GetNbrImages()>2)
                    s4->LoadImage((int) val*(s4->GetNbrImages()-1)/100+1);
		GLRenderScene();
	        SetValue(m_start_stop,0,Value);
		break;
		}

	//@V@:Case repere
		case m_repere:
		{
		if (GetValue(m_repere)) s1->Show(); else s1->Hide();
		if (!GetValue(m_start_stop)) GLRenderScene();
		break;
		}

	//@V@:Case repere plus
		case m_reperep:
		{
		repere_plus();
		if (!GetValue(m_start_stop)) GLRenderScene();
		break;
		}

	//@V@:Case repere moins
		case m_reperem:
		{
		repere_moins();
		if (!GetValue(m_start_stop)) GLRenderScene();
		break;
		}
	    
	//@V@:Case node radius smaller
		case m_noderadm:
		{
		MultNodeRadius(0.8);
		if (!GetValue(m_start_stop)) GLRenderScene();
		break;
		}

	//@V@:Case node radius larger
		case m_noderadp:
		{
		MultNodeRadius(1.25);
		if (!GetValue(m_start_stop)) GLRenderScene();
		break;
		}
	    
	//@V@:Case line radius smaller
		case m_lineradm:
		{
		MultEdgeRadius(0.8);
		if (!GetValue(m_start_stop)) GLRenderScene();
		break;
		}

	//@V@:Case line radius larger
		case m_lineradp:
		{
		MultEdgeRadius(1.25);
		if (!GetValue(m_start_stop)) GLRenderScene();
		break;
		}

	//@V@:Case vue selon x positif
		case m_vuexp:
		{
		vue_selon_x_plus();
		if (!GetValue(m_start_stop)) GLRenderScene();
		break;
		}

	//@V@:Case vue selon x négatif
		case m_vuexm:
		{
		vue_selon_x_moins();
		if (!GetValue(m_start_stop)) GLRenderScene();
		break;
		}
	
	//@V@:Case vue selon y positif
		case m_vueyp:
		{
		vue_selon_y_plus();
		if (!GetValue(m_start_stop)) GLRenderScene();
		break;
		}

	//@V@:Case vue selon y négatif
		case m_vueym:
		{
		vue_selon_y_moins();
		if (!GetValue(m_start_stop)) GLRenderScene();
		break;
		}

	//@V@:Case vue selon z positif
		case m_vuezp:
		{
		vue_selon_z_plus();
		if (!GetValue(m_start_stop)) GLRenderScene();
		break;
		}

	//@V@:Case vue selon z négatif
		case m_vuezm:
		{
		vue_selon_z_moins();
		if (!GetValue(m_start_stop)) GLRenderScene();
		break;
		}

	//@V@:Case phi_plus
		case m_phiplus:
		{
		phi_plus();
		if (!GetValue(m_start_stop)) GLRenderScene();
		break;
		}

	//@V@:Case phi_moins
		case m_phimoins:
		{
		phi_moins();
		if (!GetValue(m_start_stop)) GLRenderScene();
		break;
		}

	//@V@:Case theta_plus
		case m_thetaplus:
		{
		theta_plus();
		if (!GetValue(m_start_stop)) GLRenderScene();
		break;
		}

	//@V@:Case theta_moins
		case m_thetamoins:
		{
		theta_moins();
		if (!GetValue(m_start_stop)) GLRenderScene();
		break;
		}

	//@V@:Case gamma_plus
		case m_gammaplus:
		{
		gamma_plus();
		if (!GetValue(m_start_stop)) GLRenderScene();
		break;
		}

	//@V@:Case gamma_moins
		case m_gammamoins:
		{
		gamma_moins();
		if (!GetValue(m_start_stop)) GLRenderScene();
		break;
		}

	//@V@:Case distance plus
		//case m_distpl:
		//{
		//distance_plus();
		//GLRenderScene();
		//break;
		//}

	//@V@:Case distance plus plus
		//case m_distpp:
		//{
		//distance_plus_plus();
		//GLRenderScene();
		//break;
		//}

	//@V@:Case distance moins
		//case m_distm:
		//{
		//distance_moins();
		//GLRenderScene();
		//break;
		//}

	//@V@:Case distance moins moins
		//case m_distmm:
		//{
		//distance_moins_moins();
		//GLRenderScene();
		//break;
		//}

	//@V@:Case depl x positif
		case m_vuecentre:
		{
		vue_centre();
		if (!GetValue(m_start_stop)) GLRenderScene();
		break;
		}

//@V@:Case depl vers haut
		case m_dephaut:
		{
		depl_haut();
		if (!GetValue(m_start_stop)) GLRenderScene();
		break;
		}
//@V@:Case depl vers bas
		case m_depbas:
		{
		depl_bas();
		if (!GetValue(m_start_stop)) GLRenderScene();
		break;
		}
//@V@:Case depl vers gauche
		case m_depgauche:
		{
		depl_gauche();
		if (!GetValue(m_start_stop)) GLRenderScene();
		break;
		}
//@V@:Case depl vers droite
		case m_depdroite:
		{
		depl_droite();
		if (!GetValue(m_start_stop)) GLRenderScene();
		break;
		}

//@V@:Case multi vues
		case m_multi:
		{
		if (GetValue(m_multi)) SetMultiView(); else SetOneView();
		if (!GetValue(m_start_stop)) GLRenderScene();
		break;
		}

//@V@:Case multi vues
		case m_sldspeed:
		{
		_auxTimer->TimerStop();
		_auxTimer->TimerSet(201-2*val);
		break;
		}

//@V@:Case multi vues
		case m_loop:
		{		
		if (GetValue(m_loop)) LOOP=1; else LOOP=0;
		break;
		}

//@V@:Case auxTimer
		case cmdAuxTimer:	// Event from aux timer
	  {
	  myCanvas->TimerAnimate();
          if (s2!=NULL)
	      {
	      if (s2->GetNbrImages()>1) SetValue(m_sldprogress,
                    (int) 100*(s2->GetCurImage()-1)/(s2->GetNbrImages()-1),Value);
	      if ((!LOOP)&&(s2->GetCurImage()==s2->GetNbrImages()))
                 {
	         stop_anime();
	         SetValue(m_start_stop,0,Value);
	         }
              }
          if (s3!=NULL)
	      {
	      if (s3->GetNbrImages()>1) SetValue(m_sldprogress,
                    (int) 100*(s3->GetCurImage()-1)/(s3->GetNbrImages()-1),Value);
	      if ((!LOOP)&&(s3->GetCurImage()==s3->GetNbrImages()))
                 {
	         stop_anime();
	         SetValue(m_start_stop,0,Value);
	         }
              }
          if (s4!=NULL)
	      {
	      if (s4->GetNbrImages()>1) SetValue(m_sldprogress,
                    (int) 100*(s4->GetCurImage()-1)/(s4->GetNbrImages()-1),Value);
	      if ((!LOOP)&&(s4->GetCurImage()==s4->GetNbrImages()))
                 {
	         stop_anime();
	         SetValue(m_start_stop,0,Value);
	         }
              }
	  break;
	  }	//@V@:EndCase 
	 
//@V@:Case default
// instructions exécutées par défaut
	  	default:		// route unhandled commands up
	  {
	    vCmdWindow::WindowCommand(id, val, cType);
		break;
	  }
      }
  }
