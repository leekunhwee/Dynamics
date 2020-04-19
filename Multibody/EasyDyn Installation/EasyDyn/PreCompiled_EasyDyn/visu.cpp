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

#include <EasyDyn/visu.h>
#include <EasyDyn/spline.h>

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <fstream>
using namespace std;

mth mthvisudefault;

/*************************************************************************/

int shape::GetNbrNodes() 
{
	int nbrtotnodes=0;
 	for (int imth=0;imth<nbrmth;imth++)
	{
		nbrtotnodes+=nbrnodes[imth];
	}
	return nbrtotnodes; 
}

/*************************************************************************/

int shape::GetNbrEdges() { return nbredges; }

/*************************************************************************/

int shape::GetNbrSurf() { return nbrsurf; }

/*************************************************************************/

shape::shape(mth **mt_per_pt,int n_pt,vec **pt,int n_edges, int *enode1,int *enode2,
	int *eclr,int n_surf,int *n_nodesurf,int **num_nodesurf,int *sclr,int protection)
{
  if (protection==1)
  {
	security=1;
	mthref=mt_per_pt;
	nbrmth=n_pt;
	nbrnodes=new int[nbrmth];
	for (int imth=0;imth<nbrmth;imth++) nbrnodes[imth]=1;	
	nodecoord=new vec *[nbrmth];
	for (int imth=0; imth<nbrmth; imth++) 
	{
		nodecoord[imth]=new vec[nbrnodes[imth]];
		for (int imthnode=0;imthnode<nbrnodes[imth];imthnode++)
		{
			nodecoord[imth][imthnode]=pt[imth][imthnode];
		}
	}
	nbredges=n_edges;
	edgenode1=new int[nbredges];
	edgenode2=new int[nbredges];
	edgecolor=new int[nbredges];
	for (int iedge=0; iedge<nbrsurf; iedge++)
	{
		edgenode1[iedge]=enode1[iedge];
		edgenode2[iedge]=enode2[iedge];
		edgecolor[iedge]=eclr[iedge];
	}
	
	nbrsurf=n_surf;
	nbrnodesurf=new int[nbrsurf];
	numnodesurf=new int*[nbrsurf];
	surfcolor=new int[nbrsurf];
	for (int isurf=0; isurf<nbrsurf; isurf++)
	{
		nbrnodesurf[isurf]=n_nodesurf[isurf];
		numnodesurf[isurf]=new int[nbrnodesurf[isurf]];
		for (int ipt=0; ipt<nbrnodesurf[isurf]; ipt++)
		{
			numnodesurf[isurf][ipt]=num_nodesurf[isurf][ipt];
		}
		surfcolor[isurf]=sclr[isurf];
	}
  }
  else
  {
	security=0;
	mthref=mt_per_pt;
	nbrmth=n_pt;
	nbrnodes=new int[nbrmth];
	for (int imth=0;imth<nbrmth;imth++) nbrnodes[imth]=1;
	nodecoord=pt;
	nbredges=n_edges;
	edgenode1=enode1;
	edgenode2=enode2;
	edgecolor=eclr;
	nbrsurf=n_surf;
	nbrnodesurf=n_nodesurf;
	numnodesurf=num_nodesurf;
	surfcolor=sclr;
  }
}

/*************************************************************************/

shape::shape(mth **mt_per_pt,int n_pt,vec **pt,int protection)
{
  if (protection==1)
  {
	security=1;
	mthref=mt_per_pt;
	nbrmth=n_pt;
	nbrnodes=new int[nbrmth];
	for (int imth=0;imth<nbrmth;imth++) nbrnodes[imth]=1;	
	nodecoord=new vec *[nbrmth];
	for (int imth=0; imth<nbrmth; imth++) 
	{
		nodecoord[imth]=new vec[nbrnodes[imth]];
		for (int imthnode=0;imthnode<nbrnodes[imth];imthnode++)
		{
			nodecoord[imth][imthnode]=pt[imth][imthnode];
		}
	}
  }
  else
  {
	security=0;
	mthref=mt_per_pt;
	nbrmth=n_pt;
	nbrnodes=new int[nbrmth];
	for (int imth=0;imth<nbrmth;imth++) nbrnodes[imth]=1;
	nodecoord=pt;
  }
}

/*************************************************************************/

shape::shape(mth **mt_per_pt, istream &stream)
{
	security=1;
	int n_pt;
	int eclr;
	mthref=mt_per_pt;
	stream >> n_pt >> eclr;
	nbrmth=n_pt;
	nbrnodes=new int[nbrmth];
	for (int imth=0;imth<nbrmth;imth++) nbrnodes[imth]=1;	
	nodecoord=new vec *[nbrmth];
	for (int imth=0; imth<nbrmth; imth++)
	{
		nodecoord[imth]=new vec[nbrnodes[imth]];
		for (int imthnode=0;imthnode<nbrnodes[imth];imthnode++)
		{
			stream >> nodecoord[imth][imthnode].x >> nodecoord[imth][imthnode].y >> nodecoord[imth][imthnode].z;
		}
	}
}

/*************************************************************************/

shape::shape(int n_mth,mth **mt,int *n_pt_per_mth,vec **pt,int n_edges,
	int *enode1,int *enode2, int *eclr,int n_surf,int *n_nodesurf,int **num_nodesurf,int *sclr,int protection)
{
  if (protection==1)
  {
	security=1;
	nbrmth=n_mth;
	mthref=mt;
	nbrnodes=new int[nbrmth];
	for (int imth=0;imth<nbrmth;imth++) nbrnodes[imth]=n_pt_per_mth[imth];	
	nodecoord=new vec *[nbrmth];
	for (int imth=0; imth<nbrmth; imth++) 
	{
		nodecoord[imth]=new vec[nbrnodes[imth]];
		for (int imthnode=0;imthnode<nbrnodes[imth];imthnode++)
		{
			nodecoord[imth][imthnode]=pt[imth][imthnode];
		}
	}
	nbredges=n_edges;
	edgenode1=new int[nbredges];
	edgenode2=new int[nbredges];
	edgecolor=new int[nbredges];
	for (int iedge=0; iedge<nbrsurf; iedge++)
	{
		edgenode1[iedge]=enode1[iedge];
		edgenode2[iedge]=enode2[iedge];
		edgecolor[iedge]=eclr[iedge];
	}
	
	nbrsurf=n_surf;
	nbrnodesurf=new int[nbrsurf];
	numnodesurf=new int*[nbrsurf];
	surfcolor=new int[nbrsurf];
	for (int isurf=0; isurf<nbrsurf; isurf++)
	{
		nbrnodesurf[isurf]=n_nodesurf[isurf];
		numnodesurf[isurf]=new int[nbrnodesurf[isurf]];
		for (int ipt=0; ipt<nbrnodesurf[isurf]; ipt++)
		{
			numnodesurf[isurf][ipt]=num_nodesurf[isurf][ipt];
		}
		surfcolor[isurf]=sclr[isurf];
	}
  }
  else
  {
	security=0;
	mthref=mt;
	nbrmth=n_mth;
	nbrnodes=n_pt_per_mth;
	nodecoord=pt;
	nbredges=n_edges;
	edgenode1=enode1;
	edgenode2=enode2;
	edgecolor=eclr;
	nbrsurf=n_surf;
	nbrnodesurf=n_nodesurf;
	numnodesurf=num_nodesurf;
	surfcolor=sclr;
  }
}

/*************************************************************************/

shape::shape(int n_mth,mth **mt,int *n_pt_per_mth,vec **pt,int protection)
{
  if (protection==1)
  {
	security=1;
	mthref=mt;
	nbrmth=n_mth;
	nbrnodes=new int[nbrmth];
	for (int imth=0;imth<nbrmth;imth++) nbrnodes[imth]=n_pt_per_mth[imth];	
	nodecoord=new vec *[nbrmth];
	for (int imth=0; imth<nbrmth; imth++) 
	{
		nodecoord[imth]=new vec[nbrnodes[imth]];
		for (int imthnode=0;imthnode<nbrnodes[imth];imthnode++)
		{
			nodecoord[imth][imthnode]=pt[imth][imthnode];
		}
	}
  }
  else
  {
	security=0;
	mthref=mt;
	nbrmth=n_mth;
	nbrnodes=n_pt_per_mth;
	nodecoord=pt;
  }
}

/*************************************************************************/

shape::shape(int n_mth, mth **mt, istream &stream)
{
	security=1;
	int eclr;
	stream >> nbrmth >> eclr;
	mthref=mt;
	nbrnodes=new int[nbrmth];
	nodecoord=new vec *[nbrmth];
	for (int imth=0; imth<nbrmth; imth++) 
	{
		stream >> nbrnodes[imth];
		nodecoord[imth]=new vec[nbrnodes[imth]];
		for (int imthnode=0;imthnode<nbrnodes[imth];imthnode++)
		{
			stream >> nodecoord[imth][imthnode].x >> nodecoord[imth][imthnode].y >> nodecoord[imth][imthnode].z;
		}
	}
}

/*************************************************************************/

shape::~shape()
{
	if (security==1)
	{
		delete nbrnodes;
		for (int imth=0; imth<nbrmth; imth++) delete nodecoord[imth];
		delete nodecoord;
		delete edgenode1;
		delete edgenode2;
		delete edgecolor;
		delete nbrnodesurf;
		for (int isurf=0;isurf<nbrsurf;isurf++) delete numnodesurf[isurf];
		delete numnodesurf;
		delete surfcolor;
	}
	else
	{
		// !!!!!!! But the user has the responsability of deallocating all the pointers : nbrnodes,
		// nodecoord[imth],nodecoord,edgenode1,edgenode2,edgecolor,nbrnodesurf,numnodesurf[isurf],
		// numnodesurf and surfcolor
	}
}

/*************************************************************************/

void shape::WriteNodes(ostream &stream, int initnode)
  {
  int inode=0;
  vec globalcoord;
  for (int imth=0; imth<nbrmth; imth++)
    {
	for (int imthnode=0; imthnode<nbrnodes[imth]; imthnode++)
	{
	globalcoord=(*mthref[imth])*nodecoord[imth][imthnode];
	stream << inode+initnode << " " << globalcoord.x
		<< " " << globalcoord.y << " " << globalcoord.z << "\n";
	inode++;
	}
    }
  inode--;
  }

/*************************************************************************/

void shape::WriteNodes(mth *mthvisu, ostream &stream, int initnode,
                       int norot)
  {
  int inode=0;
  vec globalcoord;
  mth mthtmp;
  if (norot) { mthtmp.e.x=-mthvisu[0].e.x; mthtmp.e.y=-mthvisu[0].e.y;
               mthtmp.e.z=-mthvisu[0].e.z;}
  else mthtmp=mthvisu[0].inv();
  for (int imth=0; imth<nbrmth; imth++)
    {
	for (int imthnode=0; imthnode<nbrnodes[imth]; imthnode++)
	{
	globalcoord=mthtmp*(*mthref[imth])*nodecoord[imth][imthnode];
	stream << inode+initnode << " " << globalcoord.x
		<< " " << globalcoord.y << " " << globalcoord.z << "\n";
	inode++;
	}
    }
  inode--;
  }

/*************************************************************************/

void shape::WriteEdges(ostream &stream, int initnode, int initedge)
  {
  int iedge;
  for (iedge=0; iedge<nbredges; iedge++)
    stream << iedge+initedge << " " << edgenode1[iedge]+initnode << " "
	   << edgenode2[iedge]+initnode << " " << edgecolor[iedge] << "\n";
  }

/*************************************************************************/

void shape::WriteSurf(ostream &stream, int initnode, int initsurf)
  {
  int isurf, inode;
  for (isurf=0; isurf<nbrsurf; isurf++)
    {
    stream << isurf+initsurf << " " << nbrnodesurf[isurf];
    for (inode=0;inode<nbrnodesurf[isurf];inode++)
      stream << " " << numnodesurf[isurf][inode]+initnode;
    stream << " " << surfcolor[isurf] << "\n";
    }
  }

/*************************************************************************/

void shape::WriteCoord(ostream &stream)
  {
  vec globalcoord;
  for (int imth=0; imth<nbrmth; imth++)
    {
	for (int imthnode=0; imthnode<nbrnodes[imth]; imthnode++)
	{
	globalcoord=(*mthref[imth])*nodecoord[imth][imthnode];
	stream << globalcoord << "\n";
	}
    }
  }

/*************************************************************************/

void shape::WriteCoord(mth *mthvisu, ostream &stream, int norot)
  {
  vec globalcoord;
  mth mthtmp;
  if (norot) { mthtmp.e.x=-mthvisu[0].e.x; mthtmp.e.y=-mthvisu[0].e.y;
               mthtmp.e.z=-mthvisu[0].e.z;}
  else mthtmp=mthvisu[0].inv();
  for (int imth=0; imth<nbrmth; imth++)
    {
	for (int imthnode=0; imthnode<nbrnodes[imth]; imthnode++)
	{
	globalcoord=mthtmp*(*mthref[imth])*nodecoord[imth][imthnode];
	stream << globalcoord << "\n";
	}
    }
  }

/*************************************************************************/

void shape::WriteShape(ostream &stream)
{
	stream << nbrmth << endl << endl;
	for (int imth=0;imth<nbrmth;imth++)
	{
		stream << mthref[imth];
		stream << nbrnodes[imth] << " : ";
		for (int imthnode=0;imthnode<nbrnodes[imth];imthnode++)
		{
			stream << nodecoord[imth][imthnode] << " ";
		}
		stream << endl;
	}
	stream << endl;
	stream << nbredges << endl;
	for (int iedge=0;iedge<nbredges;iedge++)
	{
		stream << iedge << " " << edgenode1[iedge] << "," << edgenode2[iedge] << " " << edgecolor[iedge] << endl;
	}
	stream << endl;
	stream << nbrsurf << endl;
	for (int isurf=0;isurf<nbrsurf;isurf++)
	{
		stream << isurf << " " << nbrnodesurf[isurf] << " : ";
		for (int ipt=0;ipt<nbrnodesurf[isurf];ipt++)
		{
			stream << numnodesurf[isurf][ipt] << " ";
		}
		stream << endl;
	}
	stream << endl;
}

/*************************************************************************/

line::line(mth *mt, vec pt1, vec pt2, int eclr)
  {
  mthref=new mth *[1];
  mthref[0]=mt;
  point1=pt1;
  point2=pt2;
  ecolor=eclr;
  build();
  }

/*************************************************************************/

line::line(mth *mt, istream &stream)
  {
  mthref=new mth *[1];
  mthref[0]=mt;
  stream >> point1 >> point2 >> ecolor;
  build();
  }

/*************************************************************************/

void line::build()
  {
  nbrmth=1;   nbrnodes=new int[nbrmth];   nbrnodes[0]=2;   nbredges=1;   nbrsurf=0;

  nodecoord=new vec *[nbrmth];
  for (int imth=0; imth<nbrmth; imth++) nodecoord[imth]=new vec[nbrnodes[imth]];
  edgenode1=new int[nbredges];
  edgenode2=new int[nbredges];
  edgecolor=new int[nbredges];
  for (int iedge=0; iedge<nbredges; iedge++) edgecolor[iedge]=ecolor;

  nodecoord[0][0]=point1;
  nodecoord[0][1]=point2;

  edgenode1[0]=0; edgenode2[0]=1;
  }

/*************************************************************************/

line::~line()
  {
  delete nbrnodes;
  for (int imth=0; imth<nbrmth; imth++) delete nodecoord[imth];
  delete nodecoord;
  delete edgenode1;
  delete edgenode2;
  delete edgecolor;
  }

/*************************************************************************/

void line::WriteShape(ostream &stream)
    {
    stream << point1 << point2 << ecolor << "\n";
    };

/*************************************************************************/

triangle::triangle(mth *mt,
		  vec c1, vec c2, vec c3,
		  int eclr, int sclr)
  {
  mthref=new mth *[1];
  mthref[0]=mt;
  corner1=c1;
  corner2=c2;
  corner3=c3;
  ecolor=eclr;
  scolor=sclr;
  build();
  }

/*************************************************************************/

triangle::triangle(mth *mt, istream &stream)
  {
  mthref=new mth *[1];
  mthref[0]=mt;
  stream >> corner1 >> corner2 >> corner3 >> ecolor >> scolor;
  build();
  }

/*************************************************************************/

void triangle::build()
  {
  nbrmth=1; nbrnodes=new int[nbrmth]; nbrnodes[0]=3; nbredges=3; nbrsurf=1;

  nodecoord=new vec *[nbrmth];
  for (int imth=0; imth<nbrmth; imth++) nodecoord[imth]=new vec[nbrnodes[imth]];
  edgenode1=new int[nbredges];
  edgenode2=new int[nbredges];
  edgecolor=new int[nbredges];
  for (int iedge=0; iedge<nbredges; iedge++) edgecolor[iedge]=ecolor;
  nbrnodesurf=new int[nbrsurf];
  numnodesurf=new int*[nbrsurf];
  surfcolor=new int[nbrsurf];
  nbrnodesurf[0]=3;
  numnodesurf[0]=new int[3];
  surfcolor[0]=scolor;

  nodecoord[0][0]=corner1;
  nodecoord[0][1]=corner2;
  nodecoord[0][2]=corner3;

  edgenode1[0]=0; edgenode2[0]=1;
  edgenode1[1]=1; edgenode2[1]=2;
  edgenode1[2]=2; edgenode2[2]=0;

  numnodesurf[0][0]=0;
  numnodesurf[0][1]=1;
  numnodesurf[0][2]=2;
  }

/*************************************************************************/

triangle::~triangle()
  {
  delete nbrnodes;
  for (int imth=0; imth<nbrmth; imth++) delete nodecoord[imth];
  delete nodecoord;
  delete edgenode1;
  delete edgenode2;
  delete edgecolor;
  delete nbrnodesurf;
  for (int isurf=0;isurf<nbrsurf;isurf++) delete numnodesurf[isurf];
  delete numnodesurf;
  delete surfcolor;
  }

/*************************************************************************/

void triangle::WriteShape(ostream &stream)
  {
  stream << corner1 << corner2 << corner3;
  stream << ecolor << " " << scolor << "\n";
  }

/*************************************************************************/

triangle2::triangle2(mth *mt1, vec c1, mth* mt2, vec c2, mth* mt3, vec c3,
                     int eclr, int sclr)
  {
  mthref=new mth *[3];
  mthref[0]=mt1;
  corner1=c1;
  mthref[1]=mt2;
  corner2=c2;
  mthref[2]=mt3;
  corner3=c3;
  ecolor=eclr;
  scolor=sclr;
  build();
  }

/*************************************************************************/

triangle2::triangle2(mth *mt, istream &stream)
  {
    mthref=new mth *[3];
    mthref[0]=mt;
    mthref[1]=mt;
    mthref[2]=mt;
    stream >> corner1 >> corner2 >> corner3 >> ecolor >> scolor;
    build();
  }

/*************************************************************************/

void triangle2::build()
  {
  nbrmth=3; nbrnodes=new int[nbrmth]; 
  nbrnodes[0]=1; nbrnodes[1]=1; nbrnodes[2]=1;
  nbredges=3; nbrsurf=1;

  nodecoord=new vec *[nbrmth];
  for (int imth=0; imth<nbrmth; imth++) nodecoord[imth]=new vec[nbrnodes[imth]];
  edgenode1=new int[nbredges];
  edgenode2=new int[nbredges];
  edgecolor=new int[nbredges];
  for (int iedge=0; iedge<nbredges; iedge++) edgecolor[iedge]=ecolor;
  nbrnodesurf=new int[nbrsurf];
  numnodesurf=new int*[nbrsurf];
  surfcolor=new int[nbrsurf];
  nbrnodesurf[0]=3;
  numnodesurf[0]=new int[3];
  surfcolor[0]=scolor;

  nodecoord[0][0]=corner1;
  nodecoord[1][0]=corner2;
  nodecoord[2][0]=corner3;

  edgenode1[0]=0; edgenode2[0]=1;
  edgenode1[1]=1; edgenode2[1]=2;
  edgenode1[2]=2; edgenode2[2]=0;

  numnodesurf[0][0]=0;
  numnodesurf[0][1]=1;
  numnodesurf[0][2]=2;
  }

/*************************************************************************/

triangle2::~triangle2()
  {
  delete nbrnodes;
  for (int imth=0; imth<nbrmth; imth++) delete nodecoord[imth];
  delete nodecoord;
  delete edgenode1;
  delete edgenode2;
  delete edgecolor;
  delete nbrnodesurf;
  for (int isurf=0;isurf<nbrsurf;isurf++) delete numnodesurf[isurf];
  delete numnodesurf;
  delete surfcolor;
  }

/*************************************************************************/

void triangle2::WriteShape(ostream &stream)
  {
  stream << corner1 << corner2 << corner3;
  stream << ecolor << " " << scolor << "\n";
  }

/*************************************************************************/

box::box(mth *mt, vec c1, vec c2, int eclr, int sclr)
  {
  mthref=new mth *[1];
  mthref[0]=mt;
  corner1=c1;
  corner2=c2;
  ecolor=eclr;
  scolor=sclr;
  build();
  }

/*************************************************************************/

box::box(mth *mt, istream &stream)
  {
  mthref=new mth *[1];
  mthref[0]=mt;
  stream >> corner1 >> corner2 >> ecolor >> scolor;
  build();
  }

/*************************************************************************/

void box::build()
  {
  nbrmth=1; nbrnodes=new int[nbrmth]; nbrnodes[0]=8; nbredges=12; nbrsurf=6;

  nodecoord=new vec *[nbrmth];
  for (int imth=0; imth<nbrmth; imth++) nodecoord[imth]=new vec[nbrnodes[imth]];
  edgenode1=new int[nbredges];
  edgenode2=new int[nbredges];
  edgecolor=new int[nbredges];
  for (int iedge=0; iedge<nbredges; iedge++) edgecolor[iedge]=ecolor;
  nbrnodesurf=new int[nbrsurf];
  numnodesurf=new int*[nbrsurf];
  surfcolor=new int[nbrsurf];
  for (int isurf=0; isurf<nbrsurf; isurf++)
    {
    nbrnodesurf[isurf]=4;
    numnodesurf[isurf]=new int[4];
    surfcolor[isurf]=scolor;
    }

  nodecoord[0][0]=corner1;
  nodecoord[0][1].x=corner1.x; nodecoord[0][1].y=corner2.y; nodecoord[0][1].z=corner1.z;
  nodecoord[0][2].x=corner2.x; nodecoord[0][2].y=corner2.y; nodecoord[0][2].z=corner1.z;
  nodecoord[0][3].x=corner2.x; nodecoord[0][3].y=corner1.y; nodecoord[0][3].z=corner1.z;
  nodecoord[0][4].x=corner1.x; nodecoord[0][4].y=corner1.y; nodecoord[0][4].z=corner2.z;
  nodecoord[0][5].x=corner1.x; nodecoord[0][5].y=corner2.y; nodecoord[0][5].z=corner2.z;
  nodecoord[0][6]=corner2;
  nodecoord[0][7].x=corner2.x; nodecoord[0][7].y=corner1.y; nodecoord[0][7].z=corner2.z;

  edgenode1[0]=0; edgenode2[0]=1;
  edgenode1[1]=1; edgenode2[1]=2;
  edgenode1[2]=2; edgenode2[2]=3;
  edgenode1[3]=3; edgenode2[3]=0;
  edgenode1[4]=4; edgenode2[4]=5;
  edgenode1[5]=5; edgenode2[5]=6;
  edgenode1[6]=6; edgenode2[6]=7;
  edgenode1[7]=7; edgenode2[7]=4;
  edgenode1[8]=0; edgenode2[8]=4;
  edgenode1[9]=1; edgenode2[9]=5;
  edgenode1[10]=2; edgenode2[10]=6;
  edgenode1[11]=3; edgenode2[11]=7;

  numnodesurf[0][0]=0; numnodesurf[0][1]=1;
  numnodesurf[0][2]=2; numnodesurf[0][3]=3;
  numnodesurf[1][0]=4; numnodesurf[1][1]=5;
  numnodesurf[1][2]=6; numnodesurf[1][3]=7;
  numnodesurf[2][0]=0; numnodesurf[2][1]=1;
  numnodesurf[2][2]=5; numnodesurf[2][3]=4;
  numnodesurf[3][0]=1; numnodesurf[3][1]=2;
  numnodesurf[3][2]=6; numnodesurf[3][3]=5;
  numnodesurf[4][0]=2; numnodesurf[4][1]=3;
  numnodesurf[4][2]=7; numnodesurf[4][3]=6;
  numnodesurf[5][0]=3; numnodesurf[5][1]=0;
  numnodesurf[5][2]=4; numnodesurf[5][3]=7;
  }

/*************************************************************************/

box::~box()
  {
  delete nbrnodes;
  for (int imth=0; imth<nbrmth; imth++) delete nodecoord[imth];
  delete nodecoord;
  delete edgenode1;
  delete edgenode2;
  delete edgecolor;
  delete nbrnodesurf;
  for (int isurf=0;isurf<nbrsurf;isurf++) delete numnodesurf[isurf];
  delete numnodesurf;
  delete surfcolor;
  }

/*************************************************************************/

void box::WriteShape(ostream &stream)
  {
  stream << corner1 << corner2 << ecolor << " " << scolor << "\n";
  }

/*************************************************************************/

grspline::grspline(mth *mt, vec *p, vec *t, int nt, int ns,
		   int eclr)
  {
  mthref=new mth *[1];
  mthref[0]=mt;
  nbrtr=nt;
  nbrseg=ns;
  pt=new vec[nbrtr+1];
  tgt=new vec[nbrtr+1];
  for (int i=0; i<=nbrtr; i++) { pt[i]=p[i]; tgt[i]=t[i]; }
  ecolor=eclr;
  build();
  }

/*************************************************************************/

grspline::grspline(mth *mt, istream &stream)
  {
  mthref=new mth *[1];
  mthref[0]=mt;
  stream >> nbrtr >> nbrseg;
  pt=new vec[nbrtr+1];
  tgt=new vec[nbrtr+1];
  for (int i=0; i<=nbrtr; i++) { stream >> pt[i] >> tgt[i]; }
  stream >> ecolor;
  build();
  }

/*************************************************************************/

void grspline::build()
  {
  nbrmth=1; nbrnodes=new int[nbrmth]; nbrnodes[0]=nbrtr*nbrseg+1; nbredges=nbrtr*nbrseg; nbrsurf=0;

  nodecoord=new vec *[nbrmth];
  for (int imth=0; imth<nbrmth; imth++) nodecoord[imth]=new vec[nbrnodes[imth]];
  edgenode1=new int[nbredges];
  edgenode2=new int[nbredges];
  edgecolor=new int[nbredges];
  int iedge;
  for (iedge=0; iedge<nbredges; iedge++)
    {
    edgecolor[iedge]=ecolor;
    edgenode1[iedge]=iedge;
    edgenode2[iedge]=iedge+1;
    }

  int ipt,iseg;
  for (ipt=0; ipt<nbrtr; ipt++) for (iseg=0; iseg<nbrseg; iseg++)
    {
    double u=(1.0*iseg)/nbrseg;
    nodecoord[0][nbrseg*ipt+iseg]=pt[ipt]*N1(1-u)+pt[ipt+1]*N1(u)
			      -tgt[ipt]*N2(1-u)+tgt[ipt+1]*N2(u);
    }
  nodecoord[0][nbrnodes[0]-1]=pt[nbrtr];
  }

/*************************************************************************/

grspline::~grspline()
  {
  delete nbrnodes;
  for (int imth=0; imth<nbrmth; imth++) delete nodecoord[imth];
  delete nodecoord;
  delete edgenode1;
  delete edgenode2;
  delete edgecolor;
  }

/*************************************************************************/

void grspline::WriteShape(ostream &stream)
  {
  stream << nbrtr << " " << nbrseg << "\n";
  for (int i=0; i<=nbrtr; i++) { stream << pt[i] << tgt[i]; }
  stream << ecolor << "\n";
  }

/*************************************************************************/

frustum::frustum(mth *mt, vec ap1, vec ap2,
	 double r1, double r2, int nbs,
	 int eclr, int sclr)
  {
  mthref=new mth *[1];
  mthref[0]=mt;
  axispt1=ap1;
  axispt2=ap2;
  radius1=r1;
  radius2=r2;
  nbrsect=nbs;
  ecolor=eclr;
  scolor=sclr;
  build();
  }

/*************************************************************************/

frustum::frustum(mth *mt, istream &stream)
  {
  mthref=new mth *[1];
  mthref[0]=mt;
  stream >> axispt1 >> axispt2;
  stream >> radius1 >> radius2 >> nbrsect >> ecolor >> scolor;
  build();
  }

/*************************************************************************/

void frustum::build()
  {
  nbrmth=1; nbrnodes=new int[nbrmth]; nbrnodes[0]=2*nbrsect; nbredges=3*nbrsect; nbrsurf=nbrsect+2;

  nodecoord=new vec *[nbrmth];
  for (int imth=0; imth<nbrmth; imth++) nodecoord[imth]=new vec[nbrnodes[imth]];
  edgenode1=new int[nbredges];
  edgenode2=new int[nbredges];
  edgecolor=new int[nbredges];
  for (int iedge=0; iedge<nbredges; iedge++) edgecolor[iedge]=ecolor;
  nbrnodesurf=new int[nbrsurf];
  numnodesurf=new int*[nbrsurf];
  surfcolor=new int[nbrsurf];
  int isurf;
  for (isurf=0; isurf<nbrsurf; isurf++) surfcolor[isurf]=scolor;
  for (isurf=0; isurf<nbrsect; isurf++)
    {
    nbrnodesurf[isurf]=4;
    numnodesurf[isurf]=new int[4];
    }
  nbrnodesurf[nbrsect]=nbrsect; numnodesurf[nbrsect]=new int[nbrsect];
  nbrnodesurf[nbrsect+1]=nbrsect; numnodesurf[nbrsect+1]=new int[nbrsect];

  vec xaxis(1,0,0), yaxis(0,1,0), zaxis=axispt2-axispt1;
  if (zaxis.length()<1E-6) zaxis.put(0.0,0.0,1.0);
  if (fabs(zaxis.x) > fabs(zaxis.y))
    {
    xaxis=yaxis^zaxis; xaxis.unite();
    yaxis=zaxis^xaxis; yaxis.unite();
    }
  else
    {
    yaxis=zaxis^xaxis; yaxis.unite();
    xaxis=yaxis^zaxis; xaxis.unite();
    }

  for(int i=0; i<nbrsect; i++)
    {
    nodecoord[0][i]=axispt1+xaxis*radius1*cos(i*(6.2832/nbrsect))
			+yaxis*radius1*sin(i*(6.2832/nbrsect));
    nodecoord[0][nbrsect+i]=axispt2+xaxis*radius2*cos(i*(6.2832/nbrsect))
			+yaxis*radius2*sin(i*(6.2832/nbrsect));
    edgenode1[i]=i;
    edgenode2[i]=(i+1)%nbrsect;
    edgenode1[nbrsect+i]=nbrsect+i;
    edgenode2[nbrsect+i]=nbrsect+(i+1)%nbrsect;
    edgenode1[2*nbrsect+i]=i;
    edgenode2[2*nbrsect+i]=nbrsect+i;
    numnodesurf[i][0]=i;
    numnodesurf[i][1]=(i+1)%nbrsect;
    numnodesurf[i][2]=nbrsect+(i+1)%nbrsect;
    numnodesurf[i][3]=nbrsect+i;
    numnodesurf[nbrsect][i]=i;
    numnodesurf[nbrsect+1][i]=nbrsect+i;
    }
  }

/*************************************************************************/

frustum::~frustum()
  {
  delete nbrnodes;
  for (int imth=0; imth<nbrmth; imth++) delete nodecoord[imth];
  delete nodecoord;
  delete edgenode1;
  delete edgenode2;
  delete edgecolor;
  delete nbrnodesurf;
  for (int isurf=0;isurf<nbrsurf;isurf++) delete numnodesurf[isurf];
  delete numnodesurf;
  delete surfcolor;
  }

/*************************************************************************/

void frustum::WriteShape(ostream &stream)
  {
  stream << axispt1 << axispt2;
  stream << radius1 << " " << radius2 << " " << nbrsect << "\n";
  stream << ecolor << " " << scolor << "\n";
  }

/*************************************************************************/

tyre::tyre(mth *mt, vec ap1, vec ap2,
	 double ri, double re, int nbs,
	 int eclr, int sclr)
  {
  mthref=new mth *[1];
  mthref[0]=mt;
  axispt1=ap1;
  axispt2=ap2;
  radint=ri;
  radext=re;
  nbrsect=nbs;
  ecolor=eclr;
  scolor=sclr;
  build();
  }

/*************************************************************************/

tyre::tyre(mth *mt, istream &stream)
  {
  mthref=new mth *[1];
  mthref[0]=mt;
  stream >> axispt1 >> axispt2;
  stream >> radint >> radext >> nbrsect >> ecolor >> scolor;
  build();
  }

/*************************************************************************/

void tyre::build()
  {
  nbrmth=1; nbrnodes=new int[nbrmth]; nbrnodes[0]=4*nbrsect; nbredges=7*nbrsect; nbrsurf=3*nbrsect;

  nodecoord=new vec *[nbrmth];
  for (int imth=0; imth<nbrmth; imth++) nodecoord[imth]=new vec[nbrnodes[imth]];
  edgenode1=new int[nbredges];
  edgenode2=new int[nbredges];
  edgecolor=new int[nbredges];
  for (int iedge=0; iedge<nbredges; iedge++) edgecolor[iedge]=ecolor;
  nbrnodesurf=new int[nbrsurf];
  numnodesurf=new int*[nbrsurf];
  surfcolor=new int[nbrsurf];
  int isurf;
  for (isurf=0; isurf<nbrsurf; isurf++) surfcolor[isurf]=scolor;
  for (isurf=0; isurf<nbrsurf; isurf++)
    {
    nbrnodesurf[isurf]=4;
    numnodesurf[isurf]=new int[4];
    }

  vec xaxis(1,0,0), yaxis(0,1,0), zaxis=axispt2-axispt1;
  if (zaxis.length()<1E-6) zaxis.put(0.0,0.0,1.0);
  if (fabs(zaxis.x) > fabs(zaxis.y))
    {
    xaxis=yaxis^zaxis; xaxis.unite();
    yaxis=zaxis^xaxis; yaxis.unite();
    }
  else
    {
    yaxis=zaxis^xaxis; yaxis.unite();
    xaxis=yaxis^zaxis; xaxis.unite();
    }

  for(int i=0; i<nbrsect; i++)
    {
    nodecoord[0][i]=axispt1+xaxis*radint*cos(i*(6.2832/nbrsect))
			+yaxis*radint*sin(i*(6.2832/nbrsect));
    nodecoord[0][nbrsect+i]=axispt1+xaxis*radext*cos(i*(6.2832/nbrsect))
			+yaxis*radext*sin(i*(6.2832/nbrsect));
    nodecoord[0][2*nbrsect+i]=axispt2+xaxis*radext*cos(i*(6.2832/nbrsect))
			+yaxis*radext*sin(i*(6.2832/nbrsect));
    nodecoord[0][3*nbrsect+i]=axispt2+xaxis*radint*cos(i*(6.2832/nbrsect))
			+yaxis*radint*sin(i*(6.2832/nbrsect));
    edgenode1[i]=i;
    edgenode2[i]=(i+1)%nbrsect;
    edgenode1[nbrsect+i]=i;
    edgenode2[nbrsect+i]=nbrsect+i;
    edgenode1[2*nbrsect+i]=nbrsect+i;
    edgenode2[2*nbrsect+i]=nbrsect+(i+1)%nbrsect;
    edgenode1[3*nbrsect+i]=nbrsect+i;
    edgenode2[3*nbrsect+i]=2*nbrsect+i;
    edgenode1[4*nbrsect+i]=2*nbrsect+i;
    edgenode2[4*nbrsect+i]=2*nbrsect+(i+1)%nbrsect;
    edgenode1[5*nbrsect+i]=2*nbrsect+i;
    edgenode2[5*nbrsect+i]=3*nbrsect+i;
    edgenode1[6*nbrsect+i]=3*nbrsect+i;
    edgenode2[6*nbrsect+i]=3*nbrsect+(i+1)%nbrsect;
    numnodesurf[i][0]=i;
    numnodesurf[i][1]=(i+1)%nbrsect;
    numnodesurf[i][2]=nbrsect+(i+1)%nbrsect;
    numnodesurf[i][3]=nbrsect+i;
    numnodesurf[nbrsect+i][0]=nbrsect+i;
    numnodesurf[nbrsect+i][1]=nbrsect+(i+1)%nbrsect;
    numnodesurf[nbrsect+i][2]=2*nbrsect+(i+1)%nbrsect;
    numnodesurf[nbrsect+i][3]=2*nbrsect+i;
    numnodesurf[2*nbrsect+i][0]=2*nbrsect+i;
    numnodesurf[2*nbrsect+i][1]=2*nbrsect+(i+1)%nbrsect;
    numnodesurf[2*nbrsect+i][2]=3*nbrsect+(i+1)%nbrsect;
    numnodesurf[2*nbrsect+i][3]=3*nbrsect+i;
    }
  }

/*************************************************************************/

tyre::~tyre()
  {
  delete nbrnodes;
  for (int imth=0; imth<nbrmth; imth++) delete nodecoord[imth];
  delete nodecoord;
  delete edgenode1;
  delete edgenode2;
  delete edgecolor;
  delete nbrnodesurf;
  for (int isurf=0;isurf<nbrsurf;isurf++) delete numnodesurf[isurf];
  delete numnodesurf;
  delete surfcolor;
  }

/*************************************************************************/

void tyre::WriteShape(ostream &stream)
  {
  stream << axispt1 << axispt2;
  stream << radint << " " << radext << " " << nbrsect << "\n";
  stream << ecolor << " " << scolor << "\n";
  };

/*************************************************************************/

sphere::sphere(mth *mt, vec ce, double r, int nsl, int nbs, int eclr, int sclr)
  {
  mthref=new mth *[1];
  mthref[0]=mt;
  center=ce;
  nbrslices=nsl;
  nbrsect=nbs;
  radius=r;
  ecolor=eclr;
  scolor=sclr;
  build();
  }

/*************************************************************************/

sphere::sphere(mth *mt, istream &stream)
  {
  mthref=new mth *[1];
  mthref[0]=mt;
  stream >> center;
  stream >> radius >> nbrslices >> nbrsect >> ecolor >> scolor;
  build();
  }

/*************************************************************************/

void sphere::build()
  {
  nbrmth=1; nbrnodes=new int[nbrmth]; nbrnodes[0]=nbrsect*nbrslices+2; 
  nbredges=nbrsect*(2*nbrslices+1); nbrsurf=nbrsect*(nbrslices+1);

  nodecoord=new vec *[nbrmth];
  for (int imth=0; imth<nbrmth; imth++) nodecoord[imth]=new vec[nbrnodes[imth]];
  edgenode1=new int[nbredges];
  edgenode2=new int[nbredges];
  edgecolor=new int[nbredges];
  for (int iedge=0; iedge<nbredges; iedge++) edgecolor[iedge]=ecolor;
  nbrnodesurf=new int[nbrsurf];
  numnodesurf=new int*[nbrsurf];
  surfcolor=new int[nbrsurf];
  int isurf;
  for (isurf=0; isurf<nbrsurf; isurf++) surfcolor[isurf]=scolor;
  for (int isurf=0; isurf<nbrsect; isurf++)
      {
      nbrnodesurf[isurf]=3;
      numnodesurf[isurf]=new int[3];
      nbrnodesurf[nbrsect*nbrslices+isurf]=3;
      numnodesurf[nbrsect*nbrslices+isurf]=new int[3];
      for (int jj=1; jj<nbrslices; jj++)
        {
        nbrnodesurf[nbrsect*jj+isurf]=4;
        numnodesurf[nbrsect*jj+isurf]=new int[4];
	}
      }
  double rp,zp,th1,th2;
  nodecoord[0][0]=center+vcoord(0,0,radius);
  nodecoord[0][nbrnodes[0]-1]=center+vcoord(0,0,-radius);
  for (int isli=0; isli<nbrslices; isli++)
    {
    th1=(isli+1)*3.14159265/(nbrslices+1.0);
    rp=radius*sin(th1);
    zp=radius*cos(th1);
    for (int isect=0; isect<nbrsect; isect++)
      {  
      th2=2*isect*3.14159265/nbrsect;
      nodecoord[0][isli*nbrsect+isect+1]=center
                                         +vcoord(rp*cos(th2),rp*sin(th2),zp);
      }
    }

  for (int isect=0; isect<nbrsect; isect++)
    {  
    edgenode1[isect]=0; edgenode2[isect]=isect+1;
    edgenode1[nbrsect+isect]=isect+1; 
    edgenode2[nbrsect+isect]=(isect+1)%nbrsect+1;
    for (int isli=1; isli<nbrslices; isli++)
      {
      edgenode1[2*isli*nbrsect+isect]=(isli-1)*nbrsect+isect+1;
      edgenode2[2*isli*nbrsect+isect]=isli*nbrsect+isect+1;
      edgenode1[2*isli*nbrsect+nbrsect+isect]=(isli)*nbrsect+isect+1;
      edgenode2[2*isli*nbrsect+nbrsect+isect]=(isli)*nbrsect+1+(isect+1)%nbrsect;
      }
    edgenode1[2*nbrslices*nbrsect+isect]=(nbrslices-1)*nbrsect+isect+1;
    edgenode2[2*nbrslices*nbrsect+isect]=nbrslices*nbrsect+1;
    }

  for (int isect=0; isect<nbrsect; isect++)
    {  
    numnodesurf[isect][0]=0;
    numnodesurf[isect][1]=isect+1 ;
    numnodesurf[isect][2]=1+(isect+1) % nbrsect;
    for (int isli=1; isli<nbrslices; isli++)
      {
      numnodesurf[isli*nbrsect+isect][0]=(isli-1)*nbrsect+isect+1;
      numnodesurf[isli*nbrsect+isect][1]=isli*nbrsect+isect+1;
      numnodesurf[isli*nbrsect+isect][2]=isli*nbrsect+1+(isect+1)%nbrsect;
      numnodesurf[isli*nbrsect+isect][3]=(isli-1)*nbrsect+1+(isect+1)%nbrsect;
      }
    numnodesurf[nbrslices*nbrsect+isect][0]=(nbrslices-1)*nbrsect+isect+1;
    numnodesurf[nbrslices*nbrsect+isect][1]=nbrslices*nbrsect+1;
    numnodesurf[nbrslices*nbrsect+isect][2]=(nbrslices-1)*nbrsect+1+(isect+1)%nbrsect;
    }
  }

/*************************************************************************/

sphere::~sphere()
  {
  delete nbrnodes;
  for (int imth=0; imth<nbrmth; imth++) delete nodecoord[imth];
  delete nodecoord;
  delete edgenode1;
  delete edgenode2;
  delete edgecolor;
  delete nbrnodesurf;
  for (int isurf=0;isurf<nbrsurf;isurf++) delete numnodesurf[isurf];
  delete numnodesurf;
  delete surfcolor;
  }

/*************************************************************************/

void sphere::WriteShape(ostream &stream)
  {
  stream << center << "\n";
  stream << radius << " " << nbrslices << " " << nbrsect << "\n";
  stream << ecolor << " " << scolor << "\n";
  }

/*************************************************************************/

gear::gear(mth *mt, vec ap1, vec ap2,
	 double mod, int nd, double b0, double r,
	 int eclr, int sclr)
  {
  mthref=new mth *[1];
  mthref[0]=mt;
  axispt1=ap1;
  axispt2=ap2;
  module=mod;
  nbrdents=nd;
  beta0=b0;
  holeradius=r;
  ecolor=eclr;
  scolor=sclr;
  build();
  }

/*************************************************************************/

gear::gear(mth *mt, istream &stream)
  {
  mthref=new mth *[1];
  mthref[0]=mt;
  stream >> axispt1 >> axispt2;
  stream >> module >> nbrdents >> beta0 >> holeradius >> ecolor >> scolor;
  build();
  }

/*************************************************************************/

void gear::build()
  {
  nbrmth=1; nbrnodes=new int[nbrmth]; nbrnodes[0]=22*nbrdents; nbredges=26*nbrdents; nbrsurf=13*nbrdents;

  nodecoord=new vec *[nbrmth];
  for (int imth=0; imth<nbrmth; imth++) nodecoord[imth]=new vec[nbrnodes[imth]];
  edgenode1=new int[nbredges];
  edgenode2=new int[nbredges];
  edgecolor=new int[nbredges];
  for (int iedge=0; iedge<nbredges; iedge++) edgecolor[iedge]=ecolor;
  nbrnodesurf=new int[nbrsurf];
  numnodesurf=new int*[nbrsurf];
  surfcolor=new int[nbrsurf];
  int isurf;
  for (isurf=0; isurf<nbrsurf; isurf++) surfcolor[isurf]=scolor;
  for (isurf=0; isurf<(11*nbrdents); isurf++)
    {
    nbrnodesurf[isurf]=4;
    numnodesurf[isurf]=new int[4];
    }
  for (isurf=(11*nbrdents); isurf<(13*nbrdents); isurf++)
    {
    nbrnodesurf[isurf]=13;
    numnodesurf[isurf]=new int[13];
    }
  //nbrnodesurf[nbrsect]=nbrsect; numnodesurf[nbrsect]=new int[nbrsect];
  //nbrnodesurf[nbrsect+1]=nbrsect; numnodesurf[nbrsect+1]=new int[nbrsect];

  vec xaxis(1,0,0), yaxis(0,1,0), zaxis=axispt2-axispt1;
  if (zaxis.length()<1E-6) zaxis.put(0.0,0.0,1.0);
  if (fabs(zaxis.x) > fabs(zaxis.y))
    {
    xaxis=yaxis^zaxis; xaxis.unite();
    yaxis=zaxis^xaxis; yaxis.unite();
    }
  else
    {
    yaxis=zaxis^xaxis; yaxis.unite();
    xaxis=yaxis^zaxis; xaxis.unite();
    }

  // Calcul des caracteristiques de la roue dentee
  int i,n,curnode;
  double mypi, rb, rp, ri, re, alpha, alphap, alphai, alphae, beta, 
         betadent, theta, thetap, thetai, thetae, R;
  mypi=acos(-1.0);
  rp=nbrdents*module/2.0;
  alphap=mypi/9;
  thetap=tan(alphap)-alphap;
  rb=rp*cos(alphap);
  if (holeradius<rp) ri=rp-1.25*module; else ri=rp-module;
  if (ri<rb) ri=rb;
  alphai=acos(rb/ri);
  thetai=tan(alphai)-alphai;
  if (holeradius<rp) re=rp+module; else re=rp+1.25*module;
  alphae=acos(rb/re);
  thetae=tan(alphae)-alphae;
  // Construction des nodes
  curnode=-1;
  for (n=0; n<nbrdents; n++)
    {
    betadent=2*n*mypi/nbrdents+beta0;
    for (i=4; i>=0; i--)
      {
      alpha=alphae*(4-i)/4.0+alphai*i/4.0;
      theta=tan(alpha)-alpha;
      beta=thetap-theta+mypi/(2*nbrdents);
      R=rb/cos(alpha);
      curnode++;
      nodecoord[0][curnode]=axispt1+R*xaxis*cos(betadent-beta)
                                +R*yaxis*sin(betadent-beta);
      nodecoord[0][curnode+11*nbrdents]=axispt2+R*xaxis*cos(betadent-beta)
                                +R*yaxis*sin(betadent-beta);
      }
    for (i=0; i<=4; i++)
      {
      alpha=alphae*(4-i)/4.0+alphai*i/4.0;
      theta=tan(alpha)-alpha;
      beta=thetap-theta+3.14159265/(2*nbrdents);
      R=rb/cos(alpha);
      curnode++;
      nodecoord[0][curnode]=axispt1+R*xaxis*cos(betadent+beta)
                                +R*yaxis*sin(betadent+beta);
      nodecoord[0][curnode+11*nbrdents]=axispt2+R*xaxis*cos(betadent+beta)
                                +R*yaxis*sin(betadent+beta);
      }
    beta=thetap-thetai+mypi/(2*nbrdents);
    nodecoord[0][10*nbrdents+n]=axispt1+holeradius*xaxis*cos(betadent+beta)
                                    +holeradius*yaxis*sin(betadent+beta);
    nodecoord[0][21*nbrdents+n]=axispt2+holeradius*xaxis*cos(betadent+beta)
                                    +holeradius*yaxis*sin(betadent+beta);
    }
  // Construction des aretes
  for (i=0; i<(10*nbrdents); i++)
    {
      edgenode1[i]=i;
      edgenode2[i]=(i+1)%(10*nbrdents);
      edgenode1[i+11*nbrdents]=edgenode1[i]+11*nbrdents;
      edgenode2[i+11*nbrdents]=edgenode2[i]+11*nbrdents;
    }
  for (i=0; i<nbrdents; i++)
    {
      edgenode1[10*nbrdents+i]=10*nbrdents+i;
      edgenode2[10*nbrdents+i]=(10*nbrdents)+((i+1)%nbrdents);
      edgenode1[21*nbrdents+i]=edgenode1[10*nbrdents+i]+11*nbrdents;
      edgenode2[21*nbrdents+i]=edgenode2[10*nbrdents+i]+11*nbrdents;
    }
  for (i=0; i<(nbrdents); i++)
    {
      edgenode1[22*nbrdents+4*i]=i*10;
      edgenode2[22*nbrdents+4*i]=i*10+11*nbrdents;
      edgenode1[22*nbrdents+4*i+1]=i*10+4;
      edgenode2[22*nbrdents+4*i+1]=i*10+4+11*nbrdents;
      edgenode1[22*nbrdents+4*i+2]=i*10+5;
      edgenode2[22*nbrdents+4*i+2]=i*10+5+11*nbrdents;
      edgenode1[22*nbrdents+4*i+3]=i*10+9;
      edgenode2[22*nbrdents+4*i+3]=i*10+9+11*nbrdents;
    }
  // Construction des faces
  for (i=0; i<10*nbrdents; i++)
    {
      numnodesurf[i][0]=i;  
      numnodesurf[i][1]=(i+1)%(10*nbrdents);  
      numnodesurf[i][2]=numnodesurf[i][1]+11*nbrdents;
      numnodesurf[i][3]=numnodesurf[i][0]+11*nbrdents;
    }
  for (i=0; i<nbrdents; i++)
    {
      numnodesurf[10*nbrdents+i][0]=10*nbrdents+i;
      numnodesurf[10*nbrdents+i][1]=(10*nbrdents)+((i+1)%nbrdents);
      numnodesurf[10*nbrdents+i][2]=numnodesurf[10*nbrdents+i][1]+11*nbrdents;
      numnodesurf[10*nbrdents+i][3]=numnodesurf[10*nbrdents+i][0]+11*nbrdents;
    }
  if (holeradius<rp) for (i=0; i<nbrdents; i++)
    {
      numnodesurf[11*nbrdents+i][0]=i*10;
      numnodesurf[11*nbrdents+i][1]=i*10+1;
      numnodesurf[11*nbrdents+i][2]=i*10+2;
      numnodesurf[11*nbrdents+i][3]=i*10+3;
      numnodesurf[11*nbrdents+i][4]=i*10+4;
      numnodesurf[11*nbrdents+i][5]=i*10+5;
      numnodesurf[11*nbrdents+i][6]=i*10+6;
      numnodesurf[11*nbrdents+i][7]=i*10+7;
      numnodesurf[11*nbrdents+i][8]=i*10+8;
      numnodesurf[11*nbrdents+i][9]=i*10+9;
      numnodesurf[11*nbrdents+i][10]=(i*10+10)%(10*nbrdents);
      numnodesurf[11*nbrdents+i][11]=10*nbrdents+(i+1)%nbrdents;
      numnodesurf[11*nbrdents+i][12]=10*nbrdents+i;

      numnodesurf[12*nbrdents+i][0]=i*10+11*nbrdents;
      numnodesurf[12*nbrdents+i][1]=i*10+1+11*nbrdents;
      numnodesurf[12*nbrdents+i][2]=i*10+2+11*nbrdents;
      numnodesurf[12*nbrdents+i][3]=i*10+3+11*nbrdents;
      numnodesurf[12*nbrdents+i][4]=i*10+4+11*nbrdents;
      numnodesurf[12*nbrdents+i][5]=i*10+5+11*nbrdents;
      numnodesurf[12*nbrdents+i][6]=i*10+6+11*nbrdents;
      numnodesurf[12*nbrdents+i][7]=i*10+7+11*nbrdents;
      numnodesurf[12*nbrdents+i][8]=i*10+8+11*nbrdents;
      numnodesurf[12*nbrdents+i][9]=i*10+9+11*nbrdents;
      numnodesurf[12*nbrdents+i][10]=(i*10+10)%(10*nbrdents)+11*nbrdents;
      numnodesurf[12*nbrdents+i][11]=10*nbrdents+(i+1)%nbrdents+11*nbrdents;
      numnodesurf[12*nbrdents+i][12]=10*nbrdents+i+11*nbrdents;
    }
  if (holeradius>rp) for (i=0; i<nbrdents; i++)
    {
      numnodesurf[11*nbrdents+i][0]=i*10+5;
      numnodesurf[11*nbrdents+i][1]=i*10+6;
      numnodesurf[11*nbrdents+i][2]=i*10+7;
      numnodesurf[11*nbrdents+i][3]=i*10+8;
      numnodesurf[11*nbrdents+i][4]=i*10+9;
      numnodesurf[11*nbrdents+i][5]=(i*10+10)%(10*nbrdents);
      numnodesurf[11*nbrdents+i][6]=(i*10+11)%(10*nbrdents);
      numnodesurf[11*nbrdents+i][7]=(i*10+12)%(10*nbrdents);
      numnodesurf[11*nbrdents+i][8]=(i*10+13)%(10*nbrdents);
      numnodesurf[11*nbrdents+i][9]=(i*10+14)%(10*nbrdents);
      numnodesurf[11*nbrdents+i][10]=(i*10+15)%(10*nbrdents);
      numnodesurf[11*nbrdents+i][11]=10*nbrdents+(i+1)%nbrdents;
      numnodesurf[11*nbrdents+i][12]=10*nbrdents+i;

      numnodesurf[12*nbrdents+i][0]=i*10+5+11*nbrdents;
      numnodesurf[12*nbrdents+i][1]=i*10+6+11*nbrdents;
      numnodesurf[12*nbrdents+i][2]=i*10+7+11*nbrdents;
      numnodesurf[12*nbrdents+i][3]=i*10+8+11*nbrdents;
      numnodesurf[12*nbrdents+i][4]=i*10+9+11*nbrdents;
      numnodesurf[12*nbrdents+i][5]=(i*10+10)+11*nbrdents;
      numnodesurf[12*nbrdents+i][6]=(i*10+11)%(10*nbrdents)+11*nbrdents;
      numnodesurf[12*nbrdents+i][7]=(i*10+12)%(10*nbrdents)+11*nbrdents;
      numnodesurf[12*nbrdents+i][8]=(i*10+13)%(10*nbrdents)+11*nbrdents;
      numnodesurf[12*nbrdents+i][9]=(i*10+14)%(10*nbrdents)+11*nbrdents;
      numnodesurf[12*nbrdents+i][10]=(i*10+15)%(10*nbrdents)+11*nbrdents;
      numnodesurf[12*nbrdents+i][11]=10*nbrdents+(i+1)%nbrdents+11*nbrdents;
      numnodesurf[12*nbrdents+i][12]=10*nbrdents+i+11*nbrdents;
    }
  }

/*************************************************************************/

gear::~gear()
  {
  delete nbrnodes;
  for (int imth=0; imth<nbrmth; imth++) delete nodecoord[imth];
  delete nodecoord;
  delete edgenode1;
  delete edgenode2;
  delete edgecolor;
  delete nbrnodesurf;
  for (int isurf=0;isurf<nbrsurf;isurf++) delete numnodesurf[isurf];
  delete numnodesurf;
  delete surfcolor;
  }

/*************************************************************************/

void gear::WriteShape(ostream &stream)
  {
  stream << axispt1 << axispt2;
  stream << module << " " << nbrdents << " " << beta0 << " " 
         << holeradius << "\n";
  stream << ecolor << " " << scolor << "\n";
  }

/*************************************************************************/

habfile::habfile(mth *mt, const char *fn)
  {
  mthref=new mth *[1];
  mthref[0]=mt;
  filename=new char[strlen(fn)+1];
  strcpy(filename,fn);
  build();
  }

/*************************************************************************/

habfile::habfile(mth *mt, istream &stream)
  {
  mthref=new mth *[1];
  mthref[0]=mt;
  char fn[100];
  stream >> fn;
  filename=new char[strlen(fn)+1];
  strcpy(filename,fn);
  build();
  }

/*************************************************************************/

void habfile::build()
  {
  nbrmth=1;
  ifstream infile(filename);
  if (!infile)
    {
    cout << "Cannot open " << filename << "\n";
    exit(1);
    }
  int ecolor;
  nbrmth=1; nbrnodes=new int[nbrmth];
  infile >> nbrnodes[0] >> nbredges >> ecolor;
  nodecoord=new vec *[nbrmth];
  for (int imth=0; imth<nbrmth; imth++) nodecoord[imth]=new vec[nbrnodes[imth]];
  for (int imthnode=0; imthnode<nbrnodes[0]; imthnode++)
    infile >> nodecoord[0][imthnode].x >> nodecoord[0][imthnode].y >> nodecoord[0][imthnode].z;

  edgenode1=new int[nbredges];
  edgenode2=new int[nbredges];
  edgecolor=new int[nbredges];
  for (int iedge=0; iedge<nbredges; iedge++)
    {
    infile >> edgenode1[iedge] >> edgenode2[iedge];
    edgenode1[iedge]--; edgenode2[iedge]--;
    edgecolor[iedge]=ecolor;
    }

  int scolor;
  infile >> nbrsurf >> scolor;
  nbrnodesurf=new int[nbrsurf];
  numnodesurf=new int*[nbrsurf];
  surfcolor=new int[nbrsurf];
  for (int isurf=0; isurf<nbrsurf; isurf++)
    {
    infile >> nbrnodesurf[isurf];
    numnodesurf[isurf]=new int[nbrnodesurf[isurf]];
    surfcolor[isurf]=scolor;
    for (int ipt=0; ipt<nbrnodesurf[isurf]; ipt++)
      {
      infile >> numnodesurf[isurf][ipt];
      numnodesurf[isurf][ipt]--;
      }
    }
  infile.close();
  }

/*************************************************************************/

habfile::~habfile()
  {
  delete nbrnodes;
  for (int imth=0; imth<nbrmth; imth++) delete nodecoord[imth];
  delete nodecoord;
  delete edgenode1;
  delete edgenode2;
  delete edgecolor;
  delete nbrnodesurf;
  for (int isurf=0;isurf<nbrsurf;isurf++) delete numnodesurf[isurf];
  delete numnodesurf;
  delete surfcolor;
  }

/*************************************************************************/

void habfile::WriteShape(ostream &stream)
  {
  stream << filename << "\n";
  }

/*************************************************************************/

var_surface::var_surface(mth **mt_per_pt,int n_pt,vec **pt,int sclr,int protection):shape(mt_per_pt,n_pt,pt,protection)
  {
	nbrpt=n_pt;
	scolor=sclr;
	build();
  }

/*************************************************************************/

var_surface::var_surface(mth **mt_per_pt, istream &stream):shape(mt_per_pt, stream)
  {
	stream >> nbrpt >> scolor;
	build();
  }

/*************************************************************************/

void var_surface::build()
 {
	if (nbrpt==2)
	{
		nbredges=1; nbrsurf=0;
  		edgenode1=new int[nbredges];
 		edgenode2=new int[nbredges];
  		edgecolor=new int[nbredges];
  		edgenode1[0]=0;
		edgenode2[0]=1;
		edgecolor[0]=scolor;
	}
	else
	{
		nbredges=0; nbrsurf=1;
		nbrnodesurf=new int[nbrsurf];
		numnodesurf=new int*[nbrsurf];
		surfcolor=new int[nbrsurf];
		nbrnodesurf[0]=nbrpt;
		numnodesurf[0]=new int[nbrnodesurf[0]];
		surfcolor[0]=scolor;
		for (int ipt=0; ipt<nbrnodesurf[0]; ipt++)
		{
			numnodesurf[0][ipt]=ipt;
		}
	}
 }

/*************************************************************************/

var_surface::~var_surface()
  {
	if (security==1)
	{
		if (nbrpt==2)
		{
			delete nbrnodes;
			for (int imth=0; imth<nbrmth; imth++) delete nodecoord[imth];
			delete nodecoord;
			delete edgenode1;
			delete edgenode2;
			delete edgecolor;
		}
		else
		{
			delete nbrnodes;
			for (int imth=0; imth<nbrmth; imth++) nodecoord[imth];
			delete nodecoord;
			delete nbrnodesurf;
			for (int isurf=0;isurf<nbrsurf;isurf++) delete numnodesurf[isurf];
			delete numnodesurf;
			delete surfcolor;
		}
	}
	else
	{
		if (nbrpt==2)
		{
			delete nbrnodes;
			delete edgenode1;
			delete edgenode2;
			delete edgecolor;
		}
		else
		{
			delete nbrnodes;
			delete nbrnodesurf;
			for (int isurf=0;isurf<nbrsurf;isurf++) delete numnodesurf[isurf];
			delete numnodesurf;
			delete surfcolor;
		}
	
// !!!! But the user has the responsability of deallocating the pointers : nodecoord[imth],nodecoord
	}
  }

/*************************************************************************/

void var_surface::WriteShape(ostream &stream)
{
  stream << "Pay attention. They are relative coordinates, but not with respect to only one frame !\n" ;
  for (int imth=0; imth<nbrmth; imth++)
  {
	stream << mthref[imth] << endl;
	for (int imthnode=0; imthnode<nbrnodes[imth]; imthnode++)
	{
		stream << nodecoord[imth][imthnode] << " ";	
	}
	stream << endl;
  }
  stream << scolor << "\n";
}

/*************************************************************************/

var_path::var_path(mth **mt_per_pt,int n_pt,vec **pt,int eclr,int protection):shape(mt_per_pt,n_pt,pt,protection)
  {
	nbrpt=n_pt;
	ecolor=eclr;
	build();
  }

/*************************************************************************/

var_path::var_path(mth **mt_per_pt, istream &stream):shape(mt_per_pt, stream)
  {
	stream >> nbrpt >> ecolor;
	build();
  }

/*************************************************************************/

void var_path::build()
{
	nbredges=nbrpt-1; nbrsurf=0;
	edgenode1=new int[nbredges];
	edgenode2=new int[nbredges];
	edgecolor=new int[nbredges];
	for (int iedge=0; iedge<nbredges; iedge++)
	{
		edgenode1[iedge]=iedge;
		edgenode2[iedge]=iedge+1;
		edgecolor[iedge]=ecolor;
	}
}

/*************************************************************************/

var_path::~var_path()
  {
	if (security==1)
  	{
		delete nbrnodes;
		for (int imth=0; imth<nbrmth; imth++) delete nodecoord[imth];
		delete nodecoord;
		delete edgenode1;
		delete edgenode2;
		delete edgecolor;
	}
	else
	{
// !!!! But the user has the responsability of deallocating the pointers : nodecoord[imth],nodecoord
		delete nbrnodes;
		delete edgenode1;
		delete edgenode2;
		delete edgecolor;
	}
  }

/*************************************************************************/

void var_path::WriteShape(ostream &stream)
{
  stream << "Pay attention. They are relative coordinates, but not with respect to only one frame !\n" ;
  for (int imth=0; imth<nbrmth; imth++)
  {
	stream << *mthref[imth];
	for (int imthnode=0; imthnode<nbrnodes[imth]; imthnode++)
	{
		stream << " (" << nodecoord[imth][imthnode].x << " " << nodecoord[imth][imthnode].y << " " << nodecoord[imth][imthnode].z << ")" << endl;
	}
	stream << endl;
  }
  stream << endl;
  stream << nbredges << endl;
  for (int iedge=0;iedge<nbredges;iedge++)
  {
	stream << edgenode1[iedge] << " " << edgenode2[iedge] << endl;
  }
  stream << endl;
  stream << ecolor << "\n";
}
/*************************************************************************/

spring::spring(mth *mt1, vec en1, mth* mt2, vec en2, 
               double r, double l, int nsp, int ns, int eclr)
  {
  mthref=new mth *[2];
  mthref[0]=mt1;
  end1=en1;
  mthref[1]=mt2;
  end2=en2;
  radius=r;
  lrod=l;
  nbrspirals=nsp;
  nbrsect=ns;
  ecolor=eclr;
  build();
  }

/*************************************************************************/

spring::spring(mth *mt, istream &stream)
  {
  mthref=new mth *[2];
  mthref[0]=mt;
  mthref[1]=mt;
  stream >> end1 >> end2 ;
  stream >> radius >> lrod >> nbrspirals >> nbrsect >> ecolor;
  build();
  }

/*************************************************************************/

void spring::buildnodes()
  {
  vec xaxis(1,0,0), yaxis(0,1,0), zaxis;
  zaxis=(mthref[0]->inv()*(mthref[1][0]*end2))-end1;
  double lstep=(zaxis.length()-2*lrod)/(nbrspirals+1.0);
  if (zaxis.length()<1E-6) zaxis.put(0.0,0.0,1.0);
  zaxis.unite();
  if (fabs(zaxis.x) > fabs(zaxis.y))
    {
    xaxis=yaxis^zaxis; xaxis.unite();
    yaxis=zaxis^xaxis; yaxis.unite();
    }
  else
    {
    yaxis=zaxis^xaxis; yaxis.unite();
    xaxis=yaxis^zaxis; xaxis.unite();
    }

  nodecoord[0][0]=end1;
  nodecoord[0][1]=end1+lrod*zaxis;
  for(int inode=1; inode<nbrnodes[0]-2; inode++)
    nodecoord[0][inode+1]=end1+zaxis*(lrod+lstep*(0.5+(inode-1)*1.0/nbrsect))
                        +xaxis*radius*cos((inode-1)*(6.2832/nbrsect))
			+yaxis*radius*sin((inode-1)*(6.2832/nbrsect));
  nodecoord[0][nbrnodes[0]-1]=end1+(lrod+lstep*(nbrspirals+1))*zaxis;
  nodecoord[1][0]=end2;
  }

/*************************************************************************/

void spring::build()
  {
  nbrmth=2; nbrnodes=new int[nbrmth]; 
  nbrnodes[0]=nbrspirals*nbrsect+4; nbrnodes[1]=1;
  nbredges=nbrspirals*nbrsect+4; nbrsurf=0;

  nodecoord=new vec *[nbrmth];
  for (int imth=0; imth<nbrmth; imth++) nodecoord[imth]=new vec[nbrnodes[imth]];
  edgenode1=new int[nbredges];
  edgenode2=new int[nbredges];
  edgecolor=new int[nbredges];
  for (int iedge=0; iedge<nbredges; iedge++) edgecolor[iedge]=ecolor;

  buildnodes(); 

  for (int iedge=0; iedge<nbredges; iedge++)
    {
    edgenode1[iedge]=iedge;
    edgenode2[iedge]=iedge+1;
    }
  }

/*************************************************************************/

spring::~spring()
  {
  delete nbrnodes;
  for (int imth=0; imth<nbrmth; imth++) delete nodecoord[imth];
  delete nodecoord;
  delete edgenode1;
  delete edgenode2;
  delete edgecolor;
  }

/*************************************************************************/

void spring::WriteNodes(ostream &stream, int initnode)
  {
  buildnodes();
  shape::WriteNodes(stream, initnode);
  }

/*************************************************************************/

void spring::WriteNodes(mth *mthvisu, ostream &stream, int initnode,
                       int norot)
  {
  buildnodes();
  shape::WriteNodes(mthvisu, stream, initnode, norot);
  }

/*************************************************************************/

void spring::WriteCoord(ostream &stream)
  {
  buildnodes();
  shape::WriteCoord(stream);
  }

/*************************************************************************/

void spring::WriteCoord(mth *mthvisu, ostream &stream, int norot)
  {
  buildnodes();
  shape::WriteCoord(mthvisu, stream, norot);
  }

/*************************************************************************/

void spring::WriteShape(ostream &stream)
  {
  stream << end1 << " " << end2 << "\n";
  stream << radius << " " << lrod << " " << nbrspirals << " " 
         << nbrsect << "\n";
  stream << ecolor << "\n";
  };

/*************************************************************************/

scene::scene()
  {
  nbrshapes=0;
  nbrnodes=0;
  nbredges=0;
  nbrsurf=0;
  firstshape=0;
  curshape=0;
  mthvisu=&mthvisudefault;
  norot=0;
  }

/*************************************************************************/

void scene::AddShape(shape *newshape)
  {
  nbrshapes++;
  nbrnodes+=newshape->GetNbrNodes();
  nbredges+=newshape->GetNbrEdges();
  nbrsurf+=newshape->GetNbrSurf();
  if (nbrshapes==1)
    {
    firstshape=newshape;
    curshape=newshape;
    curshape->nextshape=0;
    }
  else
    {
    curshape->nextshape=newshape;
    curshape=newshape;
    curshape->nextshape=0;
    }
  }

/*************************************************************************/

void scene::CreateVolFile(const char *filename)
  {
  ofstream volfile(filename);
  volfile << "scene" << "\n";
  int ishape,curinitnode,curinitedge,curinitsurf;
  shape *curshapetmp;
  // Sauvegarde des noeuds
  curshapetmp=firstshape;
  curinitnode=1;
  volfile << nbrnodes << "\n";
  for (ishape=0; ishape<nbrshapes; ishape++)
  {
    curshapetmp->WriteNodes(mthvisu,volfile,curinitnode,norot);
    curinitnode+=curshapetmp->GetNbrNodes();
    curshapetmp=curshapetmp->nextshape;
    }
  // Sauvegarde des aretes
  curshapetmp=firstshape;
  curinitnode=1;
  curinitedge=1;
  volfile << nbredges << "\n";
  for (ishape=0; ishape<nbrshapes; ishape++)
    {
    curshapetmp->WriteEdges(volfile,curinitnode,curinitedge);
    curinitnode+=curshapetmp->GetNbrNodes();
    curinitedge+=curshapetmp->GetNbrEdges();
    curshapetmp=curshapetmp->nextshape;
    }
  // Sauvegarde des surfaces
  curshapetmp=firstshape;
  curinitnode=1;
  curinitsurf=1;
  volfile << nbrsurf << "\n";
  for (ishape=0; ishape<nbrshapes; ishape++)
    {
    curshapetmp->WriteSurf(volfile,curinitnode,curinitsurf);
    curinitnode+=curshapetmp->GetNbrNodes();
    curinitsurf+=curshapetmp->GetNbrSurf();
    curshapetmp=curshapetmp->nextshape;
    }
  volfile.close();
  }

/*************************************************************************/

void scene::WriteCoord(ostream &stream)
  {
  int ishape;
  shape *curshapetmp;
  // Sauvegarde des noeuds
  curshapetmp=firstshape;
  for (ishape=0; ishape<nbrshapes; ishape++)
    {
    curshapetmp->WriteCoord(mthvisu,stream,norot);
    curshapetmp=curshapetmp->nextshape;
    }
  }

/*************************************************************************/

void scene::SetVisuFrame(mth *mthtmp, int nrt)
{
mthvisu=mthtmp;
norot=nrt;
}

/*************************************************************************/


