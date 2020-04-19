
# module CAGeM.py
# python port of EasyDyn symbolic generator CAGeM, standing for
# Computer Aided Generation of Motion
# This module provides functions and classes to generate a C++ EasyDyn application

from sympy import *
import re
import time
import os
t0 = time.clock()
#os.system("pause") 

################################################################################
#                        generic routine for simplification                    #
# Option simplification or various simplification options                      #
# should be implemented here                                                   #
################################################################################

def mysimplify(expr):
    return simplify(expr)

################################################################################
#                   routine for building symbolic arrays                       #
#          as a replacement of symarray which depends on numpy                 #
################################################################################

def mysymarray(ArrayName,ArraySize):
    myarray=[]
    for item in range(0,ArraySize):
        myarray.append(symbols(ArrayName+"["+str(item)+"]"))
    return myarray


################################################################################
#       routines for exporting symbolic expressions, vectors or mth to C       #
################################################################################

def exportExpressionToCProgram(outFile,VarName,expr):
    cexpr=ccode(expr)
    #cexpr=re.sub(r'_q__(\d+)',r'q[\1]',cexpr)
    #cexpr=re.sub(r'_qd__(\d+)',r'qd[\1]',cexpr)
    #cexpr=re.sub(r'_qdd__(\d+)',r'qdd[\1]',cexpr)
    #cexpr=re.sub(r'_p__(\d+)',r'p[\1]',cexpr)
    #cexpr=re.sub(r'_pd__(\d+)',r'pd[\1]',cexpr)
    #cexpr=re.sub(r'_pdd__(\d+)',r'pdd[\1]',cexpr)
    strtmp=VarName+"="+cexpr+";\n"
    outFile.write(strtmp)

def exportSymbolicVectorToCProgram(outFile,VecName,SymVec):
    if SymVec[0]!=0:
        VarName=VecName+".x"
        exportExpressionToCProgram(outFile,VarName,SymVec[0])
    if SymVec[1]!=0:
        VarName=VecName+".y"
        exportExpressionToCProgram(outFile,VarName,SymVec[1])
    if SymVec[2]!=0:
        VarName=VecName+".z"
        exportExpressionToCProgram(outFile,VarName,SymVec[2])

def exportSymbolicMTHToCProgram(outFile,MTHName,SymMTH):
    VarName=MTHName+".R.r11"
    exportExpressionToCProgram(outFile,VarName,SymMTH[0,0])
    VarName=MTHName+".R.r12"
    exportExpressionToCProgram(outFile,VarName,SymMTH[0,1])
    VarName=MTHName+".R.r13"
    exportExpressionToCProgram(outFile,VarName,SymMTH[0,2])
    VarName=MTHName+".R.r21"
    exportExpressionToCProgram(outFile,VarName,SymMTH[1,0])
    VarName=MTHName+".R.r22"
    exportExpressionToCProgram(outFile,VarName,SymMTH[1,1])
    VarName=MTHName+".R.r23"
    exportExpressionToCProgram(outFile,VarName,SymMTH[1,2])
    VarName=MTHName+".R.r31"
    exportExpressionToCProgram(outFile,VarName,SymMTH[2,0])
    VarName=MTHName+".R.r32"
    exportExpressionToCProgram(outFile,VarName,SymMTH[2,1])
    VarName=MTHName+".R.r33"
    exportExpressionToCProgram(outFile,VarName,SymMTH[2,2])
    VarName=MTHName+".e.x"
    exportExpressionToCProgram(outFile,VarName,SymMTH[0,3])
    VarName=MTHName+".e.y"
    exportExpressionToCProgram(outFile,VarName,SymMTH[1,3])
    VarName=MTHName+".e.z"
    exportExpressionToCProgram(outFile,VarName,SymMTH[2,3])

def exportSymbolicCoroDispToCProgram(outFile,MTHName,SymMTH):
    VarName=MTHName+"="
    c1expr=ccode(SymMTH[0,3])
    c2expr=ccode(SymMTH[1,3])
    c3expr=ccode(SymMTH[2,3])
    strtmp=VarName+"Tdisp("+c1expr+","+c2expr+","+c3expr+");\n"
    outFile.write(strtmp)
################################################################################
#       routines providing elementary homogeneous transformation matrices      #
################################################################################

def Trotx(theta):
    C=cos(theta)
    S=sin(theta)
    M=Matrix([[1,0,0,0],[0,C,-S,0],[0,S,C,0],[0,0,0,1]])
    return M

def Troty(theta):
    C=cos(theta)
    S=sin(theta)
    M=Matrix([[C,0,S,0],[0,1,0,0],[-S,0,C,0],[0,0,0,1]])
    return M

def Trotz(theta):
    C=cos(theta)
    S=sin(theta)
    M=Matrix([[C,-S,0,0],[S,C,0,0],[0,0,1,0],[0,0,0,1]])
    return M

def Tdisp(dx,dy,dz):
    M=Matrix([[1,0,0,dx],[0,1,0,dy],[0,0,1,dz],[0,0,0,1]])
    return M

################################################################################
#           Definition of class FrameClass =  frame for bodies                 #
################################################################################
class FrameClass(object):
    "Class representing a frame for a body and its kinematics"
    def __init__(self):
        self.T0F=Matrix([[0,0,0,0],[0,0,0,0],[0,0,0,0]])
        self.TrefF=Matrix([[0,0,0,0],[0,0,0,0],[0,0,0,0]])
        self.TtoF=Matrix([[0,0,0,0],[0,0,0,0],[0,0,0,0]])
        self.TcoroF0=Matrix([[0,0,0,0],[0,0,0,0],[0,0,0,0]])
        self.TtoFd=Matrix([[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]]) #Derivative
        self.FlexBodyRef=[-1,-1]
        self.vF=Matrix([[0],[0],[0]])
        self.aF=Matrix([[0],[0],[0]])
        self.omega=Matrix([[0],[0],[0]])
        self.omegad=Matrix([[0],[0],[0]])
        self.vFpartialq=[]
        self.vFpartialp=[]
        self.vFpartial=[]
        self.omegaPartialq=[]
        self.omegaPartialp=[]
        self.omegaPartial=[]
        
        self.ExtraResults=0
        
        self.PosX=0
        self.PosY=0
        self.PosZ=0
        self.VelX=0
        self.VelY=0
        self.VelZ=0
        self.AccX=0
        self.AccY=0
        self.AccZ=0
        self.OmegaX=0
        self.OmegaY=0
        self.OmegaZ=0
        self.OmegadX=0
        self.OmegadY=0
        self.OmegadZ=0
        
        self.mainframe=-1
        self.frameref=-1
        
    def ReferenceFrame(self,mainframe=-1,frameref=-1):
        self.mainframe=mainframe
        self.frameref=frameref    
    def ExportExtraRes(self,OutFile,BodyName):
        if(self.PosX==1):
            VarName=BodyName+".T0F.e.x << \" \" <<"
            OutFile.write(str(VarName))
            self.ExtraResults=self.ExtraResults+1
        if(self.PosY==1):
            VarName=BodyName+".T0F.e.y << \" \" <<"
            OutFile.write(str(VarName))
            self.ExtraResults=self.ExtraResults+1
        if(self.PosZ==1):
            VarName=BodyName+".T0F.e.z << \" \" <<"
            OutFile.write(str(VarName))
            self.ExtraResults=self.ExtraResults+1
        if(self.VelX==1):
            VarName=BodyName+".vF.x << \" \" <<"
            OutFile.write(str(VarName))
            self.ExtraResults=self.ExtraResults+1
        if(self.VelY==1):
            VarName=BodyName+".vF.y << \" \" <<"
            OutFile.write(str(VarName))
            self.ExtraResults=self.ExtraResults+1
        if(self.VelZ==1):
            VarName=BodyName+".vF.z << \" \" <<"
            OutFile.write(str(VarName))
            self.ExtraResults=self.ExtraResults+1
        if(self.AccX==1):
            VarName=BodyName+".aF.x << \" \" <<"
            OutFile.write(str(VarName))
            self.ExtraResults=self.ExtraResults+1
        if(self.AccY==1):
            VarName=BodyName+".aF.y << \" \" <<"
            OutFile.write(str(VarName))
            self.ExtraResults=self.ExtraResults+1
        if(self.AccZ==1):
            VarName=BodyName+".aF.z << \" \" <<"
            OutFile.write(str(VarName))
            self.ExtraResults=self.ExtraResults+1
        if(self.OmegaX==1):
            VarName=BodyName+".omega.x << \" \" <<"
            OutFile.write(str(VarName))
            self.ExtraResults=self.ExtraResults+1
        if(self.OmegaY==1):
            VarName=BodyName+".omega.y << \" \" <<"
            OutFile.write(str(VarName))
            self.ExtraResults=self.ExtraResults+1
        if(self.OmegaZ==1):
            VarName=BodyName+".omega.z << \" \" <<"
            OutFile.write(str(VarName))
            self.ExtraResults=self.ExtraResults+1
        if(self.OmegadX==1):
            VarName=BodyName+".omegad.x << \" \" <<"
            OutFile.write(str(VarName))
            self.ExtraResults=self.ExtraResults+1
        if(self.OmegadY==1):
            VarName=BodyName+".omegad.y << \" \" <<"
            OutFile.write(str(VarName))
            self.ExtraResults=self.ExtraResults+1
        if(self.OmegadZ==1):
            VarName=BodyName+".omegad.z << \" \" <<"
            OutFile.write(str(VarName))
            self.ExtraResults=self.ExtraResults+1
    def ExportExtraHeader(self,OutFile,BodyName):
        if(self.PosX==1):
            VarName=BodyName+"ex\" << \" \" <<"
            OutFile.write(str(VarName))
        if(self.PosY==1):
            VarName=BodyName+"ey\" << \" \" <<"
            OutFile.write(str(VarName))
        if(self.PosZ==1):
            VarName=BodyName+"ez\" << \" \" <<"
            OutFile.write(str(VarName))
        if(self.VelX==1):
            VarName=BodyName+"vx\" << \" \" <<"
            OutFile.write(str(VarName))
        if(self.VelY==1):
            VarName=BodyName+"vy\" << \" \" <<"
            OutFile.write(str(VarName))
        if(self.VelZ==1):
            VarName=BodyName+"vz\" << \" \" <<"
            OutFile.write(str(VarName))
        if(self.AccX==1):
            VarName=BodyName+"ax\" << \" \" <<"
            OutFile.write(str(VarName))
        if(self.AccY==1):
            VarName=BodyName+"ay\" << \" \" <<"
            OutFile.write(str(VarName))
        if(self.AccZ==1):
            VarName=BodyName+"az\" << \" \" <<"
            OutFile.write(str(VarName))
        if(self.OmegaX==1):
            VarName=BodyName+"wx\" << \" \" <<"
            OutFile.write(str(VarName))
        if(self.OmegaY==1):
            VarName=BodyName+"wy\" << \" \" <<"
            OutFile.write(str(VarName))
        if(self.OmegaZ==1):
            VarName=BodyName+"wz\" << \" \" <<"
            OutFile.write(str(VarName))
        if(self.OmegadX==1):
            VarName=BodyName+"wdx\" << \" \" <<"
            OutFile.write(str(VarName))
        if(self.OmegadY==1):
            VarName=BodyName+"wdy\" << \" \" <<"
            OutFile.write(str(VarName))
        if(self.OmegadZ==1):
            VarName=BodyName+"wdz\" << \" \" <<"
            OutFile.write(str(VarName))
    def ExportExtraPlot(self,OutFile,ibody,iframe,nameres,nbrcolumn):
        if(self.PosX==1):
            nbrcolumn=nbrcolumn+1
            OutFile.write("set term postscript eps color \"Times-Roman\" 20  \n")
            OutFile.write("set output \"Body"+str(ibody)+"_Frame"+str(iframe)+"_ex.eps\" \nset ylabel \"Position x Body ["+str(ibody)+"] Frame ["+str(iframe)+"]\" \nplot")    
            OutFile.write(" '"+str(nameres)+".res' using 1:"+str(nbrcolumn)+" title 'Body"+str(ibody)+"_Frame"+str(iframe)+"_ex' with line ")      
            OutFile.write("\nset term pop \nreplot \npause -1 'Next plot ?'\n\n")
        if(self.PosY==1):
            nbrcolumn=nbrcolumn+1
            OutFile.write("set term postscript eps color \"Times-Roman\" 20  \n")
            OutFile.write("set output \"Body"+str(ibody)+"_Frame"+str(iframe)+"_ey.eps\" \nset ylabel \"Position y Body ["+str(ibody)+"] Frame ["+str(iframe)+"]\" \nplot")    
            OutFile.write(" '"+str(nameres)+".res' using 1:"+str(nbrcolumn)+" title 'Body"+str(ibody)+"_Frame"+str(iframe)+"_ey' with line ")      
            OutFile.write("\nset term pop \nreplot \npause -1 'Next plot ?'\n\n")
        if(self.PosZ==1):
            nbrcolumn=nbrcolumn+1
            OutFile.write("set term postscript eps color \"Times-Roman\" 20  \n")
            OutFile.write("set output \"Body"+str(ibody)+"_Frame"+str(iframe)+"_ez.eps\" \nset ylabel \"Position z Body ["+str(ibody)+"] Frame ["+str(iframe)+"]\" \nplot")    
            OutFile.write(" '"+str(nameres)+".res' using 1:"+str(nbrcolumn)+" title 'Body"+str(ibody)+"_Frame"+str(iframe)+"_ez' with line ")      
            OutFile.write("\nset term pop \nreplot \npause -1 'Next plot ?'\n\n")
        if(self.VelX==1):
            nbrcolumn=nbrcolumn+1
            OutFile.write("set term postscript eps color \"Times-Roman\" 20  \n")
            OutFile.write("set output \"Body"+str(ibody)+"_Frame"+str(iframe)+"_vx.eps\" \nset ylabel \"Velocity x Body ["+str(ibody)+"] Frame ["+str(iframe)+"]\" \nplot")    
            OutFile.write(" '"+str(nameres)+".res' using 1:"+str(nbrcolumn)+" title 'Body"+str(ibody)+"_Frame"+str(iframe)+"_vx' with line ")      
            OutFile.write("\nset term pop \nreplot \npause -1 'Next plot ?'\n\n")
        if(self.VelY==1):
            nbrcolumn=nbrcolumn+1
            OutFile.write("set term postscript eps color \"Times-Roman\" 20  \n")
            OutFile.write("set output \"Body"+str(ibody)+"_Frame"+str(iframe)+"_vy.eps\" \nset ylabel \"Velocity y Body ["+str(ibody)+"] Frame ["+str(iframe)+"]\" \nplot")    
            OutFile.write(" '"+str(nameres)+".res' using 1:"+str(nbrcolumn)+" title 'Body"+str(ibody)+"_Frame"+str(iframe)+"_vy' with line ")      
            OutFile.write("\nset term pop \nreplot \npause -1 'Next plot ?'\n\n")
        if(self.VelZ==1):
            nbrcolumn=nbrcolumn+1
            OutFile.write("set term postscript eps color \"Times-Roman\" 20  \n")
            OutFile.write("set output \"Body"+str(ibody)+"_Frame"+str(iframe)+"_vz.eps\" \nset ylabel \"Velocity z Body ["+str(ibody)+"] Frame ["+str(iframe)+"]\" \nplot")    
            OutFile.write(" '"+str(nameres)+".res' using 1:"+str(nbrcolumn)+" title 'Body"+str(ibody)+"_Frame"+str(iframe)+"_vz' with line ")      
            OutFile.write("\nset term pop \nreplot \npause -1 'Next plot ?'\n\n")
        if(self.AccX==1):
            nbrcolumn=nbrcolumn+1
            OutFile.write("set term postscript eps color \"Times-Roman\" 20  \n")
            OutFile.write("set output \"Body"+str(ibody)+"_Frame"+str(iframe)+"_ax.eps\" \nset ylabel \"Acceleration x Body ["+str(ibody)+"] Frame ["+str(iframe)+"]\" \nplot")    
            OutFile.write(" '"+str(nameres)+".res' using 1:"+str(nbrcolumn)+" title 'Body"+str(ibody)+"_Frame"+str(iframe)+"_ax' with line ")      
            OutFile.write("\nset term pop \nreplot \npause -1 'Next plot ?'\n\n")
        if(self.AccY==1):
            nbrcolumn=nbrcolumn+1
            OutFile.write("set term postscript eps color \"Times-Roman\" 20  \n")
            OutFile.write("set output \"Body"+str(ibody)+"_Frame"+str(iframe)+"_ay.eps\" \nset ylabel \"Acceleration y Body ["+str(ibody)+"] Frame ["+str(iframe)+"]\" \nplot")    
            OutFile.write(" '"+str(nameres)+".res' using 1:"+str(nbrcolumn)+" title 'Body"+str(ibody)+"_Frame"+str(iframe)+"_ay' with line ")      
            OutFile.write("\nset term pop \nreplot \npause -1 'Next plot ?'\n\n")
        if(self.AccZ==1):
            nbrcolumn=nbrcolumn+1
            OutFile.write("set term postscript eps color \"Times-Roman\" 20  \n")
            OutFile.write("set output \"Body"+str(ibody)+"_Frame"+str(iframe)+"_az.eps\" \nset ylabel \"Acceleration z Body ["+str(ibody)+"] Frame ["+str(iframe)+"]\" \nplot")    
            OutFile.write(" '"+str(nameres)+".res' using 1:"+str(nbrcolumn)+" title 'Body"+str(ibody)+"_Frame"+str(iframe)+"_az' with line ")      
            OutFile.write("\nset term pop \nreplot \npause -1 'Next plot ?'\n\n")
        if(self.OmegaX==1):
            nbrcolumn=nbrcolumn+1
            OutFile.write("set term postscript eps color \"Times-Roman\" 20  \n")
            OutFile.write("set output \"Body"+str(ibody)+"_Frame"+str(iframe)+"_wx.eps\" \nset ylabel \"Omega x Body ["+str(ibody)+"] Frame ["+str(iframe)+"]\" \nplot")    
            OutFile.write(" '"+str(nameres)+".res' using 1:"+str(nbrcolumn)+" title 'Body"+str(ibody)+"_Frame"+str(iframe)+"_wx' with line ")      
            OutFile.write("\nset term pop \nreplot \npause -1 'Next plot ?'\n\n")
        if(self.OmegaY==1):
            nbrcolumn=nbrcolumn+1
            OutFile.write("set term postscript eps color \"Times-Roman\" 20  \n")
            OutFile.write("set output \"Body"+str(ibody)+"_Frame"+str(iframe)+"_wy.eps\" \nset ylabel \"Omega y Body ["+str(ibody)+"] Frame ["+str(iframe)+"]\" \nplot")    
            OutFile.write(" '"+str(nameres)+".res' using 1:"+str(nbrcolumn)+" title 'Body"+str(ibody)+"_Frame"+str(iframe)+"_wy' with line ")      
            OutFile.write("\nset term pop \nreplot \npause -1 'Next plot ?'\n\n")
        if(self.OmegaZ==1):
            nbrcolumn=nbrcolumn+1
            OutFile.write("set term postscript eps color \"Times-Roman\" 20  \n")
            OutFile.write("set output \"Body"+str(ibody)+"_Frame"+str(iframe)+"_wz.eps\" \nset ylabel \"Omega z Body ["+str(ibody)+"] Frame ["+str(iframe)+"]\" \nplot")    
            OutFile.write(" '"+str(nameres)+".res' using 1:"+str(nbrcolumn)+" title 'Body"+str(ibody)+"_Frame"+str(iframe)+"_wz' with line ")      
            OutFile.write("\nset term pop \nreplot \npause -1 'Next plot ?'\n\n")
        if(self.OmegadX==1):
            nbrcolumn=nbrcolumn+1
            OutFile.write("set term postscript eps color \"Times-Roman\" 20  \n")
            OutFile.write("set output \"Body"+str(ibody)+"_Frame"+str(iframe)+"_wdx.eps\" \nset ylabel \"Omegad x Body ["+str(ibody)+"] Frame ["+str(iframe)+"]\" \nplot")    
            OutFile.write(" '"+str(nameres)+".res' using 1:"+str(nbrcolumn)+" title 'Body"+str(ibody)+"_Frame"+str(iframe)+"_wdx' with line ")      
            OutFile.write("\nset term pop \nreplot \npause -1 'Next plot ?'\n\n")
        if(self.OmegadY==1):
            nbrcolumn=nbrcolumn+1
            OutFile.write("set term postscript eps color \"Times-Roman\" 20  \n")
            OutFile.write("set output \"Body"+str(ibody)+"_f"+str(iframe)+"_wdy.eps\" \nset ylabel \"Omegad y Body ["+str(ibody)+"] Frame ["+str(iframe)+"]\" \nplot")    
            OutFile.write(" '"+str(nameres)+".res' using 1:"+str(nbrcolumn)+" title 'Body"+str(ibody)+"_Frame"+str(iframe)+"_wdy' with line ")      
            OutFile.write("\nset term pop \nreplot \npause -1 'Next plot ?'\n\n")
        if(self.OmegadZ==1):
            nbrcolumn=nbrcolumn+1
            OutFile.write("set term postscript eps color \"Times-Roman\" 20  \n")
            OutFile.write("set output \"Body"+str(ibody)+"_Frame"+str(iframe)+"_wdz.eps\" \nset ylabel \"Omegad z Body ["+str(ibody)+"] Frame ["+str(iframe)+"]\" \nplot")    
            OutFile.write(" '"+str(nameres)+".res' using 1:"+str(nbrcolumn)+" title 'Body"+str(ibody)+"_Frame"+str(iframe)+"_wdz' with line ")      
            OutFile.write("\nset term pop \nreplot \npause -1 'Next plot ?'\n\n")
    def PlotPositionX(self):
        self.PosX=1
    def PlotPositionY(self):
        self.PosY=1
    def PlotPositionZ(self):
        self.PosZ=1
    def PlotVelocityX(self):
        self.VelX=1
    def PlotVelocityY(self):
        self.VelY=1
    def PlotVelocityZ(self):
        self.VelZ=1
    def PlotAccelerationX(self):
        self.AccX=1
    def PlotAccelerationY(self):
        self.AccY=1
    def PlotAccelerationZ(self):
        self.AccZ=1
    def PlotOmegaX(self):
        self.OmegaX=1
    def PlotOmegaY(self):
        self.OmegaY=1
    def PlotOmegaZ(self):
        self.OmegaZ=1
    def PlotOmegadX(self):
        self.OmegadX=1
    def PlotOmegadY(self):
        self.OmegadY=1
    def PlotOmegadZ(self):
        self.OmegadZ=1
    def affiche(self):
        if (self.mainframe<0):
            print('T0F=',self.TtoF)
            print('vF=',self.vF)
            print('aF=',self.aF)
            print('omega=',self.omega)
            print('omegad=',self.omegad)
        if (self.mainframe>=0):
            print('FlexBodyRef=',self.FlexBodyRef[0])
            print('TrefF=',self.TtoF)
            print('vFrel=',self.vF)
            print('aFrel=',self.aF)
            print('omegarel=',self.omega)
            print('omegadrel=',self.omegad)
    def ExportCoroDispData(self,OutFile,FrameName):
        VarName=FrameName+".TcoroF0"
        exportSymbolicCoroDispToCProgram(OutFile,VarName,self.TcoroF0)
    def ExportKinematics(self,OutFile,FrameName):
        if (self.mainframe<0):
            VarName=FrameName+".T0F"
            exportSymbolicMTHToCProgram(OutFile,VarName,self.TtoF)
            VarName=FrameName+".vF"
            exportSymbolicVectorToCProgram(OutFile,VarName,self.vF)
            VarName=FrameName+".aF"
            exportSymbolicVectorToCProgram(OutFile,VarName,self.aF)
            VarName=FrameName+".omega"
            exportSymbolicVectorToCProgram(OutFile,VarName,self.omega)
            VarName=FrameName+".omegad"
            exportSymbolicVectorToCProgram(OutFile,VarName,self.omegad)
        if (self.mainframe>=0):
            VarName=FrameName+".TrefF"
            exportSymbolicMTHToCProgram(OutFile,VarName,self.TtoF)
            VarName=FrameName+".vFrel"
            exportSymbolicVectorToCProgram(OutFile,VarName,self.vF)
            VarName=FrameName+".aFrel"
            exportSymbolicVectorToCProgram(OutFile,VarName,self.aF)
            VarName=FrameName+".omegarel"
            exportSymbolicVectorToCProgram(OutFile,VarName,self.omega)
            VarName=FrameName+".omegadrel"
            exportSymbolicVectorToCProgram(OutFile,VarName,self.omegad)
    def ExportPartialVelocities(self,OutFile,FrameName):
        for idof in range(0,len(self.vFpartial)):
            if (self.mainframe<0):
                VarName=FrameName+".vFpartial["+str(idof)+"]"
                exportSymbolicVectorToCProgram(OutFile,VarName,self.vFpartial[idof])
                VarName=FrameName+".omegapartial["+str(idof)+"]"
                exportSymbolicVectorToCProgram(OutFile,VarName,self.omegaPartial[idof])
            if (self.mainframe>=0):
                VarName=FrameName+".vFrelpartial["+str(idof)+"]"
                exportSymbolicVectorToCProgram(OutFile,VarName,self.vFpartial[idof])
                VarName=FrameName+".omegarelpartial["+str(idof)+"]"
                exportSymbolicVectorToCProgram(OutFile,VarName,self.omegaPartial[idof])
    #   #   #    
################################################################################
#               Definition of class BodyClass = Main Frame                     #
################################################################################

class BodyClass(object):
    "Class representing a rigid body and its kinematics"
    def __init__(self,Nnodes):
        self.T0F=Matrix([[0,0,0,0],[0,0,0,0],[0,0,0,0]])
        self.TrefF=Matrix([[0,0,0,0],[0,0,0,0],[0,0,0,0]])
        self.TtoF=Matrix([[0,0,0,0],[0,0,0,0],[0,0,0,0]])
        
        self.TtoFd=Matrix([[0,0,0,0],[0,0,0,0],[0,0,0,0],[0,0,0,0]]) #Derivative

        self.vF=Matrix([[0],[0],[0]])
        self.aF=Matrix([[0],[0],[0]])
        self.omega=Matrix([[0],[0],[0]])
        self.omegad=Matrix([[0],[0],[0]])
        self.vFpartialq=[]
        self.vFpartialp=[]
        self.vFpartial=[]
        self.omegaPartialq=[]
        self.omegaPartialp=[]
        self.omegaPartial=[]
        
        self.nbrnodes = Nnodes       
        self.frame=[]
        for iframe in range(0,self.nbrnodes):
            self.frame.append(FrameClass())
        
        self.ExtraResults=0
        
        self.PosX=0
        self.PosY=0
        self.PosZ=0
        self.VelX=0
        self.VelY=0
        self.VelZ=0
        self.AccX=0
        self.AccY=0
        self.AccZ=0
        self.OmegaX=0
        self.OmegaY=0
        self.OmegaZ=0
        self.OmegadX=0
        self.OmegadY=0
        self.OmegadZ=0
        
        self.mainframe=-1
        self.frameref=-1
        
    def Set(self,name="",mass=0,Ixx=0,Iyy=0,Izz=0,Ixy=0,Ixz=0,Iyz=0,length=0,section=0,rho=0,EYoung=0,nu=0,AlphaDamp=0,BetaDamp=0):
        self.name=name
        self.mass=mass
        self.length=length
        self.section=section
        self.rho=rho
        self.EYoung=EYoung
        self.nu=nu
        self.Ixx=Ixx
        self.Iyy=Iyy 
        self.Izz=Izz
        self.Ixy=Ixy 
        self.Ixz=Ixz
        self.Iyz=Iyz   
        self.AlphaDamp=AlphaDamp
        self.BetaDamp=BetaDamp   
    def ReferenceFrame(self,mainframe=-1,frameref=-1):
        self.mainframe=mainframe
        self.frameref=frameref
    def PlotPositionX(self):
        self.PosX=1
    def PlotPositionY(self):
        self.PosY=1
    def PlotPositionZ(self):
        self.PosZ=1
    def PlotVelocityX(self):
        self.VelX=1
    def PlotVelocityY(self):
        self.VelY=1
    def PlotVelocityZ(self):
        self.VelZ=1
    def PlotAccelerationX(self):
        self.AccX=1
    def PlotAccelerationY(self):
        self.AccY=1
    def PlotAccelerationZ(self):
        self.AccZ=1
    def PlotOmegaX(self):
        self.OmegaX=1
    def PlotOmegaY(self):
        self.OmegaY=1
    def PlotOmegaZ(self):
        self.OmegaZ=1
    def PlotOmegadX(self):
        self.OmegadX=1
    def PlotOmegadY(self):
        self.OmegadY=1
    def PlotOmegadZ(self):
        self.OmegadZ=1
    def ExportExtraRes(self,OutFile,BodyName):
        if(self.PosX==1):
            VarName=BodyName+".T0G.e.x << \" \" <<"
            OutFile.write(str(VarName))
            self.ExtraResults=self.ExtraResults+1
        if(self.PosY==1):
            VarName=BodyName+".T0G.e.y << \" \" <<"
            OutFile.write(str(VarName))
            self.ExtraResults=self.ExtraResults+1
        if(self.PosZ==1):
            VarName=BodyName+".T0G.e.z << \" \" <<"
            OutFile.write(str(VarName))
            self.ExtraResults=self.ExtraResults+1
        if(self.VelX==1):
            VarName=BodyName+".vG.x << \" \" <<"
            OutFile.write(str(VarName))
            self.ExtraResults=self.ExtraResults+1
        if(self.VelY==1):
            VarName=BodyName+".vG.y << \" \" <<"
            OutFile.write(str(VarName))
            self.ExtraResults=self.ExtraResults+1
        if(self.VelZ==1):
            VarName=BodyName+".vG.z << \" \" <<"
            OutFile.write(str(VarName))
            self.ExtraResults=self.ExtraResults+1
        if(self.AccX==1):
            VarName=BodyName+".aG.x << \" \" <<"
            OutFile.write(str(VarName))
            self.ExtraResults=self.ExtraResults+1
        if(self.AccY==1):
            VarName=BodyName+".aG.y << \" \" <<"
            OutFile.write(str(VarName))
            self.ExtraResults=self.ExtraResults+1
        if(self.AccZ==1):
            VarName=BodyName+".aG.z << \" \" <<"
            OutFile.write(str(VarName))
            self.ExtraResults=self.ExtraResults+1
        if(self.OmegaX==1):
            VarName=BodyName+".omega.x << \" \" <<"
            OutFile.write(str(VarName))
            self.ExtraResults=self.ExtraResults+1
        if(self.OmegaY==1):
            VarName=BodyName+".omega.y << \" \" <<"
            OutFile.write(str(VarName))
            self.ExtraResults=self.ExtraResults+1
        if(self.OmegaZ==1):
            VarName=BodyName+".omega.z << \" \" <<"
            OutFile.write(str(VarName))
            self.ExtraResults=self.ExtraResults+1
        if(self.OmegadX==1):
            VarName=BodyName+".omegad.x << \" \" <<"
            OutFile.write(str(VarName))
            self.ExtraResults=self.ExtraResults+1
        if(self.OmegadY==1):
            VarName=BodyName+".omegad.y << \" \" <<"
            OutFile.write(str(VarName))
            self.ExtraResults=self.ExtraResults+1
        if(self.OmegadZ==1):
            VarName=BodyName+".omegad.z << \" \" <<"
            OutFile.write(str(VarName))
            self.ExtraResults=self.ExtraResults+1
    def ExportExtraResPlusPlus(self,OutFile,BodyName):
        if(self.PosX==1):
            VarName=BodyName+".T0F.e.x << \" \" <<"
            OutFile.write(str(VarName))
            self.ExtraResults=self.ExtraResults+1
        if(self.PosY==1):
            VarName=BodyName+".T0F.e.y << \" \" <<"
            OutFile.write(str(VarName))
            self.ExtraResults=self.ExtraResults+1
        if(self.PosZ==1):
            VarName=BodyName+".T0F.e.z << \" \" <<"
            OutFile.write(str(VarName))
            self.ExtraResults=self.ExtraResults+1
        if(self.VelX==1):
            VarName=BodyName+".vF.x << \" \" <<"
            OutFile.write(str(VarName))
            self.ExtraResults=self.ExtraResults+1
        if(self.VelY==1):
            VarName=BodyName+".vF.y << \" \" <<"
            OutFile.write(str(VarName))
            self.ExtraResults=self.ExtraResults+1
        if(self.VelZ==1):
            VarName=BodyName+".vF.z << \" \" <<"
            OutFile.write(str(VarName))
            self.ExtraResults=self.ExtraResults+1
        if(self.AccX==1):
            VarName=BodyName+".aF.x << \" \" <<"
            OutFile.write(str(VarName))
            self.ExtraResults=self.ExtraResults+1
        if(self.AccY==1):
            VarName=BodyName+".aF.y << \" \" <<"
            OutFile.write(str(VarName))
            self.ExtraResults=self.ExtraResults+1
        if(self.AccZ==1):
            VarName=BodyName+".aF.z << \" \" <<"
            OutFile.write(str(VarName))
            self.ExtraResults=self.ExtraResults+1
        if(self.OmegaX==1):
            VarName=BodyName+".omega.x << \" \" <<"
            OutFile.write(str(VarName))
            self.ExtraResults=self.ExtraResults+1
        if(self.OmegaY==1):
            VarName=BodyName+".omega.y << \" \" <<"
            OutFile.write(str(VarName))
            self.ExtraResults=self.ExtraResults+1
        if(self.OmegaZ==1):
            VarName=BodyName+".omega.z << \" \" <<"
            OutFile.write(str(VarName))
            self.ExtraResults=self.ExtraResults+1
        if(self.OmegadX==1):
            VarName=BodyName+".omegad.x << \" \" <<"
            OutFile.write(str(VarName))
            self.ExtraResults=self.ExtraResults+1
        if(self.OmegadY==1):
            VarName=BodyName+".omegad.y << \" \" <<"
            OutFile.write(str(VarName))
            self.ExtraResults=self.ExtraResults+1
        if(self.OmegadZ==1):
            VarName=BodyName+".omegad.z << \" \" <<"
            OutFile.write(str(VarName))
            self.ExtraResults=self.ExtraResults+1
    def ExportExtraHeader(self,OutFile,BodyName):
        if(self.PosX==1):
            VarName=BodyName+"ex\" << \" \" <<"
            OutFile.write(str(VarName))
        if(self.PosY==1):
            VarName=BodyName+"ey\" << \" \" <<"
            OutFile.write(str(VarName))
        if(self.PosZ==1):
            VarName=BodyName+"ez\" << \" \" <<"
            OutFile.write(str(VarName))
        if(self.VelX==1):
            VarName=BodyName+"vx\" << \" \" <<"
            OutFile.write(str(VarName))
        if(self.VelY==1):
            VarName=BodyName+"vy\" << \" \" <<"
            OutFile.write(str(VarName))
        if(self.VelZ==1):
            VarName=BodyName+"vz\" << \" \" <<"
            OutFile.write(str(VarName))
        if(self.AccX==1):
            VarName=BodyName+"ax\" << \" \" <<"
            OutFile.write(str(VarName))
        if(self.AccY==1):
            VarName=BodyName+"ay\" << \" \" <<"
            OutFile.write(str(VarName))
        if(self.AccZ==1):
            VarName=BodyName+"az\" << \" \" <<"
            OutFile.write(str(VarName))
        if(self.OmegaX==1):
            VarName=BodyName+"wx\" << \" \" <<"
            OutFile.write(str(VarName))
        if(self.OmegaY==1):
            VarName=BodyName+"wy\" << \" \" <<"
            OutFile.write(str(VarName))
        if(self.OmegaZ==1):
            VarName=BodyName+"wz\" << \" \" <<"
            OutFile.write(str(VarName))
        if(self.OmegadX==1):
            VarName=BodyName+"wdx\" << \" \" <<"
            OutFile.write(str(VarName))
        if(self.OmegadY==1):
            VarName=BodyName+"wdy\" << \" \" <<"
            OutFile.write(str(VarName))
        if(self.OmegadZ==1):
            VarName=BodyName+"wdz\" << \" \" <<"
            OutFile.write(str(VarName))
    def ExportExtraPlot(self,OutFile,ibody,nameres,nbrcolumn):
        if(self.PosX==1):
            nbrcolumn=nbrcolumn+1
            OutFile.write("set term postscript eps color \"Times-Roman\" 20  \n")
            OutFile.write("set output \"Body"+str(ibody)+"_ex.eps\" \nset ylabel \"Position x Body ["+str(ibody)+"]\" \nplot")    
            OutFile.write(" '"+str(nameres)+".res' using 1:"+str(nbrcolumn)+" title 'Body"+str(ibody)+"_ex' with line ")      
            OutFile.write("\nset term pop \nreplot \npause -1 'Next plot ?'\n\n")
        if(self.PosY==1):
            nbrcolumn=nbrcolumn+1
            OutFile.write("set term postscript eps color \"Times-Roman\" 20  \n")
            OutFile.write("set output \"Body"+str(ibody)+"_ey.eps\" \nset ylabel \"Position y Body ["+str(ibody)+"]\" \nplot")    
            OutFile.write(" '"+str(nameres)+".res' using 1:"+str(nbrcolumn)+" title 'Body"+str(ibody)+"_ey' with line ")      
            OutFile.write("\nset term pop \nreplot \npause -1 'Next plot ?'\n\n")
        if(self.PosZ==1):
            nbrcolumn=nbrcolumn+1
            OutFile.write("set term postscript eps color \"Times-Roman\" 20  \n")
            OutFile.write("set output \"Body"+str(ibody)+"_ez.eps\" \nset ylabel \"Position z Body ["+str(ibody)+"]\" \nplot")    
            OutFile.write(" '"+str(nameres)+".res' using 1:"+str(nbrcolumn)+" title 'Body"+str(ibody)+"_ez' with line ")      
            OutFile.write("\nset term pop \nreplot \npause -1 'Next plot ?'\n\n")
        if(self.VelX==1):
            nbrcolumn=nbrcolumn+1
            OutFile.write("set term postscript eps color \"Times-Roman\" 20  \n")
            OutFile.write("set output \"Body"+str(ibody)+"_vx.eps\" \nset ylabel \"Velocity x Body ["+str(ibody)+"]\" \nplot")    
            OutFile.write(" '"+str(nameres)+".res' using 1:"+str(nbrcolumn)+" title 'Body"+str(ibody)+"_vx' with line ")      
            OutFile.write("\nset term pop \nreplot \npause -1 'Next plot ?'\n\n")
        if(self.VelY==1):
            nbrcolumn=nbrcolumn+1
            OutFile.write("set term postscript eps color \"Times-Roman\" 20  \n")
            OutFile.write("set output \"Body"+str(ibody)+"_vy.eps\" \nset ylabel \"Velocity y Body ["+str(ibody)+"]\" \nplot")    
            OutFile.write(" '"+str(nameres)+".res' using 1:"+str(nbrcolumn)+" title 'Body"+str(ibody)+"_vy' with line ")      
            OutFile.write("\nset term pop \nreplot \npause -1 'Next plot ?'\n\n")
        if(self.VelZ==1):
            nbrcolumn=nbrcolumn+1
            OutFile.write("set term postscript eps color \"Times-Roman\" 20  \n")
            OutFile.write("set output \"Body"+str(ibody)+"_vz.eps\" \nset ylabel \"Velocity z Body ["+str(ibody)+"]\" \nplot")    
            OutFile.write(" '"+str(nameres)+".res' using 1:"+str(nbrcolumn)+" title 'Body"+str(ibody)+"_vz' with line ")      
            OutFile.write("\nset term pop \nreplot \npause -1 'Next plot ?'\n\n")
        if(self.AccX==1):
            nbrcolumn=nbrcolumn+1
            OutFile.write("set term postscript eps color \"Times-Roman\" 20  \n")
            OutFile.write("set output \"Body"+str(ibody)+"_ax.eps\" \nset ylabel \"Acceleration x Body ["+str(ibody)+"]\" \nplot")    
            OutFile.write(" '"+str(nameres)+".res' using 1:"+str(nbrcolumn)+" title 'Body"+str(ibody)+"_ax' with line ")      
            OutFile.write("\nset term pop \nreplot \npause -1 'Next plot ?'\n\n")
        if(self.AccY==1):
            nbrcolumn=nbrcolumn+1
            OutFile.write("set term postscript eps color \"Times-Roman\" 20  \n")
            OutFile.write("set output \"Body"+str(ibody)+"_ay.eps\" \nset ylabel \"Acceleration y Body ["+str(ibody)+"]\" \nplot")    
            OutFile.write(" '"+str(nameres)+".res' using 1:"+str(nbrcolumn)+" title 'Body"+str(ibody)+"_ay' with line ")      
            OutFile.write("\nset term pop \nreplot \npause -1 'Next plot ?'\n\n")
        if(self.AccZ==1):
            nbrcolumn=nbrcolumn+1
            OutFile.write("set term postscript eps color \"Times-Roman\" 20  \n")
            OutFile.write("set output \"Body"+str(ibody)+"_az.eps\" \nset ylabel \"Acceleration z Body ["+str(ibody)+"]\" \nplot")    
            OutFile.write(" '"+str(nameres)+".res' using 1:"+str(nbrcolumn)+" title 'Body"+str(ibody)+"_az' with line ")      
            OutFile.write("\nset term pop \nreplot \npause -1 'Next plot ?'\n\n")
        if(self.OmegaX==1):
            nbrcolumn=nbrcolumn+1
            OutFile.write("set term postscript eps color \"Times-Roman\" 20  \n")
            OutFile.write("set output \"Body"+str(ibody)+"_wx.eps\" \nset ylabel \"Omega x Body ["+str(ibody)+"]\" \nplot")    
            OutFile.write(" '"+str(nameres)+".res' using 1:"+str(nbrcolumn)+" title 'Body"+str(ibody)+"_wx' with line ")      
            OutFile.write("\nset term pop \nreplot \npause -1 'Next plot ?'\n\n")
        if(self.OmegaY==1):
            nbrcolumn=nbrcolumn+1
            OutFile.write("set term postscript eps color \"Times-Roman\" 20  \n")
            OutFile.write("set output \"Body"+str(ibody)+"_wy.eps\" \nset ylabel \"Omega y Body ["+str(ibody)+"]\" \nplot")    
            OutFile.write(" '"+str(nameres)+".res' using 1:"+str(nbrcolumn)+" title 'Body"+str(ibody)+"_wy' with line ")      
            OutFile.write("\nset term pop \nreplot \npause -1 'Next plot ?'\n\n")
        if(self.OmegaZ==1):
            nbrcolumn=nbrcolumn+1
            OutFile.write("set term postscript eps color \"Times-Roman\" 20  \n")
            OutFile.write("set output \"Body"+str(ibody)+"_wz.eps\" \nset ylabel \"Omega z Body ["+str(ibody)+"]\" \nplot")    
            OutFile.write(" '"+str(nameres)+".res' using 1:"+str(nbrcolumn)+" title 'Body"+str(ibody)+"_wz' with line ")      
            OutFile.write("\nset term pop \nreplot \npause -1 'Next plot ?'\n\n")
        if(self.OmegadX==1):
            nbrcolumn=nbrcolumn+1
            OutFile.write("set term postscript eps color \"Times-Roman\" 20  \n")
            OutFile.write("set output \"Body"+str(ibody)+"_wdx.eps\" \nset ylabel \"Omegad x Body ["+str(ibody)+"]\" \nplot")    
            OutFile.write(" '"+str(nameres)+".res' using 1:"+str(nbrcolumn)+" title 'Body"+str(ibody)+"_wdx' with line ")      
            OutFile.write("\nset term pop \nreplot \npause -1 'Next plot ?'\n\n")
        if(self.OmegadY==1):
            nbrcolumn=nbrcolumn+1
            OutFile.write("set term postscript eps color \"Times-Roman\" 20  \n")
            OutFile.write("set output \"Body"+str(ibody)+"_wdy.eps\" \nset ylabel \"Omegad y Body ["+str(ibody)+"]\" \nplot")    
            OutFile.write(" '"+str(nameres)+".res' using 1:"+str(nbrcolumn)+" title 'Body"+str(ibody)+"_wdy' with line ")      
            OutFile.write("\nset term pop \nreplot \npause -1 'Next plot ?'\n\n")
        if(self.OmegadZ==1):
            nbrcolumn=nbrcolumn+1
            OutFile.write("set term postscript eps color \"Times-Roman\" 20  \n")
            OutFile.write("set output \"Body"+str(ibody)+"_wdz.eps\" \nset ylabel \"Omegad z Body ["+str(ibody)+"]\" \nplot")    
            OutFile.write(" '"+str(nameres)+".res' using 1:"+str(nbrcolumn)+" title 'Body"+str(ibody)+"_wdz' with line ")      
            OutFile.write("\nset term pop \nreplot \npause -1 'Next plot ?'\n\n")
    def affiche(self):
        print('mass=',self.mass)
        print('PhiG=',self.Ixx,self.Iyy,self.Izz,self.Ixy,self.Ixz,self.Iyz)
        if (self.mainframe<0):
            print('T0F=',self.TtoF)
            print('vG=',self.vF)
            print('aG=',self.aF)
            print('omega=',self.omega)
            print('omegad=',self.omegad)
        if (self.mainframe>=0):
            print('BodyRef=',self.mainframe)
            print('TrefF=',self.TtoF)
            print('vGrel=',self.vF)
            print('aGrel=',self.aF)
            print('omegarel=',self.omega)
            print('omegadrel=',self.omegad)
    def ExportInertiaData(self,OutFile,BodyName):
        VarName=BodyName+".mass"
        exportExpressionToCProgram(OutFile,VarName,self.mass)
        VarName=BodyName+".PhiG.Ixx"
        exportExpressionToCProgram(OutFile,VarName,self.Ixx)
        VarName=BodyName+".PhiG.Iyy"
        exportExpressionToCProgram(OutFile,VarName,self.Iyy)
        VarName=BodyName+".PhiG.Izz"
        exportExpressionToCProgram(OutFile,VarName,self.Izz)
        VarName=BodyName+".PhiG.Ixy"
        exportExpressionToCProgram(OutFile,VarName,self.Ixy)
        VarName=BodyName+".PhiG.Ixz"
        exportExpressionToCProgram(OutFile,VarName,self.Ixz)
        VarName=BodyName+".PhiG.Iyz"
        exportExpressionToCProgram(OutFile,VarName,self.Iyz)
    def ExportKinematics(self,OutFile,BodyName):
        if (self.mainframe<0):
            VarName=BodyName+".T0G"
            exportSymbolicMTHToCProgram(OutFile,VarName,self.TtoF)
            VarName=BodyName+".vG"
            exportSymbolicVectorToCProgram(OutFile,VarName,self.vF)
            VarName=BodyName+".aG"
            exportSymbolicVectorToCProgram(OutFile,VarName,self.aF)
            VarName=BodyName+".omega"
            exportSymbolicVectorToCProgram(OutFile,VarName,self.omega)
            VarName=BodyName+".omegad"
            exportSymbolicVectorToCProgram(OutFile,VarName,self.omegad)
        if (self.mainframe>=0):
            VarName=BodyName+".TrefG"
            exportSymbolicMTHToCProgram(OutFile,VarName,self.TtoF)
            VarName=BodyName+".vGrel"
            exportSymbolicVectorToCProgram(OutFile,VarName,self.vF)
            VarName=BodyName+".aGrel"
            exportSymbolicVectorToCProgram(OutFile,VarName,self.aF)
            VarName=BodyName+".omegarel"
            exportSymbolicVectorToCProgram(OutFile,VarName,self.omega)
            VarName=BodyName+".omegadrel"
            exportSymbolicVectorToCProgram(OutFile,VarName,self.omegad)
    def ExportPartialVelocities(self,OutFile,BodyName):
        for idof in range(0,len(self.vFpartial)):
            if (self.mainframe<0):
                VarName=BodyName+".vGpartial["+str(idof)+"]"
                exportSymbolicVectorToCProgram(OutFile,VarName,self.vFpartial[idof])
                VarName=BodyName+".omegapartial["+str(idof)+"]"
                exportSymbolicVectorToCProgram(OutFile,VarName,self.omegaPartial[idof])
            if (self.mainframe>=0):
                VarName=BodyName+".vGrelpartial["+str(idof)+"]"
                exportSymbolicVectorToCProgram(OutFile,VarName,self.vFpartial[idof])
                VarName=BodyName+".omegarelpartial["+str(idof)+"]"
                exportSymbolicVectorToCProgram(OutFile,VarName,self.omegaPartial[idof])
    def ExportKinematicsPlusPlus(self,OutFile,BodyName):
        if (self.mainframe<0):
            VarName=BodyName+".T0F"
            exportSymbolicMTHToCProgram(OutFile,VarName,self.TtoF)
            if(self.mass>0 and self.Ixx>0 and self.Iyy>0 and self.Izz>0):
                VarName=BodyName+".vF"
                exportSymbolicVectorToCProgram(OutFile,VarName,self.vF)
                VarName=BodyName+".aF"
                exportSymbolicVectorToCProgram(OutFile,VarName,self.aF)
            VarName=BodyName+".omega"
            exportSymbolicVectorToCProgram(OutFile,VarName,self.omega)
            if(self.mass>0 and self.Ixx>0 and self.Iyy>0 and self.Izz>0):
                VarName=BodyName+".omegad"
                exportSymbolicVectorToCProgram(OutFile,VarName,self.omegad)
        if (self.mainframe>=0):
            VarName=BodyName+".TrefF"
            exportSymbolicMTHToCProgram(OutFile,VarName,self.TtoF)
            if(self.mass>0 and self.Ixx>0 and self.Iyy>0 and self.Izz>0):
                VarName=BodyName+".vFrel"
                exportSymbolicVectorToCProgram(OutFile,VarName,self.vF)
                VarName=BodyName+".aFrel"
                exportSymbolicVectorToCProgram(OutFile,VarName,self.aF)
            VarName=BodyName+".omegarel"
            exportSymbolicVectorToCProgram(OutFile,VarName,self.omega)
            if(self.mass>0 and self.Ixx>0 and self.Iyy>0 and self.Izz>0):
                VarName=BodyName+".omegadrel"
                exportSymbolicVectorToCProgram(OutFile,VarName,self.omegad)
    def ExportPartialVelocitiesPlusPlus(self,OutFile,BodyName):
        for idof in range(0,len(self.vFpartial)):
            if (self.mainframe<0):
                if(self.mass>0 and self.Ixx>0 and self.Iyy>0 and self.Izz>0):
                    VarName=BodyName+".vFpartial["+str(idof)+"]"
                    exportSymbolicVectorToCProgram(OutFile,VarName,self.vFpartial[idof])
                VarName=BodyName+".omegapartial["+str(idof)+"]"
                exportSymbolicVectorToCProgram(OutFile,VarName,self.omegaPartial[idof])
            if (self.mainframe>=0):
                if(self.mass>0 and self.Ixx>0 and self.Iyy>0 and self.Izz>0):
                    VarName=BodyName+".vFrelpartial["+str(idof)+"]"
                    exportSymbolicVectorToCProgram(OutFile,VarName,self.vFpartial[idof])
                VarName=BodyName+".omegarelpartial["+str(idof)+"]"
                exportSymbolicVectorToCProgram(OutFile,VarName,self.omegaPartial[idof])
    #   #   #   
################################################################################
#             Definition of class MBSClass = multibody system                  #
################################################################################

class MBSClass(object):
    "Class representing a multibody system"
    def __init__(self,nbrbody,nbrdof,nbrdep,ApplicationTitle,ApplicationFileName,nbrframe=0):
        
        self.nbrRigidBody=0 # Number of rigid bodies
        self.nbrFlexBeam=0  # Number of flexible beams
        self.nbrFlexBody=0  # Number of flexible bodies
        
        self.nbrbody=nbrbody
        self.nbrdof=nbrdof
        self.nbrdep=nbrdep
        self.ApplicationTitle=ApplicationTitle
        self.ApplicationFileName=ApplicationFileName  
        
        self.nbrinput=0
        self.input=[]
        self.nbroutput=0
        self.output=[]
        self.nbrexternalinput=0
        self.externalinput=[]        
            
        # Init Table for Frames
        self.nbrframeperbody=[0 for idof in range(0,self.nbrbody)] #init table for nbr frames per body
        for ibody in range(0,(self.nbrbody)):
            if(nbrframe==0):
                self.nbrframeperbody[ibody]=0   
            else:      
                self.nbrframeperbody[ibody]=nbrframe[ibody] 
        # Init Bodies
        self.body=[]
        for ibody in range(0,self.nbrbody):
            if (nbrframe==0):
                self.body.append(BodyClass(self.nbrframeperbody[ibody]))
            else:
                self.body.append(BodyClass(self.nbrframeperbody[ibody]))
        self.t=symbols('t')
        self.q=mysymarray('q',nbrdof)
        self.qd=mysymarray('qd',nbrdof)
        self.qdd=mysymarray('qdd',nbrdof)
        self.p=mysymarray('p',nbrdep)
        self.pd=mysymarray('pd',nbrdep)
        self.pdd=mysymarray('pdd',nbrdep)
        self.pexpr=[]
        self.pdexpr=[]
        self.pddexpr=[]
        self.qini=[0 for idof in range(0,self.nbrdof)] 
        self.qdini=[0 for idof in range(0,self.nbrdof)]
        self.qdoflocked=[0 for idof in range(0,self.nbrdof)]
        self.lockintegration=0 # variable to preserve doflocked for the integration in NewmarkIntegration()
        for idep in range(0,nbrdep):
            self.pexpr.append(0)
            self.pdexpr.append(0)
            self.pddexpr.append(0)
        self.gravity=[0,0,0]
        self.tfinal=2
        self.hsave=0.01
        self.hmax=0.005
        self.ForcesInAppEff=" "
        self.GlobalVar=" "
        self.MainVar=" "
        self.SaveDataVar=" " 
        self.ComputeResidualVar=" "
        self.OverWriteNbrdof=0
        
        self.StaticFlag=0
        self.PoleFlag=0
        self.PerformanceTestFlag=0
        self.EASYDYNPlusPlusFlag=0
        
        self.StartTime_Kinematics=0
        self.EndTime_Kinematics=0
        print('MBS object initialized')      
    def Getqpt(self):
        return self.q, self.p, self.pexpr, self.t     
    def PrintSystemSize(self):
        print('nbrbody=',self.nbrbody)
        print('nbrdof=',self.nbrdof)
        print('nbrdep=',self.nbrdep)
    # -------- end MBSClass.affiche
    # ------ begin MBSClass.SetGravity
    def SetGravity(self,gx,gy,gz):
        self.gravity[0]=gx
        self.gravity[1]=gy
        self.gravity[2]=gz
    # -------- end MBSClass.SetGravity
    # ------ begin MBSClass.SetIntegrationParameters
    def SetIntegrationParameters(self,tfinal,hsave,hmax,newmarkonestep=0):
        self.tfinal=tfinal
        self.hsave=hsave
        self.hmax=hmax
        self.newmarkonestep=newmarkonestep
        print('Integration Parameters ')
    def Force(self,ForcesInAppEff):
        self.ForcesInAppEff=ForcesInAppEff
    def ChangeNbrDof(self,OverWriteNbrdof):
        self.OverWriteNbrdof=OverWriteNbrdof
    def SetGlobalVar(self,GlobalVar):
        self.GlobalVar=GlobalVar
    def SetMainVar(self,MainVar):
        self.MainVar=MainVar
    def DefineSaveData(self,SaveDataVar):
        self.SaveDataVar=SaveDataVar
    def DefineComputeResidual(self,ComputeResidualVar):
        self.ComputeResidualVar=ComputeResidualVar
    def EasyDynFlags(self,STATIC=0,POLE=0,TEST=0,EASEASYDYNPlusPlus=0):
        self.StaticFlag=STATIC
        self.PoleFlag=POLE
        self.PerformanceTestFlag=TEST
        self.EASYDYNPlusPlusFlag=EASEASYDYNPlusPlus
    # -------- end MBSClass.SetIntegrationParameters
    # ------ begin MBSClass.ComputeKinematics
    def ComputeKinematics(self):
        self.StartTime_Kinematics = time.clock() # Time initialisation
        # Some verifications for bodies and frames
        for ibody in range(0,self.nbrbody):
            if self.body[ibody].mainframe<0:
                if  self.body[ibody].T0F.is_zero:
                    print("Error in body "+str(ibody)+" : mainframe=-1 and T0F not defined.")
                    time.sleep(4)
                    exit("Error ")
                self.body[ibody].TtoF=self.body[ibody].T0F
            else:
                if  self.body[ibody].TrefF.is_zero:
                    print("Error in body "+str(ibody)+" : mainframe>=0 and TrefF not defined.")
                    time.sleep(4)
                    exit("Error ")
                self.body[ibody].TtoF=self.body[ibody].TrefF
            self.body[ibody].TtoF=mysimplify(self.body[ibody].TtoF)
            for iframe in range(0,self.body[ibody].nbrnodes):
                if self.body[ibody].frame[iframe].mainframe<0:
                    if  self.body[ibody].frame[iframe].T0F.is_zero:
                        print("Error in body "+str(ibody)+" frame "+str(iframe)+": mainframe=-1 and T0F not defined.")
                        time.sleep(4)
                        exit("Error ")
                    self.body[ibody].frame[iframe].TtoF=self.body[ibody].frame[iframe].T0F
                else:
                    if  self.body[ibody].frame[iframe].TrefF.is_zero:
                        print("Error in body "+str(ibody)+" frame "+str(iframe)+" : mainframe>=0 and TrefF not defined.")
                        time.sleep(4)
                        exit("Error ")
                    self.body[ibody].frame[iframe].TtoF=self.body[ibody].frame[iframe].TrefF
                self.body[ibody].TtoF=mysimplify(self.body[ibody].TtoF)
        # Calcul de l'expression des derivees des variables dependantes
        for idep in range(0,self.nbrdep):
            print('Beginning processing of dependent variable ',idep)
            self.pexpr[idep]=mysimplify(self.pexpr[idep])
            print('Simplification of dependent variable ',idep,' finished')
            for idof in range(0,self.nbrdof):
                temp=self.pexpr[idep].diff(self.q[idof])
                temp=mysimplify(temp)
                self.pdexpr[idep]=self.pdexpr[idep]+temp*self.qd[idof]
                self.pddexpr[idep]=self.pddexpr[idep]+temp*self.qdd[idof]
                for jdof in range(0,self.nbrdof):
                    temp2=temp.diff(self.q[jdof])
                    temp2=mysimplify(temp2)
                    self.pddexpr[idep]=self.pddexpr[idep]+temp2*self.qd[idof]*self.qd[jdof]
            print('Processing of dependent variable ',idep,' completed')
        # Calcul de la cinematique des corps rigides
        for ibody in range(0,self.nbrbody):
            print('Processing kinematics of Body [',ibody,']...')
            self.body[ibody].TtoF=mysimplify(self.body[ibody].TtoF)
            # calcul des vitesses partielles directes par rapport aux parametres de configuration
            RDFinal=Matrix([[0,0,0],[0,0,0],[0,0,0]])
            for idof in range(0,self.nbrdof):
                pos=self.body[ibody].TtoF[0:3,3]
                self.body[ibody].vFpartialq.append(pos.diff(self.q[idof]))
                self.body[ibody].vFpartialq[idof]=mysimplify(self.body[ibody].vFpartialq[idof])
                self.body[ibody].vFpartial.append(self.body[ibody].vFpartialq[idof])
                R=self.body[ibody].TtoF[0:3,0:3]
                RD=R.diff(self.q[idof])

                TempRD=RD*self.qd[idof]
                RDFinal = RDFinal +TempRD  # HN Derivative of the Rotation matrix for each body
                
                RDRT=RD*R.T
                self.body[ibody].omegaPartialq.append(Matrix([RDRT[2,1],RDRT[0,2],RDRT[1,0]]))
                self.body[ibody].omegaPartialq[idof]=mysimplify(self.body[ibody].omegaPartialq[idof])
                self.body[ibody].omegaPartial.append(self.body[ibody].omegaPartialq[idof])
            
            self.body[ibody].TtoFd[0:3,0:3]=RDFinal # HN Retrieving the derivative of the rotation matrix in the derivative trfo homog matrix
            
            # calcul des vitesses partielles directes par rapport aux parametres de configuration
            for idep in range(0,self.nbrdep):
                pos=self.body[ibody].TtoF[0:3,3]
                self.body[ibody].vFpartialp.append(pos.diff(self.p[idep]))
                self.body[ibody].vFpartialp[idep]=mysimplify(self.body[ibody].vFpartialp[idep])
                R=self.body[ibody].TtoF[0:3,0:3]
                RD=R.diff(self.p[idep])
                RDRT=RD*R.T
                self.body[ibody].omegaPartialp.append(Matrix([RDRT[2,1],RDRT[0,2],RDRT[1,0]]))
                self.body[ibody].omegaPartialp[idep]=mysimplify(self.body[ibody].omegaPartialp[idep])
            # ajout des contributions venant des variables dependantes dans les vitesses partielles totales
            for idof in range(0,self.nbrdof):
                for idep in range(0,self.nbrdep):
                    temp=self.pexpr[idep].diff(self.q[idof])
                    temp=mysimplify(temp)
                    self.body[ibody].vFpartial[idof]=self.body[ibody].vFpartial[idof] +\
                                                  self.body[ibody].vFpartialp[idep]*temp
                    self.body[ibody].omegaPartial[idof]=self.body[ibody].omegaPartial[idof] +\
                                                  self.body[ibody].omegaPartialp[idep]*temp
            # Calcul des termes de vitesse en qd
            for idof in range(0,self.nbrdof):
                self.body[ibody].vF=self.body[ibody].vF+self.body[ibody].vFpartialq[idof]*self.qd[idof]
                self.body[ibody].omega=self.body[ibody].omega+self.body[ibody].omegaPartialq[idof]*self.qd[idof]
            
            self.body[ibody].TtoFd[0:3,3] =self.body[ibody].vF # HN adding vF in the derivative trfo homog matrix
            
            # Calcul des termes de vitesse en pd
            for idep in range(0,self.nbrdep):
                self.body[ibody].vF=self.body[ibody].vF+self.body[ibody].vFpartialp[idep]*self.pd[idep]
                self.body[ibody].omega=self.body[ibody].omega+self.body[ibody].omegaPartialp[idep]*self.pd[idep]
            # Calcul des termes d'acceleration en qdd
            for idof in range(0,self.nbrdof):
                self.body[ibody].aF=self.body[ibody].aF+self.body[ibody].vFpartialq[idof]*self.qdd[idof]
                self.body[ibody].omegad=self.body[ibody].omegad+self.body[ibody].omegaPartialq[idof]*self.qdd[idof]
            # Calcul des termes d'acceleration en pdd
            for idep in range(0,self.nbrdep):
                self.body[ibody].aF=self.body[ibody].aF+self.body[ibody].vFpartialp[idep]*self.pdd[idep]
                self.body[ibody].omegad=self.body[ibody].omegad+self.body[ibody].omegaPartialp[idep]*self.pdd[idep]
            # Calcul des termes d'acceleration en qd*pd
            for idof in range(0,self.nbrdof):
                for idep in range(0,self.nbrdep):
                    temp=self.body[ibody].vFpartialq[idof].diff(self.p[idep]) \
                        +self.body[ibody].vFpartialp[idep].diff(self.q[idof])
                    temp=mysimplify(temp)
                    self.body[ibody].aF=self.body[ibody].aF+temp*self.qd[idof]*self.pd[idep]
                    temp=self.body[ibody].omegaPartialq[idof].diff(self.p[idep]) \
                        +self.body[ibody].omegaPartialp[idep].diff(self.q[idof])
                    temp=mysimplify(temp)
                    self.body[ibody].omegad=self.body[ibody].omegad+temp*self.qd[idof]*self.pd[idep]
            # Calcul des termes d'acceleration en qd*qd
            for idof in range(0,self.nbrdof):
                for jdof in range(0,self.nbrdof):
                    temp=self.body[ibody].vFpartialq[idof].diff(self.q[jdof])
                    temp=mysimplify(temp)
                    self.body[ibody].aF=self.body[ibody].aF+temp*self.qd[idof]*self.qd[jdof]
                    temp=self.body[ibody].omegaPartialq[idof].diff(self.q[jdof])
                    temp=mysimplify(temp)
                    self.body[ibody].omegad=self.body[ibody].omegad+temp*self.qd[idof]*self.qd[jdof]
            # Calcul des termes d'acceleration en pd*pd
            for idep in range(0,self.nbrdep):
                for jdep in range(0,self.nbrdep):
                    temp=self.body[ibody].vFpartialp[idep].diff(self.p[jdep])
                    temp=mysimplify(temp)
                    self.body[ibody].aF=self.body[ibody].aF+temp*self.pd[idep]*self.pd[jdep]
                    temp=self.body[ibody].omegaPartialp[idep].diff(self.p[jdep])
                    temp=mysimplify(temp)
                    self.body[ibody].omegad=self.body[ibody].omegad+temp*self.pd[idep]*self.pd[jdep]
            # Ajout des termes dependant directement du temps
            # for idof in range(0,self.nbrdof):
            # vitesses dependant du temps
            pos=self.body[ibody].TtoF[0:3,3]
            temp=pos.diff(self.t)
            temp=mysimplify(temp)
            self.body[ibody].vF=self.body[ibody].vF+temp
            R=self.body[ibody].TtoF[0:3,0:3]
            RD=R.diff(self.t)
            RDRT=RD*R.T
            temp=Matrix([RDRT[2,1],RDRT[0,2],RDRT[1,0]])
            temp=mysimplify(temp)
            self.body[ibody].omega=self.body[ibody].omega+temp
            # accelerations dependant du temps
            temp=self.body[ibody].vF.diff(self.t)
            temp=mysimplify(temp)
            self.body[ibody].aF=self.body[ibody].aF+temp
            temp=self.body[ibody].omega.diff(self.t)
            temp=mysimplify(temp)
            self.body[ibody].omegad=self.body[ibody].omegad+temp
            # Kinematics for each frame
            for iframe in range(0,self.body[ibody].nbrnodes):
                print('Processing kinematics of Body [',ibody,'] Frame [',iframe,']...')
                self.body[ibody].frame[iframe].TtoF=mysimplify(self.body[ibody].frame[iframe].TtoF)
                # calcul des vitesses partielles directes par rapport aux parametres de configuration
                RDFinalFrame=Matrix([[0,0,0],[0,0,0],[0,0,0]])
                for idof in range(0,self.nbrdof):
                    pos=self.body[ibody].frame[iframe].TtoF[0:3,3]
                    self.body[ibody].frame[iframe].vFpartialq.append(pos.diff(self.q[idof]))
                    self.body[ibody].frame[iframe].vFpartialq[idof]=mysimplify(self.body[ibody].frame[iframe].vFpartialq[idof])
                    self.body[ibody].frame[iframe].vFpartial.append(self.body[ibody].frame[iframe].vFpartialq[idof])
                    R=self.body[ibody].frame[iframe].TtoF[0:3,0:3]
                    RD=R.diff(self.q[idof])
                    
                    TempRDFrame=RD*self.qd[idof]
                    RDFinalFrame=RDFinalFrame +TempRDFrame # HN Derivative of the rotation matrix for each frame
                    
                    RDRT=RD*R.T
                    self.body[ibody].frame[iframe].omegaPartialq.append(Matrix([RDRT[2,1],RDRT[0,2],RDRT[1,0]]))
                    self.body[ibody].frame[iframe].omegaPartialq[idof]=mysimplify(self.body[ibody].frame[iframe].omegaPartialq[idof])
                    self.body[ibody].frame[iframe].omegaPartial.append(self.body[ibody].frame[iframe].omegaPartialq[idof])
                    
                self.body[ibody].frame[iframe].TtoFd[0:3,0:3]=RDFinalFrame # HN Retrieving the derivative of the rotation matrix for each frame    
                    
                # calcul des vitesses partielles directes par rapport aux parametres de configuration
                for idep in range(0,self.nbrdep):
                    pos=self.body[ibody].frame[iframe].TtoF[0:3,3]
                    self.body[ibody].frame[iframe].vFpartialp.append(pos.diff(self.p[idep]))
                    self.body[ibody].frame[iframe].vFpartialp[idep]=mysimplify(self.body[ibody].frame[iframe].vFpartialp[idep])
                    R=self.body[ibody].frame[iframe].TtoF[0:3,0:3]
                    RD=R.diff(self.p[idep])
                    RDRT=RD*R.T
                    self.body[ibody].frame[iframe].omegaPartialp.append(Matrix([RDRT[2,1],RDRT[0,2],RDRT[1,0]]))
                    self.body[ibody].frame[iframe].omegaPartialp[idep]=mysimplify(self.body[ibody].frame[iframe].omegaPartialp[idep])
                # ajout des contributions venant des variables dependantes dans les vitesses partielles totales
                for idof in range(0,self.nbrdof):
                    for idep in range(0,self.nbrdep):
                        temp=self.pexpr[idep].diff(self.q[idof])
                        temp=mysimplify(temp)
                        self.body[ibody].frame[iframe].vFpartial[idof]=self.body[ibody].frame[iframe].vFpartial[idof] +\
                                                      self.body[ibody].frame[iframe].vFpartialp[idep]*temp
                        self.body[ibody].frame[iframe].omegaPartial[idof]=self.body[ibody].frame[iframe].omegaPartial[idof] +\
                                                      self.body[ibody].frame[iframe].omegaPartialp[idep]*temp
                # Calcul des termes de vitesse en qd
                for idof in range(0,self.nbrdof):
                    self.body[ibody].frame[iframe].vF=self.body[ibody].frame[iframe].vF+self.body[ibody].frame[iframe].vFpartialq[idof]*self.qd[idof]
                    self.body[ibody].frame[iframe].omega=self.body[ibody].frame[iframe].omega+self.body[ibody].frame[iframe].omegaPartialq[idof]*self.qd[idof]
                    
                    self.body[ibody].frame[iframe].TtoFd[0:3,3] =self.body[ibody].frame[iframe].vF # HN adding vF in the derivative trfo homog matrix
                    
                # Calcul des termes de vitesse en pd
                for idep in range(0,self.nbrdep):
                    self.body[ibody].frame[iframe].vF=self.body[ibody].frame[iframe].vF+self.body[ibody].frame[iframe].vFpartialp[idep]*self.pd[idep]
                    self.body[ibody].frame[iframe].omega=self.body[ibody].frame[iframe].omega+self.body[ibody].frame[iframe].omegaPartialp[idep]*self.pd[idep]
                # Calcul des termes d'acceleration en qdd
                for idof in range(0,self.nbrdof):
                    self.body[ibody].frame[iframe].aF=self.body[ibody].frame[iframe].aF+self.body[ibody].frame[iframe].vFpartialq[idof]*self.qdd[idof]
                    self.body[ibody].frame[iframe].omegad=self.body[ibody].frame[iframe].omegad+self.body[ibody].frame[iframe].omegaPartialq[idof]*self.qdd[idof]
                # Calcul des termes d'acceleration en pdd
                for idep in range(0,self.nbrdep):
                    self.body[ibody].frame[iframe].aF=self.body[ibody].frame[iframe].aF+self.body[ibody].frame[iframe].vFpartialp[idep]*self.pdd[idep]
                    self.body[ibody].frame[iframe].omegad=self.body[ibody].frame[iframe].omegad+self.body[ibody].frame[iframe].omegaPartialp[idep]*self.pdd[idep]
                # Calcul des termes d'acceleration en qd*pd
                for idof in range(0,self.nbrdof):
                    for idep in range(0,self.nbrdep):
                        temp=self.body[ibody].frame[iframe].vFpartialq[idof].diff(self.p[idep]) \
                            +self.body[ibody].frame[iframe].vFpartialp[idep].diff(self.q[idof])
                        temp=mysimplify(temp)
                        self.body[ibody].frame[iframe].aF=self.body[ibody].frame[iframe].aF+temp*self.qd[idof]*self.pd[idep]
                        temp=self.body[ibody].frame[iframe].omegaPartialq[idof].diff(self.p[idep]) \
                            +self.body[ibody].frame[iframe].omegaPartialp[idep].diff(self.q[idof])
                        temp=mysimplify(temp)
                        self.body[ibody].frame[iframe].omegad=self.body[ibody].frame[iframe].omegad+temp*self.qd[idof]*self.pd[idep]
                # Calcul des termes d'acceleration en qd*qd
                for idof in range(0,self.nbrdof):
                    for jdof in range(0,self.nbrdof):
                        temp=self.body[ibody].frame[iframe].vFpartialq[idof].diff(self.q[jdof])
                        temp=mysimplify(temp)
                        self.body[ibody].frame[iframe].aF=self.body[ibody].frame[iframe].aF+temp*self.qd[idof]*self.qd[jdof]
                        temp=self.body[ibody].frame[iframe].omegaPartialq[idof].diff(self.q[jdof])
                        temp=mysimplify(temp)
                        self.body[ibody].frame[iframe].omegad=self.body[ibody].frame[iframe].omegad+temp*self.qd[idof]*self.qd[jdof]
                # Calcul des termes d'acceleration en pd*pd
                for idep in range(0,self.nbrdep):
                    for jdep in range(0,self.nbrdep):
                        temp=self.body[ibody].frame[iframe].vFpartialp[idep].diff(self.p[jdep])
                        temp=mysimplify(temp)
                        self.body[ibody].frame[iframe].aF=self.body[ibody].frame[iframe].aF+temp*self.pd[idep]*self.pd[jdep]
                        temp=self.body[ibody].frame[iframe].omegaPartialp[idep].diff(self.p[jdep])
                        temp=mysimplify(temp)
                        self.body[ibody].frame[iframe].omegad=self.body[ibody].frame[iframe].omegad+temp*self.pd[idep]*self.pd[jdep]
                # Ajout des termes dependant directement du temps
                for idof in range(0,self.nbrdof):
                    # vitesses dependant du temps
                    pos=self.body[ibody].frame[iframe].TtoF[0:3,3]
                    temp=pos.diff(self.t)
                    temp=mysimplify(temp)
                    self.body[ibody].frame[iframe].vF=self.body[ibody].frame[iframe].vF+temp
                    R=self.body[ibody].frame[iframe].TtoF[0:3,0:3]
                    RD=R.diff(self.t)
                    RDRT=RD*R.T
                    temp=Matrix([RDRT[2,1],RDRT[0,2],RDRT[1,0]])
                    temp=mysimplify(temp)
                    self.body[ibody].frame[iframe].omega=self.body[ibody].frame[iframe].omega+temp
                    # accelerations dependant du temps
                    temp=self.body[ibody].frame[iframe].vF.diff(self.t)
                    temp=mysimplify(temp)
                    self.body[ibody].frame[iframe].aF=self.body[ibody].frame[iframe].aF+temp
                    temp=self.body[ibody].frame[iframe].omega.diff(self.t)
                    temp=mysimplify(temp)
                    self.body[ibody].frame[iframe].omegad=self.body[ibody].frame[iframe].omegad+temp    
        self.EndTime_Kinematics = round(time.clock(),1)
    # -------- end MBSClass.ComputeKinematics
    # ------ begin MBSClass.PrintKinematics
    def PrintKinematics(self):
        for ibody in range(0,self.nbrbody):
            if (self.body[ibody].mainframe<0):
                print('Body[',ibody,'].T0F=',self.body[ibody].TtoF)
                print('Body[',ibody,'].vF=',self.body[ibody].vF)
                print('Body[',ibody,'].omega=',self.body[ibody].omega)
                print('Body[',ibody,'].aF=',self.body[ibody].aF)
                print('Body[',ibody,'].omegad=',self.body[ibody].omegad)
                for idof in range(0,self.nbrdof):
                    print('Body[',ibody,'].vFpartialq[',idof,']=',self.body[ibody].vFpartialq[idof])
                    print('Body[',ibody,'].omegaPartialq[',idof,']=',self.body[ibody].omegaPartialq[idof])
                for idep in range(0,self.nbrdep):
                    print('Body[',ibody,'].vFpartialp[',idep,']=',self.body[ibody].vFpartialp[idep])
                    print('Body[',ibody,'].omegaPartialp[',idep,']=',self.body[ibody].omegaPartialp[idep])
                for idof in range(0,self.nbrdof):
                    print('Body[',ibody,'].vFpartial[',idof,']=',self.body[ibody].vFpartial[idof])
                    print('Body[',ibody,'].omegaPartial[',idof,']=',self.body[ibody].omegaPartial[idof])
            if (self.body[ibody].mainframe>=0):
                print('Body[',ibody,'].BodyRef=',self.body[ibody].mainframe)
                print('Body[',ibody,'].TrefF=',self.body[ibody].TtoF)
                print('Body[',ibody,'].vFrel=',self.body[ibody].vF)
                print('Body[',ibody,'].omegarel=',self.body[ibody].omega)
                print('Body[',ibody,'].aFrel=',self.body[ibody].aF)
                print('Body[',ibody,'].omegadrel=',self.body[ibody].omegad)
                for idof in range(0,self.nbrdof):
                    print('Body[',ibody,'].vFrelPartialq[',idof,']=',self.body[ibody].vFpartialq[idof])
                    print('Body[',ibody,'].omegarelPartialq[',idof,']=',self.body[ibody].omegaPartialq[idof])
                for idep in range(0,self.nbrdep):
                    print('Body[',ibody,'].vFrelPartialp[',idep,']=',self.body[ibody].vFpartialp[idep])
                    print('Body[',ibody,'].omegarelPartialp[',idep,']=',self.body[ibody].omegaPartialp[idep])
                for idof in range(0,self.nbrdof):
                    print('Body[',ibody,'].vFrelPartial[',idof,']=',self.body[ibody].vFpartial[idof])
                    print('Body[',ibody,'].omegarelPartial[',idof,']=',self.body[ibody].omegaPartial[idof])
            for iframe in range(0,self.body[ibody].nbrnodes):
                if (self.body[ibody].frame[iframe].mainframe<0):
                    print('Body[',ibody,'].NodeFrame[',iframe,'].T0F=',self.body[ibody].frame[iframe].TtoF)
                    print('Body[',ibody,'].NodeFrame[',iframe,'].vF=',self.body[ibody].frame[iframe].vF)
                    print('Body[',ibody,'].NodeFrame[',iframe,'].omega=',self.body[ibody].frame[iframe].omega)
                    print('Body[',ibody,'].NodeFrame[',iframe,'].aF=',self.body[ibody].frame[iframe].aF)
                    print('Body[',ibody,'].NodeFrame[',iframe,'].omegad=',self.body[ibody].frame[iframe].omegad)
                    for idof in range(0,self.nbrdof):
                        print('Body[',ibody,'].NodeFrame[',iframe,'].vFpartialq[',idof,']=',self.body[ibody].frame[iframe].vFpartialq[idof])
                        print('Body[',ibody,'].NodeFrame[',iframe,'].omegaPartialq[',idof,']=',self.body[ibody].frame[iframe].omegaPartialq[idof])
                    for idep in range(0,self.nbrdep):
                        print('Body[',ibody,'].NodeFrame[',iframe,'].vFpartialp[',idep,']=',self.body[ibody].frame[iframe].vFpartialp[idep])
                        print('Body[',ibody,'].NodeFrame[',iframe,'].omegaPartialp[',idep,']=',self.body[ibody].frame[iframe].omegaPartialp[idep])
                    for idof in range(0,self.nbrdof):
                        print('Body[',ibody,'].NodeFrame[',iframe,'].vFpartial[',idof,']=',self.body[ibody].frame[iframe].vFpartial[idof])
                        print('Body[',ibody,'].NodeFrame[',iframe,'].omegaPartial[',idof,']=',self.body[ibody].frame[iframe].omegaPartial[idof])
                if (self.body[ibody].frame[iframe].mainframe>=0):
                    print('Body[',ibody,'].NodeFrame[',iframe,'].TrefF=',self.body[ibody].frame[iframe].TtoF)
                    print('Body[',ibody,'].NodeFrame[',iframe,'].vFrel=',self.body[ibody].frame[iframe].vF)
                    print('Body[',ibody,'].NodeFrame[',iframe,'].omegarel=',self.body[ibody].frame[iframe].omega)
                    print('Body[',ibody,'].NodeFrame[',iframe,'].aFrel=',self.body[ibody].frame[iframe].aF)
                    print('Body[',ibody,'].NodeFrame[',iframe,'].omegadrel=',self.body[ibody].frame[iframe].omegad)
                    for idof in range(0,self.nbrdof):
                        print('Body[',ibody,'].NodeFrame[',iframe,'].vFrelPartialq[',idof,']=',self.body[ibody].frame[iframe].vFpartialq[idof])
                        print('Body[',ibody,'].NodeFrame[',iframe,'].omegarelPartialq[',idof,']=',self.body[ibody].frame[iframe].omegaPartialq[idof])
                    for idep in range(0,self.nbrdep):
                        print('Body[',ibody,'].NodeFrame[',iframe,'].vFrelPartialp[',idep,']=',self.body[ibody].frame[iframe].vFpartialp[idep])
                        print('Body[',ibody,'].NodeFrame[',iframe,'].omegarelPartialp[',idep,']=',self.body[ibody].frame[iframe].omegaPartialp[idep])
                    for idof in range(0,self.nbrdof):
                        print('Body[',ibody,'].NodeFrame[',iframe,'].vFrelPartial[',idof,']=',self.body[ibody].frame[iframe].vFpartial[idof])
                        print('Body[',ibody,'].NodeFrame[',iframe,'].omegarelPartial[',idof,']=',self.body[ibody].frame[iframe].omegaPartial[idof])
    # -------- end MBSClass.ComputeKinematic
    # ------ begin MBSClass.ExportEasyDynHeader
    def ExportEasyDynHeader(self,OutFile):
        OutFile.write("// File in C++ format\n")
        OutFile.write("// Application title: "+self.ApplicationTitle+"\n")
        OutFile.write("//\n")
        OutFile.write("//    copyright (C) 2009 Olivier VERLINDEN\n")
        OutFile.write("//    Service de Mecanique rationnelle, Dynamique et Vibrations\n")
        OutFile.write("//    Faculte Polytechnique de Mons\n")
        OutFile.write("//    31, Bd Dolez, 7000 MONS (Belgium)\n")
        OutFile.write("//    Olivier.Verlinden@fpms.ac.be\n")
        OutFile.write("//\n")
        OutFile.write("// This file is part of EasyDyn\n")
        OutFile.write("//\n")
        OutFile.write("// EasyDyn is free software; you can redistribute it and/or modify it under the\n")
        OutFile.write("// terms of the GNU GeneralPublic License as published by the Free Software\n")
        OutFile.write("// Foundation; either version 2, or (at your option) any later version.\n")
        OutFile.write("//\n")
        OutFile.write("// EasyDyn is distributed in the hope that it will be useful, but WITHOUT ANY\n")
        OutFile.write("// WARRANTY; without the implied warranty of MERCHANTABILITY or FITNESS FOR\n")
        OutFile.write("// A PARTICULAR PURPOSE. See the GNU General Public License for more details.\n")
        OutFile.write("//\n")
        OutFile.write("// You should have received a copy of the GNU General Public License along with\n")
        OutFile.write("// EasyDyn; see the file COPYING.  If not, write to the Free Software\n")
        OutFile.write("// Foundation, 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.\n\n")
        OutFile.write("//\n\n")
        OutFile.write("#define EASYDYNMBSMAIN  // to declare the global variable\n")
        OutFile.write("#define EASYDYNMBSADVANCED // This option implies that a specific \n")
        OutFile.write("                           // version of ComputePartialVelocities() is provided\n")
        OutFile.write("#include <stdio.h>\n")
        OutFile.write("#include <math.h>\n")
        OutFile.write("#include <EasyDyn/mbs.h>\n")
        OutFile.write("#include <EasyDyn/visu.h>\n")
        OutFile.write("#include <EasyDyn/vec.h>\n")
        OutFile.write("#include <fstream>\n\n")
        OutFile.write("using namespace std;\n")
        OutFile.write("scene thescene;\n")
        OutFile.write("ofstream VanFile;\n\n")
        if(self.GlobalVar != " "):
            OutFile.write("// Global variables\n")
            OutFile.write(self.GlobalVar)
            OutFile.write("\n")
        OutFile.write("//---------------------------------------------------\n\n")
        OutFile.write("void WriteDataHeader(ostream &OutFile)\n")
        OutFile.write("{\n")
        OutFile.write("WriteStateVariablesHeader(OutFile);\n")
        OutFile.write("OutFile <<")
        for ibody in range(0,self.nbrbody):
            BodyName=" \"body"+str(ibody)+"_"
            self.body[ibody].ExportExtraHeader(OutFile,BodyName)
        OutFile.write(" endl;\n")
        OutFile.write("}\n\n")
        OutFile.write("//----------------------------------------------------\n\n")
        OutFile.write("void SaveData(ostream &OutFile)\n")
        OutFile.write("{\n")
        
        # Define variables in SaveData procedure
        if(self.SaveDataVar!=" "):
            OutFile.write("\n")
            OutFile.write("  "+self.SaveDataVar+";\n")
            OutFile.write("\n")    
        OutFile.write("SaveStateVariables(OutFile);\n")
        OutFile.write("OutFile <<")
        for ibody in range(0,self.nbrbody):
            BodyName="   body["+str(ibody)+"]"
            self.body[ibody].ExportExtraRes(OutFile,BodyName)
        OutFile.write(" endl;\n")
        OutFile.write("thescene.WriteCoord(VanFile);\n")
        OutFile.write("}\n")
        
        
        
        OutFile.write("\n//----------------------------------------------------\n")
    # ------ begin MBSClass.ExportEasyDynHeader
    # ------ begin MBSClass.ExportEasyDynPlusPlusHeader
    def ExportEasyDynPlusPlusHeader(self,OutFile):
        OutFile.write("// File in C++ format\n")
        OutFile.write("// Application title: "+self.ApplicationTitle+"\n")
        OutFile.write("//\n")
        OutFile.write("//    copyright (C) 2009 Olivier VERLINDEN\n")
        OutFile.write("//    Service de Mecanique rationnelle, Dynamique et Vibrations\n")
        OutFile.write("//    Faculte Polytechnique de Mons\n")
        OutFile.write("//    31, Bd Dolez, 7000 MONS (Belgium)\n")
        OutFile.write("//    Olivier.Verlinden@fpms.ac.be\n")
        OutFile.write("//\n")
        OutFile.write("// This file is part of EasyDyn\n")
        OutFile.write("//\n")
        OutFile.write("// EasyDyn is free software; you can redistribute it and/or modify it under the\n")
        OutFile.write("// terms of the GNU GeneralPublic License as published by the Free Software\n")
        OutFile.write("// Foundation; either version 2, or (at your option) any later version.\n")
        OutFile.write("//\n")
        OutFile.write("// EasyDyn is distributed in the hope that it will be useful, but WITHOUT ANY\n")
        OutFile.write("// WARRANTY; without the implied warranty of MERCHANTABILITY or FITNESS FOR\n")
        OutFile.write("// A PARTICULAR PURPOSE. See the GNU General Public License for more details.\n")
        OutFile.write("//\n")
        OutFile.write("// You should have received a copy of the GNU General Public License along with\n")
        OutFile.write("// EasyDyn; see the file COPYING.  If not, write to the Free Software\n")
        OutFile.write("// Foundation, 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.\n\n")
        OutFile.write("#include <stdio.h>\n")
        OutFile.write("#include <math.h>\n")
        OutFile.write("#include <EasyDyn++/mbs.h>\n")
        OutFile.write("#include <EasyDyn++/visu.h>\n")
        OutFile.write("#include <EasyDyn++/vec.h>\n")
        OutFile.write("#include <fstream>\n\n")
        OutFile.write("#define ADVANCED // This option implies that a specific version of \n")
        OutFile.write("                 // ComputePartialVelocities() is provided \n\n")
        OutFile.write("using namespace std;\n")
        OutFile.write("scene thescene;\n")
        OutFile.write("ofstream VanFile;\n\n")
        if(self.GlobalVar != " "):
            OutFile.write("// Global variables\n")
            OutFile.write(self.GlobalVar)
            OutFile.write("\n")
        OutFile.write("//---------------------------------------------------\n\n")
        OutFile.write("//System class: "+self.ApplicationFileName+"\n")
        OutFile.write("class "+self.ApplicationFileName+": public SysMBS\n")
        OutFile.write("  {\n")
        OutFile.write("  public:\n")
        OutFile.write('    '+self.ApplicationFileName+'(): SysMBS("'+self.ApplicationFileName+'",'+str(self.nbrbody)+'/*nbrbody*/,')
        if(self.OverWriteNbrdof>0):
            OutFile.write(str(self.OverWriteNbrdof)+'/*nbrdof*/,')
        else:
            OutFile.write(str(self.nbrdof)+'/*nbrdof*/,')
        OutFile.write(str(self.nbrdep)+'/*nbrdep*/')
        if(self.nbrinput>0):
            OutFile.write(','+str(self.nbrinput)+'/*nbrinput*/')
        if(self.nbroutput>0):
            OutFile.write(','+str(self.nbroutput)+'/*nbroutput*/')
        if(self.nbrexternalinput>0):
            OutFile.write(','+str(self.nbrexternalinput)+'/*nbrexternalinput*/')
        OutFile.write(')\n')
        OutFile.write("    {\n")
        for ibody in range(0,self.nbrbody):
            BodyName="      body["+str(ibody)+"] "
            # If rigid body
            if(self.body[ibody].mass>0 and self.body[ibody].Ixx>0 and self.body[ibody].Iyy>0 and self.body[ibody].Izz>0 and self.body[ibody].Ixy>=0 and self.body[ibody].Ixz>=0 and self.body[ibody].Iyz>=0):
                OutFile.write(BodyName+"= new RigidBody("+str(self.body[ibody].mass)+"/*mass*/,"+str(self.body[ibody].Ixx)+"/*Ixx*/,"+str(self.body[ibody].Iyy)+"/*Iyy*/,"+str(self.body[ibody].Izz)+"/*Izz*/,"+str(self.body[ibody].Ixy)+"/*Ixy*/,"+str(self.body[ibody].Ixz)+"/*Ixz*/,"+str(self.body[ibody].Iyz)+"/*Iyz*/,"+str(self.nbrdof)+"/*nbrdof*/,"+str(self.body[ibody].nbrnodes)+"/*nbrnodes*/);\n")
                self.nbrRigidBody=self.nbrRigidBody+1 # Count the number of rigid bodies
            # if flexible body
            elif(self.body[ibody].name!=0 and self.body[ibody].nbrnodes>=2 and self.body[ibody].AlphaDamp>=0 and self.body[ibody].BetaDamp>=0 and self.body[ibody].length==0):
                OutFile.write(BodyName+'= new FlexibleBody("'+str(self.body[ibody].name)+'"/*FileName*/,')
                OutFile.write(str(self.nbrdof)+"/*nbrdof*/,"+str(self.body[ibody].nbrnodes)+"/*nbrnodes*/,")
                OutFile.write(str(self.body[ibody].AlphaDamp)+"/*AlphaDamp*/,"+str(self.body[ibody].BetaDamp)+"/*BetaDamp*/);\n")
                self.nbrFlexBody=self.nbrFlexBody+1 # Count the number of flexible bodies
            # if flexible beam
            elif(self.body[ibody].length>0 and self.body[ibody].section>0 and self.body[ibody].EYoung>0 and self.body[ibody].nu>0 and self.body[ibody].rho>0 and self.body[ibody].Iyy>0 and self.body[ibody].Izz>0 and self.body[ibody].AlphaDamp>=0 and self.body[ibody].BetaDamp>=0 and self.body[ibody].nbrnodes>=2):
                OutFile.write(BodyName+"= new FlexibleBeam("+str(self.body[ibody].length)+"/*length*/,"+str(self.body[ibody].section)+"/*section*/,")
                OutFile.write(str(self.body[ibody].EYoung)+"/*EYoung*/,"+str(self.body[ibody].nu)+"/*nu*/,"+str(self.body[ibody].rho)+"/*rho*/,")
                OutFile.write(str(self.body[ibody].Iyy)+"/*IGeomyy*/,"+str(self.body[ibody].Izz)+"/*IGeomzz*/,")
                OutFile.write(str(self.nbrdof)+"/*nbrdof*/,"+str(self.body[ibody].nbrnodes)+"/*nbrnodes*/,")
                OutFile.write(str(self.body[ibody].AlphaDamp)+"/*AlphaDamp*/,"+str(self.body[ibody].BetaDamp)+"/*BetaDamp*/);\n")  
                self.nbrFlexBeam=self.nbrFlexBeam+1 # Count the number of flexible beams
            else: 
                print(" ")
                print("Error in inertia properties of body ["+str(ibody)+"]")
                time.sleep(4)
                exit("Error ")
        for ibody in range(0,self.nbrbody):
            #if it is a flexible body 
            if((self.body[ibody].length>0 and self.body[ibody].section>0 and self.body[ibody].EYoung>0) or self.body[ibody].name!=0): 
        #       OutFile.write("      //Rigid disp of nodes wrt corotational frame\n")
                for iframe in range(0,self.body[ibody].nbrnodes):
                    CoroName="      body["+str(ibody)+"]"+"->NodeFrame["+str(iframe)+"]"
                    self.body[ibody].frame[iframe].ExportCoroDispData(OutFile,CoroName)
        #OutFile.write("      //Relative motion(s)\n")
        for ibody in range(0,self.nbrbody):
            if(self.body[ibody].nbrnodes==0):
                if (self.body[ibody].mainframe>=0 and self.body[ibody].frameref==-1):
                    OutFile.write("      body["+str(ibody)+"]->MainFrame.RelMotionRefFrame=&(body["\
                                  +str(self.body[ibody].mainframe)+"]->MainFrame);\n")
                elif(self.body[ibody].mainframe>=0 and self.body[ibody].frameref>=0):
                    OutFile.write("      body["+str(ibody)+"]->MainFrame.RelMotionRefFrame=&(body["\
                                  +str(self.body[ibody].mainframe)+"]->NodeFrame["+str(self.body[ibody].frameref)+"]);\n")
            else:
                if (self.body[ibody].mainframe>=0 and self.body[ibody].frameref==-1):
                    OutFile.write("      body["+str(ibody)+"]->MainFrame.RelMotionRefFrame=&(body["\
                                  +str(self.body[ibody].mainframe)+"]->MainFrame);\n")
                elif(self.body[ibody].mainframe>=0 and self.body[ibody].frameref>=0):
                    OutFile.write("      body["+str(ibody)+"]->MainFrame.RelMotionRefFrame=&(body["\
                                  +str(self.body[ibody].mainframe)+"]->NodeFrame["+str(self.body[ibody].frameref)+"]);\n")
                for iframe in range(0,self.body[ibody].nbrnodes):
                    if (self.body[ibody].frame[iframe].mainframe>=0 and self.body[ibody].frame[iframe].frameref==-1):
                        OutFile.write("      body["+str(ibody)+"]->NodeFrame["+str(iframe)+"].RelMotionRefFrame=&(body["\
                                      +str(self.body[ibody].frame[iframe].mainframe)+"]->MainFrame);\n")
                    elif(self.body[ibody].frame[iframe].mainframe>=0 and self.body[ibody].frame[iframe].frameref>=0):
                        OutFile.write("      body["+str(ibody)+"]->NodeFrame["+str(iframe)+"].RelMotionRefFrame=&(body["\
                                      +str(self.body[ibody].frame[iframe].mainframe)+"]->NodeFrame["+str(self.body[ibody].frame[iframe].frameref)+"]);\n")
        OutFile.write("      gravity.put("+ccode(self.gravity[0])+","+ccode(self.gravity[1])\
                           +","+ccode(self.gravity[2])+");\n")
        OutFile.write("    };\n")
        OutFile.write("    void SaveData(ostream &OutFile);\n")
        OutFile.write("    void ComputeMotion();\n")
        OutFile.write("    #ifdef ADVANCED\n")
        OutFile.write("    void ComputePartialVelocities();\n")
        OutFile.write("    #endif\n")
        OutFile.write("    void AddAppliedForces();\n")
        OutFile.write("\n")
        if(self.ComputeResidualVar != " "):
            OutFile.write("    void ComputeResidual();\n")
        if(self.nbrinput>0):
            OutFile.write("    void ComputeInputs();\n")
        if(self.nbroutput>0):
            OutFile.write("    void ComputeOutputs();\n")
        if(self.nbrexternalinput>0):
            OutFile.write("    void DefineExternalInputs();\n")
        OutFile.write("  };\n\n")
        #OutFile.write("void WriteDataHeader(ostream &OutFile)\n")
        #OutFile.write("{\n")
        #OutFile.write("WriteStateVariablesHeader(OutFile);\n")
        #OutFile.write("OutFile <<")
        #for ibody in range(0,self.nbrbody):
        #    BodyName=" \"body"+str(ibody)+"_"
        #    self.body[ibody].ExportExtraHeader(OutFile,BodyName)
        #    for iframe in range(0,self.body[ibody].nbrnodes):
        #        BodyName=" \"body"+str(ibody)+"_f"+str(iframe)+"_"
        #        self.body[ibody].frame[iframe].ExportExtraHeader(OutFile,BodyName)
        #OutFile.write(" endl;\n")
        #OutFile.write("}\n\n")
        OutFile.write("//----------------------------------------------------\n\n")
        OutFile.write("void "+self.ApplicationFileName+"::SaveData(ostream &OutFile)\n")
        OutFile.write("{\n")
        
        # Define variables in SaveData procedure
        if(self.SaveDataVar!=" "):
            OutFile.write("\n")
            OutFile.write("  "+self.SaveDataVar+";\n")
            OutFile.write("\n")
            
        OutFile.write("SaveStateVariables(OutFile);\n")
        OutFile.write("OutFile <<")
        for ibody in range(0,self.nbrbody):
            BodyName="   body["+str(ibody)+"]"
            self.body[ibody].ExportExtraResPlusPlus(OutFile,BodyName)
            for iframe in range(0,self.body[ibody].nbrnodes):
                BodyName=" body["+str(ibody)+"]->NodeFrame["+str(iframe)+"]"
                self.body[ibody].frame[iframe].ExportExtraRes(OutFile,BodyName)
        OutFile.write(" endl;\n")
        OutFile.write("thescene.WriteCoord(VanFile);\n")
        OutFile.write("}\n")
        OutFile.write("\n//----------------------------------------------------\n")
    # ------ begin MBSClass.ExportEasyDynPlusPlusHeader
    # ------ begin MBSClass.ExportSetInertiaData
    def ExportSetInertiaData(self,OutFile):
        OutFile.write("\nvoid SetInertiaData()\n{\n")
        for ibody in range(0,self.nbrbody):
            BodyName="body["+str(ibody)+"]"
            self.body[ibody].ExportInertiaData(OutFile,BodyName)
        OutFile.write("}\n\n//-----------------------------------------\n")
    # -------- end MBSClass.ExportSetInertiadata
    # ------ begin MBSClass.ExportComputeMotion
    def ExportComputeMotion(self,OutFile): 
        OutFile.write("\nvoid ComputeMotion()\n{\n")
        OutFile.write("//Homogenous transformation matrices of each body\n")
        OutFile.write("//Insert kinematics generated by python/sympy\n\n")
        for idep in range(0,self.nbrdep):
            VarName="p["+str(idep)+"]"
            exportExpressionToCProgram(OutFile,VarName,self.pexpr[idep])
            VarName="pd["+str(idep)+"]"
            exportExpressionToCProgram(OutFile,VarName,self.pdexpr[idep])
            VarName="pdd["+str(idep)+"]"
            exportExpressionToCProgram(OutFile,VarName,self.pddexpr[idep])
        for ibody in range(0,self.nbrbody):
            OutFile.write("//Body["+str(ibody)+"]\n")
            BodyName="body["+str(ibody)+"]"
            self.body[ibody].ExportKinematics(OutFile,BodyName)
            if(self.body[ibody].mainframe>=0):
                OutFile.write("ComposeMotion("+str(ibody)+","+str(self.body[ibody].mainframe)+");\n")
        OutFile.write("}\n\n//-----------------------------------------\n") 
    # -------- end MBSClass.ExportComputeMotion
    # ------ begin MBSClass.ExportComputeMotionPlusPlus
    def ExportComputeMotionPlusPlus(self,OutFile): 
        OutFile.write("\nvoid "+self.ApplicationFileName+"::ComputeMotion()\n{\n")
        OutFile.write("//Homogenous transformation matrices of each body\n")
        OutFile.write("//Insert kinematics generated by python/sympy\n\n")
        for idep in range(0,self.nbrdep):
            VarName="p["+str(idep)+"]"
            exportExpressionToCProgram(OutFile,VarName,self.pexpr[idep])
            VarName="pd["+str(idep)+"]"
            exportExpressionToCProgram(OutFile,VarName,self.pdexpr[idep])
            VarName="pdd["+str(idep)+"]"
            exportExpressionToCProgram(OutFile,VarName,self.pddexpr[idep])
        for ibody in range(0,self.nbrbody):
            if(self.body[ibody].nbrnodes==0):
                OutFile.write("//Body["+str(ibody)+"]\n")
                BodyName="body["+str(ibody)+"]->MainFrame"
                self.body[ibody].ExportKinematicsPlusPlus(OutFile,BodyName)
                if(self.body[ibody].mainframe>=0):
                    OutFile.write("body["+str(ibody)+"]->MainFrame.ComposeMotion();\n") 
            else:
                OutFile.write("//Body["+str(ibody)+"]\n")
                BodyName="body["+str(ibody)+"]->MainFrame"
                self.body[ibody].ExportKinematicsPlusPlus(OutFile,BodyName)
                if(self.body[ibody].mainframe>=0):
                    OutFile.write("body["+str(ibody)+"]->MainFrame.ComposeMotion();\n") 
                for iframe in range(0,self.body[ibody].nbrnodes):
                    OutFile.write("//Body["+str(ibody)+"] Node["+str(iframe)+"]\n")
                    BodyName="body["+str(ibody)+"]->NodeFrame["+str(iframe)+"]"
                    self.body[ibody].frame[iframe].ExportKinematics(OutFile,BodyName)
                    if(self.body[ibody].frame[iframe].mainframe>=0):
                        OutFile.write("body["+str(ibody)+"]->NodeFrame["+str(iframe)+"].ComposeMotion();\n") 
        OutFile.write("}\n\n//-----------------------------------------\n")
    # -------- end MBSClass.ExportComputeMotionPlusPlus
    # ------ begin MBSClass.ExportComputePartialVelocities
    def ExportComputePartialVelocities(self,OutFile):
        OutFile.write("\nvoid ComputePartialVelocities()\n{\n")
        for ibody in range(0,self.nbrbody):
            OutFile.write("//Body["+str(ibody)+"]\n")
            BodyName="body["+str(ibody)+"]"
            self.body[ibody].ExportPartialVelocities(OutFile,BodyName)
            if(self.body[ibody].mainframe>=0):
                OutFile.write("ComposePartialVelocities("+str(ibody)+","+str(self.body[ibody].mainframe)+");\n")
        OutFile.write("}\n\n//-----------------------------------------\n") 
    # ------ end MBSClass.ExportComputePartialVelocities
    # ------ begin MBSClass.ExportComputePartialVelocitiesPlusPlus
    def ExportComputePartialVelocitiesPlusPlus(self,OutFile):
        OutFile.write("#ifdef ADVANCED\n")
        OutFile.write("void "+self.ApplicationFileName+"::ComputePartialVelocities()\n{\n")
        for ibody in range(0,self.nbrbody):
            if(self.body[ibody].nbrnodes==0):
                OutFile.write("//Body["+str(ibody)+"]\n")
                BodyName="body["+str(ibody)+"]->MainFrame"
                self.body[ibody].ExportPartialVelocitiesPlusPlus(OutFile,BodyName)
                if(self.body[ibody].mainframe>=0):
                    OutFile.write("body["+str(ibody)+"]->MainFrame.ComposePartialVelocities();\n")
            else:
                OutFile.write("//Body["+str(ibody)+"]\n")
                BodyName="body["+str(ibody)+"]->MainFrame"
                self.body[ibody].ExportPartialVelocitiesPlusPlus(OutFile,BodyName)
                if(self.body[ibody].mainframe>=0):
                    OutFile.write("body["+str(ibody)+"]->MainFrame.ComposePartialVelocities();\n")
                for iframe in range(0,self.body[ibody].nbrnodes):
                    OutFile.write("//Body["+str(ibody)+"] Node["+str(iframe)+"]\n")
                    BodyName="body["+str(ibody)+"]->NodeFrame["+str(iframe)+"]"
                    self.body[ibody].frame[iframe].ExportPartialVelocities(OutFile,BodyName)
                    if(self.body[ibody].frame[iframe].mainframe>=0):
                        OutFile.write("body["+str(ibody)+"]->NodeFrame["+str(iframe)+"].ComposePartialVelocities();\n")
        OutFile.write("}\n")
        OutFile.write("#endif\n")
        OutFile.write("\n\n//-----------------------------------------\n")
    # -------- end MBSClass.ExportComputePartialVelocitiesPlusPlus
    # ------ begin MBSClass.ExportAddAppliedForces
    def ExportAddAppliedForces(self,OutFile):
        OutFile.write("\nvoid AddAppliedEfforts()\n{\n")
        OutFile.write("//Contribution of gravity\n")
        OutFile.write("vec gravity("+ccode(self.gravity[0])+","+ccode(self.gravity[1])\
                       +","+ccode(self.gravity[2])+");\n")
        OutFile.write("AddGravityForces(gravity);\n")
        # Defining input in EasyDyn classic
        if(self.nbrinput>0):
            OutFile.write("\n") 
            OutFile.write("//Input(s) defined by the user\n") 
            for iInput in range(0,self.nbrinput):
                OutFile.write("   "+self.input[iInput]+";\n")
            OutFile.write("\n")  
        if(self.ForcesInAppEff != " "): # si on a mis des forces appliquées
            OutFile.write("//Contribution of user defined forces\n")
            OutFile.write(self.ForcesInAppEff)
            #OutFile.write('#include "'+self.ApplicationFileName+'.AppEff.cpp'+'"\n')
        OutFile.write("}\n\n//-----------------------------------------\n")
    # -------- end MBSClass.ExportAddAppliedForces
    # ------ begin MBSClass.ExportAddAppliedForcesPlusPlus
    def ExportAddAppliedForcesPlusPlus(self,OutFile):
        OutFile.write("\nvoid "+self.ApplicationFileName+"::AddAppliedForces()\n{\n")
        if(self.ForcesInAppEff != " "): # si on a mis des forces appliquées
            OutFile.write("//Contribution of user defined forces\n")
            OutFile.write(self.ForcesInAppEff)
            #OutFile.write('#include "'+self.ApplicationFileName+'.AppEff.cpp'+'"\n')
        OutFile.write("}\n\n//-----------------------------------------\n")
    # -------- end MBSClass.ExportAddAppliedForcesPlusPlus
    # ------ begin MBSClass.ExportComputeResidual
    def ExportComputeResidual(self,OutFile):
        #OutFile.write("\nvoid "+self.ApplicationFileName+"::ComputeResidual()\n{\n")
        #OutFile.write("ComputeResidualmbs();\n")
        #OutFile.write("}\n\n//-----------------------------------------\n")
        OutFile.write("\nvoid ComputeResidual()\n{\n")
        if(self.ComputeResidualVar!=" "):
            OutFile.write("  "+self.ComputeResidualVar+"\n")
        else:
            OutFile.write("ComputeResidualmbs();\n")
        OutFile.write("}\n\n//-----------------------------------------\n")
    # -------- end MBSClass.ExportComputeResidual
    # ------ begin MBSClass.ExportComputeResidualPlusPlus
    def ExportComputeResidualPlusPlus(self,OutFile):
        #OutFile.write("\nvoid "+self.ApplicationFileName+"::ComputeResidual()\n{\n")
        #OutFile.write("ComputeResidualmbs();\n")
        #OutFile.write("}\n\n//-----------------------------------------\n")
        if(self.ComputeResidualVar != " "):
            OutFile.write("\nvoid "+self.ApplicationFileName+"::ComputeResidual()\n{\n")
            OutFile.write("  "+self.ComputeResidualVar+"\n")
            OutFile.write("}\n\n//-----------------------------------------\n")
    # -------- end MBSClass.ExportComputeResidualPlusPlus
    # ------ begin MBSClass ExportComputeInputs
    def ExportComputeInputs(self,OutFile):
        if(self.nbrinput > 0):
            OutFile.write("\nvoid "+self.ApplicationFileName+"::ComputeInputs()\n{\n")
            for iInput in range(0,self.nbrinput):
                OutFile.write("   "+self.input[iInput]+";\n")
            OutFile.write("}\n\n//-----------------------------------------\n")
    def ExportComputeOutputs(self,OutFile):
        if(self.nbroutput > 0):
            OutFile.write("\nvoid "+self.ApplicationFileName+"::ComputeOutputs()\n{\n")
            for iOutput in range(0,self.nbroutput):
                OutFile.write("   "+self.output[iOutput]+";\n")
            OutFile.write("}\n\n//-----------------------------------------\n")
    def ExportDefineExternalInputs(self,OutFile):
        if(self.nbrexternalinput > 0):
            OutFile.write("\nvoid "+self.ApplicationFileName+"::DefineExternalInputs()\n{\n")
            for iExternalInput in range(0,self.nbrexternalinput):
                OutFile.write("   "+self.externalinput[iExternalInput]+";\n")
            OutFile.write("}\n\n//-----------------------------------------\n")
    # ------ end MBSClass ExportComputeInputs
    # ------ begin MBSClass.Exportmain
    def Exportmain(self,OutFile):
        OutFile.write("\nint main()\n{\n")
        if(self.MainVar != " "):
            OutFile.write(self.MainVar)
            OutFile.write("\n")
        OutFile.write("// Initialization and memory allocation\n")
        if(self.OverWriteNbrdof>0):
            OutFile.write("nbrdof="+str(self.OverWriteNbrdof)+";\n")
        else:
            OutFile.write("nbrdof="+str(self.nbrdof)+";\n")
        OutFile.write("nbrdep="+str(self.nbrdep)+";\n")
        OutFile.write("nbrbody="+str(self.nbrbody)+";\n")
        if(self.nbrinput>0):
            OutFile.write("nbrinput="+str(self.nbrinput)+";\n")
        OutFile.write("application=new char["+str(len(self.ApplicationFileName)+2)+"];\n")
        OutFile.write('strcpy(application,"'+self.ApplicationFileName+'");\n')
        OutFile.write("InitEasyDynmbs();\n")
        OutFile.write("// Let's create the shapes\n")
        for ibody in range(0,self.nbrbody):
            OutFile.write("shape *s"+str(ibody)+";\n")
            OutFile.write("s"+str(ibody)+"=new box(&(body["+str(ibody)\
                           +"].T0G),vcoord(-0.1,-0.1,-0.1),vcoord(0.1,0.1,0.1),1,2);\n")
            OutFile.write("thescene.AddShape(s"+str(ibody)+");\n")
        OutFile.write("// Uncomment the following line if you want a moving observer\n")
        OutFile.write("// thescene.SetVisuFrame(&Tref);\n")
        OutFile.write("// Let's open an animation file\n")
        OutFile.write('VanFile.open("'+self.ApplicationFileName+'.van");\n')
        OutFile.write("// Initial configuration\n")
        for idof in range(0,self.nbrdof):
            if(self.qini[idof]!=0):
                OutFile.write("q["+str(idof)+"]="+str(ccode(self.qini[idof]))+";\n")
            if(self.qdini[idof]!=0):
                OutFile.write("qd["+str(idof)+"]="+str(ccode(self.qdini[idof]))+";\n")
        
        DofLockedFlag=0        
        for idof in range(0,self.nbrdof):
            if(self.qdoflocked[idof]!=0):
                DofLockedFlag=DofLockedFlag+1
        if(DofLockedFlag>0):        
            OutFile.write("// Blocked degree(s) of freedom \n")
            OutFile.write("int *doflocked;\n")
            OutFile.write("doflocked=new int[nbrdof];\n")
        for idof in range(0,self.nbrdof):
            if(self.qdoflocked[idof]!=0):
                OutFile.write("doflocked["+str(idof)+"]="+str(self.qdoflocked[idof])+";\n")
                
        OutFile.write("// Let's save the structure of the scene\n")
        OutFile.write("ComputeMotion();\n")
        OutFile.write('thescene.CreateVolFile("'+self.ApplicationFileName+'.vol");\n')
        
        if(self.PerformanceTestFlag==1):
            OutFile.write("// Let's perform the efficiency test\n")
            OutFile.write("ofstream PerfFile;\n")
            OutFile.write('PerfFile.open("'+self.ApplicationFileName+'.perf");\n')
            OutFile.write("PerformanceTest(PerfFile);\n")
            OutFile.write("PerfFile.close();\n\n")
            
        if(self.StaticFlag==1):  
            OutFile.write("// Searching for equilibrium position\n")  
            if(DofLockedFlag>0): 
                OutFile.write("StaticEquilibrium(doflocked);\n\n")
            else:
                OutFile.write("StaticEquilibrium();\n\n")
                
        if(self.PoleFlag==1):
            OutFile.write("// Let's calculate the poles\n")  
            OutFile.write('cout << "Eigen value analysis" << endl;\n')
            if(DofLockedFlag>0): 
                OutFile.write("SaveLinearizedSystem(doflocked);\n")
                OutFile.write("ComputePoles(doflocked);\n")
                OutFile.write("CreateVmoFile(thescene);\n")
            else:
                OutFile.write("SaveLinearizedSystem();\n")
                OutFile.write("ComputePoles();\n")    
                OutFile.write("CreateVmoFile(thescene);\n")
            OutFile.write('cout << "Eigen values computed !" << endl;\n\n') 
            
        OutFile.write("// Let's perform the integration !\n")
        if (self.hmax>self.hsave):
            self.hmax=self.hsave/2.0
        if(self.newmarkonestep==1):
            OutFile.write("\n// Integration with NewmarkOneStep\n")
            OutFile.write("// Opening of the result file \n")
            OutFile.write('ofstream ResFile("'+self.ApplicationFileName+'.res"); \n')
            OutFile.write("// Determining initial accelerations  \n")
            OutFile.write("t=0; \n")
            OutFile.write("double errqd; \n")
            OutFile.write("NewmarkOneStep(0,errqd); \n")
            OutFile.write("// Saving initial state \n")
            OutFile.write("SaveData(ResFile); \n")
            OutFile.write("double dt="+str(self.hsave)+"; \n")
            OutFile.write("double tfinal="+str(self.tfinal)+"; \n")
            OutFile.write("int code; \n")
            OutFile.write("while (t<tfinal) \n")
            OutFile.write("  { \n")
            OutFile.write("  code=NewmarkOneStep(dt,errqd); \n")
            OutFile.write("  cout << \"t=\" << t << endl; \n")
            OutFile.write("  if (code==1) { cout << \"No convergence\"; system(\"pause\"); exit(1); } \n")
            OutFile.write("  if (code==2) { cout << \"Numerical trouble\"; system(\"pause\"); exit(2); } \n")
            OutFile.write("  SaveData(ResFile); \n")
            OutFile.write("  }\n")
            OutFile.write("ResFile.close();\n")
        else:
            if(self.lockintegration>0 and DofLockedFlag>0):
                OutFile.write("NewmarkIntegration("+str(self.tfinal)+","\
                          +str(self.hsave)+","+str(self.hmax)+",doflocked);\n")
            else:
                OutFile.write("NewmarkIntegration("+str(self.tfinal)+","\
                          +str(self.hsave)+","+str(self.hmax)+");\n")
        OutFile.write("// The clean way to finish\n")
        OutFile.write("EndEasyDynmbs();\n")
        OutFile.write("VanFile.close();\n")
        OutFile.write("\n}\n//-----------------------------------------\n")
    # ------ end MBSClass.Exportmain
    # ------ begin MBSClass.ExportmainPlusPlus
    def ExportmainPlusPlus(self,OutFile):
        OutFile.write("\nint main()\n{\n\n")
        if(self.MainVar != " "):
            OutFile.write(self.MainVar)
            OutFile.write("\n")
        OutFile.write(self.ApplicationFileName+" "+self.ApplicationFileName+"Obj;\n\n")
        OutFile.write("// Let's create the shapes\n")
        for ibody in range(0,self.nbrbody):
            OutFile.write("shape *s"+str(ibody)+";\n")
            OutFile.write("s"+str(ibody)+"=new box(&("+self.ApplicationFileName+"Obj.body["+str(ibody)\
                           +"]->MainFrame.T0F),vcoord(-0.1,-0.1,-0.1),vcoord(0.1,0.1,0.1),1,2);\n")
            OutFile.write("thescene.AddShape(s"+str(ibody)+");\n")
            for iframe in range(0,self.body[ibody].nbrnodes):
                OutFile.write("shape *s"+str(ibody)+"f"+str(iframe)+";\n")
                OutFile.write("s"+str(ibody)+"f"+str(iframe)+"=new box(&("+self.ApplicationFileName+"Obj.body["+str(ibody)\
                               +"]->NodeFrame["+str(iframe)+"].T0F),vcoord(-0.1,-0.1,-0.1),vcoord(0.1,0.1,0.1),1,2);\n")
                OutFile.write("thescene.AddShape(s"+str(ibody)+"f"+str(iframe)+");\n")
        OutFile.write("\n// Uncomment the following line if you want a moving observer\n")
        OutFile.write("// thescene.SetVisuFrame(&Tref);\n\n")
        OutFile.write("// Let's open an animation file\n")
        OutFile.write('VanFile.open("'+self.ApplicationFileName+'.van");\n\n')
        OutFile.write("// Initial configuration of the system\n")
        for idof in range(0,self.nbrdof):
            if(self.qini[idof]!=0):
                OutFile.write(self.ApplicationFileName+"Obj.q["+str(idof)+"]="+str(ccode(self.qini[idof]))+";\n")
            if(self.qdini[idof]!=0):
                OutFile.write(self.ApplicationFileName+"Obj.qd["+str(idof)+"]="+str(ccode(self.qdini[idof]))+";\n")
                
        DofLockedFlag=0        
        for idof in range(0,self.nbrdof):
            if(self.qdoflocked[idof]!=0):
                DofLockedFlag=DofLockedFlag+1
        if(DofLockedFlag>0):        
            OutFile.write("// Blocked degree(s) of freedom \n")
            OutFile.write("int *doflocked;\n")
            OutFile.write("doflocked=new int["+str(self.nbrdof)+"];\n")
        for idof in range(0,self.nbrdof):
            if(self.qdoflocked[idof]!=0):
                OutFile.write("doflocked["+str(idof)+"]="+str(self.qdoflocked[idof])+";\n")
                
        OutFile.write("\n")
        OutFile.write("// Let's save the structure of the scene\n")
        OutFile.write(self.ApplicationFileName+"Obj.ComputeMotion();\n")
        OutFile.write('thescene.CreateVolFile("'+self.ApplicationFileName+'.vol");\n\n')
        
        if(self.PerformanceTestFlag==1):
            OutFile.write("// Let's perform the efficiency test\n")
            OutFile.write("ofstream PerfFile;\n")
            OutFile.write('PerfFile.open("'+self.ApplicationFileName+'.perf");\n')
            OutFile.write(self.ApplicationFileName+"Obj.PerformanceTest(PerfFile);\n")
            OutFile.write("PerfFile.close();\n\n")
            
        if(self.StaticFlag==1):  
            OutFile.write("// Searching for equilibrium position\n")  
            if(DofLockedFlag>0):  
                OutFile.write(self.ApplicationFileName+"Obj.StaticEquilibrium(doflocked);\n\n")
            else:
                OutFile.write(self.ApplicationFileName+"Obj.StaticEquilibrium();\n\n")
            
        if(self.PoleFlag==1):
            OutFile.write("// Let's calculate the poles\n")  
            OutFile.write('cout << "Eigen value analysis" << endl;\n')
            if(DofLockedFlag>0):  
                OutFile.write(self.ApplicationFileName+"Obj.SaveLinearizedSystem(doflocked);\n")
                OutFile.write(self.ApplicationFileName+"Obj.ComputePoles(doflocked);\n")
                OutFile.write(self.ApplicationFileName+"Obj.CreateVmoFile(thescene);\n")
            else:
                OutFile.write(self.ApplicationFileName+"Obj.SaveLinearizedSystem();\n")
                OutFile.write(self.ApplicationFileName+"Obj.ComputePoles();\n")
                OutFile.write(self.ApplicationFileName+"Obj.CreateVmoFile(thescene);\n")
            OutFile.write('cout << "Eigen values computed !" << endl;\n\n')        

        OutFile.write("// Let's perform the integration\n")
        if (self.hmax>self.hsave):
            self.hmax=self.hsave/2.0
        if(self.newmarkonestep==1):
            OutFile.write("\n// Integration with NewmarkOneStep\n")
            OutFile.write("// Opening of the result file \n")
            OutFile.write('ofstream ResFile("'+self.ApplicationFileName+'.res"); \n')
            OutFile.write("// Determining initial accelerations  \n")
            OutFile.write(self.ApplicationFileName+"Obj.t=0; \n")
            OutFile.write("double errqd; \n")
            OutFile.write(self.ApplicationFileName+"Obj.NewmarkOneStep(0,errqd); \n")
            OutFile.write("// Saving initial state \n")
            OutFile.write(self.ApplicationFileName+"Obj.SaveData(ResFile); \n")
            OutFile.write("double dt="+str(self.hsave)+"; \n")
            OutFile.write("double tfinal="+str(self.tfinal)+"; \n")
            OutFile.write("int code; \n")
            OutFile.write("while ("+self.ApplicationFileName+"Obj.t<tfinal) \n")
            OutFile.write("  { \n")
            OutFile.write("  code="+self.ApplicationFileName+"Obj.NewmarkOneStep(dt,errqd); \n")
            OutFile.write("  cout << \"t=\" << "+self.ApplicationFileName+"Obj.t << endl; \n")
            OutFile.write("  if (code==1) { cout << \"No convergence\"; system(\"pause\"); exit(1); } \n")
            OutFile.write("  if (code==2) { cout << \"Numerical trouble\"; system(\"pause\"); exit(2); } \n")
            OutFile.write("  "+self.ApplicationFileName+"Obj.SaveData(ResFile); \n")
            OutFile.write("  }\n")
            OutFile.write("ResFile.close();\n")
        else:
            if(self.lockintegration>0 and DofLockedFlag>0):
                OutFile.write(self.ApplicationFileName+"Obj.NewmarkIntegration("+str(self.tfinal)+","\
                            +str(self.hsave)+","+str(self.hmax)+",doflocked);\n")
            else:
                OutFile.write(self.ApplicationFileName+"Obj.NewmarkIntegration("+str(self.tfinal)+","\
                            +str(self.hsave)+","+str(self.hmax)+");\n")
        
        OutFile.write("\n// The clean way to finish\n")
        OutFile.write("VanFile.close();\n")
        OutFile.write("}\n//-----------------------------------------\n")
    # -------- end MBSClass.ExportmainPlusPlus
    # ------ begin ExportGnuplotScript
    def ExportGnuplotScript(self):
        #.plt file to plot graphs
        PlotFileName=self.ApplicationFileName+".plt"
        print("Creating Plt Plot File "+PlotFileName)
        PlotFile = open(PlotFileName, "w")
                
        PlotFile.write("reset\n")
        PlotFile.write("set xlabel \"Time [s]\" \n")
        PlotFile.write("set grid \n\n")
        PlotFile.write("set term postscript eps color \"Times-Roman\" 20  \n")
        PlotFile.write("set output \"figure1.eps\" \n")
        PlotFile.write("set ylabel \"displacements\" \n")
        PlotFile.write("plot")
        for idof in  range (0,self.nbrdof):
            PlotFile.write(" '"+self.ApplicationFileName+".res' using 1:"+str(idof*3+2)+" title 'q_"+str(idof)+"' with line ")
            if idof <((self.nbrdof)-1):
                PlotFile.write(",")               
        PlotFile.write("\nset term pop \n")
        PlotFile.write("replot \n")
        PlotFile.write("pause -1 'Next plot (velocity level)?' \n\n")
        PlotFile.write("set term postscript eps color \"Times-Roman\" 20  \n")
        PlotFile.write("set output \"figure2.eps\" \n")
        PlotFile.write("set ylabel \"velocities\" \n")
        PlotFile.write("plot")
        for idof in  range (0,self.nbrdof):
            PlotFile.write(" '"+self.ApplicationFileName+".res' using 1:"+str(idof*3+3)+" title 'qd_"+str(idof)+"' with line ")
            if idof <((self.nbrdof)-1):
                PlotFile.write(",")               
        PlotFile.write("\nset term pop \n")
        PlotFile.write("replot \n")
        PlotFile.write("pause -1 'Next plot (acceleration level)?' \n\n")
        PlotFile.write("set term postscript eps color \"Times-Roman\" 20  \n")
        PlotFile.write("set output \"figure3.eps\" \n")
        PlotFile.write("set ylabel \"accelerations\" \n")
        PlotFile.write("plot")
        for idof in  range (0,self.nbrdof):
            PlotFile.write(" '"+self.ApplicationFileName+".res' using 1:"+str(idof*3+4)+" title 'qdd_"+str(idof)+"' with line ")
            if idof <((self.nbrdof)-1):
                PlotFile.write(",")               
        PlotFile.write("\nset term pop \n")
        PlotFile.write("replot \n")
        PlotFile.write("pause -1 \n") 
        PlotFile.close()
        print(" ") 
    # ------ end ExportGnuplotScript
    # ------ begin ExportExtraPlot
    def ExportExtraPlot(self,OutFile):
        OutFile.write("reset\n")
        OutFile.write("set xlabel \"Time [s]\" \n")
        OutFile.write("set grid \n\n")
        nbrcolumn=self.nbrdof*3+1
        for ibody in range (0,self.nbrbody):
            self.body[ibody].ExportExtraPlot(OutFile,ibody,self.ApplicationFileName,nbrcolumn)
            nbrcolumn=nbrcolumn+self.body[ibody].ExtraResults
            for iframe in range(0,self.body[ibody].nbrnodes):
                self.body[ibody].frame[iframe].ExportExtraPlot(OutFile,ibody,iframe,self.ApplicationFileName,nbrcolumn)
                nbrcolumn=nbrcolumn+self.body[ibody].frame[iframe].ExtraResults
        OutFile.write("pause -1 \n") 
    # ------ end ExportExtraPlot
    # ------ begin ExportFR_Latex_Report 
    def ExportFR_Latex_Report(self):
        print("Processing FR Latex files") 
        #.tex file to Latex File
        TexFRFileName="Rapport_"+(self.ApplicationFileName)+".tex"
        print("Creating 'French' Tex File "+TexFRFileName)
        TexFRFile = open(TexFRFileName, "w")
        
        TexFRFile.write("\\documentclass[a4paper,11pt]{article}\n")
        TexFRFile.write("\\topmargin -0.5cm\n")
        TexFRFile.write("\\headheight 0.7cm\n")
        TexFRFile.write("\\headsep    0.8cm\n")
        TexFRFile.write("\\topskip      0cm\n")
        TexFRFile.write("\\textheight  24cm\n")
        TexFRFile.write("\\oddsidemargin  0cm\n")
        TexFRFile.write("\\evensidemargin 0cm\n")
        TexFRFile.write("\\textwidth   16cm\n")
        TexFRFile.write("\\usepackage{verbatim}\n")
        TexFRFile.write("\\usepackage{array} %>{$}m{3cm}<{$}\n")
        TexFRFile.write("\\usepackage{amssymb}\n")
        TexFRFile.write("\\usepackage{amsmath}\n")
        TexFRFile.write("\\usepackage[dvips]{graphicx}\n")
        TexFRFile.write("\\usepackage[T1]{fontenc}\n")
        TexFRFile.write("\\usepackage[francais]{babel}\n")
        TexFRFile.write("\\usepackage{fancyhdr}\n")
        TexFRFile.write("\\pagestyle{fancy}\n")
        TexFRFile.write("\\fancyhead{}\n")
        TexFRFile.write("\\fancyfoot{}\n")
        TexFRFile.write("\\fancyhead[L]{"+self.ApplicationTitle+"}\n")
        TexFRFile.write("\\fancyhead[R]{\\thepage}\n")
        TexFRFile.write("\\title{\\textsf{"+self.ApplicationTitle+"}}\n")
        TexFRFile.write("\\author{}\n")
        TexFRFile.write("\\begin{document}\n")
        TexFRFile.write("\\sloppy\n")
        TexFRFile.write("\\maketitle\n")
        TexFRFile.write("\n")
        TexFRFile.write("\\section{Donn\\'ees g\\'en\\'erales du m\\'ecanisme \\'etudi\\'e}\n")
        
        # Cas dans lequel il n'y a pas de corps flexibles => cas purement rigide
        if (self.nbrFlexBeam+self.nbrFlexBody)==0: 
            TexFRFile.write("\nLe nombre de solide(s) du syst\\`eme s'\\'el\\`eve \\`a "+str(self.nbrbody)+" (d\\'efini par la variable globale \\texttt{nbrbody}). ")
            TexFRFile.write("Par la suite, ils se d\\'efinissent par $S_j$ avec $j$ allant de 0 \\`a "+str(self.nbrbody-1)+". ")
            if self.nbrdof >1:
                TexFRFile.write("Le nombre de degr\\'es de libert\\'e s'\\'el\\`eve, quant \\`a lui, \\`a "+str(self.nbrdof)+" (\\texttt{nbrdof}). ")
                TexFRFile.write("Les param\\`etres de configuration sont donc d\\'efinis par $q_i$ ($i$ allant de 0 \\`a "+str(self.nbrdof-1)+").\\\\\n\n")
            if self.nbrdof ==1:
                TexFRFile.write("Le m\\'ecanisme analys\\'e est desmodrome (\\`A un degr\\'e de libert\\'e) et donc ")
                TexFRFile.write("un seul param\\`etre de configuration est d\\'efini par $q_0$.\\\\\n\n")
            TexFRFile.write("Les donn\\'ees inertielles, d\\'efinies dans le fichier utilisateur, se r\\'esument aux masses $m_{Si}$ et ")
            TexFRFile.write("aux tenseurs d'inertie $\\Phi_{G,Si}$ d\\'efinis au centre de gravit\\'e de chaque solide concern\\'e.\\\\\n\n")
        # Cas général: avec corps rigides, flexibles et poutres flexibles
        else:  
            TexFRFile.write("\nLe nombre de solide(s) du syst\\`eme s'\\'el\\`eve \\`a "+str(self.nbrbody)+" (d\\'efini par la variable globale \\texttt{nbrbody}): ")
            TexFRFile.write("\n\\begin{itemize}\n")
            if (self.nbrRigidBody)>0:
                TexFRFile.write("  \\item On d\\'enombre "+str(self.nbrRigidBody)+" corps rigide(s).\n")
            if (self.nbrFlexBeam)>0:
                TexFRFile.write("  \\item On d\\'enombre "+str(self.nbrFlexBeam)+" poutre(s) flexible(s).\n")
            if (self.nbrFlexBody)>0:
                TexFRFile.write("  \\item On d\\'enombre "+str(self.nbrFlexBody)+" corps flexible(s).\n")
            TexFRFile.write("\\end{itemize}\n")
            if self.nbrdof >1:
                TexFRFile.write("Le nombre de degr\\'es de libert\\'e s'\\'el\\`eve, quant \\`a lui, \\`a "+str(self.nbrdof)+" (\\texttt{nbrdof}). ")
                TexFRFile.write("Les param\\`etres de configuration sont donc d\\'efinis par $q_i$ ($i$ allant de 0 \\`a "+str(self.nbrdof-1)+").\\\\\n\n")
            if self.nbrdof ==1:
                TexFRFile.write("Le m\\'ecanisme analys\\'e est desmodrome (\\`a un degr\\'e de libert\\'e) et donc ")
                TexFRFile.write("un seul param\\`etre de configuration est d\\'efini par $q_0$.\\\\\n\n")
                
        # Cas dans lequel il n'y a pas de corps flexibles => cas purement rigide
        if (self.nbrFlexBeam+self.nbrFlexBody)==0:
            for ibody in range(0,self.nbrbody):
                TexFRFile.write("\\[ m_{S"+str(ibody)+"}="+str(latex(self.body[ibody].mass))+"\\,kg \\]\n")
            TexFRFile.write("\\\\\n")
        # Cas général: avec corps rigides, flexibles et poutres flexibles
        else:  
            if (self.nbrRigidBody)>0:
                TexFRFile.write("\\subsection*{Corps rigide(s): donn\\'ees inertielles}\n\n")
                for ibody in range(0,self.nbrbody):
                    if(self.body[ibody].mass>0 and self.body[ibody].Ixx>0 and self.body[ibody].Iyy>0 and self.body[ibody].Izz>0 and self.body[ibody].Ixy>=0 and self.body[ibody].Ixz>=0 and self.body[ibody].Iyz>=0):
                        TexFRFile.write("Corps rigide S"+str(ibody)+": ")
                        TexFRFile.write("\\[ m_{S"+str(ibody)+"}="+str(latex(self.body[ibody].mass))+"\\,kg \\]\n")
                        TexFRFile.write("\\[\\Phi_{G,S"+str(ibody)+"}="
                        +str(latex(Matrix(((self.body[ibody].Ixx,-self.body[ibody].Ixy,-self.body[ibody].Ixz),
                        (-self.body[ibody].Ixy,self.body[ibody].Iyy,-self.body[ibody].Ixy),
                        (-self.body[ibody].Ixz,-self.body[ibody].Iyz,self.body[ibody].Izz)))))+" \\quad \\textrm{, en $kg.m^2$}\\]\\\\\n")
            # si c'est une poutre flexible: afficher ses propriétés
            if (self.nbrFlexBeam)>0:
                TexFRFile.write("\\subsection*{Poutre(s) flexible(s): propri\\'et\\'e(s)}\n\n")
                for ibody in range(0,self.nbrbody):
                    # si c'est une poutre flexible: afficher ses propriétés
                    if(self.body[ibody].length>0 and self.body[ibody].section>0 and self.body[ibody].EYoung>0 and self.body[ibody].nu>0 and self.body[ibody].rho>0 and self.body[ibody].Iyy>0 and self.body[ibody].Izz>0 and self.body[ibody].AlphaDamp>=0 and self.body[ibody].BetaDamp>=0 and self.body[ibody].nbrnodes>=2):
                        TexFRFile.write("Poutre flexible S"+str(ibody)+": ")
                        TexFRFile.write("\n\\begin{itemize}\n")
                        TexFRFile.write("  \\item Longueur: "+str(self.body[ibody].length)+" m;\n")
                        TexFRFile.write("  \\item Section: "+str(self.body[ibody].section)+" $m^2$;\n")
                        TexFRFile.write("  \\item Module de Young: "+str(self.body[ibody].EYoung)+" Pa;\n")
                        TexFRFile.write("  \\item Coefficient de Poisson: "+str(self.body[ibody].nu)+";\n")
                        TexFRFile.write("  \\item Masse volumique: "+str(self.body[ibody].rho)+" $kg/m^3$;\n")
                        TexFRFile.write("  \\item Inertie g\\'eom\\'etrique $I_{yy}$: "+str(self.body[ibody].Iyy)+" $m^4$;\n")
                        TexFRFile.write("  \\item Inertie g\\'eom\\'etrique $I_{zz}$: "+str(self.body[ibody].Izz)+" $m^4$;\n")
                        TexFRFile.write("  \\item Nombre de frames: "+str(self.body[ibody].nbrnodes)+".\\\\\n")
                        TexFRFile.write("\\end{itemize}\n")
            if (self.nbrFlexBody)>0:
                TexFRFile.write("\\subsection*{Corps flexible(s): matrice(s) de masse et de raideur}\n\n")
                for ibody in range(0,self.nbrbody):
                    if(self.body[ibody].name!=0 and self.body[ibody].nbrnodes>=2 and self.body[ibody].AlphaDamp>=0 and self.body[ibody].BetaDamp>=0 and self.body[ibody].length==0):
                        TexFRFile.write("Corps flexible S"+str(ibody)+": \\\\\n")
                        # open Mass Matrix file
                        with open(str(self.body[ibody].name)+".mm") as f:
                            MatriceMM = []
                            for line in f:
                                line = line.split() # to deal with blank 
                                if line:            # lines (ie skip them)
                                    line = [i for i in line]
                                    MatriceMM.append(line) 
                                    MM=Matrix(MatriceMM)
                        TexFRFile.write("\\[M_{"+str(ibody)+"}="+str(latex(MM))+"\\]\\\\\n\n")
                        # open Stiffness Matrix file
                        with open(str(self.body[ibody].name)+".kk") as f:
                            MatriceKK = []
                            for line in f:
                                line = line.split() # to deal with blank 
                                if line:            # lines (ie skip them)
                                    line = [i for i in line]
                                    MatriceKK.append(line) 
                                    KK=Matrix(MatriceKK)
                        TexFRFile.write("\\[K_{"+str(ibody)+"}="+str(latex(KK))+"\\]\\\\\n\n")
        
        # Cas dans lequel il n'y a pas de corps flexibles => cas purement rigide
        if (self.nbrFlexBeam+self.nbrFlexBody)==0:
            for ibody in range(0,self.nbrbody):
                TexFRFile.write("\\[\\Phi_{G,S"+str(ibody)+"}="
                +str(latex(Matrix(((self.body[ibody].Ixx,-self.body[ibody].Ixy,-self.body[ibody].Ixz),
                (-self.body[ibody].Ixy,self.body[ibody].Iyy,-self.body[ibody].Ixy),
                (-self.body[ibody].Ixz,-self.body[ibody].Iyz,self.body[ibody].Izz)))))+" \\quad \\textrm{, en $kg.m^2$}\\]\\\\\n")
        TexFRFile.write("\n\\section{Param\\`etres cin\\'ematiques calcul\\'es avec Sympy}\n\n")
        TexFRFile.write("Ces param\\`etres ont \\'et\\'e calcul\\'es \\`a partir du fichier utilisateur \\texttt{"+str(self.ApplicationFileName)+".py} et ce, \n")
        TexFRFile.write("avec un temps de calcul \\emph{CPU} de "+str(round(self.EndTime_Kinematics - self.StartTime_Kinematics,3))+" seconde(s). \\\\\n")
        TexFRFile.write("\n")
        
        # Cas dans lequel il n'y a pas de corps flexibles => cas purement rigide
        if (self.nbrFlexBeam+self.nbrFlexBody)==0: 
            NbRelBodies=0 # Compute the number of relative bodies
            for ibody in range(0,self.nbrbody):
                if self.body[ibody].mainframe>=0:
                    NbRelBodies = NbRelBodies+1
            if NbRelBodies>1:
                TexFRFile.write("\\subsection*{Les mouvements relatifs} \n\n")
                TexFRFile.write("Certains mouvements de solides ont \\'et\\'e d\\'efinis de mani\\`ere relative par rapport \\`a d'autres solides.")
                TexFRFile.write(" Ces solides sont au nombre de "+str(NbRelBodies)+" et se d\\'efinissent de la mani\\'ere suivante :\n")
                TexFRFile.write("\\begin{itemize}\n")
                for ibody in range(0,self.nbrbody):
                    if self.body[ibody].mainframe>=0:
                        TexFRFile.write("  \\item Le mouvement du solide $S_{"+str(ibody)+"}$ est d\\'efini par rapport \\`a celui de $S_{"+str(self.body[ibody].mainframe)+"}$.\n")
                TexFRFile.write("\\end{itemize}\n")
            else:
                for ibody in range(0,self.nbrbody):
                    if self.body[ibody].mainframe>=0:
                        TexFRFile.write("\\noindent \\emph{Remarque} : Le mouvement du solide $S_{"+str(ibody)+"}$ est un mouvement relatif et est d\\'efini par rapport \\`a celui de $S_{"+str(self.body[ibody].mainframe)+"}$. \n") 
        else:
            Jar=0
            for ibody in range(0,self.nbrbody):
                if(self.body[ibody].mainframe>=0 or self.body[ibody].frameref>=0):
                    Jar=Jar+1
                for iframe in range(0,self.body[ibody].nbrnodes):
                    if(self.body[ibody].frame[iframe].mainframe>=0 or self.body[ibody].frame[iframe].frameref>=0):
                        Jar=Jar+1
            if Jar>0: # S'il y a des mouvements relatifs
                TexFRFile.write("\\subsection*{Les mouvements relatifs} \n\n")
                for ibody in range(0,self.nbrbody):
                    #corps purement rigide avec frames possibles
                    if(self.body[ibody].mass>0 and self.body[ibody].Ixx>0 and self.body[ibody].Iyy>0 and self.body[ibody].Izz>0 and self.body[ibody].Ixy>=0 and self.body[ibody].Ixz>=0 and self.body[ibody].Iyz>=0):
                        DisplayOK=0
                        if(self.body[ibody].mainframe>=0 or self.body[ibody].frameref>=0):
                            #Si mouvement relatif, on affiche
                            DisplayOK=1 
                        PassDisplayOK=0
                        for iframe in range(0,self.body[ibody].nbrnodes):
                            #Si mouvement relatif, on affiche
                            if(( (self.body[ibody].frame[iframe].mainframe>=0 or self.body[ibody].frame[iframe].frameref>=0) or DisplayOK==1) and PassDisplayOK==0):
                                TexFRFile.write("Corps rigide S"+str(ibody)+" avec "+str(self.body[ibody].nbrnodes)+" frames\n") 
                                PassDisplayOK=1
                    #poutre flexible avec ses frames
                    if(self.body[ibody].length>0 and self.body[ibody].section>0 and self.body[ibody].EYoung>0 and self.body[ibody].nu>0 and self.body[ibody].rho>0 and self.body[ibody].Iyy>0 and self.body[ibody].Izz>0 and self.body[ibody].AlphaDamp>=0 and self.body[ibody].BetaDamp>=0 and self.body[ibody].nbrnodes>=2):
                        DisplayOK=0
                        if(self.body[ibody].mainframe>=0 or self.body[ibody].frameref>=0):
                            #Si mouvement relatif, on affiche
                            DisplayOK=1
                        PassDisplayOK=0 
                        for iframe in range(0,self.body[ibody].nbrnodes):
                            #Si mouvement relatif, on affiche
                            if(( (self.body[ibody].frame[iframe].mainframe>=0 or self.body[ibody].frame[iframe].frameref>=0) or DisplayOK==1) and PassDisplayOK==0):
                                TexFRFile.write("Poutre flexible S"+str(ibody)+" avec "+str(self.body[ibody].nbrnodes)+" frames\n") 
                                PassDisplayOK=1
                    #corps flexible avec ses frames
                    if(self.body[ibody].name!=0 and self.body[ibody].nbrnodes>=2 and self.body[ibody].AlphaDamp>=0 and self.body[ibody].BetaDamp>=0 and self.body[ibody].length==0):
                        DisplayOK=0
                        if(self.body[ibody].mainframe>=0 or self.body[ibody].frameref>=0):
                            #Si mouvement relatif, on affiche
                            DisplayOK=1 
                        PassDisplayOK=0
                        for iframe in range(0,self.body[ibody].nbrnodes):
                            #Si mouvement relatif, on affiche
                            if(( (self.body[ibody].frame[iframe].mainframe>=0 or self.body[ibody].frame[iframe].frameref>=0) or DisplayOK==1) and PassDisplayOK==0):
                                TexFRFile.write("Corps flexible S"+str(ibody)+" avec "+str(self.body[ibody].nbrnodes)+" frames\n") 
                                PassDisplayOK=1  
    
                    if(self.body[ibody].mainframe<0):
                    # Rigid body en absolu (display nothing)
                        Jar2=0
                        for iframe in range(0,self.body[ibody].nbrnodes):
                            if(self.body[ibody].frame[iframe].mainframe>=0):
                                Jar2=Jar2+1
                        if Jar2>0:
                            TexFRFile.write("\\begin{itemize}\n")
                        for iframe in range(0,self.body[ibody].nbrnodes):
                            # Frame en absolu (display nothing)
                            #if(self.body[ibody].frame[iframe].mainframe<0 and self.body[ibody].frame[iframe].frameref<0):
                                #TexFRFile.write("\\[T_{0F,S"+str(ibody)+"F"+str(iframe)+"}="+str(latex(self.body[ibody].frame[iframe].TtoF))+"\\]\\\\\n\n")
                            # Flexible sur un solide rigide
                            if(self.body[ibody].frame[iframe].mainframe>=0 and self.body[ibody].frame[iframe].frameref<0):
                                TexFRFile.write("  \\item dont le frame F"+str(iframe)+" est relatif au frame principal du solide S"+str(self.body[ibody].frame[iframe].mainframe)+". \n")
                            # Flexible sur un solide flexible
                            elif(self.body[ibody].frame[iframe].mainframe>=0 and self.body[ibody].frame[iframe].frameref>=0):
                                TexFRFile.write("  \\item dont le frame F"+str(iframe)+" est relatif au frame F"+str(self.body[ibody].frame[iframe].frameref)+" du solide S"+str(self.body[ibody].frame[iframe].mainframe)+". \n")
                        if Jar2>0:
                            TexFRFile.write("\\end{itemize}\n")
                    # Mainframe du body est en relatif
                    else:
                        TexFRFile.write("\\begin{itemize}\n")
                        # un solide rigide par rapport à un solide rigide
                        if (self.body[ibody].frameref<0):
                            TexFRFile.write("  \\item dont le frame principal est relatif au frame principal du solide S"+str(self.body[ibody].mainframe)+". \n")    
                        # un solide rigide par rapport à un solide flexible
                        else:
                            TexFRFile.write("  \\item dont le frame principal est relatif au frame F"+str(self.body[ibody].frameref)+" du solide S"+str(self.body[ibody].mainframe)+". \n")
                        Jar2=0 
                        for iframe in range(0,self.body[ibody].nbrnodes):
                            if(self.body[ibody].frame[iframe].mainframe>=0):
                                Jar2=Jar2+1
                        PassEndItem=0
                        if Jar2==0:
                            TexFRFile.write("\\end{itemize}\n")
                            PassEndItem=1
                        for iframe in range(0,self.body[ibody].nbrnodes):
                            # Frame en absolu
                            #if(self.body[ibody].frame[iframe].mainframe<0 and self.body[ibody].frame[iframe].frameref<0):
                                #TexFRFile.write("\\[T_{0F,S"+str(ibody)+"F"+str(iframe)+"}="+str(latex(self.body[ibody].frame[iframe].TtoF))+"\\]\\\\\n\n")
                            # Flexible sur un solide rigide
                            if(self.body[ibody].frame[iframe].mainframe>=0 and self.body[ibody].frame[iframe].frameref<0):
                                TexFRFile.write("  \\item dont le frame F"+str(iframe)+" est relatif au frame principal du solide S"+str(self.body[ibody].frame[iframe].mainframe)+". \n")
                                # Flexible sur un solide flexible
                            elif(self.body[ibody].frame[iframe].mainframe>=0 and self.body[ibody].frame[iframe].frameref>=0):
                                TexFRFile.write("  \\item dont le frame F"+str(iframe)+" est relatif au frame F"+str(self.body[ibody].frame[iframe].frameref)+" du solide S"+str(self.body[ibody].frame[iframe].mainframe)+". \n")
                        if(PassEndItem==0):
                            TexFRFile.write("\\end{itemize}\n")       
        TexFRFile.write("\n")

        if self.nbrdep > 0:
            TexFRFile.write("\\subsection*{Les variables d\\'ependantes}\n\n")
            TexFRFile.write("L'utilisation de variables interm\\'ediaires s'est av\\'er\\'ee n\\'ecessaire afin de r\\'eduire les temps de calculs. ")
            TexFRFile.write("Elles sont au nombre de "+str(self.nbrdep)+" (d\\'efini par $p_i$, $i$ allant de 0 \\`a "+str(self.nbrdep-1)+") :\n")
            for idep in range(0,self.nbrdep):
                TexFRFile.write("\\[p_"+str(idep)+"="+str(latex(self.pexpr[idep]))+"\\]\n")
                TexFRFile.write("\\[\\dot{p}_"+str(idep)+"="+str(latex(self.pdexpr[idep]))+"\\]\n")
                TexFRFile.write("\\[\\ddot{p}_"+str(idep)+"="+str(latex(self.pddexpr[idep]))+"\\]\n")
        # Cas dans lequel il n'y a pas de corps flexibles => cas purement rigide
        if (self.nbrFlexBeam+self.nbrFlexBody)==0: 
            TexFRFile.write("\\subsection*{Les matrices de transformation homog\\`ene pour chaque solide}\n\n")
            for ibody in range(0,self.nbrbody):
                if self.body[ibody].mainframe<0:
                    TexFRFile.write("\\[T_{0G,S"+str(ibody)+"}="+str(latex(self.body[ibody].TtoF))+"\\]\\\\\n\n")
                else:
                    TexFRFile.write("\\[T_{refG,S"+str(ibody)+"/S"+str(self.body[ibody].mainframe)+"}="+str(latex(self.body[ibody].TtoF))+"\\]\\\\\n\n")
            TexFRFile.write("\\subsection*{Leur d\\'eriv\\'ee par rapport au temps}\n\n")
            for ibody in range(0,self.nbrbody):
                if self.body[ibody].mainframe<0:
                    TexFRFile.write("\\[\\dot{T}_{0G,S"+str(ibody)+"}="+str(latex(self.body[ibody].TtoFd))+"\\]\\\\\n\n")
                else:
                    TexFRFile.write("\\[\\dot{T}_{refG,S"+str(ibody)+"/S"+str(self.body[ibody].mainframe)+"}="+str(latex(self.body[ibody].TtoFd))+"\\]\\\\\n\n")
            TexFRFile.write("\\subsection*{Les vitesses des centres de gravit\\'e des solides}\n\n")
            for ibody in range(0,self.nbrbody):
                if self.body[ibody].mainframe<0:
                    TexFRFile.write("\\[\\vec{v}_{G,S"+str(ibody)+"}="+str(latex(self.body[ibody].vF))+"\\]\\\\\n\n")
                else:
                    TexFRFile.write("\\[\\left\\{\\vec{v}_{G,S"+str(ibody)+"/S"+str(self.body[ibody].mainframe)+"}\\right\\}_{S"+str(self.body[ibody].mainframe)+"}="+str(latex(self.body[ibody].vF))+"\\]\\\\\n\n")
            TexFRFile.write("\\subsection*{Les acc\\'el\\'erations des centres de gravit\\'e des solides}\n\n")
            for ibody in range(0,self.nbrbody):
                if self.body[ibody].mainframe<0:
                    TexFRFile.write("\\[\\vec{a}_{G,S"+str(ibody)+"}="+str(latex(self.body[ibody].aF))+"\\]\\\\\n\n")
                else:
                    TexFRFile.write("\\[\\left\\{\\vec{a}_{G,S"+str(ibody)+"/S"+str(self.body[ibody].mainframe)+"}\\right\\}_{S"+str(self.body[ibody].mainframe)+"}="+str(latex(self.body[ibody].aF))+"\\]\\\\\n\n")
            TexFRFile.write("\\subsection*{Les vitesses de rotation des solides}\n\n")
            for ibody in range(0,self.nbrbody):
                if self.body[ibody].mainframe<0:
                    TexFRFile.write("\\[\\vec{\omega}_{G,S"+str(ibody)+"}="+str(latex(self.body[ibody].omega))+"\\]\\\\\n\n")
                else:
                    TexFRFile.write("\\[\\left\\{\\vec{\omega}_{G,S"+str(ibody)+"/S"+str(self.body[ibody].mainframe)+"}\\right\\}_{S"+str(self.body[ibody].mainframe)+"}="+str(latex(self.body[ibody].omega))+"\\]\\\\\n\n")
            TexFRFile.write("\\subsection*{Les acc\\'el\\'erations de rotation des solides}\n\n")
            for ibody in range(0,self.nbrbody):
                if self.body[ibody].mainframe<0:
                    TexFRFile.write("\\[\\vec{\dot{\omega}}_{G,S"+str(ibody)+"}="+str(latex(self.body[ibody].omegad))+"\\]\\\\\n\n")
                else:
                    TexFRFile.write("\\[\\left\\{\\vec{\dot{\omega}}_{G,S"+str(ibody)+"/S"+str(self.body[ibody].mainframe)+"}\\right\\}_{S"+str(self.body[ibody].mainframe)+"}="+str(latex(self.body[ibody].omegad))+"\\]\\\\\n\n")
        else:
            TexFRFile.write("\\subsection*{Les matrices de transformation homog\\`ene pour chaque solide et chaque frame}\n\n")
            for ibody in range(0,self.nbrbody):
                #corps purement rigide avec frames possibles
                if(self.body[ibody].mass>0 and self.body[ibody].Ixx>0 and self.body[ibody].Iyy>0 and self.body[ibody].Izz>0 and self.body[ibody].Ixy>=0 and self.body[ibody].Ixz>=0 and self.body[ibody].Iyz>=0):
                    TexFRFile.write("Corps rigide S"+str(ibody)+" avec "+str(self.body[ibody].nbrnodes)+" frame(s) \\\\\n")   
                #poutre flexible avec ses frames
                if(self.body[ibody].length>0 and self.body[ibody].section>0 and self.body[ibody].EYoung>0 and self.body[ibody].nu>0 and self.body[ibody].rho>0 and self.body[ibody].Iyy>0 and self.body[ibody].Izz>0 and self.body[ibody].AlphaDamp>=0 and self.body[ibody].BetaDamp>=0 and self.body[ibody].nbrnodes>=2):
                    TexFRFile.write("Poutre flexible S"+str(ibody)+" avec "+str(self.body[ibody].nbrnodes)+" frames \\\\\n")   
                #corps flexible avec ses frames
                if(self.body[ibody].name!=0 and self.body[ibody].nbrnodes>=2 and self.body[ibody].AlphaDamp>=0 and self.body[ibody].BetaDamp>=0 and self.body[ibody].length==0):
                    TexFRFile.write("Corps flexible S"+str(ibody)+" avec "+str(self.body[ibody].nbrnodes)+" frames \\\\\n")

                if(self.body[ibody].mainframe<0):
                    # Rigid body en absolu
                    TexFRFile.write("\\[T_{0F,S"+str(ibody)+"}="+str(latex(self.body[ibody].TtoF))+"\\]\\\\\n\n")
                    for iframe in range(0,self.body[ibody].nbrnodes):
                        # Frame en absolu
                        if(self.body[ibody].frame[iframe].mainframe<0 and self.body[ibody].frame[iframe].frameref<0):
                            TexFRFile.write("\\[T_{0F,S"+str(ibody)+"F"+str(iframe)+"}="+str(latex(self.body[ibody].frame[iframe].TtoF))+"\\]\\\\\n\n")
                        # Flexible sur un solide rigide
                        elif(self.body[ibody].frame[iframe].mainframe>=0 and self.body[ibody].frame[iframe].frameref<0):
                            TexFRFile.write("\\[T_{refF,S"+str(ibody)+"F"+str(iframe)+"/S"+str(self.body[ibody].frame[iframe].mainframe)+"}="+str(latex(self.body[ibody].frame[iframe].TtoF))+"\\]\\\\\n\n")
                        # Flexible sur un solide flexible
                        elif(self.body[ibody].frame[iframe].mainframe>=0 and self.body[ibody].frame[iframe].frameref>=0):
                            TexFRFile.write("\\[T_{refF,S"+str(ibody)+"F"+str(iframe)+"/S"+str(self.body[ibody].frame[iframe].mainframe)+"F"+str(self.body[ibody].frame[iframe].frameref)+"}="+str(latex(self.body[ibody].frame[iframe].TtoF))+"\\]\\\\\n\n")
                # Mainframe du body est en relatif
                else:
                    # un solide rigide par rapport à un solide rigide
                    if (self.body[ibody].frameref<0):
                        TexFRFile.write("\\[T_{refF,S"+str(ibody)+"/S"+str(self.body[ibody].mainframe)+"}="+str(latex(self.body[ibody].TtoF))+"\\]\\\\\n\n")    
                    # un solide rigide par rapport à un solide flexible
                    else:
                        TexFRFile.write("\\[T_{refF,S"+str(ibody)+"/S"+str(self.body[ibody].mainframe)+"F"+str(self.body[ibody].frameref)+"}="+str(latex(self.body[ibody].TtoF))+"\\]\\\\\n\n")
                    for iframe in range(0,self.body[ibody].nbrnodes):
                        # Frame en absolu
                        if(self.body[ibody].frame[iframe].mainframe<0 and self.body[ibody].frame[iframe].frameref<0):
                            TexFRFile.write("\\[T_{0F,S"+str(ibody)+"F"+str(iframe)+"}="+str(latex(self.body[ibody].frame[iframe].TtoF))+"\\]\\\\\n\n")
                        # Flexible sur un solide rigide
                        elif(self.body[ibody].frame[iframe].mainframe>=0 and self.body[ibody].frame[iframe].frameref<0):
                            TexFRFile.write("\\[T_{refF,S"+str(ibody)+"F"+str(iframe)+"/S"+str(self.body[ibody].frame[iframe].mainframe)+"}="+str(latex(self.body[ibody].frame[iframe].TtoF))+"\\]\\\\\n\n")
                        # Flexible sur un solide flexible
                        elif(self.body[ibody].frame[iframe].mainframe>=0 and self.body[ibody].frame[iframe].frameref>=0):
                            TexFRFile.write("\\[T_{refF,S"+str(ibody)+"F"+str(iframe)+"/S"+str(self.body[ibody].frame[iframe].mainframe)+"F"+str(self.body[ibody].frame[iframe].frameref)+"}="+str(latex(self.body[ibody].frame[iframe].TtoF))+"\\]\\\\\n\n")
                      
            TexFRFile.write("\\subsection*{Leur d\\'eriv\\'ee par rapport au temps}\n\n")
            for ibody in range(0,self.nbrbody):
                #corps purement rigide avec frames possibles
                if(self.body[ibody].mass>0 and self.body[ibody].Ixx>0 and self.body[ibody].Iyy>0 and self.body[ibody].Izz>0 and self.body[ibody].Ixy>=0 and self.body[ibody].Ixz>=0 and self.body[ibody].Iyz>=0):
                    TexFRFile.write("Corps rigide S"+str(ibody)+" avec "+str(self.body[ibody].nbrnodes)+" frames \\\\\n")   
                #poutre flexible avec ses frames
                if(self.body[ibody].length>0 and self.body[ibody].section>0 and self.body[ibody].EYoung>0 and self.body[ibody].nu>0 and self.body[ibody].rho>0 and self.body[ibody].Iyy>0 and self.body[ibody].Izz>0 and self.body[ibody].AlphaDamp>=0 and self.body[ibody].BetaDamp>=0 and self.body[ibody].nbrnodes>=2):
                    TexFRFile.write("Poutre flexible S"+str(ibody)+" avec "+str(self.body[ibody].nbrnodes)+" frames \\\\\n")   
                #corps flexible avec ses frames
                if(self.body[ibody].name!=0 and self.body[ibody].nbrnodes>=2 and self.body[ibody].AlphaDamp>=0 and self.body[ibody].BetaDamp>=0 and self.body[ibody].length==0):
                    TexFRFile.write("Corps flexible S"+str(ibody)+" avec "+str(self.body[ibody].nbrnodes)+" frames \\\\\n")

                if(self.body[ibody].mainframe<0):
                    # Rigid body en absolu
                    TexFRFile.write("\\[\\dot{T}_{0F,S"+str(ibody)+"}="+str(latex(self.body[ibody].TtoFd))+"\\]\\\\\n\n")
                    for iframe in range(0,self.body[ibody].nbrnodes):
                        # Frame en absolu
                        if(self.body[ibody].frame[iframe].mainframe<0 and self.body[ibody].frame[iframe].frameref<0):
                            TexFRFile.write("\\[\\dot{T}_{0F,S"+str(ibody)+"F"+str(iframe)+"}="+str(latex(self.body[ibody].frame[iframe].TtoFd))+"\\]\\\\\n\n")
                        # Flexible sur un solide rigide
                        elif(self.body[ibody].frame[iframe].mainframe>=0 and self.body[ibody].frame[iframe].frameref<0):
                            TexFRFile.write("\\[\\dot{T}_{refF,S"+str(ibody)+"F"+str(iframe)+"/S"+str(self.body[ibody].frame[iframe].mainframe)+"}="+str(latex(self.body[ibody].frame[iframe].TtoFd))+"\\]\\\\\n\n")
                        # Flexible sur un solide flexible
                        elif(self.body[ibody].frame[iframe].mainframe>=0 and self.body[ibody].frame[iframe].frameref>=0):
                            TexFRFile.write("\\[\\dot{T}_{refF,S"+str(ibody)+"F"+str(iframe)+"/S"+str(self.body[ibody].frame[iframe].mainframe)+"F"+str(self.body[ibody].frame[iframe].frameref)+"}="+str(latex(self.body[ibody].frame[iframe].TtoFd))+"\\]\\\\\n\n")
                # Mainframe du body est en relatif
                else:
                    # un solide rigide par rapport à un solide rigide
                    if (self.body[ibody].frameref<0):
                        TexFRFile.write("\\[\\dot{T}_{refF,S"+str(ibody)+"/S"+str(self.body[ibody].mainframe)+"}="+str(latex(self.body[ibody].TtoFd))+"\\]\\\\\n\n")    
                    # un solide rigide par rapport à un solide flexible
                    else:
                        TexFRFile.write("\\[\\dot{T}_{refF,S"+str(ibody)+"/S"+str(self.body[ibody].mainframe)+"F"+str(self.body[ibody].frameref)+"}="+str(latex(self.body[ibody].TtoFd))+"\\]\\\\\n\n")
                    for iframe in range(0,self.body[ibody].nbrnodes):
                        # Frame en absolu
                        if(self.body[ibody].frame[iframe].mainframe<0 and self.body[ibody].frame[iframe].frameref<0):
                            TexFRFile.write("\\[\\dot{T}_{0F,S"+str(ibody)+"F"+str(iframe)+"}="+str(latex(self.body[ibody].frame[iframe].TtoFd))+"\\]\\\\\n\n")
                        # Flexible sur un solide rigide
                        elif(self.body[ibody].frame[iframe].mainframe>=0 and self.body[ibody].frame[iframe].frameref<0):
                            TexFRFile.write("\\[\\dot{T}_{refF,S"+str(ibody)+"F"+str(iframe)+"/S"+str(self.body[ibody].frame[iframe].mainframe)+"}="+str(latex(self.body[ibody].frame[iframe].TtoFd))+"\\]\\\\\n\n")
                        # Flexible sur un solide flexible
                        elif(self.body[ibody].frame[iframe].mainframe>=0 and self.body[ibody].frame[iframe].frameref>=0):
                            TexFRFile.write("\\[\\dot{T}_{refF,S"+str(ibody)+"F"+str(iframe)+"/S"+str(self.body[ibody].frame[iframe].mainframe)+"F"+str(self.body[ibody].frame[iframe].frameref)+"}="+str(latex(self.body[ibody].frame[iframe].TtoFd))+"\\]\\\\\n\n")

            TexFRFile.write("\\subsection*{Les vitesses des centres de gravit\\'e de chaque solide et de chaque frame}\n\n")
            for ibody in range(0,self.nbrbody):
                #corps purement rigide avec frames possibles
                if(self.body[ibody].mass>0 and self.body[ibody].Ixx>0 and self.body[ibody].Iyy>0 and self.body[ibody].Izz>0 and self.body[ibody].Ixy>=0 and self.body[ibody].Ixz>=0 and self.body[ibody].Iyz>=0):
                    TexFRFile.write("Corps rigide S"+str(ibody)+" avec "+str(self.body[ibody].nbrnodes)+" frames \\\\\n")   
                #poutre flexible avec ses frames
                if(self.body[ibody].length>0 and self.body[ibody].section>0 and self.body[ibody].EYoung>0 and self.body[ibody].nu>0 and self.body[ibody].rho>0 and self.body[ibody].Iyy>0 and self.body[ibody].Izz>0 and self.body[ibody].AlphaDamp>=0 and self.body[ibody].BetaDamp>=0 and self.body[ibody].nbrnodes>=2):
                    TexFRFile.write("Poutre flexible S"+str(ibody)+" avec "+str(self.body[ibody].nbrnodes)+" frames \\\\\n")   
                #corps flexible avec ses frames
                if(self.body[ibody].name!=0 and self.body[ibody].nbrnodes>=2 and self.body[ibody].AlphaDamp>=0 and self.body[ibody].BetaDamp>=0 and self.body[ibody].length==0):
                    TexFRFile.write("Corps flexible S"+str(ibody)+" avec "+str(self.body[ibody].nbrnodes)+" frames \\\\\n")

                if(self.body[ibody].mainframe<0):
                    # Rigid body en absolu
                    TexFRFile.write("\\[\\vec{v}_{0F,S"+str(ibody)+"}="+str(latex(self.body[ibody].vF))+"\\]\\\\\n\n")
                    for iframe in range(0,self.body[ibody].nbrnodes):
                        # Frame en absolu
                        if(self.body[ibody].frame[iframe].mainframe<0 and self.body[ibody].frame[iframe].frameref<0):
                            TexFRFile.write("\\[\\vec{v}_{0F,S"+str(ibody)+"F"+str(iframe)+"}="+str(latex(self.body[ibody].frame[iframe].vF))+"\\]\\\\\n\n")
                        # Flexible sur un solide rigide
                        elif(self.body[ibody].frame[iframe].mainframe>=0 and self.body[ibody].frame[iframe].frameref<0):
                            TexFRFile.write("\\[\\left\\{\\vec{v}_{refF,S"+str(ibody)+"F"+str(iframe)+"/S"+str(self.body[ibody].frame[iframe].mainframe)+"}\\right\\}_{S"+str(self.body[ibody].frame[iframe].mainframe)+"}="+str(latex(self.body[ibody].frame[iframe].vF))+"\\]\\\\\n\n")
                        # Flexible sur un solide flexible
                        elif(self.body[ibody].frame[iframe].mainframe>=0 and self.body[ibody].frame[iframe].frameref>=0):
                            TexFRFile.write("\\[\\left\\{\\vec{v}_{refF,S"+str(ibody)+"F"+str(iframe)+"/S"+str(self.body[ibody].frame[iframe].mainframe)+"F"+str(self.body[ibody].frame[iframe].frameref)+"}\\right\\}_{S"+str(self.body[ibody].frame[iframe].mainframe)+"F"+str(self.body[ibody].frame[iframe].frameref)+"}="+str(latex(self.body[ibody].frame[iframe].vF))+"\\]\\\\\n\n")
                # Mainframe du body est en relatif
                else:
                    # un solide rigide par rapport à un solide rigide
                    if (self.body[ibody].frameref<0):
                        TexFRFile.write("\\[\\left\\{\\vec{v}_{refF,S"+str(ibody)+"/S"+str(self.body[ibody].mainframe)+"}\\right\\}_{S"+str(self.body[ibody].mainframe)+"}="+str(latex(self.body[ibody].vF))+"\\]\\\\\n\n")    
                    # un solide rigide par rapport à un solide flexible
                    else:
                        TexFRFile.write("\\[\\left\\{\\vec{v}_{refF,S"+str(ibody)+"/S"+str(self.body[ibody].mainframe)+"F"+str(self.body[ibody].frameref)+"}\\right\\}_{S"+str(self.body[ibody].mainframe)+"F"+str(self.body[ibody].frameref)+"}="+str(latex(self.body[ibody].vF))+"\\]\\\\\n\n")
                    for iframe in range(0,self.body[ibody].nbrnodes):
                        # Frame en absolu
                        if(self.body[ibody].frame[iframe].mainframe<0 and self.body[ibody].frame[iframe].frameref<0):
                            TexFRFile.write("\\[\\vec{v}_{0F,S"+str(ibody)+"F"+str(iframe)+"}="+str(latex(self.body[ibody].frame[iframe].vF))+"\\]\\\\\n\n")
                        # Flexible sur un solide rigide
                        elif(self.body[ibody].frame[iframe].mainframe>=0 and self.body[ibody].frame[iframe].frameref<0):
                            TexFRFile.write("\\[\\left\\{\\vec{v}_{refF,S"+str(ibody)+"F"+str(iframe)+"/S"+str(self.body[ibody].frame[iframe].mainframe)+"}\\right\\}_{S"+str(self.body[ibody].frame[iframe].mainframe)+"}="+str(latex(self.body[ibody].frame[iframe].vF))+"\\]\\\\\n\n")
                        # Flexible sur un solide flexible
                        elif(self.body[ibody].frame[iframe].mainframe>=0 and self.body[ibody].frame[iframe].frameref>=0):
                            TexFRFile.write("\\[\\left\\{\\vec{v}_{refF,S"+str(ibody)+"F"+str(iframe)+"/S"+str(self.body[ibody].frame[iframe].mainframe)+"F"+str(self.body[ibody].frame[iframe].frameref)+"}\\right\\}_{S"+str(self.body[ibody].frame[iframe].mainframe)+"F"+str(self.body[ibody].frame[iframe].frameref)+"}="+str(latex(self.body[ibody].frame[iframe].vF))+"\\]\\\\\n\n")

            TexFRFile.write("\\subsection*{Les acc\\'el\\'erations des centres de gravit\\'e de chaque solide et de chaque frame}\n\n")
            for ibody in range(0,self.nbrbody):
                #corps purement rigide avec frames possibles
                if(self.body[ibody].mass>0 and self.body[ibody].Ixx>0 and self.body[ibody].Iyy>0 and self.body[ibody].Izz>0 and self.body[ibody].Ixy>=0 and self.body[ibody].Ixz>=0 and self.body[ibody].Iyz>=0):
                    TexFRFile.write("Corps rigide S"+str(ibody)+" avec "+str(self.body[ibody].nbrnodes)+" frames \\\\\n")   
                #poutre flexible avec ses frames
                if(self.body[ibody].length>0 and self.body[ibody].section>0 and self.body[ibody].EYoung>0 and self.body[ibody].nu>0 and self.body[ibody].rho>0 and self.body[ibody].Iyy>0 and self.body[ibody].Izz>0 and self.body[ibody].AlphaDamp>=0 and self.body[ibody].BetaDamp>=0 and self.body[ibody].nbrnodes>=2):
                    TexFRFile.write("Poutre flexible S"+str(ibody)+" avec "+str(self.body[ibody].nbrnodes)+" frames \\\\\n")   
                #corps flexible avec ses frames
                if(self.body[ibody].name!=0 and self.body[ibody].nbrnodes>=2 and self.body[ibody].AlphaDamp>=0 and self.body[ibody].BetaDamp>=0 and self.body[ibody].length==0):
                    TexFRFile.write("Corps flexible S"+str(ibody)+" avec "+str(self.body[ibody].nbrnodes)+" frames \\\\\n")

                if(self.body[ibody].mainframe<0):
                    # Rigid body en absolu
                    TexFRFile.write("\\[\\vec{a}_{0F,S"+str(ibody)+"}="+str(latex(self.body[ibody].aF))+"\\]\\\\\n\n")
                    for iframe in range(0,self.body[ibody].nbrnodes):
                        # Frame en absolu
                        if(self.body[ibody].frame[iframe].mainframe<0 and self.body[ibody].frame[iframe].frameref<0):
                            TexFRFile.write("\\[\\vec{a}_{0F,S"+str(ibody)+"F"+str(iframe)+"}="+str(latex(self.body[ibody].frame[iframe].aF))+"\\]\\\\\n\n")
                        # Flexible sur un solide rigide
                        elif(self.body[ibody].frame[iframe].mainframe>=0 and self.body[ibody].frame[iframe].frameref<0):
                            TexFRFile.write("\\[\\left\\{\\vec{a}_{refF,S"+str(ibody)+"F"+str(iframe)+"/S"+str(self.body[ibody].frame[iframe].mainframe)+"}\\right\\}_{S"+str(self.body[ibody].frame[iframe].mainframe)+"}="+str(latex(self.body[ibody].frame[iframe].aF))+"\\]\\\\\n\n")
                        # Flexible sur un solide flexible
                        elif(self.body[ibody].frame[iframe].mainframe>=0 and self.body[ibody].frame[iframe].frameref>=0):
                            TexFRFile.write("\\[\\left\\{\\vec{a}_{refF,S"+str(ibody)+"F"+str(iframe)+"/S"+str(self.body[ibody].frame[iframe].mainframe)+"F"+str(self.body[ibody].frame[iframe].frameref)+"}\\right\\}_{S"+str(self.body[ibody].frame[iframe].mainframe)+"F"+str(self.body[ibody].frame[iframe].frameref)+"}="+str(latex(self.body[ibody].frame[iframe].aF))+"\\]\\\\\n\n")
                # Mainframe du body est en relatif
                else:
                    # un solide rigide par rapport à un solide rigide
                    if (self.body[ibody].frameref<0):
                        TexFRFile.write("\\[\\left\\{\\vec{a}_{refF,S"+str(ibody)+"/S"+str(self.body[ibody].mainframe)+"}\\right\\}_{S"+str(self.body[ibody].mainframe)+"}="+str(latex(self.body[ibody].aF))+"\\]\\\\\n\n")    
                    # un solide rigide par rapport à un solide flexible
                    else:
                        TexFRFile.write("\\[\\left\\{\\vec{a}_{refF,S"+str(ibody)+"/S"+str(self.body[ibody].mainframe)+"F"+str(self.body[ibody].frameref)+"}\\right\\}_{S"+str(self.body[ibody].mainframe)+"F"+str(self.body[ibody].frameref)+"}="+str(latex(self.body[ibody].aF))+"\\]\\\\\n\n")
                    for iframe in range(0,self.body[ibody].nbrnodes):
                        # Frame en absolu
                        if(self.body[ibody].frame[iframe].mainframe<0 and self.body[ibody].frame[iframe].frameref<0):
                            TexFRFile.write("\\[\\vec{a}_{0F,S"+str(ibody)+"F"+str(iframe)+"}="+str(latex(self.body[ibody].frame[iframe].aF))+"\\]\\\\\n\n")
                        # Flexible sur un solide rigide
                        elif(self.body[ibody].frame[iframe].mainframe>=0 and self.body[ibody].frame[iframe].frameref<0):
                            TexFRFile.write("\\[\\left\\{\\vec{a}_{refF,S"+str(ibody)+"F"+str(iframe)+"/S"+str(self.body[ibody].frame[iframe].mainframe)+"}\\right\\}_{S"+str(self.body[ibody].frame[iframe].mainframe)+"}="+str(latex(self.body[ibody].frame[iframe].aF))+"\\]\\\\\n\n")
                        # Flexible sur un solide flexible
                        elif(self.body[ibody].frame[iframe].mainframe>=0 and self.body[ibody].frame[iframe].frameref>=0):
                            TexFRFile.write("\\[\\left\\{\\vec{a}_{refF,S"+str(ibody)+"F"+str(iframe)+"/S"+str(self.body[ibody].frame[iframe].mainframe)+"F"+str(self.body[ibody].frame[iframe].frameref)+"}\\right\\}_{S"+str(self.body[ibody].frame[iframe].mainframe)+"F"+str(self.body[ibody].frame[iframe].frameref)+"}="+str(latex(self.body[ibody].frame[iframe].aF))+"\\]\\\\\n\n")
            TexFRFile.write("\\subsection*{Les vitesses de rotation de chaque solide et de chaque frame}\n\n")
            for ibody in range(0,self.nbrbody):
                #corps purement rigide avec frames possibles
                if(self.body[ibody].mass>0 and self.body[ibody].Ixx>0 and self.body[ibody].Iyy>0 and self.body[ibody].Izz>0 and self.body[ibody].Ixy>=0 and self.body[ibody].Ixz>=0 and self.body[ibody].Iyz>=0):
                    TexFRFile.write("Corps rigide S"+str(ibody)+" avec "+str(self.body[ibody].nbrnodes)+" frames \\\\\n")   
                #poutre flexible avec ses frames
                if(self.body[ibody].length>0 and self.body[ibody].section>0 and self.body[ibody].EYoung>0 and self.body[ibody].nu>0 and self.body[ibody].rho>0 and self.body[ibody].Iyy>0 and self.body[ibody].Izz>0 and self.body[ibody].AlphaDamp>=0 and self.body[ibody].BetaDamp>=0 and self.body[ibody].nbrnodes>=2):
                    TexFRFile.write("Poutre flexible S"+str(ibody)+" avec "+str(self.body[ibody].nbrnodes)+" frames \\\\\n")   
                #corps flexible avec ses frames
                if(self.body[ibody].name!=0 and self.body[ibody].nbrnodes>=2 and self.body[ibody].AlphaDamp>=0 and self.body[ibody].BetaDamp>=0 and self.body[ibody].length==0):
                    TexFRFile.write("Corps flexible S"+str(ibody)+" avec "+str(self.body[ibody].nbrnodes)+" frames \\\\\n")

                if(self.body[ibody].mainframe<0):
                    # Rigid body en absolu
                    TexFRFile.write("\\[\\vec{\omega}_{0F,S"+str(ibody)+"}="+str(latex(self.body[ibody].omega))+"\\]\\\\\n\n")
                    for iframe in range(0,self.body[ibody].nbrnodes):
                        # Frame en absolu
                        if(self.body[ibody].frame[iframe].mainframe<0 and self.body[ibody].frame[iframe].frameref<0):
                            TexFRFile.write("\\[\\vec{\omega}_{0F,S"+str(ibody)+"F"+str(iframe)+"}="+str(latex(self.body[ibody].frame[iframe].omega))+"\\]\\\\\n\n")
                        # Flexible sur un solide rigide
                        elif(self.body[ibody].frame[iframe].mainframe>=0 and self.body[ibody].frame[iframe].frameref<0):
                            TexFRFile.write("\\[\\left\\{\\vec{\omega}_{refF,S"+str(ibody)+"F"+str(iframe)+"/S"+str(self.body[ibody].frame[iframe].mainframe)+"}\\right\\}_{S"+str(self.body[ibody].frame[iframe].mainframe)+"}="+str(latex(self.body[ibody].frame[iframe].omega))+"\\]\\\\\n\n")
                        # Flexible sur un solide flexible
                        elif(self.body[ibody].frame[iframe].mainframe>=0 and self.body[ibody].frame[iframe].frameref>=0):
                            TexFRFile.write("\\[\\left\\{\\vec{\omega}_{refF,S"+str(ibody)+"F"+str(iframe)+"/S"+str(self.body[ibody].frame[iframe].mainframe)+"F"+str(self.body[ibody].frame[iframe].frameref)+"}\\right\\}_{S"+str(self.body[ibody].frame[iframe].mainframe)+"F"+str(self.body[ibody].frame[iframe].frameref)+"}="+str(latex(self.body[ibody].frame[iframe].omega))+"\\]\\\\\n\n")
                # Mainframe du body est en relatif
                else:
                    # un solide rigide par rapport à un solide rigide
                    if (self.body[ibody].frameref<0):
                        TexFRFile.write("\\[\\left\\{\\vec{\omega}_{refF,S"+str(ibody)+"/S"+str(self.body[ibody].mainframe)+"}\\right\\}_{S"+str(self.body[ibody].mainframe)+"}="+str(latex(self.body[ibody].omega))+"\\]\\\\\n\n")    
                    # un solide rigide par rapport à un solide flexible
                    else:
                        TexFRFile.write("\\[\\left\\{\\vec{\omega}_{refF,S"+str(ibody)+"/S"+str(self.body[ibody].mainframe)+"F"+str(self.body[ibody].frameref)+"}\\right\\}_{S"+str(self.body[ibody].mainframe)+"F"+str(self.body[ibody].frameref)+"}="+str(latex(self.body[ibody].omega))+"\\]\\\\\n\n")
                    for iframe in range(0,self.body[ibody].nbrnodes):
                        # Frame en absolu
                        if(self.body[ibody].frame[iframe].mainframe<0 and self.body[ibody].frame[iframe].frameref<0):
                            TexFRFile.write("\\[\\vec{\omega}_{0F,S"+str(ibody)+"F"+str(iframe)+"}="+str(latex(self.body[ibody].frame[iframe].omega))+"\\]\\\\\n\n")
                        # Flexible sur un solide rigide
                        elif(self.body[ibody].frame[iframe].mainframe>=0 and self.body[ibody].frame[iframe].frameref<0):
                            TexFRFile.write("\\[\\left\\{\\vec{\omega}_{refF,S"+str(ibody)+"F"+str(iframe)+"/S"+str(self.body[ibody].frame[iframe].mainframe)+"}\\right\\}_{S"+str(self.body[ibody].frame[iframe].mainframe)+"}="+str(latex(self.body[ibody].frame[iframe].omega))+"\\]\\\\\n\n")
                        # Flexible sur un solide flexible
                        elif(self.body[ibody].frame[iframe].mainframe>=0 and self.body[ibody].frame[iframe].frameref>=0):
                            TexFRFile.write("\\[\\left\\{\\vec{\omega}_{refF,S"+str(ibody)+"F"+str(iframe)+"/S"+str(self.body[ibody].frame[iframe].mainframe)+"F"+str(self.body[ibody].frame[iframe].frameref)+"}\\right\\}_{S"+str(self.body[ibody].frame[iframe].mainframe)+"F"+str(self.body[ibody].frame[iframe].frameref)+"}="+str(latex(self.body[ibody].frame[iframe].omega))+"\\]\\\\\n\n")
            TexFRFile.write("\\subsection*{Les acc\\'el\\'erations de rotation de chaque solide et de chaque frame}\n\n")
            for ibody in range(0,self.nbrbody):
                #corps purement rigide avec frames possibles
                if(self.body[ibody].mass>0 and self.body[ibody].Ixx>0 and self.body[ibody].Iyy>0 and self.body[ibody].Izz>0 and self.body[ibody].Ixy>=0 and self.body[ibody].Ixz>=0 and self.body[ibody].Iyz>=0):
                    TexFRFile.write("Corps rigide S"+str(ibody)+" avec "+str(self.body[ibody].nbrnodes)+" frames \\\\\n")   
                #poutre flexible avec ses frames
                if(self.body[ibody].length>0 and self.body[ibody].section>0 and self.body[ibody].EYoung>0 and self.body[ibody].nu>0 and self.body[ibody].rho>0 and self.body[ibody].Iyy>0 and self.body[ibody].Izz>0 and self.body[ibody].AlphaDamp>=0 and self.body[ibody].BetaDamp>=0 and self.body[ibody].nbrnodes>=2):
                    TexFRFile.write("Poutre flexible S"+str(ibody)+" avec "+str(self.body[ibody].nbrnodes)+" frames \\\\\n")   
                #corps flexible avec ses frames
                if(self.body[ibody].name!=0 and self.body[ibody].nbrnodes>=2 and self.body[ibody].AlphaDamp>=0 and self.body[ibody].BetaDamp>=0 and self.body[ibody].length==0):
                    TexFRFile.write("Corps flexible S"+str(ibody)+" avec "+str(self.body[ibody].nbrnodes)+" frames \\\\\n")

                if(self.body[ibody].mainframe<0):
                    # Rigid body en absolu
                    TexFRFile.write("\\[\\vec{\dot{\omega}}_{0F,S"+str(ibody)+"}="+str(latex(self.body[ibody].omegad))+"\\]\\\\\n\n")
                    for iframe in range(0,self.body[ibody].nbrnodes):
                        # Frame en absolu
                        if(self.body[ibody].frame[iframe].mainframe<0 and self.body[ibody].frame[iframe].frameref<0):
                            TexFRFile.write("\\[\\vec{\dot{\omega}}_{0F,S"+str(ibody)+"F"+str(iframe)+"}="+str(latex(self.body[ibody].frame[iframe].omegad))+"\\]\\\\\n\n")
                        # Flexible sur un solide rigide
                        elif(self.body[ibody].frame[iframe].mainframe>=0 and self.body[ibody].frame[iframe].frameref<0):
                            TexFRFile.write("\\[\\left\\{\\vec{\dot{\omega}}_{refF,S"+str(ibody)+"F"+str(iframe)+"/S"+str(self.body[ibody].frame[iframe].mainframe)+"}\\right\\}_{S"+str(self.body[ibody].frame[iframe].mainframe)+"}="+str(latex(self.body[ibody].frame[iframe].omegad))+"\\]\\\\\n\n")
                        # Flexible sur un solide flexible
                        elif(self.body[ibody].frame[iframe].mainframe>=0 and self.body[ibody].frame[iframe].frameref>=0):
                            TexFRFile.write("\\[\\left\\{\\vec{\dot{\omega}}_{refF,S"+str(ibody)+"F"+str(iframe)+"/S"+str(self.body[ibody].frame[iframe].mainframe)+"F"+str(self.body[ibody].frame[iframe].frameref)+"}\\right\\}_{S"+str(self.body[ibody].frame[iframe].mainframe)+"F"+str(self.body[ibody].frame[iframe].frameref)+"}="+str(latex(self.body[ibody].frame[iframe].omegad))+"\\]\\\\\n\n")
                # Mainframe du body est en relatif
                else:
                    # un solide rigide par rapport à un solide rigide
                    if (self.body[ibody].frameref<0):
                        TexFRFile.write("\\[\\left\\{\\vec{\dot{\omega}}_{refF,S"+str(ibody)+"/S"+str(self.body[ibody].mainframe)+"}\\right\\}_{S"+str(self.body[ibody].mainframe)+"}="+str(latex(self.body[ibody].omegad))+"\\]\\\\\n\n")    
                    # un solide rigide par rapport à un solide flexible
                    else:
                        TexFRFile.write("\\[\\left\\{\\vec{\dot{\omega}}_{refF,S"+str(ibody)+"/S"+str(self.body[ibody].mainframe)+"F"+str(self.body[ibody].frameref)+"}\\right\\}_{S"+str(self.body[ibody].mainframe)+"F"+str(self.body[ibody].frameref)+"}="+str(latex(self.body[ibody].omegad))+"\\]\\\\\n\n")
                    for iframe in range(0,self.body[ibody].nbrnodes):
                        # Frame en absolu
                        if(self.body[ibody].frame[iframe].mainframe<0 and self.body[ibody].frame[iframe].frameref<0):
                            TexFRFile.write("\\[\\vec{\dot{\omega}}_{0F,S"+str(ibody)+"F"+str(iframe)+"}="+str(latex(self.body[ibody].frame[iframe].omegad))+"\\]\\\\\n\n")
                        # Flexible sur un solide rigide
                        elif(self.body[ibody].frame[iframe].mainframe>=0 and self.body[ibody].frame[iframe].frameref<0):
                            TexFRFile.write("\\[\\left\\{\\vec{\dot{\omega}}_{refF,S"+str(ibody)+"F"+str(iframe)+"/S"+str(self.body[ibody].frame[iframe].mainframe)+"}\\right\\}_{S"+str(self.body[ibody].frame[iframe].mainframe)+"}="+str(latex(self.body[ibody].frame[iframe].omegad))+"\\]\\\\\n\n")
                        # Flexible sur un solide flexible
                        elif(self.body[ibody].frame[iframe].mainframe>=0 and self.body[ibody].frame[iframe].frameref>=0):
                            TexFRFile.write("\\[\\left\\{\\vec{\dot{\omega}}_{refF,S"+str(ibody)+"F"+str(iframe)+"/S"+str(self.body[ibody].frame[iframe].mainframe)+"F"+str(self.body[ibody].frame[iframe].frameref)+"}\\right\\}_{S"+str(self.body[ibody].frame[iframe].mainframe)+"F"+str(self.body[ibody].frame[iframe].frameref)+"}="+str(latex(self.body[ibody].frame[iframe].omegad))+"\\]\\\\\n\n") 

        #if FORCES==1:
        if(self.ForcesInAppEff != " "):
            TexFRFile.write("\n\\section{D\\'efinition des forces externes}\n\n")
            TexFRFile.write("Les forces externes ont \\'et\\'e d\\'efinies dans la ")
            if(self.EASYDYNPlusPlusFlag==0):
                TexFRFile.write("proc\\'edure \\texttt{AddAppliedEfforts()}. ")
            else:
                TexFRFile.write("proc\\'edure \\texttt{AddAppliedForces()}. ")
        
        TexFRFile.write("\n\\section{Simulation du m\\'ecanisme}\n\n")
        TexFRFile.write("La simulation temporelle s'est faite via la proc\\'edure \\texttt{NewmarkIntegration}")
        TexFRFile.write("(\\emph{FinalTime} , \\emph{StepSave} , \\emph{StepMax}) du fichier \\texttt{"+str(self.ApplicationFileName)+".cpp} o\\`u :\n")
        TexFRFile.write("\\begin{itemize}\n")
        TexFRFile.write("   \\item \\emph{FinalTime} est la dur\\'ee de la simulation (="+str(self.tfinal)+" $s$),\n")
        TexFRFile.write("   \\item \\emph{StepSave} est le pas de temps dans l'int\\'egration num\\'erique (="+str(self.hsave)+" $s$),\n")
        TexFRFile.write("   \\item \\emph{StepMax} est un pas de temps limite (="+str(self.hmax)+" $s$),\n")
        TexFRFile.write("\\end{itemize}\n")
        TexFRFile.write("\\noindent les conditions initiales ")
        InitialCondFlag=0 # y a t il des conditions initiales
        for idof in range(0,self.nbrdof):
            if ((abs(self.qini[idof])>0) or (abs(self.qdini[idof])>0)):
                InitialCondFlag=InitialCondFlag+1
        if InitialCondFlag>0:
                TexFRFile.write("\\'etant ")
        for idof in range(0,self.nbrdof):
            if self.qini[idof] !=0:
                TexFRFile.write("$q_{"+str(idof)+"} = "+str(self.qini[idof])+"$, ")
            if self.qdini[idof] !=0:
                TexFRFile.write("$\dot{q}_{"+str(idof)+"} = "+str(self.qdini[idof])+"$, ")
        if InitialCondFlag==0:
            TexFRFile.write(" \\'etant toutes nulles.\\\\\n")
        else:    
            TexFRFile.write("les autres \\'etant nulles.\\\\\n")
        TexFRFile.write("\n\\section{R\\'esultats}\n\n")
        TexFRFile.write("L'affichage des \\'evolutions temporelles des diff\\'erents param\\`etres de configuration ainsi que leurs ")
        TexFRFile.write("d\\'eriv\\'ees premi\\`ere et seconde s'\\'etablit assez facilement sous \\textsf{Gnuplot} comme nous le montrent ")
        TexFRFile.write("les figures \\ref{figure1} \\`a \\ref{figure3} cr\\'e\\'ees par le code suivant:\\\\\n\n")
  
        TexFRFile.write("\\footnotesize\n")
        TexFRFile.write("\\verbatiminput{"+str(self.ApplicationFileName)+".plt}\n")
        TexFRFile.write("\\normalsize\n\n")
        
        TexFRFile.write("\\begin{figure}[h!tb]\n")
        TexFRFile.write("   \\begin{center}\n")
        TexFRFile.write("      \\includegraphics[width=0.75\\textwidth]{figure1.eps}\n")
        TexFRFile.write("   \\end{center}\n")
        TexFRFile.write("\\caption{\\'Evolution temporelle de(s) param\\`etre(s) de configuration}\n")
        TexFRFile.write("\\label{figure1}\n")
        TexFRFile.write("\\end{figure}\n\n")

        TexFRFile.write("\\begin{figure}[h!tb]\n")
        TexFRFile.write("   \\begin{center}\n")
        TexFRFile.write("      \\includegraphics[width=0.75\\textwidth]{figure2.eps}\n")
        TexFRFile.write("   \\end{center}\n")
        TexFRFile.write("\\caption{D\\'eriv\\'ee(s) premi\\`ere(s) de(s) param\\`etre(s) de configuration}\n")
        TexFRFile.write("\\label{figure2}\n")
        TexFRFile.write("\\end{figure}\n\n")

        TexFRFile.write("\\begin{figure}[h!tb]\n")
        TexFRFile.write("   \\begin{center}\n")
        TexFRFile.write("      \\includegraphics[width=0.75\\textwidth]{figure3.eps}\n")
        TexFRFile.write("   \\end{center}\n")
        TexFRFile.write("\\caption{D\\'eriv\\'ee(s) seconde(s) de(s) param\\`etre(s) de configuration}\n")
        TexFRFile.write("\\label{figure3}\n")
        TexFRFile.write("\\end{figure}\n\n")
                
        TexFRFile.write("\n\\appendix\n\n")
        TexFRFile.write("\\clearpage\n")
        TexFRFile.write("\\section{Listing du fichier Python}\n")
        TexFRFile.write("\\footnotesize\n")
        TexFRFile.write("\\verbatiminput{"+str(self.ApplicationFileName)+".py}\n\n")
        TexFRFile.write("\\end{document}\n")  
        
        TexFRFile.close()
    
        # Post-process of FR Latex Code to replace qd by \dot{}
        filedataFR = None
        with open('Rapport_'+str(self.ApplicationFileName)+'.tex', 'r') as file :
            filedataFR = file.read()
        # Replace the target string
        for idof in range(0,self.nbrdof):
            filedataFR = filedataFR.replace('q['+str(idof)+']', 'q_'+str(idof))
            filedataFR = filedataFR.replace('qd['+str(idof)+']', '\dot{q}_'+str(idof))
            filedataFR = filedataFR.replace('qdd['+str(idof)+']', '\ddot{q}_'+str(idof))
            filedataFR = filedataFR.replace('qddd['+str(idof)+']', '\dddot{q}_'+str(idof))
            filedataFR = filedataFR.replace('qdddd['+str(idof)+']', '\ddddot{q}_'+str(idof))
        for idof in range(0,self.nbrdep):
            filedataFR = filedataFR.replace('p['+str(idof)+']', 'p_'+str(idof))
            filedataFR = filedataFR.replace('pd['+str(idof)+']', '\dot{p}_'+str(idof))
            filedataFR = filedataFR.replace('pdd['+str(idof)+']', '\ddot{p}_'+str(idof))
            filedataFR = filedataFR.replace('pddd['+str(idof)+']', '\dddot{p}_'+str(idof))
            filedataFR = filedataFR.replace('pdddd['+str(idof)+']', '\ddddot{p}_'+str(idof))
        filedataFR = filedataFR.replace(' pi', ' \pi')
        # Write the file out again
        with open('Rapport_'+str(self.ApplicationFileName)+'.tex', 'w') as file:
            file.write(filedataFR)
            
    # ------ end ExportFR_Latex_Report
    # ------ begin ExportUK_Latex_Report
    def ExportUK_Latex_Report(self):
        print("Processing EN Latex files") 
        #.tex file to Latex File
        TexENFileName="Report_"+(self.ApplicationFileName)+".tex"
        print("Creating 'English' Tex File "+TexENFileName)
        TexENFile = open(TexENFileName, "w")
        
        TexENFile.write("\\documentclass[a4paper,11pt]{article}\n")
        TexENFile.write("\\topmargin -0.5cm\n")
        TexENFile.write("\\headheight 0.7cm\n")
        TexENFile.write("\\headsep    0.8cm\n")
        TexENFile.write("\\topskip      0cm\n")
        TexENFile.write("\\textheight  24cm\n")
        TexENFile.write("\\oddsidemargin  0cm\n")
        TexENFile.write("\\evensidemargin 0cm\n")
        TexENFile.write("\\textwidth   16cm\n")
        TexENFile.write("\\usepackage{verbatim}\n")
        TexENFile.write("\\usepackage{array} %>{$}m{3cm}<{$}\n")
        TexENFile.write("\\usepackage{amssymb}\n")
        TexENFile.write("\\usepackage{amsmath}\n")
        TexENFile.write("\\usepackage[dvips]{graphicx}\n")
        TexENFile.write("\\usepackage[T1]{fontenc}\n")
        TexENFile.write("\\usepackage[english]{babel}\n")
        TexENFile.write("\\usepackage{fancyhdr}\n")
        TexENFile.write("\\pagestyle{fancy}\n")
        TexENFile.write("\\fancyhead{}\n")
        TexENFile.write("\\fancyfoot{}\n")
        TexENFile.write("\\fancyhead[L]{"+self.ApplicationTitle+"}\n")
        TexENFile.write("\\fancyhead[R]{\\thepage}\n")
        TexENFile.write("\\title{\\textsf{"+self.ApplicationTitle+"}}\n")
        TexENFile.write("\\author{}\n")
        TexENFile.write("\\begin{document}\n")
        TexENFile.write("\\sloppy\n")
        TexENFile.write("\\maketitle\n")
        TexENFile.write("\n")
        TexENFile.write("\\section{General data of the studied mechanism}\n")
        
        # Cas dans lequel il n'y a pas de corps flexibles => cas purement rigide
        if (self.nbrFlexBeam+self.nbrFlexBody)==0: 
            TexENFile.write("\nThe system comprises "+str(self.nbrbody)+" bodies (defined by the global variable \\texttt{nbrbody}). ")
            TexENFile.write("Each body is called $S_j$ ($j$ from 0 to "+str(self.nbrbody-1)+"). ")
            if self.nbrdof >1:
                TexENFile.write("The number of degrees of freedom of the system is "+str(self.nbrdof)+" (\\texttt{nbrdof}). ")
                TexENFile.write("The configuration parameters are denoted by $q_i$ ($i$ from 0 to "+str(self.nbrdof-1)+").\\\\\n\n")
            if self.nbrdof ==1:
                TexENFile.write("The studied mechanism is describe by one degree of freedom and so ")
                TexENFile.write("only one configuration parameter is defined by $q_0$.\\\\\n\n")
            TexENFile.write("The inertial data, given by the user, consist of the mass $m_{Si}$ and ")
            TexENFile.write("the inertia tensor $\Phi_{G,Si}$ of each body $i$ expressed with respect to the center of gravity.\\\\\n\n")
        # Cas général: avec corps rigides, flexibles et poutres flexibles
        else:  
            TexENFile.write("\nThe system comprises "+str(self.nbrbody)+" bodies (defined by the global variable \\texttt{nbrbody}): ")
            TexENFile.write("\n\\begin{itemize}\n")
            if (self.nbrRigidBody)>0:
                if(self.nbrRigidBody==1):
                    TexENFile.write("  \\item There is "+str(self.nbrRigidBody)+" rigid body.\n")
                else:
                    TexENFile.write("  \\item There are "+str(self.nbrRigidBody)+" rigid bodies.\n")
            if (self.nbrFlexBeam)>0:
                if(self.nbrFlexBeam==1):
                    TexENFile.write("  \\item There is "+str(self.nbrFlexBeam)+" flexible beam.\n")
                else:
                    TexENFile.write("  \\item There are "+str(self.nbrFlexBeam)+" flexible beams.\n")
            if (self.nbrFlexBody)>0:
                if(self.nbrFlexBody==1):
                    TexENFile.write("  \\item There is "+str(self.nbrFlexBody)+" flexible body\n")
                else:
                    TexENFile.write("  \\item There are "+str(self.nbrFlexBody)+" flexible bodies\n")
            TexENFile.write("\\end{itemize}\n")
            if self.nbrdof >1:
                TexENFile.write("The number of degrees of freedom of the system is "+str(self.nbrdof)+" (\\texttt{nbrdof}). ")
                TexENFile.write("The configuration parameters are denoted by $q_i$ ($i$ from 0 to "+str(self.nbrdof-1)+").\\\\\n\n")
            if self.nbrdof ==1:
                TexENFile.write("The studied mechanism is describe by one degree of freedom and so ")
                TexENFile.write("only one configuration parameter is defined by $q_0$.\\\\\n\n")
                
        # Cas dans lequel il n'y a pas de corps flexibles => cas purement rigide
        if (self.nbrFlexBeam+self.nbrFlexBody)==0:
            for ibody in range(0,self.nbrbody):
                TexENFile.write("\\[ m_{S"+str(ibody)+"}="+str(latex(self.body[ibody].mass))+"\\,kg \\]\n")
            TexENFile.write("\\\\\n")
        # Cas général: avec corps rigides, flexibles et poutres flexibles
        else:  
            if (self.nbrRigidBody)>0:
                TexENFile.write("\\subsection*{Rigid bodie(s): inertia data}\n\n")
                for ibody in range(0,self.nbrbody):
                    if(self.body[ibody].mass>0 and self.body[ibody].Ixx>0 and self.body[ibody].Iyy>0 and self.body[ibody].Izz>0 and self.body[ibody].Ixy>=0 and self.body[ibody].Ixz>=0 and self.body[ibody].Iyz>=0):
                        TexENFile.write("Rigid body S"+str(ibody)+": ")
                        TexENFile.write("\\[ m_{S"+str(ibody)+"}="+str(latex(self.body[ibody].mass))+"\\,kg \\]\n")
                        TexENFile.write("\\[\\Phi_{G,S"+str(ibody)+"}="
                        +str(latex(Matrix(((self.body[ibody].Ixx,-self.body[ibody].Ixy,-self.body[ibody].Ixz),
                        (-self.body[ibody].Ixy,self.body[ibody].Iyy,-self.body[ibody].Ixy),
                        (-self.body[ibody].Ixz,-self.body[ibody].Iyz,self.body[ibody].Izz)))))+" \\quad \\textrm{, en $kg.m^2$}\\]\\\\\n")
            # si c'est une poutre flexible: afficher ses propriétés
            if (self.nbrFlexBeam)>0:
                TexENFile.write("\\subsection*{Flexible beam(s): propertie(s)}\n\n")
                for ibody in range(0,self.nbrbody):
                    # si c'est une poutre flexible: afficher ses propriétés
                    if(self.body[ibody].length>0 and self.body[ibody].section>0 and self.body[ibody].EYoung>0 and self.body[ibody].nu>0 and self.body[ibody].rho>0 and self.body[ibody].Iyy>0 and self.body[ibody].Izz>0 and self.body[ibody].AlphaDamp>=0 and self.body[ibody].BetaDamp>=0 and self.body[ibody].nbrnodes>=2):
                        TexENFile.write("Flexible beam S"+str(ibody)+": ")
                        TexENFile.write("\n\\begin{itemize}\n")
                        TexENFile.write("  \\item Length: "+str(self.body[ibody].length)+" m;\n")
                        TexENFile.write("  \\item Section: "+str(self.body[ibody].section)+" $m^2$;\n")
                        TexENFile.write("  \\item Young's modulus: "+str(self.body[ibody].EYoung)+" Pa;\n")
                        TexENFile.write("  \\item Poisson's ratio: "+str(self.body[ibody].nu)+";\n")
                        TexENFile.write("  \\item Density: "+str(self.body[ibody].rho)+" $kg/m^3$;\n")
                        TexENFile.write("  \\item Geometric moment of inertia $I_{yy}$: "+str(self.body[ibody].Iyy)+" $m^4$;\n")
                        TexENFile.write("  \\item Geometric moment of inertia $I_{zz}$: "+str(self.body[ibody].Izz)+" $m^4$;\n")
                        TexENFile.write("  \\item Number of frames: "+str(self.body[ibody].nbrnodes)+".\\\\\n")
                        TexENFile.write("\\end{itemize}\n")
            if (self.nbrFlexBody)>0:
                TexENFile.write("\\subsection*{Flexible bodi(s): mass and stiffness matrices}\n\n")
                for ibody in range(0,self.nbrbody):
                    if(self.body[ibody].name!=0 and self.body[ibody].nbrnodes>=2 and self.body[ibody].AlphaDamp>=0 and self.body[ibody].BetaDamp>=0 and self.body[ibody].length==0):
                        TexENFile.write("Flexible body S"+str(ibody)+": \\\\\n")
                        # open Mass Matrix file
                        with open(str(self.body[ibody].name)+".mm") as f:
                            MatriceMM = []
                            for line in f:
                                line = line.split() # to deal with blank 
                                if line:            # lines (ie skip them)
                                    line = [i for i in line]
                                    MatriceMM.append(line) 
                                    MM=Matrix(MatriceMM)
                        TexENFile.write("\\[M_{"+str(ibody)+"}="+str(latex(MM))+"\\]\\\\\n\n")
                        # open Stiffness Matrix file
                        with open(str(self.body[ibody].name)+".kk") as f:
                            MatriceKK = []
                            for line in f:
                                line = line.split() # to deal with blank 
                                if line:            # lines (ie skip them)
                                    line = [i for i in line]
                                    MatriceKK.append(line) 
                                    KK=Matrix(MatriceKK)
                        TexENFile.write("\\[K_{"+str(ibody)+"}="+str(latex(KK))+"\\]\\\\\n\n")
        
        # Cas dans lequel il n'y a pas de corps flexibles => cas purement rigide
        if (self.nbrFlexBeam+self.nbrFlexBody)==0:
            for ibody in range(0,self.nbrbody):
                TexENFile.write("\\[\\Phi_{G,S"+str(ibody)+"}="
                +str(latex(Matrix(((self.body[ibody].Ixx,-self.body[ibody].Ixy,-self.body[ibody].Ixz),
                (-self.body[ibody].Ixy,self.body[ibody].Iyy,-self.body[ibody].Ixy),
                (-self.body[ibody].Ixz,-self.body[ibody].Iyz,self.body[ibody].Izz)))))+" \\quad \\textrm{, en $kg.m^2$}\\]\\\\\n")
        TexENFile.write("\n\\section{Complete kinematics calculated by Sympy}\n\n")
        TexENFile.write("The following parameters have been calculated from the user's file \\texttt{"+str(self.ApplicationFileName)+".py} with \n")
        TexENFile.write("a \\emph{CPU} time of "+str(round(self.EndTime_Kinematics - self.StartTime_Kinematics,3))+" second(s). \\\\\n")
        TexENFile.write("\n")
        
        # Cas dans lequel il n'y a pas de corps flexibles => cas purement rigide
        if (self.nbrFlexBeam+self.nbrFlexBody)==0: 
            NbRelBodies=0 # Compute the number of relative bodies
            for ibody in range(0,self.nbrbody):
                if self.body[ibody].mainframe>=0:
                    NbRelBodies = NbRelBodies+1
            if NbRelBodies>1:
                TexENFile.write("\\subsection*{Relative motions} \n\n")
                TexENFile.write("The motion of some bodies has been defined as a relative motion with respect to another body.")
                TexENFile.write(" It is the case of "+str(NbRelBodies)+" bodies, for which the motion is defined in the following manner :\n")
                TexENFile.write("\\begin{itemize}\n")
                for ibody in range(0,self.nbrbody):
                    if self.body[ibody].mainframe>=0:
                        TexENFile.write("  \\item Motion of body $S_{"+str(ibody)+"}$ is given with respect to body $S_{"+str(self.body[ibody].mainframe)+"}$.\n")
                TexENFile.write("\\end{itemize}\n")
            else:
                for ibody in range(0,self.nbrbody):
                    if self.body[ibody].mainframe>=0:
                        TexENFile.write("\\noindent \\emph{Note} : The motion of body $S_{"+str(ibody)+"}$ has been defined as a relative motion with respect to body $S_{"+str(self.body[ibody].mainframe)+"}$. \n") 
        else:
            Jar=0
            for ibody in range(0,self.nbrbody):
                if(self.body[ibody].mainframe>=0 or self.body[ibody].frameref>=0):
                    Jar=Jar+1
                for iframe in range(0,self.body[ibody].nbrnodes):
                    if(self.body[ibody].frame[iframe].mainframe>=0 or self.body[ibody].frame[iframe].frameref>=0):
                        Jar=Jar+1
            if Jar>0: # S'il y a des mouvements relatifs
                TexENFile.write("\\subsection*{Relative motions} \n\n")
                for ibody in range(0,self.nbrbody):
                    #corps purement rigide avec frames possibles
                    if(self.body[ibody].mass>0 and self.body[ibody].Ixx>0 and self.body[ibody].Iyy>0 and self.body[ibody].Izz>0 and self.body[ibody].Ixy>=0 and self.body[ibody].Ixz>=0 and self.body[ibody].Iyz>=0):
                        DisplayOK=0
                        if(self.body[ibody].mainframe>=0 or self.body[ibody].frameref>=0):
                            #Si mouvement relatif, on affiche
                            DisplayOK=1 
                        PassDisplayOK=0
                        for iframe in range(0,self.body[ibody].nbrnodes):
                            #Si mouvement relatif, on affiche
                            if(( (self.body[ibody].frame[iframe].mainframe>=0 or self.body[ibody].frame[iframe].frameref>=0) or DisplayOK==1) and PassDisplayOK==0):
                                TexENFile.write("Rigid body S"+str(ibody)+" with "+str(self.body[ibody].nbrnodes)+" frames\n") 
                                PassDisplayOK=1
                    #poutre flexible avec ses frames
                    if(self.body[ibody].length>0 and self.body[ibody].section>0 and self.body[ibody].EYoung>0 and self.body[ibody].nu>0 and self.body[ibody].rho>0 and self.body[ibody].Iyy>0 and self.body[ibody].Izz>0 and self.body[ibody].AlphaDamp>=0 and self.body[ibody].BetaDamp>=0 and self.body[ibody].nbrnodes>=2):
                        DisplayOK=0
                        if(self.body[ibody].mainframe>=0 or self.body[ibody].frameref>=0):
                            #Si mouvement relatif, on affiche
                            DisplayOK=1
                        PassDisplayOK=0 
                        for iframe in range(0,self.body[ibody].nbrnodes):
                            #Si mouvement relatif, on affiche
                            if(( (self.body[ibody].frame[iframe].mainframe>=0 or self.body[ibody].frame[iframe].frameref>=0) or DisplayOK==1) and PassDisplayOK==0):
                                TexENFile.write("Flexible beam S"+str(ibody)+" with "+str(self.body[ibody].nbrnodes)+" frames\n") 
                                PassDisplayOK=1
                    #corps flexible avec ses frames
                    if(self.body[ibody].name!=0 and self.body[ibody].nbrnodes>=2 and self.body[ibody].AlphaDamp>=0 and self.body[ibody].BetaDamp>=0 and self.body[ibody].length==0):
                        DisplayOK=0
                        if(self.body[ibody].mainframe>=0 or self.body[ibody].frameref>=0):
                            #Si mouvement relatif, on affiche
                            DisplayOK=1 
                        PassDisplayOK=0
                        for iframe in range(0,self.body[ibody].nbrnodes):
                            #Si mouvement relatif, on affiche
                            if(( (self.body[ibody].frame[iframe].mainframe>=0 or self.body[ibody].frame[iframe].frameref>=0) or DisplayOK==1) and PassDisplayOK==0):
                                TexENFile.write("Flexible body S"+str(ibody)+" with "+str(self.body[ibody].nbrnodes)+" frames\n") 
                                PassDisplayOK=1  
    
                    if(self.body[ibody].mainframe<0):
                    # Rigid body en absolu (display nothing)
                        Jar2=0
                        for iframe in range(0,self.body[ibody].nbrnodes):
                            if(self.body[ibody].frame[iframe].mainframe>=0):
                                Jar2=Jar2+1
                        if Jar2>0:
                            TexENFile.write("\\begin{itemize}\n")
                        for iframe in range(0,self.body[ibody].nbrnodes):
                            # Frame en absolu (display nothing)
                            #if(self.body[ibody].frame[iframe].mainframe<0 and self.body[ibody].frame[iframe].frameref<0):
                                #TexENFile.write("\\[T_{0F,S"+str(ibody)+"F"+str(iframe)+"}="+str(latex(self.body[ibody].frame[iframe].TtoF))+"\\]\\\\\n\n")
                            # Flexible sur un solide rigide
                            if(self.body[ibody].frame[iframe].mainframe>=0 and self.body[ibody].frame[iframe].frameref<0):
                                TexENFile.write("  \\item whose frame F"+str(iframe)+" is relative to the main frame of body S"+str(self.body[ibody].frame[iframe].mainframe)+". \n")
                            # Flexible sur un solide flexible
                            elif(self.body[ibody].frame[iframe].mainframe>=0 and self.body[ibody].frame[iframe].frameref>=0):
                                TexENFile.write("  \\item whose frame F"+str(iframe)+" is relative to frame F"+str(self.body[ibody].frame[iframe].frameref)+" of body S"+str(self.body[ibody].frame[iframe].mainframe)+". \n")
                        if Jar2>0:
                            TexENFile.write("\\end{itemize}\n")
                    # Mainframe du body est en relatif
                    else:
                        TexENFile.write("\\begin{itemize}\n")
                        # un solide rigide par rapport à un solide rigide
                        if (self.body[ibody].frameref<0):
                            TexENFile.write("  \\item whose main frame is relative to the main frame of body S"+str(self.body[ibody].mainframe)+". \n")    
                        # un solide rigide par rapport à un solide flexible
                        else:
                            TexENFile.write("  \\item whose main frame is relative to frame F"+str(self.body[ibody].frameref)+" of body S"+str(self.body[ibody].mainframe)+". \n")
                        Jar2=0 
                        for iframe in range(0,self.body[ibody].nbrnodes):
                            if(self.body[ibody].frame[iframe].mainframe>=0):
                                Jar2=Jar2+1
                        PassEndItem=0
                        if Jar2==0:
                            TexENFile.write("\\end{itemize}\n")
                            PassEndItem=1
                        for iframe in range(0,self.body[ibody].nbrnodes):
                            # Frame en absolu
                            #if(self.body[ibody].frame[iframe].mainframe<0 and self.body[ibody].frame[iframe].frameref<0):
                                #TexENFile.write("\\[T_{0F,S"+str(ibody)+"F"+str(iframe)+"}="+str(latex(self.body[ibody].frame[iframe].TtoF))+"\\]\\\\\n\n")
                            # Flexible sur un solide rigide
                            if(self.body[ibody].frame[iframe].mainframe>=0 and self.body[ibody].frame[iframe].frameref<0):
                                TexENFile.write("  \\item whose frame F"+str(iframe)+" is relative to main frame of body S"+str(self.body[ibody].frame[iframe].mainframe)+". \n")
                                # Flexible sur un solide flexible
                            elif(self.body[ibody].frame[iframe].mainframe>=0 and self.body[ibody].frame[iframe].frameref>=0):
                                TexENFile.write("  \\item whose frame F"+str(iframe)+" is relative to frame F"+str(self.body[ibody].frame[iframe].frameref)+" of body S"+str(self.body[ibody].frame[iframe].mainframe)+". \n")
                        if(PassEndItem==0):
                            TexENFile.write("\\end{itemize}\n")       
        TexENFile.write("\n")

        if self.nbrdep > 0:
            TexENFile.write("\\subsection*{Dependent variables}\n\n")
            TexENFile.write("Dependent variables have been used to perform an efficient simulation. ")
            TexENFile.write(str(self.nbrdep)+" parameters have been defined, denoted by $p_i$ ($i$ from 0 to "+str(self.nbrdep-1)+") :\n")
            for idep in range(0,self.nbrdep):
                TexENFile.write("\\[p_"+str(idep)+"="+str(latex(self.pexpr[idep]))+"\\]\n")
                TexENFile.write("\\[\\dot{p}_"+str(idep)+"="+str(latex(self.pdexpr[idep]))+"\\]\n")
                TexENFile.write("\\[\\ddot{p}_"+str(idep)+"="+str(latex(self.pddexpr[idep]))+"\\]\n")
        # Cas dans lequel il n'y a pas de corps flexibles => cas purement rigide
        if (self.nbrFlexBeam+self.nbrFlexBody)==0: 
            TexENFile.write("\\subsection*{Homogeneous transformation matrix of each body}\n\n")
            for ibody in range(0,self.nbrbody):
                if self.body[ibody].mainframe<0:
                    TexENFile.write("\\[T_{0G,S"+str(ibody)+"}="+str(latex(self.body[ibody].TtoF))+"\\]\\\\\n\n")
                else:
                    TexENFile.write("\\[T_{refG,S"+str(ibody)+"/S"+str(self.body[ibody].mainframe)+"}="+str(latex(self.body[ibody].TtoF))+"\\]\\\\\n\n")
            TexENFile.write("\\subsection*{Their time derivative}\n\n")
            for ibody in range(0,self.nbrbody):
                if self.body[ibody].mainframe<0:
                    TexENFile.write("\\[\\dot{T}_{0G,S"+str(ibody)+"}="+str(latex(self.body[ibody].TtoFd))+"\\]\\\\\n\n")
                else:
                    TexENFile.write("\\[\\dot{T}_{refG,S"+str(ibody)+"/S"+str(self.body[ibody].mainframe)+"}="+str(latex(self.body[ibody].TtoFd))+"\\]\\\\\n\n")
            TexENFile.write("\\subsection*{The velocity of the center of gravity of each body}\n\n")
            for ibody in range(0,self.nbrbody):
                if self.body[ibody].mainframe<0:
                    TexENFile.write("\\[\\vec{v}_{G,S"+str(ibody)+"}="+str(latex(self.body[ibody].vF))+"\\]\\\\\n\n")
                else:
                    TexENFile.write("\\[\\left\\{\\vec{v}_{G,S"+str(ibody)+"/S"+str(self.body[ibody].mainframe)+"}\\right\\}_{S"+str(self.body[ibody].mainframe)+"}="+str(latex(self.body[ibody].vF))+"\\]\\\\\n\n")
            TexENFile.write("\\subsection*{The acceleration of the center of gravity of each body}\n\n")
            for ibody in range(0,self.nbrbody):
                if self.body[ibody].mainframe<0:
                    TexENFile.write("\\[\\vec{a}_{G,S"+str(ibody)+"}="+str(latex(self.body[ibody].aF))+"\\]\\\\\n\n")
                else:
                    TexENFile.write("\\[\\left\\{\\vec{a}_{G,S"+str(ibody)+"/S"+str(self.body[ibody].mainframe)+"}\\right\\}_{S"+str(self.body[ibody].mainframe)+"}="+str(latex(self.body[ibody].aF))+"\\]\\\\\n\n")
            TexENFile.write("\\subsection*{The rotation velocity of each body}\n\n")
            for ibody in range(0,self.nbrbody):
                if self.body[ibody].mainframe<0:
                    TexENFile.write("\\[\\vec{\omega}_{G,S"+str(ibody)+"}="+str(latex(self.body[ibody].omega))+"\\]\\\\\n\n")
                else:
                    TexENFile.write("\\[\\left\\{\\vec{\omega}_{G,S"+str(ibody)+"/S"+str(self.body[ibody].mainframe)+"}\\right\\}_{S"+str(self.body[ibody].mainframe)+"}="+str(latex(self.body[ibody].omega))+"\\]\\\\\n\n")
            TexENFile.write("\\subsection*{The rotation acceleration of each body}\n\n")
            for ibody in range(0,self.nbrbody):
                if self.body[ibody].mainframe<0:
                    TexENFile.write("\\[\\vec{\dot{\omega}}_{G,S"+str(ibody)+"}="+str(latex(self.body[ibody].omegad))+"\\]\\\\\n\n")
                else:
                    TexENFile.write("\\[\\left\\{\\vec{\dot{\omega}}_{G,S"+str(ibody)+"/S"+str(self.body[ibody].mainframe)+"}\\right\\}_{S"+str(self.body[ibody].mainframe)+"}="+str(latex(self.body[ibody].omegad))+"\\]\\\\\n\n")
        else:
            TexENFile.write("\\subsection*{Homogeneous transformation matrix of each body and each frame}\n\n")
            for ibody in range(0,self.nbrbody):
                #corps purement rigide avec frames possibles
                if(self.body[ibody].mass>0 and self.body[ibody].Ixx>0 and self.body[ibody].Iyy>0 and self.body[ibody].Izz>0 and self.body[ibody].Ixy>=0 and self.body[ibody].Ixz>=0 and self.body[ibody].Iyz>=0):
                    TexENFile.write("Rigid body S"+str(ibody)+" with "+str(self.body[ibody].nbrnodes)+" frame(s) \\\\\n")   
                #poutre flexible avec ses frames
                if(self.body[ibody].length>0 and self.body[ibody].section>0 and self.body[ibody].EYoung>0 and self.body[ibody].nu>0 and self.body[ibody].rho>0 and self.body[ibody].Iyy>0 and self.body[ibody].Izz>0 and self.body[ibody].AlphaDamp>=0 and self.body[ibody].BetaDamp>=0 and self.body[ibody].nbrnodes>=2):
                    TexENFile.write("Flexible beam S"+str(ibody)+" with "+str(self.body[ibody].nbrnodes)+" frames \\\\\n")   
                #corps flexible avec ses frames
                if(self.body[ibody].name!=0 and self.body[ibody].nbrnodes>=2 and self.body[ibody].AlphaDamp>=0 and self.body[ibody].BetaDamp>=0 and self.body[ibody].length==0):
                    TexENFile.write("Flexible body S"+str(ibody)+" with "+str(self.body[ibody].nbrnodes)+" frames \\\\\n")

                if(self.body[ibody].mainframe<0):
                    # Rigid body en absolu
                    TexENFile.write("\\[T_{0F,S"+str(ibody)+"}="+str(latex(self.body[ibody].TtoF))+"\\]\\\\\n\n")
                    for iframe in range(0,self.body[ibody].nbrnodes):
                        # Frame en absolu
                        if(self.body[ibody].frame[iframe].mainframe<0 and self.body[ibody].frame[iframe].frameref<0):
                            TexENFile.write("\\[T_{0F,S"+str(ibody)+"F"+str(iframe)+"}="+str(latex(self.body[ibody].frame[iframe].TtoF))+"\\]\\\\\n\n")
                        # Flexible sur un solide rigide
                        elif(self.body[ibody].frame[iframe].mainframe>=0 and self.body[ibody].frame[iframe].frameref<0):
                            TexENFile.write("\\[T_{refF,S"+str(ibody)+"F"+str(iframe)+"/S"+str(self.body[ibody].frame[iframe].mainframe)+"}="+str(latex(self.body[ibody].frame[iframe].TtoF))+"\\]\\\\\n\n")
                        # Flexible sur un solide flexible
                        elif(self.body[ibody].frame[iframe].mainframe>=0 and self.body[ibody].frame[iframe].frameref>=0):
                            TexENFile.write("\\[T_{refF,S"+str(ibody)+"F"+str(iframe)+"/S"+str(self.body[ibody].frame[iframe].mainframe)+"F"+str(self.body[ibody].frame[iframe].frameref)+"}="+str(latex(self.body[ibody].frame[iframe].TtoF))+"\\]\\\\\n\n")
                # Mainframe du body est en relatif
                else:
                    # un solide rigide par rapport à un solide rigide
                    if (self.body[ibody].frameref<0):
                        TexENFile.write("\\[T_{refF,S"+str(ibody)+"/S"+str(self.body[ibody].mainframe)+"}="+str(latex(self.body[ibody].TtoF))+"\\]\\\\\n\n")    
                    # un solide rigide par rapport à un solide flexible
                    else:
                        TexENFile.write("\\[T_{refF,S"+str(ibody)+"/S"+str(self.body[ibody].mainframe)+"F"+str(self.body[ibody].frameref)+"}="+str(latex(self.body[ibody].TtoF))+"\\]\\\\\n\n")
                    for iframe in range(0,self.body[ibody].nbrnodes):
                        # Frame en absolu
                        if(self.body[ibody].frame[iframe].mainframe<0 and self.body[ibody].frame[iframe].frameref<0):
                            TexENFile.write("\\[T_{0F,S"+str(ibody)+"F"+str(iframe)+"}="+str(latex(self.body[ibody].frame[iframe].TtoF))+"\\]\\\\\n\n")
                        # Flexible sur un solide rigide
                        elif(self.body[ibody].frame[iframe].mainframe>=0 and self.body[ibody].frame[iframe].frameref<0):
                            TexENFile.write("\\[T_{refF,S"+str(ibody)+"F"+str(iframe)+"/S"+str(self.body[ibody].frame[iframe].mainframe)+"}="+str(latex(self.body[ibody].frame[iframe].TtoF))+"\\]\\\\\n\n")
                        # Flexible sur un solide flexible
                        elif(self.body[ibody].frame[iframe].mainframe>=0 and self.body[ibody].frame[iframe].frameref>=0):
                            TexENFile.write("\\[T_{refF,S"+str(ibody)+"F"+str(iframe)+"/S"+str(self.body[ibody].frame[iframe].mainframe)+"F"+str(self.body[ibody].frame[iframe].frameref)+"}="+str(latex(self.body[ibody].frame[iframe].TtoF))+"\\]\\\\\n\n")
                      
            TexENFile.write("\\subsection*{Their time derviative}\n\n")
            for ibody in range(0,self.nbrbody):
                #corps purement rigide avec frames possibles
                if(self.body[ibody].mass>0 and self.body[ibody].Ixx>0 and self.body[ibody].Iyy>0 and self.body[ibody].Izz>0 and self.body[ibody].Ixy>=0 and self.body[ibody].Ixz>=0 and self.body[ibody].Iyz>=0):
                    TexENFile.write("Rigid body S"+str(ibody)+" with "+str(self.body[ibody].nbrnodes)+" frames \\\\\n")   
                #poutre flexible avec ses frames
                if(self.body[ibody].length>0 and self.body[ibody].section>0 and self.body[ibody].EYoung>0 and self.body[ibody].nu>0 and self.body[ibody].rho>0 and self.body[ibody].Iyy>0 and self.body[ibody].Izz>0 and self.body[ibody].AlphaDamp>=0 and self.body[ibody].BetaDamp>=0 and self.body[ibody].nbrnodes>=2):
                    TexENFile.write("Flexible beam S"+str(ibody)+" with "+str(self.body[ibody].nbrnodes)+" frames \\\\\n")   
                #corps flexible avec ses frames
                if(self.body[ibody].name!=0 and self.body[ibody].nbrnodes>=2 and self.body[ibody].AlphaDamp>=0 and self.body[ibody].BetaDamp>=0 and self.body[ibody].length==0):
                    TexENFile.write("Flexible body S"+str(ibody)+" with "+str(self.body[ibody].nbrnodes)+" frames \\\\\n")

                if(self.body[ibody].mainframe<0):
                    # Rigid body en absolu
                    TexENFile.write("\\[\\dot{T}_{0F,S"+str(ibody)+"}="+str(latex(self.body[ibody].TtoFd))+"\\]\\\\\n\n")
                    for iframe in range(0,self.body[ibody].nbrnodes):
                        # Frame en absolu
                        if(self.body[ibody].frame[iframe].mainframe<0 and self.body[ibody].frame[iframe].frameref<0):
                            TexENFile.write("\\[\\dot{T}_{0F,S"+str(ibody)+"F"+str(iframe)+"}="+str(latex(self.body[ibody].frame[iframe].TtoFd))+"\\]\\\\\n\n")
                        # Flexible sur un solide rigide
                        elif(self.body[ibody].frame[iframe].mainframe>=0 and self.body[ibody].frame[iframe].frameref<0):
                            TexENFile.write("\\[\\dot{T}_{refF,S"+str(ibody)+"F"+str(iframe)+"/S"+str(self.body[ibody].frame[iframe].mainframe)+"}="+str(latex(self.body[ibody].frame[iframe].TtoFd))+"\\]\\\\\n\n")
                        # Flexible sur un solide flexible
                        elif(self.body[ibody].frame[iframe].mainframe>=0 and self.body[ibody].frame[iframe].frameref>=0):
                            TexENFile.write("\\[\\dot{T}_{refF,S"+str(ibody)+"F"+str(iframe)+"/S"+str(self.body[ibody].frame[iframe].mainframe)+"F"+str(self.body[ibody].frame[iframe].frameref)+"}="+str(latex(self.body[ibody].frame[iframe].TtoFd))+"\\]\\\\\n\n")
                # Mainframe du body est en relatif
                else:
                    # un solide rigide par rapport à un solide rigide
                    if (self.body[ibody].frameref<0):
                        TexENFile.write("\\[\\dot{T}_{refF,S"+str(ibody)+"/S"+str(self.body[ibody].mainframe)+"}="+str(latex(self.body[ibody].TtoFd))+"\\]\\\\\n\n")    
                    # un solide rigide par rapport à un solide flexible
                    else:
                        TexENFile.write("\\[\\dot{T}_{refF,S"+str(ibody)+"/S"+str(self.body[ibody].mainframe)+"F"+str(self.body[ibody].frameref)+"}="+str(latex(self.body[ibody].TtoFd))+"\\]\\\\\n\n")
                    for iframe in range(0,self.body[ibody].nbrnodes):
                        # Frame en absolu
                        if(self.body[ibody].frame[iframe].mainframe<0 and self.body[ibody].frame[iframe].frameref<0):
                            TexENFile.write("\\[\\dot{T}_{0F,S"+str(ibody)+"F"+str(iframe)+"}="+str(latex(self.body[ibody].frame[iframe].TtoFd))+"\\]\\\\\n\n")
                        # Flexible sur un solide rigide
                        elif(self.body[ibody].frame[iframe].mainframe>=0 and self.body[ibody].frame[iframe].frameref<0):
                            TexENFile.write("\\[\\dot{T}_{refF,S"+str(ibody)+"F"+str(iframe)+"/S"+str(self.body[ibody].frame[iframe].mainframe)+"}="+str(latex(self.body[ibody].frame[iframe].TtoFd))+"\\]\\\\\n\n")
                        # Flexible sur un solide flexible
                        elif(self.body[ibody].frame[iframe].mainframe>=0 and self.body[ibody].frame[iframe].frameref>=0):
                            TexENFile.write("\\[\\dot{T}_{refF,S"+str(ibody)+"F"+str(iframe)+"/S"+str(self.body[ibody].frame[iframe].mainframe)+"F"+str(self.body[ibody].frame[iframe].frameref)+"}="+str(latex(self.body[ibody].frame[iframe].TtoFd))+"\\]\\\\\n\n")

            TexENFile.write("\\subsection*{The velocity of the center of gravity of each body and each frame}\n\n")
            for ibody in range(0,self.nbrbody):
                #corps purement rigide avec frames possibles
                if(self.body[ibody].mass>0 and self.body[ibody].Ixx>0 and self.body[ibody].Iyy>0 and self.body[ibody].Izz>0 and self.body[ibody].Ixy>=0 and self.body[ibody].Ixz>=0 and self.body[ibody].Iyz>=0):
                    TexENFile.write("Rigid body S"+str(ibody)+" with "+str(self.body[ibody].nbrnodes)+" frames \\\\\n")   
                #poutre flexible avec ses frames
                if(self.body[ibody].length>0 and self.body[ibody].section>0 and self.body[ibody].EYoung>0 and self.body[ibody].nu>0 and self.body[ibody].rho>0 and self.body[ibody].Iyy>0 and self.body[ibody].Izz>0 and self.body[ibody].AlphaDamp>=0 and self.body[ibody].BetaDamp>=0 and self.body[ibody].nbrnodes>=2):
                    TexENFile.write("Flexible beam S"+str(ibody)+" with "+str(self.body[ibody].nbrnodes)+" frames \\\\\n")   
                #corps flexible avec ses frames
                if(self.body[ibody].name!=0 and self.body[ibody].nbrnodes>=2 and self.body[ibody].AlphaDamp>=0 and self.body[ibody].BetaDamp>=0 and self.body[ibody].length==0):
                    TexENFile.write("Flexible body S"+str(ibody)+" with "+str(self.body[ibody].nbrnodes)+" frames \\\\\n")

                if(self.body[ibody].mainframe<0):
                    # Rigid body en absolu
                    TexENFile.write("\\[\\vec{v}_{0F,S"+str(ibody)+"}="+str(latex(self.body[ibody].vF))+"\\]\\\\\n\n")
                    for iframe in range(0,self.body[ibody].nbrnodes):
                        # Frame en absolu
                        if(self.body[ibody].frame[iframe].mainframe<0 and self.body[ibody].frame[iframe].frameref<0):
                            TexENFile.write("\\[\\vec{v}_{0F,S"+str(ibody)+"F"+str(iframe)+"}="+str(latex(self.body[ibody].frame[iframe].vF))+"\\]\\\\\n\n")
                        # Flexible sur un solide rigide
                        elif(self.body[ibody].frame[iframe].mainframe>=0 and self.body[ibody].frame[iframe].frameref<0):
                            TexENFile.write("\\[\\left\\{\\vec{v}_{refF,S"+str(ibody)+"F"+str(iframe)+"/S"+str(self.body[ibody].frame[iframe].mainframe)+"}\\right\\}_{S"+str(self.body[ibody].frame[iframe].mainframe)+"}="+str(latex(self.body[ibody].frame[iframe].vF))+"\\]\\\\\n\n")
                        # Flexible sur un solide flexible
                        elif(self.body[ibody].frame[iframe].mainframe>=0 and self.body[ibody].frame[iframe].frameref>=0):
                            TexENFile.write("\\[\\left\\{\\vec{v}_{refF,S"+str(ibody)+"F"+str(iframe)+"/S"+str(self.body[ibody].frame[iframe].mainframe)+"F"+str(self.body[ibody].frame[iframe].frameref)+"}\\right\\}_{S"+str(self.body[ibody].frame[iframe].mainframe)+"F"+str(self.body[ibody].frame[iframe].frameref)+"}="+str(latex(self.body[ibody].frame[iframe].vF))+"\\]\\\\\n\n")
                # Mainframe du body est en relatif
                else:
                    # un solide rigide par rapport à un solide rigide
                    if (self.body[ibody].frameref<0):
                        TexENFile.write("\\[\\left\\{\\vec{v}_{refF,S"+str(ibody)+"/S"+str(self.body[ibody].mainframe)+"}\\right\\}_{S"+str(self.body[ibody].mainframe)+"}="+str(latex(self.body[ibody].vF))+"\\]\\\\\n\n")    
                    # un solide rigide par rapport à un solide flexible
                    else:
                        TexENFile.write("\\[\\left\\{\\vec{v}_{refF,S"+str(ibody)+"/S"+str(self.body[ibody].mainframe)+"F"+str(self.body[ibody].frameref)+"}\\right\\}_{S"+str(self.body[ibody].mainframe)+"F"+str(self.body[ibody].frameref)+"}="+str(latex(self.body[ibody].vF))+"\\]\\\\\n\n")
                    for iframe in range(0,self.body[ibody].nbrnodes):
                        # Frame en absolu
                        if(self.body[ibody].frame[iframe].mainframe<0 and self.body[ibody].frame[iframe].frameref<0):
                            TexENFile.write("\\[\\vec{v}_{0F,S"+str(ibody)+"F"+str(iframe)+"}="+str(latex(self.body[ibody].frame[iframe].vF))+"\\]\\\\\n\n")
                        # Flexible sur un solide rigide
                        elif(self.body[ibody].frame[iframe].mainframe>=0 and self.body[ibody].frame[iframe].frameref<0):
                            TexENFile.write("\\[\\left\\{\\vec{v}_{refF,S"+str(ibody)+"F"+str(iframe)+"/S"+str(self.body[ibody].frame[iframe].mainframe)+"}\\right\\}_{S"+str(self.body[ibody].frame[iframe].mainframe)+"}="+str(latex(self.body[ibody].frame[iframe].vF))+"\\]\\\\\n\n")
                        # Flexible sur un solide flexible
                        elif(self.body[ibody].frame[iframe].mainframe>=0 and self.body[ibody].frame[iframe].frameref>=0):
                            TexENFile.write("\\[\\left\\{\\vec{v}_{refF,S"+str(ibody)+"F"+str(iframe)+"/S"+str(self.body[ibody].frame[iframe].mainframe)+"F"+str(self.body[ibody].frame[iframe].frameref)+"}\\right\\}_{S"+str(self.body[ibody].frame[iframe].mainframe)+"F"+str(self.body[ibody].frame[iframe].frameref)+"}="+str(latex(self.body[ibody].frame[iframe].vF))+"\\]\\\\\n\n")

            TexENFile.write("\\subsection*{The acceleration of the center of gravity of each body and each frame}\n\n")
            for ibody in range(0,self.nbrbody):
                #corps purement rigide avec frames possibles
                if(self.body[ibody].mass>0 and self.body[ibody].Ixx>0 and self.body[ibody].Iyy>0 and self.body[ibody].Izz>0 and self.body[ibody].Ixy>=0 and self.body[ibody].Ixz>=0 and self.body[ibody].Iyz>=0):
                    TexENFile.write("Rigid body S"+str(ibody)+" with "+str(self.body[ibody].nbrnodes)+" frames \\\\\n")   
                #poutre flexible avec ses frames
                if(self.body[ibody].length>0 and self.body[ibody].section>0 and self.body[ibody].EYoung>0 and self.body[ibody].nu>0 and self.body[ibody].rho>0 and self.body[ibody].Iyy>0 and self.body[ibody].Izz>0 and self.body[ibody].AlphaDamp>=0 and self.body[ibody].BetaDamp>=0 and self.body[ibody].nbrnodes>=2):
                    TexENFile.write("Flexible beam S"+str(ibody)+" with "+str(self.body[ibody].nbrnodes)+" frames \\\\\n")   
                #corps flexible avec ses frames
                if(self.body[ibody].name!=0 and self.body[ibody].nbrnodes>=2 and self.body[ibody].AlphaDamp>=0 and self.body[ibody].BetaDamp>=0 and self.body[ibody].length==0):
                    TexENFile.write("Flexible body S"+str(ibody)+" with "+str(self.body[ibody].nbrnodes)+" frames \\\\\n")

                if(self.body[ibody].mainframe<0):
                    # Rigid body en absolu
                    TexENFile.write("\\[\\vec{a}_{0F,S"+str(ibody)+"}="+str(latex(self.body[ibody].aF))+"\\]\\\\\n\n")
                    for iframe in range(0,self.body[ibody].nbrnodes):
                        # Frame en absolu
                        if(self.body[ibody].frame[iframe].mainframe<0 and self.body[ibody].frame[iframe].frameref<0):
                            TexENFile.write("\\[\\vec{a}_{0F,S"+str(ibody)+"F"+str(iframe)+"}="+str(latex(self.body[ibody].frame[iframe].aF))+"\\]\\\\\n\n")
                        # Flexible sur un solide rigide
                        elif(self.body[ibody].frame[iframe].mainframe>=0 and self.body[ibody].frame[iframe].frameref<0):
                            TexENFile.write("\\[\\left\\{\\vec{a}_{refF,S"+str(ibody)+"F"+str(iframe)+"/S"+str(self.body[ibody].frame[iframe].mainframe)+"}\\right\\}_{S"+str(self.body[ibody].frame[iframe].mainframe)+"}="+str(latex(self.body[ibody].frame[iframe].aF))+"\\]\\\\\n\n")
                        # Flexible sur un solide flexible
                        elif(self.body[ibody].frame[iframe].mainframe>=0 and self.body[ibody].frame[iframe].frameref>=0):
                            TexENFile.write("\\[\\left\\{\\vec{a}_{refF,S"+str(ibody)+"F"+str(iframe)+"/S"+str(self.body[ibody].frame[iframe].mainframe)+"F"+str(self.body[ibody].frame[iframe].frameref)+"}\\right\\}_{S"+str(self.body[ibody].frame[iframe].mainframe)+"F"+str(self.body[ibody].frame[iframe].frameref)+"}="+str(latex(self.body[ibody].frame[iframe].aF))+"\\]\\\\\n\n")
                # Mainframe du body est en relatif
                else:
                    # un solide rigide par rapport à un solide rigide
                    if (self.body[ibody].frameref<0):
                        TexENFile.write("\\[\\left\\{\\vec{a}_{refF,S"+str(ibody)+"/S"+str(self.body[ibody].mainframe)+"}\\right\\}_{S"+str(self.body[ibody].mainframe)+"}="+str(latex(self.body[ibody].aF))+"\\]\\\\\n\n")    
                    # un solide rigide par rapport à un solide flexible
                    else:
                        TexENFile.write("\\[\\left\\{\\vec{a}_{refF,S"+str(ibody)+"/S"+str(self.body[ibody].mainframe)+"F"+str(self.body[ibody].frameref)+"}\\right\\}_{S"+str(self.body[ibody].mainframe)+"F"+str(self.body[ibody].frameref)+"}="+str(latex(self.body[ibody].aF))+"\\]\\\\\n\n")
                    for iframe in range(0,self.body[ibody].nbrnodes):
                        # Frame en absolu
                        if(self.body[ibody].frame[iframe].mainframe<0 and self.body[ibody].frame[iframe].frameref<0):
                            TexENFile.write("\\[\\vec{a}_{0F,S"+str(ibody)+"F"+str(iframe)+"}="+str(latex(self.body[ibody].frame[iframe].aF))+"\\]\\\\\n\n")
                        # Flexible sur un solide rigide
                        elif(self.body[ibody].frame[iframe].mainframe>=0 and self.body[ibody].frame[iframe].frameref<0):
                            TexENFile.write("\\[\\left\\{\\vec{a}_{refF,S"+str(ibody)+"F"+str(iframe)+"/S"+str(self.body[ibody].frame[iframe].mainframe)+"}\\right\\}_{S"+str(self.body[ibody].frame[iframe].mainframe)+"}="+str(latex(self.body[ibody].frame[iframe].aF))+"\\]\\\\\n\n")
                        # Flexible sur un solide flexible
                        elif(self.body[ibody].frame[iframe].mainframe>=0 and self.body[ibody].frame[iframe].frameref>=0):
                            TexENFile.write("\\[\\left\\{\\vec{a}_{refF,S"+str(ibody)+"F"+str(iframe)+"/S"+str(self.body[ibody].frame[iframe].mainframe)+"F"+str(self.body[ibody].frame[iframe].frameref)+"}\\right\\}_{S"+str(self.body[ibody].frame[iframe].mainframe)+"F"+str(self.body[ibody].frame[iframe].frameref)+"}="+str(latex(self.body[ibody].frame[iframe].aF))+"\\]\\\\\n\n")
            TexENFile.write("\\subsection*{The rotation velocity of each body and each frame}\n\n")
            for ibody in range(0,self.nbrbody):
                #corps purement rigide avec frames possibles
                if(self.body[ibody].mass>0 and self.body[ibody].Ixx>0 and self.body[ibody].Iyy>0 and self.body[ibody].Izz>0 and self.body[ibody].Ixy>=0 and self.body[ibody].Ixz>=0 and self.body[ibody].Iyz>=0):
                    TexENFile.write("Rigid body S"+str(ibody)+" with "+str(self.body[ibody].nbrnodes)+" frames \\\\\n")   
                #poutre flexible avec ses frames
                if(self.body[ibody].length>0 and self.body[ibody].section>0 and self.body[ibody].EYoung>0 and self.body[ibody].nu>0 and self.body[ibody].rho>0 and self.body[ibody].Iyy>0 and self.body[ibody].Izz>0 and self.body[ibody].AlphaDamp>=0 and self.body[ibody].BetaDamp>=0 and self.body[ibody].nbrnodes>=2):
                    TexENFile.write("Flexible beam S"+str(ibody)+" with "+str(self.body[ibody].nbrnodes)+" frames \\\\\n")   
                #corps flexible avec ses frames
                if(self.body[ibody].name!=0 and self.body[ibody].nbrnodes>=2 and self.body[ibody].AlphaDamp>=0 and self.body[ibody].BetaDamp>=0 and self.body[ibody].length==0):
                    TexENFile.write("Flexible body S"+str(ibody)+" with "+str(self.body[ibody].nbrnodes)+" frames \\\\\n")

                if(self.body[ibody].mainframe<0):
                    # Rigid body en absolu
                    TexENFile.write("\\[\\vec{\omega}_{0F,S"+str(ibody)+"}="+str(latex(self.body[ibody].omega))+"\\]\\\\\n\n")
                    for iframe in range(0,self.body[ibody].nbrnodes):
                        # Frame en absolu
                        if(self.body[ibody].frame[iframe].mainframe<0 and self.body[ibody].frame[iframe].frameref<0):
                            TexENFile.write("\\[\\vec{\omega}_{0F,S"+str(ibody)+"F"+str(iframe)+"}="+str(latex(self.body[ibody].frame[iframe].omega))+"\\]\\\\\n\n")
                        # Flexible sur un solide rigide
                        elif(self.body[ibody].frame[iframe].mainframe>=0 and self.body[ibody].frame[iframe].frameref<0):
                            TexENFile.write("\\[\\left\\{\\vec{\omega}_{refF,S"+str(ibody)+"F"+str(iframe)+"/S"+str(self.body[ibody].frame[iframe].mainframe)+"}\\right\\}_{S"+str(self.body[ibody].frame[iframe].mainframe)+"}="+str(latex(self.body[ibody].frame[iframe].omega))+"\\]\\\\\n\n")
                        # Flexible sur un solide flexible
                        elif(self.body[ibody].frame[iframe].mainframe>=0 and self.body[ibody].frame[iframe].frameref>=0):
                            TexENFile.write("\\[\\left\\{\\vec{\omega}_{refF,S"+str(ibody)+"F"+str(iframe)+"/S"+str(self.body[ibody].frame[iframe].mainframe)+"F"+str(self.body[ibody].frame[iframe].frameref)+"}\\right\\}_{S"+str(self.body[ibody].frame[iframe].mainframe)+"F"+str(self.body[ibody].frame[iframe].frameref)+"}="+str(latex(self.body[ibody].frame[iframe].omega))+"\\]\\\\\n\n")
                # Mainframe du body est en relatif
                else:
                    # un solide rigide par rapport à un solide rigide
                    if (self.body[ibody].frameref<0):
                        TexENFile.write("\\[\\left\\{\\vec{\omega}_{refF,S"+str(ibody)+"/S"+str(self.body[ibody].mainframe)+"}\\right\\}_{S"+str(self.body[ibody].mainframe)+"}="+str(latex(self.body[ibody].omega))+"\\]\\\\\n\n")    
                    # un solide rigide par rapport à un solide flexible
                    else:
                        TexENFile.write("\\[\\left\\{\\vec{\omega}_{refF,S"+str(ibody)+"/S"+str(self.body[ibody].mainframe)+"F"+str(self.body[ibody].frameref)+"}\\right\\}_{S"+str(self.body[ibody].mainframe)+"F"+str(self.body[ibody].frameref)+"}="+str(latex(self.body[ibody].omega))+"\\]\\\\\n\n")
                    for iframe in range(0,self.body[ibody].nbrnodes):
                        # Frame en absolu
                        if(self.body[ibody].frame[iframe].mainframe<0 and self.body[ibody].frame[iframe].frameref<0):
                            TexENFile.write("\\[\\vec{\omega}_{0F,S"+str(ibody)+"F"+str(iframe)+"}="+str(latex(self.body[ibody].frame[iframe].omega))+"\\]\\\\\n\n")
                        # Flexible sur un solide rigide
                        elif(self.body[ibody].frame[iframe].mainframe>=0 and self.body[ibody].frame[iframe].frameref<0):
                            TexENFile.write("\\[\\left\\{\\vec{\omega}_{refF,S"+str(ibody)+"F"+str(iframe)+"/S"+str(self.body[ibody].frame[iframe].mainframe)+"}\\right\\}_{S"+str(self.body[ibody].frame[iframe].mainframe)+"}="+str(latex(self.body[ibody].frame[iframe].omega))+"\\]\\\\\n\n")
                        # Flexible sur un solide flexible
                        elif(self.body[ibody].frame[iframe].mainframe>=0 and self.body[ibody].frame[iframe].frameref>=0):
                            TexENFile.write("\\[\\left\\{\\vec{\omega}_{refF,S"+str(ibody)+"F"+str(iframe)+"/S"+str(self.body[ibody].frame[iframe].mainframe)+"F"+str(self.body[ibody].frame[iframe].frameref)+"}\\right\\}_{S"+str(self.body[ibody].frame[iframe].mainframe)+"F"+str(self.body[ibody].frame[iframe].frameref)+"}="+str(latex(self.body[ibody].frame[iframe].omega))+"\\]\\\\\n\n")
            TexENFile.write("\\subsection*{The rotation acceleration of each body and each frame}\n\n")
            for ibody in range(0,self.nbrbody):
                #corps purement rigide avec frames possibles
                if(self.body[ibody].mass>0 and self.body[ibody].Ixx>0 and self.body[ibody].Iyy>0 and self.body[ibody].Izz>0 and self.body[ibody].Ixy>=0 and self.body[ibody].Ixz>=0 and self.body[ibody].Iyz>=0):
                    TexENFile.write("Rigid body S"+str(ibody)+" with "+str(self.body[ibody].nbrnodes)+" frames \\\\\n")   
                #poutre flexible avec ses frames
                if(self.body[ibody].length>0 and self.body[ibody].section>0 and self.body[ibody].EYoung>0 and self.body[ibody].nu>0 and self.body[ibody].rho>0 and self.body[ibody].Iyy>0 and self.body[ibody].Izz>0 and self.body[ibody].AlphaDamp>=0 and self.body[ibody].BetaDamp>=0 and self.body[ibody].nbrnodes>=2):
                    TexENFile.write("Flexible beam S"+str(ibody)+" with "+str(self.body[ibody].nbrnodes)+" frames \\\\\n")   
                #corps flexible avec ses frames
                if(self.body[ibody].name!=0 and self.body[ibody].nbrnodes>=2 and self.body[ibody].AlphaDamp>=0 and self.body[ibody].BetaDamp>=0 and self.body[ibody].length==0):
                    TexENFile.write("Flexible body S"+str(ibody)+" with "+str(self.body[ibody].nbrnodes)+" frames \\\\\n")

                if(self.body[ibody].mainframe<0):
                    # Rigid body en absolu
                    TexENFile.write("\\[\\vec{\dot{\omega}}_{0F,S"+str(ibody)+"}="+str(latex(self.body[ibody].omegad))+"\\]\\\\\n\n")
                    for iframe in range(0,self.body[ibody].nbrnodes):
                        # Frame en absolu
                        if(self.body[ibody].frame[iframe].mainframe<0 and self.body[ibody].frame[iframe].frameref<0):
                            TexENFile.write("\\[\\vec{\dot{\omega}}_{0F,S"+str(ibody)+"F"+str(iframe)+"}="+str(latex(self.body[ibody].frame[iframe].omegad))+"\\]\\\\\n\n")
                        # Flexible sur un solide rigide
                        elif(self.body[ibody].frame[iframe].mainframe>=0 and self.body[ibody].frame[iframe].frameref<0):
                            TexENFile.write("\\[\\left\\{\\vec{\dot{\omega}}_{refF,S"+str(ibody)+"F"+str(iframe)+"/S"+str(self.body[ibody].frame[iframe].mainframe)+"}\\right\\}_{S"+str(self.body[ibody].frame[iframe].mainframe)+"}="+str(latex(self.body[ibody].frame[iframe].omegad))+"\\]\\\\\n\n")
                        # Flexible sur un solide flexible
                        elif(self.body[ibody].frame[iframe].mainframe>=0 and self.body[ibody].frame[iframe].frameref>=0):
                            TexENFile.write("\\[\\left\\{\\vec{\dot{\omega}}_{refF,S"+str(ibody)+"F"+str(iframe)+"/S"+str(self.body[ibody].frame[iframe].mainframe)+"F"+str(self.body[ibody].frame[iframe].frameref)+"}\\right\\}_{S"+str(self.body[ibody].frame[iframe].mainframe)+"F"+str(self.body[ibody].frame[iframe].frameref)+"}="+str(latex(self.body[ibody].frame[iframe].omegad))+"\\]\\\\\n\n")
                # Mainframe du body est en relatif
                else:
                    # un solide rigide par rapport à un solide rigide
                    if (self.body[ibody].frameref<0):
                        TexENFile.write("\\[\\left\\{\\vec{\dot{\omega}}_{refF,S"+str(ibody)+"/S"+str(self.body[ibody].mainframe)+"}\\right\\}_{S"+str(self.body[ibody].mainframe)+"}="+str(latex(self.body[ibody].omegad))+"\\]\\\\\n\n")    
                    # un solide rigide par rapport à un solide flexible
                    else:
                        TexENFile.write("\\[\\left\\{\\vec{\dot{\omega}}_{refF,S"+str(ibody)+"/S"+str(self.body[ibody].mainframe)+"F"+str(self.body[ibody].frameref)+"}\\right\\}_{S"+str(self.body[ibody].mainframe)+"F"+str(self.body[ibody].frameref)+"}="+str(latex(self.body[ibody].omegad))+"\\]\\\\\n\n")
                    for iframe in range(0,self.body[ibody].nbrnodes):
                        # Frame en absolu
                        if(self.body[ibody].frame[iframe].mainframe<0 and self.body[ibody].frame[iframe].frameref<0):
                            TexENFile.write("\\[\\vec{\dot{\omega}}_{0F,S"+str(ibody)+"F"+str(iframe)+"}="+str(latex(self.body[ibody].frame[iframe].omegad))+"\\]\\\\\n\n")
                        # Flexible sur un solide rigide
                        elif(self.body[ibody].frame[iframe].mainframe>=0 and self.body[ibody].frame[iframe].frameref<0):
                            TexENFile.write("\\[\\left\\{\\vec{\dot{\omega}}_{refF,S"+str(ibody)+"F"+str(iframe)+"/S"+str(self.body[ibody].frame[iframe].mainframe)+"}\\right\\}_{S"+str(self.body[ibody].frame[iframe].mainframe)+"}="+str(latex(self.body[ibody].frame[iframe].omegad))+"\\]\\\\\n\n")
                        # Flexible sur un solide flexible
                        elif(self.body[ibody].frame[iframe].mainframe>=0 and self.body[ibody].frame[iframe].frameref>=0):
                            TexENFile.write("\\[\\left\\{\\vec{\dot{\omega}}_{refF,S"+str(ibody)+"F"+str(iframe)+"/S"+str(self.body[ibody].frame[iframe].mainframe)+"F"+str(self.body[ibody].frame[iframe].frameref)+"}\\right\\}_{S"+str(self.body[ibody].frame[iframe].mainframe)+"F"+str(self.body[ibody].frame[iframe].frameref)+"}="+str(latex(self.body[ibody].frame[iframe].omegad))+"\\]\\\\\n\n") 

        #if FORCES==1:
        if(self.ForcesInAppEff != " "):
            TexENFile.write("\n\\section{Definition of external forces}\n\n")
            TexENFile.write("External forces have been defined in the ")
            if (self.EASYDYNPlusPlusFlag==0):
                TexENFile.write("procedure \\texttt{AddAppliedEfforts()}. ")
            else:
                TexENFile.write("procedure \\texttt{AddAppliedForces()}. ")
        
        TexENFile.write("\n\\section{Simulation}\n\n")
        TexENFile.write("The routine \\texttt{NewmarkIntegration} performs the integration of the equations of motion up to ")
        TexENFile.write("time \\emph{FinalTime} by regular time intervals equal to \\emph{StepSave} and with the maximum allowed time step \\emph{StepMax} defined in the file \\texttt{"+str(self.ApplicationFileName)+"}.cpp. ")
        TexENFile.write("The following values are used:\n")
        TexENFile.write("\\begin{itemize}\n")
        TexENFile.write("   \\item \\emph{FinalTime} is the simulation duration (="+str(self.tfinal)+" $s$),\n")
        TexENFile.write("   \\item \\emph{StepSave} is the time step in the numerical integration (="+str(self.hsave)+" $s$),\n")
        TexENFile.write("   \\item \\emph{StepMax} is the maximum allowed time step (="+str(self.hmax)+" $s$),\n")
        TexENFile.write("\\end{itemize}\n")
        TexENFile.write("\\noindent the initial conditions ")
        InitialCondFlag=0 # y a t il des conditions initiales
        for idof in range(0,self.nbrdof):
            if ((abs(self.qini[idof])>0) or (abs(self.qdini[idof])>0)):
                InitialCondFlag=InitialCondFlag+1
        if InitialCondFlag>0:
                TexENFile.write("being ")
        for idof in range(0,self.nbrdof):
            if self.qini[idof] !=0:
                TexENFile.write("$q_{"+str(idof)+"} = "+str(self.qini[idof])+"$, ")
            if self.qdini[idof] !=0:
                TexENFile.write("$\dot{q}_{"+str(idof)+"} = "+str(self.qdini[idof])+"$, ")
        if InitialCondFlag==0:
            TexENFile.write(" \\'being all zero.\\\\\n")
        else:    
            TexENFile.write("the others being equal to zero.\\\\\n")
        TexENFile.write("\n\\section{Results}\n\n")
        TexENFile.write("The time evolution of the different configuration parameters and their first and second derivatives ")
        TexENFile.write("can easily be plotted by \\textsf{Gnuplot} as seen in ")
        TexENFile.write("figures \\ref{figure1} to \\ref{figure3} with the code listed below:\\\\\n\n")
  
        TexENFile.write("\\footnotesize\n")
        TexENFile.write("\\verbatiminput{"+str(self.ApplicationFileName)+".plt}\n")
        TexENFile.write("\\normalsize\n\n")
        
        TexENFile.write("\\begin{figure}[h!tb]\n")
        TexENFile.write("   \\begin{center}\n")
        TexENFile.write("      \\includegraphics[width=0.75\\textwidth]{figure1.eps}\n")
        TexENFile.write("   \\end{center}\n")
        TexENFile.write("\\caption{Time evolution of configuration parameters}\n")
        TexENFile.write("\\label{figure1}\n")
        TexENFile.write("\\end{figure}\n\n")

        TexENFile.write("\\begin{figure}[h!tb]\n")
        TexENFile.write("   \\begin{center}\n")
        TexENFile.write("      \\includegraphics[width=0.75\\textwidth]{figure2.eps}\n")
        TexENFile.write("   \\end{center}\n")
        TexENFile.write("\\caption{Time evolution of time derivatives of configuration parameters}\n")
        TexENFile.write("\\label{figure2}\n")
        TexENFile.write("\\end{figure}\n\n")

        TexENFile.write("\\begin{figure}[h!tb]\n")
        TexENFile.write("   \\begin{center}\n")
        TexENFile.write("      \\includegraphics[width=0.75\\textwidth]{figure3.eps}\n")
        TexENFile.write("   \\end{center}\n")
        TexENFile.write("\\caption{Time evolution of second time derivatives of configuration parameters}\n")
        TexENFile.write("\\label{figure3}\n")
        TexENFile.write("\\end{figure}\n\n")
                
        TexENFile.write("\n\\appendix\n\n")
        TexENFile.write("\\clearpage\n")
        TexENFile.write("\\section{User's Python code}\n")
        TexENFile.write("\\footnotesize\n")
        TexENFile.write("\\verbatiminput{"+str(self.ApplicationFileName)+".py}\n\n")
        TexENFile.write("\\end{document}\n") 
        
        TexENFile.close()

        # Post-process of EN Latex Code to replace qd by \dot{}  
        filedataEN = None
        with open('Report_'+str(self.ApplicationFileName)+'.tex', 'r') as file :
            filedataEN = file.read()
        # Replace the target string
        for idof in range(0,self.nbrdof):
            filedataEN = filedataEN.replace('q['+str(idof)+']', 'q_'+str(idof))
            filedataEN = filedataEN.replace('qd['+str(idof)+']', '\dot{q}_'+str(idof))
            filedataEN = filedataEN.replace('qdd['+str(idof)+']', '\ddot{q}_'+str(idof))
            filedataEN = filedataEN.replace('qddd['+str(idof)+']', '\dddot{q}_'+str(idof))
            filedataEN = filedataEN.replace('qdddd['+str(idof)+']', '\ddddot{q}_'+str(idof))
        for idof in range(0,self.nbrdep):
            filedataEN = filedataEN.replace('p['+str(idof)+']', 'p_'+str(idof))
            filedataEN = filedataEN.replace('pd['+str(idof)+']', '\dot{p}_'+str(idof))
            filedataEN = filedataEN.replace('pdd['+str(idof)+']', '\ddot{p}_'+str(idof))
            filedataEN = filedataEN.replace('pddd['+str(idof)+']', '\dddot{p}_'+str(idof))
            filedataEN = filedataEN.replace('pdddd['+str(idof)+']', '\ddddot{p}_'+str(idof))
        filedataEN = filedataEN.replace(' pi', ' \pi')
        # Write the file out again
        with open('Report_'+str(self.ApplicationFileName)+'.tex', 'w') as file:
            file.write(filedataEN) 
    # ------ end ExportUK_Latex_Report 
    # ------ begin DefineInput 
    def DefineInput(self, input=" "):
        self.nbrinput=self.nbrinput+1
        self.input.append(input)
    # ------ end DefineInput
    # ------ begin DefineOutput 
    def DefineOutput(self, output=" "):
        self.nbroutput=self.nbroutput+1
        self.output.append(output)
    # ------ end DefineOutput
    # ------ begin DefineExternalInput 
    def DefineExternalInput(self, externalinput=" "):
        self.nbrexternalinput=self.nbrexternalinput+1
        self.externalinput.append(externalinput)
    # ------ end DefineExternalInput
    # ------ begin MBSClass.ExportEasyDynProgram
    def ExportEasyDynProgram(self):
        # Some verifications before exporting
        if(self.EASYDYNPlusPlusFlag==0): # EasyDyn classic
            for ibody in range(0,self.nbrbody):
                if self.body[ibody].mass<=0:
                    print("\nError in body ["+str(ibody)+"] : mass="+str(self.body[ibody].mass)+" <=0 !")
                    time.sleep(10)
                    exit("Error")
                if self.body[ibody].Ixx<=0:
                    print("\nError in body ["+str(ibody)+"] : Ixx="+str(self.body[ibody].Ixx)+" <=0 !")
                    time.sleep(10)
                    exit("Error ")
                if self.body[ibody].Iyy<=0:
                    print("\nError in body ["+str(ibody)+"] : Iyy="+str(self.body[ibody].Iyy)+" <=0 !")
                    time.sleep(10)
                    exit("Error ")
                if self.body[ibody].Izz<=0:
                    print("\nError in body ["+str(ibody)+"] : Izz="+str(self.body[ibody].Izz)+" <=0 !")
                    time.sleep(10)
                    exit("Error ")
                if self.body[ibody].Ixx>(self.body[ibody].Iyy+self.body[ibody].Izz):
                    print("\nError in body ["+str(ibody)+"] : Ixx> (Iyy+Izz) ! ")
                    time.sleep(10)
                    exit("Error ")
                if self.body[ibody].Iyy>(self.body[ibody].Ixx+self.body[ibody].Izz):
                    print("\nError in body ["+str(ibody)+"] : Iyy> (Ixx+Izz) ! ")
                    time.sleep(10)
                    exit("Error ")
                if self.body[ibody].Izz>(self.body[ibody].Ixx+self.body[ibody].Iyy):
                    print("\nError in body ["+str(ibody)+"] : Izz> (Ixx+Iyy) ! ")
                    time.sleep(10)
                    exit("Error ")
        CppFileName=self.ApplicationFileName+".cpp"
        print("Creating EasyDyn Application "+CppFileName)
        CppFile = open(CppFileName, "w")
        if(self.EASYDYNPlusPlusFlag==0): # EasyDyn classic
            self.ExportEasyDynHeader(CppFile)
        else:
            self.ExportEasyDynPlusPlusHeader(CppFile)
        if(self.EASYDYNPlusPlusFlag==0):    
            self.ExportSetInertiaData(CppFile)
        if(self.EASYDYNPlusPlusFlag==0): 
            self.ExportComputeMotion(CppFile)
        else:
            self.ExportComputeMotionPlusPlus(CppFile)
        if(self.EASYDYNPlusPlusFlag==0): 
            self.ExportComputePartialVelocities(CppFile)
        else:
            self.ExportComputePartialVelocitiesPlusPlus(CppFile)
        if(self.EASYDYNPlusPlusFlag==0): # EasyDyn classic
            self.ExportAddAppliedForces(CppFile)
        else:
            self.ExportAddAppliedForcesPlusPlus(CppFile)
        if(self.EASYDYNPlusPlusFlag==0): # EasyDyn classic
            self.ExportComputeResidual(CppFile)
        else:
            self.ExportComputeResidualPlusPlus(CppFile)
        if(self.EASYDYNPlusPlusFlag!=0):    
            self.ExportComputeInputs(CppFile)
        if(self.EASYDYNPlusPlusFlag!=0):    
            self.ExportComputeOutputs(CppFile)
        if(self.EASYDYNPlusPlusFlag!=0):    
            self.ExportDefineExternalInputs(CppFile)
        if(self.EASYDYNPlusPlusFlag==0):
            self.Exportmain(CppFile)
        else:
            self.ExportmainPlusPlus(CppFile)
        CppFile.close()
        
        #if (self.ForcesInAppEff != " " and self.OverWriteForces == 0):
        #    CppFileName=self.ApplicationFileName+".AppEff.cpp"
        #    print("Creating Applied Forces File "+CppFileName)
        #    CppFile = open(CppFileName, "w")
        #    #CppFile.write("// body[x]->MainFrame.R+=(Force)*body[x]->MainFrame.T0F.R.ux();\n")
        #    #CppFile.write("// body[x]->MainFrame.MG+=(Torque)*body[x]->MainFrame..T0F.R.uz();\n")
        #    #CppFile.write("// body[x]->NodeFrame[x].R+=(Force)*body[x]->NodeFrame[x].T0F.R.ux();\n")
        #    #CppFile.write("// body[x]->NodeFrame[x].MG+=(Torque)*body[x]->NodeFrame[x].T0F.R.uz();\n")
        #    CppFile.write(self.ForcesInAppEff)
        #    CppFile.close()
    
        #.plt file to plot extra graphs
        ExtraPlotFlag=0
        for ibody in range(0,self.nbrbody):
            if(self.body[ibody].ExtraResults>0):
                ExtraPlotFlag=ExtraPlotFlag+1
            for iframe in range(0,self.body[ibody].nbrnodes):
                if(self.body[ibody].frame[iframe].ExtraResults>0):
                    ExtraPlotFlag=ExtraPlotFlag+1
        if(ExtraPlotFlag>0):        
            CppFileName=self.ApplicationFileName+"_ExtraResults.plt"
            print("Creating Plt Plot File Extra Results "+CppFileName)
            CppFile = open(CppFileName, "w")
            self.ExportExtraPlot(CppFile)
            CppFile.close()
                                      
        print(" ")
        print("Time: "+str((time.clock() - t0))+" s")
        print("Files Completed !")
        time.sleep(0.5) # delays for 5 seconds to display info screen
    # -------- end MBSClass.ExportEasyDynProgram
