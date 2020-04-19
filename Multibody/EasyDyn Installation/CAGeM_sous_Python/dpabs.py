#!/usr/bin/python3

from cagem import *

dp3=MBSClass(nbrbody=2, nbrdof=2, nbrdep=0, ApplicationTitle="Simulation of a double pendulum",\
    ApplicationFileName="dpabs")
q,p,pexpr,t=dp3.Getqpt()

dp3.SetGravity(0,-9.81,0)

l0=1.2
l1=1.1

# Inertia characteristics
dp3.body[0].Set(mass=1.1,Ixx=1,Iyy=1,Izz=l0*l0/12*1.1)
dp3.body[1].Set(mass=0.9,Ixx=1,Iyy=1,Izz=l1*l1/12*0.9)

# Definition of the position matrices
dp3.body[0].T0F=Trotz(q[0]) * Tdisp(0,-l0/2,0)
dp3.body[1].T0F=Trotz(q[0]) * Tdisp(0,-l0/2,0) * Tdisp(0,-l0/2,0) * Trotz(q[1]) * Tdisp(0,-l1/2,0)

# Initial conditions
dp3.qini[1]=1

# Simulation parameters
dp3.SetIntegrationParameters(tfinal=5, hsave=0.01, hmax=0.005)

dp3.PrintSystemSize()

dp3.ComputeKinematics()
dp3.PrintKinematics()

# Set STATIC to 1 in case you want CAGeM to generate the code
# to search for static equilibrium before integration
STATIC=0
# Set POLE to 1 if you want to perform an eigen value analysis
POLE=0
#Set TEST to 1 if you want to perform the efficiency tests
TEST=0

dp3.EasyDynFlags(STATIC,POLE,TEST)
dp3.ExportEasyDynProgram()

# Uncomment if you want the LaTeX report in French
dp3.ExportFR_Latex_Report()
# Uncomment if you want the LaTeX report in English
dp3.ExportUK_Latex_Report()
#Uncomment to plot the evolution of position, velocity and #acceleration
dp3.ExportGnuplotScript()
