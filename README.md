# Two_Phase_Flow_Simulation

A one-dimensional structured two-phase flow solver based on the Baer-Nunziato two-fluid model

This prototype code allows to trigger simulations of a liquid-vapor flow in a one dimensional straight pipeline.

#Prerequisites

-Computer with a Linux distribution as OS
-Having the C++ compiler g++ installed
-Having the basic TeX-Live packages installed (A LateX report is automatically generated at the end of the simulation)
-Having the basic Tikz     packages installed (A LateX report is automatically generated at the end of the simulation)

#Absolute path to fulfill the appropriate library linkage

This Git repository contains an archive "Libraries.tar.gz" of the different library folders that are necessary.
After having cloned the present Git repository to your local repository, open the "makefile" file and copy the absolute path
towards the Untared libraries folder in front of the LIB_DIR variable. 
For example, on my personal computer, the "makefile" contained:

LIB_DIR=/home/iampietro/Documents/Administratif/After_PHD/Alan_Turing_Institute/Solver_Baer_Nunziato/Libraries/

#Permission to execute the RunSimulation.exe file

In a terminal, just type: 

chmod +x RunSimulation.exe

#After having executed the file "RunSimulation.exe"

1. The archives "Libraries.tar.gz", "Figures.tar.gz" and "Report.tar.gz" are untared.
  -"Figures.tar.gz" contains reference solutions obtained by using solvers published in the two-phase flow literature.
  These solutions will be compared, at final time, with the computed solution of the current simulation.
  -"Report.tar.gz" contains 



