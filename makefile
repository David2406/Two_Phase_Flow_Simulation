#Prototype de makefile
#Auteur: david iampietro

#Structure du dossier dans lequel se trouve le makefile:
#-Tous les fichiers "*.c", "*.cxx", "*.h" se trouve dans le meme dossier
#-Ce dosser contient également un dossier "Libraries" contenant les bibliothèques nécessaires au projet

CC=gcc            #Choix du compilateur C
CXX=g++           #choix du compilateur C++
RM=rm -f          #Variable de suppression de fichiers
FIND=find . -name #Variable pour trouver les fichiers
#CXXFLAGS=-g -Wall #Options de compilation pour CXX 
CXXFLAGS=-g -Wno-deprecated-declarations #Prevent Warnings related to 'deprecated functions in librairies' 
LDFLAGS=-g        #Options pour la phase de "linkage" entre les fichiers "*.o" et les eventuelles bibliothèques

#Tous les fichiers "*.cxx" , attention si par exemple main.cxx contient un " #include "Mesh.h" " alors main.cxx
#doit etre placé à droite de Mesh.cxx sur cette ligne.
SRCS=Param.cxx Mesh.cxx ThermoLaw.cxx Sol_Isentropic.cxx Solver_Isentropic.cxx main.cxx  
#Tous les fichiers ".h"
HEADERS=Param.h Mesh.h ThermoLaw.h Sol_Isentropic.h Solver_Isentropic.h 

OBJS=$(subst .cxx,.o,$(SRCS)) # Variable identique à la variable SRCS sauf que tous les fichiers ont une extension ".o" et non plus ".cxx"

#Chemin absolu menant vers le dossier contenant tous les dossiers bibliothèque
LIB_DIR=/home/iampietro/Documents/Administratif/After_PHD/Alan_Turing_Institute/Solver_Baer_Nunziato/Libraries/
#Si par exemple dans le dossier Libraries il y a un dossier bibliothèque nommé "eigen", je fourni son nom ici
LDLIBS=eigen/   

OUTPUT_DIR=./Output/

$(info library path is [${LIB_DIR}${LDLIBS}])

#Quand le makefile est lancé, il lance la target nommée "Low_Mach"

all: Low_Mach 

#La target Low_Mach applique l'instruction  "$(CXX) $(LDFLAGS) -o Low_Mach  $(OBJS)" 
#Signification de cette instruction: appel du compilateur CXX avec les options de linkage LDFLAGS pour produire un fichier
#éxécutable "Low_Mach" linké aux fichiers contenus dans la variable OBJS (ie tous les fichiers "*.o" du projet
#la dépendance de la target "Low_Mach" est $(OBJS), cela veut dire que l'instruction de cette target sera relancée si la variable dépendance OBJS est modifiée
#c'est à dire si un des fichiers ".o" a été modifié

Low_Mach: $(OBJS)
	$(CXX) $(LDFLAGS) -o Simulation.exe  $(OBJS)  

#Définition d'une règle implicit: cette règle définit pour target tous les fichiers "*.o" ainsi que tous les fichiers contenus dans la variable HEADERS, ie les headers du projet.  
#Son instruction est: $(CXX) -I $(LIB_DIR)$(LDLIBS) $(CXXFLAGS) -c  $< 
#Signification: appel du compilateur CXX en lui fournissant le chemin vers la bibliothèque LDLIBS, en l'occurence eigen; avec les options de compilation
#CXXFLAGS. -c signifie compile. Le compilateur va donc compiler les fichiers "*.cxx". Si ces derniers ont des include vers des fichiers ".h" de la bibliothèque LDLIBS, 
#l'instruction -I $(LIB_DIR)$(LDLIBS) va permettre au compilateur de les trouver. 
#Enfin $< signifie: produit un fichier du meme nom que la target, en l'occurence un fichier "*.o"

%.o: %.cxx  $(HEADERS) 
	$(CXX) -std=c++11 -I $(LIB_DIR)$(LDLIBS) $(CXXFLAGS) -c  $<

#make clean pour supprimer tous les fichiers ".o" générés

clean: 
	$(RM) *.exe *.o *.vtk *.dat *_TIME *.swp .*swn .*.swp 
	find $(OUTPUT_DIR) -name "*" -print0 | xargs -0 rm

