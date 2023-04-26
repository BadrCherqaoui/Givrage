# Compilateur utilisé
CC=g++



# On choisit comment on compile
CXX_FLAGS = $(DEBUG_FLAG)

# Le nom de l'exécutable
PROG = main_seq

# Les fichiers source à compiler
SRC = main.cc timescheme.cpp operations.cpp goutte.cpp Mesh.cpp captation.cpp

# La commande complète : compile seulement si un fichier a été modifié
main: $(SRC)
	mpic++ -std=c++11 -O2 -o main main.cc goutte.cpp operations.cpp timescheme.cpp Mesh.cpp captation.cpp

$(PROG) : $(SRC)
	g++ -std=c++11 -O2 -o main main.cc goutte.cpp operations.cpp timescheme.cpp Mesh.cpp captation.cpp

run : 
	mpirun --oversubscribe -np 2 ./main 


# Évite de devoir connaitre le nom de l'exécutable
all : $(PROG)

# Supprime l'exécutable, les fichiers binaires (.o) et les fichiers
# temporaires de sauvegarde (~)
clean :
	rm -f *.o *~ ../Affichage/*.dat  $(PROG) 