SRC=$(PWD)/src
ANL=$(PWD)/analyse
OBJ=$(PWD)/obj

CXX=g++ -std=c++11
ARMA_INC=/c11scratch/arijitd/libraries/arma/include
ARMA_LIB=/c11scratch/arijitd/libraries/arma/lib64
CXX_FLAGS=-O3 -Wall -m64 -I$(ARMA_INC) -L$(ARMA_LIB) -larmadillo -fopenmp# ,--start-group $(MKL)/libmkl_intel_lp64.a  $(MKL)/libmkl_core.a -Wl,--end-group -lpthread  -openmp -lmkl_intel_ilp64 -lmkl_intel_thread -lmkl_core -liomp5 -mcmodel large -shared-intel -traceback -ltbb
PROG=2D_Ising_MC

### make all
all: $(PROG) ana

$(PROG): $(OBJ)/main.o $(OBJ)/observables.o $(OBJ)/Ising2dMC.o $(OBJ)/parm.o
	$(CXX) $(OBJ)/main.o $(OBJ)/observables.o $(OBJ)/Ising2dMC.o $(OBJ)/parm.o $(CXX_FLAGS)

$(OBJ)/main.o: $(SRC)/main.cpp
	$(CXX) -c $(SRC)/main.cpp -o $(OBJ)/main.o $(CXX_FLAGS)

$(OBJ)/observables.o: $(SRC)/observables.cpp
	$(CXX) -c $(SRC)/observables.cpp -o $(OBJ)/observables.o $(CXX_FLAGS)

$(OBJ)/Ising2dMC.o: $(SRC)/Ising2dMC.cpp
	$(CXX) -c $(SRC)/Ising2dMC.cpp -o $(OBJ)/Ising2dMC.o $(CXX_FLAGS)

$(OBJ)/parm.o: $(SRC)/parm.cpp
	$(CXX) -c $(SRC)/parm.cpp -o $(OBJ)/parm.o $(CXX_FLAGS)

############################################

ana: $(OBJ)/analyse.o $(OBJ)/StrucFac.o $(OBJ)/parm.o 
	$(CXX) $(OBJ)/analyse.o $(OBJ)/StrucFac.o $(OBJ)/parm.o -o $(PWD)/ana.out $(CXX_FLAGS)

$(OBJ)/analyse.o: $(ANL)/analyse.cpp $(OBJ)/StrucFac.o
	$(CXX) -c $(ANL)/analyse.cpp -o $(OBJ)/analyse.o $(CXX_FLAGS)

$(OBJ)/StrucFac.o: $(ANL)/StrucFac.cpp
	$(CXX) -c $(ANL)/StrucFac.cpp -o $(OBJ)/StrucFac.o $(CXX_FLAGS)

cleanana:
	rm $(OBJ)/analyse.o $(OBJ)/StrucFac.o $(PWD)/neqana.out

clean:
	rm $(OBJ)/* $(PWD)/a.out $(PWD)/ana.out
