
all: jfactors_py.so run.exe
clean: rm obj/*.o

OPTS = -O3 -ffast-math -fPIC

BOOSTINCPATH =-I/opt/ioa/software/boost/1.55.0/include
BOOSTLIBPATH =-L/opt/ioa/software/boost/1.55.0/lib
PYTHONINCPATH=-I/opt/ioa/software/python/2.7.8/include/python2.7/
PYTHONLIBPATH=-L/opt/ioa/software/python/2.7.8/lib/
POTPATH=-I/data/jls/tact/tact/pot/inc/ -I/data/jls/tact/tact/general/coordtransforms/inc/ -I/data/jls/tact/tact/general/ -I/opt/ioa/software/gsl/1.16/
POTLIBPATH=-L/data/jls/tact/tact/pot/ -L/opt/ioa/software/gsl/1.16/lib/

AAPATH=/data/jls/tact/tact/aa/
AAINCPATH=$(AAPATH)inc/
AALIBPATH=$(AAPATH)lib/
GTESTPATH=/data/jls/gtest-1.7.0/
GTESTINC=$(GTESTPATH)/include/
GTESTLIB=$(GTESTPATH)/lib/.libs/

DFPATH=../df/inc/
EDDINVPATH=../edd_inv/inc/
DFSRCPATH=../df/src/
EDDINVSRCPATH=../edd_inv/src/

INC = $(PYTHONINCPATH) $(BOOSTINCPATH) $(POTPATH) -Iinc/ -I$(GTESTINC) -I$(EDDINVPATH) -I$(DFPATH) -I$(AAINCPATH)
LIB = $(PYTHONLIBPATH) $(BOOSTLIBPATH) $(POTLIBPATH) -L$(AALIBPATH)

OBJ = obj
SRC = src

# %.o : %.cpp
	# g++ --std=c++11 -c $< -o $@ $(INC) -fPIC -g

$(OBJ)/%.o : $(SRC)/%.cpp
	g++ --std=c++11 $(OPTS) -c $< -o $@ $(INC)

$(OBJ)/%.o : $(EDDINVSRCPATH)/%.cpp
	g++ --std=c++11 $(OPTS) -c $< -o $@ $(INC) -I$(GTESTINC)

$(OBJ)/%.o : $(DFSRCPATH)/%.cpp
	g++ --std=c++11 $(OPTS) -c $< -o $@ $(INC) -I$(GTESTINC)

OBJECTS = $(OBJ)/density.o $(OBJ)/observed_density.o $(OBJ)/two_component_model.o  $(OBJ)/analytic_results.o $(OBJ)/df_base.o $(OBJ)/moments.o $(OBJ)/edd_inv.o

test.exe: $(OBJECTS) obj/test.o
	g++ --std=c++11 -g $(OPTS) $(OBJECTS) obj/test.o $(LIB) -L$(GTESTLIB) -Wl,-R$(GTESTLIB) -lgtest -lpot_js_nt -lcuba -lboost_python -lpython2.7  -L/home/jls/work/code/general/cuba/ -lgsl -lgslcblas -lgcov -o test.exe

run.exe: $(OBJECTS) obj/run.o
	g++ --std=c++11 -g $(OBJECTS) obj/run.o $(LIB) -lpot_js_nt -lcuba -lboost_python -lpython2.7  -L/home/jls/work/code/general/cuba/ -lgsl -lgslcblas -o run.exe

jfactors_py.so: $(OBJECTS) obj/jfactors_py.o
	g++ --std=c++11 -fPIC -shared $(OBJECTS) obj/jfactors_py.o $(LIB) -laa -lpot_js_nt -lcuba -lboost_python -lpython2.7  -L/home/jls/work/code/general/cuba/ -lgsl -lgslcblas -o jfactors_py.so


# main.exe: jfactors_new.o
# 	g++ --std=c++11 -g jfactors_new.o  $(LIB) -lpot_js_nt -lcuba -lboost_python -lpython2.7  -L/home/jls/work/code/general/cuba/ -lgsl -lgslcblas -o main.exe
# jfactors.so: jfactors.o
# 	g++ --std=c++11 -fPIC -shared jfactors.o $(LIB) -lcuba -lboost_python -lpython2.7  -L/home/jls/work/code/general/cuba/ -lgsl -lgslcblas -o jfactors.so

# jfactors_new.so: jfactors_new.o
# 	g++ --std=c++11 -fPIC -shared jfactors_new.o $(LIB) -lpot_js_nt -lcuba -lboost_python -lpython2.7  -L/home/jls/work/code/general/cuba/ -lgsl -lgslcblas -o jfactors_new.so
