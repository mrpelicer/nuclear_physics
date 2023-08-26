CC = g++
RM = rm -r
LIBDIR = ../../include
MKDIR = mkdir -p
BUILDDIR = ../../include/build

WALL = -std=c++14 -Wall -lm -O3
GSL_LIB= -lgsl -lgslcblas
CERES_LIB= -lglog -I/usr/include/eigen3 -lceres
CUBA_LIB = -lcuba
EIGEN_LIB= -I/usr/include/eigen3/unsupported
BOOST_LIB= -I/usr/local/include/boost_1_82_0/
CUBA_LIB = -lcuba
PYTHON_LIB= -I/usr/include/python3.6/ -lpython3.6m

OBJS = 	$(BUILDDIR)/interpolator.o \
			$(BUILDDIR)/tov_solver.o
				
APP = pt

all: 	$(APP) 

$(APP) : $(BUILDDIR) \
	$(OBJS) \
	phase_transition.cpp
	$(CC) -o $(APP) phase_transition.cpp $(OBJS) $(GSL_LIB) $(CERES_LIB) $(EIGEN_LIB) $(BOOST_LIB) $(WALL) 

$(BUILDDIR):
	$(MKDIR) $(BUILDDIR)

$(BUILDDIR)/interpolator.o : $(LIBDIR)/interpolator.hpp \
	$(LIBDIR)/interpolator.cpp \
	$(LIBDIR)/interpolator.hpp
	$(CC) -c $(LIBDIR)/interpolator.cpp -o $(BUILDDIR)/interpolator.o $(GSL_LIB)  $(WALL) 

$(BUILDDIR)/tov_solver.o : $(LIBDIR)/tov_solver.hpp \
	$(LIBDIR)/tov_solver.cpp \
	$(LIBDIR)/tov_solver.hpp
	$(CC) -c $(LIBDIR)/tov_solver.cpp -o $(BUILDDIR)/tov_solver.o $(GSL_LIB) $(CERES_LIB)  $(WALL)

clean:
	$(RM) $(OBJS)  $(APP)