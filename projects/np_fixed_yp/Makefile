CC = g++
RM = rm -r
LIBDIR = ../../include
SOLVERDIR = ../../include/solvers
MKDIR = mkdir -p
BUILDDIR = ../../include/build

WALL = -std=c++17 -Wall -lm -O3
GSL_LIB= -lgsl -lgslcblas
CERES_LIB= -I/usr/include/eigen3 -lceres -lglog 
EIGEN_LIB= -I/usr/include/eigen3/unsupported

CFLAGS = $(WALL) $(GSL_LIB) $(CERES_LIB) $(EIGEN_LIB) $(WALL)


OBJS = 	$(BUILDDIR)/particles.o \
				$(BUILDDIR)/rmf_walecka.o \
				$(BUILDDIR)/interpolator.o \
				
APP = hmg

$(APP) : $(BUILDDIR) \
	$(OBJS) \
	homogeneous.cpp
	$(CC) -o $(APP) homogeneous.cpp $(OBJS) $(CFLAGS)

$(BUILDDIR):
	$(MKDIR) $(BUILDDIR)


$(BUILDDIR)/particles.o : $(LIBDIR)/particles.hpp \
	$(LIBDIR)/particles.cpp \
	$(LIBDIR)/particles.hpp
	$(CC) -c $(LIBDIR)/particles.cpp -o $(BUILDDIR)/particles.o $(CFLAGS)

$(BUILDDIR)/rmf_walecka.o : $(LIBDIR)/rmf_walecka.hpp \
	$(LIBDIR)/rmf_walecka.cpp \
	$(LIBDIR)/rmf_walecka.hpp
	$(CC) -c $(LIBDIR)/rmf_walecka.cpp -o $(BUILDDIR)/rmf_walecka.o $(CFLAGS)


$(BUILDDIR)/interpolator.o : $(LIBDIR)/interpolator.hpp \
	$(LIBDIR)/interpolator.cpp \
	$(LIBDIR)/interpolator.hpp
	$(CC) -c $(LIBDIR)/interpolator.cpp -o $(BUILDDIR)/interpolator.o $(CFLAGS)


$(BUILDDIR)/MesonFieldSolver.o : $(SOLVERDIR)/MesonFieldSolver.hpp \
	$(SOLVERDIR)/MesonFieldSolver.cpp \
	$(SOLVERDIR)/MesonFieldSolver.hpp
	$(CC) -c $(SOLVERDIR)/MesonFieldSolver.cpp -o $(BUILDDIR)/MesonFieldSolver.o $(CFLAGS)

clean:
	$(RM) $(OBJS) $(APP)