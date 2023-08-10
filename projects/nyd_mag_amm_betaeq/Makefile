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
BOOST_LIB= -I/usr/local/include/boost_1_76_0/
CUBA_LIB = -lcuba
PYTHON_LIB= -I/usr/include/python3.6/ -lpython3.6m

OBJS = 	$(BUILDDIR)/particles.o \
				$(BUILDDIR)/rmf_non_linear_walecka.o \
				$(BUILDDIR)/interpolator.o
#  $(BUILDDIR)/core.o \
#  $(BUILDDIR)/crust.o \
#  $(BUILDDIR)/lepton.o \
#  $(BUILDDIR)/nuclear_en.o \
#  $(BUILDDIR)/coulomb.o \
#  $(BUILDDIR)/nuclear_surface_en.o \
#  $(BUILDDIR)/nuclear_matter.o \
#  $(BUILDDIR)/empirical.o 
APP = hmg

$(APP) : $(BUILDDIR) \
	$(OBJS) \
	homogeneous_betaEq.cpp
	$(CC) -o $(APP) homogeneous_betaEq.cpp $(OBJS) $(GSL_LIB) $(CERES_LIB) $(WALL) 

$(BUILDDIR):
	$(MKDIR) $(BUILDDIR)


$(BUILDDIR)/particles.o : $(LIBDIR)/particles.h \
	$(LIBDIR)/particles.cpp \
	$(LIBDIR)/particles.h
	$(CC) -c $(LIBDIR)/particles.cpp -o $(BUILDDIR)/particles.o $(GSL_LIB) $(CERES_LIB)  $(WALL)

$(BUILDDIR)/rmf_non_linear_walecka.o : $(LIBDIR)/rmf_non_linear_walecka.h \
	$(LIBDIR)/rmf_non_linear_walecka.cpp \
	$(LIBDIR)/rmf_non_linear_walecka.h
	$(CC) -c $(LIBDIR)/rmf_non_linear_walecka.cpp -o $(BUILDDIR)/rmf_non_linear_walecka.o $(GSL_LIB) $(CERES_LIB)  $(WALL)


$(BUILDDIR)/interpolator.o : $(LIBDIR)/interpolator.h \
	$(LIBDIR)/interpolator.cpp \
	$(LIBDIR)/interpolator.h
	$(CC) -c $(LIBDIR)/interpolator.cpp -o $(BUILDDIR)/interpolator.o $(GSL_LIB)  $(WALL) 

clean:
	$(RM) $(OBJS) $(APP)