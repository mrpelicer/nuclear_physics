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

OBJS = 	$(BUILDDIR)/interpolator.o \
				$(BUILDDIR)/bag_model.o \
				$(BUILDDIR)/particles.o

OBJS_PT= $(BUILDDIR)/interpolator.o \
				 $(BUILDDIR)/bag_model.o \
				 $(BUILDDIR)/particles.o \
				 $(BUILDDIR)/rmf_non_linear_walecka.o

APP = sym
APP2= beta
APP3= ptflavor

all: 	$(APP) $(APP2) $(APP3)

$(APP) : $(BUILDDIR) \
	$(OBJS) \
	symmetricMatter.cpp
	$(CC) -o $(APP) symmetricMatter.cpp $(OBJS) $(GSL_LIB) $(CERES_LIB) $(WALL) 

$(APP2) : $(BUILDDIR) \
	$(OBJS) \
	betaEquilibrium.cpp
	$(CC) -o $(APP2) betaEquilibrium.cpp $(OBJS) $(GSL_LIB) $(CERES_LIB) $(WALL) 

$(APP3) : $(BUILDDIR) \
	$(OBJS_PT) \
	phase_diagram.cpp
	$(CC) -o $(APP3) phase_diagram.cpp $(OBJS_PT) $(GSL_LIB) $(CERES_LIB) $(WALL) 


$(BUILDDIR):
	$(MKDIR) $(BUILDDIR)

$(BUILDDIR)/interpolator.o : $(LIBDIR)/interpolator.h \
	$(LIBDIR)/interpolator.cpp \
	$(LIBDIR)/interpolator.h
	$(CC) -c $(LIBDIR)/interpolator.cpp -o $(BUILDDIR)/interpolator.o $(GSL_LIB)  $(WALL) 

$(BUILDDIR)/particles.o : $(LIBDIR)/particles.h \
	$(LIBDIR)/particles.cpp \
	$(LIBDIR)/particles.h
	$(CC) -c $(LIBDIR)/particles.cpp -o $(BUILDDIR)/particles.o $(GSL_LIB) $(CERES_LIB)  $(WALL)


$(BUILDDIR)/bag_model.o : $(LIBDIR)/bag_model.h \
	$(LIBDIR)/bag_model.cpp \
	$(LIBDIR)/particles.h
	$(CC) -c $(LIBDIR)/bag_model.cpp -o $(BUILDDIR)/bag_model.o $(GSL_LIB) $(CERES_LIB)  $(WALL)
	
$(BUILDDIR)/rmf_non_linear_walecka.o : $(LIBDIR)/rmf_non_linear_walecka.h \
	$(LIBDIR)/rmf_non_linear_walecka.cpp \
	$(LIBDIR)/rmf_non_linear_walecka.h
	$(CC) -c $(LIBDIR)/rmf_non_linear_walecka.cpp -o $(BUILDDIR)/rmf_non_linear_walecka.o $(GSL_LIB) $(CERES_LIB)  $(WALL)

clean:
	$(RM) $(OBJS) $(OBJS_PT) $(APP)