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

OBJS = 	$(BUILDDIR)/pasta.o \
				$(BUILDDIR)/particles.o \
				$(BUILDDIR)/rmf_non_linear_walecka.o \
				$(BUILDDIR)/interpolator.o

APP= cld


all: $(APP)

$(APP) : $(BUILDDIR) \
	$(OBJS) \
	pasta_liquid_drop.cpp
	$(CC) -o $(APP)  pasta_liquid_drop.cpp $(OBJS) $(GSL_LIB) $(CERES_LIB) $(EIGEN_LIB) $(CUBA_LIB) $(WALL) 
	
$(BUILDDIR):
	$(MKDIR) $(BUILDDIR)


$(BUILDDIR)/pasta.o : $(LIBDIR)/pasta.h \
	$(LIBDIR)/pasta.cpp \
	$(LIBDIR)/pasta.h
	$(CC) -c $(LIBDIR)/pasta.cpp -o $(BUILDDIR)/pasta.o $(GSL_LIB) $(CERES_LIB) $(WALL)


$(BUILDDIR)/particles.o : $(LIBDIR)/particles.h \
	$(LIBDIR)/particles.cpp \
	$(LIBDIR)/particles.h
	$(CC) -c $(LIBDIR)/particles.cpp -o $(BUILDDIR)/particles.o $(GSL_LIB) $(CERES_LIB)  $(WALL)


$(BUILDDIR)/rmf_non_linear_walecka.o : $(LIBDIR)/rmf_non_linear_walecka.h \
	$(LIBDIR)/rmf_non_linear_walecka.cpp \
	$(LIBDIR)/particles.h
	$(CC) -c $(LIBDIR)/rmf_non_linear_walecka.cpp -o $(BUILDDIR)/rmf_non_linear_walecka.o $(GSL_LIB) $(CERES_LIB)  $(WALL)


$(BUILDDIR)/interpolator.o : $(LIBDIR)/interpolator.h \
	$(LIBDIR)/interpolator.cpp \
	$(LIBDIR)/interpolator.h
	$(CC) -c $(LIBDIR)/interpolator.cpp -o $(BUILDDIR)/interpolator.o $(GSL_LIB)  $(WALL) 

clean:
	$(RM) $(OBJS) $(APP)