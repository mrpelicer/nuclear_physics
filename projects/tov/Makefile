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
		$(BUILDDIR)/tov_solver.o

APP = tov

$(APP) : $(BUILDDIR) \
	$(OBJS) \
	tov_solver.cpp
	$(CC) -o $(APP) tov_solver.cpp $(OBJS) $(GSL_LIB) $(WALL) 

$(BUILDDIR):
	$(MKDIR) $(BUILDDIR)


$(BUILDDIR)/interpolator.o : $(LIBDIR)/interpolator.h \
	$(LIBDIR)/interpolator.cpp \
	$(LIBDIR)/interpolator.h
	$(CC) -c $(LIBDIR)/interpolator.cpp -o $(BUILDDIR)/interpolator.o $(GSL_LIB)  $(WALL) 


$(BUILDDIR)/tov_solver.o : $(LIBDIR)/tov_solver.h \
	$(LIBDIR)/tov_solver.cpp \
	$(LIBDIR)/tov_solver.h
	$(CC) -c $(LIBDIR)/tov_solver.cpp -o $(BUILDDIR)/tov_solver.o $(GSL_LIB)  $(WALL) 

clean:
	$(RM) $(OBJS) $(APP)