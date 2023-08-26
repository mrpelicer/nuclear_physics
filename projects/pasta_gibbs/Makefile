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
				$(BUILDDIR)/pasta.o \
				$(BUILDDIR)/particles.o \
				$(BUILDDIR)/rmf_non_linear_walecka.o
			

APP = pasta
OUT1= ypfixed_pasta
OUT2= betaeq_pasta

$(APP) : $(BUILDDIR) \
	$(OBJS) \
	pasta_gibbs.cpp \
	pasta_betaeq.cpp
	$(CC) -o $(OUT1) pasta_gibbs.cpp $(OBJS) $(GSL_LIB) $(CERES_LIB)  $(CUBA_LIB) $(WALL) 
	$(CC) -o $(OUT2) pasta_betaeq.cpp $(OBJS) $(GSL_LIB) $(CERES_LIB)  $(CUBA_LIB) $(WALL) 

$(BUILDDIR):
	$(MKDIR) $(BUILDDIR)


$(BUILDDIR)/pasta.o : $(LIBDIR)/particles.hpp \
	$(LIBDIR)/pasta.cpp \
	$(LIBDIR)/pasta.hpp
	$(CC) -c $(LIBDIR)/pasta.cpp -o $(BUILDDIR)/pasta.o $(GSL_LIB) $(CERES_LIB) $(CUBA_LIB) $(WALL)


$(BUILDDIR)/particles.o : $(LIBDIR)/particles.hpp \
	$(LIBDIR)/particles.cpp \
	$(LIBDIR)/particles.hpp
	$(CC) -c $(LIBDIR)/particles.cpp -o $(BUILDDIR)/particles.o $(GSL_LIB) $(CERES_LIB)  $(WALL)


$(BUILDDIR)/rmf_non_linear_walecka.o : $(LIBDIR)/rmf_walecka.hpp \
	$(LIBDIR)/rmf_walecka.cpp \
	$(LIBDIR)/particles.hpp
	$(CC) -c $(LIBDIR)/rmf_walecka.cpp -o $(BUILDDIR)/rmf_non_linear_walecka.o $(GSL_LIB) $(CERES_LIB)  $(WALL)


$(BUILDDIR)/interpolator.o : $(LIBDIR)/interpolator.hpp \
	$(LIBDIR)/interpolator.cpp \
	$(LIBDIR)/interpolator.hpp
	$(CC) -c $(LIBDIR)/interpolator.cpp -o $(BUILDDIR)/interpolator.o $(GSL_LIB)  $(WALL) 

clean:
	$(RM) $(OBJS) $(APP)