CC = g++
RM = rm -r
LIBDIR = ../../include
MKDIR = mkdir -p
BUILDDIR = ../../include/build

WALL = -std=c++17 -Wall -lm -O3
GSL_LIB= -lgsl -lgslcblas
CERES_LIB= -I/usr/include/eigen3 -lceres -lglog 

OBJS = 	$(BUILDDIR)/particles.o \
				$(BUILDDIR)/rmf_non_linear_walecka.o \
				$(BUILDDIR)/interpolator.o

APP = hmg

$(APP) : $(BUILDDIR) \
	$(OBJS) \
	homogeneous.cpp
	$(CC) -o $(APP) homogeneous.cpp $(OBJS) $(GSL_LIB) $(CERES_LIB) $(WALL) 

$(BUILDDIR):
	$(MKDIR) $(BUILDDIR)

$(BUILDDIR)/%.o: $(SRCDIR)/%.cpp $(OBJDIR)/%.d | build
	$(cc) $(CXXFLAGS) -c $< -o $@

$(BUILDDIR)/particles.o : $(LIBDIR)/particles.hpp \
	$(LIBDIR)/particles.cpp \
	$(LIBDIR)/particles.hpp
	$(CC) -c $(LIBDIR)/particles.cpp -o $(BUILDDIR)/particles.o $(GSL_LIB) $(CERES_LIB) 

$(BUILDDIR)/rmf_non_linear_walecka.o : $(LIBDIR)/rmf_walecka.hpp \
	$(LIBDIR)/rmf_walecka.cpp \
	$(LIBDIR)/rmf_walecka.hpp
	$(CC) -c $(LIBDIR)/rmf_walecka.cpp -o $(BUILDDIR)/rmf_non_linear_walecka.o $(GSL_LIB) $(CERES_LIB)  $(WALL)


$(BUILDDIR)/interpolator.o : $(LIBDIR)/interpolator.hpp \
	$(LIBDIR)/interpolator.cpp \
	$(LIBDIR)/interpolator.hpp
	$(CC) -c $(LIBDIR)/interpolator.cpp -o $(BUILDDIR)/interpolator.o $(GSL_LIB)  $(WALL) 

clean:
	$(RM) $(OBJS) $(APP)
