include ../../local_settings

DEPS= 	$(BUILDDIR)/interpolator.o \
		$(BUILDDIR)/tov_solver.o

all: build tov

build:
	cd ../../.; $(MAKE)
	
tov:
	$(CC) -o tov tov_solver.cpp  $(DEPS) $(GSL_LIB) $(CERES_LIB) $(BOOST_LIB) $(WALL) 

clean:
	$(RM) tov