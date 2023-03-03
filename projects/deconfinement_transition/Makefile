include ../../local_settings

DEPS= 	$(BUILDDIR)/interpolator.o \
		$(BUILDDIR)/tov_solver.o

all: build pt

build:
	cd ../../.; $(MAKE)
	
pt:
	$(CC) -o pt phase_transition.cpp  $(DEPS) $(GSL_LIB) $(CERES_LIB) $(BOOST_LIB) $(WALL) 

clean:
	$(RM) pt