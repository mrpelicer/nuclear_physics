
include ../../local_settings

DEPS= 	$(BUILDDIR)/interpolator.o \
		$(BUILDDIR)/bag_model.o \
		$(BUILDDIR)/rmf_walecka.o \
		$(BUILDDIR)/particles.o

all: build sym beta phase_transition phase_diagram

build:
	cd ../../.; $(MAKE)
	$(MKDIR) data/
sym:
	$(CC) -o sym symmetricMatter.cpp  $(DEPS) $(GSL_LIB) $(CERES_LIB) $(WALL) 

beta:
	$(CC) -o beta betaEquilibrium.cpp $(DEPS) $(GSL_LIB) $(CERES_LIB) $(WALL) 

phase_transition: 
	$(CC) -o phase_transition phase_transition.cpp $(DEPS) $(GSL_LIB) $(CERES_LIB) $(WALL) 

phase_diagram: 
	$(CC) -o phase_diagram phase_diagram.cpp $(DEPS) $(GSL_LIB) $(CERES_LIB) $(WALL) 

clean:
	$(RM) sym beta  phase_transition phase_diagram