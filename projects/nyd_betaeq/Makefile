
include ../../local_settings

DEPS= 	$(BUILDDIR)/interpolator.o \
		$(BUILDDIR)/rmf_walecka.o \
		$(BUILDDIR)/particles.o

all: build densloop bloop

build:
	cd ../../.; $(MAKE)

densloop:
	$(CC) -o densloop hmg_betaEq.cpp  $(DEPS) $(GSL_LIB) $(CERES_LIB) $(WALL) 

bloop:
	$(CC) -o bloop hmg_betaEq_Bloop.cpp  $(DEPS) $(GSL_LIB) $(CERES_LIB) $(WALL) 

clean:
	$(RM) densloop bloop