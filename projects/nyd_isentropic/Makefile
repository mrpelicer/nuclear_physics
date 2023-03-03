
include ../../local_settings

DEPS= 	$(BUILDDIR)/interpolator.o \
		$(BUILDDIR)/rmf_walecka.o \
		$(BUILDDIR)/particles.o

all: build densloop

build:
	cd ../../.; $(MAKE)
	$(MKDIR) data/

densloop:
	$(CC) -o hmg trapped.cpp  $(DEPS) $(GSL_LIB) $(CERES_LIB) $(WALL) 

clean:
	$(RM) hmg