
include ../../local_settings

DEPS= 	$(BUILDDIR)/interpolator.o \
		$(BUILDDIR)/rmf_walecka.o \
		$(BUILDDIR)/particles.o

all: build run

build:
	cd ../../.; $(MAKE)
run:
	$(CC) -o hmg homogeneous.cpp  $(DEPS) $(GSL_LIB) $(CERES_LIB) $(WALL) 

clean:
	$(RM) hmg