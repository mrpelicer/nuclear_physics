
include ../../local_settings

DEPS= 	$(BUILDDIR)/interpolator.o \
		$(BUILDDIR)/pasta.o \
		$(BUILDDIR)/rmf_walecka.o \
		$(BUILDDIR)/particles.o

all: build densloop sizeloop xloop tloop

build:
	cd ../../.; $(MAKE)

densloop:
	$(CC) -o densloop pasta_densloop.cpp $(DEPS) $(GSL_LIB) $(CERES_LIB) $(CUBA_LIB) $(WALL)

sizeloop:
	$(CC) -o sizeloop pasta_Lloop.cpp $(DEPS) $(GSL_LIB) $(CERES_LIB) $(CUBA_LIB) $(WALL)

xloop:
	$(CC) -o xloop pasta_xloop.cpp $(DEPS) $(GSL_LIB) $(CERES_LIB) $(CUBA_LIB) $(WALL)

tloop:
	$(CC) -o tloop pasta_tloop.cpp $(DEPS) $(GSL_LIB) $(CERES_LIB) $(CUBA_LIB) $(WALL)

clean:
	$(RM) densloop sizeloop xloop tloop