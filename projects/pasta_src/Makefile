
include ../../local_settings

DEPS= 	$(BUILDDIR)/interpolator.o \
		$(BUILDDIR)/pasta.o \
		$(BUILDDIR)/rmf_walecka.o \
		$(BUILDDIR)/particles.o

all: build gibbs_yp

build:
	cd ../../.; $(MAKE)
	$(MKDIR) data/

gibbs_yp:
	$(CC) -o gibbs_yp gibbs_src_yp.cpp $(DEPS) $(GSL_LIB) $(CERES_LIB) $(CUBA_LIB) $(WALL)

clean:
	$(RM) gibbs_yp