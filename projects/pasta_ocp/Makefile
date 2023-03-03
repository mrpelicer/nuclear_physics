
include ../../local_settings

DEPS= 	$(BUILDDIR)/interpolator.o \
		$(BUILDDIR)/pasta.o \
		$(BUILDDIR)/rmf_walecka.o \
		$(BUILDDIR)/particles.o

all: build gibbs_yp gibbs_beta cld_yp

build:
	cd ../../.; $(MAKE)
	$(MKDIR) data/

gibbs_yp:
	$(CC) -o gibbs_yp gibbs_yp.cpp $(DEPS) $(GSL_LIB) $(CERES_LIB) $(CUBA_LIB) $(WALL)

gibbs_beta:
	$(CC) -o gibbs_beta gibbs_betaeq.cpp $(DEPS) $(GSL_LIB) $(CERES_LIB) $(CUBA_LIB) $(WALL)

cld_yp:
	$(CC) -o cld_yp cld_yp.cpp $(DEPS) $(GSL_LIB) $(CERES_LIB) $(CUBA_LIB) $(WALL)

clean:
	$(RM) gibbs_yp gibbs_beta cld_yp