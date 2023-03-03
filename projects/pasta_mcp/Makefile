
include ../../local_settings

DEPS= 	$(BUILDDIR)/interpolator.o \
		$(BUILDDIR)/pasta.o \
		$(BUILDDIR)/rmf_walecka.o \
		$(BUILDDIR)/particles.o

all: build mcp

build:
	cd ../../.; $(MAKE)
	$(MKDIR) data/

mcp:
	$(CC) -o mcp mcp.cpp $(DEPS) $(GSL_LIB) $(CERES_LIB) $(CUBA_LIB) $(WALL)
	
clean:
	$(RM) hmg