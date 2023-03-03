include local_settings
	
SOURCES := $(wildcard $(LIBDIR)/*.cpp)
OBJS := $(patsubst $(LIBDIR)/%.cpp, $(BUILDDIR)/%.o, $(SOURCES))

all: folder_build $(OBJS)

folder_build:	
	$(MKDIR) build/

#build the library
$(BUILDDIR)/%.o : $(LIBDIR)/%.cpp
	$(CC) -c $< -o $@ $(LIBS) 

clean:
	$(RM) $(OBJS)
