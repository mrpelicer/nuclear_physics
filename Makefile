include local_settings
	
SOURCES := $(wildcard $(LIBDIR)/*.cpp)
OBJS := $(patsubst $(LIBDIR)/%.cpp, $(BUILDDIR)/%.o, $(SOURCES))

all: $(BUILDDIR) $(OBJS)

#build the library
$(BUILDDIR)/%.o : $(LIBDIR)/%.cpp
	$(CC) -c $< -o $@ $(LIBS) 

clean:
	$(RM) $(OBJS)
