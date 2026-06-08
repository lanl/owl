
# paths
bindir = $(PWD)
objdir = ./obj
moddir = ./obj

# dependencies
object1 = module_su.o create_test.o
object2 = build_vs.o

obj1 = $(addprefix $(objdir)/, $(object1))
obj2 = $(addprefix $(objdir)/, $(object2))

# targets
exec1 = $(bindir)/create_test
exec2 = $(bindir)/build_vs

all: makedir $(exec1) $(exec2)

# options -- libflit.a must proceed
include $(HOME)/src/flit/src/Makefile.in
inc = $(base_inc) -I$(HOME)/src/flit/lib -I$(HOME)/src/rgm/lib
lflags = $(HOME)/src/rgm/lib/librgm.a $(HOME)/src/flit/lib/libflit.a $(base_lflags)
fflags = $(base_fflags)
ccflags = $(base_ccflags)

# compile
$(objdir)/%.o : ./%.f90
	$(fc) -o $@ -c $(fflags) $(inc) $<

$(objdir)/%.o : ../../../src/lib/%.f90
	$(fc) -o $@ -c $(fflags) $(inc) $<

# link
$(exec1) : $(obj1)
	$(fc) -o $@ $^ $(lflags) $(inc)

$(exec2) : $(obj2)
	$(fc) -o $@ $^ $(lflags) $(inc)

# make directory
makedir:
	-@mkdir -p $(bindir)
	-@mkdir -p $(objdir)
	-@mkdir -p $(moddir)

# clean
clean:
	-@rm -rf $(objdir)/*.o $(moddir)/*.mod
