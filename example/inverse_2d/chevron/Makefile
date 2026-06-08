
# paths
bindir = $(PWD)
objdir = ./obj
moddir = ./obj

# dependencies
object = module_su.o create_test.o

obj = $(addprefix $(objdir)/, $(object))

# targets
exec = $(bindir)/create_test

all: makedir $(exec)

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
$(exec) : $(obj)
	$(fc) -o $@ $^ $(lflags) $(inc)

# make directory
makedir:
	-@mkdir -p $(bindir)
	-@mkdir -p $(objdir)
	-@mkdir -p $(moddir)

# clean
clean:
	-@rm -rf $(objdir)/*.o $(moddir)/*.mod
