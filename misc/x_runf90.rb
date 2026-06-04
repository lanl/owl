#!/usr/bin/env ruby

if ARGV.length == 0

	puts ""
	puts " RUNF90 -- compile and run a F90 file  "
	puts ""
	puts " Usage: "
	puts "   x_runf90 f90file "
	puts " "

else

	bn = File.basename(ARGV[0], ".f90")

	system "rm -rf ./obj makefile_#{bn} exec_#{bn} "

	f = File.open('./makefile_' + bn, 'w+')
	f.puts ''
	f.puts '# paths '
	f.puts 'bindir = $(PWD)'
	f.puts 'objdir = ./obj'
	f.puts 'moddir = ./obj'
	f.puts ''
	f.puts '# dependencies'
	f.puts 'object = ' + bn + '.o'
	f.puts ''
	f.puts 'obj = $(addprefix $(objdir)/, $(object))'
	f.puts ''
	f.puts '# targets'
	f.puts 'exec = $(bindir)/exec_' + bn
	f.puts ''
	f.puts 'all: makedir $(exec)'
	f.puts ''
	f.puts '# options -- libflit.a must proceed '
	f.puts 'include $(HOME)/src/flit/src/Makefile.in'
	f.puts 'inc = $(base_inc) -I$(HOME)/src/flit/lib -I$(HOME)/src/rgm/lib -I$(HOME)/src/geof/lib '
	f.puts 'lflags = $(HOME)/src/rgm/lib/librgm.a $(HOME)/src/geof/lib/libgeof.a $(HOME)/src/flit/lib/libflit.a $(base_lflags)'
	f.puts 'fflags = $(base_fflags) '
	f.puts 'cflags = $(base_cflags)'
	f.puts 'cxxflags = $(base_cxxflags)'
	f.puts ''
	f.puts '# compile'
	f.puts '$(objdir)/%.o : ./%.f90'
	f.puts '	$(fc) -o $@ -c $(fflags) $(inc) $<'
	f.puts ''
	f.puts '# link'
	f.puts '$(exec) : $(obj)'
	f.puts '	$(fc) -o $@ $^ $(lflags) $(inc)'
	f.puts ''
	f.puts '# make directory'
	f.puts 'makedir:'
	f.puts '	-@mkdir -p $(bindir)'
	f.puts '	-@mkdir -p $(objdir)'
	f.puts '	-@mkdir -p $(moddir)'
	f.puts ''
	f.puts '# clean'
	f.puts 'clean:'
	f.puts '	-@rm -rf $(objdir)/*.o $(moddir)/*.mod '
	f.close

	system "make -f makefile_#{bn} "
	rest_args = ARGV[1..].join(" ")
	system "./exec_#{bn} #{rest_args}"
	system "rm -rf ./obj makefile_#{bn} exec_#{bn} "

end
