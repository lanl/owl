#!/usr/bin/env ruby

require "fileutils"

# ruby scripts
execs = ["x_runf90"]

for f in execs

	FileUtils.rm_f(Dir.home + "/bin/" + f)
	File.chmod(0700, Dir.pwd + "/" + f + ".rb")
	FileUtils.symlink(Dir.pwd+"/" + f + ".rb", Dir.home + "/bin/" + f)
	puts f + ' installed '

end
