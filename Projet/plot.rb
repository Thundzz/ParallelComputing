#!/usr/bin/env ruby

def squash_files(i)
	1.upto(i-1) do |k|
		system("cat sol#{k} >> sol0")
	end
end

def print_plotfile(file)
	system("gnuplot #{file}")
end

plotfile = "gnuplot.cmd" 

# k = ARGV[0].to_i
# cmd = "mpiexec -np #{k} ./main.out"
# p cmd
# system(cmd)
squash_files ARGV[0].to_i
print_plotfile plotfile