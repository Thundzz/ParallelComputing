#!/usr/bin/env ruby

def squash_files(i)
	1.upto(i-1) do |k|
		system("cat sol#{k} >> sol0")
	end
end

def print_plotfile(file)
	system("gnuplot #{file}")
end

def speedup(n)
	time = {}
	1.upto(n) do |k|
		cmd = "/usr/lib64/openmpi/bin/mpiexec -np #{k} ./main.out"
		p cmd
		system(cmd)
		File.open("time#{k}", "r") do |file|
			file.each_line do |line|
				time[k] = line.to_s.split(" ")[0].to_f
			end
		end
	end
	seq = time[1]
	for i in 1..n
		time[i] = seq/time[i]
	end
	File.open("speedup", "w") do |file|
		1.upto(n) do |k|
			file.puts("#{k} #{time[k]}")
		end
	end
end

plotfile = "gnuplot2.cmd" 

#squash_files ARGV[0].to_i
#print_plotfile plotfile
speedup ARGV[0].to_i