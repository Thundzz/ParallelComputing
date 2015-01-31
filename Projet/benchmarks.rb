#!/usr/bin/env ruby

def squash_files(i)
	1.upto(i) do |k|
		l = k-1
		system("cat sol#{l} >> sol0_#{i}")
	end
end

def speedup_gnuplot()
	input = "speedup.dat"

	File.open("speedup.plt", 'w+') do |f|  
		f.puts("set xlabel \"Nombre de processus\"\n"+
				"set ylabel \"Speedup\"\n"+
				"set term png\n"+
				"set output \"speedup.png\"\n"+
				"set key on inside bottom right\n"+
				"plot \""+input+"\" using 1:2 with linespoints title \"Speedup\"")
	end

	system("gnuplot speedup.plt")
end

def speedup(n)
	system ("rm -f speedup.dat")
	time = {}
	1.upto(n) do |k|
		cmd = "mpiexec -np #{k} ./main.out"
		p cmd
		system(cmd)
		File.open("time#{k}", "r") do |file|
			file.each_line do |line|
				time[k] = line.to_s.split(" ")[0].to_f
			end
		end
		squash_files(k)
	end
	seq = time[1]
	for i in 1..n
		time[i] = seq/time[i]
	end
	File.open("speedup.dat", "a+") do |file|
		1.upto(n) do |k|
			file.puts("#{k} #{time[k]}")
		end
	end
end

speedup ARGV[0].to_i
speedup_gnuplot()