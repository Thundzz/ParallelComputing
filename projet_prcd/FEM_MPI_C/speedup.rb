#!/usr/bin/env ruby

def gen_datafiles(max)
  for i in 1..max
    system("rm Data*")
    system("rm Sol*.plt")
    system("rm time.dat")
    gen_meshdata(i)
    system("./preprocess.exe #{i}")
    folder  = "metis/#{i}"
    system("rm -rf " + folder)
    system("mkdir " + folder)
    system("/usr/lib64/openmpi/bin/mpirun -n #{i} fem.exe")
    for j in 0..(i-1)
      fname = "Data0#{j}.In"
      solname = "Sol0#{j}.plt"
      timename = "time.dat"
      system("cp " + fname + " " + folder + "/")
      system("cp " + solname + " " + folder + "/")
      system("cp " + timename + " " + folder + "/")
    end
  end
end

def gen_meshdata(n)
  if(n == 1)
    system("cat MeshFor0.Data > Mesh.Data")
  else
    system("./data2tec.exe")
    system("./partdmesh dualformetis.dat #{n}")
    system("cp meshprogc.data Mesh.Data")
    system("cat dualformetis.dat.epart.#{n} >> Mesh.Data")
  end
end

def calc_speedup(n)
  time = []
  for i in 1..n
    File.open("metis/#{i}/time.dat") do |file|
      file.each do |line|
        time[i] = line
      end
    end
  end
  File.open("metis/speedup.dat", "w+") do |file|
    for i in 1..n
      spdUp = time[1].to_f/time[i].to_f
      file.puts("#{i} #{spdUp}")
    end
  end
end


def print_speedup
  File.open("metis/plot.plt", "w+") do |file|
      file.puts("set terminal png size 800,600")
      file.puts("set output 'metis/plot.png'")
      file.puts("plot 'metis/speedup.dat' with linespoints")
  end
  system("gnuplot metis/plot.plt")
  system("eog metis/plot.png")
end

gen_datafiles(6)
calc_speedup(6)
print_speedup()