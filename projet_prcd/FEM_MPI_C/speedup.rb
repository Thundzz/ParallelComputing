#!/usr/bin/env ruby

def gen_datafiles(max)
  system("mkdir metis")
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
    timename = "time.dat"
    system("mv " + timename + " " + folder + "/")
    for j in 0..(i-1)
      fname = "Data0#{j}.In" unless j>= 10
      fname = "Data#{j}.In" if j>= 10
      solname = "Sol0#{j}.plt"
      system("mv " + fname + " " + folder + "/")
      system("mv " + solname + " " + folder + "/")
    end
  end
end


def gen_datafiles_scotch(max)
  system("mkdir scotch")
  for i in 1..max
    system("rm Data*")
    system("rm Sol*.plt")
    system("rm time.dat")
    gen_meshdata_scotch(i)
    system("./preprocess.exe #{i}")
    folder  = "scotch/#{i}"
    system("rm -rf " + folder)
    system("mkdir " + folder)
    system("/usr/lib64/openmpi/bin/mpirun -n #{i} fem.exe")
    timename = "time.dat"
    system("mv " + timename + " " + folder + "/")
    for j in 0..(i-1)
      fname = "Data0#{j}.In" unless j>= 10
      fname = "Data#{j}.In" if j>= 10
      solname = "Sol0#{j}.plt"
      system("mv " + fname + " " + folder + "/")
      system("mv " + solname + " " + folder + "/")
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

def gen_meshdata_scotch(n)
  if(n == 1)
    system("cat MeshFor0Scotch.Data > Mesh.Data")
  else
    system("./data2tec.exe")
    system("echo cmplt #{n} | ./gmap dualforscotch.grf")
    system("./fromscotch")
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
      file.puts("plot 'metis/speedup.dat' with linespoints, x ")
  end
  system("gnuplot metis/plot.plt")
  system("eog metis/plot.png")
end


def gen_compute_and_print(n)
  gen_datafiles(n)
  calc_speedup(n)
  print_speedup()
end

def gen_compute_and_print_scotch(n)
  gen_datafiles_scotch(n)
end

gen_compute_and_print_scotch(ARGV[0].to_i)

#gen_datafile_for(33)