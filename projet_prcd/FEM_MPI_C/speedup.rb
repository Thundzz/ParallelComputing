#!/usr/bin/env ruby

def gen_datafiles(max)
  for i in 2..max
    system("rm Data*")
    system("rm Sol*")
    gen_meshdata(i)
    system("./preprocess.exe #{i}")
    folder  = "metis/#{i}"
    system("rm -rf " + folder)
    system("mkdir " + folder)
    system("/usr/lib64/openmpi/bin/mpirun -n #{i} fem.exe")
    for j in 0..(i-1)
      fname = "Data0#{j}.In"
      solname = "Sol0#{j}.plt"
      system("cp " + fname + " " + folder + "/")
      system("cp " + solname + " " + folder + "/")
    end
  end
end

def gen_meshdata(n)
  system("./data2tec.exe")
  system("./partdmesh dualformetis.dat #{n}")
  system("cp meshprogc.data Mesh.Data")
  system("cat dualformetis.dat.epart.#{n} >> Mesh.Data")
end

gen_datafiles(2)