#How to use this stuff:

## 1st Step, create the sub-meshes

* ./data2texc.exe -> meshprogc.dat -> dualformetis.dat


* metis-4.0.3/partdmesh dualformetis.dat 4
* *This gives you dualformetis.dat.epart.4*

## Then with METIS

* cp meshprogc.dat Mesh.*D*ata
* cat dualformetis dat.epart.4>> Mesh.Data
* gcc -lm -o preprocess.exe Preprocess.c
* ./preprocess.exe 4


## The output file Data00.In for example should look like:

0
27 -> with Data01
............
25 -> with Data02
............
9  -> with Data03
............


Which gives you the elements that are in common between the data sets, and who they are !!