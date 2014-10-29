/* small prog to create Mesh.Data from scotch informations */

#include <stdio.h>

main (int argc, char** argv)
{
FILE *InFile1,*InFile2,*OutFile;

int Nodes, Elements, I,J,K,L ;
float x,y;

InFile1=fopen("meshprogc.data","r");
InFile2=fopen("dualforscotch.map","r");
OutFile=fopen("Mesh.Data","w");

fscanf(InFile1,"%d\n",&Nodes);
fprintf(OutFile,"%d\n",Nodes);

for (I=0; I<Nodes; I++) {
   fscanf(InFile1,"%d %f %f\n", &K, &x,&y);
   fprintf(OutFile,"%d %f %f\n", K, x,y);
}

fscanf(InFile1,"%d\n",&Elements);
fprintf(OutFile,"%d\n",Elements);

for (I=0; I<Elements; I++) {
   fscanf(InFile1,"%d %d %d\n", &J, &K,&L);
   fprintf(OutFile,"%d %d %d\n", J,K,L);
}

fscanf(InFile2,"%d\n",&Elements);
for (I=0; I<Elements; I++) {
   fscanf(InFile2,"%d %d\n", &J, &K);
   fprintf(OutFile,"%d\n", K );
}




fclose(InFile1);
fclose(InFile2);
fclose(OutFile);


}
