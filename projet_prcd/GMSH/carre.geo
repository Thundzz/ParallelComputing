// Commentaires

lc = 1.; //longuer caracteristique

Point(1)={0,0,0,lc};
Point(2)={10,0,0,lc};
Point(3)={10,10,0,lc};
Point(4)={0,10,0,lc};

// Definition lignes

Line(1)={1,2};
Line(2)={2,3};
Line(3)={3,4};
Line(4)={4,1};

// Definition de la Surface

Line Loop(5) = {1,2,3,4};
Plane Surface(0)={5};
