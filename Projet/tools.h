#ifndef TOOLS_H
#define TOOLS_H


/** Calcule le i et j correspondant au n passé en paramètre, et les affecte
*   aux pointeurs i et j.
*/
void nloc(int *i, int *j, int n, int Nx)


/* Fonction de conditions aux bords (sin x + cos y)
 */
double f( double posx, double posy, double t)

/** Calcule le second membre
*/
void RightHandSide(int N, int Nx, int M, double dx, double dy, double Cx, double Cy,double *RHS)

#endif
