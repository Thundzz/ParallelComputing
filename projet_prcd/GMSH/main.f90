program testmesh
  use modmesh
  implicit none

  type(Mesh)::M
  type(Point),dimension(npmax)::vertex
  type(Triangle),dimension(ntmax)::tr
  type(Edge),dimension(nedmax)::ed

  integer :: i,flag
  integer :: nnode, nelem 
  real(kind=8)  :: x,y,z

  open(20,file='meshprogc.data',status='unknown')
  open(21,file='dualformetis.data',status='unknown')

  Call readmesh(M,"carre.mesh")
  vertex = M%mpoint
  tr = M%mtriangle

  nnode = M%np
  nelem = M%nt

  write(20,*) nnode

  do i = 1, nnode

    flag = vertex(i)%ref
    if (flag /= 0) flag = -1
    write(20,*) flag,vertex(i)%coor(1),vertex(i)%coor(2)

  enddo
write(20,*) nelem
write(21,*) nelem, 1  ! 1 == type d element triangle
  do i=1, nelem
    write(20,*) tr(i)%v
    write(21,*) tr(i)%v
  enddo

  close(20);close(21)

!  call savemesh(M,"out.mesh")
end program testmesh
