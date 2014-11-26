module modmesh
  implicit none
  !nombre max de points et de triangles
  integer,parameter::npmax = 1e6
  integer,parameter::nedmax = 5e5
  integer,parameter::ntmax = 6e6

  integer,dimension(0:2,1:2),parameter::iare=reshape((/1,2,0,2,0,1/),(/3,2/)) !(1,2),(2,0),(0,1)
  !type Point
  type Point
     integer::num                    !numero du point
     real(kind=8),dimension(3)::coor !coordonnee du point ; nb : en 2D la 3e coord est nulle
     integer::ref                    !tag associe au point
  end type Point
  !type Arete : sert a priori uniquement pour les aretes de frontieres
  type Edge
     integer::num               !numero de l'arete
     integer,dimension(0:1)::v    !indice des 2 sommets composant l'arete
     integer::numT              !1 tr contenant cet arete (numT = 3*k+i)
     integer::ref               !tag associe a l'arete
     logical::infront           !vrai si l'arete appartient au front, faux sinon
  end type Edge
  !type Triangle
  type Triangle
     integer::num               !numero du triangle
     integer,dimension(0:2)::v    !indice des 3 sommets composant le triangle
     integer,dimension(0:2)::ve   !indice des 3 aretes composant le triangle
     integer::ref               !tag associe au triangle
     integer,dimension(0:2)::adja !tab des voisins du triangles
  end type Triangle
  !type Mesh
  type Mesh
     integer::np,ned,nt                              !nb de point, d'arete, de triangle
     type(Point),dimension(npmax)::mpoint
     type(Triangle),dimension(ntmax)::mtriangle
     type(Edge),dimension(nedmax)::medge
  end type Mesh

contains
  subroutine readmesh(maill,filename)
    implicit none
    type(Mesh),intent(inout)::maill
    character(len=*),intent(in)::filename
    character(len=50)::ch
    real(kind=8)::dtmp
    integer::ver,i


    open(unit=1,file=filename,status="old")
    !
    read(1,*) ch,ver
    print*,"on a lu",trim(ch)
    print*,ver,trim(ch)=="MeshVersionFormatted"
    
    do while(trim(ch)/="End")
       select case(trim(ch))
       case ("MeshVersionFormatted")
          print*,"VERSION READING"
       case ("Dimension")
          print*,"DIMENSION READING"
       case("Vertices")
          read(1,*) maill%np
          print*,"NUMBER OF NODES : ",maill%np
          do i=1,maill%np
             maill%mpoint(i)%num = i
             read(1,*) maill%mpoint(i)%coor,dtmp
             maill%mpoint(i)%ref = int(dtmp)
          end do
       case("Triangles")
          read(1,*) maill%nt
          print*,"NUMBER OF TRIANGLES : ",maill%nt
          do i=1,maill%nt
             read(1,*) maill%mtriangle(i)%v,maill%mtriangle(i)%ref
          end do
       case("Edges")
          read(1,*) maill%ned
          print*,"NUMBER OF EDGES : ",maill%ned
          do i=1,maill%ned
             read(1,*) maill%medge(i)%v,maill%medge(i)%ref
          end do
       end select
 print*,"nombre de lignes ",maill%medge(maill%ned)%ref
       read(1,*) ch
      
    end do
    
    close(1)
  end subroutine readmesh

  subroutine savemesh(maill,filename)
    implicit none
    type(Mesh),intent(in)::maill
    type(Point),dimension(npmax)::vertex
    type(Triangle),dimension(ntmax)::tr
    type(Edge),dimension(nedmax)::ed
    character(len=*),intent(in)::filename
    character(len=50)::ch
    real(kind=8)::dtmp
    integer::ver,i

    vertex = maill%mpoint
    tr     = maill%mtriangle
    ed     = maill%medge

 
    open(unit=1,file=filename,status="unknown",form='formatted')

    write(1,*) "MeshVersionFormatted 2"
    write(1,*) "Dimension 3"
    write(1,*) "Vertices"
    write(1,'(I10)') maill%np
    do i=1,maill%np
       write(1,'(E19.12,2X,E19.12,2X,E19.12,2X,I3)') vertex(i)%coor,vertex(i)%ref
    end do
    if(maill%ned/=0) then
       write(1,*) "Edges"
       write(1,'(I10)') maill%ned
       do i=1,maill%ned
          write(1,*) ed(i)%v,ed(i)%ref
       end do
    end if
    if(maill%nt/=0) then
       write(1,*) "Triangles"
       write(1,'(I10)') maill%nt
       do i=1,maill%nt
          write(1,*) tr(i)%v,tr(i)%ref
       end do
    end if   

    write(1,*) "End"
  end subroutine savemesh
end module modmesh
