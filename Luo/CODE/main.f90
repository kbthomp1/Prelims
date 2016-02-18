program main

  use kinddefs,  only : dp
  use gridtools, only : readgrid_lou, basis_function, face_norm
  use solver,    only : get_lhspo

implicit none

!==============================================================================
!                        VARIABLE DECLARATIONS
!==============================================================================

integer      :: ndimn,ntype,nelem,npoin,nface,nsteps
integer      :: i,j,nnode,ielem,ipoin,ib
integer      :: ip1,ip2,ip3

real(dp)      :: uinf,vinf,phi_ib,v_dummy

real(dp), dimension(:),   allocatable :: rhspo,phi,Vx_local,Vy_local,Vxarea
real(dp), dimension(:),   allocatable :: Vyarea,area,Vx,Vy
real(dp), dimension(:),   allocatable :: Vt

real(dp), dimension(:,:), allocatable :: coord,geoel,lhspo,rface

integer, dimension(:,:), allocatable :: inpoel,bcface

character(len=100)                   :: gridfile,bc_case,tec_dataname
character(len=100)                   :: res_dataname,comp_name

!==============================================================================
!                         NAMELIST I/O 
!==============================================================================

namelist /fe_input/ gridfile,nnode,uinf,vinf,bc_case,nsteps,tec_dataname,  &
                    res_dataname,comp_name
open(11,file="fe_input.nml",status="old")
read(11,nml=fe_input)

close(11)

!==============================================================================
!                    OPEN TROUBLESHOOTING FILE
!==============================================================================

open(30,file="trouble.dat",status="replace")

! Read the grid
call readgrid_lou(ndimn,ntype,nelem,npoin,nface,nnode,inpoel,         &
                  gridfile,coord,bcface)

! assemble geometry related matricies
call basis_function(nelem,geoel,inpoel,coord)

!=============================================================================
!                         GLOBAL STIFFNESS MATRIX
!=============================================================================

allocate(lhspo(npoin,npoin))

lhspo = get_lhspo(npoin,nelem,nnode,inpoel,geoel)
 
open(34,file="A_matrix",status="replace")
open(35,file="A_matrix_diagonals",status="replace")

do j=1,npoin
  write(34,*) (lhspo(i,j),i=1,npoin)
  write(35,*) lhspo(j,j)
end do

close(34)
close(35)

! Formulate the load vector
call face_norm(rface,coord,bcface,nface,ndimn,npoin)
call rhslap(bcface,rface,uinf,vinf,rhspo)

!==============================================================================
!                        IMPOSE DIRCHLET BC
!==============================================================================

allocate(phi(npoin))

phi(1:npoin) = 0.0

select case(bc_case)

  case("channel")
    ib = bcface(1,1)
    phi_ib = 1.0

    lhspo(ib,ib) = lhspo(ib,ib)*1.0e+20
    rhspo(ib)   = lhspo(ib,ib)*phi_ib
  
  case("cylinder")
    call cylpot(coord,npoin,phi_ib,v_dummy,bcface(1,1)) 

    ib = bcface(1,1)

    lhspo(ib,ib) = lhspo(ib,ib)*1.0e+20
    rhspo(ib)   = lhspo(ib,ib)*phi_ib

  case default
    write(*,*) "ERROR - No Case Selected!"

end select

!==============================================================================
!                         SOLVE MATRIX
!==============================================================================

!call conjgrad(lhspo,rhspo,phi,npoin,nsteps)
!call pjac(lhspo,rhspo,phi,npoin,nsteps,res_dataname)
call GSM(lhspo,rhspo,phi,npoin,1.0e-05)

!do i=1,npoin
!  write(*,*) phi(i)
!end do

!==============================================================================
!                      CALCULATE VELOCITIES
!==============================================================================

allocate(Vx_local(nelem))
allocate(Vy_local(nelem))
allocate(Vxarea(npoin))
allocate(Vyarea(npoin))
allocate(area(npoin))
allocate(Vx(npoin))
allocate(Vy(npoin))
allocate(Vt(npoin))

Vx_local(1:nelem) = 0.0
Vy_local(1:nelem) = 0.0
Vxarea(1:npoin)   = 0.0
Vyarea(1:npoin)   = 0.0
area(1:npoin)     = 0.0

do ielem=1,nelem

  ip1=inpoel(1,ielem)
  ip2=inpoel(2,ielem)
  ip3=inpoel(3,ielem)

  Vx_local(ielem) = (geoel(1,ielem)*(phi(ip1)-phi(ip3))       &
                   + geoel(2,ielem)*(phi(ip2)-phi(ip3)))*geoel(5,ielem) 

  Vy_local(ielem) = (geoel(3,ielem)*(phi(ip1)-phi(ip3))       &
                   + geoel(4,ielem)*(phi(ip2)-phi(ip3)))*geoel(5,ielem) 

  Vxarea(ip1) = Vxarea(ip1)+Vx_local(ielem)
  Vxarea(ip2) = Vxarea(ip2)+Vx_local(ielem)
  Vxarea(ip3) = Vxarea(ip3)+Vx_local(ielem)

  Vyarea(ip1) = Vyarea(ip1)+Vy_local(ielem)
  Vyarea(ip2) = Vyarea(ip2)+Vy_local(ielem)
  Vyarea(ip3) = Vyarea(ip3)+Vy_local(ielem)
  
  area(ip1)  = area(ip1)+geoel(5,ielem)
  area(ip2)  = area(ip2)+geoel(5,ielem)
  area(ip3)  = area(ip3)+geoel(5,ielem)

end do

do ipoin=1,npoin
  
  Vx(ipoin) = Vxarea(ipoin)/area(ipoin)
  Vy(ipoin) = Vyarea(ipoin)/area(ipoin)

  Vt(ipoin) = sqrt(Vx(ipoin)**2 + Vy(ipoin)**2)

end do

!==============================================================================
!                        TECPLOT OUTPUT
!==============================================================================

open(21,file=tec_dataname,status='replace')

write(21,*) 'TITLE = "',trim(tec_dataname),'"'
write(21,*) 'VARIABLES = "X" "Y" "phi" "Vx" "Vy" "Vt"'
write(21,*) 'ZONE NODES=',npoin,",ELEMENTS=",nelem,",DATAPACKING=POINT,",  &
            "ZONETYPE=FETRIANGLE"
do i=1,npoin
  write(21,*) coord(1,i),coord(2,i),phi(i),Vx(i),Vy(i),Vt(i)
end do

do j=1,nelem
  write(21,*) (inpoel(i,j),i=1,3)
end do

write(30,*) "RHSPO"
do i=1,npoin
  write(30,*) rhspo(i)
end do

write(30,*)
write(30,*) "LHSPO"
do i=1,npoin
  write(30,*) lhspo(i,i)
end do

close(21)
close(30)

if(trim(bc_case)=="channel") then

open(15,file="channel_surface_lower.dat",status="replace")
open(16,file="channel_surface_upper.dat",status="replace")

do i=1,npoin
  if(bcface(3,i)==2) then
    ip1=bcface(1,i)
    ip2=bcface(2,i)
    
    write(15,*) coord(1,ip1),coord(2,ip1),Vt(ip1)
    write(15,*) coord(1,ip2),coord(2,ip2),Vt(ip2)
  
  else 
    exit
  end if
end do

do i=1,npoin
  if(bcface(3,i)==2) then
    ip1=bcface(1,i)
    ip2=bcface(2,i)
    
    write(16,*) coord(1,ip1),coord(2,ip1),Vt(ip1)
    write(16,*) coord(1,ip2),coord(2,ip2),Vt(ip2)
  end if
end do

end if

!==============================================================================
!                         SUBROUTINES
!==============================================================================

contains

!===============================LOAD VECTOR====================================
subroutine rhslap(bcface,rface,uinf,vinf,rhspo)

integer  bcface(ndimn+1,nface),iface,ip1,ip2
real(dp)  rface(ndimn,nface),uinf,vinf
real(dp)  cface

real(dp), dimension(:), allocatable :: rhspo

allocate(rhspo(npoin))

rhspo(1:npoin) = 0.0

do iface=1,nface
  if(bcface(3,iface) == 4) then
    
    ip1   = bcface(1,iface)
    ip2   = bcface(2,iface)
    
!    itype = bcface(4,iface)
!
!    roinf = uchar(1,itype)
!     uinf = uchar(2,itype)
!     vinf = uchar(3,itype)
!     pinf = uchar(4,itype)
     
    cface = 0.5*(uinf*rface(1,iface) + vinf*rface(2,iface))

    rhspo(ip1) = rhspo(ip1) + cface
    rhspo(ip2) = rhspo(ip2) + cface

  end if
end do

return
end subroutine

!======================CYLINDER POTENTIAL CALCULATION===========================
subroutine cylpot(coord,npoin,phi_ib,v_ib,ibpoin)

integer ibpoin,npoin
real(dp) coord(ndimn,npoin),phi_ib,v_ib,x,y

x=coord(1,ibpoin)
y=coord(2,ibpoin)

phi_ib = (1+(0.5)**2/(x**2+y**2))*x

v_ib   = sqrt(1+2*0.5**2/(x**2+y**2)*(2*y**2/(x**2+y**2)-1)+0.5**4/(x**2+y**2)**2)

end subroutine

!======================CONJUGATE GRADIENT SOLVER=================================
subroutine conjgrad(A,b,x,npoin,nsteps)

real(dp)                   :: A(npoin,npoin),rsold,rsnew,alpha,pAp
real(dp), dimension(npoin) :: b,x,ax,r,p,Ap

integer                   :: npoin,nsteps,n

ax(1:npoin)  = 0.0
rsold        = 0.0

do i=1,npoin

  do j=1,npoin   
    ax(i)=ax(i)+A(i,j)*x(j)
  end do

  r(i)=b(i)-ax(i)
  p(1:npoin)=r(1:npoin)
  rsold=rsold+r(i)*r(i)
end do

do n=1,nsteps

  Ap(1:npoin) = 0.0
  pAp   = 0.0
  rsnew = 0.0      

  do i=1,npoin

    do j=1,npoin
      Ap(i)=Ap(i)+A(i,j)*p(j)
    end do

    do j=1,npoin
      pAp=pAp+p(j)*Ap(j)
    end do
 
    alpha=rsold/pAp
    x(i)=x(i)+alpha*p(i);
    r(i)=r(i)-alpha*Ap(i);
    rsnew=rsnew+r(i)*r(i)
  
    if(sqrt(rsnew)<1e-10) then
      write(*,*) "Solution has converged at,",n,",iterations"
      exit
    end if

    p(i)=r(i)+rsnew/rsold*p(i)
    rsold=rsnew
  end do
end do

end subroutine

!=============================POINT JACOBI SOLVER===============================

subroutine pjac(A,b,x,npoin,nsteps,res_dataname)

real(dp)                   :: A(npoin,npoin),resnorm,dx,ad
real(dp), dimension(npoin) :: b,x,ax,r

integer                   :: npoin,nsteps,n

character(len=100)        :: res_dataname

open(33,file=res_dataname,status="replace")

x(1:npoin) = 0.0

do n=1,nsteps

  resnorm     = 0.0
  dx          = 0.0
  ax(1:npoin) = 0.0

  do i=1,npoin
 
    !write(*,*) i, A(i,i),b(i)
 
    do j=1,npoin   
      ax(i)=ax(i)+A(i,j)*x(j)
    end do
  
    r(i)    = b(i)-ax(i)
    resnorm = resnorm + r(i)*r(i)

    ad = 1/A(i,i)

    dx = r(i)*ad

    x(i)=x(i)+dx
  
  end do
    
  if(sqrt(resnorm)<1.0e-05) then
    write(*,*) "Solution has converged at,",n,"iterations"
    exit
  end if

  write(33,*) n,sqrt(resnorm)
  write(*,*)  n,sqrt(resnorm)

end do

close(33)

end subroutine

!========================SOR SOLVER==========================================
!subroutine SOR(A,b,x,jmax,nsteps)
!
!real(dp)                       :: eps,pi,anorm,anormb,omega
!real(dp), dimension(jmax)      :: b,x
!real(dp), dimension(jmax,jmax) :: A
!
!integer i,j,n,nsteps
!
!pi     = 4.D0*DATAN(1.D0)
!anormb = 0.0
!
!rjac   = 1-pi**2/(2*jmax**2)
!
!do i=2,jmax-1
!  do j=2,jmax-1
!    anormb=anormb+abs(b(i,j))
!  end do
!end do
!
!omega = 1.0
!
!do n=1,nsteps
!  anorm=0.0
!  do i=2,jmax-1
!    do j=2,jmax-1
!      if(mod(i+j,2) == mod(n,2)) then
!        resid(A

!========================GAUSS SEIDEL SOLVER=================================
      SUBROUTINE GSM(A1,B1,X01,N,TOL)

      IMPLICIT NONE

      REAL(dp)::A1(N,N),B1(N),X01(N),X(N),AX(N),r(N),NORM0,NORM,SUM1,SUM2
      REAL   ::TOL
      INTEGER::K=1,N,i,j

      OPEN(2,file="GSM.dat",status="replace")

      WRITE (2,*)'RESULT FOR GAUSS-SEIDEL METHOD'

  11  DO I=1,N
        SUM1=0.0
        SUM2=0.0
          DO J=1,N
            IF (J.LT.I) SUM1=SUM1+A1(I,J)*X(J)
            IF (J.GT.I) SUM2=SUM2+A1(I,J)*X01(J)
          END DO
        X(I)=(B1(I)-SUM1-SUM2)/A1(I,I)
      END DO

      WRITE (2,20)K,(X(I),I=1,N)
  20  FORMAT(2X,I3,3(2X,F9.6))

      NORM    = 0.0 
      ax(1:N) = 0.0

      DO i=1,N
        DO j=1,N   
          ax(i)=ax(i)+A1(i,j)*X(j)
        end do
  
        r(i)    = B1(i)-ax(i)
        NORM = NORM + r(i)*r(i)
      END DO

      IF(K==1) NORM0=sqrt(NORM)
      
      K=K+1

      NORM=sqrt(NORM)/NORM0

      write(*,*) K,NORM

      IF (NORM.LT.TOL) then
        write(*,*) 'The solution has converged at, ',k,'iterations'
        GO TO 12
      END IF

      DO I=1,N
        X01(I)=X(I)
      END DO

      GO TO 11

   12 END SUBROUTINE

end program main
