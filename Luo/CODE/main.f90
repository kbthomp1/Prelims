program main

  use kinddefs,  only : dp
  use gridtools, only : readgrid_lou, basis_function, face_norm
  use solver,    only : get_lhspo, set_bc, get_soln

implicit none

!==============================================================================
!                        VARIABLE DECLARATIONS
!==============================================================================

integer      :: ndimn,ntype,nelem,npoin,nface,nsteps
integer      :: i,j,nnode,ielem,ipoin
integer      :: ip1,ip2,ip3

real(dp)      :: uinf,vinf,tolerance

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

allocate(lhspo(npoin,npoin))

lhspo = get_lhspo(npoin,nelem,nnode,inpoel,geoel)
 
! Formulate the load vector
call face_norm(rface,coord,bcface,nface,ndimn,npoin)
call rhslap(bcface,rface,uinf,vinf,rhspo)

allocate(phi(npoin))
  
call set_bc(phi,lhspo,rhspo,npoin,bcface)

!==============================================================================
!                         SOLVE MATRIX
!==============================================================================

!call conjgrad(lhspo,rhspo,phi,npoin,nsteps)
!call pjac(lhspo,rhspo,phi,npoin,nsteps,res_dataname)

tolerance = 1.E-1_dp
call GSM(lhspo,rhspo,phi,npoin,tolerance)

!do i=1,npoin
!  write(*,*) phi(i)
!end do

allocate(Vx(npoin))
allocate(Vy(npoin))
allocate(Vt(npoin))

call get_soln(Vx,Vy,Vt,phi,geoel,inpoel,npoin,nelem)

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

!========================GAUSS SEIDEL SOLVER=================================
      SUBROUTINE GSM(A1,B1,X01,N,TOL)

      IMPLICIT NONE

      REAL(dp) :: A1(N,N),B1(N),X01(N)
      REAL(dp) :: X(N),AX(N),r(N)
      REAL(dp) :: NORM0,NORM,SUM1,SUM2
      REAL(dp) :: TOL
      INTEGER :: K=1,N,i,j

  11  DO I=1,N
        SUM1=0.0
        SUM2=0.0
          DO J=1,N
            IF (J.LT.I) SUM1=SUM1+A1(I,J)*X(J)
            IF (J.GT.I) SUM2=SUM2+A1(I,J)*X01(J)
          END DO
        X(I)=(B1(I)-SUM1-SUM2)/A1(I,I)
      END DO

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
