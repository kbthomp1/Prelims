PROGRAM main

use grid
use flowsolver
use output

implicit none

!-------------------------------------------------------------------------------
!///////////////////////////////////////////////////////////////////////////////
! Author: Kyle Thompson
!
! This program solves the Euler Equations for subsonic and supersonic flow. At
! this time, no limiter has been implemented, so supersonic flow simulation 
! should only be done with scheme=1 (first-order FV scheme) to avoid instability.
! 
! OUTPUT:
! All output is currently set to write to two files "test.dat" for the main flow
! solution, and "surface.dat" for the surface solution.  There is limited ability
! for tecplot animation.  In the namelist, animation must be set to 1, and the
! filenames tec_dataname and surface_dataname should be given, preferably in a
! separate folder, to avoid overwriting the final output.  The information is 
! written to a new zone at a frequency of iterations set by the control 
! animation_freq in the namelist. NOTE: all solutions are steady solutions, so
! the animation shown is not indicative of a transient flow, this is only meant
! to visualize information propagation in the flow solution (and because I think
! it looks cool). Also, convergence study information is written to the file
! "convergence.dat" and contains the average mesh spacing and entropy generation
!  as well as the gridfile that was run.
!
! INPUT:
! Most of the input parameters in the namelist (euler_input.nml) are 
! self-explanatory, but attention should be given to the parameters debug, pout
! and tolerance.  debug turns on the nsteps control (also specified in the
! namelist), which is the number of timesteps run.  NOTE: if debug is "On" then
! the solution will stop at nsteps without warning. DO NOT TURN ON DEBUG WHEN 
! EXPECTING A CONVERGED SOLUTION.  pout controls the amount of iterations between
! onscreen writeouts of the current L2 error norm. This number should be set to
! >500 for large grids to improve computation speed.  tolerance controls the 
! convergence tolerance to be used to exit the Runge-Kutta time integeration loop
!
!///////////////////////////////////////////////////////////////////////////////
!-------------------------------------------------------------------------------

!===============================================================================
!                           VARIABLE DECLARATIONS
!===============================================================================

real,dimension(:,:),allocatable    :: coord,del,fcent,ecent
real,dimension(:,:,:),allocatable  :: unkno,nodeunkno,rhsel,unold,gradlhs
real,dimension(:),allocatable      :: upoin,vpoin,rhopoin,ppoin,Vtpoin
real,dimension(:),allocatable      :: area,dtlocal,xgrad,ygrad

integer,dimension(:,:),allocatable :: inpoel,esuel,bcface,intfac
integer,dimension(:),allocatable   :: esup1,esup2

integer ndimn,ntype,nelem,npoin,nface,nnode,nfael,numfac,ndegr,neqn
integer i,j,ipoin,n,nsteps,nstag,ip1,ip2,icoun,scheme,pout !,ielem
integer restartcount,animation_freq,animationcount,animation

real g,vinf,rhoinf,minf,alpha,uchar(3),tolerance
real u,v,rho,p,c,e,dt,cfl,L2,L2norm,entropyL2,error,pi,meshspace

character(len=3)   :: debug
character(len=100) :: gridfile,tec_dataname,surface_dataname,tstep,residname
character(len=100) :: iter_str

!===============================================================================
!                             NAMELIST I/O
!===============================================================================

namelist /euler_input/ gridfile,neqn,vinf,rhoinf,minf,alpha,g,      &
                       ndegr,tec_dataname,cfl,nstag,nsteps,tstep,   &
                       tolerance,surface_dataname,debug,scheme,     &
                       residname,pout,animation_freq,animation


open(11,file="euler_input.nml",status="old")

read(11,nml=euler_input)

close(11)

if(animation==1) then
  !delete old output file and open new one (for animation purposes)
  open(23,file=trim(tec_dataname)//".dat",status="replace")
  open(24,file=trim(surface_dataname)//".dat",status="replace")
  write(23,*) 'TITLE = "',trim(tec_dataname),'"'
  write(23,*) 'VARIABLES = "X" "Y" "u" "v" "rho" "p"'
  write(24,*) 'TITLE = "',trim(surface_dataname),'"'
  write(24,*) 'VARIABLES = "X" "Y" "Velocity" "rho" "p"'
  close(23)
  close(24)
end if

!===============================================================================
!                             READ IN GRID FILE
!===============================================================================

call readgrid_lou(ndimn,ntype,nelem,npoin,nface,inpoel,     &
                        gridfile,coord,bcface)
nnode = ntype  !number of nodes per element
nfael = ntype  !number of faces per element

uchar(1)=vinf
uchar(2)=rhoinf
uchar(3)=minf

!define pi
pi=4.0*atan(1.0)

!convert alpha from degrees to radians
alpha=alpha*pi/180.0

allocate(unkno(ndegr,neqn,nelem+nface))
allocate(unold(ndegr,neqn,nelem+nface))
allocate(rhsel(ndegr,neqn,nelem+nface))
allocate(gradlhs(2,2,nelem+nface))
allocate(nodeunkno(ndegr,neqn,npoin))
allocate(upoin(npoin))
allocate(vpoin(npoin))
allocate(rhopoin(npoin))
allocate(Vtpoin(npoin))
allocate(ppoin(npoin))
allocate(area(nelem))
allocate(dtlocal(nelem))
allocate(xgrad(nelem))
allocate(ygrad(nelem))

!==============================================================================
!                   CONTRACT TO 1-DIM ARRAY STORAGE
!==============================================================================

call getesup(inpoel,npoin,nelem,nnode,esup1,esup2)
call getesuel(inpoel,esup1,esup2,nfael,nelem,npoin,esuel)
call getintfac(esuel,inpoel,nelem,nfael,numfac,nface,intfac)
call getnormal(intfac,coord,numfac,nface,del)
call getbcflag(esup1,esup2,inpoel,bcface,nface,    &
                       npoin,nfael,intfac)
call getcentroid(nface,numfac,nelem,intfac,inpoel,coord,fcent,ecent)
call getgradlhs(ecent,fcent,intfac,nface,numfac,nelem,gradlhs)
call getcellarea(inpoel,coord,nelem,area)

write(*,*) 'completed cell area calculation'

!==============================================================================
!                       SET INITIAL CONDITIONS
!==============================================================================

call setinit(unkno,intfac,nface,numfac,uchar,alpha,g)
write(*,*) 'completed initial condition initialization'

!Initialize gradients to zero (all values are constant)
unkno(2:ndegr,:,:)=0.0

!==============================================================================
!                 RUNGE-KUTTA TIME INTEGRATION
!==============================================================================

n     = 0
icoun = 0
error = 1.0e10

animationcount = 0

open(18,file=residname,status="replace")
write(18,*) 'VARIABLES = "iteration" "mass residual"'

!Time Loop
do while(error>tolerance)

  animationcount=animationcount+1
  icoun=icoun+1 
  n=n+1
  if(debug=="On") then
    if(n>nsteps) exit
  end if

  unold(1:ndegr,1:neqn,1:nelem) = unkno(1:ndegr,1:neqn,1:nelem)

  call gettstep(intfac,cfl,nface,numfac,unkno,area,del,nelem,dt,dtlocal,g)
  
  call explrk(neqn,nstag,nelem,numfac,nface,ndegr,dt,dtlocal,unkno,   &
                   unold,rhsel,intfac,uchar,del,ecent,fcent,alpha,    &
                   area,tstep,scheme,gradlhs,g)

  call getL2(unkno,unold,area,nelem,L2)
  
  if(n==1) L2norm=L2
  error=L2/L2norm

  write(18,*) n,error
  
  if(icoun==pout) then
    write(*,*) 'iteration:',n,'L2 error =',error
    icoun=0
  end if
  
  if(animation==1) then
    if(animationcount==animation_freq) then
      call tecoutput(area,unkno,esup1,esup2,npoin,neqn,n,nface,nelem, &
                         ndegr,inpoel,coord,bcface,tec_dataname,surface_dataname,g)
    animationcount=0
    end if 
  end if
end do 

write(*,*) 'Solution converged at n=',n,'iterations'
write(*,*) 'L2 error =',error

!==============================================================================
!                 COMPUTE ENTROPY GENERATION
!==============================================================================

call geterror(unkno,nelem,uchar,area,g,entropyL2)
call getmeshspace(area,del,nface,numfac,nelem,intfac,meshspace)

write(34,*) trim(gridfile),meshspace,entropyL2


!==============================================================================
!                 AREA WEIGHTED AVERAGE FROM CELLS TO NODES
!==============================================================================

call eltonode(area,unkno,esup1,esup2,npoin,neqn,nodeunkno)
 
!==============================================================================
!                      CONVERT TO PRIMITIVE VARIABLES
!==============================================================================

do ipoin=1,npoin
  call ctop_conversion(nodeunkno,u,v,rho,p,c,e,g,ipoin)
  upoin(ipoin)=u
  vpoin(ipoin)=v
  ppoin(ipoin)=p
  rhopoin(ipoin)=rho
  Vtpoin(ipoin)=sqrt(u**2+v**2)
end do

!==============================================================================
!                          TECPLOT MAIN OUTPUT
!==============================================================================

open(21,file="test.dat",status='replace')

write(21,*) 'TITLE = "',trim(tec_dataname),'"'
write(21,*) 'VARIABLES = "X" "Y" "u" "v" "rho" "p"'
write(21,*) 'ZONE NODES=',npoin,",ELEMENTS=",nelem,",DATAPACKING=POINT,",  &
            "ZONETYPE=FETRIANGLE"
do i=1,npoin
  write(21,*) coord(1,i),coord(2,i),upoin(i),vpoin(i),rhopoin(i),ppoin(i)
end do

do j=1,nelem
  write(21,*) (inpoel(i,j),i=1,3)
end do

!==============================================================================
!                          TECPLOT MAIN OUTPUT
!==============================================================================

open(22,file="surface.dat",status='replace')

write(22,*) 'TITLE = "',trim(tec_dataname),'"'
write(22,*) 'VARIABLES = "X" "Y" "Velocity" "rho" "p"'

do i=1,nface
  if(bcface(3,i)==2) then
    ip1=bcface(1,i)
    ip2=bcface(2,i)
    
    write(22,*) coord(1,ip1),coord(2,ip1),Vtpoin(ip1),rhopoin(ip1),ppoin(ip1)
    write(22,*) coord(1,ip2),coord(2,ip2),Vtpoin(ip2),rhopoin(ip2),ppoin(ip2)
  end if
end do

close(21)
close(22)
close(14)
close(15) 
close(16)
close(17)
close(18)

end
