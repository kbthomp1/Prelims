MODULE flowsolver
implicit none

contains

!========================RUNGE-KUTTA SUBROUTINE=========================

subroutine explrk(neqn,nstag,nelem,numfac,nface,ndegr,dt,dtlocal,unkno, &
                   unold,rhsel,intfac,uchar,del,ecent,fcent,alpha,      &
                   area,tstep,scheme,gradlhs,g)
implicit none

real,dimension(:,:,:)  :: unkno,unold,rhsel,gradlhs
real,dimension(:,:)    :: del,fcent,ecent
integer,dimension(:,:) :: intfac

real    :: dt,dtlocal(:),coeff,alpha,alfa,uchar(:),area(:),g
integer :: neqn,nstag,istag,ie,nelem,nface,numfac,ndegr,scheme

real valfa(5,5)

character(len=100) :: tstep

!array of R-K Coefficients

data valfa  /1.0    , 0.0    , 0.0   , 0.0  , 0.0,  &
             0.5    , 1.0    , 0.0   , 0.0  , 0.0,  &
             0.1919 , 0.4390 , 1.0   , 0.0  , 0.0,  &
             0.25   , 0.3333 , 0.5   , 1.0  , 0.0,  &
             0.25   , 0.1666 , 0.375 , 0.5  , 1.0   /

!Loop over stages
select case(tstep)

case("global")
do istag=1,nstag
  alfa = valfa(istag,nstag)

  call getrhs(unkno,rhsel,intfac,neqn,nface,numfac,nelem,  &
              ndegr,uchar,del,ecent,fcent,scheme,alpha,    &
              area,gradlhs,g)

  do ie=1,nelem
    coeff = alfa*dt/area(ie)
    unkno(1:ndegr,1:neqn,ie) = unold(1:ndegr,1:neqn,ie)   &
                               - coeff*rhsel(1:ndegr,1:neqn,ie)
  end do
end do

case("local")
do istag=1,nstag
  alfa = valfa(istag,nstag)

  call getrhs(unkno,rhsel,intfac,neqn,nface,numfac,nelem,   &
              ndegr,uchar,del,ecent,fcent,scheme,alpha,     &
              area,gradlhs,g)

  do ie=1,nelem
    coeff = alfa*dtlocal(ie)/area(ie)
    unkno(1:ndegr,1:neqn,ie) = unold(1:ndegr,1:neqn,ie)   &
                               - coeff*rhsel(1:ndegr,1:neqn,ie)
  end do
end do

end select

end subroutine

!=======================================================================

subroutine getrhs(unkno,rhsel,intfac,neqn,nface,numfac,nelem,  &
                  ndegr,uchar,del,ecent,fcent,scheme,alpha,    &
                  area,gradlhs,g)
implicit none

real,dimension(:,:,:) :: unkno,rhsel,gradlhs
real,dimension(:,:)   :: del,ecent,fcent
real uchar(:),alpha,area(:),g
integer intfac(:,:),neqn,nface,numfac,nelem,ndegr,scheme
integer ielem

!set ghost cell values
call getghostv(unkno,intfac,nface,nelem,uchar,alpha,del,   &
                     fcent,ecent,neqn,scheme,g)

!initialize the RHS
rhsel(1:ndegr,1:neqn,1:nelem+nface)=0.0

!get fluxes
call vlflux(intfac,ecent,fcent,numfac,nface,neqn,nelem,ndegr,  &
                  unkno,rhsel,scheme,gradlhs,g,del,area)

end subroutine

!==============================VAN LEER FVS=============================

subroutine vlflux(intfac,ecent,fcent,numfac,nface,neqn,nelem,ndegr,  &
                  unkno,rhsel,scheme,gradlhs,g,del,area)
implicit none

real g,c1,rho1,u1,v1,vn1,q1,p1,mn1,e1,gp1,gm1,g2m1
real c2,rho2,u2,v2,vn2,q2,p2,e2,mn2,nx,ny,area(:)

real,dimension(neqn)  :: flux,pflux,nflux
real,dimension(:,:,:) :: unkno,rhsel,gradlhs
real,dimension(:,:)   :: del,ecent,fcent

integer icell,jcell,iface,ieqn,neqn,numfac
integer intfac(:,:),nface,nelem,scheme,ndegr
integer i,j,ielem

!for VL fvs, u1 = rho, u2=rho*u, and u3=rho*E (THESE ARE THE CONSERVATIVE VARS)

gp1  = g+1.0
gm1  = g-1.0
g2m1 = g**2.0-1.0

if(scheme==2) then
  call getgrads(gradlhs,numfac,nface,nelem,intfac,unkno,  &
                neqn,ndegr,ecent)
else if(scheme==3) then
  call getgg_grad(unkno,area,nelem,nface,numfac,intfac,del,neqn,ndegr,g)
end if
  
do iface=1,nface+numfac
  if(iface<=nface) then
    jcell = intfac(1,iface)
    icell = intfac(2,iface)
  else
    icell = intfac(1,iface)
    jcell = intfac(2,iface)
  end if

  !compute unit normal vector
  nx=del(1,iface)/del(3,iface)
  ny=del(2,iface)/del(3,iface)

  if(scheme==1) then 
    call ctop_conversion(unkno,u1,v1,rho1,p1,c1,e1,g,icell)
    call ctop_conversion(unkno,u2,v2,rho2,p2,c2,e2,g,jcell)
  else if(scheme==2 .or. scheme==3) then
    if(iface<=nface) then
      call getleftright(unkno,fcent,ecent,neqn,nelem,nface,  &
                        u1,v1,rho1,p1,c1,e1,g,icell,iface)
      call ctop_conversion(unkno,u2,v2,rho2,p2,c2,e2,g,jcell)
    else
      call getleftright(unkno,fcent,ecent,neqn,nelem,nface,  &
                        u1,v1,rho1,p1,c1,e1,g,icell,iface)
      call getleftright(unkno,fcent,ecent,neqn,nelem,nface,  &
                        u2,v2,rho2,p2,c2,e2,g,jcell,iface)
    end if
  else
    write(*,*) 'Incorrect scheme selected'
  end if
    
  !Compute the normal velocity and mach number
  vn1 = u1*nx+v1*ny 
  vn2 = u2*nx+v2*ny
  mn1 = vn1/c1
  mn2 = vn2/c2

  !POSITIVE FLUXES => ui state
  if(mn1>=1.0) then !supersonic outflow
   
    !fluxes come from left state
    pflux(1) = vn1*rho1
    pflux(2) = vn1*rho1*u1+p1*nx
    pflux(3) = vn1*rho1*v1+p1*ny
    pflux(4) = vn1*(e1+p1)

  else if(mn1<=-1.0) then !supersonic inflow
 
    !fluxes come from right state
    pflux(1:neqn) = 0 

  else !subsonic case

    q1=u1**2.0+v1**2.0
 
    pflux(1) = 0.25*rho1*c1*(mn1+1.0)**2.0
    pflux(2) = pflux(1)*(u1+(nx*(-vn1+2.0*c1))/g)
    pflux(3) = pflux(1)*(v1+(ny*(-vn1+2.0*c1))/g)
    pflux(4) = pflux(1)*(0.5*(q1-vn1**2.0)+(gm1*vn1+2.0*c1)**2.0/(2.0*g2m1))

  end if

  !NEGATIVE FLUXES => ui+1 state
  if(mn2>=1.0) then !supersonic outflow
   
    nflux(1:4) = 0

  else if(mn2<=-1.0) then !supersonic inflow
 
    !fluxes come from right state
    nflux(1) = vn2*rho2
    nflux(2) = vn2*rho2*u2+p2*nx
    nflux(3) = vn2*rho2*v2+p2*ny
    nflux(4) = vn2*(e2+p2)

  else !subsonic case

    q2=u2**2.0+v2**2.0
 
    nflux(1) = -0.25*rho2*c2*(mn2-1.0)**2.0
    nflux(2) = nflux(1)*(u2+(nx*(-vn2-2.0*c2))/g)
    nflux(3) = nflux(1)*(v2+(ny*(-vn2-2.0*c2))/g)
    nflux(4) = nflux(1)*(0.5*(q2-vn2**2.0)+(gm1*vn2-2.0*c2)**2.0/(2.0*g2m1))
 
  end if

  !get normal fluxes
  do ieqn=1,neqn
    flux(ieqn) = (pflux(ieqn)+nflux(ieqn))*del(3,iface)
  end do

  !define face fluxes => fi+1/2=f+(ui)+f-(ui+1)
  do ieqn=1,neqn
    rhsel(1,ieqn,icell)=rhsel(1,ieqn,icell)+flux(ieqn)
    rhsel(1,ieqn,jcell)=rhsel(1,ieqn,jcell)-flux(ieqn)
  end do

end do

end subroutine

!========================================================================

subroutine getghostv(unkno,intfac,nface,nelem,uchar,alpha,del,   &
                     fcent,ecent,neqn,scheme,g)

implicit none

real,dimension(:,:,:) :: unkno
real rho,u,v,p,c,e,mn,g,gm1
real rhoinf,uinf,vinf,pinf
real nx,ny,del(:,:),uchar(:),alpha
real fcent(:,:),ecent(:,:)
!real rhog,ug,vg,pg,cg,vng,vni,vninf,vtg,ci,R1,R2,cinf

integer intfac(:,:),nelem,nface,neqn
integer iface,bcflag,ighost,iel,scheme

gm1=g-1.0

!freestream values
uinf=uchar(1)*cos(alpha)
vinf=uchar(1)*sin(alpha)
rhoinf=uchar(2)
pinf=rhoinf*uchar(1)**2.0/(g*uchar(3)**2.0)

do iface=1,nface
  ighost=intfac(1,iface)
  bcflag=intfac(5,iface)

  iel=intfac(2,iface)

  if(scheme==1) then
    call ctop_conversion(unkno,u,v,rho,p,c,e,g,iel)
  else if(scheme==2 .or. scheme==3) then
    call getleftright(unkno,fcent,ecent,neqn,nelem,nface,  &
                        u,v,rho,p,c,e,g,iel,iface)
    !write(16,*) 'BC element:',iel,u,v,rho,p
  end if

  nx=del(1,iface)/del(3,iface)
  ny=del(2,iface)/del(3,iface)

  mn = (u*nx+v*ny)/c

  select case(bcflag)
  
  case(2) !solid surface
    unkno(1,1,ighost)=rho
    unkno(1,2,ighost)=rho*(u-2.0*(u*nx+v*ny)*nx)
    unkno(1,3,ighost)=rho*(v-2.0*(u*nx+v*ny)*ny)
    unkno(1,4,ighost)=e
    !unkno(1,1,ighost)=rhoinf
    !unkno(1,2,ighost)=rhoinf*uinf
    !unkno(1,3,ighost)=rhoinf*vinf
    !unkno(1,4,ighost)=pinf/gm1+0.5*rhoinf*uinf**2.0+0.5*rhoinf*vinf**2.0

    !write(35,*) 'solid boundary ghost cell:',ighost,iel

  case(4) !inflow/outflow
    !supersonic inflow
    if(mn<=-1.0) then
      unkno(1,1,ighost)=rhoinf
      unkno(1,2,ighost)=rhoinf*uinf
      unkno(1,3,ighost)=rhoinf*vinf
      unkno(1,4,ighost)=pinf/gm1+0.5*rhoinf*uinf**2.0+0.5*rhoinf*vinf**2.0
    !supersonic outflow
    else if(mn>=1.0) then
      unkno(1,1,ighost)=rho
      unkno(1,2,ighost)=rho*u
      unkno(1,3,ighost)=rho*v
      unkno(1,4,ighost)=p/gm1+0.5*rho*u**2.0+0.5*rho*v**2.0
    !subsonic inflow 
    else if(mn>-1.0 .and. mn<0) then
      !!free stream
      !vninf = uinf*nx+vinf*ny
      !cinf  = sqrt(g*pinf/rhoinf)

      !!interior cell
      !vni   = u*nx+v*ny
      !ci    = sqrt(g*p/rho)

      !!compute riemann invariants
      !R1   = vni+2.0*ci/gm1
      !R2   = vninf+2.0*cinf/gm1 
      !vtg  = u*(-ny)+v*(nx)
      !vng  = (R1+R2)/2.0
      !cg   = gm1/2.0*(R2-vng)

      !!compute ghost cell primitive var value
      !rhog = (g/cg*p/rho**g)**(1/(1.0-g))
      !pg   = cg**2.0*rhog/g        
      !ug   = vng*nx+vtg*(-ny)
      !vg   = vng*ny+vtg*(nx)

      !!compute ghost cell state
      !unkno(1,1,ighost)=rhog
      !unkno(1,2,ighost)=rhog*ug
      !unkno(1,3,ighost)=rhog*vg
      !unkno(1,4,ighost)=pg/gm1+0.5*rhog*ug**2.0+0.5*rhog*vg**2.0

      !write(*,*) 'for iel =',iel,rhog,'cg =',cg,R2,vng

      unkno(1,1,ighost)=rhoinf
      unkno(1,2,ighost)=rhoinf*uinf
      unkno(1,3,ighost)=rhoinf*vinf
      unkno(1,4,ighost)=pinf/gm1+0.5*rhoinf*uinf**2.0+0.5*rhoinf*vinf**2.0
 
    !subsonic outflow
    else if(mn>=0 .and. mn<1.0) then
      unkno(1,1,ighost)=rhoinf
      unkno(1,2,ighost)=rhoinf*uinf
      unkno(1,3,ighost)=rhoinf*vinf
      unkno(1,4,ighost)=pinf/gm1+0.5*rhoinf*uinf**2.0+0.5*rhoinf*vinf**2.0
    else
     write(*,*) 'ERROR - BOUNDARY CONDTION FLAG NOT RECOGNIZED => mn =',mn
     write(*,*) 'iface:',iface
     write(*,*) unkno(1,1:4,iel)
 
    end if

  case default
     write(*,*) 'ERROR - BOUNDARY CONDTION FLAG NOT RECOGNIZED'

  end select
  
  
end do

end subroutine

!===============CONSERVATIVE TO PRIMITIVE VARIABLE CONVERSION============

subroutine ctop_conversion(unkno,u,v,rho,p,c,e,g,icell)
implicit none

integer,intent(in) :: icell

real,intent(in)  :: unkno(:,:,:),g
real,intent(out) :: u,v,rho,p,c,e

  rho = unkno(1,1,icell)
  u   = unkno(1,2,icell)/rho
  v   = unkno(1,3,icell)/rho
  p   = (g-1.0)*(unkno(1,4,icell)-0.5*rho*u**2.0-0.5*rho*v**2.0)
  e   = unkno(1,4,icell)
  c   = sqrt(g*p/rho)
 
end subroutine

!========================================================================

subroutine setinit(unkno,intfac,nface,numfac,uchar,alpha,g)
implicit none

real unkno(:,:,:)
real,intent(in) :: uchar(:),alpha,g
real u,v,rho,p 

integer iface,numfac,iel,ier,intfac(:,:),nface

u=uchar(1)*cos(alpha)
v=uchar(1)*sin(alpha)
rho=uchar(2)
p=rho*uchar(1)**2.0/(g*uchar(3)**2.0)

!loop over interior cells only
do iface=1,numfac
  iel=intfac(1,iface+nface)
  ier=intfac(2,iface+nface)
  
  unkno(1,1,iel)=rho
  unkno(1,2,iel)=rho*u
  unkno(1,3,iel)=rho*v
  unkno(1,4,iel)=p/(g-1.0)+0.5*rho*u**2.0+0.5*rho*v**2.0
  unkno(1,1,ier)=rho
  unkno(1,2,ier)=rho*u
  unkno(1,3,ier)=rho*v
  unkno(1,4,ier)=p/(g-1.0)+0.5*rho*u**2.0+0.5*rho*v**2.0
end do

end subroutine

!========================================================================

subroutine eltonode(area,unkno,esup1,esup2,npoin,neqn,output)
implicit none

real area(:),unkno(:,:,:),totalarea,output(:,:,:)
integer npoin,neqn,ipoin,jelem,esup1(:),esup2(:),istor

output(1,1:neqn,1:npoin)=0.0

do ipoin=1,npoin
  totalarea=0.0
  do istor=esup2(ipoin)+1,esup2(ipoin+1)
    jelem=esup1(istor)
    totalarea=totalarea+area(jelem)
    output(1,1:neqn,ipoin)=output(1,1:neqn,ipoin)+  &
                        unkno(1,1:neqn,jelem)*area(jelem)
  end do
  output(1,1:4,ipoin)=output(1,1:neqn,ipoin)/totalarea
end do
        
end subroutine

!========================================================================

subroutine gettstep(intfac,cfl,nface,numfac,unkno,area,del,nelem,dt,dtlocal,g)
implicit none

real    :: unkno(:,:,:),dtlocal(:),del(:,:),area(:),elwsp(nelem)
integer :: intfac(:,:),nface,numfac,nelem

real ui,vi,rhoi,pi,ci,ei,g,lambda
real uj,vj,rhoj,pj,cj,ej,nx,ny,dt,cfl

integer iface,iel,jel,ielem

elwsp(1:nelem)=0.0d0

do iface=1,numfac+nface
  
  nx=del(1,iface)/del(3,iface)
  ny=del(2,iface)/del(3,iface)

  if(iface<=nface) then
    iel=intfac(2,iface)
    jel=intfac(1,iface)

    call ctop_conversion(unkno,ui,vi,rhoi,pi,ci,ei,g,iel)
    
    lambda=del(3,iface)*abs(nx*ui+ny*vi+ci)
 
    elwsp(iel)=elwsp(iel)+lambda
  else
    iel=intfac(1,iface)
    jel=intfac(2,iface)
    call ctop_conversion(unkno,ui,vi,rhoi,pi,ci,ei,g,iel)
    call ctop_conversion(unkno,uj,vj,rhoj,pj,cj,ej,g,jel)
    
    lambda=del(3,iface)*(0.5)*abs(nx*ui+nx*uj+ny*vi+ny*vj+ci+cj)
    
    elwsp(iel)=elwsp(iel)+lambda
    elwsp(jel)=elwsp(jel)+lambda
  end if 

end do

do ielem=1,nelem
  dtlocal(ielem)=cfl*area(ielem)/elwsp(ielem)
end do

dt=minval(dtlocal)

end subroutine

!========================================================================

subroutine getL2(unkno,unold,area,nelem,L2)
implicit none

real,dimension(:,:,:) :: unkno,unold

real area(:),L2
integer ielem,nelem

L2=0.0

do ielem=1,nelem
  L2=L2+(unkno(1,1,ielem)-unold(1,1,ielem))**2.0*area(ielem)
end do

L2=sqrt(L2)

end subroutine

!========================================================================

subroutine getgrads(gradlhs,numfac,nface,nelem,intfac,unkno,  &
                    neqn,ndegr,ecent)
implicit none

real,dimension(:,:),intent(in)   :: ecent
real,dimension(:,:,:),intent(in) :: gradlhs
integer,intent(in) :: numfac,nface,neqn,ndegr,nelem,intfac(:,:)

real,dimension(neqn) :: ui,uj,du
real gradrhs(2,neqn,nelem),unkno(:,:,:),xi,yi,xj,yj,dx,dy,wj
integer icell,iface,jcell,i,j,ielem

!initialize to zero
gradrhs(:,:,:)=0.0
unkno(2:ndegr,1:neqn,1:nelem+nface)=0.0

!Least squares method
do iface=nface+1,numfac+nface
  icell=intfac(1,iface)
  jcell=intfac(2,iface)

  xi=ecent(1,icell)
  yi=ecent(2,icell)
  xj=ecent(1,jcell)
  yj=ecent(2,jcell)

  dx=xj-xi
  dy=yj-yi

  ui(1:neqn)=unkno(1,1:neqn,icell)
  uj(1:neqn)=unkno(1,1:neqn,jcell)

  du(1:neqn)=uj(1:neqn)-ui(1:neqn)

  !weighting factor => wj
  wj=1.0/sqrt(dx**2.0+dy**2.0)
 
  gradrhs(1,1:neqn,icell)=gradrhs(1,1:neqn,icell)+wj**2.0*dx*du(1:neqn)
  gradrhs(2,1:neqn,icell)=gradrhs(2,1:neqn,icell)+wj**2.0*dy*du(1:neqn)
  gradrhs(1,1:neqn,jcell)=gradrhs(1,1:neqn,jcell)+wj**2.0*dx*du(1:neqn)
  gradrhs(2,1:neqn,jcell)=gradrhs(2,1:neqn,jcell)+wj**2.0*dy*du(1:neqn)
end do 

do ielem=1,nelem
  do i=1,2
    unkno(2,1:neqn,ielem)=unkno(2,1:neqn,ielem)+gradrhs(i,1:neqn,ielem)*  &
                                                gradlhs(1,i,ielem)
    unkno(3,1:neqn,ielem)=unkno(3,1:neqn,ielem)+gradrhs(i,1:neqn,ielem)*  &
                                                gradlhs(2,i,ielem)
  end do
end do

end subroutine

!========================================================================

subroutine getleftright(unkno,fcent,ecent,neqn,nelem,nface,  &
                        u,v,rho,p,c,e,g,icell,iface)
implicit none

real,intent(in)  :: unkno(:,:,:),fcent(:,:),ecent(:,:),g
real,intent(out) :: u,v,rho,p,e,c

real dx,dy,tmp(1,neqn,1)
integer ieqn,neqn,icell,iface,nelem,nface

dx=fcent(1,iface)-ecent(1,icell) 
dy=fcent(2,iface)-ecent(2,icell)
!write(16,*) "element:",icell,"face:",iface,"dx:",dx,"dy:",dy

do ieqn=1,neqn
  tmp(1,ieqn,1)=unkno(1,ieqn,icell)+dx*unkno(2,ieqn,icell)   &
                                   +dy*unkno(3,ieqn,icell)
end do
  
  call ctop_conversion(tmp,u,v,rho,p,c,e,g,1)

end subroutine

!========================================================================

subroutine geterror(unkno,nelem,uchar,area,g,error)

real,intent(in)  :: unkno(:,:,:),uchar(:),area(:),g
real,intent(out) :: error

real u,v,rho,p,c,e
real pinf,rhoinf,diff
integer nelem,ielem

error=0.0

rhoinf=uchar(2)
pinf=rhoinf*uchar(1)**2.0/(g*uchar(3)**2.0)

do ielem=1,nelem
  call ctop_conversion(unkno,u,v,rho,p,c,e,g,ielem)
   
  diff=p/pinf*(rhoinf/rho)**g-1.0  
  error = error + diff**2.0*area(ielem)
end do

error=sqrt(error)

end subroutine

!========================================================================

subroutine getmeshspace(area,del,nface,numfac,nelem,intfac,meshspace)

real,intent(in)  :: area(:),del(:,:)
real,intent(out) :: meshspace

real perimeter(nelem)
integer nelem,ielem,iface,nface,numfac,intfac(:,:),iel,jel

perimeter(:)=0.0

do iface=1,nface+numfac
  iel=intfac(1,iface)
  jel=intfac(2,iface)

  if(iface<=nface) then
    perimeter(jel)=perimeter(jel)+del(3,jel)
  else
    perimeter(iel)=perimeter(iel)+del(3,iel)
    perimeter(jel)=perimeter(jel)+del(3,jel)
  end if
end do

meshspace=0.0

do ielem=1,nelem
  meshspace = meshspace+area(ielem)/perimeter(nelem)
end do

meshspace=meshspace/nelem

end subroutine 
  
!=======================================================================

subroutine getgg_grad(unkno,area,nelem,nface,numfac,intfac,del,neqn,ndegr,g)

real,intent(in)      :: del(:,:),area(:),g
real,dimension(neqn) :: ui,uj,dudx,dudy
integer,intent(in)   :: intfac(:,:),nface,nelem,neqn,ndegr,numfac

real unkno(:,:,:)
integer iface,iel,jel,ielem,idegr

unkno(2:ndegr,:,:)=0.0

do iface=1,numfac+nface
  if(iface<=nface) then
    iel=intfac(2,iface)
    jel=intfac(1,iface)
  else
    iel=intfac(1,iface)
    jel=intfac(2,iface)
  end if

  ui(1:neqn)=unkno(1,1:neqn,iel)
  uj(1:neqn)=unkno(1,1:neqn,jel)

  dudx(1:neqn)=0.5*(ui(1:neqn)+uj(1:neqn))*del(1,iface)
  dudy(1:neqn)=0.5*(ui(1:neqn)+uj(1:neqn))*del(2,iface)

  !du/dx derivatives
  unkno(2,1:neqn,iel)=unkno(2,1:neqn,iel)+dudx(1:neqn)
  unkno(2,1:neqn,jel)=unkno(2,1:neqn,jel)-dudx(1:neqn)

  !du/dy derivatives
  unkno(3,1:neqn,iel)=unkno(3,1:neqn,iel)+dudy(1:neqn)
  unkno(3,1:neqn,jel)=unkno(3,1:neqn,jel)-dudy(1:neqn)

end do

!divide by area
do ielem=1,nelem
  do idegr=2,ndegr
    unkno(idegr,1:neqn,ielem)=unkno(idegr,1:neqn,ielem)/area(ielem)
  end do
end do

unkno(2:ndegr,1:neqn,nelem+1:nelem+nface)=0.0

end subroutine

end module flowsolver
