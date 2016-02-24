Module grid

implicit none

contains

!==========================GRID FILE READ======================================
subroutine readgrid_lou(ndimn,ntype,nelem,npoin,nface,inpoel,     &
                        gridfile,coord,bcface)

implicit none

integer      :: titleline,ndimn,ntype,nelem,npoin,nface
integer      :: i,j,elem_dummy,poin_dummy,face_dummy

real, dimension(:,:),allocatable  :: coord

integer, dimension(:,:), allocatable :: inpoel,bcface

character(len=100) :: gridfile

open(12,file=gridfile,status="old")
read(12,*) titleline

do i=1,titleline
  read(12,*)
end do

read(12,*) 
read(12,*) ndimn,ntype
read(12,*)
read(12,*) nelem,npoin,nface
read(12,*)

allocate(inpoel(ntype,nelem))
allocate(coord(ndimn,npoin))
allocate(bcface(ndimn+1,nface))

do j=1,nelem
  read(12,*) elem_dummy,(inpoel(i,j),i=1,ntype)
end do

read(12,*)

do j=1,npoin
  read(12,*) poin_dummy,(coord(i,j),i=1,ndimn)
end do

do j=1,npoin+2
  read(12,*)
end do

do j=1,nface
  read(12,*) face_dummy,(bcface(i,j),i=1,ndimn+1)
end do

close(12)

end subroutine

!=============================================================================

subroutine getesup(inpoel,npoin,nelem,nnode,esup1,esup2)
implicit none

integer,dimension(:),allocatable :: esup1,esup2

integer ipoin,npoin,ielem,nelem,inode,nnode,mesup,istor
integer inpoel(:,:)

mesup = npoin*6

allocate(esup2(npoin+1))
allocate(esup1(mesup))

!initialize esup2
esup2(1:npoin+1) = 0

!1st pass
do ielem=1,nelem
  do inode=1,nnode
    esup2(inpoel(inode,ielem)+1)=esup2(inpoel(inode,ielem)+1)+1
  end do
end do

do ipoin=2,npoin+1
  esup2(ipoin)=esup2(ipoin)+esup2(ipoin-1)
end do

!2nd pass
do ielem=1,nelem
  do inode=1,nnode
    ipoin=inpoel(inode,ielem)
    istor=esup2(ipoin)+1
    esup2(ipoin)=istor
    esup1(istor)=ielem
  end do
end do

do ipoin=npoin+1,2,-1
  esup2(ipoin)=esup2(ipoin-1)
end do

esup2(1)=0

end subroutine

!=============================================================================

subroutine getpsup(inpoel,esup1,esup2,npoin,nnode,psup1,psup2)

integer,dimension(:,:),intent(in) :: inpoel
integer,dimension(:),intent(in)   :: esup1,esup2
integer lpoin(npoin),psup1(npoin*6),psup2(npoin+1)
integer istor,npoin,jpoin,ipoin,ielem,inode,iesup,nnode

lpoin(1:npoin)=0
psup2(1)=0
istor=0

do ipoin=1,npoin
  lpoin(ipoin)=ipoin
  do iesup=esup2(ipoin)+1,esup2(ipoin+1)
    ielem=esup1(iesup)
    do inode=1,nnode
      jpoin=inpoel(inode,ielem)
      if(lpoin(jpoin) /= ipoin) then
        istor=istor+1
        psup1(istor)=jpoin
        lpoin(jpoin)=ipoin
      end if
    end do
  end do
  psup2(ipoin+1)=istor
end do

end subroutine

!===========================================================================

subroutine getesuel(inpoel,esup1,esup2,nfael,nelem,npoin,esuel)
implicit none

integer,dimension(:,:) :: inpoel
integer,dimension(:)   :: esup1,esup2

integer,dimension(:,:),allocatable :: esuel

integer nfael,nelem,lpoin(npoin)
integer ifael,jfael,ielem,jelem,lhelp(2),lhelp2(2),ipoin,jpoin,npoin
integer jnofa,istor,icoun

allocate(esuel(nfael,nelem))

lpoin(1:npoin)=0
esuel(1:nfael,1:nelem)=0

do ielem=1,nelem
  do ifael=1,nfael
    if(ifael<nfael) then
      lhelp(1:2)=inpoel(ifael:ifael+1,ielem)
    else
      lhelp(1)=inpoel(ifael,ielem)
      lhelp(2)=inpoel(1,ielem)
    end if
    lpoin(lhelp(1:2))=1
    ipoin=lhelp(1)

    do istor=esup2(ipoin)+1,esup2(ipoin+1)
      jelem=esup1(istor)
      if(jelem/=ielem) then
        do jfael=1,nfael
          if(jfael<nfael) then
            lhelp2(1:2)=inpoel(jfael:jfael+1,jelem)
          else
            lhelp2(1)=inpoel(jfael,jelem)
            lhelp2(2)=inpoel(1,jelem)
          end if

          icoun=0
          do jnofa=1,2
            jpoin=lhelp2(jnofa)
            if(lpoin(jpoin)==1) then
              icoun=icoun+1
            end if
          end do
          if(icoun==2) then
            esuel(ifael,ielem)=jelem
            esuel(jfael,jelem)=ielem
          end if
        end do
      end if
    end do
    lpoin(lhelp(1:2))=0
  end do
end do

end subroutine

!=====================================================================

subroutine getintfac(esuel,inpoel,nelem,nfael,numfac,nface,intfac)

integer,dimension(:,:) :: esuel,inpoel
integer,dimension(:,:),allocatable :: intfac

integer nelem,ielem,jelem
integer ifael,nfael,numfac,iface,intbcfac(4,nface),nface
integer ier,iel,ip1,ip2

numfac=0
iface=1

do ielem=1,nelem
  do ifael=1,nfael
    iel=ielem
    ier=esuel(ifael,ielem)
    if(iel<ier) numfac=numfac+1
  end do
end do

do ielem=1,nelem
  do ifael=1,nfael
    jelem=esuel(ifael,ielem)

    if(ifael<nfael) then
      ip1=inpoel(ifael,ielem)
      ip2=inpoel(ifael+1,ielem)
    else
      ip1=inpoel(ifael,ielem)
      ip2=inpoel(1,ielem)
    end if

    if(jelem==0) then
      intbcfac(1,iface)=nelem+iface
      intbcfac(2,iface)=ielem
      intbcfac(3,iface)=ip1
      intbcfac(4,iface)=ip2
      iface=iface+1
    end if
 
  end do
end do

allocate(intfac(5,numfac+nface))
intfac(1:4,1:nface)=intbcfac(1:4,1:nface)
intfac(5,1:(numfac+nface))=0

do ielem=1,nelem
  do ifael=1,nfael
    iel=ielem
    ier=esuel(ifael,ielem)
    if(ifael<nfael) then
      ip1=inpoel(ifael,ielem)
      ip2=inpoel(ifael+1,ielem)
    else
      ip1=inpoel(ifael,ielem)
      ip2=inpoel(1,ielem)
    end if
    if(iel<ier) then
      intfac(1,iface)=iel
      intfac(2,iface)=ier
      intfac(3,iface)=ip1
      intfac(4,iface)=ip2
      iface=iface+1
    end if
  end do
end do
     
end subroutine
     
!========================================================================

subroutine getnormal(intfac,coord,numfac,nface,del)

integer,dimension(:,:),intent(in) :: intfac
real,dimension(:,:),intent(in)    :: coord

real,dimension(:,:),allocatable   :: del

real    :: x1,x2,y1,y2
integer :: iface,numfac,nface,ip1,ip2,i

allocate(del(3,numfac+nface))

!write(*,*) '******** normal vectors *********'

do iface=1,numfac+nface
  ip1=intfac(3,iface)
  ip2=intfac(4,iface)
  
  x1=coord(1,ip1)
  x2=coord(1,ip2)
  
  y1=coord(2,ip1)
  y2=coord(2,ip2)

  del(1,iface)=y2-y1
  del(2,iface)=-(x2-x1)
  del(3,iface)=sqrt(del(1,iface)**2.0+del(2,iface)**2.0)

end do

end subroutine

!========================================================================

subroutine getbcflag(esup1,esup2,inpoel,bcface,nface,    &
                       npoin,nfael,intfac)
implicit none

integer,dimension(:,:) :: bcface,inpoel,intfac
integer,dimension(:)   :: esup1,esup2

integer nface,nfael,ifael,istor,npoin,ipoin,jelem,i
integer lhelp(2),lhelp2(nfael),lpoin(npoin),bcflag,iface,ip1,ip2

lpoin(1:npoin)=0

do iface=1,nface

  lhelp(1:2) = bcface(1:2,iface)
  bcflag     = bcface(3,iface)
  lpoin(lhelp(1:2))=1
  ipoin=lhelp(1)

  do istor=esup2(ipoin)+1,esup2(ipoin+1)
      jelem=esup1(istor)
      do ifael=1,nfael
        lhelp2(ifael)=inpoel(ifael,jelem)
      end do

      lpoin(lhelp2(1:nfael))=lpoin(lhelp2(1:nfael))+1

      if(lpoin(lhelp(2))==2) then
        ip1=lhelp(1)
        ip2=lhelp(2)
        do i=1,nface
          if(intfac(2,i)==jelem) then
            if(ip1==intfac(3,i) .and. ip2==intfac(4,i)) then
              intfac(5,i)=bcflag
              exit
            else if(ip1==intfac(4,i) .and. ip2==intfac(3,i)) then
              intfac(5,i)=bcflag
              exit
            end if
          end if
        end do
        lpoin(lhelp2(1:nfael))=lpoin(lhelp2(1:nfael))-1 
      else
        lpoin(lhelp2(1:nfael))=lpoin(lhelp2(1:nfael))-1 
      end if

    end do

  lpoin(lhelp(1:2))=0

end do
  
end subroutine 

!========================================================================

subroutine getcellarea(inpoel,coord,nelem,area)

real coord(:,:),area(:),x1,x2,x3,y1,y2,y3
integer inpoel(:,:),ip1,ip2,ip3,ielem,nelem

do ielem=1,nelem
  ip1=inpoel(1,ielem)
  ip2=inpoel(2,ielem)
  ip3=inpoel(3,ielem)

  x1=coord(1,ip1)
  x2=coord(1,ip2)
  x3=coord(1,ip3)
  
  y1=coord(2,ip1)
  y2=coord(2,ip2)
  y3=coord(2,ip3)
  
  area(ielem)=0.5*sqrt((x1*y3-x2*y3-x1*y1+x2*y1)**2.0+   &
                       (y1*x3-y2*x3-y1*x1+y2*x1)**2.0)
end do

end subroutine

!========================================================================

subroutine getcentroid(nface,numfac,nelem,intfac,inpoel,coord,fcent,ecent)

integer,dimension(:,:),intent(in)           :: inpoel,intfac
real,dimension(:,:),allocatable,intent(out) :: fcent,ecent

real,intent(in)    :: coord(:,:)
integer,intent(in) :: nface,numfac,nelem

integer iface,ielem,ip1,ip2,ip3
real x1,x2,x3,y1,y2,y3

allocate(fcent(2,nface+numfac))
allocate(ecent(2,nelem))

!compute face centroid
do iface=1,nface+numfac
  ip1=intfac(3,iface)
  ip2=intfac(4,iface)
  x1=coord(1,ip1)
  y1=coord(2,ip1)
  x2=coord(1,ip2)
  y2=coord(2,ip2)
  
  !compute x-centroid of iface
  fcent(1,iface)=(x1+x2)*0.5
  !compute y-centroid of iface
  fcent(2,iface)=(y1+y2)*0.5
end do

do ielem=1,nelem
  ip1=inpoel(1,ielem)
  ip2=inpoel(2,ielem)
  ip3=inpoel(3,ielem)
  x1=coord(1,ip1)
  y1=coord(2,ip1)
  x2=coord(1,ip2)
  y2=coord(2,ip2)
  x3=coord(1,ip3)
  y3=coord(2,ip3)

  !compute x-centroid of ielem
  ecent(1,ielem)=(x1+x2+x3)/3.0
  !compute y-centroid of ielem
  ecent(2,ielem)=(y1+y2+y3)/3.0
end do

end subroutine

!========================================================================

subroutine getgradlhs(ecent,fcent,intfac,nface,numfac,nelem,gradlhs)

real,intent(in)    :: ecent(:,:),fcent(:,:)
real,intent(out)   :: gradlhs(:,:,:)
integer,intent(in) :: intfac(:,:),nface,numfac,nelem

integer :: icell,jcell,iface,i,ielem
real    :: xi,xj,yi,yj,wj,dx,dy,tmp(2,2),D


write(*,*) 'computing gradlhs'
gradlhs(:,:,:) = 0.0

do iface=1,numfac+nface
  if(iface<=nface) then
    icell=intfac(2,iface)
    jcell=intfac(1,iface)

    dx=2.0*(fcent(1,iface)-ecent(1,icell))
    dy=2.0*(fcent(2,iface)-ecent(2,icell))
  else
    icell=intfac(1,iface)
    jcell=intfac(2,iface)

    xi=ecent(1,icell)
    yi=ecent(2,icell)
    xj=ecent(1,jcell)
    yj=ecent(2,jcell)
  
    dx=xj-xi
    dy=yj-yi
  end if
  
  !weighting factor => wj
  wj=1.0/sqrt(dx**2.0+dy**2.0)

  gradlhs(1,1,icell)=gradlhs(1,1,icell)+wj**2.0*dx**2.0
  gradlhs(1,2,icell)=gradlhs(1,2,icell)+wj**2.0*dx*dy
  gradlhs(2,1,icell)=gradlhs(2,1,icell)+wj**2.0*dx*dy
  gradlhs(2,2,icell)=gradlhs(2,2,icell)+wj**2.0*dy**2.0

  gradlhs(1,1,jcell)=gradlhs(1,1,jcell)+wj**2.0*dx**2.0
  gradlhs(1,2,jcell)=gradlhs(1,2,jcell)+wj**2.0*dx*dy
  gradlhs(2,1,jcell)=gradlhs(2,1,jcell)+wj**2.0*dx*dy
  gradlhs(2,2,jcell)=gradlhs(2,2,jcell)+wj**2.0*dy**2.0
end do


do ielem=1,nelem
  D = gradlhs(1,1,ielem)*gradlhs(2,2,ielem)   &
      -gradlhs(1,2,ielem)*gradlhs(2,1,ielem)
  
  tmp(1,1) = gradlhs(2,2,ielem)
  tmp(2,2) = gradlhs(1,1,ielem)
  tmp(1,2) =-gradlhs(1,2,ielem)
  tmp(2,1) =-gradlhs(2,1,ielem)

  tmp(:,:)=tmp(:,:)/D

  !NOTE: This is the matrix inverse of the lhs
  gradlhs(:,:,ielem)=tmp(:,:)

end do

write(*,*) 'completed computing gradlhs'

end subroutine

end module grid
