Module output

use flowsolver

implicit none

contains

subroutine tecoutput(area,unkno,esup1,esup2,npoin,neqn,iteration,nface,nelem,  &
                     ndegr,inpoel,coord,bcface,tec_dataname,surface_dataname,g)
implicit none

real,intent(in)    :: unkno(:,:,:),area(:),coord(:,:),g
integer,intent(in) :: npoin,ndegr,neqn,nface,nelem,iteration
integer,intent(in) :: inpoel(:,:),bcface(:,:),esup1(:),esup2(:)

real,dimension(npoin) :: upoin,vpoin,ppoin,rhopoin,Vtpoin
real u,v,rho,p,c,e,nodeunkno(ndegr,neqn,nelem+nface)

integer ipoin,ip1,ip2,i,j

character(len=10) iter_str
character(len=*)  tec_dataname,surface_dataname

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

!write iteration # to string
write(iter_str,'(i6)') iteration
call rmblank(iter_str)

open(21,file=trim(tec_dataname)//".dat",position='append',status='old')

write(21,*) 'ZONE T="',trim(iter_str)//'iterations','" NODES=',npoin,   &
             ",ELEMENTS=",nelem,",DATAPACKING=POINT,ZONETYPE=FETRIANGLE"
do i=1,npoin
  write(21,*) coord(1,i),coord(2,i),upoin(i),vpoin(i),rhopoin(i),ppoin(i)
end do

do j=1,nelem
  write(21,*) (inpoel(i,j),i=1,3)
end do

close(21)

!==============================================================================
!                          TECPLOT MAIN OUTPUT
!==============================================================================

open(22,file=trim(surface_dataname)//".dat",position='append',status='old')

write(22,*) 'ZONE T=',trim(iter_str)

do i=1,nface
  if(bcface(3,i)==2) then
    ip1=bcface(1,i)
    ip2=bcface(2,i)
    
    write(22,*) coord(1,ip1),coord(2,ip1),Vtpoin(ip1),rhopoin(ip1),ppoin(ip1)
    write(22,*) coord(1,ip2),coord(2,ip2),Vtpoin(ip2),rhopoin(ip2),ppoin(ip2)
  end if
end do

close(22)

end subroutine

!========================================================================

subroutine rmblank(string)

character(len=*) string

string=adjustl(string)
string=trim(string)

end subroutine

end module output



