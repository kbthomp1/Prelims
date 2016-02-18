module io_helpers

  use kinddefs, only : dp

  implicit none
  private

  public :: write_tec_volume

contains

!============================ WRITE_TEC_OUTPUT ==============================80
! Write tecplot .dat files with solution in point format for volume
!============================================================================80
  subroutine write_tec_volume(tec_dataname,npoin,nelem,coord,inpoel,bcface,    &
                              phi,Vx,Vy,Vt)

    integer,                 intent(in) :: npoin, nelem
    integer, dimension(:,:), intent(in) :: inpoel, bcface

    real(dp), dimension(:),   intent(in) :: phi, Vx, Vy, Vt
    real(dp), dimension(:,:), intent(in) :: coord

    character(len=*), intent(in) :: tec_dataname

    integer :: i,j

  continue

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
    
    close(21)

  end subroutine write_tec_volume

  subroutine write_tec_surface()
  end subroutine write_tec_surface

end module io_helpers
