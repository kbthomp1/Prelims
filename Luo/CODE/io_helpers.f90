module io_helpers

  use kinddefs, only : dp

  implicit none
  private

  public :: read_namelist
  public :: write_tec_volume
  public :: write_tec_surface

contains

  subroutine read_namelist

    use namelist_data, only : uinf, vinf, gridfile, nnode, uinf, vinf, 
                              nsteps, tec_dataname, lin_solver

    namelist /fe_input/ gridfile, nnode, uinf, vinf, bc_case, nsteps,          &
                        tec_dataname, lin_solver

  continue

    open(11,file="fe_input.nml",status="old")
    read(11,nml=fe_input)
    close(11)

  end subroutine read_namelist

!============================ WRITE_TEC_VOLUME ==============================80
! Write tecplot .dat files with solution in point format for volume
!============================================================================80
  subroutine write_tec_volume(tec_dataname,npoin,nelem,coord,inpoel, phi,     &
                              Vx,Vy,Vt)

    integer,                 intent(in) :: npoin, nelem
    integer, dimension(:,:), intent(in) :: inpoel

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

!============================ WRITE_TEC_SURFACE =============================80
! Write tecplot .dat files with solution in point format for surface(s)
!============================================================================80
  subroutine write_tec_surface(npoin,bcface,coord,Vt)

    integer,                 intent(in) :: npoin
    integer, dimension(:,:), intent(in) :: bcface

    real(dp), dimension(:),   intent(in) :: Vt
    real(dp), dimension(:,:), intent(in) :: coord

    integer :: i, ip1, ip2

  continue

    open(15,file="channel_surface_lower.dat",status="replace")
    open(16,file="channel_surface_upper.dat",status="replace")
    
    write(15,*) 'TITLE = "Lower Surface"'
    write(15,*) 'VARIABLES = "X" "Y" "Vt"'
    write(16,*) 'TITLE = "Lower Surface"'
    write(16,*) 'VARIABLES = "X" "Y" "Vt"'

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

  end subroutine write_tec_surface

end module io_helpers
