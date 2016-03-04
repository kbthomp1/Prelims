module io_helpers

  use kinddefs,  only : dp
  use gridtools, only : gridtype

  implicit none
  private

  public :: read_namelist
  public :: write_tec_volume
  public :: write_tec_surface
  public :: write_restart
  public :: read_tec_volume
  public :: read_restart

contains

!============================ READ_NAMELIST ==================================80
! Set default and override with user-specified namelist options
!=============================================================================80

  subroutine read_namelist

    use namelist_data, only : uinf, vinf, gridfile, nnode, uinf, vinf, nsteps, &
                              tec_dataname, lin_solver, tolerance, dt,         &
                              read_restart, restart_file, rk_order, cfl,       &
                              output_freq, read_tec_restart, tec_output_freq,  &
                              write_restart_freq, relative_tol, absolute_tol

    namelist /fe_input/ gridfile, nnode, uinf, vinf, nsteps, tec_dataname,     &
                        lin_solver, tolerance, dt, read_restart, restart_file, &
                        cfl, rk_order, output_freq, read_tec_restart,          &
                        tec_output_freq, write_restart_freq, relative_tol,     &
                        absolute_tol

  continue

    !Grid filename
    gridfile = "gridfile"

    !Number of verticies per element (triangles, quads, etc.)
    nnode = 3

    !Timestep to explicitly advance solution in pseudo-time
    dt = 1.E-6_dp

    !X-direction freestream velocity
    uinf = 1.0_dp

    !Y-direction freestream velocity
    vinf = 0.0_dp

    !Number of iterations
    nsteps = 1

    !Tecplot output filename
    tec_dataname = "output.dat"

    !Linear solver selection
    lin_solver = "gauss seidel"

    !Tolerance criterion of residual to exit linear solve
    tolerance = 1.E-20_dp

    !Tolerance criterion of relative residual to signal convergence
    relative_tol = 1.E-20_dp

    !Tolerance criterion of absolute residual to signal convergence
    absolute_tol = 1.E-20_dp

    !Flag to read a restart file
    read_restart = .false.

    !Flag to read a tecplot output file and interpolat to dofs
    read_tec_restart = .false.

    !Restart file in tecplot point format
    restart_file = "output.dat"

    !Order for runge-kutta time integration
    rk_order = 2

    !CFL number for determining timestep
    cfl = 1.0_dp

    !Frequency the restart file is written out
    write_restart_freq = 500

    !Output error norm frequency
    output_freq = 1

    !Output tecplot file frequency
    tec_output_freq = 0

    open(11,file="fe_input.nml",status="old")
    read(11,nml=fe_input)
    close(11)

  end subroutine read_namelist

!============================ WRITE_TEC_VOLUME ==============================80
! Write tecplot .dat files with solution in point format for volume
!============================================================================80
  subroutine write_tec_volume(tec_dataname,grid,phi,Vx,Vy,Vt)

    type(gridtype), intent(in) :: grid

    real(dp), dimension(:),   intent(in) :: Vx, Vy, Vt, phi

    character(len=*), intent(in) :: tec_dataname

    integer :: i,j

  continue

    open(21,file=tec_dataname,status='replace')
    
    write(21,*) 'TITLE = "',trim(tec_dataname),'"'
    write(21,*) 'VARIABLES = "X" "Y" "phi" "Vx" "Vy" "Vt"'
    write(21,*) 'ZONE NODES=',grid%npoin,",ELEMENTS=",grid%nelem,              &
                ",DATAPACKING=POINT,","ZONETYPE=FETRIANGLE"
    do i=1,grid%npoin
      write(21,*) grid%coord(1,i),grid%coord(2,i),phi(i),Vx(i),Vy(i),Vt(i)
    end do
    
    do j=1,grid%nelem
      write(21,*) (grid%inpoel(i,j),i=1,3)
    end do
    
    close(21)

  end subroutine write_tec_volume

!============================ WRITE_TEC_SURFACE =============================80
! Write tecplot .dat files with solution in point format for surface(s)
!============================================================================80
  subroutine write_tec_surface(grid,Vt)

    type(gridtype), intent(in) :: grid

    real(dp), dimension(:),   intent(in) :: Vt

    integer :: i, ip1, ip2

  continue

    open(15,file="channel_surface_lower.dat",status="replace")
    open(16,file="channel_surface_upper.dat",status="replace")
    
    write(15,*) 'TITLE = "Lower Surface"'
    write(15,*) 'VARIABLES = "X" "Y" "Vt"'
    write(16,*) 'TITLE = "Lower Surface"'
    write(16,*) 'VARIABLES = "X" "Y" "Vt"'

    do i=1,grid%nbcface
      if(grid%bcface(3,i)==2) then
        ip1=grid%bcface(1,i)
        ip2=grid%bcface(2,i)
        
        write(15,*) grid%coord(1,ip1),grid%coord(2,ip1),Vt(ip1)
        write(15,*) grid%coord(1,ip2),grid%coord(2,ip2),Vt(ip2)
      else 
        exit
      end if
    end do
    
    do i=1,grid%nbcface
      if(grid%bcface(3,i)==2) then
        ip1=grid%bcface(1,i)
        ip2=grid%bcface(2,i)
        
        write(16,*) grid%coord(1,ip1),grid%coord(2,ip1),Vt(ip1)
        write(16,*) grid%coord(1,ip2),grid%coord(2,ip2),Vt(ip2)
      end if
    end do

    close(15)
    close(16)

  end subroutine write_tec_surface

!============================ READ_TEC_VOLUME ===============================80
! read tecplot .dat files with solution in point format for volume
!============================================================================80
  subroutine read_tec_volume(tec_dataname,grid,phi,ndof)

    use flux_functions, only : get_global_dof

    type(gridtype),            intent(in)  :: grid
    integer,                   intent(in)  :: ndof
    real(dp), dimension(ndof), intent(out) :: phi

    real(dp), dimension(grid%npoin) :: nodal_phi

    character(len=*), intent(in) :: tec_dataname

    real(dp) :: dummy

    integer :: i, dof1, dof2, dof3
    integer :: ip1, ip2, ip3, ielem

  continue

    open(21,file=tec_dataname,status='old')
    do i = 1,3
      read(21,*)
    end do
    do i=1,grid%npoin
      read(21,*) dummy,dummy,nodal_phi(i)
    end do
    close(21)

    do ielem=1,grid%nelem
      ip1=grid%inpoel(1,ielem)
      ip2=grid%inpoel(2,ielem)
      ip3=grid%inpoel(3,ielem)

      dof1=get_global_dof(ip1,ielem,grid)
      dof2=get_global_dof(ip2,ielem,grid)
      dof3=get_global_dof(ip3,ielem,grid)

      phi(dof1) = nodal_phi(ip1)
      phi(dof2) = nodal_phi(ip2)
      phi(dof3) = nodal_phi(ip3)
    end do

  end subroutine read_tec_volume

  subroutine write_restart(phi,ndof)
    integer,                   intent(in) :: ndof
    real(dp), dimension(ndof), intent(in) :: phi
    integer :: i
  continue
    open(51,file="flow.restart", status="replace")
    do i = 1, ndof
      write(51,*) phi(i)
    end do
    close(51)
  end subroutine write_restart

  subroutine read_restart(phi,ndof)
    integer,                   intent(in)  :: ndof
    real(dp), dimension(ndof), intent(out) :: phi
    integer :: i
  continue
    open(52,file="flow.restart", status="old")
    do i = 1, ndof
      read(52,*) phi(i)
    end do
    close(52)
  end subroutine read_restart

end module io_helpers
