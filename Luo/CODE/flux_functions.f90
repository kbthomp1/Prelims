module flux_functions

  use kinddefs,  only : dp
  use gridtools, only : gridtype, get_global_dof

  implicit none
  private

  public :: get_mass_matrix_inv
  public :: add_lift_primal_domain
  public :: add_lift_second_domain
  public :: add_jump_second_face
  public :: add_flux_primal_face
  public :: add_boundary_flux

  real(dp), parameter :: my_4th = 0.25_dp

  real(dp), parameter :: zero = 0.0_dp
  real(dp), parameter :: half = 0.5_dp
  real(dp), parameter :: one  = 1.0_dp
  real(dp), parameter :: two  = 2.0_dp

contains

!===================== GET_MASS_MATRIX ==================================80
! Get the inverse of the mass matrix, \int(B_i*B_j*dV)
!============================================================================80
  function get_mass_matrix_inv(D) result(m_inv)
    real(dp), intent(in) :: D  ! 2*(tri area)
    real(dp), dimension(3,3) :: m_inv
    real(dp), parameter :: my_6  = 6._dp
    real(dp), parameter :: my_18 = 18._dp
    integer :: i
  continue
    m_inv = -my_6
    do i = 1,3
      m_inv(i,i) = my_18
    end do
    m_inv = m_inv/D
  end function get_mass_matrix_inv

!===================== ADD_LIFT_PRIMAL_DOMAIN ===============================80
! Add local lift functions to primal domain integral
!============================================================================80
  subroutine add_lift_primal_domain(icell,jcell,iface,grid,lhspo)

    integer,        intent(in) :: icell, jcell, iface
    type(gridtype), intent(in) :: grid
    real(dp), dimension(:,:), intent(inout) :: lhspo

    real(dp), dimension(3) :: bx, by
    integer,  dimension(2) :: cells
    integer  :: ilift, ip, i, cell, l
    real(dp) :: cell_area

  continue

    ilift = global_lift_coord(iface,grid)
    cells(1) = icell
    cells(2) = jcell

    do l = 1,2
      cell = cells(l)
      call get_basis(bx,by,grid,cell)
      cell_area = half*grid%geoel(5,cell)
      do i=1,grid%nnode
        ip = grid%inpoel(i,cell)
        lhspo(ip,ilift) = lhspo(ip,ilift) - (bx(i) + by(i))*cell_area
      end do
    end do

  end subroutine add_lift_primal_domain

!===================== ADD_LIFT_SECOND_DOMAIN ===============================80
! Add local lift to lift equation domain integral
!============================================================================80
  subroutine add_lift_second_domain(icell,jcell,iface,grid,lhspo)

    integer,        intent(in) :: icell, jcell, iface
    type(gridtype), intent(in) :: grid

    real(dp), dimension(:,:), intent(inout) :: lhspo

    integer :: ilift

    real(dp) :: icell_area, jcell_area

  continue

    ilift = global_lift_coord(iface,grid)
    icell_area = half*grid%geoel(5,icell)
    jcell_area = half*grid%geoel(5,jcell)
    lhspo(ilift,ilift) = lhspo(ilift,ilift) + two*(icell_area+jcell_area)

  end subroutine add_lift_second_domain

!===================== ADD_LIFT_SECOND_DOMAIN ===============================80
! Add jump condition to lift equation
!============================================================================80
  subroutine add_jump_second_face(icell,jcell,iface,grid,lhspo)

    integer,        intent(in) :: icell, jcell, iface
    type(gridtype), intent(in) :: grid

    real(dp), dimension(:,:), intent(inout) :: lhspo

    integer :: ilift, ip, jp, i

    real(dp) :: nx, ny, face_area, integrated_const
    real(dp) :: eta = 4._dp

  continue

    ilift = global_lift_coord(iface,grid)
    call get_face_normal(nx,ny,face_area,grid,iface)

    ! result of integral of jump function at iface
    integrated_const = my_4th*eta*(nx+ny)*face_area

    do i=1,grid%nnode
      ip = grid%inpoel(i,icell)
      jp = grid%inpoel(i,jcell)
      lhspo(ilift,jp) = lhspo(ilift,jp) + integrated_const
      lhspo(ilift,ip) = lhspo(ilift,ip) - integrated_const
    end do

  end subroutine add_jump_second_face

!======================= ADD_FLUX_PRIMAL_FACE ===============================80
! Add flux to primal equation at iface
!============================================================================80
  subroutine add_flux_primal_face(icell,jcell,iface,grid,phi,residual)

    integer,        intent(in) :: icell, jcell, iface
    type(gridtype), intent(in) :: grid

    real(dp), dimension(:), intent(in)    :: phi
    real(dp), dimension(:), intent(inout) :: residual

    real(dp), dimension(3) :: bx, by, grad

    integer, dimension(2) :: cells
    integer :: ip, i, l
    integer :: cell, inode, idof, jdof

    real(dp) :: nx, ny, face_area, flux, phi_local

  continue

    flux = zero

    cells(1) = icell
    cells(2) = jcell
    call get_face_normal(nx,ny,face_area,grid,iface)

    left_right_average: do l = 1,2

      cell = cells(l) !index of cell to form contibutions
      ! Basis functions -> dB/dx, dB/dy
      call get_basis(bx,by,grid,cell)

      ! half of averaged gradient
      grad(1:3) = half*(bx(1:3)*nx + by(1:3)*ny)

      ! integrated coefficients of phi
      do i = 1,grid%nnode
        ip = grid%inpoel(i,icell)
        phi_local = phi(get_global_dof(ip,icell,grid))
        flux = flux + phi_local*grad(i)*half*face_area
      end do

    end do left_right_average

    ! integrated coefficient of local lift
    !lift_coef = my_4th*face_area

    loop_face_nodes: do i = 3,4
      inode = grid%intfac(i,iface)
      idof = get_global_dof(inode,icell,grid)
      jdof = get_global_dof(inode,jcell,grid)
      residual(idof) = residual(idof) + flux  
      residual(jdof) = residual(jdof) - flux  !jcell has negative normal
    end do loop_face_nodes

  end subroutine add_flux_primal_face

!======================= ADD_BOUNDARY_FLUX ==================================80
! Add flux from Neumann BC
!============================================================================80
  subroutine add_boundary_flux(residual,grid,uinf,vinf)

    type(gridtype),         intent(in)    :: grid
    real(dp),               intent(in)    :: uinf,vinf
    real(dp), dimension(:), intent(inout) :: residual

    real(dp) :: cface, nx, ny, face_area
    integer :: iface, ip1, ip2, icell

  continue

    do iface=1,grid%nface
      if(grid%intfac(5,iface) == 4) then
       
        ! map the boundary nodes to the solution node ordering
        icell = grid%intfac(2,iface) ! interior cell
        ip1 = get_global_dof(grid%intfac(3,iface),icell,grid)
        ip2 = get_global_dof(grid%intfac(4,iface),icell,grid)

        call get_face_normal(nx,ny,face_area,grid,iface)
        cface = 0.5_dp*(uinf*nx + vinf*ny)*face_area
    
        residual(ip1) = residual(ip1) - cface
        residual(ip2) = residual(ip2) - cface
    
      end if
    end do

  end subroutine add_boundary_flux

  function global_lift_coord(iface,grid) result(ilift)
    integer,        intent(in)  :: iface
    type(gridtype), intent(in)  :: grid
    integer :: ilift
  continue
    ilift = grid%npoin + iface - grid%nface
  end function global_lift_coord

  subroutine get_basis(bx,by,grid,icell)
    integer,                intent(in)  :: icell
    type(gridtype),         intent(in)  :: grid
    real(dp), dimension(3), intent(out) :: bx, by
  continue
    bx(1:2) = grid%geoel(1:2,icell)
    bx(3)   = -(bx(1) - bx(2))
    by(1:2) = grid%geoel(3:4,icell)
    by(3)   = -(by(1) - by(2))
  end subroutine get_basis

  subroutine get_face_normal(nx,ny,face_area,grid,iface)
    integer,        intent(in) :: iface
    type(gridtype), intent(in) :: grid
    real(dp),      intent(out) :: nx, ny, face_area
  continue
    !compute unit normal vector
    nx = grid%del(1,iface)/grid%del(3,iface)
    ny = grid%del(2,iface)/grid%del(3,iface)
    face_area = grid%del(3,iface)
  end subroutine get_face_normal

end module flux_functions
