module flux_functions

  use kinddefs,  only : dp
  use gridtools, only : gridtype

  implicit none
  private

  public :: get_mass_matrix
  public :: add_lift_primal_domain
  public :: add_lift_second_domain
  public :: add_jump_second_face
  public :: add_flux_primal_face

  real(dp), parameter :: my_4th = 0.25_dp

  real(dp), parameter :: zero = 0.0_dp
  real(dp), parameter :: half = 0.5_dp
  real(dp), parameter :: one  = 1.0_dp
  real(dp), parameter :: two  = 2.0_dp

contains

!===================== GET_MASS_MATRIX ==================================80
! Get the inverse of the mass matrix, \int(B_i*B_j*dV)
!============================================================================80
  function get_mass_matrix(D) result(M)
    real(dp), intent(in) :: D  ! 2*(tri area)
    real(dp), dimension(3,3) :: M
    real(dp), parameter :: my_24th = 1._dp/24._dp
    integer :: i
  continue
    M = D*my_24th
    do i = 1,3
      M(i,i) = two*D*my_24th
    end do
  end function get_mass_matrix

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
    integer :: ip, jp, i, l
    integer :: cell

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

      ! local nodal values

      ! integrated coefficients of phi
      do i = 1,grid%nnode
        phi_local = phi(grid%inpoel(i,cell))
        flux = flux + phi_local*grad(i)*half*face_area
      end do

    end do left_right_average

      ! integrated coefficient of local lift
      !lift_coef = my_4th*face_area

    do i=1,grid%nnode
      ip = grid%inpoel(i,icell)
      jp = grid%inpoel(i,jcell)
      !lhspo(ip,ilift) = lhspo(ip,ilift) + sgn(s)*lift_coef
      residual(ip) = residual(ip) + flux
      residual(jp) = residual(ip) - flux
    end do

  end subroutine add_flux_primal_face

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
