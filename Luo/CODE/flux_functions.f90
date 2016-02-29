module flux_functions

  use kinddefs,  only : dp
  use gridtools, only : gridtype

  implicit none
  private

  public :: get_inv_mass_matrix
  public :: add_lift_primal_domain
  public :: add_lift_second_domain
  public :: add_jump_second_face
  public :: add_flux_primal_face

  real(dp), parameter :: my_4th = 0.25_dp

  real(dp), parameter :: zero = 0.0_dp
  real(dp), parameter :: half = 0.5_dp
  real(dp), parameter :: one  = 1.0_dp
  real(dp), parameter :: two  = 1.0_dp

contains

!===================== GET_INV_MASS_MATRIX ==================================80
! Get the inverse of the mass matrix, \int(B_i*B_j*dV)
!============================================================================80
  function get_inv_mass_matrix(D) result(M_inv)
    real(dp), intent(in) :: D  ! 2*(tri area)
    real(dp), dimension(3,3) :: M_inv
    real(dp), parameter :: six   = 6._dp
    real(dp), parameter :: my_18 = 18._dp
    integer :: i
  continue
    M_inv = -six
    do i = 1,3
      M_inv(i,i) = my_18
    end do
    M_inv = M_inv/D
  end function get_inv_mass_matrix

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
  subroutine add_flux_primal_face(icell,jcell,iface,grid,lhspo)

    integer,        intent(in) :: icell, jcell, iface
    type(gridtype), intent(in) :: grid

    real(dp), dimension(:,:), intent(inout) :: lhspo

    real(dp), dimension(3) :: bx, by, grad
    real(dp), dimension(3) :: phi_coefs
    real(dp), dimension(2) :: sgn

    integer, dimension(2) :: cells
    integer :: ilift, ip, jp, i, j, l, s
    integer :: cell, cell_idx

    real(dp) :: nx, ny, face_area, lift_coef

  continue

    ilift = global_lift_coord(iface,grid)
    cells(1) = icell
    cells(2) = jcell
    sgn(1)   = one   ! add flux to left cell
    sgn(2)   = -one  ! subtract flux from right cell
    call get_face_normal(nx,ny,face_area,grid,iface)

    flip_sign: do s = 1,2
      cell_idx = cells(s) !index of cell to add/subtract flux
      left_right_cells: do l = 1,2

        cell = cells(l) !index of cell to form contibutions
        ! Basis functions -> dB/dx, dB/dy
        call get_basis(bx,by,grid,cell)

        ! half of averaged gradient
        grad(1:3) = half*(bx(1:3)*nx + by(1:3)*ny)

        ! integrated coefficients of phi
        phi_coefs(1:3) = grad(1:3)*half*face_area

        ! integrated coefficient of local lift
        lift_coef = my_4th*face_area

        do i=1,grid%nnode
          ip = grid%inpoel(i,cell_idx)
          lhspo(ip,ilift) = lhspo(ip,ilift) + sgn(s)*lift_coef
          do j=1,grid%nnode
            jp = grid%inpoel(j,cell)
            lhspo(ip,jp) = lhspo(ip,jp) - sgn(s)*phi_coefs(j)
          end do
        end do

      end do left_right_cells
    end do flip_sign

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
