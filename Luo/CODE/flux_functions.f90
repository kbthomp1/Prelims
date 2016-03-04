module flux_functions

  use kinddefs,  only : dp
  use gridtools, only : gridtype, get_global_dof

  implicit none
  private

  public :: get_mass_matrix_inv
  public :: add_lift_primal_domain
  public :: add_flux_primal_face
  public :: add_boundary_flux
  public :: compute_local_lift
  public :: get_global_dof
  public :: get_basis

  real(dp), parameter :: my_8th = 0.125_dp

  real(dp), parameter :: zero = 0.0_dp
  real(dp), parameter :: half = 0.5_dp
  real(dp), parameter :: one  = 1.0_dp
  real(dp), parameter :: two  = 2.0_dp

contains

!===================== COMPUTE_LOCAL_LIFT ===================================80
! Compute the local lifting function at interior faces
!============================================================================80
  subroutine compute_local_lift(phi,lift,grid,ndof)

    integer,        intent(in) :: ndof
    type(gridtype), intent(in) :: grid

    real(dp), dimension(ndof),        intent(in)  :: phi
    real(dp), dimension(grid%numfac), intent(out) :: lift

    integer :: iface, icell, jcell, i, fp1, fp2

    real(dp) :: eta, jump, iarea, jarea
    real(dp) :: nx, ny, face_area
    real(dp) :: u1_l, u2_l, u1_r, u2_r

  continue

    eta = real(grid%nnode+1 ,dp)
    !eta = 40._dp
    !eta = 0.001_dp
    interior_faces: do iface = grid%nbcface+1, grid%nface

      i = iface - grid%nbcface
      icell = grid%intfac(1,iface)  ! left cell
      jcell = grid%intfac(2,iface)  ! right cell

      call get_face_normal(nx,ny,face_area,grid,iface)

      iarea = half*grid%geoel(5,icell)
      jarea = half*grid%geoel(5,jcell)

      fp1 = grid%intfac(3,iface) ! point 1 on face
      fp2 = grid%intfac(4,iface) ! point 2 on face

      ! Get degrees on freedom on both sides of face
      u1_l = phi(get_global_dof(fp1,icell,grid))
      u2_l = phi(get_global_dof(fp2,icell,grid))
      u1_r = phi(get_global_dof(fp1,jcell,grid))
      u2_r = phi(get_global_dof(fp2,jcell,grid))

      jump = (u1_r + u2_r - u1_l - u2_l)*(nx+ny)
      !write(*,*) "CHECK: u",u1_l,u2_l, u1_r, u2_r

      lift(i) = -eta*my_8th*jump

    end do interior_faces

  end subroutine compute_local_lift

!========================= GET_MASS_MATRIX ==================================80
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
  subroutine add_lift_primal_domain(icell,jcell,iface,grid,lift,residual)

    integer,                intent(in) :: icell, jcell, iface
    type(gridtype),         intent(in) :: grid
    real(dp), dimension(:), intent(in) :: lift
    real(dp), dimension(:), intent(inout) :: residual

    real(dp), dimension(3) :: bx, by
    integer,  dimension(2) :: cells
    integer  :: ip, i, cell, l, k
    real(dp) :: face_area

  continue

    cells(1) = icell
    cells(2) = jcell
    k = iface - grid%nbcface

    do l = 1,2
      cell = cells(l)
      call get_basis(bx,by,grid,cell)
      face_area = grid%del(3,iface)
      do i=1,grid%nnode
        ip = get_global_dof(grid%inpoel(i,cell),cell,grid)
        residual(ip) = residual(ip) + (bx(i) + by(i))*lift(k)*face_area
      end do
    end do

  end subroutine add_lift_primal_domain

!======================= ADD_FLUX_PRIMAL_FACE ===============================80
! Add flux to primal equation at iface
!============================================================================80
  subroutine add_flux_primal_face(icell,jcell,iface,grid,phi,lift,residual)

    integer,        intent(in) :: icell, jcell, iface
    type(gridtype), intent(in) :: grid

    real(dp), dimension(:), intent(in)    :: phi
    real(dp), dimension(:), intent(in)    :: lift
    real(dp), dimension(:), intent(inout) :: residual

    real(dp), dimension(3) :: bx, by, grad

    integer, dimension(2) :: cells
    integer :: ip, i, l, k
    integer :: cell, inode, idof, jdof

    real(dp) :: nx, ny, face_area, flux, phi_local

  continue

    flux = zero

    cells(1) = icell
    cells(2) = jcell
    k = iface - grid%nbcface
    call get_face_normal(nx,ny,face_area,grid,iface)

    left_right_average: do l = 1,2

      cell = cells(l) !index of cell to form contibutions
      ! Basis functions -> dB/dx, dB/dy
      call get_basis(bx,by,grid,cell)

      ! half of averaged gradient
      grad(1:3) = half*(bx(1:3)*nx + by(1:3)*ny)

      ! integrated coefficients of phi
      do i = 1,grid%nnode
        ip = grid%inpoel(i,cell)
        phi_local = phi(get_global_dof(ip,cell,grid))
        flux = flux + phi_local*grad(i)*half*face_area
      end do

    end do left_right_average

    ! integrated local lifting operator from cell face
    flux = flux - lift(k)*half*(nx+ny)*face_area

    !write(*,*) "CHECK: face:",iface, icell,jcell
    !write(*,*) "CHECK: flux = ",flux

    loop_face_nodes: do i = 3,4
      inode = grid%intfac(i,iface)
      idof = get_global_dof(inode,icell,grid)
      jdof = get_global_dof(inode,jcell,grid)
      residual(idof) = residual(idof) - flux  
      residual(jdof) = residual(jdof) + flux  !jcell has negative normal
    end do loop_face_nodes

  end subroutine add_flux_primal_face

!======================= ADD_BOUNDARY_FLUX ==================================80
! Add flux from Neumann BC
!============================================================================80
  subroutine add_boundary_flux(residual,grid,uinf,vinf)

    type(gridtype),         intent(in)    :: grid
    real(dp),               intent(in)    :: uinf,vinf
    real(dp), dimension(:), intent(inout) :: residual

    real(dp) :: flux, nx, ny, face_area
    integer :: iface, ip1, ip2, icell

  continue

    ! Specify Neumann BC
    do iface=1,grid%nbcface
      if(grid%intfac(5,iface) == 4) then
       
        ! map the boundary nodes to the solution node ordering
        icell = grid%intfac(2,iface) ! interior cell
        ip1 = get_global_dof(grid%intfac(3,iface),icell,grid)
        ip2 = get_global_dof(grid%intfac(4,iface),icell,grid)

        call get_face_normal(nx,ny,face_area,grid,iface)
        flux = 0.5_dp*(uinf*nx + vinf*ny)*face_area
        !write(*,*) "CHECK: face: ",iface, icell
        !write(*,*) "CHECK: bc flux = ",flux
    
        residual(ip1) = residual(ip1) - flux
        residual(ip2) = residual(ip2) - flux
    
      end if
    end do

  end subroutine add_boundary_flux

  subroutine get_basis(bx,by,grid,icell)
    integer,                intent(in)  :: icell
    type(gridtype),         intent(in)  :: grid
    real(dp), dimension(3), intent(out) :: bx, by
  continue
    bx(1:2) = grid%geoel(1:2,icell)
    bx(3)   = -(bx(1) + bx(2))
    by(1:2) = grid%geoel(3:4,icell)
    by(3)   = -(by(1) + by(2))
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
