module linalg

  use kinddefs, only : dp

  implicit none
  private

  public :: conjgrad
  public :: jacobi
  public :: gauss_seidel

  real(dp), parameter :: zero = 0.0_dp
  real(dp), parameter :: one  = 1.0_dp

contains

!================================= CONJGRAD ===============================80
! Solve linear system by conjugate gradient method
!==========================================================================80
  subroutine conjgrad(A,b,x,npoin,nsteps)
  
  real(dp)                   :: rsold,rsnew,alpha,pAp
  real(dp), dimension(npoin) :: b,x,ax,r,p,Ap

  real(dp), dimension(npoin,npoin) :: A
  
  integer                   :: npoin,nsteps,n,i,j

  continue
  
  ax(1:npoin)  = 0.0
  rsold        = 0.0
  
  do i=1,npoin
  
    do j=1,npoin   
      ax(i)=ax(i)+A(i,j)*x(j)
    end do
  
    r(i)=b(i)-ax(i)
    p(1:npoin)=r(1:npoin)
    rsold=rsold+r(i)*r(i)
  end do
  
  do n=1,nsteps
  
    Ap(1:npoin) = 0.0
    pAp   = 0.0
    rsnew = 0.0      
  
    do i=1,npoin
  
      do j=1,npoin
        Ap(i)=Ap(i)+A(i,j)*p(j)
      end do
  
      do j=1,npoin
        pAp=pAp+p(j)*Ap(j)
      end do
   
      alpha=rsold/pAp
      x(i)=x(i)+alpha*p(i);
      r(i)=r(i)-alpha*Ap(i);
      rsnew=rsnew+r(i)*r(i)
    
      if(sqrt(rsnew)<1e-10) then
        write(*,*) "Solution has converged at,",n,",iterations"
        exit
      end if
  
      p(i)=r(i)+rsnew/rsold*p(i)
      rsold=rsnew
    end do
  end do
  
  end subroutine
  
!================================== JACOBI ==================================80
! Solve the linear system via the point jacobi method
!============================================================================80
  
  subroutine jacobi(A,b,x,ndof,nsteps,tolerance)

    integer, intent(in) :: ndof, nsteps

    real(dp), intent(in) :: tolerance

    real(dp), dimension(ndof), intent(in)  :: b
    real(dp), dimension(ndof), intent(out) :: x

    real(dp), dimension(ndof,ndof), intent(in) :: A
  
    real(dp), dimension(ndof) :: xnew

    real(dp) :: resnorm,ax,r

    integer :: n,i,j
  
  continue
  
    x = zero
    
    write(*,"(A,4x,A)") " Iteration","L2_norm"

    do n=1,nsteps
    
      resnorm = zero
   
      do i=1,ndof
        ax = zero
        do j=1,ndof   
          ax = ax + A(i,j)*x(j)
        end do
      
        ! Get residual
        r = b(i)-ax
        resnorm = resnorm + r**2

        r = b(i) - (ax-A(i,i)*x(i))
    
        xnew(i) = r/A(i,i)
      end do

      x = xnew
        
      write(*,"(i6,6x,e12.5)") n, sqrt(resnorm)
    
      if(sqrt(resnorm)<tolerance) then
        write(*,*) "Solution has converged at,",n,"iterations"
        exit
      end if
    
    end do
  
  end subroutine
  
!============================== GAUSS_SEIDEL ================================80
! Solve the linear system via the gauss seidel method
!============================================================================80
  subroutine gauss_seidel(A1,b1,x01,n,nsteps,tol)
  
    real(dp) :: A1(n,n),b1(n),x01(n)
    real(dp) :: x(n),Ax(n),r(n)
    real(dp) :: norm0,norm,sum1,sum2
    real(dp) :: tol

    integer  :: k,n,i,j,nsteps,istep

  continue

    k = 0

    write(*,"(A,4x,A)") " Iteration","L2_norm"

    outer_loop: do istep = 1,nsteps
        do i=1,n
          sum1=0.0
          sum2=0.0
            do j=1,n
              if (j.LT.i) sum1=sum1+A1(i,j)*x(j)
              if (j.GT.i) sum2=sum2+A1(i,j)*x01(j)
            end do
          x(i)=(b1(i)-sum1-sum2)/A1(i,i)
        end do
  
        norm    = 0.0 
        ax(1:n) = 0.0
  
        do i=1,n
          do j=1,n   
            ax(i)=ax(i)+A1(i,j)*x(j)
          end do
    
          r(i)    = b1(i)-ax(i)
          norm = norm + r(i)*r(i)
        end do
  
        if(k==0) norm0=sqrt(norm)
        
        k=k+1
  
        norm=sqrt(norm)/norm0
  
        write(*,"(i6,6x,e12.5)") k,norm
  
        if (norm.LT.tol) then
          write(*,*) 'The solution has converged at, ',k,'iterations'
          exit
        end if
  
        do i=1,n
          x01(i)=x(i)
        end do

      end do outer_loop
  
    end subroutine gauss_seidel

end module linalg
