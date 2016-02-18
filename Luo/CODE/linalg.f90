module linalg

  use kinddefs, only : dp

  implicit none
  private

  public :: conjgrad
  public :: jacobi
  public :: gauss_seidel

contains

!================================= CONJGRAD ===============================80
! Solve linear system by conjugate gradient method
!==========================================================================80
  subroutine conjgrad(A,b,x,npoin,nsteps)
  
  real(dp)                   :: A(npoin,npoin),rsold,rsnew,alpha,pAp
  real(dp), dimension(npoin) :: b,x,ax,r,p,Ap
  
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
  
  subroutine jacobi(A,b,x,npoin,nsteps)
  
    real(dp)                   :: A(npoin,npoin),resnorm,dx,ad
    real(dp), dimension(npoin) :: b,x,ax,r
    
    integer                   :: npoin,nsteps,n,i,j
  
  continue
  
    x(1:npoin) = 0.0
    
    do n=1,nsteps
    
      resnorm     = 0.0
      dx          = 0.0
      ax(1:npoin) = 0.0
    
      do i=1,npoin
     
        do j=1,npoin   
          ax(i)=ax(i)+A(i,j)*x(j)
        end do
      
        r(i)    = b(i)-ax(i)
        resnorm = resnorm + r(i)*r(i)
    
        ad = 1/A(i,i)
    
        dx = r(i)*ad
    
        x(i)=x(i)+dx
      
      end do
        
      if(sqrt(resnorm)<1.0e-05) then
        write(*,*) "Solution has converged at,",n,"iterations"
        exit
      end if
    
      write(*,*)  n,sqrt(resnorm)
    
    end do
  
  end subroutine
  
!============================== GAUSS_SEIDEL ================================80
! Solve the linear system via the gauss seidel method
!============================================================================80
        SUBROUTINE gauss_seidel(A1,B1,X01,N,TOL)
  
        IMPLICIT NONE
  
        REAL(dp) :: A1(N,N),B1(N),X01(N)
        REAL(dp) :: X(N),AX(N),r(N)
        REAL(dp) :: NORM0,NORM,SUM1,SUM2
        REAL(dp) :: TOL
        INTEGER :: K=1,N,i,j
  
    11  DO I=1,N
          SUM1=0.0
          SUM2=0.0
            DO J=1,N
              IF (J.LT.I) SUM1=SUM1+A1(I,J)*X(J)
              IF (J.GT.I) SUM2=SUM2+A1(I,J)*X01(J)
            END DO
          X(I)=(B1(I)-SUM1-SUM2)/A1(I,I)
        END DO
  
        NORM    = 0.0 
        ax(1:N) = 0.0
  
        DO i=1,N
          DO j=1,N   
            ax(i)=ax(i)+A1(i,j)*X(j)
          end do
    
          r(i)    = B1(i)-ax(i)
          NORM = NORM + r(i)*r(i)
        END DO
  
        IF(K==1) NORM0=sqrt(NORM)
        
        K=K+1
  
        NORM=sqrt(NORM)/NORM0
  
        write(*,*) K,NORM
  
        IF (NORM.LT.TOL) then
          write(*,*) 'The solution has converged at, ',k,'iterations'
          GO TO 12
        END IF
  
        DO I=1,N
          X01(I)=X(I)
        END DO
  
        GO TO 11
  
     12 END SUBROUTINE

end module linalg
