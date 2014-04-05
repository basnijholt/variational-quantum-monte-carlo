module functions
  implicit none
  ! module parameters:
  integer, parameter :: maxiter = 200
  real(8), parameter :: tol = 1.d-14

  private
  public avg, var, solve, linspace, cusp

contains

!----------------------------------------------------------------------------!  
real(8) function avg(x)
    real(8), intent(in) :: x(:)
    avg = sum(x) / size(x)
end function avg
!----------------------------------------------------------------------------!  

!----------------------------------------------------------------------------!  
real(8) function var(x)
    real(8), intent(in) :: x(:)
    var = sum(x*x) / size(x) - avg(x)**2
end function var
!----------------------------------------------------------------------------!  

!----------------------------------------------------------------------------!  
subroutine linspace(x_min, x_max, x_N, lin_space)
    real(8), intent(in) :: x_min, x_max
    integer, intent(in) :: x_N
    real(8), intent(out) :: lin_space(:)
    integer :: i
    do i = 1, x_N
        lin_space(i) = (x_max - x_min) * real(i-1,8) / real(x_N-1,8) + x_min
    end do
end subroutine
!----------------------------------------------------------------------------!  

real(8) function f(a, s)
    real(8), intent(in) :: a, s
    f = 1.0_8/(1.0_8+exp(-s/a))-a
end function

real(8) function fp(a, s)
    real(8), intent(in) :: a, s
    fp = -s*exp(-s/a)/((1+exp(-s/a))**2*a**2)-1
end function

!----------------------------------------------------------------------------!  
subroutine cusp(a, s)
    real(8), intent(in) :: s
    real(8), intent(out) :: a
    real(8) :: deltax
    integer :: k
    a = 0.5 ! guess

    do k=1,maxiter        
        if (abs(f(a, s)) < tol) then
            exit  ! jump out of do loop
        endif
        deltax = f(a, s)/fp(a, s)
        a = a - deltax
    enddo
    if (k > maxiter) then
        ! might not have converged
        if (abs(f(a, s)) > tol) then
            print *, '*** Warning: has not yet converged'
        endif
    endif 

end subroutine
!----------------------------------------------------------------------------!  

!----------------------------------------------------------------------------!  
subroutine solve(f, fp, x0, x)

    ! Estimate the zero of f(x) using Newton's method. 
    ! Input:
    !   f:  the function to find a root of
    !   fp: function returning the derivative f'
    !   x0: the initial guess
    ! Returns:
    !   the estimate x satisfying f(x)=0 (assumes Newton converged!) 
     
    implicit none
    real(8), intent(in) :: x0
    real(8), external :: f, fp
    real(8), intent(out) :: x

    ! Declare any local variables:
    real(8) :: deltax, fx, fxprime
    integer :: k

    ! initial guess
    x = x0

    ! Newton iteration to find a zero of f(x) 

    do k=1,maxiter

        ! evaluate function and its derivative:
        fx = f(x)
        fxprime = fp(x)

        if (abs(fx) < tol) then
            exit  ! jump out of do loop
        endif

        ! compute Newton increment x:
        deltax = fx/fxprime

        ! update x:
        x = x - deltax

    enddo

    if (k > maxiter) then
        ! might not have converged

        fx = f(x)
        if (abs(fx) > tol) then
            print *, '*** Warning: has not yet converged'
        endif
    endif 


end subroutine solve
!----------------------------------------------------------------------------!  

end module functions
