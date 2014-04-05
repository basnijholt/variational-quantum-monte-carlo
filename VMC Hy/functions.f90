module functions
  implicit none
  private
  
  public avg, var

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
real(8) function rnd() result(rnd_number)
    call RANDOM_NUMBER(rnd_number)
end function
!----------------------------------------------------------------------------!  

!----------------------------------------------------------------------------!  
real(8) function rnd_gauss() result(rand_gauss)
    real(8), parameter :: pi = 4*atan(1._8)
    rand_gauss = sqrt(-2.0*log(rnd()))*sin(2.0_8*pi*rnd())
end function
!----------------------------------------------------------------------------!  

end module functions
