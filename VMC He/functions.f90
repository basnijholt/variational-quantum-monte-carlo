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

end module functions
