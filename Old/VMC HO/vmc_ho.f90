program vmc_ho
    use randomseed
    implicit none
    integer, parameter :: steps =     100000    ! tot # of steps
    integer, parameter :: eq_steps =  10000     ! # equilibration steps
    integer, parameter :: N =         1000      ! # of walkers
    integer, parameter :: N_alpha = 21          ! # of different alphas I want to calculate in the range [alpha_min,alpha_max]
    real(8), parameter :: alpha_min = 0.45_8    ! minimum alpha
    real(8), parameter :: alpha_max = 0.55_8    ! maximum alpha
    real(8) :: E(steps), E_array(N_alpha,steps) ! energy, array with energies of all alphas
    real(8) :: alpha, alphas(N_alpha)           ! alpha and the array with all the different alphas
    real(8) :: walker(N)                        ! positions of the walkers
    integer :: i
    
    call make_alphas
    
    do i = 1, N_alpha
        alpha = alphas(i)
        call init_random_seed()         ! random seed module
        call initialize_walkers         ! initialize walkers
        call metropolis                 ! the MC steps
        E_array(i,:) = E(:)/N
        call print_all(1)               ! print alpha, MC results and analytical result
    end do
    call print_all(2)                   ! print data to files
    
    contains

! Generate initial configuration using random positions for the walkers
!----------------------------------------------------------------------------!  
subroutine initialize_walkers
    call random_number(walker)
    walker = walker - 0.5_8
end subroutine
!----------------------------------------------------------------------------!  
    
! For every walker in the configuration propose a move from x to x_new
! the acceptance is transition_probability
!----------------------------------------------------------------------------!  
subroutine metropolis
    real(8) :: rnd_uni(N), rnd_uni2(N), delta, x, x_new
    integer :: counter, step, i
    E = 0._8
    delta = 1._8
    counter = 0
    do step = 1, steps
        ! generate an array with normal distributed numbers
        call random_number(rnd_uni)
        call random_number(rnd_uni2)
        rnd_uni2 = rnd_uni2 - .5_8
        
        ! Shift all walkers to a new position
        do i = 1,N
            x = walker(i)
            x_new = walker(i) + rnd_uni2(i) * delta
            
            ! move the walker if random value is smaller then the transition probability
            ! and count the number of walker moves that are accepted
            if (rnd_uni(i) < transition_probability(x, x_new)) then
                counter = counter + 1
                walker(i) = x_new
                x = x_new
            end if
            
            ! after equilibration the local energy is calculated
            if (step > eq_steps) then
                E(step) = E(step) + E_local(x)
            end if 
            
        end do
        
        ! adapt the of the width of the normal distribution such that 50% is accepted
        if (mod(step,100) == 0) then
            delta = delta*counter/(50._8*N)
            counter = 0
        end if
    end do
    
end subroutine
!----------------------------------------------------------------------------!  


!----------------------------------------------------------------------------!  
subroutine print_all(switch)       
    integer :: switch
    character(30) :: fmtstring
    
    if (switch == 1) then
      print *, 'alpha', alpha
      print *, 'MC result, energy, var', avg(E(eq_steps:steps))/N, var(E(eq_steps:steps))/sqrt(real(N))
      print *, 'analytical energy, var', 0.5*alpha+1/(8*alpha), (1-4*alpha**2)**2/(32*alpha**2)
      print *, ''
    end if
    
    if (switch == 2) then
      write (fmtstring, '(a,I2,a)') '(',N_alpha,'F10.5)'
      open(11, file='E.dat')
      write (11, fmtstring) E_array(:,eq_steps+1:steps)
      open(12, file='data_int.dat')
      write(12, *) N, steps, eq_steps, N_alpha
      open(13, file='data_float.dat')
      write(13, *) alpha_min, alpha_max
    end if
end subroutine
!----------------------------------------------------------------------------!  

!----------------------------------------------------------------------------!  
subroutine make_alphas
    do i = 1, N_alpha
        alphas(i) = (alpha_max - alpha_min) *real(i-1,8) / real(N_alpha-1,8) + alpha_min      
    end do
end subroutine
!----------------------------------------------------------------------------!  

!----------------------------------------------------------------------------!  
real(8) function transition_probability(x, x_new) result(transition_prob)
    real(8), intent(in) :: x, x_new
    transition_prob = exp(-2*alpha*(x_new**2-x**2))
end function
!----------------------------------------------------------------------------!  

!----------------------------------------------------------------------------!  
real(8) function E_local(x) result(E_loc)
    real(8), intent(in) :: x
    E_loc = alpha + x**2*(0.5_8-2._8*alpha**2)
end function
!----------------------------------------------------------------------------! 

!----------------------------------------------------------------------------!  
real(8) function rnd() result(rnd_number)
    call RANDOM_NUMBER(rnd_number)
end function
!----------------------------------------------------------------------------!  

!----------------------------------------------------------------------------!  
real(8) function rnd_gauss() result(rand_gauss)
    real(8), parameter :: pi = 4*atan(1._8)
    rand_gauss = sqrt(-2._8*log(rnd()))*sin(2._8*pi*rnd())
end function
!----------------------------------------------------------------------------!  

!----------------------------------------------------------------------------!  
real(8) function avg(x) result(average)
  real(8), intent(in) :: x(:)
  average = sum(x) / size(x)
end function avg
!----------------------------------------------------------------------------!  

!----------------------------------------------------------------------------!  
real(8) function var(x) result(variance)
  real(8), intent(in) :: x(:)
  variance = sum(x*x) / size(x) - avg(x)**2
end function var
!----------------------------------------------------------------------------!  

end program vmc_ho
