program vmc_hy
    use randomseed
    implicit none
    integer, parameter :: steps =     100000      ! tot # of steps
    integer, parameter :: eq_steps =  10000       ! # equilibration steps
    integer, parameter :: N =         1000        ! # of walkers
    integer, parameter :: N_alpha = 21            ! # of different alphas I want to calculate in the range [alpha_min,alpha_max]
    real(8), parameter :: alpha_min = 0.5_8       ! minimum alpha
    real(8), parameter :: alpha_max = 1.5_8       ! maximum alpha
    real(8) ::  E(steps), E_array(N_alpha,steps)  ! energy, array with energies of all alphas
    real(8) :: alpha, alphas(N_alpha)             ! alpha and the array with all the different alphas
    real(8) :: walker(3,N)                        ! positions of the walkers in 3-dimensions
    integer :: i
    
    call make_alphas
        
    do i = 1, N_alpha
        alpha = alphas(i)
        call init_random_seed()         ! random seed module
        call initialize_walkers         ! initialize walkers
        call metropolis                 ! the MC steps
        E_array(i,:) = E/N        
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
    real(8) :: rnd_uni(3,N), rnd_uni2(N), delta, r, r_new, x(3), x_new(3)
    integer :: counter, step, i, xyz
    E = 0._8
    counter = 0
    delta = 1._8
    
    do step = 1, steps
        ! generate an array with uniform random distributed numbers
        call random_number(rnd_uni)
        call random_number(rnd_uni2)
        rnd_uni = rnd_uni - 0.5_8
        ! Shift all walkers to a new position
        do i = 1, N
            r = 0._8
            r_new = 0._8
            
            do xyz = 1, 3
                x(xyz) = walker(xyz,i)
                x_new(xyz) = x(xyz) + rnd_uni(xyz,i)*delta
                r = r + x(xyz)**2
                r_new = r_new + x_new(xyz)**2
            end do
            
            r = sqrt(r)
            r_new = sqrt(r_new)
            
            ! move the walker if random value is smaller then the transition probability
            ! and count the number of walker moves that are accepted
            if (rnd_uni2(i) < transition_probability(r, r_new)) then
                counter = counter + 1
                walker(:,i) = x_new(:)
                r = r_new
            end if
            
            ! after equilibration the local energy is calculated
            if (step > eq_steps) then
                E(step) = E(step) + E_local(r)
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
      print *, 'MC result, energy, var', avg(E(eq_steps:steps))/N, var(E(eq_steps:steps)/N)
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
real(8) function transition_probability(r, r_new) result(transition_prob)
    real(8) :: r, r_new
    transition_prob = exp(-2*alpha*(r_new-r))
end function
!----------------------------------------------------------------------------!  

!----------------------------------------------------------------------------!  
real(8) function E_local(r) result(E_loc)
    real(8), intent(in) :: r
    E_loc = -1/r-0.5_8*alpha*(alpha-2._8/r)
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
    rand_gauss = sqrt(-2.0*log(rnd()))*sin(2.0_8*pi*rnd())
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

end program vmc_hy
