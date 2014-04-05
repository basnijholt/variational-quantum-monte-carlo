program vmc_he
    use randomseed
    implicit none
    integer, parameter :: steps =     30000            ! tot # of steps
    integer, parameter :: eq_steps =  4000             ! # equilibration steps
    integer, parameter :: N =         300              ! # of walkers
    integer, parameter :: N_alpha = 2                   ! # of different alphas I want to calculate in the range [alpha_min,alpha_max]
    real(8), parameter :: alpha_min = 0.10_8            ! minimum alpha
    real(8), parameter :: alpha_max = 0.20_8            ! maximum alpha
    integer, parameter :: electrons = 2                 ! # of electrons in Helium
    real(8) :: E(steps), E_array(N_alpha,steps)    ! energy, array with energies of all alphas
    real(8) :: alpha, alphas(N_alpha)                   ! alpha and the array with all the different alphas
    real(8) :: walker(3, N, electrons)
    integer :: i
    
    call make_alphas
    
    do i = 1, N_alpha
        alpha = alphas(i)
        call init_random_seed()         ! random seed module
        call initialize_walkers         ! initialize walkers
        call metropolis                 ! the MC steps
        E_array(i,:) = E/N              ! 
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
    real(8), dimension(3,electrons) :: x, x_new
    real(8) :: rnd_uni(3,N,electrons), rnd_uni2(N), delta
    integer :: counter, step, i
    E = 0._8
    counter = 0
    delta = 1._8
    do step = 1, steps
        
        ! generate arrays with uniformly distributed numbers
        call random_number(rnd_uni)
        call random_number(rnd_uni2)
        rnd_uni = rnd_uni - 0.5_8

        ! Shift all walkers to a new position
        do i = 1, N             ! walkers
            x = walker(:,i,:)
            x_new = x + rnd_uni(:,i,:)*delta

            ! move the walker if random value is smaller then the transition probability
            ! and count the number of walker moves that are accepted
            if (rnd_uni2(i) < transition_probability(x, x_new)) then
                counter = counter + 1
                walker(:,i,:) = x_new
                x = x_new               
            end if
            
            ! after equilibration the local energy is calculated
            if (step > eq_steps) then
                E(step) = E(step) + E_local(x)
            end if 
        end do
        
        ! adapt the of the width of the normal distribution such that 50% is accepted
        if (mod(step,100) == 0) then
            delta = delta*counter/(50.0d0*N)
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
real(8) function E_local(x) result(E_loc)
    real(8), dimension(3,electrons), intent(in) :: x
    real(8) :: fac, fac2, r1(3), r2(3), r12(3), r1_abs, r2_abs, r12_abs
    r1 = x(:,1)
    r2 = x(:,2)
    r12 = r1-r2
    
    r1_abs = sqrt(dot_product(r1,r1))
    r2_abs = sqrt(dot_product(r2,r2))
    r12_abs = sqrt(dot_product(r12,r12))
    
    fac = 1.0d0/(1.0d0+alpha*r12_abs)
    fac2 = dot_product(r1/r1_abs-r2/r2_abs,r12)
    
    E_loc= -4.d0 + 1/r12_abs*(fac2*fac**2 -fac**3 + 1) - 0.25d0*fac**4
end function
!----------------------------------------------------------------------------! 

!----------------------------------------------------------------------------!  
real(8) function transition_probability(x, x_new) result(transition_prob_ratio)
    real(8), dimension(3,electrons), intent(in) :: x, x_new
    real(8), dimension(3) :: r1, r2, r12, r1_new, r2_new, r12_new
    real(8) :: r1_abs, r2_abs, r12_abs, transition_prob
    real(8) :: r1_abs_new, r2_abs_new, r12_abs_new, transition_prob_new
    
    r1 = x(:,1)
    r2 = x(:,2)
    r12 = r1-r2    

    r1_abs = sqrt(dot_product(r1,r1))
    r2_abs = sqrt(dot_product(r2,r2))
    r12_abs = sqrt(dot_product(r12,r12))
    transition_prob = exp(-2._8*(r1_abs+r2_abs))*exp(r12_abs*0.5_8/(1._8+alpha*r12_abs))

    r1_new = x_new(:,1)
    r2_new = x_new(:,2)
    r12_new = r1_new-r2_new   

    r1_abs_new = sqrt(dot_product(r1_new,r1_new))
    r2_abs_new = sqrt(dot_product(r2_new,r2_new))
    r12_abs_new = sqrt(dot_product(r12_new,r12_new))
    transition_prob_new = exp(-2._8*(r1_abs_new+r2_abs_new))*exp(r12_abs_new*0.5_8/(1._8+alpha*r12_abs_new))
  
    transition_prob_ratio = (transition_prob_new/transition_prob)**2
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

end program vmc_he
