program vmc_ho
    use randomseed
    use functions
    implicit none
    real(8), allocatable :: E(:), E_array(:,:), x(:), dphidalpha(:), E_dphidalpha(:)
    real(8) :: alpha, gamma
    integer :: i, steps, eq_steps, N, electrons, j

    call init_random_seed()
    call parameters
    do i = 1, 30
        call initialize_walkers         ! initialize walkers
        call metropolis                 ! the MC steps
        print '(a,1F6.3)',       'alpha:  ', alpha
        print '(a,1F9.6,1F8.5)', 'energy: ', avg(E(eq_steps:steps)), var(E(eq_steps:steps))
        alpha = alpha - 2*gamma*( avg(E_dphidalpha(eq_steps:steps)) - avg(E(eq_steps:steps)) * avg(dphidalpha(eq_steps:steps)) )
    end do

    contains

!----------------------------------------------------------------------------!      

subroutine parameters   
    alpha =     1.0_8           ! initial guess
    steps =     30000           ! tot # of steps
    eq_steps =  4000            ! # equilibration steps
    N =         400             ! # of walkers
    gamma =     0.1             ! damping constant
!     open(14, file='vmc.params')
!     read(14, *) steps
!     read(14, *) eq_steps
!     read(14, *) N
!     close(14)
    allocate(x(N), E(steps), E_array(1,eq_steps:steps), dphidalpha(steps), E_dphidalpha(steps))

end subroutine

!----------------------------------------------------------------------------!    
    
!----------------------------------------------------------------------------!  
subroutine initialize_walkers
    call random_number(x)
    x = x - 0.5_8
end subroutine
!----------------------------------------------------------------------------!  

!----------------------------------------------------------------------------!  

subroutine metropolis
    real(8), dimension(N) :: x_new, rnd_uni, rnd_uni2, transition_prob_ratio
    real(8) :: delta
    integer :: step, i, counter, mask(N)
    delta = 1._8
    counter = 0
    
    do step = 1, steps

        call random_number(rnd_uni)
        call random_number(rnd_uni2)
        rnd_uni = rnd_uni - 0.5_8
        x_new = x + rnd_uni*delta

        transition_prob_ratio = (transition_probability(x_new)/transition_probability(x))**2
        where (rnd_uni2 < transition_prob_ratio)
            mask = 1
        elsewhere
            mask = 0
        end where
        
        x = merge(x_new,x,mask==1)
        
        if (step >= eq_steps) call E_local(x, E(step), dphidalpha(step), E_dphidalpha(step))   
        
        ! adapt the of the width of the normal distribution such that 50% is accepted
        counter = counter + sum(mask)
        if (mod(step,100) == 0) then
            delta = delta*counter/(50.0d0*N)
            counter = 0
        end if
    end do
end subroutine

!----------------------------------------------------------------------------! 

!----------------------------------------------------------------------------!  

function transition_probability(x)
    real(8), dimension(N) :: transition_probability
    real(8), dimension(N), intent(in) :: x
    transition_probability =  exp(-alpha*(x**2))
end function

!----------------------------------------------------------------------------!  

!----------------------------------------------------------------------------!  

subroutine E_local(x, E, dphidalpha, E_dphidalpha)
    real(8), dimension(N), intent(in) :: x
    real(8), dimension(N) :: dphidalpha_array, E_tot
    real(8) :: E, dphidalpha, E_dphidalpha

    E_tot = alpha + x**2*(0.5_8-2._8*alpha**2)
    E = sum(E_tot)/N
    dphidalpha_array = -x**2
    dphidalpha = sum(dphidalpha_array)/N
    E_dphidalpha = sum(E_tot*dphidalpha_array)/N
end subroutine
!----------------------------------------------------------------------------! 

end program vmc_ho
