program vmc_hy
    use randomseed
    use functions
    implicit none
    real(8), allocatable :: E(:), E_array(:,:), x(:, :), dphidalpha(:), E_dphidalpha(:)
    real(8) :: alpha
    integer :: i, steps, eq_steps, N, electrons, j
 
    call init_random_seed()
    call parameters

    do i = 1, 5
        call initialize_walkers         ! initialize walkers
        call metropolis                 ! the MC steps
        print '(a,1F6.3)',       'alpha:  ', alpha
        print '(a,1F9.6,1F8.5)', 'energy: ', avg(E(eq_steps:steps)), var(E(eq_steps:steps))
        alpha = alpha - 2*( avg(E_dphidalpha(eq_steps:steps)) - avg(E(eq_steps:steps)) * avg(dphidalpha(eq_steps:steps)) )
    end do
    
    contains

!----------------------------------------------------------------------------!      
subroutine parameters   
    alpha = 8.40_8              ! initial guess
    steps =     30000           ! tot # of steps
    eq_steps =  4000            ! # equilibration steps
    N =         400             ! # of walkers
!     open(14, file='vmc.params')
!     read(14, *) steps
!     read(14, *) eq_steps
!     read(14, *) N
!     close(14)
    allocate(x(3, N), E(steps), E_array(1,eq_steps:steps), dphidalpha(steps), E_dphidalpha(steps))
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
    real(8), dimension(3,N) :: x_new, rnd_uni
    real(8) :: rnd_uni2(N), delta, transition_prob_ratio(N)
    integer :: step, i, counter
    integer :: mask(N), mask_array(3,N)
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

        mask_array = spread(mask,DIM=1,NCOPIES=3)
        x = merge(x_new,x,mask_array==1)
 
        
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
subroutine r(x, r1, r_abs)
    real(8), dimension(3,N), intent(in) :: x
    real(8), dimension(3,N) :: r1
    real(8), dimension(N) :: r_abs
    r1 = x(:,:)
    r_abs = NORM2(r1, dim=1)
end subroutine
!----------------------------------------------------------------------------!

!----------------------------------------------------------------------------!  
function transition_probability(x)
    real(8), dimension(N) :: transition_probability
    real(8), dimension(3,N), intent(in) :: x
    real(8), dimension(3,N) :: r1
    real(8), dimension(N) :: r_abs
    call r(x, r1, r_abs)
    transition_probability = exp(-alpha*(r_abs))
end function
!----------------------------------------------------------------------------!  

!----------------------------------------------------------------------------!  
subroutine E_local(x, E, dphidalpha, E_dphidalpha)
    real(8), dimension(3,N), intent(in) :: x
    real(8), dimension(3,N) :: r1
    real(8), dimension(N) :: r_abs
    real(8), dimension(N) :: dphidalpha_array, E_tot
    real(8) :: E, dphidalpha, E_dphidalpha
    call r(x, r1, r_abs)
    E_tot = -1/r_abs-0.5_8*alpha*(alpha-2._8/r_abs)
    E = sum(E_tot)/N
    dphidalpha_array = -r_abs
    dphidalpha = sum(dphidalpha_array)/N
    E_dphidalpha = sum(E_tot*dphidalpha_array)/N
end subroutine
!----------------------------------------------------------------------------! 

end program vmc_hy