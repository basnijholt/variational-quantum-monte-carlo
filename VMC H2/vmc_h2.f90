program vmc_h2
    use randomseed
    use functions
    implicit none
    real(8), allocatable :: E(:), x(:, :, :), dphidbeta(:), E_dphidbeta(:), s_vector(:), E_array(:,:), beta_array(:,:), E_all(:,:), E_all_var(:,:)
    real(8) :: beta, alpha, gamma, s, s_min, s_max, a
    integer :: steps, eq_steps, N, N_s, electrons, i, j, minimization_steps

    call parameters
    call init_random_seed()
    
!$OMP PARALLEL DO PRIVATE(i, j, x, a, s, beta, E, dphidbeta, E_dphidbeta) num_threads(2)
    do i = 1, N_s 
    s = s_vector(i)
    beta = 0.4
    call cusp(a, s)
        do j = 1, minimization_steps
            call initialize_walkers(x)                                  ! initialize walkers
            call metropolis(x, a, s, beta, E, dphidbeta, E_dphidbeta)   ! the MC steps
            call new_beta(i, j, beta, E, dphidbeta, E_dphidbeta)        ! calculate new beta for the minimization
            call save_and_print(i, j, E, beta)
        end do 
    end do
!$OMP END PARALLEL DO

    call print_all(1)
    
  contains

!----------------------------------------------------------------------------!
subroutine initialize_walkers(x)
    real(8), intent(inout) :: x(:,:,:)
    call random_number(x)
    x = x - 0.5_8
end subroutine
!----------------------------------------------------------------------------!  

!----------------------------------------------------------------------------!  
subroutine new_beta(i, j, beta, E, dphidbeta, E_dphidbeta)
    real(8), intent(inout) :: beta
    real(8), intent(in) :: E_dphidbeta(:), E(:), dphidbeta(:)
    integer, intent(in) :: i, j
    print '(a,1I2,1I2)',      'i,j:    ', i,j
    print *,         'beta:   ', beta
    print *,         'energy: ', avg(E(eq_steps:steps))
    print *,         ' '
    beta = beta - gamma*2*( avg(E_dphidbeta(eq_steps:steps)) - avg(E(eq_steps:steps)) * avg(dphidbeta(eq_steps:steps)) )
end subroutine
!----------------------------------------------------------------------------!  

!----------------------------------------------------------------------------!
subroutine parameters
    electrons = 2                 ! # of electrons in Helium
    alpha = 2._8                  ! constant
    gamma = 1.0_8                 ! damping factor

    steps = 100000                 ! # of MC steps
    eq_steps = 10000               ! # of equilibrium steps
    N = 1000                       ! # of walkers
    N_s = 25                      ! # of proton-proton distances
    s_min = 1.375_8!.7_8                  ! minimum proton-proton distance
    s_max = 1.43_8!4.5_8                 ! maximum proton-proton distance
    minimization_steps = 25       ! # of steps in which beta is optimized for a minimal energy 
!     open(14, file='vmc.params')
!     read(14, *) steps
!     read(14, *) eq_steps
!     read(14, *) N
!     read(14, *) N_s
!     read(14, *) s_min
!     read(14, *) s_max
!     read(14, *) minimization_steps
!     close(14)
    allocate(x(3, N, electrons), E(steps), dphidbeta(steps), E_dphidbeta(steps), E_array(N_s,steps), s_vector(N_s), beta_array(N_s, minimization_steps), E_all(N_s, minimization_steps), E_all_var(N_s, minimization_steps))
        
    call linspace(s_min, s_max, N_s, s_vector)
end subroutine
!----------------------------------------------------------------------------!

!----------------------------------------------------------------------------!
subroutine metropolis(x, a, s, beta, E, dphidbeta, E_dphidbeta)
    real(8), intent(out) :: E(:), dphidbeta(:), E_dphidbeta(:)
    real(8), intent(inout) :: x(:,:,:)
    real(8), intent(in) :: a, s, beta
    real(8), dimension(3,N,electrons) :: x_new, rnd_uni
    real(8) :: rnd_uni2(N), delta, transition_prob_ratio(N)
    integer :: step, counter
    integer :: mask(N), mask_array(3,N,electrons)
    delta = 1._8
    counter = 0
    do step = 1, steps
        call random_number(rnd_uni)
        call random_number(rnd_uni2)
        rnd_uni = rnd_uni - 0.5_8
        
        x_new = x + rnd_uni*delta
        
        transition_prob_ratio = (transition_probability(x_new, a, s, beta)/transition_probability(x, a, s, beta))**2
        
        where (rnd_uni2 < transition_prob_ratio)
            mask = 1
        elsewhere
            mask = 0
        end where
 
        mask_array = spread(spread(mask,DIM=1,NCOPIES=3),DIM=3,NCOPIES=2)
        x = merge(x_new,x,mask_array==1)
        counter = counter + sum(mask)

        if (step >= eq_steps) call E_local(x, a, s, beta, E(step), dphidbeta(step), E_dphidbeta(step))
        
        if (mod(step,100) == 0) then
            delta = delta*counter/(50.0d0*N)
            counter = 0
        end if
    end do
    
end subroutine
!----------------------------------------------------------------------------!  

!----------------------------------------------------------------------------!  
subroutine r(x, a, s, r1, r1_L, r1_L_abs, r1_R, r1_R_abs, r2, r2_L, r2_L_abs, r2_R, r2_R_abs, r12, r12_abs, phi_1_L, phi_1_R, phi_1, phi_2_L, phi_2_R, phi_2)
    real(8), intent(in) :: x(:,:,:), a, s
    real(8), dimension(3,N), intent(out) :: r1, r1_L, r1_R, r2, r2_L, r2_R, r12
    real(8), dimension(N), intent(out) :: phi_1, phi_1_L, phi_1_R, phi_2, phi_2_L, phi_2_R, r1_L_abs, r1_R_abs, r2_L_abs, r2_R_abs, r12_abs
    real(8) :: nuclei(3,N), i_hat(3,N)
    i_hat = 0._8
    i_hat(1,:) = 1._8
    nuclei = 0.5_8*s*i_hat

    r1 = x(:,:,1)
    r2 = x(:,:,2)
    r1_L = r1 + nuclei
    r1_R = r1 - nuclei
    r2_L = r2 + nuclei
    r2_R = r2 - nuclei
    r12 = r1 - r2

    r1_L_abs = NORM2(r1_L,dim=1)
    r1_R_abs = NORM2(r1_R,dim=1)
    r2_L_abs = NORM2(r2_L,dim=1)
    r2_R_abs = NORM2(r2_R,dim=1)
    r12_abs = NORM2(r12,dim=1)

    phi_1_L = exp(-r1_L_abs/a)
    phi_1_R = exp(-r1_R_abs/a)
    phi_1 = phi_1_L + phi_1_R
    
    phi_2_L = exp(-r2_L_abs/a)
    phi_2_R = exp(-r2_R_abs/a)
    phi_2 = phi_2_L + phi_2_R

end subroutine
!----------------------------------------------------------------------------!

!----------------------------------------------------------------------------!
subroutine E_local(x, a, s, beta, E, dphidbeta, E_dphidbeta)
    real(8), intent(in) :: x(:,:,:), a, s, beta
    real(8), intent(out) :: dphidbeta, E_dphidbeta, E
    real(8), dimension(3,N) :: r1, r1_L, r1_R, r2, r2_L, r2_R, r12
    real(8), dimension(N) :: phi_1, phi_1_L, phi_1_R, phi_2, phi_2_L, phi_2_R, r1_L_abs, r1_R_abs, r2_L_abs, r2_R_abs, r12_abs
    real(8), dimension(3,N) :: r1_L_hat, r1_R_hat, r2_L_hat, r2_R_hat, r12_hat
    real(8), dimension(N) ::  E1, E2, E3, E4, E5_1, E5_2, E5_3, E5, E6, E_tot, dphidbeta_array
    
    call r(x, a, s, r1, r1_L, r1_L_abs, r1_R, r1_R_abs, r2, r2_L, r2_L_abs, r2_R, r2_R_abs, r12, r12_abs, phi_1_L, phi_1_R, phi_1, phi_2_L, phi_2_R, phi_2)

    r1_L_hat = r1_L/spread(r1_L_abs,1,3)
    r1_R_hat = r1_R/spread(r1_R_abs,1,3)
    r2_L_hat = r2_L/spread(r2_L_abs,1,3)
    r2_R_hat = r2_R/spread(r2_R_abs,1,3)
    r12_hat = r12/spread(r12_abs,1,3)

    E1 = -1._8/a**2
    E2 = (phi_1_L/r1_L_abs+phi_1_R/r1_R_abs)/(a*phi_1)
    E3 = (phi_2_L/r2_L_abs+phi_2_R/r2_R_abs)/(a*phi_2)
    E4 = 1._8/r12_abs - 1._8/r1_L_abs - 1._8/r1_R_abs - 1._8/r2_L_abs - 1._8/r2_R_abs

    E5_1 = (phi_1_L*sum(r1_L_hat*r12_hat,dim=1) + phi_1_R*sum(r1_R_hat*r12_hat,dim=1)) / phi_1
    E5_2 = (phi_2_L*sum(r2_L_hat*r12_hat,dim=1) + phi_2_R*sum(r2_R_hat*r12_hat,dim=1)) / phi_2
    E5_3 = 1._8/(2._8*a*(1._8+beta*r12_abs)**2)

    E5 = (E5_1-E5_2)*E5_3
    E6 = - ((4._8*beta+1._8)*r12_abs+4.0_8)/(4.0_8*(1+beta*r12_abs)**4._8*r12_abs)
    
    E_tot = E1 + E2 + E3 + E4 + E5 + E6 + 1._8/s
    E = sum(E_tot)/N
    
    dphidbeta_array = -r12_abs**2/(alpha*(1+beta*r12_abs)**2)
    dphidbeta = sum(dphidbeta_array)/N
    E_dphidbeta = sum(E_tot*dphidbeta_array)/N
end subroutine
!----------------------------------------------------------------------------! 

!----------------------------------------------------------------------------!  
function transition_probability(x, a, s, beta)
    real(8), dimension(N) :: transition_probability
    real(8), intent(in) :: x(:,:,:), a, s, beta
    real(8), dimension(3,N) :: r1, r1_L, r1_R, r2, r2_L, r2_R, r12
    real(8), dimension(N) :: phi_1, phi_1_L, phi_1_R, phi_2, phi_2_L, phi_2_R, r1_L_abs, r1_R_abs, r2_L_abs, r2_R_abs, r12_abs

    call r(x, a, s, r1, r1_L, r1_L_abs, r1_R, r1_R_abs, r2, r2_L, r2_L_abs, r2_R, r2_R_abs, r12, r12_abs, phi_1_L, phi_1_R, phi_1, phi_2_L, phi_2_R, phi_2)
    transition_probability = phi_1*phi_2*exp( r12_abs/(alpha*(1._8+beta*r12_abs)) )
    
end function
!----------------------------------------------------------------------------!

!----------------------------------------------------------------------------!  
subroutine print_all(switch)       
    integer :: switch
    character(30) :: fmtstring
    
    if (switch == 1) then
      write (fmtstring, '(a,I2,a)') '(',N_s,'F13.9)'
      open(11, file='E.dat')
      write (11, fmtstring) E_array(:,eq_steps:steps)
      open(12, file='data_int.dat')
      write(12, '(1I8)') N, steps, eq_steps, N_s, minimization_steps
      open(13, file='data_float.dat')
      write(13, '(1F10.6)') s_min, s_max
      open(14, file='s.dat')
      write(14, '(1F10.6)') s_vector
      open(15, file='beta.dat')
      write(15, fmtstring) beta_array
      open(16, file='E_all.dat')
      write(16, fmtstring) E_all
      open(17, file='E_all_var.dat')
      write(17, fmtstring) E_all_var      
    end if

end subroutine
!----------------------------------------------------------------------------! 

!----------------------------------------------------------------------------!
subroutine save_and_print(i, j, E, beta)
    real(8), intent(in) :: E(:), beta
    integer, intent(in) :: i, j
    if (j == minimization_steps) E_array(i,:) = E
    beta_array(i,j) = beta
    E_all(i,j) = avg(E(eq_steps:steps))
    E_all_var(i,j) = var(E(eq_steps:steps))
end subroutine
!----------------------------------------------------------------------------!
  
end program vmc_h2
