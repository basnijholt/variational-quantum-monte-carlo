program vmc_h2
    use randomseed
    implicit none
    real(8), allocatable :: E(:), E_array(:,:), s_array(:), walker(:, :, :)
    real(8) :: beta, s_min, s_max, s, alpha, a 
    integer :: i, steps, eq_steps, N, N_s, electrons
    call parameters    
    call init_random_seed()         ! random seed module
    call make_s
    
    do i = 1, N_s
        beta = 10
        call initialize_walkers         ! initialize walkers
        call metropolis                 ! the MC steps
        E_array(i,:) = E/N              ! 
        call print_all(1)               ! print beta, MC results and analytical result to screen
    end do
    call print_all(2)                   ! print data to files
    
    contains

!----------------------------------------------------------------------------!  
subroutine initialize_walkers
    call random_number(walker)
    walker = walker - 0.5_8
end subroutine
!----------------------------------------------------------------------------!  

!----------------------------------------------------------------------------!    
subroutine parameters   
    electrons = 2                 ! # of electrons in Helium
    s = 1.3988657845_8            ! distance of nuclei (74E-12)/(0.529E-10)
    a = 0.840751247975988_8       ! solution to 1/(1+exp(-1.3988657845/a))=a
    alpha = 2._8

    open(14, file='vmc.params')
    read(14, *) steps
    read(14, *) eq_steps
    read(14, *) N
    read(14, *) N_s
    read(14, *) s_min
    read(14, *) s_max
    close(14)
    
    allocate(walker(3, N, electrons), s_array(N_s), E(steps), E_array(N_s,steps))
end subroutine
!----------------------------------------------------------------------------!  

!----------------------------------------------------------------------------!  
subroutine metropolis
    real(8), dimension(3,electrons) :: x, x_new
    real(8) :: rnd_uni(3,N,electrons), rnd_uni2(N), delta, transition_prob_ratio
    integer :: counter, step, i
    E = 0._8
    counter = 0
    delta = 1._8
    do step = 1, steps
        call random_number(rnd_uni)
        call random_number(rnd_uni2)
        rnd_uni = rnd_uni - 0.5_8

        do i = 1, N             ! walkers
            x = walker(:,i,:)
            x_new = x + rnd_uni(:,i,:)*delta
            
            transition_prob_ratio = (transition_probability(x_new)/transition_probability(x))**2
            if (rnd_uni2(i) < transition_prob_ratio) then
                counter = counter + 1
                walker(:,i,:) = x_new
                x = x_new
            end if

            if (step > eq_steps) then
                E(step) = E(step) + E_local(x)
            end if 
        end do
        
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
      print *, 'beta', beta
      print *, 'MC result, energy, var', avg(E(eq_steps:steps))/N, var(E(eq_steps:steps)/N)
      print *, ''
    end if
    
    if (switch == 2) then
      write (fmtstring, '(a,I2,a)') '(',N_s,'F10.5)'
      open(11, file='E.dat')
      write (11, fmtstring) E_array(:,eq_steps+1:steps)
      open(12, file='data_int.dat')
      write(12, *) N, steps, eq_steps, N_s
      open(13, file='data_float.dat')
      write(13, *) s_min, s_max
    end if
end subroutine
!----------------------------------------------------------------------------! 

!----------------------------------------------------------------------------!  
subroutine make_s
    do i = 1, N_s
        s_array(i) = (s_max - s_min) * real(i-1,8) / real(N_s-1,8) + s_min
    end do
end subroutine
!----------------------------------------------------------------------------!  

!----------------------------------------------------------------------------!  
subroutine r(x, r1, r1_L, r1_L_abs, r1_R, r1_R_abs, r2, r2_L, r2_L_abs, r2_R, r2_R_abs, r12, r12_abs, phi_1_L, phi_1_R, phi_1, phi_2_L, phi_2_R, phi_2)
    real(8), dimension(3,electrons), intent(in) :: x
    real(8), dimension(3) :: r1, r1_L, r1_R, r2, r2_L, r2_R, r12
    real(8) :: phi_1, phi_1_L, phi_1_R, phi_2, phi_2_L, phi_2_R
    real(8) :: r1_L_abs, r1_R_abs, r2_L_abs, r2_R_abs, r12_abs
    real(8) :: nuclei(3), i_hat(3)

    i_hat(1) = 1._8
    i_hat(2) = 0._8
    i_hat(3) = 0._8
    nuclei = 0.5_8*s*i_hat
    
    r1 = x(:,1)
    r1_L = r1 + nuclei
    r1_L_abs = NORM2(r1_L)
    r1_R = r1 - nuclei
    r1_R_abs = NORM2(r1_R)
    r2 = x(:,2)
    r2_L = r2 + nuclei
    r2_L_abs = NORM2(r2_L)
    r2_R = r2 - nuclei
    r2_R_abs = NORM2(r2_R)
    r12 = r1 - r2
    r12_abs = NORM2(r12)
    
    phi_1_L = exp(-r1_L_abs/a)
    phi_1_R = exp(-r1_R_abs/a)
    phi_1 = phi_1_L + phi_1_R
    
    phi_2_L = exp(-r2_L_abs/a)
    phi_2_R = exp(-r2_R_abs/a)
    phi_2 = phi_2_L + phi_2_R
end subroutine
!----------------------------------------------------------------------------!

!----------------------------------------------------------------------------!
real(8) function E_local(x) result(E_loc)
    real(8), dimension(3,electrons), intent(in) :: x
    real(8), dimension(3) :: r1, r1_L, r1_R, r2, r2_L, r2_R, r12
    real(8) :: phi_1, phi_1_L, phi_1_R, phi_2, phi_2_L, phi_2_R, r1_L_abs, r1_R_abs, r2_L_abs, r2_R_abs, r12_abs
    real(8) :: dot_prod1(3), dot_prod2(3), E1, E2, E3, E4
    
    call r(x, r1, r1_L, r1_L_abs, r1_R, r1_R_abs, r2, r2_L, r2_L_abs, r2_R, r2_R_abs, r12, r12_abs, phi_1_L, phi_1_R, phi_1, phi_2_L, phi_2_R, phi_2)

    ! I split up the energy in multiple terms
    E1 = -1._8/a**2 + 1._8/(a*phi_1)*(phi_1_L/r1_L_abs+phi_1_R/r1_R_abs) + 1._8/(a*phi_2)*(phi_2_L/r2_L_abs+phi_2_R/r2_R_abs)
    E2 = -(1._8/r1_L_abs + 1._8/r1_R_abs + 1._8/r2_L_abs + 1._8/r2_R_abs) + 1._8/r12_abs
    
    dot_prod1 = (phi_1_L*r1_L/r1_L_abs + phi_1_R*r1_R/r1_R_abs)/phi_1 - (phi_2_L*r2_L/r2_L_abs + phi_2_R*r2_R/r2_R_abs)/phi_2 
    dot_prod2 = (r12/r12_abs)/(2._8*a*(1._8+beta*r12_abs)**2)
    E3 = dot_product(dot_prod1,dot_prod2)
    E4 = -((4._8*beta+1._8)*r12_abs+4._8)/(4._8*(1._8+beta*r12_abs)**4*r12_abs)
    E_loc = E1 + E2 + E3 + E4
end function
!----------------------------------------------------------------------------! 

!----------------------------------------------------------------------------!  
real(8) function transition_probability(x) result(transition_prob)
    real(8), dimension(3,electrons), intent(in) :: x
    real(8), dimension(3) :: r1, r1_L, r1_R, r2, r2_L, r2_R, r12
    real(8) :: phi_1, phi_1_L, phi_1_R, phi_2, phi_2_L, phi_2_R, r1_L_abs, r1_R_abs, r2_L_abs, r2_R_abs, r12_abs
    
    call r(x, r1, r1_L, r1_L_abs, r1_R, r1_R_abs, r2, r2_L, r2_L_abs, r2_R, r2_R_abs, r12, r12_abs, phi_1_L, phi_1_R, phi_1, phi_2_L, phi_2_R, phi_2)
    
    transition_prob = phi_1*phi_2*exp( r12_abs/(alpha*(1._8+beta*r12_abs)) )
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


end program vmc_h2
