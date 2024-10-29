subroutine rte_1ray(ro, te, op, x, in, tu_1ray, x_tu1, qrad, qrad_tu1)
!!!----- def & init -----!!!
    ! rte_multiray_def
    use rte_mr_def
    implicit none
    integer :: i, j, k
    integer :: idir

    ! input data
    real(8), dimension(nx, ny, nz), intent(in) :: ro, te, op
    real(8), dimension(nx), intent(in) :: x

    ! output data
    real(8), dimension(2, lx, ly, lz), intent(out) :: in
    real(8), dimension(lx, ly, lz), intent(out) :: tu_1ray
    integer, dimension(ly, lz), intent(out) :: x_tu1
    real(8), dimension(lx, ly, lz), intent(out) :: qrad
    real(8), dimension(ly, lz), intent(out) :: qrad_tu1

    ! alpha & beta
    real(8), allocatable :: alpha_mhd(:, :, :), beta_mhd(:, :, :)
    real(8), allocatable :: alpha_(:, :, :), beta_(:, :, :)
    real(8), allocatable :: alpha(:, :, :), beta(:, :, :)

    ! delta optical depth & delta intensity
    real(8) :: log_alpha_down, alpha_down
    real(8) :: log_beta_down, beta_down
    real(8) :: log_alpha_up, alpha_up
    real(8) :: log_beta_up, beta_up
    real(8) :: lray
    real(8), allocatable :: delta_tu(:, :, :, :), delta_in(:, :, :, :)
    real(8), allocatable :: exp_delta_tu(:, :, :, :)

    real(8), allocatable :: delta_tu_1ray(:, :, :)
    real(8), allocatable :: fx(:, :, :)
    real(8) :: fx_up, fx_down


    ! memory
    allocate(alpha_mhd(nx, ny, nz))
    allocate(beta_mhd(nx, ny, nz))
    allocate(alpha_(nx, ny+2, nz+2))
    allocate(beta_(nx, ny+2, nz+2))
    allocate(alpha(lx, ly, lz))
    allocate(beta(lx, ly, lz))
    allocate(delta_tu(2, lx, ly, lz))
    allocate(delta_in(2, lx, ly, lz))
    allocate(exp_delta_tu(2, lx, ly, lz))
    allocate(delta_tu_1ray(lx, ly, lz))
    allocate(fx(lx, ly, lz))



!------------------------------------------------------------------------------------------------------------------------!
!!!----- absorption coefficient & planck function -----!!!
    print *, "# absorption coefficient & planck function #"

    ! alpha & planck function on the MHD grid
    alpha_mhd = ro * op
    beta_mhd = (sig * te*te*te*te) / pi


    ! insert
    alpha_(:, 2:ny+1, 2:nz+1) = alpha_mhd
    alpha_(:, 1, 2:nz+1)      = alpha_mhd(:, ny, :) 
    alpha_(:, ny+2, 2:nz+1)   = alpha_mhd(:, 1, :)
    alpha_(:, :, 1)           = alpha_(:, :, nz)
    alpha_(:, :, nz+2)        = alpha_(:, :, 2)

    beta_(:, 2:ny+1, 2:nz+1) = beta_mhd
    beta_(:, 1, 2:nz+1)      = beta_mhd(:, ny, :) 
    beta_(:, ny+2, 2:nz+1)   = beta_mhd(:, 1, :)
    beta_(:, :, 1)           = beta_(:, :, nz)
    beta_(:, :, nz+2)        = beta_(:, :, 2)


    ! alpha & planck function on the RAD grid
    do i = 1, lx
    do j = 1, ly
    do k = 1, lz
        alpha(i, j, k) = exp((log(alpha_(i,   j,   k)) + log(alpha_(i,   j,   k+1)) &
                         & + log(alpha_(i,   j+1, k)) + log(alpha_(i,   j+1, k+1)) &
                         & + log(alpha_(i+1, j,   k)) + log(alpha_(i+1, j,   k+1)) &
                         & + log(alpha_(i+1, j+1, k)) + log(alpha_(i+1, j+1, k+1))) &
                         & / 8)

        beta(i, j, k) = exp((log(beta_(i,   j,   k)) + log(beta_(i,   j,   k+1)) &
                         & + log(beta_(i,   j+1, k)) + log(beta_(i,   j+1, k+1)) &
                         & + log(beta_(i+1, j,   k)) + log(beta_(i+1, j,   k+1)) &
                         & + log(beta_(i+1, j+1, k)) + log(beta_(i+1, j+1, k+1))) &
                         & / 8)   
    enddo
    enddo
    enddo




!!!----- intensity -----!!!
    print *, "# intensity #"

    ! ray direction
    call ray_def()

    ! initial condision
    ! y-z plane
    do j = 1, ly
    do k = 1, lz
        in(1, lx, j, k) = 1.0e-10
        in(2, 1,  j, k) = beta(1, j, k)
    enddo
    enddo

    do idir = 1, 2
        
        print *, "idir: ", idir

        do i = ista(idir), iend(idir), step(idir)
        do j = 1, ly
        do k = 1, lz

            ! downstrean
            alpha_down     = alpha(i + shift(idir), j, k)
            log_alpha_down = log(alpha_down)

            beta_down     = beta(i + shift(idir), j, k)
            log_beta_down = log(beta_down)


            ! upstream 
            alpha_up = alpha(i + shift_s(idir), j, k)
            log_alpha_up = log(alpha_up)

            beta_up =  beta(i + shift_s(idir), j, k)
            log_beta_up = log(beta_up)


            ! ray lenght lray(m) in short characteristics
            lray = (x(i) - x(i-1)) / mux_1ray 


            ! delta optical depth & delta intensity
            ! delta_tu(idir, i, j, k) = ((alpha_down - alpha_up) * lray) / (log_alpha_down - log_alpha_up)
            delta_tu(idir, i, j, k) = rte_delta_tu(alpha_down, alpha_up, lray)
            exp_delta_tu(idir, i, j, k) = rte_exp(- delta_tu(idir, i, j, k))

            ! delta_in(idir, i, j, k) = &
            !     & ((beta_down - beta_up * exp_delta_tu(idir, i, j, k)) * delta_tu(idir, i, j, k)) &
            !     & / &
            !     & (log_beta_down - log_beta_up + delta_tu(idir, i, j, k))  
            delta_in(idir, i, j, k) = rte_delta_in(beta_down, beta_up, log_beta_down, log_beta_up, delta_tu(idir, i, j, k))
        enddo
        enddo                                          
        enddo


        ! 1ray optical depth
        tu_1ray(lx, :, :) = 0.0
        x_tu1 = lx
    
        do i = lx, 2, -1
        do j = 1, ly
        do k = 1, lz
            tu_1ray(i-1, j, k) = tu_1ray(i, j, k) + delta_tu(1, i, j, k)
    
            if((1 - tu_1ray(i, j, k)) * (1 - tu_1ray(i-1, j, k)) <= 0.0) then
                x_tu1(j, k) = i
            endif
        enddo
        enddo
        enddo


        ! intensity
        do i = ista(idir), iend(idir), step(idir)
        do j = 1, ly
        do k = 1, lz
            in(idir, i+shift(idir), j, k) =  & 
                & in(idir, i+shift_s(idir), j, k) * exp_delta_tu(idir, i, j, k) + delta_in(idir, i, j, k)
        enddo
        enddo
        enddo
    enddo



!!!----- radiative heating rate -----!!!
    print *, "# radiative heating rate #"

    ! radiative heating rate
    do i = 2, lx
    do j = 2, ly
    do k = 2, lz
        ! qrad_f
        fx_up   = (in(2, i,   j, k) - in(1, i,   j, k)) * mux_1ray * pi * 2.0
        fx_down = (in(2, i-1, j, k) - in(1, i-1, j, k)) * mux_1ray * pi * 2.0
        
        ! qrad
        qrad(i, j, k) = - (fx_up - fx_down) / ((x(i) - x(i-1)) * ro(i, j-1, k-1) * te(i, j-1, k-1))
    enddo
    enddo
    enddo


    ! qrad in the plane of tu=1
    do j = 2, ly
    do k = 2, lz
        qrad_tu1(j, k) = qrad(x_tu1(j, k), j, k)
    enddo
    enddo  


!------------------------------------------------------------------------------------------------------------------------!
  
    ! please comment 
    deallocate(alpha_mhd)
    deallocate(beta_mhd)
    deallocate(alpha_)
    deallocate(beta_)
    deallocate(alpha)
    deallocate(beta)
    deallocate(delta_tu)
    deallocate(delta_in)
    deallocate(exp_delta_tu)
    deallocate(delta_tu_1ray)
    deallocate(fx)


contains
!-------------------------------------------------------------------------------------------------------!
    !!!----- function -----!!!
    ! delta_tu
    real(8) function rte_delta_tu(alpha_down, alpha_up, lray) result(delta_tu)
        use rte_mr_def
        implicit none
        real(8), intent(in) :: alpha_down, alpha_up
        real(8), intent(in) :: lray

        real(8) :: alpha_du, alpha_dum
        real(8) :: alpha_du15, alpha_du15m
        real(8) :: delta_tu0, delta_tu1
        real(8) :: tmp

        !---------------------------------------------------------------------------!

        alpha_du = alpha_down / alpha_up - 1.d0
        ! for zero divide
        alpha_dum = sign(1.d0, alpha_du) * max(abs(alpha_du), epc) 

        alpha_du15 = 1.d0 - 0.5d0 * alpha_du
        ! for zero divide
        alpha_du15m = sign(1.d0, alpha_du15) * max(abs(alpha_du15), epc) 

        delta_tu0 = alpha_up * alpha_dum / log(alpha_dum + 1.d0)
        delta_tu1 = alpha_up / alpha_du15m

        tmp = 0.5d0 + sign(0.5d0, abs(alpha_du) - epc)
        delta_tu = (delta_tu0 * tmp + delta_tu1 * (1.d0 - tmp)) * lray

        return
    end function rte_delta_tu


    ! delta_in
    real(8) function rte_delta_in(beta_down, beta_up, log_beta_down, log_beta_up, delta_tu) result(delta_in)
        use rte_mr_def
        implicit none
        real(KIND(0.d0)), intent(in) :: beta_down, beta_up
        real(KIND(0.d0)), intent(in) :: log_beta_down, log_beta_up
        real(KIND(0.d0)), intent(in) :: delta_tu

        real(KIND(0.d0)) :: exp_delta_tu
        real(KIND(0.d0)) :: beta_upt
        real(KIND(0.d0)) :: beta_du, beta_dum
        real(KIND(0.d0)) :: log_beta_du_delta_tu, log_beta_du_delta_tu0, log_beta_du_delta_tu1
        real(KIND(0.d0)) :: beta_du15, beta_du15m
        real(KIND(0.d0)) :: delta_in0, delta_in1
        real(KIND(0.d0)) :: tmp

    !----------------------------------------------------------------------------------------| 
        
        exp_delta_tu = rte_exp(- delta_tu)

        beta_upt = beta_up * exp_delta_tu

        beta_du = beta_upt / beta_down - 1.d0
        beta_dum = sign(1.d0, beta_du) * max(abs(beta_du), epc) 

        beta_du15 = 1.d0 - 0.5d0 * beta_du
        ! for zero divide
        beta_du15m = sign(1.d0, beta_du15) * max(abs(beta_du15), epc) 

        tmp = 0.5d0 + sign(0.5d0, abs(beta_du) - epc)
        log_beta_du_delta_tu0 = log_beta_down - log_beta_up + delta_tu
        log_beta_du_delta_tu1 = 1.d0 ! value is not important

        log_beta_du_delta_tu = log_beta_du_delta_tu0 * tmp + log_beta_du_delta_tu1 * (1.d0 - tmp)

        delta_in0 = (beta_down - beta_upt) / log_beta_du_delta_tu
        delta_in1 = beta_down / beta_du15m

        delta_in = (delta_in0 * tmp + delta_in1 * (1.d0 - tmp)) * delta_tu
        
        return
    end function rte_delta_in

    ! exp(x)
    real(8) function rte_exp(x)
        implicit none
        real(8), intent(in) :: x

        real(8) :: tmp
        real(8), parameter :: xmin = -700.d0

    !---------------------------------------------------------------------------!
    
        tmp = 0.5d0 + sign(0.5d0, x - xmin)
        rte_exp = dexp(max(x, xmin)) * tmp

        return
    end function rte_exp 
    
    
end subroutine rte_1ray