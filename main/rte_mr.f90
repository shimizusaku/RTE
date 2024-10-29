subroutine rte_multiray(ro, te, se, op, x, y, z, in, tu_1ray, x_tu1, qrad, qrad_tu1)
!!!----- def & init -----!!!
    ! rte multiray def
    use rte_mr_def
    implicit none
    integer :: i, j, k
    integer :: idir, jdir, kdir
    integer :: count

    logical :: converged

    ! inout data
    real(8), dimension(nx, ny, nz), intent(in) :: ro, te, se, op
    real(8), dimension(nx), intent(in) :: x, y, z

    ! output data
    real(8), dimension(mrad, 2, 2, 2, lx, ly, lz), intent(out) :: in
    real(8), dimension(nx, ly, lz), intent(out) :: tu_1ray
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
    real(8), allocatable :: delta_tu(:, :, :, :, :, :, :), delta_in(:, :, :, :, :, :, :)
    real(8), allocatable :: exp_delta_tu(:, :, :, :, :, :, :)

    ! intensity
    real(8), allocatable :: in_up(:, :, :, :, :, :, :)
    real(8) :: irecv_j, isend_j
    real(8) :: irecv_k, isend_k
    real(8), allocatable :: err_j(:, :), err_k(:, :)
    real(8) :: max_err_j, max_err_k, max_err

    ! radiative heating rate
    real(8), allocatable :: delta_tu_1ray(:, :, :)
    real(8), allocatable :: fx(:, :, :), fy(:, :, :), fz(:, :, :), j0(:, :, :)
    real(8) :: fx_up, fx_down
    real(8) :: fy_up, fy_down
    real(8) :: fz_up, fz_down
    real(8) :: j_mhd
    ! real(8) :: tu_norm
    real(8) :: qrad_f, qrad_j
    real(8) :: tmp


    
    ! memory
    allocate(alpha_mhd(nx, ny, nz))
    allocate(beta_mhd(nx, ny, nz))
    allocate(alpha_(nx, ny+2, nz+2))
    allocate(beta_(nx, ny+2, nz+2))
    allocate(alpha(lx, ly, lz))
    allocate(beta(lx, ly, lz))
    allocate(delta_tu(mrad, 2, 2, 2, lx, ly, lz))
    allocate(delta_in(mrad, 2, 2, 2, lx, ly, lz))
    allocate(in_up(mrad, 2, 2, 2, lx, ly, lz))
    allocate(err_j(lx-1, lz))
    allocate(err_k(lx-1, ly))
    allocate(delta_tu_1ray(lx, ly, lz))
    allocate(fx(lx, ly, lz))
    allocate(fy(lx, ly, lz))
    allocate(fz(lx, ly, lz))
    allocate(j0(lx, ly, lz))
    allocate(exp_delta_tu(mrad, 2, 2, 2, lx, ly, lz))



!-------------------------------------------------------------------------------------------------------!
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
    
    


!!!----- delta opatical depth & delta intensity -----!!!
    print *, "# delta opacity & delta in #"

    ! ray direction
    call ray_def()

    ! initial condision
    ! y-z plane
    do m = 1, mrad
    do jdir = 1, 2
    do kdir = 1, 2
    do j = 1, ly
    do k = 1, lz
        in(m, 1, jdir, kdir, lx, j, k) = 1.0e-10
        in(m, 2, jdir, kdir, 1,  j, k) = beta(1, j, k)
    enddo
    enddo
    enddo
    enddo
    enddo

    ! x-z plane
    do m = 1, mrad
    do idir = 1, 2
    do kdir = 1, 2
    do i = 1, lx
    do k = 1, lz
        in(m, idir, 1, kdir, i, ly, k) = beta(i, ly, k)
        in(m, idir, 2, kdir, i, 1,  k) = beta(i, 1,  k)
    enddo
    enddo
    enddo   
    enddo
    enddo

    ! x-y plane
    do m = 1, mrad
    do idir = 1, 2
    do jdir = 1, 2
    do i = 1, lx
    do j = 1, ly        
        in(m, idir, jdir, 1, i, j, lz) = beta(i, j, lz)
        in(m, idir, jdir, 2, i, j, 1 ) = beta(i, j, 1 )
    enddo
    enddo
    enddo
    enddo    
    enddo


!!!----- intensity -----!!!
    print *, "# intensity #"

    do m = 1, mrad
    do idir = 1, 2
    do jdir = 1, 2
    do kdir = 1, 2
        
        print *, "m, idir, jdir, kdir: ", m, idir, jdir, kdir

        do i = ista(idir), iend(idir), step(idir)
        do j = jsta(jdir), jend(jdir), step(jdir)
        do k = ksta(kdir), kend(kdir), step(kdir)

            ! downstream
            alpha_down     = alpha(i + shift(idir), j + shift(jdir), k + shift(kdir))
            log_alpha_down = log(alpha_down)

            beta_down     = beta(i + shift(idir), j + shift(jdir), k + shift(kdir))
            log_beta_down = log(beta_down)


            ! upstream 
            log_alpha_up = ((log(alpha(i+di(m, idir, 1), j+dj(m, jdir, 1), k+dk(m, kdir, 1))) * dc1(m) * dc2(m) &
                        & + log(alpha(i+di(m, idir, 2), j+dj(m, jdir, 2), k+dk(m, kdir, 2))) * dd1(m) * dc2(m) &
                        & + log(alpha(i+di(m, idir, 3), j+dj(m, jdir, 3), k+dk(m, kdir, 3))) * dc1(m) * dd2(m) &
                        & + log(alpha(i+di(m, idir, 4), j+dj(m, jdir, 4), k+dk(m, kdir, 4))) * dd1(m) * dd2(m) &
                        &   ) * di1(m) * di2(m))
            alpha_up = exp(log_alpha_up)

            log_beta_up = ((log(beta(i+di(m, idir, 1), j+dj(m, jdir, 1), k+dk(m, kdir, 1))) * dc1(m) * dc2(m) &
                        & + log(beta(i+di(m, idir, 2), j+dj(m, jdir, 2), k+dk(m, kdir, 2))) * dd1(m) * dc2(m) &
                        & + log(beta(i+di(m, idir, 3), j+dj(m, jdir, 3), k+dk(m, kdir, 3))) * dc1(m) * dd2(m) &
                        & + log(beta(i+di(m, idir, 4), j+dj(m, jdir, 4), k+dk(m, kdir, 4))) * dd1(m) * dd2(m) &
                        &   ) * di1(m) * di2(m))
            beta_up = exp(log_beta_up)


            ! ray lenght lray(m) in short characteristics
            lray = 0.0 
            if(m == 1) then
                lray = lr(m) * (x(i) - x(i-1))

            elseif(m == 2) then
                lray = lr(m) * (y(i) - y(i-1))
            
            elseif(m == 3) then
                lray = lr(m) * (z(i) - z(i-1))
            
            endif


            ! delta optical depth & delta intensity
            delta_tu(m, idir, jdir, kdir, i, j, k) = rte_delta_tu(alpha_down, alpha_up, lray)
            exp_delta_tu(m, idir, jdir, kdir, i, j, k) = rte_exp(- delta_tu(m, idir, jdir, kdir, i, j, k))
            ! delta_tu(m, idir, jdir, kdir, i, j, k) = ((alpha_down - alpha_up) * lray) / (log_alpha_down - log_alpha_up)

            delta_in(m, idir, jdir, kdir, i, j, k) = & 
                & rte_delta_in(beta_down, beta_up, log_beta_down, log_beta_up, delta_tu(m, idir, jdir, kdir, i, j, k))
            ! delta_in(m, idir, jdir, kdir, i, j, k) = &
            !   & ((beta_down - beta_up * exp_delta_tu(m, idir, jdir, kdir, i, j, k)) * delta_tu(m, idir, jdir, kdir, i, j, k)) &
            !     & / &
            !     & (log_beta_down - log_beta_up + delta_tu(m, idir, jdir, kdir, i, j, k))    
        enddo
        enddo                                          
        enddo


        ! calculation
        count = 0

        if(m == 1) then
            converged = .false.

            do while (.not. converged)
                do i = ista(idir), iend(idir), step(idir)
                do j = jsta(jdir), jend(jdir), step(jdir)
                do k = ksta(kdir), kend(kdir), step(kdir)
                    in_up(m, idir, jdir, kdir, i, j, k) = & 
                      & (((in(m, idir, jdir, kdir, i+di(m, idir, 1), j+dj(m, jdir, 1), k+dk(m, kdir, 1))) * dc1(m) * dc2(m) &
                      & + (in(m, idir, jdir, kdir, i+di(m, idir, 2), j+dj(m, jdir, 2), k+dk(m, kdir, 2))) * dd1(m) * dc2(m) &
                      & + (in(m, idir, jdir, kdir, i+di(m, idir, 3), j+dj(m, jdir, 3), k+dk(m, kdir, 3))) * dc1(m) * dd2(m) &
                      & + (in(m, idir, jdir, kdir, i+di(m, idir, 4), j+dj(m, jdir, 4), k+dk(m, kdir, 4))) * dd1(m) * dd2(m) &
                      &   ) * di1(m) * di2(m))
                    
                    in(m, idir, jdir, kdir, i+shift(idir), j+shift(jdir), k+shift(kdir)) = & 
                      &   in_up(m, idir, jdir, kdir, i, j, k) * exp_delta_tu(m, idir, jdir, kdir, i, j, k) & 
                      & + delta_in(m, idir, jdir, kdir, i, j, k)
                enddo
                enddo
                enddo

                ! error between sta and end
                do i = 2 + shift(idir), lx + shift(idir)
                do j = 1, ly
                do k = 1, lz
                    irecv_j = in(m, idir, jdir, kdir, i, jsta(jdir) + shift_s(jdir), k)
                    isend_j = in(m, idir, jdir, kdir, i, jend(jdir) + shift(jdir)  , k)
                    err_j(i, k) = abs(isend_j - irecv_j) / isend_j

                    irecv_k = in(m, idir, jdir, kdir, i, j, ksta(kdir) + shift_s(kdir))
                    isend_k = in(m, idir, jdir, kdir, i, j, kend(kdir) + shift(kdir)  )
                    err_k(i, j) = abs(isend_k - irecv_k) / isend_k
                enddo
                enddo
                enddo

                ! max error
                max_err_j = maxval(err_j)
                max_err_k = maxval(err_k)
                max_err = max(max_err_j, max_err_k)
                print *, max_err_j, max_err_k

                ! check
                if(max_err < 1e-3) then
                    ! Terminate the iterative calculation
                    print *, "count: ", count
                    converged = .true.
                
                else
                    ! update of sta value by end value
                    in(m, idir, jdir, kdir, :, jsta(jdir) + shift_s(jdir), :) = & 
                      &  in(m, idir, jdir, kdir, :, jend(jdir) + shift(jdir)  , :)
                    in(m, idir, jdir, kdir, :, :, ksta(kdir) + shift_s(kdir)) = &
                      &  in(m, idir, jdir, kdir, :, :, kend(kdir) + shift(kdir)  )

                    count = count + 1
                endif
            enddo
        
        elseif(m == 2) then
            converged = .false.

            do while (.not. converged)
                do j = jsta(jdir), jend(jdir), step(jdir)
                do i = ista(idir), iend(idir), step(idir)
                do k = ksta(kdir), kend(kdir), step(kdir)
                    in_up(m, idir, jdir, kdir, i, j, k) = & 
                        & (((in(m, idir, jdir, kdir, i+di(m, idir, 1), j+dj(m, jdir, 1), k+dk(m, kdir, 1))) * dc1(m) * dc2(m) &
                        & + (in(m, idir, jdir, kdir, i+di(m, idir, 2), j+dj(m, jdir, 2), k+dk(m, kdir, 2))) * dd1(m) * dc2(m) &
                        & + (in(m, idir, jdir, kdir, i+di(m, idir, 3), j+dj(m, jdir, 3), k+dk(m, kdir, 3))) * dc1(m) * dd2(m) &
                        & + (in(m, idir, jdir, kdir, i+di(m, idir, 4), j+dj(m, jdir, 4), k+dk(m, kdir, 4))) * dd1(m) * dd2(m) &
                        &   ) * di1(m) * di2(m))
                    
                    in(m, idir, jdir, kdir, i+shift(idir), j+shift(jdir), k+shift(kdir)) = & 
                        &   in_up(m, idir, jdir, kdir, i, j, k) * exp_delta_tu(m, idir, jdir, kdir, i, j, k) & 
                        & + delta_in(m, idir, jdir, kdir, i, j, k)
                enddo
                enddo
                enddo

                ! error between sta and end
                do i = 2 + shift(idir), lx + shift(idir)
                do j = 1, ly
                do k = 1, lz
                    irecv_j = in(m, idir, jdir, kdir, i, jsta(jdir) + shift_s(jdir), k)
                    isend_j = in(m, idir, jdir, kdir, i, jend(jdir) + shift(jdir)  , k)
                    err_j(i, k) = abs(isend_j - irecv_j) / isend_j

                    irecv_k = in(m, idir, jdir, kdir, i, j, ksta(kdir) + shift_s(kdir))
                    isend_k = in(m, idir, jdir, kdir, i, j, kend(kdir) + shift(kdir)  )
                    err_k(i, j) = abs(isend_k - irecv_k) / isend_k
                enddo
                enddo
                enddo

                ! max error
                max_err_j = maxval(err_j)
                max_err_k = maxval(err_k)
                max_err = max(max_err_j, max_err_k)
                print *, max_err_j, max_err_k


                ! check
                if(max_err < 1e-3) then
                    ! Terminate the iterative calculation
                    print *, "count: ", count
                    converged = .true.
                
                else
                    ! update of sta value by end value
                    in(m, idir, jdir, kdir, :, jsta(jdir) + shift_s(jdir), :) = & 
                        &  in(m, idir, jdir, kdir, :, jend(jdir) + shift(jdir)  , :)
                    in(m, idir, jdir, kdir, :, :, ksta(kdir) + shift_s(kdir)) = &
                        &  in(m, idir, jdir, kdir, :, :, kend(kdir) + shift(kdir)  )

                    count = count + 1
                endif
            enddo



        elseif(m == 3) then
            converged = .false.

            do while (.not. converged)
                do k = ksta(kdir), kend(kdir), step(kdir) 
                do i = ista(idir), iend(idir), step(idir)
                do j = jsta(jdir), jend(jdir), step(jdir)
                    in_up(m, idir, jdir, kdir, i, j, k) = & 
                        & (((in(m, idir, jdir, kdir, i+di(m, idir, 1), j+dj(m, jdir, 1), k+dk(m, kdir, 1))) * dc1(m) * dc2(m) &
                        & + (in(m, idir, jdir, kdir, i+di(m, idir, 2), j+dj(m, jdir, 2), k+dk(m, kdir, 2))) * dd1(m) * dc2(m) &
                        & + (in(m, idir, jdir, kdir, i+di(m, idir, 3), j+dj(m, jdir, 3), k+dk(m, kdir, 3))) * dc1(m) * dd2(m) &
                        & + (in(m, idir, jdir, kdir, i+di(m, idir, 4), j+dj(m, jdir, 4), k+dk(m, kdir, 4))) * dd1(m) * dd2(m) &
                        &   ) * di1(m) * di2(m))
                    
                    in(m, idir, jdir, kdir, i+shift(idir), j+shift(jdir), k+shift(kdir)) = & 
                        &   in_up(m, idir, jdir, kdir, i, j, k) * exp_delta_tu(m, idir, jdir, kdir, i, j, k) & 
                        & + delta_in(m, idir, jdir, kdir, i, j, k)
                enddo
                enddo
                enddo

                ! error between sta and end
                do i = 2 + shift(idir), lx + shift(idir)
                do j = 1, ly
                do k = 1, lz
                    irecv_j = in(m, idir, jdir, kdir, i, jsta(jdir) + shift_s(jdir), k)
                    isend_j = in(m, idir, jdir, kdir, i, jend(jdir) + shift(jdir)  , k)
                    err_j(i, k) = abs(isend_j - irecv_j) / isend_j

                    irecv_k = in(m, idir, jdir, kdir, i, j, ksta(kdir) + shift_s(kdir))
                    isend_k = in(m, idir, jdir, kdir, i, j, kend(kdir) + shift(kdir)  )
                    err_k(i, j) = abs(isend_k - irecv_k) / isend_k
                enddo
                enddo
                enddo

                ! max error
                max_err_j = maxval(err_j)
                max_err_k = maxval(err_k)
                max_err = max(max_err_j, max_err_k)
                print *, max_err_j, max_err_k


                ! check
                if(max_err < 1e-3) then
                    ! Terminate the iterative calculation
                    print *, "count: ", count
                    converged = .true.
                
                else
                    ! update of sta value by end value
                    in(m, idir, jdir, kdir, :, jsta(jdir) + shift_s(jdir), :) = & 
                        &  in(m, idir, jdir, kdir, :, jend(jdir) + shift(jdir)  , :)
                    in(m, idir, jdir, kdir, :, :, ksta(kdir) + shift_s(kdir)) = &
                        &  in(m, idir, jdir, kdir, :, :, kend(kdir) + shift(kdir)  )

                    count = count + 1
                endif
            enddo
        endif
    enddo
    enddo
    enddo
    enddo



!!!----- radiative heating rate -----!!!
    print *, "# radiative heating rate #"

    ! optical depth 1ray
    tu_1ray(lx+1, :, :) = 0.0

    do i = lx, 1, -1
    do j = 1, ly
    do k = 1, lz
        ! downstream
        alpha_down = alpha_(i, j, k)

        ! upstream
        alpha_up   = alpha_(i+1, j, k)

        ! ray lenght lray(m) in short characteristics
        lray = (x(i+1) - x(i)) / mux_1ray 

        ! delta_optical depth 1ray & optical depth 1ray
        delta_tu_1ray(i, j, k) = rte_delta_tu(alpha_down, alpha_up, lray)
        tu_1ray(i, j, k) = tu_1ray(i+1, j, k) + delta_tu_1ray(i, j, k)
        ! delta_tu_1ray(i, j, k) = ((alpha_(i, j, k) - alpha_(i-1, j, k)) * (x(i) - x(i-1))) & 
        !                      & / (log(alpha_(i, j, k)) - log(alpha_(i-1, j, k)))
        ! tu_1ray(i-1, j, k) = tu_1ray(i, j, k) + delta_tu_1ray(i, j, k)

        ! xpoint in the plane of tu=1
        if((1 - tu_1ray(i, j, k)) * (1 - tu_1ray(i+1, j, k)) <= 0.0) then
            x_tu1(j, k) = i
        endif
    enddo
    enddo
    enddo


    ! flux
    do i = 1, lx
    do j = 1, ly
    do k = 1, lz
    do m = 1, mrad
    do idir = 1, 2
    do jdir = 1, 2
    do kdir = 1, 2
        fx(i, j, k) = fx(i, j, k) + (4.d0 * pi * omr(m) * mux(m) * step(idir) * in(m, idir, jdir, kdir, i, j, k))
        fy(i, j, k) = fy(i, j, k) + (4.d0 * pi * omr(m) * muy(m) * step(jdir) * in(m, idir, jdir, kdir, i, j, k))
        fz(i, j, k) = fz(i, j, k) + (4.d0 * pi * omr(m) * muz(m) * step(kdir) * in(m, idir, jdir, kdir, i, j, k))

        j0(i, j, k) = j0(i, j, k) + (omr(m) * in(m, idir, jdir, kdir, i, j, k))
    enddo
    enddo
    enddo
    enddo
    enddo
    enddo
    enddo


    ! radiative heating rate
    do i = 2, lx
    do j = 2, ly
    do k = 2, lz
        ! qrad_f
        fx_up   = (fx(i,   j, k) + fx(i,   j-1, k) + fx(i,   j, k-1) + fx(i,   j-1, k-1)) / 4
        fx_down = (fx(i-1, j, k) + fx(i-1, j-1, k) + fx(i-1, j, k-1) + fx(i-1, j-1, k-1)) / 4

        fy_up   = (fy(i, j,   k) + fy(i-1, j,   k) + fy(i, j,   k-1) + fy(i-1, j,   k-1)) / 4
        fy_down = (fy(i, j-1, k) + fy(i-1, j-1, k) + fy(i, j-1, k-1) + fy(i-1, j-1, k-1)) / 4

        fz_up   = (fz(i, j, k  ) + fz(i, j-1, k  ) + fz(i-1, j, k  ) + fz(i-1, j-1, k  )) / 4
        fz_down = (fz(i, j, k-1) + fz(i, j-1, k-1) + fz(i-1, j, k-1) + fz(i-1, j-1, k-1)) / 4

        qrad_f = - (fx_up - fx_down) / (x(i) - x(i-1)) &
              &  - (fy_up - fy_down) / (y(i) - y(i-1)) &
              &  - (fz_up - fz_down) / (z(i) - z(i-1))
        
        ! qrad_j
        j_mhd = (j0(i,   j, k) + j0(i,   j-1, k) + j0(i,   j, k-1) + j0(i,   j-1, k-1) & 
             & + j0(i-1, j, k) + j0(i-1, j-1, k) + j0(i-1, j, k-1) + j0(i-1, j-1, k-1)) / 8

        qrad_j = 4.d0 * pi * alpha_mhd(i, j-1, k-1) * (j_mhd - beta_mhd(i, j-1, k-1))


        ! qrad
        ! tu_norm = exp(- (tu_1ray(i, j, k) / tu0))
        ! qrad(i, j, k) = (qrad_j * tu_norm + qrad_f * (1 - tu_norm)) / (ro(i, j-1, k-1) * te(i, j-1, k-1))
        
        call judge_qrad()

        tmp = min(max((log(tu_1ray(i, j, k)) - log_tumin) / delta_log_tu + 1, 1.d0), real(nx_jf_judge))
        tmp = (tmp - dble(int(tmp))) * jk_coef + jf_judge(int(tmp))
        qrad(i, j, k) = & 
            & ((qrad_j * tmp + qrad_f * (1.d0 - tmp)) / (ro(i, j-1, k-1) * te(i, j-1, k-1))) &
            & * (0.5d0 + sign(0.5d0, tu_1ray(i, j, k) - 1.d-4)) &
            & + (se(i, j-1, k-1) / te(i, j-1, k-1)) * max(0.d0, (2000.d0 - te(i, j-1, k-1)/1.d0))
    enddo
    enddo
    enddo


    ! qrad in the plane of tu=1
    do j = 2, ly
    do k = 2, lz
        qrad_tu1(j, k) = qrad(x_tu1(j, k), j, k)
    enddo
    enddo



!--------------------------------------------------------------------------------------------------------------------------!

    ! please comment 
    deallocate(alpha_mhd)
    deallocate(beta_mhd)
    deallocate(alpha_)
    deallocate(beta_)
    deallocate(alpha)
    deallocate(beta)
    deallocate(delta_tu)
    deallocate(delta_in)
    deallocate(in_up)
    deallocate(err_j)
    deallocate(err_k)
    deallocate(delta_tu_1ray)
    deallocate(fx)
    deallocate(fy)
    deallocate(fz)
    deallocate(j0)
    deallocate(exp_delta_tu)



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

end subroutine rte_multiray