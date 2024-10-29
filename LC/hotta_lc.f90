subroutine rte_longcharacteristic(ro, te, op)
    !!!----- def & init ----- !!!
        ! rte multiray def
        use test_lc_def
        implicit none
        integer :: i, j, k
        integer :: ip, jp, kp

    
        ! inout data
        real(8), dimension(nx, ny, nz), intent(in) :: ro, te, op   
        
        ! ray 
        real(8), dimension(ly-1) :: iu, ju, ku
        real(8), dimension(ly-1) :: id, jd, kd
        
        
        ! alpha & beta
        real(8), allocatable :: alpha_mhd(:, :, :), beta_mhd(:, :, :)
        real(8), allocatable :: alpha_(:, :, :), beta_(:, :, :)
        real(8), allocatable :: alpha(:, :, :), beta(:, :, :)
    
        ! delta optical depth & delta intensity
        real(8) :: alpha_up, log_alpha_up
    
    
        ! memory
        allocate(alpha_mhd(nx, ny, nz))
        allocate(beta_mhd(nx, ny, nz))
        allocate(alpha_(nx, ny+2, nz+2))
        allocate(beta_(nx, ny+2, nz+2))
        allocate(alpha(lx, ly, lz))
        allocate(beta(lx, ly, lz))
        ! allocate(delta_tu(mrad, 2, 2, 2, lx, ly, lz))
        ! allocate(delta_in(mrad, 2, 2, 2, lx, ly, lz))
        ! allocate(exp_delta_tu(mrad, 2, 2, 2, lx, ly, lz))
    
    
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
    
        
        
    !!!----- intensity -----!!!
        print *, "# intensity #"
    
        do m = 1, 1
        do idir = 2, 2
        do jdir = 1, 1
        do kdir = 1, 1
            
        print *, "m, idir, jdir, kdir: ", m, idir, jdir, kdir

        ip = ipsta(idir)

        do jp = jpsta(jdir), jpend(jdir), pstep(jdir)
        do kp = kpsta(kdir), kpend(kdir), pstep(kdir)

            do i = 1, ncell(m)

                ! upstream
                iu(i) = ip + (i-1) * step(m, 1, idir)
                ju(i) = jp + (i-1) * step(m, 2, jdir)
                ku(i) = kp + (i-1) * step(m, 3, kdir)
                

                log_alpha_up = &
                    & ((log(alpha(int(iu(i))+di(m, idir, 1), int(ju(i))+dj(m, jdir, 1), int(ku(i))+dk(m, kdir, 1))) * dc1(m) * dc2(m) &
                    & + log(alpha(int(iu(i))+di(m, idir, 2), int(ju(i))+dj(m, jdir, 2), int(ku(i))+dk(m, kdir, 2))) * dd1(m) * dc2(m) &
                    & + log(alpha(int(iu(i))+di(m, idir, 3), int(ju(i))+dj(m, jdir, 3), int(ku(i))+dk(m, kdir, 3))) * dc1(m) * dd2(m) &
                    & + log(alpha(int(iu(i))+di(m, idir, 4), int(ju(i))+dj(m, jdir, 4), int(ku(i))+dk(m, kdir, 4))) * dd1(m) * dd2(m) &
                    &  ) * di1(m) * di2(m))

                ! downstream
                id(i) = iu(i) + step(m, 1, idir)
                jd(i) = ju(i) + step(m, 2, jdir)
                kd(i) = ku(i) + step(m, 3, kdir)


            enddo

            




            


        enddo
        enddo
    
        enddo
        enddo
        enddo
        enddo
    
    end subroutine rte_longcharacteristic