subroutine rte_longcharacteristic(ro, te, op)
!!!----- def & init ----- !!!
    ! rte multiray def
    use rte_lc_def
    implicit none
    integer :: i, j, k
    integer :: ip, jp, kp
    real(8) :: im, jm, km
    integer :: im_down

    ! inout data
    real(8), dimension(nx, ny, nz), intent(in) :: ro, te, op
    ! real(8), dimension(nx, ny, nz), intent(in) :: se
    ! real(8), dimension(nx), intent(in) :: x, y, z

    
    ! alpha & beta
    real(8), allocatable :: alpha_mhd(:, :, :), beta_mhd(:, :, :)
    real(8), allocatable :: alpha_(:, :, :), beta_(:, :, :)
    real(8), allocatable :: alpha(:, :, :), beta(:, :, :)

    ! delta optical depth & delta intensity
    real(8) :: log_alpha_down, alpha_down
    ! real(8) :: log_beta_down, beta_down
    ! real(8) :: log_alpha_up, alpha_up
    ! real(8) :: log_beta_up, beta_up
    ! real(8) :: lray
    ! real(8), allocatable :: delta_tu(:, :, :, :, :, :, :), delta_in(:, :, :, :, :, :, :)
    ! real(8), allocatable :: exp_delta_tu(:, :, :, :, :, :, :)


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

        do ip = ipsta(idir), ipend(idir), pstep(idir)
        do jp = jpsta(jdir), jpend(jdir), pstep(jdir)
        do kp = kpsta(kdir), kpend(kdir), pstep(kdir)

            ! upstream boundary point
            ista(ip, jp, kp, m, idir, jdir, kdir) = real(ip) - (step(m, 1, idir) * real(abs(ip - ipsta(idir))))
            jsta(ip, jp, kp, m, idir, jdir, kdir) = real(jp) - (step(m, 2, jdir) * real(abs(jp - jpsta(jdir))))
            ksta(ip, jp, kp, m, idir, jdir, kdir) = real(kp) - (step(m, 3, kdir) * real(abs(kp - kpsta(kdir))))

            ! downstream boundary point
            ! iend(ip, jp, kp, m, idir, jdir, kdir) = real(ip) + (step(m, 1, idir) * real(abs(ip - ipend(idir))))
            ! jend(ip, jp, kp, m, idir, jdir, kdir) = real(jp) + (step(m, 2, jdir) * real(abs(jp - jpend(jdir))))
            ! kend(ip, jp, kp, m, idir, jdir, kdir) = real(kp) + (step(m, 3, kdir) * real(abs(kp - kpend(kdir))))

            ! number of cell
            ncell(:) = [lx, ly, lz]



            do i = 1, ncell(m)

                ! move point
                im = ista(ip, jp, kp, m, idir, jdir, kdir) + (step(m, 1, idir) * (i - 1))  
                jm = jsta(ip, jp, kp, m, idir, jdir, kdir) + (step(m, 2, jdir) * (i - 1)) 
                km = ksta(ip, jp, kp, m, idir, jdir, kdir) + (step(m, 2, kdir) * (i - 1))


                ! downstream
                ! alpha_down = alpha(im + shift(m, 1, idir), jm + shift(m, 2, jdir), km + shift(m, 3, kdir))
                alpha_down = ((log(alpha(int(im+shift(m, 1, idir)) + di(m, idir, 1), int(jm+shift(m, 2, jdir)) + dj(m, jdir, 1), int(km+shift(m, 3, kdir)) + dk(m, kdir, 1))) * dc1(m) * dc2(m) &
                           & + log(alpha(i+di(m, idir, 2), j+dj(m, jdir, 2), k+dk(m, kdir, 2))) * dd1(m) * dc2(m) &
                           & + log(alpha(i+di(m, idir, 3), j+dj(m, jdir, 3), k+dk(m, kdir, 3))) * dc1(m) * dd2(m) &
                           & + log(alpha(i+di(m, idir, 4), j+dj(m, jdir, 4), k+dk(m, kdir, 4))) * dd1(m) * dd2(m) &
                           &   ) * di1(m) * di2(m))

            enddo



            


        enddo
        enddo
        enddo

    enddo
    enddo
    enddo
    enddo

end subroutine rte_longcharacteristic