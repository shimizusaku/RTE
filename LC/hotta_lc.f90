subroutine rte_longcharacteristic(ro, te, op, x, y, z, tu_1ray)
!!!----- def & init ----- !!!
    ! rte multiray def
    use test_lc_def
    implicit none
    integer :: i, j, k
    integer :: ip, jp, kp

    ! inout data
    real(8), dimension(nx, ny, nz), intent(in) :: ro, te, op   
    real(8), dimension(nx), intent(in) :: x, y, z

    ! output data
    real(8), dimension(lx, ly, lz, lz), intent(out) :: tu_1ray


    ! ray 
    real(8), dimension(ly-1) :: iu, ju, ku
    real(8), dimension(ly-1) :: id, jd, kd
    
    
    ! alpha & beta
    real(8), allocatable :: alpha_mhd(:, :, :), beta_mhd(:, :, :)
    real(8), allocatable :: alpha_(:, :, :), beta_(:, :, :)
    real(8), allocatable :: alpha(:, :, :), beta(:, :, :)

    ! delta optical depth & delta intensity
    real(8) :: alpha_up, log_alpha_up
    real(8) :: beta_up, log_beta_up
    real(8) :: alpha_down, log_alpha_down
    real(8) :: beta_down, log_beta_down
    real(8) :: lray
    real(8), allocatable :: delta_tu(:, :, :, :), exp_delta_tu(:, :, :, :)

    ! optical depth 1ray


    ! memory
    allocate(alpha_mhd(nx, ny, nz))
    allocate(beta_mhd(nx, ny, nz))
    allocate(alpha_(nx, ny+2, nz+2))
    allocate(beta_(nx, ny+2, nz+2))
    allocate(alpha(lx, ly, lz))
    allocate(beta(lx, ly, lz))
    allocate(delta_tu(lx, ly, lz, lz))
    ! allocate(delta_in(mrad, 2, 2, 2, lx, ly, lz))
    allocate(exp_delta_tu(lx, ly, lz, lz))


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
    do jdir = 2, 2
    do kdir = 2, 2
        
    print *, "m, idir, jdir, kdir: ", m, idir, jdir, kdir

    ip = ipsta(idir)

    do jp = jpsta(jdir), jpend(jdir), pstep(jdir)
        print *, "jp : ", jp

    do kp = kpsta(kdir), kpend(kdir), pstep(kdir)

        do i = 1, ncell(m)
            !!! upstream
            iu(i) = ip + (i-1) * step(m, 1, idir)
            ju(i) = jp + (i-1) * step(m, 2, jdir)
            ku(i) = kp + (i-1) * step(m, 3, kdir)

            if((iu(i) - real(ipend(idir))) > 0 ) then
                iu(i) = ipsta(idir) + (iu(i) - real(ipend(idir)))            
            endif

            if((ju(i) - real(jpend(jdir))) > 0 ) then
                ju(i) = jpsta(jdir) + (ju(i) - real(jpend(jdir)))            
            endif

            if((ku(i) - real(kpend(kdir))) > 0 ) then
                ku(i) = kpsta(kdir) + (ku(i) - real(kpend(kdir)))            
            endif


            ! alpha and beta
            duc1 = ju(i) - int(ju(i))
            dud1 = dyy - duc1
            duc2 = ku(i) - int(ku(i))
            dud2 = dzz - duc2

            log_alpha_up = &
                & ((log(alpha(int(iu(i))+di(m, idir, 1), int(ju(i))+dj(m, jdir, 1), int(ku(i))+dk(m, kdir, 1))) * dud1 * dud2 &
                & + log(alpha(int(iu(i))+di(m, idir, 2), int(ju(i))+dj(m, jdir, 2), int(ku(i))+dk(m, kdir, 2))) * duc1 * dud2 &
                & + log(alpha(int(iu(i))+di(m, idir, 3), int(ju(i))+dj(m, jdir, 3), int(ku(i))+dk(m, kdir, 3))) * dud1 * duc2 &
                & + log(alpha(int(iu(i))+di(m, idir, 4), int(ju(i))+dj(m, jdir, 4), int(ku(i))+dk(m, kdir, 4))) * duc1 * duc2 &
            &  ) * di1(m) * di2(m))
            alpha_up = exp(log_alpha_up)

            log_beta_up = &
            & ((log(beta(int(iu(i))+di(m, idir, 1), int(ju(i))+dj(m, jdir, 1), int(ku(i))+dk(m, kdir, 1))) * dud1 * dud2 &
            & + log(beta(int(iu(i))+di(m, idir, 2), int(ju(i))+dj(m, jdir, 2), int(ku(i))+dk(m, kdir, 2))) * duc1 * dud2 &
            & + log(beta(int(iu(i))+di(m, idir, 3), int(ju(i))+dj(m, jdir, 3), int(ku(i))+dk(m, kdir, 3))) * dud1 * duc2 &
            & + log(beta(int(iu(i))+di(m, idir, 4), int(ju(i))+dj(m, jdir, 4), int(ku(i))+dk(m, kdir, 4))) * duc1 * duc2 &
            &  ) * di1(m) * di2(m))
            beta_up = exp(log_beta_up)



            !!! downstream
            id(i) = iu(i) + step(m, 1, idir)
            jd(i) = ju(i) + step(m, 2, jdir)
            kd(i) = ku(i) + step(m, 3, kdir)

            if((id(i) - real(ipend(idir))) > 0 ) then
                id(i) = ipsta(idir) + (id(i) - real(ipend(idir)))
            endif

            if((jd(i) - real(jpend(jdir))) > 0 ) then
                jd(i) = jpsta(jdir) + (jd(i) - real(jpend(jdir)))                
            endif

            if((kd(i) - real(kpend(kdir))) > 0 ) then
                kd(i) = kpsta(kdir) + (kd(i) - real(kpend(kdir)))                
            endif


            ! alpha and beta
            ddc1 = jd(i) - int(jd(i))
            ddd1 = dyy - ddc1
            ddc2 = kd(i) - int(kd(i))
            ddd2 = dzz - ddc2

            log_alpha_down = &
                & ((log(alpha(int(id(i))+di(m, idir, 1), int(jd(i))+dj(m, jdir, 1), int(kd(i))+dk(m, kdir, 1))) * ddd1 * ddd2 &
                & + log(alpha(int(id(i))+di(m, idir, 2), int(jd(i))+dj(m, jdir, 2), int(kd(i))+dk(m, kdir, 2))) * ddc1 * ddd2 &
                & + log(alpha(int(id(i))+di(m, idir, 3), int(jd(i))+dj(m, jdir, 3), int(kd(i))+dk(m, kdir, 3))) * ddd1 * ddc2 &
                & + log(alpha(int(id(i))+di(m, idir, 4), int(jd(i))+dj(m, jdir, 4), int(kd(i))+dk(m, kdir, 4))) * ddc1 * ddc2 &
            &  ) * di1(m) * di2(m))
            alpha_down = exp(log_alpha_down)
            
            log_beta_down = &
            & ((log(beta(int(id(i))+di(m, idir, 1), int(jd(i))+dj(m, jdir, 1), int(kd(i))+dk(m, kdir, 1))) * ddd1 * ddd2 &
            & + log(beta(int(id(i))+di(m, idir, 2), int(jd(i))+dj(m, jdir, 2), int(kd(i))+dk(m, kdir, 2))) * ddc1 * ddd2 &
            & + log(beta(int(id(i))+di(m, idir, 3), int(jd(i))+dj(m, jdir, 3), int(kd(i))+dk(m, kdir, 3))) * ddd1 * ddc2 &
            & + log(beta(int(id(i))+di(m, idir, 4), int(jd(i))+dj(m, jdir, 4), int(kd(i))+dk(m, kdir, 4))) * ddc1 * ddc2 &
            &  ) * di1(m) * di2(m))
            beta_down = exp(log_beta_down)


            !!! ray lenght lray(m) in short characteristics
            lray = 0.d0
            if(m == 1) then
                lray = lr(m) * (x(i) - x(i-1))

            elseif(m == 2) then
                lray = lr(m) * (y(i) - y(i-1))
            
            elseif(m == 3) then
                lray = lr(m) * (z(i) - z(i-1))
            
            endif


            !!! delta optical depth & delta intensity
            delta_tu(ip, jp, kp, i) = rte_delta_tu(alpha_down, alpha_up, lray)
            exp_delta_tu(ip, jp, kp, i) = rte_exp(- delta_tu(ip, jp, kp, i))

        enddo

        



    enddo
    enddo

    enddo
    enddo
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


contains
!-------------------------------------------------------------------------------------------------------!
    !!!----- function -----!!!
    ! delta_tu
    real(8) function rte_delta_tu(alpha_down, alpha_up, lray) result(delta_tu)
        use test_lc_def
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


    ! ! delta_in
    ! real(8) function rte_delta_in(beta_down, beta_up, log_beta_down, log_beta_up, delta_tu) result(delta_in)
    !     use test_lc_def
    !     implicit none
    !     real(KIND(0.d0)), intent(in) :: beta_down, beta_up
    !     real(KIND(0.d0)), intent(in) :: log_beta_down, log_beta_up
    !     real(KIND(0.d0)), intent(in) :: delta_tu

    !     real(KIND(0.d0)) :: exp_delta_tu
    !     real(KIND(0.d0)) :: beta_upt
    !     real(KIND(0.d0)) :: beta_du, beta_dum
    !     real(KIND(0.d0)) :: log_beta_du_delta_tu, log_beta_du_delta_tu0, log_beta_du_delta_tu1
    !     real(KIND(0.d0)) :: beta_du15, beta_du15m
    !     real(KIND(0.d0)) :: delta_in0, delta_in1
    !     real(KIND(0.d0)) :: tmp

    ! !----------------------------------------------------------------------------------------| 
        
    !     exp_delta_tu = rte_exp(- delta_tu)

    !     beta_upt = beta_up * exp_delta_tu

    !     beta_du = beta_upt / beta_down - 1.d0
    !     beta_dum = sign(1.d0, beta_du) * max(abs(beta_du), epc) 

    !     beta_du15 = 1.d0 - 0.5d0 * beta_du
    !     ! for zero divide
    !     beta_du15m = sign(1.d0, beta_du15) * max(abs(beta_du15), epc) 

    !     tmp = 0.5d0 + sign(0.5d0, abs(beta_du) - epc)
    !     log_beta_du_delta_tu0 = log_beta_down - log_beta_up + delta_tu
    !     log_beta_du_delta_tu1 = 1.d0 ! value is not important

    !     log_beta_du_delta_tu = log_beta_du_delta_tu0 * tmp + log_beta_du_delta_tu1 * (1.d0 - tmp)

    !     delta_in0 = (beta_down - beta_upt) / log_beta_du_delta_tu
    !     delta_in1 = beta_down / beta_du15m

    !     delta_in = (delta_in0 * tmp + delta_in1 * (1.d0 - tmp)) * delta_tu
        
    !     return
    ! end function rte_delta_in


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

end subroutine rte_longcharacteristic