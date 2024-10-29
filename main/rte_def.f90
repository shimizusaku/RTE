module rte_mr_def
  implicit none
  ! def
  real(8), parameter :: sig = 5.67e-5
  real(8), parameter :: pi = 3.14159265358979323846
  real(8), parameter :: tu0 = 0.1d0
  real(8), parameter :: epc = 1.d-4

  
  ! grid
  integer, parameter :: nx = 256, ny = 256, nz = 256
  integer, parameter :: lx = nx - 1
  integer, parameter :: ly = ny + 1, lz = nz + 1

  ! grid span
  integer, parameter :: dxx = 1, dyy = 1, dzz = 1
  integer, parameter :: dxxi = 1 / dxx, dyyi = 1 / dyy, dzzi = 1 / dzz

  ! ray direction
  integer, dimension(2) :: step
  integer, dimension(2) :: shift
  integer, dimension(2) :: shift_s
  integer, dimension(2) :: ista
  integer, dimension(2) :: iend
  integer, dimension(2) :: jsta
  integer, dimension(2) :: jend 
  integer, dimension(2) :: ksta
  integer, dimension(2) :: kend
  real(8), parameter :: mux_1ray = 1.d0 /sqrt(3.d0)

  ! short characteristics
  integer :: m
  integer, parameter :: mrad = 3
  real(8), dimension(mrad) :: lr(mrad)
  real(8), dimension(mrad), save :: mux
  real(8), dimension(mrad), save :: muy
  real(8), dimension(mrad), save :: muz
  real(8), dimension(mrad), save :: omr

  ! value to determine the index
  real(8), dimension(mrad) :: dxc, dyc, dzc
  real(8), dimension(mrad) :: dxd, dyd, dzd
  integer, dimension(mrad, 2, 4) :: di, dj, dk
  real(8), dimension(mrad) :: dc1, dd1, dc2, dd2
  integer, dimension(mrad) :: di1, di2


  ! value to determine qrad
  integer, parameter :: nx_jf_judge = 2
  integer, dimension(nx_jf_judge) :: jf_judge
  real(8), parameter :: tumin = 1.d-2, tumax = 2.d-1
  real(8) :: log_tumin, log_tumax
  real(8) :: delta_log_tu
  real(8) :: jk_coef



contains

  subroutine ray_def()
    ! ray direction
    step = [-1, 1]
    shift = [-1, 0]
    shift_s = [0, -1]
    ista = [lx, 2]
    iend = [2, lx]
    jsta = [ly, 2]
    jend = [2, ly]
    ksta = [lz, 2]
    kend = [2, lz]

    ! short characteristics
    mux = [sqrt(7.d0 / 9.d0), 1.d0 / 3.d0,       1.d0 / 3.d0       ]
    muy = [1.0d0 / 3.0d0,     sqrt(7.d0 / 9.d0), 1.d0 / 3.d0       ]
    muz = [1.d0 / 3.d0,       1.d0 / 3.d0,       sqrt(7.d0 / 9.d0) ]
    omr = [1.d0 / 24.d0,      1.d0 / 24.d0,      1.d0 / 24.d0      ]

    do m = 1, mrad
      lr(m) = min(dxx / mux(m), dyy / muy(m), dzz / muz(m))
      
      dxc(m) = mux(m) * lr(m)
      dyc(m) = muy(m) * lr(m)
      dzc(m) = muz(m) * lr(m)
  
      dxd(m) = dxx - dxc(m)
      dyd(m) = dyy - dyc(m)
      dzd(m) = dzz - dzc(m)
    enddo

    ! value to determine the index
    do m = 1, mrad   
      
      ! y-z plane
      if (dxx / mux(m) <= dyy / muy(m) .and. dxx / mux(m) <= dzz / muz(m)) then
        dc1(m) = dyc(m)
        dd1(m) = dyd(m)
        dc2(m) = dzc(m)
        dd2(m) = dzd(m)
        di1(m) = 1.0d0 / dyy
        di2(m) = 1.0d0 / dzz
        di(m, 1, :) = 0
        di(m, 2, :) = -1
        dj(m, 1, :) = [0, -1, 0, -1]
        dk(m, 1, :) = [0, 0, -1, -1]
        dj(m, 2, :) = [-1, 0, -1, 0]
        dk(m, 2, :) = [-1, -1, 0, 0]    

      ! x-z plane
      elseif (dyy / muy(m) <= dxx / mux(m) .and. dyy / muy(m) <= dzz / muz(m)) then
        dc1(m) = dxc(m)
        dd1(m) = dxd(m)
        dc2(m) = dzc(m)
        dd2(m) = dzd(m)
        di1(m) = 1.0d0 / dxx
        di2(m) = 1.0d0 / dzz
        dj(m, 1, :) = 0
        dj(m, 2, :) = -1
        di(m, 1, :) = [0, -1, 0, -1]
        dk(m, 1, :) = [0, 0, -1, -1]
        di(m, 2, :) = [-1, 0, -1, 0]
        dk(m, 2, :) = [-1, -1, 0, 0]
    
      ! x-y plane
      elseif (dzz / muz(m) <= dxx / mux(m) .and. dzz / muz(m) <= dyy / muy(m)) then
        dc1(m) = dxc(m)
        dd1(m) = dxd(m)
        dc2(m) = dyc(m)
        dd2(m) = dyd(m)
        di1(m) = 1.0d0 / dxx
        di2(m) = 1.0d0 / dyy
        dk(m, 1, :) = 0
        dk(m, 2, :) = -1
        di(m, 1, :) = [0, -1, 0, -1]
        dj(m, 1, :) = [0, 0, -1, -1]
        di(m, 2, :) = [-1, 0, -1, 0]
        dj(m, 2, :) = [-1, -1, 0, 0]
      endif        
    enddo

  end subroutine ray_def

  subroutine judge_qrad()

    jf_judge = [1, 0]
    log_tumin = log(tumin)
    log_tumax = log(tumax)
    delta_log_tu = (log_tumax - log_tumin) / dble(nx_jf_judge - 1)
    jk_coef = delta_log_tu / (log_tumin - log_tumax)

  end subroutine judge_qrad

end module rte_mr_def