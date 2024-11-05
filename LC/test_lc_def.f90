module test_lc_def
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
    integer, dimension(2) :: pstep
    integer, dimension(2) :: pshift
    integer, dimension(2) :: pshift_s
    real(8), dimension(3, 3, 2) :: step
    real(8), parameter :: mux_1ray = 1.d0 /sqrt(3.d0)

    integer, dimension(2) :: ipsta
    integer, dimension(2) :: ipend
    integer, dimension(2) :: jpsta
    integer, dimension(2) :: jpend 
    integer, dimension(2) :: kpsta
    integer, dimension(2) :: kpend
  
    integer :: idir, jdir, kdir
    integer, dimension(3) :: ncell

  
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
    real(8) :: duc1, dud1, duc2, dud2
    real(8) :: ddc1, ddd1, ddc2, ddd2
    integer, dimension(mrad) :: di1, di2
  
  
  contains
  
    subroutine ray_def()
      ! ray direction
      pstep = [-1, 1]
      pshift = [-1, 0]
      pshift_s = [0, -1]
      ipsta = [lx-1, 1]
      ipend = [1, lx-1]
      jpsta = [ly-1, 1]
      jpend = [1, ly-1]
      kpsta = [lz-1, 1]
      kpend = [1, lz-1]

      ncell = [lx-1, ly-1, lz-1]
  
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

      enddo

  
      do m = 1, mrad
        lr(m) = min(dxx / mux(m), dyy / muy(m), dzz / muz(m))
      enddo
  
  
      do m = 1, mrad
        
        do idir = 1, 2
        do jdir = 1, 2
        do kdir = 1, 2
  
          step(m, 1, :) = [- dxc(m) * real(pstep(idir)), dxc(m) * real(pstep(idir))]
          step(m, 2, :) = [- dyc(m) * real(pstep(jdir)), dyc(m) * real(pstep(jdir))]
          step(m, 3, :) = [- dzc(m) * real(pstep(kdir)), dzc(m) * real(pstep(kdir))]
          
        enddo
        enddo
        enddo

      enddo
  
      ! value to determine the index
      do m = 1, mrad   
        
        ! y-z plane
        if (dxx / mux(m) <= dyy / muy(m) .and. dxx / mux(m) <= dzz / muz(m)) then

          di1(m) = 1.0d0 / dyy
          di2(m) = 1.0d0 / dzz
          di(m, 2, :) = 0
          dj(m, 2, :) = [0, 1, 0, 1]
          dk(m, 2, :) = [0, 0, 1, 1]    
  
        endif        
      enddo
  
    end subroutine ray_def

  
  end module test_lc_def