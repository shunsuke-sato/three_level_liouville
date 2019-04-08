module global_variables
  implicit none

! mathematical constants
  real(8),parameter :: pi = 4d0*atan(1d0)
  complex(8),parameter :: zi = (0d0, 1d0)

! physical constants
  real(8),parameter :: ev = 1d0/27.2114d0
  real(8),parameter :: fs = 1d0/0.024189d0

! model paramters
  real(8) :: Egap, Egap_23
  real(8) :: d_12, d_23

  complex(8) :: zrho_dm(3,3)
  real(8) :: Hmat(3,3)
!  complex(8) :: zHmat(3,3)

! relaxation paramters
  real(8) :: T2_12, T2_13, T2_23

! parameters for time-propagation
  integer :: nt, nt_cycle
  real(8) :: dt, Tprop
  real(8),allocatable :: tt(:)
  

! laser paraemter
  real(8) :: E0_1, omega0_1, tpulse_1
  real(8),allocatable :: Et_1(:),Et_1_dt2(:)
  real(8) :: E0_2, omega0_2, tpulse_2, tdelay
  real(8),allocatable :: Et_2(:),Et_2_dt2(:)

end module global_variables
!-------------------------------------------------------------------------------
program main
  use global_variables
  implicit none

  call input
  call initialize
  call time_propagation

end program main
!-------------------------------------------------------------------------------
subroutine input
  use global_variables
  implicit none

  Egap = 54.8d0*ev
  Egap_23 = 0.175d0*ev

  d_12 = 1d-2
  d_23 = 1d-2

  T2_12 = 1d0/(0.27d0*ev)

  Tprop = 60d0*fs
  dt = 0.1d0

  E0_1 = 1d-2
  omega0_1 = 1.55d0*ev
  tpulse_1 = 20d0*fs

  E0_2 = 1d-4
  omega0_2 = 55d0*ev
  tpulse_2 = 0.5d0*fs
  tdelay = 0d0*fs

end subroutine input
!-------------------------------------------------------------------------------
subroutine initialize
  use global_variables
  implicit none

  zrho_dm = 0d0
  zrho_dm(1,1) = 1d0

  call init_laser

end subroutine initialize
!-------------------------------------------------------------------------------
subroutine init_lalser
  use global_variables
  implicit none
  integer :: it
  real(8) :: xx, tt

  allocate(Et_1(-1:nt+1),Et_1_dt2(-1:nt+1))
  allocate(Et_2(-1:nt+1),Et_2_dt2(-1:nt+1))
  Et_1 = 0d0; Et_1_dt2 = 0d0
  Et_2 = 0d0; Et_2_dt2 = 0d0

! pump pulse
  do it = 0, nt+1
    tt = dt*it
    xx = tt - 0.5d0*tpulse_1
    if(abs(xx)<0.5d0*tpulse_1)then
      Et_1(it) = -E0_1/omega_1*cos(pi*xx/tpulse_1)**2*sin(omega_1*xx)
    end if

    tt = dt*it+0.5d0*dt
    xx = tt - 0.5d0*tpulse_1
    if(abs(xx)<0.5d0*tpulse_1)then
      Et_1_dt2(it) = -E0_1/omega_1*cos(pi*xx/tpulse_1)**2*sin(omega_1*xx)
    end if

  end do

! probe pulse
  do it = 0, nt+1
    tt = dt*it
    xx = tt - 0.5d0*tpulse_1 - tdelay
    if(abs(xx)<0.5d0*tpulse_2)then
      Et_2(it) = -E0_2/omega_1*cos(pi*xx/tpulse_2)**4*sin(omega_2*xx)
    end if

    tt = dt*it+0.5d0*dt
    xx = tt - 0.5d0*tpulse_1 - tdelay
    if(abs(xx)<0.5d0*tpulse_2)then
      Et_2_dt2(it) = -E0_2/omega_2*cos(pi*xx/tpulse_2)**4*sin(omega_2*xx)
    end if

  end do


end subroutine init_lalser
!-------------------------------------------------------------------------------
subroutine time_propagation
  use global_variables
  implicit none
  integer :: it
  real(8) :: dipole

  open(20,file='Et_dipole.out')
  
  it = 0
  dipole = 2d0*d_12*real(zrho_dm(1,2))
  write(20,"(999e26.16e3)")dt*it,Et_1(it),Et_2(it),dipole

  do it = 0, nt

    call dt_evolve(it)
    dipole = 2d0*d_12*real(zrho_dm(1,2))
    write(20,"(999e26.16e3)")dt*(it+1),Et_1(it+1),Et_2(it+1),dipole


  end do


  close(20)

end subroutine time_propagation
!-------------------------------------------------------------------------------
subroutine dt_evolve(it)
  use global_variables
  implicit none
  integer,intent(in) :: it
  integer :: irk
  real(8) :: H12, H23
  complex(8) :: zrho_rk(3,3,4)
  complex(8) :: zrho_t(3,3)

  Hmat = 0d0
  Hmat(1,1) = -0.5d0*Egap
  Hmat(2,2) =  0.5d0*Egap
  Hmat(3,3) =  0.5d0*Egap+Egap_23


! at time, t
  H23 = Et_1(it)*d_23
  H12 = Et_2(it)*d_12
  Hmat(1,2) = H12; Hmat(2,1) = H12
  Hmat(2,3) = H23; Hmat(3,2) = H23

!RK1
  irk = 1
  zrho_t = zrho_dm
  zrho_rk(:,:,irk) = -zi*(matmul(Hmat,zrho_t)-matmul(zrho_t,Hmat))
  zrho_rk(1,2,irk) = zrho_rk(1,2,irk) -zrho_dm(1,2)/T2_12
  zrho_rk(2,1,irk) = zrho_rk(2,1,irk) -zrho_dm(2,1)/T2_12

! at time, t+dt/2
  H23 = Et_1_dt2(it)*d_23
  H12 = Et_2_dt2(it)*d_12
  Hmat(1,2) = H12; Hmat(2,1) = H12
  Hmat(2,3) = H23; Hmat(3,2) = H32

!RK2
  irk = 2
  zrho_t = zrho_dm + 0.5d0*dt*zrho_rk(:,:,1)
  zrho_rk(:,:,irk) = -zi*(matmul(Hmat,zrho_t)-matmul(zrho_t,Hmat))
  zrho_rk(1,2,irk) = zrho_rk(1,2,irk) -zrho_dm(1,2)/T2_12
  zrho_rk(2,1,irk) = zrho_rk(2,1,irk) -zrho_dm(2,1)/T2_12

!RK3
  irk = 3
  zrho_t = zrho_dm + 0.5d0*dt*zrho_rk(:,:,2)
  zrho_rk(:,:,irk) = -zi*(matmul(Hmat,zrho_t)-matmul(zrho_t,Hmat))
  zrho_rk(1,2,irk) = zrho_rk(1,2,irk) -zrho_dm(1,2)/T2_12
  zrho_rk(2,1,irk) = zrho_rk(2,1,irk) -zrho_dm(2,1)/T2_12

! at time, t+dt
  H23 = Et_1(it+1)*d_23
  H12 = Et_2(it+1)*d_12
  Hmat(1,2) = H12; Hmat(2,1) = H12
  Hmat(2,3) = H23; Hmat(3,2) = H23

!RK4
  irk = 4
  zrho_t = zrho_dm + dt*zrho_rk(:,:,3)
  zrho_rk(:,:,irk) = -zi*(matmul(Hmat,zrho_t)-matmul(zrho_t,Hmat))
  zrho_rk(1,2,irk) = zrho_rk(1,2,irk) -zrho_dm(1,2)/T2_12
  zrho_rk(2,1,irk) = zrho_rk(2,1,irk) -zrho_dm(2,1)/T2_12


  zrho_dm = zrho_dm + dt/6d0*(zrho_rk(:,:,1) &
                         +2d0*zrho_rk(:,:,2) &
                         +2d0*zrho_rk(:,:,3) &
                         +    zrho_rk(:,:,4))

end subroutine dt_evolve
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
