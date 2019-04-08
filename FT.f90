program main
  implicit none
  complex(8),parameter :: zi = (0d0, 1d0)
  real(8),parameter :: pi = 4d0*atan(1d0)
  integer,parameter :: nt = 24805
  real(8),parameter :: ev = 1d0/27.2114d0
  integer,parameter :: nw = 1000
  real(8),parameter :: wi = 40d0*ev, wf = 60d0*ev
  real(8),parameter :: dw = (wf-wi)/nw
  complex(8) :: zeps0(0:nw)
  complex(8) :: zeps(0:nw)
  real(8) :: tt(0:nt),dt
  real(8) :: Et(0:nt), dipole(0:nt)
  complex(8) :: zEw, zDw, zchi, zexpw
  real(8) :: ww
  integer :: it, iw
  real(8) :: f1,f2,f3

  do iw = 0, nw
    ww = wi + dw*iw
    zeps0(iw) = 0.93d0 + zI*0.09d0*exp(-0.01d0*(ww/ev-55d0))
  end do
  zeps0 = zeps0**2

  open(20,file='Et_dipole.out')
  read(20,*)
  do it = 0, nt
    read(20,*)tt(it),f1,Et(it),dipole(it)
  end do
  close(20)
  dt = tt(1)-tt(0)


  do iw = 0, nw
    zEw = 0d0
    zDw = 0d0
    ww = wi + dw*iw
    do it = 0, nt
      zexpw = exp(zI*ww*tt(it))
      zEw = zEw + Et(it)*zexpw
      zDw = zDw + dipole(it)*zexpw
    end do
    zEw = zEw*dt
    zDw = zDw*dt
    zchi = -zDw/zEw

    zeps(iw) = zeps0(iw) + 4d0*pi*zchi

  end do
  
  open(30,file='zeps.out')
  do iw = 0, nw
    ww = wi + dw*iw
    write(30,"(999e26.16e3)")ww,zeps(iw),sqrt(zeps(iw)),zeps0(iw)
  end do
  close(30)

end program main
