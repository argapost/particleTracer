subroutine p_initialize(npx, npy, npz, Lx, Ly, Lz, px, py, pz, nprtcls)

  integer :: npx, npy, npz, nprtcls
  integer :: i
  real(4) :: px(nprtcls), py(nprtcls), pz(nprtcls)
  real(4) :: z(npz), y(npy)
  real(4) :: Ym(npz, npy), Zm(npz, npy)
  real(4) :: x1, x2, x3, Lx, Ly, Lz, epsi

  ! y = (/( i * 0.005, i=0,npy-1)/)
  ! z = (/( i * 0.005, i=0,npz-1)/)

  ! call meshgrid(y, z, Ym, Zm, npy, npz)

  ! Initialize particle position
  ! px = 0.
  ! py = pack(Ym, .true.)
  ! pz = pack(Zm, .true.)

  epsi = 1.E-09

  ! Random initialization
  do i = 1, nprtcls

    call RANDOM_NUMBER(x1)
    call RANDOM_NUMBER(x2)
    call RANDOM_NUMBER(x3)

    x1 = min(x1, 1.-epsi)
    x1 = max(x1, epsi)
    x2 = min(x2, 1.-epsi)
    x2 = max(x2, epsi)
    x3 = min(x3, 1.-epsi)
    x3 = max(x3, epsi)

    px(i) = x1*Lx
    py(i) = x2*Ly
    pz(i) = x3*Lz

  end do

end subroutine p_initialize

! Fibonacci Uniformly Distributed points on a sphere
subroutine sphere_initialize(sx, sy, sz, nprtcls)

  integer :: nsamples
  integer :: i, j
  real(4), parameter :: pi = 3.1415925
  real(4) :: phi, theta, radius, x, y, z
  real(4) :: sx(nprtcls), sy(nprtcls), sz(nprtcls)

  phi = pi*(3.0 - sqrt(5.0)) ! golden angle in radius

  j = 1
  do i = 0, nprtcls - 1

    y = 1.0 - (i/real(nprtcls - 1.0))*2.0
    radius = sqrt(1.0 - y*y)

    theta = phi*i

    x = cos(theta)*radius
    z = sin(theta)*radius

    if (y > 0.0) then
      sx(j) = x
      sy(j) = y
      sz(j) = z

      sx(j + 1) = -x
      sy(j + 1) = -y
      sz(j + 1) = -z

      j = j + 2
    end if

  end do

end subroutine sphere_initialize

subroutine meshgrid(x, y, Xm, Ym, nx, ny)
  real, intent(in) :: x(nx), y(ny)
  real, intent(out) :: Xm(ny, nx), Ym(ny, nx)
  Xm = spread(x, 1, size(y))
  Ym = spread(y, 2, size(x))
end subroutine
