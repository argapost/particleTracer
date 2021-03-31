subroutine p_initialize(npx, npy, npz, Lx, Ly, Lz, px, py, pz, nprtcls)

   integer :: npx, npy, npz, nprtcls
   integer :: i
   real(4) :: px(nprtcls), py(nprtcls), pz(nprtcls)
   real(4) :: z(npz), y(npy)
   real(4) :: Ym(npz, npy), Zm(npz, npy)
   real(4) :: x1, x2, x3, epsi

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

      px(i) = x1*(Lx/2)
      py(i) = x2*Ly
      pz(i) = x3*Lz
      if (i < 100) then
         print *, px(i), py(i), pz(i)
      end if

   end do

end subroutine p_initialize

subroutine meshgrid(x, y, Xm, Ym, nx, ny)
   real, intent(in) :: x(nx), y(ny)
   real, intent(out) :: Xm(ny, nx), Ym(ny, nx)
   Xm = spread(x, 1, size(y))
   Ym = spread(y, 2, size(x))
end subroutine
