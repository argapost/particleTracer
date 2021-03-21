program particleTracer

  use netcdf
  use MPI
  use :: interpolate_mod

  implicit none

  ! real(4), external :: interpolate

  integer, parameter :: nx=512
  integer, parameter :: ny=385
  integer, parameter :: nz=512

  integer, parameter :: nt=3100
  integer,parameter :: istep=1, istart=0
  integer :: it, itsave, timestep, save_every=20
  integer :: nt_saved, ncid_save

  real :: Lx, Ly, Lz, pi 
  real(4) :: dt, time, time_prev

  integer, parameter :: nprtcls=10000

  real(4) :: px(nprtcls), py(nprtcls), pz(nprtcls)
  real(4) :: pxs(nprtcls), pys(nprtcls), pzs(nprtcls)
  real(4) :: pxp, pyp, pzp
  real(4) :: pu(nprtcls), pv(nprtcls), pw(nprtcls)
  real(4) :: pup, pvp, pwp


  ! 3D arrays
  real(4) :: u(1:nx, 1:ny, 1:nz)
  real(4) :: v(1:nx, 1:ny, 1:nz)
  real(4) :: w(1:nx, 1:ny, 1:nz)

  real(4) :: grid_y(ny)
  integer :: ncid, varid(1)

  integer ::  ip, nb_procs, OMP_GET_NUM_THREADS

  REAL :: t1,t2
  integer :: nb_periodes_initial
  integer :: nb_periodes_final
  integer :: nb_periodes_max
  integer :: nb_periodes_sec
  integer :: nb_periodes     
  real ::  temps_elapsed  

  !=================================================================
  !                        Initialisations.
  !=================================================================

  nb_procs = 1

  !$OMP PARALLEL
  !$      nb_procs = OMP_GET_NUM_THREADS()
  !$OMP END PARALLEL


  write(*,*)'---------------------------'
  write(*,*)'Number of processors :',nb_procs
  write(*,*)'---------------------------'

  !-------------------------------------------------------
  !       initialization of all constant and parameters
  !-------------------------------------------------------

  pi = acos(-1.)

  Lx = 2.0 * pi
  Ly = 2.0
  Lz = 1.0 * pi

  it = 0
  itsave = 1
  timestep = (it*istep)+istart
  time = 0.
  time_prev = 0.

  nt_saved = ceiling(real(nt / save_every))

  ncid_save = 0

  CALL CPU_TIME(t1)
  CALL SYSTEM_CLOCK(COUNT_RATE=nb_periodes_sec, &
                    COUNT_MAX=nb_periodes_max)
  CALL SYSTEM_CLOCK(COUNT=nb_periodes_initial)

  ! Load y - grid
  print *, "Read y-grid"
  call io_check(nf90_open(path='grid_y.nc', &
                 mode=nf90_nowrite,ncid=ncid))

  call io_check(nf90_inq_varid(ncid,'grid_y',varid(1)))
  call io_check(nf90_get_var(ncid,varid(1),grid_y,count=(/ny/)))
  call io_check(nf90_close(ncid))

  !-------------------------------------------------------
  !              initialization of lagrangian particles
  !-------------------------------------------------------
 
  ! initialize particle position
  call p_initialize(1, 100, 100, Lx, Ly, Lz, px, py, pz, nprtcls)

  ! Load first time step
  call load_velocity(u, v, w, nx, ny, nz, timestep, time)
  time_prev = time

  ! Interpolate velocity at initial particle position

!$OMP PARALLEL
!$OMP DO SCHEDULE(RUNTIME)
  do ip=1,nprtcls

    pu(ip) = interpolate(u, grid_y, nx, ny, nz, Lx, Ly, Lz, px(ip), py(ip), pz(ip))
    pv(ip) = interpolate(v, grid_y, nx, ny, nz, Lx, Ly, Lz, px(ip), py(ip), pz(ip))
    pw(ip) = interpolate(w, grid_y, nx, ny, nz, Lx, Ly, Lz, px(ip), py(ip), pz(ip))

  enddo
!$OMP END DO
!$OMP END PARALLEL

  ! Save initial position and interpolated fields
  call p_save(grid_y, nx, ny, nz, Lx, Ly, Lz, pxs, pys, pzs, pu, pv, pw, nprtcls, nt_saved, itsave, timestep, time, ncid_save)

  do it=1,nt-1

    ! Load velocity from netcdf
    timestep = (it * istep) + istart
    call load_velocity(u, v, w, nx, ny, nz, timestep, time)

    dt = time - time_prev
    time_prev = time

    ! write(*,102)'time =',time,'it =',it,'dt =',dt
    print *, 'time =',time,'it =',it,'dt =',dt

!$OMP PARALLEL DEFAULT(SHARED), PRIVATE(pxp, pyp, pzp, pup, pvp, pwp)
!$OMP DO SCHEDULE(RUNTIME)
    do ip=1,nprtcls

      ! Find x(t+dt) with u(t) --- Predictor step
      pxp = px(ip) + dt * pu(ip)
      pyp = py(ip) + dt * pv(ip)
      pzp = pz(ip) + dt * pw(ip)

      ! If particle predictor position outside of the box reset inside the box
      if (pxp .gt. Lx) pxp = pxp - Lx
      if (pyp .gt. Ly) pyp = Ly
      if (pzp .gt. Lz) pzp = pzp - Lz

      if (pxp .le. 0.) pxp = pxp + Lx
      if (pyp .le. 0.) pyp = 0
      if (pzp .le. 0.) pzp = pzp + Lz

      ! Interpolate velocity at predictor location u(t+dt)
      pup = interpolate(u, grid_y, nx, ny, nz, Lx, Ly, Lz, pxp, pyp, pzp)
      pvp = interpolate(v, grid_y, nx, ny, nz, Lx, Ly, Lz, pxp, pyp, pzp)
      pwp = interpolate(w, grid_y, nx, ny, nz, Lx, Ly, Lz, pxp, pyp, pzp)

      ! Find x(t+dt) with (u(t) + u(t+dt))/2 --- Corrector step
      px(ip) = px(ip) + dt * (pu(ip) + pup) / 2
      py(ip) = py(ip) + dt * (pv(ip) + pvp) / 2
      pz(ip) = pz(ip) + dt * (pw(ip) + pwp) / 2

      ! Location of particles that can go outside of the box
      pxs(ip) = pxs(ip) + dt * (pu(ip) + pup) / 2
      ! pys(ip) = pys(ip) + dt * (pv(ip) + pvp) / 2
      pzs(ip) = pzs(ip) + dt * (pw(ip) + pwp) / 2

      ! If particle position outside of the box reset inside the box
      if (px(ip) .ge. Lx) px(ip) = px(ip) - Lx
      if (py(ip) .ge. Ly) py(ip) = Ly
      if (pz(ip) .ge. Lz) pz(ip) = pz(ip) - Lz

      if (px(ip) .le. 0.) px(ip) = px(ip) + Lx
      if (py(ip) .le. 0.) py(ip) = 0
      if (pz(ip) .le. 0.) pz(ip) = pz(ip) + Lz

      ! Interpolate velocity u(t+dt) at final location after correction
      pu(ip) = interpolate(u, grid_y, nx, ny, nz, Lx, Ly, Lz, px(ip), py(ip), pz(ip))
      pv(ip) = interpolate(v, grid_y, nx, ny, nz, Lx, Ly, Lz, px(ip), py(ip), pz(ip))
      pw(ip) = interpolate(w, grid_y, nx, ny, nz, Lx, Ly, Lz, px(ip), py(ip), pz(ip))

    enddo
!$OMP END DO
!$OMP END PARALLEL

    ! Save particles position and interpolated fields
    if (mod(it, save_every) .eq. 0) then
      itsave = itsave + 1
      print *, 'Save trajectories for itsave=', itsave
      call p_save(grid_y, nx, ny, nz, Lx, Ly, Lz, pxs, py, pzs, pu, pv, pw, nprtcls, nt_saved, itsave, timestep, time, ncid_save)
    endif

  enddo
  ! Close netcdf file
  call io_check(nf90_close(ncid_save))


  CALL CPU_TIME(t2)
  CALL SYSTEM_CLOCK(COUNT=nb_periodes_final)
  nb_periodes = nb_periodes_final - nb_periodes_initial
  IF (nb_periodes_final < nb_periodes_initial) THEN
    nb_periodes = nb_periodes + nb_periodes_max
  ENDIF
  temps_elapsed   =real(nb_periodes) / nb_periodes_sec
  write(*,*)'ELAPSED TIME :',temps_elapsed
  write(*,*)'CPU TIME :',t2-t1

  write(*,*)' '
  write(*,*)'END PROGRAM'
  write(*,*)' '

end program


subroutine io_check(status)
  
  use netcdf

  integer :: status
  
  if (status .ne. nf90_noerr) then
     print *, nf90_strerror(status)
     stop 'Problem with NetCDF'
  endif
  
  return
end subroutine io_check

