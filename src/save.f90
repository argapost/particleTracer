subroutine p_save(grid_y, nx, ny, nz, Lx, Ly, Lz, &
                  px, pz, pxs, py, pzs, pu, pv, pw, dumdy, duvdy, dvvdy, &
                  nprtcls, nspheres, nt_saved, itsave, timestep, time, output_fn, ncid_save)
  use netcdf
  use :: interpolate_mod

  integer :: nx, ny, nz, nprtcls, nspheres, nt_saved
  integer :: timestep, itsave
  integer :: ip, is
  real(4) :: Lx, Ly, Lz, pi
  real(4) :: grid_y(ny)

  real(4) :: time

  real(4) :: px(nprtcls, nspheres), py(nprtcls, nspheres), pz(nprtcls, nspheres), pxs(nprtcls, nspheres), pzs(nprtcls, nspheres)
  real(4) :: pufl(nprtcls, nspheres), pu(nprtcls, nspheres), pv(nprtcls, nspheres), pw(nprtcls, nspheres)
  real(4) :: pp(nprtcls, nspheres), peps(nprtcls, nspheres)
  real(4) :: pdudx(nprtcls, nspheres), pdvdx(nprtcls, nspheres), pdwdx(nprtcls, nspheres), pdpdx(nprtcls, nspheres)
  real(4) :: pdudy(nprtcls, nspheres), pdvdy(nprtcls, nspheres), pdwdy(nprtcls, nspheres), pdpdy(nprtcls, nspheres)
  real(4) :: pdudz(nprtcls, nspheres), pdvdz(nprtcls, nspheres), pdwdz(nprtcls, nspheres), pdpdz(nprtcls, nspheres)
  real(4) :: pdudxdx(nprtcls, nspheres), pdvdxdx(nprtcls, nspheres), pdwdxdx(nprtcls, nspheres)
  real(4) :: pdudydy(nprtcls, nspheres), pdvdydy(nprtcls, nspheres), pdwdydy(nprtcls, nspheres)
  real(4) :: pdudzdz(nprtcls, nspheres), pdvdzdz(nprtcls, nspheres), pdwdzdz(nprtcls, nspheres)
  real(4) :: pdumdy(nprtcls, nspheres), pduvdy(nprtcls, nspheres), pdvvdy(nprtcls, nspheres)
  real(4) :: pdudt(nprtcls, nspheres), pdvdt(nprtcls, nspheres), pdwdt(nprtcls, nspheres)

  ! 3D arrays
  real(4) :: eps(nx, ny, nz), p(nx, ny, nz), ufl(nx, ny, nz)
  real(4) :: dudx(nx, ny, nz), dvdx(nx, ny, nz), dwdx(nx, ny, nz), dpdx(nx, ny, nz)
  real(4) :: dudy(nx, ny, nz), dvdy(nx, ny, nz), dwdy(nx, ny, nz), dpdy(nx, ny, nz)
  real(4) :: dudz(nx, ny, nz), dvdz(nx, ny, nz), dwdz(nx, ny, nz), dpdz(nx, ny, nz)
  real(4) :: dudxdx(nx, ny, nz), dvdxdx(nx, ny, nz), dwdxdx(nx, ny, nz)
  real(4) :: dudydy(nx, ny, nz), dvdydy(nx, ny, nz), dwdydy(nx, ny, nz)
  real(4) :: dudzdz(nx, ny, nz), dvdzdz(nx, ny, nz), dwdzdz(nx, ny, nz)
  real(4) :: dumdy(ny), duvdy(ny), dvvdy(ny)
  real(4) :: dudt(nx, ny, nz), dvdt(nx, ny, nz), dwdt(nx, ny, nz)

  integer :: varid_i(27), varid_o(37), startv_o(3), countv_o(3), countv_i(3), dimid(3)
  integer :: ncid, ncid_save

  character(100) :: output_fn
  character(100) :: int2char, case_fn = "re9502pipi."
  character(100) :: data_dir = "/gpfsscratch/rech/avl/ulj39ir/Cases/TCF/Jimenez/Re950/data/"
  ! Load from NetCDF

  write (int2char, '(I5.5)') timestep
  call io_check(nf90_open(path=trim(data_dir)//'u_du_ddu_eps_p/'//trim(case_fn)//trim(int2char)//'.u_du_ddu_eps_p.nc', &
                          mode=nf90_nowrite, ncid=ncid))

  call io_check(nf90_inq_varid(ncid, 'dudx', varid_i(1)))
  call io_check(nf90_inq_varid(ncid, 'dudy', varid_i(2)))
  call io_check(nf90_inq_varid(ncid, 'dudz', varid_i(3)))

  call io_check(nf90_inq_varid(ncid, 'dvdx', varid_i(4)))
  call io_check(nf90_inq_varid(ncid, 'dvdy', varid_i(5)))
  call io_check(nf90_inq_varid(ncid, 'dvdz', varid_i(6)))

  call io_check(nf90_inq_varid(ncid, 'dwdx', varid_i(7)))
  call io_check(nf90_inq_varid(ncid, 'dwdy', varid_i(8)))
  call io_check(nf90_inq_varid(ncid, 'dwdz', varid_i(9)))

  call io_check(nf90_inq_varid(ncid, 'eps', varid_i(10)))

  call io_check(nf90_inq_varid(ncid, 'p', varid_i(11)))
  call io_check(nf90_inq_varid(ncid, 'dpdx', varid_i(12)))
  call io_check(nf90_inq_varid(ncid, 'dpdy', varid_i(13)))
  call io_check(nf90_inq_varid(ncid, 'dpdz', varid_i(14)))

  call io_check(nf90_inq_varid(ncid, 'dudxdx', varid_i(15)))
  call io_check(nf90_inq_varid(ncid, 'dudydy', varid_i(16)))
  call io_check(nf90_inq_varid(ncid, 'dudzdz', varid_i(17)))

  call io_check(nf90_inq_varid(ncid, 'dvdxdx', varid_i(18)))
  call io_check(nf90_inq_varid(ncid, 'dvdydy', varid_i(19)))
  call io_check(nf90_inq_varid(ncid, 'dvdzdz', varid_i(20)))

  call io_check(nf90_inq_varid(ncid, 'dwdxdx', varid_i(21)))
  call io_check(nf90_inq_varid(ncid, 'dwdydy', varid_i(22)))
  call io_check(nf90_inq_varid(ncid, 'dwdzdz', varid_i(23)))

  call io_check(nf90_inq_varid(ncid, 'u', varid_i(24)))

  call io_check(nf90_inq_varid(ncid, 'dudt', varid_i(25)))
  call io_check(nf90_inq_varid(ncid, 'dvdt', varid_i(26)))
  call io_check(nf90_inq_varid(ncid, 'dwdt', varid_i(27)))

  countv_i(1) = nx
  countv_i(2) = ny
  countv_i(3) = nz

  call io_check(nf90_get_var(ncid, varid_i(1), dudx, count=countv_i))
  call io_check(nf90_get_var(ncid, varid_i(2), dudy, count=countv_i))
  call io_check(nf90_get_var(ncid, varid_i(3), dudz, count=countv_i))

  call io_check(nf90_get_var(ncid, varid_i(4), dvdx, count=countv_i))
  call io_check(nf90_get_var(ncid, varid_i(5), dvdy, count=countv_i))
  call io_check(nf90_get_var(ncid, varid_i(6), dvdz, count=countv_i))

  call io_check(nf90_get_var(ncid, varid_i(7), dwdx, count=countv_i))
  call io_check(nf90_get_var(ncid, varid_i(8), dwdy, count=countv_i))
  call io_check(nf90_get_var(ncid, varid_i(9), dwdz, count=countv_i))

  call io_check(nf90_get_var(ncid, varid_i(10), eps, count=countv_i))

  call io_check(nf90_get_var(ncid, varid_i(11), p, count=countv_i))
  call io_check(nf90_get_var(ncid, varid_i(12), dpdx, count=countv_i))
  call io_check(nf90_get_var(ncid, varid_i(13), dpdy, count=countv_i))
  call io_check(nf90_get_var(ncid, varid_i(14), dpdz, count=countv_i))

  call io_check(nf90_get_var(ncid, varid_i(15), dudxdx, count=countv_i))
  call io_check(nf90_get_var(ncid, varid_i(16), dudydy, count=countv_i))
  call io_check(nf90_get_var(ncid, varid_i(17), dudzdz, count=countv_i))

  call io_check(nf90_get_var(ncid, varid_i(18), dvdxdx, count=countv_i))
  call io_check(nf90_get_var(ncid, varid_i(19), dvdydy, count=countv_i))
  call io_check(nf90_get_var(ncid, varid_i(20), dvdzdz, count=countv_i))

  call io_check(nf90_get_var(ncid, varid_i(21), dwdxdx, count=countv_i))
  call io_check(nf90_get_var(ncid, varid_i(22), dwdydy, count=countv_i))
  call io_check(nf90_get_var(ncid, varid_i(23), dwdzdz, count=countv_i))

  call io_check(nf90_get_var(ncid, varid_i(24), ufl, count=countv_i))

  call io_check(nf90_get_var(ncid, varid_i(25), dudt, count=countv_i))
  call io_check(nf90_get_var(ncid, varid_i(26), dvdt, count=countv_i))
  call io_check(nf90_get_var(ncid, varid_i(27), dwdt, count=countv_i))

  call io_check(nf90_close(ncid))

  do is = 1, nspheres
!$OMP PARALLEL
!$OMP DO SCHEDULE(RUNTIME)
    do ip = 1, nprtcls

      pufl(ip, is) = interpolate(ufl, grid_y, nx, ny, nz, Lx, Ly, Lz, px(ip, is), py(ip, is), pz(ip, is))

      pdudx(ip, is) = interpolate(dudx, grid_y, nx, ny, nz, Lx, Ly, Lz, px(ip, is), py(ip, is), pz(ip, is))
      pdudy(ip, is) = interpolate(dudy, grid_y, nx, ny, nz, Lx, Ly, Lz, px(ip, is), py(ip, is), pz(ip, is))
      pdudz(ip, is) = interpolate(dudz, grid_y, nx, ny, nz, Lx, Ly, Lz, px(ip, is), py(ip, is), pz(ip, is))

      pdvdx(ip, is) = interpolate(dvdx, grid_y, nx, ny, nz, Lx, Ly, Lz, px(ip, is), py(ip, is), pz(ip, is))
      pdvdy(ip, is) = interpolate(dvdy, grid_y, nx, ny, nz, Lx, Ly, Lz, px(ip, is), py(ip, is), pz(ip, is))
      pdvdz(ip, is) = interpolate(dvdz, grid_y, nx, ny, nz, Lx, Ly, Lz, px(ip, is), py(ip, is), pz(ip, is))

      pdwdx(ip, is) = interpolate(dwdx, grid_y, nx, ny, nz, Lx, Ly, Lz, px(ip, is), py(ip, is), pz(ip, is))
      pdwdy(ip, is) = interpolate(dwdy, grid_y, nx, ny, nz, Lx, Ly, Lz, px(ip, is), py(ip, is), pz(ip, is))
      pdwdz(ip, is) = interpolate(dwdz, grid_y, nx, ny, nz, Lx, Ly, Lz, px(ip, is), py(ip, is), pz(ip, is))

      peps(ip, is) = interpolate(eps, grid_y, nx, ny, nz, Lx, Ly, Lz, px(ip, is), py(ip, is), pz(ip, is))

      pp(ip, is) = interpolate(p, grid_y, nx, ny, nz, Lx, Ly, Lz, px(ip, is), py(ip, is), pz(ip, is))
      pdpdx(ip, is) = interpolate(dpdx, grid_y, nx, ny, nz, Lx, Ly, Lz, px(ip, is), py(ip, is), pz(ip, is))
      pdpdy(ip, is) = interpolate(dpdy, grid_y, nx, ny, nz, Lx, Ly, Lz, px(ip, is), py(ip, is), pz(ip, is))
      pdpdz(ip, is) = interpolate(dpdz, grid_y, nx, ny, nz, Lx, Ly, Lz, px(ip, is), py(ip, is), pz(ip, is))

      pdudxdx(ip, is) = interpolate(dudxdx, grid_y, nx, ny, nz, Lx, Ly, Lz, px(ip, is), py(ip, is), pz(ip, is))
      pdudydy(ip, is) = interpolate(dudydy, grid_y, nx, ny, nz, Lx, Ly, Lz, px(ip, is), py(ip, is), pz(ip, is))
      pdudzdz(ip, is) = interpolate(dudzdz, grid_y, nx, ny, nz, Lx, Ly, Lz, px(ip, is), py(ip, is), pz(ip, is))

      pdvdxdx(ip, is) = interpolate(dvdxdx, grid_y, nx, ny, nz, Lx, Ly, Lz, px(ip, is), py(ip, is), pz(ip, is))
      pdvdydy(ip, is) = interpolate(dvdydy, grid_y, nx, ny, nz, Lx, Ly, Lz, px(ip, is), py(ip, is), pz(ip, is))
      pdvdzdz(ip, is) = interpolate(dvdzdz, grid_y, nx, ny, nz, Lx, Ly, Lz, px(ip, is), py(ip, is), pz(ip, is))

      pdwdxdx(ip, is) = interpolate(dwdxdx, grid_y, nx, ny, nz, Lx, Ly, Lz, px(ip, is), py(ip, is), pz(ip, is))
      pdwdydy(ip, is) = interpolate(dwdydy, grid_y, nx, ny, nz, Lx, Ly, Lz, px(ip, is), py(ip, is), pz(ip, is))
      pdwdzdz(ip, is) = interpolate(dwdzdz, grid_y, nx, ny, nz, Lx, Ly, Lz, px(ip, is), py(ip, is), pz(ip, is))

      pdudt(ip, is) = interpolate(dudt, grid_y, nx, ny, nz, Lx, Ly, Lz, px(ip, is), py(ip, is), pz(ip, is))
      pdvdt(ip, is) = interpolate(dvdt, grid_y, nx, ny, nz, Lx, Ly, Lz, px(ip, is), py(ip, is), pz(ip, is))
      pdwdt(ip, is) = interpolate(dwdt, grid_y, nx, ny, nz, Lx, Ly, Lz, px(ip, is), py(ip, is), pz(ip, is))

      pdumdy(ip, is) = interpolate1d(dumdy, grid_y, ny, Ly, py(ip, is))
      pduvdy(ip, is) = interpolate1d(duvdy, grid_y, ny, Ly, py(ip, is))
      pdvvdy(ip, is) = interpolate1d(dvvdy, grid_y, ny, Ly, py(ip, is))

    end do
!$OMP END DO
!$OMP END PARALLEL
  end do

  ! Save to NetCDF
  if (itsave == 1) then
    call io_check(nf90_create(path=trim(data_dir)//'particles/'//trim(case_fn)//trim(output_fn)//'.nc', &
                              cmode=or(nf90_clobber, nf90_netcdf4), ncid=ncid_save))

    call io_check(nf90_def_dim(ncid_save, 'particles', nprtcls, dimid(1)))
    call io_check(nf90_def_dim(ncid_save, 'spheres', nspheres, dimid(2)))
    call io_check(nf90_def_dim(ncid_save, 't', nt_saved, dimid(3)))

    call io_check(nf90_def_var(ncid_save, 'px', nf90_float, dimid, varid_o(1)))
    call io_check(nf90_def_var(ncid_save, 'py', nf90_float, dimid, varid_o(2)))
    call io_check(nf90_def_var(ncid_save, 'pz', nf90_float, dimid, varid_o(3)))

    call io_check(nf90_def_var(ncid_save, 'pu', nf90_float, dimid, varid_o(4)))
    call io_check(nf90_def_var(ncid_save, 'pv', nf90_float, dimid, varid_o(5)))
    call io_check(nf90_def_var(ncid_save, 'pw', nf90_float, dimid, varid_o(6)))

    call io_check(nf90_def_var(ncid_save, 'pdudx', nf90_float, dimid, varid_o(7)))
    call io_check(nf90_def_var(ncid_save, 'pdudy', nf90_float, dimid, varid_o(8)))
    call io_check(nf90_def_var(ncid_save, 'pdudz', nf90_float, dimid, varid_o(9)))

    call io_check(nf90_def_var(ncid_save, 'pdvdx', nf90_float, dimid, varid_o(10)))
    call io_check(nf90_def_var(ncid_save, 'pdvdy', nf90_float, dimid, varid_o(11)))
    call io_check(nf90_def_var(ncid_save, 'pdvdz', nf90_float, dimid, varid_o(12)))

    call io_check(nf90_def_var(ncid_save, 'pdwdx', nf90_float, dimid, varid_o(13)))
    call io_check(nf90_def_var(ncid_save, 'pdwdy', nf90_float, dimid, varid_o(14)))
    call io_check(nf90_def_var(ncid_save, 'pdwdz', nf90_float, dimid, varid_o(15)))

    call io_check(nf90_def_var(ncid_save, 'pdudxdx', nf90_float, dimid, varid_o(16)))
    call io_check(nf90_def_var(ncid_save, 'pdudydy', nf90_float, dimid, varid_o(17)))
    call io_check(nf90_def_var(ncid_save, 'pdudzdz', nf90_float, dimid, varid_o(18)))

    call io_check(nf90_def_var(ncid_save, 'pdvdxdx', nf90_float, dimid, varid_o(19)))
    call io_check(nf90_def_var(ncid_save, 'pdvdydy', nf90_float, dimid, varid_o(20)))
    call io_check(nf90_def_var(ncid_save, 'pdvdzdz', nf90_float, dimid, varid_o(21)))

    call io_check(nf90_def_var(ncid_save, 'pdwdxdx', nf90_float, dimid, varid_o(22)))
    call io_check(nf90_def_var(ncid_save, 'pdwdydy', nf90_float, dimid, varid_o(23)))
    call io_check(nf90_def_var(ncid_save, 'pdwdzdz', nf90_float, dimid, varid_o(24)))

    call io_check(nf90_def_var(ncid_save, 'peps', nf90_float, dimid, varid_o(25)))

    call io_check(nf90_def_var(ncid_save, 'pp', nf90_float, dimid, varid_o(26)))
    call io_check(nf90_def_var(ncid_save, 'pdpdx', nf90_float, dimid, varid_o(27)))
    call io_check(nf90_def_var(ncid_save, 'pdpdy', nf90_float, dimid, varid_o(28)))
    call io_check(nf90_def_var(ncid_save, 'pdpdz', nf90_float, dimid, varid_o(29)))

    call io_check(nf90_def_var(ncid_save, 'pdumdy', nf90_float, dimid, varid_o(30)))
    call io_check(nf90_def_var(ncid_save, 'pduvdy', nf90_float, dimid, varid_o(31)))
    call io_check(nf90_def_var(ncid_save, 'pdvvdy', nf90_float, dimid, varid_o(32)))

    call io_check(nf90_def_var(ncid_save, 'pufl', nf90_float, dimid, varid_o(33)))

    call io_check(nf90_def_var(ncid_save, 'pdudt', nf90_float, dimid, varid_o(34)))
    call io_check(nf90_def_var(ncid_save, 'pdvdt', nf90_float, dimid, varid_o(35)))
    call io_check(nf90_def_var(ncid_save, 'pdwdt', nf90_float, dimid, varid_o(36)))

    call io_check(nf90_def_var(ncid_save, 'time', nf90_float, dimid(3), varid_o(37)))

    call io_check(nf90_enddef(ncid_save))
  end if

  startv_o(1) = 1
  startv_o(2) = 1
  startv_o(3) = itsave

  countv_o(1) = nprtcls
  countv_o(2) = nspheres
  countv_o(3) = 1

  call io_check(nf90_put_var(ncid_save, varid_o(1), pxs, startv_o, countv_o))
  call io_check(nf90_put_var(ncid_save, varid_o(2), py, startv_o, countv_o))
  call io_check(nf90_put_var(ncid_save, varid_o(3), pzs, startv_o, countv_o))

  call io_check(nf90_put_var(ncid_save, varid_o(4), pu, startv_o, countv_o))
  call io_check(nf90_put_var(ncid_save, varid_o(5), pv, startv_o, countv_o))
  call io_check(nf90_put_var(ncid_save, varid_o(6), pw, startv_o, countv_o))

  call io_check(nf90_put_var(ncid_save, varid_o(7), pdudx, startv_o, countv_o))
  call io_check(nf90_put_var(ncid_save, varid_o(8), pdudy, startv_o, countv_o))
  call io_check(nf90_put_var(ncid_save, varid_o(9), pdudz, startv_o, countv_o))

  call io_check(nf90_put_var(ncid_save, varid_o(10), pdvdx, startv_o, countv_o))
  call io_check(nf90_put_var(ncid_save, varid_o(11), pdvdy, startv_o, countv_o))
  call io_check(nf90_put_var(ncid_save, varid_o(12), pdvdz, startv_o, countv_o))

  call io_check(nf90_put_var(ncid_save, varid_o(13), pdwdx, startv_o, countv_o))
  call io_check(nf90_put_var(ncid_save, varid_o(14), pdwdy, startv_o, countv_o))
  call io_check(nf90_put_var(ncid_save, varid_o(15), pdwdz, startv_o, countv_o))

  call io_check(nf90_put_var(ncid_save, varid_o(16), pdudxdx, startv_o, countv_o))
  call io_check(nf90_put_var(ncid_save, varid_o(17), pdudydy, startv_o, countv_o))
  call io_check(nf90_put_var(ncid_save, varid_o(18), pdudzdz, startv_o, countv_o))

  call io_check(nf90_put_var(ncid_save, varid_o(19), pdvdxdx, startv_o, countv_o))
  call io_check(nf90_put_var(ncid_save, varid_o(20), pdvdydy, startv_o, countv_o))
  call io_check(nf90_put_var(ncid_save, varid_o(21), pdvdzdz, startv_o, countv_o))

  call io_check(nf90_put_var(ncid_save, varid_o(22), pdwdxdx, startv_o, countv_o))
  call io_check(nf90_put_var(ncid_save, varid_o(23), pdwdydy, startv_o, countv_o))
  call io_check(nf90_put_var(ncid_save, varid_o(24), pdwdzdz, startv_o, countv_o))

  call io_check(nf90_put_var(ncid_save, varid_o(25), peps, startv_o, countv_o))

  call io_check(nf90_put_var(ncid_save, varid_o(26), pp, startv_o, countv_o))
  call io_check(nf90_put_var(ncid_save, varid_o(27), pdpdx, startv_o, countv_o))
  call io_check(nf90_put_var(ncid_save, varid_o(28), pdpdy, startv_o, countv_o))
  call io_check(nf90_put_var(ncid_save, varid_o(29), pdpdz, startv_o, countv_o))

  call io_check(nf90_put_var(ncid_save, varid_o(30), pdumdy, startv_o, countv_o))
  call io_check(nf90_put_var(ncid_save, varid_o(31), pduvdy, startv_o, countv_o))
  call io_check(nf90_put_var(ncid_save, varid_o(32), pdvvdy, startv_o, countv_o))

  call io_check(nf90_put_var(ncid_save, varid_o(33), pufl, startv_o, countv_o))

  call io_check(nf90_put_var(ncid_save, varid_o(34), pdudt, startv_o, countv_o))
  call io_check(nf90_put_var(ncid_save, varid_o(35), pdvdt, startv_o, countv_o))
  call io_check(nf90_put_var(ncid_save, varid_o(36), pdwdt, startv_o, countv_o))

  call io_check(nf90_put_var(ncid_save, varid_o(37), time, (/startv_o(3)/)))

  ! if (itsave == nt_saved) then
  !   call io_check(nf90_close(ncid_save))
  ! endif

end subroutine p_save
