     
      subroutine read_settings(av,bcs)
                  use types
                  implicit none
                  type(t_appvars), intent(out) :: av
                  type(t_bconds), intent(out) :: bcs
                  character(len=80) :: tempname
                  read(5,*) tempname
                  av%casename = trim(tempname)
                  read(5,*) av%rgas, av%gam
                  read(5,*) av%cfl, av%sfac, av%d_max
                  read(5,*) av%nsteps
                  read(5,*) av%ni,av%nj
                  av%cp = av%rgas * av%gam / (av%gam - 1.0)
                  av%cv = av%cp / av%gam
                  av%fgam = (av%gam - 1.0) / av%gam
                  av%sfac = av%sfac * av%cfl
                  av%d_max = av%d_max * av%cfl
                  av%d_avg = 0.5 * av%d_max
                  read(5,*) bcs%pstag, bcs%tstag, bcs%alpha, bcs%rfin
                  bcs%alpha = bcs%alpha * 3.14159 / 180.0
                  bcs%rostag = bcs%pstag / (av%rgas * bcs%tstag )
                  read(5,*) bcs%p_out
                  read(5,*) av%facsec, av%fcorr
                  write(6,*) '-------------------------------------------------------'
                  write(6,*) 'Solver begins on ', av%casename, ' case'
                  write(6,*)
                  write(6,*) 'Read application variables from file'
                  write(6,*) '  rgas =', av%rgas, 'cp =', av%cp, 'cv =', av%cv
                  write(6,*) '  CFL =', av%cfl, 'sfac =', av%sfac
                  write(6,*) '  Convergence  d_max =', av%d_max
                  write(6,*) '  Mesh size  ni =', av%ni, 'nj =', av%nj
                  write(6,*) '  Inlet  pstag =', bcs%pstag, 'tstag =', bcs%tstag, &
                      'alpha = ', bcs%alpha
                  write(6,*) '  Outlet  p_out =', bcs%p_out
                  write(6,*) '  Facsec = ', av%facsec, 'Fcorr = ', av%fcorr
                  write(6,*) '-------------------------------------------------------'
                  end subroutine read_settings