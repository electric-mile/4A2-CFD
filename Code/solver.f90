program solver
use types
implicit none
type(t_appvars) :: av
type(t_bconds) :: bcs
type(t_match) :: p
type(t_geometry) :: geom
type(t_grid) :: g
real :: d_max = 1, d_avg = 1
integer :: nstep, nconv = 5, ncheck = 5, i, j
integer :: nrkut, nrkuts
real :: time
call read_settings(av,bcs)
if(av%ni /= -1) then
      call allocate_arrays(av,g,bcs)
      call read_geom(av,geom)
      call generate_mesh(geom,g)
else 
      call read_mesh(av,g,bcs,p)
end if
call calc_areas(g)
call write_output(av,g,1)
call check_mesh(g)
call flow_guess(av,g,bcs,2)
call write_output(av,g,2)
call set_timestep(av,g,bcs)
open(unit=3,file='conv_' // av%casename // '.csv')
open(unit=11,file='stopit')
write(11,*) 0; close(11);
g%corr_ro = 0.0
g%corr_roe = 0.0
g%corr_rovx = 0.0
g%corr_rovy = 0.0
time = 0.0
if (av%casename == 'tube') then
      av%nsteps = int(0.2/av%dt_total)+1
end if

nrkuts = 4
do nstep = 1, av%nsteps
      av%nstep = nstep
      
      if (time + av%dt_total > 0.2 .and. av%casename == 'tube') then
            write(6,*) 'dt before:', av%dt_total
            av%dt_total = 0.2 - time
            write(6,*) 'dt after:', av%dt_total
      end if
      time = time + av%dt_total
      !write(6,*) 'Time:', time, 'nstep:', nstep, 'dt:', av%dt_total
      g%ro_start = g%ro
      g%roe_start = g%roe
      g%rovx_start = g%rovx
      g%rovy_start = g%rovy
      do nrkut = 1, nrkuts
            av%dt = av%dt_total / (1 + nrkuts - nrkut)
            call set_secondary(av,g)
            call apply_bconds(av,g,bcs)
            call euler_iteration(av,g)
      end do
      if(mod(av%nstep,nconv) == 0) then
            call check_conv(av,g,d_avg,d_max)
      end if
      if(mod(av%nstep,ncheck) == 0) then
            call check_stop(av,g)
      end if
      if(d_max < av%d_max .and. d_avg < av%d_avg) then
            write(6,*) 'Calculation converged in', nstep,'iterations'
            exit
      end if
end do
write(6,*) 'Calculation completed after', av%nstep,'iterations'
call write_output(av,g,3)
close(3)
end program solver
