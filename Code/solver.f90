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
integer :: nrkut, nrkuts, ios
real :: x_analytical(1000), ro_analytical(1000), p_analytical(1000), vx_analytical(1000), mach_analytical(100)
real, allocatable :: mach(:,:), v(:,:), t(:,:), vx(:,:), p2(:,:)
real :: error_left = 0, error_right = 0, error_total = 0
real :: time
real :: x,y, x_a, x_b
call read_settings(av,bcs)
if(av%ni /= -1) then
      call allocate_arrays(av,g,bcs)
      call read_geom(av,geom)
      call generate_mesh(geom,g, av)
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
allocate(mach(g%ni,g%nj), v(g%ni,g%nj), t(g%ni,g%nj), vx(g%ni,g%nj), p2(g%ni,g%nj))
g%corr_ro = 0.0
g%corr_roe = 0.0
g%corr_rovx = 0.0
g%corr_rovy = 0.0
time = 0.0
if (av%casename == 'tube') then
      av%nsteps = int(0.2/av%dt_total)+1
end if

nrkuts = 1
do nstep = 1, av%nsteps
      av%nstep = nstep
      
      if (time + av%dt_total > 0.2 .and. av%casename == 'tube') then
            write(6,*) 'dt before:', av%dt_total
            av%dt_total = 0.2 - time
            write(6,*) 'dt after:', av%dt_total
      end if
      time = time + av%dt_total
      if (av%casename == 'tube') then
            write(6,*) 'Time:', time, 'nstep:', nstep, 'dt:', av%dt_total
      end if
      g%ro_start = g%ro
      g%roe_start = g%roe
      g%rovx_start = g%rovx
      g%rovy_start = g%rovy
      do nrkut = 1, nrkuts
            av%dt = av%dt_total / (1 + nrkuts - nrkut)
            call set_secondary(av,g)
            if (av%casename /= 'tube') then
                  call apply_bconds(av,g,bcs)
            end if
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
v = sqrt(g%rovx**2 + g%rovy**2)/g%ro
t = (g%roe - 0.5 * g%ro * (v**2)) / (g%ro * av%cv)
vx = g%rovx / g%ro
p2 = g%ro * av%rgas * t
mach = v / ((av%gam * av%rgas * t)**0.5)
if (av%casename == 'waves') then
      x_a= 0.0
      x_b= 149.7

      do i = 1, g%ni
            write(6,*) "x:", g%x(i,1)
            do j = 1, g%nj
                  x = g%x(i,j)
                  y = g%y(i,j)
            
                  if (x < x_a) then
                        !mach_anal(i,j) = 1.8
                        error_left = error_left + ((mach(i,j)-1.8)/1.8)**2
                  else if (x > x_b) then
                        !mach_anal(i,j) = 1.2
                        error_right = error_right + ((mach(i,j)-1.2)/1.2)**2
                  else if (x > x_a .and. x < x_b) then
                        if (y > (x*100/149.7)) then
                              !mach_anal(i,j) = 1.8
                              error_left = error_left + ((mach(i,j)-1.8)/1.8)**2
                        else if (y < ((82.31/24.3)*x-407.6863)) then
                              !mach_anal(i,j) = 1.2
                              error_right = error_right + ((mach(i,j)-1.2)/1.2)**2
                        end if
                  end if  
            end do
      end do
      error_total = error_left + error_right
      write(6,*) "Error in left boundary:", error_left
      write(6,*) "Error in right boundary:", error_right
      write(6,*) "Total error:", error_total
else if (av%casename == 'tube') then
      open(unit=10, file = '../Cases/sod.raw', status = 'old', action = 'read', iostat = ios)
      if (ios /= 0) then
            write(6,*) 'Error cannot open file'
            stop
      end if
      do i = 1, 1000
            read(10,*, iostat = ios) x_analytical(i), ro_analytical(i), p_analytical(i), vx_analytical(i)
            if (ios /= 0) then
                  write(6,*) 'Error: Reading failed at line', i 
                  stop
            end if
            write(6,*) 'Read line:',i, x_analytical(i), ro_analytical(i), p_analytical(i), vx_analytical(i)
      end do
      close(10)
      error_total = 0.0
      do i = 1, 1000
            error_total = error_total + ((g%ro(i,av%nj/2) - ro_analytical(i)))**2
            error_left = error_left + ((vx(i,av%nj/2) - vx_analytical(i)))**2
            error_right = error_right + ((p2(i,av%nj/2) - p_analytical(i)))**2
      end do
      error_total = sqrt(error_total/1000)
      error_left = sqrt(error_left/1000)
      error_right = sqrt(error_right/1000)
      write(6,*) 'avg ro Error:', error_total
      write(6,*) 'avg vx Error:', error_left
      write(6,*) 'avg p Error:', error_right
end if
call write_output(av,g,3)

!Error calculatiom ig





close(3)
end program solver
