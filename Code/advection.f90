module types
      type t_grid
            real :: cfl, dt_total, a, phi_inlet, phi_start
            integer :: ni
            real, dimension(:), allocatable :: x, phi
      end type t_grid
end module types

program advection
      use types
      use routines
      implicit none
      type(t_grid) :: g
      integer :: ni, n, nsteps
      real :: x_min, x_max, dt_total, dt, dx

      write(6,*) 'Started the advection example program'
      x_min = 0; x_max = 1; ni = 51;
      dt_total = 0.4
      g%cfl = 0.4; g%ni = ni; g%a = 1; 
      g%phi_inlet = 1; g%phi_start = 0;
      allocate(g%x(ni),g%phi(ni))
      call linspace(x_min,x_max,g%x)
      dx = g%x(2) - g%x(1)
      write(6,*) '  Grid spacing =', dx
      dt = g%cfl * dx / g%a
      write(6,*) '  Timestep =', dt
      nsteps = ceiling(dt_total / dt)
      write(6,*) '  Total number of timesteps required =', nsteps, &
              ' for runtime of dt_total =', dt_total
      g%phi = g%phi_start
      do n = 0, nsteps
            g%phi(1) = g%phi_inlet
            call upwind(g%phi,dt,dx,g%a,ni)
            if(mod(n,10) == 0) then
                  write(6,'(A,I4)') '  Step', n
                  write(6,'(*(F4.2,1X))') g%phi
            end if
      end do
      open(unit=1,file='advection_output.txt')
      write(1,*) g%x; write(1,*) g%phi;
      close(1)
end program advection

subroutine upwind(phi,dt,dx,a,ni)
      implicit none
      real, dimension(ni), intent(inout) :: phi
      real, intent(in) :: dt, dx, a
      integer, intent(in) :: ni 
      real, dimension(ni-1) :: dphi
      dphi = phi(2:ni) - phi(1:ni-1)
      phi(2:ni) = phi(2:ni) - dt * a * dphi / dx 
end subroutine upwind
