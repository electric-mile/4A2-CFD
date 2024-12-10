      
      subroutine set_timestep(av,g,bcs)

!     This subroutine sets a single value for the time step based on the 
!     stagnation speed of sound and the minimum length scale of any element

!     Explicitly declare the required variables
      use types
      implicit none
      type(t_appvars), intent(inout) :: av
      type(t_grid), intent(in) :: g
      type(t_bconds), intent(in) :: bcs
      ! astag: stagnation speed of sound
      ! v_max: maximum flow speed
      real :: astag, v_max

!     Calculate the stagnation speed of sound from the inlet stagnation
      if (bcs%tstag > 0.0) then
          astag = sqrt(av%gam * av%rgas * bcs%tstag)
      else
          write(6,*) 'Set Timestep Error: bcs%tstag must be positive.'
          astag = 0.0
      end if
!     INSERT
      astag = sqrt(av%gam * av%rgas * bcs%tstag)

!     Assume that the maximum flow speed is also equal to "astag". This will 
!     be pessimistic for subsonic flows but may be optimistic for supersonic 
!     flows. In the latter case the length of the time step as determined by 
!     may need to be reduced by improving this routine or varying the CFL number
!     INSERT
      v_max = 2.0 * astag
      if (g%l_min <= 0.0) then
            write(6,*) 'Set Timestep Error: g%l_min must be positive and non-zero.'
            return
      end if

!     Calculate the timestep using the CFL number and store it in "av%dt"
!     INSERT
      !av%dt = av%cfl * g%l_min / v_max

      ! Runge-Kutta
      av%dt_total = av%cfl * g%l_min / v_max
      !av%dt_total = av%dt

!     Print the calculated timestep and some intermediate values
!     INSERT
      write(6,*) 'Calculated timestep: ', av%dt, &
                 ', Stagnation speed of sound: ', astag, &
                 ', Maximum velocity: ', v_max, &
                 ', Minimum length scale: ', g%l_min
      end subroutine set_timestep


