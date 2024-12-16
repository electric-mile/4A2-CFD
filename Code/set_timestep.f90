      subroutine set_timestep(av,g,bcs)
            use types
            implicit none
            type(t_appvars), intent(inout) :: av
            type(t_grid), intent(in) :: g
            type(t_bconds), intent(in) :: bcs
            real :: astag, v_max
            if (bcs%tstag > 0.0) then
                astag = sqrt(av%gam * av%rgas * bcs%tstag)
            else
                write(6,*) 'Set Timestep Error: bcs%tstag must be positive.'
                astag = 0.0
            end if
            astag = sqrt(av%gam * av%rgas * bcs%tstag)
            v_max = astag
            if (g%l_min <= 0.0) then
                write(6,*) 'Set Timestep Error: g%l_min must be positive and non-zero.'
                return
            end if
            if (av%casename == 'bump' .or. av%casename == 'bend') then
                av%dt_total = av%cfl * g%l_min / (2*v_max)
            else
                av%dt_total = av%cfl * g%l_min / (20*v_max)
            end if
            write(6,*) 'Calculated timestep: ', av%dt, &
                       ', Stagnation speed of sound: ', astag, &
                       ', Maximum velocity: ', v_max, &
                       ', Minimum length scale: ', g%l_min
      end subroutine set_timestep