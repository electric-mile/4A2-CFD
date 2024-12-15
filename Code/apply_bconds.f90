      subroutine apply_bconds(av,g,bcs)

      use types
      implicit none
      type(t_appvars), intent(in) :: av
      type(t_grid), intent(inout) :: g
      type(t_bconds), intent(inout) :: bcs
      integer :: j
      real :: t_in(av%nj), v_in(av%nj), p_in(av%nj), ro_in

      if(av%nstep == 1) then
          bcs%ro = g%ro(1,:)

      else
          bcs%ro = bcs%rfin * g%ro(1,:) + (1 - bcs%rfin) * bcs%ro
      endif
      bcs%ro = min(bcs%ro,0.9999 * bcs%rostag)

      ! Check for valid stagnation temperature and pressure
      if (bcs%tstag <= 0.0) then
            write(6,*) 'Error: bcs%tstag must be positive.'
            stop
      end if
      if (bcs%pstag <= 0.0) then
            write(6,*) 'Error: bcs%pstag must be positive.'
            stop
      end if

      !  Calculate inlet conditions
      t_in(:) = bcs%tstag * (bcs%ro / bcs%rostag)**(av%gam - 1.0)
      ! p_in(:) = bcs%ro * t_in(:) * av%rgas
      v_in(:) = sqrt(2 * av%cp * (bcs%tstag - t_in(:)))


      ! t_in(:) = bcs%tstag*(bcs%ro/bcs%rostag)**(av%gam-1.0)
      ! g%p(1,:) = bcs%ro * t_in(:) * av%rgas
      if (av%casename == 'waves') then
            g%p(1,:) = bcs%p_out
      end if
      ! v_in(:) = sqrt(2*av%cp*(bcs%tstag-t_in(:)))
      g%vx(1,:) = v_in(:) * cos(bcs%alpha)
      g%rovx(1,:) = bcs%ro * g%vx(1,:)
      g%vy(1,:) = v_in(:) * sin(bcs%alpha)
      g%rovy(1,:) = bcs%ro * g%vy(1,:)
      g%roe(1,:) = bcs%ro * (av%cv * t_in(:) + 0.5 * v_in(:)**2)
      g%hstag(1,:) = av%cp * t_in(:) + (v_in(:)**2)/2
      ! g%p(size(g%ro, 1), :) = bcs%p_out
      
      end subroutine apply_bconds
