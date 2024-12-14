      subroutine set_secondary(av,g)
            use types
            implicit none
            type(t_appvars), intent(in) :: av
            type(t_grid), intent(inout) :: g
            integer :: i, j
            real, dimension(:,:), allocatable :: t, v2
            g%vx = g%rovx / g%ro
            g%vy = g%rovy / g%ro
            v2 = g%vx**2 + g%vy**2
            t = (g%roe/g%ro - (0.5 * v2)) / av%cv
            g%hstag = av%cp * t + 0.5 * v2
            g%p = g%ro * av%rgas * t
      end subroutine set_secondary
      
      