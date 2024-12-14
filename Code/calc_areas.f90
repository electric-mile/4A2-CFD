      subroutine calc_areas(g)
      use types
      use routines
      implicit none
      type(t_grid), intent(inout) :: g
      integer :: ni, nj
      integer :: i, j
      real :: dxa, dya, dxb, dyb, min_temp, a, b, c, d
      real, dimension(g%ni-1, g%nj -1, 2) :: v_a, v_b
      ni = g%ni
      nj = g%nj
      do i = 1, ni-1
            do j = 1, nj-1
                  g%area = 0.5 * abs((g%x(i+1, j+1) - g%x(i, j)) * &
                                     (g%y(i, j+1) - g%y(i+1, j)) - &
                                     (g%y(i+1, j+1) - g%y(i, j)) * &
                                     (g%x(i, j+1) - g%x(i+1, j)))
            end do
      end do
      do j = 1, nj-1
            g%ly_i(:, j) = g%x(:, j) - g%x(:, j+1)
            g%lx_i(:, j) = g%y(:, j+1) - g%y(:, j)
      end do
      do i = 1, ni-1
            g%ly_j(i, :) = g%x(i+1, :) - g%x(i, :)
            g%lx_j(i, :) = g%y(i, :) - g%y(i+1, :)
      end do
      g%l_min = 100
      do i = 1, ni-1
            do j = 1, nj-1
                  min_temp = min(hypot(g%lx_i(i, j), g%ly_i(i, j)), hypot(g%lx_j(i, j), g%ly_j(i, j)))
                  if (min_temp < g%l_min) then
                        g%l_min = sqrt((min_temp * min_temp))
                  end if
            end do
      end do
      write(6, *) 'area ', minval(g%area(:, :))
      write(6, *) 'Calculated cell areas and facet lengths'
      write(6, *) '  Overall minimum element size = ', g%l_min
      write(6, *)
      end subroutine calc_areas
