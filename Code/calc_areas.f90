      
      subroutine calc_areas(g)

!     Calculate the area of the quadrilateral cells and the lengths of the side
!     facets

!     Explicitly declare the required variables
      use types
      use routines
      implicit none
      type(t_grid), intent(inout) :: g
      integer :: ni, nj

!     Declare integers or any extra variables you need here
!     INSERT
      integer :: i, j
      real :: dxa, dya, dxb, dyb, min_temp,a,b,c,d
      real, dimension(g%ni-1, g%nj -1,2) :: v_a,v_b
!     Get the size of the mesh and store locally for convenience
      ni = g%ni; nj = g%nj;
      

!     Calculate the areas of the cells and store in g%area. The area of any
!     quadrilateral is half of the magnitude of the cross product of the two
!     vectors that form the diagonals. Check the order of your product so that
!     the values come out positive! You can do this using two nested loops in
!     the i and j-directions or in a vectorised way by indexing the coordinate
!     arrays with lists of indices
!     INSERT
      do i = 1, ni-1
            do j =1 , nj -1
                  g%area = 0.5* abs( (g%x(i+1,j+1) -g%x(i,j)) * &
                                     (g%y(i,j+1) -g%y(i+1,j)) -&
                                     (g%y(i+1, j+1) - g%y(i, j))*&
                                     (g%x(i, j+1) - g%x(i+1, j)))
            end do
      end do

      
!     Calculate the projected lengths in the x and y-directions on all of the
!     "i = const" facets and store them in g%lx_i and g%ly_i. When combined
!     together these two components define a vector that is normal to the facet,
!     pointing inwards towards the centre of the cell. This is only the case for
!     the left hand side of the cell, the vector stored in position i,j points
!     towards the centre of the i,j cell
!     INSERT
      do j=1, nj-1
            g%ly_i(:,j) = g%x(:, j) - g%x(:, j+1)
            g%lx_i(:,j) = g%y(:, j+1) - g%y(:, j)
      end do
!     Now repeat the calculation for the project lengths on the "j=const"
!     facets. 
!     INSERT
      do i=1, ni-1
            g%ly_j(i,:) = g%x(i+1, :)-g%x(i,:)
            g%lx_j(i,:) = g%y(i, :)-g%y(i+1,:)
      end do
!     Find the minimum length scale in the mesh, this is defined as the length
!     of the shortest side of all the cells. Call this length "l_min", it is used
!     to set the timestep from the CFL number. Start by calculating the lengths
!     of the i and j facets by using the intrinsic function "hypot", this avoids
!     underflow and overflow errors. Then find the overal minimum value using
!     both the "min" and "minval" functions.
!     INSERT
      g%l_min = 100
      do i=1, ni-1
            do j=1, nj-1
                  min_temp = min(hypot(g%lx_i(i,j),g%ly_i(i,j)),hypot(g%lx_j(i,j),g%ly_j(i,j)))
                  if (min_temp < g%l_min) then
                        g%l_min = sqrt((min_temp*min_temp))
                  end if
            end do
      end do
      write(6,*) 'area ', minval(g%area(:,:))
      write(6,*) 'Calculated cell areas and facet lengths'
      write(6,*) '  Overall minimum element size = ', g%l_min
      write(6,*)

      end subroutine calc_areas
