      
      subroutine check_mesh(g)

!     Check the cell area and facet length calculations before attempting to
!     solve the flow, make sure you do this for both the "bump" and "bend" test
!     cases

!     Explicitly declare the required variables
      use types
      implicit none
      type(t_grid), intent(inout) :: g
      real :: area_min, dx_error, dy_error, tol
      integer :: ni, nj
      integer :: i,j
      real :: dx_sum, dy_sum

!     Get the size of the mesh and store locally for convenience
      ni = g%ni; nj = g%nj;

!     Exact checking of floating point numbers never goes well, define a
!     small tolerance to use for a comparative operation instead
      tol = 1e-4 * g%l_min

!     Check that all of the cell areas are positive, either with the intrinsic
!     "minval" function or with nested do loops. Print the output to the screen
!     and flag negative numbers as an error with an if statement to "stop" the
!     program
!     INSERT
		area_min = minval(g%area)
		if (area_min <= 0.0) then
			write(6,*) "Check Mesh Error: -ve cell area found"
			write(6,*) "min cell area = ", area_min
			stop
		else
			write(6,*) "All cell areas +ve"
		end if


!     Next check that the sum of the edge vectors around every quadrilateral is 
!     very nearly zero in both the x and y-coordinate directions. You can
!     complete this with some elementwise addition of the arrays and use of the
!     "maxval" and "abs" intrinsic functions.
!     INSERT
	dx_error = 0.0
	dy_error = 0.0
	do i = 1, ni-1
		do j = 1, nj-1
			dx_sum = (g%x(i+1,j) - g%x(i,j)) + (g%x(i+1,j+1) - g%x(i+1,j)) + &
					 (g%x(i,j+1) - g%x(i+1,j+1)) + (g%x(i,j) - g%x(i,j+1))
			dy_sum = (g%y(i+1,j) - g%y(i,j)) + (g%y(i+1,j+1) - g%y(i+1,j)) + &
					 (g%y(i,j+1) - g%y(i+1,j+1)) + (g%y(i,j) - g%y(i,j+1))
			dx_error = max(dx_error, abs(dx_sum))
			dy_error = max(dy_error, abs(dy_sum))
			! Check that dx_error and dy_error are very small
			if (dx_error > tol .or. dy_error > tol) then
				write(6,*) "Check Mesh Edge vector sum failed"
				write(6,*) "dx_error = ", dx_error, "dy_error = ", dy_error
				stop
			end if
		end do
	end do

	! Check that dx_error and dy_error are very small
	if (dx_error > tol .or. dy_error > tol) then
		write(6,*) "Check Mesh Edge vector sum failed"
		write(6,*) "dx_error = ", dx_error, "dy_error = ", dy_error
		stop
	else
		write(6,*) "Check Mesh Passed"
	end if
	

!     It may be worthwhile to complete some other checks, the prevous call to
!     the "write_output" subroutine has written a file that you can read and
!     postprocess using the Python script plot_mesh.py. This program also has
!     access to all of the mesh parameters used within the solver that you could
!     inspect graphically.
	if (dx_error > tol .or. dy_error > tol) then
		write(6,*) "Check Mesh Edge vector sum failed"
		write(6,*) "dx_error = ", dx_error, "dy_error = ", dy_error
		stop
	else
		write(6,*) "Check Mesh Passed"
	end if

!     Print a blank line
      write(6,*)

      end subroutine check_mesh
