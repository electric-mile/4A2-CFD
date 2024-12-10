
      subroutine euler_iteration(av,g)
      

!     This subroutine calculates the fluxes into each cell and then sums them to
!     update the primary flow properties

!     Explicitly declare the required variables
      use types
      use flux_stencil
      use smooth_stencil
      implicit none
      type(t_appvars), intent(in) :: av
      type(t_grid), intent(inout) :: g
      real, dimension(g%ni,g%nj-1) :: mass_i, flux_i
      real, dimension(g%ni-1,g%nj) :: mass_j, flux_j
      integer :: i, j, ni, nj
!     Get the block size and store locally for convenience
      ni = g%ni; nj = g%nj

      ! Check for valid grid dimensions
      if (ni <= 0 .or. nj <= 0) then
            write(6,*) 'Euler Iteration Error: Grid dimensions must be positive and non-zero.'
            stop
      end if
      !write(6,*) 'rovx: ', g%rovx(:,:)
!     Setup the continuity equation by calculating the mass flow through
!     the facets in both the i and j-directions. Store these values in
!     "mass_i" and "mass_j"
!     INSERT
      ! If lx_j is -ve that is fine as its an indication that flow is leaving cell
      mass_i(:,1:nj-1) = ((g%rovx(:,1:nj-1)+g%rovx(:,2:nj))*g%lx_i(:,1:nj-1) + &
                          (g%rovy(:,1:nj-1)+g%rovy(:,2:nj))*g%ly_i(:,1:nj-1))/2
      mass_j(1:ni-1,:) = ((g%rovx(1:ni-1,:)+g%rovx(2:ni,:))*g%lx_j(1:ni-1,:) + &
                          (g%rovy(1:ni-1,:)+g%rovy(2:ni,:))*g%ly_j(1:ni-1,:))/2

!     Apply the wall boundary condition by checking that two nodes at the
!     end of a facet are both on a wall, if so then the appropriate mass
!     flow array is set to have zero flow through that facet
      where(g%wall(1:ni-1,:) .and. g%wall(2:ni,:)) mass_j = 0 
      where(g%wall(:,1:nj-1) .and. g%wall(:,2:nj)) mass_i = 0 
!     Update the density with mass fluxes by calling "sum_fluxes"
!     INSERT
      call sum_fluxes(av, mass_i, mass_j, g%area, g%ro, g%dro, g%ro_start)

      ! Print diagnostic information
      ! if (mod(av%nstep, 100) == 0) then
      !       write(6,*) 'Iteration:', av%nstep
      !       write(6,*) 'g%dro:', g%dro
      !       write(6,*) 'g%ro:', g%ro
      ! end if

      ! Check for non-positive density values
      if (any(g%ro <= 0.0)) then
            write(6,*) 'Euler Iteration Error: Non-positive density values detected.'
            !stop
      end if

!     Setup the conservation of energy equation by calculated the enthalpy flux
!     and storing the values in "flux_i" and "flux_j", you will need "mass_i"
!     and "mass_j" from before
!     INSERT
	do j = 1, nj-1
            flux_i(:,j) = mass_i(:,j)*(g%hstag(:,j) + g%hstag(:,j+1))/2
      end do
      do i = 1, ni-1
		flux_j(i,:) = mass_j(i,:)*(g%hstag(i,:) + g%hstag(i+1,:))/2
      end do

!     Update the internal energy with enthalpy fluxes
!     INSERT
      call sum_fluxes(av, flux_i, flux_j, g%area, g%roe, g%droe, g%roe_start)

!     Setup the x-momentum equation including momentum flux and pressure forces
!     INSERT
	
	do j = 1, nj-1
            flux_i(:,j) = mass_i(:,j)*(g%vx(:,j) + g%vx(:,j+1))/2 + ((g%p(:,j)+g%p(:,j+1))/2)*(g%lx_i(:,j))
      end do
      do i = 1, ni-1
            flux_j(i,:) = mass_j(i,:)*(g%vx(i,:) + g%vx(i+1,:))/2 + ((g%p(i,:)+g%p(i+1,:))/2)*(g%lx_j(i,:))
      end do
!     Update the x-momentum with momentum flux
!     INSERT
      call sum_fluxes(av, flux_i, flux_j, g%area, g%rovx, g%drovx, g%rovx_start)

!     Setup the y-momentum equation including momentum flux and pressure forces
!     INSERT
	do j = 1, nj-1
            flux_i(:,j) = mass_i(:,j)*(g%vy(:,j) + g%vy(:,j+1))/2 + ((g%p(:,j)+g%p(:,j+1))/2)*(g%ly_i(:,j))
      end do
      do i = 1, ni-1
            flux_j(i,:) = mass_j(i,:)*(g%vy(i,:) + g%vy(i+1,:))/2 + ((g%p(i,:)+g%p(i+1,:))/2)*(g%ly_j(i,:))
      end do
!     Update the y-momentum with momentum flux
!     INSERT
      call sum_fluxes(av, flux_i, flux_j, g%area, g%rovy, g%drovy, g%rovy_start)

!     Add artificial viscosity by smoothing all of the primary flow variables
      call smooth_array(av,g%ro, g%corr_ro)
      call smooth_array(av,g%roe, g%corr_roe)
      call smooth_array(av,g%rovx, g%corr_rovx)
      call smooth_array(av,g%rovy, g%corr_rovy)
!      write(6,*) '----------EULER ITERATION----------'
      end subroutine euler_iteration


