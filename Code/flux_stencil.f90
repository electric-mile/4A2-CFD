
      module flux_stencil
      
      contains

      subroutine sum_fluxes(av,flux_i,flux_j,area,prop,dcell, start)

      use types
      implicit none
      type(t_appvars), intent(in) :: av
      real, intent(in) :: area(:,:), start(:,:)
      real, intent(inout) :: flux_i(:,:), flux_j(:,:), prop(:,:)
      integer :: order
      ! Crocco method
      real, intent(inout) :: dcell(:,:)
      !real, dimension(size(dcell,1),size(dcell,2)) :: dcell_temp
      real :: dcell_temp(size(dcell,1),size(dcell,2))
      
      !Set Facsec here
      ! 0.0 for original Lax method
      ! 0.5 for Adams-bashforth method, 2nd order accurate in time
      ! >0.5 extra 2nd derivative may increase stability or require less smoothing
      real :: facsec = 0.0

      real, dimension(size(prop,1),size(prop,2)) :: dnode
      integer :: ni, nj, i, j

      ! Set order here
      ! 2 for second-order scheme
      ! 4 for fourth-order scheme
      order = 2

      ni = size(prop,1)
      nj = size(prop,2)

      do i = 1, ni-1
            do j = 1, nj-1
                  if (area(i, j) <= 1.0e-10) then
                        write(6,*) 'flux stencil Error: Zero or very small cell area detected at (', i, ',', j, ')'
                        write(6,*) 'area size: ', area(i,j)
                        stop
                  end if
            end do
      end do

      ! Select between second-order and fourth-order schemes
      if (order == 2) then
            call compute_flux_second_order(av, prop, flux_i, ni, nj)
      else if (order == 4) then
            call compute_flux_fourth_order(av, prop, flux_i, ni, nj)
      end if

      !Crocco Method
      dcell_temp = dcell
      
      dcell = (av%dt/area)*(flux_i(1:ni-1,:) - flux_i(2:ni,:) + flux_j(:,1:nj-1) - flux_j(:,2:nj))

      !Crocco Method
      dcell = (1.0 + facsec) * dcell - facsec * dcell_temp
      dcell_temp = dcell

      dnode(2:ni-1,2:nj-1) = (dcell(2:ni-1, 2:nj-1) + dcell(2:ni-1,1:nj-2) + dcell(1:ni-2,1:nj-2) + dcell(1:ni-2, 2:nj-1))/4

      do j = 2, nj-1
            dnode(1, j) = 0.5 * (dcell(1, j) + dcell(1, j-1))
            dnode(ni, j) = 0.5 * (dcell(ni-1, j) + dcell(ni-1, j-1))
      end do
      do i = 2, ni-1
            dnode(i, 1) = 0.5 * (dcell(i, 1) + dcell(i-1, 1))
            dnode(i, nj) = 0.5 * (dcell(i, nj-1) + dcell(i-1, nj-1))
      end do

      dnode(1, 1) = dcell(1, 1)
      dnode(1, nj) = dcell(1, nj-1)
      dnode(ni, 1) = dcell(ni-1, 1)
      dnode(ni, nj) = dcell(ni-1, nj-1)

      ! Runge-Kutta
      prop = start + dnode

      ! Old Code
      !prop = prop + dnode

      if (any(isnan(prop))) then
            write(6,*) 'Error: NaN values detected in the prop array after updating.'
            stop
      end if

      end subroutine sum_fluxes

      subroutine compute_flux_second_order(av, prop, flux, ni, nj)
      use types
      implicit none
      type(t_appvars), intent(in) :: av
      real, intent(in) :: prop(:,:)
      real, intent(out) :: flux(:,:)
      integer, intent(in) :: ni, nj
      integer :: i, j
      
      do i = 1, ni-1
            do j = 1, nj
                  flux(i,j) = 0.5 * (prop(i,j) + prop(i+1,j))
            end do
      end do
      end subroutine compute_flux_second_order

      subroutine compute_flux_fourth_order(av, prop, flux, ni, nj)
      use types
      implicit none
      type(t_appvars), intent(in) :: av
      real, intent(in) :: prop(:,:)
      real, intent(out) :: flux(:,:)
      integer, intent(in) :: ni, nj
      integer :: i, j
      real :: a, b, c, d, x
      
      do i = 3, ni-2
            do j = 1, nj
                  ! Fit a cubic polynomial using values at i-2, i-1, i+1, i+2
                  a = (-prop(i-2,j) + 3*prop(i-1,j) - 3*prop(i+1,j) + prop(i+2,j)) / 6.0
                  b = (prop(i-2,j) - 2*prop(i-1,j) + prop(i+1,j)) / 2.0
                  c = (-prop(i-2,j) + prop(i+1,j)) / 2.0
                  d = (prop(i-2,j) + 4*prop(i-1,j) + prop(i+1,j)) / 6.0
      
                  ! Evaluate the polynomial's integral over the interface to compute the flux
                  x = 0.5
                  flux(i,j) = a*x**3 + b*x**2 + c*x + d
            end do
      end do
      
      ! Handle boundaries with one-sided stencils
      do j = 1, nj
            flux(1,j) = prop(1,j)
            flux(2,j) = (prop(1,j) + prop(2,j)) / 2.0
            flux(ni-1,j) = (prop(ni-1,j) + prop(ni,j)) / 2.0
            flux(ni,j) = prop(ni,j)
      end do
      end subroutine compute_flux_fourth_order

      end module flux_stencil


