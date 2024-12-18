
      module flux_stencil

!     Packaging a subroutine in a module allows it to receive the data
!     conveniently as assumed shape arrays
      
      contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine sum_fluxes(av,flux_i,flux_j,area,prop,dcell, start)

!     This subroutine sums the fluxes into each cell, calculates the change in 
!     the cell property inside, distributes the change to the four nodes of the
!     cell and then adds it onto the flow property

!     Explicitly declare the required variables
      use types
      implicit none
      type(t_appvars), intent(in) :: av
      real, intent(in) :: flux_i(:,:), flux_j(:,:), area(:,:), start(:,:)
      real, intent(inout) :: prop(:,:)

      ! Old Code
      ! real, intent(out) :: dcell(:,:)
      
      ! Crocco method
      real, intent(inout) :: dcell(:,:)
      !real, dimension(size(dcell,1),size(dcell,2)) :: dcell_temp
      real :: dcell_temp(size(dcell,1),size(dcell,2))
      !Set Facsec here
      ! 0.0 for original Lax method
      ! 0.5 for Adams-bashforth method, 2nd order accurate in time
      ! >0.5 extra 2nd derivative may increase stability or require less smoothing

      ! dnode array stores the distributed changes to the nodes
      real, dimension(size(prop,1),size(prop,2)) :: dnode
      integer :: ni, nj, i, j

      ni = size(prop,1)
      nj = size(prop,2)


!     Use the finite volume method to find the change in the variables "prop"
!     over the timestep "dt", save it in the array "dcell"
!     INSERT

      do i = 1, ni-1
            do j = 1, nj-1
                  if (area(i, j) <= 1.0e-10) then
                        write(6,*) 'flux stencil Error: Zero or very small cell area detected at (', i, ',', j, ')'
                        write(6,*) 'area size: ', area(i,j)
                        stop
                  end if
            end do
      end do

      !Crocco Method
      dcell_temp = dcell
      
      dcell = (av%dt/area)*(flux_i(1:ni-1,:) - flux_i(2:ni,:) + flux_j(:,1:nj-1) - flux_j(:,2:nj))

      !Crocco Method
      dcell = (1.0 + av%facsec) * dcell - av%facsec * dcell_temp
      dcell_temp = dcell

!     Now distribute the changes equally to the four corners of each cell. Each 
!     interior grid point receives one quarter of the change from each of the 
!     four cells adjacent to it.
!     INSERT

      dnode(2:ni-1,2:nj-1) = (dcell(2:ni-1, 2:nj-1) + dcell(2:ni-1,1:nj-2) + dcell(1:ni-2,1:nj-2) + dcell(1:ni-2, 2:nj-1))/4

!     Bounding edge nodes do not have four adjacent cells and so must be treated
!     differently, they only recieve half the change from each of the two
!     adjacent cells. Distribute the changes for the "i = 1 & ni" edges as well
!     as the "j = 1 & nj" edges. 
!     INSERT

      do j = 2, nj-1
            dnode(1, j) = 0.5 * (dcell(1, j) + dcell(1, j-1))
            dnode(ni, j) = 0.5 * (dcell(ni-1, j) + dcell(ni-1, j-1))
      end do
      do i = 2, ni-1
            dnode(i, 1) = 0.5 * (dcell(i, 1) + dcell(i-1, 1))
            dnode(i, nj) = 0.5 * (dcell(i, nj-1) + dcell(i-1, nj-1))
      end do

!     Finally distribute the changes to be to the four bounding corner points, 
!     these receive the full change from the single cell of which they form one 
!     corner.
!     INSERT
      dnode(1, 1) = dcell(1, 1)
      dnode(1, nj) = dcell(1, nj-1)
      dnode(ni, 1) = dcell(ni-1, 1)
      dnode(ni, nj) = dcell(ni-1, nj-1)

!     Update the solution by adding the changes at the nodes "dnode" to the flow
!     property "prop"
!     INSERT
      !write(6,*) dnode(:,:)
      !write(6,*) ''

      ! Runge-Kutta
      prop = start + dnode

      !prop = prop + dnode
      !write(6,*) prop(:,:)


      ! Check for NaN values in the prop array
      if (any(isnan(prop))) then
            write(6,*) 'Error: NaN values detected in the prop array after updating.'
            stop
      end if

      end subroutine sum_fluxes

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      end module flux_stencil


