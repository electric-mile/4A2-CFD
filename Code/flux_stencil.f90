      module flux_stencil

      contains

      subroutine sum_fluxes(av,flux_i,flux_j,area,start,prop,dcell)

      use types
      implicit none
      type(t_appvars), intent(in) :: av
      
      real, intent(inout) :: flux_i(:,:), flux_j(:,:), area(:,:)
      real, intent(inout) :: prop(:,:)
      real, intent(inout) :: start(:,:)
      real, intent(inout) :: dcell(:,:)
      real, dimension(size(prop,1),size(prop,2)) :: dnode
      real, dimension(size(dcell,1),size(dcell,2)) :: dcell_temp
      integer :: ni, nj, i, j
      ni = size(prop,1)
      nj = size(prop,2)

      if (any(area <= 0.0)) then
            write(6,*) 'Flux Stencil: Error: Cell area must be positive and non-zero.'
            stop
      end if
      dcell_temp = dcell

      dcell(1:ni-1,1:nj-1) = (flux_i(1:ni-1,:) - flux_i(2:ni,:) + flux_j(:,1:nj-1) - flux_j(:,2:nj)) * (av%dt/area) 


      dcell = (1.0 + av%facsec) * dcell - av%facsec * dcell_temp

      dnode(2:ni-1,2:nj-1) = (dcell(1:ni-2, 1:nj-2) + dcell(2:ni-1,2:nj-1) + dcell(2:ni-1,1:nj-2) + dcell(1:ni-2, 2:nj-1))/4
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
      prop = start + dnode

      dcell(1:ni-1,1:nj-1) = (flux_i(1:ni-1,:) - flux_i(2:ni,:) + flux_j(:,1:nj-1) - flux_j(:,2:nj))*(av%dt/area)

      if (any(isnan(prop))) then
            write(6,*) 'Error: NaN values detected in the prop array after updating.'
            !stop
      end if
      end subroutine sum_fluxes
      end module flux_stencil