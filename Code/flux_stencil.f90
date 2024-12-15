      module flux_stencil
      use types
      implicit none
      contains


      subroutine sum_fluxes(av,flux_i,flux_j,area,start,prop,dcell, g)

      type(t_appvars), intent(in) :: av
      
      real, intent(inout) :: flux_i(:,:), flux_j(:,:), area(:,:)
      real, intent(inout) :: prop(:,:)
      real, intent(inout) :: start(:,:)
      real, intent(inout) :: dcell(:,:)
      type(t_grid), intent(inout) :: g
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

      ! call compute_dcell_cubic(dcell, flux_i, flux_j, area, av%dt, ni, nj, g)

      dcell = (1.0 + av%facsec) * dcell - av%facsec * dcell_temp

      dnode(2:ni-1,2:nj-1) = (dcell(1:ni-2, 1:nj-2) + dcell(2:ni-1,2:nj-1) + dcell(2:ni-1,1:nj-2) + dcell(1:ni-2, 2:nj-1))/4
      !      do j = 2, nj
      !            dnode(1, j) = 0.5 * (dcell(1, j) + dcell(1, j-1))
      !            dnode(ni, j) = 0.5 * (dcell(ni-1, j) + dcell(ni-1, j-1))
      !      end do
      !      do i = 2, ni
      !            dnode(i, 1) = 0.5 * (dcell(i, 1) + dcell(i-1, 1))
      !            dnode(i, nj) = 0.5 * (dcell(i, nj-1) + dcell(i-1, nj-1))
      !      end do

      dnode(2:ni-1,[1,nj]) = (dcell(1:ni-2, [1,nj-1]) + dcell(2:ni-1,[1,nj-1]))/2
      dnode([1,ni],2:nj-1) = (dcell([1,ni-1], 1:nj-2) + dcell([1,ni-1], 2:nj-1))/2

      dnode(1, 1) = dcell(1, 1)
      dnode(1, nj) = dcell(1, nj-1)
      dnode(ni, 1) = dcell(ni-1, 1)
      dnode(ni, nj) = dcell(ni-1, nj-1)
      prop = start + dnode
      !write(6,*) 'prop', prop

      dcell(1:ni-1,1:nj-1) = (flux_i(1:ni-1,:) - flux_i(2:ni,:) + flux_j(:,1:nj-1) - flux_j(:,2:nj))*(av%dt/area)

      if (any(isnan(prop))) then
            write(6,*) 'Error: NaN values detected in the prop array after updating.'
            !stop
      end if
      end subroutine sum_fluxes
      
      subroutine fourth_order_accuracy(g, prop, prop_i, prop_j)
            real, intent(in) :: prop(:,:)
            type(t_grid), intent(in) :: g
            real, intent(out) :: prop_i(:,:), prop_j(:,:)
            integer :: ni, nj, i, j
            real :: a,b,c,d

            ni  = size(prop,1)
            nj = size(prop,2)

            do j = 1, nj
                  prop_i(1,j) = 0.5*(3*prop(1,j) - 4*prop(2,j) + prop(3,j))
                  prop_i(2,j) = (-prop(1,j) + 6*prop(2,j) - 3*prop(3,j) - 2*prop(4,j))/8
                  do i = 3, ni - 2
                        a = prop(i-2, j)
                        b = prop(i-1, j)
                        c = prop(i+1, j)
                        d = prop(i+2, j)
                        prop_i(i,j) = (-a + 7*b + 7*c - d)/12
                  end do
                  prop_i(ni-1,j) = (-prop(ni-3,j) + 6*prop(ni-2,j) + 6*prop(ni-2,j) - 3*prop(ni-1,j) + prop(ni,j))/8
                  prop_i(ni,j) = 0.5*(3*prop(ni,j) - 4*prop(ni-1,j) + prop(ni-2,j))
            end do
            do i = 1, ni
                  prop_j(i,1) = 0.5*(3*prop(i,1) - 4*prop(i,2) + prop(i,3))
                  prop_j(i,2) = (-prop(i,1) + 6*prop(i,2) - 3*prop(i,3) - 2*prop(i,4))/8
                  do j = 3, nj - 2
                        a = prop(i, j - 2)
                        b = prop(i, j - 1)
                        c = prop(i, j + 1)
                        d = prop(i, j + 2)
                        prop_j(i, j) = (-a + 7.0 * b + 7.0 * c - d) / 12.0
                  end do
                  prop_j(i, nj - 1) = (-prop(i, nj - 3) + 6.0 * prop(i, nj - 2) - 3.0 * prop(i, nj - 1) + prop(i, nj)) / 8.0
                  prop_j(i, nj) = 0.5 * (3.0 * prop(i, nj) - 4.0 * prop(i, nj - 1) + prop(i, nj - 2))
            end do
      end subroutine fourth_order_accuracy

      end module flux_stencil