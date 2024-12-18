      module flux_stencil
      use types
      implicit none
      contains
      subroutine sum_fluxes(av,flux_i,flux_j,area,start,prop,dcell)
      !write(6,*) "sum_fluxes"
      type(t_appvars), intent(in) :: av
      
      real, intent(inout) :: flux_i(:,:), flux_j(:,:), area(:,:)
      real, intent(inout) :: prop(:,:)
      real, intent(inout) :: start(:,:)
      real, intent(inout) :: dcell(:,:)
      !real, dimension(size(prop,1),size(prop,2)) :: dnode
      !real, dimension(size(dcell,1),size(dcell,2)) :: dcell_temp
      real, allocatable :: dcell_temp(:,:), dnode(:,:)
      

      integer :: ni, nj, i, j
      ni = av%ni
      nj = av%nj
      !write(6,*) "hi"

      allocate(dcell_temp(ni, nj))
      allocate(dnode(size(prop, 1), size(prop, 2)))

      ! Debugging statements
      !write(6,*) "sum_fluxes: ni =", ni, ", nj =", nj
      !write(6,*) "dcell size", size(dcell,1), size(dcell,2)
      !write(6,*) "dcell_temp size", size(dcell_temp,1), size(dcell_temp,2)


      if (any(area <= 0.0)) then
            write(6,*) 'Flux Stencil: Error: Cell area must be positive and non-zero.'
            stop
      end if
      dcell_temp = dcell


      dcell(1:ni-1,1:nj-1) = (flux_i(1:ni-1,:) - flux_i(2:ni,:) + flux_j(:,1:nj-1)-flux_j(:,2:nj)) * (av%dt/area)

      dcell = (1+av%facsec) * dcell - av%facsec * dcell_temp

      dnode(2:ni-1,2:nj-1) = (dcell(1:ni-2, 1:nj-2)&
    + dcell(2:ni-1,2:nj-1)&
    + dcell(2:ni-1,1:nj-2)&
    + dcell(1:ni-2,2:nj-1))/4
      !      do j = 2, nj
      !            dnode(1, j) = 0.5 * (dcell(1, j) + dcell(1, j-1))
      !            dnode(ni, j) = 0.5 * (dcell(ni-1, j) + dcell(ni-1, j-1))
      !      end do
      !      do i = 2, ni
      !            dnode(i, 1) = 0.5 * (dcell(i, 1) + dcell(i-1, 1))
      !            dnode(i, nj) = 0.5 * (dcell(i, nj-1) + dcell(i-1, nj-1))
      !      end do

      dnode(2:ni-1,[1,nj]) = 0.5*(dcell(1:ni-2,[1,nj-1])+dcell(2:ni-1,[1,nj-1]))
      dnode([1,ni],2:nj-1) = 0.5*(dcell([1,ni-1],1:nj-2)+dcell([1,ni-1],2:nj-1))

      dnode(1,1) = dcell(1,1)
      dnode(1,nj) = dcell(1, nj-1)
      dnode(ni,1) = dcell(ni-1,1)
      dnode(ni,nj) = dcell(ni-1, nj-1)
      prop = start + dnode
      !write(6,*) 'prop', prop

      dcell(1:ni-1,1:nj-1) = (flux_i(1:ni-1,:) - flux_i(2:ni,:) + flux_j(:,1:nj-1) - flux_j(:,2:nj))*(av%dt/area)

      if (any(isnan(prop))) then
            write(6,*) 'Flux Stencil Error: NaN values detected in the prop array after updating.'
            stop
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
            prop_i = 0.0
            prop_j = 0.0

            do j = 1, nj-1
                  prop_i(1,j) = (-25.0 * prop(1,j) + 48.0 * prop(2,j) - 36.0 * prop(3,j) + &
                   16.0 * prop(4,j) - 3.0 * prop(5,j)) / 12.0
                   prop_i(2,j) = (-3.0 * prop(1,j) - 10.0 * prop(2,j) + 18.0 * prop(3,j) - &
                   6.0 * prop(4,j) + prop(5,j)) / 12.0
                  do i = 3, ni - 2
                        a = prop(i-2, j)
                        b = prop(i-1, j)
                        c = prop(i+1, j)
                        d = prop(i+2, j)
                        prop_i(i,j) = (-a + 7*b + 7*c - d)/12
                  end do
                  prop_i(ni-1,j) = (3.0 * prop(ni,j) + 10.0 * prop(ni-1,j) - 18.0 * prop(ni-2,j) + &
                      6.0 * prop(ni-3,j) - prop(ni-4,j)) / 12.0
                      prop_i(ni,j) = (25.0 * prop(ni,j) - 48.0 * prop(ni-1,j) + 36.0 * prop(ni-2,j) - &
                      16.0 * prop(ni-3,j) + 3.0 * prop(ni-4,j)) / 12.0
            end do
            do i = 1, ni-1
                  prop_j(i,1) = (-25.0 * prop(i,1) + 48.0 * prop(i,2) - 36.0 * prop(i,3) + &
                   16.0 * prop(i,4) - 3.0 * prop(i,5)) / 12.0
                   prop_j(i,2) = (-3.0 * prop(i,1) - 10.0 * prop(i,2) + 18.0 * prop(i,3) - &
                   6.0 * prop(i,4) + prop(i,5)) / 12.0
                  do j = 3, nj - 2
                        a = prop(i, j - 2)
                        b = prop(i, j - 1)
                        c = prop(i, j + 1)
                        d = prop(i, j + 2)
                        prop_j(i, j) = (-a + 7.0 * b + 7.0 * c - d) / 12.0
                  end do
                  prop_j(i, nj-1) = (3.0 * prop(i, nj) + 10.0 * prop(i, nj-1) - 18.0 * prop(i, nj-2) + &
                       6.0 * prop(i, nj-3) - prop(i, nj-4)) / 12.0
                       prop_j(i, nj) = (25.0 * prop(i, nj) - 48.0 * prop(i, nj-1) + 36.0 * prop(i, nj-2) - &
                       16.0 * prop(i, nj-3) + 3.0 * prop(i, nj-4)) / 12.0
            end do
      end subroutine fourth_order_accuracy

      subroutine fourth_order_flux(g,p_var,pi_var_out,pj_var_out) ! Friends code, doesnt work with mine

            use types
            implicit none
            type(t_grid), intent(inout) :: g
            real, intent(in) :: p_var(:,:)
            real, intent(out) :: pi_var_out(:,:), pj_var_out(:,:)
            real :: a1,a2,a3,a4
            integer :: i, j, ni, nj
        
            ni = g%ni; nj = g%nj
        
            a1 = -1.0/24.0
            a2 = 13.0/24.0
            a3 = 13.0/24.0
            a4 = -1.0/24.0
        
        
            do j = 2, nj-2
                pi_var_out(:,j) = a1*p_var(:,j-1)+a2*p_var(:,j)+a3*p_var(:,j+1)+a4*p_var(:,j+2)
        
            enddo
            do i = 2, ni-2
                pj_var_out(i,:) = a1*p_var(i-1,:)+a2*p_var(i,:)+a3*p_var(i+1,:)+a4*p_var(i+2,:)
        
            enddo
        
            pi_var_out(:,[1,nj-1]) = 3.0*p_var(:,[1,nj-1]) -3.0*p_var(:,[2,nj-2]) + p_var(:,[3,nj-3])
            pj_var_out([1,ni-1],:) = 3.0*p_var([1,ni-1],:) -3.0*p_var([2,ni-2],:) + p_var([3,ni-3],:)
        
        
            end subroutine fourth_order_flux
      subroutine flux_limiter(prop, sigma_i, sigma_j)
            real, intent(in) :: prop(:,:)
            real, allocatable :: r_i(:,:), r_j(:,:)
            real, intent(out) :: sigma_i(:,:), sigma_j(:,:)
            integer :: ni, nj, i, j
            integer :: limiter = 1

            ni  = size(prop,1)
            nj = size(prop,2)
            allocate(r_i(ni,nj-1), r_j(ni-1,nj))

            r_i(3:ni,1:nj-1) = (prop(2:ni-1,1:nj-1) - prop(1:ni-2,1:nj-1)) / (prop(3:ni,1:nj-1) - prop(2:ni-1,1:nj-1))
            r_j(1:ni-1,3:nj) = (prop(1:ni-1,2:nj-1) - prop(1:ni-1,1:nj-2)) / (prop(1:ni-1,3:nj) - prop(1:ni-1,2:nj-1))
            r_i([1,2],1:nj-1) = 0
            r_j(1:ni-1,[1,2]) = 0
            select case (limiter)
            case (1) ! Minmod
                sigma_i(1:ni,1:nj-1) = max(0.0, min(1.0, r_i(1:ni,1:nj-1)))
                sigma_j(1:ni-1,1:nj) = max(0.0, min(1.0, r_j(1:ni-1,1:nj)))
            case (2) ! van-Leer
                sigma_i = (r_i + abs(r_i)) / (1.0 + abs(r_i))
                sigma_j = (r_j + abs(r_j)) / (1.0 + abs(r_j))
            case (3) ! Superbee
                sigma_i = max(0.0, min(1.0, 2.0*r_i), min(2.0, r_i))
                sigma_j = max(0.0, min(1.0, 2.0*r_j), min(2.0, r_j))
            end select
            end subroutine flux_limiter


            
      end module flux_stencil