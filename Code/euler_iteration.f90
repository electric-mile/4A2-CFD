subroutine euler_iteration(av,g)

use types
use flux_stencil
use smooth_stencil
implicit none
type(t_appvars), intent(in) :: av
type(t_grid), intent(inout) :: g
real, dimension(g%ni,g%nj-1) :: mass_i, flux_i, rovx_i, rovy_i, hstag_i, vx_i, vy_i, p_i, mass_i_1st, mass_i_4th, &
                                sigma_i_rovx, sigma_i_rovy, sigma_i_hstag, sigma_i_vx, sigma_i_p, sigma_i_vy
real, dimension(g%ni-1,g%nj) :: mass_j, flux_j, rovx_j, rovy_j, hstag_j, vx_j, vy_j, p_j, mass_j_1st, mass_j_4th, &
                                sigma_j_rovx, sigma_j_rovy, sigma_j_hstag, sigma_j_vx, sigma_j_p, sigma_j_vy
real, dimension(g%ni,g%nj) :: mach, v, t

integer :: i, j, ni, nj
real :: dx
real :: dx_i, dy_j

ni = g%ni; nj = g%nj
! Debugging statements
!print *, "euler_iteration: ni =", ni, ", nj =", nj

if (ni <= 0 .or. nj <= 0) then
      write(6,*) 'Euler Iteration Error: Grid dimensions must be positive and non-zero.'
      stop
end if
if (av%order == 1) then
            mass_i(:,1:nj-1)= (((g%rovx(:,1:nj-1) + g%rovx(:,2:nj)) * g%lx_i(:,1:nj-1)) &
                        +((g%rovy(:,1:nj-1) + g%rovy(:,2:nj)) * g%ly_i(:,1:nj-1)))/2
      

            mass_j(1:ni-1,:) = (((g%rovx(1:ni-1,:) + g%rovx(2:ni,:)) * g%lx_j(1:ni-1,:)) &
                         +((g%rovy(1:ni-1,:) + g%rovy(2:ni,:)) * g%ly_j(1:ni-1,:)))/2

            if(av%casename == 'tube') then
                  where(g%wall(1:ni-1,:) .and. g%wall(2:ni,:)) mass_j = 0 
                  where(g%wall(:,1:nj-1) .and. g%wall(:,2:nj)) mass_i = 0
                  
            else
                  where(g%wall(1:ni-1,:) .and. g%wall(2:ni,:)) mass_j = 0 
                  where(g%wall(:,1:nj-1) .and. g%wall(:,2:nj)) mass_i = 0
            endif

            call sum_fluxes(av,mass_i,mass_j,g%area,g%ro_start,g%ro,g%dro)

            flux_i(:, 1:nj-1) = mass_i(:,1:nj-1) *(g%hstag(:,1:nj-1) +g%hstag(:,2:nj))/2
            flux_j(1:ni-1,:) = mass_j(1:ni-1,:) *(g%hstag(1:ni-1, :) + g%hstag(2:ni,:))/2

            call sum_fluxes(av,flux_i,flux_j,g%area,g%roe_start,g%roe,g%droe)

            flux_i(:,1:nj-1) = mass_i(:,1:nj-1)*(g%vx(:,1:nj-1) + g%vx(:,2:nj))/2.0 &
            +((g%p(:,1:nj-1)+g%p(:,2:nj))/2.0)*(g%lx_i(:,1:nj-1))

            flux_j(1:ni-1,:) = mass_j(1:ni-1,:)*(g%vx(1:ni-1,:) + g%vx(2:ni,:))/2.0 &
            +((g%p(1:ni-1,:)+g%p(2:ni,:))/2.0)*(g%lx_j(1:ni-1,:))

            call sum_fluxes(av,flux_i,flux_j,g%area,g%rovx_start,g%rovx,g%drovx)

            flux_i(:,1:nj-1) = mass_i(:,1:nj-1)*(g%vy(:,1:nj-1) + g%vy(:,2:nj))/2.0&
            +((g%p(:,1:nj-1)+g%p(:,2:nj))/2.0)*(g%ly_i(:,1:nj-1))

            flux_j(1:ni-1,:)=mass_j(1:ni-1,:)*(g%vy(1:ni-1,:) + g%vy(2:ni,:))/2.0&
            +((g%p(1:ni-1,:)+g%p(2:ni,:))/2.0)*(g%ly_j(1:ni-1,:))

            call sum_fluxes(av,flux_i,flux_j,g%area,g%rovy_start,g%rovy,g%drovy)
   
else if (av%order == 2) then
       ! 2nd order
      mass_i(:,1:nj-1) = ((g%rovx(:,1:nj-1)+g%rovx(:,2:nj))*g%lx_i(:,1:nj-1) + &
                    (g%rovy(:,1:nj-1)+g%rovy(:,2:nj))*g%ly_i(:,1:nj-1))/2
      mass_j(1:ni-1,:) = ((g%rovx(1:ni-1,:)+g%rovx(2:ni,:))*g%lx_j(1:ni-1,:) + &
                    (g%rovy(1:ni-1,:)+g%rovy(2:ni,:))*g%ly_j(1:ni-1,:))/2
      where(g%wall(1:ni-1,:) .and. g%wall(2:ni,:)) mass_j = 0 
      where(g%wall(:,1:nj-1) .and. g%wall(:,2:nj)) mass_i = 0 
      call sum_fluxes(av, mass_i, mass_j, g%area, g%ro_start, g%ro, g%dro)
      if (any(g%ro <= 0.0)) then
            write(6,*) 'Euler Iteration Error: Non-positive density values detected.'
            do i = 1, ni
            do j = 1, nj
                  if (g%ro(i,j) <= 0.0) then
                        write(6,*) 'Non-positive density at (', i, ',', j, '): ro =', g%ro(i,j)
                  end if
            end do
            end do
            stop
      end if
      do j = 1, nj-1
            flux_i(:,j) = mass_i(:,j)*(g%hstag(:,j) + g%hstag(:,j+1))/2
      end do
      do i = 1, ni-1
            flux_j(i,:) = mass_j(i,:)*(g%hstag(i,:) + g%hstag(i+1,:))/2
      end do
      call sum_fluxes(av, flux_i, flux_j, g%area, g%roe_start, g%roe, g%droe)
      do j = 1, nj-1
            flux_i(:,j) = mass_i(:,j)*(g%vx(:,j) + g%vx(:,j+1))/2 + ((g%p(:,j)+g%p(:,j+1))/2)*(g%lx_i(:,j))
      end do
      do i = 1, ni-1
            flux_j(i,:) = mass_j(i,:)*(g%vx(i,:) + g%vx(i+1,:))/2 + ((g%p(i,:)+g%p(i+1,:))/2)*(g%lx_j(i,:))
      end do
      call sum_fluxes(av, flux_i, flux_j, g%area, g%rovx_start, g%rovx, g%drovx)
      do j = 1, nj-1
            flux_i(:,j) = mass_i(:,j)*(g%vy(:,j) + g%vy(:,j+1))/2 + ((g%p(:,j)+g%p(:,j+1))/2)*(g%ly_i(:,j))
      end do
      do i = 1, ni-1
            flux_j(i,:) = mass_j(i,:)*(g%vy(i,:) + g%vy(i+1,:))/2 + ((g%p(i,:)+g%p(i+1,:))/2)*(g%ly_j(i,:))
      end do
      call sum_fluxes(av, flux_i, flux_j, g%area, g%rovy_start, g%rovy, g%drovy)
else if (av%order == 4) then ! 4th order, no flux_lim

      call fourth_order_accuracy(g, g%rovx, rovx_i, rovx_j)

      call fourth_order_accuracy(g, g%rovy, rovy_i, rovy_j)

      mass_i(:,1:nj-1) = rovx_i(:,1:nj-1)*g%lx_i(:,1:nj-1) + &
                    rovy_i(:,1:nj-1)*g%ly_i(:,1:nj-1)

      mass_j(1:ni-1,:) = rovx_j(1:ni-1,:)*g%lx_j(1:ni-1,:) + &
                    rovy_j(1:ni-1,:)*g%ly_j(1:ni-1,:)

      where(g%wall(1:ni-1,:) .and. g%wall(2:ni,:)) mass_j = 0 
      where(g%wall(:,1:nj-1) .and. g%wall(:,2:nj)) mass_i = 0 
      call sum_fluxes(av, mass_i, mass_j, g%area, g%ro_start, g%ro, g%dro)

      if (any(g%ro <= 0.0)) then
            write(6,*) 'Euler Iteration Error: Non-positive density values detected.'
            do i = 1, ni
            do j = 1, nj
                  if (g%ro(i,j) <= 0.0) then
                        write(6,*) 'Non-positive density at (', i, ',', j, '): ro =', g%ro(i,j)
                  end if
            end do
            end do
            stop
      end if

      call fourth_order_accuracy(g, g%hstag, hstag_i, hstag_j)

      do j = 1, nj-1
            flux_i(:,j) = mass_i(:,j) * hstag_i(:,j)
      end do

      do i = 1, ni-1
            flux_j(i,:) = mass_j(i,:)*hstag_j(i,:)
      end do
      call sum_fluxes(av, flux_i, flux_j, g%area, g%roe_start, g%roe, g%droe)

      call fourth_order_accuracy(g, g%vx, vx_i, vx_j)

      call fourth_order_accuracy(g, g%p, p_i, p_j)

      do j = 1, nj-1
            flux_i(:,j) = mass_i(:,j)*(vx_i(:,j)) + (p_i(:,j))*(g%lx_i(:,j))
      end do

      do i = 1, ni-1
            flux_j(i,:) = mass_j(i,:)*(vx_j(i,:)) + (p_j(i,:))*(g%lx_j(i,:))
      end do
      call sum_fluxes(av, flux_i, flux_j, g%area, g%rovx_start, g%rovx, g%drovx)
      !call fourth_order_accuracy(g, g%vy, vy_i, vy_j)
      !call fourth_order_accuracy(g, g%p, p_i, p_j)

!      do j = 1, nj-1
!            flux_i(:,j) = mass_i(:,j)*(vy_i(:,j)) + (p_i(:,j))*(g%ly_i(:,j))
!     end do
!
!      do i = 1, ni-1
!            flux_j(i,:) = mass_j(i,:)*(vy_j(i,:)) + (p_j(i,:))*(g%ly_j(i,:))
!      end do
      !write(6,*) "drovy size", size(g%drovy,1), size(g%drovy,2)

      do j = 1, nj-1
            flux_i(:,j) = mass_i(:,j)*(g%vy(:,j) + g%vy(:,j+1))/2 + ((g%p(:,j)+g%p(:,j+1))/2)*(g%ly_i(:,j))
      end do
      do i = 1, ni-1
            flux_j(i,:) = mass_j(i,:)*(g%vy(i,:) + g%vy(i+1,:))/2 + ((g%p(i,:)+g%p(i+1,:))/2)*(g%ly_j(i,:))
      end do

      call sum_fluxes(av, flux_i, flux_j, g%area, g%rovy_start, g%rovy, g%drovy) ! This is the one it fails in
      !write(6,*) "bash"
else if (av%order == 5) then  ! Friends code, 4th order, no flux_lim
      
      call fourth_order_flux(g,g%rovx,rovx_i,rovx_j)
      call fourth_order_flux (g,g%rovy,rovy_i,rovy_j)
      
      mass_i(:,1:nj-1) = rovx_i(:,1:nj-1)*g%lx_i(:,1:nj-1) + rovy_i(:,1:nj-1) * g%ly_i(:,1:nj-1)
      mass_j(1:ni-1,:) = rovx_j(1:ni-1,:)* g%lx_j(1:ni-1,:) + rovy_j(1:ni-1,:) * g%ly_j(1:ni-1,:)
      
      if(av%casename == 'tube') then
            where(g%wall(1:ni-1,:) .and. g%wall(2:ni,:)) mass_j = 0 
            where(g%wall(:,1:nj-1) .and. g%wall(:,2:nj)) mass_i = 0
!                 added 0 mass flux through the ends of the tube, assuming only occurs in the i direction                  
      else
            where(g%wall(1:ni-1,:) .and. g%wall(2:ni,:)) mass_j = 0 
            where(g%wall(:,1:nj-1) .and. g%wall(:,2:nj)) mass_i = 0
      endif

      call sum_fluxes(av,mass_i,mass_j,g%area,g%ro_start,g%ro,g%dro)

      call fourth_order_flux(g,g%hstag,hstag_i,hstag_j)

      flux_i(:, 1:nj-1) = mass_i(:,1:nj-1) * hstag_i(:,1:nj-1)
      flux_j(1:ni-1,:) = mass_j(1:ni-1,:) * hstag_j(1:ni-1,:)

      call sum_fluxes(av,flux_i,flux_j,g%area,g%roe_start,g%roe,g%droe)

      call fourth_order_flux(g,g%vx,vx_i,vx_j)
      call fourth_order_flux(g,g%p,p_i,p_j)


      flux_i(:,1:nj-1) = mass_i(:,1:nj-1)*(vx_i(:,1:nj-1))&
      +(p_i(:,1:nj-1))*(g%lx_i(:,1:nj-1))

      flux_j(1:ni-1,:) = mass_j(1:ni-1,:)*(vx_j(1:ni-1,:))&
      +(p_j(1:ni-1,:))*(g%lx_j(1:ni-1,:))

      call sum_fluxes(av,flux_i,flux_j,g%area,g%rovx_start,g%rovx,g%drovx)

      call fourth_order_flux(g,g%vy,vy_i,vy_j)
      !call fourth_order_flux(g,g%p,p_i,p_j)

      flux_i(:,1:nj-1) = mass_i(:,1:nj-1)*(vy_i(:,1:nj-1))&
      +(p_i(:,1:nj-1))*(g%ly_i(:,1:nj-1))

      flux_j(1:ni-1,:) = mass_j(1:ni-1,:)*(vy_j(1:ni-1,:))&
      +(p_j(1:ni-1,:))*(g%ly_j(1:ni-1,:))

      call sum_fluxes(av,flux_i,flux_j,g%area,g%rovy_start,g%rovy,g%drovy)
else if (av%order == 6) then  ! Friends code, 4th order, with flux_lim
      call fourth_order_flux(g,g%rovx,rovx_i,rovx_j)
      call fourth_order_flux (g,g%rovy,rovy_i,rovy_j)
      call flux_limiter(g%rovx,sigma_i_rovx,sigma_j_rovx)
      call flux_limiter(g%rovy,sigma_i_rovy,sigma_j_rovy)

      mass_i(2:ni,1:nj-1) = (g%rovx(1:ni-1,1:nj-1) + sigma_i_rovx(2:ni,1:nj-1) * &
                             (rovx_i(2:ni,1:nj-1) - g%rovx(1:ni-1,1:nj-1))) * &
                             g%lx_i(2:ni,1:nj-1) + &
                             (g%rovy(1:ni-1,1:nj-1) + sigma_i_rovy(2:ni,1:nj-1) * &
                             (rovy_i(2:ni,1:nj-1) - g%rovy(1:ni-1,1:nj-1))) * &
                             g%ly_i(2:ni,1:nj-1)
      mass_j(1:ni-1,2:nj) = (g%rovx(1:ni-1,1:nj-1) + sigma_j_rovx(1:ni-1,2:nj) * &
                             (rovx_j(1:ni-1,2:nj) - g%rovx(1:ni-1,1:nj-1))) * &
                             g%lx_j(1:ni-1,2:nj) + &
                             (g%rovy(1:ni-1,1:nj-1) + sigma_j_rovy(1:ni-1,2:nj) * &
                             (rovy_j(1:ni-1,2:nj) - g%rovy(1:ni-1,1:nj-1))) * &
                             g%ly_j(1:ni-1,2:nj)

      mass_i(1,1:nj-1) = mass_i(2,1:nj-1)
      mass_j(1:ni-1,1) = mass_j(1:ni-1,2)

      mass_i_1st(2:ni,1:nj-1) = g%rovx(1:ni-1,1:nj-1)*g%lx_i(2:ni,1:nj-1) + g%rovy(1:ni-1,1:nj-1) * g%ly_i(2:ni,1:nj-1)
      mass_j_1st(1:ni-1,2:nj) = g%rovx(1:ni-1,1:nj-1)* g%lx_j(1:ni-1,2:nj) + g%rovy(1:ni-1,1:nj-1) * g%ly_j(1:ni-1,2:nj)
      mass_i_1st(1,1:nj-1) = mass_i_1st(2,1:nj-1)
      mass_j_1st(1:ni-1,1) = mass_j_1st(1:ni-1,2)

      mass_i_4th(:,1:nj-1) = rovx_i(:,1:nj-1)*g%lx_i(:,1:nj-1) + rovy_i(:,1:nj-1) * g%ly_i(:,1:nj-1)
      mass_j_4th(1:ni-1,:) = rovx_j(1:ni-1,:)* g%lx_j(1:ni-1,:) + rovy_j(1:ni-1,:) * g%ly_j(1:ni-1,:)
      
      if(av%casename == 'tube') then
            where(g%wall(1:ni-1,:) .and. g%wall(2:ni,:)) mass_j = 0 
            where(g%wall(:,1:nj-1) .and. g%wall(:,2:nj)) mass_i = 0
!                 added 0 mass flux through the ends of the tube, assuming only occurs in the i direction                  
      else
            where(g%wall(1:ni-1,:) .and. g%wall(2:ni,:)) mass_j = 0 
            where(g%wall(:,1:nj-1) .and. g%wall(:,2:nj)) mass_i = 0
      endif

      call sum_fluxes(av,mass_i,mass_j,g%area,g%ro_start,g%ro,g%dro)

      call fourth_order_flux(g,g%hstag,hstag_i,hstag_j)
      call flux_limiter(g%hstag,sigma_i_hstag,sigma_j_hstag)

      flux_i(2:ni,1:nj-1) = (mass_i_1st(2:ni,1:nj-1) * g%hstag(1:ni-1,1:nj-1)*(1.0 - sigma_i_hstag(2:ni,1:nj-1))) + &
                            (mass_i_4th(2:ni,1:nj-1)* sigma_i_hstag(2:ni,1:nj-1) * hstag_i(2:ni,1:nj-1))
      flux_j(1:ni-1,2:nj) = (mass_j_1st(1:ni-1,2:nj) * g%hstag(1:ni-1,1:nj-1)*(1.0 - sigma_j_hstag(1:ni-1,2:nj))) + &
                              (mass_j_4th(1:ni-1,2:nj)* sigma_j_hstag(1:ni-1,2:nj) * hstag_j(1:ni-1,2:nj))

      flux_i(1,1:nj-1) = flux_i(2,1:nj-1)
      flux_j(1:ni-1,1) = flux_j(1:ni-1,2)

      !flux_i(:, 1:nj-1) = mass_i(:,1:nj-1) * hstag_i(:,1:nj-1)
      !flux_j(1:ni-1,:) = mass_j(1:ni-1,:) * hstag_j(1:ni-1,:)

      call sum_fluxes(av,flux_i,flux_j,g%area,g%roe_start,g%roe,g%droe)

      call fourth_order_flux(g,g%vx,vx_i,vx_j)
      call fourth_order_flux(g,g%p,p_i,p_j)
      call flux_limiter(g%vx,sigma_i_vx,sigma_j_vx)
      call flux_limiter(g%p,sigma_i_p,sigma_j_p)

      flux_i(2:ni,1:nj-1) = (mass_i_1st(2:ni,1:nj-1) * g%vx(1:ni-1,1:nj-1)*(1.0-sigma_i_vx(2:ni,1:nj-1))) + &
                              (mass_i_4th(2:ni,1:nj-1) * sigma_i_vx(2:ni,1:nj-1) * vx_i(2:ni,1:nj-1)) + &
                            (g%p(1:ni-1,1:nj-1) + sigma_i_p(2:ni,1:nj-1) * &
                            (p_i(2:ni,1:nj-1) - g%p(1:ni-1,1:nj-1))) * &
                            g%lx_i(2:ni,1:nj-1)
      flux_j(1:ni-1,2:nj) = (mass_j_1st(1:ni-1,2:nj) * g%vx(1:ni-1,1:nj-1)*(1.0-sigma_j_vx(1:ni-1,2:nj))) + &
                              (mass_j_4th(1:ni-1,2:nj) * sigma_j_vx(1:ni-1,2:nj) * vx_j(1:ni-1,2:nj)) + &
                            (g%p(1:ni-1,1:nj-1) + sigma_j_p(1:ni-1,2:nj) * &
                            (p_j(1:ni-1,2:nj) - g%p(1:ni-1,1:nj-1))) * &
                            g%lx_j(1:ni-1,2:nj)
      flux_i(1,1:nj-1) = flux_i(2,1:nj-1)
      flux_j(1:ni-1,1) = flux_j(1:ni-1,2)

      !flux_i(:,1:nj-1) = mass_i(:,1:nj-1)*(vx_i(:,1:nj-1))&
      !+(p_i(:,1:nj-1))*(g%lx_i(:,1:nj-1))

      !flux_j(1:ni-1,:) = mass_j(1:ni-1,:)*(vx_j(1:ni-1,:))&
      !+(p_j(1:ni-1,:))*(g%lx_j(1:ni-1,:))

      call sum_fluxes(av,flux_i,flux_j,g%area,g%rovx_start,g%rovx,g%drovx)

      call fourth_order_flux(g,g%vy,vy_i,vy_j)
      call flux_limiter(g%vy,sigma_i_vy,sigma_j_vy)
      
      flux_i(2:ni,1:nj-1) = (mass_i_1st(2:ni,1:nj-1) * g%vy(1:ni-1,1:nj-1) * (1.0 -sigma_i_vy(2:ni,1:nj-1))) + &
                              (mass_i_4th(2:ni,1:nj-1) * sigma_i_vy(2:ni,1:nj-1) * vy_i(2:ni,1:nj-1)) + &
                            (g%p(1:ni-1,1:nj-1) + sigma_i_p(2:ni,1:nj-1) * &
                            (p_i(2:ni,1:nj-1) - g%p(1:ni-1,1:nj-1))) * &
                            g%ly_i(2:ni,1:nj-1)

      flux_j(1:ni-1,2:nj) = (mass_j_1st(1:ni-1,2:nj) * g%vy(1:ni-1,1:nj-1) * (1.0 -sigma_j_vy(1:ni-1,2:nj))) + &
                              (mass_j_4th(1:ni-1,2:nj) * sigma_j_vy(1:ni-1,2:nj) * vy_j(1:ni-1,2:nj)) + &
                              (g%p(1:ni-1,1:nj-1) + sigma_j_p(1:ni-1,2:nj) * &
                              (p_j(1:ni-1,2:nj) - g%p(1:ni-1,1:nj-1))) * &
                              g%ly_j(1:ni-1,2:nj)
      flux_i(1,1:nj-1) = flux_i(2,1:nj-1)
      flux_j(1:ni-1,1) = flux_j(1:ni-1,2)

      !flux_i(:,1:nj-1) = mass_i(:,1:nj-1)*(vy_i(:,1:nj-1))&
      !+(p_i(:,1:nj-1))*(g%ly_i(:,1:nj-1))

      !flux_j(1:ni-1,:) = mass_j(1:ni-1,:)*(vy_j(1:ni-1,:))&
      !+(p_j(1:ni-1,:))*(g%ly_j(1:ni-1,:))

      call sum_fluxes(av,flux_i,flux_j,g%area,g%rovy_start,g%rovy,g%drovy)
else if (av%order == 7) then  !4th order, with flux limiter
      call fourth_order_accuracy(g,g%rovx,rovx_i,rovx_j)
      call fourth_order_accuracy(g,g%rovy,rovy_i,rovy_j)
      call flux_limiter(g%rovx,sigma_i_rovx,sigma_j_rovx)
      call flux_limiter(g%rovy,sigma_i_rovy,sigma_j_rovy)

      mass_i(2:ni,1:nj-1) = (g%rovx(1:ni-1,1:nj-1) + sigma_i_rovx(2:ni,1:nj-1) * &
                             (rovx_i(2:ni,1:nj-1) - g%rovx(1:ni-1,1:nj-1))) * &
                             g%lx_i(2:ni,1:nj-1) + &
                             (g%rovy(1:ni-1,1:nj-1) + sigma_i_rovy(2:ni,1:nj-1) * &
                             (rovy_i(2:ni,1:nj-1) - g%rovy(1:ni-1,1:nj-1))) * &
                             g%ly_i(2:ni,1:nj-1)
      mass_i(1,1:nj-1) = mass_i(2,1:nj-1)
      mass_j(1:ni-1,2:nj) = (g%rovx(1:ni-1,1:nj-1) + sigma_j_rovx(1:ni-1,2:nj) * &
                             (rovx_j(1:ni-1,2:nj) - g%rovx(1:ni-1,1:nj-1))) * &
                             g%lx_j(1:ni-1,2:nj) + &
                             (g%rovy(1:ni-1,1:nj-1) + sigma_j_rovy(1:ni-1,2:nj) * &
                             (rovy_j(1:ni-1,2:nj) - g%rovy(1:ni-1,1:nj-1))) * &
                             g%ly_j(1:ni-1,2:nj)
      mass_j(1:ni-1,1) = mass_j(1:ni-1,2)

      mass_i_1st(2:ni,1:nj-1) = g%rovx(1:ni-1,1:nj-1)*g%lx_i(2:ni,1:nj-1) + g%rovy(1:ni-1,1:nj-1) * g%ly_i(2:ni,1:nj-1)
      mass_j_1st(1:ni-1,2:nj) = g%rovx(1:ni-1,1:nj-1)* g%lx_j(1:ni-1,2:nj) + g%rovy(1:ni-1,1:nj-1) * g%ly_j(1:ni-1,2:nj)
      mass_i_1st(1,1:nj-1) = mass_i_1st(2,1:nj-1)
      mass_j_1st(1:ni-1,1) = mass_j_1st(1:ni-1,2)
      
      if(av%casename == 'tube') then
            where(g%wall(1:ni-1,:) .and. g%wall(2:ni,:)) mass_j = 0 
            where(g%wall(:,1:nj-1) .and. g%wall(:,2:nj)) mass_i = 0
!                 added 0 mass flux through the ends of the tube, assuming only occurs in the i direction                  
      else
            where(g%wall(1:ni-1,:) .and. g%wall(2:ni,:)) mass_j = 0 
            where(g%wall(:,1:nj-1) .and. g%wall(:,2:nj)) mass_i = 0
      endif

      call sum_fluxes(av,mass_i,mass_j,g%area,g%ro_start,g%ro,g%dro)

      call fourth_order_accuracy(g,g%hstag,hstag_i,hstag_j)
      call flux_limiter(g%hstag,sigma_i_hstag,sigma_j_hstag)

      flux_i(2:ni,1:nj-1) = (mass_i_1st(2:ni,1:nj-1) * g%hstag(1:ni-1,1:nj-1)*(1.0 - sigma_i_hstag(2:ni,1:nj-1))) + &
                            (mass_i_4th(2:ni,1:nj-1)* sigma_i_hstag(2:ni,1:nj-1) * hstag_i(2:ni,1:nj-1))
      flux_j(1:ni-1,2:nj) = (mass_j_1st(1:ni-1,2:nj) * g%hstag(1:ni-1,1:nj-1)*(1.0 - sigma_j_hstag(1:ni-1,2:nj))) + &
                              (mass_j_4th(1:ni-1,2:nj)* sigma_j_hstag(1:ni-1,2:nj) * hstag_j(1:ni-1,2:nj))

      !flux_i(:, 1:nj-1) = mass_i(:,1:nj-1) * hstag_i(:,1:nj-1)
      !flux_j(1:ni-1,:) = mass_j(1:ni-1,:) * hstag_j(1:ni-1,:)

      call sum_fluxes(av,flux_i,flux_j,g%area,g%roe_start,g%roe,g%droe)

      call fourth_order_accuracy(g,g%vx,vx_i,vx_j)
      call fourth_order_accuracy(g,g%p,p_i,p_j)
      call flux_limiter(g%vx,sigma_i_vx,sigma_j_vx)
      call flux_limiter(g%p,sigma_i_p,sigma_j_p)

      flux_i(2:ni,1:nj-1) = (mass_i_1st(2:ni,1:nj-1) * g%vx(1:ni-1,1:nj-1)*(1.0-sigma_i_vx(2:ni,1:nj-1))) + &
                              (mass_i_4th(2:ni,1:nj-1) * sigma_i_vx(2:ni,1:nj-1) * vx_i(2:ni,1:nj-1)) + &
                            (g%p(1:ni-1,1:nj-1) + sigma_i_p(2:ni,1:nj-1) * &
                            (p_i(2:ni,1:nj-1) - g%p(1:ni-1,1:nj-1))) * &
                            g%lx_i(2:ni,1:nj-1)
      flux_j(1:ni-1,2:nj) = (mass_j_1st(1:ni-1,2:nj) * g%vx(1:ni-1,1:nj-1)*(1.0-sigma_j_vx(1:ni-1,2:nj))) + &
                              (mass_j_4th(1:ni-1,2:nj) * sigma_j_vx(1:ni-1,2:nj) * vx_j(1:ni-1,2:nj)) + &
                            (g%p(1:ni-1,1:nj-1) + sigma_j_p(1:ni-1,2:nj) * &
                            (p_j(1:ni-1,2:nj) - g%p(1:ni-1,1:nj-1))) * &
                            g%lx_j(1:ni-1,2:nj)

      flux_i(1,1:nj-1) = flux_i(2,1:nj-1)
      flux_j(1:ni-1,1) = flux_j(1:ni-1,2)

      !flux_i(:,1:nj-1) = mass_i(:,1:nj-1)*(vx_i(:,1:nj-1))&
      !+(p_i(:,1:nj-1))*(g%lx_i(:,1:nj-1))

      !flux_j(1:ni-1,:) = mass_j(1:ni-1,:)*(vx_j(1:ni-1,:))&
      !+(p_j(1:ni-1,:))*(g%lx_j(1:ni-1,:))

      call sum_fluxes(av,flux_i,flux_j,g%area,g%rovx_start,g%rovx,g%drovx)

      call fourth_order_accuracy(g,g%vy,vy_i,vy_j)
      call flux_limiter(g%vy,sigma_i_vy,sigma_j_vy)
      
      flux_i(2:ni,1:nj-1) = (mass_i_1st(2:ni,1:nj-1) * g%vy(1:ni-1,1:nj-1) * (1.0 -sigma_i_vy(2:ni,1:nj-1))) + &
                              (mass_i_4th(2:ni,1:nj-1) * sigma_i_vy(2:ni,1:nj-1) * vy_i(2:ni,1:nj-1)) + &
                            (g%p(1:ni-1,1:nj-1) + sigma_i_p(2:ni,1:nj-1) * &
                            (p_i(2:ni,1:nj-1) - g%p(1:ni-1,1:nj-1))) * &
                            g%ly_i(2:ni,1:nj-1)

      flux_j(1:ni-1,2:nj) = (mass_j_1st(1:ni-1,2:nj) * g%vy(1:ni-1,1:nj-1) * (1.0 -sigma_j_vy(1:ni-1,2:nj))) + &
                              (mass_j_4th(1:ni-1,2:nj) * sigma_j_vy(1:ni-1,2:nj) * vy_j(1:ni-1,2:nj)) + &
                              (g%p(1:ni-1,1:nj-1) + sigma_j_p(1:ni-1,2:nj) * &
                              (p_j(1:ni-1,2:nj) - g%p(1:ni-1,1:nj-1))) * &
                              g%ly_j(1:ni-1,2:nj)
      flux_i(1,1:nj-1) = flux_i(2,1:nj-1)
      flux_j(1:ni-1,1) = flux_j(1:ni-1,2)

      !flux_i(:,1:nj-1) = mass_i(:,1:nj-1)*(vy_i(:,1:nj-1))&
      !+(p_i(:,1:nj-1))*(g%ly_i(:,1:nj-1))

      !flux_j(1:ni-1,:) = mass_j(1:ni-1,:)*(vy_j(1:ni-1,:))&
      !+(p_j(1:ni-1,:))*(g%ly_j(1:ni-1,:))

      call sum_fluxes(av,flux_i,flux_j,g%area,g%rovy_start,g%rovy,g%drovy)
end if 
call smooth_array(av,g%ro, g%corr_ro)
call smooth_array(av,g%roe, g%corr_roe)
call smooth_array(av,g%rovx, g%corr_rovx)
call smooth_array(av,g%rovy, g%corr_rovy)


if (av%casename /= 'tube') then
      v = sqrt(g%rovx**2 + g%rovy**2)/g%ro
      t = (g%roe - 0.5 * g%ro * (v**2)) / (g%ro * av%cv)
      mach = v / ((av%gam * av%rgas * t)**0.5)
      write(6,*) 'Euler Iteration: Mach Number'
      write(6,*) 'roe =', g%roe(2,2),', ro =', g%ro(2,2), 'at ', av%nstep ! roe changes first
      write(6,*) 'mach =', minval(mach(2,:)), maxval(mach(2,:)),', velocity =', v(2,2),', tstat =', t(2,2), 'at ', av%nstep ! tstat changes first
end if
end subroutine euler_iteration

