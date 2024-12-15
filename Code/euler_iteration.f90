subroutine euler_iteration(av,g)

use types
use flux_stencil
use smooth_stencil
implicit none
type(t_appvars), intent(in) :: av
type(t_grid), intent(inout) :: g
real, dimension(g%ni,g%nj-1) :: mass_i, flux_i, rovx_i, rovy_i, hstag_i, vx_i, vy_i, p_i
real, dimension(g%ni-1,g%nj) :: mass_j, flux_j, rovx_j, rovy_j, hstag_j, vx_j, vy_j, p_j
real, dimension(g%ni-1,g%nj-1) :: mach, v, t
integer :: order = 4

integer :: i, j, ni, nj
real :: dx

ni = g%ni; nj = g%nj

if (ni <= 0 .or. nj <= 0) then
      write(6,*) 'Euler Iteration Error: Grid dimensions must be positive and non-zero.'
      stop
end if
select case (order)
case (1)
      mass_i(:,1:nj-1) = ((g%rovx(:,1:nj-1)+g%rovx(:,2:nj))*g%lx_i(:,1:nj-1) + &
                    (g%rovy(:,1:nj-1)+g%rovy(:,2:nj))*g%ly_i(:,1:nj-1))/2
      mass_j(1:ni-1,:) = ((g%rovx(1:ni-1,:)+g%rovx(2:ni,:))*g%lx_j(1:ni-1,:) + &
                    (g%rovy(1:ni-1,:)+g%rovy(2:ni,:))*g%ly_j(1:ni-1,:))/2
      where(g%wall(1:ni-1,:) .and. g%wall(2:ni,:)) mass_j = 0 
      where(g%wall(:,1:nj-1) .and. g%wall(:,2:nj)) mass_i = 0 
      call sum_fluxes(av, mass_i, mass_j, g%area, g%ro_start, g%ro, g%dro, g)
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
      call sum_fluxes(av, flux_i, flux_j, g%area, g%roe_start, g%roe, g%droe, g)
      do j = 1, nj-1
            flux_i(:,j) = mass_i(:,j)*(g%vx(:,j) + g%vx(:,j+1))/2 + ((g%p(:,j)+g%p(:,j+1))/2)*(g%lx_i(:,j))
      end do
      do i = 1, ni-1
            flux_j(i,:) = mass_j(i,:)*(g%vx(i,:) + g%vx(i+1,:))/2 + ((g%p(i,:)+g%p(i+1,:))/2)*(g%lx_j(i,:))
      end do
      call sum_fluxes(av, flux_i, flux_j, g%area, g%rovx_start, g%rovx, g%drovx, g)
      do j = 1, nj-1
            flux_i(:,j) = mass_i(:,j)*(g%vy(:,j) + g%vy(:,j+1))/2 + ((g%p(:,j)+g%p(:,j+1))/2)*(g%ly_i(:,j))
      end do
      do i = 1, ni-1
            flux_j(i,:) = mass_j(i,:)*(g%vy(i,:) + g%vy(i+1,:))/2 + ((g%p(i,:)+g%p(i+1,:))/2)*(g%ly_j(i,:))
      end do
      call sum_fluxes(av, flux_i, flux_j, g%area, g%rovy_start, g%rovy, g%drovy, g)
case(2)
      mass_i(:,1:nj-1) = ((g%rovx(:,1:nj-1)+g%rovx(:,2:nj))*g%lx_i(:,1:nj-1) + &
                    (g%rovy(:,1:nj-1)+g%rovy(:,2:nj))*g%ly_i(:,1:nj-1))/2
      mass_j(1:ni-1,:) = ((g%rovx(1:ni-1,:)+g%rovx(2:ni,:))*g%lx_j(1:ni-1,:) + &
                    (g%rovy(1:ni-1,:)+g%rovy(2:ni,:))*g%ly_j(1:ni-1,:))/2
      where(g%wall(1:ni-1,:) .and. g%wall(2:ni,:)) mass_j = 0 
      where(g%wall(:,1:nj-1) .and. g%wall(:,2:nj)) mass_i = 0 
      call sum_fluxes(av, mass_i, mass_j, g%area, g%ro_start, g%ro, g%dro, g)
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
      call sum_fluxes(av, flux_i, flux_j, g%area, g%roe_start, g%roe, g%droe, g)
      do j = 1, nj-1
            flux_i(:,j) = mass_i(:,j)*(g%vx(:,j) + g%vx(:,j+1))/2 + ((g%p(:,j)+g%p(:,j+1))/2)*(g%lx_i(:,j))
      end do
      do i = 1, ni-1
            flux_j(i,:) = mass_j(i,:)*(g%vx(i,:) + g%vx(i+1,:))/2 + ((g%p(i,:)+g%p(i+1,:))/2)*(g%lx_j(i,:))
      end do
      call sum_fluxes(av, flux_i, flux_j, g%area, g%rovx_start, g%rovx, g%drovx, g)
      do j = 1, nj-1
            flux_i(:,j) = mass_i(:,j)*(g%vy(:,j) + g%vy(:,j+1))/2 + ((g%p(:,j)+g%p(:,j+1))/2)*(g%ly_i(:,j))
      end do
      do i = 1, ni-1
            flux_j(i,:) = mass_j(i,:)*(g%vy(i,:) + g%vy(i+1,:))/2 + ((g%p(i,:)+g%p(i+1,:))/2)*(g%ly_j(i,:))
      end do
      call sum_fluxes(av, flux_i, flux_j, g%area, g%rovy_start, g%rovy, g%drovy, g)
case(4)
      call fourth_order_accuracy(g, g%rovx, rovx_i, rovx_j)
      call fourth_order_accuracy(g, g%rovy, rovy_i, rovy_j)
      mass_i(:,1:nj-1) = rovx_i(:,1:nj-1)*g%lx_i(:,1:nj-1) + &
                    rovy_i(:,1:nj-1)*g%ly_i(:,1:nj-1)
      mass_j(1:ni-1,:) = rovx_j(1:ni-1,:)*g%lx_j(1:ni-1,:) + &
                    rovy_j(1:ni-1,:)*g%ly_j(1:ni-1,:)
      where(g%wall(1:ni-1,:) .and. g%wall(2:ni,:)) mass_j = 0 
      where(g%wall(:,1:nj-1) .and. g%wall(:,2:nj)) mass_i = 0 
      call sum_fluxes(av, mass_i, mass_j, g%area, g%ro_start, g%ro, g%dro, g)
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
      call sum_fluxes(av, flux_i, flux_j, g%area, g%roe_start, g%roe, g%droe, g)

      call fourth_order_accuracy(g, g%vx, vx_i, vx_j)
      call fourth_order_accuracy(g, g%p, p_i, p_j)
      do j = 1, nj-1
            flux_i(:,j) = mass_i(:,j)*(vx_i(:,j)) + (p_i(:,j))*(g%lx_i(:,j))
      end do
      do i = 1, ni-1
            flux_j(i,:) = mass_j(i,:)*(vx_j(i,:)) + (p_j(i,:))*(g%lx_j(i,:))
      end do
      call sum_fluxes(av, flux_i, flux_j, g%area, g%rovx_start, g%rovx, g%drovx, g)

      call fourth_order_accuracy(g, g%vy, vy_i, vy_j)
      call fourth_order_accuracy(g, g%p, p_i, p_j)
      do j = 1, nj-1
            flux_i(:,j) = mass_i(:,j)*(vy_i(:,j)) + (p_i(:,j))*(g%ly_i(:,j))
      end do
      do i = 1, ni-1
            flux_j(i,:) = mass_j(i,:)*(vy_j(i,:)) + (p_j(i,:))*(g%ly_j(i,:))
      end do
      call sum_fluxes(av, flux_i, flux_j, g%area, g%rovy_start, g%rovy, g%drovy, g)
end select
call smooth_array(av,g%ro, g%corr_ro)
call smooth_array(av,g%roe, g%corr_roe)
call smooth_array(av,g%rovx, g%corr_rovx)
call smooth_array(av,g%rovy, g%corr_rovy)

v = sqrt(g%rovx**2 + g%rovy**2)/g%ro
t = (g%roe - 0.5 * g%ro * (v**2)) / (g%ro * av%cv)
mach = v / ((av%gam * av%rgas * t)**0.5)

write(6,*) 'Euler Iteration: Mach Number'
write(6,*) 'roe =', g%roe(2,2),', ro =', g%ro(2,2), 'at ', av%nstep ! roe changes first
write(6,*) 'mach =', mach(2,2),', velocity =', v(2,2),', tstat =', t(2,2), 'at ', av%nstep ! tstat changes first
end subroutine euler_iteration

