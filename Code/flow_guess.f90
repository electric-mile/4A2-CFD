      subroutine flow_guess(av,g,bcs,guesstype)

      use types
      use routines
      implicit none
      type(t_appvars), intent(inout) :: av
      type(t_grid), intent(inout) :: g
      type(t_bconds), intent(in) :: bcs
      integer, intent(in) :: guesstype
      integer :: i, j, ni, nj, j_mid
      real :: t_out, v_out, ro_out, lx, ly, l
      real :: l_i(g%ni), v_guess(g%ni), ro_guess(g%ni)
      real :: t_static(g%ni), t_lim, mach, p_static(g%ni), t_guess(g%ni), mach_guess(g%ni)
      real :: l_tot, m_f_r, a_vel, mach_lim

      ni = g%ni
      nj = g%nj

      t_out = bcs%tstag * (bcs%p_out / bcs%pstag)**av%fgam
      v_out = (2 * av%cp * (bcs%tstag - t_out))**0.5
      ro_out = bcs%p_out / (av%rgas * t_out)

      if(guesstype == 1) then

            g%ro = ro_out 
            g%roe  = g%ro * (av%cv * t_out + 0.5 * v_out**2)

            j_mid = nj / 2
            do i = 1,ni-1
                  lx = g%lx_j(i,j_mid); ly = g%ly_j(i,j_mid); 
                  l = hypot(lx,ly)
                  g%rovx(i,:) = g%ro(i,:) * v_out * ly / l
                  g%rovy(i,:) = -g%ro(i,:) * v_out * lx / l
            end do

            g%rovx(ni,:) = g%rovx(ni-1,:)
            g%rovy(ni,:) = g%rovy(ni-1,:)

            write(6,*) 'Flow Guess: Crude flow guess calculated'
            write(6,*) '  At first point ro =', g%ro(1,1), 'roe =', &
                  g%roe(1,1), 'rovx =', g%rovx(1,1), 'rovy =', g%rovy(1,1)
            write(6,*)

      else if(guesstype == 2) then 
            l_i = sum(hypot(g%lx_i, g%ly_i), 2)
            m_f_r = ro_out * v_out * l_i(ni)

            if (av%casename == 'bump') then
                  mach_lim = 1.00
                  t_lim = bcs%tstag / (1.0 + (av%gam - 1.0)/2.0 * (mach_lim**2.0))
                  v_guess = v_out
                  t_static = max(t_lim, bcs%tstag - (v_guess**2)/(2.0*av%cp))
                  ro_guess = (bcs%pstag * (t_static/bcs%tstag) ** (1.0/av%fgam))/(av%rgas * t_static)
                  v_guess = m_f_r/(ro_guess*l_i)

            else
                  mach_lim = 3.0
                  ro_guess = ro_out
                  v_guess = v_out
                  t_static = t_out

                  ! mach_lim = 3.00
                  ! t_lim = bcs%tstag / (1.0 + (av%gam - 1.0)/2.0 * (mach_lim**2.0))
                  ! v_guess = v_out
                  ! t_static = max(t_lim, bcs%tstag - (v_guess**2)/(2.0*av%cp))
                  ! ro_guess = (bcs%pstag * (t_static/bcs%tstag) ** (1.0/av%fgam))/(av%rgas * t_static)
                  ! v_guess = m_f_r/(ro_guess*l_i)

            endif

            do i = 1, ni-1

                  do j = 1, nj
                        lx = g%lx_j(i,j)
                        ly = g%ly_j(i,j)
                        l = hypot(lx,ly)
                        g%ro(i, j) = ro_guess(i)
                        g%rovx(i,j) = g%ro(i, j) * v_guess(i) * ly / l
                        g%rovy(i,j) = -g%ro(i, j) * v_guess(i) * lx / l
                        g%roe(i, j) = g%ro(i,j)*(0.50*(v_guess(i)**2.0) + (av%cv * t_static(i)))
                        if (g%ro(i,j) <0.0) then
                              write(6,*) 'Negative density at position (', i, ',', j, '): ro =', g%ro(i,j)
                              stop
                        end if
                  end do
            end do

            do j = 1, nj
                  g%rovx(ni, j) = g%rovx(ni-1,j)
                  g%rovy(ni, j) = g%rovy(ni-1,j)
                  g%roe (ni, j) = g%roe (ni-1,j)
                  g%ro (ni, j) = g%ro (ni-1,j)
            end do

            write(6,*) 'Flow Guess: Improved flow guess calculated'
            write(6,*) '  At first point ro =', g%ro(1,1), 'roe =', &
                g%roe(1,1), 'rovx =', g%rovx(1,1), 'rovy =', g%rovy(1,1)
            write(6,*)
      end if

      av%ro_ref = sum(g%ro(1,:)) / nj
      av%roe_ref = sum(g%roe(1,:)) / nj
      av%rov_ref = max(sum(g%rovx(1,:)),sum(g%rovy(1,:))) / nj

      end subroutine flow_guess
