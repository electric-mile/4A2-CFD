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
      real :: l_i(g%ni), v_guess(g%ni), ro_guess(g%ni), l_temp(g%nj)
      real :: t_static(g%ni), t_lim, mach, p_static(g%ni), t_guess(g%ni), mach_guess(g%ni)
      real :: l_tot, m_f_r, a_vel, mach_lim

      ni = g%ni
      nj = g%nj

      t_out = (bcs%tstag * (bcs%p_out / bcs%pstag)**av%fgam) 
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
            do i = 1 , ni
                  do j = 1, nj - 1
                        l_temp(j) = hypot(g%lx_i(i,j), g%ly_i(i,j))
                  end do
                  l_i(i) = sum(l_temp)  
            end do
            m_f_r = ro_out * v_out * l_i(ni)

            if (av%casename == 'bump') then
                  mach_lim = 1.00
                  t_lim = bcs%tstag / (1.0 + (av%gam - 1.0)/2.0 * (mach_lim**2.0))
                  v_guess = m_f_r/(ro_out*l_i)
                  t_static = max(t_lim, bcs%tstag - (v_guess**2)/(2.0*av%cp))
                  ro_guess = (bcs%pstag * (t_static/bcs%tstag) ** (1.0/av%fgam))/(av%rgas * t_static)
                  v_guess = m_f_r/(ro_guess*l_i)
            else if (av%casename == 'bend') then
                  mach_lim = 1.0
                  t_lim = bcs%tstag / (1.0 + (av%gam - 1.0)/2.0 * (mach_lim**2.0))
                  v_guess = m_f_r/(ro_out*l_i)
                  t_static = max(t_lim, bcs%tstag - (v_guess**2)/(2.0*av%cp))
                  ro_guess = (bcs%pstag * (t_static/bcs%tstag) ** (1.0/av%fgam))/(av%rgas * t_static)
                  v_guess = m_f_r/(ro_guess*l_i)
                  ! else if (av%casename == 'tube') then
                  !       ro_guess(1:(ni+1)/2) = 1.0
                  !       t_static(1:(ni+1)/2) = (bcs%tstag * (1.0/bcs%pstag)**av%fgam)
                  !       v_guess = 0.0
                  !       ro_guess((ni+2)/2:ni) = 0.125
                  !       t_static((ni+2)/2:ni) = (bcs%tstag * (0.1/bcs%pstag)**av%fgam)      
            else
                  mach_lim = 3.0
                  ro_guess = ro_out
                  v_guess = v_out
                  t_static = t_out
            endif

            write(6,*) 'Inlet velocity is:', v_guess(1),'m/s'
            write(6,*) 'Inlet Mach No:', v_guess(1)/(sqrt(av%gam*av%rgas*T_static(1)))

            do i = 1, ni -1 
                  do j = 1, nj  
                        lx = g%lx_j(i,j); ly = g%ly_j(i,j); 
                        l = hypot(lx,ly)
                        g%ro(i,j) = ro_guess(i)
                        g%rovx(i,j) = g%ro(i,j) * v_guess(i) * (ly / l)
                        g%rovy(i,j) = -g%ro(i,j) * v_guess(i) * (lx/ l)
                        g%roe(i,j) =  ro_guess(i) * (((v_guess(i)**2)/2.0) + (av%cv * t_static(i)))
                  end do
            end do

            g%ro(ni,:) = g%ro(ni-1,:)
            g%roe(ni,:) = g%roe(ni-1,:)
            g%rovx(ni,:) = g%rovx(ni-1,:)
            g%rovy(ni,:) = g%rovy(ni-1,:)

            write(6,*) 'Flow Guess: Improved flow guess calculated'
            write(6,*) '  At first point ro =', g%ro(1,1), 'roe =', &
                g%roe(1,1), 'rovx =', g%rovx(1,1), 'rovy =', g%rovy(1,1)
            write(6,*)
      end if

      av%ro_ref = sum(g%ro(1,:)) / nj
      av%roe_ref = sum(g%roe(1,:)) / nj
      av%rov_ref = max(sum(g%rovx(1,:)),sum(g%rovy(1,:))) / nj

      end subroutine flow_guess
