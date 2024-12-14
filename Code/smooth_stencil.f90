module smooth_stencil
contains
subroutine smooth_array(av, prop, corr)
      use types
      implicit none
      type(t_appvars), intent(in) :: av
      real, intent(inout) :: prop(:,:), corr(:,:)
      real, dimension(size(prop,1),size(prop,2)) :: prop_avg, corr_total
      integer :: ni, nj, i, j
      ni = size(prop,1)
      nj = size(prop,2)
      if (ni <= 0 .or. nj <= 0) then
              write(6,*) 'Smooth Stencil Error: Grid dimensions must be positive and non-zero.'
              stop
      end if
      do i = 2, ni - 1
              do j = 2, nj - 1
                    prop_avg(i,j) = (prop(i+1,j) + prop(i,j-1) + prop(i,j+1) + prop(i-1,j))/4
              end do
      end do
      do j = 2, nj-1
              prop_avg(1,j) = (prop(1,j-1 )+prop(1,j+1)+2*prop(2,j) - prop(3,j))/3.0
              prop_avg(ni,j) = (prop(ni,j-1)+prop(ni,j+1)+2*prop(ni-1,j)-prop(ni-2,j ))/ 3.0
      end do
      do i = 2, ni-1
              prop_avg(i,1) = (prop(i-1,1)+prop(i+1,1)+2*prop(i,2) - prop(i,3))/3.0
              prop_avg(i,nj) = (prop(i-1,nj)+prop(i+1,nj)+2*prop(i,nj-1)-prop(i,nj-2))/3.0
      end do
      prop_avg([1,ni],[1,nj]) = prop([1,ni],[1,nj])
      corr_total = av%fcorr * (prop - prop_avg)
      corr = 0.99 * corr + 0.01 * corr_total
      prop = (1.0 - av%sfac) * prop + av%sfac * (prop_avg + corr) 
      if (any(isnan(prop))) then
              write(6,*) 'Smooth Stencil Error: NaN values detected in the prop array after smoothing.'
      end if
end subroutine smooth_array
end module smooth_stencil
