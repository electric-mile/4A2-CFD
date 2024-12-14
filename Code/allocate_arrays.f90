      subroutine allocate_arrays(av,g,bcs)
      use types
      implicit none
      type(t_appvars), intent(in) :: av
      type(t_grid), intent(inout) :: g
      type(t_bconds), intent(inout) :: bcs
      integer :: ni, nj
      ni = av%ni; nj = av%nj;
      g%ni = ni; g%nj = nj;
      allocate(g%wall(ni,nj))
      allocate(bcs%ro(nj),bcs%p(nj))
      allocate(g%x(ni,nj),g%y(ni,nj))
      allocate(g%area(ni-1,nj-1),g%lx_i(ni,nj-1),g%ly_i(ni,nj-1), &
            g%lx_j(ni-1,nj),g%ly_j(ni-1,nj))
      allocate(g%ro(ni,nj),g%rovx(ni,nj),g%rovy(ni,nj),g%roe(ni,nj))
      allocate(g%corr_ro(ni, nj), g%corr_roe(ni, nj), g%corr_rovx(ni, nj), g%corr_rovy(ni, nj))
      allocate(g%dro(ni-1,nj-1),g%drovx(ni-1,nj-1), &
            g%drovy(ni-1,nj-1),g%droe(ni-1,nj-1))
      allocate(g%p(ni,nj),g%hstag(ni,nj),g%vx(ni,nj),g%vy(ni,nj))
      end subroutine allocate_arrays      
