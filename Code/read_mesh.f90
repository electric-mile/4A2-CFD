subroutine read_mesh(av,g,bcs,p)
use types
implicit none
type(t_appvars), intent(inout) :: av
type(t_grid), allocatable, intent(out) :: g(:)
type(t_bconds), intent(inout) :: bcs
type(t_match), allocatable, intent(out) :: p(:)
integer :: n, ni, nj, m
open(2,file='mesh_' // av%casename // '.bin',form='unformatted', &
      access='stream')
read(2) av%nn; allocate(g(av%nn));
write(6,*) 'Read multi-block mesh from file'
write(6,*) '  Number of blocks nn =', av%nn
do n = 1,av%nn
      read(2) ni, nj
      g(n)%ni = ni; g(n)%nj = nj;
      write(6,*) '  Size of block', n, 'ni =', ni, 'nj =', nj
      allocate(g(n)%x(ni,nj),g(n)%y(ni,nj))
      read(2) g(n)%x; read(2) g(n)%y;
      allocate(g(n)%wall(ni,nj))
      read(2) g(n)%wall;
      allocate(g(n)%area(ni-1,nj-1),g(n)%lx_i(ni,nj-1),g(n)%ly_i(ni,nj-1), &
            g(n)%lx_j(ni-1,nj),g(n)%ly_j(ni-1,nj))
      allocate(g(n)%ro(ni,nj),g(n)%rovx(ni,nj),g(n)%rovy(ni,nj), &
            g(n)%roe(ni,nj))
      allocate(g(n)%dro(ni-1,nj-1),g(n)%drovx(ni-1,nj-1), &
            g(n)%drovy(ni-1,nj-1),g(n)%droe(ni-1,nj-1))
      allocate(g(n)%p(ni,nj),g(n)%hstag(ni,nj),g(n)%vx(ni,nj), &
            g(n)%vy(ni,nj))
end do
read(2) bcs%n_in, bcs%n_out
write(6,*) 'Read Mesh:  Inlet on block', bcs%n_in, 'outlet on block', bcs%n_out
allocate(bcs%ro(g(bcs%n_in)%nj),bcs%p(g(bcs%n_in)%nj))
read(2) av%nm; allocate(p(av%nm));
write(6,*) 'Read Mesh:  Number of matching patches nm =', av%nm
do m = 1,av%nm 
      read(2) p(m)%nk, p(m)%n_1, p(m)%n_2 
      write(6,*) 'Read Mesh:  Size of patch', m, 'nk =', p(m)%nk, 'joins block', &
            p(m)%n_1, 'and', p(m)%n_2
      allocate(p(m)%i_1(p(m)%nk),p(m)%j_1(p(m)%nk),p(m)%i_2(p(m)%nk), &
            p(m)%j_2(p(m)%nk))
      read(2) p(m)%i_1,p(m)%j_1,p(m)%i_2,p(m)%j_2
end do
close(2)
write(6,*)
end subroutine read_mesh
