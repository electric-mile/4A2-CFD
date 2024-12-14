subroutine check_conv(av,g,d_avg,d_max)
use types
use routines
implicit none
type(t_appvars), intent(in) :: av
type(t_grid), intent(in) :: g
real, intent(out) :: d_avg, d_max
real, dimension(g%ni-1,g%nj-1) :: dro, droe, drovx, drovy
integer :: ijx_max(2), ijy_max(2), ij_max(2), ncells
real :: dro_max, drovx_max, drovy_max, droe_max, dro_avg, drovx_avg, &
      drovy_avg, droe_avg, flow_ratio
character(len=100) :: fmt_step
ncells = size(g%dro)
dro = abs(g%dro); droe = abs(g%droe);
drovx = abs(g%drovx); drovy = abs(g%drovy);
dro_avg = sum(abs(dro)) / (ncells * av%ro_ref)
droe_avg = sum(abs(droe)) / (ncells * av%roe_ref)
drovx_avg = sum(abs(drovx)) / (ncells * av%rov_ref)
drovy_avg = sum(abs(drovy)) / (ncells * av%rov_ref)
dro_max = maxval(dro) / av%ro_ref; droe_max = maxval(droe) / av%roe_ref;
ijx_max = maxloc(drovx); ijy_max = maxloc(drovy);
drovx_max = drovx(ijx_max(1),ijx_max(2)) / av%rov_ref
drovy_max = drovy(ijy_max(1),ijy_max(2)) / av%rov_ref
if(drovx_avg > drovy_avg) then
      d_max = drovx_max; d_avg = drovx_avg; ij_max = ijx_max;
else
      d_max = drovy_max; d_avg = drovy_avg; ij_max = ijy_max;
end if
write(3,'(i13,8e15.6)') av%nstep, dro_avg, droe_avg, drovx_avg, &
      drovy_avg, dro_max, droe_max, drovx_max, drovy_max
write(6,*) 'Time step number ', av%nstep
fmt_step = '(a,e10.3,a,i4,a,i4,a,e10.3)'
write(*,fmt_step) '   d_max =', d_max, ' at i =', ij_max(1), ', j =', &
      ij_max(2), ', d_avg =', d_avg
end subroutine check_conv
