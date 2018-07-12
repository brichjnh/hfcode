subroutine time_checker(ff,titstring)
use nrtype
implicit none

integer(i4b), intent(in) ::ff
character(len=*), intent(in) :: titstring

integer(i4b) :: tmin
real(dp), save :: t0, tnow, tprec
real(dp) :: tsec, dt

if (ff.eq.0) then
    call cpu_time(t0)
    tnow=t0
    return
else if (ff.eq.1) then
    tprec=tnow
    call cpu_time(tnow)
    dt=tnow-t0
    tmin=int(dt/60._dp)
    tsec=dt-60._dp*real(tmin,dp)
    write (9,'(A,I3,A,F7.3,A)') titstring,tmin," min",tsec," sec"
    return
else
    tprec=tnow
    call cpu_time(tnow)
    dt=tnow-tprec 
    tmin=int(dt/60._dp)
    tsec=dt-60._dp*real(tmin,dp)
    write (9,'(A,I3,A,F7.3,A)') titstring,tmin," min",tsec," sec"
    return
end if

end subroutine time_checker

subroutine time_diff(ff,cyctime,tottime)
use nrtype
implicit none

integer(i4b), intent(in) ::ff
real(dp), intent(out) :: cyctime,tottime

integer(i4b) :: tmin
real(dp), save :: t0, tnow, tprec

cyctime=0._dp
tottime=0._dp
if (ff.eq.0) then
    call cpu_time(t0)
    tnow=t0
    return
else if (ff.eq.1) then
    tprec=tnow
    call cpu_time(tnow)
    cyctime=tnow-tprec
    tottime=tnow-t0
    return
end if

end subroutine time_diff


subroutine time_diff_omp(ff,titstring)
use nrtype; use omp_lib
implicit none

integer(i4b), intent(in) ::ff
character(len=*), intent(in) :: titstring
real(dp) :: cyctime,tottime

integer(i4b) :: tmin
real(dp), save :: t0, tnow, tprec

cyctime=0._dp
tottime=0._dp
if (ff.eq.0) then
    t0=omp_get_wtime()
    tnow=t0
    return
else if (ff.eq.1) then
    tprec=tnow
    tnow=omp_get_wtime()
    cyctime=tnow-tprec
    tottime=tnow-t0
    write (9,'(A,2F10.3)') titstring,cyctime,tottime
    return
end if

end subroutine time_diff_omp

