SUBROUTINE index_sort(n,arr,indarr)
USE nrtype
implicit none 

integer(i4b), intent(in) :: n
integer(i4b), intent(out) :: indarr(n)
REAL(dP), INTENT(IN) :: arr(n)

INTEGER(I4B), PARAMETER :: NN=15, NSTACK=50

REAL(dp) :: a
INTEGER(I4B) :: k,i,j,indext,jstack,l,r, istack(nstack), bla

do i=1,n
  indarr(i)=i
end do

jstack=0
l=1
r=n
do
  if (r-l < NN) then
    do j=l+1,r
      indext=indarr(j)
      a=arr(indext)
      do i=j-1,l,-1
        if (arr(indarr(i)) .ge. a) exit
        indarr(i+1)=indarr(i)
      end do
      indarr(i+1)=indext
    end do
    if (jstack .eq. 0) RETURN
    r=istack(jstack)
    l=istack(jstack-1)
    jstack=jstack-2
  else
    k=(l+r)/2
    bla=indarr(k)
    indarr(k)=indarr(l+1)
    indarr(l+1)=bla
    call icomp_xchg(indarr(l),indarr(r))
    call icomp_xchg(indarr(l+1),indarr(r))
    call icomp_xchg(indarr(l),indarr(l+1))
    i=l+1
    j=r
    indext=indarr(l+1)
    a=arr(indext)
    do
      do
        i=i+1
        if (arr(indarr(i)) .le. a) exit
      end do
      do
        j=j-1
        if (arr(indarr(j)) .ge. a) exit
      end do
      if (j < i) exit
      bla=indarr(i)
      indarr(i)=indarr(j)
      indarr(j)=bla 
    end do
    indarr(l+1)=indarr(j)
    indarr(j)=indext
    jstack=jstack+2
    if (jstack > NSTACK) then
       write (*,*) "error in schw sort"
       stop
    end if
    if (r-i+1 .ge. j-l) then
      istack(jstack)=r
      istack(jstack-1)=i
      r=j-1
    else
      istack(jstack)=j-1
      istack(jstack-1)=l
      l=i
    end if
  end if
end do

CONTAINS

SUBROUTINE icomp_xchg(i,j)
INTEGER(I4B), INTENT(INOUT) :: i,j
INTEGER(I4B) :: swp
if (arr(j) .gt. arr(i)) then
    swp=i
    i=j
    j=swp
end if
END SUBROUTINE icomp_xchg
END SUBROUTINE index_sort
