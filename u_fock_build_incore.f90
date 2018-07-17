subroutine u_fock_build_incore(fa,fb,pa,pb,pj)
use nrtype; use molprops
implicit none

real(dp), intent(in) :: pa(nb,nb), pb(nb,nb), pj(nb,nb)
real(dp), intent(inout) :: fa(nb,nb), fb(nb,nb)

integer(i4b) :: ij,kl,ijkl,ia,ja,ka,la,i,j
integer(i8b) :: ii
logical(lgt) :: sameabab, allsame, sameij, samejk, samekl
real(dp) :: fac, faca, facb

kl=1; ij=1
do i=1,naaaa
  ia=aonaaaa(i)
  fac=pj(ia,ia)*i2eaaaa(i)
  fa(ia,ia)=fa(ia,ia)+fac
  fb(ia,ia)=fb(ia,ia)+fac
  fa(ia,ia)=fa(ia,ia)-pa(ia,ia)*i2eaaaa(i)
  fb(ia,ia)=fb(ia,ia)-pb(ia,ia)*i2eaaaa(i)
end do

do i=1,naabb
  ia=revind2(1,aonaabb(i)) ; ka=revind2(2,aonaabb(i))
  fac=pj(ka,ka)*i2eaabb(i)
  fa(ia,ia)=fa(ia,ia)+fac
  fb(ia,ia)=fb(ia,ia)+fac
  fac=pj(ia,ia)*i2eaabb(i)
  fa(ka,ka)=fa(ka,ka)+fac
  fb(ka,ka)=fb(ka,ka)+fac
  fac=pa(ia,ka)*i2eaabb(i)
  fa(ia,ka)=fa(ia,ka)-fac
  fac=pb(ia,ka)*i2eaabb(i)
  fb(ia,ka)=fb(ia,ka)-fac
end do

do i=1,naabc
  ia=aonaabc(1,i)
  ka=revind2(1,aonaabc(2,i)) ; la=revind2(2,aonaabc(2,i))
  fac=2._dp*pj(ka,la)*i2eaabc(i)
  fa(ia,ia)=fa(ia,ia)+fac
  fb(ia,ia)=fb(ia,ia)+fac
  faca=pa(ia,ka)*i2eaabc(i)
  facb=pb(ia,ka)*i2eaabc(i)
  if (ia.gt.la) then
    fa(ia,la)=fa(ia,la)-faca
    fb(ia,la)=fb(ia,la)-facb
  else if (ia.eq.la) then
    fa(la,ia)=fa(la,ia)-faca*2._dp
    fb(la,ia)=fb(la,ia)-facb*2._dp
  else
    fa(la,ia)=fa(la,ia)-faca
    fb(la,ia)=fb(la,ia)-facb
  end if
  faca=pa(ia,la)*i2eaabc(i)
  facb=pb(ia,la)*i2eaabc(i)
  if (ia.gt.ka) then
    fa(ia,ka)=fa(ia,ka)-faca
    fb(ia,ka)=fb(ia,ka)-facb
  else if (ia.eq.ka) then
    fa(ia,ka)=fa(ia,ka)-faca*2._dp
    fb(ia,ka)=fb(ia,ka)-facb*2._dp
  else
    fa(ka,ia)=fa(ka,ia)-faca
    fb(ka,ia)=fb(ka,ia)-facb
  end if
  fac=pj(ia,ia)*i2eaabc(i)
  if (ka.gt.la) then
    fa(ka,la)=fa(ka,la)+fac 
    fb(ka,la)=fb(ka,la)+fac 
  else if (ka.eq.la) then
    fa(ka,la)=fa(ka,la)+fac*2._dp
    fb(ka,la)=fb(ka,la)+fac*2._dp
  else
    fa(la,ka)=fa(la,ka)+fac 
    fb(la,ka)=fb(la,ka)+fac 
  end if
end do
do i=1,nabab
  ia=revind2(1,aonabab(i)) ; ja=revind2(2,aonabab(i))
  fac=2._dp*pj(ia,ja)*i2eabab(i)
  Fa(ia,ja) = Fa(ia,ja) + fac
  Fb(ia,ja) = Fb(ia,ja) + fac
  faca=pa(ia,ja)*i2eabab(i)
  fa(ia,ja)=fa(ia,ja)-faca
  faca=pb(ia,ja)*i2eabab(i)
  fb(ia,ja)=fb(ia,ja)-facb
  faca=pa(ja,ja)*i2eabab(i)
  fa(ia,ia)=fa(ia,ia)-faca
  facb=pb(ja,ja)*i2eabab(i)
  fb(ia,ia)=fb(ia,ia)-facb
  faca=pa(ia,ia)*i2eabab(i)
  fa(ja,ja)=fa(ja,ja)-faca
  facb=pb(ia,ia)*i2eabab(i)
  fb(ja,ja)=fb(ja,ja)-facb
end do
!write (*,*) "Done",nabab,"(ab|ab) ints"
do i=1,nabcd_d
  ia=revind2(1,aonabcd_d(1,i)) ; ja=revind2(2,aonabcd_d(1,i))
  ka=revind2(1,aonabcd_d(2,i)) ; la=revind2(2,aonabcd_d(2,i))
  fac=2._dp*pj(la,ka)*i2eabcd_d(i)
  fa(ia,ja)=fa(ia,ja)+fac
  fb(ia,ja)=fb(ia,ja)+fac
  faca=pa(ja,ka)*i2eabcd_d(i)
  facb=pb(ja,ka)*i2eabcd_d(i)
  if (ia.gt.la) then
    fa(ia,la)=fa(ia,la)-faca
    fb(ia,la)=fb(ia,la)-facb
  else if (ia.eq.la) then
    fa(ia,la)=fa(ia,la)-faca*2._dp
    fb(ia,la)=fb(ia,la)-facb*2._dp
  else
    fa(la,ia)=fa(la,ia)-faca
    fb(la,ia)=fb(la,ia)-facb
  end if
  faca=pa(ja,la)*i2eabcd_d(i)
  facb=pb(ja,la)*i2eabcd_d(i)
  if (ia.gt.ka) then
    fa(ia,ka)=fa(ia,ka)-faca
    fb(ia,ka)=fb(ia,ka)-facb
  else if (ia.eq.ka) then
    fa(ia,ka)=fa(ia,ka)-faca*2._dp
    fb(ia,ka)=fb(ia,ka)-facb*2._dp
  else if (ia.lt.ka) then
    fa(ka,ia)=fa(ka,ia)-faca
    fb(ka,ia)=fb(ka,ia)-facb
  end if
  faca=pa(ia,ka)*i2eabcd_d(i)
  facb=pb(ia,ka)*i2eabcd_d(i)
  if (ja.gt.la) then
    fa(ja,la)=fa(ja,la)-faca
    fb(ja,la)=fb(ja,la)-facb
  else if (ja.eq.la) then
    fa(ja,la)=fa(ja,la)-faca*2._dp
    fb(ja,la)=fb(ja,la)-facb*2._dp
  else if (ja.lt.la) then
    fa(la,ja)=fa(la,ja)-faca
    fb(la,ja)=fb(la,ja)-facb
  end if
  faca=pa(ia,la)*i2eabcd_d(i)
  facb=pb(ia,la)*i2eabcd_d(i)
  if (ja.gt.ka) then
    fa(ja,ka)=fa(ja,ka)-faca
    fb(ja,ka)=fb(ja,ka)-facb
  else if (ja.eq.ka) then
    fa(ja,ka)=fa(ja,ka)-faca*2._dp
    fb(ja,ka)=fb(ja,ka)-facb*2._dp
  else if (ja.lt.ka) then
    fa(ka,ja)=fa(ka,ja)-faca
    fb(ka,ja)=fb(ka,ja)-facb
  end if
  fac=2._dp*pj(ia,ja)*i2eabcd_d(i)
  fa(ka,la)=fa(ka,la)+fac
  fb(ka,la)=fb(ka,la)+fac
end do
do ii=1,nabcd
  ia=revind2(1,aonabcd(1,ii)) ; ja=revind2(2,aonabcd(1,ii))
  ka=revind2(1,aonabcd(2,ii)) ; la=revind2(2,aonabcd(2,ii))
  fac=2._dp*pj(la,ka)*i2eabcd(ii)
  fa(ia,ja)=fa(ia,ja)+fac
  fb(ia,ja)=fb(ia,ja)+fac
  faca=pa(ja,ka)*i2eabcd(ii)
  facb=pb(ja,ka)*i2eabcd(ii)
  if (ia.gt.la) then
    fa(ia,la)=fa(ia,la)-faca
    fb(ia,la)=fb(ia,la)-facb
  else if (ia.eq.la) then
    fa(ia,la)=fa(ia,la)-faca*2._dp
    fb(ia,la)=fb(ia,la)-facb*2._dp
  else
    fa(la,ia)=fa(la,ia)-faca
    fb(la,ia)=fb(la,ia)-facb
  end if
  faca=pa(ja,la)*i2eabcd(ii)
  facb=pb(ja,la)*i2eabcd(ii)
  if (ia.gt.ka) then
    fa(ia,ka)=fa(ia,ka)-faca
    fb(ia,ka)=fb(ia,ka)-facb
  else if (ia.eq.ka) then
    fa(ia,ka)=fa(ia,ka)-faca*2._dp
    fb(ia,ka)=fb(ia,ka)-facb*2._dp
  else if (ia.lt.ka) then
    fa(ka,ia)=fa(ka,ia)-faca
    fb(ka,ia)=fb(ka,ia)-facb
  end if
  faca=pa(ia,ka)*i2eabcd(ii)
  facb=pb(ia,ka)*i2eabcd(ii)
  if (ja.gt.la) then
    fa(ja,la)=fa(ja,la)-faca
    fb(ja,la)=fb(ja,la)-facb
  else if (ja.eq.la) then
    fa(ja,la)=fa(ja,la)-faca*2._dp
    fb(ja,la)=fb(ja,la)-facb*2._dp
  else if (ja.lt.la) then
    fa(la,ja)=fa(la,ja)-faca
    fb(la,ja)=fb(la,ja)-facb
  end if
  faca=pa(ia,la)*i2eabcd(ii)
  facb=pb(ia,la)*i2eabcd(ii)
  if (ja.gt.ka) then
    fa(ja,ka)=fa(ja,ka)-faca
    fb(ja,ka)=fb(ja,ka)-facb
  else if (ja.eq.ka) then
    fa(ja,ka)=fa(ja,ka)-faca*2._dp
    fb(ja,ka)=fb(ja,ka)-facb*2._dp
  else if (ja.lt.ka) then
    fa(ka,ja)=fa(ka,ja)-faca
    fb(ka,ja)=fb(ka,ja)-facb
  end if
  fac=2._dp*pj(ia,ja)*i2eabcd(ii)
  fa(ka,la)=fa(ka,la)+fac
  fb(ka,la)=fb(ka,la)+fac
end do
!write (*,*) "Done",nabcd,"(ab|cd) ints"
do i=2,nb
  do j=1,i-1
    fa(j,i)=fa(i,j)
    fb(j,i)=fb(i,j)
  end do
end do

end subroutine u_fock_build_incore

