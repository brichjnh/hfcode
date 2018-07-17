subroutine save_abcd_integral(i,j,k,l,i2eval)
use nrtype; use molprops
implicit none

integer(i4b), intent(in) :: i,j,k,l
real(dp), intent(in) :: i2eval

logical(lgt) :: ieqj, keql, jeqk, ieqk, ieql, jeql

ieqj=(i.eq.j)
jeqk=(j.eq.k)
keql=(k.eq.l)
if (ieqj.and.jeqk.and.keql) then
  naaaa=naaaa+1
  i2eaaaa(naaaa)=i2eval
  aonaaaa(naaaa)=i
  return
else if (ieqj.and.keql) then
  naabb=naabb+1
  i2eaabb(naabb)=i2eval
  aonaabb(naabb)=ind2(i,k)
  return
else if (ieqj) then
  naabc=naabc+1
  i2eaabc(naabc)=i2eval
  aonaabc(:,naabc)=(/i,ind2(k,l)/)
  return
else if (keql) then
  naabc=naabc+1
  i2eaabc(naabc)=i2eval
  aonaabc(:,naabc)=(/k,ind2(i,j)/)
  return
end if
if ((jeqk.and.(i.eq.l)).or.((i.eq.k).and.(j.eq.l))) then
  nabab=nabab+1
  i2eabab(nabab)=i2eval
  aonabab(nabab)=ind2(i,j)
  return
end if
nabcd=nabcd+1
i2eabcd(nabcd)=i2eval
aonabcd(:,nabcd)=(/ind2(i,j),ind2(k,l)/)
return

end subroutine save_abcd_integral

subroutine save_abcd_integral_d(i,j,k,l,i2eval)
use nrtype; use molprops
implicit none

integer(i4b), intent(in) :: i,j,k,l
real(dp), intent(in) :: i2eval

logical(lgt) :: ieqj, keql, jeqk, ieqk, ieql, jeql

ieqj=(i.eq.j)
jeqk=(j.eq.k)
keql=(k.eq.l)
if (ieqj.and.jeqk.and.keql) then
  naaaa=naaaa+1
  i2eaaaa(naaaa)=i2eval
  aonaaaa(naaaa)=i
  return
else if (ieqj.and.keql) then
  naabb=naabb+1
  i2eaabb(naabb)=i2eval
  aonaabb(naabb)=ind2(i,k)
  return
else if (ieqj) then
  naabc=naabc+1
  i2eaabc(naabc)=i2eval
  aonaabc(:,naabc)=(/i,ind2(k,l)/)
  return
else if (keql) then
  naabc=naabc+1
  i2eaabc(naabc)=i2eval
  aonaabc(:,naabc)=(/k,ind2(i,j)/)
  return
end if
if ((jeqk.and.(i.eq.l)).or.((i.eq.k).and.(j.eq.l))) then
  nabab=nabab+1
  i2eabab(nabab)=i2eval
  aonabab(nabab)=ind2(i,j)
  return
end if
nabcd_d=nabcd_d+1
i2eabcd_d(nabcd_d)=i2eval
aonabcd_d(:,nabcd_d)=(/ind2(i,j),ind2(k,l)/)
return

end subroutine save_abcd_integral_d

