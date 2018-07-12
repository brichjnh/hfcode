subroutine fock_build_incore(f,p)
use nrtype; use molprops
implicit none

real(dp), intent(in) :: p(nb,nb)
real(dp), intent(inout) :: f(nb,nb)

integer(i4b) :: ij,kl,ijkl,ia,ja,ka,la
logical(lgt) :: sameabab, allsame, sameij, samekl
real(dp) :: fac

kl=1; ij=1
do ijkl=1,nb4
   ia=revind2(1,ij) ; ja=revind2(2,ij)
   ka=revind2(1,kl) ; la=revind2(2,kl)
   kl=kl+1
   if (kl.gt.ij) then
     kl=1
     ij=ij+1
   end if
   if (abs(i2e(ijkl)).lt.intthresh) cycle
   sameij=(ia.eq.ja)
   samekl=(ka.eq.la)
   sameabab=(((ia.eq.ka).and.(ja.eq.la)).or.((ia.eq.la).and.(ja.eq.ka)))
   allsame=(sameij.and.samekl.and.(ja.eq.ka))
   if (allsame) then
     f(ia,ia)=f(ia,ia)+.5_dp*p(ia,ia)*i2e(ijkl)
   else if (sameij.and.samekl) then
     f(ia,ia)=f(ia,ia)+p(ka,ka)*i2e(ijkl)
     fac=.5_dp*p(ia,ka)*i2e(ijkl)
     f(ia,ka)=f(ia,ka)-fac
     f(ka,ia)=f(ka,ia)-fac
     f(ka,ka)=f(ka,ka)+p(ia,ia)*i2e(ijkl)
   else if (sameij) then
     f(ia,ia)=f(ia,ia)+2._dp*p(ka,la)*i2e(ijkl)
     fac=.5_dp*p(ia,ka)*i2e(ijkl)
     f(ia,la)=f(ia,la)-fac
     f(la,ia)=f(la,ia)-fac
     fac=.5_dp*p(ia,la)*i2e(ijkl)
     f(ia,ka)=f(ia,ka)-fac
     f(ka,ia)=f(ka,ia)-fac
     fac=p(ia,ia)*i2e(ijkl)
     f(ka,la)=f(ka,la)+fac
     f(la,ka)=f(la,ka)+fac
   else if (samekl) then
     f(ka,ka)=f(ka,ka)+2._dp*p(ia,ja)*i2e(ijkl)
     fac=.5_dp*p(ka,ia)*i2e(ijkl)
     f(ka,ja)=f(ka,ja)-fac
     f(ja,ka)=f(ja,ka)-fac
     fac=.5_dp*p(ka,ja)*i2e(ijkl)
     f(ka,ia)=f(ka,ia)-fac
     f(ia,ka)=f(ia,ka)-fac
     fac=p(ka,ka)*i2e(ijkl)
     f(ia,ja)=f(ia,ja)+fac
     f(ja,ia)=f(ja,ia)+fac
   else if (sameabab) then
     fac=1.5_dp*p(ja,ia)*i2e(ijkl)
     F(ia,ja) = F(ia,ja) + fac
     F(ja,ia) = F(ja,ia) + fac
     F(ia,ia) = F(ia,ia) - .5_dp*p(ja,ja) * i2e(ijkl)
     F(ja,ja) = F(ja,ja) - .5_dp*p(ia,ia) * i2e(ijkl)
   else
     fac=2._dp*  p(la,ka) * i2e(ijkl)
     F(ia,ja) = F(ia,ja) +  fac
     F(ja,ia) = F(ja,ia) +  fac
     fac=.5_dp*p(ja,ka)*i2e(ijkl)
     F(ia,la) = F(ia,la) - fac
     F(la,ia) = F(la,ia) - fac
     fac=.5_dp*p(ja,la)*i2e(ijkl)
     F(ia,ka) = F(ia,ka) - fac
     F(ka,ia) = F(ka,ia) - fac
     fac=.5_dp*p(ia,ka)*i2e(ijkl)
     F(ja,la) = F(ja,la) - fac
     F(la,ja) = F(la,ja) - fac
     fac=.5_dp*p(ia,la)*i2e(ijkl)
     F(ja,ka) = F(ja,ka) - fac
     F(ka,ja) = F(ka,ja) - fac
     fac=2._dp*p(ja,ia)*i2e(ijkl)
     F(ka,la) = F(ka,la) +  fac
     F(la,ka) = F(la,ka) +   fac
   end if
end do

end subroutine fock_build_incore

