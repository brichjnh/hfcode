subroutine evaluate_maxdm(p,maxdm)
use nrtype ; use molprops
implicit none

real(dp), intent(in) :: p(nb,nb)
real(dp), intent(out) :: maxdm(nab)

integer(i4b) :: i, sh1, sh2, ao1, ao2, ab1

do i=1,nss
  sh1=ssprs(1,i)
  sh2=ssprs(2,i)
  ab1=ind2(sh1,sh2)
  ao1=aosh(sh1)
  ao2=aosh(sh2)
  maxdm(ab1)=abs(p(ao1,ao2))
end do

do i=1,nsp
  sh1=spprs(1,i)
  sh2=spprs(2,i)
  ab1=ind2(sh1,sh2)
  ao1=aosh(sh1)
  ao2=aosh(sh2)
  maxdm(ab1)=maxval(abs(p(ao1,ao2:ao2+2)))
end do

do i=1,nsd
  sh1=sdprs(1,i)
  sh2=sdprs(2,i)
  ab1=ind2(sh1,sh2)
  ao1=aosh(sh1)
  ao2=aosh(sh2)
  maxdm(ab1)=maxval(abs(p(ao1,ao2:ao2+4)))
end do

do i=1,nsl
  sh1=slprs(1,i)
  sh2=slprs(2,i)
  ab1=ind2(sh1,sh2)
  ao1=aosh(sh1)
  ao2=aosh(sh2)
  maxdm(ab1)=maxval(abs(p(ao1,ao2:ao2+3)))
end do

do i=1,npp
  sh1=ppprs(1,i)
  sh2=ppprs(2,i)
  ab1=ind2(sh1,sh2)
  ao1=aosh(sh1)
  ao2=aosh(sh2)
  maxdm(ab1)=maxval(abs(p(ao1:ao1+2,ao2:ao2+2)))
end do

do i=1,npd
  sh1=pdprs(1,i)
  sh2=pdprs(2,i)
  ab1=ind2(sh1,sh2)
  ao1=aosh(sh1)
  ao2=aosh(sh2)
  maxdm(ab1)=maxval(abs(p(ao1:ao1+2,ao2:ao2+4)))
end do

do i=1,npl
  sh1=plprs(1,i)
  sh2=plprs(2,i)
  ab1=ind2(sh1,sh2)
  ao1=aosh(sh1)
  ao2=aosh(sh2)
  maxdm(ab1)=maxval(abs(p(ao1:ao1+2,ao2:ao2+3)))
end do

do i=1,ndd
  sh1=ddprs(1,i)
  sh2=ddprs(2,i)
  ab1=ind2(sh1,sh2)
  ao1=aosh(sh1)
  ao2=aosh(sh2)
  maxdm(ab1)=maxval(abs(p(ao1:ao1+4,ao2:ao2+4)))
end do

do i=1,ndl
  sh1=dlprs(1,i)
  sh2=dlprs(2,i)
  ab1=ind2(sh1,sh2)
  ao1=aosh(sh1)
  ao2=aosh(sh2)
  maxdm(ab1)=maxval(abs(p(ao1:ao1+3,ao2:ao2+4)))
end do

do i=1,nll
  sh1=llprs(1,i)
  sh2=llprs(2,i)
  ab1=ind2(sh1,sh2)
  ao1=aosh(sh1)
  ao2=aosh(sh2)
  maxdm(ab1)=maxval(abs(p(ao1:ao1+3,ao2:ao2+3)))
end do

end subroutine evaluate_maxdm

