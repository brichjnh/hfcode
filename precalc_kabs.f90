subroutine precalcs_kabs()
use nrtype ; use molprops
implicit none

! This subroutine evaluates all the functions that depend on only two basis shells
!  and that are needed for evaluating TEIs.
! These are p = a+b, pm1=1/p, q = ab/(a+b), Px=(aAx + bBx)/p, Kab = exp(-q*rab**2)/p*djk(a)*djk(b)

integer(i4b) :: sh1, sh2, aind, ab, zzab, za, zb, zza
integer(i4b) :: i, jss,jsp,jsd,jsl,jpp,jpd,jpl,jdd,jdl,jll
real(dp) :: rab(3), extrathresh
character(3) :: bla

extrathresh=intthresh*0.01_dp

! First generate some functions needed for evaluating integrals

call buildindices()
call buildharmond()
call buildgammarats()
call eval_boys_func_on_grid()

! Start by allocating the arrays. They will be 2-dimensional (3 for Px):
! The first index for Px is the Cartesian coordinates
! The first index, otherwise, is 'kk' and refers to the pair of zetas. It is equal to (k1-1)*maxprim+k2
! The last 'ab'=a*(a-1)/2+b denotes the pair of shells. Only the cases where a>=b are stored

nzz=maxprim**2
allocate(pab(nzz,nab),pm1(nzz,nab),hpm1(nzz,nab),qab(nzz,nab),Px(3,nzz,nab),kab(nzz,nab))
allocate(rab2(nab),negab(nzz,nab))
allocate(ddab(nzz,nab))

! Then work out how many pairs of each type there are, and store them in integer arrays ssprs, etc.

jss=0 ; jsp=0; jsd=0; jsl=0 ; jpp=0
jpd=0 ; jpl=0 ; jdd=0 ; jdl=0 ; jll=0 
negab=.false.

do sh1=1,nsh   ! First create a list of shell pairs of different types
  do sh2=1,sh1
    if ((shtyp(sh1).eq.1).and.(shtyp(sh2).eq.1)) then
      jss=jss+1
      ssprs(:,jss)=(/sh1,sh2/)
    else if ((shtyp(sh1).eq.1).and.(shtyp(sh2).eq.2)) then
      jsp=jsp+1
      spprs(:,jsp)=(/sh1,sh2/)
    else if ((shtyp(sh1).eq.2).and.(shtyp(sh2).eq.1)) then
      jsp=jsp+1
      spprs(:,jsp)=(/sh2,sh1/)
    else if ((shtyp(sh1).eq.1).and.(shtyp(sh2).eq.3)) then
      jsd=jsd+1
      sdprs(:,jsd)=(/sh1,sh2/)
    else if ((shtyp(sh1).eq.3).and.(shtyp(sh2).eq.1)) then
      jsd=jsd+1
      sdprs(:,jsd)=(/sh2,sh1/)
    else if ((shtyp(sh1).eq.2).and.(shtyp(sh2).eq.3)) then
      jpd=jpd+1
      pdprs(:,jpd)=(/sh1,sh2/)
    else if ((shtyp(sh1).eq.3).and.(shtyp(sh2).eq.2)) then
      jpd=jpd+1
      pdprs(:,jpd)=(/sh2,sh1/)
    else if ((shtyp(sh1).eq.1).and.(shtyp(sh2).eq.12)) then
      jsl=jsl+1
      slprs(:,jsl)=(/sh1,sh2/)
    else if ((shtyp(sh1).eq.12).and.(shtyp(sh2).eq.1)) then
      jsl=jsl+1
      slprs(:,jsl)=(/sh2,sh1/)
    else if ((shtyp(sh1).eq.2).and.(shtyp(sh2).eq.2)) then
      jpp=jpp+1
      ppprs(:,jpp)=(/sh1,sh2/)
    else if ((shtyp(sh1).eq.2).and.(shtyp(sh2).eq.12)) then
      jpl=jpl+1
      plprs(:,jpl)=(/sh1,sh2/)
    else if ((shtyp(sh1).eq.12).and.(shtyp(sh2).eq.2)) then
      jpl=jpl+1
      plprs(:,jpl)=(/sh2,sh1/)
    else if ((shtyp(sh1).eq.3).and.(shtyp(sh2).eq.3)) then
      jdd=jdd+1
      ddprs(:,jdd)=(/sh1,sh2/)
    else if ((shtyp(sh1).eq.3).and.(shtyp(sh2).eq.12)) then
      jdl=jdl+1
      dlprs(:,jdl)=(/sh1,sh2/)
    else if ((shtyp(sh1).eq.12).and.(shtyp(sh2).eq.3)) then
      jdl=jdl+1
      dlprs(:,jdl)=(/sh2,sh1/)
    else if ((shtyp(sh1).eq.12).and.(shtyp(sh2).eq.12)) then
      jll=jll+1
      llprs(:,jll)=(/sh1,sh2/)
    end if
  end do
end do

! Then at last loop over pairs, and precalc various things

do i=1,nss
  sh1=ssprs(1,i)
  sh2=ssprs(2,i)
  ab=ind2(sh1,sh2)
  rab2(ab)=r2AB(atsh(sh1),atsh(sh2))
  do za=1,nzet(sh1)
    zza=(za-1)*maxprim  
    do zb=1,nzet(sh2)
      zzab=zza+zb
      pab(zzab,ab)=zeta(za,sh1)+zeta(zb,sh2)
      pm1(zzab,ab)=1._dp/pab(zzab,ab)
      hpm1(zzab,ab)=.5_dp*pm1(zzab,ab)
      qab(zzab,ab)=zeta(za,sh1)*zeta(zb,sh2)*pm1(zzab,ab)
      Px(:,zzab,ab)=(zeta(za,sh1)*xsh(:,sh1)+zeta(zb,sh2)*xsh(:,sh2))*pm1(zzab,ab)
      ddab(zzab,ab)=djk(za,sh1)*djk(zb,sh2)
      kab(zzab,ab)=rt2pi54*exp(-qab(zzab,ab)*rab2(ab))*pm1(zzab,ab)*ddab(zzab,ab)
      if (abs(kab(zzab,ab))*sqrt(pm1(zzab,ab)).lt.extrathresh) then
          negab(zzab,ab)=.true.
      end if
    end do
  end do
end do
do i=1,nsp
  sh1=spprs(1,i)
  sh2=spprs(2,i)
  ab=ind2(sh1,sh2)
  rab2(ab)=r2AB(atsh(sh1),atsh(sh2))
  do za=1,nzet(sh1)
    zza=(za-1)*maxprim  
    do zb=1,nzet(sh2)
      zzab=zza+zb
      pab(zzab,ab)=zeta(za,sh1)+zeta(zb,sh2)
      pm1(zzab,ab)=1._dp/pab(zzab,ab)
      hpm1(zzab,ab)=.5_dp*pm1(zzab,ab)
      qab(zzab,ab)=zeta(za,sh1)*zeta(zb,sh2)*pm1(zzab,ab)
      Px(:,zzab,ab)=(zeta(za,sh1)*xsh(:,sh1)+zeta(zb,sh2)*xsh(:,sh2))*pm1(zzab,ab)
      ddab(zzab,ab)=djk(za,sh1)*djk(zb,sh2)
      kab(zzab,ab)=rt2pi54*exp(-qab(zzab,ab)*rab2(ab))*pm1(zzab,ab)*ddab(zzab,ab)
      if (abs(kab(zzab,ab))*sqrt(pm1(zzab,ab)).lt.extrathresh) negab(zzab,ab)=.true.
    end do
  end do
end do
do i=1,nsd
  sh1=sdprs(1,i)
  sh2=sdprs(2,i)
  ab=ind2(sh1,sh2)
  rab2(ab)=r2AB(atsh(sh1),atsh(sh2))
  do za=1,nzet(sh1)
    zza=(za-1)*maxprim  
    do zb=1,nzet(sh2)
      zzab=zza+zb
      pab(zzab,ab)=zeta(za,sh1)+zeta(zb,sh2)
      pm1(zzab,ab)=1._dp/pab(zzab,ab)
      hpm1(zzab,ab)=.5_dp*pm1(zzab,ab)
      qab(zzab,ab)=zeta(za,sh1)*zeta(zb,sh2)*pm1(zzab,ab)
      Px(:,zzab,ab)=(zeta(za,sh1)*xsh(:,sh1)+zeta(zb,sh2)*xsh(:,sh2))*pm1(zzab,ab)
      ddab(zzab,ab)=djk(za,sh1)*djk(zb,sh2)
      kab(zzab,ab)=rt2pi54*exp(-qab(zzab,ab)*rab2(ab))*pm1(zzab,ab)*ddab(zzab,ab)
      if (abs(kab(zzab,ab))*sqrt(pm1(zzab,ab)).lt.extrathresh) negab(zzab,ab)=.true.
    end do
  end do
end do
do i=1,nsl
  sh1=slprs(1,i)
  sh2=slprs(2,i)
  ab=ind2(sh1,sh2)
  rab2(ab)=r2AB(atsh(sh1),atsh(sh2))
  do za=1,nzet(sh1)
    zza=(za-1)*maxprim  
    do zb=1,nzet(sh2)
      zzab=zza+zb
      pab(zzab,ab)=zeta(za,sh1)+zeta(zb,sh2)
      pm1(zzab,ab)=1._dp/pab(zzab,ab)
      hpm1(zzab,ab)=.5_dp*pm1(zzab,ab)
      qab(zzab,ab)=zeta(za,sh1)*zeta(zb,sh2)*pm1(zzab,ab)
      Px(:,zzab,ab)=(zeta(za,sh1)*xsh(:,sh1)+zeta(zb,sh2)*xsh(:,sh2))*pm1(zzab,ab)
      ddab(zzab,ab)=djk(za,sh1)*djk(zb,sh2)
      kab(zzab,ab)=rt2pi54*exp(-qab(zzab,ab)*rab2(ab))*pm1(zzab,ab)*ddab(zzab,ab)
      if (abs(kab(zzab,ab))*sqrt(pm1(zzab,ab)).lt.extrathresh) negab(zzab,ab)=.true.
    end do
  end do
end do
do i=1,npp
  sh1=ppprs(1,i)
  sh2=ppprs(2,i)
  ab=ind2(sh1,sh2)
  rab2(ab)=r2AB(atsh(sh1),atsh(sh2))
  do za=1,nzet(sh1)
    zza=(za-1)*maxprim  
    do zb=1,nzet(sh2)
      zzab=zza+zb
      pab(zzab,ab)=zeta(za,sh1)+zeta(zb,sh2)
      pm1(zzab,ab)=1._dp/pab(zzab,ab)
      hpm1(zzab,ab)=.5_dp*pm1(zzab,ab)
      qab(zzab,ab)=zeta(za,sh1)*zeta(zb,sh2)*pm1(zzab,ab)
      Px(:,zzab,ab)=(zeta(za,sh1)*xsh(:,sh1)+zeta(zb,sh2)*xsh(:,sh2))*pm1(zzab,ab)
      ddab(zzab,ab)=djk(za,sh1)*djk(zb,sh2)
      kab(zzab,ab)=rt2pi54*exp(-qab(zzab,ab)*rab2(ab))*pm1(zzab,ab)*ddab(zzab,ab)
      if (abs(kab(zzab,ab))*sqrt(pm1(zzab,ab)).lt.extrathresh) negab(zzab,ab)=.true.
    end do
  end do
end do
do i=1,npd
  sh1=pdprs(1,i)
  sh2=pdprs(2,i)
  ab=ind2(sh1,sh2)
  rab2(ab)=r2AB(atsh(sh1),atsh(sh2))
  do za=1,nzet(sh1)
    zza=(za-1)*maxprim  
    do zb=1,nzet(sh2)
      zzab=zza+zb
      pab(zzab,ab)=zeta(za,sh1)+zeta(zb,sh2)
      pm1(zzab,ab)=1._dp/pab(zzab,ab)
      hpm1(zzab,ab)=.5_dp*pm1(zzab,ab)
      qab(zzab,ab)=zeta(za,sh1)*zeta(zb,sh2)*pm1(zzab,ab)
      Px(:,zzab,ab)=(zeta(za,sh1)*xsh(:,sh1)+zeta(zb,sh2)*xsh(:,sh2))*pm1(zzab,ab)
      ddab(zzab,ab)=djk(za,sh1)*djk(zb,sh2)
      kab(zzab,ab)=rt2pi54*exp(-qab(zzab,ab)*rab2(ab))*pm1(zzab,ab)*ddab(zzab,ab)
      if (abs(kab(zzab,ab))*sqrt(pm1(zzab,ab)).lt.extrathresh) negab(zzab,ab)=.true.
    end do
  end do
end do
do i=1,npl
  sh1=plprs(1,i)
  sh2=plprs(2,i)
  ab=ind2(sh1,sh2)
  rab2(ab)=r2AB(atsh(sh1),atsh(sh2))
  do za=1,nzet(sh1)
    zza=(za-1)*maxprim  
    do zb=1,nzet(sh2)
      zzab=zza+zb
      pab(zzab,ab)=zeta(za,sh1)+zeta(zb,sh2)
      pm1(zzab,ab)=1._dp/pab(zzab,ab)
      hpm1(zzab,ab)=.5_dp*pm1(zzab,ab)
      qab(zzab,ab)=zeta(za,sh1)*zeta(zb,sh2)*pm1(zzab,ab)
      Px(:,zzab,ab)=(zeta(za,sh1)*xsh(:,sh1)+zeta(zb,sh2)*xsh(:,sh2))*pm1(zzab,ab)
      ddab(zzab,ab)=djk(za,sh1)*djk(zb,sh2)
      kab(zzab,ab)=rt2pi54*exp(-qab(zzab,ab)*rab2(ab))*pm1(zzab,ab)*ddab(zzab,ab)
      if (abs(kab(zzab,ab))*sqrt(pm1(zzab,ab)).lt.extrathresh) negab(zzab,ab)=.true.
    end do
  end do
end do
do i=1,ndd
  sh1=ddprs(1,i)
  sh2=ddprs(2,i)
  ab=ind2(sh1,sh2)
  rab2(ab)=r2AB(atsh(sh1),atsh(sh2))
  do za=1,nzet(sh1)
    zza=(za-1)*maxprim  
    do zb=1,nzet(sh2)
      zzab=zza+zb
      pab(zzab,ab)=zeta(za,sh1)+zeta(zb,sh2)
      pm1(zzab,ab)=1._dp/pab(zzab,ab)
      hpm1(zzab,ab)=.5_dp*pm1(zzab,ab)
      qab(zzab,ab)=zeta(za,sh1)*zeta(zb,sh2)*pm1(zzab,ab)
      Px(:,zzab,ab)=(zeta(za,sh1)*xsh(:,sh1)+zeta(zb,sh2)*xsh(:,sh2))*pm1(zzab,ab)
      ddab(zzab,ab)=djk(za,sh1)*djk(zb,sh2)
      kab(zzab,ab)=rt2pi54*exp(-qab(zzab,ab)*rab2(ab))*pm1(zzab,ab)*ddab(zzab,ab)
      if (abs(kab(zzab,ab))*sqrt(pm1(zzab,ab)).lt.extrathresh) negab(zzab,ab)=.true.
    end do
  end do
end do
do i=1,ndl
  sh1=dlprs(1,i)
  sh2=dlprs(2,i)
  ab=ind2(sh1,sh2)
  rab2(ab)=r2AB(atsh(sh1),atsh(sh2))
  do za=1,nzet(sh1)
    zza=(za-1)*maxprim  
    do zb=1,nzet(sh2)
      zzab=zza+zb
      pab(zzab,ab)=zeta(za,sh1)+zeta(zb,sh2)
      pm1(zzab,ab)=1._dp/pab(zzab,ab)
      hpm1(zzab,ab)=.5_dp*pm1(zzab,ab)
      qab(zzab,ab)=zeta(za,sh1)*zeta(zb,sh2)*pm1(zzab,ab)
      Px(:,zzab,ab)=(zeta(za,sh1)*xsh(:,sh1)+zeta(zb,sh2)*xsh(:,sh2))*pm1(zzab,ab)
      ddab(zzab,ab)=djk(za,sh1)*djk(zb,sh2)
      kab(zzab,ab)=rt2pi54*exp(-qab(zzab,ab)*rab2(ab))*pm1(zzab,ab)*ddab(zzab,ab)
      if (abs(kab(zzab,ab))*sqrt(pm1(zzab,ab)).lt.extrathresh) negab(zzab,ab)=.true.
    end do
  end do
end do
do i=1,nll
  sh1=llprs(1,i)
  sh2=llprs(2,i)
  ab=ind2(sh1,sh2)
  rab2(ab)=r2AB(atsh(sh1),atsh(sh2))
  do za=1,nzet(sh1)
    zza=(za-1)*maxprim  
    do zb=1,nzet(sh2)
      zzab=zza+zb
      pab(zzab,ab)=zeta(za,sh1)+zeta(zb,sh2)
      pm1(zzab,ab)=1._dp/pab(zzab,ab)
      hpm1(zzab,ab)=.5_dp*pm1(zzab,ab)
      qab(zzab,ab)=zeta(za,sh1)*zeta(zb,sh2)*pm1(zzab,ab)
      Px(:,zzab,ab)=(zeta(za,sh1)*xsh(:,sh1)+zeta(zb,sh2)*xsh(:,sh2))*pm1(zzab,ab)
      ddab(zzab,ab)=djk(za,sh1)*djk(zb,sh2)
      kab(zzab,ab)=rt2pi54*exp(-qab(zzab,ab)*rab2(ab))*pm1(zzab,ab)*ddab(zzab,ab)
      if (abs(kab(zzab,ab))*sqrt(pm1(zzab,ab)).lt.extrathresh) negab(zzab,ab)=.true.
    end do
  end do
end do
end subroutine precalcs_kabs

