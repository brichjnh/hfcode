subroutine cmpt2e()
use nrtype; use molprops ; use s_sp_l_terms ; use s_d_l_terms ; use s_p_d_terms
implicit none

integer(i4b) :: ab1, ab2, sh1,sh2,sh3,sh4, ao1, ao2,ao3,ao4, nschw,atc,atd, n1, n2, n3, n4
integer(i4b) :: ii, jj, i, j, k, l, ij, kl, ijkl, ja, ka, la, ierr
real(dp) :: essss, esssp(3), esspp(3,3), espsp(3,3), esppp(3,3,3), epppp(3,3,3,3)
real(dp) :: esssd(5), essdd(5,5), esdsd(5,5), esddd(5,5,5), edddd(5,5,5,5)
real(dp) :: esspd(3,5), espsd(3,5), esppd(3,3,5), esdpp(5,3,3)
real(dp) :: espdd(3,5,5), esdpd(5,3,5)
real(dp) :: epppd(3,3,3,5), eppdd(3,3,5,5), epdpd(3,5,3,5), epddd(3,5,5,5)
real(dp) :: esssl(4), eslsl(4,4), essll(4,4), eslll(4,4,4), ellll(4,4,4,4)
real(dp) :: epppl(3,3,3,4), eplpl(3,4,3,4), eppll(3,3,4,4), eplll(3,4,4,4)
real(dp) :: esspl(3,4), espsl(3,4), espll(3,4,4), eslpl(4,3,4)
real(dp) :: esppl(3,3,4), eslpp(4,3,3)
real(dp) :: essdl(5,4), esdsl(5,4), esddl(5,5,4), esldd(4,5,5)
real(dp) :: esdll(5,4,4), esldl(4,5,4)
real(dp) :: edddl(5,5,5,4), eddll(5,5,4,4), edldl(5,4,5,4), edlll(5,4,4,4)

! This subroutine manages computation of the two-=electron integrals 
! To do this, we make several classes of shell pairs and pairs of pairs, based on:
!  * their angular momentum s, p, d, l
!  * whether the bra or ket pairs are 'diagonal' i.e. involve the same shell
!  * whether the bra and ket pair and the same
! Each unique integral is then stored based on its ijkl index. Duplicates are avoided as much as possible, but clearly 
!   they will simply overlap each other, so it does not really matter
! Order of funcs, bra/kets: s, p, d, l, |ss), |sp), |sd), |sl), |pp), |pd), |pl), |dd), |dl), |ll)

! First work out how many TEIs of each type there may be, and create arrays to store them
allocate (i2e(nb4),stat=ierr)
if (ierr.ne.0) then
    write (*,*) "Allocate i2e failed"
    stop
end if
i2e=0._dp
schw=0._dp

! First calculate the 'diagonal' or 'Schwartz' integrals that have the same bra and ket - but different shells within these, i.e. (ab|ab)
! In due course I'll save these for Schwartz screening, but not in this version

! By shell pair type.
do ii=1,nss 
  sh1=ssprs(1,ii)
  ao1=aosh(sh1)
  sh2=ssprs(2,ii)
  ao2=aosh(sh2)
  ab1=ind2(sh1,sh2)
  call i2e_ssss(sh1,sh2,sh1,sh2,essss)
  schw(ab1)=sqrt(abs(essss))
  if (abs(essss).gt.intthresh) then
    ij=ind2(ao1,ao2)
    ijkl=ioff(ij)+ij
    i2e(ijkl)=essss
  end if
end do

! Now do (sp| pairs, which yield 3 (ab|ab) integrals and 3 unique (ab|ac)
do ii=1,nsp
  sh1=spprs(1,ii)
  ao1=aosh(sh1)
  sh2=spprs(2,ii)   ! First problem. sh2 can be smaller than sh1
  ao2=aosh(sh2)
  call i2e_spsp(sh1,sh2,sh1,sh2,espsp)
  ab1=ind2(sh1,sh2)
  schw(ab1)=sqrt(maxval(abs(espsp)))
  ! Here by permutation, we know that espsp(i,j) must equal espsp(j,i) so only store one of them
  do i=1,3
    ij=ind2(ao1,ao2+i-1)
    do j=1,i
      if (abs(espsp(i,j)).gt.intthresh) then
        kl=ind2(ao1,ao2+j-1)
        if (ij.gt.kl) then
          ijkl=ioff(ij)+kl
        else
          ijkl=ioff(kl)+ij
        end if
        i2e(ijkl)=espsp(i,j)
      end if
    end do
  end do
end do

! Now do (sd| pairs
do ii=1,nsd
  sh1=sdprs(1,ii)
  ao1=aosh(sh1)
  sh2=sdprs(2,ii)   ! First problem. sh2 can be smaller than sh1
  ao2=aosh(sh2)
  call i2e_sdsd(sh1,sh2,sh1,sh2,esdsd)
  ab1=ind2(sh1,sh2)
  schw(ab1)=sqrt(maxval(abs(esdsd)))
  ! Here by permutation, we know that esdsd(i,j) must equal esdsd(j,i) so only store one of them
  do i=1,5
    ij=ind2(ao1,ao2+i-1)
    do j=1,i
      if (abs(esdsd(i,j)).gt.intthresh) then
        kl=ind2(ao1,ao2+j-1)
        if (ij.gt.kl) then
          ijkl=ioff(ij)+kl
        else
          ijkl=ioff(kl)+ij
        end if
        i2e(ijkl)=esdsd(i,j)
      end if
    end do
  end do
end do

! Now do (sl| pairs
do ii=1,nsl
  sh1=slprs(1,ii)
  ao1=aosh(sh1)
  sh2=slprs(2,ii)   ! First problem. sh2 can be smaller than sh1
  ao2=aosh(sh2)
  call i2e_slsl(sh1,sh2,sh1,sh2,eslsl)
  ab1=ind2(sh1,sh2)
  schw(ab1)=sqrt(maxval(abs(eslsl)))
  do i=1,4
    ij=ind2(ao1,ao2+i-1)
    do j=1,i
      if (abs(eslsl(i,j)).gt.intthresh) then
        kl=ind2(ao1,ao2+j-1)
        if (ij.gt.kl) then
          ijkl=ioff(ij)+kl
        else
          ijkl=ioff(kl)+ij
        end if
        i2e(ijkl)=eslsl(i,j)
      end if
    end do
  end do
end do

! Now do (pp| pairs
do ii=1,npp
  sh1=ppprs(1,ii)
  ao1=aosh(sh1)
  sh2=ppprs(2,ii)
  ao2=aosh(sh2)
  call i2e_pppp(sh1,sh2,sh1,sh2,epppp)
  ab1=ind2(sh1,sh2)
  schw(ab1)=sqrt(maxval(abs(epppp)))
  do i=1,3
    do j=1,3
      ij=ind2(ao1+i-1,ao2+j-1)
      do k=1,i
        do l=1,3
          if (abs(epppp(i,j,k,l)).gt.intthresh) then
            kl=ind2(ao1+k-1,ao2+l-1)
            if (ij.gt.kl) then
              ijkl=ioff(ij)+kl
            else
              ijkl=ioff(kl)+ij
            end if
            i2e(ijkl)=epppp(i,j,k,l)
          end if
        end do
      end do
    end do
  end do
end do

! Now do (pd| pairs
do ii=1,npd
  sh1=pdprs(1,ii)
  ao1=aosh(sh1)
  sh2=pdprs(2,ii)
  ao2=aosh(sh2)
  call i2e_pdpd(sh1,sh2,sh1,sh2,epdpd)
  ab1=ind2(sh1,sh2)
  schw(ab1)=sqrt(maxval(abs(epdpd)))
  do i=1,3
    do j=1,5
      ij=ind2(ao1+i-1,ao2+j-1)
      do k=1,3
        do l=1,5
          if (abs(epdpd(i,j,k,l)).gt.intthresh) then
            kl=ind2(ao1+k-1,ao2+l-1)
            if (ij.gt.kl) then
              ijkl=ioff(ij)+kl
            else
              ijkl=ioff(kl)+ij
            end if
            i2e(ijkl)=epdpd(i,j,k,l)
          end if
        end do
      end do
    end do
  end do
end do

! Now do (pl| pairs
do ii=1,npl
  sh1=plprs(1,ii)
  ao1=aosh(sh1)
  sh2=plprs(2,ii)   ! First problem. sh2 can be smaller than sh1
  ao2=aosh(sh2)
  call i2e_plpl(sh1,sh2,sh1,sh2,eplpl)
  ab1=ind2(sh1,sh2)
  schw(ab1)=sqrt(maxval(abs(eplpl)))
  do i=1,3
    do j=1,4
      ij=ind2(ao1+i-1,ao2+j-1)
      do k=1,i
        do l=1,4
          if (abs(eplpl(i,j,k,l)).gt.intthresh) then
            kl=ind2(ao1+k-1,ao2+l-1)
            if (ij.gt.kl) then
              ijkl=ioff(ij)+kl
            else
              ijkl=ioff(kl)+ij
            end if
            i2e(ijkl)=eplpl(i,j,k,l)
          end if
        end do
      end do
    end do
  end do
end do

! Now do (dd| pairs
do ii=1,ndd
  sh1=ddprs(1,ii)
  ao1=aosh(sh1)
  sh2=ddprs(2,ii)
  ao2=aosh(sh2)
  call i2e_dddd(sh1,sh2,sh1,sh2,edddd)
  ab1=ind2(sh1,sh2)
  schw(ab1)=sqrt(maxval(abs(edddd)))
  do i=1,5
    do j=1,5
      ij=ind2(ao1+i-1,ao2+j-1)
      do k=1,i
        do l=1,5
          if (abs(edddd(i,j,k,l)).gt.intthresh) then
            kl=ind2(ao1+k-1,ao2+l-1)
            if (ij.gt.kl) then
              ijkl=ioff(ij)+kl
            else
              ijkl=ioff(kl)+ij
            end if
            i2e(ijkl)=edddd(i,j,k,l)
          end if
        end do
      end do
    end do
  end do
end do

! Now do (dl| pairs
do ii=1,ndl
  sh1=dlprs(1,ii)
  ao1=aosh(sh1)
  sh2=dlprs(2,ii)
  ao2=aosh(sh2)
  call i2e_dldl(sh1,sh2,sh1,sh2,edldl)
  ab1=ind2(sh1,sh2)
  schw(ab1)=sqrt(maxval(abs(edldl)))
  do i=1,5
    do j=1,4
      ij=ind2(ao1+i-1,ao2+j-1)
      do k=1,5
        do l=1,4
          if (abs(edldl(i,j,k,l)).gt.intthresh) then
            kl=ind2(ao1+k-1,ao2+l-1)
            if (ij.gt.kl) then
              ijkl=ioff(ij)+kl
            else
              ijkl=ioff(kl)+ij
            end if
            i2e(ijkl)=edldl(i,j,k,l)
          end if
        end do
      end do
    end do
  end do
end do

! Now do (ll'| pairs
do ii=1,nll
  sh1=llprs(1,ii)
  ao1=aosh(sh1)
  sh2=llprs(2,ii)
  ao2=aosh(sh2)
  call i2e_llll(sh1,sh2,sh1,sh2,ellll)
  ab1=ind2(sh1,sh2)
  schw(ab1)=sqrt(maxval(abs(ellll)))
  do i=1,4
    do j=1,4
      ij=ind2(ao1+i-1,ao2+j-1)
      do k=1,i
        do l=1,4
          if (abs(ellll(i,j,k,l)).gt.intthresh) then
            kl=ind2(ao1+k-1,ao2+l-1)
            if (ij.gt.kl) then
              ijkl=ioff(ij)+kl
            else
              ijkl=ioff(kl)+ij
            end if
            i2e(ijkl)=ellll(i,j,k,l)
          end if
        end do
      end do
    end do
  end do
end do

!call time_checker(-1,"Two electron integrals - diagonal terms, timing:")



! ********************************************
! Now loop over off-diagonal shell pair pairs.
! Start with (ss|ss)
nschw=0
do ii=2,nss
  sh1=ssprs(1,ii)
  ao1=aosh(sh1)
  sh2=ssprs(2,ii)
  ao2=aosh(sh2)
  ij=ind2(ao1,ao2)
  ab1=ind2(sh1,sh2)
  do jj=1,ii-1
    sh3=ssprs(1,jj)
    ao3=aosh(sh3)
    sh4=ssprs(2,jj)
    ao4=aosh(sh4)
    ab2=ind2(sh3,sh4)
    if ((schw(ab1)*schw(ab2)).lt.schwthresh) then
      nschw=nschw+1
      cycle
    end if
    call i2e_ssss(sh1,sh2,sh3,sh4,essss)
    if (abs(essss).gt.intthresh) then
      kl=ind2(ao3,ao4)
      if (ij.gt.kl) then
        ijkl=ioff(ij)+kl
      else
        ijkl=ioff(kl)+ij
      end if
      i2e(ijkl)=essss
    end if
  end do
end do
write (9,'(A,I8,A,I8)') "Number of screened (ss|ss) integral sets:",nschw," out of:",nss*(nss+1)/2
!call time_checker(-1,"Two electron integrals - (ss|ss) terms, timing:")

nschw=0
! Then (ss|sp)
do ii=1,nss
  sh1=ssprs(1,ii)
  ao1=aosh(sh1)
  sh2=ssprs(2,ii)
  ao2=aosh(sh2)
  ab1=ind2(sh1,sh2)
  ij=ind2(ao1,ao2)  
  do jj=1,nsp
    sh3=spprs(1,jj)
    ao3=aosh(sh3)
    sh4=spprs(2,jj)
    ao4=aosh(sh4)
    ab2=ind2(sh3,sh4)
    if ((schw(ab1)*schw(ab2)).lt.schwthresh) then
      nschw=nschw+1
      cycle
    end if
    call i2e_sssp(sh1,sh2,sh3,sh4,esssp)
    do i=1,3
      if (abs(esssp(i)).gt.intthresh) then
        kl=ind2(ao3,ao4+i-1)
        if (ij.gt.kl) then
          ijkl=ioff(ij)+kl
        else
          ijkl=ioff(kl)+ij
        end if
        i2e(ijkl)=esssp(i)
      end if
    end do
  end do
end do
write (9,'(A,I8,A,I8)') "Number of screened (ss|sp) integral sets:",nschw," out of:",nss*nsp
!call time_checker(-1,"Two electron integrals - (ss|sp) terms, timing:")

nschw=0
! Then (ss|sd)
do ii=1,nss
  sh1=ssprs(1,ii)
  ao1=aosh(sh1)
  sh2=ssprs(2,ii)
  ao2=aosh(sh2)
  ab1=ind2(sh1,sh2)
  ij=ind2(ao1,ao2)  
  do jj=1,nsd
    sh3=sdprs(1,jj)
    ao3=aosh(sh3)
    sh4=sdprs(2,jj)
    ao4=aosh(sh4)
    ab2=ind2(sh3,sh4)
    if ((schw(ab1)*schw(ab2)).lt.schwthresh) then
      nschw=nschw+1
      cycle
    end if
    call i2e_sssd(sh1,sh2,sh3,sh4,esssd)
    do i=1,5
      if (abs(esssd(i)).gt.intthresh) then
        kl=ind2(ao3,ao4+i-1)
        if (ij.gt.kl) then
          ijkl=ioff(ij)+kl
        else
          ijkl=ioff(kl)+ij
        end if
        i2e(ijkl)=esssd(i)
      end if
    end do
  end do
end do
write (9,'(A,I8,A,I8)') "Number of screened (ss|sd) integral sets:",nschw," out of:",nss*nsd
!call time_checker(-1,"Two electron integrals - (ss|sd) terms, timing:")

! Then (ss|sl)
nschw=0
do ii=1,nss
  sh1=ssprs(1,ii)
  ao1=aosh(sh1)
  sh2=ssprs(2,ii)
  ao2=aosh(sh2)
  ab1=ind2(sh1,sh2)
  ij=ind2(ao1,ao2)  
  do jj=1,nsl
    sh3=slprs(1,jj)
    ao3=aosh(sh3)
    sh4=slprs(2,jj)
    ao4=aosh(sh4)
    ab2=ind2(sh3,sh4)
    if ((schw(ab1)*schw(ab2)).lt.schwthresh) then
      nschw=nschw+1
      cycle
    end if
    call i2e_sssl(sh1,sh2,sh3,sh4,esssl)
    do i=1,4
      if (abs(esssl(i)).gt.intthresh) then
        kl=ind2(ao3,ao4+i-1)
        if (ij.gt.kl) then
          ijkl=ioff(ij)+kl
        else
          ijkl=ioff(kl)+ij
        end if
        i2e(ijkl)=esssl(i)
      end if
    end do
  end do
end do
write (9,'(A,I8,A,I8)') "Number of screened (ss|sl) integral sets:",nschw," out of:",nss*nsl
!call time_checker(-1,"Two electron integrals - (ss|sl) terms, timing:")

! Then (ss|pp), which can be (ab|cd) or (aa|bc) or (ab|cc) or (aa|cc)
nschw=0
do ii=1,nss
  sh1=ssprs(1,ii)
  ao1=aosh(sh1)
  sh2=ssprs(2,ii)
  ao2=aosh(sh2)
  ij=ind2(ao1,ao2)
  ab1=ind2(sh1,sh2)
  do jj=1,npp
    sh3=ppprs(1,jj)
    ao3=aosh(sh3)
    sh4=ppprs(2,jj)
    ao4=aosh(sh4)
    ab2=ind2(sh3,sh4)
    if ((schw(ab1)*schw(ab2)).lt.schwthresh) then
      nschw=nschw+1
      cycle
    end if
    call i2e_sspp(sh1,sh2,sh3,sh4,esspp)
    do i=1,3
      do j=1,3
        if (abs(esspp(i,j)).gt.intthresh) then
          kl=ind2(ao3+i-1,ao4+j-1)
          if (ij.gt.kl) then
            ijkl=ioff(ij)+kl
          else
            ijkl=ioff(kl)+ij
          end if
          i2e(ijkl)=esspp(i,j)
        end if
      end do
    end do
  end do
end do
write (9,'(A,I8,A,I8)') "Number of screened (ss|pp) integral sets:",nschw," out of:",nss*npp

! Then (ss|pd)
nschw=0
do ii=1,nss
  sh1=ssprs(1,ii)
  ao1=aosh(sh1)
  sh2=ssprs(2,ii)
  ao2=aosh(sh2)
  ij=ind2(ao1,ao2)
  ab1=ind2(sh1,sh2)
  do jj=1,npd
    sh3=pdprs(1,jj)
    ao3=aosh(sh3)
    sh4=pdprs(2,jj)
    ao4=aosh(sh4)
    ab2=ind2(sh3,sh4)
    if ((schw(ab1)*schw(ab2)).lt.schwthresh) then
      nschw=nschw+1
      cycle
    end if
    call i2e_sspd(sh1,sh2,sh3,sh4,esspd)
    do i=1,3
      do j=1,5
        if (abs(esspd(i,j)).gt.intthresh) then
          kl=ind2(ao3+i-1,ao4+j-1)
          if (ij.gt.kl) then
            ijkl=ioff(ij)+kl
          else
            ijkl=ioff(kl)+ij
          end if
          i2e(ijkl)=esspd(i,j)
        end if
      end do
    end do
  end do
end do
write (9,'(A,I8,A,I8)') "Number of screened (ss|pd) integral sets:",nschw," out of:",nss*npd

! Then (ss|pl)
nschw=0
do ii=1,nss
  sh1=ssprs(1,ii)
  ao1=aosh(sh1)
  sh2=ssprs(2,ii)
  ao2=aosh(sh2)
  ij=ind2(ao1,ao2)
  ab1=ind2(sh1,sh2)
  do jj=1,npl
    sh3=plprs(1,jj)
    ao3=aosh(sh3)
    sh4=plprs(2,jj)
    ao4=aosh(sh4)
    ab2=ind2(sh3,sh4)
    if ((schw(ab1)*schw(ab2)).lt.schwthresh) then
      nschw=nschw+1
      cycle
    end if
    call i2e_sspl(sh1,sh2,sh3,sh4,esspl)
    do i=1,3
      do j=1,4
        if (abs(esspl(i,j)).gt.intthresh) then
          kl=ind2(ao3+i-1,ao4+j-1)
          if (ij.gt.kl) then
            ijkl=ioff(ij)+kl
          else
            ijkl=ioff(kl)+ij
          end if
          i2e(ijkl)=esspl(i,j)
        end if
      end do
    end do
  end do
end do
write (9,'(A,I8,A,I8)') "Number of screened (ss|pl) integral sets:",nschw," out of:",nss*npl

! Then (ss|dd)
nschw=0
do ii=1,nss
  sh1=ssprs(1,ii)
  ao1=aosh(sh1)
  sh2=ssprs(2,ii)
  ao2=aosh(sh2)
  ij=ind2(ao1,ao2)
  ab1=ind2(sh1,sh2)
  do jj=1,ndd
    sh3=ddprs(1,jj)
    ao3=aosh(sh3)
    sh4=ddprs(2,jj)
    ao4=aosh(sh4)
    ab2=ind2(sh3,sh4)
    if ((schw(ab1)*schw(ab2)).lt.schwthresh) then
      nschw=nschw+1
      cycle
    end if
    call i2e_ssdd(sh1,sh2,sh3,sh4,essdd)
    do i=1,5
      do j=1,5
        if (abs(essdd(i,j)).gt.intthresh) then
          kl=ind2(ao3+i-1,ao4+j-1)
          if (ij.gt.kl) then
            ijkl=ioff(ij)+kl
          else
            ijkl=ioff(kl)+ij
          end if
          i2e(ijkl)=essdd(i,j)
        end if
      end do
    end do
  end do
end do
write (9,'(A,I8,A,I8)') "Number of screened (ss|dd) integral sets:",nschw," out of:",nss*ndd
!call time_checker(-1,"Two electron integrals - (ss|dd) terms, timing:")

! Then (ss|dl)
nschw=0
do ii=1,nss
  sh1=ssprs(1,ii)
  ao1=aosh(sh1)
  sh2=ssprs(2,ii)
  ao2=aosh(sh2)
  ij=ind2(ao1,ao2)
  ab1=ind2(sh1,sh2)
  do jj=1,ndl
    sh3=dlprs(1,jj)
    ao3=aosh(sh3)
    sh4=dlprs(2,jj)
    ao4=aosh(sh4)
    ab2=ind2(sh3,sh4)
    if ((schw(ab1)*schw(ab2)).lt.schwthresh) then
      nschw=nschw+1
      cycle
    end if
    call i2e_ssdl(sh1,sh2,sh3,sh4,essdl)
    do i=1,5
      do j=1,4
        if (abs(essdl(i,j)).gt.intthresh) then
          kl=ind2(ao3+i-1,ao4+j-1)
          if (ij.gt.kl) then
            ijkl=ioff(ij)+kl
          else
            ijkl=ioff(kl)+ij
          end if
          i2e(ijkl)=essdl(i,j)
        end if
      end do
    end do
  end do
end do
write (9,'(A,I8,A,I8)') "Number of screened (ss|dl) integral sets:",nschw," out of:",nss*ndl
!call time_checker(-1,"Two electron integrals - (ss|dl) terms, timing:")

! Then (ss|ll)
nschw=0
do ii=1,nss
  sh1=ssprs(1,ii)
  ao1=aosh(sh1)
  sh2=ssprs(2,ii)
  ao2=aosh(sh2)
  ij=ind2(ao1,ao2)
  ab1=ind2(sh1,sh2)
  do jj=1,nll
    sh3=llprs(1,jj)
    ao3=aosh(sh3)
    sh4=llprs(2,jj)
    ao4=aosh(sh4)
    ab2=ind2(sh3,sh4)
    if ((schw(ab1)*schw(ab2)).lt.schwthresh) then
      nschw=nschw+1
      cycle
    end if
    call i2e_ssll(sh1,sh2,sh3,sh4,essll)
    do i=1,4
      do j=1,4
        if (abs(essll(i,j)).gt.intthresh) then
          kl=ind2(ao3+i-1,ao4+j-1)
          if (ij.gt.kl) then
            ijkl=ioff(ij)+kl
          else
            ijkl=ioff(kl)+ij
          end if
          i2e(ijkl)=essll(i,j)
        end if
      end do
    end do
  end do
end do
write (9,'(A,I8,A,I8)') "Number of screened (ss|ll) integral sets:",nschw," out of:",nss*nll
!call time_checker(-1,"Two electron integrals - (ss|ll) terms, timing:")

! Then (sp|sp)
nschw=0
do ii=2,nsp
  sh1=spprs(1,ii)
  ao1=aosh(sh1)
  sh2=spprs(2,ii)
  ao2=aosh(sh2)
  ab1=ind2(sh1,sh2)
  do jj=1,ii-1
    sh3=spprs(1,jj)
    ao3=aosh(sh3)
    sh4=spprs(2,jj)
    ao4=aosh(sh4)
    ab2=ind2(sh3,sh4)
    if ((schw(ab1)*schw(ab2)).lt.schwthresh) then
      nschw=nschw+1
      cycle
    end if
    atc=atsh(sh3) ; atd=atsh(sh4)
    call i2e_spsp(sh1,sh2,sh3,sh4,espsp)
    do i=1,3
      do j=1,3
        if (abs(espsp(i,j)).gt.intthresh) then
          ij=ind2(ao1,ao2+i-1)
          kl=ind2(ao3,ao4+j-1)
          if (ij.gt.kl) then
            ijkl=ioff(ij)+kl
          else
            ijkl=ioff(kl)+ij
          end if
          i2e(ijkl)=espsp(i,j)
        end if
      end do
    end do
  end do
end do
write (9,'(A,I8,A,I8)') "Number of screened (sp|sp) integral sets:",nschw," out of:",nsp*(nsp+1)/2

! Then (sp|sd)
nschw=0
do ii=1,nsp
  sh1=spprs(1,ii)
  ao1=aosh(sh1)
  sh2=spprs(2,ii)
  ao2=aosh(sh2)
  ab1=ind2(sh1,sh2)
  do jj=1,nsd
    sh3=sdprs(1,jj)
    ao3=aosh(sh3)
    sh4=sdprs(2,jj)
    ao4=aosh(sh4)
    ab2=ind2(sh3,sh4)
    if ((schw(ab1)*schw(ab2)).lt.schwthresh) then
      nschw=nschw+1
      cycle
    end if
    call i2e_spsd(sh1,sh2,sh3,sh4,espsd)
    do i=1,3
      do j=1,5
        if (abs(espsd(i,j)).gt.intthresh) then
          ij=ind2(ao1,ao2+i-1)
          kl=ind2(ao3,ao4+j-1)
          if (ij.gt.kl) then
            ijkl=ioff(ij)+kl
          else
            ijkl=ioff(kl)+ij
          end if
          i2e(ijkl)=espsd(i,j)
        end if
      end do
    end do
  end do
end do
write (9,'(A,I8,A,I8)') "Number of screened (sp|sd) integral sets:",nschw," out of:",nsp*nsd

! Then (sp|sl)
nschw=0
do ii=1,nsp
  sh1=spprs(1,ii)
  ao1=aosh(sh1)
  sh2=spprs(2,ii)
  ao2=aosh(sh2)
  ab1=ind2(sh1,sh2)
  do jj=1,nsl
    sh3=slprs(1,jj)
    ao3=aosh(sh3)
    sh4=slprs(2,jj)
    ao4=aosh(sh4)
    ab2=ind2(sh3,sh4)
    if ((schw(ab1)*schw(ab2)).lt.schwthresh) then
      nschw=nschw+1
      cycle
    end if
    call i2e_spsl(sh1,sh2,sh3,sh4,espsl)
    do i=1,3
      do j=1,4
        if (abs(espsl(i,j)).gt.intthresh) then
          ij=ind2(ao1,ao2+i-1)
          kl=ind2(ao3,ao4+j-1)
          if (ij.gt.kl) then
            ijkl=ioff(ij)+kl
          else
            ijkl=ioff(kl)+ij
          end if
          i2e(ijkl)=espsl(i,j)
        end if
      end do
    end do
  end do
end do
write (9,'(A,I8,A,I8)') "Number of screened (sp|sl) integral sets:",nschw," out of:",nsp*nsl

! Then (sp|pp)
nschw=0
do ii=1,nsp
  sh1=spprs(1,ii)
  ao1=aosh(sh1)
  sh2=spprs(2,ii)
  ao2=aosh(sh2)
  ab1=ind2(sh1,sh2)
  do jj=1,npp
    sh3=ppprs(1,jj)
    ao3=aosh(sh3)
    sh4=ppprs(2,jj)
    ao4=aosh(sh4)
    ab2=ind2(sh3,sh4)
    if ((schw(ab1)*schw(ab2)).lt.schwthresh) then
      nschw=nschw+1
      cycle
    end if
    call i2e_sppp(sh1,sh2,sh3,sh4,esppp)
    do l=1,3
      do k=1,3
        kl=ind2(ao3+k-1,ao4+l-1)
        do j=1,3
          if (abs(esppp(j,k,l)).gt.intthresh) then 
            ij=ind2(ao1,ao2+j-1)
            if (ij.gt.kl) then
              ijkl=ioff(ij)+kl
            else
              ijkl=ioff(kl)+ij
            end if
            i2e(ijkl)=esppp(j,k,l)
          end if
        end do
      end do
    end do
  end do
end do
write (9,'(A,I8,A,I8)') "Number of screened (sp|pp) integral sets:",nschw," out of:",nsp*npp

! Then (sp|pd)
nschw=0
do ii=1,nsp
  sh1=spprs(1,ii)
  ao1=aosh(sh1)
  sh2=spprs(2,ii)
  ao2=aosh(sh2)
  ab1=ind2(sh1,sh2)
  do jj=1,npd
    sh3=pdprs(1,jj)
    ao3=aosh(sh3)
    sh4=pdprs(2,jj)
    ao4=aosh(sh4)
    ab2=ind2(sh3,sh4)
    if ((schw(ab1)*schw(ab2)).lt.schwthresh) then
      nschw=nschw+1
      cycle
    end if
    call i2e_sppd(sh1,sh2,sh3,sh4,esppd)
    do l=1,5
      do k=1,3
        kl=ind2(ao3+k-1,ao4+l-1)
        do j=1,3
          if (abs(esppd(j,k,l)).gt.intthresh) then 
            ij=ind2(ao1,ao2+j-1)
            if (ij.gt.kl) then
              ijkl=ioff(ij)+kl
            else
              ijkl=ioff(kl)+ij
            end if
            i2e(ijkl)=esppd(j,k,l)
          end if
        end do
      end do
    end do
  end do
end do
write (9,'(A,I8,A,I8)') "Number of screened (sp|pd) integral sets:",nschw," out of:",nsp*npd

! Then (sp|pl)
nschw=0
do ii=1,nsp
  sh1=spprs(1,ii)
  ao1=aosh(sh1)
  sh2=spprs(2,ii)
  ao2=aosh(sh2)
  ab1=ind2(sh1,sh2)
  do jj=1,npl
    sh3=plprs(1,jj)
    ao3=aosh(sh3)
    sh4=plprs(2,jj)
    ao4=aosh(sh4)
    ab2=ind2(sh3,sh4)
    if ((schw(ab1)*schw(ab2)).lt.schwthresh) then
      nschw=nschw+1
      cycle
    end if
    call i2e_sppl(sh1,sh2,sh3,sh4,esppl)
    do l=1,4
      do k=1,3
        kl=ind2(ao3+k-1,ao4+l-1)
        do j=1,3
          if (abs(esppl(j,k,l)).gt.intthresh) then 
            ij=ind2(ao1,ao2+j-1)
            if (ij.gt.kl) then
              ijkl=ioff(ij)+kl
            else
              ijkl=ioff(kl)+ij
            end if
            i2e(ijkl)=esppl(j,k,l)
          end if
        end do
      end do
    end do
  end do
end do
write (9,'(A,I8,A,I8)') "Number of screened (sp|pl) integral sets:",nschw," out of:",nsp*npl

! Then (sp|dd)
nschw=0
do ii=1,nsp
  sh1=spprs(1,ii)
  ao1=aosh(sh1)
  sh2=spprs(2,ii)
  ao2=aosh(sh2)
  ab1=ind2(sh1,sh2)
  do jj=1,ndd
    sh3=ddprs(1,jj)
    ao3=aosh(sh3)
    sh4=ddprs(2,jj)
    ao4=aosh(sh4)
    ab2=ind2(sh3,sh4)
    if ((schw(ab1)*schw(ab2)).lt.schwthresh) then
      nschw=nschw+1
      cycle
    end if
    call i2e_spdd(sh1,sh2,sh3,sh4,espdd)
    do l=1,5
      do k=1,5
        kl=ind2(ao3+k-1,ao4+l-1)
        do j=1,3
          if (abs(espdd(j,k,l)).gt.intthresh) then 
            ij=ind2(ao1,ao2+j-1)
            if (ij.gt.kl) then
              ijkl=ioff(ij)+kl
            else
              ijkl=ioff(kl)+ij
            end if
            i2e(ijkl)=espdd(j,k,l)
          end if
        end do
      end do
    end do
  end do
end do
write (9,'(A,I8,A,I8)') "Number of screened (sp|dd) integral sets:",nschw," out of:",nsp*ndd

! Then (sp|dl)
!nschw=0
!do ii=1,nsp
!  sh1=spprs(1,ii)
!  ao1=aosh(sh1)
!  sh2=spprs(2,ii)
!  ao2=aosh(sh2)
!  ab1=ind2(sh1,sh2)
!  do jj=1,ndl
!    sh3=dlprs(1,jj)
!    ao3=aosh(sh3)
!    sh4=dlprs(2,jj)
!    ao4=aosh(sh4)
!    ab2=ind2(sh3,sh4)
!    if ((schw(ab1)*schw(ab2)).lt.schwthresh) then
!      nschw=nschw+1
!      cycle
!    end if
!    call i2e_spdl(sh1,sh2,sh3,sh4,espdl)
!    do l=1,4
!      do k=1,5
!        kl=ind2(ao3+k-1,ao4+l-1)
!        do j=1,3
!          if (abs(espdl(j,k,l)).gt.intthresh) then 
!            ij=ind2(ao1,ao2+j-1)
!            ijkl=ind2(ij,kl)
!            i2e(ijkl)=espdl(j,k,l)
!          end if
!        end do
!      end do
!    end do
!  end do
!end do
!write (9,'(A,I8)') "Number of screened (sp|dl) integrals:",nschw

! Then (sp|ll)
nschw=0
do ii=1,nsp
  sh1=spprs(1,ii)
  ao1=aosh(sh1)
  sh2=spprs(2,ii)
  ao2=aosh(sh2)
  ab1=ind2(sh1,sh2)
  do jj=1,nll
    sh3=llprs(1,jj)
    ao3=aosh(sh3)
    sh4=llprs(2,jj)
    ao4=aosh(sh4)
    ab2=ind2(sh3,sh4)
    if ((schw(ab1)*schw(ab2)).lt.schwthresh) then
      nschw=nschw+1
      cycle
    end if
    call i2e_spll(sh1,sh2,sh3,sh4,espll)
    do l=1,4
      do k=1,4
        kl=ind2(ao3+k-1,ao4+l-1)
        do j=1,3
          if (abs(espll(j,k,l)).gt.intthresh) then 
            ij=ind2(ao1,ao2+j-1)
            if (ij.gt.kl) then
              ijkl=ioff(ij)+kl
            else
              ijkl=ioff(kl)+ij
            end if
            i2e(ijkl)=espll(j,k,l)
          end if
        end do
      end do
    end do
  end do
end do
write (9,'(A,I8,A,I8)') "Number of screened (sp|ll) integral sets:",nschw," out of:",nsp*nll

! Then (sd|sd)
nschw=0
do ii=2,nsd
  sh1=sdprs(1,ii)
  ao1=aosh(sh1)
  sh2=sdprs(2,ii)
  ao2=aosh(sh2)
  ab1=ind2(sh1,sh2)
  do jj=1,ii
    sh3=sdprs(1,jj)
    ao3=aosh(sh3)
    sh4=sdprs(2,jj)
    ao4=aosh(sh4)
    ab2=ind2(sh3,sh4)
    if ((schw(ab1)*schw(ab2)).lt.schwthresh) then
      nschw=nschw+1
      cycle
    end if
    call i2e_sdsd(sh1,sh2,sh3,sh4,esdsd)
    do l=1,5
      kl=ind2(ao3,ao4+l-1)
      do j=1,5
        if (abs(esdsd(j,l)).gt.intthresh) then 
          ij=ind2(ao1,ao2+j-1)
          if (ij.gt.kl) then
            ijkl=ioff(ij)+kl
          else
            ijkl=ioff(kl)+ij
          end if
          i2e(ijkl)=esdsd(j,l)
        end if
      end do
    end do
  end do
end do
write (9,'(A,I8,A,I8)') "Number of screened (sd|sd) integral sets:",nschw," out of:",nsd*(nsd+1)/2

! Then (sd|sl)
nschw=0
do ii=1,nsd
  sh1=sdprs(1,ii)
  ao1=aosh(sh1)
  sh2=sdprs(2,ii)
  ao2=aosh(sh2)
  ab1=ind2(sh1,sh2)
  do jj=1,nsl
    sh3=slprs(1,jj)
    ao3=aosh(sh3)
    sh4=slprs(2,jj)
    ao4=aosh(sh4)
    ab2=ind2(sh3,sh4)
    if ((schw(ab1)*schw(ab2)).lt.schwthresh) then
      nschw=nschw+1
      cycle
    end if
    call i2e_sdsl(sh1,sh2,sh3,sh4,esdsl)
    do l=1,4
      kl=ind2(ao3,ao4+l-1)
      do j=1,5
        if (abs(esdsl(j,l)).gt.intthresh) then 
          ij=ind2(ao1,ao2+j-1)
          if (ij.gt.kl) then
            ijkl=ioff(ij)+kl
          else
            ijkl=ioff(kl)+ij
          end if
          i2e(ijkl)=esdsl(j,l)
        end if
      end do
    end do
  end do
end do
write (9,'(A,I8,A,I8)') "Number of screened (sd|sl) integral sets:",nschw," out of:",nsd*nsl
!call time_checker(-1,"Two electron integrals - (sd|sl) terms, timing:")

! Then (sd|pp)
nschw=0
do ii=1,nsd
  sh1=sdprs(1,ii)
  ao1=aosh(sh1)
  sh2=sdprs(2,ii)
  ao2=aosh(sh2)
  ab1=ind2(sh1,sh2)
  do jj=1,npp
    sh3=ppprs(1,jj)
    ao3=aosh(sh3)
    sh4=ppprs(2,jj)
    ao4=aosh(sh4)
    ab2=ind2(sh3,sh4)
    if ((schw(ab1)*schw(ab2)).lt.schwthresh) then
      nschw=nschw+1
      cycle
    end if
    call i2e_sdpp(sh1,sh2,sh3,sh4,esdpp)
    do l=1,3
      do k=1,3
        kl=ind2(ao3+k-1,ao4+l-1)
        do j=1,5
          if (abs(esdpp(j,k,l)).gt.intthresh) then 
            ij=ind2(ao1,ao2+j-1)
            if (ij.gt.kl) then
              ijkl=ioff(ij)+kl
            else
              ijkl=ioff(kl)+ij
            end if
            i2e(ijkl)=esdpp(j,k,l)
          end if
        end do
      end do
    end do
  end do
end do

! Then (sd|pd)
nschw=0
do ii=1,nsd
  sh1=sdprs(1,ii)
  ao1=aosh(sh1)
  sh2=sdprs(2,ii)
  ao2=aosh(sh2)
  ab1=ind2(sh1,sh2)
  do jj=1,npd
    sh3=pdprs(1,jj)
    ao3=aosh(sh3)
    sh4=pdprs(2,jj)
    ao4=aosh(sh4)
    ab2=ind2(sh3,sh4)
    if ((schw(ab1)*schw(ab2)).lt.schwthresh) then
      nschw=nschw+1
      cycle
    end if
    call i2e_sdpd(sh1,sh2,sh3,sh4,esdpd)
    do l=1,5
      do k=1,3
        kl=ind2(ao3+k-1,ao4+l-1)
        do j=1,5
          if (abs(esdpd(j,k,l)).gt.intthresh) then 
            ij=ind2(ao1,ao2+j-1)
            if (ij.gt.kl) then
              ijkl=ioff(ij)+kl
            else
              ijkl=ioff(kl)+ij
            end if
            i2e(ijkl)=esdpd(j,k,l)
          end if
        end do
      end do
    end do
  end do
end do
write (9,'(A,I8,A,I8)') "Number of screened (sd|pd) integral sets:",nschw," out of:",nsd*npd

! Then (sd|pl)
!nschw=0
!do ii=1,nsd
!  sh1=sdprs(1,ii)
!  ao1=aosh(sh1)
!  sh2=sdprs(2,ii)
!  ao2=aosh(sh2)
!  ab1=ind2(sh1,sh2)
!  do jj=1,npl
!    sh3=plprs(1,jj)
!    ao3=aosh(sh3)
!    sh4=plprs(2,jj)
!    ao4=aosh(sh4)
!    ab2=ind2(sh3,sh4)
!    if ((schw(ab1)*schw(ab2)).lt.schwthresh) then
!      nschw=nschw+1
!      cycle
!    end if
!    call i2e_sdpl(sh1,sh2,sh3,sh4,esdpl)
!    do l=1,4
!      do k=1,3
!        kl=ind2(ao3+k-1,ao4+l-1)
!        do j=1,5
!          if (abs(esdpl(j,k,l)).gt.intthresh) then 
!            ij=ind2(ao1,ao2+j-1)
!            ijkl=ind2(ij,kl)
!            i2e(ijkl)=esdpl(j,k,l)
!          end if
!        end do
!      end do
!    end do
!  end do
!end do
!write (9,'(A,I8)') "Number of screened (sd|pl) integrals:",nschw

! Then (sd|dd)
nschw=0
do ii=1,nsd
  sh1=sdprs(1,ii)
  ao1=aosh(sh1)
  sh2=sdprs(2,ii)
  ao2=aosh(sh2)
  ab1=ind2(sh1,sh2)
  do jj=1,ndd
    sh3=ddprs(1,jj)
    ao3=aosh(sh3)
    sh4=ddprs(2,jj)
    ao4=aosh(sh4)
    ab2=ind2(sh3,sh4)
    if ((schw(ab1)*schw(ab2)).lt.schwthresh) then
      nschw=nschw+1
      cycle
    end if
    call i2e_sddd(sh1,sh2,sh3,sh4,esddd)
    do l=1,5
      do k=1,5
        kl=ind2(ao3+k-1,ao4+l-1)
        do j=1,5
          if (abs(esddd(j,k,l)).gt.intthresh) then 
            ij=ind2(ao1,ao2+j-1)
            if (ij.gt.kl) then
              ijkl=ioff(ij)+kl
            else
              ijkl=ioff(kl)+ij
            end if
            i2e(ijkl)=esddd(j,k,l)
          end if
        end do
      end do
    end do
  end do
end do
write (9,'(A,I8,A,I8)') "Number of screened (sd|dd) integral sets:",nschw," out of:",nsd*ndd
!call time_checker(-1,"Two electron integrals - (sd|dd) terms, timing:")

! Then (sd|dl)
nschw=0
do ii=1,nsd
  sh1=sdprs(1,ii)
  ao1=aosh(sh1)
  sh2=sdprs(2,ii)
  ao2=aosh(sh2)
  ab1=ind2(sh1,sh2)
  do jj=1,ndl
    sh3=dlprs(1,jj)
    ao3=aosh(sh3)
    sh4=dlprs(2,jj)
    ao4=aosh(sh4)
    ab2=ind2(sh3,sh4)
    if ((schw(ab1)*schw(ab2)).lt.schwthresh) then
      nschw=nschw+1
      cycle
    end if
    call i2e_sddl(sh1,sh2,sh3,sh4,esddl)
    do l=1,4
      do k=1,5
        kl=ind2(ao3+k-1,ao4+l-1)
        do j=1,5
          if (abs(esddl(j,k,l)).gt.intthresh) then 
            ij=ind2(ao1,ao2+j-1)
            if (ij.gt.kl) then
              ijkl=ioff(ij)+kl
            else
              ijkl=ioff(kl)+ij
            end if
            i2e(ijkl)=esddl(j,k,l)
          end if
        end do
      end do
    end do
  end do
end do
write (9,'(A,I8,A,I8)') "Number of screened (sd|dl) integral sets:",nschw," out of:",nsd*ndl
!call time_checker(-1,"Two electron integrals - (sd|dl) terms, timing:")

! Then (sd|ll)
nschw=0
do ii=1,nsd
  sh1=sdprs(1,ii)
  ao1=aosh(sh1)
  sh2=sdprs(2,ii)
  ao2=aosh(sh2)
  ab1=ind2(sh1,sh2)
  do jj=1,nll
    sh3=llprs(1,jj)
    ao3=aosh(sh3)
    sh4=llprs(2,jj)
    ao4=aosh(sh4)
    ab2=ind2(sh3,sh4)
    if ((schw(ab1)*schw(ab2)).lt.schwthresh) then
      nschw=nschw+1
      cycle
    end if
    call i2e_sdll(sh1,sh2,sh3,sh4,esdll)
    do l=1,4
      do k=1,4
        kl=ind2(ao3+k-1,ao4+l-1)
        do j=1,5
          if (abs(esdll(j,k,l)).gt.intthresh) then 
            ij=ind2(ao1,ao2+j-1)
            if (ij.gt.kl) then
              ijkl=ioff(ij)+kl
            else
              ijkl=ioff(kl)+ij
            end if
            i2e(ijkl)=esdll(j,k,l)
          end if
        end do
      end do
    end do
  end do
end do
write (9,'(A,I8,A,I8)') "Number of screened (sd|ll) integral sets:",nschw," out of:",nsd*nll
!call time_checker(-1,"Two electron integrals - (sd|ll) terms, timing:")

! Then (sl|sl)
nschw=0
n1=0;n2=0;n3=0;n4=0
do ii=2,nsl
  sh1=slprs(1,ii)
  ao1=aosh(sh1)
  sh2=slprs(2,ii)
  ao2=aosh(sh2)
  ab1=ind2(sh1,sh2)
  do jj=1,ii-1
    sh3=slprs(1,jj)
    ao3=aosh(sh3)
    sh4=slprs(2,jj)
    ao4=aosh(sh4)
    ab2=ind2(sh3,sh4)
    if ((schw(ab1)*schw(ab2)).lt.schwthresh) then
      nschw=nschw+1
      cycle
    end if
    call i2e_slsl(sh1,sh2,sh3,sh4,eslsl)
    do l=1,4
      kl=ind2(ao3,ao4+l-1)
      do j=1,4
        if (abs(eslsl(j,l)).gt.intthresh) then 
          ij=ind2(ao1,ao2+j-1)
          if (ij.gt.kl) then
            ijkl=ioff(ij)+kl
          else
            ijkl=ioff(kl)+ij
          end if
          i2e(ijkl)=eslsl(j,l)
        end if
      end do
    end do
  end do
end do
write (9,'(A,I8,A,I8)') "Number of screened (sl|sl) integral sets:",nschw," out of:",nsl*(nsl+1)/2
!call time_checker(-1,"Two electron integrals - (sl|sl) terms, timing:")

! Then (sl|pp)
nschw=0
do ii=1,nsl
  sh1=slprs(1,ii)
  ao1=aosh(sh1)
  sh2=slprs(2,ii)
  ao2=aosh(sh2)
  ab1=ind2(sh1,sh2)
  do jj=1,npp
    sh3=ppprs(1,jj)
    ao3=aosh(sh3)
    sh4=ppprs(2,jj)
    ao4=aosh(sh4)
    ab2=ind2(sh3,sh4)
    if ((schw(ab1)*schw(ab2)).lt.schwthresh) then
      nschw=nschw+1
      cycle
    end if
    call i2e_slpp(sh1,sh2,sh3,sh4,eslpp)
    do l=1,3
      do k=1,3
        kl=ind2(ao3+k-1,ao4+l-1)
        do j=1,4
          if (abs(eslpp(j,k,l)).gt.intthresh) then 
            ij=ind2(ao1,ao2+j-1)
            if (ij.gt.kl) then
              ijkl=ioff(ij)+kl
            else
              ijkl=ioff(kl)+ij
            end if
            i2e(ijkl)=eslpp(j,k,l)
          end if
        end do
      end do
    end do
  end do
end do
write (9,'(A,I8,A,I8)') "Number of screened (sl|pp) integral sets:",nschw," out of:",nsl*npp

! Then (sl|pd)
!nschw=0
!do ii=1,nsl
!  sh1=slprs(1,ii)
!  ao1=aosh(sh1)
!  sh2=slprs(2,ii)
!  ao2=aosh(sh2)
!  ab1=ind2(sh1,sh2)
!  do jj=1,npd
!    sh3=pdprs(1,jj)
!    ao3=aosh(sh3)
!    sh4=pdprs(2,jj)
!    ao4=aosh(sh4)
!    ab2=ind2(sh3,sh4)
!    if ((schw(ab1)*schw(ab2)).lt.schwthresh) then
!      nschw=nschw+1
!      cycle
!    end if
!    call i2e_slpd(sh1,sh2,sh3,sh4,eslpd)
!    do l=1,5
!      do k=1,3
!        kl=ind2(ao3+k-1,ao4+l-1)
!        do j=1,4
!          if (abs(eslpd(j,k,l)).gt.intthresh) then 
!            ij=ind2(ao1,ao2+j-1)
!            ijkl=ind2(ij,kl)
!            i2e(ijkl)=eslpd(j,k,l)
!          end if
!        end do
!      end do
!    end do
!  end do
!end do

! Then (sl|pl)
nschw=0
do ii=1,nsl
  sh1=slprs(1,ii)
  ao1=aosh(sh1)
  sh2=slprs(2,ii)
  ao2=aosh(sh2)
  ab1=ind2(sh1,sh2)
  do jj=1,npl
    sh3=plprs(1,jj)
    ao3=aosh(sh3)
    sh4=plprs(2,jj)
    ao4=aosh(sh4)
    ab2=ind2(sh3,sh4)
    if ((schw(ab1)*schw(ab2)).lt.schwthresh) then
      nschw=nschw+1
      cycle
    end if
    call i2e_slpl(sh1,sh2,sh3,sh4,eslpl)
    do l=1,4
      do k=1,3
        kl=ind2(ao3+k-1,ao4+l-1)
        do j=1,4
          if (abs(eslpl(j,k,l)).gt.intthresh) then 
            ij=ind2(ao1,ao2+j-1)
            if (ij.gt.kl) then
              ijkl=ioff(ij)+kl
            else
              ijkl=ioff(kl)+ij
            end if
            i2e(ijkl)=eslpl(j,k,l)
          end if
        end do
      end do
    end do
  end do
end do
write (9,'(A,I8,A,I8)') "Number of screened (sl|pl) integral sets:",nschw," out of:",nsl*npl

! Then (sl|dd)
nschw=0
do ii=1,nsl
  sh1=slprs(1,ii)
  ao1=aosh(sh1)
  sh2=slprs(2,ii)
  ao2=aosh(sh2)
  ab1=ind2(sh1,sh2)
  do jj=1,ndd
    sh3=ddprs(1,jj)
    ao3=aosh(sh3)
    sh4=ddprs(2,jj)
    ao4=aosh(sh4)
    ab2=ind2(sh3,sh4)
    if ((schw(ab1)*schw(ab2)).lt.schwthresh) then
      nschw=nschw+1
      cycle
    end if
    call i2e_sldd(sh1,sh2,sh3,sh4,esldd)
    do la=1,5
      do ka=1,5
        kl=ind2(ka+ao3-1,la+ao4-1)
        do ja=1,4
          if (abs(esldd(ja,ka,la)).gt.intthresh) then 
            ij=ind2(ao1,ao2+ja-1)
            if (ij.gt.kl) then
              ijkl=ioff(ij)+kl
            else
              ijkl=ioff(kl)+ij
            end if
            i2e(ijkl)=esldd(ja,ka,la)
          end if
        end do
      end do
    end do
  end do
end do
write (9,'(A,I8,A,I8)') "Number of screened (sl|dd) integral sets:",nschw," out of:",nsl*ndd
!call time_checker(-1,"Two electron integrals - (sl|dd) terms, timing:")

! Then, (sl|dl)
nschw=0
do ii=1,nsl
  sh1=slprs(1,ii)
  ao1=aosh(sh1)
  sh2=slprs(2,ii)
  ao2=aosh(sh2)
  ab1=ind2(sh1,sh2)
  do jj=1,ndl
    sh3=dlprs(1,jj)
    ao3=aosh(sh3)
    sh4=dlprs(2,jj)
    ao4=aosh(sh4)
    ab2=ind2(sh3,sh4)
    if ((schw(ab1)*schw(ab2)).lt.schwthresh) then
      nschw=nschw+1
      cycle
    end if
    call i2e_sldl(sh1,sh2,sh3,sh4,esldl)
    do la=1,4
      do ka=1,5
        kl=ind2(ka+ao3-1,la+ao4-1)
        do ja=1,4
          if (abs(esldl(ja,ka,la)).gt.intthresh) then 
            ij=ind2(ao1,ao2+ja-1)
            if (ij.gt.kl) then
              ijkl=ioff(ij)+kl
            else
              ijkl=ioff(kl)+ij
            end if
            i2e(ijkl)=esldl(ja,ka,la)
          end if
        end do
      end do
    end do
  end do
end do
write (9,'(A,I8,A,I8)') "Number of screened (sl|dl) integral sets:",nschw," out of:",nsl*ndl
!call time_checker(-1,"Two electron integrals - (sl|dl) terms, timing:")

! Then (sl|ll)
nschw=0
do ii=1,nsl
  sh1=slprs(1,ii)
  ao1=aosh(sh1)
  sh2=slprs(2,ii)
  ao2=aosh(sh2)
  ab1=ind2(sh1,sh2)
  do jj=1,nll
    sh3=llprs(1,jj)
    ao3=aosh(sh3)
    sh4=llprs(2,jj)
    ao4=aosh(sh4)
    ab2=ind2(sh3,sh4)
    if ((schw(ab1)*schw(ab2)).lt.schwthresh) then
      nschw=nschw+1
      cycle
    end if
    call i2e_slll(sh1,sh2,sh3,sh4,eslll)
    do la=1,4
      do ka=1,4
        kl=ind2(ka+ao3-1,la+ao4-1)
        do ja=1,4
          if (abs(eslll(ja,ka,la)).gt.intthresh) then 
            ij=ind2(ao1,ao2+ja-1)
            if (ij.gt.kl) then
              ijkl=ioff(ij)+kl
            else
              ijkl=ioff(kl)+ij
            end if
            i2e(ijkl)=eslll(ja,ka,la)
          end if
        end do
      end do
    end do
  end do
end do
write (9,'(A,I8,A,I8)') "Number of screened (sl|ll) integral sets:",nschw," out of:",nsl*nll
!call time_checker(-1,"Two electron integrals - (sl|ll) terms, timing:")

! Then (pp|pp)
nschw=0
do ii=2,npp
  sh1=ppprs(1,ii)
  ao1=aosh(sh1)
  sh2=ppprs(2,ii)
  ao2=aosh(sh2)
  ab1=ind2(sh1,sh2)
  do jj=1,ii-1
    sh3=ppprs(1,jj)
    ao3=aosh(sh3)
    sh4=ppprs(2,jj)
    ao4=aosh(sh4)
    ab2=ind2(sh3,sh4)
    if ((schw(ab1)*schw(ab2)).lt.schwthresh) then
      nschw=nschw+1
      cycle
    end if
    call i2e_pppp(sh1,sh2,sh3,sh4,epppp)
    do l=1,3
      do k=1,3
        kl=ind2(ao3+k-1,ao4+l-1)
        do j=1,3
          do i=1,3
            if (abs(epppp(i,j,k,l)).gt.intthresh) then
              ij=ind2(ao1+i-1,ao2+j-1)
              if (ij.gt.kl) then
                ijkl=ioff(ij)+kl
              else
                ijkl=ioff(kl)+ij
              end if
              i2e(ijkl)=epppp(i,j,k,l)
            end if
          end do
        end do
      end do
    end do
  end do
end do
write (9,'(A,I8,A,I8)') "Number of screened (pp|pp) integral sets:",nschw," out of:",npp*(npp+1)/2

! Then (pp|pd)
nschw=0
do ii=1,npp
  sh1=ppprs(1,ii)
  ao1=aosh(sh1)
  sh2=ppprs(2,ii)
  ao2=aosh(sh2)
  ab1=ind2(sh1,sh2)
  do jj=1,npd  
    sh3=pdprs(1,jj)
    ao3=aosh(sh3)
    sh4=pdprs(2,jj)
    ao4=aosh(sh4)
    ab2=ind2(sh3,sh4)
    if ((schw(ab1)*schw(ab2)).lt.schwthresh) then
      nschw=nschw+1
      cycle
    end if
    call i2e_pppd(sh1,sh2,sh3,sh4,epppd)
    do l=1,5
      do k=1,3
        kl=ind2(ao3+k-1,ao4+l-1)
        do j=1,3
          do i=1,3
            if (abs(epppd(i,j,k,l)).gt.intthresh) then
              ij=ind2(ao1+i-1,ao2+j-1)
              if (ij.gt.kl) then
                ijkl=ioff(ij)+kl
              else
                ijkl=ioff(kl)+ij
              end if
              i2e(ijkl)=epppd(i,j,k,l)
            end if
          end do
        end do
      end do
    end do
  end do
end do
write (9,'(A,I8,A,I8)') "Number of screened (pp|pd) integral sets:",nschw," out of:",npp*npd

! Then (pp|pl)
nschw=0
do ii=1,npp
  sh1=ppprs(1,ii)
  ao1=aosh(sh1)
  sh2=ppprs(2,ii)
  ao2=aosh(sh2)
  ab1=ind2(sh1,sh2)
  do jj=1,npl  
    sh3=plprs(1,jj)
    ao3=aosh(sh3)
    sh4=plprs(2,jj)
    ao4=aosh(sh4)
    ab2=ind2(sh3,sh4)
    if ((schw(ab1)*schw(ab2)).lt.schwthresh) then
      nschw=nschw+1
      cycle
    end if
    call i2e_pppl(sh1,sh2,sh3,sh4,epppl)
    do l=1,4
      do k=1,3
        kl=ind2(ao3+k-1,ao4+l-1)
        do j=1,3
          do i=1,3
            if (abs(epppl(i,j,k,l)).gt.intthresh) then
              ij=ind2(ao1+i-1,ao2+j-1)
              if (ij.gt.kl) then
                ijkl=ioff(ij)+kl
              else
                ijkl=ioff(kl)+ij
              end if
              i2e(ijkl)=epppl(i,j,k,l)
            end if
          end do
        end do
      end do
    end do
  end do
end do
write (9,'(A,I8,A,I8)') "Number of screened (pp|pl) integral sets:",nschw," out of:",npp*npl

! Then (pp|dd)
nschw=0
do ii=1,npp
  sh1=ppprs(1,ii)
  ao1=aosh(sh1)
  sh2=ppprs(2,ii)
  ao2=aosh(sh2)
  ab1=ind2(sh1,sh2)
  do jj=1,ndd  
    sh3=ddprs(1,jj)
    ao3=aosh(sh3)
    sh4=ddprs(2,jj)
    ao4=aosh(sh4)
    ab2=ind2(sh3,sh4)
    if ((schw(ab1)*schw(ab2)).lt.schwthresh) then
      nschw=nschw+1
      cycle
    end if
    call i2e_ppdd(sh1,sh2,sh3,sh4,eppdd)
    do l=1,5
      do k=1,5
        kl=ind2(ao3+k-1,ao4+l-1)
        do j=1,3
          do i=1,3
            if (abs(eppdd(i,j,k,l)).gt.intthresh) then
              ij=ind2(ao1+i-1,ao2+j-1)
              if (ij.gt.kl) then
                ijkl=ioff(ij)+kl
              else
                ijkl=ioff(kl)+ij
              end if
              i2e(ijkl)=eppdd(i,j,k,l)
            end if
          end do
        end do
      end do
    end do
  end do
end do
write (9,'(A,I8,A,I8)') "Number of screened (pp|dd) integral sets:",nschw," out of:",npp*ndd

! Then (pp|dl) which i'll skip

! Then (pp|ll)
nschw=0
do ii=1,npp
  sh1=ppprs(1,ii)
  ao1=aosh(sh1)
  sh2=ppprs(2,ii)
  ao2=aosh(sh2)
  ab1=ind2(sh1,sh2)
  do jj=1,nll  
    sh3=llprs(1,jj)
    ao3=aosh(sh3)
    sh4=llprs(2,jj)
    ao4=aosh(sh4)
    ab2=ind2(sh3,sh4)
    if ((schw(ab1)*schw(ab2)).lt.schwthresh) then
      nschw=nschw+1
      cycle
    end if
    call i2e_ppll(sh1,sh2,sh3,sh4,eppll)
    do l=1,4
      do k=1,4
        kl=ind2(ao3+k-1,ao4+l-1)
        do j=1,3
          do i=1,3
            if (abs(eppll(i,j,k,l)).gt.intthresh) then
              ij=ind2(ao1+i-1,ao2+j-1)
              if (ij.gt.kl) then
                ijkl=ioff(ij)+kl
              else
                ijkl=ioff(kl)+ij
              end if
              i2e(ijkl)=eppll(i,j,k,l)
            end if
          end do
        end do
      end do
    end do
  end do
end do
write (9,'(A,I8,A,I8)') "Number of screened (pp|ll) integral sets:",nschw," out of:",npp*nll

! Then (pd|pd)
nschw=0
do ii=2,npd
  sh1=pdprs(1,ii)
  ao1=aosh(sh1)
  sh2=pdprs(2,ii)
  ao2=aosh(sh2)
  ab1=ind2(sh1,sh2)
  do jj=1,ii-1
    sh3=pdprs(1,jj)
    ao3=aosh(sh3)
    sh4=pdprs(2,jj)
    ao4=aosh(sh4)
    ab2=ind2(sh3,sh4)
    if ((schw(ab1)*schw(ab2)).lt.schwthresh) then
      nschw=nschw+1
      cycle
    end if
    call i2e_pdpd(sh1,sh2,sh3,sh4,epdpd)
    do l=1,5
      do k=1,3
        kl=ind2(ao3+k-1,ao4+l-1)
        do j=1,5
          do i=1,3
            if (abs(epdpd(i,j,k,l)).gt.intthresh) then
              ij=ind2(ao1+i-1,ao2+j-1)
              if (ij.gt.kl) then
                ijkl=ioff(ij)+kl
              else
                ijkl=ioff(kl)+ij
              end if
              i2e(ijkl)=epdpd(i,j,k,l)
            end if
          end do
        end do
      end do
    end do
  end do
end do
write (9,'(A,I8,A,I8)') "Number of screened (pd|pd) integral sets:",nschw," out of:",npd*(npd+1)/2

! Gonna skip (pd|pl)

! Then (pd|dd)
nschw=0
do ii=1,npd
  sh1=pdprs(1,ii)
  ao1=aosh(sh1)
  sh2=pdprs(2,ii)
  ao2=aosh(sh2)
  ab1=ind2(sh1,sh2)
  do jj=1,ndd
    sh3=ddprs(1,jj)
    ao3=aosh(sh3)
    sh4=ddprs(2,jj)
    ao4=aosh(sh4)
    ab2=ind2(sh3,sh4)
    if ((schw(ab1)*schw(ab2)).lt.schwthresh) then
      nschw=nschw+1
      cycle
    end if
    call i2e_pddd(sh1,sh2,sh3,sh4,epddd)
    do l=1,5
      do k=1,5
        kl=ind2(ao3+k-1,ao4+l-1)
        do j=1,5
          do i=1,3
            if (abs(epddd(i,j,k,l)).gt.intthresh) then
              ij=ind2(ao1+i-1,ao2+j-1)
              if (ij.gt.kl) then
                ijkl=ioff(ij)+kl
              else
                ijkl=ioff(kl)+ij
              end if
              i2e(ijkl)=epddd(i,j,k,l)
            end if
          end do
        end do
      end do
    end do
  end do
end do
write (9,'(A,I8,A,I8)') "Number of screened (pd|dd) integral sets:",nschw," out of:",npd*ndd


! I will skip (pd|dl) and (pd|ll)

! Then (pl|pl)
nschw=0
do ii=2,npl
  sh1=plprs(1,ii)
  ao1=aosh(sh1)
  sh2=plprs(2,ii)
  ao2=aosh(sh2)
  ab1=ind2(sh1,sh2)
  do jj=1,ii-1
    sh3=plprs(1,jj)
    ao3=aosh(sh3)
    sh4=plprs(2,jj)
    ao4=aosh(sh4)
    ab2=ind2(sh3,sh4)
    if ((schw(ab1)*schw(ab2)).lt.schwthresh) then
      nschw=nschw+1
      cycle
    end if
    call i2e_plpl(sh1,sh2,sh3,sh4,eplpl)
    do l=1,4
      do k=1,3
        kl=ind2(ao3+k-1,ao4+l-1)
        do j=1,4
          do i=1,3
            if (abs(eplpl(i,j,k,l)).gt.intthresh) then
              ij=ind2(ao1+i-1,ao2+j-1)
              if (ij.gt.kl) then
                ijkl=ioff(ij)+kl
              else
                ijkl=ioff(kl)+ij
              end if
              i2e(ijkl)=eplpl(i,j,k,l)
            end if
          end do
        end do
      end do
    end do
  end do
end do
write (9,'(A,I8,A,I8)') "Number of screened (pl|pl) integral sets:",nschw," out of:",npl*(npl+1)/2

! Also skipping (pl|dd) and (pl|dl)

! Then (pl|ll)
nschw=0
do ii=1,npl
  sh1=plprs(1,ii)
  ao1=aosh(sh1)
  sh2=plprs(2,ii)
  ao2=aosh(sh2)
  ab1=ind2(sh1,sh2)
  do jj=1,nll  
    sh3=llprs(1,jj)
    ao3=aosh(sh3)
    sh4=llprs(2,jj)
    ao4=aosh(sh4)
    ab2=ind2(sh3,sh4)
    if ((schw(ab1)*schw(ab2)).lt.schwthresh) then
      nschw=nschw+1
      cycle
    end if
    call i2e_plll(sh1,sh2,sh3,sh4,eplll)
    do l=1,4
      do k=1,4
        kl=ind2(ao3+k-1,ao4+l-1)
        do j=1,4
          do i=1,3
            if (abs(eplll(i,j,k,l)).gt.intthresh) then
              ij=ind2(ao1+i-1,ao2+j-1)
              if (ij.gt.kl) then
                ijkl=ioff(ij)+kl
              else
                ijkl=ioff(kl)+ij
              end if
              i2e(ijkl)=eplll(i,j,k,l)
            end if
          end do
        end do
      end do
    end do
  end do
end do
write (9,'(A,I8,A,I8)') "Number of screened (pl|ll) integral sets:",nschw," out of:",npl*nll

! Then (dd|dd)
nschw=0
do ii=2,ndd
  sh1=ddprs(1,ii)
  ao1=aosh(sh1)
  sh2=ddprs(2,ii)
  ao2=aosh(sh2)
  ab1=ind2(sh1,sh2)
  do jj=1,ii-1
    sh3=ddprs(1,jj)
    ao3=aosh(sh3)
    sh4=ddprs(2,jj)
    ao4=aosh(sh4)
    ab2=ind2(sh3,sh4)
    if ((schw(ab1)*schw(ab2)).lt.schwthresh) then
      nschw=nschw+1
      cycle
    end if
    call i2e_dddd(sh1,sh2,sh3,sh4,edddd)
    do l=1,5
      do k=1,5
        kl=ind2(ao3+k-1,ao4+l-1)
        do j=1,5
          do i=1,5
            if (abs(edddd(i,j,k,l)).gt.intthresh) then
              ij=ind2(ao1+i-1,ao2+j-1)
              if (ij.gt.kl) then
                ijkl=ioff(ij)+kl
              else
                ijkl=ioff(kl)+ij
              end if
              i2e(ijkl)=edddd(i,j,k,l)
            end if
          end do
        end do
      end do
    end do
  end do
end do
write (9,'(A,I8,A,I8)') "Number of screened (dd|dd) integral sets:",nschw," out of:",ndd*(ndd+1)/2
!call time_checker(-1,"Two electron integrals - (dd|dd) terms, timing:")

! Then (dd|dl)
nschw=0
do ii=1,ndd
  sh1=ddprs(1,ii)
  ao1=aosh(sh1)
  sh2=ddprs(2,ii)
  ao2=aosh(sh2)
  ab1=ind2(sh1,sh2)
  do jj=1,ndl
    sh3=dlprs(1,jj)
    ao3=aosh(sh3)
    sh4=dlprs(2,jj)
    ao4=aosh(sh4)
    ab2=ind2(sh3,sh4)
    if ((schw(ab1)*schw(ab2)).lt.schwthresh) then
      nschw=nschw+1
      cycle
    end if
    call i2e_dddl(sh1,sh2,sh3,sh4,edddl)
    do l=1,4
      do k=1,5
        kl=ind2(ao3+k-1,ao4+l-1)
        do j=1,5
          do i=1,5
            if (abs(edddl(i,j,k,l)).gt.intthresh) then
              ij=ind2(ao1+i-1,ao2+j-1)
              if (ij.gt.kl) then
                ijkl=ioff(ij)+kl
              else
                ijkl=ioff(kl)+ij
              end if
              i2e(ijkl)=edddl(i,j,k,l)
            end if
          end do
        end do
      end do
    end do
  end do
end do
write (9,'(A,I8,A,I8)') "Number of screened (dd|dl) integral sets:",nschw," out of:",ndd*ndl
!call time_checker(-1,"Two electron integrals - (dd|dl) terms, timing:")

! Then (dd|ll)
nschw=0
do ii=1,ndd
  sh1=ddprs(1,ii)
  ao1=aosh(sh1)
  sh2=ddprs(2,ii)
  ao2=aosh(sh2)
  ab1=ind2(sh1,sh2)
  do jj=1,nll
    sh3=llprs(1,jj)
    ao3=aosh(sh3)
    sh4=llprs(2,jj)
    ao4=aosh(sh4)
    ab2=ind2(sh3,sh4)
    if ((schw(ab1)*schw(ab2)).lt.schwthresh) then
      nschw=nschw+1
      cycle
    end if
    call i2e_ddll(sh1,sh2,sh3,sh4,eddll)
    do l=1,4
      do k=1,4
        kl=ind2(ao3+k-1,ao4+l-1)
        do j=1,5
          do i=1,5
            if (abs(eddll(i,j,k,l)).gt.intthresh) then
              ij=ind2(ao1+i-1,ao2+j-1)
              if (ij.gt.kl) then
                ijkl=ioff(ij)+kl
              else
                ijkl=ioff(kl)+ij
              end if
              i2e(ijkl)=eddll(i,j,k,l)
            end if
          end do
        end do
      end do
    end do
  end do
end do
write (9,'(A,I8,A,I8)') "Number of screened (dd|ll) integral sets:",nschw," out of:",ndd*nll
!call time_checker(-1,"Two electron integrals - (dd|ll) terms, timing:")

! Then (dl|dl)
nschw=0
do ii=2,ndl
  sh1=dlprs(1,ii)
  ao1=aosh(sh1)
  sh2=dlprs(2,ii)
  ao2=aosh(sh2)
  ab1=ind2(sh1,sh2)
  do jj=1,ii-1
    sh3=dlprs(1,jj)
    ao3=aosh(sh3)
    sh4=dlprs(2,jj)
    ao4=aosh(sh4)
    ab2=ind2(sh3,sh4)
    if ((schw(ab1)*schw(ab2)).lt.schwthresh) then
      nschw=nschw+1
      cycle
    end if
    call i2e_dldl(sh1,sh2,sh3,sh4,edldl)
    do l=1,4
      do k=1,5
        kl=ind2(ao3+k-1,ao4+l-1)
        do j=1,4
          do i=1,5
            if (abs(edldl(i,j,k,l)).gt.intthresh) then
              ij=ind2(ao1+i-1,ao2+j-1)
              if (ij.gt.kl) then
                ijkl=ioff(ij)+kl
              else
                ijkl=ioff(kl)+ij
              end if
              i2e(ijkl)=edldl(i,j,k,l)
            end if
          end do
        end do
      end do
    end do
  end do
end do
write (9,'(A,I8,A,I8)') "Number of screened (dl|dl) integral sets:",nschw," out of:",ndl*(ndl+1)/2
!call time_checker(-1,"Two electron integrals - (dl|dl) terms, timing:")

! Then (dl|ll)
nschw=0
do ii=1,ndl
  sh1=dlprs(1,ii)
  ao1=aosh(sh1)
  sh2=dlprs(2,ii)
  ao2=aosh(sh2)
  ab1=ind2(sh1,sh2)
  do jj=1,nll
    sh3=llprs(1,jj)
    ao3=aosh(sh3)
    sh4=llprs(2,jj)
    ao4=aosh(sh4)
    ab2=ind2(sh3,sh4)
    if ((schw(ab1)*schw(ab2)).lt.schwthresh) then
      nschw=nschw+1
      cycle
    end if
    call i2e_dlll(sh1,sh2,sh3,sh4,edlll)
    do l=1,4
      do k=1,4
        kl=ind2(ao3+k-1,ao4+l-1)
        do j=1,4
          do i=1,5
            if (abs(edlll(i,j,k,l)).gt.intthresh) then
              ij=ind2(ao1+i-1,ao2+j-1)
              if (ij.gt.kl) then
                ijkl=ioff(ij)+kl
              else
                ijkl=ioff(kl)+ij
              end if
              i2e(ijkl)=edlll(i,j,k,l)
            end if
          end do
        end do
      end do
    end do
  end do
end do
write (9,'(A,I8,A,I8)') "Number of screened (dl|ll) integral sets:",nschw," out of:",ndl*nll
!call time_checker(-1,"Two electron integrals - (dl|ll) terms, timing:")

! Then (ll|ll)
nschw=0
do ii=2,nll
  sh1=llprs(1,ii)
  ao1=aosh(sh1)
  sh2=llprs(2,ii)
  ao2=aosh(sh2)
  ab1=ind2(sh1,sh2)
  do jj=1,ii-1
    sh3=llprs(1,jj)
    ao3=aosh(sh3)
    sh4=llprs(2,jj)
    ao4=aosh(sh4)
    ab2=ind2(sh3,sh4)
    if ((schw(ab1)*schw(ab2)).lt.schwthresh) then
      nschw=nschw+1
      cycle
    end if
    call i2e_llll(sh1,sh2,sh3,sh4,ellll)
    do l=1,4
      do k=1,4
        kl=ind2(ao3+k-1,ao4+l-1)
        do j=1,4
          do i=1,4
            if (abs(ellll(i,j,k,l)).gt.intthresh) then
              ij=ind2(ao1+i-1,ao2+j-1)
              if (ij.gt.kl) then
                ijkl=ioff(ij)+kl
              else
                ijkl=ioff(kl)+ij
              end if
              i2e(ijkl)=ellll(i,j,k,l)
            end if
          end do
        end do
      end do
    end do
  end do
end do
write (9,'(A,I8,A,I8)') "Number of screened (ll|ll) integral sets:",nschw," out of:",nll*(nll+1)/2
!call time_checker(-1,"Two electron integrals - (ll|ll) terms, timing:")


!kl=1; ij=1
!do i=1,nb4
!   ia=revind2(1,ij) ; ja=revind2(2,ij)
!   ka=revind2(1,kl) ; la=revind2(2,kl)
!   kl=kl+1
!   if (kl.gt.ij) then
!     kl=1
!     ij=ij+1
!   end if
!  if (abs(i2e(i)).lt.intthresh) then
!    write(982,'(i8,4i4,F14.9)') i,ia,ja,ka,la,abs(i2e(i))
!  else
!    write(982,'(i8,4i4,F14.9)') i,ia,ja,ka,la,i2e(i)
!  end if
!end do

end subroutine cmpt2e

