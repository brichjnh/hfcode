subroutine fock_build_non_diag(f,p,thr,nprsscr)
use nrtype; use molprops ; use omp_lib ; use s_sp_l_terms ; use s_d_l_terms ; use s_p_d_terms
implicit none

! This is the routine called to build the Fock matrix in cycles other than the first
! So the Schwarz integrals are available.
! This only does the OFF-diagonal shell pairs

real(dp), intent(in) :: p(nb,nb), thr
real(dp), intent(inout) :: f(nb,nb)
integer(i4b), intent(out) :: nprsscr

integer(i4b) :: thrn, nthr
integer(i4b) :: i, j,ij,k,l,ii,jj, ia, ja, ka, la, ab1, ab2, sh1, sh2, sh3, sh4
integer(i4b), save :: first_iter
logical(lgt) :: sameabab, allsame, sameij, samejk, samekl, bla
real(dp) :: essss, esssp(3), esspp(3,3), espsp(3,3), esppp(3,3,3), epppp(3,3,3,3)
real(dp) :: esssd(5), essdd(5,5), esdsd(5,5), esddd(5,5,5), edddd(5,5,5,5), t1, t2
real(dp) :: esspd(3,5),espsd(3,5),esppd(3,3,5),espdd(3,5,5),esdpp(5,3,3),esdpd(5,3,5),epppd(3,3,3,5)
real(dp) :: esssl(4), esspl(3,4), essll(4,4), espsl(3,4), esppl(3,3,4), espll(3,4,4)
real(dp) :: eslsl(4,4), eslpp(4,3,3), eslpl(4,3,4), eslll(4,4,4), epppl(3,3,3,4)
real(dp) :: esldl(4,5,4), esldd(4,5,5), esdll(5,4,4), esddl(5,5,4), esdsl(5,4), essdl(5,4)
real(dp) :: eppll(3,3,4,4), eplpl(3,4,3,4), eplll(3,4,4,4), ellll(4,4,4,4), edlll(5,4,4,4)
real(dp) :: eppdd(3,3,5,5), epdpd(3,5,3,5), epddd(3,5,5,5), eddll(5,5,4,4), edldl(5,4,5,4), edddl(5,5,5,4)
real(dp) :: maxdm(nab), fac, maxdmcoul, maxdmexch, maxdmval
real(dp) :: floc(nb,nb)

call evaluate_maxdm(p,maxdm)

nprsscr=0

! Start with (ss|ss)

!$OMP PARALLEL default(private) shared(maxdm,f,p,ssprs,spprs,first_iter, &
!$OMP      sdprs,slprs,ppprs,pdprs,plprs,ddprs,dlprs,llprs, &
!$OMP      schwss,schwsp,schwsd,schwsl,schwpp,schwpd,schwpl,schwdd,schwdl,schwll, &
!$OMP      orderss,ordersp,ordersd,ordersl,orderpp,orderpd,orderpl,orderdd,orderdl,orderll, &
!$OMP      ind2,thr,nss,nsp,nsd,nsl,npp,npd,npl,ndd,ndl,nll,aosh,atsh)
nthr = OMP_GET_NUM_THREADS()
thrn = OMP_GET_THREAD_NUM()
if (first_iter.ne.20) then
  !$OMP MASTER
  write (*,*) "Using this number of threads:",nthr
  write (9,*) "Using this number of threads:",nthr
  !$OMP END MASTER
  first_iter=20
end if
floc=0._dp
call time_diff_omp(0," ")
!$OMP DO SCHEDULE(DYNAMIC)
do ii=2,nss
  sh1=ssprs(1,orderss(ii))
  sh2=ssprs(2,orderss(ii))
  ia=aosh(sh1) ; ja=aosh(sh2)
  sameij=(sh1.eq.sh2)
  do jj=1,ii-1
    sh3=ssprs(1,orderss(jj))
    sh4=ssprs(2,orderss(jj))
    ka=aosh(sh3) ; la=aosh(sh4)
    samekl=(sh3.eq.sh4)
    samejk=(sh2.eq.sh3)
    maxdmcoul=max(maxdm(ind2(sh1,sh2)),maxdm(ind2(sh3,sh4)))
    maxdmexch=max(maxdm(ind2(sh1,sh3)),maxdm(ind2(sh1,sh4)),maxdm(ind2(sh2,sh3)),&
       & maxdm(ind2(sh2,sh4)))
    maxdmval=max(maxdmcoul,.25_dp*maxdmexch)
    if (schwss(orderss(ii))*schwss(orderss(jj))*maxdmval.lt.thr) then
        nprsscr=nprsscr+1
        cycle
    end if
    call i2e_ssss(sh1,sh2,sh3,sh4,essss)
    if (abs(essss).lt.thr) cycle
    if (sameij.and.samekl) then   ! (aa|bb) - note (Aa|aa) is impossible because this is off-diagonal
      floc(ia,ia)=floc(ia,ia)+p(ka,ka)*essss
      fac=.5_dp*p(ia,ka)*essss
      floc(ia,ka)=floc(ia,ka)-fac
      floc(ka,ia)=floc(ka,ia)-fac
      floc(ka,ka)=floc(ka,ka)+p(ia,ia)*essss
    else if (sameij) then         ! (aa|bc) or (aa|ab)
      floc(ia,ia)=floc(ia,ia)+2._dp*p(ka,la)*essss
      fac=.5_dp*p(ia,ka)*essss
      floc(ia,la)=floc(ia,la)-fac
      floc(la,ia)=floc(la,ia)-fac
      fac=.5_dp*p(ia,la)*essss
      floc(ia,ka)=floc(ia,ka)-fac
      floc(ka,ia)=floc(ka,ia)-fac
      fac=p(ia,ia)*essss
      floc(ka,la)=floc(ka,la)+fac
      floc(la,ka)=floc(la,ka)+fac
    else if (samekl) then   ! (ab|cc) or (ab|aa)
      floc(ka,ka)=floc(ka,ka)+2._dp*p(ia,ja)*essss
      fac=.5_dp*p(ka,ia)*essss
      floc(ka,ja)=floc(ka,ja)-fac
      floc(ja,ka)=floc(ja,ka)-fac
      fac=.5_dp*p(ka,ja)*essss
      floc(ka,ia)=floc(ka,ia)-fac
      floc(ia,ka)=floc(ia,ka)-fac
      fac=p(ka,ka)*essss
      floc(ia,ja)=floc(ia,ja)+fac
      floc(ja,ia)=floc(ja,ia)+fac
    else       !  (ab|cd) or (ab|ac)
      fac=2._dp*  p(la,ka) * essss
      floc(ia,ja) = floc(ia,ja) +  fac
      floc(ja,ia) = floc(ja,ia) +  fac
      fac=.5_dp*p(ja,ka)*essss
      floc(ia,la) = floc(ia,la) - fac
      floc(la,ia) = floc(la,ia) - fac
      fac=.5_dp*p(ja,la)*essss
      floc(ia,ka) = floc(ia,ka) - fac
      floc(ka,ia) = floc(ka,ia) - fac
      fac=.5_dp*p(ia,ka)*essss
      floc(ja,la) = floc(ja,la) - fac
      floc(la,ja) = floc(la,ja) - fac
      fac=.5_dp*p(ia,la)*essss
      floc(ja,ka) = floc(ja,ka) - fac
      floc(ka,ja) = floc(ka,ja) - fac
      fac=2._dp*p(ja,ia)*essss
      floc(ka,la) = floc(ka,la) +  fac
      floc(la,ka) = floc(la,ka) +   fac
    end if
  end do
end do
!$OMP END DO
!$OMP MASTER
if (prtlevl.ge.4) call time_diff_omp(1,"(ss|ss) integ:")
!$OMP END MASTER

!$OMP DO SCHEDULE(DYNAMIC)
! (ss|sp)
do ii=1,nss
  sh1=ssprs(1,orderss(ii))
  sh2=ssprs(2,orderss(ii))
  ab1=ind2(sh1,sh2)
  ia=aosh(sh1) ; ja=aosh(sh2)
  sameij=(ia.eq.ja)
  do jj=1,nsp
    sh3=spprs(1,ordersp(jj))
    sh4=spprs(2,ordersp(jj))
    ab2=ind2(sh3,sh4)
    ka=aosh(sh3) ; la=aosh(sh4)
    maxdmcoul=max(maxdm(ind2(sh1,sh2)),maxdm(ind2(sh3,sh4)))
    maxdmexch=max(maxdm(ind2(sh1,sh3)),maxdm(ind2(sh1,sh4)),maxdm(ind2(sh2,sh3)),&
       & maxdm(ind2(sh2,sh4)))
    maxdmval=max(maxdmcoul,.25_dp*maxdmexch)
    if (schwss(orderss(ii))*schwsp(ordersp(jj))*maxdmval.lt.thr) then
        nprsscr=nprsscr+1
        cycle
    end if
    call i2e_sssp(sh1,sh2,sh3,sh4,esssp)
    if (maxval(abs(esssp)).lt.thr) cycle
    if (sameij) then   ! (aa|bc) or (aa|ab)
      do l=1,3
        la=aosh(sh4)+l-1
        floc(ia,ia)=floc(ia,ia)+2._dp*p(ka,la)*esssp(l)
        fac=.5_dp*p(ia,ka)*esssp(l)
        floc(ia,la)=floc(ia,la)-fac
        floc(la,ia)=floc(la,ia)-fac
        fac=.5_dp*p(ia,la)*esssp(l)
        floc(ia,ka)=floc(ia,ka)-fac
        floc(ka,ia)=floc(ka,ia)-fac
        fac=p(ia,ia)*esssp(l)
        floc(ka,la)=floc(ka,la)+fac
        floc(la,ka)=floc(la,ka)+fac
      end do
    else       !  (ab|cd) or (ab|ac)
      do l=1,3
        la=aosh(sh4)+l-1
        fac=2._dp*  p(la,ka) * esssp(l)
        floc(ia,ja) = floc(ia,ja) +  fac
        floc(ja,ia) = floc(ja,ia) +  fac
        fac=.5_dp*p(ja,ka)*esssp(l)
        floc(ia,la) = floc(ia,la) - fac
        floc(la,ia) = floc(la,ia) - fac
        fac=.5_dp*p(ja,la)*esssp(l)
        floc(ia,ka) = floc(ia,ka) - fac
        floc(ka,ia) = floc(ka,ia) - fac
        fac=.5_dp*p(ia,ka)*esssp(l)
        floc(ja,la) = floc(ja,la) - fac
        floc(la,ja) = floc(la,ja) - fac
        fac=.5_dp*p(ia,la)*esssp(l)
        floc(ja,ka) = floc(ja,ka) - fac
        floc(ka,ja) = floc(ka,ja) - fac
        fac=2._dp*p(ja,ia)*esssp(l)
        floc(ka,la) = floc(ka,la) +  fac
        floc(la,ka) = floc(la,ka) +   fac
      end do
    end if
  end do
end do
!$OMP END DO
!$OMP MASTER
if (prtlevl.ge.4) call time_diff_omp(1,"(ss|sp) integ:")
!$OMP END MASTER

!$OMP DO SCHEDULE(DYNAMIC)
do ii=1,nss ! (ss|sd)
  sh1=ssprs(1,orderss(ii))
  sh2=ssprs(2,orderss(ii))
  ab1=ind2(sh1,sh2)
  ia=aosh(sh1) ; ja=aosh(sh2)
  sameij=(ia.eq.ja)
  do jj=1,nsd
    sh3=sdprs(1,ordersd(jj))
    sh4=sdprs(2,ordersd(jj))
    ab2=ind2(sh3,sh4)
    maxdmcoul=max(maxdm(ind2(sh1,sh2)),maxdm(ind2(sh3,sh4)))
    maxdmexch=max(maxdm(ind2(sh1,sh3)),maxdm(ind2(sh1,sh4)),maxdm(ind2(sh2,sh3)),&
       & maxdm(ind2(sh2,sh4)))
    maxdmval=max(maxdmcoul,.25_dp*maxdmexch)
    if (schwss(orderss(ii))*schwsd(ordersd(jj))*maxdmval.lt.thr) then
        nprsscr=nprsscr+1
        cycle
    end if
    ka=aosh(sh3) ; la=aosh(sh4)
    call i2e_sssd(sh1,sh2,sh3,sh4,esssd)
    if (maxval(abs(esssd)).lt.thr) cycle
    if (sameij) then   ! (aa|bc) or (aa|ab)
      do l=1,5
        la=aosh(sh4)+l-1
        floc(ia,ia)=floc(ia,ia)+2._dp*p(ka,la)*esssd(l)
        fac=.5_dp*p(ia,ka)*esssd(l)
        floc(ia,la)=floc(ia,la)-fac
        floc(la,ia)=floc(la,ia)-fac
        fac=.5_dp*p(ia,la)*esssd(l)
        floc(ia,ka)=floc(ia,ka)-fac
        floc(ka,ia)=floc(ka,ia)-fac
        fac=p(ia,ia)*esssd(l)
        floc(ka,la)=floc(ka,la)+fac
        floc(la,ka)=floc(la,ka)+fac
      end do
    else       !  (ab|cd) or (ab|ac)
      do l=1,5
        la=aosh(sh4)+l-1
        fac=2._dp*  p(la,ka) * esssd(l)
        floc(ia,ja) = floc(ia,ja) +  fac
        floc(ja,ia) = floc(ja,ia) +  fac
        fac=.5_dp*p(ja,ka)*esssd(l)
        floc(ia,la) = floc(ia,la) - fac
        floc(la,ia) = floc(la,ia) - fac
        fac=.5_dp*p(ja,la)*esssd(l)
        floc(ia,ka) = floc(ia,ka) - fac
        floc(ka,ia) = floc(ka,ia) - fac
        fac=.5_dp*p(ia,ka)*esssd(l)
        floc(ja,la) = floc(ja,la) - fac
        floc(la,ja) = floc(la,ja) - fac
        fac=.5_dp*p(ia,la)*esssd(l)
        floc(ja,ka) = floc(ja,ka) - fac
        floc(ka,ja) = floc(ka,ja) - fac
        fac=2._dp*p(ja,ia)*esssd(l)
        floc(ka,la) = floc(ka,la) +  fac
        floc(la,ka) = floc(la,ka) +   fac
      end do
    end if
  end do
end do
!$OMP END DO
!$OMP MASTER
if (prtlevl.ge.4) call time_diff_omp(1,"(ss|sd) integ:")
!$OMP END MASTER

!$OMP DO SCHEDULE(DYNAMIC)

! (ss|sl)

do ii=1,nss
  sh1=ssprs(1,orderss(ii))
  sh2=ssprs(2,orderss(ii))
  ab1=ind2(sh1,sh2)
  ia=aosh(sh1) ; ja=aosh(sh2)
  sameij=(ia.eq.ja)
  do jj=1,nsl
    sh3=slprs(1,ordersl(jj))
    sh4=slprs(2,ordersl(jj))
    ab2=ind2(sh3,sh4)
    ka=aosh(sh3) ; la=aosh(sh4)
    maxdmcoul=max(maxdm(ind2(sh1,sh2)),maxdm(ind2(sh3,sh4)))
    maxdmexch=max(maxdm(ind2(sh1,sh3)),maxdm(ind2(sh1,sh4)),maxdm(ind2(sh2,sh3)),&
       & maxdm(ind2(sh2,sh4)))
    maxdmval=max(maxdmcoul,.25_dp*maxdmexch)
    if (schwss(orderss(ii))*schwsl(ordersl(jj))*maxdmval.lt.thr) then
        nprsscr=nprsscr+1
        cycle
    end if
    call i2e_sssl(sh1,sh2,sh3,sh4,esssl)
    if (maxval(abs(esssl)).lt.thr) cycle
    if (sameij) then   ! (aa|bc) or (aa|ab)
      do l=1,4
        la=aosh(sh4)+l-1
        floc(ia,ia)=floc(ia,ia)+2._dp*p(ka,la)*esssl(l)
        fac=.5_dp*p(ia,ka)*esssl(l)
        floc(ia,la)=floc(ia,la)-fac
        floc(la,ia)=floc(la,ia)-fac
        fac=.5_dp*p(ia,la)*esssl(l)
        floc(ia,ka)=floc(ia,ka)-fac
        floc(ka,ia)=floc(ka,ia)-fac
        fac=p(ia,ia)*esssl(l)
        floc(ka,la)=floc(ka,la)+fac
        floc(la,ka)=floc(la,ka)+fac
      end do
    else       !  (ab|cd) or (ab|ac)
      do l=1,4
        la=aosh(sh4)+l-1
        fac=2._dp*  p(la,ka) * esssl(l)
        floc(ia,ja) = floc(ia,ja) +  fac
        floc(ja,ia) = floc(ja,ia) +  fac
        fac=.5_dp*p(ja,ka)*esssl(l)
        floc(ia,la) = floc(ia,la) - fac
        floc(la,ia) = floc(la,ia) - fac
        fac=.5_dp*p(ja,la)*esssl(l)
        floc(ia,ka) = floc(ia,ka) - fac
        floc(ka,ia) = floc(ka,ia) - fac
        fac=.5_dp*p(ia,ka)*esssl(l)
        floc(ja,la) = floc(ja,la) - fac
        floc(la,ja) = floc(la,ja) - fac
        fac=.5_dp*p(ia,la)*esssl(l)
        floc(ja,ka) = floc(ja,ka) - fac
        floc(ka,ja) = floc(ka,ja) - fac
        fac=2._dp*p(ja,ia)*esssl(l)
        floc(ka,la) = floc(ka,la) +  fac
        floc(la,ka) = floc(la,ka) +   fac
      end do
    end if
  end do
end do
!$OMP END DO
!$OMP MASTER
if (prtlevl.ge.4) call time_diff_omp(1,"(ss|sl) integ:")
!$OMP END MASTER

!$OMP DO SCHEDULE(DYNAMIC)

do ii=1,nss  ! (ss|pp)
  sh1=ssprs(1,orderss(ii))
  sh2=ssprs(2,orderss(ii))
  ab1=ind2(sh1,sh2)
  ia=aosh(sh1) ; ja=aosh(sh2)
  sameij=(ia.eq.ja)
  do jj=1,npp
    sh3=ppprs(1,orderpp(jj))
    sh4=ppprs(2,orderpp(jj))
    ab2=ind2(sh3,sh4)
    ka=aosh(sh3) ; la=aosh(sh4)
    maxdmcoul=max(maxdm(ind2(sh1,sh2)),maxdm(ind2(sh3,sh4)))
    maxdmexch=max(maxdm(ind2(sh1,sh3)),maxdm(ind2(sh1,sh4)),maxdm(ind2(sh2,sh3)),&
       & maxdm(ind2(sh2,sh4)))
    maxdmval=max(maxdmcoul,.25_dp*maxdmexch)
    if (schwss(orderss(ii))*schwpp(orderpp(jj))*maxdmval.lt.thr) then
        nprsscr=nprsscr+1
        cycle
    end if
    call i2e_sspp(sh1,sh2,sh3,sh4,esspp)
    if (maxval(abs(esspp)).lt.thr) cycle
    samekl=(sh3.eq.sh4)
    do k=1,3
      ka=aosh(sh3)+k-1
      do l=1,3
        la=aosh(sh4)+l-1
        if (sameij.and.samekl.and.(k.eq.l)) then   ! (aa|bb)
          floc(ia,ia)=floc(ia,ia)+p(ka,ka)*esspp(k,k)
          fac=.5_dp*p(ia,ka)*esspp(k,k)
          floc(ia,ka)=floc(ia,ka)-fac
          floc(ka,ia)=floc(ka,ia)-fac
          floc(ka,ka)=floc(ka,ka)+p(ia,ia)*esspp(k,k)
        else if (sameij.and.samekl.and.(k.lt.l)) then   ! These are same as the below - only include one of the two
          cycle
        else if (sameij.and.samekl.and.(k.gt.l)) then   ! (aa|bc)
          floc(ia,ia)=floc(ia,ia)+2._dp*p(ka,la)*esspp(k,l)
          fac=.5_dp*p(ia,ka)*esspp(k,l)
          floc(ia,la)=floc(ia,la)-fac
          floc(la,ia)=floc(la,ia)-fac
          fac=.5_dp*p(ia,la)*esspp(k,l)
          floc(ia,ka)=floc(ia,ka)-fac
          floc(ka,ia)=floc(ka,ia)-fac
          fac=p(ia,ia)*esspp(k,l)
          floc(ka,la)=floc(ka,la)+fac
          floc(la,ka)=floc(la,ka)+fac
        else if (sameij) then   ! (aa|bc) so actually do the same as above...
          floc(ia,ia)=floc(ia,ia)+2._dp*p(ka,la)*esspp(k,l)
          fac=.5_dp*p(ia,ka)*esspp(k,l)
          floc(ia,la)=floc(ia,la)-fac
          floc(la,ia)=floc(la,ia)-fac
          fac=.5_dp*p(ia,la)*esspp(k,l)
          floc(ia,ka)=floc(ia,ka)-fac
          floc(ka,ia)=floc(ka,ia)-fac
          fac=p(ia,ia)*esspp(k,l)
          floc(ka,la)=floc(ka,la)+fac
          floc(la,ka)=floc(la,ka)+fac
        else if (samekl.and.(k.eq.l)) then   ! (ab|cc)
          floc(ka,ka)=floc(ka,ka)+2._dp*p(ia,ja)*esspp(k,k)
          fac=.5_dp*p(ka,ia)*esspp(k,k)
          floc(ka,ja)=floc(ka,ja)-fac
          floc(ja,ka)=floc(ja,ka)-fac
          fac=.5_dp*p(ka,ja)*esspp(k,k)
          floc(ka,ia)=floc(ka,ia)-fac
          floc(ia,ka)=floc(ia,ka)-fac
          fac=p(ka,ka)*esspp(k,k)
          floc(ia,ja)=floc(ia,ja)+fac
          floc(ja,ia)=floc(ja,ia)+fac
        else if (samekl.and.(k.lt.l)) then   ! (ab|cd) but same as below case
          cycle
        else    ! This includes the k.eq.l case, with k>l, and the i.ne.j, k.ne.l case
          fac=2._dp*  p(la,ka) * esspp(k,l)
          floc(ia,ja) = floc(ia,ja) +  fac
          floc(ja,ia) = floc(ja,ia) +  fac
          fac=.5_dp*p(ja,ka)*esspp(k,l)
          floc(ia,la) = floc(ia,la) - fac
          floc(la,ia) = floc(la,ia) - fac
          fac=.5_dp*p(ja,la)*esspp(k,l)
          floc(ia,ka) = floc(ia,ka) - fac
          floc(ka,ia) = floc(ka,ia) - fac
          fac=.5_dp*p(ia,ka)*esspp(k,l)
          floc(ja,la) = floc(ja,la) - fac
          floc(la,ja) = floc(la,ja) - fac
          fac=.5_dp*p(ia,la)*esspp(k,l)
          floc(ja,ka) = floc(ja,ka) - fac
          floc(ka,ja) = floc(ka,ja) - fac
          fac=2._dp*p(ja,ia)*esspp(k,l)
          floc(ka,la) = floc(ka,la) +  fac
          floc(la,ka) = floc(la,ka) +   fac
        end if
      end do
    end do
  end do
end do
!$OMP END DO
!$OMP MASTER
if (prtlevl.ge.4) call time_diff_omp(1,"(ss|pp) integ:")
!$OMP END MASTER

!$OMP DO SCHEDULE(DYNAMIC)
 
do ii=1,nss    ! (ss|pd)
  sh1=ssprs(1,orderss(ii))
  sh2=ssprs(2,orderss(ii))
  ab1=ind2(sh1,sh2)
  ia=aosh(sh1) ; ja=aosh(sh2)
  sameij=(ia.eq.ja)
  do jj=1,npd
    sh3=pdprs(1,orderpd(jj))
    sh4=pdprs(2,orderpd(jj))
    ab2=ind2(sh3,sh4)
    maxdmcoul=max(maxdm(ind2(sh1,sh2)),maxdm(ind2(sh3,sh4)))
    maxdmexch=max(maxdm(ind2(sh1,sh3)),maxdm(ind2(sh1,sh4)),maxdm(ind2(sh2,sh3)),&
       & maxdm(ind2(sh2,sh4)))
    maxdmval=max(maxdmcoul,.25_dp*maxdmexch)
    if (schwss(orderss(ii))*schwpd(orderpd(jj))*maxdmval.lt.thr) then
        nprsscr=nprsscr+1
        cycle
    end if
    ka=aosh(sh3) ; la=aosh(sh4)
    call i2e_sspd(sh1,sh2,sh3,sh4,esspd)
    if (maxval(abs(esspd)).lt.thr) cycle
    if (sameij) then   ! (aa|bc) or (aa|ab)
      do k=1,3
        ka=aosh(sh3)+k-1
        do l=1,5
          la=aosh(sh4)+l-1
          floc(ia,ia)=floc(ia,ia)+2._dp*p(ka,la)*esspd(k,l)
          fac=.5_dp*p(ia,ka)*esspd(k,l)
          floc(ia,la)=floc(ia,la)-fac
          floc(la,ia)=floc(la,ia)-fac
          fac=.5_dp*p(ia,la)*esspd(k,l)
          floc(ia,ka)=floc(ia,ka)-fac
          floc(ka,ia)=floc(ka,ia)-fac
          fac=p(ia,ia)*esspd(k,l)
          floc(ka,la)=floc(ka,la)+fac
          floc(la,ka)=floc(la,ka)+fac
        end do
      end do
    else       !  (ab|cd) or (ab|ac)
      do k=1,3
        ka=aosh(sh3)+k-1
        do l=1,5
          la=aosh(sh4)+l-1
          fac=2._dp*  p(la,ka) * esspd(k,l)
          floc(ia,ja) = floc(ia,ja) +  fac
          floc(ja,ia) = floc(ja,ia) +  fac
          fac=.5_dp*p(ja,ka)*esspd(k,l)
          floc(ia,la) = floc(ia,la) - fac
          floc(la,ia) = floc(la,ia) - fac
          fac=.5_dp*p(ja,la)*esspd(k,l)
          floc(ia,ka) = floc(ia,ka) - fac
          floc(ka,ia) = floc(ka,ia) - fac
          fac=.5_dp*p(ia,ka)*esspd(k,l)
          floc(ja,la) = floc(ja,la) - fac
          floc(la,ja) = floc(la,ja) - fac
          fac=.5_dp*p(ia,la)*esspd(k,l)
          floc(ja,ka) = floc(ja,ka) - fac
          floc(ka,ja) = floc(ka,ja) - fac
          fac=2._dp*p(ja,ia)*esspd(k,l)
          floc(ka,la) = floc(ka,la) +  fac
          floc(la,ka) = floc(la,ka) +   fac
        end do
      end do
    end if
  end do
end do
!$OMP END DO
!$OMP MASTER
if (prtlevl.ge.4) call time_diff_omp(1,"(ss|pd) integ:")
!$OMP END MASTER

!$OMP DO SCHEDULE(DYNAMIC)

do ii=1,nss    ! (ss|pl)
  sh1=ssprs(1,orderss(ii))
  sh2=ssprs(2,orderss(ii))
  ab1=ind2(sh1,sh2)
  ia=aosh(sh1) ; ja=aosh(sh2)
  sameij=(ia.eq.ja)
  do jj=1,npl
    sh3=plprs(1,orderpl(jj))
    sh4=plprs(2,orderpl(jj))
    ab2=ind2(sh3,sh4)
    ka=aosh(sh3) ; la=aosh(sh4)
    maxdmcoul=max(maxdm(ind2(sh1,sh2)),maxdm(ind2(sh3,sh4)))
    maxdmexch=max(maxdm(ind2(sh1,sh3)),maxdm(ind2(sh1,sh4)),maxdm(ind2(sh2,sh3)),&
       & maxdm(ind2(sh2,sh4)))
    maxdmval=max(maxdmcoul,.25_dp*maxdmexch)
    if (schwss(orderss(ii))*schwpl(orderpl(jj))*maxdmval.lt.thr) then
        nprsscr=nprsscr+1
        cycle
    end if
    call i2e_sspl(sh1,sh2,sh3,sh4,esspl)
    if (maxval(abs(esspl)).lt.thr) cycle
    if (sameij) then   ! (aa|bc) or (aa|ab)
      do k=1,3
        ka=aosh(sh3)+k-1
        do l=1,4
          la=aosh(sh4)+l-1
          floc(ia,ia)=floc(ia,ia)+2._dp*p(ka,la)*esspl(k,l)
          fac=.5_dp*p(ia,ka)*esspl(k,l)
          floc(ia,la)=floc(ia,la)-fac
          floc(la,ia)=floc(la,ia)-fac
          fac=.5_dp*p(ia,la)*esspl(k,l)
          floc(ia,ka)=floc(ia,ka)-fac
          floc(ka,ia)=floc(ka,ia)-fac
          fac=p(ia,ia)*esspl(k,l)
          floc(ka,la)=floc(ka,la)+fac
          floc(la,ka)=floc(la,ka)+fac
        end do
      end do
    else       !  (ab|cd) or (ab|ac)
      do k=1,3
        ka=aosh(sh3)+k-1
        do l=1,4
          la=aosh(sh4)+l-1
          fac=2._dp*  p(la,ka) * esspl(k,l)
          floc(ia,ja) = floc(ia,ja) +  fac
          floc(ja,ia) = floc(ja,ia) +  fac
          fac=.5_dp*p(ja,ka)*esspl(k,l)
          floc(ia,la) = floc(ia,la) - fac
          floc(la,ia) = floc(la,ia) - fac
          fac=.5_dp*p(ja,la)*esspl(k,l)
          floc(ia,ka) = floc(ia,ka) - fac
          floc(ka,ia) = floc(ka,ia) - fac
          fac=.5_dp*p(ia,ka)*esspl(k,l)
          floc(ja,la) = floc(ja,la) - fac
          floc(la,ja) = floc(la,ja) - fac
          fac=.5_dp*p(ia,la)*esspl(k,l)
          floc(ja,ka) = floc(ja,ka) - fac
          floc(ka,ja) = floc(ka,ja) - fac
          fac=2._dp*p(ja,ia)*esspl(k,l)
          floc(ka,la) = floc(ka,la) +  fac
          floc(la,ka) = floc(la,ka) +   fac
        end do
      end do
    end if
  end do
end do
!$OMP END DO
!$OMP MASTER
if (prtlevl.ge.4) call time_diff_omp(1,"(ss|pl) integ:")
!$OMP END MASTER

!$OMP DO SCHEDULE(DYNAMIC)

do ii=1,nss  ! (ss|dd)
  sh1=ssprs(1,orderss(ii))
  sh2=ssprs(2,orderss(ii))
  ab1=ind2(sh1,sh2)
  ia=aosh(sh1) ; ja=aosh(sh2)
  sameij=(ia.eq.ja)
  do jj=1,ndd
    sh3=ddprs(1,orderdd(jj))
    sh4=ddprs(2,orderdd(jj))
    ab2=ind2(sh3,sh4)
    maxdmcoul=max(maxdm(ind2(sh1,sh2)),maxdm(ind2(sh3,sh4)))
    maxdmexch=max(maxdm(ind2(sh1,sh3)),maxdm(ind2(sh1,sh4)),maxdm(ind2(sh2,sh3)),&
       & maxdm(ind2(sh2,sh4)))
    maxdmval=max(maxdmcoul,.25_dp*maxdmexch)
    if (schwss(orderss(ii))*schwdd(orderdd(jj))*maxdmval.lt.thr) then
        nprsscr=nprsscr+1
        cycle
    end if
    ka=aosh(sh3) ; la=aosh(sh4)
    call i2e_ssdd(sh1,sh2,sh3,sh4,essdd)
    if (maxval(abs(essdd)).lt.thr) cycle
    samekl=(sh3.eq.sh4)
    do k=1,5
      ka=aosh(sh3)+k-1
      do l=1,5
        la=aosh(sh4)+l-1
        if (sameij.and.samekl.and.(k.eq.l)) then   ! (aa|bb)
          floc(ia,ia)=floc(ia,ia)+p(ka,ka)*essdd(k,k)
          fac=.5_dp*p(ia,ka)*essdd(k,k)
          floc(ia,ka)=floc(ia,ka)-fac
          floc(ka,ia)=floc(ka,ia)-fac
          floc(ka,ka)=floc(ka,ka)+p(ia,ia)*essdd(k,k)
        else if (sameij.and.samekl.and.(k.lt.l)) then   ! These are same as the below - only include one of the two
          cycle
        else if (sameij.and.samekl.and.(k.gt.l)) then   ! (aa|bc)
          floc(ia,ia)=floc(ia,ia)+2._dp*p(ka,la)*essdd(k,l)
          fac=.5_dp*p(ia,ka)*essdd(k,l)
          floc(ia,la)=floc(ia,la)-fac
          floc(la,ia)=floc(la,ia)-fac
          fac=.5_dp*p(ia,la)*essdd(k,l)
          floc(ia,ka)=floc(ia,ka)-fac
          floc(ka,ia)=floc(ka,ia)-fac
          fac=p(ia,ia)*essdd(k,l)
          floc(ka,la)=floc(ka,la)+fac
          floc(la,ka)=floc(la,ka)+fac
        else if (sameij) then   ! (aa|bc) so actually do the same as above...
          floc(ia,ia)=floc(ia,ia)+2._dp*p(ka,la)*essdd(k,l)
          fac=.5_dp*p(ia,ka)*essdd(k,l)
          floc(ia,la)=floc(ia,la)-fac
          floc(la,ia)=floc(la,ia)-fac
          fac=.5_dp*p(ia,la)*essdd(k,l)
          floc(ia,ka)=floc(ia,ka)-fac
          floc(ka,ia)=floc(ka,ia)-fac
          fac=p(ia,ia)*essdd(k,l)
          floc(ka,la)=floc(ka,la)+fac
          floc(la,ka)=floc(la,ka)+fac
        else if (samekl.and.(k.eq.l)) then   ! (ab|cc)
          floc(ka,ka)=floc(ka,ka)+2._dp*p(ia,ja)*essdd(k,k)
          fac=.5_dp*p(ka,ia)*essdd(k,k)
          floc(ka,ja)=floc(ka,ja)-fac
          floc(ja,ka)=floc(ja,ka)-fac
          fac=.5_dp*p(ka,ja)*essdd(k,k)
          floc(ka,ia)=floc(ka,ia)-fac
          floc(ia,ka)=floc(ia,ka)-fac
          fac=p(ka,ka)*essdd(k,k)
          floc(ia,ja)=floc(ia,ja)+fac
          floc(ja,ia)=floc(ja,ia)+fac
        else if (samekl.and.(k.lt.l)) then   ! (ab|cd) but same as below case
          cycle
        else    ! This includes the k.eq.l case, with k>l, and the i.ne.j, k.ne.l case
          fac=2._dp*  p(la,ka) * essdd(k,l)
          floc(ia,ja) = floc(ia,ja) +  fac
          floc(ja,ia) = floc(ja,ia) +  fac
          fac=.5_dp*p(ja,ka)*essdd(k,l)
          floc(ia,la) = floc(ia,la) - fac
          floc(la,ia) = floc(la,ia) - fac
          fac=.5_dp*p(ja,la)*essdd(k,l)
          floc(ia,ka) = floc(ia,ka) - fac
          floc(ka,ia) = floc(ka,ia) - fac
          fac=.5_dp*p(ia,ka)*essdd(k,l)
          floc(ja,la) = floc(ja,la) - fac
          floc(la,ja) = floc(la,ja) - fac
          fac=.5_dp*p(ia,la)*essdd(k,l)
          floc(ja,ka) = floc(ja,ka) - fac
          floc(ka,ja) = floc(ka,ja) - fac
          fac=2._dp*p(ja,ia)*essdd(k,l)
          floc(ka,la) = floc(ka,la) +  fac
          floc(la,ka) = floc(la,ka) +   fac
        end if
      end do
    end do
  end do
end do

!$OMP END DO
!$OMP MASTER
if (prtlevl.ge.4) call time_diff_omp(1,"(ss|dd) integ:")
!$OMP END MASTER

!$OMP DO SCHEDULE(DYNAMIC)

do ii=1,nss    ! (ss|dl)
  sh1=ssprs(1,orderss(ii))
  sh2=ssprs(2,orderss(ii))
  ab1=ind2(sh1,sh2)
  ia=aosh(sh1) ; ja=aosh(sh2)
  sameij=(ia.eq.ja)
  do jj=1,ndl
    sh3=dlprs(1,orderdl(jj))
    sh4=dlprs(2,orderdl(jj))
    ab2=ind2(sh3,sh4)
    ka=aosh(sh3) ; la=aosh(sh4)
    maxdmcoul=max(maxdm(ind2(sh1,sh2)),maxdm(ind2(sh3,sh4)))
    maxdmexch=max(maxdm(ind2(sh1,sh3)),maxdm(ind2(sh1,sh4)),maxdm(ind2(sh2,sh3)),&
       & maxdm(ind2(sh2,sh4)))
    maxdmval=max(maxdmcoul,.25_dp*maxdmexch)
    if (schwss(orderss(ii))*schwdl(orderdl(jj))*maxdmval.lt.thr) then
        nprsscr=nprsscr+1
        cycle
    end if
    call i2e_ssdl(sh1,sh2,sh3,sh4,essdl)
    if (maxval(abs(essdl)).lt.thr) cycle
    if (sameij) then   ! (aa|bc) or (aa|ab)
      do k=1,5
        ka=aosh(sh3)+k-1
        do l=1,4
          la=aosh(sh4)+l-1
          floc(ia,ia)=floc(ia,ia)+2._dp*p(ka,la)*essdl(k,l)
          fac=.5_dp*p(ia,ka)*essdl(k,l)
          floc(ia,la)=floc(ia,la)-fac
          floc(la,ia)=floc(la,ia)-fac
          fac=.5_dp*p(ia,la)*essdl(k,l)
          floc(ia,ka)=floc(ia,ka)-fac
          floc(ka,ia)=floc(ka,ia)-fac
          fac=p(ia,ia)*essdl(k,l)
          floc(ka,la)=floc(ka,la)+fac
          floc(la,ka)=floc(la,ka)+fac
        end do
      end do
    else       !  (ab|cd) or (ab|ac)
      do k=1,5
        ka=aosh(sh3)+k-1
        do l=1,4
          la=aosh(sh4)+l-1
          fac=2._dp*  p(la,ka) * essdl(k,l)
          floc(ia,ja) = floc(ia,ja) +  fac
          floc(ja,ia) = floc(ja,ia) +  fac
          fac=.5_dp*p(ja,ka)*essdl(k,l)
          floc(ia,la) = floc(ia,la) - fac
          floc(la,ia) = floc(la,ia) - fac
          fac=.5_dp*p(ja,la)*essdl(k,l)
          floc(ia,ka) = floc(ia,ka) - fac
          floc(ka,ia) = floc(ka,ia) - fac
          fac=.5_dp*p(ia,ka)*essdl(k,l)
          floc(ja,la) = floc(ja,la) - fac
          floc(la,ja) = floc(la,ja) - fac
          fac=.5_dp*p(ia,la)*essdl(k,l)
          floc(ja,ka) = floc(ja,ka) - fac
          floc(ka,ja) = floc(ka,ja) - fac
          fac=2._dp*p(ja,ia)*essdl(k,l)
          floc(ka,la) = floc(ka,la) +  fac
          floc(la,ka) = floc(la,ka) +   fac
        end do
      end do
    end if
  end do
end do

!$OMP END DO
!$OMP MASTER
if (prtlevl.ge.4) call time_diff_omp(1,"(ss|dl) integ:")
!$OMP END MASTER

!$OMP DO SCHEDULE(DYNAMIC)

do ii=1,nss  ! (ss|ll)
  sh1=ssprs(1,orderss(ii))
  sh2=ssprs(2,orderss(ii))
  ab1=ind2(sh1,sh2)
  ia=aosh(sh1) ; ja=aosh(sh2)
  sameij=(ia.eq.ja)
  do jj=1,nll
    sh3=llprs(1,orderll(jj))
    sh4=llprs(2,orderll(jj))
    ab2=ind2(sh3,sh4)
    ka=aosh(sh3) ; la=aosh(sh4)
    maxdmcoul=max(maxdm(ind2(sh1,sh2)),maxdm(ind2(sh3,sh4)))
    maxdmexch=max(maxdm(ind2(sh1,sh3)),maxdm(ind2(sh1,sh4)),maxdm(ind2(sh2,sh3)),&
       & maxdm(ind2(sh2,sh4)))
    maxdmval=max(maxdmcoul,.25_dp*maxdmexch)
    if (schwss(orderss(ii))*schwll(orderll(jj))*maxdmval.lt.thr) then
        nprsscr=nprsscr+1
        cycle
    end if
    call i2e_ssll(sh1,sh2,sh3,sh4,essll)
    if (maxval(abs(essll)).lt.thr) cycle
    samekl=(sh3.eq.sh4)
    do k=1,4
      ka=aosh(sh3)+k-1
      do l=1,4
        la=aosh(sh4)+l-1
        if (sameij.and.samekl.and.(k.eq.l)) then   ! (aa|bb)
          floc(ia,ia)=floc(ia,ia)+p(ka,ka)*essll(k,k)
          fac=.5_dp*p(ia,ka)*essll(k,k)
          floc(ia,ka)=floc(ia,ka)-fac
          floc(ka,ia)=floc(ka,ia)-fac
          floc(ka,ka)=floc(ka,ka)+p(ia,ia)*essll(k,k)
        else if (sameij.and.samekl.and.(k.lt.l)) then   ! These are same as the below - only include one of the two
          cycle
        else if (sameij.and.samekl.and.(k.gt.l)) then   ! (aa|bc)
          floc(ia,ia)=floc(ia,ia)+2._dp*p(ka,la)*essll(k,l)
          fac=.5_dp*p(ia,ka)*essll(k,l)
          floc(ia,la)=floc(ia,la)-fac
          floc(la,ia)=floc(la,ia)-fac
          fac=.5_dp*p(ia,la)*essll(k,l)
          floc(ia,ka)=floc(ia,ka)-fac
          floc(ka,ia)=floc(ka,ia)-fac
          fac=p(ia,ia)*essll(k,l)
          floc(ka,la)=floc(ka,la)+fac
          floc(la,ka)=floc(la,ka)+fac
        else if (sameij) then   ! (aa|bc) so actually do the same as above...
          floc(ia,ia)=floc(ia,ia)+2._dp*p(ka,la)*essll(k,l)
          fac=.5_dp*p(ia,ka)*essll(k,l)
          floc(ia,la)=floc(ia,la)-fac
          floc(la,ia)=floc(la,ia)-fac
          fac=.5_dp*p(ia,la)*essll(k,l)
          floc(ia,ka)=floc(ia,ka)-fac
          floc(ka,ia)=floc(ka,ia)-fac
          fac=p(ia,ia)*essll(k,l)
          floc(ka,la)=floc(ka,la)+fac
          floc(la,ka)=floc(la,ka)+fac
        else if (samekl.and.(k.eq.l)) then   ! (ab|cc)
          floc(ka,ka)=floc(ka,ka)+2._dp*p(ia,ja)*essll(k,k)
          fac=.5_dp*p(ka,ia)*essll(k,k)
          floc(ka,ja)=floc(ka,ja)-fac
          floc(ja,ka)=floc(ja,ka)-fac
          fac=.5_dp*p(ka,ja)*essll(k,k)
          floc(ka,ia)=floc(ka,ia)-fac
          floc(ia,ka)=floc(ia,ka)-fac
          fac=p(ka,ka)*essll(k,k)
          floc(ia,ja)=floc(ia,ja)+fac
          floc(ja,ia)=floc(ja,ia)+fac
        else if (samekl.and.(k.lt.l)) then   ! (ab|cd) but same as below case
          cycle
        else    ! This includes the k.eq.l case, with k>l, and the i.ne.j, k.ne.l case
          fac=2._dp*  p(la,ka) * essll(k,l)
          floc(ia,ja) = floc(ia,ja) +  fac
          floc(ja,ia) = floc(ja,ia) +  fac
          fac=.5_dp*p(ja,ka)*essll(k,l)
          floc(ia,la) = floc(ia,la) - fac
          floc(la,ia) = floc(la,ia) - fac
          fac=.5_dp*p(ja,la)*essll(k,l)
          floc(ia,ka) = floc(ia,ka) - fac
          floc(ka,ia) = floc(ka,ia) - fac
          fac=.5_dp*p(ia,ka)*essll(k,l)
          floc(ja,la) = floc(ja,la) - fac
          floc(la,ja) = floc(la,ja) - fac
          fac=.5_dp*p(ia,la)*essll(k,l)
          floc(ja,ka) = floc(ja,ka) - fac
          floc(ka,ja) = floc(ka,ja) - fac
          fac=2._dp*p(ja,ia)*essll(k,l)
          floc(ka,la) = floc(ka,la) +  fac
          floc(la,ka) = floc(la,ka) +   fac
        end if
      end do
    end do
  end do
end do

!$OMP END DO
!$OMP MASTER
if (prtlevl.ge.4) call time_diff_omp(1,"(ss|ll) integ:")
!$OMP END MASTER

!$OMP DO SCHEDULE(DYNAMIC)
do ii=2,nsp  !  (sp|sp)
  sh1=spprs(1,ordersp(ii))
  sh2=spprs(2,ordersp(ii))
  ab1=ind2(sh1,sh2)
  ia=aosh(sh1)
  do jj=1,ii-1
    sh3=spprs(1,ordersp(jj))
    sh4=spprs(2,ordersp(jj))
    ab2=ind2(sh3,sh4)
    ka=aosh(sh3)
    maxdmcoul=max(maxdm(ind2(sh1,sh2)),maxdm(ind2(sh3,sh4)))
    maxdmexch=max(maxdm(ind2(sh1,sh3)),maxdm(ind2(sh1,sh4)),maxdm(ind2(sh2,sh3)),&
       & maxdm(ind2(sh2,sh4)))
    maxdmval=max(maxdmcoul,.25_dp*maxdmexch)
    if (schwsp(ordersp(ii))*schwsp(ordersp(jj))*maxdmval.lt.thr) then
        nprsscr=nprsscr+1
        cycle
    end if
    call i2e_spsp(sh1,sh2,sh3,sh4,espsp)
    if (maxval(abs(espsp)).lt.thr) cycle
    ! Since the bra and ket shell pairs cannot both be the same, at most we can have
    !  (ab|ac) here, so no need to do anything special.
    do k=1,3
      ja=aosh(sh2)+k-1
      do l=1,3
        la=aosh(sh4)+l-1
        fac=2._dp*  p(la,ka) * espsp(k,l)
        floc(ia,ja) = floc(ia,ja) +  fac
        floc(ja,ia) = floc(ja,ia) +  fac
        fac=.5_dp*p(ja,ka)*espsp(k,l)
        floc(ia,la) = floc(ia,la) - fac
        floc(la,ia) = floc(la,ia) - fac
        fac=.5_dp*p(ja,la)*espsp(k,l)
        floc(ia,ka) = floc(ia,ka) - fac
        floc(ka,ia) = floc(ka,ia) - fac
        fac=.5_dp*p(ia,ka)*espsp(k,l)
        floc(ja,la) = floc(ja,la) - fac
        floc(la,ja) = floc(la,ja) - fac
        fac=.5_dp*p(ia,la)*espsp(k,l)
        floc(ja,ka) = floc(ja,ka) - fac
        floc(ka,ja) = floc(ka,ja) - fac
        fac=2._dp*p(ja,ia)*espsp(k,l)
        floc(ka,la) = floc(ka,la) +  fac
        floc(la,ka) = floc(la,ka) +   fac
      end do
    end do
  end do
end do
!$OMP END DO
!$OMP MASTER
if (prtlevl.ge.4) call time_diff_omp(1,"(sp|sp) integ:")
!$OMP END MASTER

!$OMP DO SCHEDULE(DYNAMIC)

do ii=1,nsp  !  (sp|sd)
  sh1=spprs(1,ordersp(ii))
  sh2=spprs(2,ordersp(ii))
  ab1=ind2(sh1,sh2)
  ia=aosh(sh1)
  do jj=1,nsd
    sh3=sdprs(1,ordersd(jj))
    sh4=sdprs(2,ordersd(jj))
    ab2=ind2(sh3,sh4)
    maxdmcoul=max(maxdm(ind2(sh1,sh2)),maxdm(ind2(sh3,sh4)))
    maxdmexch=max(maxdm(ind2(sh1,sh3)),maxdm(ind2(sh1,sh4)),maxdm(ind2(sh2,sh3)),&
       & maxdm(ind2(sh2,sh4)))
    maxdmval=max(maxdmcoul,.25_dp*maxdmexch)
    if (schwsp(ordersp(ii))*schwsd(ordersd(jj))*maxdmval.lt.thr) then
        nprsscr=nprsscr+1
        cycle
    end if
    ka=aosh(sh3)
    call i2e_spsd(sh1,sh2,sh3,sh4,espsd)
    if (maxval(abs(espsd)).lt.thr) cycle
    ! At most we can have (ab|ac) here, so no need to do anything special.
    do k=1,3
      ja=aosh(sh2)+k-1
      do l=1,5
        la=aosh(sh4)+l-1
        fac=2._dp*  p(la,ka) * espsd(k,l)
        floc(ia,ja) = floc(ia,ja) +  fac
        floc(ja,ia) = floc(ja,ia) +  fac
        fac=.5_dp*p(ja,ka)*espsd(k,l)
        floc(ia,la) = floc(ia,la) - fac
        floc(la,ia) = floc(la,ia) - fac
        fac=.5_dp*p(ja,la)*espsd(k,l)
        floc(ia,ka) = floc(ia,ka) - fac
        floc(ka,ia) = floc(ka,ia) - fac
        fac=.5_dp*p(ia,ka)*espsd(k,l)
        floc(ja,la) = floc(ja,la) - fac
        floc(la,ja) = floc(la,ja) - fac
        fac=.5_dp*p(ia,la)*espsd(k,l)
        floc(ja,ka) = floc(ja,ka) - fac
        floc(ka,ja) = floc(ka,ja) - fac
        fac=2._dp*p(ja,ia)*espsd(k,l)
        floc(ka,la) = floc(ka,la) +  fac
        floc(la,ka) = floc(la,ka) +   fac
      end do
    end do
  end do
end do
!$OMP END DO
!$OMP MASTER
if (prtlevl.ge.4) call time_diff_omp(1,"(sp|sd) integ:")
!$OMP END MASTER

!$OMP DO SCHEDULE(DYNAMIC)

do ii=1,nsp  !  (sp|sl)
  sh1=spprs(1,ordersp(ii))
  sh2=spprs(2,ordersp(ii))
  ab1=ind2(sh1,sh2)
  ia=aosh(sh1)
  do jj=1,nsl
    sh3=slprs(1,ordersl(jj))
    sh4=slprs(2,ordersl(jj))
    ab2=ind2(sh3,sh4)
    maxdmcoul=max(maxdm(ind2(sh1,sh2)),maxdm(ind2(sh3,sh4)))
    maxdmexch=max(maxdm(ind2(sh1,sh3)),maxdm(ind2(sh1,sh4)),maxdm(ind2(sh2,sh3)),&
       & maxdm(ind2(sh2,sh4)))
    maxdmval=max(maxdmcoul,.25_dp*maxdmexch)
    if (schwsp(ordersp(ii))*schwsl(ordersl(jj))*maxdmval.lt.thr) then
        nprsscr=nprsscr+1
        cycle
    end if
    ka=aosh(sh3)
    call i2e_spsl(sh1,sh2,sh3,sh4,espsl)
    if (maxval(abs(espsl)).lt.thr) cycle
    ! At most we can have (ab|ac) here, so no need to do anything special.
    do k=1,3
      ja=aosh(sh2)+k-1
      do l=1,4
        la=aosh(sh4)+l-1
        fac=2._dp*  p(la,ka) * espsl(k,l)
        floc(ia,ja) = floc(ia,ja) +  fac
        floc(ja,ia) = floc(ja,ia) +  fac
        fac=.5_dp*p(ja,ka)*espsl(k,l)
        floc(ia,la) = floc(ia,la) - fac
        floc(la,ia) = floc(la,ia) - fac
        fac=.5_dp*p(ja,la)*espsl(k,l)
        floc(ia,ka) = floc(ia,ka) - fac
        floc(ka,ia) = floc(ka,ia) - fac
        fac=.5_dp*p(ia,ka)*espsl(k,l)
        floc(ja,la) = floc(ja,la) - fac
        floc(la,ja) = floc(la,ja) - fac
        fac=.5_dp*p(ia,la)*espsl(k,l)
        floc(ja,ka) = floc(ja,ka) - fac
        floc(ka,ja) = floc(ka,ja) - fac
        fac=2._dp*p(ja,ia)*espsl(k,l)
        floc(ka,la) = floc(ka,la) +  fac
        floc(la,ka) = floc(la,ka) +   fac
      end do
    end do
  end do
end do
!$OMP END DO
!$OMP MASTER
if (prtlevl.ge.4) call time_diff_omp(1,"(sp|sl) integ:")
!$OMP END MASTER


!$OMP DO SCHEDULE(DYNAMIC)
do ii=1,nsp  ! (sp|pp)
  sh1=spprs(1,ordersp(ii))
  sh2=spprs(2,ordersp(ii))
  ab1=ind2(sh1,sh2)
  ia=aosh(sh1)
  do jj=1,npp
    sh3=ppprs(1,orderpp(jj))
    sh4=ppprs(2,orderpp(jj))
    ab2=ind2(sh3,sh4)
    maxdmcoul=max(maxdm(ind2(sh1,sh2)),maxdm(ind2(sh3,sh4)))
    maxdmexch=max(maxdm(ind2(sh1,sh3)),maxdm(ind2(sh1,sh4)),maxdm(ind2(sh2,sh3)),&
       & maxdm(ind2(sh2,sh4)))
    maxdmval=max(maxdmcoul,.25_dp*maxdmexch)
    if (schwsp(ordersp(ii))*schwpp(orderpp(jj))*maxdmval.lt.thr) then
       nprsscr=nprsscr+1
       cycle
    end if
    call i2e_sppp(sh1,sh2,sh3,sh4,esppp)
    if (maxval(abs(esppp)).lt.thr) cycle
    ! The only special case here is if shells sh3 and sh4 are the same.
    samekl=(sh3.eq.sh4)
    do j=1,3
      ja=aosh(sh2)+j-1
      do k=1,3
        ka=aosh(sh3)+k-1
        do l=1,3
          la=aosh(sh4)+l-1
          if (samekl.and.(k.eq.l)) then  ! when k=l these are (ab|cc)
            floc(ka,ka)=floc(ka,ka)+2._dp*p(ia,ja)*esppp(j,k,k)
            fac=.5_dp*p(ka,ia)*esppp(j,k,k)
            floc(ka,ja)=floc(ka,ja)-fac
            floc(ja,ka)=floc(ja,ka)-fac
            fac=.5_dp*p(ka,ja)*esppp(j,k,k)
            floc(ka,ia)=floc(ka,ia)-fac
            floc(ia,ka)=floc(ia,ka)-fac
            fac=p(ka,ka)*esppp(j,k,k)
            floc(ia,ja)=floc(ia,ja)+fac
            floc(ja,ia)=floc(ja,ia)+fac
          else if (samekl.and.(k.lt.l)) then
            cycle
          else
            fac=2._dp*  p(la,ka) * esppp(j,k,l)
            floc(ia,ja) = floc(ia,ja) +  fac
            floc(ja,ia) = floc(ja,ia) +  fac
            fac=.5_dp*p(ja,ka)*esppp(j,k,l)
            floc(ia,la) = floc(ia,la) - fac
            floc(la,ia) = floc(la,ia) - fac
            fac=.5_dp*p(ja,la)*esppp(j,k,l)
            floc(ia,ka) = floc(ia,ka) - fac
            floc(ka,ia) = floc(ka,ia) - fac
            fac=.5_dp*p(ia,ka)*esppp(j,k,l)
            floc(ja,la) = floc(ja,la) - fac
            floc(la,ja) = floc(la,ja) - fac
            fac=.5_dp*p(ia,la)*esppp(j,k,l)
            floc(ja,ka) = floc(ja,ka) - fac
            floc(ka,ja) = floc(ka,ja) - fac
            fac=2._dp*p(ja,ia)*esppp(j,k,l)
            floc(ka,la) = floc(ka,la) +  fac
            floc(la,ka) = floc(la,ka) +  fac
          end if
        end do
      end do
    end do
  end do
end do
!$OMP END DO
!$OMP MASTER
if (prtlevl.ge.4) call time_diff_omp(1,"(sp|pp) integ:")
!$OMP END MASTER


!$OMP DO SCHEDULE(DYNAMIC)

do ii=1,nsp  !  (sp|pd)
  sh1=spprs(1,ordersp(ii))
  sh2=spprs(2,ordersp(ii))
  ab1=ind2(sh1,sh2)
  ia=aosh(sh1)
  do jj=1,npd
    sh3=pdprs(1,orderpd(jj))
    sh4=pdprs(2,orderpd(jj))
    ab2=ind2(sh3,sh4)
    maxdmcoul=max(maxdm(ind2(sh1,sh2)),maxdm(ind2(sh3,sh4)))
    maxdmexch=max(maxdm(ind2(sh1,sh3)),maxdm(ind2(sh1,sh4)),maxdm(ind2(sh2,sh3)),&
       & maxdm(ind2(sh2,sh4)))
    maxdmval=max(maxdmcoul,.25_dp*maxdmexch)
    if (schwsp(ordersp(ii))*schwpd(orderpd(jj))*maxdmval.lt.thr) then
        nprsscr=nprsscr+1
        cycle
    end if
    ka=aosh(sh3)
    call i2e_sppd(sh1,sh2,sh3,sh4,esppd)
    if (maxval(abs(esppd)).lt.thr) cycle
    ! At most we can have (ab|ac) here, so no need to do anything special.
    do j=1,3
      ja=aosh(sh2)+j-1
      do k=1,3
        ka=aosh(sh3)+k-1
        do l=1,5
          la=aosh(sh4)+l-1
          fac=2._dp*  p(la,ka) * esppd(j,k,l)
          floc(ia,ja) = floc(ia,ja) +  fac
          floc(ja,ia) = floc(ja,ia) +  fac
          fac=.5_dp*p(ja,ka)*esppd(j,k,l)
          floc(ia,la) = floc(ia,la) - fac
          floc(la,ia) = floc(la,ia) - fac
          fac=.5_dp*p(ja,la)*esppd(j,k,l)
          floc(ia,ka) = floc(ia,ka) - fac
          floc(ka,ia) = floc(ka,ia) - fac
          fac=.5_dp*p(ia,ka)*esppd(j,k,l)
          floc(ja,la) = floc(ja,la) - fac
          floc(la,ja) = floc(la,ja) - fac
          fac=.5_dp*p(ia,la)*esppd(j,k,l)
          floc(ja,ka) = floc(ja,ka) - fac
          floc(ka,ja) = floc(ka,ja) - fac
          fac=2._dp*p(ja,ia)*esppd(j,k,l)
          floc(ka,la) = floc(ka,la) +  fac
          floc(la,ka) = floc(la,ka) +   fac
        end do
      end do
    end do
  end do
end do

!$OMP END DO
!$OMP MASTER
if (prtlevl.ge.4) call time_diff_omp(1,"(sp|pd) integ:")
!$OMP END MASTER


!$OMP DO SCHEDULE(DYNAMIC)
do ii=1,nsp  !  (sp|pl)
  sh1=spprs(1,ordersp(ii))
  sh2=spprs(2,ordersp(ii))
  ab1=ind2(sh1,sh2)
  ia=aosh(sh1)
  do jj=1,npl
    sh3=plprs(1,orderpl(jj))
    sh4=plprs(2,orderpl(jj))
    ab2=ind2(sh3,sh4)
    ka=aosh(sh3)
    maxdmcoul=max(maxdm(ind2(sh1,sh2)),maxdm(ind2(sh3,sh4)))
    maxdmexch=max(maxdm(ind2(sh1,sh3)),maxdm(ind2(sh1,sh4)),maxdm(ind2(sh2,sh3)),&
       & maxdm(ind2(sh2,sh4)))
    maxdmval=max(maxdmcoul,.25_dp*maxdmexch)
    if (schwsp(ordersp(ii))*schwpl(orderpl(jj))*maxdmval.lt.thr) then
        nprsscr=nprsscr+1
        cycle
    end if
    call i2e_sppl(sh1,sh2,sh3,sh4,esppl)
    if (maxval(abs(esppl)).lt.thr) cycle
    ! At most we can have (ab|ac) here, so no need to do anything special.
    do j=1,3
      ja=aosh(sh2)+j-1
      do k=1,3
        ka=aosh(sh3)+k-1
        do l=1,4
          la=aosh(sh4)+l-1
          fac=2._dp*  p(la,ka) * esppl(j,k,l)
          floc(ia,ja) = floc(ia,ja) +  fac
          floc(ja,ia) = floc(ja,ia) +  fac
          fac=.5_dp*p(ja,ka)*esppl(j,k,l)
          floc(ia,la) = floc(ia,la) - fac
          floc(la,ia) = floc(la,ia) - fac
          fac=.5_dp*p(ja,la)*esppl(j,k,l)
          floc(ia,ka) = floc(ia,ka) - fac
          floc(ka,ia) = floc(ka,ia) - fac
          fac=.5_dp*p(ia,ka)*esddl(j,k,l)
          floc(ja,la) = floc(ja,la) - fac
          floc(la,ja) = floc(la,ja) - fac
          fac=.5_dp*p(ia,la)*esddl(j,k,l)
          floc(ja,ka) = floc(ja,ka) - fac
          floc(ka,ja) = floc(ka,ja) - fac
          fac=2._dp*p(ja,ia)*esddl(j,k,l)
          floc(ka,la) = floc(ka,la) +  fac
          floc(la,ka) = floc(la,ka) +   fac
        end do
      end do
    end do
  end do
end do

!$OMP END DO
!$OMP MASTER
if (prtlevl.ge.4) call time_diff_omp(1,"(sd|dl) integ:")
!$OMP END MASTER


!$OMP DO SCHEDULE(DYNAMIC)
do ii=1,nsp  ! (sp|dd)
  sh1=spprs(1,ordersp(ii))
  sh2=spprs(2,ordersp(ii))
  ab1=ind2(sh1,sh2)
  ia=aosh(sh1)
  do jj=1,ndd
    sh3=ddprs(1,orderdd(jj))
    sh4=ddprs(2,orderdd(jj))
    ab2=ind2(sh3,sh4)
    maxdmcoul=max(maxdm(ind2(sh1,sh2)),maxdm(ind2(sh3,sh4)))
    maxdmexch=max(maxdm(ind2(sh1,sh3)),maxdm(ind2(sh1,sh4)),maxdm(ind2(sh2,sh3)),&
       & maxdm(ind2(sh2,sh4)))
    maxdmval=max(maxdmcoul,.25_dp*maxdmexch)
    if (schwsp(ordersp(ii))*schwdd(orderdd(jj))*maxdmval.lt.thr) then
       nprsscr=nprsscr+1
       cycle
    end if
    call i2e_spdd(sh1,sh2,sh3,sh4,espdd)
    if (maxval(abs(espdd)).lt.thr) cycle
    ! The only special case here is if shells sh3 and sh4 are the same.
    samekl=(sh3.eq.sh4)
    do j=1,3
      ja=aosh(sh2)+j-1
      do k=1,5
        ka=aosh(sh3)+k-1
        do l=1,5
          la=aosh(sh4)+l-1
          if (samekl.and.(k.eq.l)) then  ! when k=l these are (ab|cc)
            floc(ka,ka)=floc(ka,ka)+2._dp*p(ia,ja)*espdd(j,k,k)
            fac=.5_dp*p(ka,ia)*espdd(j,k,k)
            floc(ka,ja)=floc(ka,ja)-fac
            floc(ja,ka)=floc(ja,ka)-fac
            fac=.5_dp*p(ka,ja)*espdd(j,k,k)
            floc(ka,ia)=floc(ka,ia)-fac
            floc(ia,ka)=floc(ia,ka)-fac
            fac=p(ka,ka)*espdd(j,k,k)
            floc(ia,ja)=floc(ia,ja)+fac
            floc(ja,ia)=floc(ja,ia)+fac
          else if (samekl.and.(k.lt.l)) then
            cycle
          else
            fac=2._dp*  p(la,ka) * espdd(j,k,l)
            floc(ia,ja) = floc(ia,ja) +  fac
            floc(ja,ia) = floc(ja,ia) +  fac
            fac=.5_dp*p(ja,ka)*espdd(j,k,l)
            floc(ia,la) = floc(ia,la) - fac
            floc(la,ia) = floc(la,ia) - fac
            fac=.5_dp*p(ja,la)*espdd(j,k,l)
            floc(ia,ka) = floc(ia,ka) - fac
            floc(ka,ia) = floc(ka,ia) - fac
            fac=.5_dp*p(ia,ka)*espdd(j,k,l)
            floc(ja,la) = floc(ja,la) - fac
            floc(la,ja) = floc(la,ja) - fac
            fac=.5_dp*p(ia,la)*espdd(j,k,l)
            floc(ja,ka) = floc(ja,ka) - fac
            floc(ka,ja) = floc(ka,ja) - fac
            fac=2._dp*p(ja,ia)*espdd(j,k,l)
            floc(ka,la) = floc(ka,la) +  fac
            floc(la,ka) = floc(la,ka) +  fac
          end if
        end do
      end do
    end do
  end do
end do
!$OMP END DO
!$OMP MASTER
if (prtlevl.ge.4) call time_diff_omp(1,"(sp|dd) integ:")
!$OMP END MASTER


!$OMP DO SCHEDULE(DYNAMIC)

! necessary routine for (sp|dl) not yet written so skip it for now

do ii=1,nsp  ! (sp|ll)
  sh1=spprs(1,ordersp(ii))
  sh2=spprs(2,ordersp(ii))
  ab1=ind2(sh1,sh2)
  ia=aosh(sh1)
  do jj=1,nll
    sh3=llprs(1,orderll(jj))
    sh4=llprs(2,orderll(jj))
    ab2=ind2(sh3,sh4)
    maxdmcoul=max(maxdm(ind2(sh1,sh2)),maxdm(ind2(sh3,sh4)))
    maxdmexch=max(maxdm(ind2(sh1,sh3)),maxdm(ind2(sh1,sh4)),maxdm(ind2(sh2,sh3)),&
       & maxdm(ind2(sh2,sh4)))
    maxdmval=max(maxdmcoul,.25_dp*maxdmexch)
    if (schwsp(ordersp(ii))*schwll(orderll(jj))*maxdmval.lt.thr) then
       nprsscr=nprsscr+1
       cycle
    end if
    call i2e_spll(sh1,sh2,sh3,sh4,espll)
    if (maxval(abs(espll)).lt.thr) cycle
    ! The only special case here is if shells sh3 and sh4 are the same.
    samekl=(sh3.eq.sh4)
    do j=1,3
      ja=aosh(sh2)+j-1
      do k=1,4
        ka=aosh(sh3)+k-1
        do l=1,4
          la=aosh(sh4)+l-1
          if (samekl.and.(k.eq.l)) then  ! when k=l these are (ab|cc)
            floc(ka,ka)=floc(ka,ka)+2._dp*p(ia,ja)*espll(j,k,k)
            fac=.5_dp*p(ka,ia)*espll(j,k,k)
            floc(ka,ja)=floc(ka,ja)-fac
            floc(ja,ka)=floc(ja,ka)-fac
            fac=.5_dp*p(ka,ja)*espll(j,k,k)
            floc(ka,ia)=floc(ka,ia)-fac
            floc(ia,ka)=floc(ia,ka)-fac
            fac=p(ka,ka)*espll(j,k,k)
            floc(ia,ja)=floc(ia,ja)+fac
            floc(ja,ia)=floc(ja,ia)+fac
          else if (samekl.and.(k.lt.l)) then
            cycle
          else
            fac=2._dp*  p(la,ka) * espll(j,k,l)
            floc(ia,ja) = floc(ia,ja) +  fac
            floc(ja,ia) = floc(ja,ia) +  fac
            fac=.5_dp*p(ja,ka)*espll(j,k,l)
            floc(ia,la) = floc(ia,la) - fac
            floc(la,ia) = floc(la,ia) - fac
            fac=.5_dp*p(ja,la)*espll(j,k,l)
            floc(ia,ka) = floc(ia,ka) - fac
            floc(ka,ia) = floc(ka,ia) - fac
            fac=.5_dp*p(ia,ka)*espll(j,k,l)
            floc(ja,la) = floc(ja,la) - fac
            floc(la,ja) = floc(la,ja) - fac
            fac=.5_dp*p(ia,la)*espll(j,k,l)
            floc(ja,ka) = floc(ja,ka) - fac
            floc(ka,ja) = floc(ka,ja) - fac
            fac=2._dp*p(ja,ia)*espll(j,k,l)
            floc(ka,la) = floc(ka,la) +  fac
            floc(la,ka) = floc(la,ka) +  fac
          end if
        end do
      end do
    end do
  end do
end do

!$OMP END DO
!$OMP MASTER
if (prtlevl.ge.4) call time_diff_omp(1,"(sp|ll) integ:")
!$OMP END MASTER


!$OMP DO SCHEDULE(DYNAMIC)
do ii=2,nsd  !  (sd|sd)
  sh1=sdprs(1,ordersd(ii))
  sh2=sdprs(2,ordersd(ii))
  ab1=ind2(sh1,sh2)
  ia=aosh(sh1)
  do jj=1,ii-1
    sh3=sdprs(1,ordersd(jj))
    sh4=sdprs(2,ordersd(jj))
    ab2=ind2(sh3,sh4)
    ka=aosh(sh3)
    maxdmcoul=max(maxdm(ind2(sh1,sh2)),maxdm(ind2(sh3,sh4)))
    maxdmexch=max(maxdm(ind2(sh1,sh3)),maxdm(ind2(sh1,sh4)),maxdm(ind2(sh2,sh3)),&
       & maxdm(ind2(sh2,sh4)))
    maxdmval=max(maxdmcoul,.25_dp*maxdmexch)
    if (schwsd(ordersd(ii))*schwsd(ordersd(jj))*maxdmval.lt.thr) then
        nprsscr=nprsscr+1
        cycle
    end if
    call i2e_sdsd(sh1,sh2,sh3,sh4,esdsd)
    if (maxval(abs(esdsd)).lt.thr) cycle
    ! Since the bra and ket shell pairs cannot both be the same, at most we can have
    !  (ab|ac) here, so no need to do anything special.
    do k=1,5
      ja=aosh(sh2)+k-1
      do l=1,5
        la=aosh(sh4)+l-1
        fac=2._dp*  p(la,ka) * esdsd(k,l)
        floc(ia,ja) = floc(ia,ja) +  fac
        floc(ja,ia) = floc(ja,ia) +  fac
        fac=.5_dp*p(ja,ka)*esdsd(k,l)
        floc(ia,la) = floc(ia,la) - fac
        floc(la,ia) = floc(la,ia) - fac
        fac=.5_dp*p(ja,la)*esdsd(k,l)
        floc(ia,ka) = floc(ia,ka) - fac
        floc(ka,ia) = floc(ka,ia) - fac
        fac=.5_dp*p(ia,ka)*esdsd(k,l)
        floc(ja,la) = floc(ja,la) - fac
        floc(la,ja) = floc(la,ja) - fac
        fac=.5_dp*p(ia,la)*esdsd(k,l)
        floc(ja,ka) = floc(ja,ka) - fac
        floc(ka,ja) = floc(ka,ja) - fac
        fac=2._dp*p(ja,ia)*esdsd(k,l)
        floc(ka,la) = floc(ka,la) +  fac
        floc(la,ka) = floc(la,ka) +   fac
      end do
    end do
  end do
end do

!$OMP END DO
!$OMP MASTER
if (prtlevl.ge.4) call time_diff_omp(1,"(sd|sd) integ:")
!$OMP END MASTER

!$OMP DO SCHEDULE(DYNAMIC)
do ii=1,nsd  !  (sd|sl)
  sh1=sdprs(1,ordersd(ii))
  sh2=sdprs(2,ordersd(ii))
  ab1=ind2(sh1,sh2)
  ia=aosh(sh1)
  do jj=1,nsl
    sh3=slprs(1,ordersl(jj))
    sh4=slprs(2,ordersl(jj))
    ab2=ind2(sh3,sh4)
    maxdmcoul=max(maxdm(ind2(sh1,sh2)),maxdm(ind2(sh3,sh4)))
    maxdmexch=max(maxdm(ind2(sh1,sh3)),maxdm(ind2(sh1,sh4)),maxdm(ind2(sh2,sh3)),&
       & maxdm(ind2(sh2,sh4)))
    maxdmval=max(maxdmcoul,.25_dp*maxdmexch)
    if (schwsd(ordersd(ii))*schwsl(ordersl(jj))*maxdmval.lt.thr) then
        nprsscr=nprsscr+1
        cycle
    end if
    ka=aosh(sh3)
    call i2e_sdsl(sh1,sh2,sh3,sh4,esdsl)
    if (maxval(abs(esdsl)).lt.thr) cycle
    ! At most we can have (ab|ac) here, so no need to do anything special.
    do k=1,5
      ja=aosh(sh2)+k-1
      do l=1,4
        la=aosh(sh4)+l-1
        fac=2._dp*  p(la,ka) * esdsl(k,l)
        floc(ia,ja) = floc(ia,ja) +  fac
        floc(ja,ia) = floc(ja,ia) +  fac
        fac=.5_dp*p(ja,ka)*esdsl(k,l)
        floc(ia,la) = floc(ia,la) - fac
        floc(la,ia) = floc(la,ia) - fac
        fac=.5_dp*p(ja,la)*esdsl(k,l)
        floc(ia,ka) = floc(ia,ka) - fac
        floc(ka,ia) = floc(ka,ia) - fac
        fac=.5_dp*p(ia,ka)*esdsl(k,l)
        floc(ja,la) = floc(ja,la) - fac
        floc(la,ja) = floc(la,ja) - fac
        fac=.5_dp*p(ia,la)*esdsl(k,l)
        floc(ja,ka) = floc(ja,ka) - fac
        floc(ka,ja) = floc(ka,ja) - fac
        fac=2._dp*p(ja,ia)*esdsl(k,l)
        floc(ka,la) = floc(ka,la) +  fac
        floc(la,ka) = floc(la,ka) +   fac
      end do
    end do
  end do
end do

!$OMP END DO
!$OMP MASTER
if (prtlevl.ge.4) call time_diff_omp(1,"(sd|sl) integ:")
!$OMP END MASTER

!$OMP DO SCHEDULE(DYNAMIC)
do ii=1,nsd  ! (sd|pp)
  sh1=sdprs(1,ordersd(ii))
  sh2=sdprs(2,ordersd(ii))
  ab1=ind2(sh1,sh2)
  ia=aosh(sh1)
  do jj=1,npp
    sh3=ppprs(1,orderpp(jj))
    sh4=ppprs(2,orderpp(jj))
    ab2=ind2(sh3,sh4)
    maxdmcoul=max(maxdm(ind2(sh1,sh2)),maxdm(ind2(sh3,sh4)))
    maxdmexch=max(maxdm(ind2(sh1,sh3)),maxdm(ind2(sh1,sh4)),maxdm(ind2(sh2,sh3)),&
       & maxdm(ind2(sh2,sh4)))
    maxdmval=max(maxdmcoul,.25_dp*maxdmexch)
    if (schwsd(ordersd(ii))*schwpp(orderpp(jj))*maxdmval.lt.thr) then
       nprsscr=nprsscr+1
       cycle
    end if
    call i2e_sdpp(sh1,sh2,sh3,sh4,esdpp)
    if (maxval(abs(esdpp)).lt.thr) cycle
    ! The only special case here is if shells sh3 and sh4 are the same.
    samekl=(sh3.eq.sh4)
    do j=1,5
      ja=aosh(sh2)+j-1
      do k=1,3
        ka=aosh(sh3)+k-1
        do l=1,3
          la=aosh(sh4)+l-1
          if (samekl.and.(k.eq.l)) then  ! when k=l these are (ab|cc)
            floc(ka,ka)=floc(ka,ka)+2._dp*p(ia,ja)*esdpp(j,k,k)
            fac=.5_dp*p(ka,ia)*esdpp(j,k,k)
            floc(ka,ja)=floc(ka,ja)-fac
            floc(ja,ka)=floc(ja,ka)-fac
            fac=.5_dp*p(ka,ja)*esdpp(j,k,k)
            floc(ka,ia)=floc(ka,ia)-fac
            floc(ia,ka)=floc(ia,ka)-fac
            fac=p(ka,ka)*esdpp(j,k,k)
            floc(ia,ja)=floc(ia,ja)+fac
            floc(ja,ia)=floc(ja,ia)+fac
          else if (samekl.and.(k.lt.l)) then
            cycle
          else
            fac=2._dp*  p(la,ka) * esdpp(j,k,l)
            floc(ia,ja) = floc(ia,ja) +  fac
            floc(ja,ia) = floc(ja,ia) +  fac
            fac=.5_dp*p(ja,ka)*esdpp(j,k,l)
            floc(ia,la) = floc(ia,la) - fac
            floc(la,ia) = floc(la,ia) - fac
            fac=.5_dp*p(ja,la)*esdpp(j,k,l)
            floc(ia,ka) = floc(ia,ka) - fac
            floc(ka,ia) = floc(ka,ia) - fac
            fac=.5_dp*p(ia,ka)*esdpp(j,k,l)
            floc(ja,la) = floc(ja,la) - fac
            floc(la,ja) = floc(la,ja) - fac
            fac=.5_dp*p(ia,la)*esdpp(j,k,l)
            floc(ja,ka) = floc(ja,ka) - fac
            floc(ka,ja) = floc(ka,ja) - fac
            fac=2._dp*p(ja,ia)*esdpp(j,k,l)
            floc(ka,la) = floc(ka,la) +  fac
            floc(la,ka) = floc(la,ka) +  fac
          end if
        end do
      end do
    end do
  end do
end do

!$OMP END DO
!$OMP MASTER
if (prtlevl.ge.4) call time_diff_omp(1,"(sd|pp) integ:")
!$OMP END MASTER

!$OMP DO SCHEDULE(DYNAMIC)
do ii=1,nsd  !  (sd|pd)
  sh1=sdprs(1,ordersd(ii))
  sh2=sdprs(2,ordersd(ii))
  ab1=ind2(sh1,sh2)
  ia=aosh(sh1)
  do jj=1,npd
    sh3=pdprs(1,orderpd(jj))
    sh4=pdprs(2,orderpd(jj))
    ab2=ind2(sh3,sh4)
    ka=aosh(sh3)
    maxdmcoul=max(maxdm(ind2(sh1,sh2)),maxdm(ind2(sh3,sh4)))
    maxdmexch=max(maxdm(ind2(sh1,sh3)),maxdm(ind2(sh1,sh4)),maxdm(ind2(sh2,sh3)),&
       & maxdm(ind2(sh2,sh4)))
    maxdmval=max(maxdmcoul,.25_dp*maxdmexch)
    if (schwsd(ordersd(ii))*schwpd(orderpd(jj))*maxdmval.lt.thr) then
        nprsscr=nprsscr+1
        cycle
    end if
    call i2e_sdpd(sh1,sh2,sh3,sh4,esdpd)
    if (maxval(abs(esdpd)).lt.thr) cycle
    ! At most we can have (ab|ac) here, so no need to do anything special.
    do j=1,5
      ja=aosh(sh2)+j-1
      do k=1,3
        ka=aosh(sh3)+k-1
        do l=1,5
          la=aosh(sh4)+l-1
          fac=2._dp*  p(la,ka) * esdpd(j,k,l)
          floc(ia,ja) = floc(ia,ja) +  fac
          floc(ja,ia) = floc(ja,ia) +  fac
          fac=.5_dp*p(ja,ka)*esdpd(j,k,l)
          floc(ia,la) = floc(ia,la) - fac
          floc(la,ia) = floc(la,ia) - fac
          fac=.5_dp*p(ja,la)*esdpd(j,k,l)
          floc(ia,ka) = floc(ia,ka) - fac
          floc(ka,ia) = floc(ka,ia) - fac
          fac=.5_dp*p(ia,ka)*esdpd(j,k,l)
          floc(ja,la) = floc(ja,la) - fac
          floc(la,ja) = floc(la,ja) - fac
          fac=.5_dp*p(ia,la)*esdpd(j,k,l)
          floc(ja,ka) = floc(ja,ka) - fac
          floc(ka,ja) = floc(ka,ja) - fac
          fac=2._dp*p(ja,ia)*esdpd(j,k,l)
          floc(ka,la) = floc(ka,la) +  fac
          floc(la,ka) = floc(la,ka) +   fac
        end do
      end do
    end do
  end do
end do
!$OMP END DO
!$OMP MASTER
if (prtlevl.ge.4) call time_diff_omp(1,"(sd|pd) integ:")
!$OMP END MASTER

! Skipping (sd|pl) for now

!$OMP DO SCHEDULE(DYNAMIC)
do ii=1,nsd  ! (sd|dd)
  sh1=sdprs(1,ordersd(ii))
  sh2=sdprs(2,ordersd(ii))
  ab1=ind2(sh1,sh2)
  ia=aosh(sh1)
  do jj=1,ndd
    sh3=ddprs(1,orderdd(jj))
    sh4=ddprs(2,orderdd(jj))
    ab2=ind2(sh3,sh4)
    maxdmcoul=max(maxdm(ind2(sh1,sh2)),maxdm(ind2(sh3,sh4)))
    maxdmexch=max(maxdm(ind2(sh1,sh3)),maxdm(ind2(sh1,sh4)),maxdm(ind2(sh2,sh3)),&
       & maxdm(ind2(sh2,sh4)))
    maxdmval=max(maxdmcoul,.25_dp*maxdmexch)
    if (schwsd(ordersd(ii))*schwdd(orderdd(jj))*maxdmval.lt.thr) then
       nprsscr=nprsscr+1
       cycle
    end if
    call i2e_sddd(sh1,sh2,sh3,sh4,esddd)
    if (maxval(abs(esddd)).lt.thr) cycle
    ! The only special case here is if shells sh3 and sh4 are the same.
    samekl=(sh3.eq.sh4)
    do j=1,5
      ja=aosh(sh2)+j-1
      do k=1,5
        ka=aosh(sh3)+k-1
        do l=1,5
          la=aosh(sh4)+l-1
          if (samekl.and.(k.eq.l)) then  ! when k=l these are (ab|cc)
            floc(ka,ka)=floc(ka,ka)+2._dp*p(ia,ja)*esddd(j,k,k)
            fac=.5_dp*p(ka,ia)*esddd(j,k,k)
            floc(ka,ja)=floc(ka,ja)-fac
            floc(ja,ka)=floc(ja,ka)-fac
            fac=.5_dp*p(ka,ja)*esddd(j,k,k)
            floc(ka,ia)=floc(ka,ia)-fac
            floc(ia,ka)=floc(ia,ka)-fac
            fac=p(ka,ka)*esddd(j,k,k)
            floc(ia,ja)=floc(ia,ja)+fac
            floc(ja,ia)=floc(ja,ia)+fac
          else if (samekl.and.(k.lt.l)) then
            cycle
          else
            fac=2._dp*  p(la,ka) * esddd(j,k,l)
            floc(ia,ja) = floc(ia,ja) +  fac
            floc(ja,ia) = floc(ja,ia) +  fac
            fac=.5_dp*p(ja,ka)*esddd(j,k,l)
            floc(ia,la) = floc(ia,la) - fac
            floc(la,ia) = floc(la,ia) - fac
            fac=.5_dp*p(ja,la)*esddd(j,k,l)
            floc(ia,ka) = floc(ia,ka) - fac
            floc(ka,ia) = floc(ka,ia) - fac
            fac=.5_dp*p(ia,ka)*esddd(j,k,l)
            floc(ja,la) = floc(ja,la) - fac
            floc(la,ja) = floc(la,ja) - fac
            fac=.5_dp*p(ia,la)*esddd(j,k,l)
            floc(ja,ka) = floc(ja,ka) - fac
            floc(ka,ja) = floc(ka,ja) - fac
            fac=2._dp*p(ja,ia)*esddd(j,k,l)
            floc(ka,la) = floc(ka,la) +  fac
            floc(la,ka) = floc(la,ka) +  fac
          end if
        end do
      end do
    end do
  end do
end do
!$OMP END DO
!$OMP MASTER
if (prtlevl.ge.4) call time_diff_omp(1,"(sd|dd) integ:")
!$OMP END MASTER

!$OMP DO SCHEDULE(DYNAMIC)
do ii=1,nsd  !  (sd|dl)
  sh1=sdprs(1,ordersd(ii))
  sh2=sdprs(2,ordersd(ii))
  ab1=ind2(sh1,sh2)
  ia=aosh(sh1)
  do jj=1,ndl
    sh3=dlprs(1,orderdl(jj))
    sh4=dlprs(2,orderdl(jj))
    ab2=ind2(sh3,sh4)
    ka=aosh(sh3)
    maxdmcoul=max(maxdm(ind2(sh1,sh2)),maxdm(ind2(sh3,sh4)))
    maxdmexch=max(maxdm(ind2(sh1,sh3)),maxdm(ind2(sh1,sh4)),maxdm(ind2(sh2,sh3)),&
       & maxdm(ind2(sh2,sh4)))
    maxdmval=max(maxdmcoul,.25_dp*maxdmexch)
    if (schwsd(ordersd(ii))*schwdl(orderdl(jj))*maxdmval.lt.thr) then
        nprsscr=nprsscr+1
        cycle
    end if
    call i2e_sddl(sh1,sh2,sh3,sh4,esddl)
    if (maxval(abs(esddl)).lt.thr) cycle
    ! At most we can have (ab|ac) here, so no need to do anything special.
    do j=1,5
      ja=aosh(sh2)+j-1
      do k=1,5
        ka=aosh(sh3)+k-1
        do l=1,4
          la=aosh(sh4)+l-1
          fac=2._dp*  p(la,ka) * esddl(j,k,l)
          floc(ia,ja) = floc(ia,ja) +  fac
          floc(ja,ia) = floc(ja,ia) +  fac
          fac=.5_dp*p(ja,ka)*esddl(j,k,l)
          floc(ia,la) = floc(ia,la) - fac
          floc(la,ia) = floc(la,ia) - fac
          fac=.5_dp*p(ja,la)*esddl(j,k,l)
          floc(ia,ka) = floc(ia,ka) - fac
          floc(ka,ia) = floc(ka,ia) - fac
          fac=.5_dp*p(ia,ka)*esddl(j,k,l)
          floc(ja,la) = floc(ja,la) - fac
          floc(la,ja) = floc(la,ja) - fac
          fac=.5_dp*p(ia,la)*esddl(j,k,l)
          floc(ja,ka) = floc(ja,ka) - fac
          floc(ka,ja) = floc(ka,ja) - fac
          fac=2._dp*p(ja,ia)*esddl(j,k,l)
          floc(ka,la) = floc(ka,la) +  fac
          floc(la,ka) = floc(la,ka) +   fac
        end do
      end do
    end do
  end do
end do
!$OMP END DO
!$OMP MASTER
if (prtlevl.ge.4) call time_diff_omp(1,"(sd|dl) integ:")
!$OMP END MASTER

!$OMP DO SCHEDULE(DYNAMIC)
do ii=1,nsd  ! (sd|ll)
  sh1=sdprs(1,ordersd(ii))
  sh2=sdprs(2,ordersd(ii))
  ab1=ind2(sh1,sh2)
  ia=aosh(sh1)
  do jj=1,nll
    sh3=llprs(1,orderll(jj))
    sh4=llprs(2,orderll(jj))
    ab2=ind2(sh3,sh4)
    maxdmcoul=max(maxdm(ind2(sh1,sh2)),maxdm(ind2(sh3,sh4)))
    maxdmexch=max(maxdm(ind2(sh1,sh3)),maxdm(ind2(sh1,sh4)),maxdm(ind2(sh2,sh3)),&
       & maxdm(ind2(sh2,sh4)))
    maxdmval=max(maxdmcoul,.25_dp*maxdmexch)
    if (schwsd(ordersd(ii))*schwll(orderll(jj))*maxdmval.lt.thr) then
       nprsscr=nprsscr+1
       cycle
    end if
    call i2e_sdll(sh1,sh2,sh3,sh4,esdll)
    if (maxval(abs(esdll)).lt.thr) cycle
    ! The only special case here is if shells sh3 and sh4 are the same.
    samekl=(sh3.eq.sh4)
    do j=1,5
      ja=aosh(sh2)+j-1
      do k=1,4
        ka=aosh(sh3)+k-1
        do l=1,4
          la=aosh(sh4)+l-1
          if (samekl.and.(k.eq.l)) then  ! when k=l these are (ab|cc)
            floc(ka,ka)=floc(ka,ka)+2._dp*p(ia,ja)*esdll(j,k,k)
            fac=.5_dp*p(ka,ia)*esdll(j,k,k)
            floc(ka,ja)=floc(ka,ja)-fac
            floc(ja,ka)=floc(ja,ka)-fac
            fac=.5_dp*p(ka,ja)*esdll(j,k,k)
            floc(ka,ia)=floc(ka,ia)-fac
            floc(ia,ka)=floc(ia,ka)-fac
            fac=p(ka,ka)*esdll(j,k,k)
            floc(ia,ja)=floc(ia,ja)+fac
            floc(ja,ia)=floc(ja,ia)+fac
          else if (samekl.and.(k.lt.l)) then
            cycle
          else
            fac=2._dp*  p(la,ka) * esdll(j,k,l)
            floc(ia,ja) = floc(ia,ja) +  fac
            floc(ja,ia) = floc(ja,ia) +  fac
            fac=.5_dp*p(ja,ka)*esdll(j,k,l)
            floc(ia,la) = floc(ia,la) - fac
            floc(la,ia) = floc(la,ia) - fac
            fac=.5_dp*p(ja,la)*esdll(j,k,l)
            floc(ia,ka) = floc(ia,ka) - fac
            floc(ka,ia) = floc(ka,ia) - fac
            fac=.5_dp*p(ia,ka)*esdll(j,k,l)
            floc(ja,la) = floc(ja,la) - fac
            floc(la,ja) = floc(la,ja) - fac
            fac=.5_dp*p(ia,la)*esdll(j,k,l)
            floc(ja,ka) = floc(ja,ka) - fac
            floc(ka,ja) = floc(ka,ja) - fac
            fac=2._dp*p(ja,ia)*esdll(j,k,l)
            floc(ka,la) = floc(ka,la) +  fac
            floc(la,ka) = floc(la,ka) +  fac
          end if
        end do
      end do
    end do
  end do
end do
!$OMP END DO
!$OMP MASTER
if (prtlevl.ge.4) call time_diff_omp(1,"(sd|ll) integ:")
!$OMP END MASTER

!$OMP DO SCHEDULE(DYNAMIC)
do ii=2,nsl  !  (sl|sl)
  sh1=slprs(1,ordersl(ii))
  sh2=slprs(2,ordersl(ii))
  ab1=ind2(sh1,sh2)
  ia=aosh(sh1)
  do jj=1,ii-1
    sh3=slprs(1,ordersl(jj))
    sh4=slprs(2,ordersl(jj))
    ab2=ind2(sh3,sh4)
    ka=aosh(sh3)
    maxdmcoul=max(maxdm(ind2(sh1,sh2)),maxdm(ind2(sh3,sh4)))
    maxdmexch=max(maxdm(ind2(sh1,sh3)),maxdm(ind2(sh1,sh4)),maxdm(ind2(sh2,sh3)),&
       & maxdm(ind2(sh2,sh4)))
    maxdmval=max(maxdmcoul,.25_dp*maxdmexch)
    if (schwsl(ordersl(ii))*schwsl(ordersl(jj))*maxdmval.lt.thr) then
        nprsscr=nprsscr+1
        cycle
    end if
    call i2e_slsl(sh1,sh2,sh3,sh4,eslsl)
    if (maxval(abs(eslsl)).lt.thr) cycle
    ! Since the bra and ket shell pairs cannot both be the same, at most we can have
    !  (ab|ac) here, so no need to do anything special.
    do k=1,4
      ja=aosh(sh2)+k-1
      do l=1,4
        la=aosh(sh4)+l-1
        fac=2._dp*  p(la,ka) * eslsl(k,l)
        floc(ia,ja) = floc(ia,ja) +  fac
        floc(ja,ia) = floc(ja,ia) +  fac
        fac=.5_dp*p(ja,ka)*eslsl(k,l)
        floc(ia,la) = floc(ia,la) - fac
        floc(la,ia) = floc(la,ia) - fac
        fac=.5_dp*p(ja,la)*eslsl(k,l)
        floc(ia,ka) = floc(ia,ka) - fac
        floc(ka,ia) = floc(ka,ia) - fac
        fac=.5_dp*p(ia,ka)*eslsl(k,l)
        floc(ja,la) = floc(ja,la) - fac
        floc(la,ja) = floc(la,ja) - fac
        fac=.5_dp*p(ia,la)*eslsl(k,l)
        floc(ja,ka) = floc(ja,ka) - fac
        floc(ka,ja) = floc(ka,ja) - fac
        fac=2._dp*p(ja,ia)*eslsl(k,l)
        floc(ka,la) = floc(ka,la) +  fac
        floc(la,ka) = floc(la,ka) +   fac
      end do
    end do
  end do
end do
!$OMP END DO
!$OMP MASTER
if (prtlevl.ge.4) call time_diff_omp(1,"(sl|sl) integ:")
!$OMP END MASTER

!$OMP DO SCHEDULE(DYNAMIC)
do ii=1,nsl  ! (sl|pp)
  sh1=slprs(1,ordersl(ii))
  sh2=slprs(2,ordersl(ii))
  ab1=ind2(sh1,sh2)
  ia=aosh(sh1)
  do jj=1,npp
    sh3=ppprs(1,orderpp(jj))
    sh4=ppprs(2,orderpp(jj))
    ab2=ind2(sh3,sh4)
    maxdmcoul=max(maxdm(ind2(sh1,sh2)),maxdm(ind2(sh3,sh4)))
    maxdmexch=max(maxdm(ind2(sh1,sh3)),maxdm(ind2(sh1,sh4)),maxdm(ind2(sh2,sh3)),&
       & maxdm(ind2(sh2,sh4)))
    maxdmval=max(maxdmcoul,.25_dp*maxdmexch)
    if (schwsl(ordersl(ii))*schwpp(orderpp(jj))*maxdmval.lt.thr) then
       nprsscr=nprsscr+1
       cycle
    end if
    call i2e_slpp(sh1,sh2,sh3,sh4,eslpp)
    if (maxval(abs(eslpp)).lt.thr) cycle
    ! The only special case here is if shells sh3 and sh4 are the same.
    samekl=(sh3.eq.sh4)
    do j=1,4
      ja=aosh(sh2)+j-1
      do k=1,3
        ka=aosh(sh3)+k-1
        do l=1,3
          la=aosh(sh4)+l-1
          if (samekl.and.(k.eq.l)) then  ! when k=l these are (ab|cc)
            floc(ka,ka)=floc(ka,ka)+2._dp*p(ia,ja)*eslpp(j,k,k)
            fac=.5_dp*p(ka,ia)*eslpp(j,k,k)
            floc(ka,ja)=floc(ka,ja)-fac
            floc(ja,ka)=floc(ja,ka)-fac
            fac=.5_dp*p(ka,ja)*eslpp(j,k,k)
            floc(ka,ia)=floc(ka,ia)-fac
            floc(ia,ka)=floc(ia,ka)-fac
            fac=p(ka,ka)*eslpp(j,k,k)
            floc(ia,ja)=floc(ia,ja)+fac
            floc(ja,ia)=floc(ja,ia)+fac
          else if (samekl.and.(k.lt.l)) then
            cycle
          else
            fac=2._dp*  p(la,ka) * eslpp(j,k,l)
            floc(ia,ja) = floc(ia,ja) +  fac
            floc(ja,ia) = floc(ja,ia) +  fac
            fac=.5_dp*p(ja,ka)*eslpp(j,k,l)
            floc(ia,la) = floc(ia,la) - fac
            floc(la,ia) = floc(la,ia) - fac
            fac=.5_dp*p(ja,la)*eslpp(j,k,l)
            floc(ia,ka) = floc(ia,ka) - fac
            floc(ka,ia) = floc(ka,ia) - fac
            fac=.5_dp*p(ia,ka)*eslpp(j,k,l)
            floc(ja,la) = floc(ja,la) - fac
            floc(la,ja) = floc(la,ja) - fac
            fac=.5_dp*p(ia,la)*eslpp(j,k,l)
            floc(ja,ka) = floc(ja,ka) - fac
            floc(ka,ja) = floc(ka,ja) - fac
            fac=2._dp*p(ja,ia)*eslpp(j,k,l)
            floc(ka,la) = floc(ka,la) +  fac
            floc(la,ka) = floc(la,ka) +  fac
          end if
        end do
      end do
    end do
  end do
end do
!$OMP END DO
!$OMP MASTER
if (prtlevl.ge.4) call time_diff_omp(1,"(sl|pp) integ:")
!$OMP END MASTER

! Skipping (sl|pd)

!$OMP DO SCHEDULE(DYNAMIC)
do ii=1,nsl  !  (sl|pl)
  sh1=slprs(1,ordersl(ii))
  sh2=slprs(2,ordersl(ii))
  ab1=ind2(sh1,sh2)
  ia=aosh(sh1)
  do jj=1,npl 
    sh3=plprs(1,orderpl(jj))
    sh4=plprs(2,orderpl(jj))
    ab2=ind2(sh3,sh4)
    ka=aosh(sh3)
    maxdmcoul=max(maxdm(ind2(sh1,sh2)),maxdm(ind2(sh3,sh4)))
    maxdmexch=max(maxdm(ind2(sh1,sh3)),maxdm(ind2(sh1,sh4)),maxdm(ind2(sh2,sh3)),&
       & maxdm(ind2(sh2,sh4)))
    maxdmval=max(maxdmcoul,.25_dp*maxdmexch)
    if (schwsl(ordersl(ii))*schwpl(orderpl(jj))*maxdmval.lt.thr) then
        nprsscr=nprsscr+1
        cycle
    end if
    call i2e_slpl(sh1,sh2,sh3,sh4,eslpl)
    if (maxval(abs(eslpl)).lt.thr) cycle
    ! at most we can have (ab|ac) here, so no need to do anything special.
    do j=1,4
      ja=aosh(sh2)+j-1
      do k=1,3
        ka=aosh(sh3)+k-1
        do l=1,4
          la=aosh(sh4)+l-1
          fac=2._dp*  p(la,ka) * eslpl(j,k,l)
          floc(ia,ja) = floc(ia,ja) +  fac
          floc(ja,ia) = floc(ja,ia) +  fac
          fac=.5_dp*p(ja,ka)*eslpl(j,k,l)
          floc(ia,la) = floc(ia,la) - fac
          floc(la,ia) = floc(la,ia) - fac
          fac=.5_dp*p(ja,la)*eslpl(j,k,l)
          floc(ia,ka) = floc(ia,ka) - fac
          floc(ka,ia) = floc(ka,ia) - fac
          fac=.5_dp*p(ia,ka)*eslpl(j,k,l)
          floc(ja,la) = floc(ja,la) - fac
          floc(la,ja) = floc(la,ja) - fac
          fac=.5_dp*p(ia,la)*eslpl(j,k,l)
          floc(ja,ka) = floc(ja,ka) - fac
          floc(ka,ja) = floc(ka,ja) - fac
          fac=2._dp*p(ja,ia)*eslpl(j,k,l)
          floc(ka,la) = floc(ka,la) +  fac
          floc(la,ka) = floc(la,ka) +   fac
        end do
      end do
    end do
  end do
end do
!$OMP END DO
!$OMP MASTER
if (prtlevl.ge.4) call time_diff_omp(1,"(sl|pl) integ:")
!$OMP END MASTER

!$OMP DO SCHEDULE(DYNAMIC)
do ii=1,nsl  ! (sl|dd)
  sh1=slprs(1,ordersl(ii))
  sh2=slprs(2,ordersl(ii))
  ab1=ind2(sh1,sh2)
  ia=aosh(sh1)
  do jj=1,ndd
    sh3=ddprs(1,orderdd(jj))
    sh4=ddprs(2,orderdd(jj))
    ab2=ind2(sh3,sh4)
    maxdmcoul=max(maxdm(ind2(sh1,sh2)),maxdm(ind2(sh3,sh4)))
    maxdmexch=max(maxdm(ind2(sh1,sh3)),maxdm(ind2(sh1,sh4)),maxdm(ind2(sh2,sh3)),&
       & maxdm(ind2(sh2,sh4)))
    maxdmval=max(maxdmcoul,.25_dp*maxdmexch)
    if (schwsl(ordersl(ii))*schwdd(orderdd(jj))*maxdmval.lt.thr) then
       nprsscr=nprsscr+1
       cycle
    end if
    call i2e_sldd(sh1,sh2,sh3,sh4,esldd)
    if (maxval(abs(esldd)).lt.thr) cycle
    ! The only special case here is if shells sh3 and sh4 are the same.
    samekl=(sh3.eq.sh4)
    do j=1,4
      ja=aosh(sh2)+j-1
      do k=1,5
        ka=aosh(sh3)+k-1
        do l=1,5
          la=aosh(sh4)+l-1
          if (samekl.and.(k.eq.l)) then  ! when k=l these are (ab|cc)
            floc(ka,ka)=floc(ka,ka)+2._dp*p(ia,ja)*esldd(j,k,k)
            fac=.5_dp*p(ka,ia)*esldd(j,k,k)
            floc(ka,ja)=floc(ka,ja)-fac
            floc(ja,ka)=floc(ja,ka)-fac
            fac=.5_dp*p(ka,ja)*esldd(j,k,k)
            floc(ka,ia)=floc(ka,ia)-fac
            floc(ia,ka)=floc(ia,ka)-fac
            fac=p(ka,ka)*esldd(j,k,k)
            floc(ia,ja)=floc(ia,ja)+fac
            floc(ja,ia)=floc(ja,ia)+fac
          else if (samekl.and.(k.lt.l)) then
            cycle
          else
            fac=2._dp*  p(la,ka) * esldd(j,k,l)
            floc(ia,ja) = floc(ia,ja) +  fac
            floc(ja,ia) = floc(ja,ia) +  fac
            fac=.5_dp*p(ja,ka)*esldd(j,k,l)
            floc(ia,la) = floc(ia,la) - fac
            floc(la,ia) = floc(la,ia) - fac
            fac=.5_dp*p(ja,la)*esldd(j,k,l)
            floc(ia,ka) = floc(ia,ka) - fac
            floc(ka,ia) = floc(ka,ia) - fac
            fac=.5_dp*p(ia,ka)*esldd(j,k,l)
            floc(ja,la) = floc(ja,la) - fac
            floc(la,ja) = floc(la,ja) - fac
            fac=.5_dp*p(ia,la)*esldd(j,k,l)
            floc(ja,ka) = floc(ja,ka) - fac
            floc(ka,ja) = floc(ka,ja) - fac
            fac=2._dp*p(ja,ia)*esldd(j,k,l)
            floc(ka,la) = floc(ka,la) +  fac
            floc(la,ka) = floc(la,ka) +  fac
          end if
        end do
      end do
    end do
  end do
end do
!$OMP END DO
!$OMP MASTER
if (prtlevl.ge.4) call time_diff_omp(1,"(sl|dd) integ:")
!$OMP END MASTER

!$OMP DO SCHEDULE(DYNAMIC)
do ii=1,nsl  !  (sl|dl)
  sh1=slprs(1,ordersl(ii))
  sh2=slprs(2,ordersl(ii))
  ab1=ind2(sh1,sh2)
  ia=aosh(sh1)
  do jj=1,ndl
    sh3=dlprs(1,orderdl(jj))
    sh4=dlprs(2,orderdl(jj))
    ab2=ind2(sh3,sh4)
    ka=aosh(sh3)
    maxdmcoul=max(maxdm(ind2(sh1,sh2)),maxdm(ind2(sh3,sh4)))
    maxdmexch=max(maxdm(ind2(sh1,sh3)),maxdm(ind2(sh1,sh4)),maxdm(ind2(sh2,sh3)),&
       & maxdm(ind2(sh2,sh4)))
    maxdmval=max(maxdmcoul,.25_dp*maxdmexch)
    if (schwsl(ordersl(ii))*schwdl(orderdl(jj))*maxdmval.lt.thr) then
        nprsscr=nprsscr+1
        cycle
    end if
    call i2e_sldl(sh1,sh2,sh3,sh4,esldl)
    if (maxval(abs(esldl)).lt.thr) cycle
    ! at most we can have (ab|ac) here, so no need to do anything special.
    do j=1,4
      ja=aosh(sh2)+j-1
      do k=1,5
        ka=aosh(sh3)+k-1
        do l=1,4
          la=aosh(sh4)+l-1
          fac=2._dp*  p(la,ka) * esldl(j,k,l)
          floc(ia,ja) = floc(ia,ja) +  fac
          floc(ja,ia) = floc(ja,ia) +  fac
          fac=.5_dp*p(ja,ka)*esldl(j,k,l)
          floc(ia,la) = floc(ia,la) - fac
          floc(la,ia) = floc(la,ia) - fac
          fac=.5_dp*p(ja,la)*esldl(j,k,l)
          floc(ia,ka) = floc(ia,ka) - fac
          floc(ka,ia) = floc(ka,ia) - fac
          fac=.5_dp*p(ia,ka)*esldl(j,k,l)
          floc(ja,la) = floc(ja,la) - fac
          floc(la,ja) = floc(la,ja) - fac
          fac=.5_dp*p(ia,la)*esldl(j,k,l)
          floc(ja,ka) = floc(ja,ka) - fac
          floc(ka,ja) = floc(ka,ja) - fac
          fac=2._dp*p(ja,ia)*esldl(j,k,l)
          floc(ka,la) = floc(ka,la) +  fac
          floc(la,ka) = floc(la,ka) +   fac
        end do
      end do
    end do
  end do
end do
!$OMP END DO
!$OMP MASTER
if (prtlevl.ge.4) call time_diff_omp(1,"(sl|dl) integ:")
!$OMP END MASTER

!$OMP DO SCHEDULE(DYNAMIC)
do ii=1,nsl  ! (sl|ll)
  sh1=slprs(1,ordersl(ii))
  sh2=slprs(2,ordersl(ii))
  ab1=ind2(sh1,sh2)
  ia=aosh(sh1)
  do jj=1,nll
    sh3=llprs(1,orderll(jj))
    sh4=llprs(2,orderll(jj))
    ab2=ind2(sh3,sh4)
    maxdmcoul=max(maxdm(ind2(sh1,sh2)),maxdm(ind2(sh3,sh4)))
    maxdmexch=max(maxdm(ind2(sh1,sh3)),maxdm(ind2(sh1,sh4)),maxdm(ind2(sh2,sh3)),&
       & maxdm(ind2(sh2,sh4)))
    maxdmval=max(maxdmcoul,.25_dp*maxdmexch)
    if (schwsl(ordersl(ii))*schwll(orderll(jj))*maxdmval.lt.thr) then
        nprsscr=nprsscr+1
        cycle
    end if
    call i2e_slll(sh1,sh2,sh3,sh4,eslll)
    if (maxval(abs(eslll)).lt.thr) cycle
    ! The only special case here is if shells sh3 and sh4 are the same.
    samekl=(sh3.eq.sh4)
    do j=1,4
      ja=aosh(sh2)+j-1
      do k=1,4
        ka=aosh(sh3)+k-1
        do l=1,4
          la=aosh(sh4)+l-1
          if (samekl.and.(k.eq.l)) then  ! when k=l these are (ab|cc)
            floc(ka,ka)=floc(ka,ka)+2._dp*p(ia,ja)*eslll(j,k,k)
            fac=.5_dp*p(ka,ia)*eslll(j,k,k)
            floc(ka,ja)=floc(ka,ja)-fac
            floc(ja,ka)=floc(ja,ka)-fac
            fac=.5_dp*p(ka,ja)*eslll(j,k,k)
            floc(ka,ia)=floc(ka,ia)-fac
            floc(ia,ka)=floc(ia,ka)-fac
            fac=p(ka,ka)*eslll(j,k,k)
            floc(ia,ja)=floc(ia,ja)+fac
            floc(ja,ia)=floc(ja,ia)+fac
          else if (samekl.and.(k.lt.l)) then
            cycle
          else
            fac=2._dp*  p(la,ka) * eslll(j,k,l)
            floc(ia,ja) = floc(ia,ja) +  fac
            floc(ja,ia) = floc(ja,ia) +  fac
            fac=.5_dp*p(ja,ka)*eslll(j,k,l)
            floc(ia,la) = floc(ia,la) - fac
            floc(la,ia) = floc(la,ia) - fac
            fac=.5_dp*p(ja,la)*eslll(j,k,l)
            floc(ia,ka) = floc(ia,ka) - fac
            floc(ka,ia) = floc(ka,ia) - fac
            fac=.5_dp*p(ia,ka)*eslll(j,k,l)
            floc(ja,la) = floc(ja,la) - fac
            floc(la,ja) = floc(la,ja) - fac
            fac=.5_dp*p(ia,la)*eslll(j,k,l)
            floc(ja,ka) = floc(ja,ka) - fac
            floc(ka,ja) = floc(ka,ja) - fac
            fac=2._dp*p(ja,ia)*eslll(j,k,l)
            floc(ka,la) = floc(ka,la) +  fac
            floc(la,ka) = floc(la,ka) +  fac
          end if
        end do
      end do
    end do
  end do
end do
!$OMP END DO
!$OMP MASTER
if (prtlevl.ge.4) call time_diff_omp(1,"(sl|ll) integ:")
!$OMP END MASTER

!$OMP DO SCHEDULE(DYNAMIC)
do ii=2,npp  ! (pp|pp)
  sh1=ppprs(1,orderpp(ii))
  sh2=ppprs(2,orderpp(ii))
  ab1=ind2(sh1,sh2)
  sameij=(sh1.eq.sh2)
  do jj=1,ii-1
    sh3=ppprs(1,orderpp(jj))
    sh4=ppprs(2,orderpp(jj))
    ab2=ind2(sh3,sh4)
    maxdmcoul=max(maxdm(ind2(sh1,sh2)),maxdm(ind2(sh3,sh4)))
    maxdmexch=max(maxdm(ind2(sh1,sh3)),maxdm(ind2(sh1,sh4)),maxdm(ind2(sh2,sh3)),&
       & maxdm(ind2(sh2,sh4)))
    maxdmval=max(maxdmcoul,.25_dp*maxdmexch)
    if (schwpp(orderpp(ii))*schwpp(orderpp(jj))*maxdmval.lt.thr) then
        nprsscr=nprsscr+1
        cycle
    end if
    call i2e_pppp(sh1,sh2,sh3,sh4,epppp)
    if (maxval(abs(epppp)).lt.thr) cycle
    ! The special cases here are ij i=j, and/or k=l. It cannot be the case that i=j=k=l or that i=k and j=l, so (aa|aa) and (ab|ab) are ruled out
    samekl=(sh3.eq.sh4)
    do i=1,3
      ia=aosh(sh1)+i-1
      do j=1,3
        ja=aosh(sh2)+j-1
        do k=1,3
          ka=aosh(sh3)+k-1
          do l=1,3
            la=aosh(sh4)+l-1
            if (sameij.and.(i.eq.j).and.samekl.and.(k.eq.l)) then  ! (aa|bb)
              floc(ia,ia)=floc(ia,ia)+p(ka,ka)*epppp(i,i,k,k)
              fac=.5_dp*p(ia,ka)*epppp(i,i,k,k)
              floc(ia,ka)=floc(ia,ka)-fac
              floc(ka,ia)=floc(ka,ia)-fac
              floc(ka,ka)=floc(ka,ka)+p(ia,ia)*epppp(i,i,k,k)
              !if (bla) write (*,*) "carried out (aa|bb) for: ",sh1,sh2,sh3,sh4,"i,j,k,l",i,j,k,l
            else if (sameij.and.(i.eq.j).and.samekl.and.(k.lt.l)) then  ! (aa|cb)
              !if (bla) write (*,*) "done nothing (aa|cb) for: ",sh1,sh2,sh3,sh4,"i,j,k,l",i,j,k,l
              cycle
            else if (sameij.and.(i.eq.j).and.samekl.and.(k.gt.l)) then  !  (aa|bc)
              floc(ia,ia)=floc(ia,ia)+2._dp*p(ka,la)*epppp(i,j,k,l)
              fac=.5_dp*p(ia,ka)*epppp(i,j,k,l)
              floc(ia,la)=floc(ia,la)-fac
              floc(la,ia)=floc(la,ia)-fac
              fac=.5_dp*p(ia,la)*epppp(i,j,k,l)
              floc(ia,ka)=floc(ia,ka)-fac
              floc(ka,ia)=floc(ka,ia)-fac
              fac=p(ia,ia)*epppp(i,j,k,l)
              floc(ka,la)=floc(ka,la)+fac
              floc(la,ka)=floc(la,ka)+fac
              !if (bla) write (*,*) "done (aa|bc) for: ",sh1,sh2,sh3,sh4,"i,j,k,l",i,j,k,l
            else if (sameij.and.(i.lt.j).and.samekl.and.(k.eq.l)) then    !  (ba|cc)
              !if (bla) write (*,*) "done nothing (ba|cc) for: ",sh1,sh2,sh3,sh4,"i,j,k,l",i,j,k,l
              cycle
            else if (samekl.and.(k.eq.l).and.sameij.and.(i.gt.j)) then    ! (ab|cc)
              floc(ka,ka)=floc(ka,ka)+2._dp*p(ia,ja)*epppp(i,j,k,l)
              fac=.5_dp*p(ka,ia)*epppp(i,j,k,l)
              floc(ka,ja)=floc(ka,ja)-fac
              floc(ja,ka)=floc(ja,ka)-fac
              fac=.5_dp*p(ka,ja)*epppp(i,j,k,l)
              floc(ka,ia)=floc(ka,ia)-fac
              floc(ia,ka)=floc(ia,ka)-fac
              fac=p(ka,ka)*epppp(i,j,k,l)
              floc(ia,ja)=floc(ia,ja)+fac
              floc(ja,ia)=floc(ja,ia)+fac
              !if (bla) write (*,*) "done (ab|cc) for: ",sh1,sh2,sh3,sh4,"i,j,k,l",i,j,k,l
            else if (sameij.and.(i.eq.j)) then  !  (aa|bc)
              floc(ia,ia)=floc(ia,ia)+2._dp*p(ka,la)*epppp(i,j,k,l)
              fac=.5_dp*p(ia,ka)*epppp(i,j,k,l)
              floc(ia,la)=floc(ia,la)-fac
              floc(la,ia)=floc(la,ia)-fac
              fac=.5_dp*p(ia,la)*epppp(i,j,k,l)
              floc(ia,ka)=floc(ia,ka)-fac
              floc(ka,ia)=floc(ka,ia)-fac
              fac=p(ia,ia)*epppp(i,j,k,l)
              floc(ka,la)=floc(ka,la)+fac
              floc(la,ka)=floc(la,ka)+fac
              !if (bla) write (*,*) "done (aa|bc) for: ",sh1,sh2,sh3,sh4,"i,j,k,l",i,j,k,l
            else if (sameij.and.(i.lt.j)) then    !  (ba|cd)
              !if (bla) write (*,*) "done nothing (ba|cd) for: ",sh1,sh2,sh3,sh4,"i,j,k,l",i,j,k,l
              cycle
            else if (samekl.and.(k.eq.l)) then                  ! (ab|cc)
              floc(ka,ka)=floc(ka,ka)+2._dp*p(ia,ja)*epppp(i,j,k,l)
              fac=.5_dp*p(ka,ia)*epppp(i,j,k,l)
              floc(ka,ja)=floc(ka,ja)-fac
              floc(ja,ka)=floc(ja,ka)-fac
              fac=.5_dp*p(ka,ja)*epppp(i,j,k,l)
              floc(ka,ia)=floc(ka,ia)-fac
              floc(ia,ka)=floc(ia,ka)-fac
              fac=p(ka,ka)*epppp(i,j,k,l)
              floc(ia,ja)=floc(ia,ja)+fac
              floc(ja,ia)=floc(ja,ia)+fac
            else if (samekl.and.(k.lt.l)) then    !  (ba|cd)
              cycle
            else  !   (ab|cd) or maybe (ab|ac)
              fac=2._dp*  p(la,ka) * epppp(i,j,k,l)
              floc(ia,ja) = floc(ia,ja) +  fac
              floc(ja,ia) = floc(ja,ia) +  fac
              fac=.5_dp*p(ja,ka)*epppp(i,j,k,l)
              floc(ia,la) = floc(ia,la) - fac
              floc(la,ia) = floc(la,ia) - fac
              fac=.5_dp*p(ja,la)*epppp(i,j,k,l)
              floc(ia,ka) = floc(ia,ka) - fac
              floc(ka,ia) = floc(ka,ia) - fac
              fac=.5_dp*p(ia,ka)*epppp(i,j,k,l)
              floc(ja,la) = floc(ja,la) - fac
              floc(la,ja) = floc(la,ja) - fac
              fac=.5_dp*p(ia,la)*epppp(i,j,k,l)
              floc(ja,ka) = floc(ja,ka) - fac
              floc(ka,ja) = floc(ka,ja) - fac
              fac=2._dp*p(ja,ia)*epppp(i,j,k,l)
              floc(ka,la) = floc(ka,la) +  fac
              floc(la,ka) = floc(la,ka) +  fac
            end if
          end do
        end do
      end do
    end do
  end do
end do
!$OMP END DO
!$OMP MASTER
if (prtlevl.ge.4) call time_diff_omp(1,"(pp|pp) integ:")
!$OMP END MASTER

!$OMP DO SCHEDULE(DYNAMIC)
do ii=1,npp  ! (pp|pd)
  sh1=ppprs(1,orderpp(ii))
  sh2=ppprs(2,orderpp(ii))
  ab1=ind2(sh1,sh2)
  do jj=1,npd
    sh3=pdprs(1,orderpd(jj))
    sh4=pdprs(2,orderpd(jj))
    ab2=ind2(sh3,sh4)
    maxdmcoul=max(maxdm(ind2(sh1,sh2)),maxdm(ind2(sh3,sh4)))
    maxdmexch=max(maxdm(ind2(sh1,sh3)),maxdm(ind2(sh1,sh4)),maxdm(ind2(sh2,sh3)),&
       & maxdm(ind2(sh2,sh4)))
    maxdmval=max(maxdmcoul,.25_dp*maxdmexch)
    if (schwpp(orderpp(ii))*schwpd(orderpd(jj))*maxdmval.lt.thr) then
       nprsscr=nprsscr+1
       cycle
    end if
    call i2e_pppd(sh1,sh2,sh3,sh4,epppd)
    if (maxval(abs(epppd)).lt.thr) cycle
    ! The special case here is if shells sh1 and sh2 are the same.
    sameij=(sh1.eq.sh2)
    do i=1,3
      ia=aosh(sh1)+i-1
      do j=1,3
        ja=aosh(sh2)+j-1
        if (sameij.and.(i.eq.j)) then  ! when k=l these are (ab|cc)
          do k=1,3
            ka=aosh(sh3)+k-1
            do l=1,5
              la=aosh(sh4)+l-1
              floc(ia,ia)=floc(ia,ia)+2._dp*p(ka,la)*epppd(i,j,k,l)
              fac=.5_dp*p(ia,ka)*epppd(i,j,k,l)
              floc(ia,la)=floc(ia,la)-fac
              floc(la,ia)=floc(la,ia)-fac
              fac=.5_dp*p(ia,la)*epppd(i,j,k,l)
              floc(ia,ka)=floc(ia,ka)-fac
              floc(ka,ia)=floc(ka,ia)-fac
              fac=p(ia,ia)*epppd(i,j,k,l)
              floc(ka,la)=floc(ka,la)+fac
              floc(la,ka)=floc(la,ka)+fac
            end do
          end do
        else if (sameij.and.(i.lt.j)) then
          cycle
        else
          do k=1,3
            ka=aosh(sh3)+k-1
            do l=1,5
              la=aosh(sh4)+l-1
              fac=2._dp*  p(la,ka) * epppd(i,j,k,l)
              floc(ia,ja) = floc(ia,ja) +  fac
              floc(ja,ia) = floc(ja,ia) +  fac
              fac=.5_dp*p(ja,ka)*epppd(i,j,k,l)
              floc(ia,la) = floc(ia,la) - fac
              floc(la,ia) = floc(la,ia) - fac
              fac=.5_dp*p(ja,la)*epppd(i,j,k,l)
              floc(ia,ka) = floc(ia,ka) - fac
              floc(ka,ia) = floc(ka,ia) - fac
              fac=.5_dp*p(ia,ka)*epppd(i,j,k,l)
              floc(ja,la) = floc(ja,la) - fac
              floc(la,ja) = floc(la,ja) - fac
              fac=.5_dp*p(ia,la)*epppd(i,j,k,l)
              floc(ja,ka) = floc(ja,ka) - fac
              floc(ka,ja) = floc(ka,ja) - fac
              fac=2._dp*p(ja,ia)*epppd(i,j,k,l)
              floc(ka,la) = floc(ka,la) +  fac
              floc(la,ka) = floc(la,ka) +  fac
            end do
          end do
        end if
      end do
    end do
  end do
end do
!$OMP END DO
!$OMP MASTER
if (prtlevl.ge.4) call time_diff_omp(1,"(pp|pd) integ:")
!$OMP END MASTER

!$OMP DO SCHEDULE(DYNAMIC)
do ii=1,npp  ! (pp|pl)
  sh1=ppprs(1,orderpp(ii))
  sh2=ppprs(2,orderpp(ii))
  ab1=ind2(sh1,sh2)
  do jj=1,npl
    sh3=plprs(1,orderpl(jj))
    sh4=plprs(2,orderpl(jj))
    ab2=ind2(sh3,sh4)
    maxdmcoul=max(maxdm(ind2(sh1,sh2)),maxdm(ind2(sh3,sh4)))
    maxdmexch=max(maxdm(ind2(sh1,sh3)),maxdm(ind2(sh1,sh4)),maxdm(ind2(sh2,sh3)),&
       & maxdm(ind2(sh2,sh4)))
    maxdmval=max(maxdmcoul,.25_dp*maxdmexch)
    if (schwpp(orderpp(ii))*schwpl(orderpl(jj))*maxdmval.lt.thr) then
       nprsscr=nprsscr+1
       cycle
    end if
    call i2e_pppl(sh1,sh2,sh3,sh4,epppl)
    if (maxval(abs(epppl)).lt.thr) cycle
    ! The special case here is if shells sh1 and sh2 are the same.
    sameij=(sh1.eq.sh2)
    do i=1,3
      ia=aosh(sh1)+i-1
      do j=1,3
        ja=aosh(sh2)+j-1
        if (sameij.and.(i.eq.j)) then  ! when k=l these are (ab|cc)
          do k=1,3
            ka=aosh(sh3)+k-1
            do l=1,4
              la=aosh(sh4)+l-1
              floc(ia,ia)=floc(ia,ia)+2._dp*p(ka,la)*epppl(i,j,k,l)
              fac=.5_dp*p(ia,ka)*epppl(i,j,k,l)
              floc(ia,la)=floc(ia,la)-fac
              floc(la,ia)=floc(la,ia)-fac
              fac=.5_dp*p(ia,la)*epppl(i,j,k,l)
              floc(ia,ka)=floc(ia,ka)-fac
              floc(ka,ia)=floc(ka,ia)-fac
              fac=p(ia,ia)*epppl(i,j,k,l)
              floc(ka,la)=floc(ka,la)+fac
              floc(la,ka)=floc(la,ka)+fac
            end do
          end do
        else if (sameij.and.(i.lt.j)) then
          cycle
        else
          do k=1,3
            ka=aosh(sh3)+k-1
            do l=1,4
              la=aosh(sh4)+l-1
              fac=2._dp*  p(la,ka) * epppl(i,j,k,l)
              floc(ia,ja) = floc(ia,ja) +  fac
              floc(ja,ia) = floc(ja,ia) +  fac
              fac=.5_dp*p(ja,ka)*epppl(i,j,k,l)
              floc(ia,la) = floc(ia,la) - fac
              floc(la,ia) = floc(la,ia) - fac
              fac=.5_dp*p(ja,la)*epppl(i,j,k,l)
              floc(ia,ka) = floc(ia,ka) - fac
              floc(ka,ia) = floc(ka,ia) - fac
              fac=.5_dp*p(ia,ka)*epppl(i,j,k,l)
              floc(ja,la) = floc(ja,la) - fac
              floc(la,ja) = floc(la,ja) - fac
              fac=.5_dp*p(ia,la)*epppl(i,j,k,l)
              floc(ja,ka) = floc(ja,ka) - fac
              floc(ka,ja) = floc(ka,ja) - fac
              fac=2._dp*p(ja,ia)*epppl(i,j,k,l)
              floc(ka,la) = floc(ka,la) +  fac
              floc(la,ka) = floc(la,ka) +  fac
            end do
          end do
        end if
      end do
    end do
  end do
end do
!$OMP END DO
!$OMP MASTER
if (prtlevl.ge.4) call time_diff_omp(1,"(pp|pl) integ:")
!$OMP END MASTER

!$OMP DO SCHEDULE(DYNAMIC)
do ii=1,npp  ! (pp|dd)
  sh1=ppprs(1,orderpp(ii))
  sh2=ppprs(2,orderpp(ii))
  ab1=ind2(sh1,sh2)
  sameij=(sh1.eq.sh2)
  do jj=1,ndd
    sh3=ddprs(1,orderdd(jj))
    sh4=ddprs(2,orderdd(jj))
    ab2=ind2(sh3,sh4)
    maxdmcoul=max(maxdm(ind2(sh1,sh2)),maxdm(ind2(sh3,sh4)))
    maxdmexch=max(maxdm(ind2(sh1,sh3)),maxdm(ind2(sh1,sh4)),maxdm(ind2(sh2,sh3)),&
       & maxdm(ind2(sh2,sh4)))
    maxdmval=max(maxdmcoul,.25_dp*maxdmexch)
    if (schwpp(orderpp(ii))*schwdd(orderdd(jj))*maxdmval.lt.thr) then
        nprsscr=nprsscr+1
        cycle
    end if
    call i2e_ppdd(sh1,sh2,sh3,sh4,eppdd)
    if (maxval(abs(eppdd)).lt.thr) cycle
    ! The special cases here are ij i=j, and/or k=l.
    samekl=(sh3.eq.sh4)
    do i=1,3
      ia=aosh(sh1)+i-1
      do j=1,3
        ja=aosh(sh2)+j-1
        do k=1,5
          ka=aosh(sh3)+k-1
          do l=1,5
            la=aosh(sh4)+l-1
            if (sameij.and.(i.eq.j).and.samekl.and.(k.eq.l)) then  ! (aa|bb)
              floc(ia,ia)=floc(ia,ia)+p(ka,ka)*eppdd(i,i,k,k)
              fac=.5_dp*p(ia,ka)*eppdd(i,i,k,k)
              floc(ia,ka)=floc(ia,ka)-fac
              floc(ka,ia)=floc(ka,ia)-fac
              floc(ka,ka)=floc(ka,ka)+p(ia,ia)*eppdd(i,i,k,k)
            else if (sameij.and.(i.eq.j).and.samekl.and.(k.lt.l)) then  ! (aa|cb)
              cycle
            else if (sameij.and.(i.eq.j).and.samekl.and.(k.gt.l)) then  !  (aa|bc)
              floc(ia,ia)=floc(ia,ia)+2._dp*p(ka,la)*eppdd(i,j,k,l)
              fac=.5_dp*p(ia,ka)*eppdd(i,j,k,l)
              floc(ia,la)=floc(ia,la)-fac
              floc(la,ia)=floc(la,ia)-fac
              fac=.5_dp*p(ia,la)*eppdd(i,j,k,l)
              floc(ia,ka)=floc(ia,ka)-fac
              floc(ka,ia)=floc(ka,ia)-fac
              fac=p(ia,ia)*eppdd(i,j,k,l)
              floc(ka,la)=floc(ka,la)+fac
              floc(la,ka)=floc(la,ka)+fac
            else if (sameij.and.(i.lt.j).and.samekl.and.(k.eq.l)) then    !  (ba|cc)
              cycle
            else if (samekl.and.(k.eq.l).and.sameij.and.(i.gt.j)) then    ! (ab|cc)
              floc(ka,ka)=floc(ka,ka)+2._dp*p(ia,ja)*eppdd(i,j,k,l)
              fac=.5_dp*p(ka,ia)*eppdd(i,j,k,l)
              floc(ka,ja)=floc(ka,ja)-fac
              floc(ja,ka)=floc(ja,ka)-fac
              fac=.5_dp*p(ka,ja)*eppdd(i,j,k,l)
              floc(ka,ia)=floc(ka,ia)-fac
              floc(ia,ka)=floc(ia,ka)-fac
              fac=p(ka,ka)*eppdd(i,j,k,l)
              floc(ia,ja)=floc(ia,ja)+fac
              floc(ja,ia)=floc(ja,ia)+fac
            else if (sameij.and.(i.eq.j)) then  !  (aa|bc)
              floc(ia,ia)=floc(ia,ia)+2._dp*p(ka,la)*eppdd(i,j,k,l)
              fac=.5_dp*p(ia,ka)*eppdd(i,j,k,l)
              floc(ia,la)=floc(ia,la)-fac
              floc(la,ia)=floc(la,ia)-fac
              fac=.5_dp*p(ia,la)*eppdd(i,j,k,l)
              floc(ia,ka)=floc(ia,ka)-fac
              floc(ka,ia)=floc(ka,ia)-fac
              fac=p(ia,ia)*eppdd(i,j,k,l)
              floc(ka,la)=floc(ka,la)+fac
              floc(la,ka)=floc(la,ka)+fac
            else if (sameij.and.(i.lt.j)) then    !  (ba|cd)
              cycle
            else if (samekl.and.(k.eq.l)) then                  ! (ab|cc)
              floc(ka,ka)=floc(ka,ka)+2._dp*p(ia,ja)*eppdd(i,j,k,l)
              fac=.5_dp*p(ka,ia)*eppdd(i,j,k,l)
              floc(ka,ja)=floc(ka,ja)-fac
              floc(ja,ka)=floc(ja,ka)-fac
              fac=.5_dp*p(ka,ja)*eppdd(i,j,k,l)
              floc(ka,ia)=floc(ka,ia)-fac
              floc(ia,ka)=floc(ia,ka)-fac
              fac=p(ka,ka)*eppdd(i,j,k,l)
              floc(ia,ja)=floc(ia,ja)+fac
              floc(ja,ia)=floc(ja,ia)+fac
            else if (samekl.and.(k.lt.l)) then    !  (ba|cd)
              cycle
            else  !   (ab|cd) or maybe (ab|ac)
              fac=2._dp*  p(la,ka) * eppdd(i,j,k,l)
              floc(ia,ja) = floc(ia,ja) +  fac
              floc(ja,ia) = floc(ja,ia) +  fac
              fac=.5_dp*p(ja,ka)*eppdd(i,j,k,l)
              floc(ia,la) = floc(ia,la) - fac
              floc(la,ia) = floc(la,ia) - fac
              fac=.5_dp*p(ja,la)*eppdd(i,j,k,l)
              floc(ia,ka) = floc(ia,ka) - fac
              floc(ka,ia) = floc(ka,ia) - fac
              fac=.5_dp*p(ia,ka)*eppdd(i,j,k,l)
              floc(ja,la) = floc(ja,la) - fac
              floc(la,ja) = floc(la,ja) - fac
              fac=.5_dp*p(ia,la)*eppdd(i,j,k,l)
              floc(ja,ka) = floc(ja,ka) - fac
              floc(ka,ja) = floc(ka,ja) - fac
              fac=2._dp*p(ja,ia)*eppdd(i,j,k,l)
              floc(ka,la) = floc(ka,la) +  fac
              floc(la,ka) = floc(la,ka) +  fac
            end if
          end do
        end do
      end do
    end do
  end do
end do
!$OMP END DO
!$OMP MASTER
if (prtlevl.ge.4) call time_diff_omp(1,"(pp|dd) integ:")
!$OMP END MASTER

! Skipping (pp|dl)


!$OMP DO SCHEDULE(DYNAMIC)
do ii=1,npp  ! (pp|ll)
  sh1=ppprs(1,orderpp(ii))
  sh2=ppprs(2,orderpp(ii))
  ab1=ind2(sh1,sh2)
  sameij=(sh1.eq.sh2)
  do jj=1,nll 
    sh3=llprs(1,orderll(jj))
    sh4=llprs(2,orderll(jj))
    ab2=ind2(sh3,sh4)
    maxdmcoul=max(maxdm(ind2(sh1,sh2)),maxdm(ind2(sh3,sh4)))
    maxdmexch=max(maxdm(ind2(sh1,sh3)),maxdm(ind2(sh1,sh4)),maxdm(ind2(sh2,sh3)),&
       & maxdm(ind2(sh2,sh4)))
    maxdmval=max(maxdmcoul,.25_dp*maxdmexch)
    if (schwpp(orderpp(ii))*schwll(orderll(jj))*maxdmval.lt.thr) then
        nprsscr=nprsscr+1
        cycle
    end if
    call i2e_ppll(sh1,sh2,sh3,sh4,eppll)
    if (maxval(abs(eppll)).lt.thr) cycle
    ! The special cases here are ij i=j, and/or k=l.
    samekl=(sh3.eq.sh4)
    do i=1,3
      ia=aosh(sh1)+i-1
      do j=1,3
        ja=aosh(sh2)+j-1
        do k=1,4
          ka=aosh(sh3)+k-1
          do l=1,4
            la=aosh(sh4)+l-1
            if (sameij.and.(i.eq.j).and.samekl.and.(k.eq.l)) then  ! (aa|bb)
              floc(ia,ia)=floc(ia,ia)+p(ka,ka)*eppll(i,i,k,k)
              fac=.5_dp*p(ia,ka)*eppll(i,i,k,k)
              floc(ia,ka)=floc(ia,ka)-fac
              floc(ka,ia)=floc(ka,ia)-fac
              floc(ka,ka)=floc(ka,ka)+p(ia,ia)*eppll(i,i,k,k)
            else if (sameij.and.(i.eq.j).and.samekl.and.(k.lt.l)) then  ! (aa|cb)
              cycle
            else if (sameij.and.(i.eq.j).and.samekl.and.(k.gt.l)) then  !  (aa|bc)
              floc(ia,ia)=floc(ia,ia)+2._dp*p(ka,la)*eppll(i,j,k,l)
              fac=.5_dp*p(ia,ka)*eppll(i,j,k,l)
              floc(ia,la)=floc(ia,la)-fac
              floc(la,ia)=floc(la,ia)-fac
              fac=.5_dp*p(ia,la)*eppll(i,j,k,l)
              floc(ia,ka)=floc(ia,ka)-fac
              floc(ka,ia)=floc(ka,ia)-fac
              fac=p(ia,ia)*eppll(i,j,k,l)
              floc(ka,la)=floc(ka,la)+fac
              floc(la,ka)=floc(la,ka)+fac
            else if (sameij.and.(i.lt.j).and.samekl.and.(k.eq.l)) then    !  (ba|cc)
              cycle
            else if (samekl.and.(k.eq.l).and.sameij.and.(i.gt.j)) then    ! (ab|cc)
              floc(ka,ka)=floc(ka,ka)+2._dp*p(ia,ja)*eppll(i,j,k,l)
              fac=.5_dp*p(ka,ia)*eppll(i,j,k,l)
              floc(ka,ja)=floc(ka,ja)-fac
              floc(ja,ka)=floc(ja,ka)-fac
              fac=.5_dp*p(ka,ja)*eppll(i,j,k,l)
              floc(ka,ia)=floc(ka,ia)-fac
              floc(ia,ka)=floc(ia,ka)-fac
              fac=p(ka,ka)*eppll(i,j,k,l)
              floc(ia,ja)=floc(ia,ja)+fac
              floc(ja,ia)=floc(ja,ia)+fac
            else if (sameij.and.(i.eq.j)) then  !  (aa|bc)
              floc(ia,ia)=floc(ia,ia)+2._dp*p(ka,la)*eppll(i,j,k,l)
              fac=.5_dp*p(ia,ka)*eppll(i,j,k,l)
              floc(ia,la)=floc(ia,la)-fac
              floc(la,ia)=floc(la,ia)-fac
              fac=.5_dp*p(ia,la)*eppll(i,j,k,l)
              floc(ia,ka)=floc(ia,ka)-fac
              floc(ka,ia)=floc(ka,ia)-fac
              fac=p(ia,ia)*eppll(i,j,k,l)
              floc(ka,la)=floc(ka,la)+fac
              floc(la,ka)=floc(la,ka)+fac
            else if (sameij.and.(i.lt.j)) then    !  (ba|cd)
              cycle
            else if (samekl.and.(k.eq.l)) then                  ! (ab|cc)
              floc(ka,ka)=floc(ka,ka)+2._dp*p(ia,ja)*eppll(i,j,k,l)
              fac=.5_dp*p(ka,ia)*eppll(i,j,k,l)
              floc(ka,ja)=floc(ka,ja)-fac
              floc(ja,ka)=floc(ja,ka)-fac
              fac=.5_dp*p(ka,ja)*eppll(i,j,k,l)
              floc(ka,ia)=floc(ka,ia)-fac
              floc(ia,ka)=floc(ia,ka)-fac
              fac=p(ka,ka)*eppll(i,j,k,l)
              floc(ia,ja)=floc(ia,ja)+fac
              floc(ja,ia)=floc(ja,ia)+fac
            else if (samekl.and.(k.lt.l)) then    !  (ba|cd)
              cycle
            else  !   (ab|cd) or maybe (ab|ac)
              fac=2._dp*  p(la,ka) * eppll(i,j,k,l)
              floc(ia,ja) = floc(ia,ja) +  fac
              floc(ja,ia) = floc(ja,ia) +  fac
              fac=.5_dp*p(ja,ka)*eppll(i,j,k,l)
              floc(ia,la) = floc(ia,la) - fac
              floc(la,ia) = floc(la,ia) - fac
              fac=.5_dp*p(ja,la)*eppll(i,j,k,l)
              floc(ia,ka) = floc(ia,ka) - fac
              floc(ka,ia) = floc(ka,ia) - fac
              fac=.5_dp*p(ia,ka)*eppll(i,j,k,l)
              floc(ja,la) = floc(ja,la) - fac
              floc(la,ja) = floc(la,ja) - fac
              fac=.5_dp*p(ia,la)*eppll(i,j,k,l)
              floc(ja,ka) = floc(ja,ka) - fac
              floc(ka,ja) = floc(ka,ja) - fac
              fac=2._dp*p(ja,ia)*eppll(i,j,k,l)
              floc(ka,la) = floc(ka,la) +  fac
              floc(la,ka) = floc(la,ka) +  fac
            end if
          end do
        end do
      end do
    end do
  end do
end do
!$OMP END DO
!$OMP MASTER
if (prtlevl.ge.4) call time_diff_omp(1,"(pp|ll) integ:")
!$OMP END MASTER

!$OMP DO SCHEDULE(DYNAMIC)
do ii=2,npd  !  (pd|pd)
  sh1=pdprs(1,orderpd(ii))
  sh2=pdprs(2,orderpd(ii))
  ab1=ind2(sh1,sh2)
  ia=aosh(sh1)
  do jj=1,ii-1
    sh3=pdprs(1,orderpd(jj))
    sh4=pdprs(2,orderpd(jj))
    ab2=ind2(sh3,sh4)
    ka=aosh(sh3)
    maxdmcoul=max(maxdm(ind2(sh1,sh2)),maxdm(ind2(sh3,sh4)))
    maxdmexch=max(maxdm(ind2(sh1,sh3)),maxdm(ind2(sh1,sh4)),maxdm(ind2(sh2,sh3)),&
       & maxdm(ind2(sh2,sh4)))
    maxdmval=max(maxdmcoul,.25_dp*maxdmexch)
    if (schwpd(orderpd(ii))*schwpd(orderpd(jj))*maxdmval.lt.thr) then
        nprsscr=nprsscr+1
        cycle
    end if
    call i2e_pdpd(sh1,sh2,sh3,sh4,epdpd)
    if (maxval(abs(epdpd)).lt.thr) cycle
    ! Since the bra and ket shell pairs cannot both be the same, at most we can have
    !  (ab|ac) here, so no need to do anything special.
    do i=1,3
      ia=aosh(sh1)+i-1
      do j=1,5
        ja=aosh(sh2)+j-1
        do k=1,3
          ka=aosh(sh3)+k-1
          do l=1,5
            la=aosh(sh4)+l-1
            fac=2._dp*  p(la,ka) * epdpd(i,j,k,l)
            floc(ia,ja) = floc(ia,ja) +  fac
            floc(ja,ia) = floc(ja,ia) +  fac
            fac=.5_dp*p(ja,ka)*epdpd(i,j,k,l)
            floc(ia,la) = floc(ia,la) - fac
            floc(la,ia) = floc(la,ia) - fac
            fac=.5_dp*p(ja,la)*epdpd(i,j,k,l)
            floc(ia,ka) = floc(ia,ka) - fac
            floc(ka,ia) = floc(ka,ia) - fac
            fac=.5_dp*p(ia,ka)*epdpd(i,j,k,l)
            floc(ja,la) = floc(ja,la) - fac
            floc(la,ja) = floc(la,ja) - fac
            fac=.5_dp*p(ia,la)*epdpd(i,j,k,l)
            floc(ja,ka) = floc(ja,ka) - fac
            floc(ka,ja) = floc(ka,ja) - fac
            fac=2._dp*p(ja,ia)*epdpd(i,j,k,l)
            floc(ka,la) = floc(ka,la) +  fac
            floc(la,ka) = floc(la,ka) +   fac
          end do
        end do
      end do
    end do
  end do
end do
!$OMP END DO
!$OMP MASTER
if (prtlevl.ge.4) call time_diff_omp(1,"(pd|pd) integ:")
!$OMP END MASTER


! Skipping (pd|pl)

!$OMP DO SCHEDULE(DYNAMIC)
do ii=1,npd  ! (pd|dd)
  sh1=pdprs(1,orderpd(ii))
  sh2=pdprs(2,orderpd(ii))
  ab1=ind2(sh1,sh2)
  ia=aosh(sh1)
  do jj=1,ndd
    sh3=ddprs(1,orderdd(jj))
    sh4=ddprs(2,orderdd(jj))
    ab2=ind2(sh3,sh4)
    maxdmcoul=max(maxdm(ind2(sh1,sh2)),maxdm(ind2(sh3,sh4)))
    maxdmexch=max(maxdm(ind2(sh1,sh3)),maxdm(ind2(sh1,sh4)),maxdm(ind2(sh2,sh3)),&
       & maxdm(ind2(sh2,sh4)))
    maxdmval=max(maxdmcoul,.25_dp*maxdmexch)
    if (schwpd(orderpd(ii))*schwdd(orderdd(jj))*maxdmval.lt.thr) then
       nprsscr=nprsscr+1
       cycle
    end if
    call i2e_pddd(sh1,sh2,sh3,sh4,epddd)
    if (maxval(abs(epddd)).lt.thr) cycle
    ! The only special case here is if shells sh3 and sh4 are the same.
    samekl=(sh3.eq.sh4)
    do i=1,3
      ia=aosh(sh1)+i-1
      do j=1,5
        ja=aosh(sh2)+j-1
        do k=1,5
          ka=aosh(sh3)+k-1
          do l=1,5
            la=aosh(sh4)+l-1
            if (samekl.and.(k.eq.l)) then  ! when k=l these are (ab|cc)
              floc(ka,ka)=floc(ka,ka)+2._dp*p(ia,ja)*epddd(i,j,k,k)
              fac=.5_dp*p(ka,ia)*epddd(i,j,k,k)
              floc(ka,ja)=floc(ka,ja)-fac
              floc(ja,ka)=floc(ja,ka)-fac
              fac=.5_dp*p(ka,ja)*epddd(i,j,k,k)
              floc(ka,ia)=floc(ka,ia)-fac
              floc(ia,ka)=floc(ia,ka)-fac
              fac=p(ka,ka)*epddd(i,j,k,k)
              floc(ia,ja)=floc(ia,ja)+fac
              floc(ja,ia)=floc(ja,ia)+fac
            else if (samekl.and.(k.lt.l)) then
              cycle
            else
              fac=2._dp*  p(la,ka) * epddd(i,j,k,l)
              floc(ia,ja) = floc(ia,ja) +  fac
              floc(ja,ia) = floc(ja,ia) +  fac
              fac=.5_dp*p(ja,ka)*epddd(i,j,k,l)
              floc(ia,la) = floc(ia,la) - fac
              floc(la,ia) = floc(la,ia) - fac
              fac=.5_dp*p(ja,la)*epddd(i,j,k,l)
              floc(ia,ka) = floc(ia,ka) - fac
              floc(ka,ia) = floc(ka,ia) - fac
              fac=.5_dp*p(ia,ka)*epddd(i,j,k,l)
              floc(ja,la) = floc(ja,la) - fac
              floc(la,ja) = floc(la,ja) - fac
              fac=.5_dp*p(ia,la)*epddd(i,j,k,l)
              floc(ja,ka) = floc(ja,ka) - fac
              floc(ka,ja) = floc(ka,ja) - fac
              fac=2._dp*p(ja,ia)*epddd(i,j,k,l)
              floc(ka,la) = floc(ka,la) +  fac
              floc(la,ka) = floc(la,ka) +  fac
            end if
          end do
        end do
      end do
    end do
  end do
end do
!$OMP END DO
!$OMP MASTER
if (prtlevl.ge.4) call time_diff_omp(1,"(pd|dd) integ:")
!$OMP END MASTER


! Skipping (pd|dl) and (pd|ll)

!$OMP DO SCHEDULE(DYNAMIC)
do ii=2,npl  !  (pl|pl)
  sh1=plprs(1,orderpl(ii))
  sh2=plprs(2,orderpl(ii))
  ab1=ind2(sh1,sh2)
  ia=aosh(sh1)
  do jj=1,ii-1
    sh3=plprs(1,orderpl(jj))
    sh4=plprs(2,orderpl(jj))
    ab2=ind2(sh3,sh4)
    ka=aosh(sh3)
    maxdmcoul=max(maxdm(ind2(sh1,sh2)),maxdm(ind2(sh3,sh4)))
    maxdmexch=max(maxdm(ind2(sh1,sh3)),maxdm(ind2(sh1,sh4)),maxdm(ind2(sh2,sh3)),&
       & maxdm(ind2(sh2,sh4)))
    maxdmval=max(maxdmcoul,.25_dp*maxdmexch)
    if (schwpl(orderpl(ii))*schwpl(orderpl(jj))*maxdmval.lt.thr) then
        nprsscr=nprsscr+1
        cycle
    end if
    call i2e_plpl(sh1,sh2,sh3,sh4,eplpl)
    if (maxval(abs(eplpl)).lt.thr) cycle
    ! Since the bra and ket shell pairs cannot both be the same, at most we can have
    !  (ab|ac) here, so no need to do anything special.
    do i=1,3
      ia=aosh(sh1)+i-1
      do j=1,4
        ja=aosh(sh2)+j-1
        do k=1,3
          ka=aosh(sh3)+k-1
          do l=1,4
            la=aosh(sh4)+l-1
            fac=2._dp*  p(la,ka) * eplpl(i,j,k,l)
            floc(ia,ja) = floc(ia,ja) +  fac
            floc(ja,ia) = floc(ja,ia) +  fac
            fac=.5_dp*p(ja,ka)*eplpl(i,j,k,l)
            floc(ia,la) = floc(ia,la) - fac
            floc(la,ia) = floc(la,ia) - fac
            fac=.5_dp*p(ja,la)*eplpl(i,j,k,l)
            floc(ia,ka) = floc(ia,ka) - fac
            floc(ka,ia) = floc(ka,ia) - fac
            fac=.5_dp*p(ia,ka)*eplpl(i,j,k,l)
            floc(ja,la) = floc(ja,la) - fac
            floc(la,ja) = floc(la,ja) - fac
            fac=.5_dp*p(ia,la)*eplpl(i,j,k,l)
            floc(ja,ka) = floc(ja,ka) - fac
            floc(ka,ja) = floc(ka,ja) - fac
            fac=2._dp*p(ja,ia)*eplpl(i,j,k,l)
            floc(ka,la) = floc(ka,la) +  fac
            floc(la,ka) = floc(la,ka) +   fac
          end do
        end do
      end do
    end do
  end do
end do
!$OMP END DO
!$OMP MASTER
if (prtlevl.ge.4) call time_diff_omp(1,"(pl|pl) integ:")
!$OMP END MASTER


! Skipping (pl|dd) and (pl|dl)

!$OMP DO SCHEDULE(DYNAMIC)
do ii=1,npl  ! (pl|ll)
  sh1=plprs(1,orderpl(ii))
  sh2=plprs(2,orderpl(ii))
  ab1=ind2(sh1,sh2)
  ia=aosh(sh1)
  do jj=1,nll
    sh3=llprs(1,orderll(jj))
    sh4=llprs(2,orderll(jj))
    ab2=ind2(sh3,sh4)
    maxdmcoul=max(maxdm(ind2(sh1,sh2)),maxdm(ind2(sh3,sh4)))
    maxdmexch=max(maxdm(ind2(sh1,sh3)),maxdm(ind2(sh1,sh4)),maxdm(ind2(sh2,sh3)),&
       & maxdm(ind2(sh2,sh4)))
    maxdmval=max(maxdmcoul,.25_dp*maxdmexch)
    if (schwpl(orderpl(ii))*schwll(orderll(jj))*maxdmval.lt.thr) then
       nprsscr=nprsscr+1
       cycle
    end if
    call i2e_plll(sh1,sh2,sh3,sh4,eplll)
    if (maxval(abs(eplll)).lt.thr) cycle
    ! The only special case here is if shells sh3 and sh4 are the same.
    samekl=(sh3.eq.sh4)
    do i=1,3
      ia=aosh(sh1)+i-1
      do j=1,4
        ja=aosh(sh2)+j-1
        do k=1,4
          ka=aosh(sh3)+k-1
          do l=1,4
            la=aosh(sh4)+l-1
            if (samekl.and.(k.eq.l)) then  ! when k=l these are (ab|cc)
              floc(ka,ka)=floc(ka,ka)+2._dp*p(ia,ja)*eplll(i,j,k,k)
              fac=.5_dp*p(ka,ia)*eplll(i,j,k,k)
              floc(ka,ja)=floc(ka,ja)-fac
              floc(ja,ka)=floc(ja,ka)-fac
              fac=.5_dp*p(ka,ja)*eplll(i,j,k,k)
              floc(ka,ia)=floc(ka,ia)-fac
              floc(ia,ka)=floc(ia,ka)-fac
              fac=p(ka,ka)*eplll(i,j,k,k)
              floc(ia,ja)=floc(ia,ja)+fac
              floc(ja,ia)=floc(ja,ia)+fac
            else if (samekl.and.(k.lt.l)) then
              cycle
            else
              fac=2._dp*  p(la,ka) * eplll(i,j,k,l)
              floc(ia,ja) = floc(ia,ja) +  fac
              floc(ja,ia) = floc(ja,ia) +  fac
              fac=.5_dp*p(ja,ka)*eplll(i,j,k,l)
              floc(ia,la) = floc(ia,la) - fac
              floc(la,ia) = floc(la,ia) - fac
              fac=.5_dp*p(ja,la)*eplll(i,j,k,l)
              floc(ia,ka) = floc(ia,ka) - fac
              floc(ka,ia) = floc(ka,ia) - fac
              fac=.5_dp*p(ia,ka)*eplll(i,j,k,l)
              floc(ja,la) = floc(ja,la) - fac
              floc(la,ja) = floc(la,ja) - fac
              fac=.5_dp*p(ia,la)*eplll(i,j,k,l)
              floc(ja,ka) = floc(ja,ka) - fac
              floc(ka,ja) = floc(ka,ja) - fac
              fac=2._dp*p(ja,ia)*eplll(i,j,k,l)
              floc(ka,la) = floc(ka,la) +  fac
              floc(la,ka) = floc(la,ka) +  fac
            end if
          end do
        end do
      end do
    end do
  end do
end do
!$OMP END DO
!$OMP MASTER
if (prtlevl.ge.4) call time_diff_omp(1,"(pl|ll) integ:")
!$OMP END MASTER

!$OMP DO SCHEDULE(DYNAMIC)
do ii=2,ndd  ! (dd|dd)
  sh1=ddprs(1,orderdd(ii))
  sh2=ddprs(2,orderdd(ii))
  ab1=ind2(sh1,sh2)
  sameij=(sh1.eq.sh2)
  do jj=1,ii-1
    sh3=ddprs(1,orderdd(jj))
    sh4=ddprs(2,orderdd(jj))
    ab2=ind2(sh3,sh4)
    maxdmcoul=max(maxdm(ind2(sh1,sh2)),maxdm(ind2(sh3,sh4)))
    maxdmexch=max(maxdm(ind2(sh1,sh3)),maxdm(ind2(sh1,sh4)),maxdm(ind2(sh2,sh3)),&
       & maxdm(ind2(sh2,sh4)))
    maxdmval=max(maxdmcoul,.25_dp*maxdmexch)
    if (schwdd(orderdd(ii))*schwdd(orderdd(jj))*maxdmval.lt.thr) then
        nprsscr=nprsscr+1
        cycle
    end if
    call i2e_dddd(sh1,sh2,sh3,sh4,edddd)
    if (maxval(abs(edddd)).lt.thr) cycle
    ! The special cases here are ij i=j, and/or k=l. It cannot be the case that i=j=k=l or that i=k and j=l, so (aa|aa) and (ab|ab) are ruled out
    samekl=(sh3.eq.sh4)
    do i=1,5
      ia=aosh(sh1)+i-1
      do j=1,5
        ja=aosh(sh2)+j-1
        do k=1,5
          ka=aosh(sh3)+k-1
          do l=1,5
            la=aosh(sh4)+l-1
            if (sameij.and.(i.eq.j).and.samekl.and.(k.eq.l)) then  ! (aa|bb)
              floc(ia,ia)=floc(ia,ia)+p(ka,ka)*edddd(i,i,k,k)
              fac=.5_dp*p(ia,ka)*edddd(i,i,k,k)
              floc(ia,ka)=floc(ia,ka)-fac
              floc(ka,ia)=floc(ka,ia)-fac
              floc(ka,ka)=floc(ka,ka)+p(ia,ia)*edddd(i,i,k,k)
            else if (sameij.and.(i.eq.j).and.samekl.and.(k.lt.l)) then  ! (aa|cb)
              cycle
            else if (sameij.and.(i.eq.j).and.samekl.and.(k.gt.l)) then  !  (aa|bc)
              floc(ia,ia)=floc(ia,ia)+2._dp*p(ka,la)*edddd(i,j,k,l)
              fac=.5_dp*p(ia,ka)*edddd(i,j,k,l)
              floc(ia,la)=floc(ia,la)-fac
              floc(la,ia)=floc(la,ia)-fac
              fac=.5_dp*p(ia,la)*edddd(i,j,k,l)
              floc(ia,ka)=floc(ia,ka)-fac
              floc(ka,ia)=floc(ka,ia)-fac
              fac=p(ia,ia)*edddd(i,j,k,l)
              floc(ka,la)=floc(ka,la)+fac
              floc(la,ka)=floc(la,ka)+fac
            else if (sameij.and.(i.lt.j).and.samekl.and.(k.eq.l)) then    !  (ba|cc)
              cycle
            else if (samekl.and.(k.eq.l).and.sameij.and.(i.gt.j)) then    ! (ab|cc)
              floc(ka,ka)=floc(ka,ka)+2._dp*p(ia,ja)*edddd(i,j,k,l)
              fac=.5_dp*p(ka,ia)*edddd(i,j,k,l)
              floc(ka,ja)=floc(ka,ja)-fac
              floc(ja,ka)=floc(ja,ka)-fac
              fac=.5_dp*p(ka,ja)*edddd(i,j,k,l)
              floc(ka,ia)=floc(ka,ia)-fac
              floc(ia,ka)=floc(ia,ka)-fac
              fac=p(ka,ka)*edddd(i,j,k,l)
              floc(ia,ja)=floc(ia,ja)+fac
              floc(ja,ia)=floc(ja,ia)+fac
            else if (sameij.and.(i.eq.j)) then  !  (aa|bc)
              floc(ia,ia)=floc(ia,ia)+2._dp*p(ka,la)*edddd(i,j,k,l)
              fac=.5_dp*p(ia,ka)*edddd(i,j,k,l)
              floc(ia,la)=floc(ia,la)-fac
              floc(la,ia)=floc(la,ia)-fac
              fac=.5_dp*p(ia,la)*edddd(i,j,k,l)
              floc(ia,ka)=floc(ia,ka)-fac
              floc(ka,ia)=floc(ka,ia)-fac
              fac=p(ia,ia)*edddd(i,j,k,l)
              floc(ka,la)=floc(ka,la)+fac
              floc(la,ka)=floc(la,ka)+fac
            else if (sameij.and.(i.lt.j)) then    !  (ba|cd)
              cycle
            else if (samekl.and.(k.eq.l)) then                  ! (ab|cc)
              floc(ka,ka)=floc(ka,ka)+2._dp*p(ia,ja)*edddd(i,j,k,l)
              fac=.5_dp*p(ka,ia)*edddd(i,j,k,l)
              floc(ka,ja)=floc(ka,ja)-fac
              floc(ja,ka)=floc(ja,ka)-fac
              fac=.5_dp*p(ka,ja)*edddd(i,j,k,l)
              floc(ka,ia)=floc(ka,ia)-fac
              floc(ia,ka)=floc(ia,ka)-fac
              fac=p(ka,ka)*edddd(i,j,k,l)
              floc(ia,ja)=floc(ia,ja)+fac
              floc(ja,ia)=floc(ja,ia)+fac
            else if (samekl.and.(k.lt.l)) then    !  (ba|cd)
              cycle
            else  !   (ab|cd) or maybe (ab|ac)
              fac=2._dp*  p(la,ka) * edddd(i,j,k,l)
              floc(ia,ja) = floc(ia,ja) +  fac
              floc(ja,ia) = floc(ja,ia) +  fac
              fac=.5_dp*p(ja,ka)*edddd(i,j,k,l)
              floc(ia,la) = floc(ia,la) - fac
              floc(la,ia) = floc(la,ia) - fac
              fac=.5_dp*p(ja,la)*edddd(i,j,k,l)
              floc(ia,ka) = floc(ia,ka) - fac
              floc(ka,ia) = floc(ka,ia) - fac
              fac=.5_dp*p(ia,ka)*edddd(i,j,k,l)
              floc(ja,la) = floc(ja,la) - fac
              floc(la,ja) = floc(la,ja) - fac
              fac=.5_dp*p(ia,la)*edddd(i,j,k,l)
              floc(ja,ka) = floc(ja,ka) - fac
              floc(ka,ja) = floc(ka,ja) - fac
              fac=2._dp*p(ja,ia)*edddd(i,j,k,l)
              floc(ka,la) = floc(ka,la) +  fac
              floc(la,ka) = floc(la,ka) +  fac
            end if
          end do
        end do
      end do
    end do
  end do
end do
!$OMP END DO
!$OMP MASTER
if (prtlevl.ge.4) call time_diff_omp(1,"(dd|dd) integ:")
!$OMP END MASTER

!$OMP DO SCHEDULE(DYNAMIC)
do ii=1,ndd  ! (dd|dl)
  sh1=ddprs(1,orderdd(ii))
  sh2=ddprs(2,orderdd(ii))
  ab1=ind2(sh1,sh2)
  do jj=1,ndl
    sh3=dlprs(1,orderdl(jj))
    sh4=dlprs(2,orderdl(jj))
    ab2=ind2(sh3,sh4)
    maxdmcoul=max(maxdm(ind2(sh1,sh2)),maxdm(ind2(sh3,sh4)))
    maxdmexch=max(maxdm(ind2(sh1,sh3)),maxdm(ind2(sh1,sh4)),maxdm(ind2(sh2,sh3)),&
       & maxdm(ind2(sh2,sh4)))
    maxdmval=max(maxdmcoul,.25_dp*maxdmexch)
    if (schwdd(orderdd(ii))*schwdl(orderdl(jj))*maxdmval.lt.thr) then
       nprsscr=nprsscr+1
       cycle
    end if
    call i2e_dddl(sh1,sh2,sh3,sh4,edddl)
    if (maxval(abs(edddl)).lt.thr) cycle
    ! The special case here is if shells sh1 and sh2 are the same.
    sameij=(sh1.eq.sh2)
    do i=1,5
      ia=aosh(sh1)+i-1
      do j=1,5
        ja=aosh(sh2)+j-1
        if (sameij.and.(i.eq.j)) then  ! when k=l these are (ab|cc)
          do k=1,5
            ka=aosh(sh3)+k-1
            do l=1,4
              la=aosh(sh4)+l-1
              floc(ia,ia)=floc(ia,ia)+2._dp*p(ka,la)*edddl(i,j,k,l)
              fac=.5_dp*p(ia,ka)*edddl(i,j,k,l)
              floc(ia,la)=floc(ia,la)-fac
              floc(la,ia)=floc(la,ia)-fac
              fac=.5_dp*p(ia,la)*edddl(i,j,k,l)
              floc(ia,ka)=floc(ia,ka)-fac
              floc(ka,ia)=floc(ka,ia)-fac
              fac=p(ia,ia)*edddl(i,j,k,l)
              floc(ka,la)=floc(ka,la)+fac
              floc(la,ka)=floc(la,ka)+fac
            end do
          end do
        else if (sameij.and.(i.lt.j)) then
          cycle
        else
          do k=1,5
            ka=aosh(sh3)+k-1
            do l=1,4
              la=aosh(sh4)+l-1
              fac=2._dp*  p(la,ka) * edddl(i,j,k,l)
              floc(ia,ja) = floc(ia,ja) +  fac
              floc(ja,ia) = floc(ja,ia) +  fac
              fac=.5_dp*p(ja,ka)*edddl(i,j,k,l)
              floc(ia,la) = floc(ia,la) - fac
              floc(la,ia) = floc(la,ia) - fac
              fac=.5_dp*p(ja,la)*edddl(i,j,k,l)
              floc(ia,ka) = floc(ia,ka) - fac
              floc(ka,ia) = floc(ka,ia) - fac
              fac=.5_dp*p(ia,ka)*edddl(i,j,k,l)
              floc(ja,la) = floc(ja,la) - fac
              floc(la,ja) = floc(la,ja) - fac
              fac=.5_dp*p(ia,la)*edddl(i,j,k,l)
              floc(ja,ka) = floc(ja,ka) - fac
              floc(ka,ja) = floc(ka,ja) - fac
              fac=2._dp*p(ja,ia)*edddl(i,j,k,l)
              floc(ka,la) = floc(ka,la) +  fac
              floc(la,ka) = floc(la,ka) +  fac
            end do
          end do
        end if
      end do
    end do
  end do
end do
!$OMP END DO
!$OMP MASTER
if (prtlevl.ge.4) call time_diff_omp(1,"(dd|dl) integ:")
!$OMP END MASTER

!$OMP DO SCHEDULE(DYNAMIC)
do ii=1,ndd  ! (dd|ll)
  sh1=ddprs(1,orderdd(ii))
  sh2=ddprs(2,orderdd(ii))
  ab1=ind2(sh1,sh2)
  sameij=(sh1.eq.sh2)
  do jj=1,nll 
    sh3=llprs(1,orderll(jj))
    sh4=llprs(2,orderll(jj))
    ab2=ind2(sh3,sh4)
    maxdmcoul=max(maxdm(ind2(sh1,sh2)),maxdm(ind2(sh3,sh4)))
    maxdmexch=max(maxdm(ind2(sh1,sh3)),maxdm(ind2(sh1,sh4)),maxdm(ind2(sh2,sh3)),&
       & maxdm(ind2(sh2,sh4)))
    maxdmval=max(maxdmcoul,.25_dp*maxdmexch)
    if (schwdd(orderdd(ii))*schwll(orderll(jj))*maxdmval.lt.thr) then
        nprsscr=nprsscr+1
        cycle
    end if
    call i2e_ddll(sh1,sh2,sh3,sh4,eddll)
    if (maxval(abs(eddll)).lt.thr) cycle
    ! The special cases here are ij i=j, and/or k=l.
    samekl=(sh3.eq.sh4)
    do i=1,5
      ia=aosh(sh1)+i-1
      do j=1,5
        ja=aosh(sh2)+j-1
        do k=1,4
          ka=aosh(sh3)+k-1
          do l=1,4
            la=aosh(sh4)+l-1
            if (sameij.and.(i.eq.j).and.samekl.and.(k.eq.l)) then  ! (aa|bb)
              floc(ia,ia)=floc(ia,ia)+p(ka,ka)*eddll(i,i,k,k)
              fac=.5_dp*p(ia,ka)*eddll(i,i,k,k)
              floc(ia,ka)=floc(ia,ka)-fac
              floc(ka,ia)=floc(ka,ia)-fac
              floc(ka,ka)=floc(ka,ka)+p(ia,ia)*eddll(i,i,k,k)
            else if (sameij.and.(i.eq.j).and.samekl.and.(k.lt.l)) then  ! (aa|cb)
              cycle
            else if (sameij.and.(i.eq.j).and.samekl.and.(k.gt.l)) then  !  (aa|bc)
              floc(ia,ia)=floc(ia,ia)+2._dp*p(ka,la)*eddll(i,j,k,l)
              fac=.5_dp*p(ia,ka)*eddll(i,j,k,l)
              floc(ia,la)=floc(ia,la)-fac
              floc(la,ia)=floc(la,ia)-fac
              fac=.5_dp*p(ia,la)*eddll(i,j,k,l)
              floc(ia,ka)=floc(ia,ka)-fac
              floc(ka,ia)=floc(ka,ia)-fac
              fac=p(ia,ia)*eddll(i,j,k,l)
              floc(ka,la)=floc(ka,la)+fac
              floc(la,ka)=floc(la,ka)+fac
            else if (sameij.and.(i.lt.j).and.samekl.and.(k.eq.l)) then    !  (ba|cc)
              cycle
            else if (samekl.and.(k.eq.l).and.sameij.and.(i.gt.j)) then    ! (ab|cc)
              floc(ka,ka)=floc(ka,ka)+2._dp*p(ia,ja)*eddll(i,j,k,l)
              fac=.5_dp*p(ka,ia)*eddll(i,j,k,l)
              floc(ka,ja)=floc(ka,ja)-fac
              floc(ja,ka)=floc(ja,ka)-fac
              fac=.5_dp*p(ka,ja)*eddll(i,j,k,l)
              floc(ka,ia)=floc(ka,ia)-fac
              floc(ia,ka)=floc(ia,ka)-fac
              fac=p(ka,ka)*eddll(i,j,k,l)
              floc(ia,ja)=floc(ia,ja)+fac
              floc(ja,ia)=floc(ja,ia)+fac
            else if (sameij.and.(i.eq.j)) then  !  (aa|bc)
              floc(ia,ia)=floc(ia,ia)+2._dp*p(ka,la)*eddll(i,j,k,l)
              fac=.5_dp*p(ia,ka)*eddll(i,j,k,l)
              floc(ia,la)=floc(ia,la)-fac
              floc(la,ia)=floc(la,ia)-fac
              fac=.5_dp*p(ia,la)*eddll(i,j,k,l)
              floc(ia,ka)=floc(ia,ka)-fac
              floc(ka,ia)=floc(ka,ia)-fac
              fac=p(ia,ia)*eddll(i,j,k,l)
              floc(ka,la)=floc(ka,la)+fac
              floc(la,ka)=floc(la,ka)+fac
            else if (sameij.and.(i.lt.j)) then    !  (ba|cd)
              cycle
            else if (samekl.and.(k.eq.l)) then                  ! (ab|cc)
              floc(ka,ka)=floc(ka,ka)+2._dp*p(ia,ja)*eddll(i,j,k,l)
              fac=.5_dp*p(ka,ia)*eddll(i,j,k,l)
              floc(ka,ja)=floc(ka,ja)-fac
              floc(ja,ka)=floc(ja,ka)-fac
              fac=.5_dp*p(ka,ja)*eddll(i,j,k,l)
              floc(ka,ia)=floc(ka,ia)-fac
              floc(ia,ka)=floc(ia,ka)-fac
              fac=p(ka,ka)*eddll(i,j,k,l)
              floc(ia,ja)=floc(ia,ja)+fac
              floc(ja,ia)=floc(ja,ia)+fac
            else if (samekl.and.(k.lt.l)) then    !  (ba|cd)
              cycle
            else  !   (ab|cd) or maybe (ab|ac)
              fac=2._dp*  p(la,ka) * eddll(i,j,k,l)
              floc(ia,ja) = floc(ia,ja) +  fac
              floc(ja,ia) = floc(ja,ia) +  fac
              fac=.5_dp*p(ja,ka)*eddll(i,j,k,l)
              floc(ia,la) = floc(ia,la) - fac
              floc(la,ia) = floc(la,ia) - fac
              fac=.5_dp*p(ja,la)*eddll(i,j,k,l)
              floc(ia,ka) = floc(ia,ka) - fac
              floc(ka,ia) = floc(ka,ia) - fac
              fac=.5_dp*p(ia,ka)*eddll(i,j,k,l)
              floc(ja,la) = floc(ja,la) - fac
              floc(la,ja) = floc(la,ja) - fac
              fac=.5_dp*p(ia,la)*eddll(i,j,k,l)
              floc(ja,ka) = floc(ja,ka) - fac
              floc(ka,ja) = floc(ka,ja) - fac
              fac=2._dp*p(ja,ia)*eddll(i,j,k,l)
              floc(ka,la) = floc(ka,la) +  fac
              floc(la,ka) = floc(la,ka) +  fac
            end if
          end do
        end do
      end do
    end do
  end do
end do
!$OMP END DO
!$OMP MASTER
if (prtlevl.ge.4) call time_diff_omp(1,"(dd|ll) integ:")
!$OMP END MASTER

!$OMP DO SCHEDULE(DYNAMIC)
do ii=2,ndl  !  (dl|dl)
  sh1=dlprs(1,orderdl(ii))
  sh2=dlprs(2,orderdl(ii))
  ab1=ind2(sh1,sh2)
  ia=aosh(sh1)
  do jj=1,ii-1
    sh3=dlprs(1,orderdl(jj))
    sh4=dlprs(2,orderdl(jj))
    ab2=ind2(sh3,sh4)
    ka=aosh(sh3)
    maxdmcoul=max(maxdm(ind2(sh1,sh2)),maxdm(ind2(sh3,sh4)))
    maxdmexch=max(maxdm(ind2(sh1,sh3)),maxdm(ind2(sh1,sh4)),maxdm(ind2(sh2,sh3)),&
       & maxdm(ind2(sh2,sh4)))
    maxdmval=max(maxdmcoul,.25_dp*maxdmexch)
    if (schwdl(orderdl(ii))*schwdl(orderdl(jj))*maxdmval.lt.thr) then
        nprsscr=nprsscr+1
        cycle
    end if
    call i2e_dldl(sh1,sh2,sh3,sh4,edldl)
    if (maxval(abs(edldl)).lt.thr) cycle
    ! Since the bra and ket shell pairs cannot both be the same, at most we can have
    !  (ab|ac) here, so no need to do anything special.
    do i=1,5
      ia=aosh(sh1)+i-1
      do j=1,4
        ja=aosh(sh2)+j-1
        do k=1,5
          ka=aosh(sh3)+k-1
          do l=1,4
            la=aosh(sh4)+l-1
            fac=2._dp*  p(la,ka) * edldl(i,j,k,l)
            floc(ia,ja) = floc(ia,ja) +  fac
            floc(ja,ia) = floc(ja,ia) +  fac
            fac=.5_dp*p(ja,ka)*edldl(i,j,k,l)
            floc(ia,la) = floc(ia,la) - fac
            floc(la,ia) = floc(la,ia) - fac
            fac=.5_dp*p(ja,la)*edldl(i,j,k,l)
            floc(ia,ka) = floc(ia,ka) - fac
            floc(ka,ia) = floc(ka,ia) - fac
            fac=.5_dp*p(ia,ka)*edldl(i,j,k,l)
            floc(ja,la) = floc(ja,la) - fac
            floc(la,ja) = floc(la,ja) - fac
            fac=.5_dp*p(ia,la)*edldl(i,j,k,l)
            floc(ja,ka) = floc(ja,ka) - fac
            floc(ka,ja) = floc(ka,ja) - fac
            fac=2._dp*p(ja,ia)*edldl(i,j,k,l)
            floc(ka,la) = floc(ka,la) +  fac
            floc(la,ka) = floc(la,ka) +   fac
          end do
        end do
      end do
    end do
  end do
end do
!$OMP END DO
!$OMP MASTER
if (prtlevl.ge.4) call time_diff_omp(1,"(dl|dl) integ:")
!$OMP END MASTER

!$OMP DO SCHEDULE(DYNAMIC)
do ii=1,ndl  ! (dl|ll)
  sh1=dlprs(1,orderdl(ii))
  sh2=dlprs(2,orderdl(ii))
  ab1=ind2(sh1,sh2)
  ia=aosh(sh1)
  do jj=1,nll
    sh3=llprs(1,orderll(jj))
    sh4=llprs(2,orderll(jj))
    ab2=ind2(sh3,sh4)
    maxdmcoul=max(maxdm(ind2(sh1,sh2)),maxdm(ind2(sh3,sh4)))
    maxdmexch=max(maxdm(ind2(sh1,sh3)),maxdm(ind2(sh1,sh4)),maxdm(ind2(sh2,sh3)),&
       & maxdm(ind2(sh2,sh4)))
    maxdmval=max(maxdmcoul,.25_dp*maxdmexch)
    if (schwdl(orderdl(ii))*schwll(orderll(jj))*maxdmval.lt.thr) then
       nprsscr=nprsscr+1
       cycle
    end if
    call i2e_dlll(sh1,sh2,sh3,sh4,edlll)
    if (maxval(abs(edlll)).lt.thr) cycle
    ! The only special case here is if shells sh3 and sh4 are the same.
    samekl=(sh3.eq.sh4)
    do i=1,5
      ia=aosh(sh1)+i-1
      do j=1,4
        ja=aosh(sh2)+j-1
        do k=1,4
          ka=aosh(sh3)+k-1
          do l=1,4
            la=aosh(sh4)+l-1
            if (samekl.and.(k.eq.l)) then  ! when k=l these are (ab|cc)
              floc(ka,ka)=floc(ka,ka)+2._dp*p(ia,ja)*edlll(i,j,k,k)
              fac=.5_dp*p(ka,ia)*edlll(i,j,k,k)
              floc(ka,ja)=floc(ka,ja)-fac
              floc(ja,ka)=floc(ja,ka)-fac
              fac=.5_dp*p(ka,ja)*edlll(i,j,k,k)
              floc(ka,ia)=floc(ka,ia)-fac
              floc(ia,ka)=floc(ia,ka)-fac
              fac=p(ka,ka)*edlll(i,j,k,k)
              floc(ia,ja)=floc(ia,ja)+fac
              floc(ja,ia)=floc(ja,ia)+fac
            else if (samekl.and.(k.lt.l)) then
              cycle
            else
              fac=2._dp*  p(la,ka) * edlll(i,j,k,l)
              floc(ia,ja) = floc(ia,ja) +  fac
              floc(ja,ia) = floc(ja,ia) +  fac
              fac=.5_dp*p(ja,ka)*edlll(i,j,k,l)
              floc(ia,la) = floc(ia,la) - fac
              floc(la,ia) = floc(la,ia) - fac
              fac=.5_dp*p(ja,la)*edlll(i,j,k,l)
              floc(ia,ka) = floc(ia,ka) - fac
              floc(ka,ia) = floc(ka,ia) - fac
              fac=.5_dp*p(ia,ka)*edlll(i,j,k,l)
              floc(ja,la) = floc(ja,la) - fac
              floc(la,ja) = floc(la,ja) - fac
              fac=.5_dp*p(ia,la)*edlll(i,j,k,l)
              floc(ja,ka) = floc(ja,ka) - fac
              floc(ka,ja) = floc(ka,ja) - fac
              fac=2._dp*p(ja,ia)*edlll(i,j,k,l)
              floc(ka,la) = floc(ka,la) +  fac
              floc(la,ka) = floc(la,ka) +  fac
            end if
          end do
        end do
      end do
    end do
  end do
end do
!$OMP END DO
!$OMP MASTER
if (prtlevl.ge.4) call time_diff_omp(1,"(dl|ll) integ:")
!$OMP END MASTER

!$OMP DO SCHEDULE(DYNAMIC)
do ii=2,nll  ! (ll|ll)
  sh1=llprs(1,orderll(ii))
  sh2=llprs(2,orderll(ii))
  ab1=ind2(sh1,sh2)
  sameij=(sh1.eq.sh2)
  do jj=1,ii-1
    sh3=llprs(1,orderll(jj))
    sh4=llprs(2,orderll(jj))
    ab2=ind2(sh3,sh4)
    maxdmcoul=max(maxdm(ind2(sh1,sh2)),maxdm(ind2(sh3,sh4)))
    maxdmexch=max(maxdm(ind2(sh1,sh3)),maxdm(ind2(sh1,sh4)),maxdm(ind2(sh2,sh3)),&
       & maxdm(ind2(sh2,sh4)))
    maxdmval=max(maxdmcoul,.25_dp*maxdmexch)
    if (schwll(orderll(ii))*schwll(orderll(jj))*maxdmval.lt.thr) then
        nprsscr=nprsscr+1
        cycle
    end if
    call i2e_llll(sh1,sh2,sh3,sh4,ellll)
    if (maxval(abs(ellll)).lt.thr) cycle
    ! The special cases here are ij i=j, and/or k=l. It cannot be the case that i=j=k=l or that i=k and j=l, so (aa|aa) and (ab|ab) are ruled out
    samekl=(sh3.eq.sh4)
    do i=1,4
      ia=aosh(sh1)+i-1
      do j=1,4
        ja=aosh(sh2)+j-1
        do k=1,4
          ka=aosh(sh3)+k-1
          do l=1,4
            la=aosh(sh4)+l-1
            if (sameij.and.(i.eq.j).and.samekl.and.(k.eq.l)) then  ! (aa|bb)
              floc(ia,ia)=floc(ia,ia)+p(ka,ka)*ellll(i,i,k,k)
              fac=.5_dp*p(ia,ka)*ellll(i,i,k,k)
              floc(ia,ka)=floc(ia,ka)-fac
              floc(ka,ia)=floc(ka,ia)-fac
              floc(ka,ka)=floc(ka,ka)+p(ia,ia)*ellll(i,i,k,k)
            else if (sameij.and.(i.eq.j).and.samekl.and.(k.lt.l)) then  ! (aa|cb)
              cycle
            else if (sameij.and.(i.eq.j).and.samekl.and.(k.gt.l)) then  !  (aa|bc)
              floc(ia,ia)=floc(ia,ia)+2._dp*p(ka,la)*ellll(i,j,k,l)
              fac=.5_dp*p(ia,ka)*ellll(i,j,k,l)
              floc(ia,la)=floc(ia,la)-fac
              floc(la,ia)=floc(la,ia)-fac
              fac=.5_dp*p(ia,la)*ellll(i,j,k,l)
              floc(ia,ka)=floc(ia,ka)-fac
              floc(ka,ia)=floc(ka,ia)-fac
              fac=p(ia,ia)*ellll(i,j,k,l)
              floc(ka,la)=floc(ka,la)+fac
              floc(la,ka)=floc(la,ka)+fac
            else if (sameij.and.(i.lt.j).and.samekl.and.(k.eq.l)) then    !  (ba|cc)
              cycle
            else if (samekl.and.(k.eq.l).and.sameij.and.(i.gt.j)) then    ! (ab|cc)
              floc(ka,ka)=floc(ka,ka)+2._dp*p(ia,ja)*ellll(i,j,k,l)
              fac=.5_dp*p(ka,ia)*ellll(i,j,k,l)
              floc(ka,ja)=floc(ka,ja)-fac
              floc(ja,ka)=floc(ja,ka)-fac
              fac=.5_dp*p(ka,ja)*ellll(i,j,k,l)
              floc(ka,ia)=floc(ka,ia)-fac
              floc(ia,ka)=floc(ia,ka)-fac
              fac=p(ka,ka)*ellll(i,j,k,l)
              floc(ia,ja)=floc(ia,ja)+fac
              floc(ja,ia)=floc(ja,ia)+fac
            else if (sameij.and.(i.eq.j)) then  !  (aa|bc)
              floc(ia,ia)=floc(ia,ia)+2._dp*p(ka,la)*ellll(i,j,k,l)
              fac=.5_dp*p(ia,ka)*ellll(i,j,k,l)
              floc(ia,la)=floc(ia,la)-fac
              floc(la,ia)=floc(la,ia)-fac
              fac=.5_dp*p(ia,la)*ellll(i,j,k,l)
              floc(ia,ka)=floc(ia,ka)-fac
              floc(ka,ia)=floc(ka,ia)-fac
              fac=p(ia,ia)*ellll(i,j,k,l)
              floc(ka,la)=floc(ka,la)+fac
              floc(la,ka)=floc(la,ka)+fac
            else if (sameij.and.(i.lt.j)) then    !  (ba|cd)
              cycle
            else if (samekl.and.(k.eq.l)) then                  ! (ab|cc)
              floc(ka,ka)=floc(ka,ka)+2._dp*p(ia,ja)*ellll(i,j,k,l)
              fac=.5_dp*p(ka,ia)*ellll(i,j,k,l)
              floc(ka,ja)=floc(ka,ja)-fac
              floc(ja,ka)=floc(ja,ka)-fac
              fac=.5_dp*p(ka,ja)*ellll(i,j,k,l)
              floc(ka,ia)=floc(ka,ia)-fac
              floc(ia,ka)=floc(ia,ka)-fac
              fac=p(ka,ka)*ellll(i,j,k,l)
              floc(ia,ja)=floc(ia,ja)+fac
              floc(ja,ia)=floc(ja,ia)+fac
            else if (samekl.and.(k.lt.l)) then    !  (ba|cd)
              cycle
            else  !   (ab|cd) or maybe (ab|ac)
              fac=2._dp*  p(la,ka) * ellll(i,j,k,l)
              floc(ia,ja) = floc(ia,ja) +  fac
              floc(ja,ia) = floc(ja,ia) +  fac
              fac=.5_dp*p(ja,ka)*ellll(i,j,k,l)
              floc(ia,la) = floc(ia,la) - fac
              floc(la,ia) = floc(la,ia) - fac
              fac=.5_dp*p(ja,la)*ellll(i,j,k,l)
              floc(ia,ka) = floc(ia,ka) - fac
              floc(ka,ia) = floc(ka,ia) - fac
              fac=.5_dp*p(ia,ka)*ellll(i,j,k,l)
              floc(ja,la) = floc(ja,la) - fac
              floc(la,ja) = floc(la,ja) - fac
              fac=.5_dp*p(ia,la)*ellll(i,j,k,l)
              floc(ja,ka) = floc(ja,ka) - fac
              floc(ka,ja) = floc(ka,ja) - fac
              fac=2._dp*p(ja,ia)*ellll(i,j,k,l)
              floc(ka,la) = floc(ka,la) +  fac
              floc(la,ka) = floc(la,ka) +  fac
            end if
          end do
        end do
      end do
    end do
  end do
end do
!$OMP END DO
!$OMP MASTER
if (prtlevl.ge.4) call time_diff_omp(1,"(ll|ll) integ:")
!$OMP END MASTER

!$OMP CRITICAL
f=f+floc
!$OMP END CRITICAL
!$OMP END PARALLEL


return
end subroutine fock_build_non_diag

