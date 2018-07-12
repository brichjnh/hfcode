subroutine fock_build_initial_diag(f,p,thr)
use nrtype; use molprops ; use s_sp_l_terms ; use s_d_l_terms ; use s_p_d_terms
implicit none

! This is the routine called to build the first Fock matrix - so the Schwarz integrals are not yet available

real(dp), intent(in) :: p(nb,nb), thr
real(dp), intent(inout) :: f(nb,nb)

integer(i4b) :: i, j, k, l, ii, ia, ja, ka, la, sh1, sh2, ao1, ao2
real(dp) :: essss, espsp(3,3), esdsd(5,5), epppp(3,3,3,3), edddd(5,5,5,5), fac
real(dp) :: eslsl(4,4), eplpl(3,4,3,4), ellll(4,4,4,4), epdpd(3,5,3,5), edldl(5,4,5,4)

! The basic idea is this: (where a and b can denote either a basis function or, more loosely, a basis function shell)
! F(a,b) = F(a,b) +    p(d,c) * i2e(a,b,c,d)
! F(a,d) = F(a,d) - .5 p(b,c) * i2e(a,b,c,d)

! I use an ordered list of shell pairs (a,b), set up previously.
! In this first version, this list is not sorted by magnitude, only by order of the functions. I use the existing lists of 
!   s-s, s-p, s-d, ... pairs
! I loop over these pairs.
!    For each shell pair, I first identify the corresponding sub-elements of the density matrix for the Coulomb part,
!    and find the largest one (in absolute value)
!    I also find the largest one (in absolute value) relevant for the exchange part.
!    Now I do a second loop, over all shell pairs >= that one.
!        using the density matrix element, I perform the p-weighted Schwarz test, for Coulomb and Exchange
!        If one of these tests is passed, I calculate the integral, which I do by calling a generic routine, that returns
!        a (5,5,5,5) array containing the integral values. For shell sets other than (dd|dd) only some elements are assigned.
!        Using these values, I increment the relevant Fock matrix elements, and cycle

do ii=1,nss
  sh1=ssprs(1,ii)
  sh2=ssprs(2,ii)
  call i2e_ssss(sh1,sh2,sh1,sh2,essss)
  schwss(ii)=sqrt(abs(essss))
  if (abs(essss).gt.intthresh) then
    ao1=aosh(sh1)
    ao2=aosh(sh2)
    if (ao1.eq.ao2) then   ! (aa|aa)
      f(ao1,ao1)=f(ao1,ao1)+.5_dp*p(ao1,ao1)*essss
    else                   ! (ab|ab)
      fac=1.5_dp*p(ao1,ao2)*essss
      f(ao2,ao1)=f(ao2,ao1)+fac
      f(ao1,ao2)=f(ao1,ao2)+fac
      f(ao1,ao1)=f(ao1,ao1)-.5_dp*p(ao2,ao2)*essss
      f(ao2,ao2)=f(ao2,ao2)-.5_dp*p(ao1,ao1)*essss
    end if
  end if
end do
if (nss.gt.0) then
call index_sort(nss,schwss,orderss)
call find_threshold(nss,schwss,orderss,shthresh,lastss)
end if

do ii=1,nsp
  sh1=spprs(1,ii)
  sh2=spprs(2,ii)
  call i2e_spsp(sh1,sh2,sh1,sh2,espsp)
  schwsp(ii)=sqrt(maxval(abs(espsp)))
  ! The integrals of type (spi|spi) are (ab|ab), whereas (spi|spj) are (ab|ac)
  ao1=aosh(sh1)
  do i=1,3
    if (abs(espsp(i,i)).gt.intthresh) then   ! the cases (spi|spi) are abab
      ao2=aosh(sh2)+i-1
      fac=1.5_dp*p(ao1,ao2)*espsp(i,i)
      f(ao2,ao1)=f(ao2,ao1)+fac
      f(ao1,ao2)=f(ao1,ao2)+fac
      f(ao1,ao1)=f(ao1,ao1)-.5_dp*p(ao2,ao2)*espsp(i,i)
      f(ao2,ao2)=f(ao2,ao2)-.5_dp*p(ao1,ao1)*espsp(i,i)
    end if
    do j=i+1,3
      if (abs(espsp(i,j)).gt.intthresh) then   ! the cases (spi|spj) , i.ne.j, are abac
         ia=ao1; ja=aosh(sh2)+i-1 ; ka=ao1; la=aosh(sh2)+j-1
         fac=2._dp*p(la,ka) * espsp(i,j)
         F(ia,ja) = F(ia,ja) + fac
         F(ja,ia) = F(ja,ia) + fac
         fac=2._dp*p(ja,ia) * espsp(i,j)
         F(ka,la) = F(ka,la) + fac
         F(la,ka) = F(la,ka) + fac
         fac=.5_dp*p(ja,ka) * espsp(i,j)
         F(ia,la) = F(ia,la) - fac
         F(la,ia) = F(la,ia) - fac
         fac=.5_dp*p(ja,la) * espsp(i,j)
         F(ia,ka) = F(ia,ka) - fac
         F(ka,ia) = F(ka,ia) - fac
         fac=.5_dp*p(ia,ka) * espsp(i,j)
         F(ja,la) = F(ja,la) - fac
         F(la,ja) = F(la,ja) - fac
         fac=.5_dp*p(ia,la) * espsp(i,j)
         F(ja,ka) = F(ja,ka) - fac
         F(ka,ja) = F(ka,ja) - fac
      end if
    end do
  end do
end do
if (nsp.gt.0) then
call index_sort(nsp,schwsp,ordersp)
call find_threshold(nsp,schwsp,ordersp,shthresh,lastsp)
end if

do ii=1,nsd
  sh1=sdprs(1,ii)
  sh2=sdprs(2,ii)
  call i2e_sdsd(sh1,sh2,sh1,sh2,esdsd)
  schwsd(ii)=sqrt(maxval(abs(esdsd)))
  ! The integrals of type (sdi|sdi) are (ab|ab), others are not
  ao1=aosh(sh1)
  do i=1,5
    if (abs(esdsd(i,i)).gt.intthresh) then   ! the cases (spi|spi) are abab
      ao2=aosh(sh2)+i-1
      fac=1.5_dp*p(ao1,ao2)*esdsd(i,i)
      f(ao2,ao1)=f(ao2,ao1)+fac
      f(ao1,ao2)=f(ao1,ao2)+fac
      f(ao1,ao1)=f(ao1,ao1)-.5_dp*p(ao2,ao2)*esdsd(i,i)
      f(ao2,ao2)=f(ao2,ao2)-.5_dp*p(ao1,ao1)*esdsd(i,i)
    end if
    do j=i+1,5
      if (abs(esdsd(i,j)).gt.intthresh) then   ! the cases (ss'|sp) and (spi|spj) , i.ne.j, are abac
         ia=ao1; ja=aosh(sh2)+i-1 ; ka=ao1; la=aosh(sh2)+j-1
         fac=2._dp*p(la,ka) * esdsd(i,j)
         F(ia,ja) = F(ia,ja) + fac
         F(ja,ia) = F(ja,ia) + fac
         fac=2._dp*p(ja,ia) * esdsd(i,j)
         F(ka,la) = F(ka,la) + fac
         F(la,ka) = F(la,ka) + fac
         fac=.5_dp*p(ja,ka) * esdsd(i,j)
         F(ia,la) = F(ia,la) - fac
         F(la,ia) = F(la,ia) - fac
         fac=.5_dp*p(ja,la) * esdsd(i,j)
         F(ia,ka) = F(ia,ka) - fac
         F(ka,ia) = F(ka,ia) - fac
         fac=.5_dp*p(ia,ka) * esdsd(i,j)
         F(ja,la) = F(ja,la) - fac
         F(la,ja) = F(la,ja) - fac
         fac=.5_dp*p(ia,la) * esdsd(i,j)
         F(ja,ka) = F(ja,ka) - fac
         F(ka,ja) = F(ka,ja) - fac
      end if
    end do
  end do
end do
if (nsd.gt.0) then
call index_sort(nsd,schwsd,ordersd)
call find_threshold(nsd,schwsd,ordersd,shthresh,lastsd)
end if

do ii=1,nsl
  sh1=slprs(1,ii)
  sh2=slprs(2,ii)
  call i2e_slsl(sh1,sh2,sh1,sh2,eslsl)
  schwsl(ii)=sqrt(maxval(abs(eslsl)))
  ! The integrals of type (ss'|ss') and (spi|spi) are (ab|ab), whereas (ss'|sp) and (spi|spj) are (ab|ac)
  ao1=aosh(sh1)
  do i=1,4
    if (abs(eslsl(i,i)).gt.intthresh) then   ! the cases (spi|spi) are abab
      ao2=aosh(sh2)+i-1
      fac=1.5_dp*p(ao1,ao2)*eslsl(i,i)
      f(ao2,ao1)=f(ao2,ao1)+fac
      f(ao1,ao2)=f(ao1,ao2)+fac
      f(ao1,ao1)=f(ao1,ao1)-.5_dp*p(ao2,ao2)*eslsl(i,i)
      f(ao2,ao2)=f(ao2,ao2)-.5_dp*p(ao1,ao1)*eslsl(i,i)
    end if
    do j=i+1,4
      if (abs(eslsl(i,j)).gt.intthresh) then   ! the cases (ss'|sp) and (spi|spj) , i.ne.j, are abac
         ia=ao1; ja=aosh(sh2)+i-1 ; ka=ao1; la=aosh(sh2)+j-1
         fac=2._dp*p(la,ka) * eslsl(i,j)
         F(ia,ja) = F(ia,ja) + fac
         F(ja,ia) = F(ja,ia) + fac
         fac=2._dp*p(ja,ia) * eslsl(i,j)
         F(ka,la) = F(ka,la) + fac
         F(la,ka) = F(la,ka) + fac
         fac=.5_dp*p(ja,ka) * eslsl(i,j)
         F(ia,la) = F(ia,la) - fac
         F(la,ia) = F(la,ia) - fac
         fac=.5_dp*p(ja,la) * eslsl(i,j)
         F(ia,ka) = F(ia,ka) - fac
         F(ka,ia) = F(ka,ia) - fac
         fac=.5_dp*p(ia,ka) * eslsl(i,j)
         F(ja,la) = F(ja,la) - fac
         F(la,ja) = F(la,ja) - fac
         fac=.5_dp*p(ia,la) * eslsl(i,j)
         F(ja,ka) = F(ja,ka) - fac
         F(ka,ja) = F(ka,ja) - fac
      end if
    end do
  end do
end do
if (nsl.gt.0) then
call index_sort(nsl,schwsl,ordersl)
call find_threshold(nsl,schwsl,ordersl,shthresh,lastsl)
end if

do ii=1,npp
  sh1=ppprs(1,ii)
  sh2=ppprs(2,ii)
  call i2e_pppp(sh1,sh2,sh1,sh2,epppp)
  schwpp(ii)=sqrt(maxval(abs(epppp)))
  ! Lots of different types among these 81 integrals, some of which are non-unique...
  ! Two main cases: the shells are the same, or they are not. 
  ! Start with same shell.
  if (sh1.eq.sh2) then  ! MAjor if statement - either sh1.eq.sh2 or not

    ! Type 1: there are three (pipi|pipi) which are (aa|aa)
    do i=1,3
      ia=aosh(sh1)+i-1
      if (abs(epppp(i,i,i,i)).gt.intthresh) then 
        f(ia,ia)=f(ia,ia)+.5_dp*p(ia,ia)*epppp(i,i,i,i)
      end if
      ! Type 2: three unique (pipi|pjpj), i.ne.j, that are (aa|bb) (6 in total: 9)
      do j=i+1,3
        if (abs(epppp(i,i,j,j)).gt.intthresh) then
          ja=ia
          ka=aosh(sh1)+j-1; la=ka
          f(ia,ia)=f(ia,ia)+p(ka,ka)*epppp(i,i,j,j)
          fac=.5_dp*p(ia,ka)*epppp(i,i,j,j)
          f(ia,ka)=f(ia,ka)-fac
          f(ka,ia)=f(ka,ia)-fac
          f(ka,ka)=f(ka,ka)+p(ia,ia)*epppp(i,i,j,j)
        end if
        ! Type 3: three unique (pipj|pipj), i.ne.j, that are (ab|ab) (12 in total: 21)
        if (abs(epppp(i,j,i,j)).gt.intthresh) then   
          ja=aosh(sh1)+j-1 ; ka=ia; la=ja
          fac=1.5_dp*p(ia,ja)*epppp(i,j,i,j)
          f(ia,ja)=f(ia,ja)+fac
          f(ja,ia)=f(ja,ia)+fac
          f(ia,ia)=f(ia,ia)-.5_dp*p(ja,ja)*epppp(i,j,i,j)
          f(ja,ja)=f(ja,ja)-.5_dp*p(ia,ia)*epppp(i,j,i,j)
        end if
      end do
      ! Type 4: a set of (pipi|pjpk) where j<>k (but either can be same as i)
      ! These are (aa|bc) or (aa|ab), both of which have the same permutants
      do j=1,2
        do k=j+1,3
          if (abs(epppp(i,i,j,k)).gt.intthresh) then
            ja=ia ; ka=aosh(sh1)+j-1 ; la=aosh(sh1)+k-1
            f(ia,ia)=f(ia,ia)+2._dp*p(ka,la)*epppp(i,i,j,k)
            fac=.5_dp*p(ia,ka)*epppp(i,i,j,k)
            f(ia,la)=f(ia,la)-fac
            f(la,ia)=f(la,ia)-fac
            fac=.5_dp*p(ia,la)*epppp(i,i,j,k)
            f(ia,ka)=f(ia,ka)-fac
            f(ka,ia)=f(ka,ia)-fac
            fac=p(ia,ia)*epppp(i,i,j,k)
            f(ka,la)=f(ka,la)+fac
            f(la,ka)=f(la,ka)+fac
          end if
          ! Type 4b: the related (pjpk|pipi) where j<>k (but either can be same as i)
          ! These are (bc|aa) or (ab|aa), both of which have the same permutants
          ! Note these are identical to the (pipi|pjpk) in this case but easier to code like ths
          if (abs(epppp(j,k,i,i)).gt.intthresh) then
            ia=aosh(sh1)+j-1 ; ja=aosh(sh1)+k-1 ; ka=aosh(sh1)+i-1 ; la=ka
            f(ka,ka)=f(ka,ka)+2._dp*p(ia,ja)*epppp(j,k,i,i)
            fac=.5_dp*p(ka,ia)*epppp(j,k,i,i)
            f(ka,ja)=f(ka,ja)-fac
            f(ja,ka)=f(ja,ka)-fac
            fac=.5_dp*p(ka,ja)*epppp(j,k,i,i)
            f(ka,ia)=f(ka,ia)-fac
            f(ia,ka)=f(ia,ka)-fac
            fac=p(ka,ka)*epppp(j,k,i,i)
            f(ia,ja)=f(ia,ja)+fac
            f(ja,ia)=f(ja,ia)+fac
          end if
        end do
      end do
      ! Finally type 5: (pipj|pipk) where both j and k are different from i
      ! these are (ab|ac)
      if (i.eq.1) then
        j=2 ; k = 3
        ia=aosh(sh1) ; ja=aosh(sh1)+1 ; ka=ia ; la=aosh(sh1)+2
      else if (i.eq.2) then
        j=1 ; k = 3
        ia=aosh(sh1)+1 ; ja=aosh(sh1) ; ka=ia ; la=aosh(sh1)+2
      else if (i.eq.3) then
        j=1 ; k = 2
        ia=aosh(sh1)+2 ; ja=aosh(sh1) ; ka=ia ; la=aosh(sh1)+1
      end if
      if (abs(epppp(i,j,i,k)).gt.thr) then
        fac=2._dp*p(la,ka) * epppp(i,j,i,k)
        F(ia,ja) = F(ia,ja) + fac
        F(ja,ia) = F(ja,ia) + fac
        fac=2._dp*p(ja,ia) * epppp(i,j,i,k)
        F(ka,la) = F(ka,la) + fac
        F(la,ka) = F(la,ka) + fac
        fac=.5_dp*p(ja,ka) * epppp(i,j,i,k)
        F(ia,la) = F(ia,la) - fac
        F(la,ia) = F(la,ia) - fac
        fac=.5_dp*p(ja,la) * epppp(i,j,i,k)
        F(ia,ka) = F(ia,ka) - fac
        F(ka,ia) = F(ka,ia) - fac
        fac=.5_dp*p(ia,ka) * epppp(i,j,i,k)
        F(ja,la) = F(ja,la) - fac
        F(la,ja) = F(la,ja) - fac
        fac=.5_dp*p(ia,la) * epppp(i,j,i,k)
        F(ja,ka) = F(ja,ka) - fac
        F(ka,ja) = F(ka,ja) - fac
      end if
    end do
  ! Second case: the two shells are not the same
  ! This is simpler because now the highest permutation symm is (ab|ab). 9 of these.
  ! Still need to be careful, though, to only use each permutationally identical integral once
  else  ! This is the major if
    do i=1,3
      ia=aosh(sh1)+i-1; ka=ia
      do j=1,3
        if (abs(epppp(i,j,i,j)).gt.intthresh) then
          ja=aosh(sh2)+j-1; la=ja
          fac=1.5_dp*p(ia,ja)*epppp(i,j,i,j)
          f(ia,ja)=f(ia,ja)+fac
          f(ja,ia)=f(ja,ia)+fac
          f(ia,ia)=f(ia,ia)-.5_dp*p(ja,ja)*epppp(i,j,i,j)
          f(ja,ja)=f(ja,ja)-.5_dp*p(ia,ia)*epppp(i,j,i,j)
        end if
        do k=1,3
          do l=1,3
            if ((abs(epppp(i,j,k,l)).gt.intthresh).and. &
              &   ((i.lt.k).or.((i.eq.k).and.(j.lt.l)))) then
              ja=aosh(sh2)+j-1 ; ka=aosh(sh1)+k-1 ; la=aosh(sh2)+l-1 
              fac=2._dp*p(la,ka) * epppp(i,j,k,l)
              F(ia,ja) = F(ia,ja) + fac
              F(ja,ia) = F(ja,ia) + fac
              fac=2._dp*p(ja,ia) * epppp(i,j,k,l)
              F(ka,la) = F(ka,la) + fac
              F(la,ka) = F(la,ka) + fac
              fac=.5_dp*p(ja,ka) * epppp(i,j,k,l)
              F(ia,la) = F(ia,la) - fac
              F(la,ia) = F(la,ia) - fac
              fac=.5_dp*p(ja,la) * epppp(i,j,k,l)
              F(ia,ka) = F(ia,ka) - fac
              F(ka,ia) = F(ka,ia) - fac
              fac=.5_dp*p(ia,ka) * epppp(i,j,k,l)
              F(ja,la) = F(ja,la) - fac
              F(la,ja) = F(la,ja) - fac
              fac=.5_dp*p(ia,la) * epppp(i,j,k,l)
              F(ja,ka) = F(ja,ka) - fac
              F(ka,ja) = F(ka,ja) - fac
            end if
          end do
        end do
      end do
    end do
  end if  ! finished the if

end do ! the 'do' over pp shell pairs
if (npp.gt.0) then
call index_sort(npp,schwpp,orderpp)
call find_threshold(npp,schwpp,orderpp,shthresh,lastpp)
end if

do ii=1,npd
  sh1=pdprs(1,ii)
  sh2=pdprs(2,ii)
  call i2e_pdpd(sh1,sh2,sh1,sh2,epdpd)
  schwpd(ii)=sqrt(maxval(abs(epdpd)))
  ! The integrals of type (pidj|pidj) are (ab|ab), whereas all others are (ab|ac) or (ab|cd)
  do i=1,3
    ia=aosh(sh1)+i-1
    do j=1,5
      ja=aosh(sh2)+j-1
      if (abs(epdpd(i,j,i,j)).gt.intthresh) then   ! the cases (pilj|pilj) are abab
        fac=1.5_dp*p(ia,ja)*epdpd(i,j,i,j)
        f(ja,ia)=f(ja,ia)+fac
        f(ia,ja)=f(ia,ja)+fac
        f(ia,ia)=f(ia,ia)-.5_dp*p(ja,ja)*epdpd(i,j,i,j)
        f(ja,ja)=f(ja,ja)-.5_dp*p(ia,ia)*epdpd(i,j,i,j)
      end if
      do k=1,3
        ka=aosh(sh1)+k-1
        do l=1,5
          la=aosh(sh2)+l-1
          if ((abs(epdpd(i,j,k,l)).gt.intthresh).and. &
               &   ((i.lt.k).or.((i.eq.k).and.(j.lt.l)))) then   ! the cases (spi|spj) , i.ne.j, are abac
            ka=aosh(sh1)+k-1; la=aosh(sh2)+l-1
            fac=2._dp*p(la,ka) * epdpd(i,j,k,l)
            F(ia,ja) = F(ia,ja) + fac
            F(ja,ia) = F(ja,ia) + fac
            fac=2._dp*p(ja,ia) * epdpd(i,j,k,l)
            F(ka,la) = F(ka,la) + fac
            F(la,ka) = F(la,ka) + fac
            fac=.5_dp*p(ja,ka) * epdpd(i,j,k,l)
            F(ia,la) = F(ia,la) - fac
            F(la,ia) = F(la,ia) - fac
            fac=.5_dp*p(ja,la) * epdpd(i,j,k,l)
            F(ia,ka) = F(ia,ka) - fac
            F(ka,ia) = F(ka,ia) - fac
            fac=.5_dp*p(ia,ka) * epdpd(i,j,k,l)
            F(ja,la) = F(ja,la) - fac
            F(la,ja) = F(la,ja) - fac
            fac=.5_dp*p(ia,la) * epdpd(i,j,k,l)
            F(ja,ka) = F(ja,ka) - fac
            F(ka,ja) = F(ka,ja) - fac
          end if
        end do
      end do
    end do
  end do
end do
if (npd.gt.0) then
call index_sort(npd,schwpd,orderpd)
call find_threshold(npd,schwpd,orderpd,shthresh,lastpd)
end if

do ii=1,npl
  sh1=plprs(1,ii)
  sh2=plprs(2,ii)
  call i2e_plpl(sh1,sh2,sh1,sh2,eplpl)
  schwpl(ii)=sqrt(maxval(abs(eplpl)))
  ! The integrals of type (pilj|pilj) are (ab|ab), whereas all others are (ab|ac) or (ab|cd)
  do i=1,3
    ia=aosh(sh1)+i-1
    do j=1,4
      ja=aosh(sh2)+j-1
      if (abs(eplpl(i,j,i,j)).gt.intthresh) then   ! the cases (pilj|pilj) are abab
        fac=1.5_dp*p(ia,ja)*eplpl(i,j,i,j)
        f(ja,ia)=f(ja,ia)+fac
        f(ia,ja)=f(ia,ja)+fac
        f(ia,ia)=f(ia,ia)-.5_dp*p(ja,ja)*eplpl(i,j,i,j)
        f(ja,ja)=f(ja,ja)-.5_dp*p(ia,ia)*eplpl(i,j,i,j)
      end if
      do k=1,3
        ka=aosh(sh1)+k-1
        do l=1,4
          la=aosh(sh2)+l-1
          if ((abs(eplpl(i,j,k,l)).gt.intthresh).and. &
               &   ((i.lt.k).or.((i.eq.k).and.(j.lt.l)))) then   ! the cases (spi|spj) , i.ne.j, are abac
            ka=aosh(sh1)+k-1; la=aosh(sh2)+l-1
            fac=2._dp*p(la,ka) * eplpl(i,j,k,l)
            F(ia,ja) = F(ia,ja) + fac
            F(ja,ia) = F(ja,ia) + fac
            fac=2._dp*p(ja,ia) * eplpl(i,j,k,l)
            F(ka,la) = F(ka,la) + fac
            F(la,ka) = F(la,ka) + fac
            fac=.5_dp*p(ja,ka) * eplpl(i,j,k,l)
            F(ia,la) = F(ia,la) - fac
            F(la,ia) = F(la,ia) - fac
            fac=.5_dp*p(ja,la) * eplpl(i,j,k,l)
            F(ia,ka) = F(ia,ka) - fac
            F(ka,ia) = F(ka,ia) - fac
            fac=.5_dp*p(ia,ka) * eplpl(i,j,k,l)
            F(ja,la) = F(ja,la) - fac
            F(la,ja) = F(la,ja) - fac
            fac=.5_dp*p(ia,la) * eplpl(i,j,k,l)
            F(ja,ka) = F(ja,ka) - fac
            F(ka,ja) = F(ka,ja) - fac
          end if
        end do
      end do
    end do
  end do
end do
if (npl.gt.0) then
call index_sort(npl,schwpl,orderpl)
call find_threshold(npl,schwpl,orderpl,shthresh,lastpl)
end if

do ii=1,ndd
  sh1=ddprs(1,ii)
  sh2=ddprs(2,ii)
  call i2e_dddd(sh1,sh2,sh1,sh2,edddd)
  schwdd(ii)=sqrt(maxval(abs(edddd)))
  ! Lots of different types among these 81 integrals, some of which are non-unique...
  ! Two main cases: the shells are the same, or they are not. 
  ! Start with same shell.
  if (sh1.eq.sh2) then  ! MAjor if statement - either sh1.eq.sh2 or not
    ! Type 1: there are some (didi|didi) which are (aa|aa)
    do i=1,5
      ia=aosh(sh1)+i-1
      do j=1,5
        ja=aosh(sh1)+j-1
        do k=1,5
          ka=aosh(sh1)+k-1
          do l=1,5
            la=aosh(sh1)+l-1
            if (abs(edddd(i,j,k,l)).gt.intthresh) then 
              if ((i.eq.j).and.(j.eq.k).and.(k.eq.l)) then  ! (aa|aa), 5 of these
                f(ia,ia)=f(ia,ia)+.5_dp*p(ia,ia)*edddd(i,i,i,i)
              else if ((i.eq.j).and.(k.eq.l).and.(j.lt.k)) then  ! (aa|bb)
                f(ia,ia)=f(ia,ia)+p(ka,ka)*edddd(i,i,k,k)
                fac=.5_dp*p(ia,ka)*edddd(i,i,k,k)
                f(ia,ka)=f(ia,ka)-fac
                f(ka,ia)=f(ka,ia)-fac
                f(ka,ka)=f(ka,ka)+p(ia,ia)*edddd(i,i,k,k)
!              else if ((i.eq.j).and.(k.eq.l).and.(j.gt.k)) then   ! (bb|aa)
!                cycle
              else if ((i.eq.k).and.(j.eq.l).and.(i.lt.j)) then   ! (ab|ab)
                fac=1.5_dp*p(ia,ja)*edddd(i,j,i,j)
                f(ia,ja)=f(ia,ja)+fac
                f(ja,ia)=f(ja,ia)+fac
                f(ia,ia)=f(ia,ia)-.5_dp*p(ja,ja)*edddd(i,j,i,j)
                f(ja,ja)=f(ja,ja)-.5_dp*p(ia,ia)*edddd(i,j,i,j)
!              else if ((i.eq.l).and.(j.eq.k).and.(i.lt.j)) then  ! (ab|ba)
!                cycle
!              else if ((i.eq.l).and.(j.eq.k).and.(i.gt.j)) then  ! (ba|ab)
!                cycle
!              else if ((i.eq.k).and.(j.eq.l).and.(i.gt.j)) then  ! (ba|ba)
!                cycle
              else if ((i.eq.j).and.(k.lt.l)) then       !  (aa|bc)
                f(ia,ia)=f(ia,ia)+2._dp*p(ka,la)*edddd(i,i,k,l)
                fac=.5_dp*p(ia,ka)*edddd(i,i,k,l)
                f(ia,la)=f(ia,la)-fac
                f(la,ia)=f(la,ia)-fac
                fac=.5_dp*p(ia,la)*edddd(i,i,k,l)
                f(ia,ka)=f(ia,ka)-fac
                f(ka,ia)=f(ka,ia)-fac
                fac=p(ia,ia)*edddd(i,i,k,l)
                f(ka,la)=f(ka,la)+fac
                f(la,ka)=f(la,ka)+fac
!              else if ((i.eq.j).and.(k.gt.l)) then       !  (aa|cb)
!                cycle
!              else if ((i.lt.j).and.(k.eq.l)) then        ! (ab|cc)
!                cycle
!              else if ((i.gt.j).and.(k.eq.l)) then        ! (ba|cc)
!                cycle
              else if ((i.eq.k).and.(j.lt.l)) then    ! (ab|ac)
                fac=2._dp*p(la,ka) * edddd(i,j,k,l)
                F(ia,ja) = F(ia,ja) + fac
                F(ja,ia) = F(ja,ia) + fac
                fac=2._dp*p(ja,ia) * edddd(i,j,k,l)
                F(ka,la) = F(ka,la) + fac
                F(la,ka) = F(la,ka) + fac
                fac=.5_dp*p(ja,ka) * edddd(i,j,k,l)
                F(ia,la) = F(ia,la) - fac
                F(la,ia) = F(la,ia) - fac
                fac=.5_dp*p(ja,la) * edddd(i,j,k,l)
                F(ia,ka) = F(ia,ka) - fac
                F(ka,ia) = F(ka,ia) - fac
                fac=.5_dp*p(ia,ka) * edddd(i,j,k,l)
                F(ja,la) = F(ja,la) - fac
                F(la,ja) = F(la,ja) - fac
                fac=.5_dp*p(ia,la) * edddd(i,j,k,l)
                F(ja,ka) = F(ja,ka) - fac
                F(ka,ja) = F(ka,ja) - fac
              else if (((i.lt.j).and.(j.lt.k).and.(k.lt.l)).or. &    ! (ab|cd)
                    &  ((i.lt.k).and.(k.lt.j).and.(j.lt.l)).or. &    ! (ac|bd)
                    &  ((i.lt.l).and.(l.lt.k).and.(k.lt.j))) then    ! (ad|cb)
                fac=2._dp*p(la,ka) * edddd(i,j,k,l)
                F(ia,ja) = F(ia,ja) + fac
                F(ja,ia) = F(ja,ia) + fac
                fac=2._dp*p(ja,ia) * edddd(i,j,k,l)
                F(ka,la) = F(ka,la) + fac
                F(la,ka) = F(la,ka) + fac
                fac=.5_dp*p(ja,ka) * edddd(i,j,k,l)
                F(ia,la) = F(ia,la) - fac
                F(la,ia) = F(la,ia) - fac
                fac=.5_dp*p(ja,la) * edddd(i,j,k,l)
                F(ia,ka) = F(ia,ka) - fac
                F(ka,ia) = F(ka,ia) - fac
                fac=.5_dp*p(ia,ka) * edddd(i,j,k,l)
                F(ja,la) = F(ja,la) - fac
                F(la,ja) = F(la,ja) - fac
                fac=.5_dp*p(ia,la) * edddd(i,j,k,l)
                F(ja,ka) = F(ja,ka) - fac
                F(ka,ja) = F(ka,ja) - fac
!              else if ((i.eq.l).and.(j.lt.k)) then    ! (ab|ca)
!                cycle
!              else if ((i.eq.k).and.(j.gt.l)) then    ! (ac|ab)
!                cycle
!              else if ((i.eq.l).and.(j.gt.k)) then    ! (ac|ba)
!                cycle
!              else if ((j.eq.l).and.(i.lt.k)) then    ! (ba|ca)
!                cycle
!              else if ((j.eq.k).and.(i.lt.l)) then    ! (ba|ac)
!                cycle
!              else if ((j.eq.l).and.(i.gt.k)) then    ! (ca|ba)
!                cycle
!              else if ((j.eq.k).and.(i.gt.l)) then    ! (ca|ab)
!                cycle
!              else    ! Other cases which are 8-fold degenerate - just do the normal thing and all eight should survive
!                f(ia,ja)=f(ia,ja)+p(ka,la)*edddd(i,j,k,l)
!                f(ia,la)=f(ia,la)-.5_dp*p(ja,ka)*edddd(i,j,k,l)
              end if
            end if
          end do
        end do
      end do
    end do
  ! Second case: the two shells are not the same
  ! This is simpler because now the highest permutation symm is (ab|ab). 16 of these.
  ! Still need to be careful, though, to only use each permutationally identical integral once
  else  ! This is the major if
    do i=1,5
      ia=aosh(sh1)+i-1
      do j=1,5
        ja=aosh(sh2)+j-1
        do k=1,5
          ka=aosh(sh1)+k-1
          do l=1,5
            la=aosh(sh2)+l-1
            !if (abs(edddd(i,j,k,l)).gt.intthresh) then 
              if ((i.eq.k).and.(j.eq.l)) then   ! (ab|ab)
                fac=1.5_dp*p(ia,ja)*edddd(i,j,i,j)
                f(ia,ja)=f(ia,ja)+fac
                f(ja,ia)=f(ja,ia)+fac
                f(ia,ia)=f(ia,ia)-.5_dp*p(ja,ja)*edddd(i,j,i,j)
                f(ja,ja)=f(ja,ja)-.5_dp*p(ia,ia)*edddd(i,j,i,j)
              else if ((i.lt.k).or.((i.eq.k).and.(j.lt.l))) then
                fac=2._dp*p(la,ka) * edddd(i,j,k,l)
                F(ia,ja) = F(ia,ja) + fac
                F(ja,ia) = F(ja,ia) + fac
                fac=2._dp*p(ja,ia) * edddd(i,j,k,l)
                F(ka,la) = F(ka,la) + fac
                F(la,ka) = F(la,ka) + fac
                fac=.5_dp*p(ja,ka) * edddd(i,j,k,l)
                F(ia,la) = F(ia,la) - fac
                F(la,ia) = F(la,ia) - fac
                fac=.5_dp*p(ja,la) * edddd(i,j,k,l)
                F(ia,ka) = F(ia,ka) - fac
                F(ka,ia) = F(ka,ia) - fac
                fac=.5_dp*p(ia,ka) * edddd(i,j,k,l)
                F(ja,la) = F(ja,la) - fac
                F(la,ja) = F(la,ja) - fac
                fac=.5_dp*p(ia,la) * edddd(i,j,k,l)
                F(ja,ka) = F(ja,ka) - fac
                F(ka,ja) = F(ka,ja) - fac
              end if
            !end if
          end do
        end do
      end do
    end do
  end if

end do ! the 'do' over dd shell pairs
if (ndd.gt.0) then
call index_sort(ndd,schwdd,orderdd)
call find_threshold(ndd,schwdd,orderdd,shthresh,lastdd)
end if

!return

do ii=1,ndl
  sh1=dlprs(1,ii)
  sh2=dlprs(2,ii)
  call i2e_dldl(sh1,sh2,sh1,sh2,edldl)
  schwdl(ii)=sqrt(maxval(abs(edldl)))
  ! The integrals of type (pilj|pilj) are (ab|ab), whereas all others are (ab|ac) or (ab|cd)
  do i=1,5
    ia=aosh(sh1)+i-1
    do j=1,4
      ja=aosh(sh2)+j-1
      if (abs(edldl(i,j,i,j)).gt.intthresh) then   ! the cases (pilj|pilj) are abab
        fac=1.5_dp*p(ia,ja)*edldl(i,j,i,j)
        f(ja,ia)=f(ja,ia)+fac
        f(ia,ja)=f(ia,ja)+fac
        f(ia,ia)=f(ia,ia)-.5_dp*p(ja,ja)*edldl(i,j,i,j)
        f(ja,ja)=f(ja,ja)-.5_dp*p(ia,ia)*edldl(i,j,i,j)
      end if
      do k=1,5
        ka=aosh(sh1)+k-1
        do l=1,4
          la=aosh(sh2)+l-1
          if ((abs(edldl(i,j,k,l)).gt.intthresh).and. &
               &   ((i.lt.k).or.((i.eq.k).and.(j.lt.l)))) then   ! the cases (spi|spj) , i.ne.j, are abac
            ka=aosh(sh1)+k-1; la=aosh(sh2)+l-1
            fac=2._dp*p(la,ka) * edldl(i,j,k,l)
            F(ia,ja) = F(ia,ja) + fac
            F(ja,ia) = F(ja,ia) + fac
            fac=2._dp*p(ja,ia) * edldl(i,j,k,l)
            F(ka,la) = F(ka,la) + fac
            F(la,ka) = F(la,ka) + fac
            fac=.5_dp*p(ja,ka) * edldl(i,j,k,l)
            F(ia,la) = F(ia,la) - fac
            F(la,ia) = F(la,ia) - fac
            fac=.5_dp*p(ja,la) * edldl(i,j,k,l)
            F(ia,ka) = F(ia,ka) - fac
            F(ka,ia) = F(ka,ia) - fac
            fac=.5_dp*p(ia,ka) * edldl(i,j,k,l)
            F(ja,la) = F(ja,la) - fac
            F(la,ja) = F(la,ja) - fac
            fac=.5_dp*p(ia,la) * edldl(i,j,k,l)
            F(ja,ka) = F(ja,ka) - fac
            F(ka,ja) = F(ka,ja) - fac
          end if
        end do
      end do
    end do
  end do
end do
if (ndl.gt.0) then
call index_sort(ndl,schwdl,orderdl)
call find_threshold(ndl,schwdl,orderdl,shthresh,lastdl)
end if


do ii=1,nll
  sh1=llprs(1,ii)
  sh2=llprs(2,ii)
  call i2e_llll(sh1,sh2,sh1,sh2,ellll)
  schwll(ii)=sqrt(maxval(abs(ellll)))
  ! Lots of different types among these 81 integrals, some of which are non-unique...
  ! Two main cases: the shells are the same, or they are not. 
  ! Start with same shell.
  if (sh1.eq.sh2) then  ! MAjor if statement - either sh1.eq.sh2 or not

    ! Type 1: there is one (ss|ss) and three (pipi|pipi) which are (aa|aa)
    do i=1,4
      ia=aosh(sh1)+i-1
      do j=1,4
        ja=aosh(sh1)+j-1
        do k=1,4
          ka=aosh(sh1)+k-1
          do l=1,4
            la=aosh(sh1)+l-1
            if (abs(ellll(i,j,k,l)).gt.intthresh) then 
              if ((i.eq.j).and.(j.eq.k).and.(k.eq.l)) then  ! (aa|aa)
                f(ia,ia)=f(ia,ia)+.5_dp*p(ia,ia)*ellll(i,i,i,i)
              else if ((i.eq.j).and.(k.eq.l).and.(j.lt.k)) then  ! (aa|bb)
                f(ia,ia)=f(ia,ia)+p(ka,ka)*ellll(i,i,k,k)
                fac=.5_dp*p(ia,ka)*ellll(i,i,k,k)
                f(ia,ka)=f(ia,ka)-fac
                f(ka,ia)=f(ka,ia)-fac
                f(ka,ka)=f(ka,ka)+p(ia,ia)*ellll(i,i,k,k)
              else if ((i.eq.j).and.(k.eq.l).and.(j.gt.k)) then   ! (bb|aa)
                cycle
              else if ((i.eq.k).and.(j.eq.l).and.(i.lt.j)) then   ! (ab|ab)
                fac=1.5_dp*p(ia,ja)*ellll(i,j,i,j)
                f(ia,ja)=f(ia,ja)+fac
                f(ja,ia)=f(ja,ia)+fac
                f(ia,ia)=f(ia,ia)-.5_dp*p(ja,ja)*ellll(i,j,i,j)
                f(ja,ja)=f(ja,ja)-.5_dp*p(ia,ia)*ellll(i,j,i,j)
              else if ((i.eq.l).and.(j.eq.k).and.(i.lt.j)) then
                cycle
              else if ((i.eq.l).and.(j.eq.k).and.(i.gt.j)) then
                cycle
              else if ((i.eq.k).and.(j.eq.l).and.(i.gt.j)) then
                cycle
              else if ((i.eq.j).and.(k.lt.l)) then       !  (aa|bc)
                f(ia,ia)=f(ia,ia)+2._dp*p(ka,la)*ellll(i,i,k,l)
                fac=.5_dp*p(ia,ka)*ellll(i,i,k,l)
                f(ia,la)=f(ia,la)-fac
                f(la,ia)=f(la,ia)-fac
                fac=.5_dp*p(ia,la)*ellll(i,i,k,l)
                f(ia,ka)=f(ia,ka)-fac
                f(ka,ia)=f(ka,ia)-fac
                fac=p(ia,ia)*ellll(i,i,k,l)
                f(ka,la)=f(ka,la)+fac
                f(la,ka)=f(la,ka)+fac
              else if ((i.eq.j).and.(k.gt.l)) then
                cycle
              else if ((i.lt.j).and.(k.eq.l)) then        ! (ab|cc)
                cycle
!                f(ka,ka)=f(ka,ka)+2._dp*p(ia,ja)*ellll(i,j,k,l)
!                fac=.5_dp*p(ka,ia)*ellll(i,j,k,l)
!                f(ka,ja)=f(ka,ja)-fac
!                f(ja,ka)=f(ja,ka)-fac
!                fac=.5_dp*p(ka,ja)*ellll(i,j,k,l)
!                f(ka,ia)=f(ka,ia)-fac
!                f(ia,ka)=f(ia,ka)-fac
!                fac=p(ka,ka)*ellll(i,j,k,l)
!                f(ia,ja)=f(ia,ja)+fac
!                f(ja,ia)=f(ja,ia)+fac
              else if ((i.gt.j).and.(k.eq.l)) then
                cycle
              else    ! Other cases which are 8-fold degenerate - just do the normal thing
                f(ia,ja)=f(ia,ja)+p(ka,la)*ellll(i,j,k,l)
                f(ia,la)=f(ia,la)-.5_dp*p(ja,ka)*ellll(i,j,k,l)
              end if
            end if
          end do
        end do
      end do
    end do
  ! Second case: the two shells are not the same
  ! This is simpler because now the highest permutation symm is (ab|ab). 16 of these.
  ! Still need to be careful, though, to only use each permutationally identical integral once
  else  ! This is the major if
    do i=1,4
      ia=aosh(sh1)+i-1
      do j=1,4
        ja=aosh(sh2)+j-1
        do k=1,4
          ka=aosh(sh1)+k-1
          do l=1,4
            la=aosh(sh2)+l-1
            if (abs(ellll(i,j,k,l)).gt.intthresh) then 
              if ((i.eq.k).and.(j.eq.l)) then   ! (ab|ab)
                fac=1.5_dp*p(ia,ja)*ellll(i,j,i,j)
                f(ia,ja)=f(ia,ja)+fac
                f(ja,ia)=f(ja,ia)+fac
                f(ia,ia)=f(ia,ia)-.5_dp*p(ja,ja)*ellll(i,j,i,j)
                f(ja,ja)=f(ja,ja)-.5_dp*p(ia,ia)*ellll(i,j,i,j)
              else if ((i.lt.k).or.((i.eq.k).and.(j.lt.l))) then
                fac=2._dp*p(la,ka) * ellll(i,j,k,l)
                F(ia,ja) = F(ia,ja) + fac
                F(ja,ia) = F(ja,ia) + fac
                fac=2._dp*p(ja,ia) * ellll(i,j,k,l)
                F(ka,la) = F(ka,la) + fac
                F(la,ka) = F(la,ka) + fac
                fac=.5_dp*p(ja,ka) * ellll(i,j,k,l)
                F(ia,la) = F(ia,la) - fac
                F(la,ia) = F(la,ia) - fac
                fac=.5_dp*p(ja,la) * ellll(i,j,k,l)
                F(ia,ka) = F(ia,ka) - fac
                F(ka,ia) = F(ka,ia) - fac
                fac=.5_dp*p(ia,ka) * ellll(i,j,k,l)
                F(ja,la) = F(ja,la) - fac
                F(la,ja) = F(la,ja) - fac
                fac=.5_dp*p(ia,la) * ellll(i,j,k,l)
                F(ja,ka) = F(ja,ka) - fac
                F(ka,ja) = F(ka,ja) - fac
              end if
            end if
          end do
        end do
      end do
    end do
  end if

end do ! the 'do' over ll shell pairs
if (nll.gt.0) then
call index_sort(nll,schwll,orderll)
call find_threshold(nll,schwll,orderll,shthresh,lastll)
end if

! Still needed: pd, dd, and dl shells. But those are for later

return
end subroutine fock_build_initial_diag

