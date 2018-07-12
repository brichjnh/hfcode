subroutine cmpt1e()
use nrtype; use molprops ; use one_electron_terms
implicit none

! This subroutine is the driver routine for calculating all the 1-electron integrals
! It calls a variety of routines that calculate S, T or V matrix elements

real(dp) :: sss,tss,vss, ssp(3),tsp(3),vsp(3)
!real(dp) :: tmpmat(7,7)
real(dp) :: ssl(4),tsl(4),vsl(4)
real(dp) :: sll(4,4),tll(4,4),vll(4,4)
real(dp) :: spl(3,4),tpl(3,4),vpl(3,4)
real(dp) :: spp(3,3),tpp(3,3),vpp(3,3)
real(dp) :: ssd(5),tsd(5),vsd(5)
real(dp) :: spd(3,5),tpd(3,5),vpd(3,5)
real(dp) :: sdd(5,5),tdd(5,5),vdd(5,5)
real(dp) :: sdl(5,4),tdl(5,4),vdl(5,4)
integer(i4b) :: sh1, sh2, ao1, ao2
integer(i4b) :: i, ii,jj,ab

! For convenience, we will store the 1-electron overlap, kinetic and nuclear
!   Coulomb matrices in their full nbasis * nbasis form
! At this stage, though, we're not going to do anything particularly clever to speed things up
! For example 1-centre and 2-centre matrix elements are treated exactly the same way.

allocate (s1e(nb,nb),t1e(nb,nb),v1e(nb,nb))
s1e = 0.d0
t1e = 0.d0
v1e = 0.d0

do ab=1,nss 
  sh1=ssprs(1,ab)
  ao1=aosh(sh1)
  sh2=ssprs(2,ab)
  ao2=aosh(sh2)
  call s_s_1e_terms(sh1,sh2,sss,tss,vss)
  s1e(ao1,ao2)=sss
  s1e(ao2,ao1)=sss
  t1e(ao1,ao2)=tss
  t1e(ao2,ao1)=tss
  v1e(ao1,ao2)=vss
  v1e(ao2,ao1)=vss
end do
do ab=1,nsp
  sh1=spprs(1,ab)
  ao1=aosh(sh1)
  sh2=spprs(2,ab)
  ao2=aosh(sh2)
  call s_p_1e_terms(sh1,sh2,ssp,tsp,vsp)
  do ii=1,3
    s1e(ao1,ao2+ii-1)=ssp(ii)
    s1e(ao2+ii-1,ao1)=ssp(ii)
    t1e(ao1,ao2+ii-1)=tsp(ii)
    t1e(ao2+ii-1,ao1)=tsp(ii)
    v1e(ao1,ao2+ii-1)=vsp(ii)
    v1e(ao2+ii-1,ao1)=vsp(ii)
  end do
end do
do ab=1,nsl
  sh1=slprs(1,ab)
  ao1=aosh(sh1)
  sh2=slprs(2,ab)
  ao2=aosh(sh2)
  call s_l_1e_terms(sh1,sh2,ssl,tsl,vsl)
  do ii=1,4
    s1e(ao1,ao2+ii-1)=ssl(ii)
    s1e(ao2+ii-1,ao1)=ssl(ii)
    t1e(ao1,ao2+ii-1)=tsl(ii)
    t1e(ao2+ii-1,ao1)=tsl(ii)
    v1e(ao1,ao2+ii-1)=vsl(ii)
    v1e(ao2+ii-1,ao1)=vsl(ii)
  end do
end do
do ab=1,npp
  sh1=ppprs(1,ab)
  ao1=aosh(sh1)
  sh2=ppprs(2,ab)
  ao2=aosh(sh2)
  call p_p_1e_terms(sh1,sh2,spp,tpp,vpp)
  do ii=1,3
    do jj=1,3
      s1e(ao1+ii-1,ao2+jj-1)=spp(ii,jj)
      s1e(ao2+jj-1,ao1+ii-1)=spp(ii,jj)
      t1e(ao1+ii-1,ao2+jj-1)=tpp(ii,jj)
      t1e(ao2+jj-1,ao1+ii-1)=tpp(ii,jj)
      v1e(ao1+ii-1,ao2+jj-1)=vpp(ii,jj)
      v1e(ao2+jj-1,ao1+ii-1)=vpp(ii,jj)
    end do
  end do
end do
do ab=1,npl
  sh1=plprs(1,ab)
  ao1=aosh(sh1)
  sh2=plprs(2,ab)
  ao2=aosh(sh2)
  call p_l_1e_terms(sh1,sh2,spl,tpl,vpl)
  do ii=1,3
    do jj=1,4
      s1e(ao1+ii-1,ao2+jj-1)=spl(ii,jj)
      s1e(ao2+jj-1,ao1+ii-1)=spl(ii,jj)
      t1e(ao1+ii-1,ao2+jj-1)=tpl(ii,jj)
      t1e(ao2+jj-1,ao1+ii-1)=tpl(ii,jj)
      v1e(ao1+ii-1,ao2+jj-1)=vpl(ii,jj)
      v1e(ao2+jj-1,ao1+ii-1)=vpl(ii,jj)
    end do
  end do
end do
do ab=1,nll
  sh1=llprs(1,ab)
  ao1=aosh(sh1)
  sh2=llprs(2,ab)
  ao2=aosh(sh2)
  call l_l_1e_terms(sh1,sh2,sll,tll,vll)
  do ii=1,4
    do jj=1,4
      s1e(ao1+ii-1,ao2+jj-1)=sll(ii,jj)
      s1e(ao2+jj-1,ao1+ii-1)=sll(ii,jj)
      t1e(ao1+ii-1,ao2+jj-1)=tll(ii,jj)
      t1e(ao2+jj-1,ao1+ii-1)=tll(ii,jj)
      v1e(ao1+ii-1,ao2+jj-1)=vll(ii,jj)
      v1e(ao2+jj-1,ao1+ii-1)=vll(ii,jj)
    end do
  end do
end do
do ab=1,nsd
  sh1=sdprs(1,ab)
  ao1=aosh(sh1)
  sh2=sdprs(2,ab)
  ao2=aosh(sh2)
  call s_d_1e_terms(sh1,sh2,ssd,tsd,vsd)
  do ii=1,5
    s1e(ao1,ao2+ii-1)=ssd(ii)
    s1e(ao2+ii-1,ao1)=ssd(ii)
    t1e(ao1,ao2+ii-1)=tsd(ii)
    t1e(ao2+ii-1,ao1)=tsd(ii)
    v1e(ao1,ao2+ii-1)=vsd(ii)
    v1e(ao2+ii-1,ao1)=vsd(ii)
  end do
end do
do ab=1,npd
  sh1=pdprs(1,ab)
  ao1=aosh(sh1)
  sh2=pdprs(2,ab)
  ao2=aosh(sh2)
  call p_d_1e_terms(sh1,sh2,spd,tpd,vpd)
  do ii=1,5
    do jj=1,3
      s1e(ao1+jj-1,ao2+ii-1)=spd(jj,ii)
      s1e(ao2+ii-1,ao1+jj-1)=spd(jj,ii)
      t1e(ao1+jj-1,ao2+ii-1)=tpd(jj,ii)
      t1e(ao2+ii-1,ao1+jj-1)=tpd(jj,ii)
      v1e(ao1+jj-1,ao2+ii-1)=vpd(jj,ii)
      v1e(ao2+ii-1,ao1+jj-1)=vpd(jj,ii)
    end do
  end do
end do
do ab=1,ndd
  sh1=ddprs(1,ab)
  ao1=aosh(sh1)
  sh2=ddprs(2,ab)
  ao2=aosh(sh2)
  call d_d_1e_terms(sh1,sh2,sdd,tdd,vdd)
  do ii=1,5
    do jj=1,5
      s1e(ao1+jj-1,ao2+ii-1)=sdd(jj,ii)
      s1e(ao2+ii-1,ao1+jj-1)=sdd(jj,ii)
      t1e(ao1+jj-1,ao2+ii-1)=tdd(jj,ii)
      t1e(ao2+ii-1,ao1+jj-1)=tdd(jj,ii)
      v1e(ao1+jj-1,ao2+ii-1)=vdd(jj,ii)
      v1e(ao2+ii-1,ao1+jj-1)=vdd(jj,ii)
    end do
  end do
end do
do ab=1,ndl
  sh1=dlprs(1,ab)
  ao1=aosh(sh1)
  sh2=dlprs(2,ab)
  ao2=aosh(sh2)
  call d_l_1e_terms(sh1,sh2,sdl,tdl,vdl)
  do ii=1,4
    do jj=1,5
      s1e(ao1+jj-1,ao2+ii-1)=sdl(jj,ii)
      s1e(ao2+ii-1,ao1+jj-1)=sdl(jj,ii)
      t1e(ao1+jj-1,ao2+ii-1)=tdl(jj,ii)
      t1e(ao2+ii-1,ao1+jj-1)=tdl(jj,ii)
      v1e(ao1+jj-1,ao2+ii-1)=vdl(jj,ii)
      v1e(ao2+ii-1,ao1+jj-1)=vdl(jj,ii)
    end do
  end do
end do

if (prtlevl.ge.3) then
  write (9,*) ""
  write (9,*) "One-electron matrix elements"
  write (9,*) "Overlap matrix elements"
  do i=1,nb
     if (nb.le.10) then
       write (9,'(I4,10F13.7)') i,s1e(i,:)
     else
       write (9,'(I4,10F13.7)') i,s1e(i,1:10)
       write (9,'(4x,10F13.7)') s1e(i,11:)
     end if
  end do
  write (9,*) ""
  write (9,*) "Kinetic matrix elements"
  do i=1,nb
     if (nb.le.10) then
       write (9,'(I4,10F13.7)') i,t1e(i,:)
     else
       write (9,'(I4,10F13.7)') i,t1e(i,1:10)
       write (9,'(4x,10F13.7)') t1e(i,11:)
     end if
  end do
  write (9,*) ""
  write (9,*) "Potential matrix elements"
  do i=1,nb
     if (nb.le.10) then
       write (9,'(I4,10F13.7)') i,v1e(i,:)
     else
       write (9,'(I4,10F13.7)') i,v1e(i,1:10)
       write (9,'(4x,10F13.7)') v1e(i,11:)
     end if
  end do
  write (9,*) "End of 1-e matrix elements"
end if

!oo(1:7)=(/1,2,6,7,3,4,5 /)
!do i=1,nb
!  do j=1,nb
!    tmpmat(i,j)=s1e(oo(i),oo(j))
!  end do
!end do
!  write (9,*) "reshuffled overlap matrix elements"
!  do i=1,nb
!     write (9,'(10F13.7)') tmpmat(i,:)
!  end do
!do i=1,nb
!  do j=1,nb
!    tmpmat(i,j)=t1e(oo(i),oo(j))
!  end do
!end do
!  write (9,*) "reshuffled kin matrix elements"
!  do i=1,nb
!     write (9,'(10F13.7)') tmpmat(i,:)
!  end do
!do i=1,nb
!  do j=1,nb
!    tmpmat(i,j)=v1e(oo(i),oo(j))
!  end do
!end do
!  write (9,*) "reshuffled pot matrix elements"
!  do i=1,nb
!     write (9,'(10F13.7)') tmpmat(i,:)
!  end do
!stop

end subroutine cmpt1e

