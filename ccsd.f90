subroutine ccsd()
use nrtype ; use molprops
implicit none

integer(i4b) :: nm, no, nv
integer(i4b) :: i,j,k,l,ii,jj,kk,ll,m,n,a,b,e,f
integer(i4b) :: iter
real(dp) :: fac1, fac2, fac3, mp2check, ecc
real(dp), allocatable :: tei(:,:,:,:), oe(:), t1(:,:), t1n(:,:)
real(dp), allocatable :: t2(:,:,:,:), t2n(:,:,:,:), tau(:,:,:,:)
real(dp), allocatable :: frdii(:,:), frdie(:,:), frdee(:,:)
real(dp), allocatable :: wrdii(:,:,:,:), wrdie(:,:,:,:), wrdee(:,:,:,:),t2wrdie(:,:,:,:)
real(dp), allocatable :: wrdee1d(:),wrdie1d(:)
real(dp), allocatable :: tfrd(:,:), tfrd2(:,:)
real(dp), allocatable :: teioooo(:), teiovvv(:), teioovv(:), teiooov(:), teivvvv(:),teiovvo(:)
integer(i4b), allocatable :: indoooo(:,:), indovvv(:,:), indoovv(:,:), indooov(:,:), indvvvv(:,:),indovvo(:,:)
integer(i4b), allocatable :: indwrdee(:,:), indwrdie(:,:)
integer(i4b) :: maxoooo, maxovvv, maxoovv, maxooov, maxvvvv,maxovvo
integer(i4b) :: noooo, novvv, noovv, nooov, nvvvv,novvo
integer(i4b) :: nwrdee,nwrdie, allocok
real(dp), parameter :: thresh=1.d-9

! In this new version, we will use the fact that the Fock matrix in the spin orbitals is diagonal.
! Start by defining some sizes and allocating some arrays

nm=2*(nb-ncore)  ! total number spin orbitals (without the core)
no=2*(nocc-ncore)  ! number occupied spin orbitals (without the core)
nv=nm-no       ! number vacant spin orbitals
allocate(oe(nm))
allocate(t1(no,nv),t2(no,no,nv,nv))
allocate(t1n(no,nv),t2n(no,no,nv,nv))
allocate(tau(no,no,nv,nv))
allocate(tfrd(nv,nv),tfrd2(no,no))
allocate(frdii(no,no),frdie(no,nv),frdee(nv,nv))
allocate(wrdii(no,no,no,no),wrdie(no,nv,nv,no),wrdee(nv,nv,nv,nv),t2wrdie(no,nv,nv,no))

! Then set up an array with orbital energies
oe(1:nm-1:2)=moen(ncore+1:nb)
oe(2:nm:2)=moen(ncore+1:nb)

! Then construct a 4-D array with all of the integrals <ij||kl> in terms of spin orbitals.
! The orbitals are taken to be a\psi_1, b\psi1, a\psi_2, ...
! <ij||kl>=<ij|kl>-<ij|lk>=[ik|jl]-[il|jk]
!  and [ik|jl] is only non-zero if both i and k and j and l have same spin
!  and [il|jk] is only non-zero if both i and l and j and k have same spin
! Note we only need these terms for non-core orbitals!!

allocate(tei(nm,nm,nm,nm),STAT=allocok)
if (allocok.ne.0) then
   write (9,*) "Problem allocating the spin-orbital array of TEIs in ccsd"
   stop
end if
tei=0._dp
do i=ncore+1,nb
   ii=2*(i-ncore)-1
   do j=ncore+1,nb
      jj=2*(j-ncore)-1
      do k=ncore+1,nb
         kk=2*(k-ncore)-1
         do l=ncore+1,nb
            ll=2*(l-ncore)-1
            fac1=m2e(i,k,j,l)
            fac2=m2e(i,l,j,k)
            tei(ii,jj,kk,ll)=fac1-fac2
            tei(ii,jj+1,kk,ll+1)=fac1
            tei(ii,jj+1,kk+1,ll)=-fac2
            tei(ii+1,jj,kk,ll+1)=-fac2
            tei(ii+1,jj,kk+1,ll)=fac1
            tei(ii+1,jj+1,kk+1,ll+1)=fac1-fac2
         end do
      end do
   end do
end do

! Once this is done, we no longer need m2e
deallocate(m2e)

! Now build bespoke 1-dimensional arrays with only the non-zero integrals
maxoooo=no**4/2
allocate(indoooo(4,maxoooo),teioooo(maxoooo))
noooo=0
do j=1,no
   do i=1,no
      do n=1,no
         do m=1,no
            fac1=tei(m,n,i,j)
            if (abs(fac1).gt.thresh) then
               noooo=noooo+1
               indoooo(:,noooo) = (/m,n,i,j/)
               teioooo(noooo) = fac1
            end if
         end do
      end do
   end do
end do
write (*,*) "number of oooo non-zero ints:",noooo

maxovvv=no*nv**3/2
allocate(indovvv(4,maxovvv),teiovvv(maxovvv))
novvv=0
do e=1,nv
   do f=1,nv
      do a=1,nv
         do m=1,no
            fac1=tei(m,a+no,f+no,e+no)
            if (abs(fac1).gt.thresh) then
               novvv=novvv+1
               indovvv(:,novvv) = (/m,a,f,e/)
               teiovvv(novvv) = fac1
            end if
         end do
      end do
   end do
end do
write (*,*) "number of ovvv non-zero ints:",novvv

maxoovv=no**2*nv**2/2
allocate(indoovv(4,maxoovv),teioovv(maxoovv))
noovv=0
do f=1,nv
   do e=1,nv
      do n=1,no
         do m=1,no
            fac1=tei(m,n,e+no,f+no)
            if (abs(fac1).gt.thresh) then
               noovv=noovv+1
               indoovv(:,noovv) = (/m,n,e,f/)
               teioovv(noovv) = fac1
            end if
         end do
      end do
   end do
end do
write (*,*) "number of oovv non-zero ints:",noovv

maxovvo=no**2*nv**2
allocate(indovvo(4,maxovvo),teiovvo(maxovvo))
novvo=0
do j=1,no
   do e=1,nv
      do b=1,nv
         do m=1,no
            fac1=tei(m,b+no,e+no,j)
            if (abs(fac1).gt.thresh) then
               novvo=novvo+1
               indovvo(:,novvo) = (/m,b,e,j/)
               teiovvo(novvo) = fac1
            end if
         end do
      end do
   end do
end do
write (*,*) "number of ovvo non-zero ints:",novvo

maxooov=no**3*nv/2
allocate(indooov(4,maxooov),teiooov(maxooov))
nooov=0
do e=1,nv
   do i=1,no
      do n=1,no
         do m=1,no
            fac1=tei(m,n,i,e+no)
            if (abs(fac1).gt.thresh) then
               nooov=nooov+1
               indooov(:,nooov) = (/m,n,i,e/)
               teiooov(nooov) = fac1
            end if
         end do
      end do
   end do
end do
write (*,*) "number of ooov non-zero ints:",nooov

maxvvvv=nv**4
allocate(indvvvv(4,maxvvvv),teivvvv(maxvvvv))
nvvvv=0
do f=1,nv
   do e=1,nv
      do b=1,nv
         do a=1,nv
            fac1=tei(a+no,b+no,e+no,f+no)
            if (abs(fac1).gt.thresh) then
               nvvvv=nvvvv+1
               indvvvv(:,nvvvv) = (/a,b,e,f/)
               teivvvv(nvvvv) = fac1
            end if
         end do
      end do
   end do
end do
write (*,*) "number of vvvv non-zero ints:",nvvvv

! Now assign initial amplitudes, and check mp2 energy

t1=0._dp
mp2check=0._dp
do b=1,nv
   do a=1,nv
      do j=1,no
         do i=1,no
            fac1=tei(i,j,a+no,b+no)/(oe(i)+oe(j)-oe(a+no)-oe(b+no))
            mp2check=mp2check+.25_dp*fac1*tei(i,j,a+no,b+no)
            t2(i,j,a,b)=fac1
         end do
      end do
   end do
end do

! And now the whole array tei is no longer needed
deallocate(tei)
allocate(indwrdee(4,maxvvvv),wrdee1d(maxvvvv))
allocate(wrdie1d(maxovvo),indwrdie(4,maxovvo))

write (9,'(A,2F18.10)') "MP2 energy in CC code:",mp2check,mp2check+etot

! Now start CC proper. In a loop in due course, just one iteration for the moment.

write (9,*) "CCSD iterations. Cycle number, ecorrel, etot:"
do iter=1,20

!  Start by building tau tilde, eq. (9), needed in eqs. (3) and (4)
   do concurrent (i=1:no,j=1:no,a=1:nv,b=1:nv)
      tau(i,j,a,b)=t2(i,j,a,b)+.5_dp*(t1(i,a)*t1(j,b)-t1(i,b)*t1(j,a))
   end do

! Assemble frdee following eq. (3) using compct integrals
   frdee=0._dp
   do ii=1,novvv
      m=indovvv(1,ii)
      a=indovvv(2,ii)
      f=indovvv(3,ii)
      e=indovvv(4,ii)
      frdee(a,e)=frdee(a,e)+t1(m,f)*teiovvv(ii)
   end do
   do ii=1,noovv
      m=indoovv(1,ii)
      n=indoovv(2,ii)
      e=indoovv(3,ii)
      f=indoovv(4,ii)
      frdee(:,e)=frdee(:,e)-.5_dp*tau(m,n,:,f)*teioovv(ii)
   end do

! Now the same for eq. (4). Here too the first two terms vanish
   frdii=0._dp
   do ii=1,nooov
      m=indooov(1,ii)
      n=indooov(2,ii)
      i=indooov(3,ii)
      e=indooov(4,ii)
      frdii(m,i)=frdii(m,i)+t1(n,e)*teiooov(ii)
   end do
   do ii=1,noovv
      m=indoovv(1,ii)
      n=indoovv(2,ii)
      e=indoovv(3,ii)
      f=indoovv(4,ii)
      frdii(m,:)=frdii(m,:)+.5_dp*tau(:,n,e,f)*teioovv(ii)
   end do
 
! Now eq (5)
   frdie=0._dp
   do ii=1,noovv
      m=indoovv(1,ii)
      n=indoovv(2,ii)
      e=indoovv(3,ii)
      f=indoovv(4,ii)
      frdie(m,e)=frdie(m,e)+t1(n,f)*teioovv(ii)
   end do

!  Now buid tau, eq. (10), needed in eqs. (6), (7) and (2)
   do b=1,nv
      do a=1,nv
         do j=1,no
            do i=1,no
               tau(i,j,a,b)=t2(i,j,a,b)+t1(i,a)*t1(j,b)-t1(i,b)*t1(j,a)
            end do
         end do
      end do
   end do

! eq (6) using compact approach - in the last term, the factor 1/4 will be swithced to 1/2 as explained in appendix
   wrdii=0._dp
   do ii=1,noooo
      m=indoooo(1,ii)
      n=indoooo(2,ii)
      i=indoooo(3,ii)
      j=indoooo(4,ii)
      wrdii(m,n,i,j)=teioooo(ii)
   end do
   do ii=1,nooov
      m=indooov(1,ii)
      n=indooov(2,ii)
      i=indooov(3,ii)
      e=indooov(4,ii)
      do j=1,no
         wrdii(m,n,i,j)=wrdii(m,n,i,j)+t1(j,e)*teiooov(ii)
         wrdii(m,n,j,i)=wrdii(m,n,j,i)-t1(j,e)*teiooov(ii)
      end do
   end do
   do i=1,no
      do j=1,no
      do ii=1,noovv
         m=indoovv(1,ii)
         n=indoovv(2,ii)
         e=indoovv(3,ii)
         f=indoovv(4,ii)
            wrdii(m,n,i,j)=wrdii(m,n,i,j)+.5_dp*tau(i,j,e,f)*teioovv(ii)
         end do
      end do
   end do

! eq (7) using compact approach - the last term will be dropped as explained in appendix
! Also note that we use here the fact that <am||ef> = <ma||fe>
   wrdee=0._dp
   do ii=1,nvvvv
      a=indvvvv(1,ii)
      b=indvvvv(2,ii)
      e=indvvvv(3,ii)
      f=indvvvv(4,ii)
      wrdee(a,b,e,f)=teivvvv(ii)
   end do
   do ii=1,novvv
      m=indovvv(1,ii)
      a=indovvv(2,ii)
      f=indovvv(3,ii)
      e=indovvv(4,ii)
      do b=1,nv
         wrdee(a,b,e,f)=wrdee(a,b,e,f)-t1(m,b)*teiovvv(ii)
         wrdee(b,a,e,f)=wrdee(b,a,e,f)+t1(m,b)*teiovvv(ii)
      end do
   end do
! Then rearrange to 1d array and keep only the non-zero ones
   nwrdee=0
   do f=1,nv
      do e=1,nv
         do b=1,nv
            do a=1,nv
               fac1=wrdee(a,b,e,f)
               if (abs(fac1).gt.thresh) then
                  nwrdee=nwrdee+1
                  indwrdee(:,nwrdee) = (/a,b,e,f/)
                  wrdee1d(nwrdee) = fac1
               end if
            end do
         end do
      end do
   end do


! eq (8) compact
!   Here we use the fact that <mn||ej> = -<mn||je>
   wrdie=0._dp
   do ii=1,novvo
      m=indovvo(1,ii)
      b=indovvo(2,ii)
      e=indovvo(3,ii)
      j=indovvo(4,ii)
      wrdie(m,b,e,j)=teiovvo(ii)
   end do
   do ii=1,novvv
      m=indovvv(1,ii)
      b=indovvv(2,ii)
      e=indovvv(3,ii)
      f=indovvv(4,ii)
      do j=1,no
         wrdie(m,b,e,j)=wrdie(m,b,e,j)+t1(j,f)*teiovvv(ii)
      end do
   end do
   do ii=1,nooov
      m=indooov(1,ii)
      n=indooov(2,ii)
      j=indooov(3,ii)
      e=indooov(4,ii)
      do b=1,nv
         wrdie(m,b,e,j)=wrdie(m,b,e,j)+t1(n,b)*teiooov(ii)
      end do
   end do
   do ii=1,noovv
      m=indoovv(1,ii)
      n=indoovv(2,ii)
      e=indoovv(3,ii)
      f=indoovv(4,ii)
      do b=1,nv
         do j=1,no
            wrdie(m,b,e,j)=wrdie(m,b,e,j)-(.5_dp*t2(j,n,f,b)+t1(j,f)*t1(n,b))*teioovv(ii)
         end do
      end do
   end do
!  Again we will re-sort these integrals so as only to keep the non-zero ones
   nwrdie=0
   do m=1,no
      do b=1,nv
         do e=1,nv
            do j=1,no
               fac1=wrdie(m,b,e,j)
               if (abs(fac1).gt.thresh) then
                  nwrdie=nwrdie+1
                  indwrdie(:,nwrdie) = (/m,b,e,j/)
                  wrdie1d(nwrdie) = fac1
               end if
            end do
         end do
      end do
   end do

! Now update the amplitudes, starting with singles, eq. (1)
! Again we are going to use the sub-blocks of integrals in 1-d arrays, requiring some re-writing.
! Note that again we can ignore the off-diagonal Fock matrix element, the first term
   do concurrent (i=1:no,a=1:nv)
      t1n(i,a)=sum(t1(i,:)*frdee(a,:))-sum(t1(:,a)*frdii(:,i))+sum(t2(i,:,a,:)*frdie(:,:))
   end do
!  For the 4th term, we need to address TEI's of the type <na||if>>.
!    These are equal to -<na||fi> which are homologous to the ovvo's <mb||ej> we already have
   do ii=1,novvo
      n=indovvo(1,ii)
      a=indovvo(2,ii)
      f=indovvo(3,ii)
      i=indovvo(4,ii)
      t1n(i,a)=t1n(i,a)+t1(n,f)*teiovvo(ii)
   end do
!  For the 5th term, we need to address TEI's of the type <ma||ef>>.
!     These are just ovvv's that we used already (albeit with different labels)
   do ii=1,novvv
      m=indovvv(1,ii)
      a=indovvv(2,ii)
      e=indovvv(3,ii)
      f=indovvv(4,ii)
      do i=1,no
         t1n(i,a)=t1n(i,a)-.5_dp*t2(i,m,e,f)*teiovvv(ii)
      end do
   end do
!  For the 6th term, we need to address TEI's of the type <nm||ei>>.
!      these are the same as the ooov's <mn||ie> used earlier, as <nm||ei>=<mn||ie>
   do ii=1,nooov
      m=indooov(1,ii)
      n=indooov(2,ii)
      i=indooov(3,ii)
      e=indooov(4,ii)
      do a=1,nv
         t1n(i,a)=t1n(i,a)-.5_dp*t2(m,n,a,e)*teiooov(ii)
      end do
   end do
   do concurrent (i=1:no,a=1:nv)
      t1n(i,a)=t1n(i,a)/(oe(i)-oe(a+no))
   end do

! Now update the doubles amplitudes, eq. (2)
! First get some sums that are re-used for multiple i, j, a or b
!  The t1 amplitudes are often zipped in a particular waywith F_round so zip them up
   do e=1,nv
      do b=1,nv
         tfrd(b,e)=frdee(b,e)-.5_dp*sum(t1(:,b)*frdie(:,e))
      end do
   end do
   do j=1,no
      do m=1,no
         tfrd2(m,j)=frdii(m,j)+.5_dp*sum(t1(j,:)*frdie(m,:))
      end do
   end do
if ((iter.eq.1).or.(iter.eq.2)) then
end if
!  First zero all items and copy the OOVV terms
   t2n=0._dp
   do ii=1,noovv
      i=indoovv(1,ii)
      j=indoovv(2,ii)
      a=indoovv(3,ii)
      b=indoovv(4,ii)
      t2n(i,j,a,b) =teioovv(ii)
   end do

!  This does the 5th term
   do ii=1,nwrdee
      a=indwrdee(1,ii)
      b=indwrdee(2,ii)
      e=indwrdee(3,ii)
      f=indwrdee(4,ii)
      do j=1,no
         do i=1,no
            t2n(i,j,a,b)=t2n(i,j,a,b)+.5_dp*tau(i,j,e,f)*wrdee1d(ii)
         end do
      end do
   end do

! The sixth term needs to be done in two goes and is stored in a temp array
   t2wrdie=0._dp
   do ii=1,nwrdie
      m=indwrdie(1,ii)
      b=indwrdie(2,ii)
      e=indwrdie(3,ii)
      j=indwrdie(4,ii)
      do i=1,no
         do a=1,nv
            t2wrdie(i,a,b,j)=t2wrdie(i,a,b,j)+t2(i,m,a,e)*wrdie1d(ii)
         end do
      end do
   end do
   do ii=1,novvo
      m=indovvo(1,ii)
      b=indovvo(2,ii)
      e=indovvo(3,ii)
      j=indovvo(4,ii)
      do i=1,no
         do a=1,nv
            t2wrdie(i,a,b,j)=t2wrdie(i,a,b,j)-t1(i,e)*t1(m,a)*teiovvo(ii)
         end do
      end do
   end do

! The seventh term requires TEIs <ab||ej> which are of the form VVVO
!   We already have ovvv <ma||fe> = <ef||am>
   do ii=1,novvv
      a=indovvv(4,ii)
      b=indovvv(3,ii)
      e=indovvv(2,ii)
      j=indovvv(1,ii)
      do i=1,no
         t2n(i,j,a,b)=t2n(i,j,a,b)+t1(i,e)*teiovvv(ii)
         t2n(j,i,a,b)=t2n(j,i,a,b)-t1(i,e)*teiovvv(ii)
      end do
   end do

! The eighth and final term requires TEIs <mb||ij> which are of the form ovoo
!   We already have ooov <mn||ie> = <ie||mn>
   do a=1,nv
      do ii=1,nooov
         m=indooov(3,ii)
         b=indooov(4,ii)
         i=indooov(1,ii)
         j=indooov(2,ii)
         fac1=t1(m,a)*teiooov(ii)
         t2n(i,j,a,b)=t2n(i,j,a,b)-fac1
         t2n(i,j,b,a)=t2n(i,j,b,a)+fac1
      end do
   end do

!  Then spread all the terms that depend on low-dimension arrays
!  These are the second, and third terms
!  The fourth term depends on six indices, but is O4V2, so relatively low, and will also be done the old way
!  Finally add the accumulated sixth terms, and divide by orbital energy difference
   do b=1,nv
      do a=1,nv
         do j=1,no
            do i=1,no
               fac1=sum(t2(i,j,a,:)*tfrd(b,:))-sum(t2(i,j,b,:)*tfrd(a,:)) &
                & -(sum(t2(i,:,a,b)*tfrd2(:,j))-sum(t2(j,:,a,b)*tfrd2(:,i)))
               fac2=sum(tau(:,:,a,b)*wrdii(:,:,i,j))
               fac3=t2wrdie(i,a,b,j)-t2wrdie(i,b,a,j)-t2wrdie(j,a,b,i)+t2wrdie(j,b,a,i)
               t2n(i,j,a,b)=t2n(i,j,a,b)+fac1+.5_dp*fac2+fac3
               t2n(i,j,a,b)=t2n(i,j,a,b)/(oe(i)+oe(j)-oe(a+no)-oe(b+no))
            end do
         end do
      end do
   end do

   t1=t1n
   t2=t2n

! Now compute the CCSD energy
   ecc=0._dp
   do ii=1,noovv
      i=indoovv(1,ii)
      j=indoovv(2,ii)
      a=indoovv(3,ii)
      b=indoovv(4,ii)
      fac1=teioovv(ii)*(.25_dp*t2(i,j,a,b)+.5_dp*t1(i,a)*t1(j,b))
      ecc=ecc+fac1
   end do

   write (9,'(I3,2F18.12)') iter,ecc,ecc+etot
end do
deallocate(oe,t1,t2,t1n,t2n,tau)
deallocate(frdii,frdie,frdee,wrdii,wrdie,wrdee)
deallocate(tfrd,tfrd2,t2wrdie)
deallocate(indoooo,teioooo,indovvv,teiovvv,indoovv,teioovv)
deallocate(indovvo,teiovvo,indooov,teiooov,indvvvv,teivvvv)
deallocate(indwrdee,wrdee1d,wrdie1d,indwrdie)
call time_checker(-1,"ccsd iteratiions finished , timing:")

end subroutine ccsd

