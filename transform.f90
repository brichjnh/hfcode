subroutine transform1e()
use nrtype; use molprops
implicit none

real(dp) :: hcore(nb,nb)

allocate (m1e(nb,nb))
hcore=t1e+v1e

m1e=0._dp
m1e=matmul(cij,matmul(hcore,transpose(cij)))

write (9,*) "Transformed 1e integrals"
!do i=1,nb
!write(9,'(5F12.6)') m1e(i,:)
!end do

end subroutine transform1e


subroutine transform2e()
use nrtype; use molprops
implicit none

! This transforms the 2-electron integrals from the Crawford-like 'ijkl' one-dimensional array to a 4-dimensional
!  array of 2-electron integrals over spin-molecular orbitals.
! It is done in several steps!!
! And is currently inefficient in both memory and speed.

real(dp), allocatable :: tmat1(:,:,:,:), tmat2(:,:,:,:), tmat3(:,:,:,:)
integer(i4b) :: i,j,k,l,mu,nu,lbd,sig, ijkl, ij, kl, ia,ja,ka,la,allocok

! First allocate the array i2e
allocate (i2e4i(nb,nb,nb,nb),STAT=allocok)
if (allocok.ne.0) then
  write (9,*) "Problem allocating memory for AO basis 2-electron integrals in Transform"
  STOP
end if

i2e4i=0._dp
kl=1; ij=1
do ijkl=1,nb4
   ia=revind2(1,ij) ; ja=revind2(2,ij)
   ka=revind2(1,kl) ; la=revind2(2,kl)
   kl=kl+1
   if (kl.gt.ij) then
      kl=1
      ij=ij+1
   end if
   i2e4i(ia,ja,ka,la)=i2e(ijkl)
   i2e4i(ja,ia,ka,la)=i2e(ijkl)
   i2e4i(ia,ja,la,ka)=i2e(ijkl)
   i2e4i(ja,ia,la,ka)=i2e(ijkl)
   i2e4i(ka,la,ia,ja)=i2e(ijkl)
   i2e4i(la,ka,ia,ja)=i2e(ijkl)
   i2e4i(ka,la,ja,ia)=i2e(ijkl)
   i2e4i(la,ka,ja,ia)=i2e(ijkl)
end do

! Once we have i2e4i, we no longer need the 1-d array of integrals
deallocate(i2e)

! Then allocate the MO-basis 4-index integral array m2e
allocate (m2e(nb,nb,nb,nb),STAT=allocok)
if (allocok.ne.0) then
  write (9,*) "Problem allocating memory for MO basis 2-electron integrals in Transform"
  STOP
end if

! Then allocate the temporary 4-index arrays for the transform
allocate (tmat1(nb,nb,nb,nb),STAT=allocok)
if (allocok.ne.0) then
  write (9,*) "Problem allocating memory for temporary array tmat1 in Transform"
  STOP
end if
allocate (tmat2(nb,nb,nb,nb),STAT=allocok)
if (allocok.ne.0) then
  write (9,*) "Problem allocating memory for temporary array tmat2 in Transform"
  STOP
end if
allocate (tmat3(nb,nb,nb,nb),STAT=allocok)
if (allocok.ne.0) then
  write (9,*) "Problem allocating memory for temporary array tmat3 in Transform"
  STOP
end if

tmat1=0._dp
tmat2=0._dp
tmat3=0._dp
m2e=0._dp
do i=1,nb
   do mu=1,nb
      tmat1(i,:,:,:)=tmat1(i,:,:,:)+cij(i,mu)*i2e4i(mu,:,:,:)
   end do
   do j=1,nb
      do nu=1,nb
         tmat2(i,j,:,:)=tmat2(i,j,:,:)+cij(j,nu)*tmat1(i,nu,:,:)
      end do
      do k=1,nb
         do lbd=1,nb
            tmat3(i,j,k,:) = tmat3(i,j,k,:)+cij(k,lbd)*tmat2(i,j,lbd,:)
         end do
         do l=1,nb
            do sig=1,nb
                m2e(i,j,k,l)=m2e(i,j,k,l)+cij(l,sig)*tmat3(i,j,k,sig)
            end do
         end do
      end do
   end do
end do

! Now we can deallocate the temporary arrays as well as the AO basis array

deallocate(i2e4i,tmat1,tmat2,tmat3)
write (9,*) ""
write (9,*) "Calculated Two electron transformed integrals"

end subroutine transform2e

subroutine calcmoenergies()
use nrtype ; use molprops
implicit none

integer(i4b) :: i, j

allocate(moen(nb))

!m2e=0._dp
do i=1,nb
   moen(i)=m1e(i,i) 
   do j=1,nocc
      moen(i)=moen(i)+2._dp*m2e(i,i,j,j)-m2e(i,j,i,j)
   end do
end do

end subroutine calcmoenergies
