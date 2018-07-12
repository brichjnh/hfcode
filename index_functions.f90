subroutine buildindices()
use nrtype ; use molprops
implicit none

! This subroutine builds several index arrays, in particular 'minus1' as well as some related one.
! For a given index(column) and row(Cartesian direction), the array contains the index
! number for the corresponding subtraction of one exponent. The null function has index 0, s has 1, and so on
! Currently the array "minus1" goes up to 165 (which is a 'k' function, z**8 - this is needed for
!  building (dd|dd) integrals when using the Giese/Hamilton/Schaefer approach.

! The array minus1s is a smaller version of minus 1, it only goes up to element 45 (z**4, enough to build all (pp|pp)
! The array 'plus1' is the reverse: for a given function, it tells you the index of the f. you get if you *add* 1 to 
!   the x, y or z index.

! Some additional arrays are used here to build these 4 arrays

! Indexps is the one important such array. It is 2-dimensional:
!   in one dimension, it contains three exponents, the other gives the sequence order of that exponent set
! Example of an entry: the row (:,10) contains the elements (1,10)=0, (2,10)=0, (3,10)=2. This is because exponent
!  x^0 y^0 z^2 is the 10th exponent. Likewise, 

! The other array is facind, which depends on three indices, impolicitly the x, y and z exponent. It comntains
! the index number for the appropriate entry. E.g. facind(0,0,2)=10

! Another array made here is one that maps integrals to storage of different permutation symmetry (aa|aa), (aa|ab), (aa|bb), (ab|ab), (aa|bc) and (ab|cd)

integer(i4b) :: i, j, l, m, n, cter, ll, lm1, mm1, nm1, ij, kl, im1, jj
integer(i4b) :: facind(0:labcdmax,0:labcdmax,0:labcdmax)

! The offset function, needs to be calculated up to nb*(nb+1)/2
allocate(ioff(nb2))
ioff(1)=0
do i=2,nb2
  ioff(i)=ioff(i-1)+i-1
end do

! Also the 2-shell indices a la Crawford made here, and its reverse
allocate(ind2(nb,nb))
do i=1,nb
  do j=1,i
    ind2(i,j)=ioff(i)+j
    ind2(j,i)=ind2(i,j)
  end do
end do
allocate (revind2(2,nb2))
! Set up an array to step back from ij or kl to i and j
ij=1
kl=1
do i=1,nb*(nb+1)/2
  revind2(:,i)=(/ij,i-ioff(ij)/)
  kl=kl+1
  if (kl.gt.ij) then
    kl=1
    ij=ij+1
  end if
end do


! Now start on minus 1 and co.
cter=1
indexps(:,1) = (/0,0,0/)
facind(0,0,0)=1
do ll=1,labcdmax
  do i=0,ll
    l=ll-i
    do n=0, i
      m=i-n
      cter=cter+1
      facind(l,m,n)=cter
      indexps(:,cter)=(/l,m,n/)
    end do
  end do
end do
plus1(:,0)=0
minus1(:,0)=0

! Now work out the plus1 and minus1 addresses
do i=1,cter
  l=indexps(1,i)
  m=indexps(2,i)
  n=indexps(3,i)
  if (i.le.plu1max) then
    plus1(1,i)=facind(l+1,m,n)
    plus1(2,i)=facind(l,m+1,n)
    plus1(3,i)=facind(l,m,n+1)
  end if
  lm1=l-1
  if (lm1.lt.0) then
    minus1(1,i)=0
  else
    minus1(1,i)=facind(lm1,m,n)
  end if
  mm1=m-1
  if (mm1.lt.0) then
    minus1(2,i)=0
  else
    minus1(2,i)=facind(l,mm1,n)
  end if
  nm1=n-1
  if (nm1.lt.0) then
    minus1(3,i)=0
  else
    minus1(3,i)=facind(l,m,nm1)
  end if
end do

min1b=0
outer: do i=2,min1max
  do jj=1,3
    im1=minus1(jj,i)
    if (im1.ne.0) then
      min1b(1,i)=jj
      min1b(2,i)=im1
      min1b(3,i)=indexps(jj,im1)
      cycle outer
    end if
  end do
end do outer

minus1s=minus1(:,0:35)
plus1s=plus1(:,0:20)
! For checking purposes, print this all out
!open (unit=12,file="Summary_of_indices.txt")
!write (12,*) "Summary of the Indices used for +1/-1"
!write (12,*) "functions (l,m,n) and their overall index:"
!do i=1,cter
!  l=indexps(1,i)
!  m=indexps(2,i)
!  n=indexps(3,i)
!  write (12,'("(",I1,",",I1,",",I1,")",I4)') l,m,n,facind(l,m,n)
!end do
 
!write (12,*) ""
!write (12,*) "Minus one and Plus 1 correspondences:"
!do i=1,cter
!  if (i.le.plu1max) then
!    write (12,'(I5," (",I1,",",I1,",",I1,")",X,3I4,4X,3I4)') i,indexps(:,i), minus1(:,i), plus1(:,i)
!  else 
!    write (12,'(I5," (",I1,",",I1,",",I1,")",X,3I4,2X)') i,indexps(:,i), minus1(:,i)
!  end if
!end do
!close(12)

end subroutine buildindices


subroutine buildharmond()
use nrtype ; use molprops
implicit none

! now contract to form spherical harmonics. Rem: Cartesians:
!     0    1    2    3    4    5    6    7    8    9   10
!     /    s    x    y    z   xx   xy   xz   yy   yz   zz
! Spherical harmonics:
!     1    2    3    4    5
!    z2   xz   yz   xy  x2-y2

matharmond=0._dp
matharmond(1,6)=1._dp
matharmond(1,1)=-.5_dp
matharmond(1,4)=-.5_dp
matharmond(2,3)=sqrt(3._dp)
matharmond(3,5)=sqrt(3._dp)
matharmond(4,2)=sqrt(3._dp)
matharmond(5,1)=sqrt(3._dp)/2._dp
matharmond(5,4)=-sqrt(3._dp)/2._dp

end subroutine buildharmond


subroutine buildgammarats()
use nrtype ; use molprops
implicit none

integer(i4b) :: i, j

allocate(gammrats(0:maxgamma,0:maxboys+3),boyszeroarg(0:maxboys),boyslrfac(maxboys),boysmrfac(maxboys))

do i=0,maxgamma
   do j=0,maxboys+3
       gammrats(i,j)=gamma(real(j,dp)+.5_dp)/gamma(real(i+j,dp)+1.5_dp)
   end do
end do

do i=0,maxboys
   boyszeroarg(i)=1._dp/real(1+2*i,dp)
end do
do i=1,maxboys
   boyslrfac(i)=real(i,dp)-0.5_dp
   boysmrfac(i)=2._dp/real(2*(i-1)+1,dp)
end do

end subroutine buildgammarats

subroutine eval_boys_func_on_grid()
use nrtype ; use molprops
implicit none

integer(i4b) :: i, j, jn
real(dp) :: x,xpow,bf(0:maxboys+3), efac
!real(dp) :: x,xpow,bf(0:maxboys+3), efac, bla(0:4), blub

allocate(boysfuncargs(0:ngamma),boysfuncvals(0:3,0:maxboys,0:ngamma))

do i=0,ngamma
  x=(real(i,dp)+0.5_dp)*twodelta
  boysfuncargs(i)=x
  bf=0._dp
  xpow=1._dp
  efac=exp(-x)/2._dp
  do jn=0,maxgamma
    bf(:)=bf(:)+xpow*gammrats(jn,:)
    xpow=xpow*x
  end do
  bf=bf*efac
  do j=0,maxboys
    boysfuncvals(0,j,i)=bf(j)+x*bf(j+1)+.5_dp*x**2*bf(j+2)+x**3*bf(j+3)/6._dp
    boysfuncvals(1,j,i)=bf(j+1)+x*bf(j+2)+.5_dp*x**2*bf(j+3)
    boysfuncvals(2,j,i)=.5_dp*bf(j+2)+.5_dp*x*bf(j+3)
    boysfuncvals(3,j,i)=bf(j+3)/6._dp
  end do
!  blub=boysfuncvals(0,0,i)-x*(boysfuncvals(1,0,i)-x*(boysfuncvals(2,0,i)-x*boysfuncvals(3,0,i)))
!  write (*,'(4F14.8)') x,bf(0),blub,bf(0)-blub
end do

end subroutine eval_boys_func_on_grid


