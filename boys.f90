module boys
use nrtype; use molprops
implicit none

contains

pure subroutine boys_0(x,f0)
use nrtype; use molprops
implicit none

real(dp), intent(in) :: x
real(dp), intent(out) :: f0

integer(i4b) :: j
real(dp), parameter :: vsmall = 1.d-15
real(dp), parameter :: boys_longrange  = 28._dp
real(dp), parameter :: boys_medrange  = 2._dp
real(dp), parameter :: hsqrtpi=sqrt(pi)/2._dp
real(dp) :: sqrtx

if (x.lt.vsmall) then
  f0=1._dp
  return
else if (x.gt.boys_longrange) then
  f0=hsqrtpi/sqrt(x)
  return
else if ((x.ge.boys_medrange).and.(x.le.boys_longrange)) then
  sqrtx=sqrt(x)
  f0=hsqrtpi*erf(sqrtx)/sqrtx
  return
else
  j=int(x/twodelta)
!  call repulsion_f(0,x,f0)
  f0=boysfuncvals(0,0,j)-x*(boysfuncvals(1,0,j)-x*(boysfuncvals(2,0,j)-x*boysfuncvals(3,0,j)))
  return
end if

end subroutine boys_0


pure subroutine boys_1(x,fn)
use nrtype; use molprops
implicit none

real(dp), intent(in) :: x
real(dp), intent(out) :: fn(0:1)

integer(i4b) :: j
real(dp) :: sqrtx
real(dp), parameter :: vsmall = 1.d-15
real(dp), parameter :: boys_longrange  = 28._dp
real(dp), parameter :: boys_medrange  = 4._dp
real(dp), parameter :: hsqrtpi=sqrt(pi)/2._dp

if (x.lt.vsmall) then
  fn=boyszeroarg(0:1)
  return
else if (x.gt.boys_longrange) then
  fn(0)=hsqrtpi/sqrt(x)
  fn(1)=fn(0)*.5_dp/x
  return
else if ((x.ge.boys_medrange).and.(x.le.boys_longrange)) then
  sqrtx=sqrt(x)
  fn(0)=hsqrtpi*erf(sqrtx)/sqrtx
  fn(1)=(fn(0)-exp(-x))/(2._dp*x)
  return
else
  j=int(x/twodelta)
  fn=boysfuncvals(0,0:1,j)-x*(boysfuncvals(1,0:1,j)-x*(boysfuncvals(2,0:1,j)-x*boysfuncvals(3,0:1,j)))
end if

end subroutine boys_1

pure subroutine boys_2(x,fn)
use nrtype; use molprops
implicit none

real(dp), intent(in) :: x
real(dp), intent(out) :: fn(0:2)

integer(i4b) :: j
real(dp) :: efac, sqrtx
real(dp), parameter :: vsmall = 1.d-15
real(dp), parameter :: boys_longrange  = 28._dp
real(dp), parameter :: boys_medrange  = 8._dp
real(dp), parameter :: hsqrtpi=sqrt(pi)/2._dp

if (x.lt.vsmall) then
  fn=boyszeroarg(0:2)
  return
else if (x.gt.boys_longrange) then
  fn(0)=hsqrtpi/sqrt(x)
  fn(1)=fn(0)*.5_dp/x
  fn(2)=fn(1)*1.5_dp/x
  return
else if ((x.ge.boys_medrange).and.(x.le.boys_longrange)) then
  sqrtx=sqrt(x)
  fn(0)=hsqrtpi*erf(sqrtx)/sqrtx
  !efac=exp(-x)
  efac=.5_dp*exp(-x)
  fn(1)=(.5_dp*fn(0)-efac)/x
  fn(2)=(1.5_dp*fn(1)-efac)/x
  return
else
  j=int(x/twodelta)
  efac=.5_dp*exp(-x)
  fn(2)=boysfuncvals(0,2,j)-x*(boysfuncvals(1,2,j)-x*(boysfuncvals(2,2,j)-x*boysfuncvals(3,2,j)))
  fn(1)=(x*fn(2)+efac)*2._dp/3._dp
  fn(0)=(x*fn(1)+efac)*2._dp
end if

end subroutine boys_2

pure subroutine boys_n(n,x,fn)
use nrtype; use molprops
implicit none

integer(i4b), intent(in) :: n
real(dp), intent(in) :: x
real(dp), intent(out) :: fn(0:n)

integer(i4b) :: i,j
real(dp) :: efac, sqrtx
real(dp), parameter :: vsmall = 1.d-15
real(dp), parameter :: boys_longrange  = 28._dp
real(dp), parameter :: boys_medrange  = 8._dp
real(dp), parameter :: hsqrtpi=sqrt(pi)/2._dp

if (x.lt.vsmall) then
  fn=boyszeroarg(0:n)
  return
else if (x.gt.boys_longrange) then
  fn(0)=hsqrtpi/sqrt(x)
  do i=1,n
    fn(i)=fn(i-1)*boyslrfac(i)/x
  end do
  return
else if ((x.ge.boys_medrange).and.(x.le.boys_longrange)) then
  sqrtx=sqrt(x)
  fn(0)=hsqrtpi*erf(sqrtx)/sqrtx
  !efac=exp(-x)
  efac=.5_dp*exp(-x)
  do i=0,n-1
    !fn(i+1)=(real((2*i+1),dp)*fn(i)-efac)/(2._dp*x)
    fn(i+1)=(boyslrfac(i+1)*fn(i)-efac)/x
  end do
  return
else
  j=int(x/twodelta)
  efac=.5_dp*exp(-x)
  fn(n)=boysfuncvals(0,n,j)-x*(boysfuncvals(1,n,j)-x*(boysfuncvals(2,n,j)-x*boysfuncvals(3,n,j)))
  do i=n,1,-1
    fn(i-1)=(x*fn(i)+efac)*boysmrfac(i)
    !fn(i-1)=(2._dp*x*fn(i)+2._dp*efac)/real(2*(i-1)+1,dp)
  end do
end if

end subroutine boys_n

end module boys

