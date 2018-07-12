subroutine read_basfile(basfile)
use nrtype; use molprops
implicit none

! In this new version, the basis set is read in by shell

character(len=30), intent(in) :: basfile

character(len=30) :: bla
integer(i4b) :: i,j,k,lbas,nz, nbfsel, ii, jj,kk, dmat

open(unit=8,file=basfile)
read(8,*) bla,nbasats,dmat
if (nbasats.gt.maxelbas) then
    write (*,*) "Basis File with more atoms than maxelbas:",nbasats,maxelbas
    stop
end if

allocate(basels(nbasats),bastypes(nbasats,maxfuncel))
allocate(basnshell(nbasats),basnspd(nbasats,4))
allocate(basnzet(nbasats,maxfuncel))
allocate(baszeta(nbasats,maxfuncel,maxprim),basdjk(nbasats,maxfuncel,maxprim))
allocate(basdjkl(nbasats,maxfuncel,maxprim))
allocate(basdmat(maxfuncat,maxfuncat,nbasats))

! Hard-wire some Hucke theory parameters: (all except H, He, C and O are presently very rough guesses;
!   the H, He, C and O values come from M06 atomic calcs)

uhs=reshape((/-.3,.8,-0.66,1.1,-1.2,-.05,-2.,-.1,-3.,-.15,-4.,-.2,-10.33,-0.43,-14.,-.5,-19.3,-0.775,-25.,-.9,-30.,-1. /),(/2,10/))
uhp=(/.8,.6,.1,-0.05,-.1,-.165,-.20,-.298,-.4,-.5/)

basnspd=0
do i=1,nbasats
   read(8,*) basnshell(i)
   read(8,*) basels(i)
   call to_lower(basels(i))
   do j=1,basnshell(i)
      read(8,*) bastypes(i,j), basnzet(i,j)
      call to_lower(bastypes(i,j))
      if (bastypes(i,j).eq."s") then
         lbas=0
         basnspd(i,1)=basnspd(i,1)+1
      else if (bastypes(i,j).eq."p") then
         lbas=1
         basnspd(i,2)=basnspd(i,2)+1
      else if (bastypes(i,j).eq."d") then
         lbas=2
         basnspd(i,3)=basnspd(i,3)+1
      else if (bastypes(i,j).eq."l") then
         lbas=12
         basnspd(i,4)=basnspd(i,4)+1
      end if
      if ((lbas.eq.0).or.(lbas.eq.1).or.(lbas.eq.2)) then
         do k=1,basnzet(i,j)
             read(8,*) baszeta(i,j,k), basdjk(i,j,k)
         end do
         nz=basnzet(i,j)
         call normalize(nz,baszeta(i,j,1:nz),lbas,basdjk(i,j,1:nz))
      else if (lbas.eq.12) then
         do k=1,basnzet(i,j)
             read(8,*) baszeta(i,j,k), basdjk(i,j,k), basdjkl(i,j,k)
         end do
         nz=basnzet(i,j)
         call normalize(nz,baszeta(i,j,1:nz),0,basdjk(i,j,1:nz))
         call normalize(nz,baszeta(i,j,1:nz),1,basdjkl(i,j,1:nz))
      end if
   end do
   read(8,*) bla
end do
!return
if (dmat.eq.0) then
  close(8)
  return
end if
read(8,*) bla
do i=1,nbasats
   nbfsel=basnspd(i,1)+3*basnspd(i,2)+basnspd(i,3)*5+basnspd(i,4)*4
   k=int(nbfsel/7)
   read (8,*) bla
   do ii=1,nbfsel
     do jj=1,k
       kk=(jj-1)*7
       read(8,*) basdmat(ii,kk+1:kk+7,i)
     end do
     kk=k*7
     if (kk.lt.nbfsel) then
       read(8,*) basdmat(ii,kk+1:nbfsel,i)
     end if
   end do
   read (8,*) bla
end do

close(8)
return

end subroutine read_basfile



subroutine normalize(n,zetas,l,coefs)
use nrtype
implicit none

! Most basis sets are reported with coefficients such that the overall contracted basis set is normalized,
! but the contraction coefficients do not include the normalization constants for the primitive gaussians
! i.e. you get CGF_j = \sum_k c_jk NPGF_k, where CGF_j is a contracted gaussian function, c_jk is the reported coefficient,
! and NPGF_k is a *normalized* primitive Gaussian function.
! For the program, it is more useful to have coefficients that include the normalization factor, so that one has:
! CGF_j = \sum_k d_jk PGF_k, where PGF_k is not normalized, and is given simply by (x-Ax)^ax (y-Ay)^ay (z-Az)^az exp(-zeta |r-rA|^2)
! This routine takes as input the c_jk, which are assumed to relate to NPGFs, and returns in ytheir place the d_jks, after normalizing.

integer(i4b), intent(in) :: n, l
real(dp), intent(in) :: zetas(n)
real(dp), intent(inout) :: coefs(n)

integer(i4b) :: i, j
real(dp) :: a, b, p, ofac, nfac
real(dp) :: r0(3), gcoef(n)

r0=0._dp

gcoef=coefs

if (l.eq.0) then
   do i=1,n
      a=zetas(i)
      ofac=(2._dp*a/PI)**0.75_dp
      gcoef(i)=gcoef(i)*ofac
   end do
else if (l.eq.1) then
   do i=1,n
      a=zetas(i)
      ofac=(2._dp/PI)**0.75_dp*2._dp*a**1.25_dp
      gcoef(i)=gcoef(i)*ofac
   end do
else if (l.eq.2) then
   do i=1,n
      a=zetas(i)
      ofac=(2._dp/PI)**0.75_dp*4._dp*a**1.75_dp/sqrt(3._dp)
      gcoef(i)=gcoef(i)*ofac
   end do
end if
nfac=0._dp
if (l.eq.0) then
   do i=1,n
      a=zetas(i)
      do j=1,n
         b=zetas(j)
         p=a+b
         nfac=nfac+pi32/sqrt(p**3)*gcoef(i)*gcoef(j)
      end do
   end do
else if (l.eq.1) then
   do i=1,n
      a=zetas(i)
      do j=1,n
         b=zetas(j)
         p=a+b
         nfac=nfac+pi32/(2._dp*p**2*sqrt(p))*gcoef(i)*gcoef(j)
      end do
   end do
else if (l.eq.2) then
! Rem by O-S, (d|d) for the same centre is 1/2p(s|d)+2/2p(p|p)
! (s|d) = 1/2p(s|s) and (p|p) = 1/2p(s|s)
! So (d|d)=1/2p(1/2p(s|s)+2/2p(s|s)) = 3 x (1/2p)^2 (s|s)
   do i=1,n
      a=zetas(i)
      do j=1,n
         b=zetas(j)
         p=a+b
         nfac=nfac+pi32*3._dp/(4._dp*p**3*sqrt(p))*gcoef(i)*gcoef(j)
      end do
   end do
end if
gcoef=gcoef/sqrt(nfac)
coefs=gcoef
return

end subroutine normalize

