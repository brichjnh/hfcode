subroutine incore_hartreefock()
use nrtype; use molprops
implicit none

integer(i4b), parameter :: maxiter=50, ndiis=5
real(dp) :: f(nb,nb),fprim(nb,nb), cmat(nb,nb), p(nb,nb), pold2(nb,nb),pold(nb,nb), mprod,invseigs(nb,nb)
real(dp) :: fold(nb,nb,ndiis),fdiis(nb,nb), hcore(nb,nb), sminhalf(nb,nb), pinit(nb,nb), tmat(nb,nb)
real(dp) :: resmax,resid(nb,nb),resido(nb,nb,ndiis), bmat(ndiis+1,ndiis+1),vecb(ndiis+1),solb(ndiis+1)
real(dp) :: eigs(nb), eeold, deltae, deltap
real(dp), parameter :: conv=1.d-5, conve=1.d-8
integer(i4b) :: i, j, k, jk, mu, nu
logical(lgt) :: diis

! Then obtain the matrix S^-1/2 that orthogonalizes the basis
call wrapper_dsyev(nb,s1e,eigs,p)
invseigs=0._dp
do i=1,nb
    invseigs(i,i)=1._dp/sqrt(eigs(i))
end do
sminhalf=matmul(p,matmul(invseigs,transpose(p)))

! Build the guess for the Fock matrix - either core hamiltonian or Huckel
if (densup) then
  p=0._dp
  k=0
  do i=1,natom
     p(k+1:k+nbfsat(i),k+1:k+nbfsat(i))=basdmat(1:nbfsat(i),1:nbfsat(i),bsref(i))
     k=k+nbfsat(i)
  end do
  write (9,*) "Built a guess density matrix by superimposing atom densities."
  call wrapper_dsyev(nb,p,eigs,cmat)
  cmat=matmul(sminhalf,cmat)
  write (9,*) "eigenvalues of density-based density matrix:"
  write (9,'(10F7.3)') eigs
  pinit=p
  pold=p
else
  hcore=t1e+v1e
  f=matmul(sminhalf,matmul(hcore,transpose(sminhalf)))
  call wrapper_dsyev(nb,f,eigs,cmat)
  cmat=matmul(sminhalf,cmat)
  if (prtlevl.ge.2) then
    write (9,*) ""
    write (9,*) "Eigenvalues for core or Huckel Hamiltonian:"
    write (9,'(5F12.5)') eigs
  end if
  if (prtlevl.ge.3) then
    write (9,*) "Eigenvectors for core or Huckel Hamiltonian:"
    do i=1,nb
      if (nb.le.7) then
        write (9,'(i3,7F12.6)') i,cmat(:,i)
      else
        write (9,'(i3,7F12.6)') i,cmat(1:7,i)
        write (9,'(3x,7F12.6)') cmat(8:,i)
      end if
    end do
  end if
  p=2._dp*matmul(cmat(:,1:nocc),transpose(cmat(:,1:nocc)))
  pold=p
end if


eelec=0._dp
! Finally assign hcore 'properly':
hcore=t1e+v1e

! Set up some things for DIIS
resido=1000._dp
fold=0._dp
diis=.false.

write (9,*) ""
write (9,*) "Starting HF cycles: cycle, etot, deltae, deltap"
do i=1,maxiter
  f=hcore
  call fock_build_incore(f,p)
  fprim=matmul(transpose(sminhalf),matmul(f,sminhalf))
  eeold=eelec
  tmat=matmul(f,matmul(p,s1e))-matmul(s1e,matmul(p,f))
  resid=matmul(transpose(sminhalf),matmul(tmat,sminhalf))
  do j=1,ndiis-1
      resido(:,:,ndiis+1-j)=resido(:,:,ndiis-j)
      fold(:,:,ndiis+1-j)=fold(:,:,ndiis-j)
  end do
  resido(:,:,1)=resid
  resmax=maxval(abs(resid))
  fold(:,:,1)=fprim
  if ((i.ge.5).and.(resmax.lt.0.1).and.(.not.diis)) then
     diis=.true.
     write (9,*) "Starting DIIS"
  end if
  if (diis) then
!    build b matrix
     bmat=0._dp
     bmat(ndiis+1,:)=-1._dp; bmat(:,ndiis+1)=-1._dp; bmat(ndiis+1,ndiis+1)=0._dp
     vecb=0._dp ; vecb(ndiis+1)=-1._dp
     do j=1,ndiis
        do k=1,j
            mprod=0._dp
            do jk=1,nb
                mprod=mprod+sum(resido(:,jk,j)*resido(:,jk,k))
            end do
            bmat(j,k)=mprod
            bmat(k,j)=bmat(j,k)
        end do
     end do
     call wrapper_dsysv(ndiis+1,bmat,vecb,solb)
     fdiis=0._dp
     do j=1,ndiis
         fdiis=fdiis+solb(j)*fold(:,:,j)
     end do
     call wrapper_dsyev(nb,fdiis,eigs,cmat)
     cmat=matmul(sminhalf,cmat)
  else
     call wrapper_dsyev(nb,fprim,eigs,cmat)
     cmat=matmul(sminhalf,cmat)
  end if

  eelec=0._dp
  do nu = 1, nb
      do mu = 1, nb
          eelec=eelec+p(nu,mu)*(hcore(mu,nu)+F(mu,nu))
      end do
  end do
  eelec=eelec/2._dp
  deltae=eelec-eeold
  deltap=maxval(abs(p-pold))
  etot=eelec+unuc
  write (9,'(I3,3F20.10)') i,etot,deltae,deltap
  if (abs(deltae).lt.conve) then
      if ((deltap.lt.conv).and.(i.gt.1)) then
           write (9,*) "SCF Converged"
           exit
      end if
  end if
  pold2=pold
  pold=p
  p=2._dp*matmul(cmat(:,1:nocc),transpose(cmat(:,1:nocc)))
!    if ((i.le.4).and.(.not.diis)) then
!        p=.5_dp*p+.5_dp*pold
!    else if ((i.gt.2).and.(i.le.6)) then
!        p=.5_dp*p+.33_dp*pold+.17_dp*pold2
!    end if
end do
write (9,*) "Final Eigenvalues:"
write (9,'(5F12.5)') eigs
if (prtlevl.ge.3) then
  write (9,*) "Final Eigenvectors:"
  do i=1,nb
  write (9,'(7F12.6)') cmat(i,:)
  end do
end if
allocate(cij(nb,nb))
cij=transpose(cmat)
!write (9,*) "Final Density matrix:"
!do i=1,nb
!  write (9,'(10F12.5)') p(i,:)
!end do
!write (9,*) "Difference Guess and Final Density matrix:"
!do i=1,nb
!  write (9,'(10F12.5)') p(i,:)-pinit(i,:)
!end do
return

end subroutine incore_hartreefock

