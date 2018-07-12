subroutine direct_hartreefock()
use nrtype; use molprops
implicit none


! My 'best protocol' for HF convergence is implemented here, and can be summarized as follows:
!  * Use atomic densities or Core Hamiltonian for initial guess.
!  * For the next two cycles, use weighting: the F matrix is assembled using the average of the previous two D matrices
!  * DIIS is turned on once the change in D matrix drops below a given threshold; DIIS uses previous 4 F matrices (so 5 in total)
!  * A user-specified level shift can be used - it is turned off once DIIS is initiated
!  * Istart using a low threshold for integrals, then switch to tight threshold at the same time as starting DIIS

integer(i4b), parameter :: maxiter=25, ndiis=8
real(dp) :: f(nb,nb),fob(nb,nb), cmat(nb,nb), cmatob(nb,nb), p(nb,nb), pold(nb,nb), mprod,invseigs(nb,nb)
real(dp) :: fold(nb,nb,ndiis),fdiis(nb,nb), hcore(nb,nb), sminhalf(nb,nb), delp(nb,nb), delf(nb,nb), faobold(nb,nb)
real(dp) :: resmax,resid(nb,nb),resido(nb,nb,ndiis),eo(ndiis),shftmat(nb,nb)
real(dp), allocatable :: vecb(:),solb(:),bmat(:,:)
real(dp) :: eigs(nb), eeold, deltae, deltap, thresh,cyctime,tottime
real(dp), parameter :: conv=1.d-5, conve=5.d-8, lowthresh=5.d-7, fdifffac=0.1_dp
integer(i4b) :: i,j,k,jk,ii, mu, nu, nusediis,nfull
logical(lgt) :: diis, fullfock, uselevshft

! First obtain the matrix S^-1/2 that orthogonalizes the basis
call wrapper_dsyev(nb,s1e,eigs,p)
invseigs=0._dp
do i=1,nb
    invseigs(i,i)=1._dp/sqrt(eigs(i))
end do
sminhalf=matmul(p,matmul(invseigs,transpose(p)))

! Set up some things for DIIS among others
resido=1000._dp
fold=0._dp
eo=0._dp
diis=.false.
fullfock=.true.
hcore=t1e+v1e
eelec=0._dp
thresh=lowthresh
! Start timer
call time_diff(0,cyctime,tottime)
if (levshft) then
   shftmat=0._dp
   do k=nocc+1,nb
      shftmat(k,k)=lshft
   end do
end if

! Build the guess for the Fock matrix - either core hamiltonian or superposition of densities
if (densup) then
  p=0._dp
  k=0
  do i=1,natom
     p(k+1:k+nbfsat(i),k+1:k+nbfsat(i))=basdmat(1:nbfsat(i),1:nbfsat(i),bsref(i))
     k=k+nbfsat(i)
  end do
  write (9,*) "Built a guess density matrix by superimposing atom densities."
  pold=p
  uselevshft=.false.
else
  fob=matmul(sminhalf,matmul(hcore,transpose(sminhalf)))
  call wrapper_dsyev(nb,fob,eigs,cmatob)
  cmat=matmul(sminhalf,cmatob)
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
  ! calc initial energy
  eelec=0._dp
  do nu = 1, nb
    do mu = 1, nb
      eelec=eelec+p(mu,nu)*hcore(mu,nu)
    end do
  end do
  write (9,'(A,2F20.10)') "Electronic & total Energy from Hcore:",eelec,eelec+unuc
  uselevshft=.true.
end if

nusediis=0

write (9,*) ""
write (9,*) "Starting HF cycles: cycle, etot, deltae, deltap, nshpr screened, tstep, ttotal"
write (9,'(A,I15)') "Total number of shell quadruples:",nab*(nab+1)/2
do i=1,maxiter
    f=hcore
    eeold=eelec
    if (fullfock.and.(i.eq.1)) then
       call fock_build_initial_diag(f,p,thresh)
       call fock_build_non_diag(f,p,thresh)
       !fullfock=.false.
       nfull=0
    else if (fullfock) then
       call fock_build_diag(f,p,thresh)
       call fock_build_non_diag(f,p,thresh)
       fullfock=.false.
       nfull=0
    else
       delp=p-pold
       delf=0._dp
       call fock_build_diag(delf,delp,thresh)
       call fock_build_non_diag(delf,delp,thresh*fdifffac)
       f=faobold+delf
       nfull=nfull+1
       if (nfull.eq.15) then
         fullfock=.true.
       end if
    end if
! get energy for new fock matrix
    eelec=0._dp
    do nu = 1, nb
        do mu = 1, nb
            eelec=eelec+p(mu,nu)*(hcore(mu,nu)+F(mu,nu))
        end do
    end do
    eelec=eelec/2._dp
    deltae=eelec-eeold
    etot=eelec+unuc
    faobold=f
    fob=matmul(sminhalf,matmul(f,transpose(sminhalf)))
    invseigs=matmul(f,matmul(p,s1e))-matmul(s1e,matmul(p,f))
    resid=matmul(transpose(sminhalf),matmul(invseigs,sminhalf))
    do j=1,ndiis-1
        resido(:,:,ndiis+1-j)=resido(:,:,ndiis-j)
        fold(:,:,ndiis+1-j)=fold(:,:,ndiis-j)
        eo(ndiis+1-j)=eo(ndiis-j)
    end do
    resido(:,:,1)=resid
    resmax=maxval(abs(resid))
    fold(:,:,1)=fob
    eo(1)=eelec
    if ((i.ge.2).and.(resmax.lt.0.3).and.(.not.diis)) then
       diis=.true.
    !   fullfock=.true.
       uselevshft=.false.
       write (9,*) "Starting DIIS"
       write (*,*) "Starting DIIS"
    end if
    if ((i.ge.7).and.(thresh.gt.intthresh).and.(resmax.lt.0.005)) then
       write (9,*) "Rebuilding Fock Matrix, lowering threshold for ints, turning off level shift, resetting DIIS"
       write (*,*) "Rebuilding Fock Matrix, lowering threshold for ints, turning off level shift, resetting DIIS"
       eo=2000.
       nusediis=0
       fullfock=.true.
       uselevshft=.false.
       thresh=intthresh
    end if
    if (diis) then
       nusediis=nusediis+1
       if (nusediis.gt.ndiis) nusediis=ndiis
       if (nusediis.gt.1) then
!        build b matrix
         allocate(vecb(nusediis+1),solb(nusediis+1),bmat(nusediis+1,nusediis+1))
         bmat=0._dp
         bmat(nusediis+1,:)=-1._dp; bmat(:,nusediis+1)=-1._dp; bmat(nusediis+1,nusediis+1)=0._dp
         vecb=0._dp ; vecb(nusediis+1)=-1._dp
         do j=1,nusediis
           do k=1,j
             mprod=0._dp
             do jk=1,nb
               mprod=mprod+sum(resido(:,jk,j)*resido(:,jk,k))
             end do
             bmat(j,k)=mprod
             bmat(k,j)=bmat(j,k)
           end do
         end do
         ii=minloc(eo,1)
         do j=1,nusediis
           if (j.eq.ii) cycle
           bmat(j,j)=bmat(j,j)*1.05
         end do
         call wrapper_dsysv(nusediis+1,bmat,vecb,solb)
!write (*,'(A,8f12.5)') "diis coeffs",solb(1:nusediis)
         if (maxval(abs(solb(1:nusediis))).lt.1.5_dp) then
           fdiis=0._dp
           do j=1,nusediis
               fdiis=fdiis+solb(j)*fold(:,:,j)
           end do
           deallocate(vecb,solb,bmat)
           fob=fdiis
         else
           nusediis=0
           deallocate(vecb,solb,bmat)
         end if
       end if
    end if
    if (uselevshft.and.levshft) then
!write (*,*) "applying level shift"
       fob=fob+matmul(cmatob,matmul(shftmat,transpose(cmatob)))
    end if
    call wrapper_dsyev(nb,fob,eigs,cmatob)
    cmat=matmul(sminhalf,cmatob)

    deltap=maxval(abs(p-pold))
    call time_diff(1,cyctime,tottime)
    write (9,'(I3,3F20.10,2F12.2)') i,etot,deltae,deltap,cyctime,tottime
    write (*,'(I3,3F20.10,2F12.2)') i,etot,deltae,deltap,cyctime,tottime
    if (abs(deltae).lt.conve) then
        if ((deltap.lt.conv).and.(i.gt.1)) then
             write (9,*) "SCF Converged"
             exit
        end if
    end if
    pold=p
    p=2._dp*matmul(cmat(:,1:nocc),transpose(cmat(:,1:nocc)))
!    if (i.le.2) then
!        p=.5_dp*p+.5_dp*pold
!        uselevshft=.true.
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

return
end subroutine direct_hartreefock
