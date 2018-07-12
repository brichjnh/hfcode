module molprops
use nrtype
implicit none

! Global program settings, limits, etc.

save

integer(i4b), parameter :: maxat=200        ! Maximum total number of atoms
! integer(i4b), parameter :: maxb=800        ! Maximum total number of basis functions
! integer(i4b), parameter :: maxsh=200        ! Maximum total number of shells
integer(i4b), parameter :: maxelbas=18     ! maximum number of element in a basis file
integer(i4b), parameter :: maxfuncel=10    ! maximum number of basis function shells for a given atom
integer(i4b), parameter :: maxfuncat=40    ! maximum number of basis functions for a given atom
integer(i4b), parameter :: maxgamma=50     ! maximum value of gamma function used for evaluating Boys function
integer(i4b), parameter :: ngamma=10000    ! number of points at which the Boys function is calculated on a grid
integer(i4b), parameter :: maxboys=9       ! maximum index of the Boys function needed. For 1st derivs = 4*maxL + 1
integer(i4b), parameter :: maxprim=10      ! maximum number of primitives for a given function shell
integer(i4b) :: prtlevl                    ! printlevel: 1=normal, 2=high, 3=ultrahigh
integer(i4b) :: nprtsh
integer(i4b),allocatable :: prsh(:,:)                   ! For testing purposes, printout this shell worth's of integras

! Huckel theory parameters
real(dp) :: uhs(2,10)   ! Huckel energies for s orbitals for H - Ne
real(dp) :: uhp(10)     ! Huckecl energies for p orbitals
real(dp), parameter :: uhnv=0.8_dp   ! Energy level shift for higher-energy s or p orbitals


! Some parameters relating to index incrementing/decrementing.
integer(i4b), parameter :: labcdmax=8          ! maximum number of primitives for a given function shell - yielding the values below as n. entries for minus1 and plus1
integer(i4b), parameter :: min1max=(labcdmax+1)*(labcdmax+2)*(labcdmax+3)/6, plu1max=(labcdmax+1)*labcdmax*(labcdmax+2)/6
integer(i4b) :: min1b(3,0:min1max)
integer(i4b) :: minus1(3,0:min1max),plus1(3,0:plu1max)     ! Tables showing how one ang momentum index converts to another - up and down
integer(i4b) :: indexps(3,0:min1max)                       ! Table showing how many ang momentum of each type are present for a given entry
integer(i4b) :: minus1s(3,0:35),plus1s(3,0:20)             ! Smaller versions of the above cross-reference Tables

integer(i4b), allocatable :: ind2(:,:)    ! compisite index, pre-computed
integer(i4b), allocatable :: ioff(:)      ! compisite index, pre-computed
integer(i4b), allocatable :: revind2(:,:) ! To get the elementary indices back


! General molecule variables

character(len=40) :: title                   ! title of job - not really used
integer(i4b) :: natom, nb, nb2, nb4, nocc, ncore       ! Total number of atoms, basis funcs, bfunc pairs, occupied orbs, core orbs
integer(i4b) :: sz, nocca, noccb                       ! Open shell: na-nb, na, nb

! Overall things like coordinates, MO coefs, energy, etc.

character(len=2), allocatable :: atlabels(:)

real(dp) :: unuc, eelec, etot, emp2, lshft  ! nuclear, electronic, total and mp2 energies, level shift
logical(lgt) :: mp2energy, ccenergy, densup, uhf, levshft, directscf

real(dp), allocatable :: atcoords(:,:)     ! atomic coordinates in bohr
real(dp), allocatable :: angcoords(:,:)    ! atomic coordinates in Angstrom
real(dp), allocatable :: nuccharges(:)     ! nuclear charges
integer(i4b), allocatable :: atznumber(:)  ! atomic number

real(dp), allocatable :: cij(:,:)          ! MO coefficients
real(dp), allocatable :: moen(:)           ! MO energies

! Variables relative to the basis set family

integer(i4b) :: nbasats  ! Number of atoms included in the basis list

character(len=2), allocatable :: basels(:)      ! labels of atoms in basis - index is atom number
character(len=1), allocatable :: bastypes(:,:)  ! s, p, d etc. - indexed per atom and shell

integer(i4b), allocatable :: basnshell(:)    ! Total number of basis function shells for a given element
integer(i4b), allocatable :: basnspd(:,:)    ! number of s, p, d, l basis function shells for a given element
integer(i4b), allocatable :: basnzet(:,:)    ! n. primitives in a given basis function shell, by atom & shell

real(dp), allocatable :: baszeta(:,:,:)  ! zetas for a given atom shell, per atom, then shell, then zetas
real(dp), allocatable :: basdjk(:,:,:)   ! contraction coefficients d_jk for s shells, per atom, shell, zetas
real(dp), allocatable :: basdjkl(:,:,:)  ! additional contraction coefficients d_jk for l shells, per atom, shell, zetas

real(dp), allocatable :: basdmat(:,:,:)  ! density matrices

! Now basis set variables for the specific molecular system
! These variables are of two types: ones needed to carry out the integrals - grouped by shell
! and ones needed to do book-keeping (say if I ever implement Mulliken population or the like) - grouped by function

integer(i4b), allocatable :: shtyp(:)      ! Type of shell: s=0, p=1, d=2, l=12
integer(i4b), allocatable :: nzet(:)       ! number of zetas for a given (contracted) shell of basis functions
integer(i4b), allocatable :: atsh(:)       ! atom on which a given shell of basis functions is centred
integer(i4b), allocatable :: nshat(:)      ! For each atom, gives the total number of shells that are based on it
integer(i4b), allocatable :: nshtat(:,:)   ! For each atom, gives the number of s, p, d shells that are based on it
integer(i4b), allocatable :: shat(:,:)     ! For each atom, gives the shells that are based on it
integer(i4b), allocatable :: nbfsat(:)     ! For each atom, gives the number of basis functions on it
integer(i4b), allocatable :: bsref(:)      ! For each atom, gives the basis set order reference
integer(i4b), allocatable :: shtyat(:,:,:) ! For each atom, gives the shells, by type, that are based on it
integer(i4b), allocatable :: aosh(:)       ! For each shell, gives the basis function reference for the *first* func in teh shell
integer(i4b), allocatable :: ssh(:)        ! Shell numbers for s shells
integer(i4b), allocatable :: psh(:)        ! Shell numbers for p shells
integer(i4b), allocatable :: dsh(:)        ! Shell numbers for d shells
integer(i4b), allocatable :: lsh(:)        ! Shell numbers for l shells
real(dp), allocatable :: xsh(:,:)          ! coordinates of atom on which a shell of basis funcs is positioned
real(dp), allocatable :: zeta(:,:)         ! zeta values for a given shell of basis functions
real(dp), allocatable :: djk(:,:)          ! "full" (including normalization fac) contraction coefs for a shell
real(dp), allocatable :: djkl(:,:)         ! "full" (including normalization fac) relative contraction coefs for p part of l shell
real(dp), allocatable :: gammrats(:,:)     ! ratios of Gamma functions used to evaluate Boys function
real(dp), allocatable :: boyszeroarg(:)    ! Values of Boys function for zero argument
real(dp), allocatable :: boysfuncvals(:,:,:) ! Values used for Ishida-like interpolation of Boys function
real(dp), allocatable :: boysfuncargs(:)   ! Values at which the Boys function is computed and stored 
real(dp), allocatable :: boyslrfac(:)      ! Factors used for upwards recursion of Boys functions at long range
real(dp), allocatable :: boysmrfac(:)      ! Factors used for downwards recursion of Boys functions at medium range
real(dp), parameter :: twodelta=1.d-3      ! Spacing of Boys function grid
real(dp) :: matharmond(5,6)                ! matrix showing link between Cartesian and sph harm d funcs
real(dp), parameter :: dnorm(6)=(/1._dp,sqrt(3._dp),sqrt(3._dp),1._dp,sqrt(3._dp),1._dp/)  ! Different normalization

integer(i4b), allocatable :: atao(:)    ! atom on which a given basis function is centred
real(dp), allocatable :: aohuck(:)      ! Huckel estimate of orbital energy
integer(i4b), allocatable :: shao(:)    ! cross-reference of the basis shell number for a given basis function
integer(i4b), allocatable :: aotyp(:)   ! code for the type of basis function: s: 1; px,y,z: 2-4; dxy,...: 5-9
integer(i4b) :: nsh, nshs, nshp, nshd, nshl  ! Total number of shells, and of s, p, d

! Finally some variables relating to integrals

integer(i4b) :: nab, nzz                 ! Number of shell pairs and max number primitive pairs
integer(i4b) :: nss, nsp, nsd, nsl, npp  ! No. of s,s, s,p, s,d, s,l, p,p shell pairs
integer(i4b) :: npd, npl, ndd, ndl, nll  ! No. of p,d, p,l, d,d, d,l, l,l shell pairs
integer(i4b), allocatable :: ssprs(:,:)  ! List containing shell number for s,s, shell pairs
integer(i4b), allocatable :: spprs(:,:)  ! List containing shell number for s,p shell pairs
integer(i4b), allocatable :: sdprs(:,:)  ! List containing shell number for s,d shell pairs
integer(i4b), allocatable :: slprs(:,:)  ! List containing shell number for s,l shell pairs
integer(i4b), allocatable :: ppprs(:,:)  ! List containing shell number for p,p shell pairs
integer(i4b), allocatable :: pdprs(:,:)  ! List containing shell number for p,d shell pairs
integer(i4b), allocatable :: plprs(:,:)  ! List containing shell number for p,l shell pairs
integer(i4b), allocatable :: ddprs(:,:)  ! List containing shell number for d,d shell pairs
integer(i4b), allocatable :: dlprs(:,:)  ! List containing shell number for d,l shell pairs
integer(i4b), allocatable :: llprs(:,:)  ! List containing shell number for l,l shell pairs
real(dp), allocatable :: pab(:,:)        ! p = zeta_a + zeta_b for all shair pairs
real(dp), allocatable :: pm1(:,:)        ! 1/p
real(dp), allocatable :: hpm1(:,:)       ! 1/2p
real(dp), allocatable :: qab(:,:)        ! zeta_a*zeta_b/pab
real(dp), allocatable :: rab2(:)         ! square of distance from centre of shell A to centre of shell B
real(dp), allocatable :: r2AB(:,:)       ! square of distance from atom A to atom B
real(dp), allocatable :: Px(:,:,:)       ! distribution centre coords (zeta_a*Ax + zeta_b*Bx)/pab

logical(lgt), allocatable :: negab(:,:)  ! set to .true. if kab/sqrt(p) is less than intthresh
real(dp), allocatable :: ddab(:,:)       ! djka * djkb
real(dp), allocatable :: kab(:,:)        ! fac * exp(-q*rAB^2) - includes contraction coefficients product
real(dp), allocatable :: schw(:)         ! square root of each 'diagonal' (ab|ab) two-electron integral, for Schwartz screening
real(dp), allocatable :: schwss(:), schwsp(:), schwsd(:), schwsl(:), schwpp(:), schwpd(:)
real(dp), allocatable :: schwpl(:), schwdd(:), schwdl(:), schwll(:)
integer(i4b), allocatable :: orderss(:), ordersp(:), ordersd(:), ordersl(:), orderpp(:), orderpd(:)
integer(i4b) :: lastss, lastsp, lastsd, lastsl, lastpp, lastpd, lastpl, lastdd, lastdl, lastll
integer(i4b), allocatable :: orderpl(:), orderdd(:), orderdl(:), orderll(:)   ! integer arrays storing the order of size of Schwartz terms
real(dp), allocatable :: s1e(:,:)        ! One electron integrals
real(dp), allocatable :: t1e(:,:)        ! Kinetic energy integrals
real(dp), allocatable :: v1e(:,:)        ! Nuclear Coulomb integrals

integer(i4b) :: naaaa,naaab,naabb,nabab,naabc,nabcd
integer(i4b) :: maxaaaa,maxaaab,maxaabb,maxabab,maxaabc,maxabcd
integer(i4b), allocatable :: indaaaa(:)
integer(i4b), allocatable :: indaaab(:,:)
integer(i4b), allocatable :: indaabb(:,:)
integer(i4b), allocatable :: indabab(:,:)
integer(i4b), allocatable :: indaabc(:,:)
integer(i4b), allocatable :: indabcd(:,:)
real(dp), parameter :: intthresh=1.d-10
real(dp), parameter :: shthresh=intthresh/10._dp
real(dp), parameter :: schwthresh=1.d-10
real(dp), allocatable :: i2e(:)
real(dp), allocatable :: i2e4i(:,:,:,:)
real(dp), allocatable :: m1e(:,:)
real(dp), allocatable :: m2e(:,:,:,:)

end module molprops


