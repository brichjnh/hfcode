subroutine spread_basis()
use nrtype; use molprops
implicit none

! re-jig the order of the shells. Make all s shells come first, then p, then d (with perhaps l also)
! From basnspd we know how many of each shell type are on each atom. The first block works out the n. of atoms.

integer(i4b) :: i, ii, j, k, atmap
integer(i4b) :: ish, ishs, ishp, ishd, iao, jshs, jshp, jshd, ishl, jshl
character(len=1) :: labsh

! First work out how many shells of each type there are
!  Also how many basis functions there are in total
!  And need to set up basis set atom labels
nshs=0
nshp=0
nshd=0
nshl=0
do i=1,natom
   do j=1,nbasats
      if (trim(basels(j)).eq.trim(atlabels(i))) then
          atmap=j  ! atom i by order is atom atmap in the basis set list
          nshs=nshs+basnspd(atmap,1)
          nshp=nshp+basnspd(atmap,2)
          nshd=nshd+basnspd(atmap,3)
          nshl=nshl+basnspd(atmap,4)
          exit
      end if
   end do
end do
nb=nshs+3*nshp+5*nshd+4*nshl   ! Total number basis funcs
nb2=nb*(nb+1)/2                ! Total number basis func pairs
nb4=nb2*(nb2+1)/2              ! Total number unique basis function quadruples
nsh=nshs+nshp+nshd+nshl        ! Total number shells
nss=nshs*(nshs+1)/2            ! Number of s-s shell pairs
nsp=nshs*nshp                  ! Number of s-p shell pairs
nsd=nshs*nshd
nsl=nshs*nshl
npp=nshp*(nshp+1)/2
npd=nshp*nshd
npl=nshp*nshl
ndd=nshd*(nshd+1)/2
ndl=nshd*nshl
nll=nshl*(nshl+1)/2
nab=nsh*(nsh+1)/2

allocate(nzet(nsh),xsh(3,nsh),atsh(nsh),shtyp(nsh))
allocate(zeta(maxprim,nsh),djk(maxprim,nsh),djkl(maxprim,nsh),aosh(nsh))
allocate(nshat(natom),shat(natom,maxfuncel),nshtat(natom,4),shtyat(natom,maxfuncel,4))
allocate(atao(nb),aotyp(nb),aohuck(nb),nbfsat(natom),bsref(natom))
allocate(ssh(nshs),psh(nshp),dsh(nshd),lsh(nshl),shao(nb))
allocate(ssprs(2,nss),spprs(2,nsp),sdprs(2,nsd),slprs(2,nsl),ppprs(2,npp))
allocate(pdprs(2,npd),plprs(2,npl),ddprs(2,ndd),dlprs(2,ndl),llprs(2,nll))
allocate(schw(nab))
allocate(schwss(nss),schwsp(nsp),schwsd(nsd),schwsl(nsl),schwpp(npp))
allocate(schwpd(npd),schwpl(npl),schwdd(ndd),schwdl(ndl),schwll(nll))
allocate(orderss(nss),ordersp(nsp),ordersd(nsd),ordersl(nsl),orderpp(npp))
allocate(orderpd(npd),orderpl(npl),orderdd(ndd),orderdl(ndl),orderll(nll))

ish=0 ! counter for the basis shell array
ishs=0 ! counter for the s basis shell array
ishp=0 ! counter for the p basis shell array
ishd=0 ! counter for the d basis shell array
ishl=0 ! counter for the d basis shell array
iao=0 ! counter for the basis function array
do i=1,natom
  ! First map the atom in the structure onto the atom in the basis set
  do j=1,nbasats
     if (trim(basels(j)).eq.trim(atlabels(i))) then
         atmap=j  ! atom i by order is atom atmap in the basis set list
         bsref(i)=j
         exit
     end if
  end do
  ! Now set up the different shells for that atom
  nshat(i)=basnshell(atmap)
  nshtat(i,:)=basnspd(atmap,:)
  nbfsat(i)=nshtat(i,1)+nshtat(i,2)*3+nshtat(i,3)*5+nshtat(i,4)*4
  jshs=0
  jshp=0
  jshd=0
  jshl=0
  do j=1,nshat(i)
    ish=ish+1
    xsh(:,ish)=atcoords(i,:)
    zeta(:,ish)=baszeta(atmap,j,:)
    djk(:,ish)=basdjk(atmap,j,:)
    nzet(ish)=basnzet(atmap,j)
    atsh(ish)=i
    shat(i,j)=ish
    aosh(ish)=iao+1
    if (bastypes(atmap,j).eq."s") then
      jshs=jshs+1
      shtyat(i,jshs,1)=ish
      shtyp(ish)=1
      ishs=ishs+1
      ssh(ishs)=ish
      iao=iao+1
      if (jshs.le.2) then
        aohuck(iao)=uhs(jshs,atznumber(i))
      else
        aohuck(iao)=uhnv
      end if
      atao(iao)=i
      aotyp(iao)=1
      shao(iao)=ish
    else if (bastypes(atmap,j).eq."p") then
      jshp=jshp+1
      shtyat(i,jshp,2)=ish
      shtyp(ish)=2
      ishp=ishp+1
      psh(ishp)=ish
      iao=iao+3
      if (jshp.le.1) then
        aohuck(iao-2:iao)=uhp(atznumber(i))
      else
        aohuck(iao-2:iao)=uhnv
      end if
      atao(iao-2:iao)=i
      aotyp(iao-2:iao)=2
      shao(iao-2:iao)=ish
    else if (bastypes(atmap,j).eq."d") then
      jshd=jshd+1
      shtyat(i,jshd,3)=ish
      shtyp(ish)=3
      ishd=ishd+1
      dsh(ishd)=ish
      iao=iao+5
      aohuck(iao-4:iao)=uhnv
      atao(iao-4:iao)=i
      aotyp(iao-4:iao)=3
      shao(iao-4:iao)=ish
    else if (bastypes(atmap,j).eq."l") then
      jshl=jshl+1
      shtyat(i,jshl,4)=ish
      shtyp(ish)=12
      ishl=ishl+1
      lsh(ishl)=ish
      iao=iao+4
      if (jshl.le.1) then
        aohuck(iao-3)=uhs(2,atznumber(i))
        aohuck(iao-2:iao)=uhp(atznumber(i))
      else
        aohuck(iao-3:iao)=uhnv
      end if
      atao(iao-3:iao)=i
      aotyp(iao-3:iao)=12
      shao(iao-3:iao)=ish
      do ii=1,nzet(ish)
        djkl(ii,ish)=basdjkl(atmap,j,ii)/basdjk(atmap,j,ii)
      end do
    end if
  end do
end do

write (9,*) ""
write (9,*) "Overview of basis shells used"
write (9,*) "Shell No, BF no., Atom No, At label, F Type, coords; zetas; gcoefs"
do i=1,nsh
    ii=atsh(i)
    k=nzet(i)
    if (shtyp(i).eq.1) labsh="S"
    if (shtyp(i).eq.2) labsh="P"
    if (shtyp(i).eq.3) labsh="D"
    if (shtyp(i).eq.12) labsh="L"
    write (9,'(2i3,I4,X,A4,2X,A,2X,3F12.6)') i,ii,aosh(i),atlabels(ii),labsh,xsh(:,i)
    write (9,'(10f12.6)') (zeta(j,i),j=1,k)
    write (9,'(10f12.6)') (djk(j,i),j=1,k)
    if (shtyp(i).eq.12) write (9,'(10f12.6)') (djkl(j,i),j=1,k)
end do
write (9,*) "Total basis set properties:"
write (9,'(A,4I3)') "Number of S, P, D, L shells:",nshs,nshp,nshd,nshl
write (9,'(A,2I4)') "Total number of shells & BFs:",nsh,nb
if ((nshd.gt.0).and.(nshl.gt.0).and.(nshp.gt.0)) then
   write (9,*) "Mixed P, D and L shells not supported. ERROR"
   stop
end if

end subroutine spread_basis



