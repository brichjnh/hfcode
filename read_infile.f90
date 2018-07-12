subroutine read_infile(infile)
use nrtype; use molprops
implicit none

character(len=30), intent(in) :: infile
real(dp), parameter :: bohr=1.889725989_dp
character(len=10) :: bla
character(len=1) :: domp2,docc,douhf,dodensup,dolevshft, dodirect
integer(i4b) :: i,j

open(unit=8,file=infile)
read(8,*) natom
if (natom.gt.maxat) then
    write (*,*) "File with more atoms than maxat:",natom,maxat
    stop
end if
read(8,*) title
read(8,*) bla,dodensup
call to_lower(dodensup)
if (dodensup.eq.'y') then
    densup=.true.
else
    densup=.false.
end if
read (8,*) bla,dodirect
call to_lower(dodirect)
if (dodirect.eq.'y') then
    directscf=.true.
else
    directscf=.false.
end if
read (8,*) bla,dolevshft,lshft
call to_lower(dolevshft)
if (dolevshft.eq.'y') then
    levshft=.true.
else
    levshft=.false.
end if
read(8,*) bla,douhf,sz
call to_lower(douhf)
if (douhf.eq.'y') then
    uhf=.true.
else
    uhf=.false.
end if
read(8,*) bla,domp2
call to_lower(domp2)
if (domp2.eq.'y') then
    mp2energy=.true.
else
    mp2energy=.false.
end if
read(8,*) bla,docc
call to_lower(docc)
if (docc.eq.'y') then
    ccenergy=.true.
else
    ccenergy=.false.
end if
read(8,*) bla,prtlevl
read(8,*) bla,ncore
allocate(atlabels(natom),nuccharges(natom),atznumber(natom),atcoords(natom,3))
allocate(angcoords(natom,3))
write (9,*) "Atom labels, znuc, coords:"
do i=1,natom
   read(8,*) atlabels(i),angcoords(i,:)
   call to_lower(atlabels(i))
   call getnuccharge(atlabels(i),nuccharges(i),atznumber(i))
   write (9,'(a3,f5.1,3f12.6)') atlabels(i), nuccharges(i), angcoords(i,:)
end do
close(8)
atcoords=angcoords*bohr
allocate(r2AB(natom,natom))
r2AB=0._dp
do i=1,natom-1
  do j=1+1,natom
    r2AB(i,j)=sum((atcoords(i,:)-atcoords(j,:))**2)
    r2AB(j,i)=r2AB(i,j)
  end do
end do
if (uhf) then
  nocca=(sum(atznumber)+sz)/2
  noccb=nocca-sz
else
  nocc=sum(atznumber)/2
end if

! Do some checks for non-implemented functionality
if (uhf.and.directscf) then
  write (*,*) "Direct SCF not yet implemented with Direct SCF"
  STOP
end if
if ((mp2energy.or.ccenergy).and.directscf) then
  write (*,*) "Correlation methods not yet compatible with Direct SCF"
  STOP
end if

return

end subroutine read_infile

subroutine to_lower(str)
use nrtype
implicit none

character(*), intent(in out) :: str
integer(i4b) :: i

do i = 1, len(str)
    select case(str(i:i))
      case("A":"Z")
        str(i:i) = achar(iachar(str(i:i))+32)
    end select
end do  
end subroutine To_Lower

subroutine getnuccharge(lab,z,zint)
use nrtype 
implicit none

character(len=2), intent(in) :: lab
integer(i4b), intent(out) :: zint
real(dp), intent(out) :: z

integer(i4b), parameter :: nelement=36
integer(i4b) :: i
character(len=2) :: pertab(nelement)

pertab = (/ "h ","he","li","be","b ","c ","n ","o ","f ","ne", &
         &  "na","mg","al","si","p ","s ","cl","ar", &
         &  "k ","ca","sc","ti","v ","cr","mn","fe","co","ni","cu","zn", &
         &  "ga","ge","as","se","br","kr" /)

z=0
do i=1,nelement
   if (trim(lab).eq.trim(pertab(i))) then
       zint=i
       z=real(i,dp)
       exit
   end if
end do

if (z.eq.0) then
    write (*,*) "element not recognized!"
    stop
end if

end subroutine getnuccharge

