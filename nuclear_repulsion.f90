subroutine nuclear_repulsion()
use nrtype; use molprops
implicit none

integer(i4b) :: i, j
real(dp) :: rnm(3), r

unuc=0._dp
do i=1, natom
   do j = i+1, natom
       rnm=atcoords(i,:)-atcoords(j,:)
       r=sqrt(sum(rnm*rnm))
       unuc=unuc+nuccharges(i)*nuccharges(j)/r
   end do
end do
write (9,'(A,F15.7)') "Nuclear repulsion energy:",unuc
write (9,'(A,2I6)') "Number of occ. orbs/ number of electrons:",nocc,2*nocc
return

end subroutine nuclear_repulsion
