subroutine mp2()
use nrtype ; use molprops
implicit none

integer(i4b) :: i,j,k,l
real(dp) :: fac1, fac2

emp2=0._dp

do i=ncore+1,nocc
   do j=ncore+1,nocc
      do k=nocc+1,nb
          do l=nocc+1,nb
              fac1=moen(i)+moen(j)-moen(k)-moen(l)
              fac2=2._dp*m2e(i,k,j,l)*m2e(k,i,l,j)-m2e(i,k,j,l)*m2e(k,j,l,i)
              emp2=emp2+fac2/fac1
          end do
      end do
   end do
end do

write (9,*) ""
write (9,*) "MP2 energy correction: Delta EMP2 and EMP2"
write (9,'(2F18.10)') emp2,etot+emp2

end subroutine mp2
