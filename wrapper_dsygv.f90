subroutine wrapper_dsygv(n,a,b,eigs,eigv)
use nrtype
implicit none

integer(i4b), intent(in) :: n
real(dp), intent(in) :: a(n,n), b(n,n)
real(dp), intent(out) :: eigs(n),eigv(n,n)

integer(i4b) :: itype, lda, lwork, info
character(len=1) :: jobz,uplo
real(dp) :: mat1(n,n), mat2(n,n), work(n*8)

itype=1   ! Solve generalized problem AX = lBX
jobz='V'  ! get eigenvectors
uplo='U'  ! upper triangle given on input (actually both are...)
lwork=8*n 
mat1=a
mat2=b
lda=n

call dsygv(itype,jobz,uplo,n,mat1,lda,mat2,lda,eigs,work,lwork,info)
eigv=mat1

return

end subroutine wrapper_dsygv


subroutine wrapper_dsyev(n,a,eigs,eigv)
use nrtype
implicit none

integer(i4b), intent(in) :: n
real(dp), intent(in) :: a(n,n)
real(dp), intent(out) :: eigs(n),eigv(n,n)

integer(i4b) :: lda, lwork, info
character(len=1) :: jobz,uplo
real(dp) :: mat(n,n), work(n*8)

jobz='V'  ! get eigenvectors
uplo='U'  ! upper triangle given on input (actually both are...)
lwork=8*n 
mat=a
lda=n

call dsyev(jobz,uplo,n,mat,lda,eigs,work,lwork,info)
eigv=mat

return

end subroutine wrapper_dsyev


subroutine wrapper_dsysv(n,a,b,x)
use nrtype
implicit none

integer(i4b), intent(in) :: n
real(dp), intent(in) :: a(n,n), b(n)
real(dp), intent(out) :: x(n)

integer(i4b) :: lda, lwork, info, ipiv(n)
character(len=1) :: uplo
real(dp) :: mat(n,n), bmat(n,1), work(8*n)

uplo='U'  ! upper triangle given on input (actually both are...)
lwork=8*n 
mat=a
lda=n
bmat(:,1)=b

call dsysv(uplo,n,1,mat,lda,ipiv,bmat,n,work,lwork,info)
x=bmat(:,1)

return

end subroutine wrapper_dsysv


