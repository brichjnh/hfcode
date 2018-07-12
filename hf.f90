program hf 
use nrtype; use molprops
implicit none

real(dp), parameter :: eps=1.d-12
character(len=30) :: arg(4), infile, basfile, outfile
integer(i4b) :: i

do i=1,3
   CALL get_command_argument(i, arg(i))
   IF (LEN_TRIM(arg(i)) == 0) then
      write (*,*) "You need to input infile basfile outfile."
      stop
   end if
end do

infile=arg(1)
basfile=arg(2)
outfile=arg(3)
open(unit=9,file=outfile)

call time_checker(0," ")

call read_infile(infile)
call nuclear_repulsion()
call read_basfile(basfile)
call spread_basis()

call time_checker(-1,"Starting Hartree-Fock program. Preambles:")

call precalcs_kabs()
call time_checker(-1,"Precalculation of various shell pair properties, timing:")
call cmpt1e()
call time_checker(-1,"One-electron integrals timing:")

if (directscf) then
   call direct_hartreefock()
   call time_checker(-1,"Direct Hartree-Fock timing:")
else
   call cmpt2e()
   call time_checker(-1,"Two-electron integrals timing:")
   if (uhf) then
     call u_incore_hartreefock()
     call time_checker(-1,"Unrestricted Hartree-Fock timing:")
   else
     call incore_hartreefock()   
     call time_checker(-1,"Restricted Conventional Hartree-Fock timing:")
     if (mp2energy.or.ccenergy) then
       call transform1e()
       call time_checker(-1,"One electron integral transform timing:")
       call transform2e()
       call time_checker(-1,"Two electron integral transform timing:")
       call calcmoenergies()
       call mp2()
       call time_checker(-1,"MP2 timing:")
       if (ccenergy) then
          call ccsd()
          call time_checker(-1,"CCSD timing:")
       end if
     end if
   end if
end if
call dealloc_molec_stuff()

call time_checker(1,"Final overall timing: ")
close(9)
end program hf

