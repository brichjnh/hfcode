subroutine dealloc_hf_integs()
use nrtype; use molprops
implicit none

deallocate(s1e,t1e,v1e)
deallocate(pab,pm1,hpm1,qab,Px,kab,rab2,schw,ddab)

end subroutine dealloc_hf_integs

subroutine dealloc_molec_stuff()
use nrtype; use molprops
implicit none


deallocate(nuccharges,atznumber,atcoords,angcoords)
deallocate(basels,bastypes,basnshell,basnspd,basnzet,baszeta,basdjk,basdjkl)
deallocate(nzet,xsh,atsh,shtyp,zeta,djk,djkl,aosh,nshat,shat,nshtat,shtyat)
deallocate(atao,aotyp,aohuck,ssh,psh,dsh,lsh,shao)
deallocate(ssprs,spprs,sdprs,slprs,ppprs,pdprs,plprs,ddprs,dlprs,llprs)
deallocate(ioff,ind2,revind2)
deallocate(schwss,schwsp,schwsd,schwsl,schwpp)
deallocate(schwpd,schwpl,schwdd,schwdl,schwll)
deallocate(orderss,ordersp,ordersd,ordersl,orderpp)
deallocate(orderpd,orderpl,orderdd,orderdl,orderll)
deallocate(gammrats)
deallocate(boysfuncargs,boysfuncvals)

end subroutine dealloc_molec_stuff


subroutine dealloc_correl()
use nrtype; use molprops
implicit none

deallocate(moen)
deallocate(m1e,m2e)

end subroutine dealloc_correl


