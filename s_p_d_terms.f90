module s_p_d_terms
use nrtype ; use molprops
implicit none

contains

pure subroutine i2e_sssd(sh1,sh2,sh3,sh4,esssd)
use nrtype ; use molprops ; use boys
implicit none

! This version exploits the horizontal relation, eq. 18 of HG & P 1988

integer(i4b), intent(in) :: sh1,sh2,sh3,sh4
real(dp), intent(out) :: esssd(5)

integer(i4b) :: ab, cd, za, zb, zc, zd, zza, zzc, zzab, zzcd, m, lp1, l, lm1, jj, jrec
real(dp) :: ooppq, qpq, poppq, rpq2, fac1,afac, wx(3)
real(dp) :: bfac(0:2), os(10,0:2), osc(10), rdq(3), rqw(3)

! First assemble the [ss|sx] (x=s,p,d) integrals and contract them

osc=0._dp

ab=ind2(sh1,sh2)
cd=ind2(sh3,sh4)
do za=1,nzet(sh1)
  zza=maxprim*(za-1)
  do zb=1,nzet(sh2)
    zzab=zza+zb
    if (negab(zzab,ab)) cycle
    do zc=1,nzet(sh3)
      zzc=maxprim*(zc-1)
      do zd=1,nzet(sh4)
        zzcd=zzc+zd
        if (negab(zzcd,cd)) cycle
        ooppq=1._dp/(pab(zzab,ab)+pab(zzcd,cd))
        poppq=pab(zzab,ab)*ooppq
        qpq=pab(zzcd,cd)*poppq
        rpq2=sum((Px(:,zzab,ab)-Px(:,zzcd,cd))**2)
        Wx=(pab(zzab,ab)*Px(:,zzab,ab)+pab(zzcd,cd)*Px(:,zzcd,cd))*ooppq
        rdq=Px(:,zzcd,cd)-xsh(:,sh4)
        rqw=Wx-Px(:,zzcd,cd)
        ! Build (ss|ss)^m, m=0,1,2
        !call repulsion_f_vector(2,qpq*rpq2,bfac)
        call boys_2(qpq*rpq2,bfac)
        os(1,:)=kab(zzab,ab)*kab(zzcd,cd)*bfac*sqrt(ooppq)
        ! OS recursion: Build (ss|sp)^m, m=0,1
        do m=0,1
          os(2:4,m)=rdq(:)*os(1,m)+rqw(:)*os(1,m+1)
        end do
        afac=(os(1,0)-poppq*os(1,1))*hpm1(zzcd,cd)
        do lp1=5,10
          jrec=min1b(1,lp1)
          l=min1b(2,lp1)
          os(lp1,0)=os(l,0)*rdq(jrec)+os(l,1)*rqw(jrec)
          lm1=minus1(jrec,l)
          if (lm1.eq.0) cycle
          os(lp1,0)=os(lp1,0)+afac
        end do
        osc=osc+os(:,0)
      end do
    end do
  end do
end do


! now contract to form spherical harmonics. Rem: Cartesians:
!     0    1    2    3    4    5    6    7    8    9   10
!     /    s    x    y    z   xx   xy   xz   yy   yz   zz
! Spherical harmonics:
!     1    2    3    4    5
!    z2   xz   yz   xy  x2-y2


do jj=1,5
  esssd(jj)=sum(matharmond(jj,:)*osc(5:10))
end do

return

end subroutine i2e_sssd



pure subroutine i2e_ssdd(sh1,sh2,sh3,sh4,essdd)
use nrtype ; use molprops ; use boys
implicit none

! This version exploits the horizontal relation, eq. 18 of HG & P 1988

integer(i4b), intent(in) :: sh1,sh2,sh3,sh4
real(dp), intent(out) :: essdd(5,5)

integer(i4b) :: ab, cd, za, zb, zc, zd, zza, zzc, zzab, zzcd, m, lp1, l, lm1
integer(i4b) :: kp1, k, jj, kk, jrec
real(dp) :: ooppq, qpq, ooq, oop, oo2q, poppq, rpq2, fac1,afac(0:2), wx(3), matl(6,5)
real(dp) :: bfac(0:4), os(35,0:4), osc(35), rcq(3), rdq(3), rqw(3), hgp(10,35), rcd(3)

! First assemble the [ss|sx] (x=s,p,d) integrals and contract them

osc=0._dp

ab=ind2(sh1,sh2)
cd=ind2(sh3,sh4)
do za=1,nzet(sh1)
  zza=maxprim*(za-1)
  do zb=1,nzet(sh2)
    zzab=zza+zb
    if (negab(zzab,ab)) cycle
    do zc=1,nzet(sh3)
      zzc=maxprim*(zc-1)
      do zd=1,nzet(sh4)
        zzcd=zzc+zd
        if (negab(zzcd,cd)) cycle
        ooppq=1._dp/(pab(zzab,ab)+pab(zzcd,cd))
        poppq=pab(zzab,ab)*ooppq
        qpq=pab(zzcd,cd)*poppq
        rpq2=sum((Px(:,zzab,ab)-Px(:,zzcd,cd))**2)
        Wx=(pab(zzab,ab)*Px(:,zzab,ab)+pab(zzcd,cd)*Px(:,zzcd,cd))*ooppq
        rcq=Px(:,zzcd,cd)-xsh(:,sh3)
        rdq=Px(:,zzcd,cd)-xsh(:,sh4)
        rqw=Wx-Px(:,zzcd,cd)
        ! Build (ss|ss)^m, m=0,1,2
        !call repulsion_f_vector(4,qpq*rpq2,bfac)
        call boys_n(4,qpq*rpq2,bfac)
        os(1,:)=kab(zzab,ab)*kab(zzcd,cd)*bfac*sqrt(ooppq)
        ! OS recursion: Build (ss|sp)^m, m=0,1
        do m=0,3
          os(2:4,m)=rdq(:)*os(1,m)+rqw(:)*os(1,m+1)
        end do
        afac(0:2)=(os(1,0:2)-poppq*os(1,1:3))*hpm1(zzcd,cd)
        do lp1=5,10
          jrec=min1b(1,lp1)
          l=min1b(2,lp1)
          os(lp1,0:2)=os(l,0:2)*rdq(jrec)+os(l,1:3)*rqw(jrec)
          lm1=minus1(jrec,l)
          if (lm1.eq.0) cycle
          os(lp1,0:2)=os(lp1,0:2)+afac(0:2)
        end do
        do lp1=11,20
          jrec=min1b(1,lp1)
          l=min1b(2,lp1)
          os(lp1,0:1)=os(l,0:1)*rdq(jrec)+os(l,1:2)*rqw(jrec)
          lm1=minus1(jrec,l)
          if (lm1.eq.0) cycle
          os(lp1,0:1)=os(lp1,0:1)+indexps(jrec,l)*hpm1(zzcd,cd)*(os(lm1,0:1)-poppq*os(lm1,1:2))
        end do
        do lp1=21,35
          jrec=min1b(1,lp1)
          l=min1b(2,lp1)
          os(lp1,0)=os(l,0)*rdq(jrec)+os(l,1)*rqw(jrec)
          lm1=minus1(jrec,l)
          if (lm1.eq.0) cycle
          os(lp1,0:1)=os(lp1,0:1)+indexps(jrec,l)*hpm1(zzcd,cd)*(os(lm1,0)-poppq*os(lm1,1))
        end do
        osc=osc+os(:,0)
      end do
    end do
  end do
end do


! Now we have (ss|ss), (ss|sp), (ss|sd), (ss|sf) and (ss|sg) integrals
! We use:
! (ss|dijdkl)=(ss|pifjkl)+rcd(j)*(ss|pidkl)  for this we need (ss|pf) and (ss|pd)
! (ss|pifjkl)=(ss|sgijkl)+rcd(i)*(ss|sfjkl)      for this I need (ss|sg) and (ss|sf)
! (ss|pidkl)=(ss|sfikl)+rcd(i)*(ss|sdkl)      for this I need (ss|sf) and (ss|sd)

! Start by moving the (ss|sd), (ss|sf) and (ss|sg) integrals
hgp(1,5:35)=osc(5:35)
rcd=xsh(:,sh4)-xsh(:,sh3)
! Now build the (ss|pf) and (ss|pd)
do l=5,20   ! covering (ss|sd) and (ss|sf), with lp1 also covering (ss|sg)
  do jj=1,3
    kp1=plus1(jj,1)
    lp1=plus1(jj,l)
    hgp(kp1,l)=hgp(1,lp1)+rcd(jj)*hgp(1,l)
  end do
end do
! Now build the (ss|dd) 
do k=2,4
  do l=5,10
    do jj=1,3
      kp1=plus1(jj,k)
      lp1=plus1(jj,l)
      hgp(kp1,l)=hgp(k,lp1)+rcd(jj)*hgp(k,l)
    end do
  end do
end do

! Finally do contraction to form spherical harmonics

do jj=1,6
  do kk=1,5
    matl(jj,kk)=sum(matharmond(kk,:)*hgp(jj+4,5:10))
  end do
end do
do jj=1,5
  do kk=1,5
    essdd(jj,kk)=sum(matharmond(jj,:)*matl(:,kk))
  end do
end do

end subroutine i2e_ssdd





pure subroutine i2e_sdsd(sh1,sh2,sh3,sh4,esdsd)
use nrtype ; use molprops ; use boys
implicit none

integer(i4b), intent(in) :: sh1,sh2,sh3,sh4
real(dp), intent(out) :: esdsd(5,5)

integer(i4b) :: j,l,jm1,lm1, lp1, jj, jp1, kk, jrec
integer(i4b) :: ab, cd, za, zb, zc, zd, zza, zzc, zzab, zzcd, m
real(dp) :: ooppq, qpq, ooq, oop, oo2q, poppq, qoppq, rpq2, fac1, wx(3)
real(dp) :: bfac(0:4), os(10,10,0:4), rbp(3), rpw(3), rdq(3), rqw(3)
real(dp) :: afac(0:2), matl(6,5)
real(dp) :: osc(6,6)

! This routine uses O-S recursion

osc=0._dp

ab=ind2(sh1,sh2)
cd=ind2(sh3,sh4)
do za=1,nzet(sh1)
   zza=maxprim*(za-1)
   do zb=1,nzet(sh2)
      zzab=zza+zb
      if (negab(zzab,ab)) cycle
      rbp=Px(:,zzab,ab)-xsh(:,sh2)
      do zc=1,nzet(sh3)
         zzc=maxprim*(zc-1)
         do zd=1,nzet(sh4)
            zzcd=zzc+zd
            if (negab(zzcd,cd)) cycle
            ooppq=1._dp/(pab(zzab,ab)+pab(zzcd,cd))
            poppq=pab(zzab,ab)*ooppq
            qoppq=pab(zzcd,cd)*ooppq
            qpq=pab(zzab,ab)*qoppq
            rpq2=sum((Px(:,zzab,ab)-Px(:,zzcd,cd))**2)
            Wx=(pab(zzab,ab)*Px(:,zzab,ab)+pab(zzcd,cd)*Px(:,zzcd,cd))*ooppq
            rdq=Px(:,zzcd,cd)-xsh(:,sh4)
            rpw=Wx-Px(:,zzab,ab)
            rqw=Wx-Px(:,zzcd,cd)
            ! Build (ss|ss)^m, m=0,1,2,3,4
            !call repulsion_f_vector(4,qpq*rpq2,bfac)
            call boys_n(4,qpq*rpq2,bfac)
            os(1,1,:)=kab(zzab,ab)*kab(zzcd,cd)*bfac*sqrt(ooppq)
            ! OS recursion: Build (ss|sp)^m, m=0,3
            do m=0,3
              os(1,2:4,m)=rdq(:)*os(1,1,m)+rqw(:)*os(1,1,m+1)
            end do
            ! OS recursion: Build (ss|sd)^m, m=0,1,2
            afac(0:2)=(os(1,1,0:2)-poppq*os(1,1,1:3))*hpm1(zzab,ab)
            do lp1=5,10
              jrec=min1b(1,lp1)
              l=min1b(2,lp1)
              os(1,lp1,0:2)=os(1,l,0:2)*rdq(jrec)+os(1,l,1:3)*rqw(jrec)
              lm1=minus1(jrec,l)
              if (lm1.eq.0) cycle
              os(1,lp1,0:2)=os(1,lp1,0:2)+afac(0:2)
            end do
            ! OS recursion: Build (sp|sd)^m, m=0,1
            do jj=1,3
              j=jj+1
              os(j,1:10,0:1)=os(1,1:10,0:1)*rbp(jj)+os(1,1:10,1:2)*rpw(jj)
              do l=2,10
                lm1=minus1(jj,l)
                if (lm1.eq.0) cycle
                os(j,l,0:1)=os(j,l,0:1)+indexps(jj,l)*os(1,lm1,1:2)*.5_dp*ooppq
              end do
            end do
            ! OS recursion: Build (sd|sd)^m, m=0
            do jp1=5,10
              jrec=min1b(1,jp1)
              j=min1b(2,jp1)
              os(jp1,1:10,0)=os(j,1:10,0)*rbp(jrec)+os(j,1:10,1)*rpw(jrec)
              jm1=minus1s(jrec,j)
              if (jm1.ne.0) then
                os(jp1,1:10,0)=os(jp1,1:10,0)+hpm1(zzab,ab)*(os(jm1,:,0)-qoppq*os(jm1,:,1))
              end if
              do l=5,10
                lm1=minus1(jrec,l)
                if (lm1.eq.0) cycle
                os(jp1,l,0)=os(jp1,l,0)+os(j,lm1,1)*.5_dp*ooppq*indexps(jrec,l)
              end do
            end do
            osc=osc+os(5:10,5:10,0)
         end do
      end do
   end do
end do

! Finally do contraction to form spherical harmonics
! This is not very efficient but very simple

do jj=1,6
  do kk=1,5
    matl(jj,kk)=sum(matharmond(kk,:)*osc(jj,:))
  end do
end do
do jj=1,5
  do kk=1,5
    esdsd(jj,kk)=sum(matharmond(jj,:)*matl(:,kk))
  end do
end do

end subroutine i2e_sdsd



pure subroutine i2e_sddd(sh1,sh2,sh3,sh4,esddd)
use nrtype ; use molprops ; use boys
implicit none

integer(i4b), intent(in) :: sh1,sh2,sh3,sh4
real(dp), intent(out) :: esddd(5,5,5)

integer(i4b) :: j,l,jm1,lm1,lp1, k, kp1, bigL, jj, jp1, li, lf, jrec
integer(i4b) :: ab, cd, za, zb, zc, zd, zza, zzc, zzab, zzcd, m
real(dp) :: ooppq, qpq, qop, poppq, rpq2
real(dp) :: bfac(0:6), osv(84,0:6),oset(10,84), rdq(3), rqw(3)
real(dp) :: osc(6,35), hgp(6,10,35), sdddc(6,6,6), tmpl(6,6,5), tmpkl(6,5,5)
real(dp) :: afac(0:4)
real(dp) :: wx(3), rab(3), rcd(3), delv(3)

! Use the VRR to make [ss|si] (l=6), then ET to make [sd|sg], then HRR

osc=0._dp
oset=0._dp
hgp=0._dp

ab=ind2(sh1,sh2)
cd=ind2(sh3,sh4)
rcd=xsh(:,sh4)-xsh(:,sh3)
rab=xsh(:,sh2)-xsh(:,sh1)
do za=1,nzet(sh1)
  zza=maxprim*(za-1)
  do zb=1,nzet(sh2)
    zzab=zza+zb
    if (negab(zzab,ab)) cycle
    do zc=1,nzet(sh3)
      zzc=maxprim*(zc-1)
      do zd=1,nzet(sh4)
        zzcd=zzc+zd
        if (negab(zzcd,cd)) cycle
        ooppq=1._dp/(pab(zzab,ab)+pab(zzcd,cd))
        poppq=pab(zzab,ab)*ooppq
        qpq=pab(zzcd,cd)*poppq
        rpq2=sum((Px(:,zzab,ab)-Px(:,zzcd,cd))**2)
        Wx=(pab(zzab,ab)*Px(:,zzab,ab)+pab(zzcd,cd)*Px(:,zzcd,cd))*ooppq
        rdq=Px(:,zzcd,cd)-xsh(:,sh4)
        rqw=Wx-Px(:,zzcd,cd)
        ! Build (ss|ss)^m, m=0,1,2,3,4,5,6
        !call repulsion_f_vector(6,qpq*rpq2,bfac)
        call boys_n(6,qpq*rpq2,bfac)
        osv(1,:)=kab(zzab,ab)*kab(zzcd,cd)*bfac*sqrt(ooppq)
        ! OS recursion: Build (ss|sp)^m, m=0-5
        do m=0,5
          osv(2:4,m)=rdq(:)*osv(1,m)+rqw(:)*osv(1,m+1)
        end do
        ! OS recursion: Build (ss|sd)^m, m=0-4
        afac(0:4)=(osv(1,0:4)-poppq*osv(1,1:5))*hpm1(zzcd,cd)
        do lp1=5,10
          jrec=min1b(1,lp1)
          l=min1b(2,lp1)
          osv(lp1,0:4)=osv(l,0:4)*rdq(jrec)+osv(l,1:5)*rqw(jrec)
          lm1=minus1(jrec,l)
          if (lm1.eq.0) cycle
          osv(lp1,0:4)=osv(lp1,0:4)+afac(0:4)
        end do
        ! With my new notation, I can write a loop all the way up to (ss|si)
        lf=10
        do bigL=3,6
          li=lf ; lf=li+(bigL+1)*(bigL+2)/2  ! These go 10-20, 20-35, 35-56 and 56-84
          m=6-bigL
          do lp1=li+1,lf
            jrec=min1b(1,lp1)
            l=min1b(2,lp1)
            osv(lp1,0:m)=osv(l,0:m)*rdq(jrec)+osv(l,1:m+1)*rqw(jrec)
            lm1=minus1(jrec,l)
            if (lm1.eq.0) cycle
            osv(lp1,0:m)=osv(lp1,0:m)+real(indexps(jrec,l),dp)*hpm1(zzcd,cd)* &
                   & (osv(lm1,0:m)-poppq*osv(lm1,1:m+1))
          end do
        end do
        qop=pab(zzcd,cd)*pm1(zzab,ab)
        delv=pm1(zzab,ab)*(zeta(za,sh1)*rab(:)+zeta(zc,sh3)*rcd(:))
        oset(1,:)=osv(:,0)
        ! Now build up the bra terms. For the HR, I need [sd|sd], [sd|sf] and [sd|sg]
        ! Equation (45) in Giese's HSERILib document becomes here:
        ! [sd|sg] = del_i* [sp|sg] -qop* [sp|sh] + n_i(g)*oo2p* [sp|sf] + n_i(p)*oo2p* [ss|sg]
        ! [sd|sf] = del_i* [sp|sf] -qop* [sp|sg] + n_i(f)*oo2p* [sp|sd] + n_i(p)*oo2p* [ss|sf]
        ! [sd|sd] = del_i* [sp|sd] -qop* [sp|sf] + n_i(d)*oo2p* [sp|sp] + n_i(p)*oo2p* [ss|sd]
        !   These in turn require that I have [sp|sp], [sp|sd], [sp|sf], [sp|sg] and [sp|sh]
        ! [sp|sp] = del_i*[ss|sp]-qop*[ss|sd]+n_i(p)*oo2p*[ss|ss]
        ! [sp|sd] = del_i*[ss|sd]-qop*[ss|sf]+n_i(d)*oo2p*[ss|sp]
        ! [sp|sf] = del_i*[ss|sf]-qop*[ss|sg]+n_i(f)*oo2p*[ss|sd]
        ! [sp|sg] = del_i*[ss|sg]-qop*[ss|sh]+n_i(g)*oo2p*[ss|sf]
        ! [sp|sh] = del_i*[ss|sh]-qop*[ss|si]+n_i(h)*oo2p*[ss|sg]
        !     to first make [sp|sp], [sp|sd], [sp|sf], [sp|sg] and [sp|sh]
        !   So I need to use [ss|ss], ...  [ss|si] from O-S VR
        do l=1,56
          do jj=1,3
            lm1=minus1(jj,l)
            lp1=plus1(jj,l)
            jp1=jj+1
            oset(jp1,l)=-delv(jj)*oset(1,l)-qop*oset(1,lp1)
            if (lm1.ne.0) then
              oset(jp1,l)=oset(jp1,l)+hpm1(zzab,ab)*indexps(jj,l)*oset(1,lm1)
            end if
          end do
        end do
        ! Then we assemble [sd|sd], [sd|sf] and [sd|sg] (g means up to 35)
        do l=1,35
          do jp1=5,10
            jrec=min1b(1,jp1)
            j=min1b(2,jp1)
            lp1=plus1(jrec,l) ! This will from 11 to 56
            jm1=minus1s(jrec,j) ! This will be either 0 or 1
            lm1=minus1(jrec,l) ! And this will go from 2 to 20
            oset(jp1,l)=-delv(jrec)*oset(j,l)-qop*oset(j,lp1)
            if (jm1.ne.0) then
              oset(jp1,l)=oset(jp1,l)+indexps(jrec,j)*hpm1(zzab,ab)*oset(jm1,l)
            end if
            if (lm1.ne.0) then
              oset(jp1,l)=oset(jp1,l)+indexps(jrec,l)*hpm1(zzab,ab)*oset(j,lm1)
            end if
          end do
        end do
        osc=osc+oset(5:10,1:35)
      end do
    end do
  end do
end do

! So now we have contracted (sd|ss), (sd|sp), ...(sd|sg) integrals
! Recursion eqs. in ket:
! (sd|dd) = (sd|pf) + (Di-Ci) * (sd|pd)
! (sd|pf) = (sd|sg) + (Di-Ci) * (sd|sf)
! (sd|pd) = (sd|sf) + (Di-Ci) * (sd|sd)
!   So the integrals I use as building blocks are (sd|sd), (sd|sf), and (sd|sg)

! Start by moving the (sg|sg) (ss|ss), (ss|sp), ...(ss|sg), (sp|ss) ... (sg,sg) integrals

hgp(:,1,1:35) = osc(:,1:35)
! Now build the (sd|pf) and (sd|pd)
do l=5,20
  do jj=1,3
    kp1=jj+1
    lp1=plus1(jj,l)
    hgp(:,kp1,l)=hgp(:,1,lp1)+rcd(jj)*hgp(:,1,l)
  end do
end do
! Now build the (sd|dd) 
do kp1=5,10
  do l=5,10
    do jj=1,3
      k=minus1s(jj,kp1)
      if (k.eq.0) cycle
      jrec=jj
      exit
    end do
    lp1=plus1(jrec,l)
    hgp(:,kp1,l)=hgp(:,k,lp1)+rcd(jj)*hgp(:,k,l)
  end do
end do


! Finally do contraction to form spherical harmonics, using temporary arrays tmp(6,6,6) and tmp2
! This is not very efficient but very simple
!First take only the Cartesians:
sdddc=hgp(:,5:10,5:10)

do l=1,5
  do k=1,6
    do j=1,6
      tmpl(j,k,l)=sum(matharmond(l,:)*sdddc(j,k,:))
    end do
  end do
end do
do l=1,5
  do k=1,5
    do j=1,6
      tmpkl(j,k,l)=sum(matharmond(k,:)*tmpl(j,:,l))
    end do
  end do
end do
do l=1,5
  do k=1,5
    do j=1,5
      esddd(j,k,l)=sum(matharmond(j,:)*tmpkl(:,k,l))
    end do
  end do
end do

return
end subroutine i2e_sddd






pure subroutine i2e_dddd(sh1,sh2,sh3,sh4,edddd)
use nrtype ; use molprops ; use boys
implicit none

integer(i4b), intent(in) :: sh1,sh2,sh3,sh4
real(dp), intent(out) :: edddd(5,5,5,5)

integer(i4b) :: j,l,jm1,lm1,li,lf, bigL, i, ip1, jj, jp1, k, kp1 , lp1
integer(i4b) :: ab, cd, za, zb, zc, zd, zza, zzc, zzab, zzcd, m, jrec
real(dp) :: ooppq, qpq, ooq, oop, qop, delv(3), poppq, rpq2, fac1
real(dp) :: rab(3), rcd(3), wx(3)
real(dp) :: bfac(0:8), osv(165,0:8),oset(35,165), rbp(3), rpw(3), rdq(3), rqw(3)
real(dp) :: osc(35,35), hgpk(5:35,10,5:35), hgpb(10,5:35,6,6), ddddc(6,6,6,6)
real(dp) :: afac(0:7), tmpl(6,6,6,5), tmpkl(6,6,5,5), tmpjkl(6,5,5,5)

! This is my first try at doing the Giese thing - first build [ss|sk] (k has L = 8)
!   then do electron transfer, then contract, then use the HGP HR

osc=0._dp

ab=ind2(sh1,sh2)
cd=ind2(sh3,sh4)
rab=xsh(:,sh2)-xsh(:,sh1)
rcd=xsh(:,sh4)-xsh(:,sh3)
do za=1,nzet(sh1)
  zza=maxprim*(za-1)
  do zb=1,nzet(sh2)
    zzab=zza+zb
    if (negab(zzab,ab)) cycle
    rbp=Px(:,zzab,ab)-xsh(:,sh2)
    do zc=1,nzet(sh3)
      zzc=maxprim*(zc-1)
      do zd=1,nzet(sh4)
        zzcd=zzc+zd
        if (negab(zzcd,cd)) cycle
        ooppq=1._dp/(pab(zzab,ab)+pab(zzcd,cd))
        poppq=pab(zzab,ab)*ooppq
        qpq=pab(zzcd,cd)*poppq
        rpq2=sum((Px(:,zzab,ab)-Px(:,zzcd,cd))**2)
        Wx=(pab(zzab,ab)*Px(:,zzab,ab)+pab(zzcd,cd)*Px(:,zzcd,cd))*ooppq
        rdq=Px(:,zzcd,cd)-xsh(:,sh4)
        rqw=Wx-Px(:,zzcd,cd)
        ! Build (ss|ss)^m, m=0,1,2,3,4,5,6,7,8
        !call repulsion_f_vector(8,qpq*rpq2,bfac)
        call boys_n(8,qpq*rpq2,bfac)
        osv(1,:)=kab(zzab,ab)*kab(zzcd,cd)*bfac*sqrt(ooppq)
        ! OS recursion: Build (ss|sp)^m, m=0-7
        do m=0,7
          osv(2:4,m)=rdq(:)*osv(1,m)+rqw(:)*osv(1,m+1)
        end do
        ! OS recursion: Build (ss|sd)^m, m=0-6
        afac(0:6)=(osv(1,0:6)-poppq*osv(1,1:7))*hpm1(zzcd,cd)
        do lp1=5,10
          jrec=min1b(1,lp1)
          l=min1b(2,lp1)
          osv(lp1,0:6)=osv(l,0:6)*rdq(jrec)+osv(l,1:7)*rqw(jrec)
          lm1=minus1(jrec,l)
          if (lm1.eq.0) cycle
          osv(lp1,0:6)=osv(lp1,0:6)+afac(0:6)
        end do
        ! With my new notation, I can write a loop all the way up to (ss|sk)
        ! OS recursion: Build (ss|sf)
        lf=10
        do bigL=3,8
          li=lf ; lf=li+(bigL+1)*(bigL+2)/2
          m=8-bigL
          do lp1=li+1,lf
            jrec=min1b(1,lp1)
            l=min1b(2,lp1)
            osv(lp1,0:m)=osv(l,0:m)*rdq(jrec)+osv(l,1:m+1)*rqw(jrec)
            lm1=minus1(jrec,l)
            if (lm1.eq.0) cycle
            osv(lp1,0:m)=osv(lp1,0:m)+real(indexps(jrec,l),dp)*hpm1(zzcd,cd)* &
                   & (osv(lm1,0:m)-poppq*osv(lm1,1:m+1))
          end do
        end do
        ! So now we have [ss|ss], [ss|sp] ... [ss|sk]
        ! We no longer need to track m, but do need to handle the angular momentum on b, 'k'.
        ! Start by moving to a different array, and computing some needed quantities
        qop=pab(zzcd,cd)*pm1(zzab,ab)
        delv=-pm1(zzab,ab)*(zeta(za,sh1)*rab(:)+zeta(zc,sh3)*rcd(:))
        oset(1,:)=osv(:,0)
        ! Now build up the bra terms. For the HR, I need (sd|sd), (sd|sf), (sd|sg), (sf|sd), (sf|sf), (sf|sg), (sg|sd), (sg|sf), (sg|sg) 
        ! Equation (45) in Giese's HSERILib document becomes here:
        ! [sg|sg] = del_i* [sf|sg] -qop* [sf|sh] + n_i(g)*oo2p* [sf|sf] + n_i(p)*oo2p* [sd|sg]
        !   basically reguiring [sf|sf], [sf|sg], [sf|sh]
        ! Following eq. (45), Giese says we need [sp|sp] to [sp|sj], [sd|sd] to [sd|si], [sf|sf] to [sf|sh]
        ! This assembles [sp|ss] to [sp|sj]
        do l=1,120
          do jj=1,3
            lm1=minus1(jj,l)
            lp1=plus1(jj,l)
            jp1=jj+1
            oset(jp1,l)=delv(jj)*oset(1,l)-qop*oset(1,lp1)
            if (lm1.ne.0) then
              oset(jp1,l)=oset(jp1,l)+hpm1(zzab,ab)*indexps(jj,l)*oset(1,lm1)
            end if
          end do
        end do
        ! This assembles [sd|sp] to [sd|si]. It requires [sp|ss]
        do l=1,84
          do jp1=5,10
            jrec=min1b(1,jp1)
            j=min1b(2,jp1)
            lp1=plus1(jrec,l)
            jm1=minus1(jrec,j)
            lm1=minus1(jrec,l)
            oset(jp1,l)=delv(jrec)*oset(j,l)-qop*oset(j,lp1)
            if (jm1.ne.0) then
              oset(jp1,l)=oset(jp1,l)+indexps(jrec,j)*hpm1(zzab,ab)*oset(jm1,l)
            end if
            if (lm1.ne.0) then
              oset(jp1,l)=oset(jp1,l)+indexps(jrec,l)*hpm1(zzab,ab)*oset(j,lm1)
            end if
          end do
        end do
        ! This assembles [sf|sd] to [sf|sh]. This requires among others [sd|sp]
        do l=1,56
          do jp1=11,20
            jrec=min1b(1,jp1)
            j=min1b(2,jp1)
            lp1=plus1(jrec,l)
            jm1=minus1s(jrec,j)
            lm1=minus1(jrec,l)
            oset(jp1,l)=delv(jrec)*oset(j,l)-qop*oset(j,lp1)
            if (jm1.ne.0) then
              oset(jp1,l)=oset(jp1,l)+indexps(jrec,j)*hpm1(zzab,ab)*oset(jm1,l)
            end if
            if (lm1.ne.0) then
              oset(jp1,l)=oset(jp1,l)+indexps(jrec,l)*hpm1(zzab,ab)*oset(j,lm1)
            end if
          end do
        end do
        ! Finally we assemble [sg|sf] and [sg|sg] (h means up to 35). This requires among others [sf,sd]
        do l=1,35
          do jp1=21,35
            jrec=min1b(1,jp1)
            j=min1b(2,jp1)
            lp1=plus1(jrec,l)
            jm1=minus1(jrec,j)
            lm1=minus1(jrec,l)
            oset(jp1,l)=delv(jrec)*oset(j,l)-qop*oset(j,lp1)
            if (jm1.ne.0) then
              oset(jp1,l)=oset(jp1,l)+indexps(jrec,j)*hpm1(zzab,ab)*oset(jm1,l)
            end if
            if (lm1.ne.0) then
              oset(jp1,l)=oset(jp1,l)+indexps(jrec,l)*hpm1(zzab,ab)*oset(j,lm1)
            end if
          end do
        end do
        osc=osc+oset(:,1:35)
      end do
    end do
  end do
end do

! So now we have (sp|sp) to (sp|sj), (sd|sd) to (sd|si), (sf|sf) to (sf|sh) and (sg|sg)

! For each of the bra and ket terms, the HRR becomes: (here for ket):
! |dd) = |pf) + (Di-Ci) * |pd)
! |pf) = |sg) + (Di-Ci) * |sf)
! |pd) = |sf) + (Di-Ci) * |sd)
! We need as input (sd|sd)*, (sd|sf)*, (sd|sg)*, (sf|sd), (sf|sf), (sf|sg), (sg|sd), (sg|sf), (sg|sg)

hgpk(5:35,1,5:35) = osc(5:35,5:35)
rcd=xsh(:,sh4)-xsh(:,sh3)
rab=xsh(:,sh2)-xsh(:,sh1)
! Now build the (sx|pf) and (sx|pd)
do l=5,20
  do jj=1,3
    kp1=1+jj
    lp1=plus1(jj,l)
    hgpk(5:35,kp1,l)=hgpk(5:35,1,lp1)+rcd(jj)*hgpk(5:35,1,l)
  end do
end do
! Now build the (sx|dd) 
do kp1=5,10
  do l=5,10
    jrec=min1b(1,kp1)
    k=min1b(2,kp1)
    lp1=plus1(jrec,l)
    hgpk(5:35,kp1,l)=hgpk(5:35,k,lp1)+rcd(jrec)*hgpk(5:35,k,l)
  end do
end do
! Now move the resulting (sd|dd), (sf|dd) and (sg|dd) ints
hgpb(1,5:35,:,:)=hgpk(5:35,5:10,5:10)
! Now build the (pf|dd) and (pd|dd)
do j=5,20
  do jj=1,3
    ip1=1+jj
    jp1=plus1(jj,j)
    hgpb(ip1,j,:,:)=hgpb(1,jp1,:,:)+rab(jj)*hgpb(1,j,:,:) 
  end do
end do
! Now build the (dd|dd) 
do ip1=5,10
  do j=5,10
    jrec=min1b(1,ip1)
    i=min1b(2,ip1)
    jp1=plus1(jrec,j)
    hgpb(ip1,j,:,:)=hgpb(i,jp1,:,:)+rab(jrec)*hgpb(i,j,:,:) 
  end do
end do


! Finally do contraction to form spherical harmonics, using temporary arrays tmp(6,6,6,6) and tmp2
! This is not very efficient but very simple
!First take only the Cartesians:
ddddc=hgpb(5:10,5:10,:,:)

! Then recontract
do l=1,5
  do k=1,6
    do j=1,6
      do i=1,6
        tmpl(i,j,k,l)=sum(matharmond(l,:)*ddddc(i,j,k,:))
      end do
    end do
  end do
end do
do l=1,5
  do k=1,5
    do j=1,6
      do i=1,6
        tmpkl(i,j,k,l)=sum(matharmond(k,:)*tmpl(i,j,:,l))
      end do
    end do
  end do
end do
do l=1,5
  do k=1,5
    do j=1,5
      do i=1,6
        tmpjkl(i,j,k,l)=sum(matharmond(j,:)*tmpkl(i,:,k,l))
      end do
    end do
  end do
end do
do l=1,5
  do k=1,5
    do j=1,5
      do i=1,5
        edddd(i,j,k,l)=sum(matharmond(i,:)*tmpjkl(:,j,k,l))
      end do
    end do
  end do
end do

return
end subroutine i2e_dddd


pure subroutine i2e_sspd(sh1,sh2,sh3,sh4,esspd)
use nrtype ; use molprops ; use boys
implicit none

! This version exploits the horizontal relation, eq. 18 of HG & P 1988

integer(i4b), intent(in) :: sh1,sh2,sh3,sh4
real(dp), intent(out) :: esspd(3,5)

integer(i4b) :: ab, cd, za, zb, zc, zd, zza, zzc, zzab, zzcd, m, lp1, l, lm1
integer(i4b) :: kp1, k, jj, kk, jrec
real(dp) :: ooppq, qpq, poppq, rpq2, fac1,afac(0:1), wx(3)
real(dp) :: bfac(0:3), os(20,0:3), osc(20), rdq(3), rqw(3), hgp(4,20), rcd(3)

! First assemble the [ss|sx] (x=s,p,d) integrals and contract them

osc=0._dp

ab=ind2(sh1,sh2)
cd=ind2(sh3,sh4)
do za=1,nzet(sh1)
  zza=maxprim*(za-1)
  do zb=1,nzet(sh2)
    zzab=zza+zb
    if (negab(zzab,ab)) cycle
    do zc=1,nzet(sh3)
      zzc=maxprim*(zc-1)
      do zd=1,nzet(sh4)
        zzcd=zzc+zd
        if (negab(zzcd,cd)) cycle
        ooppq=1._dp/(pab(zzab,ab)+pab(zzcd,cd))
        poppq=pab(zzab,ab)*ooppq
        qpq=pab(zzcd,cd)*poppq
        rpq2=sum((Px(:,zzab,ab)-Px(:,zzcd,cd))**2)
        Wx=(pab(zzab,ab)*Px(:,zzab,ab)+pab(zzcd,cd)*Px(:,zzcd,cd))*ooppq
        rdq=Px(:,zzcd,cd)-xsh(:,sh4)
        rqw=Wx-Px(:,zzcd,cd)
        ! Build (ss|ss)^m, m=0,1,2
        !call repulsion_f_vector(3,qpq*rpq2,bfac)
        call boys_n(3,qpq*rpq2,bfac)
        os(1,:)=kab(zzab,ab)*kab(zzcd,cd)*bfac*sqrt(ooppq)
        ! OS recursion: Build (ss|sp)^m, m=0,1
        do m=0,2
          os(2:4,m)=rdq(:)*os(1,m)+rqw(:)*os(1,m+1)
        end do
        afac(0:1)=(os(1,0:1)-poppq*os(1,1:2))*hpm1(zzcd,cd)
        do lp1=5,10
          jrec=min1b(1,lp1)
          l=min1b(2,lp1)
          os(lp1,0:1)=os(l,0:1)*rdq(jrec)+os(l,1:2)*rqw(jrec)
          lm1=minus1(jrec,l)
          if (lm1.eq.0) cycle
          os(lp1,0:1)=os(lp1,0:1)+afac(0:1)
        end do
        do lp1=11,20
          jrec=min1b(1,lp1)
          l=min1b(2,lp1)
          os(lp1,0)=os(l,0)*rdq(jrec)+os(l,1)*rqw(jrec)
          lm1=minus1(jrec,l)
          if (lm1.eq.0) cycle
          os(lp1,0)=os(lp1,0)+indexps(jrec,l)*hpm1(zzcd,cd)*(os(lm1,0)-poppq*os(lm1,1))
        end do
        osc=osc+os(:,0)
      end do
    end do
  end do
end do


! Now we have (ss|ss), (ss|sp), (ss|sd), (ss|sf) integrals
! We use:
! (ss|pd) = (ss|sf) + rcd*(ss|sd)

! Start by moving the (ss|sd), (ss|sf) and (ss|sg) integrals
hgp(1,:)=osc
rcd=xsh(:,sh4)-xsh(:,sh3)
do l=5,10   ! covering (ss|sd) and (ss|sf), with lp1 also covering (ss|sg)
  do jj=1,3
    kp1=jj+1
    lp1=plus1(jj,l)
    hgp(kp1,l)=hgp(1,lp1)+rcd(jj)*hgp(1,l)
  end do
end do

! Finally do contraction to form spherical harmonics

do jj=1,3
do kk=1,5
  esspd(jj,kk)=sum(matharmond(kk,:)*hgp(jj+1,5:10))
end do
end do

end subroutine i2e_sspd





pure subroutine i2e_spsd(sh1,sh2,sh3,sh4,espsd)
use nrtype ; use molprops ; use boys
implicit none

integer(i4b), intent(in) :: sh1,sh2,sh3,sh4
real(dp), intent(out) :: espsd(3,5)

integer(i4b) :: j,l,jm1,lm1, lp1, jj, jp1, kk, jrec
integer(i4b) :: ab, cd, za, zb, zc, zd, zza, zzc, zzab, zzcd, m
real(dp) :: ooppq, qpq, poppq, rpq2, fac1, wx(3)
real(dp) :: bfac(0:3), os(4,10,0:3), rbp(3), rpw(3), rdq(3), rqw(3)
real(dp) :: afac(0:1)
real(dp) :: osc(3,6)

! This routine uses O-S recursion

osc=0._dp

ab=ind2(sh1,sh2)
cd=ind2(sh3,sh4)
do za=1,nzet(sh1)
   zza=maxprim*(za-1)
   do zb=1,nzet(sh2)
      zzab=zza+zb
      if (negab(zzab,ab)) cycle
      rbp=Px(:,zzab,ab)-xsh(:,sh2)
      do zc=1,nzet(sh3)
         zzc=maxprim*(zc-1)
         do zd=1,nzet(sh4)
            zzcd=zzc+zd
            if (negab(zzcd,cd)) cycle
            ooppq=1._dp/(pab(zzab,ab)+pab(zzcd,cd))
            poppq=pab(zzab,ab)*ooppq
            qpq=pab(zzcd,cd)*poppq
            rpq2=sum((Px(:,zzab,ab)-Px(:,zzcd,cd))**2)
            Wx=(pab(zzab,ab)*Px(:,zzab,ab)+pab(zzcd,cd)*Px(:,zzcd,cd))*ooppq
            rdq=Px(:,zzcd,cd)-xsh(:,sh4)
            rpw=Wx-Px(:,zzab,ab)
            rqw=Wx-Px(:,zzcd,cd)
            ! Build (ss|ss)^m, m=0,1,2,3,4
            !call repulsion_f_vector(3,qpq*rpq2,bfac)
            call boys_n(3,qpq*rpq2,bfac)
            os(1,1,:)=kab(zzab,ab)*kab(zzcd,cd)*bfac*sqrt(ooppq)
            ! OS recursion: Build (ss|sp)^m, m=0,3
            do m=0,2
              os(1,2:4,m)=rdq(:)*os(1,1,m)+rqw(:)*os(1,1,m+1)
            end do
            ! OS recursion: Build (ss|sd)^m, m=0,1,2
            afac(0:1)=(os(1,1,0:1)-poppq*os(1,1,1:2))*hpm1(zzab,ab)
            do lp1=5,10
              jrec=min1b(1,lp1)
              l=min1b(2,lp1)
              os(1,lp1,0:1)=os(1,l,0:1)*rdq(jrec)+os(1,l,1:2)*rqw(jrec)
              lm1=minus1(jrec,l)
              if (lm1.eq.0) cycle
              os(1,lp1,0:1)=os(1,lp1,0:1)+afac(0:1)
            end do
            ! OS recursion: Build (sp|sd)^m, m=0,1
            do jj=1,3
              j=jj+1
              os(j,1:10,0)=os(1,1:10,0)*rbp(jj)+os(1,1:10,1)*rpw(jj)
              do l=2,10
                lm1=minus1(jj,l)
                if (lm1.eq.0) cycle
                os(j,l,0)=os(j,l,0)+indexps(jj,l)*os(1,lm1,1)*.5_dp*ooppq
              end do
            end do
            osc=osc+os(2:4,5:10,0)
         end do
      end do
   end do
end do

! Finally do contraction to form spherical harmonics
! This is not very efficient but very simple

do jj=1,3
do kk=1,5
  espsd(jj,kk)=sum(matharmond(kk,:)*osc(jj,:))
end do
end do

end subroutine i2e_spsd



pure subroutine i2e_sppd(sh1,sh2,sh3,sh4,esppd)
use nrtype ; use molprops ; use boys
implicit none

integer(i4b), intent(in) :: sh1,sh2,sh3,sh4
real(dp), intent(out) :: esppd(3,3,5)

integer(i4b) :: j,l,jm1,lm1,lp1, k, kp1, bigL, jj, jp1, li, lf, jrec
integer(i4b) :: ab, cd, za, zb, zc, zd, zza, zzc, zzab, zzcd, m
real(dp) :: ooppq, qpq, qop, poppq, rpq2
real(dp) :: bfac(0:4), osv(35,0:4),oset(4,35), rdq(3), rqw(3)
real(dp) :: osc(3,20), hgp(3,4,20), sppdc(3,3,6)
real(dp) :: afac(0:2)
real(dp) :: wx(3), rab(3), rcd(3), delv(3)

! Use the VRR to make [ss|si] (l=6), then ET to make [sd|sg], then HRR

osc=0._dp
oset=0._dp
hgp=0._dp

ab=ind2(sh1,sh2)
cd=ind2(sh3,sh4)
rcd=xsh(:,sh4)-xsh(:,sh3)
rab=xsh(:,sh2)-xsh(:,sh1)
do za=1,nzet(sh1)
  zza=maxprim*(za-1)
  do zb=1,nzet(sh2)
    zzab=zza+zb
    if (negab(zzab,ab)) cycle
    do zc=1,nzet(sh3)
      zzc=maxprim*(zc-1)
      do zd=1,nzet(sh4)
        zzcd=zzc+zd
        if (negab(zzcd,cd)) cycle
        ooppq=1._dp/(pab(zzab,ab)+pab(zzcd,cd))
        poppq=pab(zzab,ab)*ooppq
        qpq=pab(zzcd,cd)*poppq
        rpq2=sum((Px(:,zzab,ab)-Px(:,zzcd,cd))**2)
        Wx=(pab(zzab,ab)*Px(:,zzab,ab)+pab(zzcd,cd)*Px(:,zzcd,cd))*ooppq
        rdq=Px(:,zzcd,cd)-xsh(:,sh4)
        rqw=Wx-Px(:,zzcd,cd)
        !call repulsion_f_vector(4,qpq*rpq2,bfac)
        call boys_n(4,qpq*rpq2,bfac)
        osv(1,:)=kab(zzab,ab)*kab(zzcd,cd)*bfac*sqrt(ooppq)
        do m=0,3
          osv(2:4,m)=rdq(:)*osv(1,m)+rqw(:)*osv(1,m+1)
        end do
        ! OS recursion: Build (ss|sd)^m, m=0-4
        afac(0:2)=(osv(1,0:2)-poppq*osv(1,1:3))*hpm1(zzcd,cd)
        do lp1=5,10
          jrec=min1b(1,lp1)
          l=min1b(2,lp1)
          osv(lp1,0:2)=osv(l,0:2)*rdq(jrec)+osv(l,1:3)*rqw(jrec)
          lm1=minus1(jrec,l)
          if (lm1.eq.0) cycle
          osv(lp1,0:2)=osv(lp1,0:2)+afac(0:2)
        end do
        ! With my new notation, I can write a loop all the way up to (ss|si)
        lf=10
        do bigL=3,4
          li=lf ; lf=li+(bigL+1)*(bigL+2)/2  ! These go 10-20, 20-35, 35-56
          m=4-bigL
          do lp1=li+1,lf
            do jj=1,3
              l=minus1(jj,lp1)
              if(l.eq.0) cycle
              jrec=jj
              exit
            end do
            osv(lp1,0:m)=osv(l,0:m)*rdq(jrec)+osv(l,1:m+1)*rqw(jrec)
            lm1=minus1(jrec,l)
            if (lm1.eq.0) cycle
            osv(lp1,0:m)=osv(lp1,0:m)+real(indexps(jrec,l),dp)*hpm1(zzcd,cd)* &
                   & (osv(lm1,0:m)-poppq*osv(lm1,1:m+1))
          end do
        end do
        qop=pab(zzcd,cd)*pm1(zzab,ab)
        delv=pm1(zzab,ab)*(zeta(za,sh1)*rab(:)+zeta(zc,sh3)*rcd(:))
        oset(1,:)=osv(:,0)
        do l=1,20
          do jj=1,3
            lm1=minus1(jj,l)
            lp1=plus1(jj,l)
            jp1=jj+1
            oset(jp1,l)=-delv(jj)*oset(1,l)-qop*oset(1,lp1)
            if (lm1.ne.0) then
              oset(jp1,l)=oset(jp1,l)+hpm1(zzab,ab)*indexps(jj,l)*oset(1,lm1)
            end if
          end do
        end do
        osc=osc+oset(2:4,1:20)
      end do
    end do
  end do
end do

hgp(:,1,1:20) = osc(:,1:20)
! Build the (sp|pd)
do l=5,10
  do jj=1,3
    kp1=jj+1
    lp1=plus1(jj,l)
    hgp(:,kp1,l)=hgp(:,1,lp1)+rcd(jj)*hgp(:,1,l)
  end do
end do


! Finally do contraction to form spherical harmonics, using temporary arrays tmp(6,6,6) and tmp2
! This is not very efficient but very simple
!First take only the Cartesians:
sppdc=hgp(:,2:4,5:10)

do l=1,5
  do k=1,3
    do j=1,3
      esppd(j,k,l)=sum(matharmond(l,:)*sppdc(j,k,:))
    end do
  end do
end do

return
end subroutine i2e_sppd


pure subroutine i2e_sdpp(sh1,sh2,sh3,sh4,esdpp)
use nrtype ; use molprops ; use boys
implicit none

integer(i4b), intent(in) :: sh1,sh2,sh3,sh4
real(dp), intent(out) :: esdpp(5,3,3)

integer(i4b) :: j,l,jm1,ip1, k, kp1, bigL, jj, jp1, li, lf, lp1, lm1, i
integer(i4b) :: ab, cd, za, zb, zc, zd, zza, zzc, zzab, zzcd, m, jrec
real(dp) :: ooppq, qpq, qop, poppq, rpq2
real(dp) :: bfac(0:4), osv(35,0:4),oset(10,35), rdq(3), rqw(3)
real(dp) :: osc(10,10), hgp(6,4,10), sdppc(6,3,3)
real(dp) :: afac(0:2)
real(dp) :: wx(3), rab(3), rcd(3), delv(3)

! Use the VRR to make [ss|si] (l=6), then ET to make [sd|sg], then HRR

osc=0._dp
oset=0._dp
hgp=0._dp

ab=ind2(sh1,sh2)
cd=ind2(sh3,sh4)
rcd=xsh(:,sh4)-xsh(:,sh3)
rab=xsh(:,sh2)-xsh(:,sh1)
do za=1,nzet(sh1)
  zza=maxprim*(za-1)
  do zb=1,nzet(sh2)
    zzab=zza+zb
    if (negab(zzab,ab)) cycle
    do zc=1,nzet(sh3)
      zzc=maxprim*(zc-1)
      do zd=1,nzet(sh4)
        zzcd=zzc+zd
        if (negab(zzcd,cd)) cycle
        ooppq=1._dp/(pab(zzab,ab)+pab(zzcd,cd))
        poppq=pab(zzab,ab)*ooppq
        qpq=pab(zzcd,cd)*poppq
        rpq2=sum((Px(:,zzab,ab)-Px(:,zzcd,cd))**2)
        Wx=(pab(zzab,ab)*Px(:,zzab,ab)+pab(zzcd,cd)*Px(:,zzcd,cd))*ooppq
        rdq=Px(:,zzcd,cd)-xsh(:,sh4)
        rqw=Wx-Px(:,zzcd,cd)
        !call repulsion_f_vector(4,qpq*rpq2,bfac)
        call boys_n(4,qpq*rpq2,bfac)
        osv(1,:)=kab(zzab,ab)*kab(zzcd,cd)*bfac*sqrt(ooppq)
        do m=0,3
          osv(2:4,m)=rdq(:)*osv(1,m)+rqw(:)*osv(1,m+1)
        end do
        ! OS recursion: Build (ss|sd)^m, m=0-4
        afac(0:2)=(osv(1,0:2)-poppq*osv(1,1:3))*hpm1(zzcd,cd)
        do lp1=5,10
          do jj=1,3
            l=minus1(jj,lp1)
            if (l.eq.0) cycle
            jrec=jj
            exit
          end do
          osv(lp1,0:2)=osv(l,0:2)*rdq(jrec)+osv(l,1:3)*rqw(jrec)
          lm1=minus1(jrec,l)
          if (lm1.eq.0) cycle
          osv(lp1,0:2)=osv(lp1,0:2)+afac(0:2)
        end do
        ! With my new notation, I can write a loop all the way up to (ss|si)
        lf=10
        do bigL=3,4
          li=lf ; lf=li+(bigL+1)*(bigL+2)/2  ! These go 10-20, 20-35, 35-56
          m=4-bigL
          do lp1=li+1,lf
            do jj=1,3
              l=minus1(jj,lp1)
              if(l.eq.0) cycle
              jrec=jj
              exit
            end do
            osv(lp1,0:m)=osv(l,0:m)*rdq(jrec)+osv(l,1:m+1)*rqw(jrec)
            lm1=minus1(jrec,l)
            if (lm1.eq.0) cycle
            osv(lp1,0:m)=osv(lp1,0:m)+real(indexps(jrec,l),dp)*hpm1(zzcd,cd)* &
                   & (osv(lm1,0:m)-poppq*osv(lm1,1:m+1))
          end do
        end do
        qop=pab(zzcd,cd)*pm1(zzab,ab)
        delv=pm1(zzab,ab)*(zeta(za,sh1)*rab(:)+zeta(zc,sh3)*rcd(:))
        oset(1,:)=osv(:,0)
        do l=1,20
          do jj=1,3
            lm1=minus1(jj,l)
            lp1=plus1(jj,l)
            jp1=jj+1
            oset(jp1,l)=-delv(jj)*oset(1,l)-qop*oset(1,lp1)
            if (lm1.ne.0) then
              oset(jp1,l)=oset(jp1,l)+hpm1(zzab,ab)*indexps(jj,l)*oset(1,lm1)
            end if
          end do
        end do
        do l=1,10
          do jp1=5,10
            do jj=1,3
              j=minus1s(jj,jp1) ! This will range from 2 to 4
              if (j.eq.0) cycle
              jrec=jj
              exit
            end do
            lp1=plus1(jrec,l) ! This will from 11 to 56
            jm1=minus1s(jrec,j) ! This will be either 0 or 1
            lm1=minus1(jrec,l) ! And this will go from 2 to 20
            oset(jp1,l)=-delv(jrec)*oset(j,l)-qop*oset(j,lp1)
            if (jm1.ne.0) then
              oset(jp1,l)=oset(jp1,l)+indexps(jrec,j)*hpm1(zzab,ab)*oset(jm1,l)
            end if
            if (lm1.ne.0) then
              oset(jp1,l)=oset(jp1,l)+indexps(jrec,l)*hpm1(zzab,ab)*oset(j,lm1)
            end if
          end do
        end do
        osc=osc+oset(:,1:10)
      end do
    end do
  end do
end do

hgp(:,1,:) = osc(5:10,:)
! Build the (pp|sd)
do j=2,4
  do jj=1,3
    ip1=jj+1
    jp1=plus1(jj,j)
    hgp(:,ip1,j)=hgp(:,1,jp1)+rcd(jj)*hgp(:,1,j)
  end do
end do


! Finally do contraction to form spherical harmonics, using temporary arrays tmp(6,6,6) and tmp2
! This is not very efficient but very simple
!First take only the Cartesians:
sdppc=hgp(:,2:4,2:4)

do l=1,3
  do k=1,3
    do j=1,5
      esdpp(j,k,l)=sum(matharmond(j,:)*sdppc(:,k,l))
    end do
  end do
end do

return
end subroutine i2e_sdpp


pure subroutine i2e_spdd(sh1,sh2,sh3,sh4,espdd)
use nrtype ; use molprops ; use boys
implicit none

integer(i4b), intent(in) :: sh1,sh2,sh3,sh4
real(dp), intent(out) :: espdd(3,5,5)

integer(i4b) :: j,l,jm1,lm1,lp1, k, kp1, bigL, jj, jp1, li, lf, jrec
integer(i4b) :: ab, cd, za, zb, zc, zd, zza, zzc, zzab, zzcd, m
real(dp) :: ooppq, qpq, qop, poppq, rpq2
real(dp) :: bfac(0:5), osv(56,0:5),oset(4,56), rdq(3), rqw(3)
real(dp) :: osc(3,35), hgp(3,10,35), spddc(3,6,6), tmpl(3,6,5)
real(dp) :: afac(0:3)
real(dp) :: wx(3), rab(3), rcd(3), delv(3)

! Use the VRR to make [ss|si] (l=6), then ET to make [sd|sg], then HRR

osc=0._dp
oset=0._dp
hgp=0._dp

ab=ind2(sh1,sh2)
cd=ind2(sh3,sh4)
rcd=xsh(:,sh4)-xsh(:,sh3)
rab=xsh(:,sh2)-xsh(:,sh1)
do za=1,nzet(sh1)
  zza=maxprim*(za-1)
  do zb=1,nzet(sh2)
    zzab=zza+zb
    if (negab(zzab,ab)) cycle
    do zc=1,nzet(sh3)
      zzc=maxprim*(zc-1)
      do zd=1,nzet(sh4)
        zzcd=zzc+zd
        if (negab(zzcd,cd)) cycle
        ooppq=1._dp/(pab(zzab,ab)+pab(zzcd,cd))
        poppq=pab(zzab,ab)*ooppq
        qpq=pab(zzcd,cd)*poppq
        rpq2=sum((Px(:,zzab,ab)-Px(:,zzcd,cd))**2)
        Wx=(pab(zzab,ab)*Px(:,zzab,ab)+pab(zzcd,cd)*Px(:,zzcd,cd))*ooppq
        rdq=Px(:,zzcd,cd)-xsh(:,sh4)
        rqw=Wx-Px(:,zzcd,cd)
        !call repulsion_f_vector(5,qpq*rpq2,bfac)
        call boys_n(5,qpq*rpq2,bfac)
        osv(1,:)=kab(zzab,ab)*kab(zzcd,cd)*bfac*sqrt(ooppq)
        do m=0,4
          osv(2:4,m)=rdq(:)*osv(1,m)+rqw(:)*osv(1,m+1)
        end do
        ! OS recursion: Build (ss|sd)^m, m=0-4
        afac(0:3)=(osv(1,0:3)-poppq*osv(1,1:4))*hpm1(zzcd,cd)
        do lp1=5,10
          do jj=1,3
            l=minus1(jj,lp1)
            if (l.eq.0) cycle
            jrec=jj
            exit
          end do
          osv(lp1,0:3)=osv(l,0:3)*rdq(jrec)+osv(l,1:4)*rqw(jrec)
          lm1=minus1(jrec,l)
          if (lm1.eq.0) cycle
          osv(lp1,0:3)=osv(lp1,0:3)+afac(0:3)
        end do
        ! With my new notation, I can write a loop all the way up to (ss|si)
        lf=10
        do bigL=3,5
          li=lf ; lf=li+(bigL+1)*(bigL+2)/2  ! These go 10-20, 20-35, 35-56
          m=5-bigL
          do lp1=li+1,lf
            do jj=1,3
              l=minus1(jj,lp1)
              if(l.eq.0) cycle
              jrec=jj
              exit
            end do
            osv(lp1,0:m)=osv(l,0:m)*rdq(jrec)+osv(l,1:m+1)*rqw(jrec)
            lm1=minus1(jrec,l)
            if (lm1.eq.0) cycle
            osv(lp1,0:m)=osv(lp1,0:m)+real(indexps(jrec,l),dp)*.5_dp*pm1(zzcd,cd)* &
                   & (osv(lm1,0:m)-poppq*osv(lm1,1:m+1))
          end do
        end do
        qop=pab(zzcd,cd)*pm1(zzab,ab)
        delv=pm1(zzab,ab)*(zeta(za,sh1)*rab(:)+zeta(zc,sh3)*rcd(:))
        oset(1,:)=osv(:,0)
        do l=1,35
          do jj=1,3
            lm1=minus1(jj,l)
            lp1=plus1(jj,l)
            jp1=jj+1
            oset(jp1,l)=-delv(jj)*oset(1,l)-qop*oset(1,lp1)
            if (lm1.ne.0) then
              oset(jp1,l)=oset(jp1,l)+hpm1(zzab,ab)*indexps(jj,l)*oset(1,lm1)
            end if
          end do
        end do
        osc=osc+oset(2:4,1:35)
      end do
    end do
  end do
end do

hgp(:,1,1:35) = osc(:,1:35)
! Now build the (sd|pf) and (sd|pd)
do l=5,20
  do jj=1,3
    kp1=jj+1
    lp1=plus1(jj,l)
    hgp(:,kp1,l)=hgp(:,1,lp1)+rcd(jj)*hgp(:,1,l)
  end do
end do
! Now build the (sd|dd) 
do kp1=5,10
  do l=5,10
    do jj=1,3
      k=minus1s(jj,kp1)
      if (k.eq.0) cycle
      jrec=jj
      exit
    end do
    lp1=plus1(jrec,l)
    hgp(:,kp1,l)=hgp(:,k,lp1)+rcd(jrec)*hgp(:,k,l)
  end do
end do


! Finally do contraction to form spherical harmonics, using temporary arrays tmp(6,6,6) and tmp2
! This is not very efficient but very simple
!First take only the Cartesians:
spddc=hgp(:,5:10,5:10)

do l=1,5
  do k=1,6
    do j=1,3
      tmpl(j,k,l)=sum(matharmond(l,:)*spddc(j,k,:))
    end do
  end do
end do
do l=1,5
  do k=1,5
    do j=1,3
      espdd(j,k,l)=sum(matharmond(k,:)*tmpl(j,:,l))
    end do
  end do
end do

return
end subroutine i2e_spdd





pure subroutine i2e_sdpd(sh1,sh2,sh3,sh4,esdpd)
use nrtype ; use molprops ; use boys
implicit none

integer(i4b), intent(in) :: sh1,sh2,sh3,sh4
real(dp), intent(out) :: esdpd(5,3,5)

integer(i4b) :: j,l,jm1,lm1,lp1, k, kp1, bigL, jj, jp1, li, lf, jrec
integer(i4b) :: ab, cd, za, zb, zc, zd, zza, zzc, zzab, zzcd, m
real(dp) :: ooppq, qpq, qop, poppq, rpq2
real(dp) :: bfac(0:5), osv(56,0:5),oset(10,56), rdq(3), rqw(3)
real(dp) :: osc(6,20), hgp(6,4,20), sdpdc(6,3,6), tmpl(6,3,5)
real(dp) :: afac(0:3)
real(dp) :: wx(3), rab(3), rcd(3), delv(3)

! Use the VRR to make [ss|si] (l=6), then ET to make [sd|sg], then HRR

osc=0._dp
oset=0._dp
hgp=0._dp

ab=ind2(sh1,sh2)
cd=ind2(sh3,sh4)
rcd=xsh(:,sh4)-xsh(:,sh3)
rab=xsh(:,sh2)-xsh(:,sh1)
do za=1,nzet(sh1)
  zza=maxprim*(za-1)
  do zb=1,nzet(sh2)
    zzab=zza+zb
    if (negab(zzab,ab)) cycle
    do zc=1,nzet(sh3)
      zzc=maxprim*(zc-1)
      do zd=1,nzet(sh4)
        zzcd=zzc+zd
        if (negab(zzcd,cd)) cycle
        ooppq=1._dp/(pab(zzab,ab)+pab(zzcd,cd))
        poppq=pab(zzab,ab)*ooppq
        qpq=pab(zzcd,cd)*poppq
        rpq2=sum((Px(:,zzab,ab)-Px(:,zzcd,cd))**2)
        Wx=(pab(zzab,ab)*Px(:,zzab,ab)+pab(zzcd,cd)*Px(:,zzcd,cd))*ooppq
        rdq=Px(:,zzcd,cd)-xsh(:,sh4)
        rqw=Wx-Px(:,zzcd,cd)
        ! Build (ss|ss)^m, m=0,1,2,3,4,5,6
        !call repulsion_f_vector(5,qpq*rpq2,bfac)
        call boys_n(5,qpq*rpq2,bfac)
        osv(1,:)=kab(zzab,ab)*kab(zzcd,cd)*bfac*sqrt(ooppq)
        ! OS recursion: Build (ss|sp)^m, m=0-5
        do m=0,4
          osv(2:4,m)=rdq(:)*osv(1,m)+rqw(:)*osv(1,m+1)
        end do
        ! OS recursion: Build (ss|sd)^m, m=0-4
        afac(0:3)=(osv(1,0:3)-poppq*osv(1,1:4))*hpm1(zzcd,cd)
        do lp1=5,10
          do jj=1,3
            l=minus1(jj,lp1)
            if (l.eq.0) cycle
            jrec=jj
            exit
          end do
          osv(lp1,0:3)=osv(l,0:3)*rdq(jrec)+osv(l,1:4)*rqw(jrec)
          lm1=minus1(jrec,l)
          if (lm1.eq.0) cycle
          osv(lp1,0:3)=osv(lp1,0:3)+afac(0:3)
        end do
        ! With my new notation, I can write a loop all the way up to (ss|si)
        lf=10
        do bigL=3,5
          li=lf ; lf=li+(bigL+1)*(bigL+2)/2  ! These go 10-20, 20-35, 35-56 and 56-84
          m=5-bigL
          do lp1=li+1,lf
            do jj=1,3
              l=minus1(jj,lp1)
              if(l.eq.0) cycle
              jrec=jj
              exit
            end do
            osv(lp1,0:m)=osv(l,0:m)*rdq(jrec)+osv(l,1:m+1)*rqw(jrec)
            lm1=minus1(jrec,l)
            if (lm1.eq.0) cycle
            osv(lp1,0:m)=osv(lp1,0:m)+real(indexps(jrec,l),dp)*hpm1(zzcd,cd)* &
                   & (osv(lm1,0:m)-poppq*osv(lm1,1:m+1))
          end do
        end do
        qop=pab(zzcd,cd)*pm1(zzab,ab)
        delv=pm1(zzab,ab)*(zeta(za,sh1)*rab(:)+zeta(zc,sh3)*rcd(:))
        oset(1,:)=osv(:,0)
        do l=1,35
          do jj=1,3
            lm1=minus1(jj,l)
            lp1=plus1(jj,l)
            jp1=jj+1
            oset(jp1,l)=-delv(jj)*oset(1,l)-qop*oset(1,lp1)
            if (lm1.ne.0) then
              oset(jp1,l)=oset(jp1,l)+hpm1(zzab,ab)*indexps(jj,l)*oset(1,lm1)
            end if
          end do
        end do
        ! Then we assemble [sd|sd], [sd|sf] and [sd|sg] (g means up to 35)
        do l=1,20
          do jp1=5,10
            do jj=1,3
              j=minus1s(jj,jp1) ! This will range from 2 to 4
              if (j.eq.0) cycle
              jrec=jj
              exit
            end do
            lp1=plus1(jrec,l) ! This will from 11 to 56
            jm1=minus1s(jrec,j) ! This will be either 0 or 1
            lm1=minus1(jrec,l) ! And this will go from 2 to 20
            oset(jp1,l)=-delv(jrec)*oset(j,l)-qop*oset(j,lp1)
            if (jm1.ne.0) then
              oset(jp1,l)=oset(jp1,l)+indexps(jrec,j)*hpm1(zzab,ab)*oset(jm1,l)
            end if
            if (lm1.ne.0) then
              oset(jp1,l)=oset(jp1,l)+indexps(jrec,l)*hpm1(zzab,ab)*oset(j,lm1)
            end if
          end do
        end do
        osc=osc+oset(5:10,1:20)
      end do
    end do
  end do
end do

hgp(:,1,1:20) = osc(:,1:20)
! Now build the (sd|pf) and (sd|pd)
do l=5,10
  do jj=1,3
    kp1=jj+1
    lp1=plus1(jj,l)
    hgp(:,kp1,l)=hgp(:,1,lp1)+rcd(jj)*hgp(:,1,l)
  end do
end do

! Finally do contraction to form spherical harmonics, using temporary arrays tmp(6,6,6) and tmp2
! This is not very efficient but very simple
!First take only the Cartesians:
sdpdc=hgp(:,2:4,5:10)

do l=1,5
  do k=1,3
    do j=1,6
      tmpl(j,k,l)=sum(matharmond(l,:)*sdpdc(j,k,:))
    end do
  end do
end do
do l=1,5
  do k=1,3
    do j=1,5
      esdpd(j,k,l)=sum(matharmond(j,:)*tmpl(:,k,l))
    end do
  end do
end do

return
end subroutine i2e_sdpd


pure subroutine i2e_pppd(sh1,sh2,sh3,sh4,epppd)
use nrtype ; use molprops ; use boys
implicit none

integer(i4b), intent(in) :: sh1,sh2,sh3,sh4
real(dp), intent(out) :: epppd(3,3,3,5)

integer(i4b) :: j,l,jm1,lm1,li,lf, bigL, i, ip1, jj, jp1, k, kp1 , lp1, jrec
integer(i4b) :: ab, cd, za, zb, zc, zd, zza, zzc, zzab, zzcd, m
real(dp) :: ooppq, qpq, qop, delv(3), poppq, rpq2, fac1
real(dp) :: rab(3), rcd(3), wx(3)
real(dp) :: bfac(0:5), osv(56,0:5),oset(10,56), rdq(3), rqw(3)
real(dp) :: osc(10,20), hgpk(10,4,20), hgpb(4,10,3,6), pppdc(3,3,3,6)
real(dp) :: afac(0:3)

! This is my first try at doing the Giese thing - first build [ss|sk] (k has L = 8)
!   then do electron transfer, then contract, then use the HGP HR

osc=0._dp

ab=ind2(sh1,sh2)
cd=ind2(sh3,sh4)
rab=xsh(:,sh2)-xsh(:,sh1)
rcd=xsh(:,sh4)-xsh(:,sh3)
do za=1,nzet(sh1)
  zza=maxprim*(za-1)
  do zb=1,nzet(sh2)
    zzab=zza+zb
    if (negab(zzab,ab)) cycle
    do zc=1,nzet(sh3)
      zzc=maxprim*(zc-1)
      do zd=1,nzet(sh4)
        zzcd=zzc+zd
        if (negab(zzcd,cd)) cycle
        ooppq=1._dp/(pab(zzab,ab)+pab(zzcd,cd))
        poppq=pab(zzab,ab)*ooppq
        qpq=pab(zzcd,cd)*poppq
        rpq2=sum((Px(:,zzab,ab)-Px(:,zzcd,cd))**2)
        Wx=(pab(zzab,ab)*Px(:,zzab,ab)+pab(zzcd,cd)*Px(:,zzcd,cd))*ooppq
        rdq=Px(:,zzcd,cd)-xsh(:,sh4)
        rqw=Wx-Px(:,zzcd,cd)
        ! Build (ss|ss)^m, m=0,1,2,3,4,5,6,7,8
        call boys_n(5,qpq*rpq2,bfac)
        osv(1,:)=kab(zzab,ab)*kab(zzcd,cd)*bfac*sqrt(ooppq)
        ! OS recursion: Build (ss|sp)^m, m=0-7
        do m=0,4
          osv(2:4,m)=rdq(:)*osv(1,m)+rqw(:)*osv(1,m+1)
        end do
        ! OS recursion: Build (ss|sd)^m, m=0-6
        afac(0:3)=(osv(1,0:3)-poppq*osv(1,1:4))*hpm1(zzcd,cd)
        do lp1=5,10
          do jj=1,3
            l=minus1(jj,lp1)
            if (l.eq.0) cycle
            jrec=jj
            exit
          end do
          osv(lp1,0:3)=osv(l,0:3)*rdq(jrec)+osv(l,1:4)*rqw(jrec)
          lm1=minus1(jrec,l)
          if (lm1.eq.0) cycle
          osv(lp1,0:3)=osv(lp1,0:3)+afac(0:3)
        end do
        ! With my new notation, I can write a loop all the way up to (ss|sk)
        ! OS recursion: Build (ss|sf)
        lf=10
        do bigL=3,5
          li=lf ; lf=li+(bigL+1)*(bigL+2)/2
          m=5-bigL
          do lp1=li+1,lf
            do jj=1,3
              l=minus1(jj,lp1)
              if(l.eq.0) cycle
              jrec=jj
              exit
            end do
            osv(lp1,0:m)=osv(l,0:m)*rdq(jrec)+osv(l,1:m+1)*rqw(jrec)
            lm1=minus1(jrec,l)
            if (lm1.eq.0) cycle
            osv(lp1,0:m)=osv(lp1,0:m)+real(indexps(jrec,l),dp)*hpm1(zzcd,cd)* &
                   & (osv(lm1,0:m)-poppq*osv(lm1,1:m+1))
          end do
        end do
        ! So now we have [ss|ss], [ss|sp] ... [ss|sk]
        ! We no longer need to track m, but do need to handle the angular momentum on b, 'k'.
        ! Start by moving to a different array, and computing some needed quantities
        qop=pab(zzcd,cd)*pm1(zzab,ab)
        delv=-pm1(zzab,ab)*(zeta(za,sh1)*rab(:)+zeta(zc,sh3)*rcd(:))
        oset(1,:)=osv(:,0)
        ! Now build up the bra terms. For the HR, I need (sd|sd), (sd|sf), (sd|sg), (sf|sd), (sf|sf), (sf|sg), (sg|sd), (sg|sf), (sg|sg) 
        ! Equation (45) in Giese's HSERILib document becomes here:
        ! [sg|sg] = del_i* [sf|sg] -qop* [sf|sh] + n_i(g)*oo2p* [sf|sf] + n_i(p)*oo2p* [sd|sg]
        !   basically reguiring [sf|sf], [sf|sg], [sf|sh]
        ! Following eq. (45), Giese says we need [sp|sp] to [sp|sj], [sd|sd] to [sd|si], [sf|sf] to [sf|sh]
        ! This assembles [sp|ss] to [sp|sj]
        do l=1,35
          do jj=1,3
            lm1=minus1(jj,l)
            lp1=plus1(jj,l)
            jp1=jj+1
            oset(jp1,l)=delv(jj)*oset(1,l)-qop*oset(1,lp1)
            if (lm1.ne.0) then
              oset(jp1,l)=oset(jp1,l)+hpm1(zzab,ab)*indexps(jj,l)*oset(1,lm1)
            end if
          end do
        end do
        ! This assembles [sd|sp] to [sd|si]. It requires [sp|ss]
        do l=1,20
          do jp1=5,10
            do jj=1,3
              j=minus1s(jj,jp1)
              if (j.eq.0) cycle
              jrec=jj
              exit
            end do
            lp1=plus1(jrec,l)
            jm1=minus1(jrec,j)
            lm1=minus1(jrec,l)
            oset(jp1,l)=delv(jrec)*oset(j,l)-qop*oset(j,lp1)
            if (jm1.ne.0) then
              oset(jp1,l)=oset(jp1,l)+indexps(jrec,j)*hpm1(zzab,ab)*oset(jm1,l)
            end if
            if (lm1.ne.0) then
              oset(jp1,l)=oset(jp1,l)+indexps(jrec,l)*hpm1(zzab,ab)*oset(j,lm1)
            end if
          end do
        end do
        osc=osc+oset(:,1:20)
      end do
    end do
  end do
end do

! For each of the bra and ket terms, the HRR becomes: (here for ket):
! |dd) = |pf) + (Di-Ci) * |pd)
! |pf) = |sg) + (Di-Ci) * |sf)
! |pd) = |sf) + (Di-Ci) * |sd)
! We need as input (sd|sd)*, (sd|sf)*, (sd|sg)*, (sf|sd), (sf|sf), (sf|sg), (sg|sd), (sg|sf), (sg|sg)

hgpk(:,1,:) = osc
do l=5,10
  do jj=1,3
    kp1=1+jj
    lp1=plus1(jj,l)
    hgpk(:,kp1,l)=hgpk(:,1,lp1)+rcd(jj)*hgpk(:,1,l)
  end do
end do
hgpb(1,:,:,:)=hgpk(:,2:4,5:10)
do j=2,4
  do jj=1,3
    ip1=1+jj
    jp1=plus1(jj,j)
    hgpb(ip1,j,:,:)=hgpb(1,jp1,:,:)+rab(jj)*hgpb(1,j,:,:) 
  end do
end do

! Finally do contraction to form spherical harmonics, using temporary arrays tmp(6,6,6,6) and tmp2
! This is not very efficient but very simple
!First take only the Cartesians:
pppdc=hgpb(2:4,2:4,:,:)

! Then recontract
do l=1,5
  do k=1,3
    do j=1,3
      do i=1,3
        epppd(i,j,k,l)=sum(matharmond(l,:)*pppdc(i,j,k,:))
      end do
    end do
  end do
end do

return
end subroutine i2e_pppd




pure subroutine i2e_ppdd(sh1,sh2,sh3,sh4,eppdd)
use nrtype ; use molprops ; use boys
implicit none

integer(i4b), intent(in) :: sh1,sh2,sh3,sh4
real(dp), intent(out) :: eppdd(3,3,5,5)

integer(i4b) :: j,l,jm1,lm1,li,lf, bigL, i, ip1, jj, jp1, k, kp1 , lp1
integer(i4b) :: ab, cd, za, zb, zc, zd, zza, zzc, zzab, zzcd, m, jrec
real(dp) :: ooppq, qpq, qop, delv(3), poppq, rpq2, fac1
real(dp) :: rab(3), rcd(3), wx(3)
real(dp) :: bfac(0:6), osv(84,0:6),oset(10,84), rdq(3), rqw(3)
real(dp) :: osc(10,35), hgpk(10,10,35), hgpb(4,10,6,6), ppddc(3,3,6,6)
real(dp) :: afac(0:4), tmpl(3,3,6,5)

osc=0._dp

ab=ind2(sh1,sh2)
cd=ind2(sh3,sh4)
rab=xsh(:,sh2)-xsh(:,sh1)
rcd=xsh(:,sh4)-xsh(:,sh3)
do za=1,nzet(sh1)
  zza=maxprim*(za-1)
  do zb=1,nzet(sh2)
    zzab=zza+zb
    if (negab(zzab,ab)) cycle
    do zc=1,nzet(sh3)
      zzc=maxprim*(zc-1)
      do zd=1,nzet(sh4)
        zzcd=zzc+zd
        if (negab(zzcd,cd)) cycle
        ooppq=1._dp/(pab(zzab,ab)+pab(zzcd,cd))
        poppq=pab(zzab,ab)*ooppq
        qpq=pab(zzcd,cd)*poppq
        rpq2=sum((Px(:,zzab,ab)-Px(:,zzcd,cd))**2)
        Wx=(pab(zzab,ab)*Px(:,zzab,ab)+pab(zzcd,cd)*Px(:,zzcd,cd))*ooppq
        rdq=Px(:,zzcd,cd)-xsh(:,sh4)
        rqw=Wx-Px(:,zzcd,cd)
        ! Build (ss|ss)^m, m=0,1,2,3,4,5,6,7,8
        call boys_n(6,qpq*rpq2,bfac)
        osv(1,:)=kab(zzab,ab)*kab(zzcd,cd)*bfac*sqrt(ooppq)
        ! OS recursion: Build (ss|sp)^m, m=0-7
        do m=0,5
          osv(2:4,m)=rdq(:)*osv(1,m)+rqw(:)*osv(1,m+1)
        end do
        ! OS recursion: Build (ss|sd)^m, m=0-6
        afac(0:4)=(osv(1,0:4)-poppq*osv(1,1:5))*hpm1(zzcd,cd)
        do lp1=5,10
          do jj=1,3
            l=minus1(jj,lp1)
            if (l.eq.0) cycle
            jrec=jj
            exit
          end do
          osv(lp1,0:4)=osv(l,0:4)*rdq(jrec)+osv(l,1:5)*rqw(jrec)
          lm1=minus1(jrec,l)
          if (lm1.eq.0) cycle
          osv(lp1,0:4)=osv(lp1,0:4)+afac(0:4)
        end do
        ! With my new notation, I can write a loop all the way up to (ss|sk)
        ! OS recursion: Build (ss|sf)
        lf=10
        do bigL=3,6
          li=lf ; lf=li+(bigL+1)*(bigL+2)/2
          m=6-bigL
          do lp1=li+1,lf
            do jj=1,3
              l=minus1(jj,lp1)
              if(l.eq.0) cycle
              jrec=jj
              exit
            end do
            osv(lp1,0:m)=osv(l,0:m)*rdq(jrec)+osv(l,1:m+1)*rqw(jrec)
            lm1=minus1(jrec,l)
            if (lm1.eq.0) cycle
            osv(lp1,0:m)=osv(lp1,0:m)+real(indexps(jrec,l),dp)*hpm1(zzcd,cd)* &
                   & (osv(lm1,0:m)-poppq*osv(lm1,1:m+1))
          end do
        end do
        ! So now we have [ss|ss], [ss|sp] ... [ss|sk]
        ! We no longer need to track m, but do need to handle the angular momentum on b, 'k'.
        ! Start by moving to a different array, and computing some needed quantities
        qop=pab(zzcd,cd)*pm1(zzab,ab)
        delv=-pm1(zzab,ab)*(zeta(za,sh1)*rab(:)+zeta(zc,sh3)*rcd(:))
        oset(1,:)=osv(:,0)
        ! Now build up the bra terms. For the HR, I need (sd|sd), (sd|sf), (sd|sg), (sf|sd), (sf|sf), (sf|sg), (sg|sd), (sg|sf), (sg|sg) 
        ! Equation (45) in Giese's HSERILib document becomes here:
        ! [sg|sg] = del_i* [sf|sg] -qop* [sf|sh] + n_i(g)*oo2p* [sf|sf] + n_i(p)*oo2p* [sd|sg]
        !   basically reguiring [sf|sf], [sf|sg], [sf|sh]
        ! Following eq. (45), Giese says we need [sp|sp] to [sp|sj], [sd|sd] to [sd|si], [sf|sf] to [sf|sh]
        ! This assembles [sp|ss] to [sp|sj]
        do l=1,56
          do jj=1,3
            lm1=minus1(jj,l)
            lp1=plus1(jj,l)
            jp1=jj+1
            oset(jp1,l)=delv(jj)*oset(1,l)-qop*oset(1,lp1)
            if (lm1.ne.0) then
              oset(jp1,l)=oset(jp1,l)+hpm1(zzab,ab)*indexps(jj,l)*oset(1,lm1)
            end if
          end do
        end do
        ! This assembles [sd|sp] to [sd|si]. It requires [sp|ss]
        do l=1,35
          do jp1=5,10
            do jj=1,3
              j=minus1s(jj,jp1)
              if (j.eq.0) cycle
              jrec=jj
              exit
            end do
            lp1=plus1(jrec,l)
            jm1=minus1(jrec,j)
            lm1=minus1(jrec,l)
            oset(jp1,l)=delv(jrec)*oset(j,l)-qop*oset(j,lp1)
            if (jm1.ne.0) then
              oset(jp1,l)=oset(jp1,l)+indexps(jrec,j)*hpm1(zzab,ab)*oset(jm1,l)
            end if
            if (lm1.ne.0) then
              oset(jp1,l)=oset(jp1,l)+indexps(jrec,l)*hpm1(zzab,ab)*oset(j,lm1)
            end if
          end do
        end do
        osc=osc+oset(:,1:35)
      end do
    end do
  end do
end do

! So now we have (sd|sg) to (sp|sj), (sd|sd) to (sd|si), (sf|sf) to (sf|sh) and (sg|sg)

! For each of the bra and ket terms, the HRR becomes: (here for ket):
! |dd) = |pf) + (Di-Ci) * |pd)
! |pf) = |sg) + (Di-Ci) * |sf)
! |pd) = |sf) + (Di-Ci) * |sd)
! We need as input (sd|sd)*, (sd|sf)*, (sd|sg)*, (sf|sd), (sf|sf), (sf|sg), (sg|sd), (sg|sf), (sg|sg)

hgpk(:,1,:) = osc(:,:)
! Now build the (sx|pf) and (sx|pd)
do l=5,20
  do jj=1,3
    kp1=1+jj
    lp1=plus1(jj,l)
    hgpk(:,kp1,l)=hgpk(:,1,lp1)+rcd(jj)*hgpk(:,1,l)
  end do
end do
! Now build the (sx|dd) 
do kp1=5,10
  do l=5,10
    do jj=1,3
      k=minus1(jj,kp1)
      if (k.eq.0) cycle
      jrec=jj
      exit
    end do
    lp1=plus1(jrec,l)
    hgpk(:,kp1,l)=hgpk(:,k,lp1)+rcd(jrec)*hgpk(:,k,l)
  end do
end do
! Now move the resulting (sd|dd), (sf|dd) and (sg|dd) ints
hgpb(1,:,:,:)=hgpk(:,5:10,5:10)
! Now build the (pf|dd) and (pd|dd)
do j=2,4
  do jj=1,3
    ip1=1+jj
    jp1=plus1(jj,j)
    hgpb(ip1,j,:,:)=hgpb(1,jp1,:,:)+rab(jj)*hgpb(1,j,:,:) 
  end do
end do

! Finally do contraction to form spherical harmonics, using temporary arrays tmp(6,6,6,6) and tmp2
! This is not very efficient but very simple
!First take only the Cartesians:
ppddc=hgpb(2:4,2:4,:,:)

! Then recontract
do l=1,5
  do k=1,6
    do j=1,3
      do i=1,3
        tmpl(i,j,k,l)=sum(matharmond(l,:)*ppddc(i,j,k,:))
      end do
    end do
  end do
end do
do l=1,5
  do k=1,5
    do j=1,3
      do i=1,3
        eppdd(i,j,k,l)=sum(matharmond(k,:)*tmpl(i,j,:,l))
      end do
    end do
  end do
end do

return
end subroutine i2e_ppdd



pure subroutine i2e_pdpd(sh1,sh2,sh3,sh4,epdpd)
use nrtype ; use molprops ; use boys
implicit none

integer(i4b), intent(in) :: sh1,sh2,sh3,sh4
real(dp), intent(out) :: epdpd(3,5,3,5)

integer(i4b) :: j,l,jm1,lm1,li,lf, bigL, i, ip1, jj, jp1, k, kp1 , lp1
integer(i4b) :: ab, cd, za, zb, zc, zd, zza, zzc, zzab, zzcd, m, jrec
real(dp) :: ooppq, qpq, qop, delv(3), poppq, rpq2, fac1
real(dp) :: rab(3), rcd(3), wx(3)
real(dp) :: bfac(0:6), osv(84,0:6),oset(20,84), rdq(3), rqw(3)
real(dp) :: osc(20,20), hgpk(20,4,20), hgpb(4,20,3,6), pdpdc(3,6,3,6)
real(dp) :: afac(0:4), tmpl(3,6,3,5)

! This is my first try at doing the Giese thing - first build [ss|sk] (k has L = 8)
!   then do electron transfer, then contract, then use the HGP HR

osc=0._dp

ab=ind2(sh1,sh2)
cd=ind2(sh3,sh4)
rab=xsh(:,sh2)-xsh(:,sh1)
rcd=xsh(:,sh4)-xsh(:,sh3)
do za=1,nzet(sh1)
  zza=maxprim*(za-1)
  do zb=1,nzet(sh2)
    zzab=zza+zb
    if (negab(zzab,ab)) cycle
    do zc=1,nzet(sh3)
      zzc=maxprim*(zc-1)
      do zd=1,nzet(sh4)
        zzcd=zzc+zd
        if (negab(zzcd,cd)) cycle
        ooppq=1._dp/(pab(zzab,ab)+pab(zzcd,cd))
        poppq=pab(zzab,ab)*ooppq
        qpq=pab(zzcd,cd)*poppq
        rpq2=sum((Px(:,zzab,ab)-Px(:,zzcd,cd))**2)
        Wx=(pab(zzab,ab)*Px(:,zzab,ab)+pab(zzcd,cd)*Px(:,zzcd,cd))*ooppq
        rdq=Px(:,zzcd,cd)-xsh(:,sh4)
        rqw=Wx-Px(:,zzcd,cd)
        ! Build (ss|ss)^m, m=0,1,2,3,4,5,6,7,8
        call boys_n(6,qpq*rpq2,bfac)
        osv(1,:)=kab(zzab,ab)*kab(zzcd,cd)*bfac*sqrt(ooppq)
        ! OS recursion: Build (ss|sp)^m, m=0-7
        do m=0,5
          osv(2:4,m)=rdq(:)*osv(1,m)+rqw(:)*osv(1,m+1)
        end do
        ! OS recursion: Build (ss|sd)^m, m=0-6
        afac(0:4)=(osv(1,0:4)-poppq*osv(1,1:5))*hpm1(zzcd,cd)
        do lp1=5,10
          do jj=1,3
            l=minus1(jj,lp1)
            if (l.eq.0) cycle
            jrec=jj 
            exit
          end do
          osv(lp1,0:4)=osv(l,0:4)*rdq(jrec)+osv(l,1:5)*rqw(jrec)
          lm1=minus1(jrec,l)
          if (lm1.eq.0) cycle
          osv(lp1,0:4)=osv(lp1,0:4)+afac(0:4)
        end do
        ! With my new notation, I can write a loop all the way up to (ss|sk)
        ! OS recursion: Build (ss|sf)
        lf=10
        do bigL=3,6
          li=lf ; lf=li+(bigL+1)*(bigL+2)/2
          m=6-bigL
          do lp1=li+1,lf
            do jj=1,3
              l=minus1(jj,lp1)
              if(l.eq.0) cycle
              jrec=jj
              exit
            end do
            osv(lp1,0:m)=osv(l,0:m)*rdq(jrec)+osv(l,1:m+1)*rqw(jrec)
            lm1=minus1(jrec,l)
            if (lm1.eq.0) cycle
            osv(lp1,0:m)=osv(lp1,0:m)+real(indexps(jrec,l),dp)*hpm1(zzcd,cd)* &
                   & (osv(lm1,0:m)-poppq*osv(lm1,1:m+1))
          end do
        end do
        qop=pab(zzcd,cd)*pm1(zzab,ab)
        delv=-pm1(zzab,ab)*(zeta(za,sh1)*rab(:)+zeta(zc,sh3)*rcd(:))
        oset(1,:)=osv(:,0)
        ! Now build up the bra terms. For the HR, I need (sd|sd), (sd|sf), (sf|sd), (sf|sf)
        do l=1,56 ! First make [sp|ss] to [sp|sh], using [ss|ss] to [ss|si]
          do jj=1,3
            lm1=minus1(jj,l)
            lp1=plus1(jj,l)
            jp1=jj+1
            oset(jp1,l)=delv(jj)*oset(1,l)-qop*oset(1,lp1)
            if (lm1.ne.0) then
              oset(jp1,l)=oset(jp1,l)+hpm1(zzab,ab)*indexps(jj,l)*oset(1,lm1)
            end if
          end do
        end do
        ! This assembles [sd|ss] to [sd|sg]. It requires [sp|ss] - so this makes two of the targets, [sd|sd] and [sd|sf]
        do l=1,35
          do jp1=5,10
            do jj=1,3
              j=minus1(jj,jp1)
              if (j.eq.0) cycle
              jrec=jj
              exit
            end do
            lp1=plus1(jrec,l)
            jm1=minus1(jrec,j)
            lm1=minus1(jrec,l)
            oset(jp1,l)=delv(jrec)*oset(j,l)-qop*oset(j,lp1) ! This uses from [sp|ss] to [sp|sh]
            if (jm1.ne.0) then
              oset(jp1,l)=oset(jp1,l)+indexps(jrec,j)*hpm1(zzab,ab)*oset(jm1,l) ! And this uses [ss|ss] to [ss|sg]
            end if
            if (lm1.ne.0) then
              oset(jp1,l)=oset(jp1,l)+indexps(jrec,l)*hpm1(zzab,ab)*oset(j,lm1) ! And this uses [sp|ss] to [ss|sf]
            end if
          end do
        end do
        ! This assembles [sf|sd] to [sf|sf]
        do l=5,20
          do jp1=11,20
            do jj=1,3
              j=minus1(jj,jp1)
              if (j.eq.0) cycle
              jrec=jj
              exit
            end do
            lp1=plus1(jrec,l)
            jm1=minus1s(jrec,j)
            lm1=minus1(jrec,l)
            oset(jp1,l)=delv(jrec)*oset(j,l)-qop*oset(j,lp1)
            if (jm1.ne.0) then
              oset(jp1,l)=oset(jp1,l)+indexps(jrec,j)*hpm1(zzab,ab)*oset(jm1,l)
            end if
            if (lm1.ne.0) then
              oset(jp1,l)=oset(jp1,l)+indexps(jrec,l)*hpm1(zzab,ab)*oset(j,lm1)
            end if
          end do
        end do
        osc=osc+oset(:,1:20)
      end do
    end do
  end do
end do

! We need as input for HRR:
! (sf|sf), (sd|sf), (sf|sd), (sd|sd)

hgpk(:,1,:) = osc
! We have (xx|pd) = (xx|sf) + rcd*(xx|sd) and same for bra
do l=5,10
  do jj=1,3
    kp1=1+jj
    lp1=plus1(jj,l)
    hgpk(:,kp1,l)=hgpk(:,1,lp1)+rcd(jj)*hgpk(:,1,l)
  end do
end do
hgpb(1,:,:,:)=hgpk(:,2:4,5:10)
! We have (xx|pd) = (xx|sf) + rcd*(xx|sd) and same for bra
do j=5,10
  do jj=1,3
    ip1=1+jj
    jp1=plus1(jj,j)
    hgpb(ip1,j,:,:)=hgpb(1,jp1,:,:)+rab(jj)*hgpb(1,j,:,:) 
  end do
end do

! Finally do contraction to form spherical harmonics, using temporary arrays tmp(6,6,6,6) and tmp2
! This is not very efficient but very simple
!First take only the Cartesians:
pdpdc=hgpb(2:4,5:10,:,:)

! Then recontract
do l=1,5
  do k=1,3
    do j=1,6
      do i=1,3
        tmpl(i,j,k,l)=sum(matharmond(l,:)*pdpdc(i,j,k,:))
      end do
    end do
  end do
end do
do l=1,5
  do k=1,3
    do j=1,5
      do i=1,3
        epdpd(i,j,k,l)=sum(matharmond(j,:)*tmpl(i,:,k,l))
      end do
    end do
  end do
end do

return
end subroutine i2e_pdpd


pure subroutine i2e_pddd(sh1,sh2,sh3,sh4,epddd)
use nrtype ; use molprops ; use boys
implicit none

integer(i4b), intent(in) :: sh1,sh2,sh3,sh4
real(dp), intent(out) :: epddd(3,5,5,5)

integer(i4b) :: j,l,jm1,lm1,li,lf, bigL, i, ip1, jj, jp1, k, kp1 , lp1
integer(i4b) :: ab, cd, za, zb, zc, zd, zza, zzc, zzab, zzcd, m, jrec
real(dp) :: ooppq, qpq, qop, delv(3), poppq, rpq2, fac1
real(dp) :: rab(3), rcd(3), wx(3)
real(dp) :: bfac(0:7), osv(120,0:7),oset(20,120), rdq(3), rqw(3)
real(dp) :: osc(20,35), hgpk(20,10,35), hgpb(4,20,6,6), pdddc(3,6,6,6)
real(dp) :: afac(0:5), tmpl(3,6,6,5), tmpkl(3,6,5,5)

! This is my first try at doing the Giese thing - first build [ss|sk] (k has L = 8)
!   then do electron transfer, then contract, then use the HGP HR

osc=0._dp

ab=ind2(sh1,sh2)
cd=ind2(sh3,sh4)
rab=xsh(:,sh2)-xsh(:,sh1)
rcd=xsh(:,sh4)-xsh(:,sh3)
do za=1,nzet(sh1)
  zza=maxprim*(za-1)
  do zb=1,nzet(sh2)
    zzab=zza+zb
    if (negab(zzab,ab)) cycle
    do zc=1,nzet(sh3)
      zzc=maxprim*(zc-1)
      do zd=1,nzet(sh4)
        zzcd=zzc+zd
        if (negab(zzcd,cd)) cycle
        ooppq=1._dp/(pab(zzab,ab)+pab(zzcd,cd))
        poppq=pab(zzab,ab)*ooppq
        qpq=pab(zzcd,cd)*poppq
        rpq2=sum((Px(:,zzab,ab)-Px(:,zzcd,cd))**2)
        Wx=(pab(zzab,ab)*Px(:,zzab,ab)+pab(zzcd,cd)*Px(:,zzcd,cd))*ooppq
        rdq=Px(:,zzcd,cd)-xsh(:,sh4)
        rqw=Wx-Px(:,zzcd,cd)
        ! Build (ss|ss)^m, m=0,1,2,3,4,5,6,7,8
        call boys_n(7,qpq*rpq2,bfac)
        osv(1,:)=kab(zzab,ab)*kab(zzcd,cd)*bfac*sqrt(ooppq)
        ! OS recursion: Build (ss|sp)^m, m=0-7
        do m=0,6
          osv(2:4,m)=rdq(:)*osv(1,m)+rqw(:)*osv(1,m+1)
        end do
        ! OS recursion: Build (ss|sd)^m, m=0-6
        afac(0:5)=(osv(1,0:5)-poppq*osv(1,1:6))*hpm1(zzcd,cd)
        do lp1=5,10
          do jj=1,3
            l=minus1(jj,lp1)
            if (l.eq.0) cycle
            jrec=jj
            exit
          end do
          osv(lp1,0:5)=osv(l,0:5)*rdq(jrec)+osv(l,1:6)*rqw(jrec)
          lm1=minus1(jrec,l)
          if (lm1.eq.0) cycle
          osv(lp1,0:5)=osv(lp1,0:5)+afac(0:5)
        end do
        ! With my new notation, I can write a loop all the way up to (ss|sk)
        ! OS recursion: Build (ss|sf)
        lf=10
        do bigL=3,7
          li=lf ; lf=li+(bigL+1)*(bigL+2)/2
          m=7-bigL
          do lp1=li+1,lf
            do jj=1,3
              l=minus1(jj,lp1)
              if(l.eq.0) cycle
              jrec=jj
              exit
            end do
            osv(lp1,0:m)=osv(l,0:m)*rdq(jrec)+osv(l,1:m+1)*rqw(jrec)
            lm1=minus1(jrec,l)
            if (lm1.eq.0) cycle
            osv(lp1,0:m)=osv(lp1,0:m)+real(indexps(jrec,l),dp)*hpm1(zzcd,cd)* &
                   & (osv(lm1,0:m)-poppq*osv(lm1,1:m+1))
          end do
        end do
        qop=pab(zzcd,cd)*pm1(zzab,ab)
        delv=-pm1(zzab,ab)*(zeta(za,sh1)*rab(:)+zeta(zc,sh3)*rcd(:))
        oset(1,:)=osv(:,0)
        do l=1,84
          do jj=1,3
            lm1=minus1(jj,l)
            lp1=plus1(jj,l)
            jp1=jj+1
            oset(jp1,l)=delv(jj)*oset(1,l)-qop*oset(1,lp1)
            if (lm1.ne.0) then
              oset(jp1,l)=oset(jp1,l)+hpm1(zzab,ab)*indexps(jj,l)*oset(1,lm1)
            end if
          end do
        end do
        ! This assembles [sd|sp] to [sd|si]. It requires [sp|ss]
        do l=1,56
          do jp1=5,10
            do jj=1,3
              j=minus1s(jj,jp1)
              if (j.eq.0) cycle
              jrec=jj
              exit
            end do
            lp1=plus1(jrec,l)
            jm1=minus1(jrec,j)
            lm1=minus1(jrec,l)
            oset(jp1,l)=delv(jrec)*oset(j,l)-qop*oset(j,lp1)
            if (jm1.ne.0) then
              oset(jp1,l)=oset(jp1,l)+indexps(jrec,j)*hpm1(zzab,ab)*oset(jm1,l)
            end if
            if (lm1.ne.0) then
              oset(jp1,l)=oset(jp1,l)+indexps(jrec,l)*hpm1(zzab,ab)*oset(j,lm1)
            end if
          end do
        end do
        ! This assembles [sf|sd] to [sf|sh]. This requires among others [sd|sp]
        do l=1,35
          do jp1=11,20
            do jj=1,3
              j=minus1s(jj,jp1)
              if (j.eq.0) cycle
              jrec=jj
              exit
            end do
            lp1=plus1(jrec,l)
            jm1=minus1s(jrec,j)
            lm1=minus1(jrec,l)
            oset(jp1,l)=delv(jrec)*oset(j,l)-qop*oset(j,lp1)
            if (jm1.ne.0) then
              oset(jp1,l)=oset(jp1,l)+indexps(jrec,j)*hpm1(zzab,ab)*oset(jm1,l)
            end if
            if (lm1.ne.0) then
              oset(jp1,l)=oset(jp1,l)+indexps(jrec,l)*hpm1(zzab,ab)*oset(j,lm1)
            end if
          end do
        end do
        osc=osc+oset(:,1:35)
      end do
    end do
  end do
end do

hgpk(:,1,:) = osc
! Now build the (sx|pf) and (sx|pd)
do l=5,20
  do jj=1,3
    kp1=1+jj
    lp1=plus1(jj,l)
    hgpk(:,kp1,l)=hgpk(:,1,lp1)+rcd(jj)*hgpk(:,1,l)
  end do
end do
! Now build the (sx|dd) 
do kp1=5,10
  do l=5,10
    do jj=1,3
      k=minus1(jj,kp1)
      if (k.eq.0) cycle
      jrec=jj 
      exit
    end do
    lp1=plus1(jrec,l)
    hgpk(:,kp1,l)=hgpk(:,k,lp1)+rcd(jrec)*hgpk(:,k,l)
  end do
end do
! Now move the resulting (sd|dd), (sf|dd) and (sg|dd) ints
hgpb(1,:,:,:)=hgpk(:,5:10,5:10)
! Now build the (pf|dd) and (pd|dd)
do j=5,10
  do jj=1,3
    ip1=1+jj
    jp1=plus1(jj,j)
    hgpb(ip1,j,:,:)=hgpb(1,jp1,:,:)+rab(jj)*hgpb(1,j,:,:) 
  end do
end do

! Finally do contraction to form spherical harmonics, using temporary arrays tmp(6,6,6,6) and tmp2
! This is not very efficient but very simple
!First take only the Cartesians:
pdddc=hgpb(2:4,5:10,:,:)

! Then recontract
do l=1,5
  do k=1,6
    do j=1,6
      do i=1,3
        tmpl(i,j,k,l)=sum(matharmond(l,:)*pdddc(i,j,k,:))
      end do
    end do
  end do
end do
do l=1,5
  do k=1,5
    do j=1,6
      do i=1,3
        tmpkl(i,j,k,l)=sum(matharmond(k,:)*tmpl(i,j,:,l))
      end do
    end do
  end do
end do
do l=1,5
  do k=1,5
    do j=1,5
      do i=1,3
        epddd(i,j,k,l)=sum(matharmond(j,:)*tmpkl(i,:,k,l))
      end do
    end do
  end do
end do

return
end subroutine i2e_pddd




end module s_p_d_terms








