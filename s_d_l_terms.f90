module s_d_l_terms
use nrtype
implicit none

contains

pure subroutine i2e_dlll(sh1,sh2,sh3,sh4,edlll)
use nrtype ; use molprops ; use boys
implicit none

integer(i4b), intent(in) :: sh1,sh2,sh3,sh4
real(dp), intent(out) :: edlll(5,4,4,4)

integer(i4b) :: j,l,jm1,lm1,li,lf, bigL, i, ip1, jj, jp1, k, kp1 , lp1
integer(i4b) :: ab, cd, za, zb, zc, zd, zza, zzc, zzab, zzcd, m, jrec
real(dp) :: ooppq, qpq, poq, delv(3),qoppq, rpq2, fac1
real(dp) :: rba(3), rdc(3), wx(3)
real(dp) :: bfac(0:5), osv(56,0:5),oset(56,10), rap(3), rpw(3), osc
real(dp) :: oshb(20,4,10), oshk(6,4,10,4), dlllc(6,4,4,4), oshbc(6,4,10)
real(dp) :: afac(0:3)

! This is my first try at doing the Giese thing - first build [ss|sk] (k has L = 8)
!   then do electron transfer, then contract, then use the HGP HR

oshb=0._dp
oshk=0._dp
oset=0._dp
dlllc=0._dp

ab=ind2(sh1,sh2)
cd=ind2(sh3,sh4)
rba=xsh(:,sh1)-xsh(:,sh2)
rdc=xsh(:,sh3)-xsh(:,sh4)
do zc=1,nzet(sh3)
  zzc=maxprim*(zc-1)
  do zd=1,nzet(sh4)
    zzcd=zzc+zd
    if (negab(zzcd,cd)) cycle
    oshbc=0._dp
    do za=1,nzet(sh1)
      zza=maxprim*(za-1)
      do zb=1,nzet(sh2)
        zzab=zza+zb
        if (negab(zzab,ab)) cycle
        ooppq=1._dp/(pab(zzab,ab)+pab(zzcd,cd))
        qoppq=pab(zzcd,cd)*ooppq
        qpq=pab(zzab,ab)*qoppq
        rpq2=sum((Px(:,zzab,ab)-Px(:,zzcd,cd))**2)
        Wx=(pab(zzab,ab)*Px(:,zzab,ab)+pab(zzcd,cd)*Px(:,zzcd,cd))*ooppq
        rap=Px(:,zzab,ab)-xsh(:,sh1)
        rpw=Wx-Px(:,zzab,ab)
        call boys_n(5,qpq*rpq2,bfac)
        osv(1,:)=kab(zzab,ab)*kab(zzcd,cd)*bfac*sqrt(ooppq)
        ! OS recursion: Build (ps|ss)^m, m=0-4
        do m=0,4
          osv(2:4,m)=rap(:)*osv(1,m)+rpw(:)*osv(1,m+1)
        end do
        ! OS recursion: Build (ds|ss)^m, m=0-3
        afac(0:3)=hpm1(zzab,ab)*(osv(1,0:3)-qoppq*osv(1,1:4))
        do lp1=5,10
          do jj=1,3
            l=minus1(jj,lp1)
            if (l.eq.0) cycle
            jrec=jj
            exit
          end do
          osv(lp1,0:3)=osv(l,0:3)*rap(jrec)+osv(l,1:4)*rpw(jrec)
          lm1=minus1(jj,l)
          if (lm1.eq.0) cycle
          osv(lp1,0:3)=osv(lp1,0:3)+afac(0:3)
        end do
        ! With my new notation, I can write a loop all the way up to (ss|sk)
        ! OS recursion: Build (hs|ss)
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
            osv(lp1,0:m)=osv(l,0:m)*rap(jrec)+osv(l,1:m+1)*rpw(jrec)
            lm1=minus1(jrec,l)
            if (lm1.eq.0) cycle
            osv(lp1,0:m)=osv(lp1,0:m)+real(indexps(jrec,l),dp)*hpm1(zzab,ab)* &
                   & (osv(lm1,0:m)-qoppq*osv(lm1,1:m+1))
          end do
        end do
        poq=pab(zzab,ab)*pm1(zzcd,cd)
        delv=-pm1(zzcd,cd)*(zeta(zb,sh2)*rba(:)+zeta(zd,sh4)*rdc(:))
        oset(:,1)=osv(:,0)
        do l=1,35
          do jj=1,3
            lm1=minus1(jj,l)
            lp1=plus1(jj,l)
            jp1=jj+1
            oset(l,jp1)=delv(jj)*oset(l,1)-poq*oset(lp1,1)
            if (lm1.ne.0) then
              oset(l,jp1)=oset(l,jp1)+hpm1(zzcd,cd)*indexps(jj,l)*oset(lm1,1)
            end if
          end do
        end do
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
            oset(l,jp1)=delv(jrec)*oset(l,j)-poq*oset(lp1,j)
            if (jm1.ne.0) then
              oset(l,jp1)=oset(l,jp1)+indexps(jrec,j)*hpm1(zzcd,cd)*oset(l,jm1)
            end if
            if (lm1.ne.0) then
              oset(l,jp1)=oset(l,jp1)+indexps(jrec,l)*hpm1(zzcd,cd)*oset(lm1,j)
            end if
          end do
        end do
        oshb(:,1,:)=oset(1:20,:)
        do l=5,10
          do jj=1,3
            kp1=1+jj
            lp1=plus1(jj,l)
            oshb(l,kp1,:)=oshb(lp1,1,:)+rba(jj)*oshb(l,1,:)
          end do
        end do
        oshbc(:,1,:)=oshbc(:,1,:)+oshb(5:10,1,:)
        oshbc(:,2:4,:)=oshbc(:,2:4,:)+oshb(5:10,2:4,:)*djkl(zb,sh2)
      end do
    end do
    oshk(:,:,:,1)=oshbc
    do j=1,4
      do jj=1,3
        ip1=1+jj
        jp1=plus1(jj,j)
        oshk(:,:,j,ip1)=oshk(:,:,jp1,1)+rdc(jj)*oshk(:,:,j,1) 
      end do
    end do
    dlllc(:,:,1,1)=dlllc(:,:,1,1)+oshk(:,:,1,1)
    dlllc(:,:,2:4,1)=dlllc(:,:,2:4,1)+oshk(:,:,2:4,1)*djkl(zc,sh3)
    dlllc(:,:,1,2:4)=dlllc(:,:,1,2:4)+oshk(:,:,1,2:4)*djkl(zd,sh4)
    dlllc(:,:,2:4,2:4)=dlllc(:,:,2:4,2:4)+oshk(:,:,2:4,2:4)*djkl(zc,sh3)*djkl(zd,sh4)
  end do
end do

! Then recontract
do l=1,4
  do k=1,4
    do j=1,4
      do i=1,5
        edlll(i,j,k,l)=sum(matharmond(i,:)*dlllc(:,j,k,l))
      end do
    end do
  end do
end do

return
end subroutine i2e_dlll




pure subroutine i2e_ddll(sh1,sh2,sh3,sh4,eddll)
use nrtype ; use molprops ; use boys
implicit none

integer(i4b), intent(in) :: sh1,sh2,sh3,sh4
real(dp), intent(out) :: eddll(5,5,4,4)

integer(i4b) :: j,l,jm1,lm1,li,lf, bigL, i, ip1, jj, jp1, k, kp1 , lp1
integer(i4b) :: ab, cd, za, zb, zc, zd, zza, zzc, zzab, zzcd, m, jrec
real(dp) :: ooppq, qpq, poq, delv(3), qoppq, rpq2, fac1
real(dp) :: rba(3), rdc(3), wx(3)
real(dp) :: bfac(0:6), osv(84,0:6),oset(84,10), rap(3), rpw(3)
real(dp) :: oshb(35,10,10), oshk(6,6,10,4), ddllc(6,6,4,4)
real(dp) :: afac(0:4), tmpj(6,5,4,4), lcd

ddllc=0._dp

ab=ind2(sh1,sh2)
cd=ind2(sh3,sh4)
rba=xsh(:,sh1)-xsh(:,sh2)
rdc=xsh(:,sh3)-xsh(:,sh4)
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
        qoppq=pab(zzcd,cd)*ooppq
        qpq=pab(zzab,ab)*qoppq
        rpq2=sum((Px(:,zzab,ab)-Px(:,zzcd,cd))**2)
        Wx=(pab(zzab,ab)*Px(:,zzab,ab)+pab(zzcd,cd)*Px(:,zzcd,cd))*ooppq
        rap=Px(:,zzab,ab)-xsh(:,sh1)
        rpw=Wx-Px(:,zzab,ab)
        call boys_n(6,qpq*rpq2,bfac)
        osv(1,:)=kab(zzab,ab)*kab(zzcd,cd)*bfac*sqrt(ooppq)
        do m=0,5
          osv(2:4,m)=rap(:)*osv(1,m)+rpw(:)*osv(1,m+1)
        end do
        ! OS recursion: Build (ss|sd)^m, m=0-6
        afac(0:4)=hpm1(zzab,ab)*(osv(1,0:4)-qoppq*osv(1,1:5))
        do lp1=5,10
          do jj=1,3
            l=minus1(jj,lp1)
            if (l.eq.0) cycle
            jrec=jj
            exit
          end do
          osv(lp1,0:4)=osv(l,0:4)*rap(jrec)+osv(l,1:5)*rpw(jrec)
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
            osv(lp1,0:m)=osv(l,0:m)*rap(jrec)+osv(l,1:m+1)*rpw(jrec)
            lm1=minus1(jrec,l)
            if (lm1.eq.0) cycle
            osv(lp1,0:m)=osv(lp1,0:m)+real(indexps(jrec,l),dp)*hpm1(zzab,ab)* &
                   & (osv(lm1,0:m)-qoppq*osv(lm1,1:m+1))
          end do
        end do
        poq=pab(zzab,ab)*pm1(zzcd,cd)
        delv=-pm1(zzcd,cd)*(zeta(zb,sh2)*rba(:)+zeta(zd,sh4)*rdc(:))
        oset(:,1)=osv(:,0)
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
            oset(l,jp1)=delv(jj)*oset(l,1)-poq*oset(lp1,1)
            if (lm1.ne.0) then
              oset(l,jp1)=oset(l,jp1)+hpm1(zzcd,cd)*indexps(jj,l)*oset(lm1,1)
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
            oset(l,jp1)=delv(jrec)*oset(l,j)-poq*oset(lp1,j)
            if (jm1.ne.0) then
              oset(l,jp1)=oset(l,jp1)+indexps(jrec,j)*hpm1(zzcd,cd)*oset(l,jm1)
            end if
            if (lm1.ne.0) then
              oset(l,jp1)=oset(l,jp1)+indexps(jrec,l)*hpm1(zzcd,cd)*oset(lm1,j)
            end if
          end do
        end do
        oshb(:,1,:)=oset(1:35,:)
        do l=5,20
          do jj=1,3
            kp1=1+jj
            lp1=plus1(jj,l)
            oshb(l,kp1,:)=oshb(lp1,1,:)+rba(jj)*oshb(l,1,:)
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
            oshb(l,kp1,:)=oshb(lp1,k,:)+rba(jrec)*oshb(l,k,:)
          end do
        end do
        oshk(:,:,:,1)=oshb(5:10,5:10,:)
        do j=1,4
          do jj=1,3
            ip1=1+jj
            jp1=plus1(jj,j)
            oshk(:,:,j,ip1)=oshk(:,:,jp1,1)+rdc(jj)*oshk(:,:,j,1) 
          end do
        end do
        lcd=djkl(zc,sh3)*djkl(zd,sh4)
        ddllc(:,:,1,1)=ddllc(:,:,1,1)+oshk(:,:,1,1)
        ddllc(:,:,1,2:4)=ddllc(:,:,1,2:4)+oshk(:,:,1,2:4)*djkl(zd,sh4)
        ddllc(:,:,2:4,1)=ddllc(:,:,2:4,1)+oshk(:,:,2:4,1)*djkl(zc,sh3)
        ddllc(:,:,2:4,2:4)=ddllc(:,:,2:4,2:4)+oshk(:,:,2:4,2:4)*lcd
      end do
    end do
  end do
end do


! Then recontract
do l=1,4
  do k=1,4
    do j=1,5
      do i=1,6
        tmpj(i,j,k,l)=sum(matharmond(j,:)*ddllc(i,:,k,l))
      end do
    end do
  end do
end do
do l=1,4
  do k=1,4
    do j=1,5
      do i=1,5
        eddll(i,j,k,l)=sum(matharmond(i,:)*tmpj(:,j,k,l))
      end do
    end do
  end do
end do

return
end subroutine i2e_ddll


pure subroutine i2e_dldl(sh1,sh2,sh3,sh4,edldl)
use nrtype ; use molprops ; use boys
implicit none

integer(i4b), intent(in) :: sh1,sh2,sh3,sh4
real(dp), intent(out) :: edldl(5,4,5,4)

integer(i4b) :: j,l,jm1,lm1,li,lf, bigL, i, ip1, jj, jp1, k, kp1 , lp1
integer(i4b) :: ab, cd, za, zb, zc, zd, zza, zzc, zzab, zzcd, m, jrec
real(dp) :: ooppq, qpq, poq, delv(3), qoppq, rpq2, fac1
real(dp) :: rba(3), rdc(3), wx(3)
real(dp) :: bfac(0:6), osv(84,0:6),oset(84,20), rap(3), rpw(3)
real(dp) :: osc(20,20), oshb(20,4,20), oshk(6,4,20,4), dldlc(6,4,6,4)
real(dp) :: afac(0:4), tmpk(6,4,5,4)

! This is my first try at doing the Giese thing - first build [ss|sk] (k has L = 8)
!   then do electron transfer, then contract, then use the HGP HR

dldlc=0._dp

ab=ind2(sh1,sh2)
cd=ind2(sh3,sh4)
rba=xsh(:,sh1)-xsh(:,sh2)
rdc=xsh(:,sh3)-xsh(:,sh4)
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
        qoppq=pab(zzcd,cd)*ooppq
        qpq=pab(zzab,ab)*qoppq
        rpq2=sum((Px(:,zzab,ab)-Px(:,zzcd,cd))**2)
        Wx=(pab(zzab,ab)*Px(:,zzab,ab)+pab(zzcd,cd)*Px(:,zzcd,cd))*ooppq
        rap=Px(:,zzab,ab)-xsh(:,sh1)
        rpw=Wx-Px(:,zzab,ab)
        ! Build (ss|ss)^m, m=0,1,2,3,4,5,6,7,8
        call boys_n(6,qpq*rpq2,bfac)
        osv(1,:)=kab(zzab,ab)*kab(zzcd,cd)*bfac*sqrt(ooppq)
        ! OS recursion: Build (ss|sp)^m, m=0-7
        do m=0,5
          osv(2:4,m)=rap(:)*osv(1,m)+rpw(:)*osv(1,m+1)
        end do
        ! OS recursion: Build (ss|sd)^m, m=0-6
        afac(0:4)=hpm1(zzab,ab)*(osv(1,0:4)-qoppq*osv(1,1:5))
        do lp1=5,10
          do jj=1,3
            l=minus1(jj,lp1)
            if (l.eq.0) cycle
            jrec=jj 
            exit
          end do
          osv(lp1,0:4)=osv(l,0:4)*rap(jrec)+osv(l,1:5)*rpw(jrec)
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
            osv(lp1,0:m)=osv(l,0:m)*rap(jrec)+osv(l,1:m+1)*rpw(jrec)
            lm1=minus1(jrec,l)
            if (lm1.eq.0) cycle
            osv(lp1,0:m)=osv(lp1,0:m)+real(indexps(jrec,l),dp)*hpm1(zzab,ab)* &
                   & (osv(lm1,0:m)-qoppq*osv(lm1,1:m+1))
          end do
        end do
        poq=pab(zzab,ab)*pm1(zzcd,cd)
        delv=-pm1(zzcd,cd)*(zeta(zb,sh2)*rba(:)+zeta(zd,sh4)*rdc(:))
        oset(:,1)=osv(:,0)
        ! Now build up the bra terms. For the HR, I need (sd|sd), (sd|sf), (sf|sd), (sf|sf)
        do l=1,56 ! First make [sp|ss] to [sp|sh], using [ss|ss] to [ss|si]
          do jj=1,3
            lm1=minus1(jj,l)
            lp1=plus1(jj,l)
            jp1=jj+1
            oset(l,jp1)=delv(jj)*oset(l,1)-poq*oset(lp1,1)
            if (lm1.ne.0) then
              oset(l,jp1)=oset(l,jp1)+hpm1(zzcd,cd)*indexps(jj,l)*oset(lm1,1)
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
            oset(l,jp1)=delv(jrec)*oset(l,j)-poq*oset(lp1,j) ! This uses from [sp|ss] to [sp|sh]
            if (jm1.ne.0) then
              oset(l,jp1)=oset(l,jp1)+indexps(jrec,j)*hpm1(zzcd,cd)*oset(l,jm1) ! And this uses [ss|ss] to [ss|sg]
            end if
            if (lm1.ne.0) then
              oset(l,jp1)=oset(l,jp1)+indexps(jrec,l)*hpm1(zzcd,cd)*oset(lm1,j) ! And this uses [sp|ss] to [ss|sf]
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
            lm1=minus1s(jrec,l)
            oset(l,jp1)=delv(jrec)*oset(l,j)-poq*oset(lp1,j)
            if (jm1.ne.0) then
              oset(l,jp1)=oset(l,jp1)+indexps(jrec,j)*hpm1(zzcd,cd)*oset(l,jm1)
            end if
            if (lm1.ne.0) then
              oset(l,jp1)=oset(l,jp1)+indexps(jrec,l)*hpm1(zzcd,cd)*oset(lm1,j)
            end if
          end do
        end do
        oshb(:,1,:)=oset(1:20,:)
        do l=5,10
          do jj=1,3
            kp1=1+jj
            lp1=plus1(jj,l)
            oshb(l,kp1,:)=oshb(lp1,1,:)+rba(jj)*oshb(l,1,:)
          end do
        end do
        oshk(:,:,:,1)=oshb(5:10,:,:)
        do j=5,10
          do jj=1,3
            ip1=1+jj
            jp1=plus1(jj,j)
            oshk(:,:,j,ip1)=oshk(:,:,jp1,1)+rdc(jj)*oshk(:,:,j,1) 
          end do
        end do
        dldlc(:,1,:,1)=dldlc(:,1,:,1)+oshk(:,1,5:10,1)
        dldlc(:,1,:,2:4)=dldlc(:,1,:,2:4)+oshk(:,1,5:10,2:4)*djkl(zd,sh4)
        dldlc(:,2:4,:,1)=dldlc(:,2:4,:,1)+oshk(:,2:4,5:10,1)*djkl(zb,sh2)
        dldlc(:,2:4,:,2:4)=dldlc(:,2:4,:,2:4)+oshk(:,2:4,5:10,2:4)*djkl(zb,sh2)*djkl(zd,sh4)
      end do
    end do
  end do
end do

! Then recontract
do l=1,4
  do k=1,5
    do j=1,4
      do i=1,6
        tmpk(i,j,k,l)=sum(matharmond(k,:)*dldlc(i,j,:,l))
      end do
    end do
  end do
end do
do l=1,4
  do k=1,5
    do j=1,4
      do i=1,5
        edldl(i,j,k,l)=sum(matharmond(i,:)*tmpk(:,j,k,l))
      end do
    end do
  end do
end do

return
end subroutine i2e_dldl


pure subroutine i2e_dddl(sh1,sh2,sh3,sh4,edddl)
use nrtype ; use molprops ; use boys
implicit none

integer(i4b), intent(in) :: sh1,sh2,sh3,sh4
real(dp), intent(out) :: edddl(5,5,5,4)

integer(i4b) :: j,l,jm1,lm1,li,lf, bigL, i, ip1, jj, jp1, k, kp1 , lp1
integer(i4b) :: ab, cd, za, zb, zc, zd, zza, zzc, zzab, zzcd, m, jrec
real(dp) :: ooppq, qpq, poq, delv(3), qoppq, rpq2, fac1
real(dp) :: rba(3), rdc(3), wx(3)
real(dp) :: bfac(0:7), osv(120,0:7),oset(120,20), rap(3), rpw(3), osetc(35,20)
real(dp) :: hgpb(35,10,6,4), hgpkc(35,6,4), hgpk(35,20,4), dddlc(6,6,6,4)
real(dp) :: afac(0:5), tmpjk(6,5,5,4), tmpk(6,6,5,4)

! This is my first try at doing the Giese thing - first build [ss|sk] (k has L = 8)
!   then do electron transfer, then contract, then use the HGP HR

hgpkc=0._dp

ab=ind2(sh1,sh2)
cd=ind2(sh3,sh4)
rba=xsh(:,sh1)-xsh(:,sh2)
rdc=xsh(:,sh3)-xsh(:,sh4)
! Inverse the loop so I do the contraction for ket part only inside two loops
do zc=1,nzet(sh3)
  zzc=maxprim*(zc-1)
  do zd=1,nzet(sh4)
    zzcd=zzc+zd
    if (negab(zzcd,cd)) cycle
    osetc=0._dp
    do za=1,nzet(sh1)
      zza=maxprim*(za-1)
      do zb=1,nzet(sh2)
        zzab=zza+zb
        if (negab(zzab,ab)) cycle
        ooppq=1._dp/(pab(zzab,ab)+pab(zzcd,cd))
        qoppq=pab(zzcd,cd)*ooppq
        qpq=pab(zzab,ab)*qoppq
        rpq2=sum((Px(:,zzab,ab)-Px(:,zzcd,cd))**2)
        Wx=(pab(zzab,ab)*Px(:,zzab,ab)+pab(zzcd,cd)*Px(:,zzcd,cd))*ooppq
        rap=Px(:,zzab,ab)-xsh(:,sh1)
        rpw=Wx-Px(:,zzab,ab)
        ! Build (ss|ss)^m, m=0,1,2,3,4,5,6,7,8
        call boys_n(7,qpq*rpq2,bfac)
        osv(1,:)=kab(zzab,ab)*kab(zzcd,cd)*bfac*sqrt(ooppq)
        ! OS recursion: Build (ss|sp)^m, m=0-7
        do m=0,6
          osv(2:4,m)=rap(:)*osv(1,m)+rpw(:)*osv(1,m+1)
        end do
        ! OS recursion: Build (ss|sd)^m, m=0-6
        afac(0:5)=hpm1(zzab,ab)*(osv(1,0:5)-qoppq*osv(1,1:6))
        do lp1=5,10
          do jj=1,3
            l=minus1(jj,lp1)
            if (l.eq.0) cycle
            jrec=jj
            exit
          end do
          osv(lp1,0:5)=osv(l,0:5)*rap(jrec)+osv(l,1:6)*rpw(jrec)
          lm1=minus1(jrec,l)
          if (lm1.eq.0) cycle
          osv(lp1,0:5)=osv(lp1,0:5)+afac(0:5)
        end do
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
            osv(lp1,0:m)=osv(l,0:m)*rap(jrec)+osv(l,1:m+1)*rpw(jrec)
            lm1=minus1(jrec,l)
            if (lm1.eq.0) cycle
            osv(lp1,0:m)=osv(lp1,0:m)+real(indexps(jrec,l),dp)*hpm1(zzab,ab)* &
                   & (osv(lm1,0:m)-qoppq*osv(lm1,1:m+1))
          end do
        end do
        poq=pab(zzab,ab)*pm1(zzcd,cd)
        delv=-pm1(zzcd,cd)*(zeta(zb,sh2)*rba(:)+zeta(zd,sh4)*rdc(:))
        oset(:,1)=osv(:,0)
        do l=1,84
          do jj=1,3
            lm1=minus1(jj,l)
            lp1=plus1(jj,l)
            jp1=jj+1
            oset(l,jp1)=delv(jj)*oset(l,1)-poq*oset(lp1,1)
            if (lm1.ne.0) then
              oset(l,jp1)=oset(l,jp1)+hpm1(zzcd,cd)*indexps(jj,l)*oset(lm1,1)
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
            oset(l,jp1)=delv(jrec)*oset(l,j)-poq*oset(lp1,j)
            if (jm1.ne.0) then
              oset(l,jp1)=oset(l,jp1)+indexps(jrec,j)*hpm1(zzcd,cd)*oset(l,jm1)
            end if
            if (lm1.ne.0) then
              oset(l,jp1)=oset(l,jp1)+indexps(jrec,l)*hpm1(zzcd,cd)*oset(lm1,j)
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
            oset(l,jp1)=delv(jrec)*oset(l,j)-poq*oset(lp1,j)
            if (jm1.ne.0) then
              oset(l,jp1)=oset(l,jp1)+indexps(jrec,j)*hpm1(zzcd,cd)*oset(l,jm1)
            end if
            if (lm1.ne.0) then
              oset(l,jp1)=oset(l,jp1)+indexps(jrec,l)*hpm1(zzcd,cd)*oset(lm1,j)
            end if
          end do
        end do
        osetc=osetc+oset(1:35,:)
      end do
    end do
    hgpk(:,:,1)=osetc
    do l=5,10
      do jj=1,3
        kp1=1+jj
        lp1=plus1(jj,l)
        hgpk(:,l,kp1)=hgpk(:,lp1,1)+rdc(jj)*hgpk(:,l,1)
      end do
    end do
    hgpkc(:,:,1)=hgpkc(:,:,1)+hgpk(:,5:10,1)
    hgpkc(:,:,2:4)=hgpkc(:,:,2:4)+hgpk(:,5:10,2:4)*djkl(zd,sh4)
  end do
end do

hgpb(:,1,:,:)=hgpkc
do l=5,20
  do jj=1,3
    kp1=1+jj
    lp1=plus1(jj,l)
    hgpb(l,kp1,:,:)=hgpb(lp1,1,:,:)+rba(jj)*hgpb(l,1,:,:)
  end do
end do
do kp1=5,10
  do l=5,10
    do jj=1,3
      k=minus1(jj,kp1)
      if (k.eq.0) cycle
      jrec=jj 
      exit
    end do
    lp1=plus1(jrec,l)
    hgpb(l,kp1,:,:)=hgpb(lp1,k,:,:)+rba(jrec)*hgpb(l,k,:,:)
  end do
end do

dddlc=hgpb(5:10,5:10,:,:)

! Then recontract
do l=1,4
  do k=1,5
    do j=1,6
      do i=1,6
        tmpk(i,j,k,l)=sum(matharmond(k,:)*dddlc(i,j,:,l))
      end do
    end do
  end do
end do
do l=1,4
  do k=1,5
    do j=1,5
      do i=1,6
        tmpjk(i,j,k,l)=sum(matharmond(j,:)*tmpk(i,:,k,l))
      end do
    end do
  end do
end do
do l=1,4
  do k=1,5
    do j=1,5
      do i=1,5
        edddl(i,j,k,l)=sum(matharmond(i,:)*tmpjk(:,j,k,l))
      end do
    end do
  end do
end do

return
end subroutine i2e_dddl




pure subroutine i2e_ssdl(sh1,sh2,sh3,sh4,essdl)
use nrtype ; use molprops ; use boys
implicit none

! This version exploits the horizontal relation, eq. 18 of HG & P 1988

integer(i4b), intent(in) :: sh1,sh2,sh3,sh4
real(dp), intent(out) :: essdl(5,4)

integer(i4b) :: ab, cd, za, zb, zc, zd, zza, zzc, zzab, zzcd, m, lp1, l, lm1
integer(i4b) :: kp1, k, jj, kk, jrec, km1
real(dp) :: ooppq, qpq, oo2q, poppq, rpq2, fac1,afac(0:1), wx(3)
real(dp) :: bfac(0:3), os(10,4,0:3), osc(6,4), rdq(3), rqw(3), rcq(3)

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
        call boys_n(3,qpq*rpq2,bfac)
        os(1,1,:)=kab(zzab,ab)*kab(zzcd,cd)*bfac*sqrt(ooppq)
        ! OS recursion: Build (ss|sp)^m, m=0,1
        do m=0,2
          os(1,2:4,m)=(rdq(:)*os(1,1,m)+rqw(:)*os(1,1,m+1))*djkl(zd,sh4)
        end do
        afac(0:1)=djkl(zd,sh4)*.5_dp*pm1(zzcd,cd)*(os(1,1,0:1)-poppq*os(1,1,1:2))
        do jj=1,3
          kp1=jj+1
          os(kp1,:,0:1)=rcq(jj)*os(1,:,0:1)+rqw(jj)*os(1,:,1:2)
          os(kp1,kp1,0:1)=os(kp1,kp1,0:1)+afac(0:1)
        end do
        do kp1=5,10
          do jj=1,3
            k=minus1(jj,kp1)
            if(k.eq.0) cycle
            jrec=jj
            exit
          end do
          os(kp1,:,0)=os(k,:,0)*rcq(jrec)+os(k,:,1)*rqw(jrec)
          km1=minus1(jrec,k)
          if (km1.ne.0) then
            os(kp1,:,0)=os(kp1,:,0)+.5_dp*pm1(zzcd,cd)*(os(km1,:,0)-poppq*os(km1,:,1))
          end if
          os(kp1,jrec+1,0)=os(kp1,jrec+1,0)+djkl(zd,sh4)*.5_dp*pm1(zzcd,cd)*(os(k,1,0)-poppq*os(k,1,1))
        end do
        osc=osc+os(5:10,:,0)
      end do
    end do
  end do
end do


! Finally do contraction to form spherical harmonics

do jj=1,5
do kk=1,4
  essdl(jj,kk)=sum(matharmond(jj,:)*osc(:,kk))
end do
end do

end subroutine i2e_ssdl





pure subroutine i2e_sdsl(sh1,sh2,sh3,sh4,esdsl)
use nrtype ; use molprops ; use boys
implicit none

integer(i4b), intent(in) :: sh1,sh2,sh3,sh4
real(dp), intent(out) :: esdsl(5,4)

integer(i4b) :: j,l,jm1,lm1, lp1, jj, jp1, kk, jrec
integer(i4b) :: ab, cd, za, zb, zc, zd, zza, zzc, zzab, zzcd, m
real(dp) :: ooppq, qpq, oo2q, qoppq, rpq2, fac1, wx(3)
real(dp) :: bfac(0:3), os(10,4,0:3), rbp(3), rpw(3), rdq(3), rqw(3)
real(dp) :: afac(0:1)
real(dp) :: ketc(6,4),osc(6,4)

! This routine uses O-S recursion

osc=0._dp

ab=ind2(sh1,sh2)
cd=ind2(sh3,sh4)
do zc=1,nzet(sh3)
   zzc=maxprim*(zc-1)
   do zd=1,nzet(sh4)
      zzcd=zzc+zd
      if (negab(zzcd,cd)) cycle
      rdq=(Px(:,zzcd,cd)-xsh(:,sh4))
      ketc=0._dp
      do za=1,nzet(sh1)
         zza=maxprim*(za-1)
         do zb=1,nzet(sh2)
            zzab=zza+zb
            if (negab(zzab,ab)) cycle
            ooppq=1._dp/(pab(zzab,ab)+pab(zzcd,cd))
            qoppq=pab(zzcd,cd)*ooppq
            qpq=pab(zzab,ab)*qoppq
            rpq2=sum((Px(:,zzab,ab)-Px(:,zzcd,cd))**2)
            Wx=(pab(zzab,ab)*Px(:,zzab,ab)+pab(zzcd,cd)*Px(:,zzcd,cd))*ooppq
            rbp=Px(:,zzab,ab)-xsh(:,sh2)
            rpw=Wx-Px(:,zzab,ab)
            rqw=(Wx-Px(:,zzcd,cd))
            ! Build (ss|ss)^m, m=0,1,2,3,4
            call boys_n(3,qpq*rpq2,bfac)
            os(1,1,:)=kab(zzab,ab)*kab(zzcd,cd)*bfac*sqrt(ooppq)
            ! OS recursion: Build (ss|sp)^m, m=0,3
            do m=0,2
              os(1,2:4,m)=(rdq(:)*os(1,1,m)+rqw(:)*os(1,1,m+1))
            end do
            afac=.5_dp*ooppq*os(1,1,1:2)
            do jj=1,3
              jp1=jj+1
              os(jp1,:,0:1)=os(1,:,0:1)*rbp(jj)+os(1,:,1:2)*rpw(jj)
              os(jp1,jp1,0:1)=os(jp1,jp1,0:1)+afac
            end do
            ! OS recursion: Build (sd|sl)
            do jp1=5,10
              do jj=1,3
                j=minus1s(jj,jp1)
                if (j.eq.0) cycle
                jrec=jj
                exit
              end do
              os(jp1,:,0)=rbp(jrec)*os(j,:,0)+rpw(jrec)*os(j,:,1)
              jm1=minus1s(jrec,j)
              if (jm1.ne.0) then
                os(jp1,:,0)=os(jp1,:,0)+hpm1(zzab,ab)*(os(jm1,:,0)-qoppq*os(jm1,:,1))
              end if
              l=jrec+1
              os(jp1,l,0)=os(jp1,l,0)+.5_dp*ooppq*os(j,1,1)
            end do
            ketc=ketc+os(5:10,:,0)
         end do
      end do
      osc(:,1)=osc(:,1)+ketc(:,1)
      osc(:,2:4)=osc(:,2:4)+ketc(:,2:4)*djkl(zd,sh4)
   end do
end do

! Finally do contraction to form spherical harmonics
! This is not very efficient but very simple

do jj=1,5
do kk=1,4
  esdsl(jj,kk)=sum(matharmond(jj,:)*osc(:,kk))
end do
end do

end subroutine i2e_sdsl



pure subroutine i2e_sldl(sh1,sh2,sh3,sh4,esldl)
use nrtype ; use molprops ; use boys
implicit none

integer(i4b), intent(in) :: sh1,sh2,sh3,sh4
real(dp), intent(out) :: esldl(4,5,4)

integer(i4b) :: j,l,jm1,lm1,lp1, k, kp1, bigL, jj, jp1, li, lf, jrec
integer(i4b) :: ab, cd, za, zb, zc, zd, zza, zzc, zzab, zzcd, m
real(dp) :: ooppq, qpq, qop, poppq, rpq2
real(dp) :: bfac(0:4), osv(35,0:4),oset(4,35),osh(4,20,4), rcq(3), rqw(3)
real(dp) :: osc(4,6,4)
real(dp) :: afac(0:2)
real(dp) :: wx(3), rab(3), rdc(3), delv(3)

! Use the VRR to make [ss|si] (l=6), then ET to make [sd|sg], then HRR

osc=0._dp
oset=0._dp
osh=0._dp

ab=ind2(sh1,sh2)
cd=ind2(sh3,sh4)
rdc=xsh(:,sh3)-xsh(:,sh4)
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
        rcq=Px(:,zzcd,cd)-xsh(:,sh3)
        rqw=Wx-Px(:,zzcd,cd)
        call boys_n(4,qpq*rpq2,bfac)
        osv(1,:)=kab(zzab,ab)*kab(zzcd,cd)*bfac*sqrt(ooppq)
        do m=0,3
          osv(2:4,m)=rcq(:)*osv(1,m)+rqw(:)*osv(1,m+1)
        end do
        ! OS recursion: Build (ss|ds)^m, m=0-2
        afac(0:2)=(osv(1,0:2)-poppq*osv(1,1:3))*hpm1(zzcd,cd)
        do lp1=5,10
          do jj=1,3
            l=minus1(jj,lp1)
            if (l.eq.0) cycle
            jrec=jj
            exit
          end do
          osv(lp1,0:2)=osv(l,0:2)*rcq(jrec)+osv(l,1:3)*rqw(jrec)
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
            osv(lp1,0:m)=osv(l,0:m)*rcq(jrec)+osv(l,1:m+1)*rqw(jrec)
            lm1=minus1(jrec,l)
            if (lm1.eq.0) cycle
            osv(lp1,0:m)=osv(lp1,0:m)+real(indexps(jrec,l),dp)*hpm1(zzcd,cd)* &
                   & (osv(lm1,0:m)-poppq*osv(lm1,1:m+1))
          end do
        end do
        qop=pab(zzcd,cd)*pm1(zzab,ab)
        delv=pm1(zzab,ab)*(zeta(za,sh1)*rab(:)+zeta(zd,sh4)*rdc(:))
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
            oset(jp1,l)=oset(jp1,l)*djkl(zb,sh2)
          end do
        end do
        ! Now use HRR, (sx|dp)=(sx|fs)+rdc*(sx|ds)
        osh(:,:,1)=oset(:,1:20)
        do k=5,10
          do jj=1,3
            lp1=jj+1
            kp1=plus1(jj,k)
            osh(:,k,lp1)=djkl(zd,sh4)*(osh(:,kp1,1)+rdc(jj)*osh(:,k,1))
          end do
        end do
        osc=osc+osh(:,5:10,:)
      end do
    end do
  end do
end do

do l=1,4
  do k=1,5
    do j=1,4
      esldl(j,k,l)=sum(matharmond(k,:)*osc(j,:,l))
    end do
  end do
end do

return
end subroutine i2e_sldl


pure subroutine i2e_sdll(sh1,sh2,sh3,sh4,esdll)
use nrtype ; use molprops ; use boys
implicit none

integer(i4b), intent(in) :: sh1,sh2,sh3,sh4
real(dp), intent(out) :: esdll(5,4,4)

integer(i4b) :: j,l,jm1,ip1, k, kp1, bigL, jj, jp1, li, lf, lp1, lm1, i
integer(i4b) :: ab, cd, za, zb, zc, zd, zza, zzc, zzab, zzcd, m, jrec
real(dp) :: ooppq, qpq, qop, poppq, rpq2
real(dp) :: bfac(0:4), osv(35,0:4),oset(10,35), rdq(3), rqw(3)
real(dp) :: osh(6,4,10), osc(6,4,4)
real(dp) :: afac(0:2)
real(dp) :: wx(3), rab(3), rcd(3), delv(3), ggll

! Use the VRR to make [ss|si] (l=6), then ET to make [sd|sg], then HRR

osh=0._dp
oset=0._dp
osc=0._dp

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
        call boys_n(4,qpq*rpq2,bfac)
        osv(1,:)=kab(zzab,ab)*kab(zzcd,cd)*bfac*sqrt(ooppq)
        do m=0,3
          osv(2:4,m)=rdq(:)*osv(1,m)+rqw(:)*osv(1,m+1)
        end do
        ! OS recursion: Build (ss|sd)^m, m=0-4
        afac(0:2)=.5_dp*(osv(1,0:2)-poppq*osv(1,1:3))*pm1(zzcd,cd)
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
            osv(lp1,0:m)=osv(lp1,0:m)+real(indexps(jrec,l),dp)*.5_dp*pm1(zzcd,cd)* &
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
        osh(:,1,:)=oset(5:10,1:10)
        do j=1,4
          do jj=1,3
            ip1=jj+1
            jp1=plus1(jj,j)
            osh(:,ip1,j)=osh(:,1,jp1)+rcd(jj)*osh(:,1,j)
          end do
        end do
        ggll=djkl(zc,sh3)*djkl(zd,sh4)
        osc(:,1,1)=osc(:,1,1)+osh(:,1,1)
        osc(:,2:4,1)=osc(:,2:4,1)+osh(:,2:4,1)*djkl(zc,sh3)
        osc(:,1,2:4)=osc(:,1,2:4)+osh(:,1,2:4)*djkl(zd,sh4)
        osc(:,2:4,2:4)=osc(:,2:4,2:4)+osh(:,2:4,2:4)*ggll
      end do
    end do
  end do
end do


! Finally do contraction to form spherical harmonics, using temporary arrays tmp(6,6,6) and tmp2
! This is not very efficient but very simple
!First take only the Cartesians:

do l=1,4
  do k=1,4
    do j=1,5
      esdll(j,k,l)=sum(matharmond(j,:)*osc(:,k,l))
    end do
  end do
end do

return
end subroutine i2e_sdll



pure subroutine i2e_sldd(sh1,sh2,sh3,sh4,esldd)
use nrtype ; use molprops ; use boys
implicit none

integer(i4b), intent(in) :: sh1,sh2,sh3,sh4
real(dp), intent(out) :: esldd(4,5,5)

integer(i4b) :: j,l,jm1,lm1,lp1, k, kp1, bigL, jj, jp1, li, lf, jrec
integer(i4b) :: ab, cd, za, zb, zc, zd, zza, zzc, zzab, zzcd, m
real(dp) :: ooppq, qpq, qop, poppq, rpq2
real(dp) :: bfac(0:5), osv(56,0:5),oset(4,56), rdq(3), rqw(3)
real(dp) :: osc(4,35), hgp(4,10,35), slddc(4,6,6), tmpl(4,6,5)
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
        call boys_n(5,qpq*rpq2,bfac)
        osv(1,:)=kab(zzab,ab)*kab(zzcd,cd)*bfac*sqrt(ooppq)
        do m=0,4
          osv(2:4,m)=rdq(:)*osv(1,m)+rqw(:)*osv(1,m+1)
        end do
        ! OS recursion: Build (ss|sd)^m, m=0-4
        afac(0:3)=.5_dp*(osv(1,0:3)-poppq*osv(1,1:4))*pm1(zzcd,cd)
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
        osc(1,:)=osc(1,:)+oset(1,1:35)
        osc(2:4,:)=osc(2:4,:)+oset(2:4,1:35)*djkl(zb,sh2)
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
slddc=hgp(:,5:10,5:10)

do l=1,5
  do k=1,6
    do j=1,4
      tmpl(j,k,l)=sum(matharmond(l,:)*slddc(j,k,:))
    end do
  end do
end do
do l=1,5
  do k=1,5
    do j=1,4
      esldd(j,k,l)=sum(matharmond(k,:)*tmpl(j,:,l))
    end do
  end do
end do

return
end subroutine i2e_sldd





pure subroutine i2e_sddl(sh1,sh2,sh3,sh4,esddl)
use nrtype ; use molprops ; use boys
implicit none

integer(i4b), intent(in) :: sh1,sh2,sh3,sh4
real(dp), intent(out) :: esddl(5,5,4)

integer(i4b) :: j,l,jm1,lm1,lp1, k, kp1, bigL, jj, jp1, li, lf, jrec
integer(i4b) :: ab, cd, za, zb, zc, zd, zza, zzc, zzab, zzcd, m
real(dp) :: ooppq, qpq, qop, poppq, rpq2
real(dp) :: bfac(0:5), osv(56,0:5),oset(10,56), rcq(3), rqw(3)
real(dp) :: osh(6,20,4), sddlc(6,6,4), tmpl(6,5,4)
real(dp) :: afac(0:3)
real(dp) :: wx(3), rab(3), rdc(3), delv(3)

! Use the VRR to make [ss|si] (l=6), then ET to make [sd|sg], then HRR

sddlc=0._dp
oset=0._dp
osh=0._dp

ab=ind2(sh1,sh2)
cd=ind2(sh3,sh4)
rdc=xsh(:,sh3)-xsh(:,sh4)
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
        rcq=Px(:,zzcd,cd)-xsh(:,sh3)
        rqw=Wx-Px(:,zzcd,cd)
        ! Build (ss|ss)^m, m=0,1,2,3,4,5,6
        call boys_n(5,qpq*rpq2,bfac)
        osv(1,:)=kab(zzab,ab)*kab(zzcd,cd)*bfac*sqrt(ooppq)
        ! OS recursion: Build (ss|sp)^m, m=0-5
        do m=0,4
          osv(2:4,m)=rcq(:)*osv(1,m)+rqw(:)*osv(1,m+1)
        end do
        ! OS recursion: Build (ss|sd)^m, m=0-4
        afac(0:3)=.5_dp*(osv(1,0:3)-poppq*osv(1,1:4))*pm1(zzcd,cd)
        do lp1=5,10
          do jj=1,3
            l=minus1(jj,lp1)
            if (l.eq.0) cycle
            jrec=jj
            exit
          end do
          osv(lp1,0:3)=osv(l,0:3)*rcq(jrec)+osv(l,1:4)*rqw(jrec)
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
            osv(lp1,0:m)=osv(l,0:m)*rcq(jrec)+osv(l,1:m+1)*rqw(jrec)
            lm1=minus1(jrec,l)
            if (lm1.eq.0) cycle
            osv(lp1,0:m)=osv(lp1,0:m)+real(indexps(jrec,l),dp)*.5_dp*pm1(zzcd,cd)* &
                   & (osv(lm1,0:m)-poppq*osv(lm1,1:m+1))
          end do
        end do
        qop=pab(zzcd,cd)*pm1(zzab,ab)
        delv=pm1(zzab,ab)*(zeta(za,sh1)*rab(:)+zeta(zd,sh4)*rdc(:))
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
        osh(:,:,1)=oset(5:10,1:20)
        do k=5,10
          do jj=1,3
            lp1=jj+1
            kp1=plus1(jj,k)
            osh(:,k,lp1)=osh(:,kp1,1)+rdc(jj)*osh(:,k,1)
          end do
        end do
        sddlc(:,:,1)=sddlc(:,:,1)+osh(:,5:10,1)
        sddlc(:,:,2:4)=sddlc(:,:,2:4)+osh(:,5:10,2:4)*djkl(zd,sh4)
      end do
    end do
  end do
end do


do l=1,4
  do k=1,5
    do j=1,6
      tmpl(j,k,l)=sum(matharmond(k,:)*sddlc(j,:,l))
    end do
  end do
end do
do l=1,4
  do k=1,5
    do j=1,5
      esddl(j,k,l)=sum(matharmond(j,:)*tmpl(:,k,l))
    end do
  end do
end do

return
end subroutine i2e_sddl


pure subroutine i2e_sssl(sh1,sh2,sh3,sh4,esssl)
use nrtype ; use molprops ; use boys
implicit none

! Here I calculate (ss|sl), (ss|ll), (sl|sl), (sl|ll) and (ll|ll) shell integrals.
! Can't be bothered to work out how to do head-gordon-pople or Giese-like ET shift, so everything is just O-S

integer(i4b), intent(in) :: sh1,sh2,sh3,sh4
real(dp), intent(out) :: esssl(4)

integer(i4b) :: ab, cd, za, zb, zc, zd, zza, zzc, zzab, zzcd
real(dp) :: ooppq, qpq, poppq, rpq2, vec(3)
real(dp) :: bfac(0:1), os(0:1), rdq(3), rqw(3), wx(3), osc(4)

ab=ind2(sh1,sh2)
cd=ind2(sh3,sh4)
esssl=0._dp
do zd=1,nzet(sh4)
  osc=0._dp
  do zc=1,nzet(sh3)
    zzc=maxprim*(zc-1)
    zzcd=zzc+zd
    if (negab(zzcd,cd)) cycle
    rdq=(Px(:,zzcd,cd)-xsh(:,sh4))
    vec=pab(zzcd,cd)*Px(:,zzcd,cd)
    do za=1,nzet(sh1)
      zza=maxprim*(za-1)
      do zb=1,nzet(sh2)
        zzab=zza+zb
        if (negab(zzab,ab)) cycle
        ooppq=1._dp/(pab(zzab,ab)+pab(zzcd,cd))
        qpq=pab(zzab,ab)*pab(zzcd,cd)*ooppq
        !rpq2=sum((Px(:,zzab,ab)-Px(:,zzcd,cd))**2)
        rpq2=(Px(1,zzab,ab)-Px(1,zzcd,cd))**2+(Px(2,zzab,ab)-Px(2,zzcd,cd))**2+&
               (Px(3,zzab,ab)-Px(3,zzcd,cd))**2
        call boys_1(qpq*rpq2,bfac)
        ! Build (ss|ss)^m, m=0,1
        os(0:1)=kab(zzab,ab)*kab(zzcd,cd)*bfac(0:1)*sqrt(ooppq)
        Wx=(pab(zzab,ab)*Px(:,zzab,ab)+vec)*ooppq
        rqw=(Wx-Px(:,zzcd,cd))
        ! Build (ss|sp) part
        osc(1)=osc(1)+os(0)
        osc(2:4)=osc(2:4)+rdq(:)*os(0)+rqw(:)*os(1)
      end do
    end do
  end do
  esssl(1)=esssl(1)+osc(1)
  esssl(2:4)=esssl(2:4)+osc(2:4)*djkl(zd,sh4)
end do
return

end subroutine i2e_sssl



pure subroutine i2e_ssll(sh1,sh2,sh3,sh4,essll)
use nrtype ; use molprops ; use boys
implicit none

! This version uses standard OS VRR

integer(i4b), intent(in) :: sh1,sh2,sh3,sh4
real(dp), intent(out) :: essll(4,4)

integer(i4b) :: ab, cd, za, zb, zc, zd, zza, zzc, zzab, zzcd, m
integer(i4b) :: k, l, km1, lm1, jj
real(dp) :: ooppq,poppq, rpq2, gg, fac, qpq
real(dp) :: bfac(0:2), os(4,4,0:2), osc(4,4), rcq(3), rdq(3), rqw(3), wx(3)

! First assemble the [ss|sx] (x=s,p,d) integrals and contract them

essll=0._dp

ab=ind2(sh1,sh2)
cd=ind2(sh3,sh4)
do zc=1,nzet(sh3)
  zzc=maxprim*(zc-1)
  do zd=1,nzet(sh4)
    zzcd=zzc+zd
    if (negab(zzcd,cd)) cycle
    rcq=(Px(:,zzcd,cd)-xsh(:,sh3))
    rdq=(Px(:,zzcd,cd)-xsh(:,sh4))
    osc=0._dp
    do za=1,nzet(sh1)
      zza=maxprim*(za-1)
      do zb=1,nzet(sh2)
        zzab=zza+zb
        if (negab(zzab,ab)) cycle
        ooppq=1._dp/(pab(zzab,ab)+pab(zzcd,cd))
        poppq=pab(zzab,ab)*ooppq
        qpq=pab(zzcd,cd)*poppq
        rpq2=sum((Px(:,zzab,ab)-Px(:,zzcd,cd))**2)
        call boys_2(qpq*rpq2,bfac)
        Wx=(pab(zzab,ab)*Px(:,zzab,ab)+pab(zzcd,cd)*Px(:,zzcd,cd))*ooppq
        rqw=(Wx-Px(:,zzcd,cd))
        ! Build (ss|ss)^m, m=0,1,2
        os(1,1,0:2)=kab(zzab,ab)*kab(zzcd,cd)*bfac(0:2)*sqrt(ooppq)
        ! OS recursion: Build (ss|sp)^m, m=0,1
        do m=0,1
          os(1,2:4,m)=rdq(:)*os(1,1,m)+rqw(:)*os(1,1,m+1)
          os(2:4,1,m)=rcq(:)*os(1,1,m)+rqw(:)*os(1,1,m+1)
        end do
        ! OS recursion: Build (ss|pp) parts
        fac=(os(1,1,0)-poppq*os(1,1,1))*hpm1(zzcd,cd)
        do jj=1,3
          k=1+jj
          os(k,2:4,0)=rcq(jj)*os(1,2:4,0)+rqw(jj)*os(1,2:4,1)
          os(k,k,0)=os(k,k,0)+fac
        end do
        osc=osc+os(:,:,0)
      end do
    end do
    essll(1,1)=essll(1,1)+osc(1,1)
    essll(1,2:4)=essll(1,2:4)+osc(1,2:4)*djkl(zd,sh4)
    essll(2:4,1)=essll(2:4,1)+osc(2:4,1)*djkl(zc,sh3)
    gg=djkl(zc,sh3)*djkl(zd,sh4)
    essll(2:4,2:4)=essll(2:4,2:4)+osc(2:4,2:4)*gg
  end do
end do

end subroutine i2e_ssll




pure subroutine i2e_slsl(sh1,sh2,sh3,sh4,eslsl)
use nrtype ; use molprops ; use boys
implicit none

integer(i4b), intent(in) :: sh1,sh2,sh3,sh4
real(dp), intent(out) :: eslsl(4,4)

integer(i4b) :: ab, cd, za, zb, zc, zd, zza, zzc, zzab, zzcd, m
integer(i4b) :: j, jm1, l, jj
real(dp) :: ooppq, qpq, ooq, oop, poppq, qoppq, rpq2, fac1, wx(3)
real(dp) :: sterm(0:2), pterml(3,0:1), ptermj(3)
real(dp) :: bfac(0:2), osc(4,4), os(4,4), rbp(3), rpw(3), rdq(3), rqw(3), gg

! Simply assemble the [sx|sy] (x,y=s,p) integrals and contract them

eslsl=0._dp

ab=ind2(sh1,sh2)
cd=ind2(sh3,sh4)
do zb=1,nzet(sh2)
  do zd=1,nzet(sh4)
    osc=0._dp
    do za=1,nzet(sh1)
      zzab=maxprim*(za-1)+zb
      rbp=(Px(:,zzab,ab)-xsh(:,sh2))
      do zc=1,nzet(sh3)
        zzcd=maxprim*(zc-1)+zd
        ooppq=1._dp/(pab(zzab,ab)+pab(zzcd,cd))
        qoppq=pab(zzcd,cd)*ooppq
        qpq=pab(zzab,ab)*qoppq
        rpq2=sum((Px(:,zzab,ab)-Px(:,zzcd,cd))**2)
        Wx=(pab(zzab,ab)*Px(:,zzab,ab)+pab(zzcd,cd)*Px(:,zzcd,cd))*ooppq
        rdq=(Px(:,zzcd,cd)-xsh(:,sh4))
        rpw=(Wx-Px(:,zzab,ab))
        rqw=(Wx-Px(:,zzcd,cd))
        ! Build (ss|ss)^m, m=0,1,2
        call boys_2(qpq*rpq2,bfac)
        sterm=kab(zzab,ab)*kab(zzcd,cd)*bfac*sqrt(ooppq)
        os(1,1)=sterm(0)
        ! OS recursion: Build (ss|sp)^m, m=0,1
        do m=0,1
          pterml(:,m)=(rdq(:)*sterm(m)+rqw(:)*sterm(m+1))
        end do
        os(1,2:4)=pterml(:,0)
        os(2:4,1)=(rbp(:)*sterm(0)+rpw(:)*sterm(1))
        ! OS recursion: Build (sp|sp)
        fac1=.5_dp*sterm(1)*ooppq
        do j=1,3
          os(1+j,2:4)=(rbp(j)*pterml(:,0)+rpw(j)*pterml(:,1))
          os(j+1,j+1)=os(j+1,j+1)+fac1
        end do
        osc=osc+os(:,:)
      end do
    end do
    gg=djkl(zb,sh2)*djkl(zd,sh4)
    eslsl(1,1)=eslsl(1,1)+osc(1,1)
    eslsl(1,2:4)=eslsl(1,2:4)+osc(1,2:4)*djkl(zd,sh4)
    eslsl(2:4,1)=eslsl(2:4,1)+osc(2:4,1)*djkl(zb,sh2)
    eslsl(2:4,2:4)=eslsl(2:4,2:4)+osc(2:4,2:4)*gg
  end do
end do

end subroutine i2e_slsl


pure subroutine i2e_slll(sh1,sh2,sh3,sh4,eslll)
use nrtype ; use molprops ; use boys
implicit none

! This version simply uses O-S recursion

integer(i4b), intent(in) :: sh1,sh2,sh3,sh4
real(dp), intent(out) :: eslll(4,4,4)

integer(i4b) :: ab, cd, za, zb, zc, zd, zza, zzc, zzab, zzcd, m, jj, k
real(dp) :: ooppq, qpq, poppq, qoppq, rpq2
real(dp) :: ggjl, ggjk, ggkl, fackl(0:1), facjl(0:1), facjk(0:1), vec(3)
real(dp) :: sterm(0:3), ptermj(3,0:2),ptermk(3,0:2), pterml(3,0:2)

real(dp) :: bfac(0:3), os(4,4,4,0:1), osc(4,4,4), rpw(3), rbp(3), rcq(3),rdq(3), rqw(3), wx(3)

eslll=0._dp

ab=ind2(sh1,sh2)
cd=ind2(sh3,sh4)
do zc=1,nzet(sh3)
  zzc=maxprim*(zc-1)
  do zd=1,nzet(sh4)
    zzcd=zzc+zd
    if (negab(zzcd,cd)) cycle
    rcq=Px(:,zzcd,cd)-xsh(:,sh3)
    rdq=Px(:,zzcd,cd)-xsh(:,sh4)
    vec=pab(zzcd,cd)*Px(:,zzcd,cd)
    osc=0._dp
    do za=1,nzet(sh1)
      zza=maxprim*(za-1)
      do zb=1,nzet(sh2)
        zzab=zza+zb    ! 1
        if (negab(zzab,ab)) cycle
        rbp=Px(:,zzab,ab)-xsh(:,sh2)      ! 3
        ooppq=1._dp/(pab(zzab,ab)+pab(zzcd,cd))  ! 2
        qoppq=pab(zzcd,cd)*ooppq ! 1
        poppq=pab(zzab,ab)*ooppq ! 1
        qpq=pab(zzab,ab)*qoppq ! 1
        rpq2=sum((Px(:,zzab,ab)-Px(:,zzcd,cd))**2) ! 8
        Wx=(pab(zzab,ab)*Px(:,zzab,ab)+vec)*ooppq  ! 12
        rpw=Wx-Px(:,zzab,ab) ! 3
        rqw=Wx-Px(:,zzcd,cd) ! 3
        ! Build (ss|ss)^m, m=0,1,2,3
        call boys_n(3,qpq*rpq2,bfac)  ! 4xN
        sterm(:)=kab(zzab,ab)*kab(zzcd,cd)*bfac*sqrt(ooppq)  !   1 element out of 64 done
        os(1,1,1,0)=sterm(0)
        ! OS recursion: Build (ss|sp)^m, m=0,2
        do m=0,2
          pterml(:,m)=(rdq(:)*sterm(m)+rqw(:)*sterm(m+1))                !  4 / 64
          ptermk(:,m)=(rcq(:)*sterm(m)+rqw(:)*sterm(m+1))                !  7 / 64
          ptermj(:,m)=(rbp(:)*sterm(m)+rpw(:)*sterm(m+1))   ! 10 / 64
        end do
        os(1,1,2:4,0)=pterml(:,0)
        os(1,2:4,1,0)=ptermk(:,0)
        os(2:4,1,1,0)=ptermj(:,0)
        ! OS recursion: Build (ss|pp), (sp|sp) and (sp|ps)
        fackl=(sterm(0:1)-poppq*sterm(1:2))*hpm1(zzcd,cd)
        facjl=sterm(1:2)*.5_dp*ooppq
        facjk=sterm(1:2)*.5_dp*ooppq
        do jj=1,3
          k=jj+1
          os(1,k,2:4,0:1)=(rcq(jj)*pterml(:,0:1)+rqw(jj)*pterml(:,1:2))               ! 19
          os(1,k,k,0:1)=os(1,k,k,0:1)+fackl
          os(k,1,2:4,0:1)=(rbp(jj)*pterml(:,0:1)+rpw(jj)*pterml(:,1:2))  ! 28
          os(k,1,k,0:1)=os(k,1,k,0:1)+facjl
          os(k,2:4,1,0:1)=(rbp(jj)*ptermk(:,0:1)+rpw(jj)*ptermk(:,1:2))  ! 37
          os(k,k,1,0:1)=os(k,k,1,0:1)+facjk
        end do
        ! Last step: build the (sp|pp) bit
        do jj=1,3
          k=jj+1
          os(k,2:4,2:4,0)=(os(1,2:4,2:4,0)*rbp(jj)+os(1,2:4,2:4,1)*rpw(jj))  ! 64
          os(k,2:4,k,0)=os(k,2:4,k,0)+.5_dp*ooppq*ptermk(:,1)
          os(k,k,2:4,0)=os(k,k,2:4,0)+.5_dp*ooppq*pterml(:,1)
        end do
        osc(1,:,:)=osc(1,:,:)+os(1,:,:,0)
        osc(2:4,:,:)=osc(2:4,:,:)+os(2:4,:,:,0)*djkl(zb,sh2)
      end do
    end do
    ggkl=djkl(zc,sh3)*djkl(zd,sh4)
    eslll(:,1,1)=eslll(:,1,1)+osc(:,1,1)
    eslll(:,1,2:4)=eslll(:,1,2:4)+djkl(zd,sh4)*osc(:,1,2:4)
    eslll(:,2:4,1)=eslll(:,2:4,1)+djkl(zc,sh3)*osc(:,2:4,1)
    eslll(:,2:4,2:4)=eslll(:,2:4,2:4)+ggkl*osc(:,2:4,2:4)
  end do
end do

end subroutine i2e_slll



pure subroutine i2e_llll(sh1,sh2,sh3,sh4,ellll)
use nrtype ; use molprops ; use boys
implicit none

! This version uses only O-S recursion

integer(i4b), intent(in) :: sh1,sh2,sh3,sh4
real(dp), intent(out) :: ellll(4,4,4,4)

integer(i4b) :: ab, cd, za, zb, zc, zd, zza, zzc, zzab, zzcd, m, jj, k
real(dp) :: ooppq, hooppq, qpq, poppq, qoppq, rpq2, fac1(3)
real(dp) :: bfac(0:4), osc(4,4,4,4),os(4,4,4,4,0:4), rap(3), rpw(3), rbp(3), rcq(3), rdq(3), rqw(3)
real(dp) :: ggkl, ggij, wx(3)

! First assemble the [sx|sy] (x=s,p; y=s,p,d) integrals and contract them
! We use the matrix os(j1,j2,l1,l2,m) depending on the two angular momenta directions (0 = none,
! 1 = x, 2 = y, 3 = z) at b, the two at d, and the Boys order m.

ellll=0._dp

ab=ind2(sh1,sh2)
cd=ind2(sh3,sh4)
do za=1,nzet(sh1)
  zza=maxprim*(za-1)
  do zb=1,nzet(sh2)
    zzab=zza+zb
    if (negab(zzab,ab)) cycle
    rap=Px(:,zzab,ab)-xsh(:,sh1)
    rbp=Px(:,zzab,ab)-xsh(:,sh2)
    osc=0._dp
    do zc=1,nzet(sh3)
      zzc=maxprim*(zc-1)
      do zd=1,nzet(sh4)
        zzcd=zzc+zd
        if (negab(zzcd,cd)) cycle
        ooppq=1._dp/(pab(zzab,ab)+pab(zzcd,cd))
        hooppq=.5_dp*ooppq
        qoppq=pab(zzcd,cd)*ooppq
        poppq=pab(zzab,ab)*ooppq
        qpq=pab(zzab,ab)*qoppq
        rpq2=sum((Px(:,zzab,ab)-Px(:,zzcd,cd))**2)
        Wx=(pab(zzab,ab)*Px(:,zzab,ab)+pab(zzcd,cd)*Px(:,zzcd,cd))*ooppq
        rcq=Px(:,zzcd,cd)-xsh(:,sh3)
        rdq=Px(:,zzcd,cd)-xsh(:,sh4)
        rpw=Wx-Px(:,zzab,ab)
        rqw=Wx-Px(:,zzcd,cd)
        ! Build (ss|ss)^m, m=0,1,2,3,4
        call boys_n(4,qpq*rpq2,bfac)
        os(1,1,1,1,:)=kab(zzab,ab)*kab(zzcd,cd)*bfac*sqrt(ooppq)
        ! OS recursion: Build (ss|sp)^m, m=0,3
        do m=0,3
          os(1,1,1,2:4,m)=(rdq(:)*os(1,1,1,1,m)+rqw(:)*os(1,1,1,1,m+1))
          os(1,1,2:4,1,m)=(rcq(:)*os(1,1,1,1,m)+rqw(:)*os(1,1,1,1,m+1))
          os(1,2:4,1,1,m)=(rbp(:)*os(1,1,1,1,m)+rpw(:)*os(1,1,1,1,m+1))
          os(2:4,1,1,1,m)=(rap(:)*os(1,1,1,1,m)+rpw(:)*os(1,1,1,1,m+1))
        end do
        ! OS recursion: Build (sp|sp)^m, m=0,1,2 and other s,s,p,p terms
        fac1=hooppq*os(1,1,1,1,1:3)
        do jj=1,3
          k=jj+1
          os(1,k,1,2:4,0:2)=(rbp(jj)*os(1,1,1,2:4,0:2)+rpw(jj)*os(1,1,1,2:4,1:3))
          os(1,k,1,k,0:2)=os(1,k,1,k,0:2)+fac1
          os(k,1,1,2:4,0:2)=(rap(jj)*os(1,1,1,2:4,0:2)+rpw(jj)*os(1,1,1,2:4,1:3))
          os(k,1,1,k,0:2)=os(k,1,1,k,0:2)+fac1
          os(1,k,2:4,1,0:2)=(rbp(jj)*os(1,1,2:4,1,0:2)+rpw(jj)*os(1,1,2:4,1,1:3))
          os(1,k,k,1,0:2)=os(1,k,k,1,0:2)+fac1
          os(k,1,2:4,1,0:2)=(rap(jj)*os(1,1,2:4,1,0:2)+rpw(jj)*os(1,1,2:4,1,1:3))
          os(k,1,k,1,0:2)=os(k,1,k,1,0:2)+fac1

          os(1,1,k,2:4,0:2)=(rcq(jj)*os(1,1,1,2:4,0:2)+rqw(jj)*os(1,1,1,2:4,1:3))
          os(1,1,k,k,0:2)=os(1,1,k,k,0:2)+(os(1,1,1,1,0:2)-poppq*os(1,1,1,1,1:3))*hpm1(zzcd,cd)

          os(k,2:4,1,1,0:2)=(rap(jj)*os(1,2:4,1,1,0:2)+rpw(jj)*os(1,2:4,1,1,1:3))
          os(k,k,1,1,0:2)=os(k,k,1,1,0:2)+(os(1,1,1,1,0:2)-qoppq*os(1,1,1,1,1:3))*hpm1(zzab,ab)
        end do
        ! OS recursion: Build 4 different (sp|pp)^m, m=0,1 parts. Only need m=1 for pne of them
        do jj=1,3
          k=jj+1
          os(1,k,2:4,2:4,0:1)=(os(1,1,2:4,2:4,0:1)*rbp(jj)+os(1,1,2:4,2:4,1:2)*rpw(jj))
          os(1,k,2:4,k,0:1)=os(1,k,2:4,k,0:1)+hooppq*os(1,1,2:4,1,1:2)
          os(1,k,k,2:4,0:1)=os(1,k,k,2:4,0:1)+hooppq*os(1,1,1,2:4,1:2)
          os(k,1,2:4,2:4,0)=(os(1,1,2:4,2:4,0)*rap(jj)+os(1,1,2:4,2:4,1)*rpw(jj))
          os(k,1,2:4,k,0)=os(k,1,2:4,k,0)+hooppq*os(1,1,2:4,1,1)
          os(k,1,k,2:4,0)=os(k,1,k,2:4,0)+hooppq*os(1,1,1,2:4,1)
          os(2:4,2:4,k,1,0)=(os(2:4,2:4,1,1,0)*rcq(jj)+os(2:4,2:4,1,1,1)*rqw(jj))
          os(k,2:4,k,1,0)=os(k,2:4,k,1,0)+hooppq*os(1,2:4,1,1,1)
          os(2:4,k,k,1,0)=os(2:4,k,k,1,0)+hooppq*os(2:4,1,1,1,1)
          os(2:4,2:4,1,k,0)=(os(2:4,2:4,1,1,0)*rdq(jj)+os(2:4,2:4,1,1,1)*rqw(jj))
          os(k,2:4,1,k,0)=os(k,2:4,1,k,0)+hooppq*os(1,2:4,1,1,1)
          os(2:4,k,1,k,0)=os(2:4,k,1,k,0)+hooppq*os(2:4,1,1,1,1)
        end do
        ! Final recursion: The (pp|pp) part
        do jj=1,3
          k=jj+1
          os(k,2:4,2:4,2:4,0)=(os(1,2:4,2:4,2:4,0)*rap(jj)+os(1,2:4,2:4,2:4,1)*rpw(jj))
          os(k,k,2:4,2:4,0)=os(k,k,2:4,2:4,0)+(os(1,1,2:4,2:4,0)-qoppq*os(1,1,2:4,2:4,1))*hpm1(zzab,ab)
          os(k,2:4,k,2:4,0)=os(k,2:4,k,2:4,0)+hooppq*os(1,2:4,1,2:4,1)
          os(k,2:4,2:4,k,0)=os(k,2:4,2:4,k,0)+hooppq*os(1,2:4,2:4,1,1)
        end do
        ggkl=djkl(zc,sh3)*djkl(zd,sh4)
        osc(:,:,1,1)=osc(:,:,1,1)+os(:,:,1,1,0)
        osc(:,:,1,2:4)=osc(:,:,1,2:4)+os(:,:,1,2:4,0)*djkl(zd,sh4)
        osc(:,:,2:4,1)=osc(:,:,2:4,1)+os(:,:,2:4,1,0)*djkl(zc,sh3)
        osc(:,:,2:4,2:4)=osc(:,:,2:4,2:4)+os(:,:,2:4,2:4,0)*ggkl
      end do
    end do
    ggij=djkl(za,sh1)*djkl(zb,sh2)
    ellll(1,1,:,:)=ellll(1,1,:,:)+osc(1,1,:,:)
    ellll(1,2:4,:,:)=ellll(1,2:4,:,:)+osc(1,2:4,:,:)*djkl(zb,sh2)
    ellll(2:4,1,:,:)=ellll(2:4,1,:,:)+osc(2:4,1,:,:)*djkl(za,sh1)
    ellll(2:4,2:4,:,:)=ellll(2:4,2:4,:,:)+osc(2:4,2:4,:,:)*ggij
  end do
end do

end subroutine i2e_llll



end module s_d_l_terms

