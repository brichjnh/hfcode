module one_electron_terms
use nrtype
implicit none

contains

pure subroutine s_s_1e_terms(sh1,sh2,s,t,v)
use nrtype ; use molprops ; use boys
implicit none

integer(i4b), intent(in) :: sh1,sh2
real(dp), intent(out) :: s, t, v

integer(i4b) :: nuc, ab, za, zb, zza, zzab
real(dp) :: sterm, tterm, vterm
real(dp) :: fac1, fac2, fac3, fac5
real(dp) :: Cx(3),zc, rpc2, qr2

sterm=0._dp
tterm=0._dp
vterm=0._dp
ab=ind2(sh1,sh2)
do za =1,nzet(sh1)
  zza=maxprim*(za-1)
  do zb = 1, nzet(sh2)
    zzab=zza+zb
    fac1=pi32/pab(zzab,ab)**1.5_dp*ddab(zzab,ab)
    qr2=qab(zzab,ab)*rab2(ab)
    fac2=exp(-qr2)
    sterm=sterm+fac1*fac2
    fac3=3._dp-2._dp*qr2
    tterm=tterm+qab(zzab,ab)*fac3*fac1*fac2

    fac3=-2._dp*PI*pm1(zzab,ab)*fac2*ddab(zzab,ab)
    do nuc =1,natom
       Cx=atcoords(nuc,:)
       zc=nuccharges(nuc)
       rpc2=sum((Cx-px(:,zzab,ab))**2)
       call boys_0(pab(zzab,ab)*rpc2,fac5)
       vterm=vterm+fac3*zc*fac5
    end do
  end do
end do
s=sterm
t=tterm
v=vterm

end subroutine s_s_1e_terms


pure subroutine s_p_1e_terms(sh1,sh2,s,t,v)
use nrtype ; use molprops ; use boys
implicit none

! Here the p function is assumed to be on the second atom
! For overlap, we use the simple recurrence in 1 dimension: S_sp=S_ss * R_ap

integer(i4b), intent(in) :: sh1, sh2
real(dp), intent(out) :: s(3), t(3), v(3)

integer(i4b) :: nuc, ab, za, zb, zza, zzab
real(dp) :: sss, ssp(3),tss, rbp(3), vssm(0:1)
real(dp) :: fac1, fac2, fac3, bfuncs(0:1)
real(dp) :: Cx(3), zc, rcp(3), rcp2, efac, qr2,vterm(3)

ab=ind2(sh1,sh2) 
s=0._dp
t=0._dp
v=0._dp
do za =1,nzet(sh1)
  zza=(za-1)*maxprim
  do zb = 1, nzet(sh2)
    zzab=zza+zb
    ! Start with the s-s overlap
    fac1=pi32/pab(zzab,ab)**1.5_dp        
    qr2=qab(zzab,ab)*rab2(ab)
    efac=exp(-qr2)*ddab(zzab,ab)
    sss=efac*fac1

    ! Then use recursion to get the s-p overlaps from the ss overlaps
    rbp=Px(:,zzab,ab)-xsh(:,sh2)
    ssp(:)=sss*rbp(:)
    s=s+ssp

    ! Then use recursion to get s-p kinetic integrals from s-s and s-p overlaps and s-s kinetic
    fac2=qab(zzab,ab)*(3._dp-2._dp*qr2)
    tss=fac2*sss
    t=t+rbp*tss+2._dp*qab(zzab,ab)*ssp

    ! Then use recursion to get s-p nuclear Coulomb integrals from s-s m=1 and m=0 integrals
    fac3=-2._dp*sqrt(pab(zzab,ab)/pi)*sss
    vterm=0._dp
    do nuc=1,natom
      Cx=atcoords(nuc,:)
      zc=nuccharges(nuc)
      rcp=Px(:,zzab,ab)-Cx
      rcp2=sum(rcp**2)
      call boys_1(pab(zzab,ab)*rcp2,bfuncs)
      vssm=fac3*bfuncs
      vterm=vterm+zc*(rbp*vssm(0)-rcp*vssm(1))
    end do
    v=v+vterm
  end do
end do

end subroutine s_p_1e_terms

pure subroutine s_l_1e_terms(sh1,sh2,s,t,v)
use nrtype ; use molprops ; use boys
implicit none

! Here the l shell is assumed to be on the second atom
! (s|s) (s|T|s) and (s|V|s) integrals output first, then px, py, pz
! For overlap, we use the simple recurrence in 1 dimension: S_sp=S_ss * R_ap

integer(i4b), intent(in) :: sh1, sh2
real(dp), intent(out) :: s(4), t(4), v(4)

integer(i4b) :: nuc, ab, za, zb, zza, zzab
real(dp) :: sss, ssp(3),tss, rbp(3), vssm(0:1)
real(dp) :: fac1, fac2, fac3, bfuncs(0:1), gg
real(dp) :: Cx(3),zc, rcp(3), rcp2, efac, qr2,vterm(4)

ab=ind2(sh1,sh2) 
s=0._dp
t=0._dp
v=0._dp
do za =1,nzet(sh1)
  zza=(za-1)*maxprim
  do zb = 1, nzet(sh2)
    zzab=zza+zb
    gg=djkl(zb,sh2)
    ! Start with the s-s overlap
    fac1=pi32/pab(zzab,ab)**1.5_dp        
    qr2=qab(zzab,ab)*rab2(ab)
    efac=exp(-qr2)*ddab(zzab,ab)
    sss=efac*fac1
    s(1)=s(1)+sss

    ! Then use recursion to get the s-p overlaps from the ss overlaps
    rbp=Px(:,zzab,ab)-xsh(:,sh2)
    ssp(:)=sss*rbp(:)*gg
    s(2:4)=s(2:4)+ssp

    ! Then use recursion to get s-p kinetic integrals from s-s and s-p overlaps and s-s kinetic
    fac2=qab(zzab,ab)*(3._dp-2._dp*qr2)
    tss=fac2*sss
    t(1)=t(1)+tss
    t(2:4)=t(2:4)+(gg*rbp*tss+2._dp*qab(zzab,ab)*ssp)

    ! Then use recursion to get s-p nuclear Coulomb integrals from s-s m=1 and m=0 integrals
    fac3=-2._dp*sqrt(pab(zzab,ab)/pi)*sss
    vterm=0._dp
    do nuc=1,natom
      Cx=atcoords(nuc,:)
      zc=nuccharges(nuc)
      rcp=Px(:,zzab,ab)-Cx
      rcp2=sum(rcp**2)
      call boys_1(pab(zzab,ab)*rcp2,bfuncs)
      vssm=fac3*bfuncs
      vterm(1)=vterm(1)+zc*vssm(0)
      vterm(2:4)=vterm(2:4)+zc*gg*(rbp*vssm(0)-rcp*vssm(1))
    end do
    v=v+vterm
  end do
end do

end subroutine s_l_1e_terms


pure subroutine p_p_1e_terms(sh1,sh2,s,t,v)
use nrtype ; use molprops ; use boys
implicit none

integer(i4b), intent(in) :: sh1, sh2
real(dp), intent(out) :: s(3,3), t(3,3), v(3,3)

integer(i4b) :: nuc, ab, za, zb, zza, zzab, ii, jj
real(dp) :: fac1, fac2, fac3
real(dp) :: Cx(3),zc, rcp(3), rcp2, efac, qr2
real(dp) :: sss, ssp(3),tss,tsp(3), spp(3,3), vsp0(3),vsp1(3),vpp(3,3)
real(dp) :: rbp(3), rap(3), vssm(0:2)

s=0._dp
t=0._dp
v=0._dp
ab=ind2(sh1,sh2)
do za =1,nzet(sh1)
  zza=(za-1)*maxprim
  do zb = 1, nzet(sh2)
    zzab=zza+zb
    fac1=pi32/pab(zzab,ab)**1.5_dp      
    rap=px(:,zzab,ab)-xsh(:,sh1)
    rbp=px(:,zzab,ab)-xsh(:,sh2)
    qr2=qab(zzab,ab)*rab2(ab)
    efac=exp(-qr2)*ddab(zzab,ab)

    ! Use recursion to get the p-p overlaps from the ss and sp overlaps
    ! First recurse over ket, then bra
    sss=efac*fac1
    ssp=sss*rbp
    fac1=hpm1(zzab,ab)*sss
    do ii=1,3
      do jj=1,3
        spp(ii,jj)=ssp(jj)*rap(ii)
      end do
      spp(ii,ii)=spp(ii,ii)+fac1
    end do
    s=s+spp

    ! Then use recursion to get p-p kinetic integrals from s-s and s-p overlaps and s-s kinetic
    ! First recurse over ket, then bra
    fac2=qab(zzab,ab)*(3._dp-2._dp*qr2)
    tss=fac2*sss
    tsp=rbp*tss+2._dp*qab(zzab,ab)*ssp
    fac3=fac2*fac1
    do ii=1,3
      do jj=1,3
        t(ii,jj)=t(ii,jj)+tsp(jj)*rap(ii)
      end do
      t(ii,ii)=t(ii,ii)+fac3
    end do
    t=t+2._dp*qab(zzab,ab)*spp

    ! Finally use recursion to get p-p nuclear Coulomb integrals
    vpp=0._dp
    do nuc =1,natom
      Cx=atcoords(nuc,:)
      zc=nuccharges(nuc)
      rcp=px(:,zzab,ab)-Cx
      rcp2=sum(rcp**2)
      call boys_2(pab(zzab,ab)*rcp2,vssm)
      vsp0=rbp*vssm(0)-rcp*vssm(1)
      vsp1=rbp*vssm(1)-rcp*vssm(2)
      do ii=1,3
        do jj=1,3
          vpp(ii,jj)=vpp(ii,jj)+zc*(rap(ii)*vsp0(jj)-rcp(ii)*vsp1(jj))
        end do
        vpp(ii,ii)=vpp(ii,ii)+zc*hpm1(zzab,ab)*(vssm(0)-vssm(1))
      end do
    end do
    fac3=-2._dp*PI*pm1(zzab,ab)*efac
    v=v+vpp*fac3
   end do
end do

end subroutine p_p_1e_terms


pure subroutine p_l_1e_terms(sh1,sh2,so,to,vo)
use nrtype ; use molprops ; use boys
implicit none

integer(i4b), intent(in) :: sh1, sh2
real(dp), intent(out) :: so(3,4), to(3,4), vo(3,4)

integer(i4b) :: nuc, ab, za, zb, zza, zzab, ii, jj
real(dp) :: fac1, fac2, fac3
real(dp) :: Cx(3),zc, rcp(3), rcp2, efac, qr2
real(dp) :: sss, ssp(3),tss,tsp(3), spp(3,3), vsp0(3),vsp1(3),vmat(4,4)
real(dp) :: tpp(3,3), sps(3), tps(3)
real(dp) :: rbp(3), rap(3), vssm(0:2), s(4,4), t(4,4), v(4,4)

s=0._dp
t=0._dp
v=0._dp
ab=ind2(sh1,sh2)
do za =1,nzet(sh1)
  zza=(za-1)*maxprim
  do zb = 1, nzet(sh2)
    zzab=zza+zb
    fac1=pi32/pab(zzab,ab)**1.5_dp      
    rap=px(:,zzab,ab)-xsh(:,sh1)
    rbp=px(:,zzab,ab)-xsh(:,sh2)
    qr2=qab(zzab,ab)*rab2(ab)
    efac=exp(-qr2)*ddab(zzab,ab)

    ! Use recursion to get the p-p overlaps from the ss and sp overlaps
    ! First recurse over ket, then bra
    sss=efac*fac1
    s(1,1)=s(1,1)+sss
    sps=sss*rap
    ssp=sss*rbp
    s(2:4,1)=s(2:4,1)+sps
    s(1,2:4)=s(1,2:4)+ssp*djkl(zb,sh2)
    fac1=hpm1(zzab,ab)*sss
    do ii=1,3
      do jj=1,3
        spp(ii,jj)=ssp(jj)*rap(ii)
      end do
      spp(ii,ii)=spp(ii,ii)+fac1
    end do
    s(2:4,2:4)=s(2:4,2:4)+spp*djkl(zb,sh2)

    ! Then use recursion to get p-p kinetic integrals from s-s and s-p overlaps and s-s kinetic
    ! First recurse over ket, then bra
    fac2=qab(zzab,ab)*(3._dp-2._dp*qr2)
    tss=fac2*sss
    t(1,1)=t(1,1)+tss
    tps=rap*tss+2._dp*qab(zzab,ab)*sps
    t(2:4,1)=t(2:4,1)+tps
    tsp=rbp*tss+2._dp*qab(zzab,ab)*ssp
    t(1,2:4)=t(1,2:4)+tsp*djkl(zb,sh2)
    fac3=fac2*fac1
    do ii=1,3
      do jj=1,3
        tpp(ii,jj)=tsp(jj)*rap(ii)
      end do
      tpp(ii,ii)=tpp(ii,ii)+fac3
    end do
    tpp=tpp+2._dp*qab(zzab,ab)*spp
    t(2:4,2:4)=t(2:4,2:4)+tpp*djkl(zb,sh2)

    ! Finally use recursion to get p-p nuclear Coulomb integrals
    vmat=0._dp
    do nuc =1,natom
      Cx=atcoords(nuc,:)
      zc=nuccharges(nuc)
      rcp=px(:,zzab,ab)-Cx
      rcp2=sum(rcp**2)
      call boys_2(pab(zzab,ab)*rcp2,vssm)
      vssm=vssm*zc
      vmat(1,1)=vmat(1,1)+vssm(0)
      vsp0=rbp*vssm(0)-rcp*vssm(1)
      vsp1=rbp*vssm(1)-rcp*vssm(2)
      vmat(1,2:4)=vmat(1,2:4)+vsp0*djkl(zb,sh2)
      vmat(2:4,1)=vmat(2:4,1)+(rap*vssm(0)-rcp*vssm(1))
      do ii=1,3
        do jj=1,3
          vmat(1+ii,1+jj)=vmat(1+ii,1+jj)+(rap(ii)*vsp0(jj)-rcp(ii)*vsp1(jj))*djkl(zb,sh2)
        end do
        vmat(1+ii,1+ii)=vmat(1+ii,1+ii)+hpm1(zzab,ab)*(vssm(0)-vssm(1))*djkl(zb,sh2)
      end do
    end do
    fac3=-2._dp*PI*pm1(zzab,ab)*efac
  v=v+vmat*fac3
  end do
end do
so=s(2:4,:)
to=t(2:4,:)
vo=v(2:4,:)

end subroutine p_l_1e_terms



pure subroutine l_l_1e_terms(sh1,sh2,s,t,v)
use nrtype ; use molprops ; use boys
implicit none

integer(i4b), intent(in) :: sh1, sh2
real(dp), intent(out) :: s(4,4), t(4,4), v(4,4)

integer(i4b) :: nuc, ab, za, zb, zza, zzab, ii, jj
real(dp) :: fac1, fac2, fac3
real(dp) :: Cx(3),zc, rcp(3), rcp2, efac, qr2
real(dp) :: sss, ssp(3),tss,tsp(3), spp(3,3), vsp0(3),vsp1(3),vmat(4,4)
real(dp) :: tpp(3,3), sps(3), tps(3)
real(dp) :: rbp(3), rap(3), vssm(0:2), gg

s=0._dp
t=0._dp
v=0._dp
ab=ind2(sh1,sh2)
do za =1,nzet(sh1)
  zza=(za-1)*maxprim
  do zb = 1, nzet(sh2)
    zzab=zza+zb
    gg=djkl(za,sh1)*djkl(zb,sh2)
    fac1=pi32/pab(zzab,ab)**1.5_dp      
    rap=px(:,zzab,ab)-xsh(:,sh1)
    rbp=px(:,zzab,ab)-xsh(:,sh2)
    qr2=qab(zzab,ab)*rab2(ab)
    efac=exp(-qr2)*ddab(zzab,ab)

    ! Use recursion to get the p-p overlaps from the ss and sp overlaps
    ! First recurse over ket, then bra
    sss=efac*fac1
    s(1,1)=s(1,1)+sss
    sps=sss*rap
    ssp=sss*rbp
    s(2:4,1)=s(2:4,1)+sps*djkl(za,sh1)
    s(1,2:4)=s(1,2:4)+ssp*djkl(zb,sh2)
    fac1=hpm1(zzab,ab)*sss
    do ii=1,3
      do jj=1,3
        spp(ii,jj)=ssp(jj)*rap(ii)
      end do
      spp(ii,ii)=spp(ii,ii)+fac1
    end do
    s(2:4,2:4)=s(2:4,2:4)+spp*gg

    ! Then use recursion to get p-p kinetic integrals from s-s and s-p overlaps and s-s kinetic
    ! First recurse over ket, then bra
    fac2=qab(zzab,ab)*(3._dp-2._dp*qr2)
    tss=fac2*sss
    t(1,1)=t(1,1)+tss
    tps=rap*tss+2._dp*qab(zzab,ab)*sps
    t(2:4,1)=t(2:4,1)+tps*djkl(za,sh1)
    tsp=rbp*tss+2._dp*qab(zzab,ab)*ssp
    t(1,2:4)=t(1,2:4)+tsp*djkl(zb,sh2)
    fac3=fac2*fac1
    do ii=1,3
      do jj=1,3
        tpp(ii,jj)=tsp(jj)*rap(ii)
      end do
      tpp(ii,ii)=tpp(ii,ii)+fac3
    end do
    tpp=tpp+2._dp*qab(zzab,ab)*spp
    t(2:4,2:4)=t(2:4,2:4)+tpp*gg

    ! Finally use recursion to get p-p nuclear Coulomb integrals
    vmat=0._dp
    do nuc =1,natom
      Cx=atcoords(nuc,:)
      zc=nuccharges(nuc)
      rcp=px(:,zzab,ab)-Cx
      rcp2=sum(rcp**2)
      call boys_2(pab(zzab,ab)*rcp2,vssm)
      vssm=vssm*zc
      vmat(1,1)=vmat(1,1)+vssm(0)
      vsp0=rbp*vssm(0)-rcp*vssm(1)
      vsp1=rbp*vssm(1)-rcp*vssm(2)
      vmat(1,2:4)=vmat(1,2:4)+vsp0*djkl(zb,sh2)
      vmat(2:4,1)=vmat(2:4,1)+(rap*vssm(0)-rcp*vssm(1))*djkl(za,sh1)
      do ii=1,3
        do jj=1,3
          vmat(1+ii,1+jj)=vmat(1+ii,1+jj)+(rap(ii)*vsp0(jj)-rcp(ii)*vsp1(jj))*gg
        end do
        vmat(1+ii,1+ii)=vmat(1+ii,1+ii)+hpm1(zzab,ab)*(vssm(0)-vssm(1))*gg
      end do
    end do
    fac3=-2._dp*PI*pm1(zzab,ab)*efac
  v=v+vmat*fac3
  end do
end do

end subroutine l_l_1e_terms



pure subroutine s_d_1e_terms(sh1,sh2,s,t,v)
use nrtype ; use molprops ; use boys
implicit none

! This is the first place where I use d functions. I store these in 1-d arrays, Cartesian & sph harm
! This is the first routine where I use the generalized indices and 'minus 1' Table:
! !   0    1    2    3    4    5    6    7    8    9   10
!     /    s    x    y    z   xx   xy   xz   yy   yz   zz

! Conventional order, Cartesian: xx, xy, xz, yy, yz, zz
integer(i4b), intent(in) :: sh1, sh2
real(dp), intent(out) :: s(5), t(5), v(5)

integer(i4b) :: nuc, ab, za, zb, zza, zzab, ii, jj, j, jm1, jp1, m
real(dp) :: fac1, fac2, fac3
real(dp) :: Cx(3),zc, rcp(3), rcp2, efac, qr2, q, zaop
real(dp) :: smat(10), tmat(10),vmat(10,0:2), csmat(6),ctmat(6),cvmat(6)
real(dp) :: rbp(3), bfac(0:2)

csmat=0._dp
ctmat=0._dp
cvmat=0._dp

ab=ind2(sh1,sh2) 
do za =1,nzet(sh1)
  zza=(za-1)*maxprim
  do zb = 1, nzet(sh2)
    zzab=zza+zb
    ! Start as usual with the S-S overlap
    fac1=pi32/pab(zzab,ab)**1.5_dp      
    qr2=qab(zzab,ab)*rab2(ab)
    efac=exp(-qr2)*ddab(zzab,ab)
    smat(1)=efac*fac1

    ! Use recursion to get the s-d overlaps from the ss and sp overlaps
    rbp=px(:,zzab,ab)-xsh(:,sh2)
    smat(2:4)=smat(1)*rbp(:)
    fac2=hpm1(zzab,ab)*smat(1)
    do jp1=5,10
      do jj=1,3 
        j=minus1s(jj,jp1)
        if (j.eq.0) cycle
        smat(jp1)=smat(j)*rbp(jj)
        jm1=minus1s(jj,j)
        if (jm1.eq.0) cycle
        smat(jp1)=smat(jp1)+fac2
      end do
    end do
    csmat=csmat+smat(5:10)

    ! Now use recursion for the s-d kinetic integrals
    fac3=qab(zzab,ab)*(3._dp-2._dp*qr2)
    tmat(1)=fac3*smat(1)
    tmat(2:4)=rbp(:)*tmat(1)+2._dp*qab(zzab,ab)*smat(2:4)
    zaop=zeta(za,sh1)*pm1(zzab,ab)
    q=zaop*zeta(zb,sh2)
!    do jp1=5,10
!      tmat(jp1)=2._dp*q*smat(jp1)
!      do jj=1,3 
!        j=minus1s(jj,jp1)
!        jm1=minus1s(jj,j)
!        if (j.eq.0) cycle
!        tmat(jp1)=tmat(jp1)+rbp(jj)*tmat(j)
!        if (jm1.eq.0) cycle
!        tmat(jp1)=tmat(jp1)+tmat(1)/(2._dp*pab(zzab,ab))-zaop*smat(1)
!      end do
!    end do
    tmat(5)=rbp(1)*tmat(2)+hpm1(zzab,ab)*tmat(1)+2._dp*q*(smat(5)-smat(1)/(2._dp*zeta(zb,sh2)))
    tmat(6)=rbp(2)*tmat(2)+2._dp*q*smat(6)
    tmat(7)=rbp(3)*tmat(2)+2._dp*q*smat(7)
    tmat(8)=rbp(2)*tmat(3)+hpm1(zzab,ab)*tmat(1)+2._dp*q*(smat(8)-smat(1)/(2._dp*zeta(zb,sh2)))
    tmat(9)=rbp(3)*tmat(3)+2._dp*q*smat(9)
    tmat(10)=rbp(3)*tmat(4)+hpm1(zzab,ab)*tmat(1)+2._dp*q*(smat(10)-smat(1)/(2._dp*zeta(zb,sh2)))
    ctmat=ctmat+tmat(5:10)

    ! Finally use recursion to get s-d nuclear Coulomb integrals

    fac1=2._dp*sqrt(pab(zzab,ab)/PI)
    do nuc =1,natom
      vmat=0._dp
      Cx=atcoords(nuc,:)
      zc=nuccharges(nuc)
      rcp=px(:,zzab,ab)-Cx
      rcp2=sum(rcp**2)
      call boys_2(pab(zzab,ab)*rcp2,bfac)
      vmat(1,:)=-zc*fac1*smat(1)*bfac(:)
      do m=0,1
        vmat(2:4,m)=rbp(:)*vmat(1,m)-rcp(:)*vmat(1,m+1)
      end do
      fac2=(vmat(1,0)-vmat(1,1))*hpm1(zzab,ab)
      do jp1=5,10
        do jj=1,3 
          j=minus1s(jj,jp1)
          if (j.eq.0) cycle
          vmat(jp1,0)=rbp(jj)*vmat(j,0)-rcp(jj)*vmat(j,1)
          jm1=minus1s(jj,j)
          if (jm1.eq.0) cycle
          vmat(jp1,0)=vmat(jp1,0)+fac2
        end do
      end do
      cvmat=cvmat+vmat(5:10,0)
    end do

  end do
end do

! now contract to form spherical harmonics. Rem: Cartesians:
!     0    1    2    3    4    5    6    7    8    9   10
!     /    s    x    y    z   xx   xy   xz   yy   yz   zz
! Spherical harmonics: 1=z2,2=xz,3=yz,4=xy,5=x2-y2
do ii=1,5
  s(ii)=sum(csmat(:)*matharmond(ii,:)) 
  t(ii)=sum(ctmat(:)*matharmond(ii,:)) 
  v(ii)=sum(cvmat(:)*matharmond(ii,:)) 
end do

end subroutine s_d_1e_terms


pure subroutine p_d_1e_terms(sh1,sh2,s,t,v)
use nrtype ; use molprops ; use boys
implicit none

! This is the second place where I use d functions. I store these in 1-d arrays, Cartesian & sph harm
! This is the first routine where I use the generalized indices and 'minus 1' Table:
! !   0    1    2    3    4    5    6    7    8    9   10
!     /    s    x    y    z   xx   xy   xz   yy   yz   zz
! Note that this routine still needs checking

! Conventional order, Cartesian: xx, xy, xz, yy, yz, zz
integer(i4b), intent(in) :: sh1, sh2
real(dp), intent(out) :: s(3,5), t(3,5), v(3,5)

integer(i4b) :: nuc, ab, za, zb, zza, zzab, ii, jj, jp1, j, jm1, ip1,i,im1, m
real(dp) :: fac1, fac2, fac3, zaop, q
real(dp) :: Cx(3),zc, rcp(3), rcp2, efac, qr2
real(dp) :: smat(4,10), tmat(4,10),vmat(4,10,0:3), csmat(3,6),ctmat(3,6),cvmat(3,6)
real(dp) :: rbp(3), bfac(0:3), rap(3)

csmat=0._dp
ctmat=0._dp
cvmat=0._dp

ab=ind2(sh1,sh2) 
do za =1,nzet(sh1)
  zza=(za-1)*maxprim
  do zb = 1, nzet(sh2)
    zzab=zza+zb
    ! Start as usual with the S-S overlap
    fac1=pi32/pab(zzab,ab)**1.5_dp
    qr2=qab(zzab,ab)*rab2(ab)
    efac=exp(-qr2)*ddab(zzab,ab)
    smat(1,1)=efac*fac1

    ! Use recursion to get the s-d overlaps from the ss and sp overlaps
    rbp=px(:,zzab,ab)-xsh(:,sh2)
    smat(1,2:4)=smat(1,1)*rbp(:)
    fac2=hpm1(zzab,ab)*smat(1,1)
    do jp1=5,10
      do jj=1,3 
        j=minus1s(jj,jp1)
        if (j.eq.0) cycle
        smat(1,jp1)=smat(1,j)*rbp(jj)
        jm1=minus1s(jj,j)
        if (jm1.eq.0) cycle
        smat(1,jp1)=smat(1,jp1)+fac2
      end do
    end do
    rap=px(:,zzab,ab)-xsh(:,sh1)
    do jj=1,3
      smat(1+jj,:)=rap(jj)*smat(1,:)
      do i=2,10
        im1=minus1s(jj,i)
        if (im1.eq.0) cycle
        smat(1+jj,i)=smat(1+jj,i)+smat(1,im1)*hpm1(zzab,ab)*indexps(jj,i)
      end do
    end do
    csmat=csmat+smat(2:4,5:10)

    ! Now use recursion for the s-d kinetic integrals
    fac3=qab(zzab,ab)*(3._dp-2._dp*qr2)
    tmat(1,1)=fac3*smat(1,1)
    tmat(1,2:4)=rbp(:)*tmat(1,1)+2._dp*qab(zzab,ab)*smat(1,2:4)
    zaop=zeta(za,sh1)*pm1(zzab,ab)
    q=zaop*zeta(zb,sh2)
    tmat(1,5)=rbp(1)*tmat(1,2)+hpm1(zzab,ab)*tmat(1,1)+2._dp*q*(smat(1,5)-smat(1,1)/(2._dp*zeta(zb,sh2)))
    tmat(1,6)=rbp(2)*tmat(1,2)+2._dp*q*smat(1,6)
    tmat(1,7)=rbp(3)*tmat(1,2)+2._dp*q*smat(1,7)
    tmat(1,8)=rbp(2)*tmat(1,3)+hpm1(zzab,ab)*tmat(1,1)+2._dp*q*(smat(1,8)-smat(1,1)/(2._dp*zeta(zb,sh2)))
    tmat(1,9)=rbp(3)*tmat(1,3)+2._dp*q*smat(1,9)
    tmat(1,10)=rbp(3)*tmat(1,4)+hpm1(zzab,ab)*tmat(1,1)+2._dp*q*(smat(1,10)-smat(1,1)/(2._dp*zeta(zb,sh2)))
    do j=2,10
      tmat(2:4,j)=rap(:)*tmat(1,j)+2._dp*q*smat(2:4,j)
      do jj=1,3
        ip1=jj+1
        jm1=minus1s(jj,j)
        if (jm1.eq.0) cycle
        tmat(ip1,j)=tmat(ip1,j)+hpm1(zzab,ab)*tmat(1,jm1)*indexps(jj,j)
      end do
    end do
    ctmat=ctmat+tmat(2:4,5:10)

    ! Finally use recursion to get p-d nuclear Coulomb integrals

    fac1=2._dp*sqrt(pab(zzab,ab)/PI)
    do nuc =1,natom
      vmat=0._dp
      Cx=atcoords(nuc,:)
      zc=nuccharges(nuc)
      rcp=px(:,zzab,ab)-Cx
      rcp2=sum(rcp**2)
      call boys_n(3,pab(zzab,ab)*rcp2,bfac)
      vmat(1,1,:)=-zc*fac1*smat(1,1)*bfac(:)
      do m=0,2
        vmat(1,2:4,m)=rbp(:)*vmat(1,1,m)-rcp(:)*vmat(1,1,m+1)
      end do
      do jp1=5,10
        do jj=1,3 
          j=minus1s(jj,jp1)
          if (j.eq.0) cycle
          vmat(1,jp1,0:1)=rbp(jj)*vmat(1,j,0:1)-rcp(jj)*vmat(1,j,1:2)
          jm1=minus1s(jj,j)
          if (jm1.eq.0) cycle
          vmat(1,jp1,0:1)=vmat(1,jp1,0:1)+hpm1(zzab,ab)*(vmat(1,1,0:1)-vmat(1,1,1:2))
        end do
      end do
      do jj=1,3
        ip1=jj+1
        do j=2,10
          vmat(ip1,j,0)=rap(jj)*vmat(1,j,0)-rcp(jj)*vmat(1,j,1)
          jm1=minus1s(jj,j)
          if (jm1.eq.0) cycle
          vmat(ip1,j,0)=vmat(ip1,j,0)+hpm1(zzab,ab)*indexps(jj,j)*(vmat(1,jm1,0)-vmat(1,jm1,1))
        end do
      end do
      cvmat=cvmat+vmat(2:4,5:10,0)
    end do

  end do
end do

! now contract to form spherical harmonics. Rem: Cartesians:
!     0    1    2    3    4    5    6    7    8    9   10
!     /    s    x    y    z   xx   xy   xz   yy   yz   zz
! Spherical harmonics: 1=z2,2=xz,3=yz,4=xy,5=x2-y2

do ii=1,3
  do jj=1,5
    s(ii,jj)=sum(csmat(ii,:)*matharmond(jj,:))
    t(ii,jj)=sum(ctmat(ii,:)*matharmond(jj,:))
    v(ii,jj)=sum(cvmat(ii,:)*matharmond(jj,:))
  end do
end do

return
end subroutine p_d_1e_terms


pure subroutine d_d_1e_terms(sh1,sh2,s,t,v)
use nrtype ; use molprops ; use boys
implicit none

! This is the second place where I use d functions. I store these in 1-d arrays, Cartesian & sph harm
! This is the first routine where I use the generalized indices and 'minus 1' Table:
! !   0    1    2    3    4    5    6    7    8    9   10
!     /    s    x    y    z   xx   xy   xz   yy   yz   zz
! Note that this routine still needs checking

! Conventional order, Cartesian: xx, xy, xz, yy, yz, zz
integer(i4b), intent(in) :: sh1, sh2
real(dp), intent(out) :: s(5,5), t(5,5), v(5,5)

integer(i4b) :: nuc, ab, za, zb, zza, zzab, ii, jj, jp1, j, jm1, ip1, i, im1, m, jrec
real(dp) :: fac1
real(dp) :: Cx(3),zc, rcp(3), rcp2, efac, qr2, q
real(dp) :: smat(10,10), tmat(10,10),vmat(10,10,0:4), csmat(6,6),ctmat(6,6),cvmat(6,6)
real(dp) :: rbp(3), bfac(0:4), rap(3), smatl(6,5), tmatl(6,5), vmatl(6,5)

csmat=0._dp
ctmat=0._dp
cvmat=0._dp

ab=ind2(sh1,sh2) 
do za =1,nzet(sh1)
  zza=(za-1)*maxprim
  do zb = 1, nzet(sh2)
    zzab=zza+zb

    rap=px(:,zzab,ab)-xsh(:,sh1)
    rbp=px(:,zzab,ab)-xsh(:,sh2)
    q=qab(zzab,ab)
    qr2=q*rab2(ab)
    fac1=(pi*pm1(zzab,ab))**(1.5_dp)
    efac=exp(-qr2)*ddab(zzab,ab)

    ! Overlap integrals
    smat(1,1)=fac1*efac
    smat(1,2:4)=rbp(:)*smat(1,1)
    smat(2:4,1)=rap(:)*smat(1,1)
    do jp1=5,10
      do jj=1,3
        j=minus1s(jj,jp1)
        if (j.eq.0) cycle
        jrec=jj
        exit
      end do
      smat(1,jp1)=rbp(jrec)*smat(1,j)
      jm1=minus1s(jrec,j)
      if (jm1.eq.0) cycle
      smat(1,jp1)=smat(1,jp1)+hpm1(zzab,ab)*smat(1,jm1)
    end do
    do jj=1,3
      ip1=jj+1
      smat(ip1,:)=rap(jj)*smat(1,:)
      do j=2,10
        jm1=minus1s(jj,j)
        if (jm1.eq.0) cycle
        smat(ip1,j)=smat(ip1,j)+hpm1(zzab,ab)*smat(1,jm1)*indexps(jj,j)
      end do
    end do
    do ip1=5,10
      do jj=1,3
        i=minus1s(jj,ip1)
        if (i.eq.0) cycle
        jrec=jj
        exit
      end do
      smat(ip1,:)=rap(jrec)*smat(i,:)
      im1=minus1s(jrec,i)
      if (im1.ne.0) then
        smat(ip1,:)=smat(ip1,:)+hpm1(zzab,ab)*smat(im1,:)
      end if
      do j=2,10
        jm1=minus1s(jrec,j)
        if (jm1.eq.0) cycle
        smat(ip1,j)=smat(ip1,j)+hpm1(zzab,ab)*indexps(jrec,j)*smat(i,jm1)
      end do
    end do
    csmat=csmat+smat(5:10,5:10)

    fac1=2._dp*sqrt(pab(zzab,ab)/PI)
    do nuc =1,natom
      vmat=0._dp
      Cx=atcoords(nuc,:)
      zc=nuccharges(nuc)
      rcp=px(:,zzab,ab)-Cx
      rcp2=sum(rcp**2)
      call boys_n(4,pab(zzab,ab)*rcp2,bfac)
      ! (s|V|s)
      vmat(1,1,:)=-zc*fac1*smat(1,1)*bfac(:)
      ! (s|V|p)
      do m=0,3
        vmat(1,2:4,m)=rbp(:)*vmat(1,1,m)-rcp(:)*vmat(1,1,m+1)
      end do
      ! (s|V|d)
      do jp1=5,10
        do jj=1,3
          j=minus1s(jj,jp1)
          if (j.eq.0) cycle
          jrec=jj
          exit
        end do
        vmat(1,jp1,0:2)=rbp(jrec)*vmat(1,j,0:2)-rcp(jrec)*vmat(1,j,1:3)
        jm1=minus1s(jrec,j)
        if (jm1.eq.0) cycle
        vmat(1,jp1,0:2)=vmat(1,jp1,0:2)+hpm1(zzab,ab)*(vmat(1,jm1,0:2)-vmat(1,jm1,1:3))
      end do
      do jj=1,3
        ip1=jj+1
        do m=0,1
          vmat(ip1,:,m)=rap(jj)*vmat(1,:,m)-rcp(jj)*vmat(1,:,m+1)
        end do
        do j=2,10
          jm1=minus1s(jj,j)
          if (jm1.eq.0) cycle
          vmat(ip1,j,0:1)=vmat(ip1,j,0:1)+hpm1(zzab,ab)*(vmat(1,jm1,0:1)-vmat(1,jm1,1:2))*indexps(jj,j)
        end do
      end do
      do ip1=5,10
        do jj=1,3
          i=minus1s(jj,ip1)
          if (i.eq.0) cycle
          jrec=jj
          exit
        end do
        vmat(ip1,:,0)=rap(jrec)*vmat(i,:,0)-rcp(jrec)*vmat(i,:,1)
        im1=minus1s(jrec,i)
        if (im1.ne.0) then
          vmat(ip1,:,0)=vmat(ip1,:,0)+hpm1(zzab,ab)*(vmat(im1,:,0)-vmat(im1,:,1))
        end if
        do j=2,10
          jm1=minus1s(jrec,j)
          if (jm1.eq.0) cycle
          vmat(ip1,j,0)=vmat(ip1,j,0)+hpm1(zzab,ab)*indexps(jrec,j)*(vmat(i,jm1,0)-vmat(i,jm1,1))
        end do
      end do
      cvmat=cvmat+vmat(5:10,5:10,0)
    end do

    tmat=0._dp
    fac1=q*(3._dp-2._dp*qr2)
    tmat(1,1)=fac1*smat(1,1)
    tmat(1,2:4)=2._dp*q*smat(1,2:4)+rbp(:)*tmat(1,1)
    tmat(2:4,1)=2._dp*q*smat(2:4,1)+rap(:)*tmat(1,1)
    tmat(1,5)=rbp(1)*tmat(1,2)+hpm1(zzab,ab)*tmat(1,1)+2._dp*q*smat(1,5)-q/zeta(zb,sh2)*smat(1,1)
    tmat(1,6)=rbp(2)*tmat(1,2)+2._dp*q*smat(1,6)
    tmat(1,7)=rbp(3)*tmat(1,2)+2._dp*q*smat(1,7)
    tmat(1,8)=rbp(2)*tmat(1,3)+hpm1(zzab,ab)*tmat(1,1)+2._dp*q*smat(1,8)-q/zeta(zb,sh2)*smat(1,1)
    tmat(1,9)=rbp(3)*tmat(1,3)+2._dp*q*smat(1,9)
    tmat(1,10)=rbp(3)*tmat(1,4)+hpm1(zzab,ab)*tmat(1,1)+2._dp*q*smat(1,10)-q/zeta(zb,sh2)*smat(1,1)
    do j=2,10
      tmat(2:4,j)=rap(:)*tmat(1,j)+2._dp*q*smat(2:4,j)
      do jj=1,3
        ip1=jj+1
        jm1=minus1s(jj,j)
        if (jm1.eq.0) cycle
        tmat(ip1,j)=tmat(ip1,j)+hpm1(zzab,ab)*tmat(1,jm1)*indexps(jj,j)
      end do
    end do
    do ip1=5,10
      tmat(ip1,:)=2._dp*q*smat(ip1,:)
      do jj=1,3
        i=minus1s(jj,ip1)
        if (i.eq.0) cycle
        jrec=jj
        exit
      end do
      tmat(ip1,:)=tmat(ip1,:)+rap(jrec)*tmat(i,:)
      im1=minus1s(jrec,i)
      if (im1.ne.0) then
        tmat(ip1,:)=tmat(ip1,:)+hpm1(zzab,ab)*tmat(im1,:)-q/zeta(za,sh1)*smat(im1,:)
      end if
      do j=5,10
        jm1=minus1s(jrec,j)
        if (jm1.ne.0) then
          tmat(ip1,j)=tmat(ip1,j)+hpm1(zzab,ab)*indexps(jrec,j)*tmat(i,jm1)
        end if
      end do
    end do
    ctmat=ctmat+tmat(5:10,5:10)

  end do
end do

! now contract to form spherical harmonics. Rem: Cartesians:
!     0    1    2    3    4    5    6    7    8    9   10
!     /    s    x    y    z   xx   xy   xz   yy   yz   zz
! Spherical harmonics: 1=z2,2=xz,3=yz,4=xy,5=x2-y2

do ii=1,6
do jj=1,5
  smatl(ii,jj)=sum(csmat(ii,:)*matharmond(jj,:))
  tmatl(ii,jj)=sum(ctmat(ii,:)*matharmond(jj,:))
  vmatl(ii,jj)=sum(cvmat(ii,:)*matharmond(jj,:))
end do
end do
do ii=1,5
do jj=1,5
  s(ii,jj)=sum(smatl(:,jj)*matharmond(ii,:))
  t(ii,jj)=sum(tmatl(:,jj)*matharmond(ii,:))
  v(ii,jj)=sum(vmatl(:,jj)*matharmond(ii,:))
end do
end do

end subroutine d_d_1e_terms


pure subroutine d_l_1e_terms(sh1,sh2,s,t,v)
use nrtype ; use molprops ; use boys
implicit none

! This is the second place where I use d functions. I store these in 1-d arrays, Cartesian & sph harm
! This is the first routine where I use the generalized indices and 'minus 1' Table:
! !   0    1    2    3    4    5    6    7    8    9   10
!     /    s    x    y    z   xx   xy   xz   yy   yz   zz
! Note that this routine still needs checking

! Conventional order, Cartesian: xx, xy, xz, yy, yz, zz
integer(i4b), intent(in) :: sh1, sh2
real(dp), intent(out) :: s(5,4), t(5,4), v(5,4)

integer(i4b) :: nuc, ab, za, zb, zza, zzab, ii, jj, jp1, ip1,i,im1, m, jrec
real(dp) :: fac1, fac2, fac3, zbop, q
real(dp) :: Cx(3),zc, rcp(3), rcp2, efac, qr2
real(dp) :: smat(10,4), tmat(10,4),vmat(10,4,0:3), csmat(6,4),ctmat(6,4),cvmat(6,4)
real(dp) :: rbp(3), bfac(0:3), rap(3)

csmat=0._dp
ctmat=0._dp
cvmat=0._dp

ab=ind2(sh1,sh2)
do za =1,nzet(sh1)
  zza=(za-1)*maxprim
  do zb = 1, nzet(sh2)
    zzab=zza+zb

    ! Start as usual with the S-S overlap
    fac1=pi32/pab(zzab,ab)**1.5_dp
    qr2=qab(zzab,ab)*rab2(ab)
    efac=exp(-qr2)*ddab(zzab,ab)
    smat(1,1)=efac*fac1

    ! Use recursion to get the s-d overlaps from the ss and sp overlaps
    rap=px(:,zzab,ab)-xsh(:,sh1)
    rbp=px(:,zzab,ab)-xsh(:,sh2)
    smat(2:4,1)=smat(1,1)*rap(:)
    fac2=.5_dp*pm1(zzab,ab)*smat(1,1)
    do ip1=5,10
      do jj=1,3
        i=minus1s(jj,ip1)
        if (i.eq.0) cycle
        jrec=jj
        exit
      end do
      smat(ip1,1)=smat(i,1)*rap(jrec)
      im1=minus1s(jrec,i)
      if (im1.eq.0) cycle
      smat(ip1,1)=smat(ip1,1)+fac2
    end do
    do jj=1,3
      smat(:,1+jj)=rbp(jj)*smat(:,1)*djkl(zb,sh2)
      do i=5,10
        im1=minus1s(jj,i)
        if (im1.eq.0) cycle
        smat(i,1+jj)=smat(i,1+jj)+.5_dp*djkl(zb,sh2)*smat(im1,1)*pm1(zzab,ab)*indexps(jj,i)
      end do
    end do
    csmat=csmat+smat(5:10,:)

    ! Now use recursion for the d-s kinetic integrals
    fac3=qab(zzab,ab)*(3._dp-2._dp*qr2)
    tmat(1,1)=fac3*smat(1,1)
    tmat(2:4,1)=rap(:)*tmat(1,1)+2._dp*qab(zzab,ab)*smat(2:4,1)
    zbop=zeta(zb,sh2)*pm1(zzab,ab)
    q=zbop*zeta(za,sh1)
    do ip1=5,10
      do jj=1,3
        i=minus1s(jj,ip1)
        if (i.eq.0) cycle
        jrec=jj
        exit
      end do
      tmat(ip1,1)=rap(jrec)*tmat(i,1)+2._dp*q*smat(ip1,1)
      im1=minus1s(jrec,i)
      if (im1.eq.0) cycle
      tmat(ip1,1)=tmat(ip1,1)+hpm1(zzab,ab)*tmat(im1,1)-zbop*smat(im1,1)
    end do
    do i=5,10
      tmat(i,2:4)=rbp(:)*tmat(i,1)*djkl(zb,sh2)+2._dp*q*smat(i,2:4)
    end do
    do i=5,10
      do jj=1,3
        jp1=jj+1
        im1=minus1s(jj,i)
        if (im1.eq.0) cycle
        tmat(i,jp1)=tmat(i,jp1)+hpm1(zzab,ab)*tmat(im1,1)*indexps(jj,i)*djkl(zb,sh2)
      end do
    end do
    ctmat=ctmat+tmat(5:10,:)
    ! Finally use recursion to get p-d nuclear Coulomb integrals

    fac1=2._dp*sqrt(pab(zzab,ab)/PI)
    do nuc =1,natom
      vmat=0._dp
      Cx=atcoords(nuc,:)
      zc=nuccharges(nuc)
      rcp=px(:,zzab,ab)-Cx
      rcp2=sum(rcp**2)
      call boys_n(3,pab(zzab,ab)*rcp2,bfac)
      vmat(1,1,:)=-zc*fac1*smat(1,1)*bfac(:)
      do m=0,2
        vmat(2:4,1,m)=rap(:)*vmat(1,1,m)-rcp(:)*vmat(1,1,m+1)
      end do
      do ip1=5,10
        do jj=1,3
          i=minus1s(jj,ip1)
          if (i.eq.0) cycle
          jrec=jj
          exit
        end do
        vmat(ip1,1,0:1)=rap(jrec)*vmat(i,1,0:1)-rcp(jrec)*vmat(i,1,1:2)
        im1=minus1s(jrec,i)
        if (im1.eq.0) cycle
        vmat(ip1,1,0:1)=vmat(ip1,1,0:1)+pm1(zzab,ab)*.5_dp*(vmat(1,1,0:1)-vmat(1,1,1:2))
      end do
      do jj=1,3
        jp1=jj+1
        do i=5,10
          vmat(i,jp1,0)=(rbp(jj)*vmat(i,1,0)-rcp(jj)*vmat(i,1,1))*djkl(zb,sh2)
          im1=minus1s(jj,i)
          if (im1.eq.0) cycle
          vmat(i,jp1,0)=vmat(i,jp1,0)+.5_dp*pm1(zzab,ab)*indexps(jj,i)*(vmat(im1,1,0)-vmat(im1,1,1))*djkl(zb,sh2)
        end do
      end do
      cvmat=cvmat+vmat(5:10,:,0)
    end do

  end do
end do

! now contract to form spherical harmonics. Rem: Cartesians:
!     0    1    2    3    4    5    6    7    8    9   10
!     /    s    x    y    z   xx   xy   xz   yy   yz   zz
! Spherical harmonics: 1=z2,2=xz,3=yz,4=xy,5=x2-y2

do ii=1,5
  do jj=1,4
    s(ii,jj)=sum(csmat(:,jj)*matharmond(ii,:))
    t(ii,jj)=sum(ctmat(:,jj)*matharmond(ii,:))
    v(ii,jj)=sum(cvmat(:,jj)*matharmond(ii,:))
  end do
end do

return
end subroutine d_l_1e_terms

end module one_electron_terms
