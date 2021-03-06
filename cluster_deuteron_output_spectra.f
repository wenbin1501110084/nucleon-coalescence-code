       program cluster
c       ---------------------------------------------------------------
c       cluster production from coalaescenc of nucleons from Hydrodynamics+UrQMD
C       when you use this code, please cite these papers:
C         W. Zhao, L.Zhu, H.Zheng, C.M.Ko and H.Song, Phys. Rev. C 98, no. 5, 054905 (2018).
C         W. Zhao, C.Shen, C.M.Ko, Q.Liu and H.Song,  Phys. Rev. C 102, no. 4, 044912 (2020).   
C         W.~Zhao, K.~j.~Sun, C.~M.~Ko and X.~Luo,[arXiv:2105.14204 [nucl-th]].       
c       ---------------------------------------------------------------
       implicit double precision (a-h,o-z)
        implicit integer*4 (i-n)
        CHARACTER*200 NAME1

       double precision  mpi,mn,mnb,mp,mpb,ml,mlb,mr,
     *                    md,mdb,mt,mtb,mh,mhb,mht,mhtb,mhh,mhhb,
     *                    npt,nbpt,lpt,lbpt,nv2,nbv2,lv2,lbv2
       
        parameter (ndimp=3000,ndim=50000)
        dimension ppix(ndimp),ppiy(ndimp),ppiz(ndimp),
     *            xpi(ndimp),ypi(ndimp),zpi(ndimp),tpi(ndimp),
     *            pnx(ndim),pny(ndim),pnz(ndim),
     *            xn(ndim),yn(ndim),zn(ndim),tn(ndim),
     *                 pnbx(ndim),pnby(ndim),pnbz(ndim),
     *            xnb(ndim),ynb(ndim),znb(ndim),tnb(ndim),
     *            ppx(ndim),ppy(ndim),ppz(ndim),
     *            xp(ndim),yp(ndim),zp(ndim),tp(ndim),
     *            ppbx(ndim),ppby(ndim),ppbz(ndim),
     *            xpb(ndim),ypb(ndim),zpb(ndim),tpb(ndim),
     *            plx(ndim),ply(ndim),plz(ndim),
     *            xl(ndim),yl(ndim),zl(ndim),tl(ndim),
     *            plbx(ndim),plby(ndim),plbz(ndim),
     *            xlb(ndim),ylb(ndim),zlb(ndim),tlb(ndim),
     *            dpx(ndimp),dpy(ndimp),dpz(ndimp), 
     *            dbpx(ndimp),dbpy(ndimp),dbpz(ndimp),  
     *            tpx(ndimp),tpy(ndimp),tpz(ndimp),     
     *            tbpx(ndimp),tbpy(ndimp),tbpz(ndimp),
     *            hpx(ndimp),hpy(ndimp),hpz(ndimp),
     *            hbpx(ndimp),hbpy(ndimp),hbpz(ndimp)
     
        parameter (nbin=31)
       dimension ptbin(nbin),pipt(nbin),npt(nbin),nbpt(nbin),ppt(nbin),
     *            pbpt(nbin),lpt(nbin),lbpt(nbin),piv2(nbin),nv2(nbin),
     *            nbv2(nbin),pv2(nbin),pbv2(nbin),lv2(nbin),lbv2(nbin),
     *            rpt(nbin),dpt(nbin),dbpt(nbin),tpt(nbin),tbpt(nbin),
     *            hpt(nbin),hbpt(nbin),htpt(nbin),htbpt(nbin),
     *            hhpt(nbin),hhbpt(nbin),rv2(nbin),dv2(nbin),dbv2(nbin),
     *            tv2(nbin),tbv2(nbin),hv2(nbin),hbv2(nbin),htv2(nbin),
     *            htbv2(nbin),hhv2(nbin),hhbv2(nbin),dndyneg(15),
     *            anti_three_p(nbin-1),three_p(nbin-1),dndypos(15),
     *            p_spectra(nbin-1),antip_spectra(nbin-1),dndyantip(15),
     *            dndyp(15)

       
       Common/const/hbc2,cut
       
       !nevent is the number of event in one input file. 
       parameter (ndisk=1) 
        open (unit=88, file='input.txt', status='unknown')
       !.....Read input parameters ....... 
        READ(88,*) nevent
        READ(88,*) nmix
        READ(88,*) rd
        READ(88,*) rt
        
        open (unit=200,file="deuteron_spectra.txt",status="unknown")  
        open (unit=201,file="deuteron_dndy.txt",status="unknown")  

        call GETARG(1,NAME1)
        open (unit=100,file=NAME1,status="unknown")

       
        cut=10.
       bin =0.2
       mpi=0.14
       mn=0.94
       mnb=0.94
       mp=0.938
       mpb=0.938
       ml=1.116
       mlb=1.116
       mr=0.77
       md=2.*mp
       mdb=md
        mt=3.*mp
        mtb=mt
        mh=3.*mp
        mhb=mh
        mht=mt+ml-mn
        mhtb=mht
        mhh=mh+ml-mn
        mhhb=mhh
       hbc=0.197
       hbc2=hbc*hbc
       gr=3./4.
       gd=3./4.
       gdb=3./4.
       gt=1./4.
        gtb=1./4.
        gh=1./4.
        ghb=1./4.
        ght=1./4.
        ghtb=1./4.
        ghh=1./4.
        ghhb=1./4.
       sigmar=1.

        omegad=3./4./(mn*rd**2)
        sigmad=1./sqrt(mn*omegad)
        omegat=1./(mn*rt**2)
        sigmat1=1./sqrt(mn*omegat)
        sigmat2=sigmat1
        omegah=1./(mn*1.9661**2)
        sigmah1=1./sqrt(mn*omegah)
        sigmah2=sigmah1
        sigmaht1=sigmat1
        sigmaht2=sigmat2
        sigmahh1=sigmah1
        sigmahh2=sigmah2

c       --------------------------------------------------------------------
c       read in proton, neutron, antiproton, antineutron, lambda and
c       antilambda in midrapidity
c       --------------------------------------------------------------------
       rn=0.
       dn=0.
       dbn=0.
       trn=0.
       tbn=0.
       hn=0.
       hbn=0.
       htn=0.
       htbn=0.
       hhn=0.
       hhbn=0.
       do i=1,nbin
        ptbin(i)=(i-1.0)*bin
       pipt(i)=0.
       piv2(i)=0.
        npt(i)=0.
       nbpt(i)=0.
       ppt(i)=0.
       pbpt(i)=0.
       lpt(i)=0.
       lbpt(i)=0.
       nv2(i)=0.
       nbv2(i)=0.
       pv2(i)=0.
       pbv2(i)=0.
       lv2(i)=0.
       lbv2(i)=0.
       rpt(i)=0.
        dpt(i)=0.
       dbpt(i)=0.
       tpt(i)=0.
       tbpt(i)=0.
       hpt(i)=0.
       hbpt(i)=0.
       htpt(i)=0.
        htbpt(i)=0.
        hhpt(i)=0.
        hhbpt(i)=0.
       rv2(i)=0.
       dv2(i)=0.
       dbv2(i)=0.
       tv2(i)=0.
       tbv2(i)=0.
       hv2(i)=0.
        hbv2(i)=0.
        htv2(i)=0.
       htbv2(i)=0.
       hhv2(i)=0.
       hhbv2(i)=0.
       enddo
       
       do i=1,nbin-1
         three_p(i)=0.0
         anti_three_p(i)=0.0
         p_spectra(i)=0.0
         antip_spectra(i)=0.0
       enddo

       do j =1,15
          dndypos(j)=0.0
          dndyneg(j)=0.0
          dndyantip(j)=0.0
          dndyp(j)=0.0
      enddo

        npis=0
        nns=0
        nnbs=0
        nps=0
        npbs=0
        nls=0
        nlbs=0
                                 
       do idisk=31,30+ndisk                       ! start to circulate I
        iread=nevent/nmix                         ! start to circulate II  !iread=32000/80
       do ie=1,iread
c        npi=0   
       nn=0
       nnb=0
       np=0
       npb=0
       nl=0
       nlb=0
        do i=1,nmix                               ! start to circulate III
        read (100,*) nparticle
       !print *, nparticle
        do j=1,nparticle                          ! start to circulate IV
       read (100,*) nid,mid,xm,px,py,pz,t,x,y,z
       if(nid.eq.0)then
          nparticle=nparticle+1
          continue
       endif
        e=sqrt(px**2+py**2+pz**2+xm**2)
       yy=0.5*dlog((e+pz)/(e-pz))
       if (abs(yy).le. 3.0) then

       if (nid.eq. 2112) then 
       nn=nn+1
       pnx(nn)=px
       pny(nn)=py
       pnz(nn)=pz
       xn(nn)=x
       yn(nn)=y
       zn(nn)=z
       tn(nn)=t
       else
       endif   

       if (nid.eq. -2112) then
        nnb=nnb+1
        pnbx(nnb)=px
        pnby(nnb)=py
        pnbz(nnb)=pz
        xnb(nnb)=x
        ynb(nnb)=y
        znb(nnb)=z
        tnb(nnb)=t
        else
        endif                    ! end if II
       
       
       else
       endif

       if (nid.eq. 2212) then
       np=np+1
       ppx(np)=px
       ppy(np)=py
        ppz(np)=pz
        xp(np)=x
        yp(np)=y
        zp(np)=z
        tp(np)=t
        ptpp=sqrt(px**2+py**2)             
        epp=sqrt(ptpp**2+pz**2+md**2)
        ypp=0.5*dlog((e+pz)/(e-pz))
        do mm=0,13
          yu=mm*1.0*0.6-3.60+0.3
          yd=yu - 0.6
          if(ypp.gt.yd.and.ypp.le.yu)then
             dndyp(mm+1)=dndyp(mm+1)+1.0/0.6
          else
          endif
       enddo 
       
        if (abs(ypp).le. 0.50) then
          do m=1,nbin-1
          if (ptpp.gt.(ptbin(m)).and.ptpp.le.(ptbin(m+1))) then   !11
           p_spectra(m)=p_spectra(m)+1.0
          else                                                          !11
          endif                                                         !11
          enddo
        else
        endif                    ! end if II
       
        else
       endif

       if (nid.eq. -2212) then
       npb=npb+1
       ppbx(npb)=px
       ppby(npb)=py
        ppbz(npb)=pz
        xpb(npb)=x
        ypb(npb)=y
        zpb(npb)=z
        tpb(npb)=t     
        ptpp=sqrt(px**2+py**2)             
        epp=sqrt(ptpp**2+pz**2+md**2)
        ypp=0.5*dlog((e+pz)/(e-pz))
        do mm=0,13
          yu=mm*1.0*0.6-3.60+0.3
          yd=yu-0.6
          if(ypp.gt.yd.and.ypp.le.yu)then
             dndyantip(mm+1)=dndyantip(mm+1)+1.0/0.6
          else
          endif
       enddo 
        if (abs(ypp).le. 0.50) then
          do m=1,nbin-1
          if (ptpp.gt.(ptbin(m)).and.ptpp.le.(ptbin(m+1))) then   !11
           antip_spectra(m)=antip_spectra(m)+1.0
          else                                                          !11
          endif                                                         !11
          enddo
            
       else
       endif

c       if (nid.eq. 3122) then
c       nl=nl+1
c       plx(nl)=px
c       ply(nl)=py
c       plz(nl)=pz
c       xl(nl)=x
c       yl(nl)=y
c       zl(nl)=z
c       tl(nl)=t
c       else
c       endif      

c        if (nid.eq. -3122) then
c        nlb=nlb+1
c       plbx(nlb)=px
c       plby(nlb)=py
c       plbz(nlb)=pz
c       xlb(nlb)=x
c       ylb(nlb)=y
c       zlb(nlb)=z
c       tlb(nlb)=t
c        else
c       endif 
       else
       endif
       enddo                          ! end to circulate III
        enddo                          ! end to circulate IV

 
c        npis=npis+npi
        nns=nns+nn
        nnbs=nnbs+nnb
        nps=nps+np
        npbs=npbs+npb
c       ------------------------
c       rho meson production
c       -----------------------
c       call cluster2(5,gr,ndimp,sigmar,npi,ppix,ppiy,ppiz,mpi,tpi,xpi,
c     *                ypi,zpi,npi,ppix,ppiy,ppiz,mpi,tpi,xpi,ypi,zpi,
c     *                mr,rn,nbin,bin,ptbin,rpt,rv2)       
       
c       ------------------------
c       deuteron production   
c       ------------------------
       call cluster2(-5,gd,ndim,sigmad,nn,pnx,pny,pnz,mn,tn,xn,yn,zn,np,
     *                ppx,ppy,ppz,mp,tp,xp,yp,zp,md,dn,nbin,bin,
     *                ptbin,three_p,dndypos)          ! dpx,dpy,dpz added 

c       -------------------------       
c       antideuteron production
c       -------------------------
       call cluster2(-5,gdb,ndim,sigmad,nnb,pnbx,pnby,pnbz,mnb,tnb,
     *                xnb,ynb,znb,npb,ppbx,ppby,ppbz,mpb,tpb,xpb,ypb,
     *            zpb,mdb,dbn,nbin,bin,ptbin,anti_three_p,dndyneg)           ! dbpx,dbpy,dbpz added 



!c       -------------------------           
!c       triton production
!c       ------------------------
!        call cluster3(5,ndim,gt,sigmat1,sigmat2,nn,pnx,pny,pnz,mn,
!     *                tn,xn,yn,zn,nn,pnx,pny,pnz,mn,tn,xn,yn,zn,
!     *                np,ppx,ppy,ppz,mp,tp,xp,yp,zp,mt,trn,nbin,
!     *                bin,ptbin,three_p,1)           ! tpx,tpy,tpz added 
!        call cluster3(5,ndim,gh,sigmah1,sigmah2,np,ppx,ppy,ppz,mp,tp,xp,
!     *                yp,zp,np,ppx,ppy,ppz,mp,tp,xp,yp,zp,nn,pnx,pny,
!     *                pnz,mn,tn,xn,yn,zn,mh,hn,nbin,bin,ptbin,hpt
!     *                )           ! hpx,hpy,hpz added 


!c       -------------------------           
!c       antitriton production
!c       ------------------------
!        call cluster3(5,ndim,gtb,sigmat1,sigmat2,nnb,pnbx,pnby,
!     *                pnbz,mnb,tnb,xnb,ynb,znb,nnb,pnbx,pnby,pnbz,mnb,
!     *                tnb,xnb,ynb,znb,npb,ppbx,ppby,ppbz,mpb,tpb,xpb,
!     *       ypb,zpb,mtb,trbn,nbin,bin,ptbin,anti_three_p,-1)           ! tbpx,tbpy,tbpz added 

c       -----------------------------
c       helium3 production
!c       -----------------------------
!        call cluster3(5,ndim,gh,sigmah1,sigmah2,np,ppx,ppy,ppz,mp,tp,xp,
!     *                yp,zp,np,ppx,ppy,ppz,mp,tp,xp,yp,zp,nn,pnx,pny,
!     *                pnz,mn,tn,xn,yn,zn,mh,hn,nbin,bin,ptbin,hpt
!     *                )           ! hpx,hpy,hpz added 

!c       -----------------------------
!c       antihelium3 production
!c       -----------------------------
!        call cluster3(5,ndim,ghb,sigmah1,sigmah2,npb,ppbx,ppby,ppbz,mpb,
!     *                tpb,xpb,ypb,zpb,npb,ppbx,ppby,ppbz,mpb,tpb,xpb,
!     *                ypb,zpb,nnb,pnbx,pnby,pnbz,mnb,tnb,xnb,ynb,znb,
!     *                mhb,hbn,nbin,bin,ptbin,hbpt)           ! hbpx,hbpy,hbpz added 

c       -----------------------------
c       hypertriton  production
c       -----------------------------
c        call cluster3(-5,ndim,ght,sigmaht1,sigmaht2,nn,pnx,pny,pnz,mn,
c     *                tn,xn,yn,zn,np,ppx,ppy,ppz,mp,tp,xp,yp,zp,nl,plx,
c     *                ply,plz,ml,tl,xl,yl,zl,mht,htn,nbin,bin,ptbin,
c     *                htpt)

c       -----------------------------
c       antihypertriton  production
c       -----------------------------
c        call cluster3(-5,ndim,ghtb,sigmaht1,sigmaht2,nnb,pnbx,pnby,pnbz,
c     *                mnb,tnb,xnb,ynb,znb,npb,ppbx,ppby,ppbz,mpb,tpb,
c     *                xpb,ypb,zpb,nlb,plbx,plby,plbz,mlb,tlb,xlb,ylb,
c     *                zlb,mhtb,htbn,nbin,bin,ptbin,htbpt)

c       -----------------------------
c       hyper helium3 production
c       -----------------------------
c        call cluster3(5,ndim,ghh,sigmahh1,sigmahh2,np,ppx,ppy,ppz,mp,tp,
c     *                xp,yp,zp,np,ppx,ppy,ppz,mp,tp,xp,yp,zp,nl,plx,ply,
c     *                plz,ml,tl,xl,yl,zl,mhh,hhn,nbin,bin,ptbin,
c     *                hhpt)

c       -----------------------------
c       anti hyper helium3 production
c       -----------------------------
c        call cluster3(5,ndim,ghhb,sigmahh1,sigmahh2,npb,ppbx,ppby,ppbz,
c     *                mpb,tpb,xpb,ypb,zpb,npb,ppbx,ppby,ppbz,mpb,tpb,
c     *                xpb,ypb,zpb,nlb,plbx,plby,plbz,mlb,tlb,xlb,ylb,
c     *                zlb,mhhb,hhbn,nbin,bin,ptbin,hhbpt)
       
c       else
c          endif

       enddo                          ! end to circulate II
       enddo                          ! end to circulate I

!        write (20,*) '#------------pt spectra-------------------------'
        write(200,*)"#pT, dN/dpTdy of p, anti-p, deuteron,anti-deuteron"
        do i=1,nbin-1
        write(200,*)i*0.2-0.1, p_spectra(i)/nevent*1.0,
     *               antip_spectra(i)/nevent*1.0,
     *               three_p(i)/(nevent*nmix*1.0*bin),
     *               anti_three_p(i)/(nevent*nmix*1.0*bin)

        enddo
        write(201,*)"#y, dN/dy of p, anti-p, deuteron,    anti-deuteron"
        do j=0,12

           write(201,*)(j*1.0)*0.6-3.60, dndyp(j+1)/(nevent*1.0),
     *                 dndyantip(j+1)/(nevent*1.0), 
     *                 dndypos(j+1)/(nevent*nmix*1.0),
     *                 dndyneg(j+1)/(nevent*nmix*1.0)
        enddo
       end

       subroutine cluster2(icase,g,ndim,sigma,n1,p1x,p1y,p1z,m1,t1,x1,
     *                      y1,z1,n2,p2x,p2y,p2z,m2,t2,x2,y2,z2,md,dn,
     *                      nbin,bin,ptbin,pt,dndy)                       ! cpx(ndim),cpy(ndim),cpz(ndim) added 
c       ------------------------------------------
c       evaluate two-body coalescence      
c       ------------------------------------------
        implicit double precision (a-h,o-z)
        implicit integer*4 (i-n)
       double precision m1,m2,md 
       dimension p1x(ndim),p1y(ndim),p1z(ndim),t1(ndim),x1(ndim),
     *            y1(ndim),z1(ndim),p2x(ndim),p2y(ndim),p2z(ndim),
     *            t2(ndim),x2(ndim),y2(ndim),z2(ndim),ptbin(nbin),
     *            pt(nbin-1),dndy(15)                   ! cpx(ndim),cpy(ndim),cpz(ndim) added 
       Common/const/hbc2,cut
       sigma2=sigma*sigma

       do i=1,n1
       e1=sqrt(p1x(i)**2+p1y(i)**2+p1z(i)**2+m1**2)   
       ii=1
       if (icase.gt.0) ii=i+1
        if (ii.lt.n2) then
       do j=ii,n2   
        e2=sqrt(p2x(j)**2+p2y(j)**2+p2z(j)**2+m2**2)
        dx=sqrt(((x1(i)-x2(j))**2+(y1(i)-y2(j))**2+(z1(i)-z2(j))**2)/2.)
        if (dx.lt.cut*sigma) then                  ! start if I
        dp=sqrt(((m2*p1x(i)-m1*p2x(j))**2+(m2*p1y(i)-m1*p2y(j))**2
     *     +(m2*p1z(i)-m1*p2z(j))**2)/(m1+m2)**2*2./hbc2)
        if (dp.lt.cut/sigma) then                   ! start if II
       call wigner2(g,sigma2,e1,p1x(i),p1y(i),p1z(i),m1,t1(i),x1(i),
     *               y1(i),z1(i),e2,p2x(j),p2y(j),p2z(j),m2,t2(j),
     *               x2(j),y2(j),z2(j),w,px,py,pz)

       t=sqrt(px**2+py**2)             
       e=sqrt(t**2+pz**2+md**2)
       y=0.5*dlog((e+pz)/(e-pz))
       do mm=0,13
          yu=mm*1.0*0.6-3.60+0.3
          yd=yu-0.6
          if(y.gt.yd.and.y.le.yu)then
             dndy(mm+1)=dndy(mm+1)+w/0.6
          else
          endif
       enddo 
!       if (abs(y).le. 1.0) then
       if (abs(y).le. 0.30) then
       do m=1,nbin-1
       if (t.gt.(ptbin(m)).and.t.le.(ptbin(m+1))) then   !11
        pt(m)=pt(m)+w/0.60
       else                                                          !11
       endif                                                         !11
       enddo
       else
       endif                    ! end if II
       
       
       else
       endif                    ! end if II
       else
       endif   
        enddo
       else
       endif   
        enddo
       return 
       end

       subroutine wigner2(g,sigma2,e1,p1x,p1y,p1z,m1,t1,x1,y1,z1,e2,p2x,
     *                     p2y,p2z,m2,t2,x2,y2,z2,w,px,py,pz)
c       ----------------------------------
c       evaluate two-body wigner function
c       ----------------------------------
        implicit double precision (a-h,o-z)
        implicit integer*4 (i-n)
       double precision m1,m2
       Common/const/hbc2,cut
       betax=(p1x+p2x)/(e1+e2)  
        betay=(p1y+p2y)/(e1+e2)
        betaz=(p1z+p2z)/(e1+e2) 
       call lorentzboost(betax,betay,betaz,t1,x1,
     *                    y1,z1,t1p,x1p,y1p,z1p)   
        call lorentzboost(betax,betay,betaz,e1,p1x,
     *                    p1y,p1z,e1p,p1xp,p1yp,p1zp)
        call lorentzboost(betax,betay,betaz,t2,x2,
     *                    y2,z2,t2p,x2p,y2p,z2p)
        call lorentzboost(betax,betay,betaz,e2,p2x,
     *                    p2y,p2z,e2p,p2xp,p2yp,p2zp)
        if (t1p.gt.t2p) then
       x2p=x2p+(t1p-t2p)*p2xp/e2p
       y2p=y2p+(t1p-t2p)*p2yp/e2p
       z2p=z2p+(t1p-t2p)*p2zp/e2p
        else
       endif   
       if (t2p.gt.t1p) then
        x1p=x1p+(t2p-t1p)*p1xp/e1p
             y1p=y1p+(t2p-t1p)*p1yp/e1p
       z1p=z1p+(t2p-t1p)*p1zp/e1p
       else
       endif   
        w=g*8.*exp(-((x1p-x2p)**2+(y1p-y2p)**2+(z1p-z2p)**2)/2./sigma2)
     *    *exp(-sigma2*((m2*p1xp-m1*p2xp)**2+(m2*p1yp-m1*p2yp)**2
     *    +(m2*p1zp-m1*p2zp)**2)/(m1+m2)**2*2./hbc2)
        px=p1x+p2x
        py=p1y+p2y
        pz=p1z+p2z
       return
       end

       subroutine cluster3(icase,ndim,g,sigma1,sigma2,n1,p1x,
     *                      p1y,p1z,m1,t1,x1,y1,z1,n2,p2x,p2y,p2z,m2,t2,
     *                      x2,y2,z2,n3,p3x,p3y,p3z,m3,t3,x3,y3,z3,
     *                      mt,tn,nbin,bin,ptbin,pt,nnnn) ! cpx(ndim),cpy(ndim),cpz(ndim) added 
c       ------------------------------------------
c       evaluate three-body coalescence      
c       ------------------------------------------
        implicit double precision (a-h,o-z)
        implicit integer*4 (i-n)
       double precision m1,m2,m3,mt 
       dimension p1x(ndim),p1y(ndim),p1z(ndim),t1(ndim),x1(ndim),
     *            y1(ndim),z1(ndim),p2x(ndim),p2y(ndim),p2z(ndim),
     *            t2(ndim),x2(ndim),y2(ndim),z2(ndim),p3x(ndim),
     *            p3y(ndim),p3z(ndim),t3(ndim),x3(ndim),y3(ndim),
     *            z3(ndim),ptbin(nbin),pt(nbin-1)
!     *            cpx(ndim),cpy(ndim),cpz(ndim)                      ! cpx(ndim),cpy(ndim),cpz(ndim) added 
       Common/const/hbc2,cut
       sigma12=sigma1*sigma1
       sigma22=sigma2*sigma2
       do i=1,n1
       ii=1
       if (icase.gt. 0) ii=i+1
        if (ii.lt.n2) then
       e1=sqrt(p1x(i)**2+p1y(i)**2+p1z(i)**2+m1**2)   
        do j=ii,n2    
        e2=sqrt(p2x(j)**2+p2y(j)**2+p2z(j)**2+m2**2)
       do k=1,n3
        e3=sqrt(p3x(k)**2+p3y(k)**2+p3z(k)**2+m3**2)
        dx=sqrt(((x1(i)-x2(i))**2+(y1(i)-y2(j))**2+(z1(i)-z2(j))**2)/2.)
        if (dx.lt.cut*sigma1) then 
        dp=sqrt(((m2*p1x(i)-m1*p2x(j))**2+(m2*p1y(i)-m1*p2y(j))**2
     *      +(m2*p1z(i)-m1*p2z(j))**2)/(m1+m2)**2*2./hbc2)
        if (dp.lt.cut/sigma1) then
        dx=sqrt(((m1*x1(i)+m2*x2(j))/(m1+m2)-x3(k))**2
     *      +((m1*y1(i)+m2*y2(j))/(m1+m2)-y3(k))**2
     *      +((m1*z1(i)+m2*z2(j))/(m1+m2)-z3(k))**2*2./3.)
        if (dx.lt.cut*sigma2) then   
        dp=sqrt(((m3*(p1x(i)+p2x(j))-(m1+m2)*p3x(k))**2
     *      +(m3*(p1y(i)+p2y(j))-(m1+m2)*p3y(k))**2
     *      +(m3*(p1z(i)+p2z(j))-(m1+m2)*p3z(k))**2)/(m1+m2+m3)**2
     *      *3./2./hbc2)
        if (dp.lt.cut/sigma2) then        
        call wigner3(g,sigma12,sigma22,e1,p1x(i),p1y(i),p1z(i),m1,
     *               t1(i),x1(i),y1(i),z1(i),e2,p2x(j),p2y(j),p2z(j),
     *               m2,t2(j),x2(j),y2(j),z2(j),e3,p3x(k),p3y(k),
     *               p3z(k),m3,t3(k),x3(k),y3(k),z3(k),w,px,py,pz)




       t=sqrt(px**2+py**2)             
       e=sqrt(t**2+pz**2+mt**2)
       y=0.5*dlog((e+pz)/(e-pz)) 

       if (abs(y).le. 0.50) then
       do m=1,nbin-1
       if (t.gt.(ptbin(m)).and.t.le.(ptbin(m+1))) then   !11
        pt(m)=pt(m)+w
!!!        v2(m)=v2(m)+(px**2-py**2)/t**2*w        
       else                                                          !11
       endif                                                         !11
       enddo

       else
       endif
       else
       endif   
        else
        endif
       else
       endif
       else
        endif
        enddo
       enddo
       else
       endif   
       enddo
       return 
       end

       subroutine wigner3(g,sigma12,sigma22,e1,p1x,p1y,p1z,m1,t1,x1,y1,
     *                     z1,e2,p2x,p2y,p2z,m2,t2,x2,y2,z2,e3,p3x,p3y,
     *                     p3z,m3,t3,x3,y3,z3,w,px,py,pz)
c       ----------------------------------
c       evaluate three-body wigner function
c       ----------------------------------
        implicit double precision (a-h,o-z)
        implicit integer*4 (i-n)
       double precision m1,m2,m3
       Common/const/hbc2,cut
       betax=(p1x+p2x+p3x)/(e1+e2+e3)  
        betay=(p1y+p2y+p3y)/(e1+e2+e3)
        betaz=(p1z+p2z+p3z)/(e1+e2+e3) 
       call lorentzboost(betax,betay,betaz,t1,x1,
     *                    y1,z1,t1p,x1p,y1p,z1p)   
        call lorentzboost(betax,betay,betaz,e1,p1x,
     *                    p1y,p1z,e1p,p1xp,p1yp,p1zp)
        call lorentzboost(betax,betay,betaz,t2,x2,
     *                    y2,z2,t2p,x2p,y2p,z2p)
        call lorentzboost(betax,betay,betaz,e2,p2x,
     *                    p2y,p2z,e2p,p2xp,p2yp,p2zp)
       call lorentzboost(betax,betay,betaz,t3,x3,
     *                    y3,z3,t3p,x3p,y3p,z3p)   
        call lorentzboost(betax,betay,betaz,e3,p3x,
     *                    p3y,p3z,e3p,p3xp,p3yp,p3zp)       
       if (t1p.gt.t2p.and.t1p.gt.t3p) then
        x2p=x2p+(t1p-t2p)*p2xp/e2p
             y2p=y2p+(t1p-t2p)*p2yp/e2p
       z2p=z2p+(t1p-t2p)*p2zp/e2p
        x3p=x3p+(t1p-t3p)*p3xp/e3p
             y3p=y3p+(t1p-t3p)*p3yp/e3p
       z3p=z3p+(t1p-t3p)*p3zp/e3p
       else
       endif   
        if (t2p.gt.t1p.and.t2p.ge.t3p) then
       x1p=x1p+(t2p-t1p)*p1xp/e1p
       y1p=y1p+(t2p-t1p)*p1yp/e1p
       z1p=z1p+(t2p-t1p)*p1zp/e1p
        x3p=x3p+(t2p-t3p)*p3xp/e3p
             y3p=y3p+(t2p-t3p)*p3yp/e3p
       z3p=z3p+(t2p-t3p)*p3zp/e3p
        else
       endif   
        if (t3p.gt.t1p.and.t3p.ge.t2p) then
       x1p=x1p+(t3p-t1p)*p1xp/e1p
       y1p=y1p+(t3p-t1p)*p1yp/e1p
       z1p=z1p+(t3p-t1p)*p1zp/e1p
        x2p=x2p+(t3p-t2p)*p2xp/e2p
             y2p=y2p+(t3p-t2p)*p2yp/e2p
       z2p=z2p+(t3p-t2p)*p2zp/e2p
        else
        endif   
        w=g*8.*exp(-((x1p-x2p)**2+(y1p-y2p)**2+(z1p-z2p)**2)/2./sigma12)
     *    *exp(-sigma12*((m2*p1xp-m1*p2xp)**2+(m2*p1yp-m1*p2yp)**2
     *    +(m2*p1zp-m1*p2zp)**2)/(m1+m2)**2*2./hbc2)
     *    *8.*exp(-(((m1*x1p+m2*x2p)/(m1+m2)-x3p)**2
     *    +((m1*y1p+m2*y2p)/(m1+m2)-y3p)**2
     *    +((m1*z1p+m2*z2p)/(m1+m2)-z3p)**2)*2./3./sigma22)
     *    *exp(-sigma22*((m3*(p1xp+p2xp)-(m1+m2)*p3xp)**2
     *    +(m3*(p1yp+p2yp)-(m1+m2)*p3yp)**2
     *    +(m3*(p1zp+p2zp)-(m1+m2)*p3zp)**2)/(m1+m2+m3)**2*3./2./hbc2)
        px=p1x+p2x+p3x
        py=p1y+p2y+p3y
        pz=p1z+p2z+p3z
       return
       end

       subroutine lorentzboost(betax,betay,betaz,t,x,y,z,tp,xp,yp,zp)
c       -------------------------
c       loretnz transformation
c       -------------------------
        implicit double precision (a-h,o-z) 
            implicit integer*4 (i-n)
       tp=t
       xp=x
       yp=y
       zp=z
       beta2 = betax**2+betay**2+betaz**2
       if (beta2.gt.1.0e-5) then 
       gamma = 1./sqrt(1.-beta2)
       xlam00= gamma
       xlam01 = -gamma*betax
       xlam02 = -gamma*betay
       xlam03 = -gamma*betaz
       xlam10 = xlam01
       xlam11 = 1.+(gamma-1.)*betax**2/beta2
       xlam12 = (gamma-1.)*betax*betay/beta2
       xlam13 = (gamma-1.)*betax*betaz/beta2
       xlam20 = xlam02
       xlam21 = xlam12
       xlam22 = 1.+(gamma-1.)*betay**2/beta2
       xlam23 = (gamma-1.)*betay*betaz/beta2
       xlam30 = xlam03
       xlam31 = xlam13
       xlam32 = xlam23
       xlam33 = 1.+(gamma-1.)*betaz**2/beta2
       tp = t*xlam00+x*xlam01+y*xlam02+z*xlam03
       xp = t*xlam10+x*xlam11+y*xlam12+z*xlam13
       yp = t*xlam20+x*xlam21+y*xlam22+z*xlam23
       zp = t*xlam30+x*xlam31+y*xlam32+z*xlam33
       else
       endif   
       return
       end

        subroutine binpt(ndim,n,px,py,pz,m,nbin,bin,ptbin,pt,v2) 
    
c       -----------------------------------------   
c       bin transverse momentum and ellptic flow
c       -----------------------------------------

        implicit double precision (a-h,o-z)
        implicit integer*4 (i-n)
        double precision m

        dimension px(ndim),py(ndim),pz(ndim),ptbin(nbin),
     *            pt(nbin),v2(nbin)

        do i=1,n
        t=sqrt(px(i)**2+py(i)**2)
        if (t.le.ptbin(nbin)+bin/2.) then
        do j=1,nbin
        if (t.ge.(ptbin(j)-bin/2.).and.t.lt.(ptbin(j)+bin/2.)) then
        pt(j)=pt(j)+1.
        v2(j)=v2(j)+(px(i)**2-py(i)**2)/amax1(t**2,0.000001)
        else
        endif
        enddo
        else
        endif
        enddo
        return
        end

