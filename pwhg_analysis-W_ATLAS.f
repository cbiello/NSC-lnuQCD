c  The next subroutines open some histograms and prepare them 
c      to receive data 
c  You can substitute these  with your favourite ones

      subroutine init_hist
      implicit none
      include  'LesHouches.h'
      include 'pwhg_math.h'
      character * 6 whcprg      
      common/cwhcprg/whcprg
      real *8 dmt
      real *8 mtmin, mtmax

      call inihists
         
      dmt=10d0

      mtmin=200d0
      mtmax=3000d0
      
cccccccccccccccccccccccccc
c
c Define the histos:
c
c rate and transverse mass
      call bookupeqbins('total',1d0,0d0,1d0)
      call bookupeqbins('rate',1d0,0d0,1d0)
      call bookupeqbins('massT',dmt,mtmin,mtmax)
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
      end

      subroutine analysis(dsig0)
      implicit none
      real * 8 dsig0
      real * 8, allocatable, save :: dsig(:)
      include 'hepevt.h'
      include 'pwhg_math.h' 
      include  'LesHouches.h'
      include 'pwhg_weights.h'
      include 'pwhg_rwl.h'
      logical ini
      data ini/.true./
      save ini

      real * 8 binsize(700)
      common/pwhghistcommon/binsize
      integer ihep
      character * 6 whcprg      
      common/cwhcprg/whcprg
      data whcprg/'NLO   '/

      logical write_cuts
      integer i, mu

      integer nl, nnu
      integer il, inu
      real *8 ptl,yl,etal
      real *8 ptnu,ynu,etanu
      real *8 ylnu,etalnu,philnu,rlnu
      real *8 massT
      real *8 ppairlnu(4)

      real *8 ptlep_cut, ptmiss_cut
      real *8 ptlep_geom
      real *8 ylep_cut
      real *8 mT_min, mT_max
      
c================================================

ccccccccccccccccc
      ptlep_cut = 65d0
      ptmiss_cut = 85d0
      ptlep_geom = 35d0
      ylep_cut   = 2.4d0
      mT_min = 200.0d0
      mT_max = 5000.0d0
ccccccccccccccccc
      
      write_cuts=.false.     
      
      if(ini) then

         write(*,*) ''
         write(*,*) '*****************************'
         write(*,*) '** weights_num     = ',weights_num
         write(*,*) '** rwl_num_weights = ',rwl_num_weights
         write(*,*) '** rwl_num_groups = ',rwl_num_groups
         write(*,*) '*****************************'
         write(*,*) ''

         write_cuts=.true.

         if(weights_num.eq.0.and.rwl_num_weights.eq.0) then
            call setupmulti(1)
         else if(weights_num.ne.0.and.rwl_num_weights.eq.0) then
            call setupmulti(weights_num)
         else if(weights_num.eq.0.and.rwl_num_weights.ne.0) then
            call setupmulti(rwl_num_weights)
         else
            call setupmulti(rwl_num_weights)
         endif

         if(.not. allocated(dsig)) then
            allocate(dsig(max(10,rwl_num_weights+1)))
         endif
          
         ini=.false.
      endif

      dsig(:)=0d0

      if(weights_num.eq.0.and.rwl_num_weights.eq.0) then
         dsig(1)=dsig0
      else if(weights_num.ne.0.and.rwl_num_weights.eq.0) then
         dsig(1:weights_num)=weights_val(1:weights_num)
      else if(weights_num.eq.0.and.rwl_num_weights.ne.0) then
         dsig(1:rwl_num_weights)=rwl_weights(1:rwl_num_weights)
      else
         dsig(1:rwl_num_weights)=rwl_weights(1:rwl_num_weights)
      endif

      if(sum(abs(dsig)).eq.0) return
   
      if(write_cuts) then
         write(*,*) '********************************************'
         write(*,*) '********************************************'
         write(*,*) '                ANALYSIS CUTS               '
         write(*,*) '********************************************'
         write(*,*) '********************************************'
         write(*,*) ''
         write(*,*) 'ATLAS W:'
         write(*,*) 'cuts on leptons'
         write(*,*) 'ptlep_cut= ', ptlep_cut
         write(*,*) 'ptmiss_cut= ', ptmiss_cut
         write(*,*) 'ptlep_geom= ', ptlep_geom
         write(*,*) 'ylep_cut= ', ylep_cut
         write(*,*) ''
         write(*,*) '*******************************************'
         write(*,*) '*******************************************'
      endif

      
      call filld('total',0.5d0,dsig)      
      
c Parton level analysis
      if(whcprg.eq.'NLO'.or.whcprg.eq.'LHE') then


c$$$         do ihep=1,nhep
c$$$            print*, 'ihep= ', ihep, '-> idhep= ',idhep(ihep) 
c$$$         enddo

         nl=0
         nnu=0
         do ihep=3,nhep
            if(abs(idhep(ihep)).eq.11) then
               il=ihep
               nl=nl+1
            endif
            if(abs(idhep(ihep)).eq.12) then
               inu=ihep
               nnu=nnu+1
            endif
         enddo

c$$$         if(nl.ne.1) then
c$$$            write(*,*) "Error in pwhg_analysis: ",nl," charged leps"
c$$$            call exit(1)
c$$$         endif
c$$$         if(nnu.ne.1) then
c$$$            write(*,*) "Error in pwhg_analysis: ",nnu," charged leps"
c$$$            call exit(1)
c$$$         endif


         if(nhep.gt.4) then
!            print*, 'implement the photon isolation'
         endif
         
c Hadron level analysis
      else

         print*, "This code is developed for FO studies,"
         print*, "we should not be here."
         stop

      endif                     !parton or hadron level analysis

ccccc CHECK OF INFINITE WEIGHTS
      if(rwl_num_weights.eq.0) then
         if(dsig0+1 .eq. dsig0) then
             write(*,*) "INF weight. DISCARDING EVENT, weight = ",dsig0
             return
         endif
      else
         do i=1,rwl_num_weights
           if(dsig(i)+1 .eq. dsig(i)) then
             write(*,*) "INF weight. DISCARDING EVENT, i, weight = ",i, dsig(i)
             return
           endif
         enddo
      endif 
cccccccccccccccccccccccccccccccccc

      return

ccccc KINEMATIC RECONSTRUCTION
      call ptyeta(phep(1,il),ptl,yl,etal)
      call ptyeta(phep(1,inu),ptnu,ynu,etanu)

      call getdydetadphidr(phep(1,il),phep(1,inu),
     %    ylnu,etalnu,philnu,rlnu)

      massT = sqrt(2.0d0*ptl*ptnu*(1.0d0-cos(philnu)))

      do mu=1,4
         ppairlnu(mu)=phep(mu,il)+phep(mu,inu)
      enddo

      if(massT.le.mT_min .or. massT.ge.mT_max) return
      if(ptl .lt. ptlep_cut) return
      if(ptnu .lt. ptmiss_cut) return
      if(ptl*ptnu < ptlep_geom**2) return
      if (abs(yl).ge.ylep_cut) return
      if (abs(ynu).ge.ylep_cut) return

      call filld('rate',0.5d0,dsig)
      call filld('massT',massT,dsig)
      
c$$$c bottoms:
c$$$      call ptyeta(phep(1,ib),ptb,yb,etab)
c$$$      call ptyeta(phep(1,ibbar),ptbbar,ybbar,etabbar)
c$$$
c$$$      call getdydetadphidr(phep(1,ib),phep(1,ibbar),
c$$$     %     ybbbar,etabbbar,phibbbar,rbbbar)
c$$$
c$$$      call getdydetadphidrr(phep(1,ib),phep(1,ibbar),
c$$$     %     retabbbar)
c$$$
c$$$c Higgs:
c$$$      call ptyeta(phep(1,ihiggs),pth,yh,etah)
c$$$
c$$$      mhiggs = dsqrt(abs(phep(4,ihiggs)**2-phep(1,ihiggs)**2-
c$$$     &                   phep(2,ihiggs)**2-phep(3,ihiggs)**2))
c$$$
c$$$      call getdydetadphidr(phep(1,ib),phep(1,ihiggs),
c$$$     %     ybh,etabh,phibh,rbh)
c$$$      call getdydetadphidr(phep(1,ibbar),phep(1,ihiggs),
c$$$     %     ybbarh,etabbarh,phibbarh,rbbarh)
c$$$
c$$$      call getdydetadphidrr(phep(1,ib),phep(1,ihiggs),
c$$$     %     retabh)
c$$$      call getdydetadphidrr(phep(1,ibbar),phep(1,ihiggs),
c$$$     %     retabbarh)
c$$$
c$$$c mass of the pair
c$$$      do mu=1,4
c$$$         ppairbbbar(mu)=phep(mu,ib)+phep(mu,ibbar)
c$$$      enddo
c$$$      mbbbar=dsqrt(abs(ppairbbbar(4)**2
c$$$     1            -ppairbbbar(1)**2-ppairbbbar(2)**2-ppairbbbar(3)**2))
c$$$      ptbbbar=dsqrt(abs(ppairbbbar(1)**2+ppairbbbar(2)**2))
c$$$
c$$$c bbbar+Higgs system:
c$$$      pbbh(1:4) = ppairbbbar(1:4)+phep(1:4,ihiggs)
c$$$      ptbbh = dsqrt(abs(pbbh(1)**2+pbbh(2)**2))
c$$$      mbbh = dsqrt(abs(pbbh(4)**2-pbbh(1)**2-pbbh(2)**2-pbbh(3)**2))
c$$$
c$$$c bbbar vs Higgs:
c$$$      call getdydetadphidr(ppairbbbar(1),phep(1,ihiggs),
c$$$     %     ybbh,etabbh,phibbh,rbbh)
c$$$
c$$$      call getdydetadphidrr(ppairbbbar(1),phep(1,ihiggs),
c$$$     %     retabbh)
c$$$
c$$$
c$$$cccccccccccccccccccccccccccccc
c$$$c
c$$$c inclusive analysis:
c$$$c
c$$$c xsec:
c$$$      call filld('xsec',0.5d0,dsig)
c$$$
c$$$
c$$$      if(pth.ge.10) call filld('xsec-ptHge10',0.5d0,dsig)
c$$$      if(pth.ge.15) call filld('xsec-ptHge15',0.5d0,dsig)
c$$$      if(pth.ge.20) call filld('xsec-ptHge20',0.5d0,dsig)
c$$$      
c$$$c Higgs:
c$$$      call filld('pt_Higgs',pth,dsig)
c$$$      call filld('ptzoom_Higgs',pth,dsig)
c$$$      call filld('eta_Higgs',etah,dsig)
c$$$      call filld('y_Higgs',yh,dsig)
c$$$      call filld('mass_Higgs',mhiggs,dsig)
c$$$
c$$$c H-b pair:
c$$$      call filld('deta_Higgs_b',etabh,dsig)
c$$$      call filld('dy_Higgs_b',ybh,dsig)
c$$$      call filld('dphi_Higgs_b',phibh,dsig)
c$$$      call filld('Reta_Higgs_b',retabh,dsig)
c$$$      call filld('Ry_Higgs_b',rbh,dsig)
c$$$
c$$$c H-bbar pair:
c$$$      call filld('deta_Higgs_bbar',etabbarh,dsig)
c$$$      call filld('dy_Higgs_bbar',ybbarh,dsig)
c$$$      call filld('dphi_Higgs_bbar',phibbarh,dsig)
c$$$      call filld('Reta_Higgs_bbar',retabbarh,dsig)
c$$$      call filld('Ry_Higgs_bbar',rbbarh,dsig)
c$$$
c$$$c Higgs-b-bbar system:
c$$$      call filld('pt_Higgsbbbar',ptbbh,dsig)
c$$$      call filld('ptzoom_Higgsbbbar',ptbbh,dsig)
c$$$      call filld('mass_Higgsbbbar',mbbh,dsig)
c$$$
c$$$c bottom and anti-bottom:
c$$$c pt 
c$$$      call filld('pt_b',ptb,dsig)
c$$$      call filld('ptzoom_b',ptb,dsig)
c$$$      call filld('pt_bbar',ptbbar,dsig)
c$$$      call filld('ptzoom_bbar',ptbbar,dsig)
c$$$
c$$$c eta:
c$$$      call filld('eta_b',etab,dsig)
c$$$      call filld('eta_bbar',etabbar,dsig)
c$$$
c$$$c y:
c$$$      call filld('y_b',yb,dsig)
c$$$      call filld('y_bbar',ybbar,dsig)
c$$$
c$$$c bbbar pair:
c$$$      call filld('mass_bbbar',mbbbar,dsig)
c$$$      call filld('pt_bbbar',ptbbbar,dsig)
c$$$      call filld('ptzoom_bbbar',ptbbbar,dsig)
c$$$c
c$$$c b,bbar separation:
c$$$      call filld('deta_b_bbar',etabbbar,dsig)
c$$$      call filld('dy_b_bbar',ybbbar,dsig)
c$$$      call filld('dphi_b_bbar',phibbbar,dsig)
c$$$      call filld('Reta_b_bbar',retabbbar,dsig)
c$$$      call filld('Ry_b_bbar',rbbbar,dsig)
c$$$
c$$$
c$$$c bbbar, Higgs separation:
c$$$      call filld('Reta_bbbar_Higgs',retabbh,dsig)
c$$$      call filld('Ry_bbbar_Higgs',rbbh,dsig)
   
   
      end

cccccccccccccccccccccccccccc

      subroutine ptyeta(p,pt,y,eta)
      implicit none
      real * 8 p(1:4),pt,y,eta
      real * 8 pp,tiny
      parameter (tiny=1d-12)
      pt=sqrt(p(1)**2+p(2)**2)
      y=log((p(4)+p(3))/(p(4)-p(3)))/2
      pp=sqrt(pt**2+p(3)**2)*(1+tiny)
      eta=log((pp+p(3))/(pp-p(3)))/2
      end
 
      subroutine getdydetadphidr(p1,p2,dy,deta,dphi,dr)
      implicit none
      include 'pwhg_math.h' 
      real * 8 p1(1:4),p2(1:4),dy,deta,dphi,dr
      real * 8 y1,eta1,pt1,mass1,phi1
      real * 8 y2,eta2,pt2,mass2,phi2
      call ptyeta(p1,pt1,y1,eta1)
      call ptyeta(p2,pt2,y2,eta2)
      dy=y1-y2
      deta=eta1-eta2
      phi1=atan2(p1(1),p1(2))
      phi2=atan2(p2(1),p2(2))
      dphi=abs(phi1-phi2)
      dphi=min(dphi,2d0*pi-dphi)
      dr=sqrt(dy**2+dphi**2)
      end

      subroutine getdydetadphidrr(p1,p2,dr)
      implicit none
      include 'pwhg_math.h' 
      real * 8 p1(1:4),p2(1:4),dy,deta,dphi,dr
      real * 8 y1,eta1,pt1,mass1,phi1
      real * 8 y2,eta2,pt2,mass2,phi2

      call ptyeta(p1,pt1,y1,eta1)
      call ptyeta(p2,pt2,y2,eta2)

      deta=eta1-eta2
      phi1=atan2(p1(1),p1(2))
      phi2=atan2(p2(1),p2(2))
      dphi=abs(phi1-phi2)
      dphi=min(dphi,2d0*pi-dphi)
      dr=sqrt(deta**2+dphi**2)
      end
