      subroutine init_couplings
      implicit none
      include 'pwhg_par.h'
      include 'PhysPars.h'
      include "pwhg_flg.h" 
      real * 8 powheginput
      external powheginput
c Avoid multiple calls to this subroutine. The parameter file is opened
c but never closed ...
      logical called
      data called/.false./
      save called
      if(called) then
         return
      else
         called=.true.
      endif

      ! there can be issues, especially when having Born cuts, when the 
      ! first 100 times the amplitude is zero, then it will set the whole
      ! amplitude falsely to zero; the following flag prevents this
      flg_cachediscardzeroamps = .false.

c     For V+J, the use of flg_doublefsr is recommended.
c     We make it the default, at least when minnlo is used.
      flg_doublefsr=.true.
      if(powheginput("#doublefsr").eq.0) flg_doublefsr=.false.      
      
      if(powheginput("#verytinypars").eq.1) then
         par_isrtinycsi = 1d-12
         par_isrtinyy = 1d-12
         par_fsrtinycsi = 1d-12
         par_fsrtinyy = 1d-12
      endif

      call param_readin()
      call tophys()

*********************************************************
* Print out of all the relevant couplings
* so that they are in a log file
*********************************************************
      write(*,*) '**************************************************'
      write(*,*) '* init_couplings.f writeout'
      write(*,*) '**************************************************'
      write(*,*) 'alpha=', alpha
      write(*,*) 'gfermi=', gfermi
      write(*,*) 'alfas=', alfas
      write(*,*) 'zmass=', zmass
      write(*,*) 'tmass=', tmass
      write(*,*) 'tmass_phsp=', ph_tmass_phsp
      write(*,*) 'lmass=', lmass
      write(*,*) 'cmass=', cmass
      write(*,*) 'bmass=', bmass
      write(*,*) 'lmass=', lmass
      write(*,*) 'wmass=', wmass
      write(*,*) 'zwidth=', zwidth
      write(*,*) 'wwidth=', wwidth
      write(*,*) 'twidth=', twidth
      write(*,*) 'hmass=', hmass
      write(*,*) 'hwidth=', hwidth
      write(*,*) 'Vud=', Vud
      write(*,*) 'Vus=', Vus
      write(*,*) 'Vub=', Vub
      write(*,*) 'Vcd=', Vcd
      write(*,*) 'Vcs=', Vcs
      write(*,*) 'Vcb=', Vcb
      write(*,*) 'Vtd=', Vtd
      write(*,*) 'Vts=', Vts
      write(*,*) 'Vtb=', Vtb
      write(*,*) '**************************************************'
      end

      subroutine param_readin
        implicit none
        include 'PhysPars.h'
        include 'pwhg_st.h'
        include 'pwhg_flg.h'

        double precision, parameter :: Pi = 3.14159265358979323846d0
        integer             :: ewscheme
        real * 8, external  :: powheginput
        complex * 16         :: cwmass2, czmass2


        zmass=powheginput('#zmass')
        if(zmass<0d0) zmass = 91.1876d0
        tmass=powheginput('#tmass')
        if(tmass<0d0) tmass=173.2d0
        ph_tmass_phsp=powheginput('#tmass_phsp')
        if(ph_tmass_phsp<0d0) ph_tmass_phsp=tmass
        wmass=powheginput('#wmass')
        if(wmass<0d0) wmass = 80.385d0
        bmass=powheginput('#bmass')
        if(bmass<0d0) bmass = 0d0
        hmass=powheginput('#hmass')
        if(hmass<0d0) hmass = 125d0

        zwidth=powheginput('#zwidth')
        if(zwidth<0d0) zwidth=2.4952d0
        wwidth=powheginput('#wwidth')
        if(wwidth<0d0) wwidth=2.0854d0
        twidth = powheginput('#twidth')
        if(twidth<0d0) twidth=1.442620d0
        hwidth=powheginput('#hwidth')
        if(hwidth<0d0) hwidth=0.407d-2


        ewscheme=powheginput('#ewscheme')
        if(ewscheme<1.or.ewscheme>2) ewscheme = 2

        if (ewscheme.eq.1) then  ! MW, MZ, Gmu scheme
          gfermi=powheginput('#gfermi')
          if(gfermi<0d0) gfermi = 1.16639d-5
          if (powheginput("#complexGFermi").eq.0) then
            alpha=sqrt(2d0)/pi*gfermi*abs(wmass**2*(1d0-(wmass/zmass)**2))
          else
            cwmass2=DCMPLX(wmass**2,-wwidth*wmass)
            czmass2=DCMPLX(zmass**2,-zwidth*zmass)
            alpha=sqrt(2d0)/pi*gfermi*abs(cwmass2*(1d0-cwmass2/czmass2))
          endif
        else  ! MW, MZ, alpha scheme
          alpha=powheginput('#alpha')
          if(alpha<0d0) alpha = 1/132.50698d0
          if (powheginput("#complexGFermi").eq.0) then
            gfermi=alpha/sqrt(2d0)*pi/(1d0-(wmass/zmass)**2)/wmass**2
          else
            cwmass2=DCMPLX(wmass**2,-wwidth*wmass)
            czmass2=DCMPLX(zmass**2,-zwidth*zmass)
            gfermi=alpha/sqrt(2d0)*pi/abs((1d0-cwmass2/czmass2)*cwmass2)
          endif
        end if

        alfas = 0.119d0
        lmass = 0d0
        cmass = 0d0

c     POWHEG CKM matrix
c
c        d     s     b
c    u
c    c
c    t

!       Vud=0.97428d0
!       Vus=0.2253d0
!       Vub=0.00347d0
!       Vcd=0.2252d0
!       Vcs=0.97345d0
!       Vcb=0.0410d0
!       Vtd=0.00862d0
!       Vts=0.0403d0
!       Vtb=0.999152d0

        Vud=1d0
        Vus=1d-10
        Vub=1d-10
        Vcd=1d-10
        Vcs=1d0
        Vcb=1d-10
        Vtd=1d-10
        Vts=1d-10
        Vtb=1d0

      end subroutine param_readin


      subroutine tophys
        implicit none
        include 'PhysPars.h'
        include 'pwhg_math.h'
        include 'pwhg_physpar.h'
        real * 8 e_em,g_weak
        real * 8 powheginput
        external powheginput
        integer j
        real * 8 tmass_phsp
        logical massive_leptons

        ph_alphaem=alpha
        ph_gfermi=gfermi

        ph_Zmass = zmass
        ph_Wmass = wmass
        ph_Hmass = hmass
        ph_tmass = tmass
        ph_bmass = bmass

        ph_Zwidth = zwidth
        ph_Wwidth = wwidth
        ph_Hwidth = hwidth
        ph_twidth = twidth

        ph_CKM(1,1)=Vud
        ph_CKM(1,2)=Vus
        ph_CKM(1,3)=Vub
        ph_CKM(2,1)=Vcd
        ph_CKM(2,2)=Vcs
        ph_CKM(2,3)=Vcb
        ph_CKM(3,1)=Vtd
        ph_CKM(3,2)=Vts
        ph_CKM(3,3)=Vtb

c set up masses and widths for resonance damping factors
        physpar_pdgmasses(25) = ph_Hmass
        physpar_pdgmasses(24) = ph_Wmass
        physpar_pdgmasses(23) = ph_Zmass
        physpar_pdgmasses(5) = ph_bmass
        physpar_pdgmasses(6) = ph_tmass
        physpar_pdgmasses(22) = 0
        physpar_pdgwidths(25) = ph_Hwidth
        physpar_pdgwidths(24) = ph_Wwidth
        physpar_pdgwidths(23) = ph_Zwidth
        physpar_pdgwidths(6) = ph_twidth
        physpar_pdgwidths(22) = 0

        do j=1,physpar_npdg
           physpar_pdgmasses(-j) = physpar_pdgmasses(j)
           physpar_pdgwidths(-j) = physpar_pdgwidths(j)
        enddo

        physpar_phspmasses = physpar_pdgmasses
        physpar_phspwidths = physpar_pdgwidths

        physpar_phspmasses(6) = ph_tmass_phsp
        physpar_phspmasses(-6) = ph_tmass_phsp

!   Provide charm and bottom (and lepton) masses for momenta reshuffling in event generation
!   NOTE: it overrides masses (in particular bmass) only for the event generation
        physpar_mq(4)= 1.40d0
        if(powheginput('#c_mass').gt.0) physpar_mq(4)=powheginput('#c_mass') 
        physpar_mq(5)= 4.92d0 
        if(powheginput('#b_mass').gt.0) physpar_mq(5)=powheginput('#b_mass') 
        
        massive_leptons=.false.
        if(powheginput('#massive_leptons').eq.1) massive_leptons=.true.
        if(massive_leptons)then
           physpar_ml(1)=0.511d-3
           if(powheginput('#e_mass').gt.0) physpar_ml(1)=powheginput('#e_mass') 
           physpar_ml(2)=0.1057d0 
           if(powheginput('#mu_mass').gt.0) physpar_ml(2)=powheginput('#mu_mass') 
           physpar_ml(3)=1.777d0 
           if(powheginput('#tau_mass').gt.0) physpar_ml(3)=powheginput('#tau_mass') 
        endif
        
      end subroutine tophys


