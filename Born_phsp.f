      subroutine born_phsp(xborn)
      implicit none
      include 'nlegborn.h'
      include 'pwhg_kn.h'
      include 'pwhg_flst.h'
      include 'pwhg_physpar.h'
      include 'PhysPars.h'
      real * 8  :: p(0:3,nlegborn)
      real * 8  :: cmp(0:3,nlegborn)
      real * 8  :: masses(nlegborn)
      real * 8  :: xborn(*)
      real * 8  :: powheginput
      integer   :: iborn,j
      logical, save :: ini=.true., nores
      real * 8  :: q2min

      if(ini) then
         if(powheginput("#nores") == 1) then
            nores = .true.
         else
            nores = .false.
         endif
         ini = .false.
      endif

      if(nores) then
c provide Born phase space for resonance unaware integration

      else
c generate Born phase space from resoance information
         do iborn=1,flst_nborn
            if(flst_bornresgroup(iborn).eq.flst_ibornresgroup) exit
         enddo

         flst_ibornlength = flst_bornlength(iborn)
         flst_ireallength = flst_ibornlength + 1

         do j=1,flst_ibornlength
            kn_masses(j) = physpar_phspmasses(flst_born(j,iborn))
         enddo

         call genphasespace(xborn,flst_ibornlength,flst_born(:,iborn),
     1        flst_bornres(:,iborn),kn_beams,kn_jacborn,
     1        kn_xb1,kn_xb2,kn_sborn,kn_cmpborn,kn_pborn)


c Needed to avoid the photon pole in the Vcharge=0 run
         q2min=powheginput("qmin")**2         
         if(kn_sborn .lt. q2min) kn_jacborn=0d0
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
         
         kn_minmass = 0d0
         
      endif
      end

      subroutine born_suppression(fact)
      implicit none
      real * 8 fact
      fact=1d0
      end

      subroutine rmn_suppression(fact)
      implicit none
      real * 8 fact
      fact=1d0
      end


      subroutine regular_suppression(fact)
      implicit none
      real * 8 fact
      call rmn_suppression(fact)
      end


      subroutine global_suppression(c,fact)
      implicit none
      character * 1 c
      real * 8 fact
      fact=1d0
      end



      subroutine set_fac_ren_scales(muf,mur)
      implicit none
      include 'PhysPars.h'
      include 'nlegborn.h'
      include 'pwhg_flst.h'
      include 'pwhg_flg.h'
      include 'pwhg_kn.h'
      real * 8, intent(out) ::  muf,mur
      logical, save        :: ini=.true., nores
      real *8              :: mu0
      real *8              :: dotp,powheginput
      external dotp,powheginput
      logical, save        :: fixedscale
      if (ini) then
         if(powheginput("#fixedscale").le.0) then
            fixedscale=.false.
         else
            fixedscale=.true.
         endif
         write(*,*) '*****************************************'
         if(fixedscale) then
            write(*,*) ' mu0 = 91.1876 '
         else
            write(*,*) ' mu0= dynamic scale'
         endif
         write(*,*) ' with muf=mur=mu0'
         if(powheginput("#nores") == 1) then
            nores = .true.
         else
            nores = .false.
         endif
         ini=.false.
      endif
      if(fixedscale) then ! fixedscale: mu0 = ph_zmass
        mu0 = 91.1876d0
      else
c     set dynamical scale
         print*, 'dynamical scale still to be set'
         stop

      endif
      muf = mu0
      mur = mu0

      contains
      subroutine no_scales
      write(*,*) "Error in scale setting, exiting ..."
      call pwhg_exit(-1)
      end subroutine no_scales
      end

