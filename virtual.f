      subroutine setvirtual(p,vflav,virtual)
c Wrapper subroutine to call OL Virtual
      use openloops, only: set_parameter
      implicit none
      include 'nlegborn.h'
      include 'pwhg_st.h'
      include 'pwhg_math.h'
      include 'PhysPars.h'
      integer, parameter :: nlegs=nlegbornexternal
      real * 8, intent(in)  :: p(0:3,nlegs)
      integer,  intent(in)  :: vflav(nlegs)
      real * 8, intent(out) :: virtual
      real * 8 :: bornn, virtualCat,s, shift
      real * 8 :: nested1, nested2, nestedtotal
      real * 8 :: auxjk(nlegs,nlegs), auxmunu(0:3,0:3,nlegs)
      
      s= 2.0d0 * ( p(0,1)*p(0,2) - ( p(1,1)*p(1,2) + p(2,1)*p(2,2)
     &     + p(3,1)*p(3,2) ) )

c$$$c      print*, 's= ', s
c$$$c      st_muren2 = s      
c$$$c     Poles
c$$$c      print*, 'polm1:', 4.0d0/3.0d0*(-3d0/2d0+log(s/st_muren2))/pi
c$$$c      print*, 'polm2:', -4.0d0/3.0d0/pi
c$$$c      shift = 4.0d0/3.0d0 * 1d0/12.0d0/pi * ( 7.0d0*pi**2 +
c$$$c     & 18d0*log(s/st_muren2) - 6d0*log(s/st_muren2)**2 )
c$$$c      diff= st_alpha/2.0d0/pi*(virtual-virtualCat)/bornn/st_alpha 
c$$$c     print*, 'diff=', st_alpha/2.0d0/pi*(virtual-virtualCat)/bornn/st_alpha
c$$$c      print*, 'explicit shift= ', diff
c$$$c      print*, 'predicted shift= ', shift
c$$$c      print*, 'diff= ', (diff-shift)/4d0*3d0*12.0d0*pi
c$$$c      print*, 'b= ', bornn
c$$$c      print*, 'v/b= ', virtualCat/(2.0d0*4.0d0/3.d0*bornn)
      
c      call openloops_virtual(p,vflav,virtual,bornn)
c      print*, 'virtual= ', virtual

      call openloops_born(p,vflav,bornn,auxjk,auxmunu)
      !FLVfin
      nested1 = -8.0d0*4.0d0/3.0d0*st_alpha/2.0d0/pi*bornn
      !Finite part of eps expansion of FLVdiv
      nested2 = 4.0d0/3.0d0/12.0d0/pi*st_alpha*( 7.0d0*pi**2 +
     &     18d0*log(s/st_muren2) - 6d0*log(s/st_muren2)**2 )*bornn
      nestedtotal = nested 1 + nested2
      !Remove of finite piece due to the nested normalisation e^epsG/Gamma[1+eps]
      shift = bornn*st_alpha*(4.0d0/3.0d0)*(pi**2/12.0d0/pi)
      virtual = (nestedtotal-shift)/(st_alpha/2.0d0/pi)
      
c      print*, 'nestedtotal= ', (nestedtotal-shift)/(st_alpha/2.0d0/pi)      
      
      end
