!$Id$
module omega

   use precision_mod
   use truncation, only: n_r_max, lm_max, l_max, minc
   use radial_functions, only: r_CMB, r_ICB, i_costf_init, d_costf_init
   use blocking, only: lm2
   use logic, only: lVerbose
   use output_data, only: tag
   use plms_theta, only: plm_theta
   use cosine_transform, only: costf1
   use const, only: one, two, half
 
   implicit none

   private

   public :: outOmega

contains

   subroutine outOmega(z,omega_IC)
      !-----------------------------------------------------------------------
      !   Output of axisymmetric zonal flow omega(s) into field omega.TAG,
      !   where s is the cylindrical radius. This is done for the southern
      !   and norther hemispheres at z=+-(r_icb+0.5)
      !-----------------------------------------------------------------------

      !-- Input variables:
      complex(cp), intent(in) :: z(lm_max,n_r_max)
      real(cp),    intent(in) :: omega_IC
              
      !-- Local stuff:
      real(cp) :: dzVpLMr(l_max+1,n_r_max)

      integer :: nR,lm,l ! counter
      integer :: nNS     ! index for NHS and SHS

      integer, parameter :: nSmax=300
      integer :: nS
      real(cp) ::  sZ,zZ,dsZ
      real(cp) :: rZ,thetaZ
      real(cp) :: VpS,omega(2)
      real(cp) :: workA(lm_max,n_r_max) ! work array

      character(len=64) :: fileName

      if ( lVerbose ) write(*,*) '! Starting outOmega!'

      fileName='omega.'//tag
      open(99, file=fileName, status='unknown')

      dsZ=r_CMB/real(nSmax-1,kind=cp)

      !--- Transform to lm-space for all radial grid points:
      do nR=1,n_r_max
         dzVpLMr(1,nR)=0.0_cp
         do l=1,l_max
            lm=lm2(l,0)
            dzVpLMr(l+1,nR)=real(z(lm,nR))
         end do
      end do

      !---- Transform the contributions to cheb space for z-integral:
      call costf1(dzVpLMr,l_max+1,1,l_max+1,workA,i_costf_init,d_costf_init)

      sZ=0.0_cp
      outer: do nS=1,nSmax
         sZ=sZ+dsZ

         inner: do nNS=1,2  !  North and south hemisphere !

            if ( nNS == 1 ) then  ! south hemisphere !
               zZ=r_ICB+half
            else
               zZ=-(r_ICB+half)
            end if
            rZ    =sqrt(zZ*zZ+sZ*sZ)
            thetaZ=atan2(sZ,zZ)
            if ( rZ > r_CMB ) exit outer

              !------ Get the function values for (sZ,zCy)
            VpS=lnPAS2tr(dzVpLMr,l_max+1,r_ICB,r_CMB, &
                         l_max,minc,n_r_max,thetaZ,rZ)
            omega(nNS)=VpS/(rZ*sin(thetaZ))/omega_IC

         end do inner ! Loop over north and south hemisphere

         write(99,*) sZ,omega(1),omega(2)

      end do outer ! Loop over s

      close(99)

      if ( lVerbose ) write(*,*) '! End of outOmega!'

   end subroutine outOmega
!-----------------------------------------------------------------------------
   real(cp) function lnPAS2tr(f,lmMax,a,b,lMax,minc,nChebMax,theta,r)

      !-- Input variables
      integer,  intent(in) :: lmMax
      real(cp), intent(in) :: f(lmMax,*)
      real(cp), intent(in) :: a,b
      integer,  intent(in) :: lMax
      integer,  intent(in) :: minc
      integer,  intent(in) :: nChebMax
      real(cp), intent(in) :: theta
      real(cp), intent(in) :: r

      !-- Local variables
      real(cp) :: plm(lmMax),dthetaPlm(lmMax)
      integer :: nCheb
      real(cp) :: cheb(nChebMax)
      real(cp) :: ftr
      integer :: l
      real(cp) :: x,chebNorm

      !--- Map r to cheb intervall [-1,1]:
      !    and calculate the cheb polynomia:
      !    Note: the factor chebNorm is needed
      !    for renormalisation. Its not needed if one used
      !    costf1 for the back transform.
      x=two*(r-half*(a+b))/(b-a)
      chebNorm=sqrt(two/real(nChebMax-1,kind=cp))
      cheb(1)=one*chebNorm
      cheb(2)=x*chebNorm
      do nCheb=3,nChebMax
         cheb(nCheb)=two*x*cheb(nCheb-1)-cheb(nCheb-2)
      end do
      cheb(1)       =half*cheb(1)
      cheb(nChebMax)=half*cheb(nChebMax)

      !--- Calculate Legendres:
      call plm_theta(theta,lMax,0,minc,plm,dthetaPlm,lmMax,2)

      !--- Loop to add all contribution functions:
      ftr=0.0_cp
      do l=0,lMax
         if ( l+1 > lmMax ) then
            write(*,*) '! lmMax too small in ln2rt!'
            write(*,*) '! lm',l+1,lMax
            stop
         end if
         do nCheb=1,nChebMax
            ftr=ftr+f(l+1,nCheb)*dthetaPlm(l+1)*cheb(nCheb)
         end do
      end do

      lnPAS2tr=-ftr/(r*sin(theta))

   end function lnPAS2tr
!----------------------------------------------------------------------
end module omega
