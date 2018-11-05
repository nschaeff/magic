module preCalculations

   use constants
   use num_param
   use output_data
   use precision_mod
   use truncation, only: n_r_max, l_max, minc, n_r_ic_max, nalias, &
       &                 n_cheb_ic_max, m_max, n_cheb_max,         &
       &                 lm_max, n_phi_max, n_theta_max
   use init_fields, only: bots, tops, s_bot, s_top, n_s_bounds,    &
       &                  l_reset_t, topxi, botxi, xi_bot, xi_top, &
       &                  n_xi_bounds, omega_diff
   use parallel_mod, only: rank
   use logic, only: l_mag, l_cond_ic, l_non_rot, l_mag_LF, l_newmap,   &
       &            l_anel, l_heat, l_time_hits,  l_anelastic_liquid,  &
       &            l_cmb_field, l_storeTpot, l_storeVpot, l_storeBpot,&
       &            l_save_out, l_TO, l_TOmovie, l_r_field, l_movie,   &
       &            l_LCR, l_dt_cmb_field, l_storePot, l_non_adia,     &
       &            l_temperature_diff, l_chemical_conv, l_probe,      &
       &            l_precession, l_diff_prec
   use radial_functions, only: rscheme_oc, temp0, r_CMB, ogrun,            &
       &                       r_surface, visc, r, r_ICB, dLtemp0,         &
       &                       beta, rho0, rgrav, dbeta, alpha0,           &
       &                       dentropy0, sigma, lambda, dLkappa, kappa,   &
       &                       dLvisc, dLlambda, divKtemp0, radial,        &
       &                       transportProperties
   use physical_parameters, only: nVarEps, pr, prmag, ra, rascaled, ek,    &
       &                          ekscaled, opr, opm, o_sr, radratio,      &
       &                          sigma_ratio, CorFac, LFfac, BuoFac,      &
       &                          PolInd, nVarCond, nVarDiff, nVarVisc,    &
       &                          rho_ratio_ic, rho_ratio_ma, epsc, epsc0, &
       &                          ktops, kbots, interior_model, r_LCR,     &
       &                          n_r_LCR, mode, tmagcon, oek,             &
       &                          ktopxi, kbotxi, epscxi, epscxi0, sc, osc,&
       &                          ChemFac, raxi, Po, prec_angle, diff_prec_angle
   use horizontal_data, only: horizontal
   use integration, only: rInt_R
   use useful, only: logWrite, abortRun
   use special, only: l_curr, fac_loop

   implicit none

   private

   public :: preCalc, preCalcTimes, writeInfo

contains

   subroutine preCalc
      !
      !  Purpose of this subroutine is to initialize the calc values,     
      !  arrays, constants that are used all over the code.               
      !  The stuff is stored in the common blocks.                        
      !  MPI: This is called by every processors.                         
      !
    
      !---- Local variables:
      real(cp) :: c1,help,facIH
      real(cp) :: delmin,sr_top,si_top,sr_bot,si_bot
      real(cp) :: xir_top,xii_top,xir_bot,xii_bot
      real(cp) :: topconduc, botconduc
      integer :: n,n_r,l,m,l_bot,m_bot,l_top,m_top
      integer :: fileHandle
      character(len=76) :: fileName
      character(len=80) :: message
      real(cp) :: mom(n_r_max)
       
      !-- Determine scales depending on n_tScale,n_lScale :
      if ( n_tScale == 0 ) then
         !----- Viscous time scale:
         tScale=one
      else if ( n_tScale == 1 ) then
         !----- Magnetic time scale:
         tScale=one/prmag
      else if ( n_tScale == 2 ) then
         !----- Thermal time scale:
         tScale=one/pr
      else if ( n_tScale == 3 ) then
         !----- Rotational time scale:
         tScale=one/ek  ! or ekScaled ? (not defined yet...)
         if ( rank==0 ) then
            print*, 'Warning: rotational timescale, be sure to set dtmax large enough !'
         end if
      end if
      if ( n_lScale == 0 ) then
         !----- Outer Core:
         lScale=one
      else if ( n_lScale == 1 ) then
         !----- Total Core:
         lScale=(one-radratio)
      end if
    
      !---- Scale according to scdIFf:
      vScale  =lScale/tScale
      pScale  =tScale**2/lScale**2
      eScale  =vScale*vScale/enscale
      raScaled=ra/lScale**3
      ekScaled=ek*lScale**2
    
      if ( l_cond_ic ) O_sr=one/sigma_ratio
    
      opr=one/pr
      if ( l_mag ) then
         opm=one/prmag
      else
         opm=0.0_cp
      end if

      if ( l_chemical_conv ) osc=one/sc

      ! Note: CorFac is the factor in front of the Coriolis force. In the scaling
      !       used here (viscous time scale) this is simply the inverse Ekman num:
      !          CorFac=1/Ek
      !       In the non-rotating case there is no Coriolis force and
      !          CorFac=0
      !       LFfac is the factor the Lorentz-force in the Navier-Stokes
      !       equation is multiplied with. This depends on the magnetic
      !       scale chosen. In the standart rotating case the B-scale
      !       is sqrt(\rho \lambda \mu \Omega) where \rho is core density,
      !       \lambda is magnetic diffusivity, \mu is magnetic permeability
      !       and \Omega is rotation rate. In this scaling the square of the
      !       dimensionless magnetic field is identical with the Elasser number:
      !          Els=B^2 / (\rho\lambda\mu\Omega)
      !       The factor to be used in front of the LF is then
      !          LFfac=1/(Ek Pm).
      !       This can also be used to bring kinetic and magnetic energy to
      !       the same scale since Ek=Em/(Ek Pm).
      !       If the system is non rotating we chose the magnetic scale
      !       sqrt(\rho \mu) \nu / L where \nu is the viscous diffusivity
      !       and L=ro-ri is the chosen length scale. Then LFfac=1.
      !       Note that for kinematic dynamos the magnetic field strength
      !       is arbitrary. We nevertheless still use the same LFfac.
    
      if ( l_non_rot ) then
         CorFac=0.0_cp
         if ( l_mag .or. l_mag_LF ) then
            LFfac=one
         else
            LFfac=0.0_cp
         end if
      else
         CorFac=one/ekScaled
         if ( l_mag .or. l_mag_LF ) then
            LFfac=one/(ekScaled*prmag)
         else
            LFfac=0.0_cp
         end if
      end if

      oek = CorFac

      if ( l_precession ) CorFac = CorFac*(one+po*cos(prec_angle))
    
      ! Note: BuoFac is the factor used in front of the buoyancy force.
      !       In the scaling used here its
      !          BuoFac=Ra/Pr
      !       where Ra= ( g0 \alpha \delta T L^3)/(\nu \kappa)
      !       with g0 the CMB gravity, \alpha the thermal expansivity,
      !       \delta T the temperature scale, and \kappa the thermal
      !       diffusivity
      if ( l_heat ) then
         BuoFac=raScaled/pr
      else
         BuoFac=0.0_cp
      end if

      if ( l_chemical_conv ) then
         ChemFac=raxi/sc
      else
         ChemFac=0.0_cp
      end if
    
      dtStart=dtStart/tScale
      dtMax  =dtMax/tScale
      if ( .not. l_non_rot ) dtMax=min(dtMax,intfac*ekScaled)
      dtMin  =dtMax/1.0e6_cp
    
      !-- Calculate radial functions for all threads (chebs,r,.....):

      call radial()

      if ( ( l_newmap ) .and. (rank == 0) ) then
         fileName='rNM.'//tag
         open(newunit=fileHandle, file=fileName, status='unknown')
         do n_r=1,n_r_max
            write(fileHandle,'(I4,4ES16.8)') n_r, r(n_r)-r_icb,   &
            &                                rscheme_oc%drx(n_r), &
            &                                rscheme_oc%ddrx(n_r),&
            &                                rscheme_oc%dddrx(n_r)
         end do
         close(fileHandle)
      end if
    
      call transportProperties

      if ( ( l_anel .or. l_non_adia ) .and. ( rank == 0 ) ) then
         ! Write the equilibrium setup in anel.tag
         fileName='anel.'//tag
         open(newunit=fileHandle, file=fileName, status='unknown')
         write(fileHandle,'(11a15)') 'radius', 'temp0', 'rho0', 'beta',         &
         &                          'dbeta', 'grav', 'ds0/dr', 'div(k grad T)', &
         &                          'alpha', 'ogrun', 'dLtemp0'
         do n_r=1,n_r_max
            write(fileHandle,'(11ES16.8)') r(n_r), temp0(n_r), rho0(n_r),    &
            &                             beta(n_r), dbeta(n_r), rgrav(n_r), &
            &                             dentropy0(n_r), divKtemp0(n_r),    &
            &                             alpha0(n_r), ogrun(n_r), dLtemp0(n_r)
         end do
         close(fileHandle)
      end if
    
      !-- Write radial profiles
      if ( l_mag .and. nVarCond > 0 ) then
         fileName='varCond.'//tag
         open(newunit=fileHandle, file=fileName, status='unknown')
         write(fileHandle,'(4a15)') 'radius', 'sigma', 'lambda', 'dLlambda'
         do n_r=n_r_max,1,-1
            write(fileHandle,'(4ES16.8)') r(n_r),sigma(n_r),lambda(n_r), &
            &                             dLlambda(n_r)
         end do
         close(fileHandle)
      end if
    
      if ( ( l_heat .and. nVarDiff > 0  .or. nVarVisc > 0) .and. ( rank == 0 ) ) then
         fileName='varDiff.'//tag
         open(newunit=fileHandle, file=fileName, status='unknown')
         write(fileHandle,'(5a15)') 'radius', 'conductivity', 'kappa', &
         &                          'dLkappa', 'Prandtl'
         do n_r=n_r_max,1,-1
            write(fileHandle,'(5ES16.8)') r(n_r),kappa(n_r)*rho0(n_r), &
            &                             kappa(n_r),dLkappa(n_r),     &
            &                             pr*visc(n_r)/kappa(n_r)
         end do
         close(fileHandle)
      end if
    
      if ( ( nVarVisc > 0 ) .and. (rank == 0) ) then
         fileName='varVisc.'//tag
         open(newunit=fileHandle, file=fileName, status='unknown')
         write(fileHandle,'(7a15)') 'radius', 'dynVisc', 'kinVisc', &
         &                          'dLvisc', 'Ekman', 'Prandtl', 'Pm'
         if ( l_mag ) then
            do n_r=n_r_max,1,-1
               write(fileHandle,'(7ES16.8)') r(n_r),visc(n_r)*rho0(n_r), &
               &                             visc(n_r),dLvisc(n_r),      &
               &                             ek*visc(n_r),               &
               &                             pr*visc(n_r)/kappa(n_r),    &
               &                             prmag*visc(n_r)/lambda(n_r)
            end do
         else
            do n_r=n_r_max,1,-1
               write(fileHandle,'(7ES16.8)') r(n_r),visc(n_r)*rho0(n_r), &
               &                             visc(n_r), dLvisc(n_r),     &
               &                             ek*visc(n_r),               &
               &                             pr*visc(n_r)/kappa(n_r),    &
               &                             prmag
            end do
         end if
         close(fileHandle)
      end if
    
      l_LCR=.false.
      n_r_LCR=0
      do n_r=1,n_r_max
         if ( r_LCR <= r(n_r)/r_CMB ) then
             l_LCR=.true.
             n_r_LCR=n_r
         end if
      end do
      if ( n_r_LCR == 1 ) then
         l_LCR=.false.
         n_r_LCR=0
      end if
    
      !-- Compute some constants:
      vol_ic=four*third*pi*r_icb**3             ! Inner core volume
      vol_oc=four*third*pi*(r_cmb**3-r_icb**3)  ! Outer core volume
      surf_cmb=four*pi*r_cmb**2                ! Surface of CMB
    
      !-- Initialize everything that has to do with the horizontal representation
      !   on all threads:
      call horizontal
    
      !-- Computation of the average density (useful to compute Re and Rm)
      if ( l_anel ) then
         mom(:)=r(:)**2 * rho0(:)
         mass=four*pi/vol_oc*rInt_R(mom,r,rscheme_oc)
      else
         mass=one
      end if
    
      !-- Calculate auxiliary arrays containing effective Courant grid intervals:
      c1=one/real(l_max*(l_max+1),kind=cp)
      delxh2(1)      =c1*r_cmb**2
      delxh2(n_r_max)=c1*r_icb**2
      delxr2(1)      =(r(1)-r(2))**2
      delxr2(n_r_max)=(r(n_r_max-1)-r(n_r_max))**2
      do n_r=2,n_r_max-1
         delxh2(n_r)=c1*r(n_r)**2
         delmin=min((r(n_r-1)-r(n_r)),(r(n_r)-r(n_r+1)))
         delxr2(n_r)=delmin*delmin
      end do

      !Differential precession of mantle and outer core

      if ( l_diff_prec ) then
         omega_diff = oek * (cos(diff_prec_angle) - 1)
      else
         omega_diff = 0.0_cp
      end if

      !-- Constants used for rotating core or mantle:
      y10_norm=half*sqrt(three/pi)  ! y10=y10_norm * cos(theta)
      y11_norm=half*sqrt(1.5_cp/pi) ! y11=y11_norm * sin(theta)
    
      !----- Proportionality factor of (l=1,m=0) toroidal velocity potential
      !      and inner core rotation rate:
      c_z10_omega_ic=y10_norm/(r(n_r_max)*r(n_r_max))/rho0(n_r_max)
    
      !----- Proportionality factor of (l=1,m=0) toroidal velocity potential
      !      and mantle rotation rate:
      c_z10_omega_ma=y10_norm/(r(1)*r(1))/rho0(1)
    
      !----- Inner-core normalized moment of inertia:
      c_moi_ic=8.0_cp*pi/15.0_cp*r_icb**5*rho_ratio_ic*rho0(n_r_max)
    
      !----- Outer-core normalized moment of inertia:
      ! _moi_oc=8.0_cp*pi/15.0_cp*(r_cmb**5-r_icb**5) ! rho=cst
      do n_r=1,n_r_max
         mom(n_r)=r(n_r)**4 * rho0(n_r)
      end do
      c_moi_oc=8.0_cp*third*pi*rInt_R(mom,r,rscheme_oc)
    
      !----- Mantle normalized moment of inertia:
      c_moi_ma=8.0_cp*pi/15.0_cp*(r_surface**5-r_cmb**5)*rho_ratio_ma
    
      !----- IC normalised moment of inertia / r_icb**4 * 3/(8 pi)
      c_dt_z10_ic=0.2_cp*r_icb*rho_ratio_ic*rho0(n_r_max)
    
      !----- Mantle normalised moment of inertia / r_cmb**4 * 3/(8 pi)
      c_dt_z10_ma=0.2_cp*r_cmb*rho_ratio_ma * ( (r_surface/r_cmb)**5 - one )
    
      !----- Proportionality factor for ic lorentz_torque as used in
      !      ic torque-equation (z10):
      c_lorentz_ic=0.25_cp*sqrt(three/pi)/(r(n_r_max)*r(n_r_max))
    
      !----- Proportionality factor for mantle lorentz_torque as used in
      !      mantle torque-equation (z10):
      c_lorentz_ma=0.25_cp*sqrt(three/pi)/(r(1)*r(1))
    
      !-- Set thermal boundary conditions for fixed temp. on both boundaries:
      !----- Extract tops and bots
      if ( l_heat ) then
    
         !-- Renormalisation of heat flow coeffs and heat source.
         !      This has been copied from the Christensen and
         !      means the we use a different normalisation of
         !      input tops,bots and epsc0 than in the code.
         !      Spherical harmonics of the input tops, bots and epsc0
         !      are normalised so that the integral of Y*Y(c.c.)
         !      (c.c. stands for conjugate complex)
         !      over a spherical surface of radius 1 gives 4 Pi.
         !      In the code we use totally normalized spherical harmonic,
         !      i.e. the integral of Y*Y(c.c.) over a spherical surface
         !      of radius 1 is unity
         epsc=epsc0*sq4pi
    
         do m=0,m_max,minc
            do l=m,l_max
               bots(l,m)=zero
               tops(l,m)=zero
               do n=1,n_s_bounds
                  l_bot =int(s_bot(4*n-3))
                  m_bot =int(s_bot(4*n-2))
                  sr_bot=s_bot(4*n-1)
                  si_bot=s_bot(4*n)
                  l_top =int(s_top(4*n-3))
                  m_top =int(s_top(4*n-2))
                  sr_top=s_top(4*n-1)
                  si_top=s_top(4*n)
                  if ( l_bot == l .and. m_bot == m .and. &
                       cmplx(sr_bot,si_bot,kind=cp) /= zero ) then
                     if ( m == 0 ) si_bot=0.0_cp
                     bots(l,m)=sq4pi*cmplx(sr_bot,si_bot,kind=cp)
                     if ( kbots == 2 ) bots(l,m)=bots(l,m)*lScale
                  end if
                  if ( l_top == l .and. m_top == m .and. &
                       cmplx(sr_top,si_top,kind=cp) /= zero ) then
                     if ( m == 0 ) si_top=0.0_cp
                     tops(l,m)=sq4pi*cmplx(sr_top,si_top,kind=cp)
                     if ( ktops == 2 ) tops(l,m)=tops(l,m)*lScale
                  end if
               end do
            end do
         end do
    
         if ( nVarEps==0 ) then
            facIH=vol_oc
         else if ( nVarEps==1 ) then
            facIH=mass*vol_oc
         else if ( nVarEps == 2 ) then
            mom(:)=r(:)**2*rho0(:)*temp0(:)
            facIH=four*pi*rInt_R(mom,r,rscheme_oc)
         end if

         if ( l_temperature_diff .or. l_anelastic_liquid ) then
            topconduc = rho0(1)*kappa(1)
            botconduc = rho0(n_r_max)*kappa(n_r_max)
         else
            topconduc = rho0(1)*kappa(1)*temp0(1)
            botconduc = rho0(n_r_max)*kappa(n_r_max)*temp0(n_r_max)
         end if
    
         if ( ktops == 1 .and. kbots == 1 ) then ! Fixed entropy
    
            tops(0,0)=-r_icb**2/(r_icb**2+r_cmb**2)*sq4pi
            bots(0,0)= r_cmb**2/(r_icb**2+r_cmb**2)*sq4pi

         else if ( ktops == 3 .and. kbots == 3 ) then ! Fixed temperature contrast

            tops(0,0)=-r_icb**2/(r_icb**2+r_cmb**2)*sq4pi
            bots(0,0)= r_cmb**2/(r_icb**2+r_cmb**2)*sq4pi
    
         else if ( (ktops==2 .and. kbots==2) .or. (ktops == 4 .and. kbots==4) ) then
    
            if ( real(bots(0,0)) > 0.0_cp ) then
               write(*,*)
               write(*,*) '! NOTE: you have supplied'
               write(*,*) '! s_bot(l=0,m=0)>0 which '
               write(*,*) '! means there is a heat '
               write(*,*) '! flux into the inner core.'
               write(*,*) '! This is unrealistic!'
               write(*,*) '! Use s_bot(l=0,m=0)<0 !'
               call abortRun('Stop run in preCalc')
            end if
    
            !--- |epsc0|=1 signifies that the heat production rate is used as
            !    temperature scale! I make sure here that the heat flux through
            !    inner and outer boundary, bots and tops, balance the sources.
            if ( abs(epsc0) == one ) then
    
               !--- Make sure that all the heat comes out at CMB
               if ( tops(0,0) == 0.0_cp ) then
    
                  !--- Compensate by flux from ICB:
                  !    all over the core :
                  if ( epsc0 >= 0 ) then
                     write(*,*) '! NOTE: when the flux through the '
                     write(*,*) '! outer boundary is zero, sinks in'
                     write(*,*) '! the outer core need to balance  '
                     write(*,*) '! the flux from the ICB. Thus we  '
                     write(*,*) '! need epsc<0 !                   '
                     call abortRun('Stop run in preCalc')
                  end if
                  bots(0,0)=epsc*pr*facIH/(four*pi*r_icb**2 * botconduc )
                  call logWrite( &
                       '! CMB heat flux set to balance volume sources!')
    
               else if ( tops(0,0) /= 0.0_cp .and. bots(0,0) == 0.0_cp ) then
    
                  !--- Correct tops to balance inner sources/sinks:
                  if ( epsc0 <= 0 .and. bots(0,0) == 0.0_cp  ) then
                     write(*,*) '! NOTE: when the flux through the '
                     write(*,*) '! ICB is zero we need sources in  '
                     write(*,*) '! the outer core which means      '
                     write(*,*) '! epsc0>0.                        '
                     call abortRun('Stop run in preCalc')
                  end if
                  if ( abs(real(tops(0,0))) == sq4pi ) &
                       call logWrite('! You intend to use the CMB flux as buoy. scale??')
                  tops(0,0)=-facIH*epsc*pr/(four*pi*r_cmb**2) + &
                            radratio**2*botconduc*bots(0,0)
                  if ( tops(0,0) /= help ) call logWrite( &
                       '!!!! WARNING: CMB heat flux corrected !!!!')
    
               else if ( tops(0,0) /= 0.0_cp .and. bots(0,0) /= 0.0_cp ) then
    
                  !--- Correct tops to balance inner sources/sinks:
                  if ( epsc0 <= 0 .and. bots(0,0) == 0.0_cp  ) then
                     write(*,*) '! NOTE: when the flux through the '
                     write(*,*) '! ICB is zero we need sources in  '
                     write(*,*) '! the outer core which means      '
                     write(*,*) '! epsc0>0.                        '
                     call abortRun('Stop run in preCalc')
                  end if
                  help=four*pi*opr*facIH *            &
                       (r_icb**2*real(bots(0,0))*botconduc - &
                        r_cmb**2*real(tops(0,0))*topconduc)
                  if ( help /= epsc ) then
                     write(*,*) '! NOTE: when flux BC through the '
                     write(*,*) '! ICB and CMB are used the sources '
                     write(*,*) '! have to balance the total flux.'
                     call abortRun('Stop run in preCalc')
                  end if
    
               end if
    
            else if ( epsc0 == 0.0_cp .and. ( tops(0,0) /= 0.0_cp .or. &
                                            bots(0,0) /= 0.0_cp ) ) then
               !--- Correct epsc0 to balance the difference between
               !    flux through the inner and outer boundary:
               if ( index(interior_model,'PNS_0V2S') /= 0 .and. &
                    ktops==2 .and. kbots==2 ) then
                  call logWrite('! FIXING bots = f(tops) / epsc==0')
                  bots(0,0) = r_cmb**2*real(tops(0,0))*topconduc/(r_icb**2*botconduc)
                  write(message, '(''! real(bots(0,0))       ='',ES16.6)') real(bots(0,0))
                  call logWrite(message)
                  write(message, '(''! real(bots(0,0))/sq4pi ='',ES16.6)') real(bots(0,0))/sq4pi
                  call logWrite(message)
               end if
               epsc=four*pi/pr/facIH *          &
                    (r_icb**2*real(bots(0,0))*botconduc - &
                     r_cmb**2*real(tops(0,0))*topconduc)
               call logWrite( &
                    '! Sources introduced to balance surface heat flux!')
               write(message,'(''!      epsc0='',ES16.6)') epsc/sq4pi
               call logWrite(message)
            else if ( epsc0 /= 0.0_cp .and. tops(0,0) /= 0.0_cp .and. &
                                          bots(0,0) /= 0.0_cp ) then
               help=four*pi/pr/facIH *          &
                    (r_icb**2*real(bots(0,0))*botconduc - &
                     r_cmb**2*real(tops(0,0))*topconduc)
               if ( help /= epsc ) then
                  write(*,*) '! NOTE: when flux BC through the '
                  write(*,*) '! ICB and/or CMB is used the sources '
                  write(*,*) '! have to balance it.'
                  call abortRun('Stop run in preCalc')
               end if
            end if
    
         end if
    
         if ( l_non_adia ) then
            bots(0,0)=0.0_cp
            tops(0,0)=0.0_cp
            epsc=0.0_cp
            !epsc=four*pi/pr/epsS/vol_oc *          &
            !     (r_icb**2*dtemp0(n_r_max)*rho0(n_r_max)*kappa(n_r_max) - &
            !      r_cmb**2*dtemp0(1)*rho0(1)*kappa(1))*sq4pi
         end if
         if ( ktops == 1 ) then
            write(message,'(''! Const. entropy at outer boundary S ='',ES16.6)') &
                  real(tops(0,0))/sq4pi
            call logWrite(message)
         else if ( ktops == 2 ) then
            help=surf_cmb*topconduc*real(tops(0,0))/sq4pi
            write(message,'(''! Const. total outer boundary entropy flux    ='',ES16.6)') help
            call logWrite(message)
         else if ( ktops == 3 ) then
            write(message,'(''! Const. temp. at outer boundary S ='',ES16.6)') &
                  real(tops(0,0))/sq4pi
            call logWrite(message)
         else if ( ktops == 4 ) then
            help=surf_cmb*topconduc*real(tops(0,0))/sq4pi
            write(message,'(''! Const. total outer boundary temp. flux    ='',ES16.6)') help
            call logWrite(message)
         end if
         if ( kbots == 1 ) then
            write(message,'(''! Const. entropy at inner boundary S ='',ES16.6)') real(bots(0,0))/sq4pi
            call logWrite(message)
         else if ( kbots == 2 ) then
            help=surf_cmb*radratio**2*botconduc*real(bots(0,0))/sq4pi
            write(message, '(''! Const. total inner boundary entropy flux    ='',ES16.6)') help
            call logWrite(message)
         else if ( kbots == 3 ) then
            write(message,'(''! Const. temp. at inner boundary S ='',ES16.6)') real(bots(0,0))/sq4pi
            call logWrite(message)
         else if ( kbots == 4 ) then
            help=surf_cmb*radratio**2*botconduc*real(bots(0,0))/sq4pi
            write(message, '(''! Const. total inner boundary temp. flux    ='',ES16.6)') help
            call logWrite(message)
         end if
         help=facIH*pr*epsc/sq4pi
         write(message,'(''! Total vol. buoy. source ='',ES16.6)') help
         call logWrite(message)
    
      end if

      !-- Set  boundary conditions for chemical composition
      if ( l_chemical_conv ) then
         epscxi=epscxi0*sq4pi
    
         do m=0,m_max,minc
            do l=m,l_max
               botxi(l,m)=zero
               topxi(l,m)=zero
               do n=1,n_xi_bounds
                  l_bot =int(xi_bot(4*n-3))
                  m_bot =int(xi_bot(4*n-2))
                  xir_bot=xi_bot(4*n-1)
                  xii_bot=xi_bot(4*n)
                  l_top =int(xi_top(4*n-3))
                  m_top =int(xi_top(4*n-2))
                  xir_top=xi_top(4*n-1)
                  xii_top=xi_top(4*n)
                  if ( l_bot == l .and. m_bot == m .and. &
                       cmplx(xir_bot,xii_bot,kind=cp) /= zero ) then
                     if ( m == 0 ) xii_bot=0.0_cp
                     botxi(l,m)=sq4pi*cmplx(xir_bot,xii_bot,kind=cp)
                     if ( kbotxi == 2 ) botxi(l,m)=botxi(l,m)*lScale
                  end if
                  if ( l_top == l .and. m_top == m .and. &
                       cmplx(xir_top,xii_top,kind=cp) /= zero ) then
                     if ( m == 0 ) xii_top=0.0_cp
                     topxi(l,m)=sq4pi*cmplx(xir_top,xii_top,kind=cp)
                     if ( ktopxi == 2 ) topxi(l,m)=topxi(l,m)*lScale
                  end if
               end do
            end do
         end do
    
         facIH=vol_oc
         topconduc = rho0(1)
         botconduc = rho0(n_r_max)
    
         if ( ktopxi == 1 .and. kbotxi == 1 ) then ! Fixed chemical comp
    
            topxi(0,0)=-r_icb**2/(r_icb**2+r_cmb**2)*sq4pi
            botxi(0,0)= r_cmb**2/(r_icb**2+r_cmb**2)*sq4pi

         else if ( (ktopxi==2 .and. kbotxi==2) ) then
    
            if ( real(botxi(0,0)) > 0.0_cp ) then
               write(*,*)
               write(*,*) '! NOTE: you have supplied'
               write(*,*) '! xi_bot(l=0,m=0)>0 which '
               write(*,*) '! means there is a composition '
               write(*,*) '! flux into the inner core.'
               write(*,*) '! This is unrealistic!'
               write(*,*) '! Use xi_bot(l=0,m=0)<0 !'
               call abortRun('Stop run in preCalc')
            end if
    
            if ( abs(epscxi0) == one ) then
    
               if ( topxi(0,0) == 0.0_cp ) then
    
                  !--- Compensate by flux from ICB:
                  !    all over the core :
                  if ( epscxi0 >= 0 ) then
                     write(*,*) '! NOTE: when the flux through the '
                     write(*,*) '! outer boundary is zero, sinks in'
                     write(*,*) '! the outer core need to balance  '
                     write(*,*) '! the flux from the ICB. Thus we  '
                     write(*,*) '! need epscxi<0 !                   '
                     call abortRun('Stop run in preCalc')
                  end if
                  botxi(0,0)=epscxi*sc*facIH/(four*pi*r_icb**2*botconduc)
                  call logWrite( &
                       '! CMB heat flux set to balance volume sources!')
    
               else if ( topxi(0,0) /= 0.0_cp .and. botxi(0,0) == 0.0_cp ) then
    
                  !--- Correct topxi to balance inner sources/sinks:
                  if ( epscxi0 <= 0 .and. botxi(0,0) == 0.0_cp  ) then
                     write(*,*) '! NOTE: when the flux through the '
                     write(*,*) '! ICB is zero we need sources in  '
                     write(*,*) '! the outer core which means      '
                     write(*,*) '! epscxi0>0.                        '
                     call abortRun('Stop run in preCalc')
                  end if
                  if ( abs(real(topxi(0,0))) == sq4pi ) &
                       call logWrite('! You intend to use the CMB flux as buoy. scale??')
                  topxi(0,0)=-facIH*epscxi*sc/(four*pi*r_cmb**2) + &
                       radratio**2*botxi(0,0)*botconduc
                  if ( topxi(0,0) /= help ) call logWrite( &
                       '!!!! WARNING: CMB composition flux corrected !!!!')
    
               else if ( topxi(0,0) /= 0.0_cp .and. botxi(0,0) /= 0.0_cp ) then
    
                  !--- Correct tops to balance inner sources/sinks:
                  if ( epscxi0 <= 0 .and. botxi(0,0) == 0.0_cp  ) then
                     write(*,*) '! NOTE: when the flux through the '
                     write(*,*) '! ICB is zero we need sources in  '
                     write(*,*) '! the outer core which means      '
                     write(*,*) '! epscxi0>0.                      '
                     call abortRun('Stop run in preCalc')
                  end if
                  help=four*pi/sc/facIH *            &
                       (r_icb**2*real(botxi(0,0))*botconduc - &
                        r_cmb**2*real(topxi(0,0))*topconduc)
                  if ( help /= epscxi ) then
                     write(*,*) '! NOTE: when flux BC through the '
                     write(*,*) '! ICB and CMB are used the sources '
                     write(*,*) '! have to balance the total flux.'
                     call abortRun('Stop run in preCalc')
                  end if
    
               end if
    
            else if ( epscxi0 == 0.0_cp .and. ( topxi(0,0) /= 0.0_cp .or. &
                                            botxi(0,0) /= 0.0_cp ) ) then
               !--- Correct epscxi0 to balance the difference between
               !    flux through the inner and outer boundary:
               epscxi=four*pi/sc/facIH *                     &
                      (r_icb**2*real(botxi(0,0))*botconduc - &
                       r_cmb**2*real(topxi(0,0))*topconduc)
               call logWrite( &
                    '! Sources introduced to balance surface Fickian flux!')
               write(message,'(''!      epscxi0='',ES16.6)') epscxi/sq4pi
               call logWrite(message)
            else if ( epscxi0 /= 0.0_cp .and. topxi(0,0) /= 0.0_cp .and. &
                                          botxi(0,0) /= 0.0_cp ) then
               help=four*pi/sc/facIH *                     &
                    (r_icb**2*real(botxi(0,0))*botconduc - &
                     r_cmb**2*real(topxi(0,0))*topconduc)
               if ( help /= epscxi ) then
                  write(*,*) '! NOTE: when flux BC through the '
                  write(*,*) '! ICB and/or CMB is used the sources '
                  write(*,*) '! have to balance it.'
                  call abortRun('Stop run in preCalc')
               end if
            end if
    
         end if
    
         if ( ktopxi == 1 ) then
            write(message,'(''! Constant comp. at CMB T ='',ES16.6)') &
                  real(topxi(0,0))/sq4pi
            call logWrite(message)
         else if ( ktopxi == 2 ) then
            help=surf_cmb*topconduc*real(topxi(0,0))/sq4pi
            write(message,'(''! Const. total CMB comp. flux    ='',ES16.6)') help
            call logWrite(message)
         end if
         if ( kbotxi == 1 ) then
            write(message,'(''! Constant comp. at ICB T ='',ES16.6)') real(botxi(0,0))/sq4pi
            call logWrite(message)
         else if ( kbotxi == 2 ) then
            help=surf_cmb*radratio**2*botconduc*real(botxi(0,0))/sq4pi
            write(message, '(''! Const. total ICB comp. flux    ='',ES16.6)') help
            call logWrite(message)
         end if
         help=facIH*sc*epscxi/sq4pi
         write(message,'(''! Total vol. comp. source ='',ES16.6)') help
         call logWrite(message)
    
      end if


      !--- Compute fac_loop for current carrying loop

      if ( l_curr ) then

         allocate(fac_loop(l_max))
              
         do l=1,l_max
               
            fac_loop(l)=0.0_cp

            if (mod(l,2)/=0) then
               if(l==1) then
                  fac_loop(l)= half
               else
                  fac_loop(l)= -fac_loop(l-2)*0.64_cp*real(l,kind=cp)/ &
                  &            real(l-1,kind=cp)
               end if
            end if
               
         end do
          
      end if
       

   end subroutine preCalc
!-------------------------------------------------------------------------------
   subroutine preCalcTimes(time,n_time_step)
      !
      !  Precalc. after time, time and dthas been read from startfile.
      !

      !-- Output variables
      real(cp), intent(out) ::  time
      integer,  intent(out) :: n_time_step

      !-- Local variables:
      logical :: l_time
      integer :: n

      !----- Set time step:
      if ( l_reset_t ) then
         time=0.0_cp
         n_time_step=0
      end if

      tmagcon=tmagcon+time

      !-- Get output times:
      l_time_hits=.false.

      call get_hit_times(t_graph,n_time_hits,n_t_graph,l_time, &
           &             t_graph_start,t_graph_stop,dt_graph,  &
           &             n_graphs,n_graph_step,'graph',time,tScale)
      l_time_hits=l_time_hits .or. l_time

      call get_hit_times(t_rst,n_time_hits,n_t_rst,l_time, &
           &             t_rst_start,t_rst_stop,dt_rst,    &
           &             n_rsts,n_rst_step,'rst',time,tScale)
      l_time_hits=l_time_hits .or. l_time

      call get_hit_times(t_log,n_time_hits,n_t_log,l_time, &
           &             t_log_start,t_log_stop,dt_log,    &
           &             n_logs,n_log_step,'log',time,tScale)
      l_time_hits=l_time_hits .or. l_time

      call get_hit_times(t_spec,n_time_hits,n_t_spec,l_time, &
           &             t_spec_start,t_spec_stop,dt_spec,   &
           &             n_specs,n_spec_step,'spec',time,tScale)
      l_time_hits=l_time_hits .or. l_time

      if ( l_probe ) then
         l_probe=.false.
         call get_hit_times(t_probe,n_time_hits,n_t_probe,l_time, &
              &             t_probe_start,t_probe_stop,dt_probe,  &
              &             n_probe_out,n_probe_step,'probe',time,tScale)
         if ( n_probe_out > 0 .or. n_probe_step > 0 .or. l_time ) l_probe= .true.
         l_time_hits=l_time_hits .or. l_time
      end if

      if ( l_cmb_field ) then
         l_cmb_field=.false.
         call get_hit_times(t_cmb,n_time_hits,n_t_cmb,l_time, &
                               t_cmb_start,t_cmb_stop,dt_cmb, &
                         n_cmbs,n_cmb_step,'cmb',time,tScale)
         if ( n_cmbs > 0 .or. n_cmb_step > 0 .or. l_time ) l_cmb_field= .true. 
         l_time_hits=l_time_hits .or. l_time
      end if
      l_dt_cmb_field=l_dt_cmb_field .and. l_cmb_field

      if ( l_r_field ) then
         l_r_field=.false.
         call get_hit_times(t_r_field,n_time_hits,n_t_r_field,l_time, &
                           t_r_field_start,t_r_field_stop,dt_r_field, &
                          n_r_fields,n_r_field_step,'r_field',time,tScale)
         if ( n_r_fields > 0 .or. n_r_field_step > 0 .or. l_time ) l_r_field= .true. 
         l_time_hits=l_time_hits .or. l_time
      end iF

      if ( l_movie ) then
         call get_hit_times(t_movie,n_time_hits,n_t_movie,l_time, &
                            t_movie_start,t_movie_stop,dt_movie,  &
                  n_movie_frames,n_movie_step,'movie',time,tScale)
          l_time_hits=l_time_hits .or. l_time
      end if

      if ( l_TO ) then
         if ( n_TOs == 0 .and. n_t_TO == 0 ) n_TO_step=max(3,n_TO_step)
         call get_hit_times(t_TO,n_time_hits,n_t_TO,l_time, &
                                t_TO_start,t_TO_stop,dt_TO, &
                          n_TOs,n_TO_step,'TO',time,tScale)
         l_time_hits=l_time_hits.or.l_time
         if ( n_TOZs == 0 .and. n_t_TOZ == 0 ) n_TOZs=3
         call get_hit_times(t_TOZ,n_time_hits,n_t_TOZ,l_time, &
                               t_TOZ_start,t_TOZ_stop,dt_TOZ, &
                         n_TOZs,n_TOZ_step,'TOZ',time,tScale)
         l_time_hits=l_time_hits .or. l_time
      end if

      if ( l_TOmovie ) then
         call get_hit_times(t_TOmovie,n_time_hits,n_t_TOmovie,l_time, &
                           t_TOmovie_start,t_TOmovie_stop,dt_TOmovie, &
               n_TOmovie_frames,n_TOmovie_step,'TOmovie',time,tScale)
         l_time_hits=l_time_hits .or. l_time
      end if

      if ( l_storePot ) then
         l_storeVpot   =.true.
         n_Vpot_step   =n_pot_step
         n_Vpots       =n_pots
         t_Vpot_start  =t_pot_start
         t_Vpot_stop   =t_pot_stop
         dt_Vpot       =dt_pot
         l_storeBpot   =.true.
         n_Bpot_step   =n_pot_step
         n_Bpots       =n_pots
         t_Bpot_start  =t_pot_start
         t_Bpot_stop   =t_pot_stop
         dt_Bpot       =dt_pot
         l_storeTpot   =.true.
         n_Tpot_step   =n_pot_step
         n_Tpots       =n_pots
         t_Tpot_start  =t_pot_start
         t_Tpot_stop   =t_pot_stop
         dt_Tpot       =dt_pot
         do n=1,n_time_hits
            t_Bpot(n)=t_pot(n)
            t_Vpot(n)=t_pot(n)
            t_Tpot(n)=t_pot(n)
         end do
      end if

      if ( l_storeBpot ) then
         call get_hit_times(t_Bpot,n_time_hits,n_t_Bpot,l_time, &
                              t_Bpot_start,t_Bpot_stop,dt_Bpot, &
                        n_Bpots,n_Bpot_step,'Bpot',time,tScale)
         l_time_hits=l_time_hits .or. l_time
      end if

      if ( l_storeVpot ) then
         call get_hit_times(t_Vpot,n_time_hits,n_t_Vpot,l_time, &
                              t_Vpot_start,t_Vpot_stop,dt_Vpot, &
                        n_Vpots,n_Vpot_step,'Vpot',time,tScale)
         l_time_hits=l_time_hits .or. l_time
      end if

      if ( l_storeTpot ) then
         call get_hit_times(t_Tpot,n_time_hits,n_t_Tpot,l_time, &
                              t_Tpot_start,t_Tpot_stop,dt_Tpot, &
                        n_Tpots,n_Tpot_step,'Tpot',time,tScale)
         l_time_hits=l_time_hits .or. l_time
      end if

   end subroutine preCalcTimes
!-------------------------------------------------------------------------------
   subroutine get_hit_times(t,n_t_max,n_t,l_t,t_start,t_stop,dt, &
                             n_tot,n_step,string,time,tScale)
      !
      ! This subroutine checks whether any specific times t(*) are given
      ! on input. If so, it returns their number n_r and sets l_t        
      ! to true. If not, t(*) may also be defined by giving a time step  
      ! dt or a number n_tot of desired output times and t_stop>t_start. 
      !

      !-- Input variables:
      integer,          intent(in) :: n_t_max    ! Dimension of t(*)
      real(cp),         intent(in) :: time       ! Time of start file
      real(cp),         intent(in) :: tScale
      character(len=*), intent(in) :: string
      integer,  intent(inout) :: n_tot       ! No. of output (times) if no times defined
      integer,  intent(inout) :: n_step      ! Ouput step in no. of time steps
      real(cp), intent(inout) :: t(n_t_max)  ! Times for output
      real(cp), intent(inout) :: t_start     ! Starting time for output
      real(cp), intent(inout) :: t_stop      ! Stop time for output
      real(cp), intent(inout) :: dt          ! Time step for output

      !-- Output variables
      integer, intent(out) :: n_t        ! No. of output times
      logical, intent(out) :: l_t        ! =.true. if output times are defined

      !-- Local variables:
      integer :: n         ! Counter


      t_start=t_start/tScale
      t_stop =t_stop/tScale
      dt     =dt/tScale

      !-- Check whether any time is given explicitly:
      l_t=.false.
      n_t=0
      do n=1,n_t_max
         if ( t(n) >= 0.0_cp ) then
            t(n)=t(n)/tScale
            l_t=.true.
            n_t=n_t+1
         end if
      end do

      !-- Check times should be constructed:
      if ( t_start < time ) t_start=time
      if ( .not. l_t .and. ( dt > 0.0_cp .or. &
         ( n_tot > 0 .and. t_stop > t_start ) ) ) then

         if ( n_tot > 0 .and. dt > 0.0_cp ) then
            n_t  =n_tot
            n_tot=0
         else if ( dt > 0.0_cp ) then
            if ( t_stop > t_start ) then
               n_t=int((t_stop-t_start)/dt)+1
            else
               n_t=n_t_max
            end if
         else if ( n_tot > 0 ) then
            n_t=n_tot
            n_tot=0
            dt=(t_stop-t_start)/real(n_t-1,kind=cp)
         end if
         if ( n_t > n_t_max ) then
            write(*,*) '! Sorry, maximum no. of times for'
            write(*,*) '! output ',string
            write(*,*) '! is:',n_t_max
            write(*,*) '! Increase n_time_hits in c_output.f!'
            call abortRun('Stop run in get_hit_times')
         end if

         l_t=.true.
         if ( t_start == time ) then
            n_t=n_t-1
            t(1)=t_start+dt
         else
            t(1)=t_start
         end if

         do n=2,n_t
            t(n)=t(n-1)+dt
         end do

      end if


      if ( n_tot /= 0 .AND. n_step /= 0 ) then
         write(*,*)
         write(*,*) '! You have to either provide the total'
         write(*,*) '! number or the step for output:'
         write(*,'(A,2(A,I10))') string, "n_tot = ",n_tot,", n_step = ",n_step
         write(*,*) '! I set the step width to zero!'
         n_step=0
      end if

      if ( l_t ) then
         t_start=t(1)
         t_stop =t(n_t)
         dt     =t(2)-t(1)
      end if
            
   end subroutine get_hit_times
!------------------------------------------------------------------------------
   subroutine writeInfo(n_out)
      !
      !  Purpose of this subroutine is to write the namelist to           
      !  file unit n_out. This file has to be open before calling this    
      !  routine.                                                         
      !

      !-- Input variable:
      integer, intent(in) :: n_out

      if ( rank == 0 ) then

         !-- Output of mode:
         write(n_out,*)
         if ( mode == 0 ) then
            write(n_out,*) '! Self consistent dynamo integration.'
         else if ( mode == 1 ) then
            write(n_out,*) '! Convection integration.'
         else if ( mode == 2 ) then
            write(n_out,*) '! Kinematic dynamo integration.'
         else if ( mode == 3 ) then
            write(n_out,*) '! Magnetic decay modes.'
         else if ( mode == 4 ) then
            write(n_out,*) '! Magneto convection.'
         else if ( mode == 5 ) then
            write(n_out,*) '! Linear onset of convection.'
         else if ( mode == 6 ) then
            write(n_out,*) '! Self consistent dynamo integration without LF.'
         else if ( mode == 7 ) then
            write(n_out,*) '! Super-rotating IC, no convection, no dynamo.'
         else if ( mode == 8 ) then
            write(n_out,*) '! Super-rotating IC, no convection, dynamo.'
         else if ( mode == 9 ) then
            write(n_out,*) '! Super-rotating IC, no convection, dynamo, no LF.'
         else if ( mode == 10 ) then
            write(n_out,*) '! Super-rotating IC, no advection, no convection, no dynamo.'
         else if ( mode == 11 ) then
            write(n_out,*) '! Viscous flow, no inertia, no rotation, no dynamo.'
         end if
      end if
      if ( mode > 11 ) then
         call abortRun('Mode > 11 not implemented !')
      end if

      if (rank == 0) then
         !-- Output of name lists:
         write(n_out, '('' ! Normalized OC moment of inertia:'',ES14.6)') c_moi_oc
         write(n_out, '('' ! Normalized IC moment of inertia:'',ES14.6)') c_moi_ic
         write(n_out, '('' ! Normalized MA moment of inertia:'',ES14.6)') c_moi_ma
         write(n_out, '('' ! Normalized IC volume           :'',ES14.6)') vol_ic
         write(n_out, '('' ! Normalized OC volume           :'',ES14.6)') vol_oc
         write(n_out, '('' ! Normalized IC surface          :'',ES14.6)')  &
                       surf_cmb*radratio**2
         write(n_out, '('' ! Normalized OC surface          :'',ES14.6)') surf_cmb
         write(n_out,*)
         write(n_out,*) '! Grid parameters:'
         write(n_out,'(''  n_r_max      ='',i6, &
              &   '' = number of radial grid points'')') n_r_max
         write(n_out,'(''  n_cheb_max   ='',i6)') n_cheb_max
         write(n_out,'(''  max cheb deg.='',i6)') n_cheb_max-1
         write(n_out,'(''  n_phi_max    ='',i6, &
              &   '' = no of longitude grid points'')') n_phi_max
         write(n_out,'(''  n_theta_max  ='',i6, &
              &   '' = no of latitude grid points'')') n_theta_max
         write(n_out,'(''  n_r_ic_max   ='',i6, &
              &   '' = number of radial grid points in IC'')') n_r_ic_max
         write(n_out,'(''  n_cheb_ic_max='',i6)') n_cheb_ic_max-1
         write(n_out,'(''  max cheb deg ='',i6)') 2*(n_cheb_ic_max-1)
         write(n_out,'(''  l_max        ='',i6, '' = max degree of Plm'')') l_max
         write(n_out,'(''  m_max        ='',i6, '' = max oder of Plm'')') m_max
         write(n_out,'(''  lm_max       ='',i6, '' = no of l/m combinations'')') lm_max
         write(n_out,'(''  minc         ='',i6, '' = longitude symmetry wave no'')') minc
         write(n_out,'(''  nalias       ='',i6, &
              &   '' = spher. harm. deal. factor '')') nalias

      end if

   end subroutine writeInfo
!------------------------------------------------------------------------------
end module preCalculations
