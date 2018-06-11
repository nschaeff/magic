module nonlinear_bcs

   use precision_mod
   use geometry, only: nrp, lmP_max, n_phi_max, l_axi,            &
       &                 l_theta, u_theta, n_lmP, n_lm, &
       &                 n_r_cmb, n_r_icb
   use radial_data, only: 
   use radial_functions, only: r_cmb, r_icb, rho0
   use blocking, only: lm2l, lm2m, lm2lmP, lmP2lmPS, lmP2lmPA, nfs, &
       &               sizeThetaB
   use physical_parameters, only: sigma_ratio, conductance_ma, prmag
   use horizontal_data, only: dTheta1S_loc, dTheta1A_loc, dPhi_loc, O_sin_theta, &
       &                      dLh_loc, sn2, cosTheta
   use fft, only: fft_thetab
   use legendre_grid_to_spec, only: legTF2
   use constants, only: two
#ifdef WITH_SHTNS
   use shtns, only: spat_to_SH, spat_to_SH_dist
#endif
   use useful, only: abortRun
   use distributed_theta, only: dist_map

   implicit none

   private

   public :: get_br_v_bcs, get_b_nl_bcs, v_rigid_boundary

contains

   subroutine get_br_v_bcs(br,vt,vp,omega,O_r_E_2,O_rho, &
        &                  br_vt_lm,br_vp_lm)
      !
      !  Purpose of this subroutine is to calculate the nonlinear term    
      !  of the magnetic boundary condition for a conducting mantle or    
      !  inner core in space (r,lm).                                      
      !  Calculation is performed for the theta block:                    
      !
      !  .. code-block:: fortran
      !
      !     n_theta_min<=n_theta<=n_theta_min+n_theta_block-1        
      !
      !  On input br, vt and vp are given on all phi points and           
      !  thetas in the specific block.                                    
      !  On output the contribution of these grid points to all           
      !  degree and orders is stored in br_vt_lm and br_vp_lm.            
      !  Output is [r/sin(theta)*Br*U]=[(0,br_vt_lm,br_vp_lm)]           
      !
    
      !-- input:
      real(cp), intent(in) :: br(n_phi_max,l_theta:u_theta)      ! r**2 * B_r
      real(cp), intent(in) :: vt(n_phi_max,l_theta:u_theta)      ! r*sin(theta) U_theta
      real(cp), intent(in) :: vp(n_phi_max,l_theta:u_theta)      ! r*sin(theta) U_phi
      real(cp), intent(in) :: omega          ! rotation rate of mantle or IC
      real(cp), intent(in) :: O_r_E_2        ! 1/r**2
      real(cp), intent(in) :: O_rho          ! 1/rho0 (anelastic)
    
      !-- Output variables:
      ! br*vt/(sin(theta)**2*r**2)
      complex(cp), intent(inout) :: br_vt_lm(n_lmP)
      ! br*(vp/(sin(theta)**2*r**2)-omega_ma)
      complex(cp), intent(inout) :: br_vp_lm(n_lmP)
    
      !-- Local variables:
      integer :: n_theta     ! number of theta position
      integer :: n_theta_rel ! number of theta position in block
      integer :: n_phi       ! number of longitude
      real(cp) :: br_vt(n_phi_max,l_theta:u_theta)
      real(cp) :: br_vp(n_phi_max,l_theta:u_theta)
      real(cp) :: fac          ! 1/( r**2 sin(theta)**2 )
    
      ! 20180328 Lago
      ! This used to take the thetaBlock into account. I'm deleting it 
      ! because it is extremely confusing right now.
      ! This used to have OMP stuff in it. I'm deleting it because it needs 
      ! to be reworked anyway!
    
      do n_theta_rel=l_theta,u_theta
         fac=O_sin_theta(n_theta)*O_sin_theta(n_theta)*O_r_E_2*O_rho
         do n_phi=1,n_phi_max
            br_vt(n_phi,n_theta)= fac*br(n_phi,n_theta)*vt(n_phi,n_theta)
    
            br_vp(n_phi,n_theta)= br(n_phi,n_theta) * &
                                     ( fac*vp(n_phi,n_theta) - omega )
         end do
      end do
    
      !-- Fourier transform phi 2 m (real 2 complex!)
#ifdef WITH_SHTNS
      call spat_to_SH_dist(br_vt, br_vt_lm)
      call spat_to_SH_dist(br_vp, br_vp_lm)
#else
      if ( .not. l_axi ) then
         call fft_thetab(br_vt, -1)
         call fft_thetab(br_vp, -1)
      end if
    
      !-- Legendre transform contribution of thetas in block:
      call legTF2(l_theta,br_vt_lm,br_vp_lm,br_vt,br_vp)
#endif
    
   end subroutine get_br_v_bcs
   
!----------------------------------------------------------------------------
   subroutine get_b_nl_bcs(bc,br_vt_lm,br_vp_lm,b_nl_bc,aj_nl_bc)
      !
      !  Purpose of this subroutine is to calculate the nonlinear term    
      !  of the magnetic boundary condition for a conducting mantle in    
      !  physical space (phi,theta), assuming that the conductance        
      !  of the mantle is much smaller than that of the core.             
      !  Calculation is performed for the theta block:                    
      !
      !  .. code-block:: fortran
      !
      !      n_theta_min<=n_theta<=n_theta_min+n_theta_block-1        
      !
      ! This function has been θ-parallelized, but not fully tested yet
      ! Please delete this comment once it is certain that it works as 
      ! expected - Lago 20180502
         
      !-- Input variables:
      character(len=3), intent(in) :: bc                 ! Distinguishes 'CMB' and 'ICB'
      complex(cp),      intent(in) :: br_vt_lm(n_lmP)  ! [br*vt/(r**2*sin(theta)**2)]
      complex(cp),      intent(in) :: br_vp_lm(n_lmP)  ! [br*vp/(r**2*sin(theta)**2)

      !-- Output variables:
      complex(cp), intent(out) :: b_nl_bc(n_lm)  ! nonlinear bc for b
      complex(cp), intent(out) :: aj_nl_bc(n_lm) ! nonlinear bc for aj

      !-- Local variables:
      integer :: l,m       ! degree and order
      integer :: lm        ! position of degree and order
      integer :: lmP       ! same as lm but for l running to l_max+1
      integer :: lmPS,lmPA ! lmP for l-1 and l+1
      real(cp) :: fac

      !write(*,"(2A)") "In get_b_nl_bcs with bc=",bc

      if ( bc == 'CMB' ) then

         fac=conductance_ma*prmag

         if (dist_map%lm2(0,0) > 0)  b_nl_bc(dist_map%lm2(0,0)) = (1.0_cp,1.0_cp) 
         if (dist_map%lm2(0,0) > 0) aj_nl_bc(dist_map%lm2(0,0)) = (1.0_cp,1.0_cp) 
         
         do lm=1,n_lm
            l   =dist_map%lm2l(lm)
            m   =dist_map%lm2m(lm)
            if ((l==0) .and. (m==0)) cycle

            lmP =dist_map%lm2lmP(lm)
            lmPS=dist_map%lmP2lmPS(lmP)
            lmPA=dist_map%lmP2lmPA(lmP)
            if ( l > m ) then
               b_nl_bc(lm)= fac/dLH_loc(lm) * (  dTheta1S_loc(lm)*br_vt_lm(lmPS)  &
                                           - dTheta1A_loc(lm)*br_vt_lm(lmPA)  &
                                           +     dPhi_loc(lm)*br_vp_lm(lmP)   )
               aj_nl_bc(lm)=-fac/dLH_loc(lm) * ( dTheta1S_loc(lm)*br_vp_lm(lmPS)  &
                                           - dTheta1A_loc(lm)*br_vp_lm(lmPA)  &
                                           -     dPhi_loc(lm)*br_vt_lm(lmP)    )
            else if ( l == m ) then
               b_nl_bc(lm)= fac/dLH_loc(lm) * ( - dTheta1A_loc(lm)*br_vt_lm(lmPA)  &
                                                + dPhi_loc(lm)*br_vp_lm(lmP)   )
               aj_nl_bc(lm)=-fac/dLH_loc(lm) * ( - dTheta1A_loc(lm)*br_vp_lm(lmPA) &
                                                 - dPhi_loc(lm)*br_vt_lm(lmP)   )
            end if
         end do

      else if ( bc == 'ICB' ) then

         fac=sigma_ratio*prmag
         
         ! if m=0 is in this rank, aj_nl_bc(1) = 1.0
         if (dist_map%lm2(0,0) > 0) aj_nl_bc(dist_map%lm2(0,0)) = (1.0_cp,1.0_cp)  
         
         do lm=1,n_lm
            l   =dist_map%lm2l(lm)
            m   =dist_map%lm2m(lm)
            if ((l==0) .and. (m==0)) cycle

            lmP =dist_map%lm2lmP(lm)
            lmPS=dist_map%lmP2lmPS(lmP)
            lmPA=dist_map%lmP2lmPA(lmP)
            if ( l > m ) then
               aj_nl_bc(lm)=-fac/dLH_loc(lm) * ( dTheta1S_loc(lm)*br_vp_lm(lmPS)   &
                                           - dTheta1A_loc(lm)*br_vp_lm(lmPA)   &
                                           -     dPhi_loc(lm)*br_vt_lm(lmP)    )
            else if ( l == m ) then
               aj_nl_bc(lm)=-fac/dLH_loc(lm) * (- dTheta1A_loc(lm)*br_vp_lm(lmPA) &
                                                - dPhi_loc(lm)*br_vt_lm(lmP)    )
            end if
         end do

      else
         call abortRun('Wrong input of bc into get_b_nl_bcs')
      end if
               
   end subroutine get_b_nl_bcs
  
!-------------------------------------------------------------------------
   subroutine v_rigid_boundary(nR,omega,lDeriv,vrr,vtr,vpr,      &
            &                  cvrr,dvrdtr,dvrdpr,dvtdpr,dvpdpr)
!@>details Distributed version of the v_rigid_boundary; the difference is 
!> basically the size of the data and the indexes. I've removed nThetaStart
!> argument though, 'cause I didn't understand what was its purpose 
!> (it is always called with nThetaStart=1 anyway)
!
!@>author Rafael Lago, MPCDF, December 2017
!-------------------------------------------------------------------------------

      !-- Input of variables:
      integer,  intent(in) :: nR            ! no of radial grid point
      logical,  intent(in) :: lDeriv        ! derivatives required ?
              
      !-- Input of boundary rotation rate
      real(cp), intent(in) :: omega

      !-- output:
      real(cp), intent(out) :: vrr(n_phi_max,l_theta:u_theta)
      real(cp), intent(out) :: vpr(n_phi_max,l_theta:u_theta)
      real(cp), intent(out) :: vtr(n_phi_max,l_theta:u_theta)
      real(cp), intent(out) :: cvrr(n_phi_max,l_theta:u_theta)
      real(cp), intent(out) :: dvrdtr(n_phi_max,l_theta:u_theta)
      real(cp), intent(out) :: dvrdpr(n_phi_max,l_theta:u_theta)
      real(cp), intent(out) :: dvtdpr(n_phi_max,l_theta:u_theta)
      real(cp), intent(out) :: dvpdpr(n_phi_max,l_theta:u_theta)

      !-- Local variables:
      real(cp) :: r2
      integer :: nThetaCalc,nThetaNHS
      integer :: i, j

      if ( nR == n_r_cmb ) then
         r2=r_cmb*r_cmb
      else if ( nR == n_r_icb ) then
         r2=r_icb*r_icb
      else
         write(*,*)
         write(*,*) '! v_rigid boundary called for a grid'
         write(*,*) '! points which is not a boundary !  '
         return
      end if

      nThetaCalc=l_theta-1
      do j=l_theta,u_theta
         nThetaCalc=nThetaCalc+1
         nThetaNHS =(nThetaCalc+1)/2 ! northern hemisphere=odd n_theta
         do i=1,n_phi_max
            vrr(i,j)=0.0_cp
            vtr(i,j)=0.0_cp
            vpr(i,j)=r2*rho0(nR)*sn2(nThetaNHS)*omega
            if ( lDeriv ) then
               cvrr(i,j)  =r2*rho0(nR)*two*cosTheta(nThetaCalc)*omega
               dvrdtr(i,j)=0.0_cp
               dvrdpr(i,j)=0.0_cp
               dvtdpr(i,j)=0.0_cp
               dvpdpr(i,j)=0.0_cp
            end if
         end do
      end do

   end subroutine v_rigid_boundary
!-------------------------------------------------------------------------
end module nonlinear_bcs
