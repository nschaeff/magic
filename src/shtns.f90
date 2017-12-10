#include "perflib_preproc.cpp"
module shtns

   use precision_mod, only: cp
   use constants, only: ci
   use truncation
   use horizontal_data, only: dLh, dLH_loc, D_m_loc, gauss, theta_ord, D_m, &
                              O_sin_theta_E2
   use radial_functions, only: r
   use parallel_mod
   use fft, only: init_fft_phi, finalize_fft_phi, fft_phi_loc
   
   use libspina, only: spprintmatrix
   
   implicit none

   include "shtns.f"

   private

   public :: init_shtns, scal_to_spat, scal_to_grad_spat, pol_to_grad_spat, &
             torpol_to_spat, pol_to_curlr_spat, torpol_to_curl_spat,        &
             torpol_to_dphspat, spat_to_SH,&
             spat_to_SH_dist, sh_to_spat_dist, &
             sph_to_spat_dist, torpol_to_spat_dist, &
             pol_to_curlr_spat_dist, pol_to_grad_spat_dist, &
             torpol_to_dphspat_dist, torpol_to_curl_spat_dist
             
   public :: finalize_shtns

contains

   subroutine init_shtns()
   
   !----------------------------------------------------------------------------
   ! Lago: This is important! minc here is sort of the opposite of 
   ! mres in SHTns. In MagIC, the number of m modes stored is 
   ! n_m_max=m_max/minc+1. In SHTns, it is simply m_max.
   ! Conversely, there are m_max m modes in MagIC, though not all of them
   ! are effectivelly stored/computed. In SHTns, the number of m modes
   ! is m_max*mres+1. This makes things very confusing specially because lm2 
   ! and lmP2 arrays (and their relatives) store m_max m points, and sets those 
   ! which are not multiple of minc to 0.
   !
   ! In other words, this:
   ! call shtns_lmidx(lm_idx, l_idx, m_idx/minc)
   ! returns the same as 
   ! lm_idx = lm2(l_idx, m_idx)
   ! provided that m_idx is a multiple of minc.
   !----------------------------------------------------------------------------
      integer :: norm

      if ( rank == 0 ) then
         call shtns_verbose(1)
      end if

      call shtns_use_threads(0)

      norm = SHT_ORTHONORMAL + SHT_NO_CS_PHASE

      call shtns_set_size(l_max, m_max/minc, minc, norm)
      call shtns_precompute(SHT_GAUSS, SHT_PHI_CONTIGUOUS, &
                            1.e-10_cp, n_theta_max, n_phi_max)
      call shtns_save_cfg(0)


      call shtns_set_size(l_max+1, m_max/minc, minc, norm)
      call shtns_precompute(SHT_QUICK_INIT, SHT_PHI_CONTIGUOUS, &
                            1.e-10_cp, n_theta_max, n_phi_max)
      call shtns_save_cfg(1)
      
      if ( rank == 0 ) then
         call shtns_verbose(0)
      end if
      
      call shtns_load_cfg(0)
      
      ! Might want to add a condition for initializing this only if necessary
      call init_fft_phi
      
   end subroutine
!-------------------------------------------------------------------------------
   subroutine finalize_shtns()
!@>author Rafael Lago, MPCDF, July 2017
!-------------------------------------------------------------------------------
      call finalize_fft_phi
   end subroutine finalize_shtns
!-------------------------------------------------------------------------------
   subroutine scal_to_spat(Slm, fieldc)
      ! transform a spherical harmonic field into grid space
      complex(cp), intent(in) :: Slm(lm_max)
      real(cp), intent(out) :: fieldc(n_phi_max, n_theta_max)

      call shtns_SH_to_spat(Slm, fieldc)

   end subroutine scal_to_spat
!-------------------------------------------------------------------------------
   subroutine scal_to_grad_spat(Slm, gradtc, gradpc)
      ! transform a scalar spherical harmonic field into it's gradient
      ! on the grid
      complex(cp), intent(in) :: Slm(lm_max)
      real(cp), intent(out) :: gradtc(n_phi_max, n_theta_max)
      real(cp), intent(out) :: gradpc(n_phi_max, n_theta_max)

      call shtns_sph_to_spat(Slm, gradtc, gradpc)

   end subroutine scal_to_grad_spat
!-------------------------------------------------------------------------------
   subroutine pol_to_grad_spat(Slm, gradtc, gradpc)

      complex(cp), intent(in) :: Slm(lm_max)
      real(cp), intent(out) :: gradtc(n_phi_max, n_theta_max)
      real(cp), intent(out) :: gradpc(n_phi_max, n_theta_max)

      ! local
      complex(cp) :: Qlm(lm_max)
      integer :: lm

      do lm = 1, lm_max
         Qlm(lm) = dLh(lm) * Slm(lm)
      end do

      call shtns_sph_to_spat(Qlm, gradtc, gradpc)

   end subroutine pol_to_grad_spat
!-------------------------------------------------------------------------------
   subroutine torpol_to_spat(Wlm, dWlm, Zlm, vrc, vtc, vpc)
      complex(cp), intent(in) :: Wlm(lm_max), dWlm(lm_max), Zlm(lm_max)
      real(cp), intent(out) :: vrc(n_phi_max, n_theta_max)
      real(cp), intent(out) :: vtc(n_phi_max, n_theta_max)
      real(cp), intent(out) :: vpc(n_phi_max, n_theta_max)

      ! local
      complex(cp) :: Qlm(lm_max)
      integer :: lm

      do lm = 1, lm_max
         Qlm(lm) = dLh(lm) * Wlm(lm)
      end do

      call shtns_qst_to_spat(Qlm, dWlm, Zlm, vrc, vtc, vpc)

   end subroutine torpol_to_spat
!-------------------------------------------------------------------------------
   subroutine torpol_to_dphspat(dWlm, Zlm, dvtdp, dvpdp)
      !
      ! Computes horizontal phi derivative of a toroidal/poloidal field
      !
      complex(cp), intent(in) :: dWlm(lm_max), Zlm(lm_max)
      real(cp), intent(out) :: dvtdp(n_phi_max, n_theta_max)
      real(cp), intent(out) :: dvpdp(n_phi_max, n_theta_max)

      ! local
      complex(cp) :: Slm(lm_max), Tlm(lm_max)
      integer :: lm, it, ip
      real(cp) :: m

      do lm = 1, lm_max
         m = D_m(lm)
         Slm(lm) = ci*m*dWlm(lm)
         Tlm(lm) = ci*m*Zlm(lm)
      end do

      call shtns_sphtor_to_spat(Slm, Tlm, dvtdp, dvpdp)

      do it=1, n_theta_max
         do ip=1, n_phi_max
            dvtdp(ip, it) = dvtdp(ip, it) * O_sin_theta_E2(it)
            dvpdp(ip, it) = dvpdp(ip, it) * O_sin_theta_E2(it)
         end do
      end do

   end subroutine torpol_to_dphspat
!-------------------------------------------------------------------------------
   subroutine pol_to_curlr_spat(Qlm, cvrc)
      complex(cp), intent(in) :: Qlm(lm_max)
      real(cp), intent(out) :: cvrc(n_phi_max, n_theta_max)

      ! local
      complex(cp) :: dQlm(lm_max)
      integer :: lm


      do lm = 1, lm_max
         dQlm(lm) = dLh(lm) * Qlm(lm)
      end do

      call shtns_SH_to_spat(dQlm, cvrc)

   end subroutine pol_to_curlr_spat
!-------------------------------------------------------------------------------
   subroutine torpol_to_curl_spat(Blm, ddBlm, Jlm, dJlm, nR, &
                                 cvrc, cvtc, cvpc)
      complex(cp), intent(in) :: Blm(lm_max), ddBlm(lm_max)
      complex(cp), intent(in) :: Jlm(lm_max), dJlm(lm_max)
      integer, intent(in) :: nR
      real(cp), intent(out) :: cvrc(n_phi_max, n_theta_max)
      real(cp), intent(out) :: cvtc(n_phi_max, n_theta_max)
      real(cp), intent(out) :: cvpc(n_phi_max, n_theta_max)

      ! local
      complex(cp) :: Qlm(lm_max), Tlm(lm_max)
      integer :: lm

      do lm = 1, lm_max
         Qlm(lm) = dLh(lm) * Jlm(lm)
         Tlm(lm) = 1/r(nR)**2 * dLh(lm) * Blm(lm) - ddBlm(lm)
      end do

      call shtns_qst_to_spat(Qlm, dJlm, Tlm, cvrc, cvtc, cvpc)

   end subroutine torpol_to_curl_spat
!-------------------------------------------------------------------------------
   subroutine spat_to_SH(f, fLM)

      real(cp), intent(in) :: f(n_phi_max, n_theta_max)
      complex(cp), intent(out) :: fLM(lm_max)  ! THIS IS WRONG Should be lmP_max! LAGOTEST (Thomas already fixed this in master branch)

      call shtns_load_cfg(1)
      call shtns_spat_to_sh(f, fLM)
      call shtns_load_cfg(0)

   end subroutine spat_to_SH
   
!-------------------------------------------------------------------------------
! 
! Distributed version:
!@>TODO Update the details and notation of these functions to properly describe 
!> what their purpose
!
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
   subroutine spat_to_SH_dist(Vr_loc, QlmP_loc)
!@>details transform the spherical harmonic coefficients Qlm into its spatial 
!>representation Vr.
!> 
!@>author Rafael Lago (MPCDF) July 2017
!-------------------------------------------------------------------------------
      real(cp),     intent(inout) :: Vr_loc(n_phi_max,n_theta_loc)
      complex(cp),  intent(out)   :: QlmP_loc(lmP_loc)
      
      complex(cp) ::  QlP_loc(n_theta_max,n_m_loc)
      complex(cp) ::  transpose_loc(n_m_max,n_theta_loc)
      complex(cp) ::  Fr_loc(n_phi_max/2+1,n_theta_loc)
      
      integer :: m_idx, lm_s, lm_e, i
      
      call shtns_load_cfg(1) ! l_max + 1
      
      !@>TODO: The FFT must be performed for an array with the dimensions of 
      !> Fr_loc which may end up paded with zeroes.
      !> Is there any way to tell MKL to perform a "truncated" FFT?
      !@>TODO: Terrible performance here!
      call fft_phi_loc(Vr_loc, Fr_loc, 1)
      transpose_loc(1:n_m_max,1:n_theta_loc) = Fr_loc(1:n_m_max,1:n_theta_loc)
      
      call transpose_m_theta(transpose_loc, QlP_loc)
      
      ! Now do the Legendre transform using the new function in a loop
      do i = 1, n_m_loc
        m_idx = lmP_dist(coord_theta, i, 1)/minc
        lm_s  = lmP_dist(coord_theta, i, 3)
        lm_e  = lmP_dist(coord_theta, i, 4)
        call shtns_spat_to_sh_ml(m_idx,QlP_loc(:,i),QlmP_loc(lm_s:lm_e),l_max+1)
      end do
      
      call shtns_load_cfg(0) ! l_max

   end subroutine spat_to_SH_dist
!-------------------------------------------------------------------------------
   subroutine sh_to_spat_dist(Qlm_loc, Vr_loc)
!@>details transform the spherical harmonic coefficients Qlm into its spatial 
!> representation Vr
!
!@>author Rafael Lago (MPCDF) December 2017
!-------------------------------------------------------------------------------
      complex(cp),  intent(inout)   :: Qlm_loc(lm_loc)
      real(cp),     intent(out)     :: Vr_loc(n_phi_max, n_theta_loc)
      
      complex(cp) :: Ql_loc(n_theta_max,n_m_loc)
      complex(cp) :: transposed_loc(n_m_max,n_theta_loc)
      complex(cp) :: F_loc(n_phi_max/2+1,n_theta_loc)
      
      integer :: i, lm_s, lm_e, m_idx
      
      do i = 1, n_m_loc
        m_idx = lm_dist(coord_theta, i, 1)/minc
        lm_s  = lm_dist(coord_theta, i, 3)
        lm_e  = lm_dist(coord_theta, i, 4)
        call shtns_sh_to_spat_ml(m_idx, Qlm_loc(lm_s:lm_e), Ql_loc(:,i), l_max)
      end do
      
      call transpose_theta_m(Ql_loc, transposed_loc)
      
      !@>TODO: The FFT must be performed for an array with the dimensions of 
      !> F_loc which may end up paded with zeroes.
      !> Is there any way to tell MKL to perform a "truncated" FFT?
      !@>TODO: Terrible performance here!
      F_loc = 0.0
      F_loc(1:n_m_max,1:n_theta_loc) = transposed_loc
      
      call fft_phi_loc(Vr_loc, F_loc, -1)
      
   end subroutine sh_to_spat_dist
!-------------------------------------------------------------------------------
   subroutine sph_to_spat_dist(Slm_loc, gVt_loc, gVp_loc)
!@>details transform spheroidal spherical harmonic coefficients Slm to the 
!> spatial theta and phi components (Vt,Vp), effectively computing the 
!> gradient of S. This is the distributed counterpart of [scal_to_grad_spat]
!
!@>author Rafael Lago (MPCDF) November 2017
!-------------------------------------------------------------------------------
      complex(cp),  intent(inout)   :: Slm_loc(lm_loc)
      real(cp),     intent(out)     :: gVt_loc(n_phi_max, n_theta_loc)
      real(cp),     intent(out)     :: gVp_loc(n_phi_max, n_theta_loc)
      
      complex(cp) :: Sl_t(n_theta_max,n_m_loc), Sl_p(n_theta_max,n_m_loc)
      complex(cp) :: tranposed_loc(n_m_max,n_theta_loc)
      complex(cp) :: F_loc(n_phi_max/2+1,n_theta_loc)
      
      integer :: i, lm_s, lm_e, m_idx
      
      
      do i = 1, n_m_loc
        m_idx = lm_dist(coord_theta, i, 1)/minc
        lm_s  = lm_dist(coord_theta, i, 3)
        lm_e  = lm_dist(coord_theta, i, 4)
        call shtns_sph_to_spat_ml(m_idx, Slm_loc(lm_s:lm_e), &
                                         Sl_t(:,i), Sl_p(:,i), l_max)
      end do
      
      !>@TODO write the rest of this function such that both the poloidal and 
      !> the toroidal part can be transposed and Fourier-transformed simultanouesly
      call transpose_theta_m(Sl_t, tranposed_loc)
      F_loc = 0.0
      F_loc(1:n_m_max,1:n_theta_loc) = tranposed_loc
      call fft_phi_loc(gVt_loc, F_loc, -1)
      !--------
      call transpose_theta_m(Sl_p, tranposed_loc)
      F_loc = 0.0
      F_loc(1:n_m_max,1:n_theta_loc) = tranposed_loc
      call fft_phi_loc(gVp_loc, F_loc, -1)
            
   end subroutine
!-------------------------------------------------------------------------------
   subroutine qst_to_spat_dist(Qlm_loc, Slm_loc, Tlm_loc, &
                                  Vr_loc, Vt_loc, Vp_loc)
!@>details 3D vector transform from spherical coordinates to radial-spheroidal-
!> toroidal spectral components. 
!
!@>author Rafael Lago (MPCDF) December 2017
!-------------------------------------------------------------------------------
      complex(cp),  intent(inout) :: Qlm_loc(lm_loc)
      complex(cp),  intent(inout) :: Slm_loc(lm_loc)
      complex(cp),  intent(inout) :: Tlm_loc(lm_loc)
      real(cp),     intent(out)   :: Vr_loc(n_phi_max, n_theta_loc)
      real(cp),     intent(out)   :: Vt_loc(n_phi_max, n_theta_loc)
      real(cp),     intent(out)   :: Vp_loc(n_phi_max, n_theta_loc)
      
      complex(cp) :: Ql_loc(n_theta_max,n_m_loc)
      complex(cp) :: Sl_loc(n_theta_max,n_m_loc)
      complex(cp) :: Tl_loc(n_theta_max,n_m_loc)
      complex(cp) :: transposed_loc(n_m_max,n_theta_loc)
      complex(cp) :: F_loc(n_phi_max/2+1,n_theta_loc)
      
      integer :: i, lm_s, lm_e, m_idx
      
      do i = 1, n_m_loc
        m_idx = lm_dist(coord_theta, i, 1)/minc
        lm_s  = lm_dist(coord_theta, i, 3)
        lm_e  = lm_dist(coord_theta, i, 4)
        call shtns_qst_to_spat_ml(m_idx, Qlm_loc(lm_s:lm_e), Slm_loc(lm_s:lm_e), &
                                  Tlm_loc(lm_s:lm_e), Ql_loc(:,i), &
                                  Sl_loc(:,i), Tl_loc(:,i), l_max)
      end do
      
      !>@TODO write the rest of this function such that both the poloidal and 
      !> the toroidal part can be transposed and Fourier-transformed simultanouesly
      call transpose_theta_m(Ql_loc, transposed_loc)
      F_loc = 0.0
      F_loc(1:n_m_max,1:n_theta_loc) = transposed_loc
      call fft_phi_loc(Vr_loc, F_loc, -1)
      !--------
      call transpose_theta_m(Sl_loc, transposed_loc)
      F_loc = 0.0
      F_loc(1:n_m_max,1:n_theta_loc) = transposed_loc
      call fft_phi_loc(Vt_loc, F_loc, -1)
      !--------
      call transpose_theta_m(Tl_loc, transposed_loc)
      F_loc = 0.0
      F_loc(1:n_m_max,1:n_theta_loc) = transposed_loc
      call fft_phi_loc(Vp_loc, F_loc, -1)
      
   end subroutine qst_to_spat_dist
!-------------------------------------------------------------------------------
   subroutine pol_to_curlr_spat_dist(Plm_loc, cVr_loc)
!@>details computes the spherical harmonic coefficients Qlm throught the 
!> poloidal scalar Plm, and then transforms it into the curl of its spatial 
!> representation cVr
!> 
!@>author Rafael Lago (MPCDF) August 2017
!-------------------------------------------------------------------------------
      complex(cp),  intent(inout)   :: Plm_loc(lm_loc)
      real(cp),     intent(out)     :: cVr_loc(n_phi_max, n_theta_loc)
      
      complex(cp) :: Qlm_loc(lm_loc)
      complex(cp) :: Ql_loc(n_theta_max,n_m_loc)
      complex(cp) :: transpose_loc(n_m_max,n_theta_loc)
      complex(cp) :: F_loc(n_phi_max/2+1,n_theta_loc)
      
      integer :: i, lm_s, lm_e, m_idx
      
      Qlm_loc(1:lm_loc) = dLh_loc(1:lm_loc) * Plm_loc(1:lm_loc)
      
      call sh_to_spat_dist(Qlm_loc,cVr_loc)
      
   end subroutine pol_to_curlr_spat_dist
!-------------------------------------------------------------------------------
   subroutine pol_to_grad_spat_dist(Plm_loc, gVt_loc, gVp_loc)
!@>details Computes Qlm (spheroidal spherical harmonic coefficients) from Plm
!> (poloidal scalar) and then transforms them into the spatial theta and phi 
!> components (Vt,Vp), effectively computing the gradient of S
!> 
!@>author Rafael Lago (MPCDF) November 2017
!-------------------------------------------------------------------------------
      complex(cp),  intent(inout)   :: Plm_loc(lm_max)
      real(cp),     intent(out)     :: gVt_loc(n_phi_max, n_theta_max)
      real(cp),     intent(out)     :: gVp_loc(n_phi_max, n_theta_max)
      
      complex(cp) :: Qlm_loc(lm_loc)
      
      Qlm_loc(1:lm_loc) = dLh_loc(1:lm_loc) * Plm_loc(1:lm_loc)
      
      call sph_to_spat_dist(Qlm_loc, gVt_loc, gVp_loc)
            
   end subroutine pol_to_grad_spat_dist
!-------------------------------------------------------------------------------
   subroutine torpol_to_spat_dist(Plm_loc, Slm_loc, Tlm_loc, &
                                  Vr_loc, Vt_loc, Vp_loc)
!@>details computes the spherical harmonic coefficients Qlm throught the 
!> poloidal scalar Plm, then transforms the 3D vector from spherical 
!> coordinates to radial-spheroidal-toroidal spectral components. 
!
!@>author Rafael Lago (MPCDF) December 2017
!-------------------------------------------------------------------------------
      complex(cp),  intent(inout) :: Plm_loc(lm_loc)
      complex(cp),  intent(inout) :: Slm_loc(lm_loc)
      complex(cp),  intent(inout) :: Tlm_loc(lm_loc)
      real(cp),     intent(out)   :: Vr_loc(n_phi_max, n_theta_loc)
      real(cp),     intent(out)   :: Vt_loc(n_phi_max, n_theta_loc)
      real(cp),     intent(out)   :: Vp_loc(n_phi_max, n_theta_loc)
      
      complex(cp) :: Qlm_loc(lm_loc)
      
      Qlm_loc(1:lm_loc) = dLH_loc(1:lm_loc)*Plm_loc(1:lm_loc)
      
      call qst_to_spat_dist(Qlm_loc,Slm_loc,Tlm_loc,Vr_loc,Vt_loc,Vp_loc)
      
   end subroutine torpol_to_spat_dist
!-------------------------------------------------------------------------------
   subroutine torpol_to_dphspat_dist(dWlm_loc, Zlm_loc, dVt_loc, dVp_loc)
!@>details transform spheroidal-toroidal spherical harmonic coefficients &
!> (Slm,Tlm) to the derivative of the spatial theta and phi components (Vt,Vp).
!> 
!@>author Rafael Lago (MPCDF) November 2017
!-------------------------------------------------------------------------------
      complex(cp),  intent(inout)   :: dWlm_loc(lm_loc)
      complex(cp),  intent(inout)   :: Zlm_loc(lm_loc)
      real(cp),     intent(out)     :: dVt_loc(n_phi_max, n_theta_max)
      real(cp),     intent(out)     :: dVp_loc(n_phi_max, n_theta_max)
      
      complex(cp) :: Slm_loc(lm_loc)
      complex(cp) :: Tlm_loc(lm_loc)
      complex(cp) :: Sl_loc(n_theta_max,n_m_loc)
      complex(cp) :: Tl_loc(n_theta_max,n_m_loc)
      complex(cp) :: transpose_loc(n_m_max,n_theta_loc)
      complex(cp) :: F_loc(n_phi_max/2+1,n_theta_loc)
      
      integer  :: i, j, lm_s, lm_e, m_idx
      complex(cp) :: m
      
      !>@TODO do we need those temporary arrays?
      ! This is D_m_loc is just m_idx...
      do i = 1, lm_loc
         m = ci*D_m_loc(i)
         Slm_loc(i) = m*dWlm_loc(i)
         Tlm_loc(i) = m* Zlm_loc(i)
      end do
      
      do i = 1, n_m_loc
        m_idx = lm_dist(coord_theta, i, 1)/minc
        lm_s  = lm_dist(coord_theta, i, 3)
        lm_e  = lm_dist(coord_theta, i, 4)
        call shtns_sphtor_to_spat_ml(m_idx, Slm_loc(lm_s:lm_e), Tlm_loc(lm_s:lm_e), &
                                     Sl_loc(:,i), Tl_loc(:,i), l_max)
      end do
      
      !>@TODO write the rest of this function such that both the poloidal and 
      !> the toroidal part can be transposed and Fourier-transformed simultanouesly
      call transpose_theta_m(Sl_loc, transpose_loc)
      F_loc = 0.0
      F_loc(1:n_m_max,1:n_theta_loc) = transpose_loc
      call fft_phi_loc(dVt_loc, F_loc, -1)
      do j=1, n_theta_loc
         do i=1, n_phi_max
            dVt_loc(i, j) = dVt_loc(i, j) * O_sin_theta_E2(n_theta_beg+j-1)
         end do
      end do
      !--------
      call transpose_theta_m(Tl_loc, transpose_loc)
      F_loc = 0.0
      F_loc(1:n_m_max,1:n_theta_loc) = transpose_loc
      call fft_phi_loc(dVp_loc, F_loc, -1)
      
      do j=1, n_theta_loc
         do i=1, n_phi_max
            dVp_loc(i, j) = dVp_loc(i, j) * O_sin_theta_E2(n_theta_beg+j-1)
         end do
      end do
            
   end subroutine torpol_to_dphspat_dist
!-------------------------------------------------------------------------------
   subroutine torpol_to_curl_spat_dist(Blm_loc, ddBlm_loc, Jlm_loc, Slm_loc, &
                                       nR, cVr_loc, cVt_loc, cVp_loc)
!@>details 3D vector transform from spherical coordinates to the curl of the 
!> radial-spheroidal-toroidal spectral components. 
!> 
!@>TODO: Review the notation here
!> 
!@>author Rafael Lago (MPCDF) November 2017
!-------------------------------------------------------------------------------
      complex(cp),  intent(inout) :: Blm_loc(lm_loc)
      complex(cp),  intent(inout) :: ddBlm_loc(lm_loc)
      complex(cp),  intent(inout) :: Jlm_loc(lm_loc)
      complex(cp),  intent(inout) :: Slm_loc(lm_loc)
      integer,      intent(in)    :: nR
      real(cp),     intent(out)   :: cVr_loc(n_phi_max, n_theta_max)
      real(cp),     intent(out)   :: cVt_loc(n_phi_max, n_theta_max)
      real(cp),     intent(out)   :: cVp_loc(n_phi_max, n_theta_max)
      
      complex(cp) :: Qlm_loc(lm_loc)
      complex(cp) :: Tlm_loc(lm_loc)
      
      Qlm_loc(1:lm_loc) = dLH_loc(1:lm_loc)*Jlm_loc(1:lm_loc)
      Tlm_loc(1:lm_loc) = 1/r(nR)**2 * dLH_loc(1:lm_loc) * Blm_loc(1:lm_loc) &
                        - ddBlm_loc(1:lm_loc) 
      
      call qst_to_spat_dist(Qlm_loc,Slm_loc,Tlm_loc,cVr_loc,cVt_loc,cVp_loc)
      
   end subroutine torpol_to_curl_spat_dist

end module shtns

