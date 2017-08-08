module shtns

   use precision_mod, only: cp
   use constants, only: ci
   use truncation, only: m_max, l_max, n_theta_max, n_phi_max, &
                         minc, lm_max, n_m_max, nrp, lmP_max
   use horizontal_data, only: dLh, gauss, theta_ord, D_m, O_sin_theta_E2
   use radial_functions, only: r
   use parallel_mod
   use fft, only: init_fft_phi, finalize_fft_phi, fft_phi;
   use blocking, only: lm2, lmP2

   implicit none

   include "shtns.f"

   private

   public :: init_shtns, scal_to_spat, scal_to_grad_spat, pol_to_grad_spat, &
             torpol_to_spat, pol_to_curlr_spat, torpol_to_curl_spat,        &
             torpol_to_dphspat, spat_to_SH, spat_to_SH_parallel
   public :: finalize_shtns

contains

   subroutine init_shtns()

      integer :: nlm
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
!------------------------------------------------------------------------------
   subroutine finalize_shtns
!@>author Rafael Lago, MPCDF, July 2017
!------------------------------------------------------------------------------

      call finalize_fft_phi

   end subroutine finalize_shtns
!------------------------------------------------------------------------------
   subroutine scal_to_spat(Slm, fieldc)
      ! transform a spherical harmonic field into grid space
      complex(cp), intent(in) :: Slm(lm_max)
      real(cp), intent(out) :: fieldc(n_phi_max, n_theta_max)

      call shtns_SH_to_spat(Slm, fieldc)

   end subroutine scal_to_spat
!------------------------------------------------------------------------------
   subroutine scal_to_grad_spat(Slm, gradtc, gradpc)
      ! transform a scalar spherical harmonic field into it's gradient
      ! on the grid
      complex(cp), intent(in) :: Slm(lm_max)
      real(cp), intent(out) :: gradtc(n_phi_max, n_theta_max)
      real(cp), intent(out) :: gradpc(n_phi_max, n_theta_max)

      call shtns_sph_to_spat(Slm, gradtc, gradpc)

   end subroutine scal_to_grad_spat
!------------------------------------------------------------------------------
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
!------------------------------------------------------------------------------
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
!------------------------------------------------------------------------------
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
!------------------------------------------------------------------------------
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
!------------------------------------------------------------------------------
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
!------------------------------------------------------------------------------
   subroutine spat_to_SH(f, fLM)

      real(cp), intent(in) :: f(n_phi_max, n_theta_max)
      complex(cp), intent(out) :: fLM(lm_max)

      call shtns_load_cfg(1)
      call shtns_spat_to_sh(f, fLM)
      call shtns_load_cfg(0)

   end subroutine spat_to_SH

!------------------------------------------------------------------------------
   subroutine spat_to_SH_parallel(f, fLM, name, passed)
!@>details Computes the FFT and Legendre transform for distributed θ and m.
!> It will take the f in (φ,θ) space, and then apply FFT nθ times, to obtain 
!> \hat f in (m,θ) space. 
!> Following, it will redistribute the data such that every m has access to all 
!> of its respective θ points - obtaining \bar f in (θ,m)
!> Then it will loop over m, applying the Legendre transform to obtain the 
!> field fLM in (l,m) coordinates.
!> 
!> This uses SHtns.
!> 
!
!>@bug if you pass f as a real argument converted to complex, e.g.
!> call spat_to_SH_ml(1, cmplx(f,kind=cp), fLM, 10)
!> then C will do something silly with the memory addresses and overwrite
!> m_idx (which is, theoretically, read-only). These on-the-fly conversion 
!> are discouraged anyway.
!
!@>author Rafael Lago (MPCDF) July 2017
!------------------------------------------------------------------------------

      real(cp),     intent(inout) :: f(n_phi_max, n_theta_max)
      complex(cp),  intent(inout) :: fLM(lmP_max)
      character(len=*), intent(in) :: name
      logical :: passed
      
      complex(cp) :: f_phi_theta(n_phi_max,n_theta_max)
      complex(cp) :: f_m_theta(m_max+1,n_theta_max)
      complex(cp) :: f_theta_m(n_theta_max, m_max+1)
      complex(cp) :: new_fLM(lmP_max)
      integer :: m_idx, lm_s, lm_e
      real    :: norm_diff, norm_fLM, norm_new, cosfLM
      
      
      call shtns_load_cfg(1)
      f_phi_theta = cmplx(f, 0.0_cp, kind=cp)
      f_theta_m = 0.0_cp
      f_m_theta = 0.0_cp
      new_fLM = 0.0_cp
      
      call shtns_spat_to_sh(f, fLM)
      
      call fft_phi(f_phi_theta, f_m_theta, 1)
      f_theta_m = transpose(f_m_theta)
      
      ! Now do it using the new function in a loop
      do m_idx = 0, m_max, minc
        lm_s = lmP2(m_idx  ,m_idx)
        lm_e = lmP2(l_max+1,m_idx)
        call shtns_spat_to_sh_ml(m_idx, f_theta_m(:,m_idx+1), new_fLM(lm_s:lm_e), l_max+1) !l_max - m_idx + 1)
      end do
      
      
      ! These norms, cosine and prints are only for the debugging version.
      norm_fLM = NORM2(real(fLM))
      norm_new = NORM2(real(new_fLM))
      if (norm_new == 0.0_cp) norm_new = 1.0_cp
      if (norm_fLM == 0.0_cp) norm_fLM = 1.0_cp
      norm_diff = NORM2(real(new_fLM) - real(fLM))/real(lmP_max)
      cosfLM = dot_product(real(new_fLM),real(fLM))/(norm_fLM*norm_new)
      
      if (norm_diff > 1e-10) then
        print *, "Conversion for  ", name, " doesn't match! Norm: ", norm_diff, " cos:", cosfLM
        print *, fLM(1:10)
        print *, "~~~~~~~"
        print *, new_fLM(1:10)
        passed = .false.
      else
        if (rank == 0 ) then
          print *, name, "'s Norm:", norm_diff, " cos:", cosfLM
        end if
        passed = .true.
      end if
      
      call shtns_load_cfg(0)

   end subroutine spat_to_SH_parallel
!------------------------------------------------------------------------------
end module shtns
