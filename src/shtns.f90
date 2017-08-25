module shtns

   use precision_mod, only: cp
   use constants, only: ci
   use truncation
   use horizontal_data, only: dLh, theta_ord, D_m, O_sin_theta_E2
   use radial_functions, only: or2
   use parallel_mod
   use fft, only: init_fft_phi, finalize_fft_phi, fft_phi, fft_phi_dist
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
<<<<<<< HEAD

=======
>>>>>>> spat_to_SH_parallel function implemented and working (tested with dynamo_bench)
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
           &                1.e-10_cp, n_theta_max, n_phi_max)
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
              &                  cvrc, cvtc, cvpc)
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
         Tlm(lm) = or2(nR) * dLh(lm) * Blm(lm) - ddBlm(lm)
      end do

      call shtns_qst_to_spat(Qlm, dJlm, Tlm, cvrc, cvtc, cvpc)

   end subroutine torpol_to_curl_spat
!------------------------------------------------------------------------------
   subroutine spat_to_SH(f, fLM)

      real(cp), intent(in) :: f(n_phi_max, n_theta_max)
      complex(cp), intent(out) :: fLM(lmP_max)

      call shtns_load_cfg(1)
      call shtns_spat_to_sh(f, fLM)
      call shtns_load_cfg(0)

   end subroutine spat_to_SH

!------------------------------------------------------------------------------
   subroutine spat_to_SH_parallel(f, fLM)
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
      complex(cp),  intent(out)   :: fLM(lmP_max)
      
      complex(cp) :: f_phi_theta_loc(n_phi_max,n_theta_loc)
      complex(cp) :: f_m_theta_loc(m_max+1,n_theta_loc)
      complex(cp) :: f_theta_m_loc(n_theta_max, n_m_loc)
      
      complex(cp) :: dist_fLM(lmP_loc)
      integer :: m_idx, lm_s, lm_e, i, ierr
      
      
      call shtns_load_cfg(1)
      f_phi_theta_loc = cmplx(f(:,n_theta_beg:n_theta_end), 0.0_cp, kind=cp)
      
      !>@TODO have f_phi_theta_loc as REAL, and use REAL2CMLX transform from mkl
      call fft_phi_dist(f_phi_theta_loc, f_m_theta_loc, 1)
      call transpose_m_theta(f_m_theta_loc,f_theta_m_loc)
      
      ! Now do it using the new function in a loop
      do i = 1, n_m_loc
        m_idx = lmP_dist(coord_theta, i, 1)
        lm_s  = lmP_dist(coord_theta, i, 3)
        lm_e  = lmP_dist(coord_theta, i, 4)
        call shtns_spat_to_sh_ml(m_idx, f_theta_m_loc(:,i), dist_fLM(lm_s:lm_e), l_max+1)
      end do
      
      call gather_fLM(dist_fLM, fLM)
      
      call shtns_load_cfg(0)

   end subroutine spat_to_SH_parallel
!------------------------------------------------------------------------------
<<<<<<< HEAD
   subroutine spat_to_sphertor(f,g,fLM,gLM)

      real(cp), intent(in) :: f(n_phi_max,n_theta_max)
      real(cp), intent(in) :: g(n_phi_max,n_theta_max)
      complex(cp), intent(out) :: fLM(lmP_max)
      complex(cp), intent(out) :: gLM(lmP_max)

      call shtns_load_cfg(1)
      call shtns_spat_to_sphtor(f,g,fLM,gLM)
      call shtns_load_cfg(0)

   end subroutine spat_to_sphertor
!------------------------------------------------------------------------------
=======
   subroutine gather_fLM(fLM_local, fLM_global)
!@>details Gathers the fLM which was computed in a distributed fashion.
!> Mostly used for debugging. This function may not be performant at all.
!
!@>TODO this with mpi_allgatherv
!@>author Rafael Lago (MPCDF) August 2017
!------------------------------------------------------------------------------
      complex(cp),  intent(inout) :: fLM_local(lmP_loc)
      complex(cp),  intent(out)   :: fLM_global(lmP_max)
      
      complex(cp) ::  buffer(lmP_max)
      integer :: i, j, m_idx, lm_s_local, lm_e_local, lm_s_global, lm_e_global
      integer :: pos, ilen, Rq(n_procs_theta), ierr
      
      ! buffer will receive all messages, but they are ordered by ranks,
      ! not by m.
      
!       buffer = 0.0
      pos = 1
      do i=0,n_procs_theta-1
         ilen = sum(lmP_dist(i,:,2))
         if (coord_theta == i) buffer(pos:pos+ilen-1) = fLM_local(1:lmP_loc)
         CALL MPI_IBCAST(buffer(pos:pos+ilen-1), ilen, MPI_DOUBLE_COMPLEX, i, comm_theta, Rq(i+1), ierr)
         pos = pos + ilen
      end do
      
      CALL MPI_WAITALL(n_procs_theta, Rq, MPI_STATUSES_IGNORE, ierr)
      
      ! This basically re-orders the buffer 
      pos = 0
      do i=0,n_procs_theta-1
         do j = 1, n_m_ext
            m_idx = lmP_dist(i, j, 1)
            if (m_idx < 0) exit
            lm_s_local  = pos + lmP_dist(i, j, 3)
            lm_e_local  = pos + lmP_dist(i, j, 4)
            lm_s_global = lmP2(m_idx  ,m_idx)
            lm_e_global = lmP2(l_max+1,m_idx)
            fLM_global(lm_s_global:lm_e_global) = buffer(lm_s_local:lm_e_local)
         end do
         pos = pos + sum(lmP_dist(i,:,2))
      end do
      
   end subroutine gather_fLM
!------------------------------------------------------------------------------
   subroutine transpose_m_theta(f_m_theta, f_theta_m)
!@>details Does the transposition using alltoallv.
!
!@>TODO this with mpi_type to stride the data
!@>author Rafael Lago (MPCDF) August 2017
!------------------------------------------------------------------------------
      complex(cp), intent(inout) :: f_m_theta(m_max+1,n_theta_loc)
      complex(cp), intent(inout) :: f_theta_m(n_theta_max, n_m_loc)
      
      complex(cp) :: sendbuf((m_max+1)*n_theta_loc)
      complex(cp) :: recvbuf(n_m_loc, n_theta_max)
      
      integer :: sendcount(0:n_procs_theta-1)
      integer :: recvcount(0:n_procs_theta-1)
      integer :: senddispl(0:n_procs_theta-1)
      integer :: recvdispl(0:n_procs_theta-1)
      integer :: i, j, k, m_idx, pos, n_theta
      
      pos = 1
      do i=0,n_procs_theta-1
         ! Copy each m which belongs to the i-th rank into the send buffer
         ! column-wise. That will simplify a lot things later
         !
         !>@TODO check performance of this; implementing this with mpi_type
         !> striding the data will probably be faster
         senddispl(i) = pos-1
         do k=1,n_theta_loc
            do j=1,n_m_ext
               if (lmP_dist(i,j,1) < 0) exit
               m_idx = lmP_dist(i,j,1)
               sendbuf(pos) = f_m_theta(m_idx+1,k)
               pos = pos + 1
            end do
         end do
         sendcount(i) = pos - senddispl(i) - 1
         n_theta = n_theta_dist(i,2) - n_theta_dist(i,1) + 1
         recvdispl(i) = i*n_m_loc*n_theta
         recvcount(i) = n_m_loc*n_theta
      end do
      
      call MPI_ALLTOALLV(sendbuf, sendcount, senddispl, MPI_DOUBLE_COMPLEX, &
                         recvbuf, recvcount, recvdispl, MPI_DOUBLE_COMPLEX, &
                         comm_theta, i)
      f_theta_m = transpose(recvbuf)
      
   end subroutine transpose_m_theta
>>>>>>> spat_to_SH_parallel function implemented and working (tested with dynamo_bench)
end module shtns

