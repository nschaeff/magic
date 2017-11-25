module shtns

   use precision_mod, only: cp
   use constants, only: ci
   use truncation
   use horizontal_data, only: dLh, gauss, theta_ord, D_m, O_sin_theta_E2
   use radial_functions, only: r
   use parallel_mod
   use fft, only: init_fft_phi, finalize_fft_phi, fft_phi_loc
   use blocking, only: lm2, lmP2
   
   implicit none

   include "shtns.f"

   private

   public :: init_shtns, scal_to_spat, scal_to_grad_spat, pol_to_grad_spat, &
             torpol_to_spat, pol_to_curlr_spat, torpol_to_curl_spat,        &
             torpol_to_dphspat, spat_to_SH,&
             spat_to_SH_parallel, scal_to_spat_parallel
   public :: finalize_shtns

contains

   subroutine init_shtns()
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
   subroutine spat_to_SH_parallel(f, fLM)
!@>details Computes the FFT and Legendre transform for distributed θ and m.
!> It will take the f(φ,θ) and apply FFT nθ times, to obtain g(m,θ). 
!> Following, it will redistribute the data such that every m has access to all 
!> of its respective θ points - obtaining \bar gT(θ,m)
!> Then it will loop over m, applying the Legendre transform to obtain the 
!> field fLM in (l,m) coordinates.
!> 
!> This uses SHtns.
!>
!> Notation: 
!> ---------
!>       f: function in the real domain
!>       g: f in Fourier domain
!>      gT: the transpose of g
!>     fLM: f in LM-domain. This is the same as gT after applying the Legendre 
!>          polynomials
!>    _loc: local version of the supracited variables. In the LM-domain, the 
!>    ms are distributed amonst n_procs_theta MPI ranks. In the other domains, 
!>    the θs are distributed amonsts the n_procs_theta MPI ranks.
!
!@>author Rafael Lago (MPCDF) July 2017
!------------------------------------------------------------------------------
      real(cp),     intent(inout) :: f(n_phi_max, n_theta_max)
      complex(cp),  intent(out)   :: fLM(lmP_max)
      
      real(cp) :: f_loc(n_phi_max,n_theta_loc)
      complex(cp) :: gT_loc(n_theta_max,n_m_loc)
      complex(cp) ::  g_loc(m_max+1,n_theta_loc)
      complex(cp) ::  tmp_loc(n_phi_max/2+1,n_theta_loc)
      
      complex(cp) :: fLM_loc(lmP_loc)
      integer :: m_idx, lm_s, lm_e, i
      
      
      call shtns_load_cfg(1)
      
      ! Copies global array's data into local array
      !>@TODO receive f_loc as argument instead of f
      f_loc = f(:,n_theta_beg:n_theta_end)
      
      !@>TODO: The FFT must be performed for an array with the dimensions of tmp_loc
      !> which may end up paded with zeroes.
      !> Is there any way to tell MKL to perform a "truncated" FFT?
      !@>TODO: Terrible performance here!
      call fft_phi_loc(f_loc, tmp_loc, 1)
      g_loc = tmp_loc(1:m_max+1,1:n_theta_loc)
      
      call transpose_m_theta(g_loc, gT_loc)
      
      ! Now do the Legendre transform using the new function in a loop
      do i = 1, n_m_loc
        m_idx = lmP_dist(coord_theta, i, 1)
        lm_s  = lmP_dist(coord_theta, i, 3)
        lm_e  = lmP_dist(coord_theta, i, 4)
        call shtns_spat_to_sh_ml(m_idx, gT_loc(:,i), fLM_loc(lm_s:lm_e), l_max+1)
      end do
      
      ! Gathers the transform computed in every process in this process
      !>@TODO receive return fLM_loc instead of fLM
      call gather_fLM(fLM_loc, fLM)
      
      call shtns_load_cfg(0)

   end subroutine spat_to_SH_parallel
!------------------------------------------------------------------------------
   subroutine scal_to_spat_parallel(fLM, f)
!@>details Does the reverse of the routine above.
!> 
!> Notice that here we use [lm_dist] and NOT lmP_dist. We assume that there 
!> are (m_max+1)*(l_max) modes.
!> 
!@>author Rafael Lago (MPCDF) August 2017
!------------------------------------------------------------------------------
      complex(cp),  intent(inout)   :: fLM(lm_max)
      real(cp),     intent(out)     :: f(n_phi_max, n_theta_max)
      
      complex(cp) :: fLM_loc(lm_loc)
      complex(cp) :: gT_loc(n_theta_max,n_m_loc)
      complex(cp) ::  g_loc(m_max+1,n_theta_loc)
      complex(cp) ::  tmp_loc(n_phi_max/2+1,n_theta_loc)
      real(cp)    ::  f_loc(n_phi_max, n_theta_loc)
      
      integer :: i, lm_s, lm_e, lm_gs, lm_ge, m_idx
      
      ! Copies fLM into fLM_loc
      !>@TODO receive fLM_loc as argument instead of fLM
      do i = 1, n_m_loc
        m_idx = lm_dist(coord_theta, i, 1)
        lm_s  = lm_dist(coord_theta, i, 3)
        lm_e  = lm_dist(coord_theta, i, 4)
        call shtns_lmidx(lm_gs, m_idx, m_idx)
        call shtns_lmidx(lm_ge, l_max, m_idx)
        fLM_loc(lm_s:lm_e) = fLM(lm_gs:lm_ge)
      end do
      
      do i = 1, n_m_loc
        m_idx = lm_dist(coord_theta, i, 1)
        lm_s  = lm_dist(coord_theta, i, 3)
        lm_e  = lm_dist(coord_theta, i, 4)
        call shtns_sh_to_spat_ml(m_idx, fLM_loc(lm_s:lm_e), gT_loc(:,i), l_max)
      end do
      
      call transpose_theta_m(gT_loc, g_loc)
      
      !@>TODO: The FFT must be performed for an array with the dimensions of tmp_loc
      !> which may end up paded with zeroes.
      !> Is there any way to tell MKL to perform a "truncated" FFT?
      !@>TODO: Terrible performance here!
      tmp_loc = 0.0
      tmp_loc(1:m_max+1,1:n_theta_loc) = g_loc
      
      call fft_phi_loc(f_loc, tmp_loc, -1)
      
      call gather_f(f_loc, f)
            
   end subroutine

!------------------------------------------------------------------------------
   subroutine gather_fLM(fLM_local, fLM_global)
!@>details Gathers the fLM which was computed in a distributed fashion.
!> Mostly used for debugging. This function may not be performant at all.
!> 
!> Notice that here we use [lmP_dist] and NOT lm_dist. We assume that there 
!> are (m_max+1)*(l_max+1) modes.
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
   subroutine gather_f(f_local, f_global)
!@>details Gathers the fLM which was computed in a distributed fashion.
!> Mostly used for debugging. This function may not be performant at all.
!
!@>TODO this with mpi_allgatherv
!@>author Rafael Lago (MPCDF) August 2017
!------------------------------------------------------------------------------
      real(cp),  intent(inout) :: f_local(n_phi_max, n_theta_loc)
      real(cp),  intent(out)   :: f_global(n_phi_max, n_theta_max)
      
      integer :: i, nt, ierr
      integer :: Rq(n_procs_theta) 
      
      ! Copies local content to f_global
      f_global = 0.0
      f_global(:,n_theta_beg:n_theta_end) = f_local
      
      do i=0,n_procs_theta-1
         nt = n_theta_dist(i,2) - n_theta_dist(i,1) + 1
         CALL MPI_IBCAST(f_global(1:n_phi_max,n_theta_dist(i,1):n_theta_dist(i,2)), &
                         n_phi_max*nt, MPI_DOUBLE_PRECISION, i, comm_theta, Rq(i+1), ierr)
      end do
      
      CALL MPI_WAITALL(n_procs_theta, Rq, MPI_STATUSES_IGNORE, ierr)
      
   end subroutine gather_f
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
!------------------------------------------------------------------------------
   subroutine transpose_theta_m(f_theta_m, f_m_theta)
!@>details Does the transposition using alltoallv.
!
!@>TODO this with mpi_type to stride the data
!@>author Rafael Lago (MPCDF) August 2017
!------------------------------------------------------------------------------
      complex(cp), intent(inout) :: f_theta_m(n_theta_max, n_m_loc)
      complex(cp), intent(inout) :: f_m_theta(m_max+1,n_theta_loc)
      
      complex(cp) :: sendbuf(n_m_loc*n_theta_max)
      complex(cp) :: recvbuf(n_theta_loc,m_max+1)
      
      integer :: sendcount(0:n_procs_theta-1)
      integer :: recvcount(0:n_procs_theta-1)
      integer :: senddispl(0:n_procs_theta-1)
      integer :: recvdispl(0:n_procs_theta-1)
      integer :: i, j, pos, n_theta, s_theta, e_theta
      integer :: m_idx(n_procs_theta*n_m_ext)
      
      recvcount = 0
      pos = 1
      do i=0,n_procs_theta-1
         ! Copy each theta chunk so that the send buffer is contiguous
         !
         !>@TODO check performance of this; implementing this with mpi_type
         !> striding the data will probably be faster
         senddispl(i) = pos-1
         s_theta = n_theta_dist(i,1)
         e_theta = n_theta_dist(i,2)
         n_theta = e_theta - s_theta + 1
         do j=1, n_m_loc
            sendbuf(pos:pos + n_theta - 1) = f_theta_m(s_theta:e_theta,j)
            pos = pos + n_theta
         end do
         sendcount(i) = pos - senddispl(i) - 1
         
         recvdispl(i) = sum(recvcount)
         recvcount(i) = (n_m_ext-1) * n_theta
         if (lmP_dist(i,n_m_ext,1) > -1) then
            recvcount(i) = recvcount(i) + n_theta
         end if
      end do
      
      call MPI_ALLTOALLV(sendbuf, sendcount, senddispl, MPI_DOUBLE_COMPLEX, &
                         recvbuf, recvcount, recvdispl, MPI_DOUBLE_COMPLEX, &
                         comm_theta, i)
      
      ! Now we reorder the receiver buffer. If the m distribution looks like:
      ! rank 0: 0, 4,  8, 12, 16
      ! rank 1: 1, 5,  9, 13
      ! rank 2: 2, 6, 10, 14
      ! rank 3: 3, 7, 11, 15
      ! then the columns of recvbuf are ordered as 0,4,8,12,16,1,5,9,13(...)
      ! and so forth. m_idx will contain this ordering (+1):
      m_idx = reshape(transpose(lmP_dist(:,:,1)),(/n_procs_theta*n_m_ext/)) + 1
      j = 1
      do i = 1, n_procs_theta*n_m_ext
         if (m_idx(i) == 0) cycle
         f_m_theta(m_idx(i),:) = recvbuf(:,j)
         j = j + 1
      end do
      
   end subroutine transpose_theta_m
end module shtns

