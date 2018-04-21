module radial_data
   !
   ! This module defines the MPI decomposition in the radial direction.
   !

   use truncation, only: n_r_max
   use parallel_mod, only: coord_r, n_procs_r, nR_per_rank, nR_on_last_rank
   use logic, only: l_mag, lVerbose, l_finite_diff, l_chemical_conv, &
        &           l_double_curl, l_TP_form
 
   implicit none
 
   private
 
   integer, public :: nRstart,nRstop,nRstartMag,nRstopMag
   integer, public :: nRstartChe,nRstopChe,nRstartTP,nRstopTP
   integer, public :: nRstartDC,nRstopDC
   integer, public :: n_r_cmb,n_r_icb
 
   public :: initialize_radial_data

contains

   subroutine initialize_radial_data
      !
      ! This subroutine is used to set up the MPI decomposition in the
      ! radial direction
      !

      !-- Local variable:
      integer :: nR_remaining

      n_r_cmb=1
      n_r_icb=n_r_max

#ifdef WITH_MPI
      if ( .not. l_finite_diff ) then ! Cheb grid are restriced to odd numbers for now
         nR_per_rank = (n_r_max-1)/n_procs_r
         nRstart = n_r_cmb + coord_r*nR_per_rank
         nRstop  = n_r_cmb + (coord_r+1)*nR_per_rank - 1

         if ( coord_r == n_procs_r-1 ) then
            ! add the last point to the last process, which now has nR_per_rank+1
            ! radial points
            nRstop = nRstop+1
         end if
         nR_on_last_rank = nR_per_rank+1
      else ! In FD, any grid is allowed
         nR_per_rank = n_r_max/n_procs_r
         nRstart = n_r_cmb + coord_r*nR_per_rank
         nRstop  = n_r_cmb + (coord_r+1)*nR_per_rank - 1

         nR_remaining = n_r_max-(n_r_cmb + n_procs_r*nR_per_rank - 1)
         if ( coord_r == n_procs_r-1 ) then
            nRstop = nRstop+nR_remaining
         end if
         nR_on_last_rank = nR_per_rank+nR_remaining
      end if
#else
      nR_per_rank = n_r_max
      nR_on_last_rank = n_r_max
      nRstart = n_r_cmb
      nRstop  = n_r_icb
#endif
      nRstartMag = nRstart
      nRstartChe = nRstart
      nRstartTP  = nRstart
      nRstartDC  = nRstart
      nRstopMag  = nRstop
      nRstopChe  = nRstop
      nRstopTP   = nRstop
      nRstopDC   = nRstop
      if ( .not. l_mag           ) nRstartMag = 1
      if ( .not. l_chemical_conv ) nRstartChe = 1
      if ( .not. l_TP_form       ) nRstartTP  = 1
      if ( .not. l_double_curl   ) nRstartDC  = 1
      if ( .not. l_mag           ) nRstopMag  = 1
      if ( .not. l_chemical_conv ) nRstopChe  = 1
      if ( .not. l_TP_form       ) nRstopTP   = 1
      if ( .not. l_double_curl   ) nRstopDC   = 1

      if ( lVerbose ) then
         write(*,"(4(A,I4))") "On coord_r ",coord_r," nR is in (", &
               nRstart,",",nRstop,"), nR_per_rank is ",nR_per_rank
      end if

   end subroutine initialize_radial_data
!------------------------------------------------------------------------------
end module radial_data

