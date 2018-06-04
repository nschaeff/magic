module radial_data
   !
   ! This module defines the MPI decomposition in the radial direction.
   !

   use truncation, only: n_r_max
   use parallel_mod, only: coord_r, n_ranks_r, nR_per_rank, nR_on_last_rank
   use logic, only: l_mag, lVerbose, l_finite_diff, l_chemical_conv, &
        &           l_double_curl, l_TP_form
 
   implicit none
 
   private
 
   integer, public :: l_r,u_r,l_r_Mag,u_r_Mag
   integer, public :: l_r_Che,u_r_Che,l_r_TP,u_r_TP
   integer, public :: l_r_DC,u_r_DC
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
         nR_per_rank = (n_r_max-1)/n_ranks_r
         l_r = n_r_cmb + coord_r*nR_per_rank
         u_r  = n_r_cmb + (coord_r+1)*nR_per_rank - 1

         if ( coord_r == n_ranks_r-1 ) then
            ! add the last point to the last process, which now has nR_per_rank+1
            ! radial points
            u_r = u_r+1
         end if
         nR_on_last_rank = nR_per_rank+1
      else ! In FD, any grid is allowed
         nR_per_rank = n_r_max/n_ranks_r
         l_r = n_r_cmb + coord_r*nR_per_rank
         u_r  = n_r_cmb + (coord_r+1)*nR_per_rank - 1

         nR_remaining = n_r_max-(n_r_cmb + n_ranks_r*nR_per_rank - 1)
         if ( coord_r == n_ranks_r-1 ) then
            u_r = u_r+nR_remaining
         end if
         nR_on_last_rank = nR_per_rank+nR_remaining
      end if
#else
      nR_per_rank = n_r_max
      nR_on_last_rank = n_r_max
      l_r = n_r_cmb
      u_r  = n_r_icb
#endif
      l_r_Mag = l_r
      l_r_Che = l_r
      l_r_TP  = l_r
      l_r_DC  = l_r
      u_r_Mag  = u_r
      u_r_Che  = u_r
      u_r_TP   = u_r
      u_r_DC   = u_r
      if ( .not. l_mag           ) l_r_Mag = 1
      if ( .not. l_chemical_conv ) l_r_Che = 1
      if ( .not. l_TP_form       ) l_r_TP  = 1
      if ( .not. l_double_curl   ) l_r_DC  = 1
      if ( .not. l_mag           ) u_r_Mag  = 1
      if ( .not. l_chemical_conv ) u_r_Che  = 1
      if ( .not. l_TP_form       ) u_r_TP   = 1
      if ( .not. l_double_curl   ) u_r_DC   = 1

      if ( lVerbose ) then
         write(*,"(4(A,I4))") "On coord_r ",coord_r," nR is in (", &
               l_r,",",u_r,"), nR_per_rank is ",nR_per_rank
      end if

   end subroutine initialize_radial_data
!------------------------------------------------------------------------------
end module radial_data

