module distribute_mod
   !
   ! This module is supposed to replace 
   !  - blocking.f90
   !  - truncation.f90
   !  - distributed_theta.f90 (newly introduced by me, soon to be destroyed)
   ! 
   ! This takes the input parameters and decides how to split them. If the MPI
   ! partition was imposed via parameter file, it will follow whatever is there
   ! instead.
   !  
   !-- Author: Rafael Lago, MPCDF, May 2018
   !

   use precision_mod, only: cp
!    use logic, only: l_finite_diff, lVerbose, l_chemical_conv, l_double_curl, &
!                     l_TP_form
   use parallel_mod

   implicit none

   public
   
   integer :: n_r_max       ! number of radial grid points
   integer :: n_cheb_max    ! max degree-1 of cheb polynomia
   integer :: n_phi_tot     ! number of longitude grid points
   integer :: n_r_ic_max    ! number of grid points in inner core
   integer :: n_cheb_ic_max ! number of chebs in inner core
   
   integer :: n_r_cmb       ! ??? 
   integer :: n_r_icb       ! ??? 
   n_theta_max
   
   !-- Distribution of the points amonst the ranks
   !   Contiguous directions are stored in a 2D array
   !   Discontiguous directions are stored in a 3D array
   integer, allocatable :: gs_t_dist(:,:)
   integer, allocatable :: gs_r_dist(:,:)
   integer, allocatable :: lm_m_dist(:,:,:)
   integer, allocatable :: lm_r_dist(:,:)
   integer, allocatable :: ml_m_dist(:,:)
   integer, allocatable :: ml_l_dist(:,:,:)
   
   !------------------------------------------------------------------------------
   subroutine initialize_truncation_new
      !  
      !  This subroutine will 
      !  
      
      !-- I don't know the purpose of these two variables, 
      !   but I'll keep them around. - Lago
      n_r_cmb = 1
      n_r_icb = n_r_max
      
      call distribute_gs
      
      
      
   end subroutine initialize_truncation_new
   
   !------------------------------------------------------------------------------
   subroutine distribute_grid_r
      !
      ! This subroutine is used to set up the MPI decomposition in the
      ! radial direction
      !
      !-- Local variable:
      integer :: perRank

      perRank
      n_ranks_grid_theta
      
      
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

   end subroutine distribute_grid_r
!------------------------------------------------------------------------------
   
end module distribute_mod