module distributed_theta
!@>details This module is the distributed counterpart of blocking.f90. I think 
!>         that it would be best to merge some of parallel.f90 and 
!>         truncation.f90 here too, but it is hard due to circular 
!>         dependencies.
!@>author Rafael Lago, MPCDF, April 2018
!-------------------------------------------------------------------------------

   use precision_mod
   use parallel_mod, only: coord_theta, rank, n_procs, n_procs_r, rank2cart_theta, rank2cart_r, nR_per_rank
   use truncation, only: lm_dist, lmP_dist, n_m_loc, lmP_max, lm_max, l_max, lm_loc, lmP_loc, n_r_max, n_theta_dist
   use mem_alloc, only: memWrite, bytes_allocated
   use logic, only: l_save_out, l_finite_diff
   use output_data, only: n_log_file, log_file
   use LMmapping, only: mappings, allocate_mappings, deallocate_mappings
   use useful, only: abortRun
   use radial_data, only: n_r_cmb
   use mpi  
   
   implicit none
 
   private
   
   public :: initialize_distributed_theta, finalize_distributed_theta, &
             get_distributed_lm_mapping
   
 
   type(mappings), public :: dist_map

contains

   subroutine initialize_distributed_theta
      
      integer :: i, iRstart, iRstop, iRremaining, iThetastart, iThetastop, ierr
      integer(lip) :: local_bytes_used
      

      local_bytes_used = bytes_allocated
      call allocate_mappings(dist_map,l_max,lm_loc,lmP_loc)
      call get_distributed_lm_mapping(dist_map)
      
      if ( rank == 0 ) then
         if ( l_save_out ) then
            open(newunit=n_log_file, file=log_file, status='unknown', position='append')
         end if

         iRremaining = n_r_max-(n_r_cmb + n_procs_r*nR_per_rank - 1)
         if ( .not. l_finite_diff ) iRremaining = 1
         write(*,*) '!-- Distributed (θ,r) per rank:'
         do i=0,n_procs-1
            iRstart = n_r_cmb + rank2cart_r(i)*nR_per_rank
            iRstop  = n_r_cmb + (rank2cart_r(i)+1)*nR_per_rank - 1
            iThetastart = n_theta_dist(rank2cart_theta(i),1)
            iThetastop  = n_theta_dist(rank2cart_theta(i),2)
            if ( rank2cart_r(i) == n_procs_r-1 ) iRstop = iRstop + iRremaining
            write(*,'(A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A)') ' ! ',i,' = (',iThetastart,':',&
                  iThetastop,',',iRstart,':', iRstop, &
                  '), ',iThetastop-iThetastart+1,'x',iRstop-iRstart+1,' points'
         end do

         if ( l_save_out ) close(n_log_file)
      end if
      
      local_bytes_used = bytes_allocated-local_bytes_used
      call memWrite('distributed_theta.f90', local_bytes_used)

   end subroutine initialize_distributed_theta
!------------------------------------------------------------------------
   subroutine finalize_distributed_theta

      call deallocate_mappings(dist_map)

   end subroutine finalize_distributed_theta
!-------------------------------------------------------------------------------
   subroutine get_distributed_lm_mapping(map)
!@>details This is the distributed counterpart of get_standard_lm_blocking.
!>         Notice that this does not use the ThetaBlocking paradigma, but 
!>         I kept in this module for the sake
!@>author Rafael Lago, MPCDF, April 2018
!-------------------------------------------------------------------------------
      type(mappings), intent(inout) :: map
      
      integer, parameter :: Invalid_Idx = -1   ! Will be used for invalid entries e.g. l=1, m=2
      integer, parameter :: Remote_Idx = -2    ! Will be used for valid entries which are in another rank
      
      ! Local variables
      integer :: m,l,lm,lmP,i
      
      map%lm2  = Remote_Idx
      map%lmP2 = Remote_Idx
      map%lm2l = Invalid_Idx
      map%lm2m = Invalid_Idx
      map%lmP2l  = Invalid_Idx
      map%lmP2m  = Invalid_Idx
      map%lm2lmP = Invalid_Idx
      map%lmP2lm = Invalid_Idx
      map%l2lmAS = Invalid_Idx
      
      !@>TODO lm_dist contains lm_s, lm_e as well as lm_e-lm_s+1. Maybe it would be best to 
      !@>remove it from lm_dist and keep it all centralized here.
      
      map%lm2mc = -5 ! Lago@180405: I have no idea what this one is for. It needs to be completely implemented
      
      lm  = 0
      lmP = 0
      do i = 1, n_m_loc
        m  = lm_dist(coord_theta, i, 1)
        map%lm2 (:,m) = Invalid_Idx
        map%lmP2(:,m) = Invalid_Idx
        
        do l=m,map%l_max
          lm  = lm + 1  ! lm_dist (coord_theta, i, 3)
          lmP = lmP + 1 ! lmP_dist(coord_theta, i, 3)
          map%lm2 (l,m) = lm
          map%lm2l(lm)  = l
          map%lm2m(lm)  = m
          
!           map%lm2mc(lm) = mc
          if ( m == 0 ) map%l2lmAS(l)=lm
          
          map%lmP2(l,m)  = lmP
          map%lmP2l(lmP) = l
          map%lmP2m(lmP) = m

          map%lm2lmP(lm)  = lmP
          map%lmP2lm(lmP) = lm
        end do
        lmP = lmP + 1
        map%lmP2(map%l_max+1,m) = lmP
        map%lmP2l(lmP) = map%l_max+1
        map%lmP2m(lmP) = m
        map%lmP2lm(lmP)= Invalid_Idx
        
      end do
      if ( lm /= lm_loc ) then
         write(*,"(2(A,I6))") 'Wrong lm=',lm," != lm_loc = ", lm_loc
         call abortRun('Stop run in distribute_theta')
      end if
      if ( lmP /= lmP_loc ) then
         write(*,"(3(A,I6))") 'Wrong lmP=',lm," != lmP_loc = ", lmP_loc
         call abortRun('Stop run in distribute_theta')
      end if

      map%lm2lmS   = Invalid_Idx
      map%lm2lmA   = Invalid_Idx
      map%lmP2lmPS = Invalid_Idx
      map%lmP2lmPA = Invalid_Idx
      
      do lm=1,lm_loc
         l=map%lm2l(lm)
         m=map%lm2m(lm)
         if ( l > 0 .and. l > m ) then
            map%lm2lmS(lm)=map%lm2(l-1,m)
         else
            map%lm2lmS(lm)=-1
         end if
         if ( l < map%l_max ) then
            map%lm2lmA(lm)=map%lm2(l+1,m)
         else
            map%lm2lmA(lm)=-1
         end if
      end do
      do lmP=1,lmP_loc
         l=map%lmP2l(lmP)
         m=map%lmP2m(lmP)
         if ( l > 0 .and. l > m ) then
            map%lmP2lmPS(lmP)=map%lmP2(l-1,m)
         else
            map%lmP2lmPS(lmP)=-1
         end if
         if ( l < map%l_max+1 ) then
            map%lmP2lmPA(lmP)=map%lmP2(l+1,m)
         else
            map%lmP2lmPA(lmP)=-1
         end if
      end do
      
   end subroutine get_distributed_lm_mapping
!-------------------------------------------------------------------------------
end module distributed_theta
