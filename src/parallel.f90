module parallel_mod
   !
   !  This module contains information about the MPI partitioning and 
   !  its communicators.
   !  
   ! There are 3 spaces in MagIC
   !    - Grid Space (φ, θ, r)
   !      - φ:    contiguous, local
   !      - θ:    contiguous, distributed
   !      - r:    contiguous, distributed
   !    - LM-Space (l, m, r) radia Loop
   !      - l:    contiguous, local
   !      - m: discontiguous, distributed
   !      - r:    contiguous, distributed
   !    - ML-Space (r, m, l) MLLoop
   !      - r:    contiguous, local
   !      - m:    contiguous, distributed
   !      - l: discontiguous, distributed
   !  
   ! The available MPI ranks are split into a 2D cartesian grid. 
   ! Notice that the cartesian for (l,m,r) *hast* to be the same 
   ! as (φ,θ,r). This is because the Legendre Transform (which 
   ! converts from (φ,θ) -> (l,m)) is called inside of the radial
   ! loop for each r independently. It would be very inefficient 
   ! (and messy) to reorganize the r's in the middle of the radial loop.
   ! 
   ! The communicators are:
   ! 
   ! - comm_gs: cartesian grid for (φ, θ, r)
   !   - comm_r: only the r direction of comm_gs
   !   - comm_theta: only the θ direction of comm_gs
   !   
   ! - comm_lm: copy of comm_gs
   !   - comm_r: same as above
   !   - comm_m: copy of comm_theta
   ! 
   ! - comm_ml: cartesian grid for (r, m, l) used in the MLloop
   !   - comm_ml_m: only the m direction of comm_ml_space
   !   - comm_ml_l: only the l direction of comm_ml_space
   ! 
   ! The duplications kinda make the code more readable, however.
   ! 
   !   Author: Rafael Lago (MPCDF) May 2018
   ! 

   use MPI
   use logic, only: l_save_out, lVerbose

   implicit none

   !   MPI_COMM_WORLD
   integer :: rank, n_ranks
   
   !   Grid Space
   integer ::    comm_gs
   integer ::    comm_r,    comm_theta
   integer :: n_ranks_r, n_ranks_theta
   integer ::   coord_r,   coord_theta
   integer, allocatable :: gs2rank(:,:), &
                                      rank2theta(:), &
                                      rank2r(:)
   
   !   LM-Space (radial Loop)
   !   Those will be just copies of variables above for theta
   integer ::    comm_lm
   integer ::    comm_m
   integer :: n_ranks_m
   integer ::   coord_m
   integer, allocatable :: lm2rank(:,:), &
                                      rank2m(:)
   
   !   ML-Space (ML Loop)
   integer ::    comm_ml
   integer ::    comm_ml_l,    comm_ml_m
   integer :: n_ranks_ml_l, n_ranks_ml_m
   integer ::   coord_ml_l,   coord_ml_m
   integer, allocatable :: ml2rank(:,:), &
                                      rank2ml_m(:), &
                                      rank2ml_l(:)
   
   !   Others (might be deprecated)
   integer :: nR_per_rank,nR_on_last_rank
   integer :: nThreads
   integer :: nLMBs_per_rank
   integer :: rank_with_l1m0
   integer :: chunksize
   integer :: ierr

contains
   
   !------------------------------------------------------------------------------
   subroutine initialize_mpi_world
      !
      !-- Get number (name) of processor
      !
      call mpi_comm_rank(MPI_COMM_WORLD,rank,    ierr)
      call mpi_comm_size(MPI_COMM_WORLD,n_ranks, ierr)

      nThreads  = 1
      chunksize = 16

      if (rank .ne. 0) l_save_out = .false.
      if (rank .ne. 0) lVerbose   = .false.
      
   end subroutine initialize_mpi_world
   
   !------------------------------------------------------------------------------
   subroutine initialize_mpi_decomposition
      !
      !   Author: Rafael Lago (MPCDF) May 2018
      ! 
      integer :: dims(2), coords(2), i, irank
      logical :: periods(2)
      
      if (rank == 0) write(*,*) ' '
      call optimize_decomposition_simple
      if (rank == 0) then
         write(*,*) '! MPI Domain Decomposition Info: '
         write(*,'(A,I0,A,I0)') ' ! Grid Space (θ,r): ', n_ranks_theta, "x", n_ranks_r
         write(*,'(A,I0,A,I0)') ' !   LM Space (m,r): ', n_ranks_m, "x", n_ranks_r
         write(*,'(A,I0,A,I0)') ' !   ML Space (l,m): ', n_ranks_ml_l, "x", n_ranks_ml_m
         write(*,'(A,I0)')      ' !      Total Ranks: ', n_ranks
      end if
      call check_decomposition
      
      !-- Grid Space
      !
      ! All arrays will be allocated as (θ,r), meaning that neighbouring ranks 
      ! should have contiguous θ. For that end, I need to flip the dimensions 
      ! when creating communicators using MPI because it uses row major instead 
      ! of Fortran's column major.
      dims    = (/n_ranks_r, n_ranks_theta/)
      periods = (/.true., .false./)
            
      call MPI_Cart_Create(MPI_COMM_WORLD, 2, dims, periods, .true., comm_gs, ierr)
      call check_MPI_error(ierr)
      
      call MPI_Comm_Rank(comm_gs, irank, ierr) 
      call check_MPI_error(ierr)
      
      call MPI_Cart_Coords(comm_gs, irank, 2, coords, ierr)
      call check_MPI_error(ierr)
      
      call MPI_Comm_Split(comm_gs, coords(2), irank, comm_r, ierr)
      call check_MPI_error(ierr)
      call MPI_Comm_Rank(comm_r, coord_r, ierr) 
      call check_MPI_error(ierr)
      
      call MPI_Comm_Split(comm_gs, coords(1), irank, comm_theta, ierr)
      call MPI_Comm_Rank(comm_theta, coord_theta, ierr) 
      call check_MPI_error(ierr)
      
      allocate(rank2theta(0:n_ranks-1))
      allocate(rank2r(0:n_ranks-1))
      allocate(gs2rank(0:n_ranks_theta-1,0:n_ranks_r-1))
      
      do i=0,n_ranks-1
         call mpi_cart_coords(comm_gs, i, 2, coords, ierr)
         rank2theta(i) = coords(2)
         rank2r(i)     = coords(1)
         gs2rank(coords(2),coords(1)) = i
      end do
      
      !-- LM Space (radial Loop)
      !   This is just a copy
      !   
      !   PS: is it a problem that I'm using a periodic communicator for 
      !       lm-space as well?
      comm_lm   = comm_gs
      comm_m    = comm_theta
      n_ranks_m = n_ranks_theta
      coord_m   = coord_theta
      
      allocate(rank2m(0:n_ranks-1))
      allocate(lm2rank(0:n_ranks_m-1,0:n_ranks_r-1))
      
      lm2rank   = gs2rank
      rank2m    = rank2theta
      
      !-- ML Space (ML Loop)
      !   
      dims    = (/n_ranks_ml_m, n_ranks_ml_l/)
      periods = (/.false., .false./)
            
      call MPI_Cart_Create(MPI_COMM_WORLD, 2, dims, periods, .true., comm_ml, ierr)
      call check_MPI_error(ierr)
      
      call MPI_Comm_Rank(comm_ml, irank, ierr) 
      call check_MPI_error(ierr)
      
      call MPI_Cart_Coords(comm_ml, irank, 2, coords, ierr)
      call check_MPI_error(ierr)
      
      call MPI_Comm_Split(comm_ml, coords(2), irank, comm_ml_m, ierr)
      call check_MPI_error(ierr)
      call MPI_Comm_Rank(comm_ml_m, coord_ml_m, ierr) 
      call check_MPI_error(ierr)
      
      call MPI_Comm_Split(comm_ml, coords(1), irank, comm_ml_l, ierr)
      call MPI_Comm_Rank(comm_ml_l, coord_ml_l, ierr) 
      call check_MPI_error(ierr)
      
      allocate(rank2ml_m(0:n_ranks-1))
      allocate(rank2ml_l(0:n_ranks-1))
      allocate(ml2rank(0:n_ranks_ml_l-1,0:n_ranks_ml_m-1))
      
      do i=0,n_ranks-1
         call mpi_cart_coords(comm_ml, i, 2, coords, ierr)
         rank2ml_l(i) = coords(2)
         rank2ml_m(i) = coords(1)
         ml2rank(coords(2),coords(1)) = i
      end do
      
   end subroutine initialize_mpi_decomposition
   
   !------------------------------------------------------------------------------
   subroutine check_MPI_error(code)

      integer, intent(in) :: code
      character(len=MPI_MAX_ERROR_STRING) :: error_str
      integer :: ierr, strlen

      if (code /= MPI_SUCCESS) then
          call MPI_Error_string(code, error_str, strlen, ierr)
          write(*, '(A, A)') 'MPI error: ', trim(error_str)
          call MPI_Abort(MPI_COMM_WORLD, code, ierr)
      endif

   end subroutine check_MPI_error
   
   !------------------------------------------------------------------------------
   subroutine finalize_mpi_decomposition

      call MPI_Comm_Free(comm_gs,   ierr) 
      call MPI_Comm_Free(comm_r, ierr) 
      call MPI_Comm_Free(comm_theta, ierr) 
      
      call MPI_Comm_Free(comm_ml,   ierr)
      call MPI_Comm_Free(comm_ml_m, ierr)
      call MPI_Comm_Free(comm_ml_l, ierr)
      
      deallocate(rank2theta, rank2r, gs2rank)
      deallocate(rank2m, lm2rank)
      deallocate(rank2ml_m, rank2ml_l, ml2rank)
   end subroutine finalize_mpi_decomposition

   !------------------------------------------------------------------------------
   subroutine optimize_decomposition_simple
      !   
      !   This is a *very* simple function to split all n_ranks into a 2D 
      !   grid. In the future we might want to make this function much more 
      !   sophisticated.
      !
      !   Author: Rafael Lago (MPCDF) May 2018
      ! 
      real     :: log2
      integer  :: nPow2
      
      !-- If the grid dimensions were not given
      if (n_ranks_theta == 0 .and. n_ranks_r == 0) then
         if (rank==0) write(*,*) '! Automatic splitting of MPI ranks for Grid Space...'
         
         log2 = log(real(n_ranks)) / log(2.0) ! Gets log2(n_ranks)
         n_ranks_theta = 2**FLOOR(log2/2)
         n_ranks_r = n_ranks/n_ranks_theta
      else if (n_ranks_theta == 0) then
         n_ranks_theta = n_ranks/n_ranks_r
      else if (n_ranks_r == 0) then
         n_ranks_r = n_ranks/n_ranks_theta
      end if
      
      n_ranks_m = n_ranks_theta
         
      n_ranks_ml_m = n_ranks_r
      n_ranks_ml_l = n_ranks_theta
      
   end subroutine optimize_decomposition_simple
   
   !------------------------------------------------------------------------------
   subroutine check_decomposition

      if (n_ranks_r * n_ranks_theta .NE. n_ranks) then
        print *, '! Invalid parallelization partition! Make sure that n_ranks_r'//&
                 '* n_ranks_theta equals the number of available processes!'
        stop
      end if
      if (n_ranks_r * n_ranks_m .NE. n_ranks) then
        print *, '! Invalid parallelization partition! Make sure that n_ranks_r'//&
                 '* n_ranks_m equals the number of available processes!'
        stop
      end if
      if (n_ranks_ml_l * n_ranks_ml_m .NE. n_ranks) then
        print *, '! Invalid parallelization partition! Make sure that n_ranks_ml_l'//&
                 '* n_ranks_ml_m equals the number of available processes!'
        stop
      end if

   end subroutine check_decomposition

!-------------------------------------------------------------------------------
end module parallel_mod
