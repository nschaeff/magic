module parallel_mod
   !
   !  This module contains the blocking information
   !

#ifdef WITH_MPI
   use MPI
#endif
   use omp_lib
   use logic, only: l_save_out, lVerbose

   implicit none

   integer :: nThreads
   integer :: rank, n_procs
   integer :: n_procs_r, n_procs_theta
   integer :: comm_cart, comm_r, comm_theta
   integer :: coord_r, coord_m, coord_theta
   integer :: nR_per_rank,nR_on_last_rank
   integer :: nLMBs_per_rank
   integer :: rank_with_l1m0
   integer :: chunksize
   integer :: ierr
   integer, allocatable, protected :: cart2rank(:,:), rank2cart_theta(:), rank2cart_r(:)

contains

   subroutine parallel

      !--- Get number (name) of processor
#ifdef WITH_MPI
      call mpi_comm_rank(MPI_COMM_WORLD,rank,ierr)
      call mpi_comm_size(MPI_COMM_WORLD,n_procs, ierr)
      !write(*,"(A,I3,A,I3)") "Running MPI coord_r no. ",coord_r," out of ",n_procs
#else
      coord_r    = 0
      n_procs = 1
      ierr    = 0
#endif

#ifdef WITHOMP
      nThreads = omp_get_max_threads()
#else
      nThreads = 1
#endif

      chunksize=16

   end subroutine parallel
!------------------------------------------------------------------------------
   subroutine initialize_cartesian
#ifdef WITH_MPI
      integer :: dims(2), coords(2), colour, i
      logical :: periods(2)
      dims    = (/n_procs_theta, n_procs_r/)
      periods = (/.true., .false./)
      
      
      call MPI_Cart_Create(MPI_COMM_WORLD, 2, dims, periods, .true., comm_cart, ierr)
      call check_MPI_error(ierr)
      
      ! Overwrite "rank" ! Make sure later that 0 in MPI_COMM_WORLD is the same 0 as in this comm!
      call MPI_Comm_Rank(comm_cart, rank, ierr) 
      call check_MPI_error(ierr)
      
      call MPI_Cart_Coords(comm_cart,rank,2,coords,ierr)
      call check_MPI_error(ierr)
      
      call MPI_Comm_Split(comm_cart, coords(1), rank, comm_r, ierr)
      call check_MPI_error(ierr)
      call MPI_Comm_Rank(comm_r, coord_r, ierr) 
      call check_MPI_error(ierr)
      
      call MPI_Comm_Split(comm_cart, coords(2), rank, comm_theta, ierr)
      call MPI_Comm_Rank(comm_theta, coord_theta, ierr) 
      call check_MPI_error(ierr)
      
      if (rank .ne. 0) l_save_out = .false.
      if (rank .ne. 0) lVerbose   = .false.
      
      allocate(rank2cart_theta(0:n_procs-1))
      allocate(rank2cart_r(0:n_procs-1))
      allocate(cart2rank(0:n_procs_theta-1,0:n_procs_r-1))
      
      do i=0,n_procs-1
         call mpi_cart_coords(comm_cart, i, 2, coords, ierr)
         rank2cart_theta(i) = coords(1)
         rank2cart_r(i)     = coords(2)
         cart2rank(coords(1),coords(2)) = i
      end do
      
#endif WITH_MPI
   end subroutine initialize_cartesian
!------------------------------------------------------------------------------
   subroutine check_MPI_error(code)

      integer, intent(in) :: code
#ifdef WITH_MPI
      character(len=MPI_MAX_ERROR_STRING) :: error_str
      integer :: ierr, strlen
#endif

#ifdef WITH_MPI
      if (code /= MPI_SUCCESS) then
          call MPI_Error_string(code, error_str, strlen, ierr)
          write(*, '(A, A)') 'MPI error: ', trim(error_str)
          call MPI_Abort(MPI_COMM_WORLD, code, ierr)
      endif
#else
      if (code /= 0) then
          write(*, '(A, I4)') 'Error code: ', code
          stop
      endif
#endif

   end subroutine check_MPI_error
!------------------------------------------------------------------------------
   subroutine finalize_cartesian
#ifdef WITH_MPI
      
      ! Overwrite "rank" ! Make sure later that 0 in MPI_COMM_WORLD is the same 0 as in this comm!
      call MPI_Comm_Free(comm_cart, ierr) 
      call MPI_Comm_Free(comm_r, ierr) 
      call MPI_Comm_Free(comm_theta, ierr) 
      
      deallocate(rank2cart_theta, rank2cart_r, cart2rank)
#endif WITH_MPI
   end subroutine finalize_cartesian
!------------------------------------------------------------------------------
end module parallel_mod
