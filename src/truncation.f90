module truncation
   !
   ! This module defines the grid points and the truncation 
   !

   use precision_mod, only: cp
   use logic, only: l_finite_diff, lVerbose
   use useful, only: abortRun
   use parallel_mod, only: rank, coord_r, coord_theta, n_procs_theta
   use mpi

   implicit none

   integer :: n_r_max       ! number of radial grid points
   integer :: n_cheb_max    ! max degree-1 of cheb polynomia
   integer :: n_phi_tot     ! number of longitude grid points
   integer :: n_r_ic_max    ! number of grid points in inner core
   integer :: n_cheb_ic_max ! number of chebs in inner core
   integer :: minc          ! basic wavenumber, longitude symmetry  
   integer :: nalias        ! controls dealiasing in latitude
   logical :: l_axi         ! logical for axisymmetric calculations
   character(len=72) :: radial_scheme ! radial scheme (either Cheybev of FD)
   real(cp) :: fd_stretch    ! regular intervals over irregular intervals
   real(cp) :: fd_ratio      ! drMin over drMax (only when FD are used)
   integer :: fd_order       ! Finite difference order (for now only 2 and 4 are safe)
   integer :: fd_order_bound ! Finite difference order on the  boundaries
 
   !-- Derived quantities:
   integer :: n_phi_max   ! absolute number of phi grid-points
   integer :: n_theta_max ! number of theta grid-points
   integer :: n_theta_axi ! number of theta grid-points (axisymmetric models)
   integer :: l_max       ! max degree of Plms
   integer :: m_max       ! max order of Plms
   integer :: n_m_max     ! max number of ms (different oders)
   integer :: lm_max      ! number of l/m combinations
   integer :: lmP_max     ! number of l/m combination if l runs to l_max+1
   integer :: lm_max_real ! number of l/m combination for real representation (cos/sin)
   integer :: nrp,ncp     ! dimension of phi points in for real/complex arrays
   integer :: n_r_tot     ! total number of radial grid points
 
   !--- Now quantities for magnetic fields:
   !    Set lMag=0 if you want to save this memory (see c_fields)!
   integer :: lMagMem       ! Memory for magnetic field calculation
   integer :: n_r_maxMag    ! Number of radial points to calculate magnetic field
   integer :: n_r_ic_maxMag ! Number of radial points to calculate IC magnetic field
   integer :: n_r_totMag    ! n_r_maxMag + n_r_ic_maxMag
   integer :: l_maxMag      ! Max. degree for magnetic field calculation
   integer :: lm_maxMag     ! Max. number of l/m combinations for magnetic field calculation
 
   !-- Movie memory control:
   integer :: lMovieMem      ! Memory for movies
   integer :: ldtBMem        ! Memory for movie output
   integer :: lm_max_dtB     ! Number of l/m combinations for movie output
   integer :: n_r_max_dtB    ! Number of radial points for movie output
   integer :: n_r_ic_max_dtB ! Number of IC radial points for movie output
   integer :: lmP_max_dtB    ! Number of l/m combinations for movie output if l runs to l_max+1
 
   !--- Memory control for stress output:
   integer :: lStressMem     ! Memory for stress output
   integer :: n_r_maxStr     ! Number of radial points for stress output
   integer :: n_theta_maxStr ! Number of theta points for stress output
   integer :: n_phi_maxStr   ! Number of phi points for stress output
   
   !--- Distributed domain - for n_procs_theta > 1
   ! lmP_dist - (n_procs_theta, 2)
   ! n_theta_dist(i,1): first theta point in the i-th rank
   ! n_theta_dist(i,2): last theta ppoint in the i-th rank
   ! n_theta_beg: shortcut to n_theta_dist(coord_theta,1)
   ! n_theta_end: shortcut to n_theta_dist(coord_theta,2)
   ! n_theta_loc: number of theta points in the local rank
   !-------------------------------------------------------------
   ! lmP_dist - (n_procs_theta, n_m_ext, 4)
   ! lmP_dist(i,j,1): value of the j-th "m" in the i-th rank
   ! lmP_dist(i,j,2): length of the j-th row in fLM
   ! lmP_dist(i,j,3): where the j-th row begins in fLM
   ! lmP_dist(i,j,4): where the j-th row ends in fLM
   !-------------------------------------------------------------
   ! n_m_ext: number of m points in the local rank, +1. The extra
   !          point comes in when m_max + 1 is not divisible by the 
   !          number of ranks. If the extra point is "not in use"
   !          in the j-th rank, then lmP_dist(j,n_m_ext,:) = -1,
   !          and lmP_dist(j,n_m_ext,2) = 0
   !          This is useful mostly for looping over lmP_dist
   ! n_m_loc: number of points in the local rank. This coincides with
   !          n_m_ext in the ranks with extra points.
   ! 
   ! lmP_loc: shortcut to sum(lmP_dist(coord_theta,:,2))
   ! 
   ! - Lago
   integer :: n_theta_beg, n_theta_end, n_theta_loc
   integer :: n_m_ext, n_m_loc, lmP_loc
   integer, allocatable :: n_theta_dist(:,:), lmP_dist(:,:,:)
!    type(mpi_datatype) :: m_theta_types(:)
 
contains

   subroutine initialize_truncation

      integer :: n_r_maxML,n_r_ic_maxML,n_r_totML,l_maxML,lm_maxML
      integer :: lm_max_dL,lmP_max_dL,n_r_max_dL,n_r_ic_max_dL
      integer :: n_r_maxSL,n_theta_maxSL,n_phi_maxSL

      if ( .not. l_axi ) then
         ! absolute number of phi grid-points
         n_phi_max=n_phi_tot/minc

         ! number of theta grid-points
         n_theta_max=n_phi_tot/2

         ! max degree and order of Plms
         l_max=(nalias*n_theta_max)/30 
         m_max=(l_max/minc)*minc
      else
         n_theta_max=n_theta_axi
         n_phi_max  =1
         n_phi_tot  =1
         minc       =1
         l_max      =(nalias*n_theta_max)/30 
         m_max      =0
      end if

      ! max number of ms (different oders)
      n_m_max=m_max/minc+1

      ! number of l/m combinations
      lm_max=m_max*(l_max+1)/minc - m_max*(m_max-minc)/(2*minc)+(l_max+1-m_max)
      ! number of l/m combination if l runs to l_max+1
      lmP_max=lm_max+n_m_max

      ! number of l/m combination 
      ! for real representation (cos/sin)
      lm_max_real=2*lm_max

#if WITH_SHTNS
      nrp=n_phi_max
#else
      nrp=n_phi_max+2
#endif
      ncp=nrp/2

      ! total number of radial grid points
      n_r_tot=n_r_max+n_r_ic_max


      !--- Now quantities for magnetic fields:
      !    Set lMag=0 if you want to save this memory (see c_fields)!
      n_r_maxML     = lMagMem*n_r_max
      n_r_ic_maxML  = lMagMem*n_r_ic_max
      n_r_totML     = lMagMem*n_r_tot
      l_maxML       = lMagMem*l_max
      lm_maxML      = lMagMem*lm_max
      n_r_maxMag    = max(1,n_r_maxML)
      n_r_ic_maxMag = max(1,n_r_ic_maxML)
      n_r_totMag    = max(1,n_r_totML)
      l_maxMag      = max(1,l_maxML)
      lm_maxMag     = max(1,lm_maxML)

      !-- Movie memory control:
      lm_max_dL    =ldtBMem*lm_max
      lmP_max_dL   =ldtBMem*lmP_max
      n_r_max_dL   =ldtBMem*n_r_max
      n_r_ic_max_dL=ldtBMem*n_r_ic_max
      lm_max_dtB    =max(lm_max_DL,1) 
      lmP_max_dtB   =max(lmP_max_DL,1)
      n_r_max_dtB   =max(n_r_max_DL,1)
      n_r_ic_max_dtB=max(n_r_ic_max_DL,1)

      !--- Memory control for stress output:
      n_r_maxSL     =lStressMem*n_r_max
      n_theta_maxSL =lStressMem*n_theta_max
      n_phi_maxSL   =lStressMem*n_phi_max
      n_r_maxStr    =max(n_r_maxSL,1)
      n_theta_maxStr=max(n_theta_maxSL,1)
      n_phi_maxStr  =max(n_phi_maxSL,1)
   end subroutine initialize_truncation

!--------------------------------------------------------------------------------   
   subroutine distribute_truncation(lmP2)
!@>details Divides the number of points in theta direction evenly amongst the 
!> ranks; if the number of point is not round, the last first receive one extra 
!> point. 
!
!@>TODO load imbalance; load imbalance everywhere. I won't mind that now, because
!> the ordering will be much more complicated anyway, I'll worry about that later
!
!@>TODO deallocate the arrays allocated here
!
!@>author Rafael Lago, MPCDF, July 2017
!-----------------------------------------------------------------------
      integer, intent(in) :: lmP2(0:l_max+1,0:l_max+1)
      integer :: i, j, itmp, pos, loc
      
      
      loc = n_theta_max/n_procs_theta
      itmp = n_theta_max - loc*n_procs_theta
      allocate(n_theta_dist(0:n_procs_theta-1,2))
      do i=0,n_procs_theta-1
         n_theta_dist(i,1) = loc*i + min(i,itmp) + 1
         n_theta_dist(i,2) = loc*(i+1) + min(i+1,itmp)
      end do
      n_theta_beg = n_theta_dist(coord_theta,1)
      n_theta_end = n_theta_dist(coord_theta,2)
      n_theta_loc = n_theta_dist(coord_theta,2) - n_theta_dist(coord_theta,1) + 1
      
      
      ! --------------------------------------------------------------------------
      ! This will try to create a balanced workload, by assigning to each rank
      ! non-consecutive m points. For instance, if we have m_max=16 and 
      ! 4 processes in the theta direction, we will have:
      ! rank 0: 0, 4,  8, 12, 16
      ! rank 1: 1, 5,  9, 13
      ! rank 2: 2, 6, 10, 14
      ! rank 3: 3, 7, 11, 15
      ! 
      !>@TODO none of this will work in case minc =/= 1
      ! -------------------------------------------------------------------------- 
      n_m_ext = ceiling(real(m_max+1)/real(n_procs_theta))
      allocate(lmP_dist(0:n_procs_theta-1, n_m_ext, 4))
      lmP_dist        = -1
      lmP_dist(:,:,2) =  0
      do i=0, n_procs_theta-1
         j    = 1
         pos  = 1
         itmp = i
         do while (itmp <= m_max)
            !-------------------------------------------------------------
            ! lmP_dist(i,j,1): value of the j-th m in the i-th rank
            ! lmP_dist(i,j,2): length of the m-th row in fLM
            ! lmP_dist(i,j,3): where the m-th row begins in fLM
            ! lmP_dist(i,j,4): where the m-th row ends in fLM
            !-------------------------------------------------------------
            lmP_dist(i,j,1) = itmp         
            !>@TODO The line below can probably be replaced by l_max+2-itmp or something like that
            !>      but I better ask Thomas before, just to make sure
            lmP_dist(i,j,2) = lmP2(l_max+1,itmp) - lmP2(itmp, itmp) + 1
            lmP_dist(i,j,3) = pos
            lmP_dist(i,j,4) = pos + lmP_dist(i,j,2) - 1
            pos = lmP_dist(i,j,4) + 1
            itmp = itmp + n_procs_theta
            j = j + 1
         end do
      end do
      
      n_m_loc = n_m_ext - 1
      if (lmP_dist(coord_theta,n_m_ext,1) > -1) n_m_loc = n_m_ext
      lmP_loc = sum(lmP_dist(coord_theta,:,2))
      
      if (rank == 0) then
         print "('θ partition in rank ', I3, ': ', I5, I5, I5, ' points')", &
               0, n_theta_dist(0,1), n_theta_dist(0,2), n_theta_dist(0,2) - n_theta_dist(0,1) + 1
         do i=1, n_procs_theta-1
            print "('               rank ', I3, ': ', I5, I5, I5, ' points')", &
                  i, n_theta_dist(i,1), n_theta_dist(i,2), n_theta_dist(i,2) - n_theta_dist(i,1) + 1
         end do

         do i=0, n_procs_theta-1
            do j=1, n_m_ext
               print "('m partition in rank ', I3, ': ', I5, I5, I5, I5)", i, lmP_dist(i,j,1:4)
            end do
            print "('length in rank ', I3, ' is ', I5)", i, sum(lmP_dist(i,:,2))
         end do
         
      end if
      
   end subroutine distribute_truncation

!--------------------------------------------------------------------------------
   subroutine checkTruncation
      !  This function checks truncations and writes it
      !  into STDOUT and the log-file.                                 
      !  MPI: called only by the processor responsible for output !  

      if ( minc < 1 ) then
         call abortRun('! Wave number minc should be > 0!')
      end if
      if ( mod(n_phi_tot,minc) /= 0 .and. (.not. l_axi) ) then
         call abortRun('! Number of longitude grid points n_phi_tot must be a multiple of minc')
      end if
      if ( mod(n_phi_max,4) /= 0 .and. (.not. l_axi) ) then
         call abortRun('! Number of longitude grid points n_phi_tot/minc must be a multiple of 4')
      end if
      if ( mod(n_phi_tot,16) /= 0 .and. (.not. l_axi) ) then
         call abortRun('! Number of longitude grid points n_phi_tot must be a multiple of 16')
      end if
      if ( n_phi_max/2 <= 2 .and. (.not. l_axi) ) then
         call abortRun('! Number of longitude grid points n_phi_tot must be larger than 2*minc')
      end if

      !-- Checking radial grid:
      if ( .not. l_finite_diff ) then
         if ( mod(n_r_max-1,4) /= 0 ) then
            call abortRun('! Number n_r_max-1 should be a multiple of 4')
         end if
      end if

      if ( n_theta_max <= 2 ) then
         call abortRun('! Number of latitude grid points n_theta_max must be larger than 2')
      end if
      if ( mod(n_theta_max,4) /= 0 ) then
         call abortRun('! Number n_theta_max of latitude grid points be a multiple must be a multiple of 4')
      end if
      if ( n_cheb_max > n_r_max ) then
         call abortRun('! n_cheb_max should be <= n_r_max!')
      end if
      if ( n_cheb_max < 1 ) then
         call abortRun('! n_cheb_max should be > 1!')
      end if
      if ( (n_phi_max+1)*n_theta_max > lm_max_real*(n_r_max+2) ) then
         call abortRun('! (n_phi_max+1)*n_theta_max > lm_max_real*(n_r_max+2) !')
      end if

   end subroutine checkTruncation
!-----------------------------------------------------------------------------
end module truncation
