module truncation
   !
   ! This module defines the grid points and the truncation 
   !

   use precision_mod, only: cp
   use logic, only: l_finite_diff, lVerbose, l_chemical_conv, l_double_curl, &
       &            l_TP_form
   use useful, only: abortRun
   use parallel_mod
   use mpi
   
   implicit none

   integer :: n_r_max       ! number of radial grid points
   integer :: n_cheb_max    ! max degree-1 of cheb polynomia
   integer :: n_phi_tot     ! number of longitude grid points
   integer :: n_r_ic_max    ! number of grid points in inner core
   integer :: n_cheb_ic_max ! number of chebs in inner core
   !----------------------------------------------------------------------------
   ! Lago: This is important! minc here is sort of the opposite of 
   ! mres in SHTns. In MagIC, the number of m modes stored is 
   ! m_max/minc+1. In SHTns, it is simply m_max.
   ! Conversely, there are m_max m modes in MagIC, though not all of them
   ! are effectivelly stored/computed. In SHTns, the number of m modes
   ! is m_max*minc+1 (therein called mres). This makes things very confusing
   ! specially because lm2 and lmP2 arrays (and their relatives) store
   ! m_max m points, and sets those which are not multiple of minc to 0.
   !
   ! In other words, this:
   ! call shtns_lmidx(lm_idx, l_idx, m_idx/minc)
   ! returns the same as 
   ! lm_idx = lm2(l_idx, m_idx)
   !----------------------------------------------------------------------------
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
   integer :: n_m_max     ! max number of ms (different oders) = m_max/minc+1
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
   ! n_theta_dist(i,1): first theta point in the i-th rank
   ! n_theta_dist(i,2): last theta ppoint in the i-th rank
   ! n_theta_beg: shortcut to n_theta_dist(coord_theta,1)
   ! n_theta_end: shortcut to n_theta_dist(coord_theta,2)
   ! n_theta_loc: number of theta points in the local rank
   !-------------------------------------------------------------
   ! lm_dist - (n_procs_theta, n_m_ext, 4)
   ! lm_dist(i,j,1): value of the j-th "m" in the i-th rank
   ! lm_dist(i,j,2): length of the j-th row in Flm_loc
   ! lm_dist(i,j,3): where the j-th row begins in Flm_loc
   ! lm_dist(i,j,4): where the j-th row ends in Flm_loc
   !-------------------------------------------------------------
   ! lmP_dist(i,j,k) : same as above, but for l_max+1 instead of l_max
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
   integer :: n_m_ext, n_m_loc, lmP_loc, lm_loc
   integer :: lm_locMag
   integer :: lm_locChe
   integer :: lm_locTP
   integer :: lm_locDC
   integer, allocatable :: n_theta_dist(:,:), lmP_dist(:,:,:), lm_dist(:,:,:)
   
   !@>TODO: 
   ! lmP2_ptr is here because I cannot use blocking module here (it causes 
   ! circular dependency). The lmP2 is required, at the moment, for 
   ! gather_Flm only.  I can try to figure out a better way later, but as 
   ! of now, I just copy the pointer during distribute_truncation(). 
   ! - Lago 
   integer, private, pointer :: lmP2_ptr(:,:) , lm2_ptr(:,:) 
                              
 
   interface slice_Flm
      module procedure slice_Flm_cmplx, slice_Flm_real
   end interface slice_Flm
   
   interface slice_FlmP
      module procedure slice_FlmP_cmplx, slice_FlmP_real
   end interface slice_FlmP
 
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

!------------------------------------------------------------------------------- 
   subroutine distribute_truncation(lmP2, lm2)
!@>details Divides the number of points in theta direction evenly amongst the 
!> ranks; if the number of point is not round, the last first receive one extra 
!> point. 
!
!@>TODO load imbalance; load imbalance everywhere. I won't mind that now, 
!> because the ordering will be much more complicated anyway, I'll worry about 
!> that later
!
!@>TODO deallocate the arrays allocated here
!
!@>author Rafael Lago, MPCDF, July 2017
!-------------------------------------------------------------------------------
      integer, target, intent(in) :: lmP2(0:l_max+1,0:l_max+1), lm2(0:l_max,0:l_max)
      integer :: i, j, itmp, loc
      
      lmP2_ptr => lmP2
      lm2_ptr => lm2
      
      loc = n_theta_max/n_procs_theta
      itmp = n_theta_max - loc*n_procs_theta
      allocate(n_theta_dist(0:n_procs_theta-1,2))
      do i=0,n_procs_theta-1
         n_theta_dist(i,1) = loc*i + min(i,itmp) + 1
         n_theta_dist(i,2) = loc*(i+1) + min(i+1,itmp)
      end do
      n_theta_beg = n_theta_dist(coord_theta,1)
      n_theta_end = n_theta_dist(coord_theta,2)
      n_theta_loc = n_theta_dist(coord_theta,2) - n_theta_dist(coord_theta,1)+1
      
      n_m_ext = ceiling(real(m_max+1)/real(minc*n_procs_theta))
      
      allocate(lmP_dist(0:n_procs_theta-1, n_m_ext, 4))
      allocate(lm_dist(0:n_procs_theta-1, n_m_ext, 4))
      
      call distribute_snake(lmP_dist, l_max+1, m_max, minc)
      call distribute_snake(lm_dist,  l_max  , m_max, minc)
      
      n_m_loc = n_m_ext - 1
      if (lmP_dist(coord_theta,n_m_ext,1) > -1) n_m_loc = n_m_ext
      
      lmP_loc = sum(lmP_dist(coord_theta,:,2))
      lm_loc  = sum(lm_dist( coord_theta,:,2))
      lm_locMag = lm_loc
      lm_locChe = lm_loc
      lm_locTP  = lm_loc
      lm_locDC  = lm_loc
      if ( lm_maxMag==1 ) lm_locMag = 1
      if ( .not. l_chemical_conv ) lm_locChe = 1
      if ( .not. l_TP_form       ) lm_locTP  = 1
      if ( .not. l_double_curl   ) lm_locDC  = 1
      
      if (rank == 0) then
         print "('θ partition in rank ', I3, ': ', I5, I5, I5, ' points')", &
               0, n_theta_dist(0,1), n_theta_dist(0,2), n_theta_dist(0,2) - &
               n_theta_dist(0,1) + 1
         do i=1, n_procs_theta-1
            print "('               rank ', I3, ': ', I5, I5, I5, ' points')", &
                  i, n_theta_dist(i,1), n_theta_dist(i,2), n_theta_dist(i,2) - &
                  n_theta_dist(i,1) + 1
         end do

         do i=0, n_procs_theta-1
            do j=1, n_m_ext
               print "('lmP partition in rank ', I3, ': ', I5, I5, I5, I5)", i,&
               lmP_dist(i,j,1:4)
            end do
            print "('length in rank ', I3, ' is ', I5)", i, sum(lmP_dist(i,:,2))
         end do
      end if
      
   end subroutine distribute_truncation
   
!-------------------------------------------------------------------------------   
   subroutine distribute_round_robin(idx_dist, l_max, m_max, minc)
!@>details This will try to create a balanced workload, by assigning to each 
!> rank  non-consecutive m points. For instance, if we have m_max=16, minc=1 and 
!> n_procs_theta=4, we will have:
!> 
!> rank 0: 0, 4,  8, 12, 16
!> rank 1: 1, 5,  9, 13
!> rank 2: 2, 6, 10, 14
!> rank 3: 3, 7, 11, 15
!> 
!> This is not optimal, but seem to be a decent compromise.
! 
!@>author Rafael Lago, MPCDF, December 2017
!-------------------------------------------------------------------------------
      integer, intent(inout) :: idx_dist(0:n_procs_theta-1, n_m_ext, 4)
      integer, intent(in)    :: l_max, m_max, minc
      
      integer :: m_idx, proc_idx, row_idx
      
      idx_dist        = -minc
      idx_dist(:,:,2) =  0
      idx_dist(:,1,3) =  1
      
      m_idx=0
      row_idx = 1
      do while (m_idx <= m_max)
         do proc_idx=0, n_procs_theta-1
            idx_dist(proc_idx, row_idx, 1) = m_idx
            idx_dist(proc_idx, row_idx, 2) = l_max - m_idx + 1
            if (row_idx > 1) idx_dist(proc_idx, row_idx, 3) = &
                             idx_dist(proc_idx, row_idx-1, 4) + 1
            idx_dist(proc_idx, row_idx, 4) = idx_dist(proc_idx, row_idx, 3) &
                                           + l_max - m_idx
            
            m_idx = m_idx + minc
            if (m_idx > m_max) exit
         end do
         row_idx = row_idx + 1
      end do
     
   end subroutine distribute_round_robin
   
!-------------------------------------------------------------------------------   
   subroutine distribute_snake(idx_dist, l_max, m_max, minc)
!@>details Same as above but for Snake Ordering.
!@>author Rafael Lago, MPCDF, December 2017
!-------------------------------------------------------------------------------
      integer, intent(inout) :: idx_dist(0:n_procs_theta-1, n_m_ext, 4)
      integer, intent(in)    :: l_max, m_max, minc
      
      integer :: m_idx, proc_idx, row_idx
      
      idx_dist        = -minc
      idx_dist(:,:,2) =  0
      idx_dist(:,1,3) =  1
      
      m_idx=0
      row_idx = 1
      do while (m_idx <= m_max)
         do proc_idx=0, n_procs_theta-1
            idx_dist(proc_idx, row_idx, 1) = m_idx
            idx_dist(proc_idx, row_idx, 2) = l_max - m_idx + 1
            if (row_idx > 1) idx_dist(proc_idx, row_idx, 3) = &
                             idx_dist(proc_idx, row_idx-1, 4) + 1
            idx_dist(proc_idx, row_idx, 4) = idx_dist(proc_idx, row_idx, 3) &
                                           + l_max - m_idx
            m_idx = m_idx + minc
            if (m_idx > m_max) exit
         end do
         row_idx = row_idx + 1
         do proc_idx=n_procs_theta-1,0,-1
            if (m_idx > m_max) exit
            idx_dist(proc_idx, row_idx, 1) = m_idx
            idx_dist(proc_idx, row_idx, 2) = l_max - m_idx + 1
            if (row_idx > 1) idx_dist(proc_idx, row_idx, 3) = &
                             idx_dist(proc_idx, row_idx-1, 4) + 1
            idx_dist(proc_idx, row_idx, 4) = idx_dist(proc_idx, row_idx, 3) &
                                           + l_max - m_idx
            
            m_idx = m_idx + minc
         end do
         row_idx = row_idx + 1
      end do
      
     
   end subroutine distribute_snake
   
!-------------------------------------------------------------------------------
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

!-------------------------------------------------------------------------------
   subroutine slice_f(f_global, f_local)
!@>details Copies the relevant part of f_global into f_local
!@>author Rafael Lago (MPCDF) August 2017
!-------------------------------------------------------------------------------
      real(cp),  intent(in)  :: f_global(n_phi_max, n_theta_max)
      real(cp),  intent(out) :: f_local(n_phi_max, n_theta_loc)
      
      f_local = f_global(:,n_theta_beg:n_theta_end)
   end subroutine slice_f
   
!-------------------------------------------------------------------------------
   subroutine slice_Flm_cmplx(Flm_global, Flm_local)
!@>details Copies the relevant part of Flm_global into Flm_local
!@>author Rafael Lago (MPCDF) August 2017
!-------------------------------------------------------------------------------
      complex(cp),  intent(in)  :: Flm_global(lm_max)
      complex(cp),  intent(out) :: Flm_local(lm_loc)
      
      integer :: i, lm_s, lm_e, lm_gs, lm_ge, m_idx
      
      do i = 1, n_m_loc
        m_idx = lm_dist(coord_theta, i, 1)
        lm_s  = lm_dist(coord_theta, i, 3)
        lm_e  = lm_dist(coord_theta, i, 4)
        lm_gs = lm2_ptr(m_idx, m_idx)
        lm_ge = lm2_ptr(l_max, m_idx)
        Flm_local(lm_s:lm_e) = Flm_global(lm_gs:lm_ge)
      end do
   end subroutine slice_Flm_cmplx
   
!-------------------------------------------------------------------------------
   subroutine slice_Flm_real(Flm_global, Flm_local)
!@>details Copies the relevant part of Flm_global into Flm_local
!@>author Rafael Lago (MPCDF) August 2017
!-------------------------------------------------------------------------------
      real(cp),  intent(in)  :: Flm_global(lm_max)
      real(cp),  intent(out) :: Flm_local(lm_loc)
      
      integer :: i, lm_s, lm_e, lm_gs, lm_ge, m_idx
      
      do i = 1, n_m_loc
        m_idx = lm_dist(coord_theta, i, 1)
        lm_s  = lm_dist(coord_theta, i, 3)
        lm_e  = lm_dist(coord_theta, i, 4)
        lm_gs = lm2_ptr(m_idx, m_idx)
        lm_ge = lm2_ptr(l_max, m_idx)
        Flm_local(lm_s:lm_e) = Flm_global(lm_gs:lm_ge)
      end do
   end subroutine slice_Flm_real

!-------------------------------------------------------------------------------
   subroutine slice_FlmP_cmplx(FlmP_global, FlmP_local)
!@>details Copies the relevant part of FlmP_global into FlmP_local
!@>author Rafael Lago (MPCDF) August 2017
!-------------------------------------------------------------------------------
      complex(cp),  intent(in)  :: FlmP_global(lmP_max)
      complex(cp),  intent(out) :: FlmP_local(lmP_loc)
      
      integer :: i, lm_s, lm_e, lm_gs, lm_ge, m_idx
      
      call shtns_load_cfg(1) ! l_max + 1
      
      do i = 1, n_m_loc
        m_idx = lmP_dist(coord_theta, i, 1)
        lm_s  = lmP_dist(coord_theta, i, 3)
        lm_e  = lmP_dist(coord_theta, i, 4)
        lm_gs = lmP2_ptr(m_idx,   m_idx)
        lm_ge = lmP2_ptr(l_max+1, m_idx)
        FlmP_local(lm_s:lm_e) = FlmP_global(lm_gs:lm_ge)
      end do
      
      call shtns_load_cfg(0) ! l_max
      
   end subroutine slice_FlmP_cmplx
   
!-------------------------------------------------------------------------------
   subroutine slice_FlmP_real(FlmP_global, FlmP_local)
!@>details Copies the relevant part of FlmP_global into FlmP_local
!@>author Rafael Lago (MPCDF) August 2017
!-------------------------------------------------------------------------------
      real(cp),  intent(in)  :: FlmP_global(lmP_max)
      real(cp),  intent(out) :: FlmP_local(lmP_loc)
      
      integer :: i, lm_s, lm_e, lm_gs, lm_ge, m_idx
      
      call shtns_load_cfg(1) ! l_max + 1
      
      do i = 1, n_m_loc
        m_idx = lm_dist(coord_theta, i, 1)
        lm_s  = lm_dist(coord_theta, i, 3)
        lm_e  = lm_dist(coord_theta, i, 4)
        lm_gs = lmP2_ptr(m_idx,   m_idx)
        lm_ge = lmP2_ptr(l_max+1, m_idx)
        FlmP_local(lm_s:lm_e) = FlmP_global(lm_gs:lm_ge)
      end do
      
      call shtns_load_cfg(0) ! l_max
      
   end subroutine slice_FlmP_real

!-------------------------------------------------------------------------------
   subroutine gather_FlmP(FlmP_local, FlmP_global)
!@>details Gathers the FlmP which was computed in a distributed fashion.
!> Mostly used for debugging. This function may not be performant at all.
!> 
!@>TODO this with mpi_allgatherv
!@>author Rafael Lago (MPCDF) August 2017
!-------------------------------------------------------------------------------
      complex(cp),  intent(in)  :: FlmP_local(lmP_loc)
      complex(cp),  intent(out) :: FlmP_global(lmP_max)
      
      complex(cp) ::  buffer(lmP_max)
      integer :: i, j, m_idx, lm_s_local, lm_e_local, lm_s_global, lm_e_global
      integer :: pos, ilen, Rq(n_procs_theta), ierr
      
      ! buffer will receive all messages, but they are ordered by ranks,
      ! not by m.
      
      pos = 1
      do i=0,n_procs_theta-1
         ilen = sum(lmP_dist(i,:,2))
         if (coord_theta == i) buffer(pos:pos+ilen-1) = FlmP_local(1:lmP_loc)
         CALL MPI_IBCAST(buffer(pos:pos+ilen-1), ilen, MPI_DOUBLE_COMPLEX, i, &
                         comm_theta, Rq(i+1), ierr)
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
            lm_s_global = lmP2_ptr(m_idx  ,m_idx)
            lm_e_global = lmP2_ptr(l_max+1,m_idx)
            FlmP_global(lm_s_global:lm_e_global) = buffer(lm_s_local:lm_e_local)
         end do
         pos = pos + sum(lmP_dist(i,:,2))
      end do
      
   end subroutine gather_FlmP
   
!-------------------------------------------------------------------------------
   subroutine gather_Flm(Flm_local, Flm_global)
!@>details Gathers the Flm which was computed in a distributed fashion.
!> Mostly used for debugging. This function may not be performant at all.
!> 
!@>TODO this with mpi_allgatherv
!@>author Rafael Lago (MPCDF) August 2017
!-------------------------------------------------------------------------------
      complex(cp),  intent(in)  :: Flm_local(lm_loc)
      complex(cp),  intent(out) :: Flm_global(lm_max)
      
      complex(cp) ::  buffer(lm_max)
      integer :: i, j, m_idx, lm_s_local, lm_e_local, lm_s_global, lm_e_global
      integer :: pos, ilen, Rq(n_procs_theta), ierr
      
      ! buffer will receive all messages, but they are ordered by ranks,
      ! not by m.
      
      pos = 1
      do i=0,n_procs_theta-1
         ilen = sum(lm_dist(i,:,2))
         if (coord_theta == i) buffer(pos:pos+ilen-1) = Flm_local(1:lm_loc)
         CALL MPI_IBCAST(buffer(pos:pos+ilen-1), ilen, MPI_DOUBLE_COMPLEX, i, &
                         comm_theta, Rq(i+1), ierr)
         pos = pos + ilen
      end do
      
      CALL MPI_WAITALL(n_procs_theta, Rq, MPI_STATUSES_IGNORE, ierr)
      
      pos = 0
      do i=0,n_procs_theta-1
         do j = 1, n_m_ext
            m_idx = lm_dist(i, j, 1)
            if (m_idx < 0) exit
            lm_s_local  = pos + lm_dist(i, j, 3)
            lm_e_local  = pos + lm_dist(i, j, 4)
            lm_s_global = lm2_ptr(m_idx , m_idx)
            lm_e_global = lm2_ptr(l_max , m_idx)
            Flm_global(lm_s_global:lm_e_global) = buffer(lm_s_local:lm_e_local)
         end do
         pos = pos + sum(lm_dist(i,:,2))
      end do
      
   end subroutine gather_Flm
   
!-------------------------------------------------------------------------------
   subroutine gather_f(f_local, f_global)
!@>details Gathers the Flm which was computed in a distributed fashion.
!> Mostly used for debugging. This function may not be performant at all.
!
!@>TODO this with mpi_allgatherv
!@>author Rafael Lago (MPCDF) August 2017
!-------------------------------------------------------------------------------
      real(cp),  intent(inout) :: f_local(n_phi_max, n_theta_loc)
      real(cp),  intent(out)   :: f_global(n_phi_max, n_theta_max)
      
      integer :: i, nt, ierr
      integer :: Rq(n_procs_theta) 
      
      ! Copies local content to f_global
      f_global = 0.0
      f_global(:,n_theta_beg:n_theta_end) = f_local
      
      do i=0,n_procs_theta-1
         nt = n_theta_dist(i,2) - n_theta_dist(i,1) + 1
         CALL MPI_IBCAST(f_global(:,n_theta_dist(i,1):n_theta_dist(i,2)),   &
                         n_phi_max*nt, MPI_DOUBLE_PRECISION, i, comm_theta, &
                         Rq(i+1), ierr)
      end do
      
      CALL MPI_WAITALL(n_procs_theta, Rq, MPI_STATUSES_IGNORE, ierr)
      
   end subroutine gather_f
   
!-------------------------------------------------------------------------------
   subroutine transpose_m_theta(f_m_theta, f_theta_m)
!@>details Does the transposition using alltoallv.
!
!@>TODO this with mpi_type to stride the data
!@>author Rafael Lago (MPCDF) August 2017
!-------------------------------------------------------------------------------
      complex(cp), intent(inout) :: f_m_theta(n_m_max,n_theta_loc)
      complex(cp), intent(inout) :: f_theta_m(n_theta_max, n_m_loc)
      
      complex(cp) :: sendbuf(n_m_max*n_theta_loc)
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
               m_idx = lmP_dist(i,j,1)/minc
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
   
!-------------------------------------------------------------------------------
   subroutine transpose_theta_m(f_theta_m, f_m_theta)
!@>details Does the transposition using alltoallv.
!
!@>TODO this with mpi_type to stride the data
!@>author Rafael Lago (MPCDF) August 2017
!-------------------------------------------------------------------------------
      complex(cp), intent(inout) :: f_theta_m(n_theta_max, n_m_loc)
      complex(cp), intent(inout) :: f_m_theta(n_m_max,n_theta_loc)
      
      complex(cp) :: sendbuf(n_m_loc*n_theta_max)
      complex(cp) :: recvbuf(n_theta_loc,n_m_max)
      
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
      m_idx = reshape(transpose(lmP_dist(:,:,1)),(/n_procs_theta*n_m_ext/))/minc + 1
      j = 1
      do i = 1, n_procs_theta*n_m_ext
         if (m_idx(i) < 1) cycle
         f_m_theta(m_idx(i),:) = recvbuf(:,j)
         j = j + 1
      end do
      
   end subroutine transpose_theta_m
!-------------------------------------------------------------------------------
end module truncation
