module LMmapping
   !
   !   This module contains the mapping arrays, from the 2-dimensional (l,m) 
   !   coordinates onto the 1-dimensional (lm) coordinate. 
   !   

   use precision_mod
   use parallel_mod
   use geometry
   use mem_alloc, only: bytes_allocated, memWrite
   use useful, only: abortRun
   use mpi

   implicit none
 
   public
 
   !-- Mappings for (0:l_max, 0:m_max)
   !   
   !   lm2:  (l,m) -> (lm)       <misleading name>
   !   lm2l: (lm)  -> (l,.)
   !   lm2m: (lm)  -> (.,m)
   !   
   !   lmP*: same as above, for (0:l_max+1, 0:m_max)
   !   
   type, public :: mappings
      integer :: l_max
      integer :: m_max
      integer :: lm_max   ! deprecated; to be later replaced by %n_lm
      integer :: lmP_max  ! deprecated; to be later replaced by %n_lmP
      integer :: n_lm
      integer :: n_lmP
      integer, allocatable :: lm2(:,:),lm2l(:),lm2m(:)
      integer, allocatable :: lm2mc(:),l2lmAS(:)      
      integer, allocatable :: lm2lmS(:),lm2lmA(:)     
                                                     
      integer, allocatable :: lmP2(:,:),lmP2l(:),lmP2m(:)
      integer, allocatable :: lmP2lmPS(:),lmP2lmPA(:) 
                                                     
      integer, allocatable :: lm2lmP(:),lmP2lm(:)     
 
   end type mappings
   
   !
   ! n_ml_loc: number of ml tuples in this rank. This means that mlo goes
   !   from 1 to n_ml_loc in this rank. Each mlo corresponds to a (m,l) tuple.
   ! 
   ! dist_ml2(rank,m,l): the mlo index of tuplet (m,l) in rank. If this 
   ! tuple is not stored in "rank", then its value is -1.
   ! dist_ml2l(rank,mlo): the value of l for the index mlo in rank "rank"
   ! dist_ml2m(rank,mlo): the value of m for the index mlo in rank "rank"
   !
   ! ml2(m,l): a pointer to dist_ml2(coord_mlo,m,l)
   ! ml2l(mlo): a pointer to dist_ml2l(coord_mlo,mlo)
   ! ml2m(mlo): a pointer to dist_ml2m(coord_mlo,mlo)
   ! 
   ! li2lo(n_lo_loc): the list of all l's found in the current rank
   ! mi2mo(n_mo_loc): the list of all m's found in the current rank
   ! Notice that not all (m,l) pairs are present in this rank; for instance
   ! (2,5) and (1,6) might be in this rank, but (2,6) and (1,5) are not!
   ! 
   ! With the structure above you can loop over all local l and m as:
   !  do li=1,n_lo_loc
   !   l = li2lo(li)
   !   do mi=1,n_mo_loc
   !     m = mi2mo(mi)
   !     if (ml2%(m,l)<0) cycle
   !     *do stuff*
   !   end do
   !  end do
   ! 
   ! or as:
   !  do mlo_idx=1,n_mlo_loc
   !    l = ml2l(mlo_idx)
   !    m = ml2m(mlo_idx)
   !    /*do stuff*/
   !  end do
   ! It depends on the purpose of the loop, I guess.
   ! 
   !--------------------------------------------------------------------------
   type, public :: ml_mappings
      integer, allocatable :: i2l(:), i2m(:), i2ml(:)
      integer, allocatable :: ml2i(:,:), ml2coord(:,:)
      integer, allocatable :: milj2i(:,:), milj2m(:,:)
      integer, allocatable :: n_mi(:), lj2l(:)
      integer :: n_li
   end type ml_mappings
 
   !-- ????
   !   
   type, public :: subblocks_mappings
      integer :: nLMBs,l_max,m_max,sizeLMB2max
      integer, allocatable :: nLMBs2(:)
      integer, allocatable :: sizeLMB2(:,:)
      integer, allocatable :: lm22lm(:,:,:)
      integer, allocatable :: lm22l(:,:,:)
      integer, allocatable :: lm22m(:,:,:)
   end type subblocks_mappings
 
   type(mappings) :: map_dist_st, map_glbl_st
   type(ml_mappings) :: map_mlo
   
contains
   
   !----------------------------------------------------------------------------
   subroutine initialize_mapping
      !   
      !   Author: Rafael Lago, MPCDF, April 2018
      !   
      
      integer :: i, iRstart, iRstop, iRremaining, iThetastart, iThetastop
      integer(lip) :: local_bytes_used
      
      local_bytes_used = bytes_allocated
!       call allocate_mappings(map_glbl_st,l_max,lm_max,lmP_max,l_axi)
!       call allocate_mappings(lo_map,l_max,lm_max,lmP_max,l_axi)

      call allocate_mappings(map_glbl_st, lm_max, lmP_max, l_axi)
      call allocate_mappings(map_dist_st, n_lm_loc,n_lmP_loc,l_axi)
      call allocate_ml_mappings(map_mlo)
      
      ! (/0:n_m_max-1/)*minc: an array containing all m points
      call set_lmmapping_default(map_glbl_st, (/0:n_m_max-1/)*minc ) 
      call set_lmmapping_default(map_dist_st,   dist_m(coord_m,1:))
      
      call set_mlmapping(map_mlo)
      
      call print_mapping(map_dist_st, 'dist_st')
      call print_mapping(map_glbl_st, 'glbl_st')
      
      local_bytes_used = bytes_allocated-local_bytes_used
      call memWrite('LMmapping.f90', local_bytes_used)
   end subroutine initialize_mapping
   
   !----------------------------------------------------------------------------
   subroutine finalize_mapping

      call deallocate_mappings(map_glbl_st)
      call deallocate_mappings(map_dist_st)

   end subroutine finalize_mapping
   
   !----------------------------------------------------------------------------
   subroutine allocate_mappings(self, n_lm_len, n_lmP_len, l_axi)
      !   
      !   Allocates the mapping objects. 
      !   
      !   The l_max is taken from geometry module. It is supposed to represent 
      !   the global number of l points.
      !   n_lm_len is the number of local lm points. Similarly, n_lmP_len 
      !   is the same quantity, for l_max+1.
      !   
      !   So, if you're creating the map_dist_st, pass (n_lm_loc, n_lmP_len).
      !   If you're creating map_glbl_st, pass (n_lm_max, n_lmP_max).
      !   
      !   Author: Rafael Lago, MPCDF, April 2018
      !   

      type(mappings), intent(inout) :: self
      integer, intent(in) :: n_lm_len, n_lmP_len
      logical, intent(in) :: l_axi

      self%l_max = l_max
      if ( .not. l_axi ) then
         self%m_max = l_max
      else
         self%m_max = 0
      end if
      self%lm_max   = n_lm_len   ! deprecated; to be later replaced by %n_lm. Only used in blocking.f90
      self%lmP_max  = n_lmP_len  ! deprecated; to be later replaced by %n_lmP. Only used in blocking.f90
      self%n_lm  = n_lm_len
      self%n_lmP = n_lmP_len

      allocate( self%lm2(0:l_max,0:l_max),self%lm2l(n_lm_len),self%lm2m(n_lm_len))
      allocate( self%lm2mc(n_lm_len),self%l2lmAS(0:l_max) )
      allocate( self%lm2lmS(n_lm_len),self%lm2lmA(n_lm_len) )
      bytes_allocated = bytes_allocated + &
                        ((l_max+1)*(l_max+1)+5*n_lm_len+l_max+1)*SIZEOF_INTEGER

      allocate( self%lmP2(0:l_max+1,0:l_max+1),self%lmP2l(n_lmP_len) )
      allocate( self%lmP2m(n_lmP_len) )
      allocate( self%lmP2lmPS(n_lmP_len),self%lmP2lmPA(n_lmP_len) )
      allocate( self%lm2lmP(n_lm_len),self%lmP2lm(n_lmP_len) )
      bytes_allocated = bytes_allocated + &
                        ((l_max+2)*(l_max+2)+5*n_lmP_len+n_lm_len)*SIZEOF_INTEGER

   end subroutine allocate_mappings
   
   !----------------------------------------------------------------------------
   subroutine allocate_ml_mappings(self)
      !   
      !   Allocates the mapping objects. 
      !   
      !   The l_max is taken from geometry module. It is supposed to represent 
      !   the global number of l points.
      !   
      !   Author: Rafael Lago, MPCDF, November 2018
      !   
      type(ml_mappings), intent(inout) :: self
      
      allocate( self%ml2coord(0:l_max, 0:l_max) )
      
      ! These two will point to their target in set_mlmapping
      allocate( self%i2m(n_mlo_loc) )
      allocate( self%i2l(n_mlo_loc) )
      allocate( self%i2ml(n_mlo_loc) )
      
      allocate( self%ml2i(0:l_max, 0:l_max) )
      
!       allocate( self%li2l(n_lo_loc) )
!       allocate( self%mi2m(n_mo_loc) )
!       allocate( self%mi2m(n_mo_loc) )
      allocate( self%milj2i(n_mo_loc,n_lo_loc) )
      allocate( self%milj2m(n_mo_loc,n_lo_loc) )
      allocate( self%lj2l(n_lo_loc) )
      allocate( self%n_mi(n_lo_loc) )
      
!       bytes_allocated = bytes_allocated + &
!         (2*(l_max+1)*(l_max+1)+2*(l_max+1))*SIZEOF_INTEGER
                        
   end subroutine allocate_ml_mappings
   
   !----------------------------------------------------------------------------
   subroutine deallocate_mappings(self)
      type(mappings) :: self
      deallocate( self%lm2, self%lm2l, self%lm2m, self%lm2mc, self%l2lmAS )
      deallocate( self%lm2lmS, self%lm2lmA, self%lmP2, self%lmP2l )
      deallocate( self%lmP2m, self%lmP2lmPS, self%lmP2lmPA, self%lm2lmP )
      deallocate( self%lmP2lm )
   end subroutine deallocate_mappings
   
   !----------------------------------------------------------------------------
   subroutine deallocate_ml_mappings(self)
      type(ml_mappings) :: self
!       deallocate( self%ml2coord, self%ml2, self%ml2m, self%ml2l )
   end subroutine deallocate_ml_mappings
   
   !----------------------------------------------------------------------------
   subroutine allocate_subblocks_mappings(self,map,nLMBs,in_l_max,lmStartB, &
                                          lmStopB,l_axi)

      !-- Input variables
      type(subblocks_mappings)   :: self
      type(mappings), intent(in) :: map
      integer,        intent(in) :: nLMBs, in_l_max
      integer,        intent(in) :: lmStartB(nLMBs), lmStopB(nLMBs)
      logical,        intent(in) :: l_axi

      !-- Local variables
      integer :: nLMB,lm1,l1,max_size_of_subblock,lmStart,lmStop
      integer :: counter(0:in_l_max)

      self%nLMBs = nLMBs
      self%l_max = in_l_max

      if ( .not. l_axi ) then
         self%m_max = in_l_max
      else
         self%m_max = 0
      end if

      ! now determine the maximal size of a subblock (was sizeLMB2max parameter
      ! in former versions).
      max_size_of_subblock=0
      do nLMB=1,nLMBs
         lmStart=lmStartB(nLMB)
         lmStop =lmStopB(nLMB)
         counter=0
         do lm1=lmStart,lmStop
            l1=map%lm2l(lm1)
            counter(l1) = counter(l1) + 1
         end do
         if (maxval(counter) > max_size_of_subblock) then
            max_size_of_subblock=maxval(counter)
         end if
      end do
      self%sizeLMB2max = max_size_of_subblock

      allocate( self%nLMBs2(nLMBs),self%sizeLMB2(in_l_max+1,nLMBs) )
      allocate( self%lm22lm(self%sizeLMB2max,in_l_max+1,nLMBs) )
      allocate( self%lm22l(self%sizeLMB2max,in_l_max+1,nLMBs) )
      allocate( self%lm22m(self%sizeLMB2max,in_l_max+1,nLMBs) )
      bytes_allocated = bytes_allocated + &
                        (nLMBs+(in_l_max+1)*nLMBs+ &
                        3*(in_l_max+1)*nLMBS*self%sizeLMB2max)*SIZEOF_INTEGER

   end subroutine allocate_subblocks_mappings
   
   !----------------------------------------------------------------------------
   subroutine deallocate_subblocks_mappings(self)

      type(subblocks_mappings) :: self

      deallocate( self%nLMBs2, self%sizeLMB2, self%lm22lm, self%lm22l )
      deallocate( self%lm22m )

   end subroutine deallocate_subblocks_mappings
   
   !----------------------------------------------------------------------------
   !
   !   
   !                                 Set Methods
   !
   !
   !----------------------------------------------------------------------------
   subroutine set_lmmapping_default(map, m_arr)
      !   
      !   This function will place (all l_max) l's sequentially in increasing 
      !   order, e.g. for l_max=10 and m_arr=(/3,5/) the lm2 mapping will
      !   be:
      !   
      !   [(3,3), (4,3), ... (10,3), (5,5), (6,5), ... (10,5)]
      !   
      !   It makes no assumption about the lengths m_arr.
      !   
      !-- Remark: if you pass m_arr = (/0,1,2,...,m_max/) this function does 
      !      exactly the same distribution as the global lm map for minc=1.
      !   
      !   
      !-- Remark: this function can more or less easily be adapted to support 
      !      intervals of l, so, for instance, from l_l to u_l. Just keep this 
      !      in mind in case we need it in the future.
      !
      !   Author: Rafael Lago, MPCDF, April 2018
      !
      !-- TODO: double-check if lm2mc is correctly computed (or if it is needed
      !   at all anymore. I don't think it is.)
      !   
      type(mappings), intent(inout) :: map
      integer,        intent(in)    :: m_arr(:)
      
      integer, parameter :: Invalid_Idx = -1
      
      !-- Local variables
      integer :: m, l, m_idx, lm, lmP
      
      map%lm2      = Invalid_Idx
      map%lmP2     = Invalid_Idx
      map%lm2l     = Invalid_Idx
      map%lm2m     = Invalid_Idx
      map%lmP2l    = Invalid_Idx
      map%lmP2m    = Invalid_Idx
      map%lm2lmP   = Invalid_Idx
      map%lmP2lm   = Invalid_Idx
      map%l2lmAS   = Invalid_Idx
      map%lm2lmS   = Invalid_Idx
      map%lm2lmA   = Invalid_Idx
      map%lmP2lmPS = Invalid_Idx
      map%lmP2lmPA = Invalid_Idx
      map%lm2mc    = Invalid_Idx
      
      lm  = 0
      lmP = 0
      do m_idx = lbound(m_arr,1), ubound(m_arr,1)
         m  = m_arr(m_idx)
         if (m < 0) cycle
         
         do l=m,map%l_max
            lm  = lm  + 1 
            lmP = lmP + 1
            map%lm2 (l,m) = lm
            map%lm2l(lm)  = l
            map%lm2m(lm)  = m
            
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
      end do
      
      if ( lm /= map%n_lm ) then
         write(*,"(2(A,I6))") 'Wrong lm=',lm," != n_lm_loc = ", map%n_lm
         call abortRun('Stop run in distribute_theta')
      end if
      if ( lmP /= map%n_lmP ) then
         write(*,"(3(A,I6))") 'Wrong lmP=',lm," != n_lmP_loc = ", map%n_lmP
         call abortRun('Stop run in distribute_theta')
      end if

      do lm=1,map%n_lm
         l=map%lm2l(lm) 
         m=map%lm2m(lm)
         if ( l > 0 .and. l > m ) map%lm2lmS(lm)=map%lm2(l-1,m)
         if ( l < map%l_max )     map%lm2lmA(lm)=map%lm2(l+1,m)
         map%lm2mc(lm) =map%lm2m(lm)/minc + 1
      end do
      do lmP=1,map%n_lmP
         l=map%lmP2l(lmP)
         m=map%lmP2m(lmP)
         if ( l > 0 .and. l > m ) map%lmP2lmPS(lmP)=map%lmP2(l-1,m)
         if ( l < map%l_max+1 )   map%lmP2lmPA(lmP)=map%lmP2(l+1,m)
      end do
      
   end subroutine set_lmmapping_default
   
   !----------------------------------------------------------------------------
   subroutine set_mlmapping(map)
      !   
      !   This function assumes that m_arr is of size (n,2), where n is the 
      !   number of tuples in this mapping
      !   
      !   It makes no assumption about the lengths m_arr.
      !   
      !   Author: Rafael Lago, MPCDF, November 2018
      !
      !-- TODO: double-check if lm2mc is correctly computed (or if it is needed
      !   at all anymore. I don't think it is.)
      !   
      type(ml_mappings), intent(inout) :: map
      integer, parameter :: Invalid_Idx = -1
      
      !-- Local variables
      integer :: m, l, mi, lj, i, j, ml, irank, mlo_idx
      integer :: l_counter, m_counter
      
      map%ml2coord = Invalid_Idx
      map%i2ml     = Invalid_Idx
      map%ml2i     = Invalid_Idx
      map%i2m      = Invalid_Idx
      map%i2l      = Invalid_Idx

      ! Which ranks contain the specified (m,l) tuplet
      do irank=0,n_ranks_mlo-1
         do i=1,n_mlo_array
            m = dist_mlo(irank,i,1)
            l = dist_mlo(irank,i,2)
            if(m>=0 .and. l>=0) map%ml2coord(m,l) = irank
         end do
      end do
      
      ! Maps all local m,l tuplets into a global array of size l_max,l_max
      ! The tuples which do not belong to this rank are marked with Invalid_Idx
      do i = 1,n_mlo_loc
         m  = dist_mlo(coord_mlo,i,1)
         l  = dist_mlo(coord_mlo,i,2)
         
         if (m < 0) cycle
         if (l < 0) cycle
         
         map%i2m(i) = m
         map%i2l(i) = l
         
         ml = map_glbl_st%lm2(l,m)
         map%i2ml(i)   = ml        ! perhaps incorporate this in dist_mlo?
         map%ml2i(m,l) = i
! ! ! ! ! ! ! !          print *, "pairing: ", l, m, i
      end do
      
!       l_counter = 1
!       m_counter = 1
!       do i=0,l_max
!          if (any(dist_mlo(coord_mlo,:,1)==i)) then
!             map%mi2m(m_counter) = i
!             m_counter = m_counter + 1
!          end if
!          if (any(dist_mlo(coord_mlo,:,2)==i)) then
!             map%li2l(l_counter) = i
!             l_counter = l_counter + 1
!          end if
!       end do
!       map%n_li = l_counter - 1
!       map%n_mi = m_counter - 1
      
      
      ! Now this is complicated
      ! I'll map each (mi,lj) tuple into their atual (m,l) tuple as 
      ! well as their index i and other shenenigans. The thing is, the 
      ! j-th l and the i-th m are not the same as the (m,l)!!
!       map%lj2l   => map%helper(0,:,1)
!       map%n_mi   => map%helper(0,:,2)
!       map%milj2m => map%helper(:,:,1)
!       map%milj2i => map%helper(:,:,2)
!       map%mi2m => map%helper(:,0,1)
      lj = 0
      do l=0,l_max
         if (.not. any(dist_mlo(coord_mlo,:,2)==l)) cycle
         lj = lj + 1
         map%lj2l(lj) = l

         mi = 0
         do m=0,l_max
         
            if (map%ml2i(m,l)<0) cycle
            mi = mi + 1
            
            map%milj2m(mi,lj) = m
            map%milj2i(mi,lj) = map%ml2i(m,l)
         end do
         map%n_mi(lj) = mi
      end do
      
   end subroutine set_mlmapping
   
   !----------------------------------------------------------------------------
   subroutine print_mapping(map,name)
      !  
      !   Author: Rafael Lago, MPCDF, June 2018
      !  
      type(mappings), intent(in) :: map
      character(*), intent(in) :: name
      integer :: irank, ierr, count_l, count_m, i
      
      count_m = 0
      count_l = 0
      do i=0,map%l_max
         if (any(map%lm2(i,:) >= 0)) count_m = count_m + 1
      end do
      do i=0,map%m_max
         if (any(map%lm2(:,i) >= 0)) count_l = count_l + 1
      end do
      
      
      if (coord_m == 0 .and. coord_r == 0) then
         print "(' ! ',A,' mapping in rank ', I0, ': ', I0,' l-pts and ',I0,' m-pts (', I0, ' pts)')", name, &
            0, count_m, count_l, map%n_lm
      end if
      call mpi_barrier(mpi_comm_world,ierr)
      do irank=1,n_ranks_m-1
         if (coord_m == irank .and. coord_r == 0) then
            print "(' !                rank ', I0, ': ', I0,' l-pts and ',I0,' m-pts (', I0, ' pts)')", &
               irank, count_m, count_l, map%n_lm
         end if
         call mpi_barrier(mpi_comm_world,ierr)
      end do 
   end subroutine print_mapping
   
end module LMmapping
