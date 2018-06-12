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
 
   !-- ????
   !   
   !   
   type, public :: subblocks_mappings
      integer :: nLMBs,l_max,m_max,sizeLMB2max
      integer, allocatable :: nLMBs2(:)
      integer, allocatable :: sizeLMB2(:,:)
      integer, allocatable :: lm22lm(:,:,:)
      integer, allocatable :: lm22l(:,:,:)
      integer, allocatable :: lm22m(:,:,:)
   end type subblocks_mappings
 
   type(mappings) :: dist_map, lo_dist_map
   
contains
   
   !----------------------------------------------------------------------------
   subroutine initialize_mapping
      !   
      !   Author: Rafael Lago, MPCDF, April 2018
      !   
      
      integer :: i, iRstart, iRstop, iRremaining, iThetastart, iThetastop
      integer(lip) :: local_bytes_used
      

      local_bytes_used = bytes_allocated
!       call allocate_mappings(st_map,l_max,lm_max,lmP_max,l_axi)
!       call allocate_mappings(lo_map,l_max,lm_max,lmP_max,l_axi)

      call allocate_mappings(dist_map,   n_lm,n_lmP,l_axi)
      call allocate_mappings(lo_dist_map,n_lm,n_lmP,l_axi)
      
      call set_lmmapping_default(dist_map,    dist_m(coord_m,1:))
      call set_lmmapping_default(lo_dist_map, dist_m(coord_m,1:))
      
      call print_mapping(dist_map, 'dist')

      local_bytes_used = bytes_allocated-local_bytes_used
      call memWrite('LMmapping.f90', local_bytes_used)
   end subroutine initialize_mapping
   
   !----------------------------------------------------------------------------
   subroutine finalize_mapping

      call deallocate_mappings(dist_map)

   end subroutine finalize_mapping
   
   !----------------------------------------------------------------------------
   subroutine allocate_mappings(self, in_n_lm, in_n_lmP, l_axi)
      !   
      !   Allocates the mapping objects. 
      !   
      !   The l_max is taken from geometry module. It is supposed to represent 
      !   the global number of l points, not local. 
      !   in_n_lm is the number of local lm points. Similarly, in_n_lmP is the 
      !   same quantity, for l_max+1.
      !   
      !   Author: Rafael Lago, MPCDF, April 2018
      !   

      type(mappings), intent(inout) :: self
      integer, intent(in) :: in_n_lm, in_n_lmP
      logical, intent(in) :: l_axi

      self%l_max = l_max
      if ( .not. l_axi ) then
         self%m_max = l_max
      else
         self%m_max = 0
      end if
      self%lm_max   = in_n_lm   ! deprecated; to be later replaced by %in_n_lm
      self%lmP_max  = in_n_lmP  ! deprecated; to be later replaced by %in_n_lmP
      self%n_lm  = in_n_lm
      self%n_lmP = in_n_lmP

      allocate( self%lm2(0:l_max,0:l_max),self%lm2l(in_n_lm),self%lm2m(in_n_lm))
      allocate( self%lm2mc(in_n_lm),self%l2lmAS(0:l_max) )
      allocate( self%lm2lmS(in_n_lm),self%lm2lmA(in_n_lm) )
      bytes_allocated = bytes_allocated + &
                        ((l_max+1)*(l_max+1)+5*in_n_lm+l_max+1)*SIZEOF_INTEGER

      allocate( self%lmP2(0:l_max+1,0:l_max+1),self%lmP2l(in_n_lmP) )
      allocate( self%lmP2m(in_n_lmP) )
      allocate( self%lmP2lmPS(in_n_lmP),self%lmP2lmPA(in_n_lmP) )
      allocate( self%lm2lmP(in_n_lm),self%lmP2lm(in_n_lmP) )
      bytes_allocated = bytes_allocated + &
                        ((l_max+2)*(l_max+2)+5*in_n_lmP+in_n_lm)*SIZEOF_INTEGER

   end subroutine allocate_mappings
   
   !----------------------------------------------------------------------------
   subroutine deallocate_mappings(self)

      type(mappings) :: self

      deallocate( self%lm2, self%lm2l, self%lm2m, self%lm2mc, self%l2lmAS )
      deallocate( self%lm2lmS, self%lm2lmA, self%lmP2, self%lmP2l )
      deallocate( self%lmP2m, self%lmP2lmPS, self%lmP2lmPA, self%lm2lmP )
      deallocate( self%lmP2lm )

   end subroutine deallocate_mappings
   
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
      !   Remark: if you pass m_arr = (/0,1,2,...,m_max/) this function does 
      !      exactly the same distribution as the global lm map for minc=1.
      !   
      !   Author: Rafael Lago, MPCDF, April 2018
      !   
      !-- TODO: this function can more or less easily be adapted to support 
      !      intervals of l, so, for instance, from l_l to u_l. Just keep this 
      !      in mind in case we need it in the future
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
      
      map%lm2mc = -5 ! Lago@180405: I have no idea what this one is for. It needs to be completely implemented
      
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
      end do
      
      if ( lm /= map%n_lm ) then
         write(*,"(2(A,I6))") 'Wrong lm=',lm," != n_lm = ", map%n_lm
         call abortRun('Stop run in distribute_theta')
      end if
      if ( lmP /= map%n_lmP ) then
         write(*,"(3(A,I6))") 'Wrong lmP=',lm," != n_lmP = ", map%n_lmP
         call abortRun('Stop run in distribute_theta')
      end if

      do lm=1,map%n_lm
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
      do lmP=1,map%n_lmP
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
      
   end subroutine set_lmmapping_default
   
   
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
      do i=0,map%m_max
         if (any(map%lm2(:,i) >= 0)) count_m = count_m + 1
      end do
      do i=0,map%l_max
         if (any(map%lm2(i,:) >= 0)) count_l = count_l + 1
      end do
      
      
      if (coord_m == 0 .and. coord_r == 0) then
         print "(' ! ',A,' mapping in rank ', I0, ': ', I0,' l pts and ',I0,'m pts (', I0, ' pts)')", name, &
            0, count_m, count_l, map%n_lm
      end if
      
      call mpi_barrier(mpi_comm_world,ierr)
      do irank=1,n_ranks_m-1
         if (coord_m == irank .and. coord_r == 0) then
            print "(' !                rank ', I0, ': ', I0,' l pts and ',I0,'m pts (', I0, ' pts)')", &
               irank, count_m, count_l, map%n_lm
         end if
         call mpi_barrier(mpi_comm_world,ierr)
      end do 
   end subroutine print_mapping
   
end module LMmapping
