!$Id$
!#define OLD_THETA_BLOCKING 1
module blocking
   !--------------------------------------------------
   !  Common block containing blocking information
   !--------------------------------------------------

   use precision_mod
   use parallel_mod, only: nThreads, rank, n_procs, nLMBs_per_rank, &
                           rank_with_l1m0
   use truncation, only: lmP_max, lm_max, l_max, nrp, n_theta_max, &
                         minc, n_r_max
   use logic, only: l_save_out
   use output_data, only: nLF, log_file
   use LMmapping, only: mappings, allocate_mappings,  &
                        allocate_subblocks_mappings,  &
                        subblocks_mappings
   use useful, only: logWrite
 
   implicit none
 
   private
 
   !------------------------------------------------------------------------
   !  Dividing loops over all spherical harmonics into blocks that
   !  contain approx. nChunk harmonics . The number of blocks, nLMBs,   
   !  is a multiple of nThreadsUse (the number of processors used).
   !  Parameter nChunk controls the number (and size) of blocks.
   !  When nThreadUse is large, the size of the blocks may be 
   !  considerably smaller than the chosen nChunk,
   !  since nLMBs must be a multiple of nThreadsUse!
 
   !integer,PARAMETER :: nChunk=512
   !integer :: nThreadsMax
   ! nthreads > 1
   integer, public, pointer :: lm2(:,:),lm2l(:),lm2m(:)
   integer, public, pointer :: lm2mc(:),l2lmAS(:)
   integer, public, pointer :: lm2lmS(:),lm2lmA(:)
 
   integer, public, pointer :: lmP2(:,:),lmP2l(:)
   integer, public, pointer :: lmP2lmPS(:),lmP2lmPA(:)
 
   integer, public, pointer :: lm2lmP(:),lmP2lm(:)
 
   
   type(mappings), public, target :: st_map
   type(mappings), public, target :: lo_map
   !TYPE(mappings),TARGET :: sn_map
 
   !integer :: nLMBsMax
   integer, public :: nLMBs,sizeLMB
 
   integer, public, allocatable :: lmStartB(:),lmStopB(:)
   !integer,PARAMETER :: sizeLMB2max=201
 
   integer, public, pointer :: nLMBs2(:),sizeLMB2(:,:)
   integer, public, pointer :: lm22lm(:,:,:)
   integer, public, pointer :: lm22l(:,:,:)
   integer, public, pointer :: lm22m(:,:,:)
 
   type(subblocks_mappings), public, target :: st_sub_map, lo_sub_map,sn_sub_map
 
   integer, public :: sizeRB
 
 
   !------------------------------------------------------------------------
   !  Following divides loops over points in theta-direction (index ic) into
   !  blocks. Enhances performance by trying to decrease memory access
   !  but is not relevant for SMP parallel processing.
   !  It is tried to divide the theta loop into parts whose data
   !  fit into cache. Thus ideal block sizes (nfs) highly depend on the
   !  used computer.
   !  The value nfs used here has been determined experientally
   !  by Uli Christensen for an IBM SP2. It represents an upper bound.
   !  The real block size sizeThetaB used in the code is determined in
   !  s_prep.f by decreasing the size starting with nfs,
   !  until n_theta_max is a multiple of sizeThetaB. The maximum number of
   !  blocks is limited to 8 here (see s_prep.f).
   !  To find a new blocking for a different computer I suggest to
   !  vary the number of blocks nThetaBs (see s_prep.f) for a fixed
   !  resolution. If the ideal no. of blocks nThetaBsI has been found
   !  this determins the ideal block size:
   !         sizeThetaBI=((n_theta_max-1)/nThetaBsI+1)*n_phi_max
   !  This can than be rescaled for different resolutions n_phi_max:
   !         nfs=sizeThetaBI/n_phi_max+1
   !  The resulting block number can be kept down by multiplying this
   !  with an integer number nBDown, Uli has used K=8 here.
   !  For the memory access this is equivalent to inceasing the
   !  number of blocks.
   !  I dont know exactly why Uli has the additional 16. It may just
   !  decrease the block size to be on the safe side, which is a good
   !  idea, since you dont want to end up with blocks that are just
   !  a little bit larger than the cache.
   !  Thus it seems a good idea to use
   !        nfs=sizeThetaBI/(n_phi_max+nBSave)+1)*nBDown
   !
 
   integer, public :: nfs
   integer, parameter :: sizeThetaBI=284,nBSave=16,nBDown=8
   integer, public :: cacheblock_size_in_B=4096
 
   integer, public :: nThetaBs,sizeThetaB
 
   interface get_theta_blocking
      module procedure get_theta_blocking_cache,get_theta_blocking_OpenMP
   end interface get_theta_blocking
 
   public :: initialize_blocking, get_theta_blocking

contains

   subroutine initialize_blocking

      real(cp) :: load
      integer :: iLoad
      integer :: n
      integer :: LMB_with_l1m0,l1m0,irank

      logical,PARAMETER :: DEBUG_OUTPUT=.false.
      integer :: lm,l,m

      character(len=255) :: message


      call allocate_mappings(st_map,l_max,lm_max,lmP_max)
      call allocate_mappings(lo_map,l_max,lm_max,lmP_max)
      !call allocate_mappings(sn_map,l_max,lm_max,lmP_max)

      if ( ( rank == 0 ) .and. l_save_out ) then
         open(nLF, file=log_file, status='unknown', position='append')
      end if

      if ( rank == 0 ) then
         write(message,*) '! Number of ranks I will use:',n_procs
         call logWrite(message)
      end if

      nLMBs = n_procs
      nLMBs_per_rank = nLMBs/n_procs
      !nThreadsMax = 1
      !PRINT*,"nLMBs first = ",nLMBs
      sizeLMB=(lm_max-1)/nLMBs+1
      !nLMBs = nLMBs - (nLMBs*sizeLMB-lm_max)/sizeLMB

      if ( nLMBs*sizeLMB > lm_max ) then
         write(message,*) '! Uneven load balancing in LM blocks!'
         call logWrite(message)
         load=real(lm_max-(nLMBs-1)*sizeLMB,cp)/sizeLMB
         write(message,*) '! Load percentage of last block:',load*100.0_cp
         call logWrite(message)
         iLoad=int(load)
         if ( iLoad >= 1 ) then
            write(*,*) '! No. of redundant blocks:',iLoad
            nLMBs=nLMBs-iLoad
         end if
      end if
      allocate( lmStartB(nLMBs),lmStopB(nLMBs) )

      !--- Get radial blocking
      if ( mod(n_r_max-1,n_procs) /= 0 ) then
         if ( rank == 0 ) then
            write(*,*) 'Number of MPI ranks has to be multiple of n_r_max-1!'
            write(*,*) 'n_procs :',n_procs
            write(*,*) 'n_r_max-1:',n_r_max-1
         end if
         stop
      end if
      sizeRB=(n_r_max-1)/n_procs

      !--- Calculate lm and ml blocking:
      do n=1,nLMBs
         lmStartB(n)=(n-1)*sizeLMB+1
         lmStopB(n) =min(n*sizeLMB,lm_max)
         if ( lmStopB(n) == lm_max ) exit
      end do

      call get_standard_lm_blocking(st_map,minc)
      !call get_standard_lm_blocking(lo_map,minc)
      if (n_procs <= l_max/2) then
         !better load balancing, but only up to l_max/2
         call get_snake_lm_blocking(lo_map,minc)
         write(message,*) "Using snake ordering."
      else
         call get_lorder_lm_blocking(lo_map,minc)
         write(message,*) "Using lorder ordering."
      end if
      call logWrite(message)

      do n=1,nLMBs
         write(message,*) n,lmStartB(n),lmStopB(n),lmStopB(n)-lmStartB(n)+1
         call logWrite(message)
         if ( lmStopB(n) == lm_max ) exit
      end do

      !--- Get the block (rank+1) with the l1m0 mode
      l1m0 = lo_map%lm2(1,0)
      do n=1,nLMBs
         if ( (l1m0 >= lmStartB(n)) .and. (l1m0 <= lmStopB(n)) ) then
            LMB_with_l1m0=n
            exit
         end if
      end do

      ! which rank does have the LMB with LMB_with_l1m0?
      do irank=0,n_procs-1
         if ((LMB_with_l1m0-1 >= irank*nLMBs_per_rank).and.&
              &(LMB_with_l1m0-1 <= (irank+1)*nLMBs_per_rank-1)) then
            rank_with_l1m0 = irank
         end if
      end do
      write(message,"(2(A,I4))") "rank no ",rank_with_l1m0, &
            " has l1m0 in block ",LMB_with_l1m0
      call logWrite(message)
         
      if (DEBUG_OUTPUT) then
         ! output the lm -> l,m mapping
         if (rank == 0) then
            do lm=1,lm_max
               l=lo_map%lm2l(lm)
               m=lo_map%lm2m(lm)
               write(*,"(A,I5,2(A,I3))") "lm = ",lm," --> l=",l,", m=",m
            end do
         end if
      end if
      
      ! set the standard ordering as default
      lm2(0:,0:) => st_map%lm2
      lm2l(1:lm_max) => st_map%lm2l
      lm2m(1:lm_max) => st_map%lm2m
      lm2mc(1:lm_max)=> st_map%lm2mc
      l2lmAS(0:l_max)=> st_map%l2lmAS
      lm2lmS(1:lm_max) => st_map%lm2lmS
      lm2lmA(1:lm_max) => st_map%lm2lmA
      lmP2(0:,0:) => st_map%lmP2
      lmP2l(1:lmP_max) => st_map%lmP2l
      lmP2lmPS(1:lmP_max) => st_map%lmP2lmPS
      lmP2lmPA(1:lmP_max) => st_map%lmP2lmPA
      lm2lmP(1:lm_max) => st_map%lm2lmP
      lmP2lm(1:lmP_max) => st_map%lmP2lm

      call allocate_subblocks_mappings(st_sub_map,st_map,nLMBs,l_max,lmStartB,lmStopB)
      call allocate_subblocks_mappings(lo_sub_map,lo_map,nLMBs,l_max,lmStartB,lmStopB)
      !call allocate_subblocks_mappings(sn_sub_map,sn_map,nLMBs,l_max,lmStartB,lmStopB)

      !--- Getting lm sub-blocks:
      !PRINT*," ---------------- Making the standard subblocks -------------- "
      call get_subblocks(st_map, st_sub_map) !nLMBs,lm2l,lm2m, nLMBs2,sizeLMB2,lm22lm,lm22l,lm22m)
      !PRINT*," ---------------- Making the lorder subblocks ---------------- "
      call get_subblocks(lo_map, lo_sub_map)
      !PRINT*," ---------------- Making the snake order subblocks ----------- "
      !call get_subblocks(sn_map, sn_sub_map)

      ! default mapping
      nLMBs2(1:nLMBs) => st_sub_map%nLMBs2
      sizeLMB2(1:,1:) => st_sub_map%sizeLMB2
      lm22lm(1:,1:,1:) => st_sub_map%lm22lm
      lm22l(1:,1:,1:) => st_sub_map%lm22l
      lm22m(1:,1:,1:) => st_sub_map%lm22m

      !-- Calculate blocking parameters for blocking loops over theta:

      if (nThreads == 1) then
#ifdef OLD_THETA_BLOCKING    
         nfs=(sizeThetaBI/(n_phi_tot+nBSave)+1) * nBDown
         sizeThetaB=min(n_theta_max,nfs)
         nThetaBs  =n_theta_max/sizeThetaB
         if ( nThetaBs*sizeThetaB /= n_theta_max ) then
            write(*,*)
            write(*,*) '! n_theta_max is not multiple of nfs!'
            write(*,*) '! n_theta_max    =',n_theta_max
            write(*,*) '! nfs            =',nfs
            write(*,*) '! n_theta_max/nfs=',n_theta_max/nfs
            write(*,*) '! Please decrease sizeThetaBI or nBDown in m_blocking.F90!'
            stop
         end if
#else
         call get_theta_blocking_cache(n_theta_max,nrp,cacheblock_size_in_B, &
                                       nThetaBs,sizeThetaB)
         nfs=sizeThetaB
#endif
      else
         call get_theta_blocking_OpenMP(n_theta_max,nThreads, nThetaBs,sizeThetaB)
         nfs=sizeThetaB
      end if

      if ( rank == 0 ) then
         write(*,*) '!-- Blocking information:'
         write(*,*)
         write(*,*) '!    Number of LM-blocks:',nLMBs
         write(*,*) '!    Size   of LM-blocks:',sizeLMB
         write(*,*) '!               nThreads:',nThreads
         write(*,*)
         write(*,*) '! Number of theta blocks:',nThetaBs
         write(*,*) '!   size of theta blocks:',sizeThetaB
         write(*,*) '!       ideal size (nfs):',nfs
         write(nLF,*) '!-- Blocking information:'
         write(nLF,*)
         write(nLF,*) '!    Number of LM-blocks:',nLMBs
         write(nLF,*) '!    Size   of LM-blocks:',sizeLMB
         write(nLF,*) '!               nThreads:',nThreads
         write(nLF,*)
         write(nLF,*) '! Number of theta blocks:',nThetaBs
         write(nLF,*) '!   size of theta blocks:',sizeThetaB
         write(nLF,*) '!       ideal size (nfs):',nfs

         if ( l_save_out ) close(nLF)
      end if

   end subroutine initialize_blocking
!------------------------------------------------------------------------
   subroutine get_subblocks(map,sub_map) 

      !-- Input variables:
      type(mappings),           intent(in) :: map
      type(subblocks_mappings), intent(inout) :: sub_map
      
      !-- Local variables:
      integer :: number_of_blocks
      logical :: lStop
      integer :: size
      integer :: nB2,n,n2,n3
      integer :: help,help1(lm_max),help2(lm_max),help3(lm_max)
      integer :: lm,l,m
      integer :: check(0:l_max,0:l_max)

      logical, parameter :: DEBUG_OUTPUT=.false.

      number_of_blocks=sub_map%nLMBs
      
      check = 0
      lStop=.false.
      size=0
      nB2=0
      do n=1,nLMBs
         sub_map%nLMBs2(n)=1
         lm=lmStartB(n)
         !PRINT*,n,": lm = ",lm
         !------ Start first sub-block:
         sub_map%sizeLMB2(1,n) =1
         sub_map%lm22lm(1,1,n) =lm
         sub_map%lm22l(1,1,n)  =map%lm2l(lm)
         sub_map%lm22m(1,1,n)  =map%lm2m(lm)
         do lm=lmStartB(n)+1,lmStopB(n)
            !write(*,"(4X,A,I4)") "lm = ",lm
            do n2=1,sub_map%nLMBs2(n)
               !write(*,"(8X,A,I4)") "n2 = ",n2
               if ( sub_map%lm22l(1,n2,n) == map%lm2l(lm) ) then
                  !------ Add to old block
                  sub_map%sizeLMB2(n2,n)=sub_map%sizeLMB2(n2,n)+1
                  go to 20
               end if
            end do
            !------ Start new l-block:
            n2 = sub_map%nLMBs2(n)+1
            sub_map%nLMBs2(n)     =n2
            sub_map%sizeLMB2(n2,n)=1
20          continue
            sub_map%lm22lm(sub_map%sizeLMB2(n2,n),n2,n)=lm
            sub_map%lm22l(sub_map%sizeLMB2(n2,n),n2,n) =map%lm2l(lm)
            sub_map%lm22m(sub_map%sizeLMB2(n2,n),n2,n) =map%lm2m(lm)
         end do

         !------ Resort:
         if ( sub_map%nLMBs2(n) > 1 ) then
            do n2=1,sub_map%nLMBs2(n)
               do n3=n2+1,sub_map%nLMBs2(n)
                  if  ( sub_map%lm22m(1,n2,n) > sub_map%lm22m(1,n3,n) ) then
                     help=sub_map%sizeLMB2(n2,n)
                     do lm=1,help
                        help1(lm)=sub_map%lm22l(lm,n2,n)
                        help2(lm)=sub_map%lm22m(lm,n2,n)
                        help3(lm)=sub_map%lm22lm(lm,n2,n)
                     end do
                     sub_map%sizeLMB2(n2,n)=sub_map%sizeLMB2(n3,n)
                     do lm=1,sub_map%sizeLMB2(n2,n)
                        sub_map%lm22l(lm,n2,n) =sub_map%lm22l(lm,n3,n)
                        sub_map%lm22m(lm,n2,n) =sub_map%lm22m(lm,n3,n)
                        sub_map%lm22lm(lm,n2,n)=sub_map%lm22lm(lm,n3,n)
                     end do
                     sub_map%sizeLMB2(n3,n)=help
                     do lm=1,help
                        sub_map%lm22l(lm,n3,n) =help1(lm)
                        sub_map%lm22m(lm,n3,n) =help2(lm)
                        sub_map%lm22lm(lm,n3,n)=help3(lm)
                     end do
                  end if
               end do
            end do
         end if

         nB2=nB2+sub_map%nLMBs2(n)
         do n2=1,sub_map%nLMBs2(n)
            if ( sub_map%sizeLMB2(n2,n) > sub_map%sizeLMB2max ) then
               lStop=.true.
               size=max(size,sub_map%sizeLMB2(n2,n))
            end if
            do n3=1,sub_map%sizeLMB2(n2,n)
               l=sub_map%lm22l(n3,n2,n)
               m=sub_map%lm22m(n3,n2,n)
               check(l,m)=check(l,m)+1
            end do
            !             write(99,*) n,n2,sub_map%sizeLMB2(n2,n)
         end do
         if (DEBUG_OUTPUT) then
            if (rank == 0) then
               write(*,"(4X,2(A,I4))") "Subblocks of Block ",n,"/",nLMBs
               do n2=1,sub_map%nLMBs2(n)
                  write(*,"(8X,3(A,I4))") "subblock no. ",n2,", of ",&
                       & sub_map%nLMBs2(n)," with size ",sub_map%sizeLMB2(n2,n)
                  do n3=1,sub_map%sizeLMB2(n2,n)
                     write(*,"(10X,A,I4,A,I6,2I4)") "local lm is ",n3,      &
                          &" translates into global lm,l,m : ",             &
                          & sub_map%lm22lm(n3,n2,n),sub_map%lm22l(n3,n2,n), &
                          & sub_map%lm22m(n3,n2,n)
                  end do
               end do
            end if
         end if

      end do

      if ( lStop ) then
         write(*,*) '! Increase sizeLMB2max in m_blocking.F90!'
         write(*,*) '! to at least:',size
         stop
      end if

      do m=0,l_max,minc
         do l=m,l_max
            if ( check(l,m) == 0 ) then
               write(*,*) 'Warning, forgotten l,m:',l,m,map%lm2(l,m)
               stop
            else if ( check(l,m) > 1 ) then
               write(*,*) 'Warning, too much l,m:',l,m,check(l,m)
               stop
            end if
         end do
      end do
   end subroutine get_subblocks
!------------------------------------------------------------------------
   subroutine get_standard_lm_blocking(map,minc)

      type(mappings), intent(inout) :: map
      integer,        intent(in) :: minc
      
      ! Local variables
      integer :: m,l,lm,lmP,mc
      integer :: lmP2m(lmP_max)
      
      do m=0,map%l_max
         do l=m,map%l_max
            map%lm2(l,m)  =-1
            map%lmP2(l,m) =-1
            !check(l,m)=0
         end do
         l=map%l_max+1
         map%lmP2(l,m)=-1
      end do

      lm =0
      lmP=0
      mc =0
      do m=0,map%l_max,minc
         mc=mc+1
         !m2mc(m)=mc
         do l=m,map%l_max
            lm         =lm+1
            map%lm2l(lm)   =l
            map%lm2m(lm)   =m
            map%lm2mc(lm)  =mc
            map%lm2(l,m)   =lm
            if ( m == 0 ) map%l2lmAS(l)=lm
            lmP        =lmP+1
            map%lmP2l(lmP) = l
            lmP2m(lmP) = m
            map%lmP2(l,m)  =lmP
            !if ( m == 0 ) l2lmPAS(l)=lmP
            map%lm2lmP(lm) =lmP
            map%lmP2lm(lmP)=lm
         end do
         l=map%l_max+1    ! Extra l for lmP
         lmP=lmP+1
         map%lmP2l(lmP) =l
         lmP2m(lmP) = m
         map%lmP2(l,m)  =lmP
         !if ( m == 0 ) l2lmPAS(l)=lmP
         map%lmP2lm(lmP)=-1
      end do
      if ( lm /= map%lm_max ) then
         write(*,"(2(A,I6))") 'Wrong lm=',lm," != map%lm_max = ",map%lm_max
         stop
      end if
      if ( lmP /= map%lmP_max ) then
         write(*,*) 'Wrong lmP!'
         stop
      end if
      do lm=1,map%lm_max
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
      do lmP=1,map%lmP_max
         l=map%lmP2l(lmP)
         m=lmP2m(lmP)
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

   end subroutine get_standard_lm_blocking
!------------------------------------------------------------------------
   subroutine get_lorder_lm_blocking(map,minc)

      type(mappings), intent(inout) :: map
      integer,        intent(in) :: minc
      
      ! Local variables
      integer :: m,l,lm,lmP,mc
      integer :: lmP2m(lmP_max)
      
      do m=0,map%l_max
         do l=m,map%l_max
            map%lm2(l,m)  =-1
            map%lmP2(l,m) =-1
            !check(l,m)=0
         end do
         l=map%l_max+1
         map%lmP2(l,m)=-1
      end do

      lm =0
      lmP=0
      do l=0,map%l_max
         mc =0
         ! set l2lmAS for m==0
         map%l2lmAS(l)=lm
         do m=0,l,minc
            mc=mc+1

            lm         =lm+1
            map%lm2l(lm)   =l
            map%lm2m(lm)   =m
            map%lm2mc(lm)  =mc
            map%lm2(l,m)   =lm

            lmP        =lmP+1
            map%lmP2l(lmP) = l
            lmP2m(lmP) = m
            map%lmP2(l,m)  =lmP
            !if ( m == 0 ) l2lmPAS(l)=lmP
            map%lm2lmP(lm) =lmP
            map%lmP2lm(lmP)=lm
         end do
      end do
      l=map%l_max+1    ! Extra l for lmP
      mc =0
      do m=0,map%l_max,minc
         mc=mc+1

         lmP=lmP+1
         map%lmP2l(lmP) =l
         lmP2m(lmP) = m
         map%lmP2(l,m)  =lmP
         map%lmP2lm(lmP)=-1
      end do

      if ( lm /= map%lm_max ) then
         write(*,"(2(A,I6))") 'get_lorder_lm_blocking: Wrong lm = ',lm, &
                              " != map%lm_max = ",map%lm_max
         stop
      end if
      if ( lmP /= map%lmP_max ) then
         write(*,*) 'Wrong lmP!'
         stop
      end if
      do lm=1,map%lm_max
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
      do lmP=1,map%lmP_max
         l=map%lmP2l(lmP)
         m=lmP2m(lmP)
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

   end subroutine get_lorder_lm_blocking
!------------------------------------------------------------------------
   subroutine get_snake_lm_blocking(map,minc)

      type(mappings), intent(inout) :: map
      integer,        intent(in) :: minc
      
      ! Local variables
      integer :: l,proc,lm,m,i_l,lmP,mc
      logical :: Ascending
      integer :: lmP2m(map%lmP_max)
      integer :: l_list(0:n_procs-1,map%l_max+1)
      integer :: l_counter(0:n_procs-1)
      integer :: temp_l_counter,l0proc,pc,src_proc,temp
      integer :: temp_l_list(map%l_max+1)

      logical, parameter :: DEBUG_OUTPUT=.false.
      ! First we loop over all l values and jump for each
      ! new l value to the next process in a snake like fashion.
      proc=0
      Ascending=.true.
      l_counter=1
      do l=map%l_max,0,-1
         ! this l block is distributed to the actual proc
         l_list(proc,l_counter(proc))=l
         !write(*,"(A,3I3)") "l,l_list,l_counter=",l,l_list(proc,l_counter(proc)),l_counter(proc)
         l_counter(proc) = l_counter(proc)+1
         if (l == 0) l0proc=proc
         ! now determine on which proc to put the next l value
         if (Ascending) then
            if (proc < n_procs-1) then
               proc=proc+1
            else if (proc == n_procs-1) then
               Ascending=.false.
            end if
         else
            if (proc > 0) then
               proc=proc-1
            else if (proc == 0) then
               Ascending=.true.
            end if
         end if
      end do

      if (DEBUG_OUTPUT) then
         do proc=0,n_procs-1
            if (proc == l0proc) then
               write(*,"(A,I4,A)") "==== proc ",proc," has l=0 ===="
            else
               write(*,"(A,I4,A)") "---- proc ",proc," ----"
            end if
            do i_l=1,l_counter(proc)-1
               write(*,"(I4,2X)",advance="no") l_list(proc,i_l)
            end do
            write(*,"(A)") ""
         end do
      end if

      ! Now distribution is as equal as possible. We rotate the distribution
      ! now to have the l0proc as first process. 
      if (l0proc /= 0) then
         temp_l_list=l_list(0,:)
         temp_l_counter=l_counter(0)
         pc = 0
         do while (.true.)
            src_proc=modulo(l0proc+pc,n_procs)
            if (src_proc /= 0) then
               l_list(pc,:)=l_list(src_proc,:)
               l_counter(pc)=l_counter(src_proc)
            else
               l_list(pc,:)=temp_l_list
               l_counter(pc)=temp_l_counter
               EXIT
            end if
            ! now we can overwrite src_proc
            pc=src_proc
         end do
         l0proc=0
      end if

      ! Last step in preparation is to put the l=0 on process 0
      ! as the first l in the list
      do i_l=1,l_counter(0)-1
         if (l_list(0,i_l) == 0) then
            !write(*,"(A,I3)") "i_l = ",i_l
            temp=l_list(0,1)
            l_list(0,1)=l_list(0,i_l)
            l_list(0,i_l)=temp
            EXIT
         end if
      end do

      if (DEBUG_OUTPUT) then
         write(*,"(A)") "Ordering after the l0proc reordering:"
         do proc=0,n_procs-1
            if (proc == l0proc) then
               write(*,"(A,I4,A)") "==== proc ",proc," has l=0 ===="
            else
               write(*,"(A,I4,A)") "---- proc ",proc," ----"
            end if
            do i_l=1,l_counter(proc)-1
               write(*,"(I4,2X)",advance="no") l_list(proc,i_l)
            end do
            write(*,"(A)") ""
         end do
      end if

      lm=1
      lmP=1
      do proc=0,n_procs-1
         lmStartB(proc+1)=lm
         do i_l=1,l_counter(proc)-1
            l=l_list(proc,i_l)
            mc = 0
            !write(*,"(3I3)") i_l,proc,l
            do m=0,l,minc
               mc = mc+1
               map%lm2(l,m)=lm
               map%lm2l(lm)=l
               map%lm2m(lm)=m
               map%lm2mc(lm)=mc

               map%lmP2(l,m)=lmP
               map%lmP2l(lmP)=l
               lmP2m(lmP) = m
               map%lm2lmP(lm)=lmP
               map%lmP2lm(lmP)=lm

               lm = lm+1
               lmP = lmP+1
            end do
         end do
         lmStopB(proc+1)=lm-1
         write(*,"(I3,2(A,I6))") proc,": lmStartB=",lmStartB(proc+1), &
                                 ", lmStopB=",lmStopB(proc+1)
      end do
      
      if ( lm-1 /= map%lm_max ) then
         write(*,"(2(A,I6))") 'get_snake_lm_blocking: Wrong lm-1 = ',lm-1,&
              & " != map%lm_max = ",map%lm_max
         stop
      end if

      l=map%l_max+1    ! Extra l for lmP
      mc =0
      do m=0,map%l_max,minc
         mc=mc+1

         map%lmP2l(lmP) =l
         lmP2m(lmP) = m
         map%lmP2(l,m)  =lmP
         map%lmP2lm(lmP)=-1
         lmP=lmP+1
      end do

      if ( lmP-1 /= map%lmP_max ) then
         write(*,*) 'Wrong lmP!'
         stop
      end if

      do lm=1,map%lm_max
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
      do lmP=1,map%lmP_max
         l=map%lmP2l(lmP)
         m=lmP2m(lmP)
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

   end subroutine get_snake_lm_blocking
!------------------------------------------------------------------------
   subroutine get_theta_blocking_cache(n_theta_max,nrp,      &
                                       cacheblock_size_in_B, &
                                       nThetaBs, sizeThetaB)

      integer, intent(in) :: n_theta_max,nrp,cacheblock_size_in_B
      integer, intent(out) :: nThetaBs, sizeThetaB
    
      integer :: best_s,s,memory_size,min_s
      
      best_s=0
      min_s = 0
      ! The size of the theta blocks must be dividable by 4
      ! due to the algorithms in the legTF routines.
      do s=4,n_theta_max,4
         if ( modulo(n_theta_max,s)==0 ) then
            ! candidate found
            if (min_s == 0) min_s=s
            nThetaBs=n_theta_max/s
            memory_size=s*nrp*8
            if ( cacheblock_size_in_b/real(memory_size)  >=  1.0 ) then
               best_s=s
            else if ( cacheblock_size_in_B/memory_size  ==  0 ) then
               exit
            end if
         end if
      end do
      if ( best_s /= 0 ) then
         sizeThetaB=best_s
      else
         sizeThetaB=min_s
      end if
      nThetaBs=n_theta_max/sizeThetaB

   end subroutine get_theta_blocking_cache
!------------------------------------------------------------------------
   subroutine get_theta_blocking_OpenMP(n_theta_max, nThreads, nThetaBs, sizeThetaB)
      !
      !  This routine determines the number of theta blocks and the
      !  blocksize with respect to the number of threads.
      !
      integer, intent(in) :: n_theta_max,nThreads
      integer, intent(out) :: nThetaBs, sizeThetaB
    
      integer :: best_s,s,min_s
      
      best_s=0
      min_s = 0
      ! The size of the theta blocks must be dividable by 4
      ! due to the algorithms in the legTF routines.
      do s=4,n_theta_max,4
         if ( modulo(n_theta_max,s)==0 ) then
            ! candidate found
            if (min_s == 0) min_s=s
            nThetaBs=n_theta_max/s
            !write(*,"(3(A,I3))") "Testing s=",s,", nThreads=",nThreads,", &
            !     &              nThetaBs = ",nThetaBs
    
            if ( modulo(nThetaBs,nThreads) == 0 ) then
               best_s=s
            elseif (nThetaBs/nThreads  ==  0) then
               EXIT
            end if
         end if
      end do
      if (best_s /= 0) then
         sizeThetaB=best_s
      else
         sizeThetaB=min_s
      end if
      nThetaBs=n_theta_max/sizeThetaB

   end subroutine get_theta_blocking_OpenMP
!------------------------------------------------------------------------
end module blocking
