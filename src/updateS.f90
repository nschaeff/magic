#include "perflib_preproc.cpp"
module updateS_mod

   use omp_lib
   use precision_mod
   use mem_alloc, only: bytes_allocated
   use truncation, only: n_r_max, lm_max, l_max
   use radial_data, only: n_r_cmb, n_r_icb
   use radial_functions, only: orho1, or1, or2, beta, dentropy0, rscheme_oc,  &
       &                       kappa, dLkappa, dLtemp0, temp0
   use physical_parameters, only: opr, kbots, ktops
   use num_param, only: alpha
   use init_fields, only: tops,bots
   use blocking, only: nLMBs,st_map,lo_map,lo_sub_map,lmStartB,lmStopB
   use horizontal_data, only: dLh,hdif_S
   use logic, only: l_update_s, l_anelastic_liquid
   use LMLoop_data, only: llm, ulm
   use parallel_mod, only: rank,chunksize
   use algebra, only: prepare_mat, solve_mat
   use radial_der, only: get_ddr, get_dr
   use fields, only:  work_LMloc
   use constants, only: zero, one, two
   use useful, only: abortRun

   implicit none

   private

   !-- Local variables
   complex(cp), allocatable :: rhs1(:,:,:)
   real(cp), allocatable :: s0Mat(:,:)     ! for l=m=0  
   real(cp), allocatable :: sMat(:,:,:)
   integer, allocatable :: s0Pivot(:)
   integer, allocatable :: sPivot(:,:)
#ifdef WITH_PRECOND_S
   real(cp), allocatable :: sMat_fac(:,:)
#endif
#ifdef WITH_PRECOND_S0
   real(cp), allocatable :: s0Mat_fac(:)
#endif
   logical, public, allocatable :: lSmat(:)

   integer :: maxThreads

   public :: initialize_updateS, updateS, updateS_ala, finalize_updateS

contains

   subroutine initialize_updateS

      allocate( s0Mat(n_r_max,n_r_max) )      ! for l=m=0  
      allocate( sMat(n_r_max,n_r_max,l_max) )
      bytes_allocated = bytes_allocated+(n_r_max*n_r_max*(1+l_max))* &
      &                 SIZEOF_DEF_REAL
      allocate( s0Pivot(n_r_max) )
      allocate( sPivot(n_r_max,l_max) )
      bytes_allocated = bytes_allocated+(n_r_max+n_r_max*l_max)*SIZEOF_INTEGER
#ifdef WITH_PRECOND_S
      allocate(sMat_fac(n_r_max,l_max))
      bytes_allocated = bytes_allocated+n_r_max*l_max*SIZEOF_DEF_REAL
#endif
#ifdef WITH_PRECOND_S0
      allocate(s0Mat_fac(n_r_max))
      bytes_allocated = bytes_allocated+n_r_max*SIZEOF_DEF_REAL
#endif
      allocate( lSmat(0:l_max) )
      bytes_allocated = bytes_allocated+(l_max+1)*SIZEOF_LOGICAL

#ifdef WITHOMP
      maxThreads=omp_get_max_threads()
#else
      maxThreads=1
#endif
      allocate( rhs1(n_r_max,lo_sub_map%sizeLMB2max,0:maxThreads-1) )
      bytes_allocated = bytes_allocated + n_r_max*lo_sub_map%sizeLMB2max*&
      &                 maxThreads*SIZEOF_DEF_COMPLEX

   end subroutine initialize_updateS
!------------------------------------------------------------------------------
   subroutine finalize_updateS

      deallocate( s0Mat, sMat, s0Pivot, sPivot, lSmat )
#ifdef WITH_PRECOND_S
      deallocate( sMat_fac )
#endif
#ifdef WITH_PRECOND_S0
      deallocate( s0Mat_fac )
#endif
      deallocate( rhs1 )
  
   end subroutine finalize_updateS
!------------------------------------------------------------------------------
   subroutine updateS(s,ds,w,dVSrLM,dsdt,dsdtLast,w1,coex,dt,nLMB)
      !
      !  updates the entropy field s and its radial derivatives
      !  adds explicit part to time derivatives of s
      !

      !-- Input of variables:
      real(cp),    intent(in) :: w1        ! weight for time step !
      real(cp),    intent(in) :: coex      ! factor depending on alpha
      real(cp),    intent(in) :: dt        ! time step
      integer,     intent(in) :: nLMB
      complex(cp), intent(in) :: w(llm:ulm,n_r_max)
      complex(cp), intent(inout) :: dVSrLM(llm:ulm,n_r_max)

      !-- Input/output of scalar fields:
      complex(cp), intent(inout) :: s(llm:ulm,n_r_max)
      complex(cp), intent(inout) :: dsdt(llm:ulm,n_r_max)
      complex(cp), intent(inout) :: dsdtLast(llm:ulm,n_r_max)
      !-- Output: udpated s,ds,dsdtLast
      complex(cp), intent(out) :: ds(llm:ulm,n_r_max)

      !-- Local variables:
      real(cp) :: w2            ! weight of second time step
      real(cp) :: O_dt
      integer :: l1,m1              ! degree and order
      integer :: lm1,lmB,lm         ! position of (l,m) in array
      integer :: lmStart,lmStop
      integer :: nLMB2
      integer :: nR                 ! counts radial grid points
      integer :: n_r_out             ! counts cheb modes
      real(cp) ::  rhs(n_r_max) ! real RHS for l=m=0

      integer, pointer :: nLMBs2(:),lm2l(:),lm2m(:)
      integer, pointer :: sizeLMB2(:,:),lm2(:,:)
      integer, pointer :: lm22lm(:,:,:),lm22l(:,:,:),lm22m(:,:,:)

      integer :: threadid,nThreads,iThread,all_lms,per_thread,start_lm,stop_lm
      integer :: iChunk,nChunks,size_of_last_chunk,lmB0

      if ( .not. l_update_s ) return

      nLMBs2(1:nLMBs) => lo_sub_map%nLMBs2
      sizeLMB2(1:,1:) => lo_sub_map%sizeLMB2
      lm22lm(1:,1:,1:) => lo_sub_map%lm22lm
      lm22l(1:,1:,1:) => lo_sub_map%lm22l
      lm22m(1:,1:,1:) => lo_sub_map%lm22m
      lm2(0:,0:) => lo_map%lm2
      lm2l(1:lm_max) => lo_map%lm2l
      lm2m(1:lm_max) => lo_map%lm2m

      lmStart     =lmStartB(nLMB)
      lmStop      =lmStopB(nLMB)
      w2  =one-w1
      O_dt=one/dt


      !PERFON('upS_fin')
      !$OMP PARALLEL  &
      !$OMP private(iThread,start_lm,stop_lm,nR,lm)       &
      !$OMP shared(all_lms,per_thread,lmStart,lmStop)     &
      !$OMP shared(dVSrLM,dsdt,orho1,or2,dentropy0,w,dLh) &
      !$OMP shared(n_r_max,rscheme_oc,work_LMloc,nThreads,llm,ulm)
      !$OMP SINGLE
#ifdef WITHOMP
      nThreads=omp_get_num_threads()
#else
      nThreads=1
#endif
      !-- Get radial derivatives of s: work_LMloc,dsdtLast used as work arrays
      all_lms=lmStop-lmStart+1
      per_thread=all_lms/nThreads
      !$OMP END SINGLE
      !$OMP BARRIER
      !$OMP DO
      do iThread=0,nThreads-1
         start_lm=lmStart+iThread*per_thread
         stop_lm = start_lm+per_thread-1
         if (iThread == nThreads-1) stop_lm=lmStop

         !--- Finish calculation of dsdt:
         call get_dr( dVSrLM,work_LMloc,ulm-llm+1,start_lm-llm+1,  &
              &       stop_lm-llm+1,n_r_max,rscheme_oc, nocopy=.true. )
      end do
      !$OMP end do

      !$OMP DO
      do nR=1,n_r_max
         do lm=lmStart,lmStop
            dsdt(lm,nR)=orho1(nR)*(dsdt(lm,nR)-or2(nR)*work_LMloc(lm,nR)- &
            &           dLh(st_map%lm2(lm2l(lm),lm2m(lm)))*or2(nR)*       &
            &           dentropy0(nR)*w(lm,nR))
         end do
      end do
      !$OMP end do
      !$OMP END PARALLEL
      !PERFOFF

      ! one subblock is linked to one l value and needs therefore once the matrix
      !$OMP PARALLEL default(shared)
      !$OMP SINGLE
      do nLMB2=1,nLMBs2(nLMB)
         ! this inner loop is in principle over the m values which belong to the
         ! l value
         !$OMP TASK default(shared) &
         !$OMP firstprivate(nLMB2) &
         !$OMP private(lm,lm1,l1,m1,lmB,threadid) &
         !$OMP private(nChunks,size_of_last_chunk,iChunk)
         nChunks = (sizeLMB2(nLMB2,nLMB)+chunksize-1)/chunksize
         size_of_last_chunk = chunksize + (sizeLMB2(nLMB2,nLMB)-nChunks*chunksize)

         ! This task treats one l given by l1
         l1=lm22l(1,nLMB2,nLMB)
         !write(*,"(3(A,I3),A)") "Launching task for nLMB2=",nLMB2," (l=",l1,") and scheduling ",nChunks," subtasks."

         if ( l1 == 0 ) then
            if ( .not. lSmat(l1) ) then
#ifdef WITH_PRECOND_S0
               call get_s0Mat(dt,s0Mat,s0Pivot,s0Mat_fac)
#else
               call get_s0Mat(dt,s0Mat,s0Pivot)
#endif
               lSmat(l1)=.true.
            end if
         else
            if ( .not. lSmat(l1) ) then
#ifdef WITH_PRECOND_S
               call get_sMat(dt,l1,hdif_S(st_map%lm2(l1,0)), &
                    &        sMat(1,1,l1),sPivot(1,l1),sMat_fac(1,l1))
#else
               call get_sMat(dt,l1,hdif_S(st_map%lm2(l1,0)), &
                    &        sMat(1,1,l1),sPivot(1,l1))
#endif
               lSmat(l1)=.true.
               !write(*,"(A,I3,ES22.14)") "sMat: ",l1,SUM( sMat(:,:,l1) )
            end if
          end if

         do iChunk=1,nChunks
            !$OMP TASK default(shared) &
            !$OMP firstprivate(iChunk) &
            !$OMP private(lmB0,lmB,lm,lm1,m1,nR,n_r_out) &
            !$OMP private(threadid)
#ifdef WITHOMP
            threadid = omp_get_thread_num()
#else
            threadid = 0
#endif
            lmB0=(iChunk-1)*chunksize
            lmB=lmB0

            do lm=lmB0+1,min(iChunk*chunksize,sizeLMB2(nLMB2,nLMB))
               !do lm=1,sizeLMB2(nLMB2,nLMB)
               lm1=lm22lm(lm,nLMB2,nLMB)
               !l1 =lm22l(lm,nLMB2,nLMB)
               m1 =lm22m(lm,nLMB2,nLMB)

               if ( l1 == 0 ) then
                  rhs(1)=      real(tops(0,0))
                  rhs(n_r_max)=real(bots(0,0))
                  do nR=2,n_r_max-1
                     rhs(nR)=real(s(lm1,nR))*O_dt+w1*real(dsdt(lm1,nR))  + &
                     &       w2*real(dsdtLast(lm1,nR))
                  end do

#ifdef WITH_PRECOND_S0
                  rhs = s0Mat_fac*rhs
#endif

                  call solve_mat(s0Mat,n_r_max,n_r_max,s0Pivot,rhs)

               else ! l1  /=  0
                  lmB=lmB+1

                  rhs1(1,lmB,threadid)=      tops(l1,m1)
                  rhs1(n_r_max,lmB,threadid)=bots(l1,m1)
#ifdef WITH_PRECOND_S
                  rhs1(1,lmB,threadid)=      sMat_fac(1,l1)*rhs1(1,lmB,threadid)
                  rhs1(n_r_max,lmB,threadid)=sMat_fac(1,l1)*rhs1(n_r_max,lmB,threadid)
#endif
                  do nR=2,n_r_max-1
                     rhs1(nR,lmB,threadid)=s(lm1,nR)*O_dt+w1*dsdt(lm1,nR)  &
                     &                     +w2*dsdtLast(lm1,nR)
#ifdef WITH_PRECOND_S
                     rhs1(nR,lmB,threadid) = sMat_fac(nR,l1)*rhs1(nR,lmB,threadid)
#endif
                  end do
               end if
            end do
            !PERFOFF

            !PERFON('upS_sol')
            if ( lmB  >  lmB0 ) then
               call solve_mat(sMat(:,:,l1),n_r_max,n_r_max, &
                    &         sPivot(:,l1),rhs1(:,lmB0+1:lmB,threadid),lmB-lmB0)
            end if
            !PERFOFF

            lmB=lmB0
            !PERFON('upS_af')
            do lm=lmB0+1,min(iChunk*chunksize,sizeLMB2(nLMB2,nLMB))
             !do lm=1,sizeLMB2(nLMB2,nLMB)
               lm1=lm22lm(lm,nLMB2,nLMB)
               !l1 =lm22l(lm,nLMB2,nLMB)
               m1 =lm22m(lm,nLMB2,nLMB)
               if ( l1 == 0 ) then
                  do n_r_out=1,rscheme_oc%n_max
                     s(lm1,n_r_out)=rhs(n_r_out)
                  end do
               else
                  lmB=lmB+1
                  if ( m1 > 0 ) then
                     do n_r_out=1,rscheme_oc%n_max
                        s(lm1,n_r_out)=rhs1(n_r_out,lmB,threadid)
                     end do
                  else
                     do n_r_out=1,rscheme_oc%n_max
                        s(lm1,n_r_out)= cmplx(real(rhs1(n_r_out,lmB,threadid)), &
                        &                    0.0_cp,kind=cp)
                     end do
                  end if
               end if
            end do
            !PERFOFF
            !$OMP END TASK
         end do
         !$OMP END TASK
      end do     ! loop over lm blocks
      !$OMP END SINGLE
      !$OMP END PARALLEL

      !write(*,"(A,2ES22.12)") "s after = ",SUM(s)
      !-- set cheb modes > rscheme_oc%n_max to zero (dealiazing)
      do n_r_out=rscheme_oc%n_max+1,n_r_max
         do lm1=lmStart,lmStop
            s(lm1,n_r_out)=zero
         end do
      end do

      !PERFON('upS_drv')
      all_lms=lmStop-lmStart+1
#ifdef WITHOMP
      if (all_lms < maxThreads) then
         call omp_set_num_threads(all_lms)
         per_thread=1
      else
         per_thread=all_lms/omp_get_max_threads()
      end if
#else
      per_thread=all_lms
#endif
      !$OMP PARALLEL &
      !$OMP private(iThread,start_lm,stop_lm) &
      !$OMP shared(per_thread,lmStart,lmStop,nThreads) &
      !$OMP shared(s,ds,dsdtLast,rscheme_oc) &
      !$OMP shared(n_r_max,work_LMloc,llm,ulm) &
      !$OMP shared(n_r_cmb,n_r_icb,dsdt,coex,opr,hdif_S) &
      !$OMP shared(st_map,lm2l,lm2m,kappa,beta,dLtemp0,or1) &
      !$OMP shared(dLkappa,dLh,or2)
      !$OMP DO
      do iThread=0,nThreads-1
         start_lm=lmStart+iThread*per_thread
         stop_lm = start_lm+per_thread-1
         if (iThread == nThreads-1) stop_lm=lmStop
         call get_ddr(s, ds, work_LMloc, ulm-llm+1, start_lm-llm+1, &
              &       stop_lm-llm+1, n_r_max, rscheme_oc, l_dct_in=.false.)
         call rscheme_oc%costf1(s,ulm-llm+1,start_lm-llm+1,stop_lm-llm+1)
      end do
      !$OMP end do

      !-- Calculate explicit time step part:
      !$OMP do private(nR,lm1)
      do nR=n_r_cmb+1,n_r_icb-1
         do lm1=lmStart,lmStop
            dsdtLast(lm1,nR)=                              dsdt(lm1,nR) &
            &      - coex*opr*hdif_S(st_map%lm2(lm2l(lm1),lm2m(lm1))) * &
            &        kappa(nR) *                   ( work_LMloc(lm1,nR) &
            &        + ( beta(nR)+dLtemp0(nR)+two*or1(nR)+dLkappa(nR) ) &
            &                                              * ds(lm1,nR) &
            &        - dLh(st_map%lm2(lm2l(lm1),lm2m(lm1)))*or2(nR)     &
            &                                              *  s(lm1,nR) )
         end do
      end do
      !$OMP end do
      !$OMP END PARALLEL
#ifdef WITHOMP
      call omp_set_num_threads(maxThreads)
#endif
      !PERFOFF

   end subroutine updateS
!------------------------------------------------------------------------------
   subroutine updateS_ala(s,ds,w,dVSrLM,dsdt,dsdtLast,w1,coex,dt,nLMB)
      !
      !  updates the entropy field s and its radial derivatives
      !  adds explicit part to time derivatives of s
      !

      !-- Input of variables:
      real(cp),    intent(in) :: w1        ! weight for time step !
      real(cp),    intent(in) :: coex      ! factor depending on alpha
      real(cp),    intent(in) :: dt        ! time step
      integer,     intent(in) :: nLMB
      complex(cp), intent(in) :: w(llm:ulm,n_r_max)
      complex(cp), intent(inout) :: dVSrLM(llm:ulm,n_r_max)

      !-- Input/output of scalar fields:
      complex(cp), intent(inout) :: s(llm:ulm,n_r_max)
      complex(cp), intent(inout) :: dsdt(llm:ulm,n_r_max)
      complex(cp), intent(inout) :: dsdtLast(llm:ulm,n_r_max)
      !-- Output: udpated s,ds,dsdtLast
      complex(cp), intent(out) :: ds(llm:ulm,n_r_max)

      !-- Local variables:
      real(cp) :: w2            ! weight of second time step
      real(cp) :: O_dt
      integer :: l1,m1              ! degree and order
      integer :: lm1,lmB,lm         ! position of (l,m) in array
      integer :: lmStart,lmStop
      integer :: nLMB2
      integer :: nR                 ! counts radial grid points
      integer :: n_r_out             ! counts cheb modes
      real(cp) ::  rhs(n_r_max) ! real RHS for l=m=0

      integer, pointer :: nLMBs2(:),lm2l(:),lm2m(:)
      integer, pointer :: sizeLMB2(:,:),lm2(:,:)
      integer, pointer :: lm22lm(:,:,:),lm22l(:,:,:),lm22m(:,:,:)

      integer :: threadid,nThreads,iThread,all_lms,per_thread,start_lm,stop_lm
      integer :: iChunk,nChunks,size_of_last_chunk,lmB0

      if ( .not. l_update_s ) return

      nLMBs2(1:nLMBs) => lo_sub_map%nLMBs2
      sizeLMB2(1:,1:) => lo_sub_map%sizeLMB2
      lm22lm(1:,1:,1:) => lo_sub_map%lm22lm
      lm22l(1:,1:,1:) => lo_sub_map%lm22l
      lm22m(1:,1:,1:) => lo_sub_map%lm22m
      lm2(0:,0:) => lo_map%lm2
      lm2l(1:lm_max) => lo_map%lm2l
      lm2m(1:lm_max) => lo_map%lm2m


      lmStart     =lmStartB(nLMB)
      lmStop      =lmStopB(nLMB)
      w2  =one-w1
      O_dt=one/dt


      !PERFON('upS_fin')
      !$OMP PARALLEL  &
      !$OMP private(iThread,start_lm,stop_lm,nR,lm)        &
      !$OMP shared(all_lms,per_thread)                     &
      !$OMP shared(dVSrLM,rscheme_oc,dsdt,orho1,dentropy0) &
      !$OMP shared(dLtemp0,w,or2,lmStart,lmStop,temp0,dLh) &
      !$OMP shared(n_r_max,work_LMloc,nThreads,llm,ulm)
      !$OMP SINGLE
#ifdef WITHOMP
      nThreads=omp_get_num_threads()
#else
      nThreads=1
#endif
      !-- Get radial derivatives of s: work_LMloc,dsdtLast used as work arrays
      all_lms=lmStop-lmStart+1
      per_thread=all_lms/nThreads
      !$OMP END SINGLE
      !$OMP BARRIER
      !$OMP DO
      do iThread=0,nThreads-1
         start_lm=lmStart+iThread*per_thread
         stop_lm = start_lm+per_thread-1
         if (iThread == nThreads-1) stop_lm=lmStop

         !--- Finish calculation of dsdt:
         call get_dr( dVSrLM,work_LMloc,ulm-llm+1,start_lm-llm+1,    &
              &       stop_lm-llm+1,n_r_max,rscheme_oc, nocopy=.true. )
      end do
      !$OMP end do

      !$OMP DO
      do nR=1,n_r_max
         do lm=lmStart,lmStop
            dsdt(lm,nR)=           orho1(nR)*        dsdt(lm,nR) - & 
            &        or2(nR)*orho1(nR)*        work_LMloc(lm,nR) + &
            &        or2(nR)*orho1(nR)*dLtemp0(nR)*dVSrLM(lm,nR) - &
            &        dLh(st_map%lm2(lm2l(lm),lm2m(lm)))*or2(nR)*   &
            &        orho1(nR)*temp0(nR)*dentropy0(nR)*w(lm,nR)

         end do
      end do
      !$OMP end do
      !$OMP END PARALLEL
      !PERFOFF

      ! one subblock is linked to one l value and needs therefore once the matrix
      !$OMP PARALLEL default(shared)
      !$OMP SINGLE
      do nLMB2=1,nLMBs2(nLMB)
         ! this inner loop is in principle over the m values which belong to the
         ! l value
         !$OMP TASK default(shared) &
         !$OMP firstprivate(nLMB2) &
         !$OMP private(lm,lm1,l1,m1,lmB,threadid) &
         !$OMP private(nChunks,size_of_last_chunk,iChunk)
         nChunks = (sizeLMB2(nLMB2,nLMB)+chunksize-1)/chunksize
         size_of_last_chunk = chunksize + (sizeLMB2(nLMB2,nLMB)-nChunks*chunksize)

         ! This task treats one l given by l1
         l1=lm22l(1,nLMB2,nLMB)
         !write(*,"(3(A,I3),A)") "Launching task for nLMB2=",nLMB2," (l=",l1,") and scheduling ",nChunks," subtasks."

         if ( l1 == 0 ) then
            if ( .not. lSmat(l1) ) then
#ifdef WITH_PRECOND_S0
               call get_s0Mat(dt,s0Mat,s0Pivot,s0Mat_fac)
#else
               call get_s0Mat(dt,s0Mat,s0Pivot)
#endif
               lSmat(l1)=.true.
            end if
         else
            if ( .not. lSmat(l1) ) then
#ifdef WITH_PRECOND_S
               call get_sMat(dt,l1,hdif_S(st_map%lm2(l1,0)), &
                    &        sMat(1,1,l1),sPivot(1,l1),sMat_fac(1,l1))
#else
               call get_sMat(dt,l1,hdif_S(st_map%lm2(l1,0)), &
                    &        sMat(1,1,l1),sPivot(1,l1))
#endif
               lSmat(l1)=.true.
             !write(*,"(A,I3,ES22.14)") "sMat: ",l1,SUM( sMat(:,:,l1) )
            end if
         end if

         do iChunk=1,nChunks
            !$OMP TASK default(shared) &
            !$OMP firstprivate(iChunk) &
            !$OMP private(lmB0,lmB,lm,lm1,m1,nR,n_r_out) &
            !$OMP private(threadid)
#ifdef WITHOMP
            threadid = omp_get_thread_num()
#else
            threadid = 0
#endif
            lmB0=(iChunk-1)*chunksize
            lmB=lmB0

            do lm=lmB0+1,min(iChunk*chunksize,sizeLMB2(nLMB2,nLMB))
               !do lm=1,sizeLMB2(nLMB2,nLMB)
               lm1=lm22lm(lm,nLMB2,nLMB)
               !l1 =lm22l(lm,nLMB2,nLMB)
               m1 =lm22m(lm,nLMB2,nLMB)

               if ( l1 == 0 ) then
                  rhs(1)=      real(tops(0,0))
                  rhs(n_r_max)=real(bots(0,0))
                  do nR=2,n_r_max-1
                     rhs(nR)=real(s(lm1,nR))*O_dt+w1*real(dsdt(lm1,nR)) + &
                     &       w2*real(dsdtLast(lm1,nR))
                  end do

#ifdef WITH_PRECOND_S0
                  rhs = s0Mat_fac*rhs
#endif

                  call solve_mat(s0Mat,n_r_max,n_r_max,s0Pivot,rhs)

               else ! l1  /=  0
                  lmB=lmB+1

                  rhs1(1,lmB,threadid)=      tops(l1,m1)
                  rhs1(n_r_max,lmB,threadid)=bots(l1,m1)
#ifdef WITH_PRECOND_S
                  rhs1(1,lmB,threadid)=      sMat_fac(1,l1)*rhs1(1,lmB,threadid)
                  rhs1(n_r_max,lmB,threadid)=sMat_fac(1,l1)*rhs1(n_r_max,lmB,threadid)
#endif
                  do nR=2,n_r_max-1
                     rhs1(nR,lmB,threadid)=s(lm1,nR)*O_dt + w1*dsdt(lm1,nR)  &
                     &                     + w2*dsdtLast(lm1,nR)
#ifdef WITH_PRECOND_S
                     rhs1(nR,lmB,threadid) = sMat_fac(nR,l1)*rhs1(nR,lmB,threadid)
#endif
                  end do
               end if
            end do
            !PERFOFF

            !PERFON('upS_sol')
            if ( lmB  >  lmB0 ) then
               call solve_mat(sMat(:,:,l1),n_r_max,n_r_max, &
                    &         sPivot(:,l1),rhs1(:,lmB0+1:lmB,threadid),lmB-lmB0)
            end if
            !PERFOFF

            lmB=lmB0
            !PERFON('upS_af')
            do lm=lmB0+1,min(iChunk*chunksize,sizeLMB2(nLMB2,nLMB))
             !do lm=1,sizeLMB2(nLMB2,nLMB)
               lm1=lm22lm(lm,nLMB2,nLMB)
               !l1 =lm22l(lm,nLMB2,nLMB)
               m1 =lm22m(lm,nLMB2,nLMB)
               if ( l1 == 0 ) then
                  do n_r_out=1,rscheme_oc%n_max
                     s(lm1,n_r_out)=rhs(n_r_out)
                  end do
               else
                  lmB=lmB+1
                  if ( m1 > 0 ) then
                     do n_r_out=1,rscheme_oc%n_max
                        s(lm1,n_r_out)=rhs1(n_r_out,lmB,threadid)
                     end do
                  else
                     do n_r_out=1,rscheme_oc%n_max
                        s(lm1,n_r_out)= cmplx(real(rhs1(n_r_out,lmB,threadid)), &
                        &                    0.0_cp,kind=cp)
                     end do
                  end if
               end if
            end do
            !PERFOFF
            !$OMP END TASK
         end do
         !$OMP END TASK
      end do     ! loop over lm blocks
      !$OMP END SINGLE
      !$OMP END PARALLEL

      !write(*,"(A,2ES22.12)") "s after = ",SUM(s)
      !-- set cheb modes > rscheme_oc%n_max to zero (dealiazing)
      do n_r_out=rscheme_oc%n_max+1,n_r_max
         do lm1=lmStart,lmStop
            s(lm1,n_r_out)=zero
         end do
      end do

      !PERFON('upS_drv')
      all_lms=lmStop-lmStart+1
#ifdef WITHOMP
      if (all_lms < maxThreads) then
         call omp_set_num_threads(all_lms)
         per_thread=1
      else
         per_thread=all_lms/omp_get_max_threads()
      end if
#else
      per_thread=all_lms
#endif
      !$OMP PARALLEL &
      !$OMP private(iThread,start_lm,stop_lm) &
      !$OMP shared(per_thread,nThreads) &
      !$OMP shared(s,ds,w,dsdtLast,rscheme_oc) &
      !$OMP shared(n_r_max,work_LMloc,llm,ulm) &
      !$OMP shared(n_r_cmb,n_r_icb,lmStart,lmStop,dsdt,coex,opr,hdif_S) &
      !$OMP shared(st_map,lm2l,lm2m,kappa,beta,or1,dLkappa,dLh,or2) &
      !$OMP shared(orho1)
      !$OMP DO
      do iThread=0,nThreads-1
         start_lm=lmStart+iThread*per_thread
         stop_lm = start_lm+per_thread-1
         if (iThread == nThreads-1) stop_lm=lmStop
         call get_ddr(s, ds, work_LMloc, ulm-llm+1, start_lm-llm+1, &
              &       stop_lm-llm+1, n_r_max, rscheme_oc, l_dct_in=.false.)
         call rscheme_oc%costf1(s,ulm-llm+1,start_lm-llm+1,stop_lm-llm+1)
      end do
      !$OMP end do

      !-- Calculate explicit time step part:
      !$OMP do private(nR,lm1)
      do nR=n_r_cmb+1,n_r_icb-1
         do lm1=lmStart,lmStop
           dsdtLast(lm1,nR)=dsdt(lm1,nR) &
           &      - coex*opr*hdif_S(st_map%lm2(lm2l(lm1),lm2m(lm1)))*kappa(nR) * &
           &        (                                         work_LMloc(lm1,nR) &
           &                + ( beta(nR)+two*or1(nR)+dLkappa(nR) ) *  ds(lm1,nR) &
           &          - dLh(st_map%lm2(lm2l(lm1),lm2m(lm1)))*or2(nR)*  s(lm1,nR) ) 
         end do
      end do
      !$OMP end do
      !$OMP END PARALLEL
#ifdef WITHOMP
      call omp_set_num_threads(maxThreads)
#endif
      !PERFOFF

   end subroutine updateS_ala
!-------------------------------------------------------------------------------
#ifdef WITH_PRECOND_S0
   subroutine get_s0Mat(dt,sMat,sPivot,sMat_fac)
#else
   subroutine get_s0Mat(dt,sMat,sPivot)
#endif
      !
      !  Purpose of this subroutine is to contruct the time step matrix   
      !  sMat0                                                            
      !

      !-- Input variables
      real(cp), intent(in) :: dt

      !-- Output variables
      real(cp), intent(out) :: sMat(n_r_max,n_r_max)
      integer,  intent(out) :: sPivot(n_r_max)
#ifdef WITH_PRECOND_S0
      real(cp), intent(out) :: sMat_fac(n_r_max)
#endif

      !-- Local variables:
      integer :: info,nR_out,nR
      real(cp) :: O_dt

      O_dt=one/dt
    
      !----- Boundary condition:
      do nR_out=1,rscheme_oc%n_max
    
         if ( ktops == 1 ) then
            !--------- Constant entropy at CMB:
            sMat(1,nR_out)=rscheme_oc%rnorm*rscheme_oc%rMat(1,nR_out)
         else
            !--------- Constant flux at CMB:
            sMat(1,nR_out)=rscheme_oc%rnorm*rscheme_oc%drMat(1,nR_out)
         end if
         if ( kbots == 1 ) then
            !--------- Constant entropy at ICB:
            sMat(n_r_max,nR_out)=rscheme_oc%rnorm* &
            &                    rscheme_oc%rMat(n_r_max,nR_out)
         else
            !--------- Constant flux at ICB:
            sMat(n_r_max,nR_out)=rscheme_oc%rnorm* &
            &                    rscheme_oc%drMat(n_r_max,nR_out)
         end if
      end do
      if ( rscheme_oc%n_max < n_r_max ) then ! fill with zeros !
         do nR_out=rscheme_oc%n_max+1,n_r_max
            sMat(1,nR_out)      =0.0_cp
            sMat(n_r_max,nR_out)=0.0_cp
         end do
      end if
    
      if ( l_anelastic_liquid ) then
         do nR_out=1,n_r_max
            do nR=2,n_r_max-1
               sMat(nR,nR_out)= rscheme_oc%rnorm * (                          &
               &                            O_dt*rscheme_oc%rMat(nR,nR_out) - & 
               &      alpha*opr*kappa(nR)*(    rscheme_oc%d2rMat(nR,nR_out) + &
               & (beta(nR)+two*or1(nR)+dLkappa(nR))*                          &
               &                                rscheme_oc%drMat(nR,nR_out) ) )
            end do
         end do
      else
         do nR_out=1,n_r_max
            do nR=2,n_r_max-1
               sMat(nR,nR_out)= rscheme_oc%rnorm * (                         &
               &                           O_dt*rscheme_oc%rMat(nR,nR_out) - & 
               &     alpha*opr*kappa(nR)*(    rscheme_oc%d2rMat(nR,nR_out) + &
               & (beta(nR)+dLtemp0(nR)+two*or1(nR)+dLkappa(nR))*             &
               &                               rscheme_oc%drMat(nR,nR_out) ) )
            end do
         end do
      end if
    
      !----- Factors for highest and lowest cheb mode:
      do nR=1,n_r_max
         sMat(nR,1)      =rscheme_oc%boundary_fac*sMat(nR,1)
         sMat(nR,n_r_max)=rscheme_oc%boundary_fac*sMat(nR,n_r_max)
      end do
    
#ifdef WITH_PRECOND_S0
      ! compute the linesum of each line
      do nR=1,n_r_max
         sMat_fac(nR)=one/maxval(abs(sMat(nR,:)))
      end do
      ! now divide each line by the linesum to regularize the matrix
      do nr=1,n_r_max
         sMat(nR,:) = sMat(nR,:)*sMat_fac(nR)
      end do
#endif
    
      !---- LU decomposition:
      call prepare_mat(sMat,n_r_max,n_r_max,sPivot,info)
      if ( info /= 0 ) then
         call abortRun('! Singular matrix sMat0!')
      end if

   end subroutine get_s0Mat
!-----------------------------------------------------------------------------
#ifdef WITH_PRECOND_S
   subroutine get_Smat(dt,l,hdif,sMat,sPivot,sMat_fac)
#else
   subroutine get_Smat(dt,l,hdif,sMat,sPivot)
#endif
      !
      !  Purpose of this subroutine is to contruct the time step matricies
      !  sMat(i,j) and s0mat for the entropy equation.                    
      !
      
      !-- Input variables
      real(cp), intent(in) :: dt
      real(cp), intent(in) :: hdif
      integer,  intent(in) :: l

      !-- Output variables
      real(cp), intent(out) :: sMat(n_r_max,n_r_max)
      integer,  intent(out) :: sPivot(n_r_max)
#ifdef WITH_PRECOND_S
      real(cp),intent(out) :: sMat_fac(n_r_max)
#endif

      !-- Local variables:
      integer :: info,nR_out,nR
      real(cp) :: O_dt,dLh

#ifdef MATRIX_CHECK
      integer :: i,j
      real(cp) :: rcond
      integer ::ipiv(n_r_max),iwork(n_r_max)
      real(cp) :: work(4*n_r_max),anorm,linesum
      real(cp) :: temp_Mat(n_r_max,n_r_max)
      integer,save :: counter=0
      integer :: filehandle
      character(len=100) :: filename
#endif

      O_dt=one/dt

      dLh=real(l*(l+1),kind=cp)

      !----- Boundary coditions:
      do nR_out=1,rscheme_oc%n_max
         if ( ktops == 1 ) then
            sMat(1,nR_out)=rscheme_oc%rnorm*rscheme_oc%rMat(1,nR_out)
         else
            sMat(1,nR_out)=rscheme_oc%rnorm*rscheme_oc%drMat(1,nR_out)
         end if
         if ( kbots == 1 ) then
            sMat(n_r_max,nR_out)=rscheme_oc%rnorm* &
            &                    rscheme_oc%rMat(n_r_max,nR_out)
         else
            sMat(n_r_max,nR_out)=rscheme_oc%rnorm* &
            &                    rscheme_oc%drMat(n_r_max,nR_out)
         end if
      end do
      if ( rscheme_oc%n_max < n_r_max ) then ! fill with zeros !
         do nR_out=rscheme_oc%n_max+1,n_r_max
            sMat(1,nR_out)      =0.0_cp
            sMat(n_r_max,nR_out)=0.0_cp
         end do
      end if

      !----- Other points:
      if ( l_anelastic_liquid ) then
         do nR_out=1,n_r_max
            do nR=2,n_r_max-1
               sMat(nR,nR_out)= rscheme_oc%rnorm * (                        &
               &                          O_dt*rscheme_oc%rMat(nR,nR_out) - &
               & alpha*opr*hdif*kappa(nR)*(  rscheme_oc%d2rMat(nR,nR_out) + &
               &( beta(nR)+two*or1(nR)+dLkappa(nR) )*                       &
               &                              rscheme_oc%drMat(nR,nR_out) - &
               &      dLh*or2(nR)*             rscheme_oc%rMat(nR,nR_out) ) )
            end do
         end do
      else
         do nR_out=1,n_r_max
            do nR=2,n_r_max-1
               sMat(nR,nR_out)= rscheme_oc%rnorm * (                        &
               &                          O_dt*rscheme_oc%rMat(nR,nR_out) - &
               & alpha*opr*hdif*kappa(nR)*(  rscheme_oc%d2rMat(nR,nR_out) + &
               & ( beta(nR)+dLtemp0(nR)+                                    &
               &   two*or1(nR)+dLkappa(nR) )* rscheme_oc%drMat(nR,nR_out) - &
               &      dLh*or2(nR)*             rscheme_oc%rMat(nR,nR_out) ) )
            end do
         end do
      end if

      !----- Factor for highest and lowest cheb:
      do nR=1,n_r_max
         sMat(nR,1)      =rscheme_oc%boundary_fac*sMat(nR,1)
         sMat(nR,n_r_max)=rscheme_oc%boundary_fac*sMat(nR,n_r_max)
      end do

#ifdef WITH_PRECOND_S
      ! compute the linesum of each line
      do nR=1,n_r_max
         sMat_fac(nR)=one/maxval(abs(sMat(nR,:)))
      end do
      ! now divide each line by the linesum to regularize the matrix
      do nr=1,n_r_max
         sMat(nR,:) = sMat(nR,:)*sMat_fac(nR)
      end do
#endif

#ifdef MATRIX_CHECK
      ! copy the sMat to a temporary variable for modification
      write(filename,"(A,I3.3,A,I3.3,A)") "sMat_",l,"_",counter,".dat"
      open(newunit=filehandle,file=trim(filename))
      counter= counter+1

      do i=1,n_r_max
         do j=1,n_r_max
            write(filehandle,"(2ES20.12,1X)",advance="no") sMat(i,j)
         end do
         write(filehandle,"(A)") ""
      end do
      close(filehandle)
      temp_Mat=sMat
      anorm = 0.0_cp
      do i=1,n_r_max
         linesum = 0.0_cp
         do j=1,n_r_max
            linesum = linesum + abs(temp_Mat(i,j))
         end do
         if (linesum  >  anorm) anorm=linesum
      end do
      !write(*,"(A,ES20.12)") "anorm = ",anorm
      ! LU factorization
      call dgetrf(n_r_max,n_r_max,temp_Mat,n_r_max,ipiv,info)
      ! estimate the condition number
      call dgecon('I',n_r_max,temp_Mat,n_r_max,anorm,rcond,work,iwork,info)
      write(*,"(A,I3,A,ES11.3)") "inverse condition number of sMat for l=",l," is ",rcond
#endif

!----- LU decomposition:
      call prepare_mat(sMat,n_r_max,n_r_max,sPivot,info)
      if ( info /= 0 ) then
         call abortRun('Singular matrix sMat!')
      end if
            
   end subroutine get_Smat
!-----------------------------------------------------------------------------
end module updateS_mod
