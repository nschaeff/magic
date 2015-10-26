#include "mkl_trig_transforms.f90"

module cosine_transform_odd
 
   use constants, only: half
   use precision_mod, only: cp
   use mkl_trig_transforms

   implicit none

   private

   type, public :: costf_odd_t
      integer :: nRad
      real(cp) :: preFactor
      integer, allocatable  :: i_costf_init(:)
      real(cp), allocatable :: d_costf_init(:)
      type(DFTI_DESCRIPTOR), pointer :: r2c_handle
   contains
      procedure, private :: costf1_complex
      procedure, private :: costf1_complex_1d
      procedure, private :: costf1_real
      procedure, private :: costf1_real_1d
      procedure :: initialize
      procedure :: finalize
      generic :: costf1 => costf1_complex,costf1_real,costf1_real_1d,costf1_complex_1d
   end type costf_odd_t

contains

   subroutine initialize(this,n, ni, nd)

      class(costf_odd_t) :: this

      !-- Input variables
      integer, intent(in) :: n
      integer, intent(in) :: ni
      integer, intent(in) :: nd

      !-- Local variables:
      integer :: stat
      real(cp) :: fac
      real(cp) :: work(n)

      this%nRad = n
      allocate(this%i_costf_init(128))
      allocate(this%d_costf_init(1:5*(n-1)/2+2))

      call d_init_trig_transform(n-1,MKL_COSINE_TRANSFORM,this%i_costf_init, &
                                 this%d_costf_init,stat)
      call d_commit_trig_transform(work,this%r2c_handle,this%i_costf_init, &
                                   this%d_costf_init,stat)

      this%preFactor = sqrt(half*real(this%nRad-1,cp))
      !stat = DftiSetValue(this%r2c_handle, DFTI_FORWARD_SCALE, fac)
      !stat = DftiCommitDescriptor(this%r2c_handle)

   end subroutine initialize
!------------------------------------------------------------------------------
   subroutine finalize(this)

      class(costf_odd_t) :: this

      integer :: stat

      call free_trig_transform(this%r2c_handle, this%i_costf_init, stat)

      deallocate( this%d_costf_init )

   end subroutine finalize
!------------------------------------------------------------------------------
   subroutine costf1_real_1d(this,f,work)

      class(costf_odd_t) :: this

      !-- Input variables:
      real(cp), intent(inout) :: f(this%nRad)
      real(cp) :: work(*)

      !-- Local variables:
      integer :: stat

      call d_forward_trig_transform(f,this%r2c_handle,this%i_costf_init, &
                                    this%d_costf_init,stat)
      f(:) = this%preFactor*f(:)

   end subroutine costf1_real_1d
!------------------------------------------------------------------------------
   subroutine costf1_complex_1d(this,f,work)

      class(costf_odd_t) :: this

      !-- Input variables:
      complex(cp), intent(inout) :: f(this%nRad)
      complex(cp) :: work(this%nRad)

      !-- Local variables:
      real(cp) :: work_real(this%nRad), work_imag(this%nRad)
      integer :: stat

      work_real(:) = real(f(:))
      work_imag(:) = aimag(f(:))

      call d_forward_trig_transform(work_real,this%r2c_handle,this%i_costf_init, &
                                    this%d_costf_init,stat)

      call d_forward_trig_transform(work_imag,this%r2c_handle,this%i_costf_init, &
                                    this%d_costf_init,stat)

      f(:) = this%preFactor*cmplx(work_real, work_imag, kind=cp)

   end subroutine costf1_complex_1d
!------------------------------------------------------------------------------
   subroutine costf1_complex(this,f,n_f_max,n_f_start,n_f_stop,work)

      class(costf_odd_t) :: this

      !-- Input variables:
      integer,  intent(in) :: n_f_start,n_f_stop ! columns to be transformed
      integer,  intent(in) :: n_f_max
      complex(cp), intent(inout) :: f(n_f_max,this%nRad)
      complex(cp) :: work(n_f_max,this%nRad)

      !-- Local variables:
      real(cp) :: work_real(this%nRad)
      real(cp) :: work_imag(this%nRad)
      integer :: stat, n_f

      do n_f=n_f_start,n_f_stop
         work_real(:) = real(f(n_f,:))
         work_imag(:) = aimag(f(n_f,:))
         call d_forward_trig_transform(work_real,this%r2c_handle,this%i_costf_init, &
                                       this%d_costf_init,stat)

         call d_forward_trig_transform(work_imag,this%r2c_handle,this%i_costf_init, &
                                       this%d_costf_init,stat)
         f(n_f,:)=this%preFactor*cmplx(work_real,work_imag,kind=cp)
      end do

   end subroutine costf1_complex
!------------------------------------------------------------------------------
   subroutine costf1_real(this,f,n_f_max,n_f_start,n_f_stop,work)

      class(costf_odd_t) :: this

      !-- Input variables:
      integer,  intent(in) :: n_f_start,n_f_stop ! columns to be transformed
      integer,  intent(in) :: n_f_max
      real(cp), intent(inout) :: f(n_f_max,this%nRad)
      real(cp) :: work(n_f_max,*)

      !-- Local variables:
      integer :: stat, n_f

      do n_f=n_f_start,n_f_stop
         call d_forward_trig_transform(f(n_f,:),this%r2c_handle,this%i_costf_init, &
                                      this%d_costf_init,stat)
         f(n_f,:)=this%preFactor*f(n_f,:)
      end do

   end subroutine costf1_real
!------------------------------------------------------------------------------
end module cosine_transform_odd
