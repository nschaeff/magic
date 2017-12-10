!-------------------------------------------------------------------------------
module arrays_dist
!@>details Distributed implementation of grid_space_arrays_t and nonlinear_lm_t
!> for coord_theta>1.
!
!@>author Rafael Lago, MPCDF, December 2017
!-------------------------------------------------------------------------------
   use grid_space_arrays_mod, only: grid_space_arrays_t
   use mem_alloc, only: bytes_allocated
   use truncation
   use nonlinear_lm_mod, only: nonlinear_lm_t
   use logic, only: l_conv_nl, l_heat_nl, l_mag_nl, l_anel, l_mag_LF, &
       &            l_RMS, l_chemical_conv, l_TP_form, l_HT, l_heat, l_conv, &
       &            l_mag_kin
       
   use fields, only: s_Rloc,ds_Rloc, z_Rloc,dz_Rloc, p_Rloc,dp_Rloc, &
       &             b_Rloc,db_Rloc,ddb_Rloc, aj_Rloc,dj_Rloc,       &
       &             w_Rloc,dw_Rloc,ddw_Rloc, xi_Rloc,               &
       &             s_dist,ds_dist, z_dist,dz_dist, p_dist,dp_dist, &
       &             b_dist,db_dist,ddb_dist, aj_dist,dj_dist,       &
       &             w_dist,dw_dist,ddw_dist, xi_dist
   use precision_mod
   use radial_data, only: nRstart, nRstop

   implicit none

   private
   
   type, public, extends(grid_space_arrays_t) :: grid_space_arrays_dist_t
   contains
      procedure :: initialize => initialize_grid_space             ! override
      procedure :: slice_all  => slice_all_grid_space
      procedure :: gather_all  => gather_all_grid_space
   end type grid_space_arrays_dist_t
   
   type, public, extends(nonlinear_lm_t) :: nonlinear_lm_dist_t
   contains
      procedure :: initialize => initialize_nonlinear_lm           ! override
      procedure :: slice_all  => slice_all_nonlinear_lm
      procedure :: gather_all  => gather_all_nonlinear_lm
   end type nonlinear_lm_dist_t

contains

!-------------------------------------------------------------------------------
   subroutine gather_all_grid_space(this, gsa_glb)
!@>author Rafael Lago, MPCDF, December 2017
!-------------------------------------------------------------------------------
      class(grid_space_arrays_dist_t) :: this
      class(grid_space_arrays_t)      :: gsa_glb
      
      call gather_f(this%Advr, gsa_glb%Advr)
      call gather_f(this%Advt, gsa_glb%Advt)
      call gather_f(this%Advp, gsa_glb%Advp)
      call gather_f(this%LFr, gsa_glb%LFr)
      call gather_f(this%LFt, gsa_glb%LFt)
      call gather_f(this%LFp, gsa_glb%LFp)
      call gather_f(this%VxBr, gsa_glb%VxBr)
      call gather_f(this%VxBt, gsa_glb%VxBt)
      call gather_f(this%VxBp, gsa_glb%VxBp)
      call gather_f(this%VSr, gsa_glb%VSr)
      call gather_f(this%VSt, gsa_glb%VSt)
      call gather_f(this%VSp, gsa_glb%VSp)
      call gather_f(this%ViscHeat, gsa_glb%ViscHeat)
      call gather_f(this%OhmLoss, gsa_glb%OhmLoss)

      if ( l_TP_form ) then
         call gather_f(this%VPr, gsa_glb%VPr)
      end if

      if ( l_chemical_conv ) then
         call gather_f(this%VXir, gsa_glb%VXir)
         call gather_f(this%VXit, gsa_glb%VXit)
         call gather_f(this%VXip, gsa_glb%VXip)
      end if

      !----- Fields calculated from these help arrays by legtf:
      call gather_f(this%vrc, gsa_glb%vrc)
      call gather_f(this%vtc, gsa_glb%vtc)
      call gather_f(this%vpc, gsa_glb%vpc)
      call gather_f(this%dvrdrc, gsa_glb%dvrdrc)
      call gather_f(this%dvtdrc, gsa_glb%dvtdrc)
      call gather_f(this%dvpdrc, gsa_glb%dvpdrc)
      call gather_f(this%cvrc, gsa_glb%cvrc)
      call gather_f(this%dvrdtc, gsa_glb%dvrdtc)
      call gather_f(this%dvrdpc, gsa_glb%dvrdpc)
      call gather_f(this%dvtdpc, gsa_glb%dvtdpc)
      call gather_f(this%dvpdpc, gsa_glb%dvpdpc)
      call gather_f(this%brc, gsa_glb%brc)
      call gather_f(this%btc, gsa_glb%btc)
      call gather_f(this%bpc, gsa_glb%bpc)
      call gather_f(this%cbrc, gsa_glb%cbrc)
      call gather_f(this%cbtc, gsa_glb%cbtc)
      call gather_f(this%cbpc, gsa_glb%cbpc)
      call gather_f(this%sc, gsa_glb%sc)
      call gather_f(this%drSc, gsa_glb%drSc)
      call gather_f(this%pc, gsa_glb%pc)
      call gather_f(this%dsdtc, gsa_glb%dsdtc)
      call gather_f(this%dsdpc, gsa_glb%dsdpc)

      if ( l_chemical_conv ) then
         call gather_f(this%xic, gsa_glb%xic)
      else
         gsa_glb%xic(1,1) = this%xic(1,1)
      end if

      !-- RMS Calculations
      if ( l_RMS ) then
         call gather_f(this%Advt2, gsa_glb%Advt2)
         call gather_f(this%Advp2, gsa_glb%Advp2)
         call gather_f(this%LFt2, gsa_glb%LFt2)
         call gather_f(this%LFp2, gsa_glb%LFp2)
         call gather_f(this%CFt2, gsa_glb%CFt2)
         call gather_f(this%CFp2, gsa_glb%CFp2)
         call gather_f(this%p1, gsa_glb%p1)
         call gather_f(this%p2, gsa_glb%p2)
      end if
      
   end subroutine gather_all_grid_space

!-------------------------------------------------------------------------------
   subroutine gather_all_nonlinear_lm(this, nl_lm_glb)
!@>author Rafael Lago, MPCDF, December 2017
!-------------------------------------------------------------------------------
      class(nonlinear_lm_dist_t) :: this
      class(nonlinear_lm_t)      :: nl_lm_glb
      
      call gather_FlmP(this%AdvrLM, nl_lm_glb%AdvrLM)   
      call gather_FlmP(this%AdvtLM, nl_lm_glb%AdvtLM)   
      call gather_FlmP(this%AdvpLM, nl_lm_glb%AdvpLM)   
      call gather_FlmP(this%LFrLM, nl_lm_glb%LFrLM)    
      call gather_FlmP(this%LFtLM, nl_lm_glb%LFtLM)    
      call gather_FlmP(this%LFpLM, nl_lm_glb%LFpLM)    
      call gather_FlmP(this%VxBrLM, nl_lm_glb%VxBrLM)   
      call gather_FlmP(this%VxBtLM, nl_lm_glb%VxBtLM)   
      call gather_FlmP(this%VxBpLM, nl_lm_glb%VxBpLM)   
      call gather_FlmP(this%VSrLM, nl_lm_glb%VSrLM)    
      call gather_FlmP(this%VStLM, nl_lm_glb%VStLM)    
      call gather_FlmP(this%VSpLM, nl_lm_glb%VSpLM)    
      call gather_FlmP(this%ViscHeatLM, nl_lm_glb%ViscHeatLM)
      call gather_FlmP(this%OhmLossLM, nl_lm_glb%OhmLossLM)

      if ( l_TP_form ) then
         call gather_FlmP(this%VPrLM, nl_lm_glb%VPrLM)    
      end if
      
      if ( l_chemical_conv ) then
         call gather_FlmP(this%VXirLM, nl_lm_glb%VXirLM)    
         call gather_FlmP(this%VXitLM, nl_lm_glb%VXitLM)    
         call gather_FlmP(this%VXipLM, nl_lm_glb%VXipLM)    
      end if

      !-- RMS calculations
      if ( l_RMS ) then
         call gather_FlmP(this%Advt2LM, nl_lm_glb%Advt2LM)
         call gather_FlmP(this%Advp2LM, nl_lm_glb%Advp2LM)
         call gather_FlmP(this%LFt2LM, nl_lm_glb%LFt2LM)
         call gather_FlmP(this%LFp2LM, nl_lm_glb%LFp2LM)
         call gather_FlmP(this%CFt2LM, nl_lm_glb%CFt2LM)
         call gather_FlmP(this%CFp2LM, nl_lm_glb%CFp2LM)
         call gather_FlmP(this%p1LM, nl_lm_glb%p1LM)
         call gather_FlmP(this%p2LM, nl_lm_glb%p2LM)
      end if
      
   end subroutine gather_all_nonlinear_lm

!-------------------------------------------------------------------------------
   subroutine slice_all_grid_space(this, gsa_glb)
!@>author Rafael Lago, MPCDF, December 2017
!-------------------------------------------------------------------------------
      class(grid_space_arrays_dist_t) :: this
      class(grid_space_arrays_t)      :: gsa_glb
      
      call slice_f(gsa_glb%Advr, this%Advr)
      call slice_f(gsa_glb%Advt, this%Advt)
      call slice_f(gsa_glb%Advp, this%Advp)
      call slice_f(gsa_glb%LFr, this%LFr)
      call slice_f(gsa_glb%LFt, this%LFt)
      call slice_f(gsa_glb%LFp, this%LFp)
      call slice_f(gsa_glb%VxBr, this%VxBr)
      call slice_f(gsa_glb%VxBt, this%VxBt)
      call slice_f(gsa_glb%VxBp, this%VxBp)
      call slice_f(gsa_glb%VSr, this%VSr)
      call slice_f(gsa_glb%VSt, this%VSt)
      call slice_f(gsa_glb%VSp, this%VSp)
      call slice_f(gsa_glb%ViscHeat, this%ViscHeat)
      call slice_f(gsa_glb%OhmLoss, this%OhmLoss)

      if ( l_TP_form ) then
         call slice_f(gsa_glb%VPr, this%VPr)
      end if

      if ( l_chemical_conv ) then
         call slice_f(gsa_glb%VXir, this%VXir)
         call slice_f(gsa_glb%VXit, this%VXit)
         call slice_f(gsa_glb%VXip, this%VXip)
      end if

      !----- Fields calculated from these help arrays by legtf:
      call slice_f(gsa_glb%vrc, this%vrc )
      call slice_f(gsa_glb%vtc, this%vtc)
      call slice_f(gsa_glb%vpc, this%vpc)
      call slice_f(gsa_glb%dvrdrc, this%dvrdrc)
      call slice_f(gsa_glb%dvtdrc, this%dvtdrc)
      call slice_f(gsa_glb%dvpdrc, this%dvpdrc)
      call slice_f(gsa_glb%cvrc, this%cvrc)
      call slice_f(gsa_glb%dvrdtc, this%dvrdtc)
      call slice_f(gsa_glb%dvrdpc, this%dvrdpc)
      call slice_f(gsa_glb%dvtdpc, this%dvtdpc)
      call slice_f(gsa_glb%dvpdpc, this%dvpdpc)
      call slice_f(gsa_glb%brc, this%brc)
      call slice_f(gsa_glb%btc, this%btc)
      call slice_f(gsa_glb%bpc, this%bpc)
      call slice_f(gsa_glb%cbrc, this%cbrc)
      call slice_f(gsa_glb%cbtc, this%cbtc)
      call slice_f(gsa_glb%cbpc, this%cbpc)
      call slice_f(gsa_glb%sc, this%sc)
      call slice_f(gsa_glb%drSc, this%drSc)
      call slice_f(gsa_glb%pc, this%pc)
      call slice_f(gsa_glb%dsdtc, this%dsdtc)
      call slice_f(gsa_glb%dsdpc, this%dsdpc)

      if ( l_chemical_conv ) then
         call slice_f(gsa_glb%xic, this%xic)
      else
         this%xic(1,1) = gsa_glb%xic(1,1)
      end if

      !-- RMS Calculations
      if ( l_RMS ) then
         call slice_f(gsa_glb%Advt2, this%Advt2)
         call slice_f(gsa_glb%Advp2, this%Advp2)
         call slice_f(gsa_glb%LFt2, this%LFt2)
         call slice_f(gsa_glb%LFp2, this%LFp2)
         call slice_f(gsa_glb%CFt2, this%CFt2)
         call slice_f(gsa_glb%CFp2, this%CFp2)
         call slice_f(gsa_glb%p1, this%p1)
         call slice_f(gsa_glb%p2, this%p2)
      end if
      
   end subroutine slice_all_grid_space

!-------------------------------------------------------------------------------
   subroutine slice_all_nonlinear_lm(this, nl_lm_glb)
!@>author Rafael Lago, MPCDF, December 2017
!-------------------------------------------------------------------------------
      class(nonlinear_lm_dist_t) :: this
      class(nonlinear_lm_t)      :: nl_lm_glb
      
      call slice_FlmP(nl_lm_glb%AdvrLM, this%AdvrLM)   
      call slice_FlmP(nl_lm_glb%AdvtLM, this%AdvtLM)   
      call slice_FlmP(nl_lm_glb%AdvpLM, this%AdvpLM)   
      call slice_FlmP(nl_lm_glb%LFrLM, this%LFrLM)    
      call slice_FlmP(nl_lm_glb%LFtLM, this%LFtLM)    
      call slice_FlmP(nl_lm_glb%LFpLM, this%LFpLM)    
      call slice_FlmP(nl_lm_glb%VxBrLM, this%VxBrLM)   
      call slice_FlmP(nl_lm_glb%VxBtLM, this%VxBtLM)   
      call slice_FlmP(nl_lm_glb%VxBpLM, this%VxBpLM)   
      call slice_FlmP(nl_lm_glb%VSrLM, this%VSrLM)    
      call slice_FlmP(nl_lm_glb%VStLM, this%VStLM)    
      call slice_FlmP(nl_lm_glb%VSpLM, this%VSpLM)    
      call slice_FlmP(nl_lm_glb%ViscHeatLM, this%ViscHeatLM)
      call slice_FlmP(nl_lm_glb%OhmLossLM, this%OhmLossLM)

      if ( l_TP_form ) then
         call slice_FlmP(nl_lm_glb%VPrLM, this%VPrLM)    
      end if
      
      if ( l_chemical_conv ) then
         call slice_FlmP(nl_lm_glb%VXirLM, this%VXirLM)    
         call slice_FlmP(nl_lm_glb%VXitLM, this%VXitLM)    
         call slice_FlmP(nl_lm_glb%VXipLM, this%VXipLM)    
      end if

      !-- RMS calculations
      if ( l_RMS ) then
         call slice_FlmP(nl_lm_glb%Advt2LM, this%Advt2LM)
         call slice_FlmP(nl_lm_glb%Advp2LM, this%Advp2LM)
         call slice_FlmP(nl_lm_glb%LFt2LM, this%LFt2LM)
         call slice_FlmP(nl_lm_glb%LFp2LM, this%LFp2LM)
         call slice_FlmP(nl_lm_glb%CFt2LM, this%CFt2LM)
         call slice_FlmP(nl_lm_glb%CFp2LM, this%CFp2LM)
         call slice_FlmP(nl_lm_glb%p1LM, this%p1LM)
         call slice_FlmP(nl_lm_glb%p2LM, this%p2LM)
      end if
      
   end subroutine slice_all_nonlinear_lm

!-------------------------------------------------------------------------------   
   subroutine initialize_grid_space(this)
!@>author Rafael Lago, MPCDF, December 2017
!-------------------------------------------------------------------------------

      class(grid_space_arrays_dist_t) :: this

      allocate( this%Advr(n_phi_max,n_theta_beg:n_theta_end) )
      allocate( this%Advt(n_phi_max,n_theta_beg:n_theta_end) )
      allocate( this%Advp(n_phi_max,n_theta_beg:n_theta_end) )
      allocate( this%LFr(n_phi_max,n_theta_beg:n_theta_end) )
      allocate( this%LFt(n_phi_max,n_theta_beg:n_theta_end) )
      allocate( this%LFp(n_phi_max,n_theta_beg:n_theta_end) )
      allocate( this%VxBr(n_phi_max,n_theta_beg:n_theta_end) )
      allocate( this%VxBt(n_phi_max,n_theta_beg:n_theta_end) )
      allocate( this%VxBp(n_phi_max,n_theta_beg:n_theta_end) )
      allocate( this%VSr(n_phi_max,n_theta_beg:n_theta_end) )
      allocate( this%VSt(n_phi_max,n_theta_beg:n_theta_end) )
      allocate( this%VSp(n_phi_max,n_theta_beg:n_theta_end) )
      allocate( this%ViscHeat(n_phi_max,n_theta_beg:n_theta_end) )
      allocate( this%OhmLoss(n_phi_max,n_theta_beg:n_theta_end) )
      bytes_allocated=bytes_allocated + 14*n_phi_max*n_theta_loc*SIZEOF_DEF_REAL

      if ( l_TP_form ) then
         allocate( this%VPr(n_phi_max,n_theta_beg:n_theta_end) )
         bytes_allocated=bytes_allocated + n_phi_max*n_theta_loc*SIZEOF_DEF_REAL
      end if

      if ( l_chemical_conv ) then
         allocate( this%VXir(n_phi_max,n_theta_beg:n_theta_end) )
         allocate( this%VXit(n_phi_max,n_theta_beg:n_theta_end) )
         allocate( this%VXip(n_phi_max,n_theta_beg:n_theta_end) )
         bytes_allocated=bytes_allocated + 3*n_phi_max*n_theta_loc*SIZEOF_DEF_REAL
      end if

      !----- Fields calculated from these help arrays by legtf:
      allocate( this%vrc(n_phi_max,n_theta_beg:n_theta_end) ) 
      allocate( this%vtc(n_phi_max,n_theta_beg:n_theta_end) )
      allocate( this%vpc(n_phi_max,n_theta_beg:n_theta_end) )
      allocate( this%dvrdrc(n_phi_max,n_theta_beg:n_theta_end) )
      allocate( this%dvtdrc(n_phi_max,n_theta_beg:n_theta_end) )
      allocate( this%dvpdrc(n_phi_max,n_theta_beg:n_theta_end) )
      allocate( this%cvrc(n_phi_max,n_theta_beg:n_theta_end) )
      allocate( this%dvrdtc(n_phi_max,n_theta_beg:n_theta_end) )
      allocate( this%dvrdpc(n_phi_max,n_theta_beg:n_theta_end) )
      allocate( this%dvtdpc(n_phi_max,n_theta_beg:n_theta_end) )
      allocate( this%dvpdpc(n_phi_max,n_theta_beg:n_theta_end) )
      allocate( this%brc(n_phi_max,n_theta_beg:n_theta_end) )
      allocate( this%btc(n_phi_max,n_theta_beg:n_theta_end) )
      allocate( this%bpc(n_phi_max,n_theta_beg:n_theta_end) )
      this%btc=1.0e50_cp
      this%bpc=1.0e50_cp
      allocate( this%cbrc(n_phi_max,n_theta_beg:n_theta_end) )
      allocate( this%cbtc(n_phi_max,n_theta_beg:n_theta_end) )
      allocate( this%cbpc(n_phi_max,n_theta_beg:n_theta_end) )
      allocate( this%sc(n_phi_max,n_theta_beg:n_theta_end) )
      allocate( this%drSc(n_phi_max,n_theta_beg:n_theta_end) )
      allocate( this%pc(n_phi_max,n_theta_beg:n_theta_end) )
      allocate( this%dsdtc(n_phi_max,n_theta_beg:n_theta_end) )
      allocate( this%dsdpc(n_phi_max,n_theta_beg:n_theta_end) )
      bytes_allocated=bytes_allocated + 22*n_phi_max*n_theta_loc*SIZEOF_DEF_REAL

      if ( l_chemical_conv ) then
         allocate( this%xic(n_phi_max,n_theta_beg:n_theta_end) )
         bytes_allocated=bytes_allocated + n_phi_max*n_theta_loc*SIZEOF_DEF_REAL
      else
         allocate( this%xic(1,1) )
      end if

      !-- RMS Calculations
      if ( l_RMS ) then
         allocate ( this%Advt2(n_phi_max,n_theta_beg:n_theta_end) )
         allocate ( this%Advp2(n_phi_max,n_theta_beg:n_theta_end) )
         allocate ( this%LFt2(n_phi_max,n_theta_beg:n_theta_end) )
         allocate ( this%LFp2(n_phi_max,n_theta_beg:n_theta_end) )
         allocate ( this%CFt2(n_phi_max,n_theta_beg:n_theta_end) )
         allocate ( this%CFp2(n_phi_max,n_theta_beg:n_theta_end) )
         allocate ( this%p1(n_phi_max,n_theta_beg:n_theta_end) )
         allocate ( this%p2(n_phi_max,n_theta_beg:n_theta_end) )
         bytes_allocated=bytes_allocated + 8*n_phi_max*n_theta_loc*SIZEOF_DEF_REAL
      end if
      !write(*,"(A,I15,A)") "grid_space_arrays: allocated ",bytes_allocated,"B."

   end subroutine initialize_grid_space
   
!-------------------------------------------------------------------------------
   subroutine initialize_nonlinear_lm(this)
!@>author Rafael Lago, MPCDF, December 2017
!-------------------------------------------------------------------------------

      class(nonlinear_lm_dist_t) :: this

      allocate( this%AdvrLM(lmP_loc) )   
      allocate( this%AdvtLM(lmP_loc) )   
      allocate( this%AdvpLM(lmP_loc) )   
      allocate( this%LFrLM(lmP_loc) )    
      allocate( this%LFtLM(lmP_loc) )    
      allocate( this%LFpLM(lmP_loc) )    
      allocate( this%VxBrLM(lmP_loc) )   
      allocate( this%VxBtLM(lmP_loc) )   
      allocate( this%VxBpLM(lmP_loc) )   
      allocate( this%VSrLM(lmP_loc) )    
      allocate( this%VStLM(lmP_loc) )    
      allocate( this%VSpLM(lmP_loc) )    
      allocate( this%ViscHeatLM(lmP_loc) )
      allocate( this%OhmLossLM(lmP_loc) )
      bytes_allocated = bytes_allocated + 14*lmP_loc*SIZEOF_DEF_COMPLEX

      if ( l_TP_form ) then
         allocate( this%VPrLM(lmP_loc) )    
         bytes_allocated = bytes_allocated + lmP_loc*SIZEOF_DEF_COMPLEX
      end if
      
      if ( l_chemical_conv ) then
         allocate( this%VXirLM(lmP_loc) )    
         allocate( this%VXitLM(lmP_loc) )    
         allocate( this%VXipLM(lmP_loc) )    
         bytes_allocated = bytes_allocated + 3*lmP_loc*SIZEOF_DEF_COMPLEX
      end if

      !-- RMS calculations
      if ( l_RMS ) then
         allocate( this%Advt2LM(lmP_loc) )
         allocate( this%Advp2LM(lmP_loc) )
         allocate( this%LFt2LM(lmP_loc) )
         allocate( this%LFp2LM(lmP_loc) )
         allocate( this%CFt2LM(lmP_loc) )
         allocate( this%CFp2LM(lmP_loc) )
         allocate( this%p1LM(lmP_loc) )
         allocate( this%p2LM(lmP_loc) )
         bytes_allocated = bytes_allocated + 8*lmP_loc*SIZEOF_DEF_COMPLEX
      end if

   end subroutine initialize_nonlinear_lm


end module arrays_dist
