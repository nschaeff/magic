module fields
   !
   !  This module contains the potential fields and their radial
   !  derivatives
   !
   use precision_mod
   use mem_alloc, only: bytes_allocated
   use truncation, only: lm_max, n_r_max, lm_maxMag, n_r_maxMag, &
       &                 n_r_ic_maxMag, lm_loc, lm_locMag
   use logic, only: l_chemical_conv
   use LMLoop_data, only: llm, ulm, llmMag, ulmMag
   use radial_data, only: nRstart, nRstop
   use parallel_mod, only: coord_r, n_procs_theta
 
   implicit none

   private
 
   !-- Velocity potentials:
   complex(cp), public, allocatable, target :: flow_LMloc_container(:,:,:)
   complex(cp), public, allocatable, target :: flow_dist_container(:,:,:)
   complex(cp), public, pointer :: w_LMloc(:,:),dw_LMloc(:,:),ddw_LMloc(:,:)
   complex(cp), public, pointer :: w_dist(:,:), dw_dist(:,:), ddw_dist(:,:)
 
   complex(cp), public, pointer :: z_LMloc(:,:),dz_LMloc(:,:)
   complex(cp), public, pointer :: z_dist(:,:), dz_dist(:,:)
 
   !-- Entropy:
   complex(cp), public, allocatable, target :: s_LMloc_container(:,:,:)
   complex(cp), public, allocatable, target :: s_dist_container(:,:,:)
   complex(cp), public, pointer :: s_LMloc(:,:), ds_LMloc(:,:)
   complex(cp), public, pointer :: s_dist(:,:), ds_dist(:,:)
 
   !-- Chemical composition:
   complex(cp), public, allocatable, target :: xi_LMloc_container(:,:,:)
   complex(cp), public, allocatable, target :: xi_dist_container(:,:,:)
   complex(cp), public, pointer :: xi_LMloc(:,:), dxi_LMloc(:,:)
   complex(cp), public, pointer :: xi_dist(:,:), dxi_dist(:,:)

   !-- Pressure:
   complex(cp), public, pointer :: p_LMloc(:,:), dp_LMloc(:,:)
   complex(cp), public, pointer :: p_dist(:,:), dp_dist(:,:)
 
   !-- Magnetic field potentials:
   complex(cp), public, allocatable :: b(:,:)
   complex(cp), public, allocatable, target :: field_LMloc_container(:,:,:)
   complex(cp), public, allocatable, target :: field_dist_container(:,:,:)
   complex(cp), public, pointer :: b_LMloc(:,:), db_LMloc(:,:), ddb_LMloc(:,:)
   complex(cp), public, pointer :: b_dist(:,:), db_dist(:,:), ddb_dist(:,:)
   complex(cp), public, pointer :: aj_LMloc(:,:), dj_LMloc(:,:), ddj_LMloc(:,:)
   complex(cp), public, pointer :: aj_dist(:,:), dj_dist(:,:)
 
   !-- Magnetic field potentials in inner core:
   !   NOTE: n_r-dimension may be smaller once CHEBFT is addopted
   !         for even chebs
   complex(cp), public, allocatable :: b_ic(:,:)  
   complex(cp), public, allocatable :: db_ic(:,:)
   complex(cp), public, allocatable :: ddb_ic(:,:)
   complex(cp), public, allocatable :: aj_ic(:,:) 
   complex(cp), public, allocatable :: dj_ic(:,:)
   complex(cp), public, allocatable :: ddj_ic(:,:)
   complex(cp), public, allocatable :: b_ic_LMloc(:,:)  
   complex(cp), public, allocatable :: db_ic_LMloc(:,:)
   complex(cp), public, allocatable :: ddb_ic_LMloc(:,:)
   complex(cp), public, allocatable :: aj_ic_LMloc(:,:) 
   complex(cp), public, allocatable :: dj_ic_LMloc(:,:)
   complex(cp), public, allocatable :: ddj_ic_LMloc(:,:)

   complex(cp), public, allocatable :: work_LMloc(:,:) ! Needed in update routines
   
   !-- Rotation rates:
   real(cp), public :: omega_ic,omega_ma

   public :: initialize_fields, finalize_fields

contains

   subroutine initialize_fields

      !-- Velocity potentials:
      if ( coord_r == 0 ) then
         allocate( b(lm_maxMag,n_r_maxMag) )
         bytes_allocated = bytes_allocated +  &
                           lm_maxMag*n_r_maxMag*SIZEOF_DEF_COMPLEX
         allocate( b_ic(lm_maxMag,n_r_ic_maxMag) )  
         allocate( db_ic(lm_maxMag,n_r_ic_maxMag) )
         allocate( ddb_ic(lm_maxMag,n_r_ic_maxMag) )
         allocate( aj_ic(lm_maxMag,n_r_ic_maxMag) ) 
         allocate( dj_ic(lm_maxMag,n_r_ic_maxMag) )
         allocate( ddj_ic(lm_maxMag,n_r_ic_maxMag) )
         bytes_allocated = bytes_allocated + &
                           6*lm_maxMag*n_r_ic_maxMag*SIZEOF_DEF_COMPLEX

      else
         allocate( b(1,n_r_maxMag) )
         bytes_allocated = bytes_allocated + n_r_maxMag*SIZEOF_DEF_COMPLEX
         allocate( b_ic(1,n_r_ic_maxMag) )  
         allocate( db_ic(1,n_r_ic_maxMag) )
         allocate( ddb_ic(1,n_r_ic_maxMag) )
         allocate( aj_ic(1,n_r_ic_maxMag) ) 
         allocate( dj_ic(1,n_r_ic_maxMag) )
         allocate( ddj_ic(1,n_r_ic_maxMag) )
         bytes_allocated = bytes_allocated + 6*n_r_ic_maxMag*SIZEOF_DEF_COMPLEX
      end if
      allocate( flow_LMloc_container(llm:ulm,n_r_max,1:7) )
      w_LMloc(llm:ulm,1:n_r_max)   => flow_LMloc_container(:,:,1)
      dw_LMloc(llm:ulm,1:n_r_max)  => flow_LMloc_container(:,:,2)
      ddw_LMloc(llm:ulm,1:n_r_max) => flow_LMloc_container(:,:,3)
      z_LMloc(llm:ulm,1:n_r_max)   => flow_LMloc_container(:,:,4)
      dz_LMloc(llm:ulm,1:n_r_max)  => flow_LMloc_container(:,:,5)
      p_LMloc(llm:ulm,1:n_r_max)   => flow_LMloc_container(:,:,6)
      dp_LMloc(llm:ulm,1:n_r_max)  => flow_LMloc_container(:,:,7)

      allocate( flow_dist_container(lm_loc,nRstart:nRstop,1:7) )
      w_dist(1:lm_loc,nRstart:nRstop)   => flow_dist_container(:,:,1)
      dw_dist(1:lm_loc,nRstart:nRstop)  => flow_dist_container(:,:,2)
      ddw_dist(1:lm_loc,nRstart:nRstop) => flow_dist_container(:,:,3)
      z_dist(1:lm_loc,nRstart:nRstop)   => flow_dist_container(:,:,4)
      dz_dist(1:lm_loc,nRstart:nRstop)  => flow_dist_container(:,:,5)
      p_dist(1:lm_loc,nRstart:nRstop)   => flow_dist_container(:,:,6)
      dp_dist(1:lm_loc,nRstart:nRstop)  => flow_dist_container(:,:,7)

      !-- Entropy:
      allocate( s_LMloc_container(llm:ulm,n_r_max,1:2) )
      s_LMloc(llm:ulm,1:n_r_max)  => s_LMloc_container(:,:,1)
      ds_LMloc(llm:ulm,1:n_r_max) => s_LMloc_container(:,:,2)
      allocate( s_dist_container(lm_loc,nRstart:nRstop,1:2) )
      s_dist(1:lm_loc,nRstart:nRstop)   => s_dist_container(:,:,1)
      ds_dist(1:lm_loc,nRstart:nRstop)  => s_dist_container(:,:,2)

      bytes_allocated = bytes_allocated + &
                        9*(ulm-llm+1)*n_r_max*SIZEOF_DEF_COMPLEX
      bytes_allocated = bytes_allocated + &
                        9*lm_loc*(nRstop-nRstart+1)*SIZEOF_DEF_COMPLEX

      !-- Chemical composition:
      if ( l_chemical_conv ) then
         allocate( xi_LMloc_container(llm:ulm,n_r_max,1:2) )
         xi_LMloc(llm:ulm,1:n_r_max)  => xi_LMloc_container(:,:,1)
         dxi_LMloc(llm:ulm,1:n_r_max) => xi_LMloc_container(:,:,2)
         allocate( xi_dist_container(lm_loc,nRstart:nRstop,1:2) )
         xi_dist(1:lm_loc,nRstart:nRstop)   => xi_dist_container(:,:,1)
         dxi_dist(1:lm_loc,nRstart:nRstop)  => xi_dist_container(:,:,2)
         bytes_allocated = bytes_allocated + &
                           2*(ulm-llm+1)*n_r_max*SIZEOF_DEF_COMPLEX
         bytes_allocated = bytes_allocated + &
                           2*lm_loc*(nRstop-nRstart+1)*SIZEOF_DEF_COMPLEX
      else
         allocate( xi_LMloc_container(1,1,2) ) ! For debugging
         xi_LMloc(1:1,1:1)  => xi_LMloc_container(:,:,1)
         dxi_LMloc(1:1,1:1) => xi_LMloc_container(:,:,2)
         allocate( xi_dist_container(1,1,2) )
         xi_dist(1:1,1:1)   => xi_dist_container(:,:,1)
         dxi_dist(1:1,1:1)  => xi_dist_container(:,:,2)
      end if

      !-- Magnetic field potentials:
      allocate( field_LMloc_container(llmMag:ulmMag,n_r_maxMag,1:6) )
      b_LMloc(llmMag:ulmMag,1:n_r_maxMag)   => field_LMloc_container(:,:,1)
      db_LMloc(llmMag:ulmMag,1:n_r_maxMag)  => field_LMloc_container(:,:,2)
      ddb_LMloc(llmMag:ulmMag,1:n_r_maxMag) => field_LMloc_container(:,:,3)
      aj_LMloc(llmMag:ulmMag,1:n_r_maxMag)  => field_LMloc_container(:,:,4)
      dj_LMloc(llmMag:ulmMag,1:n_r_maxMag)  => field_LMloc_container(:,:,5)
      ddj_LMloc(llmMag:ulmMag,1:n_r_maxMag) => field_LMloc_container(:,:,6)
      bytes_allocated = bytes_allocated + &
                        6*(ulmMag-llmMag+1)*n_r_maxMag*SIZEOF_DEF_COMPLEX

      allocate( field_dist_container(lm_locMag,nRstart:nRstop,1:5) )
      b_dist(1:lm_locMag,nRstart:nRstop)   => field_dist_container(:,:,1)
      db_dist(1:lm_locMag,nRstart:nRstop)  => field_dist_container(:,:,2)
      ddb_dist(1:lm_locMag,nRstart:nRstop) => field_dist_container(:,:,3)
      aj_dist(1:lm_locMag,nRstart:nRstop)  => field_dist_container(:,:,4)
      dj_dist(1:lm_locMag,nRstart:nRstop)  => field_dist_container(:,:,5)
      bytes_allocated = bytes_allocated + &
                        5*lm_locMag*(nRstop-nRstart+1)*SIZEOF_DEF_COMPLEX

      !-- Magnetic field potentials in inner core:
      !   NOTE: n_r-dimension may be smaller once CHEBFT is adopted
      !         for even chebs
      allocate( b_ic_LMloc(llmMag:ulmMag,n_r_ic_maxMag) )  
      allocate( db_ic_LMloc(llmMag:ulmMag,n_r_ic_maxMag) )
      allocate( ddb_ic_LMloc(llmMag:ulmMag,n_r_ic_maxMag) )
      allocate( aj_ic_LMloc(llmMag:ulmMag,n_r_ic_maxMag) ) 
      allocate( dj_ic_LMloc(llmMag:ulmMag,n_r_ic_maxMag) )
      allocate( ddj_ic_LMloc(llmMag:ulmMag,n_r_ic_maxMag) )
      bytes_allocated = bytes_allocated + &
                        6*(ulmMag-llmMag+1)*n_r_ic_maxMag*SIZEOF_DEF_COMPLEX
      
      allocate( work_LMloc(llm:ulm,1:n_r_max) )
      bytes_allocated = bytes_allocated + (ulm-llm+1)*n_r_max*SIZEOF_DEF_COMPLEX

   end subroutine initialize_fields
!----------------------------------------------------------------------------
   subroutine finalize_fields

      deallocate( b, b_ic, db_ic, ddb_ic )
      deallocate( aj_ic, dj_ic, ddj_ic, flow_LMloc_container )
      deallocate( flow_dist_container, s_LMloc_container, s_dist_container )
      deallocate( field_LMloc_container, field_dist_container )
      deallocate( b_ic_LMloc )
      deallocate( db_ic_LMloc )
      deallocate( ddb_ic_LMloc )
      deallocate( aj_ic_LMloc )
      deallocate( dj_ic_LMloc, ddj_ic_LMloc )
      deallocate( xi_LMloc_container, xi_dist_container )
      deallocate( work_LMloc )
      
      nullify( xi_dist )
      nullify( dxi_dist)
      nullify( s_dist  )
      nullify( ds_dist )
      nullify( z_dist  )
      nullify( dz_dist )
      nullify( p_dist  )
      nullify( dp_dist )
      nullify( b_dist  )
      nullify( db_dist )
      nullify( ddb_dist)
      nullify( aj_dist )
      nullify( dj_dist )
      nullify( w_dist  )
      nullify( dw_dist )
      nullify( ddw_dist)

   end subroutine finalize_fields
!----------------------------------------------------------------------------
end module fields
