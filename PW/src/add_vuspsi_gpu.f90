!
! Copyright (C) 2001-2015 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE add_vuspsi_gpu( lda, n, m, hpsi_d )
  !----------------------------------------------------------------------------
  !! This routine applies the Ultra-Soft Hamiltonian to a
  !! vector psi and puts the result in hpsi.
  !! It requires the products of psi with all beta functions
  !! in array becp(nkb,m) (calculated by calbec).
  !
  USE kinds,           ONLY: DP
  USE ions_base,       ONLY: nat, ntyp => nsp, ityp
  USE lsda_mod,        ONLY: current_spin
  USE control_flags,   ONLY: gamma_only
  USE noncollin_module
  USE uspp,            ONLY : ofsbeta, nkb, vkb
  USE uspp_param,      ONLY : nh, nhm
#if defined(__OPENMP_GPU)
  USE uspp,            ONLY : deeq, deeq_nc
  USE becmod,          ONLY : bec_type, becp
#else
  USE uspp,            ONLY : deeq_d, deeq_nc_d
  USE becmod_gpum,     ONLY : bec_type_d, becp_d
#endif
  USE becmod_gpum,     ONLY : using_becp_r_d, using_becp_k_d, using_becp_nc_d
  !
  IMPLICIT NONE
  !
  ! ... I/O variables
  !
  INTEGER, INTENT(IN) :: lda
  !! leading dimension of arrays psi, spsi
  INTEGER, INTENT(IN) :: n
  !! true dimension of psi, spsi
  INTEGER, INTENT(IN) :: m
  !! number of states psi
  COMPLEX(DP), INTENT(INOUT) :: hpsi_d(lda*npol,m)
  !! V_US|psi> is added to hpsi
  !
  ! ... here the local variables
  !
#if defined(__CUDA) || defined(__OPENMP_GPU)
#if defined(__CUDA)
  attributes(DEVICE) :: hpsi_d
#endif
  !
  !
  CALL start_clock_gpu( 'add_vuspsi' )
  !
  IF ( gamma_only ) THEN
     !
     CALL add_vuspsi_gamma_gpu()
     !
  ELSE IF ( noncolin) THEN
     !
     CALL add_vuspsi_nc_gpu ()
     !
  ELSE
     !
     CALL add_vuspsi_k_gpu()
     !
  END IF
  !
  CALL stop_clock_gpu( 'add_vuspsi' )
  !
  RETURN
  !
  CONTAINS
     !
     !-----------------------------------------------------------------------
     SUBROUTINE add_vuspsi_gamma_gpu()
       !-----------------------------------------------------------------------
       !! See comments inside
       !
       USE mp,              ONLY : mp_get_comm_null, mp_circular_shift_left
       USE devxlib_buffers, ONLY : dev_buf=>gpu_buffer
       USE distools,        ONLY : ldim_block, gind_block
       !
#if defined(__CUDA)
       USE cudafor
       USE cublas
#elif defined(__OPENMP_GPU)
#if defined(MKL_ILP64)
       USE onemkl_blas_omp_offload_ilp64_no_array_check
#else
       USE onemkl_blas_omp_offload_lp64_no_array_check
#endif
#endif
       !
       IMPLICIT NONE
       REAL(DP), POINTER :: ps_d (:,:)
       INTEGER :: ierr
       INTEGER :: nproc, mype, m_loc, m_begin, ibnd_loc, icyc, icur_blk, m_max
       ! counters
       INTEGER :: jkb, ikb, ih, jh, na, nt, ibnd
       !
#if defined(__CUDA)
       attributes(device) :: ps_d
#endif
       !
       IF ( nkb == 0 ) RETURN
       !
       CALL using_becp_r_d(0)
       !
#if defined(__OPENMP_GPU)
       ASSOCIATE(becp_d=>becp, deeq_d=>deeq, r=>becp%r)
#endif
       IF( becp_d%comm == mp_get_comm_null() ) THEN
          nproc   = 1
          mype    = 0
          m_loc   = m
          m_begin = 1
          m_max   = m
       ELSE
          !
          ! becp(l,i) = <beta_l|psi_i>, with vkb(n,l)=|beta_l>
          ! in this case becp(l,i) are distributed (index i is)
          !
          nproc   = becp_d%nproc
          mype    = becp_d%mype
          m_loc   = becp_d%nbnd_loc
          m_begin = becp_d%ibnd_begin
#if defined(__OPENMP_GPU)
          m_max   = SIZE( becp_d%r, 2 )
#else
          m_max   = SIZE( becp_d%r_d, 2 )
#endif
          IF( ( m_begin + m_loc - 1 ) > m ) m_loc = m - m_begin + 1
       END IF
       !
       CALL dev_buf%lock_buffer(ps_d, (/ nkb, m_max /), ierr ) !ALLOCATE (ps_d (nkb,m_max), STAT=ierr )
       !
       IF( ierr /= 0 ) &
          CALL errore( ' add_vuspsi_gamma ', ' cannot allocate ps_d ', ABS(ierr) )
       !
#if defined(__OPENMP_GPU)
       !$omp target teams loop collapse(2)
       do jh=1,m_max
          do ih=1,nkb
             ps_d(ih,jh) = 0.D0
          enddo
       enddo
#else
       ps_d(:,:) = 0.D0
#endif
       !
       !   In becp=<vkb_i|psi_j> terms corresponding to atom na of type nt
       !   run from index i=ofsbeta(na)+1 to i=ofsbeta(na)+nh(nt)
       !
       DO nt = 1, ntyp
          !
          IF ( nh(nt) == 0 ) CYCLE
          DO na = 1, nat
             !
             IF ( ityp(na) == nt ) THEN
                !
                ! Next operation computes ps(l',i) = \sum_m deeq(l,m) becp(m',i)
                ! (l'=l+ijkb0, m'=m+ijkb0, indices run from 1 to nh(nt))
                !
                IF ( m_loc > 0 ) THEN
                  !$omp target variant dispatch use_device_ptr(deeq_d,r)
                  CALL DGEMM('N', 'N', nh(nt), m_loc, nh(nt), 1.0_dp, &
                           deeq_d(1,1,na,current_spin), nhm, &
#if defined(__OPENMP_GPU)
                           r(ofsbeta(na)+1,1), nkb, 0.0_dp, &
#else
                           becp_d%r_d(ofsbeta(na)+1,1), nkb, 0.0_dp, &
#endif
                               ps_d(ofsbeta(na)+1,1), nkb )
                  !$omp end target variant dispatch
                END IF
                !
             END IF
             !
          END DO
          !
       END DO
       !
       IF( becp_d%comm == mp_get_comm_null() ) THEN
          !
          ! Normal case: hpsi(n,i) = \sum_l beta(n,l) ps(l,i)
          ! (l runs from 1 to nkb)
          !
!$acc data present(vkb(:,:))
!$acc host_data use_device(vkb)
          !$omp target variant dispatch use_device_ptr(vkb,hpsi_d)
          CALL DGEMM( 'N', 'N', ( 2 * n ), m, nkb, 1.D0, vkb, &
                   ( 2 * lda ), ps_d, nkb, 1.D0, hpsi_d, ( 2 * lda ) )
          !$omp end target variant dispatch
!$acc end host_data
!$acc end data
       ELSE
          !
          ! parallel block multiplication of vkb and ps
          !
          icur_blk = mype
          !
          DO icyc = 0, nproc - 1

             !$omp target update from(becp%nbnd)
             m_loc   = ldim_block( becp_d%nbnd, nproc, icur_blk )
             m_begin = gind_block( 1,  becp_d%nbnd, nproc, icur_blk )

             IF( ( m_begin + m_loc - 1 ) > m ) m_loc = m - m_begin + 1

             IF( m_loc > 0 ) THEN
!$acc data present(vkb(:,:))
!$acc host_data use_device(vkb)
                !$omp target variant dispatch use_device_ptr(vkb,hpsi_d)
                CALL DGEMM( 'N', 'N', ( 2 * n ), m_loc, nkb, 1.D0, vkb, &
                   ( 2 * lda ), ps_d, nkb, 1.D0, hpsi_d( 1, m_begin ), ( 2 * lda ) )
                !$omp end target variant dispatch
!$acc end host_data
!$acc end data
             ENDIF

             ! block rotation
             !
             CALL mp_circular_shift_left( ps_d, icyc, becp_d%comm )

             icur_blk = icur_blk + 1
             IF( icur_blk == nproc ) icur_blk = 0

          ENDDO
       ENDIF
       !
       CALL dev_buf%release_buffer(ps_d, ierr) ! DEALLOCATE (ps_d)
#if defined(__OPENMP_GPU)
       ENDASSOCIATE
#endif
       !
       RETURN
       !
     END SUBROUTINE add_vuspsi_gamma_gpu
     !
     !-----------------------------------------------------------------------
     SUBROUTINE add_vuspsi_k_gpu()
       !-----------------------------------------------------------------------
       !! See add_vuspsi_gamma for comments
       !
#if defined(__CUDA)
       USE cudafor
       USE cublas
#elif defined(__OPENMP_GPU)
#if defined(MKL_ILP64)
       USE onemkl_blas_omp_offload_ilp64_no_array_check
#else
       USE onemkl_blas_omp_offload_lp64_no_array_check
#endif
#endif
       USE devxlib_buffers, ONLY : dev_buf=>gpu_buffer
       use omp_lib
       use iso_c_binding
       !
       IMPLICIT NONE
       COMPLEX(DP), POINTER :: ps_d (:,:), deeaux_d (:,:)
       INTEGER :: ierr
       ! counters
       INTEGER :: i, j, k, jkb, ikb, ih, jh, na, nt, ibnd, nhnt
       COMPLEX(DP) :: temp1, temp2

#if defined(__CUDA)
       ATTRIBUTES( DEVICE ) :: ps_d, deeaux_d
#endif
       !
       IF ( nkb == 0 ) RETURN
       !
       CALL using_becp_k_d(0)
       !
       CALL dev_buf%lock_buffer(ps_d, (/ nkb,m /), ierr ) ! ALLOCATE (ps_d (nkb,m), STAT=ierr )
       IF( ierr /= 0 ) &
          CALL errore( ' add_vuspsi_k ', ' cannot allocate ps_d ', ABS( ierr ) )
       !
       CALL dev_buf%lock_buffer(deeaux_d, (/ nhm, nhm /), ierr ) !ALLOCATE ( deeaux_d(nhm, nhm) )
       IF( ierr /= 0 ) &
          CALL errore( ' add_vuspsi_k ', ' cannot allocate deeaux_d ', ABS( ierr ) )
       !
       DO nt = 1, ntyp
          !
          IF ( nh(nt) == 0 ) CYCLE
          !
          nhnt = nh(nt)
          !
          DO na = 1, nat
             !
             IF ( ityp(na) == nt ) THEN
                !
                ! deeq is real: copy it into a complex variable to perform
                ! a zgemm - simple but sub-optimal solution
                !
                !deeaux_d(:,:) = CMPLX(deeq_d(1:nh(nt),1:nh(nt),na,current_spin), 0.0_dp, KIND=dp )
                !
!$cuf kernel do(2) <<<*,*>>>
                !$omp target teams distribute parallel do collapse(2)
                DO j = 1, nhnt
                   DO k = 1, nhnt
#if defined(__OPENMP_GPU)
                      deeaux_d(k,j) = CMPLX(deeq(k,j,na,current_spin), 0.0_dp, KIND=DP )
#else
                      deeaux_d(k,j) = CMPLX(deeq_d(k,j,na,current_spin), 0.0_dp, KIND=DP )
#endif
                   END DO
                END DO
                !
#if defined(__OPENMP_GPU)
                associate(k => becp%k)
#endif
                !$omp target variant dispatch use_device_ptr(k, ofsbeta)
                CALL ZGEMM('N','N', nh(nt), m, nh(nt), (1.0_dp,0.0_dp), &
#if defined(__OPENMP_GPU)
                           deeaux_d, nhm, k(ofsbeta(na)+1,1), nkb, &
#else
                           deeaux_d, nhm, becp_d%k_d(ofsbeta(na)+1,1), nkb, &
#endif
                          (0.0_dp, 0.0_dp), ps_d(ofsbeta(na)+1,1), nkb )
                !$omp end target variant dispatch
#if defined(__OPENMP_GPU)
                endassociate
#endif
                !
             END IF
             !
          END DO
          !
       END DO
       CALL dev_buf%release_buffer(deeaux_d, ierr) ! DEALLOCATE (deeaux_d)
       !
!$acc data present(vkb(:,:))
!$acc host_data use_device(vkb)
       !$omp target variant dispatch use_device_ptr(vkb,hpsi_d)
       CALL ZGEMM( 'N', 'N', n, m, nkb, ( 1.D0, 0.D0 ) , vkb, &
                   lda, ps_d, nkb, ( 1.D0, 0.D0 ) , hpsi_d, lda )
       !$omp end target variant dispatch
!$acc end host_data
!$acc end data
       !
       CALL dev_buf%release_buffer(ps_d, ierr) !DEALLOCATE (ps_d)
       !
       RETURN
       !
     END SUBROUTINE add_vuspsi_k_gpu
     !
     !-----------------------------------------------------------------------
     SUBROUTINE add_vuspsi_nc_gpu()
       !-----------------------------------------------------------------------
       !! See add_vuspsi_k for comments
       !
#if defined(__CUDA)
       USE cudafor
       USE cublas
#elif defined(__OPENMP_GPU)
#if defined(MKL_ILP64)
       USE onemkl_blas_omp_offload_ilp64_no_array_check
#else
       USE onemkl_blas_omp_offload_lp64_no_array_check
#endif
#endif
       USE devxlib_buffers, ONLY : dev_buf=>gpu_buffer
#if !defined(__OPENMP_GPU)
       USE becmod_gpum,   ONLY : becp_d
#else
       USE becmod,        ONLY : becp
#endif
       IMPLICIT NONE
       COMPLEX(DP), POINTER :: ps_d (:,:,:)
       INTEGER :: ierr
       ! counters
       INTEGER :: na, nt, ibnd

#if defined(__CUDA)
       ATTRIBUTES( DEVICE ) :: ps_d
#endif
       !
       IF ( nkb == 0 ) RETURN
#if defined(__OPENMP_GPU)
       ASSOCIATE(becp_d=>becp, deeq_nc_d=>deeq, nc => becp%nc)
#endif
       !
       CALL using_becp_nc_d(0)
       !
       ! ALLOCATE (ps_d( nkb, npol, m), STAT=ierr )
       CALL dev_buf%lock_buffer(ps_d, (/ nkb, npol, m /), ierr )
       !
       IF( ierr /= 0 ) &
          CALL errore( ' add_vuspsi_nc ', ' error allocating ps_d ', ABS( ierr ) )
       !
       !  OPTIMIZE HERE: possibly streamline
       !
       DO nt = 1, ntyp
          !
          IF ( nh(nt) == 0 ) CYCLE
          DO na = 1, nat
             !
             IF ( ityp(na) == nt ) THEN
                !
                !$omp target variant dispatch use_device_ptr(deeq_nc_d,nc)
                CALL ZGEMM('N','N', nh(nt), m, nh(nt), (1.0_dp,0.0_dp), &
#if defined (__OPENMP_GPU)
                           deeq_nc_d(1,1,na,1), nhm, nc(ofsbeta(na)+1,1,1), 2*nkb, &
#else
                           deeq_nc_d(1,1,na,1), nhm, becp_d%nc_d(ofsbeta(na)+1,1,1), 2*nkb, &
#endif
                          (0.0_dp, 0.0_dp), ps_d(ofsbeta(na)+1,1,1), 2*nkb )

                CALL ZGEMM('N','N', nh(nt), m, nh(nt), (1.0_dp,0.0_dp), &
#if defined (__OPENMP_GPU)
                           deeq_nc_d(1,1,na,2), nhm, nc(ofsbeta(na)+1,2,1), 2*nkb, &
#else
                           deeq_nc_d(1,1,na,2), nhm, becp_d%nc_d(ofsbeta(na)+1,2,1), 2*nkb, &
#endif
                          (1.0_dp, 0.0_dp), ps_d(ofsbeta(na)+1,1,1), 2*nkb )

                CALL ZGEMM('N','N', nh(nt), m, nh(nt), (1.0_dp,0.0_dp), &
#if defined (__OPENMP_GPU)
                           deeq_nc_d(1,1,na,3), nhm, nc(ofsbeta(na)+1,1,1), 2*nkb, &
#else
                           deeq_nc_d(1,1,na,3), nhm, becp_d%nc_d(ofsbeta(na)+1,1,1), 2*nkb, &
#endif
                          (0.0_dp, 0.0_dp), ps_d(ofsbeta(na)+1,2,1), 2*nkb )

                CALL ZGEMM('N','N', nh(nt), m, nh(nt), (1.0_dp,0.0_dp), &
#if defined (__OPENMP_GPU)
                           deeq_nc_d(1,1,na,4), nhm, nc(ofsbeta(na)+1,2,1), 2*nkb, &
#else
                           deeq_nc_d(1,1,na,4), nhm, becp_d%nc_d(ofsbeta(na)+1,2,1), 2*nkb, &
#endif
                          (1.0_dp, 0.0_dp), ps_d(ofsbeta(na)+1,2,1), 2*nkb )
                !$omp end target variant dispatch

!                DO ibnd = 1, m
!                   !
!                   DO jh = 1, nh(nt)
!                      !
!!$cuf kernel do(1) <<<*,*>>>
!                      DO ih = 1, nh(nt)
!                         !
!                         ikb = ijkb0 + ih
!                         jkb = ijkb0 + jh
!                         becpup_jkb = becp_nc_d(jkb,1,ibnd)
!                         becpdn_jkb = becp_nc_d(jkb,2,ibnd)
!                         !
!                         ps_d(ikb,1,ibnd) = ps_d(ikb,1,ibnd) +   &
!                              deeq_nc_d(ih,jh,na,1)*becpup_jkb + &
!                              deeq_nc_d(ih,jh,na,2)*becpdn_jkb
!                         ps_d(ikb,2,ibnd) = ps_d(ikb,2,ibnd) +   &
!                              deeq_nc_d(ih,jh,na,3)*becpup_jkb + &
!                              deeq_nc_d(ih,jh,na,4)*becpdn_jkb
!                         !
!                      END DO
!                      !
!                   END DO
!                   !
!                END DO
                !
             END IF
             !
          END DO
          !
       END DO
       !
!$acc data present(vkb(:,:))
!$acc host_data use_device(vkb)
       !$omp target variant dispatch use_device_ptr(vkb,hpsi_d)
       call ZGEMM ('N', 'N', n, m*npol, nkb, ( 1.D0, 0.D0 ) , vkb, &
                   lda, ps_d, nkb, ( 1.D0, 0.D0 ) , hpsi_d, lda )
       !$omp end target variant dispatch
!$acc end host_data
!$acc end data
       !
       CALL dev_buf%release_buffer(ps_d, ierr ) ! DEALLOCATE (ps_d)
#if defined(__OPENMP_GPU)
       ENDASSOCIATE
#endif
       !
       RETURN
       !
     END SUBROUTINE add_vuspsi_nc_gpu
#endif
     !
END SUBROUTINE add_vuspsi_gpu
