!
! Copyright (C) 2001-2015 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE sum_band_gpu()
  !----------------------------------------------------------------------------
  !! Calculates the symmetrized charge density and related quantities.
  !! Also computes the occupations and the sum of occupied eigenvalues.
  !
#if defined(__CUDA)
  USE cudafor
#endif
  USE kinds,                ONLY : DP
  USE ener,                 ONLY : eband
  USE control_flags,        ONLY : diago_full_acc, gamma_only, lxdm, tqr
  USE cell_base,            ONLY : at, bg, omega, tpiba
  USE ions_base,            ONLY : nat, ntyp => nsp, ityp
  USE fft_base,             ONLY : dfftp, dffts
  USE fft_interfaces,       ONLY : fwfft, invfft
  USE gvect,                ONLY : ngm, g
  USE gvecs,                ONLY : doublegrid
  USE klist,                ONLY : nks, nkstot, wk, xk, ngk, igk_k
#if !defined(__OPENMP_GPU)
  USE klist,                ONLY : igk_k_d
#endif
  USE ldaU,                 ONLY : lda_plus_u, lda_plus_u_kind, is_hubbard_back
  USE lsda_mod,             ONLY : lsda, nspin, current_spin, isk
  USE scf,                  ONLY : rho, rhoz_or_updw
  USE symme,                ONLY : sym_rho
  USE io_files,             ONLY : iunwfc, nwordwfc
  USE buffers,              ONLY : get_buffer
  USE uspp,                 ONLY : nkb, vkb, becsum, ebecsum, nhtol, nhtoj, indv, okvan
#if !defined(__OPENMP_GPU)
  USE uspp,                 ONLY : becsum_d, ebecsum_d, vkb_d
#endif
  USE uspp_param,           ONLY : upf, nh, nhm
  USE wavefunctions,        ONLY : evc, psic
  USE noncollin_module,     ONLY : noncolin, npol, nspin_mag, domag
  USE wvfct,                ONLY : nbnd, npwx, wg, et, btype
  USE mp_pools,             ONLY : inter_pool_comm
  USE mp_bands,             ONLY : inter_bgrp_comm, intra_bgrp_comm, nbgrp
  USE mp,                   ONLY : mp_sum
  USE xc_lib,               ONLY : xclib_dft_is
  USE paw_symmetry,         ONLY : PAW_symmetrize
  USE paw_variables,        ONLY : okpaw
  USE becmod,               ONLY : allocate_bec_type, deallocate_bec_type, &
                                   becp
  USE gcscf_module,         ONLY : lgcscf, gcscf_calc_nelec
#if !defined(__OPENMP_GPU)
  USE wavefunctions_gpum,   ONLY : evc_d
#endif
  USE wavefunctions_gpum,   ONLY : using_evc, using_evc_d
  USE wvfct_gpum,           ONLY : using_et
  USE becmod_subs_gpum,     ONLY : using_becp_auto
  !
  IMPLICIT NONE
  !
  ! ... local variables
  !
  INTEGER :: ir,   &! counter on 3D r points
             is,   &! counter on spin polarizations
             ig,   &! counter on g vectors
             ibnd, &! counter on bands
             ik,   &! counter on k points
             nt,   &! counter on atomic types
             npol_,&! auxiliary dimension for noncolin case
             ibnd_start, ibnd_end, this_bgrp_nbnd ! first, last and number of band in this bgrp
  REAL (DP), ALLOCATABLE :: kplusg (:)
  !
  !
  CALL start_clock_gpu( 'sum_band' )
  !
  if ( nhm > 0 ) then
     becsum(:,:,:) = 0.D0
     if (tqr) ebecsum(:,:,:) = 0.D0
#if !defined(__OPENMP_GPU)
     becsum_d(:,:,:) = 0.D0
     if (tqr) ebecsum_d(:,:,:) = 0.D0
#else
     !$omp target update to(becsum)
     if (tqr) then
        !$omp target update to(ebecsum)
     endif
#endif
  end if
  rho%of_r(:,:)      = 0.D0
  rho%of_g(:,:)      = (0.D0, 0.D0)
  if ( xclib_dft_is('meta') .OR. lxdm ) then
     rho%kin_r(:,:)      = 0.D0
     rho%kin_g(:,:)      = (0.D0, 0.D0)
  end if
  !
  ! ... calculates weights of Kohn-Sham orbitals used in calculation of rho
  !
  CALL start_clock_gpu( 'sum_band:weights' )
  CALL weights ( )
  CALL stop_clock_gpu( 'sum_band:weights' )
  !
  ! ... btype, used in diagonalization, is set here: a band is considered empty
  ! ... and computed with low accuracy only when its occupation is < 0.01, and
  ! ... only if option diago_full_acc is false; otherwise, use full accuracy
  !
  btype(:,:) = 1
  IF ( .NOT. diago_full_acc ) THEN
     !
     FORALL( ik = 1:nks, wk(ik) > 0.D0 )
        WHERE( wg(:,ik) / wk(ik) < 0.01D0 ) btype(:,ik) = 0
     END FORALL
     !
  END IF
  !
  ! ... Needed for LDA+U: compute occupations of Hubbard states
  !
  IF (lda_plus_u) THEN
    IF (lda_plus_u_kind.EQ.0) THEN
       !
       CALL new_ns(rho%ns)
       !
       DO nt = 1, ntyp
          IF (is_hubbard_back(nt)) CALL new_nsb(rho%nsb)
       ENDDO
       !
    ELSEIF (lda_plus_u_kind.EQ.1) THEN
       !
       IF (noncolin) THEN
          CALL new_ns_nc(rho%ns_nc)
       ELSE
          CALL new_ns(rho%ns)
       ENDIF
       !
    ELSEIF (lda_plus_u_kind.EQ.2) THEN
       !
       CALL new_nsg()
       !
    ENDIF
  ENDIF
  !
  ! ... for band parallelization: set band computed by this processor
  !
  call divide ( inter_bgrp_comm, nbnd, ibnd_start, ibnd_end )
  this_bgrp_nbnd = ibnd_end - ibnd_start + 1
  !
  ! ... Allocate (and later deallocate) arrays needed in specific cases
  !
  IF ( okvan ) CALL allocate_bec_type (nkb, this_bgrp_nbnd, becp, intra_bgrp_comm)
  IF ( okvan ) CALL using_becp_auto(2)
  IF (xclib_dft_is('meta') .OR. lxdm) ALLOCATE (kplusg(npwx))
  !
  ! ... specialized routines are called to sum at Gamma or for each k point
  ! ... the contribution of the wavefunctions to the charge
  ! ... The band energy contribution eband is computed together with the charge
  !
  eband         = 0.D0
  !
  CALL start_clock_gpu( 'sum_band:loop' )
  IF ( gamma_only ) THEN
     !
     CALL sum_band_gamma_gpu()
     !
  ELSE
     !
     CALL sum_band_k_gpu()
     !
  END IF
  CALL stop_clock_gpu( 'sum_band:loop' )
  CALL mp_sum( eband, inter_pool_comm )
  CALL mp_sum( eband, inter_bgrp_comm )
  !
  IF (xclib_dft_is('meta') .OR. lxdm) DEALLOCATE (kplusg)
  IF ( okvan ) CALL deallocate_bec_type ( becp )
  IF ( okvan ) CALL using_becp_auto(2)
  !
  ! ... sum charge density over pools (distributed k-points) and bands
  !
  CALL mp_sum( rho%of_r, inter_pool_comm )
  CALL mp_sum( rho%of_r, inter_bgrp_comm )
  IF ( noncolin .AND. .NOT. domag ) rho%of_r(:,2:4)=0.D0
  !
  ! ... bring the unsymmetrized rho(r) to G-space (use psic as work array)
  !
  DO is = 1, nspin
     psic(1:dffts%nnr) = rho%of_r(1:dffts%nnr,is)
     psic(dffts%nnr+1:) = 0.0_dp
     CALL fwfft (1, psic, dffts)
     rho%of_g(1:dffts%ngm,is) = psic(dffts%nl(1:dffts%ngm))
     rho%of_g(dffts%ngm+1:,is) = (0.0_dp,0.0_dp)
  END DO

  IF( okvan )  THEN
     !
     ! ... becsum is summed over bands (if bgrp_parallelization is done)
     ! ... and over k-points (but it is not symmetrized)
     !
     ! use host copy to do the comunication. This avoids going back an forth GPU data
     ! becsum=becsum_d     not needed
     ! since becsum is already uptodate, see sum_band*gpu
     !
     CALL mp_sum(becsum, inter_bgrp_comm )
     CALL mp_sum(becsum, inter_pool_comm )
#if !defined(__OPENMP_GPU)
     becsum_d=becsum
#else
     !$omp target update to(becsum)
#endif
     !
     ! ... same for ebecsum, a correction to becsum (?) in real space
     !
     IF (tqr) THEN
        !ebecsum=ebecsum_d   not needed as above
        CALL mp_sum(ebecsum, inter_pool_comm )
        CALL mp_sum(ebecsum, inter_bgrp_comm )
#if !defined(__OPENMP_GPU)
        ebecsum_d=ebecsum
#else
        !$omp target update to(ebecsum)
#endif
     ENDIF
     !
     ! ... PAW: symmetrize becsum and store it
     ! ... FIXME: the same should be done for USPP as well
     !
     IF ( okpaw ) THEN
        rho%bec(:,:,:) = becsum(:,:,:)
        CALL PAW_symmetrize(rho%bec)
     END IF
     !
     ! ... Here we add the (unsymmetrized) Ultrasoft contribution to the charge
     !
     CALL addusdens_gpu(rho%of_g(:,:))
     !
  ENDIF
  !
  ! ... symmetrize rho(G)
  !
  CALL start_clock_gpu( 'sum_band:sym_rho' )
  CALL sym_rho ( nspin_mag, rho%of_g )
  !
  ! ... synchronize rho%of_r to the calculated rho%of_g (use psic as work array)
  !
  DO is = 1, nspin_mag
     psic(:) = ( 0.D0, 0.D0 )
     psic(dfftp%nl(:)) = rho%of_g(:,is)
     IF ( gamma_only ) psic(dfftp%nlm(:)) = CONJG( rho%of_g(:,is) )
     CALL invfft (1, psic, dfftp)
     rho%of_r(:,is) = psic(:)
  END DO
  !
  ! ... rho_kin(r): sum over bands, k-points, bring to G-space, symmetrize,
  ! ... synchronize with rho_kin(G)
  !
  IF ( xclib_dft_is('meta') .OR. lxdm) THEN
     !
     CALL mp_sum( rho%kin_r, inter_pool_comm )
     CALL mp_sum( rho%kin_r, inter_bgrp_comm )
     DO is = 1, nspin
        psic(1:dffts%nnr) = rho%kin_r(1:dffts%nnr,is)
        psic(dffts%nnr+1:) = 0.0_dp
        CALL fwfft (1, psic, dffts)
        rho%kin_g(1:dffts%ngm,is) = psic(dffts%nl(1:dffts%ngm))
     END DO
     !
     IF (.NOT. gamma_only) CALL sym_rho( nspin, rho%kin_g )
     !
     DO is = 1, nspin
        psic(:) = ( 0.D0, 0.D0 )
        psic(dfftp%nl(:)) = rho%kin_g(:,is)
        IF ( gamma_only ) psic(dfftp%nlm(:)) = CONJG( rho%kin_g(:,is) )
        CALL invfft (1, psic, dfftp)
        rho%kin_r(:,is) = psic(:)
     END DO
     !
  END IF
  CALL stop_clock_gpu( 'sum_band:sym_rho' )
  !
  ! ... if LSDA rho%of_r and rho%of_g are converted from (up,dw) to
  ! ... (up+dw,up-dw) format.
  !
  IF ( nspin == 2 ) CALL rhoz_or_updw( rho, 'r_and_g', '->rhoz' )
  !
  ! ... sum number of electrons, for GC-SCF
  !
  IF ( lgcscf ) CALL gcscf_calc_nelec()
  !
  CALL stop_clock( 'sum_band' )
  !
  RETURN
  !
  CONTAINS
     !
     ! ... internal procedures
     !
     !-----------------------------------------------------------------------
     SUBROUTINE sum_band_gamma_gpu()
       !-----------------------------------------------------------------------
       !! \(\texttt{sum_band}\) - part for gamma version.
       !
       USE becmod,        ONLY : becp
       USE mp_bands,      ONLY : me_bgrp
       USE mp,            ONLY : mp_sum, mp_get_comm_null
       USE fft_helper_subroutines
       USE uspp_init,     ONLY : init_us_2
#if defined(__OPENMP_GPU)
       USE omp_lib
#endif
       !
       IMPLICIT NONE
       !
       ! ... local variables
       !
       REAL(DP) :: w1, w2
         ! weights
       INTEGER  :: npw, idx, ioff, ioff_tg, nxyp, incr, v_siz, i, j
       COMPLEX(DP), ALLOCATABLE :: tg_psi_d(:)
#if !defined(__OPENMP_GPU)
       COMPLEX(DP), ALLOCATABLE :: psic_d(:)
       REAL(DP),    ALLOCATABLE :: tg_rho_d(:), tg_rho_h(:)
       REAL(DP),    ALLOCATABLE :: rho_d(:,:)
       INTEGER,     POINTER     :: dffts_nl_d(:), dffts_nlm_d(:)
#else
       REAL(DP),    ALLOCATABLE :: tg_rho(:)
#endif
       LOGICAL :: use_tg
       INTEGER :: right_nnr, right_nr3, right_inc, ntgrp, ierr
#if defined(__CUDA)
       attributes(device) :: psic_d, tg_psi_d, tg_rho_d, rho_d
       attributes(device) :: dffts_nl_d, dffts_nlm_d
       attributes(pinned) :: tg_rho_h
#endif
       !
       CALL using_evc_d(0); CALL using_et(0)
#if !defined(__OPENMP_GPU)
       dffts_nl_d  => dffts%nl_d
       dffts_nlm_d => dffts%nlm_d
#else
       !omp_device = omp_get_default_device()
       associate(psic_d=>psic)
#endif
       !
       ! ... here we sum for each k point the contribution
       ! ... of the wavefunctions to the charge
       !
       use_tg = ( dffts%has_task_groups ) .AND. ( .NOT. (xclib_dft_is('meta') .OR. lxdm) )
       !
       incr = 2

       IF( use_tg ) THEN
          !
          v_siz = dffts%nnr_tg
          !
#if defined(__OPENMP_GPU)
          !call omp_target_alloc_f(fptr_dev=tg_psi_d, dimensions=[v_siz], omp_dev=omp_device)
          !$omp allocate allocator(omp_target_device_mem_alloc)
          allocate(tg_psi_d(v_siz))
          ALLOCATE( tg_rho( v_siz ) )
          !$omp target enter data map(alloc:tg_rho)
#else
          ALLOCATE( tg_psi_d( v_siz ) )
          ALLOCATE( tg_rho_d( v_siz ) )
          ALLOCATE( tg_rho_h( v_siz ) )
#endif
          !
          incr  = 2 *  fftx_ntgrp(dffts)
          !
       ELSE
#if defined(__OPENMP_GPU)
          !$omp target teams distribute parallel do collapse(2)
          do j=1,nspin
             do i=1,dfftp%nnr
                rho%of_r(i,j) = 0.0_DP
             enddo
          enddo
#else
          ALLOCATE( rho_d, MOLD=rho%of_r ) ! OPTIMIZE HERE, use buffers (and batched FFT)
          ALLOCATE(psic_d(dfftp%nnr))
          rho_d = 0.0_DP
#endif
       END IF
       !
       k_loop: DO ik = 1, nks
          !
          IF ( use_tg ) THEN
#if defined(__OPENMP_GPU)
             !$omp target teams distribute parallel do
             do i=1,v_siz
                tg_rho(i) = 0.0_DP
             enddo
#else
             tg_rho_d = 0.0_DP
#endif
          ENDIF
          !
          IF ( lsda ) current_spin = isk(ik)
          !
          npw = ngk(ik)
          !
          CALL start_clock_gpu( 'sum_band:buffer' )
          IF ( nks > 1 ) &
             CALL get_buffer ( evc, nwordwfc, iunwfc, ik )
          IF ( nks > 1 ) CALL using_evc(2) ! get_buffer(evc, ...) evc is updated (intent out)
          IF ( nks > 1 ) CALL using_evc_d(0) ! sync on the GPU
          !
          CALL stop_clock_gpu( 'sum_band:buffer' )
          !
          CALL start_clock_gpu( 'sum_band:init_us_2' )
          !
          IF ( nkb > 0 ) CALL init_us_2( npw, igk_k(1,ik), xk(1,ik), vkb )
          !
          CALL stop_clock_gpu( 'sum_band:init_us_2' )
          !
          ! ... here we compute the band energy: the sum of the eigenvalues
          !
          DO ibnd = ibnd_start, ibnd_end
             !
             ! ... the sum of eband and demet is the integral for
             ! ... e < ef of e n(e) which reduces for degauss=0 to the sum of
             ! ... the eigenvalues.
             !
             eband = eband + et(ibnd,ik) * wg(ibnd,ik)
             !
          END DO
          !
          DO ibnd = ibnd_start, ibnd_end, incr
             !
             IF( use_tg ) THEN
                !
#if defined(__OPENMP_GPU)
                !$omp target teams distribute parallel do
                do i=1,v_siz
                   tg_psi_d(i) = 0.0_DP
                enddo
#else
                tg_psi_d(:) = ( 0.D0, 0.D0 )
#endif
                ioff   = 0
                !
                CALL tg_get_nnr( dffts, right_nnr )
                ntgrp = fftx_ntgrp(dffts)
                !
                DO idx = 1, 2*ntgrp, 2
                   !
                   ! ... 2*ntgrp ffts at the same time
                   !
                   IF( idx + ibnd - 1 < ibnd_end ) THEN
#if defined(__OPENMP_GPU)
                      !$omp target teams distribute parallel do
                      DO j = 1, npw
                         tg_psi_d(dffts%nl (j)+ioff)=      evc(j,idx+ibnd-1)+&
                              (0.0d0,1.d0) * evc(j,idx+ibnd)
                         tg_psi_d(dffts%nlm(j)+ioff)=CONJG(evc(j,idx+ibnd-1) -&
                              (0.0d0,1.d0) * evc(j,idx+ibnd) )
                      END DO
                      !$omp end target teams distribute parallel do
#else
!$cuf kernel do(1) <<<*,*>>>
                      DO j = 1, npw
                         tg_psi_d(dffts_nl_d (j)+ioff)=     evc_d(j,idx+ibnd-1)+&
                              (0.0d0,1.d0) * evc_d(j,idx+ibnd)
                         tg_psi_d(dffts_nlm_d(j)+ioff)=CONJG(evc_d(j,idx+ibnd-1) -&
                              (0.0d0,1.d0) * evc_d(j,idx+ibnd) )
                      END DO
#endif
                   ELSE IF( idx + ibnd - 1 == ibnd_end ) THEN
#if defined(__OPENMP_GPU)
                      !$omp target teams distribute parallel do
                      DO j = 1, npw
                         tg_psi_d(dffts%nl (j)+ioff)=       evc(j,idx+ibnd-1)
                         tg_psi_d(dffts%nlm(j)+ioff)=CONJG( evc(j,idx+ibnd-1) )
                      END DO
                      !$omp end target teams distribute parallel do
#else
!$cuf kernel do(1) <<<*,*>>>
                      DO j = 1, npw
                         tg_psi_d(dffts_nl_d (j)+ioff)=       evc_d(j,idx+ibnd-1)
                         tg_psi_d(dffts_nlm_d(j)+ioff)=CONJG( evc_d(j,idx+ibnd-1) )
                      END DO
#endif
                   END IF

                   ioff = ioff + right_nnr

                END DO
                !
                !$omp dispatch
                CALL invfft (3, tg_psi_d, dffts )
                !
                ! Now the first proc of the group holds the first two bands
                ! of the 2*ntgrp bands that we are processing at the same time,
                ! the second proc. holds the third and fourth band
                ! and so on
                !
                ! Compute the proper factor for each band
                !
                idx = fftx_tgpe(dffts) + 1
                !
                ! Remember two bands are packed in a single array :
                ! proc 0 has bands ibnd   and ibnd+1
                ! proc 1 has bands ibnd+2 and ibnd+3
                ! ....
                !
                idx = 2 * idx - 1
                !
                IF( idx + ibnd - 1 < ibnd_end ) THEN
                   w1 = wg( idx + ibnd - 1, ik) / omega
                   w2 = wg( idx + ibnd    , ik) / omega
                ELSE IF( idx + ibnd - 1 == ibnd_end ) THEN
                   w1 = wg( idx + ibnd - 1, ik) / omega
                   w2 = w1
                ELSE
                   w1 = 0.0d0
                   w2 = w1
                END IF
                !
                CALL tg_get_group_nr3( dffts, right_nr3 )
                !
#if defined(__OPENMP_GPU)
                CALL get_rho_gamma_gpu(tg_rho,   dffts%nr1x*dffts%nr2x*right_nr3, w1, w2, tg_psi_d)
#else
                CALL get_rho_gamma_gpu(tg_rho_d, dffts%nr1x*dffts%nr2x*right_nr3, w1, w2, tg_psi_d)
#endif
                !
             ELSE
                !
#if defined(__OPENMP_GPU)
                !$omp target teams distribute parallel do
                do i=1,dfftp%nnr
                   psic_d(i) = (0.D0, 0.D0)
                enddo
#else
                psic_d(:) = ( 0.D0, 0.D0 )
#endif
                !
                IF ( ibnd < ibnd_end ) THEN
                   !
                   ! ... two ffts at the same time
                   !
#if defined(__OPENMP_GPU)
!$omp target teams distribute parallel do
                   DO j=1,npw
                      psic_d(dffts%nl_d(j))  = evc(j,ibnd) + &
                                              ( 0.D0, 1.D0 ) * evc(j,ibnd+1)
                      psic_d(dffts%nlm_d(j)) = CONJG( evc(j,ibnd) - &
                                              ( 0.D0, 1.D0 ) * evc(j,ibnd+1) )
                   END DO
                   !$omp end target teams distribute parallel do
#else
                   !$cuf kernel do(1)
                   DO j=1,npw
                      psic_d(dffts_nl_d(j))  = evc_d(j,ibnd) + &
                                              ( 0.D0, 1.D0 ) * evc_d(j,ibnd+1)
                      psic_d(dffts_nlm_d(j)) = CONJG( evc_d(j,ibnd) - &
                                              ( 0.D0, 1.D0 ) * evc_d(j,ibnd+1) )
                   END DO
#endif
                   !
                ELSE
                   !
#if defined(__OPENMP_GPU)
                   !$omp target teams distribute parallel do
                   DO j=1,npw
                      psic_d(dffts%nl (j)) = evc(j,ibnd)
                      psic_d(dffts%nlm(j)) = CONJG( evc(j,ibnd) )
                   END DO
                   !$omp end target teams distribute parallel do
#else
                   !$cuf kernel do(1)
                   DO j=1,npw
                      psic_d(dffts_nl_d (j))  = evc_d(j,ibnd)
                      psic_d(dffts_nlm_d(j)) = CONJG( evc_d(j,ibnd) )
                   END DO
#endif
                   !
                END IF
                !
                !$omp dispatch
                CALL invfft (2, psic_d, dffts)
                !
                w1 = wg(ibnd,ik) / omega
                !
                ! ... increment the charge density ...
                !
                IF ( ibnd < ibnd_end ) THEN
                   !
                   ! ... two ffts at the same time
                   !
                   w2 = wg(ibnd+1,ik) / omega
                   !
                ELSE
                   !
                   w2 = w1
                   !
                END IF
                !
#if defined(__OPENMP_GPU)
                CALL get_rho_gamma_gpu(rho%of_r(:,current_spin), dffts%nnr, w1, w2, psic_d(:))
#else
                CALL get_rho_gamma_gpu(rho_d   (:,current_spin), dffts%nnr, w1, w2, psic_d(:))
#endif
                !
             END IF
             !
             IF (xclib_dft_is('meta') .OR. lxdm) THEN
                CALL using_evc(0)
                DO j=1,3
                   psic(:) = ( 0.D0, 0.D0 )
                   !
                   kplusg (1:npw) = (xk(j,ik)+g(j,1:npw)) * tpiba

                   IF ( ibnd < ibnd_end ) THEN
                      ! ... two ffts at the same time
                      psic(dffts%nl (1:npw))=CMPLX(0d0, kplusg(1:npw),kind=DP) * &
                                            ( evc(1:npw,ibnd) + &
                                            ( 0.D0, 1.D0 ) * evc(1:npw,ibnd+1) )
                      psic(dffts%nlm(1:npw)) = CMPLX(0d0, -kplusg(1:npw),kind=DP) * &
                                       CONJG( evc(1:npw,ibnd) - &
                                            ( 0.D0, 1.D0 ) * evc(1:npw,ibnd+1) )
                   ELSE
                      psic(dffts%nl(1:npw)) = CMPLX(0d0, kplusg(1:npw),kind=DP) * &
                                              evc(1:npw,ibnd)
                      psic(dffts%nlm(1:npw)) = CMPLX(0d0, -kplusg(1:npw),kind=DP) * &
                                       CONJG( evc(1:npw,ibnd) )
                   END IF
                   !
                   CALL invfft (2, psic, dffts)
                   !
                   ! ... increment the kinetic energy density ...
                   !
                   DO ir = 1, dffts%nnr
                      rho%kin_r(ir,current_spin) = &
                                           rho%kin_r(ir,current_spin) + &
                                           w1 *  DBLE( psic(ir) )**2 + &
                                           w2 * AIMAG( psic(ir) )**2
                   END DO
                   !
                END DO
             END IF
             !
             !
          END DO
          !
          IF( use_tg ) THEN
#if defined(__OPENMP_GPU)
             !$omp target update from(tg_rho)
             CALL tg_reduce_rho( rho%of_r, tg_rho,   current_spin, dffts )
#else
             tg_rho_h = tg_rho_d
             CALL tg_reduce_rho( rho%of_r, tg_rho_h, current_spin, dffts )
#endif
             !
          END IF
          !
          ! ... If we have a US pseudopotential we compute here the becsum term
          !
          IF ( okvan ) CALL sum_bec_gpu ( ik, current_spin, ibnd_start,ibnd_end,this_bgrp_nbnd )
          !
       END DO k_loop
       !
       IF( .not. use_tg ) THEN
#if defined(__OPENMP_GPU)
          !$omp target update from(rho%of_r)
#else
          rho%of_r = rho_d
#endif
       END IF
       !
       ! ... with distributed <beta|psi>, sum over bands
       !
       IF ( okvan .AND. becp%comm /= mp_get_comm_null() .AND. nhm>0) THEN
          !becsum=becsum_d      not needed, since already updated in sum_bec_gpu
          CALL mp_sum( becsum, becp%comm )
#if !defined(__OPENMP_GPU)
          becsum_d=becsum
#else
          !$omp target update to(becsum)
#endif
       ENDIF
       IF ( okvan .AND. becp%comm /= mp_get_comm_null() .AND. tqr .AND. nhm>0) THEN
          !ebecsum=ebecsum_d    as above
          CALL mp_sum( ebecsum, becp%comm )
#if !defined(__OPENMP_GPU)
          ebecsum_d=ebecsum
#else
          !$omp target update to(ebecsum)
#endif
       ENDIF
       !
#if defined(__OPENMP_GPU)
       endassociate
#endif
       IF( use_tg ) THEN
#if defined(__OPENMP_GPU)
          !call omp_target_free_f(fptr_dev=tg_psi_d, omp_dev=omp_device)
          deallocate(tg_psi_d)
          !$omp target exit data map(delete:tg_rho)
          DEALLOCATE( tg_rho )
#else
          DEALLOCATE( tg_psi_d )
          DEALLOCATE( tg_rho_d )
          DEALLOCATE( tg_rho_h )
#endif
       ELSE
#if defined(__OPENMP_GPU)
          !call omp_target_free_f(fptr_dev=psic_d, omp_dev=omp_device)
          deallocate(psic)
#else
          DEALLOCATE(rho_d)
          DEALLOCATE(psic_d)
#endif
       END IF
       !
       RETURN
       !
     END SUBROUTINE sum_band_gamma_gpu
     !
     !
     !-----------------------------------------------------------------------
     SUBROUTINE sum_band_k_gpu()
       !-----------------------------------------------------------------------
       !! \(\texttt{sum_band}\) - part for k-points version
       !
#if defined(__OPENMP_GPU)
       USE omp_lib
       USE wavefunctions,      ONLY : psic_nc
#else
       USE wavefunctions_gpum, ONLY : psic_nc_d
#endif
       USE mp_bands,           ONLY : me_bgrp
       USE mp,                 ONLY : mp_sum, mp_get_comm_null
       USE control_flags,      ONLY : many_fft
       USE fft_helper_subroutines
       USE uspp_init,          ONLY : init_us_2
       !
       IMPLICIT NONE
       !
       ! ... local variables
       !
       REAL(DP) :: w1
       ! weights
       INTEGER :: npw, ipol, na, np
       !
       INTEGER  :: idx, ioff, ioff_tg, nxyp, incr, v_siz
       COMPLEX(DP), ALLOCATABLE :: tg_psi_d(:), tg_psi_nc_d(:,:)
       COMPLEX(DP), ALLOCATABLE :: psic_d(:)
#if !defined(__OPENMP_GPU)
       REAL(DP),    ALLOCATABLE :: tg_rho_d(:), tg_rho_nc_d(:,:)
       REAL(DP),    ALLOCATABLE :: tg_rho_h(:), tg_rho_nc_h(:,:)
       REAL(DP),    ALLOCATABLE :: rho_d(:,:)
       INTEGER,     POINTER     :: dffts_nl_d(:)
#else
       REAL(DP),    ALLOCATABLE :: tg_rho(:), tg_rho_nc(:,:)
       !INTEGER                  :: omp_device
#endif
       LOGICAL  :: use_tg
       INTEGER :: nnr, right_nnr, right_nr3, right_inc, ntgrp, ierr
       INTEGER :: i, j, group_size
       !
#if defined(__CUDA)
       attributes(device) :: psic_d, tg_psi_d, tg_rho_d, tg_psi_nc_d, tg_rho_nc_d
       attributes(device) :: rho_d, dffts_nl_d
       attributes(pinned) :: tg_rho_h, tg_rho_nc_h
#endif
       !
       CALL using_evc(0); CALL using_evc_d(0); CALL using_et(0)
#if !defined(__OPENMP_GPU)
       dffts_nl_d => dffts%nl_d
#else
       !omp_device = omp_get_default_device()
#endif
       !
       !
       ! ... here we sum for each k point the contribution
       ! ... of the wavefunctions to the charge
       !
       use_tg = ( dffts%has_task_groups ) .AND. ( .NOT. (xclib_dft_is('meta') .OR. lxdm) )
       !
       incr = 1
       nnr  = dffts%nnr
       !
       IF( use_tg ) THEN
          !
          v_siz = dffts%nnr_tg
          !
#if defined(__OPENMP_GPU)
          IF (noncolin) THEN
             ALLOCATE( tg_rho_nc( v_siz, nspin_mag ) )
             !$omp target enter data map(alloc:tg_rho_nc)
             !$omp allocate allocator(omp_target_device_mem_alloc)
             ALLOCATE( tg_psi_nc_d( v_siz, npol ) )
          ELSE
             ALLOCATE( tg_rho( v_siz ) )
             !$omp target enter data map(alloc:tg_rho)
             !$omp allocate allocator(omp_target_device_mem_alloc)
             ALLOCATE( tg_psi_d( v_siz ) )
          ENDIF
#else
          IF (noncolin) THEN
             ALLOCATE( tg_psi_nc_d( v_siz, npol ) )
             ALLOCATE( tg_rho_nc_d( v_siz, nspin_mag ) )
             ALLOCATE( tg_rho_nc_h( v_siz, nspin_mag ) )
          ELSE
             ALLOCATE( tg_psi_d( v_siz ) )
             ALLOCATE( tg_rho_d( v_siz ) )
             ALLOCATE( tg_rho_h( v_siz ) )
          ENDIF
#endif
          !
          incr  = fftx_ntgrp(dffts)
          !
       ELSE
          IF (noncolin .or. (xclib_dft_is('meta') .OR. lxdm)) THEN
             !$omp allocate allocator(omp_target_device_mem_alloc)
             ALLOCATE(psic_d(dffts%nnr))
             incr  = 1
          ELSE
             !$omp allocate allocator(omp_target_device_mem_alloc)
             ALLOCATE(psic_d(dffts%nnr * many_fft))
             incr  = many_fft
          END IF
          ! This is used as reduction variable on the device
#if defined(__OPENMP_GPU)
          !$omp target teams distribute parallel do collapse(2)
          do j=1,nspin
             do i=1,dfftp%nnr
                rho%of_r(i,j) = 0.0_DP
             enddo
          enddo
#else
          ALLOCATE(rho_d, MOLD=rho%of_r) ! OPTIMIZE HERE, use buffers!
          rho_d = 0.0_DP
#endif
       END IF
       !
#if defined(__OPENMP_GPU)
       associate(dffts_nl_d=>dffts%nl, psic_nc_d=>psic_nc, evc_d=>evc, &
                 igk_k_d=>igk_k, vkb_d=>vkb, rho_d=>rho%of_r)
#endif
       !
       k_loop: DO ik = 1, nks
          !
          IF( use_tg ) THEN
            IF (noncolin) THEN
#if defined(__OPENMP_GPU)
               !$omp target teams distribute parallel do collapse(2)
               do j=1,nspin_mag
                  do i=1,v_siz
                     tg_rho_nc(i,j) = 0.0_DP
                  enddo
               enddo
#else
               tg_rho_nc_d = 0.0_DP
#endif
            ELSE
#if defined(__OPENMP_GPU)
               !$omp target teams distribute parallel do
               do i=1,v_siz
                  tg_rho(i) = 0.0_DP
               enddo
#else
               tg_rho_d = 0.0_DP
#endif
            ENDIF
          ENDIF

          IF ( lsda ) current_spin = isk(ik)
          npw = ngk (ik)
          !
          CALL start_clock_gpu( 'sum_band:buffer' )
          IF ( nks > 1 ) &
             CALL get_buffer ( evc, nwordwfc, iunwfc, ik )
          IF ( nks > 1 ) CALL using_evc(2)
          IF ( nks > 1 ) CALL using_evc_d(0)  ! sync evc on GPU, OPTIMIZE (use async here)
          CALL stop_clock_gpu( 'sum_band:buffer' )
          !
          CALL start_clock_gpu( 'sum_band:init_us_2' )
          !
          IF ( nkb > 0 ) CALL init_us_2( npw, igk_k(1,ik), xk(1,ik), vkb, .true. )
          !
          CALL stop_clock_gpu( 'sum_band:init_us_2' )
          !
          ! ... here we compute the band energy: the sum of the eigenvalues
          !
          DO ibnd = ibnd_start, ibnd_end, incr
             !
             !IF( use_tg ) THEN
             DO idx = 1, incr
                IF( idx + ibnd - 1 <= ibnd_end ) eband = eband + et( idx + ibnd - 1, ik ) * wg( idx + ibnd - 1, ik )
             END DO
             !ELSE
             !   eband = eband + et( ibnd, ik ) * wg( ibnd, ik )
             !END IF
             !
             ! ... the sum of eband and demet is the integral for e < ef of
             ! ... e n(e) which reduces for degauss=0 to the sum of the
             ! ... eigenvalues
             w1 = wg(ibnd,ik) / omega
             !
             IF (noncolin) THEN
                IF( use_tg ) THEN
                   !
#if defined(__OPENMP_GPU)
                   !$omp target teams distribute parallel do collapse(2)
                   do j=1,npol
                      do i=1,v_siz
                         tg_psi_nc_d(i,j) = (0.D0, 0.D0)
                      enddo
                   enddo
#else
                   tg_psi_nc_d = ( 0.D0, 0.D0 )
#endif
                   !
                   CALL tg_get_nnr( dffts, right_nnr )
                   ntgrp = fftx_ntgrp( dffts )
                   !
                   ioff   = 0
                   !
                   DO idx = 1, ntgrp
                      !
                      ! ... ntgrp ffts at the same time
                      !
                      IF( idx + ibnd - 1 <= ibnd_end ) THEN
                         !$cuf kernel do(1)
                         !$omp target teams distribute parallel do
                         DO j = 1, npw
                            tg_psi_nc_d( dffts_nl_d(igk_k_d(j,ik) ) + ioff, 1 ) = &
                                                       evc_d( j, idx+ibnd-1 )
                            tg_psi_nc_d( dffts_nl_d(igk_k_d(j,ik) ) + ioff, 2 ) = &
                                                       evc_d( j+npwx, idx+ibnd-1 )
                         END DO
                      END IF

                      ioff = ioff + right_nnr

                   END DO
                   !
                   !$omp dispatch
                   CALL invfft (3, tg_psi_nc_d(:,1), dffts)
                   !$omp dispatch
                   CALL invfft (3, tg_psi_nc_d(:,2), dffts)
                   !
                   ! Now the first proc of the group holds the first band
                   ! of the ntgrp bands that we are processing at the same time,
                   ! the second proc. holds the second and so on
                   !
                   ! Compute the proper factor for each band
                   !
                   idx = fftx_tgpe(dffts) + 1
                   !
                   ! Remember
                   ! proc 0 has bands ibnd
                   ! proc 1 has bands ibnd+1
                   ! ....
                   !
                   IF( idx + ibnd - 1 <= ibnd_end ) THEN
                      w1 = wg( idx + ibnd - 1, ik) / omega
                   ELSE
                      w1 = 0.0d0
                   END IF
                   !
                   CALL tg_get_group_nr3( dffts, right_nr3 )
                   !
                   ! OPTIMIZE HERE : this is a sum of all densities in first spin channel
                   DO ipol=1,npol
#if defined(__OPENMP_GPU)
                      CALL get_rho_gpu(tg_rho_nc  (:,1), dffts%nr1x * dffts%nr2x* right_nr3, w1, tg_psi_nc_d(:,ipol))
#else
                      CALL get_rho_gpu(tg_rho_nc_d(:,1), dffts%nr1x * dffts%nr2x* right_nr3, w1, tg_psi_nc_d(:,ipol))
#endif
                   ENDDO
                   !
#if defined(__OPENMP_GPU)
                   IF (domag) CALL get_rho_domag_gpu(tg_rho_nc  (:,:), dffts%nr1x*dffts%nr2x*dffts%my_nr3p, w1, tg_psi_nc_d(:,:))
#else
                   IF (domag) CALL get_rho_domag_gpu(tg_rho_nc_d(:,:), dffts%nr1x*dffts%nr2x*dffts%my_nr3p, w1, tg_psi_nc_d(:,:))
#endif
                   !
                ELSE
!
!     Noncollinear case without task groups
!
#if defined(__OPENMP_GPU)
                   !$omp target teams distribute parallel do collapse(2)
                   do j=1, npol
                      do i=1, dfftp%nnr
                         psic_nc(i,j) = (0.D0,0.D0)
                      enddo
                   enddo
                   !$omp end target teams distribute parallel do
#else
                   psic_nc_d = (0.D0,0.D0)
#endif

!$cuf kernel do(1)
!$omp target teams distribute parallel do
                   DO j = 1, npw
                      psic_nc_d( dffts_nl_d(igk_k_d(j,ik) ), 1 ) = &
                                                 evc_d( j, ibnd )
                      psic_nc_d( dffts_nl_d(igk_k_d(j,ik) ), 2 ) = &
                                                 evc_d( j+npwx, ibnd )
                   END DO
                   !
                   !$omp dispatch
                   CALL invfft (2, psic_nc_d(:,1), dffts)
                   !$omp dispatch
                   CALL invfft (2, psic_nc_d(:,2), dffts)
                   !
                   ! increment the charge density ...
                   !
                   DO ipol=1,npol
                      CALL get_rho_gpu(rho_d(:,1), dffts%nnr, w1, psic_nc_d(:,ipol))
                   END DO
                   !
                   ! In this case, calculate also the three
                   ! components of the magnetization (stored in rho%of_r(ir,2-4))
                   !
                   IF (domag) THEN
                      CALL get_rho_domag_gpu(rho_d(1:,1:), dffts%nnr, w1, psic_nc_d(1:,1:))
                   ELSE
                      !$omp target teams distribute parallel do collapse(2)
                      do j=2,4
                         do i=1, dfftp%nnr
                            rho_d(i,j)=0.0_DP  ! OPTIMIZE HERE: this memset can be avoided
                         enddo
                      enddo
                      !$omp end target teams distribute parallel do
                   END IF
                   !
                END IF
                !
             ELSE
                !
                IF( use_tg ) THEN
                   !
!$cuf kernel do(1)
!$omp target teams distribute parallel do
                   DO j = 1, SIZE( tg_psi_d )
                      tg_psi_d(j) = ( 0.D0, 0.D0 )
                   END DO
                   !
                   ioff   = 0
                   !
                   CALL tg_get_nnr( dffts, right_nnr )
                   ntgrp = fftx_ntgrp( dffts )
                   !
                   DO idx = 1, ntgrp
                      !
                      ! ... ntgrp ffts at the same time
                      !
                      IF( idx + ibnd - 1 <= ibnd_end ) THEN
!$cuf kernel do(1)
!$omp target teams distribute parallel do
                         DO j = 1, npw
                            tg_psi_d( dffts_nl_d(igk_k_d(j,ik))+ioff ) = evc_d(j,idx+ibnd-1)
                         END DO
                      END IF

                      ioff = ioff + right_nnr

                   END DO
                   !
                   !$omp dispatch
                   CALL invfft (3, tg_psi_d, dffts)
                   !
                   ! Now the first proc of the group holds the first band
                   ! of the ntgrp bands that we are processing at the same time,
                   ! the second proc. holds the second and so on
                   !
                   ! Compute the proper factor for each band
                   !
                   idx = fftx_tgpe(dffts) + 1
                   !
                   ! Remember
                   ! proc 0 has bands ibnd
                   ! proc 1 has bands ibnd+1
                   ! ....
                   !
                   IF( idx + ibnd - 1 <= ibnd_end ) THEN
                      w1 = wg( idx + ibnd - 1, ik) / omega
                   ELSE
                      w1 = 0.0d0
                   END IF
                   !
                   CALL tg_get_group_nr3( dffts, right_nr3 )
                   !
#if defined(__OPENMP_GPU)
                   CALL get_rho_gpu(tg_rho  , dffts%nr1x * dffts%nr2x * right_nr3, w1, tg_psi_d)
#else
                   CALL get_rho_gpu(tg_rho_d, dffts%nr1x * dffts%nr2x * right_nr3, w1, tg_psi_d)
#endif
                   !
                ELSE IF (many_fft > 1 .and. (.not. (xclib_dft_is('meta') .OR. lxdm))) THEN
                   !
                   !!! == OPTIMIZE HERE == (setting to 0 and setting elements!)
                   group_size = MIN(many_fft, ibnd_end - (ibnd -1))
#if defined(__OPENMP_GPU)
                   !$omp target teams distribute parallel do
                   do i=1,nnr*group_size
                      psic_d(i) = (0.d0, 0.d0)
                   enddo
#else
                   psic_d(1: nnr*group_size) = (0.d0, 0.d0)
#endif
                   !
                   !$cuf kernel do(2) <<<*,*>>>
                   !$omp target teams distribute parallel do collapse(2)
                   DO i = 0, group_size-1
                      DO j = 1, npw
                         psic_d(dffts_nl_d(igk_k_d(j,ik))+i*nnr) = evc_d(j,ibnd+i)
                      END DO
                   END DO
                   !
                   !$omp dispatch
                   CALL invfft (2, psic_d, dffts, howmany=group_size)
                   !
                   ! ... increment the charge density ...
                   !
                   DO i = 0, group_size - 1
                     w1 = wg(ibnd+i,ik) / omega
                     CALL get_rho_gpu(rho_d(:,current_spin), nnr, w1, psic_d(i*nnr+1:))
                   ENDDO
                ELSE
                   !
#if defined(__OPENMP_GPU)
                   !$omp target teams distribute parallel do
                   do i = 1, dffts%nnr
                      psic_d(i) = ( 0.D0, 0.D0 )
                   enddo
#else
                   psic_d(:) = ( 0.D0, 0.D0 )
#endif
                   !
                   !$cuf kernel do(1)
                   !$omp target teams distribute parallel do
                   DO j = 1, npw
                      psic_d(dffts_nl_d(igk_k_d(j,ik))) = evc_d(j,ibnd)
                   END DO
                   !
                   !$omp dispatch
                   CALL invfft (2, psic_d, dffts)
                   !
                   ! ... increment the charge density ...
                   !
                   CALL get_rho_gpu(rho_d(:,current_spin), dffts%nnr, w1, psic_d(:))

                END IF
                !
                IF (xclib_dft_is('meta') .OR. lxdm) THEN
                   DO j=1,3
                      psic(:) = ( 0.D0, 0.D0 )
                      !
                      kplusg (1:npw) = (xk(j,ik)+g(j,igk_k(1:npw,ik))) * tpiba
                      psic(dffts%nl(igk_k(1:npw,ik)))=CMPLX(0d0,kplusg(1:npw),kind=DP) * &
                                              evc(1:npw,ibnd)
                      !
                      CALL invfft (2, psic, dffts)
                      !
                      ! ... increment the kinetic energy density ...
                      !
                      CALL get_rho(rho%kin_r(:,current_spin), dffts%nnr, w1, psic)
                   END DO
                END IF
                !
             END IF
             !
          END DO
          !
          IF( use_tg ) THEN
             !
             ! reduce the charge across task group
             !
#if defined(__OPENMP_GPU)
             IF (noncolin) then
                !$omp target update from(tg_rho_nc)
             else
                !$omp target update from(tg_rho)
             endif
             CALL tg_reduce_rho( rho%of_r, tg_rho_nc, tg_rho, current_spin, noncolin, domag, dffts )
#else
             IF (noncolin)       tg_rho_nc_h = tg_rho_nc_d
             IF (.not. noncolin) tg_rho_h    = tg_rho_d
             CALL tg_reduce_rho( rho%of_r, tg_rho_nc_h, tg_rho_h, current_spin, noncolin, domag, dffts )
#endif
             !
          END IF
          !
          ! ... If we have a US pseudopotential we compute here the becsum term
          !
          IF ( okvan ) CALL sum_bec_gpu ( ik, current_spin, ibnd_start,ibnd_end,this_bgrp_nbnd )
          !
       END DO k_loop
       !
       IF (.not. use_tg ) THEN
#if defined(__OPENMP_GPU)
          !$omp target update from(rho%of_r)
#else
          rho%of_r = rho_d
#endif
       END IF
       !
       ! ... with distributed <beta|psi>, sum over bands
       !
       IF ( okvan .AND. becp%comm /= mp_get_comm_null() .AND. nhm>0 ) THEN
          !becsum=becsum_d     not needed, since already updated in sum_bec_gpu
          CALL mp_sum( becsum, becp%comm )
#if defined(__OPENMP_GPU)
          !$omp target update to(becsum)
#else
          becsum_d=becsum
#endif
       ENDIF
       IF ( okvan .AND. becp%comm /= mp_get_comm_null() .AND. tqr .AND. nhm> 0) THEN
          !ebecsum=ebecsum_d    not needed as above
          CALL mp_sum( ebecsum, becp%comm )
#if defined(__OPENMP_GPU)
          !$omp target update to(ebecsum)
#else
          ebecsum_d=ebecsum
#endif
       ENDIF
       !
#if defined(__OPENMP_GPU)
       endassociate
       IF( use_tg ) THEN
          IF (noncolin) THEN
             DEALLOCATE( tg_psi_nc_d )
             !$omp target exit data map(delete:tg_rho_nc)
             DEALLOCATE( tg_rho_nc )
          ELSE
             DEALLOCATE( tg_psi_d )
             !$omp target exit data map(delete:tg_rho)
             DEALLOCATE( tg_rho )
          END IF
       ELSE
          DEALLOCATE(psic_d) ! OPTIMIZE HERE, use buffers!
       END IF
#else
       IF( use_tg ) THEN
          IF (noncolin) THEN
             DEALLOCATE( tg_psi_nc_d )
             DEALLOCATE( tg_rho_nc_d )
             DEALLOCATE( tg_rho_nc_h )
          ELSE
             DEALLOCATE( tg_psi_d )
             DEALLOCATE( tg_rho_d )
             DEALLOCATE( tg_rho_h )
          END IF
       ELSE
          DEALLOCATE(rho_d) ! OPTIMIZE HERE, use buffers!
          DEALLOCATE(psic_d) ! OPTIMIZE HERE, use buffers!
       END IF
#endif
       !
       RETURN
       !
     END SUBROUTINE sum_band_k_gpu
     !
     !
     SUBROUTINE get_rho_gpu(rho_loc_d, nrxxs_loc, w1_loc, psic_loc_d)

        IMPLICIT NONE

        INTEGER :: nrxxs_loc
        REAL(DP) :: rho_loc_d(:)
        REAL(DP) :: w1_loc
        COMPLEX(DP) :: psic_loc_d(:)
#if defined(__CUDA)
        attributes(device) :: rho_loc_d, psic_loc_d
#endif
        INTEGER :: ir

        !$cuf kernel do(1)
        !$omp target teams distribute parallel do
        DO ir = 1, nrxxs_loc
           !
           rho_loc_d(ir) = rho_loc_d(ir) + &
                         w1_loc * ( DBLE( psic_loc_d(ir) )**2 + &
                                   AIMAG( psic_loc_d(ir) )**2 )
        END DO
        !
     END SUBROUTINE get_rho_gpu
     !
     SUBROUTINE get_rho(rho_loc_h, nrxxs_loc, w1_loc, psic_loc_h)

        IMPLICIT NONE

        INTEGER :: nrxxs_loc
        REAL(DP) :: rho_loc_h(nrxxs_loc)
        REAL(DP) :: w1_loc
        COMPLEX(DP) :: psic_loc_h(nrxxs_loc)
        INTEGER :: ir

        DO ir = 1, nrxxs_loc
           !
           rho_loc_h(ir) = rho_loc_h(ir) + &
                         w1_loc * ( DBLE( psic_loc_h(ir) )**2 + &
                                   AIMAG( psic_loc_h(ir) )**2 )
           !
        END DO

     END SUBROUTINE get_rho

     SUBROUTINE get_rho_gamma_gpu(rho_loc_d, nrxxs_loc, w1_loc, w2_loc, psic_loc_d)

        IMPLICIT NONE

        INTEGER :: nrxxs_loc
        REAL(DP) :: rho_loc_d(nrxxs_loc)
        REAL(DP) :: w1_loc, w2_loc
        COMPLEX(DP) :: psic_loc_d(nrxxs_loc)
#if defined(__CUDA)
        attributes(device) :: rho_loc_d, psic_loc_d
#endif
        INTEGER :: ir

        !$cuf kernel do(1)
        !$omp target teams distribute parallel do
        DO ir = 1, nrxxs_loc
           !
           rho_loc_d(ir) = rho_loc_d(ir) + &
                         w1_loc * DBLE( psic_loc_d(ir) )**2 + &
                         w2_loc * AIMAG( psic_loc_d(ir) )**2
           !
        END DO

     END SUBROUTINE get_rho_gamma_gpu

     SUBROUTINE get_rho_domag_gpu(rho_loc_d, nrxxs_loc, w1_loc, psic_loc_d)

        IMPLICIT NONE

        INTEGER :: nrxxs_loc
        REAL(DP) :: rho_loc_d(:, :)
        REAL(DP) :: w1_loc
        COMPLEX(DP) :: psic_loc_d(:, :)
#if defined(__CUDA)
        attributes(device) :: rho_loc_d, psic_loc_d
#endif
        INTEGER :: ir

        !$cuf kernel do(1)
        !$omp target teams distribute parallel do
        DO ir = 1, nrxxs_loc
           !
           rho_loc_d(ir,2) = rho_loc_d(ir,2) + w1_loc*2.D0* &
                          (DBLE(psic_loc_d(ir,1))* DBLE(psic_loc_d(ir,2)) + &
                          AIMAG(psic_loc_d(ir,1))*AIMAG(psic_loc_d(ir,2)))

           rho_loc_d(ir,3) = rho_loc_d(ir,3) + w1_loc*2.D0* &
                          (DBLE(psic_loc_d(ir,1))*AIMAG(psic_loc_d(ir,2)) - &
                           DBLE(psic_loc_d(ir,2))*AIMAG(psic_loc_d(ir,1)))

           rho_loc_d(ir,4) = rho_loc_d(ir,4) + w1_loc* &
                          (DBLE(psic_loc_d(ir,1))**2+AIMAG(psic_loc_d(ir,1))**2 &
                          -DBLE(psic_loc_d(ir,2))**2-AIMAG(psic_loc_d(ir,2))**2)
           !
        END DO

     END SUBROUTINE get_rho_domag_gpu

END SUBROUTINE sum_band_gpu

!----------------------------------------------------------------------------
SUBROUTINE sum_bec_gpu ( ik, current_spin, ibnd_start, ibnd_end, this_bgrp_nbnd )
  !----------------------------------------------------------------------------
  !
  !! This routine computes the sum over bands:
  !
  !! \[ \sum_i \langle\psi_i|\beta_l\rangle w_i \langle\beta_m|\psi_i\rangle \]
  !
  !! for point "ik" and, for LSDA, spin "current_spin".
  !! Calls calbec to compute \(\text{"becp"}=\langle \beta_m|\psi_i \rangle\).
  !! Output is accumulated (unsymmetrized) into "becsum", module "uspp".
  !
  !! Routine used in sum_band (if okvan) and in compute_becsum, called by hinit1 (if okpaw).
  !
#if defined(__CUDA)
  USE cudafor
  USE cublas
#else
#define cublasZgemm zgemm
#define cublasDgemm dgemm
#endif
  USE kinds,              ONLY : DP
  USE becmod,             ONLY : becp, calbec, allocate_bec_type
  USE control_flags,      ONLY : gamma_only, tqr
  USE ions_base,          ONLY : nat, ntyp => nsp, ityp
  USE uspp,               ONLY : nkb, becsum, ebecsum, ofsbeta, vkb
#if !defined(__OPENMP_GPU)
  USE uspp,               ONLY : becsum_d, ebecsum_d, ofsbeta_d
#endif
  USE uspp_param,         ONLY : upf, nh, nhm
  USE wvfct,              ONLY : nbnd, wg, et, current_k
  USE klist,              ONLY : ngk, nkstot
  USE noncollin_module,   ONLY : noncolin, npol
  USE wavefunctions,      ONLY : evc
  USE realus,             ONLY : real_space, &
                                 invfft_orbital_gamma, calbec_rs_gamma, &
                                 invfft_orbital_k, calbec_rs_k
  USE us_exx,             ONLY : store_becxx0
  USE mp_bands,           ONLY : nbgrp,inter_bgrp_comm
  USE mp,                 ONLY : mp_sum
  USE wavefunctions_gpum, ONLY : using_evc, using_evc_d
  USE wvfct_gpum,         ONLY : using_et, using_et_d, using_wg_d
#if !defined(__OPENMP_GPU)
  USE wavefunctions_gpum, ONLY : evc_d
  USE wvfct_gpum,         ONLY : et_d, wg_d
  USE becmod_gpum,        ONLY : becp_d
#endif
  USE becmod_subs_gpum,   ONLY : calbec_gpu, using_becp_auto, using_becp_d_auto
  !
  ! Used to avoid unnecessary memcopy
  USE xc_lib,             ONLY : xclib_dft_is
#if defined(__OPENMP_GPU)
  USE omp_lib
  USE onemkl_blas_omp_offload_no_array_check
#endif
  !
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: ik, current_spin, ibnd_start, ibnd_end, this_bgrp_nbnd
  !
  COMPLEX(DP), ALLOCATABLE :: auxk1_d(:,:), auxk2_d(:,:), aux_nc_d(:,:)
  REAL(DP), ALLOCATABLE    :: auxg_d(:,:), aux_gk_d(:,:), aux_egk_d(:,:)
#if defined(__CUDA)
  attributes(DEVICE) :: auxk1_d, auxk2_d, aux_nc_d
  attributes(DEVICE) :: auxg_d, aux_gk_d, aux_egk_d
#endif
  INTEGER :: ibnd, kbnd, ibnd_loc, nbnd_loc, ibnd_begin  ! counters on bands
  INTEGER :: npw, ikb, jkb, ih, jh, ijh, na, np, is, js, nhnt
  ! counters on beta functions, atoms, atom types, spin, and auxiliary vars
  !
#if !defined(__OPENMP_GPU)
  REAL(DP), POINTER :: becp_d_r_d(:,:)
  COMPLEX(DP), POINTER :: becp_d_k_d(:,:), becp_d_nc_d(:,:,:)
#if defined(__CUDA)
  attributes(DEVICE) :: becp_d_r_d, becp_d_k_d, becp_d_nc_d
#endif
#endif
  !
  CALL using_wg_d(0)
  !
  CALL start_clock_gpu( 'sum_band:calbec' )
  npw = ngk(ik)
  IF ( .NOT. real_space ) THEN
     CALL using_evc_d(0)
     CALL using_becp_d_auto(2)
     ! calbec computes becp = <vkb_i|psi_j>
#if defined(__OPENMP_GPU)
     CALL calbec_gpu( npw, vkb, evc, becp)
#else
!$acc data present(vkb(:,:))
!$acc host_data use_device(vkb)
     CALL calbec_gpu( npw, vkb, evc_d(:,ibnd_start:ibnd_end), becp_d )
!$acc end host_data
!$acc end data
#endif
  ELSE
     CALL using_evc(0)
     CALL using_becp_auto(2)
     if (gamma_only) then
        do ibnd = ibnd_start, ibnd_end, 2
           call invfft_orbital_gamma(evc,ibnd,ibnd_end)
           call calbec_rs_gamma(ibnd,ibnd_end,becp%r)
        enddo
        call mp_sum(becp%r,inter_bgrp_comm)
     else
        current_k = ik
        becp%k = (0.d0,0.d0)
        do ibnd = ibnd_start, ibnd_end
           call invfft_orbital_k(evc,ibnd,ibnd_end)
           call calbec_rs_k(ibnd,ibnd_end)
        enddo
       call mp_sum(becp%k,inter_bgrp_comm)
     endif
  ENDIF
  CALL stop_clock_gpu( 'sum_band:calbec' )
  !
  ! In the EXX case with ultrasoft or PAW, a copy of becp will be
  ! saved in a global variable to be rotated later
  IF(xclib_dft_is('hybrid')) THEN       ! This if condition is not present in the CPU code!! Add it?
     CALL using_becp_auto(0)            ! this is very important to save useless memory copies from GPU to CPU
     CALL store_becxx0(ik, becp)
  ENDIF
  !
  CALL start_clock_gpu( 'sum_band:becsum' )
  CALL using_becp_d_auto(0)
  !
  DO np = 1, ntyp
     !
     IF ( upf(np)%tvanp ) THEN
        !
        ! allocate work space used to perform GEMM operations
        !
        IF ( gamma_only ) THEN
#if defined(__OPENMP_GPU)
           !$omp target map(from:nbnd_loc)
           nbnd_loc = becp%nbnd_loc
           !$omp end target
           !call omp_target_alloc_f(fptr_dev=auxg_d, dimensions=[nbnd_loc,nh(np)], omp_dev=omp_device)
           !$omp allocate allocator(omp_target_device_mem_alloc)
           allocate(auxg_d(nbnd_loc,nh(np)))
#else
           nbnd_loc = becp_d%nbnd_loc
           ALLOCATE( auxg_d( nbnd_loc, nh(np) ) )
#endif
        ELSE
#if defined(__OPENMP_GPU)
           !call omp_target_alloc_f(fptr_dev=auxk1_d, dimensions=[ibnd_end-ibnd_start+1,nh(np)*npol], omp_dev=omp_device)
           !call omp_target_alloc_f(fptr_dev=auxk2_d, dimensions=[ibnd_end-ibnd_start+1,nh(np)*npol], omp_dev=omp_device)
           !$omp allocate allocator(omp_target_device_mem_alloc)
#endif
           ALLOCATE( auxk1_d( ibnd_start:ibnd_end, nh(np)*npol ), &
                     auxk2_d( ibnd_start:ibnd_end, nh(np)*npol ) )
        END IF
        IF ( noncolin ) THEN
#if defined(__OPENMP_GPU)
           !call omp_target_alloc_f(fptr_dev=aux_nc_d, dimensions=[nh(np)*npol,nh(np)*npol], omp_dev=omp_device)
           !$omp allocate allocator(omp_target_device_mem_alloc)
#endif
           ALLOCATE ( aux_nc_d( nh(np)*npol,nh(np)*npol ) )
        ELSE
#if defined(__OPENMP_GPU)
           !call omp_target_alloc_f(fptr_dev=aux_gk_d, dimensions=[nh(np),nh(np)], omp_dev=omp_device)
           !$omp allocate allocator(omp_target_device_mem_alloc)
#endif
           ALLOCATE ( aux_gk_d( nh(np),nh(np) ) )
           if (tqr) THEN
#if defined(__OPENMP_GPU)
              !call omp_target_alloc_f(fptr_dev=aux_egk_d, dimensions=[nh(np),nh(np)], omp_dev=omp_device)
              !$omp allocate allocator(omp_target_device_mem_alloc)
#endif
              ALLOCATE ( aux_egk_d( nh(np),nh(np) ) )
           ENDIF
        END IF
        !
        !   In becp=<vkb_i|psi_j> terms corresponding to atom na of type nt
        !   run from index i=ofsbeta(na)+1 to i=ofsbeta(na)+nh(nt)
        !
        nhnt = nh(np)
        DO na = 1, nat
           !
           IF (ityp(na)==np) THEN
              !
              ! sum over bands: \sum_i <psi_i|beta_l><beta_m|psi_i> w_i
              ! copy into aux1, aux2 the needed data to perform a GEMM
              !
              IF ( noncolin ) THEN
                 !
#if !defined(__OPENMP_GPU)
                 becp_d_nc_d => becp_d%nc_d
#endif
                 !$cuf kernel do(2)
                 !$omp target teams distribute parallel do collapse(2)
                 DO is = 1, npol
                    DO ih = 1, nhnt
#if defined(__OPENMP_GPU)
                       ikb = ofsbeta(na) + ih
#else
                       ikb = ofsbeta_d(na) + ih
#endif
                       DO kbnd = 1, this_bgrp_nbnd
                          ibnd = ibnd_start + kbnd -1
#if defined(__OPENMP_GPU)
                          auxk1_d(ibnd,ih+(is-1)*nhnt)= becp%nc(ikb,is,kbnd)
                          auxk2_d(ibnd,ih+(is-1)*nhnt)= wg(ibnd,ik) * &
                                                        becp%nc(ikb,is,kbnd)
#else
                          auxk1_d(ibnd,ih+(is-1)*nhnt)= becp_d_nc_d(ikb,is,kbnd)
                          auxk2_d(ibnd,ih+(is-1)*nhnt)= wg_d(ibnd,ik) * &
                                                        becp_d_nc_d(ikb,is,kbnd)
#endif
                       END DO
                    END DO
                 END DO
                 !
                 !$omp target variant dispatch use_device_ptr(auxk1_d,auxk2_d,aux_nc_d)
                 CALL cublasZgemm ( 'C', 'N', npol*nhnt, npol*nhnt, this_bgrp_nbnd, &
                      (1.0_dp,0.0_dp), auxk1_d, this_bgrp_nbnd, auxk2_d, this_bgrp_nbnd, &
                      (0.0_dp,0.0_dp), aux_nc_d, npol*nhnt )
                 !$omp end target variant dispatch
                 !
              ELSE IF ( gamma_only ) THEN
                 !
#if !defined(__OPENMP_GPU)
                 becp_d_r_d => becp_d%r_d
                 ibnd_begin = becp_d%ibnd_begin
#endif
                 !$cuf kernel do(2)
                 !$omp target teams distribute parallel do collapse(2)
                 DO ih = 1, nhnt
                    DO ibnd_loc = 1, nbnd_loc
                       ibnd = (ibnd_start -1) + ibnd_loc + ibnd_begin - 1
#if defined(__OPENMP_GPU)
                       ikb = ofsbeta(na) + ih
                       auxg_d(ibnd_loc,ih) = becp%r(ikb,ibnd_loc) * wg(ibnd,ik)
#else
                       ikb = ofsbeta_d(na) + ih
                       auxg_d(ibnd_loc,ih) = becp_d_r_d(ikb,ibnd_loc) * wg_d(ibnd,ik)
#endif
                    END DO
                 END DO
#if defined(__OPENMP_GPU)
                 associate(r => becp%r)
                    !$omp target variant dispatch use_device_ptr(r,auxg_d,aux_gk_d)
                    CALL dgemm ( 'N', 'N', nhnt, nhnt, nbnd_loc, &
                         1.0_dp, r(ofsbeta(na)+1,1), nkb,    &
                         auxg_d, nbnd_loc, 0.0_dp, aux_gk_d, nhnt )
                    !$omp end target variant dispatch
                 endassociate
#else
                 CALL cublasDgemm ( 'N', 'N', nhnt, nhnt, nbnd_loc, &
                      1.0_dp, becp_d%r_d(ofsbeta(na)+1,1), nkb,    &
                      auxg_d, nbnd_loc, 0.0_dp, aux_gk_d, nhnt )
#endif
                 !
                 if (tqr) then
                   CALL using_et_d(0)
                   !$cuf kernel do(1)
                   !$omp target teams distribute parallel do
                   DO ih = 1, nhnt
#if defined(__OPENMP_GPU)
                      ikb = ofsbeta(na) + ih
#else
                      ikb = ofsbeta_d(na) + ih
#endif
                      DO ibnd_loc = 1, nbnd_loc
#if defined(__OPENMP_GPU)
                      auxg_d(ibnd_loc,ih) = et(ibnd_loc,ik) * auxg_d(ibnd_loc,ih)
#else
                      auxg_d(ibnd_loc,ih) = et_d(ibnd_loc,ik) * auxg_d(ibnd_loc,ih)
#endif
                      END DO
                   END DO

#if defined(__OPENMP_GPU)
                  associate(r => becp%r)
                     !$omp target variant dispatch use_device_ptr(r,auxg_d,aux_egk_d)
                     CALL dgemm ( 'N', 'N', nhnt, nhnt, nbnd_loc, &
                          1.0_dp, r(ofsbeta(na)+1,1), nkb,    &
                          auxg_d, nbnd_loc, 0.0_dp, aux_egk_d, nhnt )
                     !$omp end target variant dispatch
                  endassociate
#else
                  CALL cublasDgemm ( 'N', 'N', nhnt, nhnt, nbnd_loc, &
                       1.0_dp, becp_d%r_d(ofsbeta(na)+1,1), nkb,    &
                       auxg_d, nbnd_loc, 0.0_dp, aux_egk_d, nhnt )
#endif
                  end if
                 !
              ELSE
                 !
#if !defined(__OPENMP_GPU)
                 becp_d_k_d => becp_d%k_d
#endif
                 !$cuf kernel do(2) <<<*,*>>>
                 !$omp target teams distribute parallel do collapse(2)
                 DO ih = 1, nhnt
                    DO ibnd = ibnd_start, ibnd_end
#if defined(__OPENMP_GPU)
                       ikb = ofsbeta(na) + ih
                       auxk1_d(ibnd,ih) = becp%k(ikb,kbnd)
                       auxk2_d(ibnd,ih) = wg(ibnd,ik)*becp%k(ikb,kbnd)
#else
                       ikb = ofsbeta_d(na) + ih
                       auxk1_d(ibnd,ih) = becp_d_k_d(ikb,kbnd)
                       auxk2_d(ibnd,ih) = wg_d(ibnd,ik)*becp_d_k_d(ikb,kbnd)
#endif
                    END DO
                 END DO
                 !
                 ! only the real part is computed
                 !
#if defined(__OPENMP_GPU)
                 !$omp target variant dispatch use_device_ptr(auxk1_d,auxk2_d,aux_gk_d)
                 CALL dgemm ( 'C', 'N', nhnt, nhnt, 2*this_bgrp_nbnd, &
                      1.0_dp, auxk1_d, 2*this_bgrp_nbnd, auxk2_d, 2*this_bgrp_nbnd, &
                      0.0_dp, aux_gk_d, nhnt )
                 !$omp end target variant dispatch
#else
                 CALL cublasDgemm ( 'C', 'N', nhnt, nhnt, 2*this_bgrp_nbnd, &
                      1.0_dp, auxk1_d, 2*this_bgrp_nbnd, auxk2_d, 2*this_bgrp_nbnd, &
                      0.0_dp, aux_gk_d, nhnt )
#endif
                 !
                 if (tqr) then
                   CALL using_et_d(0)
                   !$cuf kernel do(2)
                   !$omp target teams distribute parallel do collapse(2)
                   DO ih = 1, nhnt
                      DO ibnd = ibnd_start, ibnd_end
#if defined(__OPENMP_GPU)
                         ikb = ofsbeta(na) + ih
                         auxk2_d(ibnd,ih) = et(ibnd,ik)*auxk2_d(ibnd,ih)
#else
                         ikb = ofsbeta_d(na) + ih
                         auxk2_d(ibnd,ih) = et_d(ibnd,ik)*auxk2_d(ibnd,ih)
#endif
                      END DO
                   END DO

                   !$omp target variant dispatch use_device_ptr(auxk1_d,auxk2_d,aux_egk_d)
                   CALL cublasDgemm ( 'C', 'N', nhnt, nhnt, 2*this_bgrp_nbnd, &
                        1.0_dp, auxk1_d, 2*this_bgrp_nbnd, auxk2_d, 2*this_bgrp_nbnd, &
                        0.0_dp, aux_egk_d, nhnt )
                   !$omp end target variant dispatch
                 end if

              END IF
              !
              ! copy output from GEMM into desired format
              !
              IF (noncolin .AND. .NOT. upf(np)%has_so) THEN
#if defined(__OPENMP_GPU)
                 CALL add_becsum_nc_gpu (na, np, aux_nc_d, becsum )
#else
                 CALL add_becsum_nc_gpu (na, np, aux_nc_d, becsum_d )
#endif
              ELSE IF (noncolin .AND. upf(np)%has_so) THEN
#if defined(__OPENMP_GPU)
                 CALL add_becsum_so_gpu (na, np, aux_nc_d, becsum )
#else
                 CALL add_becsum_so_gpu (na, np, aux_nc_d, becsum_d )
#endif
              ELSE
                 !
                 !$cuf kernel do(2) <<<*,*>>>
                 !$omp target teams distribute parallel do collapse(2)
                 DO ih = 1, nhnt
                    DO jh = 1, nhnt
                       ijh = jh + ((ih-1)*(2*nhnt-ih))/2  ! or use  ijtoh_d(ih,jh,np) ?  OPTIMIZE !!
                       !
                       ! nondiagonal terms summed and collapsed into a
                       ! single index (matrix is symmetric wrt (ih,jh))
                       !
#if defined(__OPENMP_GPU)
                       IF ( jh == ih ) THEN
                          becsum(ijh,na,current_spin) = &
                               becsum(ijh,na,current_spin) + aux_gk_d (ih,jh)
                          if (tqr) ebecsum(ijh,na,current_spin) = &
                             ebecsum(ijh,na,current_spin) + aux_egk_d (ih,jh)
                       ELSE IF ( jh > ih ) THEN
                          becsum(ijh,na,current_spin) = &
                               becsum(ijh,na,current_spin) + aux_gk_d(ih,jh)*2.0_dp
                          if (tqr) ebecsum(ijh,na,current_spin) = &
                             ebecsum(ijh,na,current_spin) + aux_egk_d(ih,jh)*2.0_dp
                       END IF
#else
                       IF ( jh == ih ) THEN
                          becsum_d(ijh,na,current_spin) = &
                               becsum_d(ijh,na,current_spin) + aux_gk_d (ih,jh)
                          if (tqr) ebecsum_d(ijh,na,current_spin) = &
                             ebecsum_d(ijh,na,current_spin) + aux_egk_d (ih,jh)
                       ELSE IF ( jh > ih ) THEN
                          becsum_d(ijh,na,current_spin) = &
                               becsum_d(ijh,na,current_spin) + aux_gk_d(ih,jh)*2.0_dp
                          if (tqr) ebecsum_d(ijh,na,current_spin) = &
                             ebecsum_d(ijh,na,current_spin) + aux_egk_d(ih,jh)*2.0_dp
                       END IF
#endif
                    END DO
                 END DO
                 !
              END IF
           END IF
           !
        END DO
        !
#if defined(__OPENMP_GPU)
        IF ( noncolin ) THEN
           !call omp_target_free_f(fptr_dev=aux_nc_d, omp_dev=omp_device)
           deallocate(aux_nc_d)
        ELSE
           !call omp_target_free_f(fptr_dev=aux_gk_d, omp_dev=omp_device)
           !if (tqr) call omp_target_free_f(fptr_dev=aux_egk_d, omp_dev=omp_device)
           deallocate(aux_gk_d)
           if (tqr) deallocate(aux_egk_d)
        END IF
        IF ( gamma_only ) THEN
           !call omp_target_free_f(fptr_dev=auxg_d, omp_dev=omp_device)
           deallocate(auxg_d)
        ELSE
           !call omp_target_free_f(fptr_dev=auxk1_d, omp_dev=omp_device)
           !call omp_target_free_f(fptr_dev=auxk2_d, omp_dev=omp_device)
           deallocate(auxk1_d)
           deallocate(auxk2_d)
        END IF
#else
        IF ( noncolin ) THEN
           DEALLOCATE ( aux_nc_d )
        ELSE
           DEALLOCATE ( aux_gk_d  )
           if (tqr) DEALLOCATE ( aux_egk_d  )
        END IF
        IF ( gamma_only ) THEN
           DEALLOCATE( auxg_d )
        ELSE
           DEALLOCATE( auxk2_d, auxk1_d )
        END IF
#endif
        !
     END IF
     !
  END DO
  !
  ! sync
  if (nhm > 0) then
#if defined(__OPENMP_GPU)
     !$omp target update from(becsum)
     if (tqr) then
        !$omp target update from(ebecsum)
     endif
#else
     becsum=becsum_d
     if (tqr) ebecsum=ebecsum_d
#endif
  endif
  !
  CALL stop_clock_gpu( 'sum_band:becsum' )
  !
END SUBROUTINE sum_bec_gpu
!
!----------------------------------------------------------------------------
SUBROUTINE add_becsum_nc_gpu ( na, np, becsum_nc_d, becsum_d )
!----------------------------------------------------------------------------
  !! This routine multiplies \(\text{becsum_nc}\) by the identity and the
  !! Pauli matrices, saves it in \(\text{becsum}\) for the calculation of
  !! augmentation charge and magnetization.
  !
#if defined(__CUDA)
  USE cudafor
#endif
  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat, ntyp => nsp, ityp
  USE uspp_param,           ONLY : nh, nhm
#if defined(__OPENMP_GPU)
  USE uspp,                 ONLY : ijtoh
#else
  USE uspp,                 ONLY : ijtoh_d
#endif
  USE lsda_mod,             ONLY : nspin
  USE noncollin_module,     ONLY : npol, nspin_mag, domag
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: na, np
  COMPLEX(DP), INTENT(IN) :: becsum_nc_d(nh(np),npol,nh(np),npol)
  REAL(DP), INTENT(INOUT) :: becsum_d(nhm*(nhm+1)/2,nat,nspin_mag)
#if defined(__CUDA)
  attributes(DEVICE) :: becsum_nc_d, becsum_d
#endif
  !
  ! ... local variables
  !
  INTEGER :: ih, jh, ijh, ipol, jpol, nhnp
  REAL(DP) :: fac
  !
  nhnp = nh(np)

  !$cuf kernel do(2) <<<*,*>>>
  !$omp target teams distribute parallel do collapse(2)
  DO ih = 1, nhnp
     DO jh = 1, nhnp
        IF ( jh >= ih ) THEN
           !ijh = jh + ((ih-1)*(2*nhnp-ih))/2  is this faster? Does it matter?
#if defined(__OPENMP_GPU)
           ijh=ijtoh(ih,jh,np)
#else
           ijh=ijtoh_d(ih,jh,np)
#endif
           IF ( ih == jh ) THEN
              fac = 1.0_dp
           ELSE
              fac = 2.0_dp
           END IF
           becsum_d(ijh,na,1)= becsum_d(ijh,na,1) + fac * &
                   DBLE( becsum_nc_d(ih,1,jh,1) + becsum_nc_d(ih,2,jh,2) )
           IF (domag) THEN
              becsum_d(ijh,na,2)= becsum_d(ijh,na,2) + fac *  &
                   DBLE( becsum_nc_d(ih,1,jh,2) + becsum_nc_d(ih,2,jh,1) )
              becsum_d(ijh,na,3)= becsum_d(ijh,na,3) + fac * DBLE( (0.d0,-1.d0)* &
                  (becsum_nc_d(ih,1,jh,2) - becsum_nc_d(ih,2,jh,1)) )
              becsum_d(ijh,na,4)= becsum_d(ijh,na,4) + fac * &
                   DBLE( becsum_nc_d(ih,1,jh,1) - becsum_nc_d(ih,2,jh,2) )
           END IF
        END IF
     END DO
  END DO

END SUBROUTINE add_becsum_nc_gpu
!
!----------------------------------------------------------------------------
SUBROUTINE add_becsum_so_gpu( na, np, becsum_nc_d, becsum_d )
  !----------------------------------------------------------------------------
  !! This routine multiplies \(\text{becsum_nc}\) by the identity and the Pauli
  !! matrices, rotates it as appropriate for the spin-orbit case, saves it in
  !! \(\text{becsum}\) for the calculation of augmentation charge and magnetization.
  !
#if defined(__CUDA)
  USE cudafor
#endif
  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat, ntyp => nsp, ityp
  USE uspp_param,           ONLY : nh, nhm
  USE noncollin_module,     ONLY : npol, nspin_mag, domag
#if defined(__OPENMP_GPU)
  USE uspp,                 ONLY : ijtoh, nhtol, nhtoj, indv
  USE upf_spinorb,          ONLY : fcoef
#else
  USE uspp,                 ONLY : ijtoh_d, nhtol_d, nhtoj_d, indv_d
  USE upf_spinorb,          ONLY : fcoef_d
#endif
  !
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: na, np
  COMPLEX(DP), INTENT(IN) :: becsum_nc_d(nh(np),npol,nh(np),npol)
  REAL(DP), INTENT(INOUT) :: becsum_d(nhm*(nhm+1)/2,nat,nspin_mag)
  !
  ! ... local variables
  !
  INTEGER :: ih, jh, lh, kh, ijh, is1, is2, nhnt
  COMPLEX(DP) :: fac

#if defined(__CUDA)
  attributes (DEVICE) :: becsum_nc_d, becsum_d
#elif defined(__OPENMP_GPU)
  associate(nhtoj_d => nhtoj, nhtol_d => nhtol, ijtoh_d => ijtoh, indv_d => indv, fcoef_d => fcoef)
#endif
  !
  nhnt = nh(np)
  !
  !$cuf kernel do(1)
  !$omp target teams distribute parallel do
  DO ih = 1, nhnt
     DO jh = 1, nhnt
        ijh=ijtoh_d(ih,jh,np)
        DO kh = 1, nhnt
           IF ( (nhtol_d(kh,np)==nhtol_d(ih,np)).AND. &
                (ABS(nhtoj_d(kh,np)-nhtoj_d(ih,np))<1.d8).AND. &
                (indv_d(kh,np)==indv_d(ih,np)) ) THEN ! same_lj(kh,ih,np)
              DO lh=1,nhnt
                 IF ( (nhtol_d(lh,np)==nhtol_d(jh,np)).AND. &
                      (ABS(nhtoj_d(lh,np)-nhtoj_d(jh,np))<1.d8).AND. &
                      (indv_d(lh,np)==indv_d(jh,np)) ) THEN   !same_lj(lh,jh,np)) THEN
                    DO is1=1,npol
                       DO is2=1,npol
                          fac=becsum_nc_d(kh,is1,lh,is2)
                          becsum_d(ijh,na,1)=becsum_d(ijh,na,1) + DBLE( fac * &
                               (fcoef_d(kh,ih,is1,1,np)*fcoef_d(jh,lh,1,is2,np) + &
                                fcoef_d(kh,ih,is1,2,np)*fcoef_d(jh,lh,2,is2,np)  ) )
                          IF (domag) THEN
                            becsum_d(ijh,na,2)=becsum_d(ijh,na,2) + DBLE( fac * &
                                (fcoef_d(kh,ih,is1,1,np)*fcoef_d(jh,lh,2,is2,np) +&
                                 fcoef_d(kh,ih,is1,2,np)*fcoef_d(jh,lh,1,is2,np)  ) )
                            becsum_d(ijh,na,3)=becsum_d(ijh,na,3) + DBLE( fac*(0.d0,-1.d0)*&
                               (fcoef_d(kh,ih,is1,1,np)*fcoef_d(jh,lh,2,is2,np) - &
                                fcoef_d(kh,ih,is1,2,np)*fcoef_d(jh,lh,1,is2,np)  ))
                            becsum_d(ijh,na,4)=becsum_d(ijh,na,4) + DBLE(fac * &
                               (fcoef_d(kh,ih,is1,1,np)*fcoef_d(jh,lh,1,is2,np) - &
                                fcoef_d(kh,ih,is1,2,np)*fcoef_d(jh,lh,2,is2,np)  ) )
                        END IF
                     END DO
                  END DO
               END IF
            END DO
         END IF
      END DO
   END DO
END DO

#if defined(__OPENMP_GPU)
  endassociate
#endif

END SUBROUTINE add_becsum_so_gpu

