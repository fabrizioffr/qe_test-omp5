!
! Copyright (C) 2021 Quantum ESPRESSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
SUBROUTINE gen_us_dj_gpu_ &
     ( npw, npwx, igk_d, xk, nat, tau, ityp, ntyp, tpiba, &
       omega, nr1, nr2, nr3, eigts1_d, eigts2_d, eigts3_d, mill_d, g_d, dvkb_d )
  !----------------------------------------------------------------------
  !! Calculates the beta function pseudopotentials with
  !! the derivative of the Bessel functions - GFPU version
  !
  ! AF: more gpu-resident variables can be used, avoiding local GPU-alloc
  !     and host2dev transfers
  !
  USE upf_kinds,       ONLY : dp
  USE upf_const,       ONLY : tpi
  USE uspp_data,       ONLY : nqx, tab, dq
#if defined(__OPENMP_GPU)
  USE omp_lib
  USE uspp,            ONLY : nkb, indv, nhtol, nhtolm
#else
  USE uspp,            ONLY : nkb, indv_d, nhtol_d, nhtolm_d
  USE uspp_data,       ONLY : tab_d
#endif
  USE uspp_param,      ONLY : upf, lmaxkb, nbetam, nh, nhm
  USE devxlib_buffers, ONLY : gpu_buffer
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: npw
  !! number ok plane waves
  INTEGER, INTENT(IN) :: npwx
  !! max number ok plane waves across k-points
  INTEGER, INTENT(IN) :: igk_d(npw)
  !! indices of plane waves k+G
  REAL(dp), INTENT(IN) :: xk(3)
  !! k-point
  INTEGER, INTENT(IN) :: nat
  !! number of atoms
  INTEGER, INTENT(IN) :: ityp(nat)
  !! index of type per atom
  INTEGER, INTENT(IN) :: ntyp
  !! number of atomic types
  REAL(DP), INTENT(IN) :: tau(3,nat)
  !! atomic positions (cc alat units)
  REAL(DP), INTENT(IN) :: tpiba
  !! rec.lattice units 2pi/a
  REAL(DP), INTENT(IN) :: omega
  !! cell volume
  INTEGER, INTENT(IN) :: nr1,nr2,nr3
  !! fft dims (dense grid)
  COMPLEX(DP), INTENT(IN) :: eigts1_d(-nr1:nr1,nat)
  !! structure factor 1
  COMPLEX(DP), INTENT(IN) :: eigts2_d(-nr2:nr2,nat)
  !! structure factor 2
  COMPLEX(DP), INTENT(IN) :: eigts3_d(-nr3:nr3,nat)
  !! structure factor 3
  INTEGER, INTENT(IN) :: mill_d(3,*)
  !! miller index map
  REAL(DP), INTENT(IN) :: g_d(3,*)
  !! g vectors (2pi/a units)
  COMPLEX(DP), INTENT(OUT) :: dvkb_d(npwx, nkb)
  !! the beta function pseudopotential
  !
  ! ... local variables
  !
  INTEGER  :: na, nt, nb, ih, l, lm, ikb, iig, i0, i1, i2, &
              i3, ig, nbm, iq, mil1, mil2, mil3,      &
              ikb_t, nht, ina, nas(nat), ierr(3)
  REAL(DP) :: px, ux, vx, wx, arg, u_ipol, xk1, xk2, xk3, qt
  COMPLEX(DP) :: pref
  INTEGER,  ALLOCATABLE :: ih_d(:), na_d(:)
  REAL(DP), ALLOCATABLE :: q(:), djl(:,:,:)
  !
  REAL(DP), POINTER :: gk_d(:,:), tau_d(:,:),  ylm_d(:,:)
  COMPLEX(DP), ALLOCATABLE :: phase_d(:), sk_d(:,:)
#if !defined(__OPENMP_GPU)
  REAL(DP), POINTER :: q_d(:), djl_d(:,:,:)
  INTEGER,  ALLOCATABLE :: ityp_d(:), nas_d(:)
#else
  INTEGER :: omp_device
#endif
  !
#if defined(__CUDA)
  attributes(DEVICE) :: igk_d, mill_d, eigts1_d, eigts2_d, eigts3_d, g_d
  attributes(DEVICE) :: gk_d, q_d, sk_d, djl_d, ylm_d, &
                        ityp_d, phase_d, ih_d, na_d, tau_d, nas_d
  attributes(DEVICE) :: dvkb_d
#endif
  !
#if defined(__CUDA) || defined(__OPENMP_GPU)
  IF (nkb == 0) RETURN
  !
  CALL gpu_buffer%lock_buffer( ylm_d, (/ npw,(lmaxkb+1)**2 /), ierr(1) )
  CALL gpu_buffer%lock_buffer( gk_d,  (/ 3,npw /), ierr(2) )
#if defined(__OPENMP_GPU)
  ALLOCATE( q(npw), djl(npw,nbetam,ntyp) )
  !$omp target enter data map(alloc:q, djl)
#else
  CALL gpu_buffer%lock_buffer( djl_d, (/ npw,nbetam,ntyp /), ierr(3) )
  ALLOCATE( q_d(npw) )
#endif
  !
  xk1 = xk(1)
  xk2 = xk(2)
  xk3 = xk(3)
  !
  !$cuf kernel do (1) <<<*,*>>>
  !$omp target teams distribute parallel do
  DO ig = 1, npw
     iig = igk_d(ig)
     gk_d(1,ig) = xk1 + g_d(1,iig)
     gk_d(2,ig) = xk2 + g_d(2,iig)
     gk_d(3,ig) = xk3 + g_d(3,iig)
#if defined(__OPENMP_GPU)
     q(ig) = gk_d(1,ig)**2 +  gk_d(2,ig)**2 + gk_d(3,ig)**2
#else
     q_d(ig) = gk_d(1,ig)**2 +  gk_d(2,ig)**2 + gk_d(3,ig)**2
#endif
  ENDDO
  !
#if defined(__OPENMP_GPU)
  CALL ylmr2_gpu( (lmaxkb+1)**2, npw, gk_d, q, ylm_d )
#else
  CALL ylmr2_gpu( (lmaxkb+1)**2, npw, gk_d, q_d, ylm_d )
#endif
  !
  DO nt = 1, ntyp
      nbm = upf(nt)%nbeta
      !$cuf kernel do (2) <<<*,*>>>
      !$omp target teams distribute parallel do collapse(2)
      DO nb = 1, nbm
         DO ig = 1, npw
#if defined(__OPENMP_GPU)
            qt = SQRT(q(ig)) * tpiba
#else
            qt = SQRT(q_d(ig)) * tpiba
#endif
            px = qt/dq - DBLE(INT(qt/dq))
            ux = 1._DP - px
            vx = 2._DP - px
            wx = 3._DP - px
            i0 = INT(qt/dq) + 1
            i1 = i0 + 1
            i2 = i0 + 2
            i3 = i0 + 3
#if defined(__OPENMP_GPU)
            djl  (ig,nb,nt) = (tab  (i0,nb,nt) * (-vx*wx-ux*wx-ux*vx)/6._DP + &
                               tab  (i1,nb,nt) * (+vx*wx-px*wx-px*vx)/2._DP - &
                               tab  (i2,nb,nt) * (+ux*wx-px*wx-px*ux)/2._DP + &
                               tab  (i3,nb,nt) * (+ux*vx-px*vx-px*ux)/6._DP)/dq
#else
            djl_d(ig,nb,nt) = (tab_d(i0,nb,nt) * (-vx*wx-ux*wx-ux*vx)/6._DP + &
                               tab_d(i1,nb,nt) * (+vx*wx-px*wx-px*vx)/2._DP - &
                               tab_d(i2,nb,nt) * (+ux*wx-px*wx-px*ux)/2._DP + &
                               tab_d(i3,nb,nt) * (+ux*vx-px*vx-px*ux)/6._DP)/dq
#endif
          ENDDO
        ENDDO
  ENDDO
  !
#if defined(__OPENMP_GPU)
  !$omp target exit data map(delete:q)
  DEALLOCATE( q )
#else
  DEALLOCATE( q_d )
#endif
  !
#if defined(__OPENMP_GPU)
  !$omp target enter data map(to:ityp, nas, tau)
#else
  ALLOCATE( ityp_d(nat), nas_d(nat) )
  ityp_d = ityp
  ALLOCATE( tau_d(3,nat) )
  tau_d  = tau
#endif
  !$omp allocate allocator(omp_target_device_mem_alloc)
  ALLOCATE( phase_d(nat) )
  !
  ina = 0
  DO nt = 1, ntyp
  DO na = 1, nat
    IF ( ityp(na) == nt ) THEN
      ina = ina + 1
      nas(ina) = na
    ENDIF
  ENDDO
  ENDDO
#if !defined(__OPENMP_GPU)
  nas_d = nas
#endif
  !
  !$cuf kernel do (1) <<<*,*>>>
  !$omp target teams distribute parallel do
  DO ina = 1, nat
#if defined(__OPENMP_GPU)
     na = nas(ina)
     arg = (xk1 * tau(1,na) + xk2 * tau(2,na) &
          + xk3 * tau(3,na) ) * tpi
#else
     na = nas_d(ina)
     arg = (xk1 * tau_d(1,na) + xk2 * tau_d(2,na) &
          + xk3 * tau_d(3,na) ) * tpi
#endif
     phase_d(na) = CMPLX( COS(arg), -SIN(arg), KIND=DP )
  ENDDO
  !
#if defined(__OPENMP_GPU)
  !$omp target exit data map(release:tau)
#else
  DEALLOCATE( tau_d )
#endif
  !
  !$omp allocate allocator (omp_target_device_mem_alloc)
  ALLOCATE( sk_d(npw,nat) )
  !
  !$cuf kernel do (2) <<<*,*>>>
  !$omp target teams distribute parallel do collapse(2)
  DO ina = 1, nat
    DO ig = 1, npw
      !
#if defined(__OPENMP_GPU)
      na = nas(ina)
#else
      na = nas_d(ina)
#endif
      iig = igk_d(ig)
      mil1 = mill_d(1,iig)
      mil2 = mill_d(2,iig)
      mil3 = mill_d(3,iig)
      sk_d(ig,na) = eigts1_d(mil1,na) * &
                    eigts2_d(mil2,na) * &
                    eigts3_d(mil3,na) * phase_d(na)
    ENDDO
  ENDDO
  !
  !
  DEALLOCATE( phase_d )
  !
  !$omp allocate allocator(omp_target_device_mem_alloc)
  ALLOCATE( ih_d(nat*nhm), na_d(nat*nhm) )
  !
  ikb_t = 0
  DO ina = 1, nat
    na = nas(ina)
    nht = nh(ityp(na))
    !$cuf kernel do (1) <<<*,*>>>
    !$omp target teams distribute parallel do
    DO ih = 1, nht
       ih_d(ikb_t+ih) = ih
       na_d(ikb_t+ih) = na
    ENDDO
    ikb_t = ikb_t + nht
  ENDDO
  !
  !$cuf kernel do (2) <<<*,*>>>
  !$omp target teams distribute parallel do collapse(2)
  DO ikb = 1, ikb_t
    DO ig = 1, npw
      ih = ih_d(ikb)
      na = na_d(ikb)
#if defined(__OPENMP_GPU)
      nt = ityp(na)
      nb = indv(ih,nt)
      lm = nhtolm(ih,nt)
      l  = nhtol(ih,nt)
      do i0=1,l-1
         pref = pref*(0._DP,-1._DP)
      enddo
      dvkb_d(ig,ikb) = CMPLX(djl  (ig,nb,nt)) * sk_d(ig,na) * &
#else
      nt = ityp_d(na)
      nb = indv_d(ih,nt)
      lm = nhtolm_d(ih,nt)
      l  = nhtol_d(ih,nt)
      pref = (0._DP,-1._DP)**l
      !
      dvkb_d(ig,ikb) = CMPLX(djl_d(ig,nb,nt)) * sk_d(ig,na) * &
#endif
                       CMPLX(ylm_d(ig,lm))  * pref
    ENDDO
  ENDDO
  !
  DEALLOCATE( sk_d )
  !
  IF (ikb_t /= nkb) CALL upf_error( 'gen_us_dj', 'unexpected error', 1 )
  !
  CALL gpu_buffer%release_buffer( ylm_d, ierr(1) )
  CALL gpu_buffer%release_buffer( gk_d, ierr(2) )
#if defined(__OPENMP_GPU)
  !$omp target exit data map(delete:djl)
  DEALLOCATE( djl )
#else
  CALL gpu_buffer%release_buffer( djl_d, ierr(3) )
#endif
  !
#if defined(__OPENMP_GPU)
  deallocate(ih_d, na_d)
  !$omp target exit data map(release:nas, indv, nhtol, nhtolm)
#else
  DEALLOCATE( ih_d, na_d, nas_d )
#endif
#endif
  !
  RETURN
  !
END SUBROUTINE gen_us_dj_gpu_
