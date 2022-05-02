!
! Copyright (C) 2021 Quantum ESPRESSSO Foundation
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
SUBROUTINE gen_us_dy_gpu_ ( npw, npwx, igk_d, xk, nat, tau, ityp, ntyp, &
                tpiba, omega, nr1, nr2, nr3, eigts1_d, eigts2_d, eigts3_d, &
                mill_d, g_d, u, dvkb_d )
  !----------------------------------------------------------------------
  !! Calculates the kleinman-bylander pseudopotentials with the
  !! derivative of the spherical harmonics projected on vector u
  !
  ! AF: more extensive use of GPU-resident vars possible
  !
  USE upf_kinds,       ONLY : dp
  USE upf_const,       ONLY : tpi
  USE uspp_data,       ONLY : nqx, tab, dq
#if defined(__OPENMP_GPU)
  USE omp_lib
  USE dmr
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
  REAL(DP), INTENT(IN) :: u(3)
  !! input: projection vector
  COMPLEX(DP), INTENT(OUT) :: dvkb_d(npwx, nkb)
  !! the beta function pseudopotential
  !
  ! ... local variables
  !
  INTEGER :: na, nt, nb, ih, l, lm, ikb, iig, ipol, i0, i1, i2, &
             i3, ig, nbm, iq, mil1, mil2, mil3, ikb_t,     &
             nht, ina, lmx2
  INTEGER :: nas(nat), ierr(4)
  !
  INTEGER, ALLOCATABLE :: ityp_d(:), ih_d(:), na_d(:), nas_d(:)
  !
  REAL(DP), ALLOCATABLE :: q(:), dylm(:,:), vkb0(:,:,:)
  !
  REAL(DP), POINTER :: gk_d(:,:), dylm_u_d(:,:), dylm_d(:,:,:)
  ! dylm = d Y_lm/dr_i in cartesian axes
  ! dylm_u as above projected on u
  COMPLEX(DP), ALLOCATABLE :: phase_d(:), sk_d(:,:)
#if !defined(__OPENMP_GPU)
  REAL(DP), ALLOCATABLE :: tau_d(:,:), q_d(:)
  REAL(DP), POINTER :: vkb0_d(:,:,:)
#else
  INTEGER :: omp_device
#endif
  !
  REAL(DP) :: px, ux, vx, wx, arg, u_ipol1, u_ipol2, u_ipol3, xk1, xk2, xk3
  COMPLEX(DP) :: pref
  !
#if defined(__CUDA)
  attributes(DEVICE) :: igk_d, mill_d, eigts1_d, eigts2_d, eigts3_d, g_d
  attributes(DEVICE) :: gk_d, q_d, sk_d, vkb0_d, dylm_u_d, dylm_d, &
                        ityp_d, phase_d, ih_d, na_d, tau_d, nas_d
#endif
  !
#if defined(__CUDA) || defined(__OPENMP_GPU)
#if defined(__OPENMP_GPU)
  omp_device = omp_get_default_device()
  call omp_target_init(array=dvkb_d, val=(0._DP,0._DP))
#else
  dvkb_d = (0._DP,0._DP)
#endif
  !
  IF (lmaxkb <= 0) RETURN
  !
  lmx2 = (lmaxkb+1)**2
  !
  CALL gpu_buffer%lock_buffer( dylm_u_d, (/ npw,lmx2 /), ierr(1) )
  CALL gpu_buffer%lock_buffer( gk_d, (/ 3,npw /), ierr(2) )
#if defined(__OPENMP_GPU)
  ALLOCATE( q(npw), vkb0(npw,nbetam,ntyp))
  !$omp target enter data map(alloc:q, vkb0)
#else
  CALL gpu_buffer%lock_buffer( vkb0_d, (/ npw,nbetam,ntyp /), ierr(3) )
  ALLOCATE( q_d(npw) )
#endif
  IF (ANY(ierr /= 0)) CALL upf_error( 'gen_us_dy_gpu', 'cannot allocate buffers', ABS(ierr) )
  !
  xk1 = xk(1)
  xk2 = xk(2)
  xk3 = xk(3)
  !
  !$cuf kernel do <<<*,*>>>
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
  CALL gpu_buffer%lock_buffer( dylm_d, (/npw,lmx2,3/), ierr(4) )
  DO ipol = 1, 3
#if defined(__OPENMP_GPU)
     CALL dylmr2_gpu( lmx2, npw, gk_d, q, dylm_d(:,:,ipol), ipol )
#else
     CALL dylmr2_gpu( lmx2, npw, gk_d, q_d, dylm_d(:,:,ipol), ipol )
#endif
  ENDDO
  !
  u_ipol1 = u(1) ; u_ipol2 = u(2) ; u_ipol3 = u(3)
  !
  !$cuf kernel do (2) <<<*,*>>>
  !$omp target teams distribute parallel do collapse(2)
  DO lm = 1, lmx2
    DO ig = 1, npw
      dylm_u_d(ig,lm) = u_ipol1*dylm_d(ig,lm,1) + &
                        u_ipol2*dylm_d(ig,lm,2) + &
                        u_ipol3*dylm_d(ig,lm,3)
    ENDDO
  ENDDO
  CALL gpu_buffer%release_buffer( dylm_d, ierr(4) )
  !
  !
  !$cuf kernel do (1) <<<*,*>>>
  !$omp target teams distribute parallel do
  DO ig = 1, npw
#if defined(__OPENMP_GPU)
     q  (ig) = SQRT(q  (ig)) * tpiba
#else
     q_d(ig) = SQRT(q_d(ig)) * tpiba
#endif
  ENDDO
  !
  !
  DO nt = 1, ntyp
     nbm = upf(nt)%nbeta
     !$cuf kernel do (2) <<<*,*>>>
     !$omp target teams distribute parallel do collapse(2)
     DO nb = 1, nbm
        DO ig = 1, npw
#if defined(__OPENMP_GPU)
           px = q  (ig)/dq - DBLE(INT(q  (ig)/dq))
           i0 = INT(q(ig)/dq) + 1
#else
           px = q_d(ig)/dq - DBLE(INT(q_d(ig)/dq))
           i0 = INT(q_d(ig)/dq) + 1
#endif
           ux = 1._DP - px
           vx = 2._DP - px
           wx = 3._DP - px
#if defined(__OPENMP_GPU)
           vkb0  (ig,nb,nt) = tab  (i0,nb,nt) * ux * vx * wx / 6._DP + &
                              tab  (i1,nb,nt) * px * vx * wx / 2._DP - &
                              tab  (i2,nb,nt) * px * ux * wx / 2._DP + &
                              tab  (i3,nb,nt) * px * ux * vx / 6._DP
#else
           vkb0_d(ig,nb,nt) = tab_d(i0,nb,nt) * ux * vx * wx / 6._DP + &
                              tab_d(i1,nb,nt) * px * vx * wx / 2._DP - &
                              tab_d(i2,nb,nt) * px * ux * wx / 2._DP + &
                              tab_d(i3,nb,nt) * px * ux * vx / 6._DP
#endif
       ENDDO
    ENDDO
    !
  ENDDO
#if defined(__OPENMP_GPU)
  !
  !$omp target exit data map(delete:q)
  DEALLOCATE( q )
  !
  !$omp target enter data map(to:ityp, nas, tau)
#else
  !
  DEALLOCATE( q_d )
  !
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
  DO ina = 1, nat
#if defined(__OPENMP_GPU)
     na = nas(ina)
     arg = (xk1 * tau  (1,na) + xk2 * tau  (2,na) &
          + xk3 * tau  (3,na) ) * tpi
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
  !$omp allocate allocator(omp_target_device_mem_alloc)
  ALLOCATE( sk_d(npw,nat) )
  !
  !$cuf kernel do (2) <<<*,*>>>
  !$omp target teams distribute parallel do collapse(2)
  DO ina = 1, nat
    DO ig = 1, npw
      !
      na = nas_d(ina)
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
      !
      dvkb_d(ig,ikb) = CMPLX(vkb0(ig,nb,nt)) * sk_d(ig,na) * &
#else
      nt = ityp_d(na)
      nb = indv_d(ih,nt)
      lm = nhtolm_d(ih,nt)
      l  = nhtol_d(ih,nt)
      pref = (0._DP, -1._DP)**l
      !
      dvkb_d(ig,ikb) = CMPLX(vkb0_d(ig,nb,nt)) * sk_d(ig,na) * &
#endif
                       CMPLX(dylm_u_d(ig,lm))  * pref / CMPLX(tpiba)
    ENDDO
  ENDDO
  !
  DEALLOCATE( sk_d )
  !
  IF (ikb_t /= nkb) CALL upf_error( 'gen_us_dy', 'unexpected error', 1 )
  !
  CALL gpu_buffer%release_buffer( dylm_u_d, ierr(1) )
  CALL gpu_buffer%release_buffer( gk_d, ierr(2) )
#if defined(__OPENMP_GPU)
  !$omp target exit data map(release:nas, indv, nhtol, nhtolm)
  !$omp target exit data map(delete:vkb0)
  DEALLOCATE(vkb0)
#else
  CALL gpu_buffer%release_buffer( vkb0_d, ierr(3) )
  DEALLOCATE( ih_d, na_d, nas_d )
#endif
  !
#endif
  !
  RETURN
  !
END SUBROUTINE gen_us_dy_gpu_
