!
! Copyright (C) 2001-2015 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
MODULE dfunct_gpum
!
CONTAINS
!
!-------------------------------------------------------------------------
SUBROUTINE newq_gpu(vr,deeq_d,skip_vltot)
  !----------------------------------------------------------------------
  !! This routine computes the integral of the perturbed potential with
  !! the Q function
  !
#if defined(__CUDA)
  USE cudafor
  USE cublas
#else
#if defined(__OPENMP_GPU)
  USE omp_lib
#endif
#define cublasZgemm Zgemm
#define cublasDGEMM Dgemm
#define cudaDGER    Dger
#endif
  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat, ntyp => nsp, ityp
  USE cell_base,            ONLY : omega, tpiba
  USE fft_base,             ONLY : dfftp
  USE fft_interfaces,       ONLY : fwfft
  USE gvect,                ONLY : g, gg, ngm, gstart, mill, eigts1, eigts2, eigts3
#if !defined(__OPENMP_GPU)
  USE gvect,                ONLY : g_d, gg_d, mill_d, eigts1_d, eigts2_d, eigts3_d
#endif
  USE lsda_mod,             ONLY : nspin
  USE scf,                  ONLY : vltot
  USE uspp_param,           ONLY : upf, lmaxq, nh, nhm
  USE control_flags,        ONLY : gamma_only
  USE wavefunctions,        ONLY : psic
#if !defined(__OPENMP_GPU)
  USE wavefunctions_gpum,   ONLY : psic_d
#endif
  USE noncollin_module,     ONLY : nspin_mag
  USE mp_bands,             ONLY : intra_bgrp_comm
  USE mp_pools,             ONLY : inter_pool_comm
  USE mp,                   ONLY : mp_sum
#if defined(__CUDA) || defined(__OPENMP_GPU)
  USE devxlib_buffers,      ONLY : buffer=>gpu_buffer
#endif
#if defined(__OPENMP_GPU)
  USE mp,                   ONLY : mp_sum_mapped
  !USE dmr,                  ONLY : omp_target_alloc_f, omp_target_free_f
#endif
  IMPLICIT NONE
  !
  !
  ! Input: potential , output: contribution to integral
  REAL(kind=dp), intent(in)  :: vr(dfftp%nnr,nspin)
  REAL(kind=dp), intent(out) :: deeq_d( nhm, nhm, nat, nspin )
#if defined(__CUDA)
  attributes(DEVICE) :: deeq_d
#endif
  LOGICAL, intent(in) :: skip_vltot !If .false. vltot is added to vr when necessary
  ! INTERNAL
  INTEGER :: ngm_s, ngm_e, ngm_l
  ! starting/ending indices, local number of G-vectors
  INTEGER :: ig, nt, ih, jh, na, is, ijh, nij, nb, nab, nhnt, ierr
  ! counters on g vectors, atom type, beta functions x 2,
  !   atoms, spin, aux, aux, beta func x2 (again)
  COMPLEX(DP), ALLOCATABLE :: vaux_d(:,:), aux_d(:,:), qgm_d(:,:)
    ! work space
  REAL(DP), ALLOCATABLE :: ylmk0_d(:,:), qmod_d(:), deeaux_d(:,:)
    ! spherical harmonics, modulus of G
  REAL(DP) :: fact
  !
#if !defined(__OPENMP_GPU)
  INTEGER, POINTER :: dfftp_nl_d(:)
  ! workaround for cuf kernel limitations
  !
#if defined(__CUDA)
  attributes(DEVICE) :: vaux_d, aux_d, qgm_d, ylmk0_d, qmod_d, deeaux_d, dfftp_nl_d
#endif
#endif
  ! variable to map index of atoms of the same type
#if !defined(__OPENMP_GPU)
  INTEGER, ALLOCATABLE :: na_to_nab_h(:)
  INTEGER, POINTER     :: na_to_nab_d(:)
#if defined(__CUDA)
  attributes(DEVICE) :: na_to_nab_d
#endif
#else
  INTEGER, ALLOCATABLE :: na_to_nab(:)
  INTEGER :: omp_device
#endif
  !
#if defined(__CUDA) || defined(__OPENMP_GPU)
#if !defined(__OPENMP_GPU)
  ALLOCATE(na_to_nab_h(nat))
  CALL buffer%lock_buffer(na_to_nab_d, nat, ierr)
#else
  ALLOCATE(na_to_nab(nat))
  !$omp target enter data map(alloc:na_to_nab)
#endif
  !
  IF ( gamma_only ) THEN
     fact = 2.0_dp
  ELSE
     fact = 1.0_dp
  ENDIF
  !
#if defined(__OPENMP_GPU)
  deeq_d(:,:,:,:) = 0.D0
#else
  !$omp target teams distribute parallel do collapse(4)
  REAL(kind=dp), intent(out) :: deeq_d( nhm, nhm, nat, nspin )
  do l=1,nspin
     do k=1,nat
        do j=1,nhm
           do i=1,nhm
              deeq_d(i,j,k,l) = 0.0D0
           enddo
        enddo
     enddo
  enddo
#endif
  !
  ! With k-point parallelization, distribute G-vectors across processors
  ! ngm_s = index of first G-vector for this processor
  ! ngm_e = index of last  G-vector for this processor
  ! ngm_l = local number of G-vectors
  !
  CALL divide (inter_pool_comm, ngm, ngm_s, ngm_e)
  ngm_l = ngm_e-ngm_s+1
  ! for the extraordinary unlikely case of more processors than G-vectors
  !
  IF ( ngm_l > 0 ) THEN
#if !defined(__OPENMP_GPU)
     ALLOCATE( vaux_d(ngm_l,nspin_mag), qmod_d(ngm_l), ylmk0_d( ngm_l, lmaxq*lmaxq ) )
#else
     !omp_device = call omp_get_default_device()
     !call omp_target_alloc_f(fptr_dev=vaux_d,  dimensions=[ngm_l,nspin_mag],    omp_dev=omp_device)
     !call omp_target_alloc_f(fptr_dev=qmod_d,  dimensions=[ngm_l],              omp_dev=omp_device)
     !call omp_target_alloc_f(fptr_dev=ylmk0_d, dimensions=[ngm_l, lmaxq*lmaxq], omp_dev=omp_device)
     !$omp allocate allocator(omp_target_device_mem_alloc)
     ALLOCATE( vaux_d(ngm_l,nspin_mag), qmod_d(ngm_l), ylmk0_d( ngm_l, lmaxq*lmaxq ) )
#endif
     !
#if !defined(__OPENMP_GPU)
     CALL ylmr2_gpu (lmaxq * lmaxq, ngm_l, g_d(1,ngm_s), gg_d(ngm_s), ylmk0_d)
     !$cuf kernel do
     DO ig = 1, ngm_l
        qmod_d (ig) = SQRT(gg_d(ngm_s+ig-1))*tpiba
     ENDDO
#else
     CALL ylmr2_gpu (lmaxq * lmaxq, ngm_l, g(1,ngm_s), gg(ngm_s), ylmk0_d)
!$omp target teams distribute parallel do
     DO ig = 1, ngm_l
        qmod_d (ig) = SQRT(gg(ngm_s+ig-1))*tpiba
     ENDDO
#endif
  END IF
  !
  ! ... fourier transform of the total effective potential
  !
#if !defined(__OPENMP_GPU)
  dfftp_nl_d => dfftp%nl_d
#endif
  DO is = 1, nspin_mag
     !
     IF ( (nspin_mag == 4 .AND. is /= 1) .or. skip_vltot ) THEN

#if !defined(__OPENMP_GPU)
        psic_d(1:dfftp%nnr) = vr(1:dfftp%nnr,is)
#else
        psic(1:dfftp%nnr) = vr(1:dfftp%nnr,is)
#endif

     ELSE
        !$omp parallel do default(shared) private(ig)
        do ig=1,dfftp%nnr
           psic(ig) = vltot(ig) + vr(ig,is)
        end do
        !$omp end parallel do
#if !defined(__OPENMP_GPU)
        psic_d(1:dfftp%nnr) = psic(1:dfftp%nnr)
#else
        !$omp target update to(psic)
#endif
     END IF
#if !defined(__OPENMP_GPU)
     CALL fwfft (1, psic_d, dfftp)
     !
     !$cuf kernel do
     do ig=1,ngm_l
        vaux_d(ig, is) = psic_d(dfftp_nl_d(ngm_s+ig-1))
     end do
#else
     CALL fwfft (1, psic, dfftp)
!$omp target teams distribute parallel do
     do ig=1,ngm_l
        vaux_d(ig, is) = psic(dfftp%nl(ngm_s+ig-1))
     end do
!$omp end target teams distribute parallel do
#endif
     !
  END DO

  DO nt = 1, ntyp
     !
     IF ( upf(nt)%tvanp ) THEN
        !
        ! count max number of atoms of type nt, create mapping table
        !
        nab = 0
        DO na = 1, nat
           IF ( ityp(na) == nt ) nab = nab + 1
           IF ( ityp(na) == nt ) THEN
#if !defined(__OPENMP_GPU)
              na_to_nab_h(na) = nab
#else
              na_to_nab(na)   = nab
#endif
           ELSE
#if !defined(__OPENMP_GPU)
              na_to_nab_h(na) = -1
#else
              na_to_nab(na)   = -1
#endif
           END IF
        END DO
        IF ( nab == 0 ) CYCLE ! No atoms for this type (?!?)
        !
#if !defined(__OPENMP_GPU)
        na_to_nab_d(1:nat) = na_to_nab_h(1:nat)
#else
        !$omp target update to(na_to_nab)
#endif
        !
        ! nij = max number of (ih,jh) pairs per atom type nt
        !
        nhnt = nh(nt)
        nij = nh(nt)*(nh(nt)+1)/2
#if !defined(__OPENMP_GPU)
        ALLOCATE ( qgm_d(ngm_l,nij) )
#else
        !omp_device = call omp_get_default_device()
        !call omp_target_alloc_f(fptr_dev=qgm_d, dimensions=[ngm_l,nij], omp_dev=omp_device)
        !$omp allocate allocator(omp_target_device_mem_alloc)
        ALLOCATE ( qgm_d(ngm_l,nij) )
#endif
        !
        ! ... Compute and store Q(G) for this atomic species
        ! ... (without structure factor)
        !
        ijh = 0
        DO ih = 1, nhnt
           DO jh = ih, nhnt
              ijh = ijh + 1
              CALL qvan2_gpu ( ngm_l, ih, jh, nt, qmod_d, qgm_d(1,ijh), ylmk0_d )
           END DO
        END DO
        !
#if !defined(__OPENMP_GPU)
        ALLOCATE ( aux_d (ngm_l, nab ), deeaux_d(nij, nab) )
#else
        !omp_device = call omp_get_default_device()
        !call omp_target_alloc_f(fptr_dev=aux_d,    dimensions=[ngm_l,nab], omp_dev=omp_device)
        !call omp_target_alloc_f(fptr_dev=deeaux_d, dimensions=[nij,nab],   omp_dev=omp_device)
        !$omp allocate allocator(omp_target_device_mem_alloc)
        ALLOCATE ( aux_d (ngm_l, nab ), deeaux_d(nij, nab) )
#endif
        !
        ! ... Compute and store V(G) times the structure factor e^(-iG*tau)
        !
        DO is = 1, nspin_mag
#if !defined(__OPENMP_GPU)
           !$cuf kernel do(2)
           DO na = 1, nat
              DO ig=1,ngm_l
                 nb = na_to_nab_d(na)
                 IF (nb > 0) &
                    aux_d(ig, nb) = vaux_d(ig,is) * CONJG ( &
                      eigts1_d(mill_d(1,ngm_s+ig-1),na) * &
                      eigts2_d(mill_d(2,ngm_s+ig-1),na) * &
                      eigts3_d(mill_d(3,ngm_s+ig-1),na) )
              END DO
           END DO
#else
           !$omp target teams distribute parallel do collapse(2)
           DO na = 1, nat
              DO ig=1,ngm_l
                 nb = na_to_nab(na)
                 IF (nb > 0) &
                    aux_d(ig, nb) = vaux_d(ig,is) * CONJG ( &
                      eigts1(mill(1,ngm_s+ig-1),na) * &
                      eigts2(mill(2,ngm_s+ig-1),na) * &
                      eigts3(mill(3,ngm_s+ig-1),na) )
              END DO
           END DO
#endif
           !
           ! ... here we compute the integral Q*V for all atoms of this kind
           !
#if defined(__OPENMP_GPU)
           !$omp target variant dispatch use_device_ptr(qgm_d,aux_d,deeaux_d)
           CALL DGEMM( 'C', 'N', nij, nab, 2*ngm_l, fact, qgm_d, 2*ngm_l, aux_d, &
                    2*ngm_l, 0.0_dp, deeaux_d, nij )
           IF ( gamma_only .AND. gstart == 2 ) &
                CALL DGER(nij, nab,-1.0_dp, qgm_d, 2*ngm_l,aux_d,2*ngm_l,deeaux_d,nij)
           !$omp end target variant dispatch
#else
           CALL cublasDGEMM( 'C', 'N', nij, nab, 2*ngm_l, fact, qgm_d, 2*ngm_l, aux_d, &
                    2*ngm_l, 0.0_dp, deeaux_d, nij )
           IF ( gamma_only .AND. gstart == 2 ) &
                CALL cudaDGER(nij, nab,-1.0_dp, qgm_d, 2*ngm_l,aux_d,2*ngm_l,deeaux_d,nij)
#endif
           !
           nhnt = nh(nt)
#if defined(__OPENMP_GPU)
           ASSOCIATE(na_to_nab_d => na_to_nab)
#endif
           !$cuf kernel do(3)
           !$omp target teams distribute parallel do collapse(3)
           DO na = 1, nat
              DO ih = 1, nhnt
                 DO jh = 1, nhnt
                    nb = na_to_nab_d(na)
                    IF (nb > 0) THEN
                       ijh = jh + ((ih-1)*(2*nhnt-ih))/2
                       IF (jh >= ih) deeq_d(ih,jh,na,is) = omega * deeaux_d(ijh,nb)
                       IF (jh > ih) deeq_d(jh,ih,na,is) = deeq_d(ih,jh,na,is)
                    END IF
                 END DO
              END DO
           END DO
#if defined(__OPENMP_GPU)
           ENDASSOCIATE
#endif
           !
        END DO
        !
#if !defined(__OPENMP_GPU)
        DEALLOCATE ( deeaux_d, aux_d, qgm_d )
#else
        !call omp_target_free_f(fptr_dev=deeaux_d, omp_dev=omp_device)
        !call omp_target_free_f(fptr_dev=aux_d,    omp_dev=omp_device)
        !call omp_target_free_f(fptr_dev=qgm_d,    omp_dev=omp_device)
        DEALLOCATE ( deeaux_d, aux_d, qgm_d )
#endif
        !
     END IF
     !
  END DO
  !
#if !defined(__OPENMP_GPU)
  DEALLOCATE( qmod_d, ylmk0_d, vaux_d )
  DEALLOCATE(na_to_nab_h); CALL buffer%release_buffer(na_to_nab_d, ierr)
#else
  !call omp_target_free_f(fptr_dev=qmod_d,  omp_dev=omp_device)
  !call omp_target_free_f(fptr_dev=ylmk0_d, omp_dev=omp_device)
  !call omp_target_free_f(fptr_dev=vaux_d,  omp_dev=omp_device)
  DEALLOCATE( qmod_d, ylmk0_d, vaux_d )
  !$omp target exit data map(delete:na_to_nab)
  DEALLOCATE(na_to_nab)
#endif
  !
  ! REPLACE THIS WITH THE NEW allgather with type or use CPU variable! OPTIMIZE HERE
  CALL mp_sum( deeq_d( :, :, :, 1:nspin_mag ), inter_pool_comm )
  CALL mp_sum( deeq_d( :, :, :, 1:nspin_mag ), intra_bgrp_comm )
#endif
  !
END SUBROUTINE newq_gpu
  !
!----------------------------------------------------------------------------
SUBROUTINE newd_gpu( )
  !----------------------------------------------------------------------------
  !! This routine computes the integral of the effective potential with
  !! the Q function and adds it to the bare ionic D term which is used
  !! to compute the non-local term in the US scheme.
  !
#if defined(__CUDA)
  use cudafor
  use cublas
#elif defined(__OPENMP_GPU)
  USE omp_lib
#endif
  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat, ntyp => nsp, ityp
  USE lsda_mod,             ONLY : nspin
  USE uspp,                 ONLY : deeq, okvan, deeq_nc
#if !defined(__OPENMP_GPU)
  USE uspp,                 ONLY : deeq_d, deeq_nc_d, dvan_d, dvan_so_d
#else
  USE uspp,                 ONLY : dvan, dvan_so
#endif
  USE uspp_param,           ONLY : upf, lmaxq, nh, nhm
  USE noncollin_module,     ONLY : noncolin, domag, nspin_mag, lspinorb
  USE uspp,                 ONLY : nhtol, nhtolm
  USE scf,                  ONLY : v
  USE realus,               ONLY : newq_r
  USE control_flags,        ONLY : tqr
  USE ldaU,                 ONLY : lda_plus_U, U_projection
  USE devxlib_buffers,      ONLY : buffer=>gpu_buffer
  !
  IMPLICIT NONE
  !
  INTEGER :: ig, nt, ih, jh, na, is, nht, nb, mb, ierr
  ! counters on g vectors, atom type, beta functions x 2,
  !   atoms, spin, aux, aux, beta func x2 (again)
  REAL(kind=dp), allocatable :: deeq_h( :,:,:,: )
#if defined(__OPENMP_GPU)
  INTEGER :: omp_default, omp_host
#endif
#if !defined(__OPENMP_GPU)
  INTEGER, POINTER :: ityp_d(:)
#if defined(__CUDA)
  attributes(DEVICE) :: ityp_d
#endif
#endif
  !
#if defined(__CUDA) || defined(__OPENMP_GPU)
  IF ( .NOT. okvan ) THEN
     !
     ! ... no ultrasoft potentials: use bare coefficients for projectors
     !
#if !defined(__OPENMP_GPU)
     CALL buffer%lock_buffer(ityp_d, nat, ierr)
     ityp_d(1:nat)=ityp(1:nat)
#else
     !$omp target data map(to:ityp)
#endif
     DO nt = 1, ntyp
        !
        nht = nh(nt)
        !
        IF ( lspinorb ) THEN
           !
#if !defined(__OPENMP_GPU)
           !$cuf kernel do(4)
           DO is =  1, nspin
              DO na = 1, nat
                 DO jh = 1, nht
                    DO ih = 1, nht
                       IF ( ityp_d(na) == nt ) deeq_nc_d(ih,jh,na,is) = dvan_so_d(ih,jh,is,nt)
                    END DO
                 END DO
              END DO
           END DO
#else
!$omp target teams distribute parallel do collapse(4)
           DO is =  1, nspin
              DO na = 1, nat
                 DO jh = 1, nht
                    DO ih = 1, nht
                       IF ( ityp(na) == nt ) deeq_nc(ih,jh,na,is) = dvan_so(ih,jh,is,nt)
                    END DO
                 END DO
              END DO
           END DO
!$omp end target teams distribute parallel do
#endif
           !
        ELSE IF ( noncolin ) THEN
           !
#if !defined(__OPENMP_GPU)
           !$cuf kernel do(3)
           DO na = 1, nat
              DO jh = 1, nht
                 DO ih = 1, nht
                    IF ( ityp_d(na) == nt ) THEN
                       deeq_nc_d(ih,jh,na,1) = dvan_d(ih,jh,nt)
                       deeq_nc_d(ih,jh,na,2) = ( 0.D0, 0.D0 )
                       deeq_nc_d(ih,jh,na,3) = ( 0.D0, 0.D0 )
                       deeq_nc_d(ih,jh,na,4) = dvan_d(ih,jh,nt)
                    END IF
                 END DO
              END DO
           END DO
#else
!$omp target teams distribute parallel do collapse(3)
           DO na = 1, nat
              DO jh = 1, nht
                 DO ih = 1, nht
                    IF ( ityp(na) == nt ) THEN
                       deeq_nc(ih,jh,na,1) = dvan(ih,jh,nt)
                       deeq_nc(ih,jh,na,2) = ( 0.D0, 0.D0 )
                       deeq_nc(ih,jh,na,3) = ( 0.D0, 0.D0 )
                       deeq_nc(ih,jh,na,4) = dvan(ih,jh,nt)
                    END IF
                 END DO
              END DO
           END DO
!$omp end target teams distribute parallel do
#endif
           !
        ELSE
           !
           if ( nht > 0 ) THEN
#if !defined(__OPENMP_GPU)
              !$cuf kernel do(4)
              DO is = 1, nspin
                 DO na = 1, nat
                    DO jh = 1, nht
                       DO ih = 1, nht
                          !
                          IF ( ityp_d(na) == nt ) deeq_d(ih,jh,na,is) = dvan_d(ih,jh,nt)
                          !
                       END DO
                    END DO
                 END DO
              END DO
              !
#else
              !$omp target teams distribute parallel do collapse(4)
              DO is = 1, nspin
                 DO na = 1, nat
                    DO jh = 1, nht
                       DO ih = 1, nht
                          !
                          IF ( ityp(na) == nt ) deeq(ih,jh,na,is) = dvan(ih,jh,nt)
                          !
                       END DO
                    END DO
                 END DO
              END DO
#endif
           end if
           !
        END IF
        !
     END DO
     !
     ! ... early return
     !
#if !defined(__OPENMP_GPU)
     CALL buffer%release_buffer(ityp_d, ierr)
#else
     !$omp end target data
#endif
     !
     ! ... sync with CPU
     if (noncolin) then
#if defined(__OPENMP_GPU)
        !$omp target update from(deeq_nc)
#else
        deeq_nc=deeq_nc_d
#endif
     else
#if defined(__OPENMP_GPU)
        !$omp target update from(deeq)
#else
        deeq=deeq_d
#endif
     endif
     !
     RETURN
     !
  END IF
  !
  CALL start_clock_gpu( 'newd' )
  allocate(deeq_h( nhm, nhm, nat, nspin ))
  !
  ! move atom type info to GPU
#if !defined(__OPENMP_GPU)
     CALL buffer%lock_buffer(ityp_d, nat, ierr)
     ityp_d(1:nat)=ityp(1:nat)
#else
     !$omp target data map(to:ityp)
#endif
  !
  IF (tqr) THEN
     CALL newq_r(v%of_r,deeq,.false.)
#if defined(__OPENMP_GPU)
     !$omp target update to(deeq)
#else
     deeq_d=deeq
#endif
  ELSE
#if !defined(__OPENMP_GPU)
     CALL newq_gpu(v%of_r,deeq_d,.false.)
#else
     CALL newq_gpu(v%of_r,deeq,.false.)
#endif
  END IF
  !
#if !defined(__OPENMP_GPU)
  IF (noncolin) call add_paw_to_deeq_gpu(deeq_d)
#else
  IF (noncolin) call add_paw_to_deeq_gpu(deeq)
#endif
  !
  types : &
  DO nt = 1, ntyp
     !
     if_noncolin:&
     IF ( noncolin ) THEN
        !
        IF (upf(nt)%has_so) THEN
           !
           CALL newd_so_gpu(nt)
           !
        ELSE
           !
           CALL newd_nc_gpu(nt)
           !
        END IF
        !
     ELSE if_noncolin
        !
        nht = nh(nt)
        !$cuf kernel do(4)
        !$omp target teams distribute parallel do collapse(4)
        DO is = 1, nspin
           DO na = 1, nat
              DO ih = 1, nht
                 DO jh = 1, nht
#if !defined(__OPENMP_GPU)
                    IF ( ityp_d(na) == nt ) THEN
                       deeq_d(ih,jh,na,is) = deeq_d(ih,jh,na,is) + dvan_d(ih,jh,nt)
                    END IF
#else
                    IF ( ityp(na) == nt ) THEN
                       deeq(ih,jh,na,is) = deeq(ih,jh,na,is) + dvan(ih,jh,nt)
                    END IF
#endif
                 END DO
              END DO
           END DO
        END DO
        !
     END IF if_noncolin
     !
  END DO types
  !
#if !defined(__OPENMP_GPU)
  IF (.NOT.noncolin) CALL add_paw_to_deeq_gpu(deeq_d)
#else
  IF (.NOT.noncolin) CALL add_paw_to_deeq_gpu(deeq)
#endif
  !
#if !defined(__OPENMP_GPU)
  IF (lda_plus_U .AND. (U_projection == 'pseudo')) CALL add_vhub_to_deeq_gpu(deeq_d)
#else
  IF (lda_plus_U .AND. (U_projection == 'pseudo')) CALL add_vhub_to_deeq_gpu(deeq)
#endif
  !
#if !defined(__OPENMP_GPU)
  CALL buffer%release_buffer(ityp_d, ierr)
#else
  !$omp end target data
#endif
  CALL stop_clock_gpu( 'newd' )
  !
  if (noncolin) then
#if defined(__OPENMP_GPU)
     !$omp target update from(deeq_nc)
#else
     deeq_nc=deeq_nc_d
#endif
  else
#if defined(__OPENMP_GPU)
     !$omp target update from(deeq)
#else
     deeq=deeq_d
#endif
  endif
  !
  RETURN
  !
  CONTAINS
    !
    !------------------------------------------------------------------------
    SUBROUTINE newd_so_gpu(nt)
      !------------------------------------------------------------------------
      !
#if defined(__OPENMP_GPU)
      USE upf_spinorb,   ONLY : fcoef
#else
      USE upf_spinorb,   ONLY : fcoef_d
#endif
      USE ions_base,     ONLY : nat
      !
      IMPLICIT NONE
      !
      INTEGER :: nt

      INTEGER :: ijs, is1, is2, kh, lh, nhnt, ih, jh, na
      !
      nhnt = nh(nt)
      ijs = 0
      !
      DO is1 = 1, 2
         !
         DO is2 =1, 2
            !
            ijs = ijs + 1
            !
            IF (domag) THEN
#if !defined(__OPENMP_GPU)
               !$cuf kernel do(3)
               DO na = 1, nat
                  !
                  DO ih = 1, nhnt
                     !
                     DO jh = 1, nhnt
                        !
                        IF ( ityp_d(na) == nt ) THEN
                           !
                           deeq_nc_d(ih,jh,na,ijs) = dvan_so_d(ih,jh,ijs,nt)
                           !
                           DO kh = 1, nhnt
                              !
                              DO lh = 1, nhnt
                                 !
                                 deeq_nc_d(ih,jh,na,ijs) = deeq_nc_d(ih,jh,na,ijs) +   &
                                      deeq_d (kh,lh,na,1)*            &
                                   (fcoef_d(ih,kh,is1,1,nt)*fcoef_d(lh,jh,1,is2,nt)  + &
                                   fcoef_d(ih,kh,is1,2,nt)*fcoef_d(lh,jh,2,is2,nt)) + &
                                   deeq_d (kh,lh,na,2)*            &
                                   (fcoef_d(ih,kh,is1,1,nt)*fcoef_d(lh,jh,2,is2,nt)  + &
                                   fcoef_d(ih,kh,is1,2,nt)*fcoef_d(lh,jh,1,is2,nt)) + &
                                   (0.D0,-1.D0)*deeq_d (kh,lh,na,3)*            &
                                   (fcoef_d(ih,kh,is1,1,nt)*fcoef_d(lh,jh,2,is2,nt)  - &
                                   fcoef_d(ih,kh,is1,2,nt)*fcoef_d(lh,jh,1,is2,nt)) + &
                                   deeq_d (kh,lh,na,4)*            &
                                   (fcoef_d(ih,kh,is1,1,nt)*fcoef_d(lh,jh,1,is2,nt)  - &
                                   fcoef_d(ih,kh,is1,2,nt)*fcoef_d(lh,jh,2,is2,nt))
                                 !
                              END DO
                              !
                           END DO
                           !
                        END IF
                        !
                     END DO
                  END DO
                  !
               END DO
#else
               !$omp target teams distribute parallel do collapse(3)
               DO na = 1, nat
                  DO ih = 1, nhnt
                     DO jh = 1, nhnt
                        IF ( ityp(na) == nt ) THEN
                           deeq_nc(ih,jh,na,ijs) = dvan_so(ih,jh,ijs,nt)
                           DO kh = 1, nhnt
                              DO lh = 1, nhnt
                                 deeq_nc(ih,jh,na,ijs) = deeq_nc(ih,jh,na,ijs) +   &
                                      deeq (kh,lh,na,1)*            &
                                   (fcoef(ih,kh,is1,1,nt)*fcoef(lh,jh,1,is2,nt)  + &
                                   fcoef(ih,kh,is1,2,nt)*fcoef(lh,jh,2,is2,nt)) + &
                                   deeq (kh,lh,na,2)*            &
                                   (fcoef(ih,kh,is1,1,nt)*fcoef(lh,jh,2,is2,nt)  + &
                                   fcoef(ih,kh,is1,2,nt)*fcoef(lh,jh,1,is2,nt)) + &
                                   (0.D0,-1.D0)*deeq (kh,lh,na,3)*            &
                                   (fcoef(ih,kh,is1,1,nt)*fcoef(lh,jh,2,is2,nt)  - &
                                   fcoef(ih,kh,is1,2,nt)*fcoef(lh,jh,1,is2,nt)) + &
                                   deeq (kh,lh,na,4)*            &
                                   (fcoef(ih,kh,is1,1,nt)*fcoef(lh,jh,1,is2,nt)  - &
                                   fcoef(ih,kh,is1,2,nt)*fcoef(lh,jh,2,is2,nt))
                              END DO
                           END DO
                        END IF
                     END DO
                  END DO
               END DO
               !$omp end target teams distribute parallel do
#endif
               !
            ELSE
               !
#if !defined(__OPENMP_GPU)
               !$cuf kernel do(3) <<<*,*>>>
               DO na = 1, nat
                  !
                  DO ih = 1, nhnt
                     !
                     DO jh = 1, nhnt
                        !
                        IF ( ityp_d(na) == nt ) THEN
                           !
                           deeq_nc_d(ih,jh,na,ijs) = dvan_so_d(ih,jh,ijs,nt)
                           !
                           DO kh = 1, nhnt
                              !
                              DO lh = 1, nhnt
                                 !
                                 deeq_nc_d(ih,jh,na,ijs) = deeq_nc_d(ih,jh,na,ijs) +   &
                                      deeq_d (kh,lh,na,1)*            &
                                   (fcoef_d(ih,kh,is1,1,nt)*fcoef_d(lh,jh,1,is2,nt)  + &
                                   fcoef_d(ih,kh,is1,2,nt)*fcoef_d(lh,jh,2,is2,nt) )
                                 !
                              END DO
                              !
                           END DO
                           !
                        END IF
                        !
                     END DO
                     !
                  END DO
                  !
               END DO
#else
               !$omp target teams distribute parallel do collapse(3)
               DO na = 1, nat
                  DO ih = 1, nhnt
                     DO jh = 1, nhnt
                        IF ( ityp(na) == nt ) THEN
                           deeq_nc(ih,jh,na,ijs) = dvan_so(ih,jh,ijs,nt)
                           DO kh = 1, nhnt
                              DO lh = 1, nhnt
                                 deeq_nc(ih,jh,na,ijs) = deeq_nc(ih,jh,na,ijs) +   &
                                      deeq (kh,lh,na,1)*            &
                                   (fcoef(ih,kh,is1,1,nt)*fcoef(lh,jh,1,is2,nt)  + &
                                   fcoef(ih,kh,is1,2,nt)*fcoef(lh,jh,2,is2,nt) )
                              END DO
                           END DO
                        END IF
                     END DO
                  END DO
               END DO
               !$omp end target teams distribute parallel do
#endif
               !
            END IF
            !
         END DO
         !
      END DO
      !
    RETURN
      !
    END SUBROUTINE newd_so_gpu
    !
    !------------------------------------------------------------------------
    SUBROUTINE newd_nc_gpu(nt)
      !------------------------------------------------------------------------
      !
      USE ions_base,     ONLY : nat
      !
      IMPLICIT NONE
      !
      INTEGER, INTENT(IN) :: nt
      INTEGER :: nhnt, na, ih, jh
      !
      nhnt = nh(nt)
      !
#if !defined(__OPENMP_GPU)
      !$cuf kernel do(3)
      DO na = 1, nat
         !
         DO ih = 1, nhnt
            !
            DO jh = 1, nhnt
               !
               IF ( ityp_d(na) == nt ) THEN
                  !
                  IF (lspinorb) THEN
                     deeq_nc_d(ih,jh,na,1) = dvan_so_d(ih,jh,1,nt) + &
                                           deeq_d(ih,jh,na,1) + deeq_d(ih,jh,na,4)
                     !
                     deeq_nc_d(ih,jh,na,4) = dvan_so_d(ih,jh,4,nt) + &
                                           deeq_d(ih,jh,na,1) - deeq_d(ih,jh,na,4)
                     !
                  ELSE
                     deeq_nc_d(ih,jh,na,1) = dvan_d(ih,jh,nt) + &
                                           deeq_d(ih,jh,na,1) + deeq_d(ih,jh,na,4)
                     !
                     deeq_nc_d(ih,jh,na,4) = dvan_d(ih,jh,nt) + &
                                           deeq_d(ih,jh,na,1) - deeq_d(ih,jh,na,4)
                     !
                  END IF
                  deeq_nc_d(ih,jh,na,2) = deeq_d(ih,jh,na,2) - &
                                        ( 0.D0, 1.D0 ) * deeq_d(ih,jh,na,3)
                  !
                  deeq_nc_d(ih,jh,na,3) = deeq_d(ih,jh,na,2) + &
                                        ( 0.D0, 1.D0 ) * deeq_d(ih,jh,na,3)
                  !
               END IF
               !
            END DO
            !
         END DO
         !
      END DO
#else
      !$omp target teams distribute parallel do collapse(3)
      DO na = 1, nat
         DO ih = 1, nhnt
            DO jh = 1, nhnt
               IF ( ityp(na) == nt ) THEN
                  IF (lspinorb) THEN
                     deeq_nc(ih,jh,na,1) = dvan_so(ih,jh,1,nt) + &
                                           deeq(ih,jh,na,1) + deeq(ih,jh,na,4)
                     deeq_nc(ih,jh,na,4) = dvan_so(ih,jh,4,nt) + &
                                           deeq(ih,jh,na,1) - deeq(ih,jh,na,4)
                  ELSE
                     deeq_nc(ih,jh,na,1) = dvan(ih,jh,nt) + &
                                           deeq(ih,jh,na,1) + deeq(ih,jh,na,4)
                     deeq_nc(ih,jh,na,4) = dvan(ih,jh,nt) + &
                                           deeq(ih,jh,na,1) - deeq(ih,jh,na,4)
                  END IF
                  deeq_nc(ih,jh,na,2) = deeq(ih,jh,na,2) - &
                                        ( 0.D0, 1.D0 ) * deeq(ih,jh,na,3)
                  deeq_nc(ih,jh,na,3) = deeq(ih,jh,na,2) + &
                                        ( 0.D0, 1.D0 ) * deeq(ih,jh,na,3)
               END IF
            END DO
         END DO
      END DO
      !$omp end target teams distribute parallel do
#endif
      !
    RETURN
    END SUBROUTINE newd_nc_gpu
#endif
    !
END SUBROUTINE newd_gpu

END MODULE dfunct_gpum
