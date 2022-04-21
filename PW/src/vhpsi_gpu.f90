!
! Copyright (C) 2001-2020 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#if !defined(__CUDA)
#define cublasDgemm dgemm
#define cublasZgemm zgemm
#endif
!-----------------------------------------------------------------------
SUBROUTINE vhpsi_gpu( ldap, np, mps, psip_d, hpsi_d )
  !-----------------------------------------------------------------------
  !! This routine computes the Hubbard potential applied to the electronic
  !! structure of the current k-point. The result is added to hpsi.
  !
  USE kinds,         ONLY : DP
  USE becmod,        ONLY : bec_type, calbec, allocate_bec_type, &
                            deallocate_bec_type
  USE ldaU,          ONLY : Hubbard_lmax, Hubbard_l, is_Hubbard,   &
                            nwfcU, wfcU, offsetU, lda_plus_u_kind, &
                            is_hubbard_back, Hubbard_l_back, offsetU_back, &
                            backall, offsetU_back1
  USE lsda_mod,      ONLY : current_spin
  USE scf,           ONLY : v
  USE ions_base,     ONLY : nat, ntyp => nsp, ityp
  USE control_flags, ONLY : gamma_only
  USE mp,            ONLY : mp_sum
  !
#if defined(__OPENMP_GPU)
  USE becmod,           ONLY : bec_type
#else
  USE becmod_gpum,      ONLY : bec_type_d
#endif
  USE becmod_subs_gpum, ONLY : allocate_bec_type_gpu, deallocate_bec_type_gpu, &
                               calbec_gpu
  USE devxlib_memcpy,   ONLY : dev_memcpy_h2d => devxlib_memcpy_h2d
  USE devxlib_memset,   ONLY : dev_memset => devxlib_memory_set
  !
#if defined(__CUDA)
  USE cudafor
  USE cublas
#elif defined(__OPENMP_GPU)
  USE onemkl_blas_omp_offload_no_array_check
#endif
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: ldap
  !! leading dimension of arrays psip, hpsi
  INTEGER, INTENT(IN) :: np
  !! true dimension of psip, hpsi
  INTEGER, INTENT(IN) :: mps
  !! number of states psip
  COMPLEX(DP), INTENT(IN) :: psip_d(ldap,mps)
  !! the wavefunction
  COMPLEX(DP), INTENT(INOUT) :: hpsi_d(ldap,mps)
  !! Hamiltonian dot psi
  !
  ! ... local variables
  !
  INTEGER :: dimwf1
#if defined(__OPENMP_GPU)
  TYPE(bec_type), TARGET :: proj_d
#else
  TYPE(bec_type_d), TARGET :: proj_d
#endif
  COMPLEX(DP), ALLOCATABLE :: wfcU_d(:,:)
  !
#if defined(__CUDA)
  attributes(DEVICE) :: wfcU_d, psip_d, hpsi_d
#endif
  !
  CALL start_clock_gpu( 'vhpsi' )
  !
  ! Offset of atomic wavefunctions initialized in setup and stored in offsetU
  !
  ! Allocate the array proj
  CALL allocate_bec_type_gpu( nwfcU, mps, proj_d )
  !
  dimwf1 = SIZE(wfcU(:,1))
#if defined(__OPENMP_GPU)
  !$omp target enter data map(to:wfcU)
#else
  ALLOCATE( wfcU_d(dimwf1,nwfcU) )
  CALL dev_memcpy(wfcU_d, wfcU)
#endif
  !
  ! proj = <wfcU|psip>
#if defined(__OPENMP_GPU)
  CALL calbec_gpu( np, wfcU, psip_d, proj_d )
#else
  CALL calbec_gpu( np, wfcU_d, psip_d, proj_d )
#endif
  !
  IF ( lda_plus_u_kind.EQ.0 .OR. lda_plus_u_kind.EQ.1 ) THEN
     CALL vhpsi_U_gpu()  ! DFT+U
  ELSEIF ( lda_plus_u_kind.EQ.2 ) THEN
     CALL errore('vhpsi', 'DFT+U+V case not implemented for GPU', 1 )
  ENDIF
  !
  CALL deallocate_bec_type_gpu( proj_d )
  !
#if defined(__OPENMP_GPU)
  !$omp target exit data map(delete:wfcU)
#else
  DEALLOCATE( wfcU_d )
#endif
  !
  CALL stop_clock_gpu( 'vhpsi' )
  !
  RETURN
  !
CONTAINS
  !
SUBROUTINE vhpsi_U_gpu()
  !
  ! This routine applies the Hubbard potential with U_I
  ! to the KS wave functions.
  !
  USE ldaU, ONLY : ldim_back, ldmx_b, Hubbard_l1_back
#if defined(__OPENMP_GPU)
  USE omp_lib
  !USE dmr,  ONLY : omp_target_alloc_f, omp_target_free_f
#endif
  !
  IMPLICIT NONE
  !
  INTEGER :: na, nt, ldim, ldim0, ldimax, ldimaxt
  !
  REAL(DP),    ALLOCATABLE :: rtemp_d(:,:), vns_d(:,:,:),  vnsb_d(:,:,:)
  COMPLEX(DP), ALLOCATABLE :: ctemp_d(:,:), vaux_d(:,:,:), vauxb_d(:,:,:)
  !
#if defined(__CUDA)
  attributes(DEVICE) :: vns_d, vnsb_d, rtemp_d, ctemp_d, vaux_d, vauxb_d
#endif

#if defined(__OPENMP_GPU)
  INTEGER :: i, j, k!, omp_device

  !omp_device = omp_get_Default_device()
#endif
  !
  ldimax = 2*Hubbard_lmax+1
  ldimaxt = MAX(ldimax, ldmx_b)
  !
  IF ( ANY(is_hubbard(:)) .OR. ANY(is_hubbard_back(:)) ) THEN
    IF (gamma_only) THEN
#if defined(__OPENMP_GPU)
      !call omp_target_alloc_f(fptr_dev=rtemp_d, dimensions[ldimaxt,mps], omp_dev=omp_device)
      !$omp allocate allocator(omp_target_device_mem_alloc)
      ALLOCATE( rtemp_d(ldimaxt,mps) )
#else
      ALLOCATE( rtemp_d(ldimaxt,mps) )
#endif
      IF (ANY(is_hubbard(:))) THEN
        !$omp allocate allocator(omp_target_device_mem_alloc)
        ALLOCATE( vns_d(ldimax,ldimax,nat) )
#if defined(__OPENMP_GPU)
        !$omp target teams distribute parallel do collapse(3) map(to:v%ns)
        do k=1,ldimax
           do j=1,ldimax
              do i=1,nat
                 vns_d(i,j,k) = v%ns(i,j,current_spin,k)
              enddo
           enddo
        enddo
#else
        vns_d = v%ns(:,:,current_spin,:)
#endif
      ENDIF
      IF (ANY(is_hubbard_back(:))) THEN
        !$omp allocate allocator(omp_target_device_mem_alloc)
        ALLOCATE( vnsb_d(ldmx_b,ldmx_b,nat) )
#if defined(__OPENMP_GPU)
        !$omp target teams distribute parallel do collapse(3) map(to:v%nsb)
        do k=1,nat
           do j=1,ldmx_b
              do i=1,ldmx_b
                 vnsb_d(i,j,k) = v%nsb(i,j,current_spin,k)
              enddo
           enddo
        enddo
#else
        vnsb_d = v%nsb(:,:,current_spin,:)
#endif
      ENDIF
    ELSE
#if defined(__OPENMP_GPU)
      !call omp_target_alloc_f(fptr_dev=ctemp_d, dimensions[ldimaxt,mps], omp_dev=omp_device)
      !$omp allocate allocator(omp_target_device_mem_alloc)
      ALLOCATE( ctemp_d(ldimaxt,mps) )
#else
      ALLOCATE( ctemp_d(ldimaxt,mps) )
#endif
      IF (ANY(is_hubbard(:))) THEN
#if defined(__OPENMP_GPU)
        !call omp_target_alloc_f(fptr_dev=vaux_d, dimensions=[ldimax,ldimax,nat], omp_dev=omp_device)
        !$omp allocate allocator(omp_target_device_mem_alloc)
        allocate(vaux_d(ldimax,ldimax,nat))
        do k=1, nat
           !$omp target teams distribute parallel do collapse(2) map(to:v%ns)
           do j=1, ldimax
              do i=1, ldimax
                 vaux_d(i,j,k) = CMPLX(v%ns(i,j,current_spin,k))
              enddo
           enddo
           !$omp end target teams distribute parallel do
        enddo
#else
        ALLOCATE( vaux_d(ldimax,ldimax,nat) )
        vaux_d = CMPLX(v%ns(:,:,current_spin,:))
#endif
      ENDIF
      IF (ANY(is_hubbard_back(:))) THEN
#if defined(__OPENMP_GPU)
        !call omp_target_alloc_f(fptr_dev=vauxb_d, dimensions=[ldmx_b,ldmx_b,nat], omp_dev=omp_device)
        !$omp allocate allocator(omp_target_device_mem_alloc)
        allocate(vauxb_d(ldmx_b,ldmx_b,nat))
        do k=1, nat
           !$omp target teams distribute parallel do collapse(2) map(to:v%nsb)
           do j=1, ldmx_b
              do i=1, ldmx_b
                 vauxb_d(i,j,k) = CMPLX(v%nsb(i,j,current_spin,k))
              enddo
           enddo
           !$omp end target teams distribute parallel do
        enddo
#else
        ALLOCATE( vauxb_d(ldmx_b,ldmx_b,nat) )
        vauxb_d = CMPLX(v%nsb(:,:,current_spin,:))
#endif
      ENDIF
    ENDIF
  ENDIF
  !
  DO nt = 1, ntyp
     !
     ! Compute the action of the Hubbard potential on the KS wave functions:
     ! V_Hub |psip > = \sum v%ns |wfcU> <wfcU|psip>
     ! where v%ns = U ( delta/2 - rho%ns ) is computed in v_of_rho
     !
     IF ( is_hubbard(nt) ) THEN
        !
        ldim = 2*Hubbard_l(nt) + 1
        !
        DO na = 1, nat
           IF ( nt == ityp(na) ) THEN
              IF (gamma_only) THEN
                 !
#if defined(__OPENMP_GPU)
                 associate(r => proj_d%r)
                    !$omp target variant dispatch use_device_ptr(vns_d,r,rtemp_d)
                    CALL cublasDgemm( 'N','N', ldim,mps,ldim, 1.0_dp, &
                         vns_d(1,1,na), ldimax, &
                         r(offsetU(na)+1,1), nwfcU, 0.0_dp, rtemp_d, ldimaxt )
                    !$omp end target variant dispatch
                 endassociate
                 !
                 !$omp target variant dispatch use_device_ptr(wfcU,rtemp_d,hpsi_d)
                 CALL cublasDgemm( 'N','N', 2*np, mps, ldim, 1.0_dp, &
                      wfcU(1,offsetU(na)+1), 2*ldap, rtemp_d, ldimaxt, &
                      1.0_dp, hpsi_d, 2*ldap )
                 !$omp end target variant dispatch
#else
                 CALL cublasDgemm( 'N','N', ldim,mps,ldim, 1.0_dp, &
                      vns_d(1,1,na), ldimax, &
                      proj_d%r_d(offsetU(na)+1,1), nwfcU, 0.0_dp, rtemp_d, ldimaxt )
                 !
                 CALL cublasDgemm( 'N','N', 2*np, mps, ldim, 1.0_dp, &
                      wfcU_d(1,offsetU(na)+1), 2*ldap, rtemp_d, ldimaxt, &
                      1.0_dp, hpsi_d, 2*ldap )
#endif
                 !
              ELSE
                 !
#if defined(__OPENMP_GPU)
                 associate(k => proj_d%k)
                    !$omp target variant dispatch use_device_ptr(vaux_d,k,ctemp_d)
                    CALL cublasZgemm( 'N', 'N', ldim, mps, ldim, (1.0_dp,0.0_dp), &
                         vaux_d(:,:,na), ldimax, k(offsetU(na)+1,1), nwfcU, &
                         (0.0_dp,0.0_dp), ctemp_d, ldimaxt )
                    !$omp end target variant dispatch
                 endassociate
                 !
                 !$omp target variant dispatch use_device_ptr(wfcU,ctemp_d,hpsi_d)
                 CALL cublasZgemm( 'N', 'N', np, mps, ldim, (1.0_dp,0.0_dp), &
                      wfcU(1,offsetU(na)+1), ldap, ctemp_d, ldimaxt, &
                      (1.0_dp,0.0_dp), hpsi_d, ldap)
                 !$omp end target variant dispatch
#else
                 CALL cublasZgemm( 'N', 'N', ldim, mps, ldim, (1.0_dp,0.0_dp), &
                      vaux_d(:,:,na), ldimax, proj_d%k_d(offsetU(na)+1,1), nwfcU, &
                      (0.0_dp,0.0_dp), ctemp_d, ldimaxt )
                 !
                 CALL cublasZgemm( 'N', 'N', np, mps, ldim, (1.0_dp,0.0_dp), &
                      wfcU_d(1,offsetU(na)+1), ldap, ctemp_d, ldimaxt, &
                      (1.0_dp,0.0_dp), hpsi_d, ldap)
#endif
              ENDIF
           ENDIF
        ENDDO
        !
     ENDIF
     !
     ! If the background is used then compute extra
     ! contribution to the Hubbard potential
     !
     IF ( is_hubbard_back(nt) ) THEN
        !
        ldim = ldim_back(nt)
        !
        DO na = 1, nat
           IF ( nt == ityp(na) ) THEN
              !
              IF (gamma_only) THEN
                 !
                 ldim = 2*Hubbard_l_back(nt)+1
                 !
#if defined(__OPENMP_GPU)
                 associate(r => proj_d%r)
                    !$omp target variant dispatch use_device_ptr(vnsb_d,r,rtemp_d)
                    CALL cublasDgemm( 'N','N', ldim,mps,ldim, 1.0_dp, &
                         vnsb_d(1,1,na),ldmx_b, &
                         r(offsetU_back(na)+1,1), &
                         nwfcU, 0.0_dp, rtemp_d, ldimaxt )
                    !$omp end target variant dispatch
                 endassociate
                 !
                 !$omp target variant dispatch use_device_ptr(WfcU,rtemp_d,hpsi_d)
                 CALL cublasDgemm( 'N','N', 2*np, mps, ldim, 1.0_dp, &
                      wfcU(1,offsetU_back(na)+1), 2*ldap, rtemp_d, &
                      ldimaxt, 1.0_dp, hpsi_d, 2*ldap )
                 !$omp end target variant dispatch
#else
                 CALL cublasDgemm( 'N','N', ldim,mps,ldim, 1.0_dp, &
                      vnsb_d(1,1,na),ldmx_b, &
                      proj_d%r_d(offsetU_back(na)+1,1), &
                      nwfcU, 0.0_dp, rtemp_d, ldimaxt )
                 !
                 CALL cublasDgemm( 'N','N', 2*np, mps, ldim, 1.0_dp, &
                      wfcU_d(1,offsetU_back(na)+1), 2*ldap, rtemp_d, &
                      ldimaxt, 1.0_dp, hpsi_d, 2*ldap )
#endif
                 !
                 IF (backall(nt)) THEN
                    !
                    ldim  = 2*Hubbard_l1_back(nt)+1
                    ldim0 = 2*Hubbard_l_back(nt)+1
                    !
#if defined(__OPENMP_GPU)
                    associate(r => proj_d%r)
                       !$omp target variant dispatch use_device_ptr(vnsb_d,r,rtemp_d)
                       CALL cublasDgemm( 'N', 'N', ldim,mps,ldim, 1.0_dp,     &
                            vnsb_d(ldim0+1,ldim0+1,na),                       &
                            ldim_back(nt), r(offsetU_back1(na)+1,1), &
                            nwfcU, 0.0_dp, rtemp_d, ldimaxt )
                       !$omp end target variant dispatch
                    endassociate
                    !
                    !$omp target variant dispatch use_device_ptr(WfcU,rtemp_d,hpsi_d)
                    CALL cublasDgemm( 'N', 'N', 2*np, mps, ldim, 1.0_dp, &
                         wfcU(1,offsetU_back1(na)+1), 2*ldap, rtemp_d, &
                         ldimaxt, 1.0_dp, hpsi_d, 2*ldap )
                    !$omp end target variant dispatch
#else
                    CALL cublasDgemm( 'N', 'N', ldim,mps,ldim, 1.0_dp,     &
                         vnsb_d(ldim0+1,ldim0+1,na),                       &
                         ldim_back(nt), proj_d%r_d(offsetU_back1(na)+1,1), &
                         nwfcU, 0.0_dp, rtemp_d, ldimaxt )
                    !
                    CALL cublasDgemm( 'N', 'N', 2*np, mps, ldim, 1.0_dp, &
                         wfcU_d(1,offsetU_back1(na)+1), 2*ldap, rtemp_d, &
                         ldimaxt, 1.0_dp, hpsi_d, 2*ldap )
#endif
                 ENDIF
                 !
              ELSE
                 !
                 ldim = 2*Hubbard_l_back(nt)+1
                 !
#if defined(__OPENMP_GPU)
                 associate(k => proj_d%k)
                    !$omp target variant dispatch use_device_ptr(vauxb_d,k,ctemp_d)
                    CALL cublasZgemm( 'N', 'N', ldim,mps,ldim, (1.0_dp,0.0_dp),   &
                         vauxb_d(:,:,na), ldmx_b, k(offsetU_back(na)+1,1), &
                         nwfcU, (0.0_dp,0.0_dp), ctemp_d, ldimaxt )
                    !$omp end target variant dispatch
                 endassociate
                 !
                 !$omp target variant dispatch use_device_ptr(wfcU,ctemp_d,hpsi_d)
                 CALL cublasZgemm( 'N', 'N', np, mps, ldim, (1.0_dp,0.0_dp), &
                      wfcU_d(1,offsetU_back(na)+1), ldap, ctemp_d,           &
                      ldimaxt, (1.0_dp,0.0_dp), hpsi_d, ldap )
                 !$omp end target variant dispatch
#else
                 CALL cublasZgemm( 'N', 'N', ldim,mps,ldim, (1.0_dp,0.0_dp),     &
                      vauxb_d(:,:,na), ldmx_b, proj_d%k_d(offsetU_back(na)+1,1), &
                      nwfcU, (0.0_dp,0.0_dp), ctemp_d, ldimaxt )
                 !
                 CALL cublasZgemm( 'N', 'N', np, mps, ldim, (1.0_dp,0.0_dp), &
                      wfcU_d(1,offsetU_back(na)+1), ldap, ctemp_d,           &
                      ldimaxt, (1.0_dp,0.0_dp), hpsi_d, ldap )
#endif
                 IF (backall(nt)) THEN
                    !
                    ldim  = 2*Hubbard_l1_back(nt)+1
                    ldim0 = 2*Hubbard_l_back(nt)+1
                    !
#if defined(__OPENMP_GPU)
                    associate(k => proj_d%k)
                       !$omp target variant dispatch use_device_ptr(vauxb_d,k,ctemp_d)
                       CALL cublasZgemm( 'N', 'N', ldim,mps,ldim,(1.0_dp,0.0_dp), &
                            vauxb_d(ldim0+1,ldim0+1,na),ldmx_b,                   &
                            k(offsetU_back1(na)+1,1), nwfcU,               &
                            (0.0_dp,0.0_dp), ctemp_d, ldimaxt )
                       !$omp end target variant dispatch
                    endassociate
                    !
                    !$omp target variant dispatch use_device_ptr(wfcU,ctemp_d,hpsi_d)
                    CALL cublasZgemm( 'N', 'N', np, mps, ldim, (1.0_dp,0.0_dp), &
                         wfcU(1,offsetU_back1(na)+1), ldap, ctemp_d,          &
                         ldimaxt, (1.0_dp,0.0_dp), hpsi_d, ldap )
                    !$omp end target variant dispatch
#else
                    CALL cublasZgemm( 'N', 'N', ldim,mps,ldim,(1.0_dp,0.0_dp), &
                         vauxb_d(ldim0+1,ldim0+1,na),ldmx_b,                   &
                         proj_d%k_d(offsetU_back1(na)+1,1), nwfcU,             &
                         (0.0_dp,0.0_dp), ctemp_d, ldimaxt )
                    !
                    CALL cublasZgemm( 'N', 'N', np, mps, ldim, (1.0_dp,0.0_dp), &
                         wfcU_d(1,offsetU_back1(na)+1), ldap, ctemp_d,          &
                         ldimaxt, (1.0_dp,0.0_dp), hpsi_d, ldap )
#endif
                 ENDIF
                 !
              ENDIF
           ENDIF
        ENDDO
        !
     ENDIF
     !
  ENDDO
  !
  IF ( ANY(is_hubbard(:)) .OR. ANY(is_hubbard_back(:)) ) THEN
    IF (gamma_only) THEN
      DEALLOCATE( rtemp_d )
      IF (ANY(is_hubbard(:))) DEALLOCATE( vns_d )
      IF (ANY(is_hubbard_back(:))) DEALLOCATE( vnsb_d )
    ELSE
#if defined(__OPENMP_GPU)
      !call omp_target_free_f(fptr_dev=ctemp_d, omp_dev=omp_device)
      !IF (ANY(is_hubbard(:)))      call omp_target_free_f( fptr_dev=vaux_d, omp_dev=omp_device )
      !IF (ANY(is_hubbard_back(:))) call omp_target_free_f( fptr_dev=vauxb_d, omp_dev=omp_device )
      deallocate(ctemp_d)
      IF (ANY(is_hubbard(:)))      deallocate(vaux_d)
      IF (ANY(is_hubbard_back(:))) deallocate(vauxb_d)
#else
      DEALLOCATE( ctemp_d )
      IF (ANY(is_hubbard(:))) DEALLOCATE( vaux_d )
      IF (ANY(is_hubbard_back(:))) DEALLOCATE( vauxb_d )
#endif
    ENDIF
  ENDIF
  !
  !
  RETURN
  !
END SUBROUTINE vhpsi_U_gpu
!
END SUBROUTINE vhpsi_gpu
!-------------------------------------------------------------------------
