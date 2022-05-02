!
! Copyright (C) 2001-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine force_corr_gpu (forcescc)
  !-----------------------------------------------------------------------
  !   This routine calculates the force term vanishing at full
  !     self-consistency. It follows the suggestion of Chan-Bohnen-Ho
  !     (PRB 47, 4771 (1993)). The true charge density is approximated
  !     by means of a free atom superposition.
  !     (alessio f.)
  ! Uses superposition of atomic charges contained in the array rho_at
  ! and read from pseudopotential files
  !
  USE kinds,                ONLY : DP
  USE constants,            ONLY : tpi
  USE atom,                 ONLY : msh, rgrid
  USE uspp_param,           ONLY : upf
  USE ions_base,            ONLY : nat, ntyp => nsp, ityp, tau
  USE cell_base,            ONLY : tpiba
  USE fft_base,             ONLY : dfftp
  USE fft_interfaces,       ONLY : fwfft
#if defined(__OPENMP_GPU)
  USE gvect,                ONLY : ngm, gstart, g, gl, ngl, igtongl
#else
  USE gvect,                ONLY : ngm, gstart, g, ngl, gl_d, igtongl_d
#endif
  USE lsda_mod,             ONLY : nspin
  USE scf,                  ONLY : vnew
  USE control_flags,        ONLY : gamma_only
  USE mp_bands,             ONLY : intra_bgrp_comm
  USE mp,                   ONLY : mp_sum
  USE devxlib_buffers,      ONLY : dev_buf => gpu_buffer
#if !defined(__OPENMP_GPU)
  USE gvect_gpum,           ONLY : g_d
#endif
  !
  USE simpsn_gpum,          ONLY : simpsn_gpu_dev
#if defined(__CUDA)
  USE cudafor
  USE cublas
#endif
  !
  implicit none
  !
  real(DP) :: forcescc (3, nat)
#if defined(__CUDA) || defined(__OPENMP_GPU)
  !
  real(DP), pointer, contiguous ::   rhocgnt_d(:),  aux_d (:,:)
#if !defined(__OPENMP_GPU)
  real(DP), pointer, contiguous ::   r_d(:), rab_d(:), rhoat_d(:), tau_d(:,:)
#endif
  complex(DP),pointer      ::   psic_d(:)
  integer, pointer      :: nl_d (:)
  real(DP), allocatable :: tmp(:)

  ! work space
  real(DP) ::  gx, arg, fact, forcesccx, forcesccy, forcesccz
  ! temp factors
  integer :: ir, isup, isdw, ig, nt, na, ndm, ierr, msh_nt, igg, glblock,igl, ierrs(7)
#if defined(__CUDA)
  ATTRIBUTES(DEVICE)  :: rhocgnt_d, aux_d, psic_d,r_d, rab_d, rhoat_d, nl_d, tau_d
#endif
#if defined(__OPENMP_GPU)
  INTEGER :: i
#endif

  ! counters
  !
  ! vnew is V_out - V_in, psic is the temp space
  !
#if !defined(__OPENMP_GPU)
  nl_d => dfftp%nl_d
#endif
  CALL dev_buf%lock_buffer(psic_d, size(vnew%of_r,1), ierrs(1) )
  if (nspin == 1 .or. nspin == 4) then
#if defined(__OPENMP_GPU)
     !$omp target teams distribute parallel do map(to:vnew%of_r)
     do i=1, size(vnew%of_r,1)
        psic_d(i) = vnew%of_r(i, 1)
     enddo
#else
     psic_d(:) = vnew%of_r (:, 1)
#endif
  else
     isup = 1
     isdw = 2
#if defined(__OPENMP_GPU)
     !$omp target teams distribute parallel do
     do i = 1, size(vnew%of_r,1)
        psic_d(i) = (vnew%of_r(i,isup) + vnew%of_r(i,isdw)) * 0.5d0
     enddo
     !$omp end target teams distribute parallel do
#else
     psic_d(:) = (vnew%of_r (:, isup) + vnew%of_r (:, isdw)) * 0.5d0
#endif
  end if
  !
  ndm = MAXVAL ( msh(1:ntyp) )
  CALL dev_buf%lock_buffer ( aux_d, [ndm, ngl], ierrs(2) )
  CALL dev_buf%lock_buffer ( rhocgnt_d, ngl, ierrs(3) )
#if !defined(__OPENMP_GPU)
  CALL dev_buf%lock_buffer ( r_d, ndm, ierrs(4))
  CALL dev_buf%lock_buffer ( rab_d, ndm, ierrs(5) )
  CALL dev_buf%lock_buffer ( rhoat_d, ndm, ierrs(6))
  !
  CALL dev_buf%lock_buffer ( tau_d, [3,nat], ierrs(7))
  !
  IF (ANY(ierrs /= 0)) CALL errore('force_corr_gpu', 'cannot allocate buffers', -1)
#endif
  !
#if defined(__OPENMP_GPU)
  !$omp target data map(to:tau, rgrid(nt)%r, rgrid(nt)%rab, upf(nt)%rho_at, dfftp%nl)
#else
  tau_d(1:3,1:nat)=tau(1:3,1:nat)
#endif
  !
  !$omp dispatch
  CALL fwfft (1, psic_d, dfftp)
  if (gamma_only) then
     fact = 2.d0
  else
     fact = 1.d0
  end if
  !
  do nt = 1, ntyp
     !
     ! Here we compute the G.ne.0 term
     !
     msh_nt = msh(nt)
#if !defined(__OPENMP_GPU)
     r_d (1:msh_nt) =  rgrid(nt)%r(1:msh_nt)
     rab_d(1:msh_nt) = rgrid(nt)%rab(1:msh_nt)
     rhoat_d (1:msh_nt) =  upf(nt)%rho_at (1:msh_nt)
#endif
     !
     !$cuf kernel do(1)
     !$omp target teams distribute parallel do
     do ig = 1, ngl
#if defined(__OPENMP_GPU)
        gx = sqrt (gl (ig) ) * tpiba
#else
        gx = sqrt (gl_d (ig) ) * tpiba
#endif
        do ir = 1, msh_nt
#if defined(__OPENMP_GPU)
           if (rgrid(nt)%r(ir) .lt.1.0d-8) then
              aux_d (ir,ig) = upf(nt)%rho_at (ir)
#else
           if (r_d(ir) .lt.1.0d-8) then
              aux_d (ir,ig) = rhoat_d (ir)
#endif
           else
#if defined(__OPENMP_GPU)
              aux_d (ir,ig) = upf(nt)%rho_at(ir) * sin(gx*rgrid(nt)%r(ir)) / (rgrid(nt)%r(ir)*gx)
#else
              aux_d (ir,ig) = rhoat_d (ir) * sin(gx*r_d(ir)) / (r_d(ir)*gx)
#endif
           endif
        enddo

#if defined(__OPENMP_GPU)
        call simpsn_gpu_dev(msh_nt, aux_d(:,ig), rgrid(nt)%rab, rhocgnt_d(ig) )
#else
        call simpsn_gpu_dev(msh_nt, aux_d(:,ig), rab_d, rhocgnt_d(ig) )
#endif

     enddo
     !
     do na = 1, nat
        if (nt.eq.ityp (na) ) then
           forcescc (1:3, na) = 0.0_DP
           forcesccx = 0.0d0
           forcesccy = 0.0d0
           forcesccz = 0.0d0
           !$cuf kernel do(1)
           !$omp target teams distribute parallel do reduction(+:forcesccx) &
           !$omp &                                   reduction(+:forcesccy) &
           !$omp &                                   reduction(+:forcesccz) &
           !$omp & map(tofrom:forcesccx, forcesccy, forcesccz)
           do ig = gstart, ngm
#if defined(__OPENMP_GPU)
              arg = (g(1, ig) * tau(1, na) + g(2, ig) * tau(2, na) &
                   + g(3, ig) * tau(3, na) ) * tpi
              forcesccx = forcesccx+ REAL(fact * &
                      rhocgnt_d (igtongl(ig) ) * CMPLX(sin(arg),cos(arg),kind=DP) * &
                      g(1,ig) * tpiba * CONJG(psic_d(dfftp%nl(ig))))
              forcesccy = forcesccy + REAL(fact * &
                      rhocgnt_d (igtongl(ig) ) * CMPLX(sin(arg),cos(arg),kind=DP) * &
                      g(2,ig) * tpiba * CONJG(psic_d(dfftp%nl(ig))))
              forcesccz = forcesccz + REAL(fact * &
                      rhocgnt_d (igtongl(ig) ) * CMPLX(sin(arg),cos(arg),kind=DP) * &
                      g(3,ig) * tpiba * CONJG(psic_d(dfftp%nl(ig))))
#else
              arg = (g_d(1, ig) * tau_d (1, na) + g_d (2, ig) * tau_d (2, na) &
                   + g_d (3, ig) * tau_d (3, na) ) * tpi
              forcesccx = forcesccx+ REAL(fact * &
                      rhocgnt_d (igtongl_d(ig) ) * CMPLX(sin(arg),cos(arg),kind=DP) * &
                      g_d(1,ig) * tpiba * CONJG(psic_d(nl_d(ig))))
              forcesccy = forcesccy + REAL(fact * &
                      rhocgnt_d (igtongl_d(ig) ) * CMPLX(sin(arg),cos(arg),kind=DP) * &
                      g_d(2,ig) * tpiba * CONJG(psic_d(nl_d(ig))))
              forcesccz = forcesccz + REAL(fact * &
                      rhocgnt_d (igtongl_d(ig) ) * CMPLX(sin(arg),cos(arg),kind=DP) * &
                      g_d(3,ig) * tpiba * CONJG(psic_d(nl_d(ig))))
#endif
           enddo
           forcescc(1, na) = forcesccx
           forcescc(2, na) = forcesccy
           forcescc(3, na) = forcesccz
        endif
     enddo
  enddo
  !
#if defined(__OPENMP_GPU)
  !$omp end target data
#else
  call dev_buf%release_buffer ( r_d, ierr )
  call dev_buf%release_buffer ( rab_d, ierr )
  call dev_buf%release_buffer ( rhoat_d, ierr )
  call dev_buf%release_buffer ( tau_d, ierr )
#endif
  call dev_buf%release_buffer ( psic_d, ierr )
  !
  call dev_buf%release_buffer ( aux_d, ierr )
  !
  call dev_buf%release_buffer ( rhocgnt_d, ierr )
  !
  call mp_sum(  forcescc, intra_bgrp_comm )
#endif
  !
  return
end subroutine force_corr_gpu

