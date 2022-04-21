! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
MODULE test_hpsi_io
   USE kinds, ONLY: dp
   IMPLICIT NONE
   PRIVATE
   PUBLIC write_data_serial, read_data_serial

   CONTAINS

   SUBROUTINE write_data_serial(filename, npw, nbnd, data)
      implicit none
      character(len=*),intent(in) :: filename
      integer, intent(in)         :: npw, nbnd
      complex(dp),intent(in)      :: data(:,:)
      integer :: unit = 855
      !
      print *, trim(filename)
      open (unit = unit, file=trim(filename), form = 'unformatted', status='unknown')
      write (unit) npw, nbnd
      write (unit) data(1:npw, 1:nbnd)
      close (unit)
   END SUBROUTINE write_data_serial

   SUBROUTINE read_data_serial(filename, npw, nbnd, data )
      implicit none
      character(len=*),intent(in) :: filename
      integer,intent(in) :: npw, nbnd
      complex(dp) :: data(:,:)
      !
      integer :: unit = 856
      integer :: npw_, nbnd_
      open (unit=unit, file=trim(filename), form='unformatted', status='old')
      read(unit) npw_, nbnd_
      if ( npw_ /= npw .or. nbnd_ /= nbnd) call errore('read_data wrong dims:', trim(filename),1)
      read(unit) data(1:npw_, 1:nbnd_)
   END SUBROUTINE read_data_serial
END MODULE test_hpsi_io

!-----------------------------------------------------------------------
PROGRAM test_hpsi
   !-----------------------------------------------------------------------
   !
   ! Sample code, showing how to read QE data and re-use QE variables.
   ! This simple code
   ! 1. reads the data directory of QE, then
   ! 2. fills the hamiltonian matrix and diagonalizes it for each k-point
   !    (conventional, not iterative diagonalization)
   ! BEWARE: don't try this for large systems! the hamiltonian matrix is
   ! large, having Npw^2 elements, and diagonalization is slow, O(Npw^3)
   !
   ! Input: namelist &inputpp [outdir=...] [prefix=...] / as in QE input
   ! (default values as in QE).
   !
   USE io_global,  ONLY : ionode, ionode_id
   USE io_files,   ONLY : tmp_dir, prefix
   USE mp_global,  ONLY : mp_startup
   USE mp_images,  ONLY : intra_image_comm, root_image
   USE mp,         ONLY : mp_bcast
   USE environment,ONLY : environment_start, environment_end
   !
   IMPLICIT NONE
   !
   LOGICAL :: needwf = .true.
   LOGICAL :: write_ref = .false.
   INTEGER :: ios
   INTEGER :: ik_start = 1 , ik_stop = 10
   CHARACTER(LEN=256) :: outdir
   CHARACTER(LEN=256), EXTERNAL :: trimcheck
   !
   NAMELIST / inputpp / outdir, prefix, write_ref, ik_start, ik_stop
   !
   ! initialise environment
   !
   CALL mp_startup ( )
   CALL environment_start ( 'TEST_HPSI' )
   !
   IF ( ionode )  CALL input_from_file ( )
   !
   !   set default values for variables in namelist
   !
   prefix = 'pwscf'
   CALL get_environment_variable( 'ESPRESSO_TMPDIR', outdir )
   IF ( trim( outdir ) == ' ' ) outdir = './'
   !
   IF ( ionode )  THEN
      !
      !     reading the namelist inputpp
      !
      READ (5, inputpp, iostat = ios)
      !
      tmp_dir = trimcheck ( outdir )
      !
   ENDIF
   !
   CALL mp_bcast (ios, ionode_id, intra_image_comm)
   IF ( ios /= 0) CALL errore ('postproc', 'reading inputpp namelist', abs(ios))
   !
   ! ... Broadcast variables
   !
   CALL mp_bcast( tmp_dir, ionode_id, intra_image_comm )
   CALL mp_bcast( prefix, ionode_id, intra_image_comm )
   CALL mp_bcast( write_ref, ionode_id, intra_image_comm)
   !
   !   Read xml file, allocate and initialize general variables
   !
   CALL read_file_new ( needwf )
   !
   CALL run_tests (ik_start, ik_stop, write_ref )
   !
   CALL environment_end ( 'TEST_HPSI' )
   !
   CALL stop_pp()
   !
END PROGRAM test_hpsi

!-----------------------------------------------------------------------
SUBROUTINE run_tests (ik_in, ik_end, write_ref )
!-----------------------------------------------------------------------
   !
   USE constants,      ONLY : rytoev
   USE kinds,          ONLY : DP
   USE becmod,         ONLY : becp, calbec, allocate_bec_type
   USE fft_base,       ONLY : dfftp
   USE klist,          ONLY : xk, nks, nkstot, igk_k, ngk
   USE lsda_mod,       ONLY : nspin, isk, current_spin
   USE io_files,       ONLY : restart_dir
   USE scf,            ONLY : vrs, vltot, v, kedtau
   USE gvecs,          ONLY : doublegrid
   USE uspp,           ONLY : nkb, vkb
   USE uspp_init,      ONLY : init_us_2
   USE wvfct,          ONLY : npwx, nbnd, current_k
   USE mp_bands,       ONLY : me_bgrp, root_bgrp, intra_bgrp_comm
   USE wavefunctions,  ONLY : evc, psic
   USE pw_restart_new, ONLY : read_collected_wfc
   USE mp,             ONLY : mp_sum
   USE test_hpsi_io
   !
   IMPLICIT NONE
   !
   INTEGER, INTENT(IN)  :: ik_in, ik_end
   LOGICAL, INTENT(IN)  :: write_ref
   !
   INCLUDE 'laxlib.fh'
   !
   COMPLEX(DP), ALLOCATABLE :: aux(:,:), aux_check(:,:)
   COMPLEX(DP), ALLOCATABLE :: hc(:,:), sc(:,:), vc(:,:)
   REAL(DP),    ALLOCATABLE :: en(:)
   INTEGER :: ik, npw, ik_start, ik_stop, i, ibnd
   CHARACTER(LEN=320) ::  filename
   LOGICAL             :: ionode = .TRUE.
   COMPLEX(DP)         :: res
   CHARACTER(LEN=6), EXTERNAL :: int_to_char
   !
   call ik_check(ik_in, ik_start, nkstot, 1, write_ref)
   call ik_check(ik_end, ik_stop, nkstot, nkstot, write_ref)
   print *, ik_start, ik_stop

   ALLOCATE( aux(npwx, nbnd ) )
   ALLOCATE( hc( nbnd, nbnd) )
   ALLOCATE( sc( nbnd, nbnd) )
   ALLOCATE( vc( nbnd, nbnd) )
   ALLOCATE( en( nbnd ) )
   CALL allocate_bec_type(nkb, nbnd, becp )
   if (use_gpu) CALL allocate_bec_type_gpu(nkb, nbnd, becp)
   CALL set_vrs(vrs,vltot,v%of_r,kedtau,v%kin_r,dfftp%nnr,nspin,doublegrid)

   DO ik = ik_start, ik_stop
      !
      CALL read_collected_wfc( restart_dir() , ik, evc )
      !
      npw = ngk(ik)
      current_k=ik
      current_spin  = isk(ik)
      !
      CALL init_us_2(npw, igk_k(1,ik), xk (1, ik), vkb)
      CALL calbec( npw, vkb, evc, becp)
      CALL g2_kin(ik)
      !
      CALL h_psi( npwx, npw, nbnd, evc, aux )
      filename = trim(restart_dir())//"hpsi_"//trim(int_to_char(ik))//".dat"
      if (write_ref ) then
         call write_data_serial(filename, npw, nbnd, aux)
      else
         allocate (aux_check(npw, nbnd) )
         call read_data_serial (filename, npw, nbnd, aux_check)
         res = (0.d0, 0.d0)
         do ibnd = 1, nbnd
            res = res + dot_product(aux_check(:,ibnd) -aux(:,ibnd),aux_check(:,ibnd) - aux(:,ibnd))
         end do

         if (ionode) then
            print '("check on hpsi, sum of residuals for k=: ", I5, 2F16.8)',  ik, res
         end if
      end if
      !
      CALL calbec ( npw, evc, aux, hc )
      filename = trim(restart_dir())//"hc_"//trim(int_to_char(ik))//".dat"
      if (write_ref) then
         call write_data_serial(filename, nbnd, nbnd, hc)
      else
         call read_data_serial(filename, nbnd, nbnd, vc )
         res = (0.d0, 0.d0)
         do i = 1, size(hc,2)
            res= res + DOT_PRODUCT(hc(:,i)-vc(:,i), hc(:,i) - vc(:,i))
         end do
         if (ionode) then
            print '("check on calbec, sum of residuals on hc for k=: ", I5, 2F16.8)', ik, res
         end if
      end if

      filename = trim(restart_dir())//"spsi_"//trim(int_to_char(ik))//".dat"
      CALL s_psi( npwx, npw, nbnd, evc, aux )
      if (write_ref) then
         call write_data_serial(filename, npw, nbnd, aux)
      else
         call read_data_serial (filename, npw, nbnd, aux_check)
         res = (0.d0, 0.d0)
         do ibnd = 1, nbnd
            res = res + dot_product(aux_check(:,ibnd) - aux(:,ibnd), aux_check(:,ibnd) - aux(:,ibnd))
         end do
         if (ionode) then
            print '("check on spsi, sum of residuals for k=: ", I5, 2F16.8)',  ik, res
         end if
      end if
      !
      CALL calbec ( npw, evc, aux, sc )
      !
      CALL diaghg( nbnd, nbnd, hc, sc, nbnd, en, vc, me_bgrp, &
          root_bgrp, intra_bgrp_comm )
       !
      print '(/,12x,"k =",3f7.4," (",i6," PWs)  bands (eV):",/)', xk(:,ik),npw
      print '(8f9.4)', en(:)*rytoev
      !
   END DO
   !
   contains
      !
      SUBROUTINE ik_check(in, out, ik_max, ik_default, lwrite)
         implicit none
         integer, intent(in)   :: in, ik_max, ik_default
         integer, intent(out)  :: out
         logical, intent(in)   :: lwrite
         if (lwrite) then
            print '("Writing references, ik from input neglected")'
            out = ik_default
         else if (in <=0) then
            out = ik_default
         else if (in > 0 .and. in <= nkstot ) then
            out = in
         else if (in > nkstot) then
            print '("Max ik in data set is ",I5)', nkstot
            call errore ('test_hpsi','wrong ik_in', in)
         end if
      END SUBROUTINE ik_check
END SUBROUTINE run_tests
