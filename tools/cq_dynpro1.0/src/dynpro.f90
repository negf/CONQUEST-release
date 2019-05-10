!#/bin/csh
!
! Copyright (C) 2002-2008 Quantum-Espresso group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!#include "f_defs.h"

PROGRAM dynpro_postproc

  USE kinds,         ONLY : DP
  USE constants,     ONLY : bohr => BOHR_RADIUS_ANGS
  USE constants,     ONLY : ry   => AUTOEV
  USE constants,     ONLY : ps   => AU_PS
  USE io_files,      ONLY : outdir


  !USE io_global,     ONLY : io_global_start
  !USE mp_global,     ONLY : mp_global_start, mp_global_end
  !USE mp_global,     ONLY : mp_startup, mp_global_end
  !USE mp,            ONLY : mp_start

  USE dynpro_mod
  USE dynpro_inp
  USE dynpro_sub,   ONLY : dynpro_sub_file, dynpro_sub_atom
  USE dynpro_set,   ONLY : dynpro_read_input !, dynpro_read_xml
  USE dynpro_sph,   ONLY : dynpro_sphere   
  !USE dynpro_rax,   ONLY : tensor, tensor_acf
  USE dynpro_cqc,   ONLY : dynpro_read_cqinp, dynpro_read_cqCoordForce

  
  IMPLICIT NONE

  INTEGER                       :: cunit, punit, funit, vunit, eunit

  !==---- CENTER cell --------==
  REAL(DP),        ALLOCATABLE  :: ct_abc(:,:), ct_xyz(:,:)

  !==---- RDF calculation ----==
  INTEGER,         ALLOCATABLE  :: event(:,:,:)
  INTEGER,         ALLOCATABLE  :: dN(:,:,:,:)
  INTEGER                       :: rd_count
  CHARACTER(len=2)              :: rd_lab    
  REAL(DP)                      :: vol
  
  !------ Sampling -----------==
  REAL(DP),        ALLOCATABLE  :: sp_xyz(:,:), sp_vel(:,:), xyz_(:,:)
  REAL(DP),        ALLOCATABLE  :: sp_for(:,:), sp_cel(:,:), sp_evp(:,:)

  INTEGER                       :: i, j, k, l, m, n, o, p, tmp_xml

  
  !==---------------------------------------------------------------------==

#if defined __INTEL
  ! ... Intel compilers v .ge.8 allocate a lot of stack space
  ! ... Stack limit is often small, thus causing SIGSEGV and crash
  CALL remove_stack_limit ( )
#endif  
      
  !  see cprstart.f90 for the meaning of the following 4 calls
  !CALL mp_start()
  !CALL mp_env( nproc, mpime, world )
  !CALL mp_global_start( root, mpime, world, nproc )
  !CALL io_global_start( mpime, root )
  !CALL get_env( 'ESPRESSO_TMPDIR', outdir )

  !  initialize mpi
  !CALL mp_start( nproc, mpime, world )
  !CALL mp_global_start( root, mpime, world, nproc )
  !CALL io_global_start( mpime, root )
  !  initialize mpi 
  !CALL mp_startup  ( )                                                             

  CALL get_env( 'ESPRESSO_TMPDIR', outdir )

  IF ( TRIM( outdir ) == ' ' ) outdir = './tmp/'

  !==---------------------------------------------------------------------==
  print*, 'call dynpro_read_cqinp() '
  CALL dynpro_read_cqinp( )
  
  print*, 'call dynpro_read_input() '  
  CALL dynpro_read_input( )

  print*, 'call dynpro_read_cqCoordForce() '  
  if (lcq_coordforce) CALL dynpro_read_cqCoordForce( )


  !==---------------------------------------------------------------------==
  !<<< variables settings >>>
  

  !==---------------------------------------------------------------------==
!!$  CALL dynpro_read_xml()
!!$
!!$  filepp =
!!$  ALLOCATE(ityp_xml(nat_xml))
!!$  ALLOCATE(label_xml(nat_xml))
!!$  ALLOCATE(atm_xml(nsp_xml))
!!$  ityp_xml = 0 ; label_xml = 'XX' ; atm_xml = 'XX'
!!$
!!$  ALLOCATE(stau_xml( 3, nat_xml))
!!$  ALLOCATE(svel_xml( 3, nat_xml))
!!$  ALLOCATE(force_xml(3, nat_xml))
!!$  stau_xml = 0.0d0 ; svel_xml = 0.0d0 ; force_xml = 0.0d0
!!$  tmp_xml = MODULO((i+j+k+l+m), 5)
!!$  n_xml = i
!!$


  print*, 'call check_file() '
  CALL check_file(filecel,3      ,1,i)
  CALL check_file(filepos,nat_xml,1,j)
  CALL check_file(filevel,nat_xml,1,k)
  CALL check_file(filefor,nat_xml,1,l)
  CALL check_file(fileevp,1      ,0,m)

  tmp_xml = MODULO((i+j+k+l+m), 5)

  n_xml   = i

  nat_xml = nat_cqinp
  nsp_xml = nsp_cqinp

  ALLOCATE(ityp_xml(nat_xml))
  ALLOCATE(label_xml(nat_xml))
  ALLOCATE(atm_xml(nsp_xml))
  
  ityp_xml  = ityp_cqcoord
  label_xml = label_cqinp
  atm_xml   = atm_cqinp

  print*, '    n steps = ',   n_xml  
  print*, '    n atoms = ', nat_xml
  print*, '    species = ', nsp_xml
  print*, '    atoms   = ', atm_xml
  print*, '    label   = ', label_xml
  print*, '    ityp    = ', ityp_xml
  
  !==---------------------------------------------------------------------==

  !<<< variables settings >>>
  iteration = n_xml          !---XML---!
  nat       = nat_xml        !---XML---!
  nsp       = nsp_xml        !---XML---!
  
  ALLOCATE(xyz(3,nat*iteration))
  ALLOCATE(abc(3,nat*iteration))
  ALLOCATE(vel(3,nat*iteration))
  ALLOCATE(for(3,nat*iteration))
  ALLOCATE(cel(3,3*iteration))
  ALLOCATE(evp(8,iteration))  
  ALLOCATE(ct_abc(3,nat*iteration))
  ALLOCATE(ct_xyz(3,nat*iteration))

  print*, 'call dynpro_files( ) '  
  punit = 10 
  CALL read_file(filepos, punit, 1, 3, nat, iteration, xyz)
  vunit = 11
  CALL read_file(filevel, vunit, 1, 3, nat, iteration, vel)
  funit = 12 
  CALL read_file(filefor, funit, 1, 3, nat, iteration, for)
  cunit = 13 
  CALL read_file(filecel, cunit, 1, 3, 3,  iteration, cel)
  eunit = 14
  CALL read_file(fileevp, eunit, 0, 8, 1,  iteration, evp)

  !==---------------------------------------------------------------------==
  !<<< sampling function >>>
  IF (lsample) THEN
     
     print*, 'call check_sampling( ) '  
     CALL check_sampling(sp_start, sp_end, evp, iteration, iprint, &
          sp_n, start, end, interval, time_step, time_dt)

     ALLOCATE(sp_xyz(3,nat*iteration))
     ALLOCATE(sp_vel(3,nat*iteration))
     ALLOCATE(sp_for(3,nat*iteration))
     ALLOCATE(sp_cel(3,3*iteration))
     ALLOCATE(sp_evp(8,iteration))  
     
     sp_xyz = xyz
     sp_vel = vel
     sp_for = for
     sp_cel = cel 
     sp_evp = evp
     
     DEALLOCATE(xyz, abc, vel, for, cel, evp)
     DEALLOCATE(ct_abc,ct_xyz)     
 
     ALLOCATE(xyz(3,nat*sp_n))
     ALLOCATE(abc(3,nat*sp_n))
     ALLOCATE(vel(3,nat*sp_n))
     ALLOCATE(for(3,nat*sp_n))
     ALLOCATE(cel(3,3*sp_n))
     ALLOCATE(evp(8,sp_n))
     ALLOCATE(ct_abc(3,nat*sp_n))
     ALLOCATE(ct_xyz(3,nat*sp_n))

     print*, 'call sampling( ) '  
     CALL sampling(sp_xyz, 3, nat, iteration, sp_n, start, end, & 
          interval, iprint, xyz)
     CALL sampling(sp_cel, 3, 3, iteration, sp_n, start, end, & 
          interval, iprint, cel)
     CALL sampling(sp_vel, 3, nat, iteration, sp_n, start, end, & 
          interval, iprint, vel)
     CALL sampling(sp_for, 3, nat, iteration, sp_n, start, end, & 
          interval, iprint, for)
     CALL sampling(sp_evp, 8, 1, iteration, sp_n, start, end, & 
          interval, iprint, evp)

     iteration = sp_n
     
     DEALLOCATE(sp_xyz)
     DEALLOCATE(sp_vel)
     DEALLOCATE(sp_for)
     DEALLOCATE(sp_cel)
     DEALLOCATE(sp_evp)  

     ELSE
        interval = 1
  END IF   

  !==---------------------------------------------------------------------==
  ! Add distance / angle analysis here

  !==---------------------------------------------------------------------==

  !<<< variables settings >>>
  IF (rmax == 0.0d0) THEN
     rmax = cel(1,1)
  END IF

  IF( .not. lau_unit) THEN
     cel  = cel*bohr 
     xyz  = xyz*bohr
     for  = for*ry/bohr
     rmax = rmax*bohr
     rmin = rmin*bohr
     vel  = vel*bohr/ps
  END IF

  print*, ps
  print*, bohr

  ALLOCATE(ityp(nat))
  ALLOCATE(ityp_old(nat))

  ALLOCATE(label(nat))
  ALLOCATE(label_old(nat))

 
  ALLOCATE(atm(nsp))
  ALLOCATE(tmp_atm(nsp))
  ALLOCATE(na(nsp))
  ALLOCATE(atomic_number(nsp))
  ALLOCATE(syst_ID_spe(nsp,2))
  ALLOCATE(syst_ID_lab(nsp)) 

  !print*, n_xml+3, 'toto'
  
  ityp    = ityp_xml     !---XML---!
  atm     = atm_xml      !---XML---!

  print*, 'dynpro: atm'
  print*, atm
  !
  print*, 'dynpro: ityp'
  print*, ityp
  !
  print*, 'call tricky( ) '
  !
  CALL tricky(nsp, nat, atm, ityp, label, syst_ID_spe, syst_ID_lab)

  print*, 'dynpro before tricky_sort: label_xml'
  print*, label_xml
  !
  print*, 'dynpro  before tricky_sort: ityp_xml'
  print*, ityp_xml
  !
  print*, 'call tricky_sort( ) '
  !
  CALL tricky_sort(nsp, nat, iteration, xyz, vel, ityp_xml, ityp, syst_ID_spe, &
       label_xml, syst_ID_lab, .true.)

  print*, 'dynpro: label'
  print*, label
  !
  print*, 'dynpro: ityp'
  print*, ityp
  
  ALLOCATE(xyz_(3,nat))

  xyz_(:,:) = xyz(:,1:nat)
  
  CALL dynpro_sub_file(xyz_)

  DEALLOCATE(xyz_)

  
  IF (lsub) THEN
     print*, 'call dynpro_sub_atom( ) '  
     CALL dynpro_sub_atom()
     
  END IF

  !print*, n_xml+6

  
  !==---------------------------------------------------------------------==  
  print*, 'call xyz_to_abc( ) '  
  CALL xyz_to_abc(nat, iteration, xyz, abc, cel)

  IF (lcenter) THEN
     print*, 'entering CENTER '            
     print*, '     call center_cell( ) '       
     CALL center_cell(nat, iteration, ct_atom, cel, abc, ct_abc, ct_xyz)
     print*, '     call write_xyz( ) '       
     CALL write_xyz(nat, iteration, label, ct_xyz, evp, filemol, 21)
     print*, 'leaving CENTER '   
  ELSE
     print*, ' call write_xyz( ) '       
     CALL write_xyz(nat, iteration, label, xyz, evp, filemol, 21)
     print*, ' leaving CENTER '     
  END IF

  
  print*, 'call write_evp( ) '       
  CALL write_evp(evp, iteration, time_step, iprint, lgrace, filener)

  !==---------------------------------------------------------------------==

  IF (lrdf) THEN
     
     print*, 'entering RDF '       
     print*, rmin, rmax, nr

     rd_count = 0
     k        = 1
     DO i = 1, nsp     
        IF (rd_atom == syst_ID_spe(i,1)) THEN
           rd_count = syst_ID_spe(i,2)
           rd_lab   = syst_ID_lab(i)
        END IF
     END DO
     
     ALLOCATE(dN(iteration,rd_count,nsp,nr))
     ALLOCATE(event(iteration,rd_count,nsp))

     print*, '   call dNab '   
     CALL dNab (nsp, nat,  iteration, nr, rmin, rmax,  & 
          ityp, syst_ID_spe, cel, abc, &
          rd_atom, rd_count, dN, event, vol)
     print*, '   call RDF '        
     CALL RDF (dN, event, iteration, rd_atom, rd_count, rd_lab, &
          nat, nsp, nr, rmin, rmax, vol, & 
          syst_ID_spe, syst_ID_lab, filerdf,  lau_unit)
     DEALLOCATE(dN)
     DEALLOCATE(event)
     print*, 'leaving RDF '     
  END IF


  !==---------------------------------------------------------------------==

  IF (lvacf) THEN
     print*, 'entering VACF ' 
     unit_vac = 777
     unit_psd = 888      
     !
     ALLOCATE(C_nsp(nsp,iteration))
     ALLOCATE(C_tot(iteration))  
     !
     print*, '   call VACF ' 
     !
     CALL VACF(nat, iteration, vel, nsp, syst_ID_spe, syst_ID_lab, &
          C_nsp, C_tot, filevac, unit_vac, time_step, interval, iprint)
     !
     DEALLOCATE(C_nsp)
     DEALLOCATE(C_tot)
     !
     print*, 'leaving VACF ' 
     !
  END IF
  !
  IF (lpsd) THEN
     !
     print*, '   call PSD '
     print*, '      call check_power2'
     !
     CALL check_power2(iteration)
     !
     CALL PSD(0, nsp, lspe, iteration, syst_ID_lab, syst_ID_spe, &
          filevac, loverlap, n_seg, win, unit_vac, m_cof, &
          time_step, interval, iprint, unit_psd, filepsd, qcfactor, temp)
     !  
  END IF
  !    
  !DEALLOCATE(C_nsp)
  !DEALLOCATE(C_tot)
  !
 
!!$  IF (lpsd) THEN
!!$     print*, '   call PSD '
!!$     
!!$     print*, '      call check_power2'
!!$     
!!$     CALL check_power2(iteration)
!!$     
!!$     CALL PSD(0, C_tot, nsp, lspe, iteration, syst_ID_lab, syst_ID_spe, &
!!$          filevac, loverlap, n_seg, win, unit_vac, m_cof, &
!!$          time_step, interval, iprint, unit_psd, filepsd, qcfactor, temp)
!!$
!!$  END IF

  !==---------------------------------------------------------------------==

  IF (lsphere) THEN
     print*, 'entering SPHERE ' 
     CALL dynpro_sphere()
     print*, 'leaving SPHERE ' 
  END IF
  
  !==---------------------------------------------------------------------==
  
  !IF (lnmr .and. lsample) THEN
  !   CALL tensor()
  !   CALL tensor_acf()
  !END IF

  !==---------------------------------------------------------------------==
  DEALLOCATE(label_xml)

  DEALLOCATE(ityp_xml)
  DEALLOCATE(atm_xml)
  
  DEALLOCATE(xyz)
  DEALLOCATE(abc)
  DEALLOCATE(vel)
  DEALLOCATE(for)
  DEALLOCATE(cel)
  DEALLOCATE(evp)  
  DEALLOCATE(ct_abc)
  DEALLOCATE(ct_xyz)

  DEALLOCATE(ityp)
  DEALLOCATE(label)
  DEALLOCATE(atm)
  DEALLOCATE(tmp_atm)
  DEALLOCATE(na)
  DEALLOCATE(atomic_number)
  DEALLOCATE(syst_ID_spe)
  DEALLOCATE(syst_ID_lab)   

  !CALL mp_global_end()

END PROGRAM dynpro_postproc



