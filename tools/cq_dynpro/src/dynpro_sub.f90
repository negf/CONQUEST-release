MODULE dynpro_sub
  
  USE kinds,   ONLY : DP
  USE dynpro_inp
  USE constants,     ONLY : bohr => BOHR_RADIUS_ANGS
  
  IMPLICIT NONE
  
  SAVE

CONTAINS
  
  SUBROUTINE dynpro_sub_file(xyz_)
    
    !USE dynpro_inp
    USE constants,     ONLY : bohr => BOHR_RADIUS_ANGS
    
    
    IMPLICIT NONE

    REAL(DP), INTENT(in)  :: xyz_(3,nat)
    
    INTEGER :: i, j

    !ALLOCATE(xyz_xml(3,nat_xml))
    
    !xyz_xml = MATMUL(ht_xml,stau_xml)
    
    open(911, file=filesub_out, status='unknown', position='rewind')
    
    WRITE(911,*) nat
    WRITE(911,*) sub_nsp 

    !DO i = 1, nat 
    !   WRITE(911,'(a2,x,3f12.4,2x,i2,x,3i)')  &
    !        label(i), (xyz_xml(j,i)*bohr, j=1,3), ityp(i), i     
    !END DO

    DO i = 1, nat 
       WRITE(911,*)  &
            label(i), (xyz_(j,i), j=1,3), ityp(i), i     
    END DO

    !DEALLOCATE(xyz_xml)
    
    close(911)
    
  END SUBROUTINE dynpro_sub_file
  
  
  SUBROUTINE dynpro_sub_atom()

    USE dynpro_inp
    USE dynpro_mod,   ONLY : permute, very_tricky,  tricky_sort
    
    IMPLICIT NONE
     
    INTEGER :: i, j

    !ALLOCATE(ityp_old(nat))
    !ALLOCATE(label_old(nat))
    
    ityp_old  = ityp
    !nsp_old   = nsp
    label_old = label

    DEALLOCATE(ityp,label,atm,tmp_atm,na,atomic_number,syst_ID_spe,syst_ID_lab)

    ALLOCATE(ityp(nat))
    ALLOCATE(label(nat))
    ALLOCATE(xyz_sub(3,nat))
    
    
    OPEN(912, file=filesub_in, status='old', position='rewind')

    READ(912,*)
    READ(912,*) sub_nsp

    !DO i = 1, nat
    !   READ(912,'(a2,x,3f12.4,2x,i2,x,3i)') &
    !        label(i), (xyz_sub(j,i), j=1,3), ityp(i)     
    !END DO

    DO i = 1, nat
       READ(912,*) &
            label(i), (xyz_sub(j,i), j=1,3), ityp(i)
       print*, label(i), (xyz_sub(j,i), j=1,3), ityp(i)
    END DO


    CLOSE(912)
    !---------!

    ityp_old  = ityp
    !label_old = label
    
    ALLOCATE(atm(nsp+sub_nsp))
    ALLOCATE(tmp_atm(nsp+sub_nsp))
    ALLOCATE(na(nsp+sub_nsp))
    ALLOCATE(atomic_number(nsp+sub_nsp))
    ALLOCATE(syst_ID_spe(nsp+sub_nsp,2))
    ALLOCATE(syst_ID_lab(nsp+sub_nsp))
    
    atm(1:nsp) = atm_xml(1:nsp)


    print*, ' call very_tricky()', perm_n

    CALL VERY_TRICKY(sub_nsp, nsp+sub_nsp, nat, atm, ityp, label_old, label, &
         syst_ID_spe, syst_ID_lab, perm_xyz_tmp, perm_n)

    
    ALLOCATE(perm_xyz(perm_n))

    
    perm_xyz(:) = perm_xyz_tmp(1:perm_n)

    nsp = nsp + sub_nsp
    
    print*, ' call permute()'
    print*,  perm_xyz(:) 
    print*, syst_ID_spe
    print*, syst_ID_lab
    print*, label
    print*,
    print*, nsp
    print*,
    print*, ityp_old
    print*
    print*, ityp

    !CALL PERMUTE(nat, iteration, xyz, perm_xyz, perm_n)
    !CALL PERMUTE(nat, iteration, vel, perm_xyz, perm_n)
    !CALL PERMUTE(nat, iteration, for, perm_xyz, perm_n)

    !nsp = nsp + sub_nsp

    call tricky_sort(nsp, nat, iteration, xyz, vel, ityp_old, ityp, syst_ID_spe, &
         label_old, syst_ID_lab, .false.)

    ct_atom  = ct_atom  - perm_n
    sph_atom = sph_atom - perm_n
    
    DEALLOCATE(xyz_sub, perm_xyz, label_old, ityp_old)


    !print*, 'I am here'

    
  END SUBROUTINE dynpro_sub_atom
  
END MODULE dynpro_sub
