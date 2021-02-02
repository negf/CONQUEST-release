module kubo_greenwood
  use datatypes

  implicit none
  
  real(double), parameter :: kb = 8.617333262145d-5 / 27.21138602d0
  complex(double_cplx), target, allocatable, save :: ngradm_summed(:,:), sigma(:,:)
  integer :: mat_dS2(3)
  
  contains
    
    subroutine dealloc_matdS2()
      use mult_module, only: free_temp_matrix
      implicit none
      
      integer :: i
      
      do i = 3, 1, -1
        call free_temp_matrix(mat_dS2(i))
      end do
      
      
    end subroutine dealloc_matdS2
  
    subroutine get_conductivity(w, expH, Ef, wtk, kp, kps, spin)
    
      use numbers
      use datatypes
      use GenComms, only: cq_abort, inode, ionode, myid,io_lun, gsum
      use PAO_grid_transform_module, only: single_PAO_to_grid, single_PAO_to_grad      
      use matrix_data, only: Hrange 
      use global_module, only: sf, atomf
      use mult_module, only: SF_to_AtomF_transform, matrix_scale,allocate_temp_matrix, &
       free_temp_matrix, matrix_product_trace, matS
      use calc_matrix_elements_module, only: act_on_vectors_new
      use functions_on_grid, only: atomfns, allocate_temp_fn_on_grid, free_temp_fn_on_grid, &
        gridfunctions, fn_on_grid       
      use set_bucket_module, only: rem_bucket, atomf_H_atomf_rem 
      use ScalapackFormat, only: matrix_size
      use calc_matrix_elements_module, only: get_matrix_elements_new
      use set_bucket_module,           only: rem_bucket, atomf_H_atomf_rem ,atomf_atomf_rem
      use negf_output_module 
      use MPI
      
      
      implicit none
    
      complex(double_cplx), intent(in) :: expH(:,:,:)
      real(double), intent(in) :: wtk,w(:,:,:), kps(3), Ef(2)
      integer, intent(in) :: spin, kp
      
      character(1),parameter :: xyz(3)=["x","y","z"]
      character(256) :: file_ds2
      integer :: phi_n_fn(2), grad_phi_m_fn(2), atom_grad_fn, ierr
      integer :: i1, i2, ii, nn, n_min, n_max, id1 ,id2, maxijk(3),jj
      integer :: matMOcoeff(2)
      complex(double_cplx) :: ngradm_tmp1(3,3), ngradm_tmp2(3,3), ngradm_tmp(3,3)   
      real(double) :: ngradm_re, ngradm_im, tmp_re, tmp_im, weight_nm, fac
        



      if (.not.allocated(sigma)) then
      
        atom_grad_fn = allocate_temp_fn_on_grid(atomf) ! gradient of paos
        do id1 = 1, 3
          mat_dS2(id1) = allocate_temp_matrix(Hrange,0,sf,sf) ! right-side overlap derivatives                
          call matrix_scale(0d0, mat_dS2(id1))     
          gridfunctions(atom_grad_fn)%griddata = 0d0    
          call single_PAO_to_grad(id1, atom_grad_fn)
          call get_matrix_elements_new(inode-1, rem_bucket(atomf_atomf_rem), &
            mat_dS2(id1), atomfns,atom_grad_fn)

        end do
        call free_temp_fn_on_grid(atom_grad_fn)
      
!~         allocate(ngradm_summed(matrix_size,matrix_size), sigma(3,3))
        allocate(sigma(3,3))
        ngradm_summed=0d0        
        sigma=0d0
        do id1 = 1, 3
          file_ds2="dS2_"//xyz(id1)
          call write_mat(mat_dS2(id1),Hrange,maxijk,file_ds2)
        end do
      end if 
        
              
      matMOcoeff(1) = allocate_temp_matrix(Hrange,0,sf,sf) ! real part
      matMOcoeff(2) = allocate_temp_matrix(Hrange,0,sf,sf) ! imaginary part      
            
      write(io_lun, fmt='(A,i8,3e24.12)') "Calculate Kubo-Greenwood conductivity for k-point ",kp,kps
!~       call get_n_min_max(w(:, kp, spin), Ef(spin), 300d0, 1d-12, n_min, n_max)
!~       nn = (n_max - n_min + 1 + 1) * (n_max - n_min + 1) * 0.5
      
      nn = (matrix_size + 1) * matrix_size * 0.5
      
      ii = 0
      ngradm_tmp = 0d0
      do i1 = 1, -1 !matrix_size ! n_min, n_max
        do i2 = 1, i1
          ii = ii+1
          if (mod(ii, nint(nn * 0.01d0)).eq.0) write(io_lun,fmt='(2i8)', advance="no") nint(ii / (nn * 0.01d0) )
          weight_nm = get_weight_nm(w(i1, kp, spin), w(i2, kp, spin), Ef(spin), 300d0, 1d-3)

                    
          call matrix_scale(0d0, matMOcoeff(1))
          call matrix_scale(0d0, matMOcoeff(2))

          call build_nm(Hrange, kps, expH(:,:,spin), i1, i2, matMOcoeff)          
          
          do id1 = 1, 3    
            do id2 = id1, id1
              tmp_re = matrix_product_trace(mat_dS2(id1), matMOcoeff(1))
              tmp_im = matrix_product_trace(mat_dS2(id1), matMOcoeff(2))
              ngradm_tmp1(id1,id2) = cmplx(tmp_re, tmp_im)
            end do
          end do
    
          call matrix_scale(0d0, matMOcoeff(1))
          call matrix_scale(0d0, matMOcoeff(2))
              
          call build_nm(Hrange, kps, expH(:,:,spin), i2, i1, matMOcoeff)
              
          do id1 = 1, 3    
            do id2 = id1, id1
              tmp_re = matrix_product_trace(mat_dS2(id1), matMOcoeff(1))
              tmp_im = matrix_product_trace(mat_dS2(id1), matMOcoeff(2))
              ngradm_tmp2(id1,id2) = cmplx(tmp_re, tmp_im)
            end do
          end do   
          

    
          if (i1.eq.i2) then
            fac=1d0
          else
            fac=2d0
          end if
          
          ngradm_tmp = ngradm_tmp + ngradm_tmp1 * ngradm_tmp2 * weight_nm * fac
                        
!~               ngradm_summed(i1, i2) = ngradm_summed(i1,i2) + ngradm_tmp 
!~               ngradm_summed(i2, i1) = ngradm_summed(i1,i2) 
    
        end do
      end do      
      
      sigma = sigma + ngradm_tmp * wtk
      
      write(io_lun,*)
      
      do id1 = 1, 3
        do id2 = id1, id1
          write(io_lun,fmt='(1A,X,1A,X,2e24.12)' ) xyz(id1),xyz(id2), ngradm_tmp(id1,id2)
        end do
      end do
      write(io_lun, fmt='(A)') " finished" 
      
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      call free_temp_matrix(matMOcoeff(2))
      call free_temp_matrix(matMOcoeff(1))
      
    end subroutine get_conductivity

    subroutine get_conductivity_direct(w, expH, wtk, kp, kps, spin)
    
      use numbers
      use datatypes
      use GenComms, only: cq_abort, inode, ionode, myid,io_lun, gsum
      use PAO_grid_transform_module, only: single_PAO_to_grid, single_PAO_to_grad      
      use matrix_data, only: Hrange
      use global_module, only: sf, atomf
      use mult_module, only: SF_to_AtomF_transform, matrix_scale,allocate_temp_matrix, &
       free_temp_matrix,matH  
      use calc_matrix_elements_module, only: act_on_vectors_new
      use functions_on_grid, only: atomfns, allocate_temp_fn_on_grid, free_temp_fn_on_grid, &
        gridfunctions, fn_on_grid       
      use set_bucket_module, only: rem_bucket, atomf_H_atomf_rem 
      use ScalapackFormat, only: matrix_size
      use calc_matrix_elements_module, only: get_matrix_elements_new
      use set_bucket_module,           only: rem_bucket, atomf_H_atomf_rem ,atomf_atomf_rem      
      
      implicit none
    
      complex(double_cplx), intent(in) :: expH(:,:,:)
      real(double), intent(in) :: wtk,w(:,:,:), kps(3)
      integer, intent(in) :: spin, kp
      
      integer :: matMOcoeff(2), phi_n_fn(2), grad_phi_m_fn(2), atom_grad_fn,i1,i2,mat_S,ijk(3)
      complex(double_cplx) :: ngradm_tmp1,ngradm_tmp2,ngradm_tmp      
    
      if (.not.allocated(ngradm_summed)) then
        allocate(ngradm_summed(matrix_size,matrix_size))
        ngradm_summed=0d0        
      end if 
    
      matMOcoeff(1) = allocate_temp_matrix(Hrange,0,sf,sf) ! real part
      matMOcoeff(2) = allocate_temp_matrix(Hrange,0,sf,sf) ! imaginary part

      
      grad_phi_m_fn(1) = allocate_temp_fn_on_grid(atomf)
      grad_phi_m_fn(2) = allocate_temp_fn_on_grid(atomf)
   
      atom_grad_fn = allocate_temp_fn_on_grid(atomf)
      gridfunctions(atom_grad_fn)%griddata = 0d0
      call single_PAO_to_grad(3, atom_grad_fn)
      
      
      do i1=1, matrix_size
        do i2=1, matrix_size
        write(io_lun,fmt='(i8)', advance="no" ) matrix_size*matrix_size-((i2-1)*matrix_size+i1)
        
        gridfunctions(grad_phi_m_fn(1))%griddata = 0d0
        gridfunctions(grad_phi_m_fn(2))%griddata = 0d0
        call matrix_scale(0d0,matMOcoeff(1))
        call matrix_scale(0d0,matMOcoeff(2))
      
        call build_nm(Hrange, kps, expH(:,:,spin),i1,i2,matMOcoeff)
        
        call act_on_vectors_new(inode-1, rem_bucket(atomf_H_atomf_rem), &
          matMOcoeff(1), grad_phi_m_fn(1), atom_grad_fn)
        call act_on_vectors_new(inode-1, rem_bucket(atomf_H_atomf_rem), &
          matMOcoeff(2), grad_phi_m_fn(2), atom_grad_fn)
      
        call get_n_grad_m(ngradm_tmp1, grad_phi_m_fn, atomfns)
      
        gridfunctions(grad_phi_m_fn(1))%griddata = 0d0
        gridfunctions(grad_phi_m_fn(2))%griddata = 0d0
        call matrix_scale(0d0,matMOcoeff(1))
        call matrix_scale(0d0,matMOcoeff(2))      
        
         call build_nm(Hrange, kps, expH(:,:,spin),i2,i1,matMOcoeff)
          
        call act_on_vectors_new(inode-1, rem_bucket(atomf_H_atomf_rem), &
            matMOcoeff(1), grad_phi_m_fn(1), atom_grad_fn)
        call act_on_vectors_new(inode-1, rem_bucket(atomf_H_atomf_rem), &
            matMOcoeff(2), grad_phi_m_fn(2), atom_grad_fn)
        
        call get_n_grad_m(ngradm_tmp2, grad_phi_m_fn, atomfns)
        
        ngradm_tmp = ngradm_tmp1*ngradm_tmp2
        ngradm_summed(i1,i2) = ngradm_summed(i1,i2) + ngradm_tmp * wtk  
      end do
    end do
    write(io_lun,*) " finished"
      
!~       write(io_lun,fmt='(7e16.8,i6)') ngradm_tmp,ngradm_tmp1,ngradm_tmp2, wtk, kp
!~       write(io_lun,*) ngradm
      
      
      call free_temp_fn_on_grid(grad_phi_m_fn(2))
      call free_temp_fn_on_grid(grad_phi_m_fn(1))

      call free_temp_matrix(matMOcoeff(2))
      call free_temp_matrix(matMOcoeff(1))
      
    end subroutine get_conductivity_direct

    subroutine get_mo_matrix(mat_range, localEig,iband,mat_coeff,kps,fac)

      !use maxima_module, only: mx_nponn, mx_at_prim
      use numbers
      use matrix_module,   only: matrix, matrix_halo
      use group_module,    only: parts
      use primary_module,  only: bundle
      use cover_module,    only: BCS_parts
      use ScalapackFormat, only: CC_to_SC,maxrow,maxcol,proc_block,    &
           SC_to_refx,SC_to_refy, block_size_r,  &
           block_size_c, blocks_c, proc_start,   &
           matrix_size
      use global_module,   only: numprocs, iprint_DM, id_glob,         &
           ni_in_cell, x_atom_cell, y_atom_cell, &
           z_atom_cell, max_wf, out_wf,io_lun
      use mpi
      use GenBlas,         only: dot
      use GenComms,        only: myid,inode,cq_abort
      use mult_module,     only: store_matrix_value_pos, matrix_pos,matrix_index
      use matrix_data,     only: mat, halo
      use species_module,  only: nsf_species
      

      implicit none

  type Krecv_data
     integer                               :: orbs
     integer,      dimension(:),   pointer :: ints
     !integer,     dimension(:),   pointer :: ndimi
     integer,      dimension(:),   pointer :: ndimj
     integer,      dimension(:,:), pointer :: prim_atom
     integer,      dimension(:,:), pointer :: locj
     real(double), dimension(:,:), pointer :: dx, dy, dz
  end type Krecv_data

      ! Passed variables


      integer :: mat_range, mat_coeff(2)
      real(double_cplx) :: fac,kps(3)
      complex(double_cplx), dimension(:,:), intent(in) :: localEig
     
      
      integer :: part, memb, neigh, ist, prim_atom, owning_proc, locatom
      integer :: Col_FSC_part, Col_FSC_seq,  FSC_atom
      integer :: i, j, k, proc
          
      integer :: maxloc, maxint, maxsend, curr, gcspart, CC, orb_count

      integer :: ierr, atom, inter, prim,  col_sup, send_proc, recv_proc
      integer, dimension(:,:), allocatable :: atom_list

      integer, dimension(:), allocatable :: current_loc_atoms, &
           LocalAtom, prim_orbs
      integer, dimension(MPI_STATUS_SIZE) :: mpi_stat

      complex(double_cplx) :: zsum, zphase
      real(double) :: dx, dy, dz, phase

      logical :: flag, flag_write_out
      integer :: FSCpart, ipart, iband, iwf, len_occ,sf1,ah,wheremat2
      type(Krecv_data), dimension(:), allocatable :: recv_info
      complex(double_cplx), dimension(:,:), allocatable :: RecvBuffer, &
         SendBuffer
      integer, dimension(:,:), allocatable :: ints, &
         send_prim, send_info, send_orbs, send_off         


    send_proc = myid
    recv_proc = myid
    maxloc = 0
    maxint = 0
    maxsend = 0
    do i=1,numprocs
       if(current_loc_atoms(i)>maxloc) maxloc = current_loc_atoms(i)
    end do
   ! Allocate recv_info
    allocate( &
         prim_orbs(bundle%mx_iprim),STAT=ierr)
    if(ierr/=0) call cq_abort('buildK: Error allocating send_info !',ierr)
    send_info = 0
    send_orbs = 0
    send_off = 0
    prim_orbs = 0
    orb_count = 0
    do j=1,bundle%n_prim
       prim_orbs(j) = orb_count
       orb_count = orb_count + nsf_species(bundle%species(j))
    end do
    allocate(recv_info(numprocs),STAT=ierr)
    if(ierr/=0) call cq_abort('buildK: Error allocating recv_info !',ierr)
    do i=1,numprocs
       ! Build the list of which primary set atoms we send to which processor
       curr = 0
       orb_count = 0
       do j=1,bundle%n_prim
          orb_count = orb_count + nsf_species(bundle%species(j))
       end do
       allocate(recv_info(i)%ints(maxloc),recv_info(i)%ndimj(maxloc), &
            recv_info(i)%prim_atom(maxint,maxloc), recv_info(i)%locj(maxint,maxloc), &
            recv_info(i)%dx(maxint,maxloc),recv_info(i)%dy(maxint,maxloc), recv_info(i)%dz(maxint,maxloc),STAT=ierr)
       recv_info(i)%orbs = 0
       recv_info(i)%ints = 0
       !recv_info(i)%ndimi = 0
       recv_info(i)%ndimj = 0
       recv_info(i)%prim_atom = 0
       recv_info(i)%locj = 0
       recv_info(i)%dx = 0.0_double
       recv_info(i)%dy = 0.0_double
       recv_info(i)%dz = 0.0_double
       if(ierr/=0) call cq_abort('buildK: Error allocating recv_info !',ierr)
    end do

    do part = 1,bundle%groups_on_node ! Loop over primary set partitions

       if(bundle%nm_nodgroup(part)>0) then ! If there are atoms in partition
          CC = parts%ngnode(parts%inode_beg(myid+1)+part-1)
          do memb = 1,bundle%nm_nodgroup(part) ! Loop over atoms

             prim_atom = bundle%nm_nodbeg(part)+memb-1
             do neigh = 1, mat(part,mat_range)%n_nab(memb) ! Loop over neighbours of atom

                ist = mat(part,mat_range)%i_acc(memb)+neigh-1
                ! Establish FSC number of neighbour
                Col_FSC_part = BCS_parts%lab_cell(mat(part,mat_range)%i_part(ist))
                Col_FSC_seq  = mat(part,mat_range)%i_seq(ist)
                ! Now use i_cc2node(Col_FSC_part) to establish which processor owns this atom
                owning_proc = parts%i_cc2node(Col_FSC_part)
                ! Find Fundamental Simulation Cell atom
                FSC_atom = id_glob(parts%icell_beg(Col_FSC_part)+Col_FSC_seq-1)
                ! Work out a map from primary atom + FSC + identifier to distance and position in data_Matrix
                locatom = LocalAtom(FSC_atom) ! Which atom in the list on the remote proc is this ?


                recv_info(owning_proc)%ints(locatom) = recv_info(owning_proc)%ints(locatom) + 1

                gcspart = BCS_parts%icover_ibeg(mat(part,mat_range)%i_part(ist))+mat(part,mat_range)%i_seq(ist)-1
                !recv_info(owning_proc)%ndimi(locatom) = mat(part,range)%ndimi(memb)
                recv_info(owning_proc)%ndimj(locatom) = mat(part,mat_range)%ndimj(ist)
                if(recv_info(owning_proc)%ints(locatom)==1) &
                     recv_info(owning_proc)%orbs = recv_info(owning_proc)%orbs + mat(part,mat_range)%ndimj(ist)
                recv_info(owning_proc)%prim_atom(recv_info(owning_proc)%ints(locatom),locatom) = prim_atom
                recv_info(owning_proc)%locj(recv_info(owning_proc)%ints(locatom),locatom) = halo(mat_range)%i_halo(gcspart)
                ! Build the distances between atoms - needed for phases
                FSCpart = BCS_parts%lab_cell(mat(part,mat_range)%i_part(ist))!gcspart)
                recv_info(owning_proc)%dx(recv_info(owning_proc)%ints(locatom),locatom) = &
                     BCS_parts%xcover(gcspart)-bundle%xprim(prim_atom)
                recv_info(owning_proc)%dy(recv_info(owning_proc)%ints(locatom),locatom) = &
                     BCS_parts%ycover(gcspart)-bundle%yprim(prim_atom)
                recv_info(owning_proc)%dz(recv_info(owning_proc)%ints(locatom),locatom) = &
                     BCS_parts%zcover(gcspart)-bundle%zprim(prim_atom)
             end do ! End do neigh=1,mat%n_nab
          end do ! End do memb =1,nm_nodgroup
       end if ! End if nm_nodgroup > 0
    end do ! End do part=1,groups_on_node    
    
      do part = 1,bundle%groups_on_node ! Loop over primary set partitions
         if(bundle%nm_nodgroup(part)>0) then ! If there are atoms in partition
            CC = parts%ngnode(parts%inode_beg(myid+1)+part-1)
            do memb = 1,bundle%nm_nodgroup(part) ! Loop over atoms
               prim_atom = bundle%nm_nodbeg(part)+memb-1

                  neigh = 1 ! only need neighbor 1 which should be the atom itself
               
                  ist = mat(part,mat_range)%i_acc(memb)+neigh-1

                  Col_FSC_part = BCS_parts%lab_cell(mat(part,mat_range)%i_part(ist))
                  Col_FSC_seq  = mat(part,mat_range)%i_seq(ist)

                  locatom =  memb!LocalAtom(FSC_atom) ! Which atom in the list on the remote proc is this ?
                  
                  gcspart = BCS_parts%icover_ibeg(mat(part,mat_range)%i_part(ist))+mat(part,mat_range)%i_seq(ist)-1
                phase = kps(1)*recv_info(recv_proc+1)%dx(inter,locatom) + kps(2)*recv_info(recv_proc+1)%dy(inter,locatom) + &
                     kps(3)*recv_info(recv_proc+1)%dz(inter,locatom)                
                  dx=kps(1)*recv_info(recv_proc+1)%dx(inter,locatom)
                  dy=kps(1)*recv_info(recv_proc+1)%dy(inter,locatom)
                  dz=kps(1)*recv_info(recv_proc+1)%dz(inter,locatom)
                  
                  phase = kps(1)*dx + kps(2)*dy + kps(3)*dz
                  if (fac.eq.0)  then
                    zphase=1d0
                  else
                    zphase=zexp(dcmplx(0d0,fac*phase))
                  end if
!~                   write(0,fmt='(9e24.12)') zphase,phase,dx,dy,dz,kps(1),kps(2),kps(3)
                  Ah = matrix_index(mat_coeff(1))
                  sf1 = mat(part,Ah)%ndimi(locatom)
                  do col_sup = 1, sf1
                    wheremat2 = (col_sup - 1) * sf1 + col_sup - 1 + mat(part,Ah)%onsite(locatom) 
                    if (fac.eq.0) then   
                      zsum = conjg(localEig(iband,prim_orbs(prim_atom)+col_sup))
                    else
                      zsum = (localEig(iband,prim_orbs(prim_atom)+col_sup))*zphase
                    end if
!~                     write(0,*) "zsum ",zsum
                    call store_matrix_value_pos(mat_coeff(1),whereMat2,real(zsum,double))  ! assume here real expansion coefficients, i.e. gamma point                    
                    call store_matrix_value_pos(mat_coeff(2),whereMat2,aimag(zsum))  ! assume here real expansion coefficients, i.e. gamma point                    
                  end do  ! End loop over NSF                                                                        
            end do ! End do memb =1,nm_nodgroup
         end if ! End if nm_nodgroup > 0
      end do ! End do part=1,groups_on_node
      
      deallocate(prim_orbs,STAT=ierr)
      if(ierr.ne.0) call cq_abort('get_mo_matrix: Error allocating send_info !',ierr)

    end subroutine get_mo_matrix
    
  subroutine build_nm(range, kps, localEig,iband,jband,matBand)

    !use maxima_module, only: mx_nponn, mx_at_prim
    use numbers
    use matrix_module,   only: matrix, matrix_halo
    use group_module,    only: parts
    use primary_module,  only: bundle
    use cover_module,    only: BCS_parts
    use ScalapackFormat, only: CC_to_SC,maxrow,maxcol,proc_block,    &
         SC_to_refx,SC_to_refy, block_size_r,  &
         block_size_c, blocks_c, proc_start,   &
         matrix_size
    use global_module,   only: numprocs, iprint_DM, id_glob,         &
         ni_in_cell, x_atom_cell, y_atom_cell, &
         z_atom_cell, max_wf, out_wf
    use mpi
    use GenBlas,         only: dot
    use GenComms,        only: myid, io_lun, cq_abort
    use mult_module,     only: store_matrix_value_pos, matrix_pos
    use matrix_data,     only: mat, halo
    use species_module,  only: nsf_species

    implicit none


  type Krecv_data
     integer                               :: orbs
     integer,      dimension(:),   pointer :: ints
     !integer,     dimension(:),   pointer :: ndimi
     integer,      dimension(:),   pointer :: ndimj
     integer,      dimension(:,:), pointer :: prim_atom
     integer,      dimension(:,:), pointer :: locj
     real(double), dimension(:,:), pointer :: dx, dy, dz
  end type Krecv_data

    ! Passed variables
    integer :: range
    real(double) :: kps(3)
    complex(double_cplx), dimension(:,:), intent(in) :: localEig
    integer :: matBand(2), iband, jband
    ! Local variables
    type(Krecv_data), dimension(:), allocatable :: recv_info
    integer :: part, memb, neigh, ist, prim_atom, owning_proc, locatom
    integer :: Row_FSC_part, Row_FSC_seq, Col_FSC_part, Col_FSC_seq, &
         FSC_atom
    integer :: SCblockr, SCblockc, SCrowc, i, j, k, proc, stat, &
         supfn_r, supfn_c
    integer :: maxloc, maxint, maxsend, curr, gcspart, CC, orb_count
    integer :: len, send_size, recv_size, send_proc, recv_proc, nsf1, &
         sendtag, recvtag
    integer :: req1, req2, ierr, atom, inter, prim, wheremat, row_sup,&
         col_sup
    integer, dimension(:,:), allocatable :: ints, atom_list, &
         send_prim, send_info, send_orbs, send_off
    integer, dimension(:), allocatable :: current_loc_atoms, &
         LocalAtom, num_send, norb_send, send_FSC, recv_to_FSC, &
         mapchunk, prim_orbs
    integer, dimension(MPI_STATUS_SIZE) :: mpi_stat
    real(double) :: phase, rfac, ifac, rcc, icc, rsum,zphase
    complex(double_cplx) :: zsum
    complex(double_cplx), dimension(:,:), allocatable :: RecvBuffer, &
         SendBuffer
    logical :: flag, flag_write_out
    integer :: FSCpart, ipart, iwf, len_occ
    ! for spin polarisation
    real(double) :: occ_correction

    

    ! get occ_correction
    occ_correction = one

    ! Allocate data and zero arrays
    allocate(ints(numprocs,bundle%mx_iprim),&
         current_loc_atoms(numprocs),atom_list(numprocs,bundle%&
         mx_iprim), LocalAtom(ni_in_cell),send_prim(numprocs,bundle%&
         mx_iprim), num_send(numprocs),norb_send(numprocs),STAT=stat)
    if(stat/=0) call cq_abort('buildK: Error allocating ints, &
         &current_loc_atoms and atom_list !',stat)
    ints = 0
    current_loc_atoms = 0
    atom_list = 0
    LocalAtom = 0
    send_prim = 0
    num_send = 0
    norb_send = 0
    
    ! Step one - work out which processors we need to exchange data with
    do part = 1,bundle%groups_on_node ! Loop over primary set partitions
       
       if(bundle%nm_nodgroup(part)>0) then ! If there are atoms in partition
          CC = parts%ngnode(parts%inode_beg(myid+1)+part-1)
          do memb = 1,bundle%nm_nodgroup(part) ! Loop over atoms
             
             prim_atom = bundle%nm_nodbeg(part)+memb-1
             do neigh = 1, mat(part,range)%n_nab(memb) ! Loop over neighbours of atom
                if(iprint_DM>=5.AND.myid==0) write(io_lun,3) myid,neigh
                ist = mat(part,range)%i_acc(memb)+neigh-1
                ! Establish FSC number of neighbour
                Col_FSC_part = BCS_parts%lab_cell(mat(part,range)%i_part(ist))
                Col_FSC_seq  = mat(part,range)%i_seq(ist)
                ! Now use i_cc2node(Col_FSC_part) to establish which processor owns this atom
                owning_proc = parts%i_cc2node(Col_FSC_part)
                ! Find Fundamental Simulation Cell atom
                FSC_atom = id_glob(parts%icell_beg(Col_FSC_part)+Col_FSC_seq-1)
                ! Find if we have seen this before
                flag = .false.
                if(iprint_DM>=5.AND.myid==0) write(io_lun,*) 'prim, neigh, FSC: ',prim_atom, neigh, FSC_atom
                if(iprint_DM>=5.AND.myid==0) write(io_lun,*) 'curr_loc_atoms: ',current_loc_atoms(owning_proc)
                if(current_loc_atoms(owning_proc)>0) then
                   do i=1,current_loc_atoms(owning_proc)
                      if(atom_list(owning_proc,i)==FSC_atom) then
                         if(iprint_DM>=5.AND.myid==0) write(io_lun,*) 'Loc atom: ',i, LocalAtom(FSC_atom)
                         ints(owning_proc,LocalAtom(FSC_atom)) = ints(owning_proc,LocalAtom(FSC_atom)) + 1
                         send_prim(owning_proc,prim_atom) = nsf_species(bundle%species(prim_atom))
                         flag = .true.
                         exit
                      end if
                   end do
                end if ! current_loc_atoms(owning_proc)>0
                if(flag) then
                   cycle
                end if
                ! Record
                current_loc_atoms(owning_proc) = current_loc_atoms(owning_proc) + 1
                atom_list(owning_proc,current_loc_atoms(owning_proc)) = FSC_atom
                LocalAtom(FSC_atom) = current_loc_atoms(owning_proc)
                ints(owning_proc,LocalAtom(FSC_atom)) = ints(owning_proc,LocalAtom(FSC_atom)) + 1
                send_prim(owning_proc,prim_atom) = nsf_species(bundle%species(prim_atom))
             end do ! End do neigh=1,mat%n_nab
          end do ! End do memb =1,nm_nodgroup
       end if ! End if nm_nodgroup > 0
    end do ! End do part=1,groups_on_node
    ! Find max value of current_loc_atoms and interactions
    
    if(iprint_DM>3.AND.myid==0) write(io_lun,*) 'buildK: Stage two'
    maxloc = 0
    maxint = 0
    maxsend = 0
    do i=1,numprocs
       if(iprint_DM>=4.AND.myid==0) write(io_lun,*) myid,' Curr loc atoms: ',i,current_loc_atoms(i)
       if(current_loc_atoms(i)>maxloc) maxloc = current_loc_atoms(i)
       do j=1,bundle%mx_iprim ! Needs to be mx_iprim because goes over primary atoms on REMOTE processors
          if(ints(i,j)>maxint) maxint = ints(i,j)
          if(send_prim(i,j)>0) num_send(i) = num_send(i) + 1
          norb_send(i) = norb_send(i) + send_prim(i,j)
          if(iprint_DM>=5.AND.myid==0) write(io_lun,4) myid,j,send_prim(i,j),num_send(i)
       end do
       if(num_send(i)>maxsend) maxsend = num_send(i)
    end do
    if(iprint_DM>=4.AND.myid==0) write(io_lun,*) myid,' Maxima: ',maxloc, maxint, maxsend
    ! Allocate recv_info
    allocate(send_info(numprocs,maxsend),send_orbs(numprocs,maxsend),send_off(numprocs,maxsend), &
         prim_orbs(bundle%mx_iprim),STAT=stat)
    if(stat/=0) call cq_abort('buildK: Error allocating send_info !',stat)
    send_info = 0
    send_orbs = 0
    send_off = 0
    prim_orbs = 0
    orb_count = 0
    do j=1,bundle%n_prim
       prim_orbs(j) = orb_count
       orb_count = orb_count + nsf_species(bundle%species(j))
    end do
    allocate(recv_info(numprocs),STAT=stat)
    if(stat/=0) call cq_abort('buildK: Error allocating recv_info !',stat)
    do i=1,numprocs
       ! Build the list of which primary set atoms we send to which processor
       curr = 0
       orb_count = 0
       do j=1,bundle%n_prim
          if(send_prim(i,j)>0) then
             curr = curr+1
             send_info(i,curr)=j
             send_off(i,curr)=orb_count
             send_orbs(i,curr)=nsf_species(bundle%species(j))
          end if
          orb_count = orb_count + nsf_species(bundle%species(j))
       end do
       allocate(recv_info(i)%ints(maxloc),recv_info(i)%ndimj(maxloc), &
            recv_info(i)%prim_atom(maxint,maxloc), recv_info(i)%locj(maxint,maxloc), &
            recv_info(i)%dx(maxint,maxloc),recv_info(i)%dy(maxint,maxloc), recv_info(i)%dz(maxint,maxloc),STAT=stat)
       recv_info(i)%orbs = 0
       recv_info(i)%ints = 0
       !recv_info(i)%ndimi = 0
       recv_info(i)%ndimj = 0
       recv_info(i)%prim_atom = 0
       recv_info(i)%locj = 0
       recv_info(i)%dx = 0.0_double
       recv_info(i)%dy = 0.0_double
       recv_info(i)%dz = 0.0_double
       if(stat/=0) call cq_abort('buildK: Error allocating recv_info !',stat)
    end do
    
    do part = 1,bundle%groups_on_node ! Loop over primary set partitions
       if(iprint_DM>=5.AND.myid==0) write(io_lun,1) myid,part
       if(bundle%nm_nodgroup(part)>0) then ! If there are atoms in partition
          CC = parts%ngnode(parts%inode_beg(myid+1)+part-1)
          do memb = 1,bundle%nm_nodgroup(part) ! Loop over atoms
             if(iprint_DM>=5.AND.myid==0) write(io_lun,2) myid,memb
             prim_atom = bundle%nm_nodbeg(part)+memb-1
             do neigh = 1, mat(part,range)%n_nab(memb) ! Loop over neighbours of atom
                if(iprint_DM>=5.AND.myid==0) write(io_lun,3) myid,neigh
                ist = mat(part,range)%i_acc(memb)+neigh-1
                ! Establish FSC number of neighbour
                Col_FSC_part = BCS_parts%lab_cell(mat(part,range)%i_part(ist))
                Col_FSC_seq  = mat(part,range)%i_seq(ist)
                ! Now use i_cc2node(Col_FSC_part) to establish which processor owns this atom
                owning_proc = parts%i_cc2node(Col_FSC_part)
                ! Find Fundamental Simulation Cell atom
                FSC_atom = id_glob(parts%icell_beg(Col_FSC_part)+Col_FSC_seq-1)
                ! Work out a map from primary atom + FSC + identifier to distance and position in data_Matrix
                locatom = LocalAtom(FSC_atom) ! Which atom in the list on the remote proc is this ?
                if(iprint_DM>=5.AND.myid==0) write(io_lun,*) myid,' own, FSC, loc: ',owning_proc, FSC_atom, locatom, &
                     recv_info(owning_proc)%ints(locatom)
                recv_info(owning_proc)%ints(locatom) = recv_info(owning_proc)%ints(locatom) + 1
                if(iprint_DM>=5.AND.myid==0) write(io_lun,*) myid,' ints: ',recv_info(owning_proc)%ints(locatom)
                gcspart = BCS_parts%icover_ibeg(mat(part,range)%i_part(ist))+mat(part,range)%i_seq(ist)-1
                !recv_info(owning_proc)%ndimi(locatom) = mat(part,range)%ndimi(memb)
                recv_info(owning_proc)%ndimj(locatom) = mat(part,range)%ndimj(ist)
                if(recv_info(owning_proc)%ints(locatom)==1) &
                     recv_info(owning_proc)%orbs = recv_info(owning_proc)%orbs + mat(part,range)%ndimj(ist)
                recv_info(owning_proc)%prim_atom(recv_info(owning_proc)%ints(locatom),locatom) = prim_atom
                recv_info(owning_proc)%locj(recv_info(owning_proc)%ints(locatom),locatom) = halo(range)%i_halo(gcspart)
                ! Build the distances between atoms - needed for phases
                FSCpart = BCS_parts%lab_cell(mat(part,range)%i_part(ist))!gcspart)
                recv_info(owning_proc)%dx(recv_info(owning_proc)%ints(locatom),locatom) = &
                     BCS_parts%xcover(gcspart)-bundle%xprim(prim_atom)
                recv_info(owning_proc)%dy(recv_info(owning_proc)%ints(locatom),locatom) = &
                     BCS_parts%ycover(gcspart)-bundle%yprim(prim_atom)
                recv_info(owning_proc)%dz(recv_info(owning_proc)%ints(locatom),locatom) = &
                     BCS_parts%zcover(gcspart)-bundle%zprim(prim_atom)
!~                  recv_info(owning_proc)%dx(recv_info(owning_proc)%ints(locatom),locatom) = BCS_parts%xcover(gcspart)- &
!~                       x_atom_cell(parts%icell_beg(FSCpart)+mat(part,range)%i_seq(ist)-1)
!~                  recv_info(owning_proc)%dy(recv_info(owning_proc)%ints(locatom),locatom) = BCS_parts%ycover(gcspart)- &
!~                       y_atom_cell(parts%icell_beg(FSCpart)+mat(part,range)%i_seq(ist)-1)
!~                  recv_info(owning_proc)%dz(recv_info(owning_proc)%ints(locatom),locatom) = BCS_parts%zcover(gcspart)- &
!~                       z_atom_cell(parts%icell_beg(FSCpart)+mat(part,range)%i_seq(ist)-1)
                     
             end do ! End do neigh=1,mat%n_nab
          end do ! End do memb =1,nm_nodgroup
       end if ! End if nm_nodgroup > 0
    end do ! End do part=1,groups_on_node
    ! Work out length
    ! NB we send all bands up to last occupied one (for deltaSCF this will include some empty)
    len = 2 !matrix_size !max(iband,jband)
    
    ! Step three - loop over processors, send and recv data and build K
    allocate(send_fsc(bundle%mx_iprim),recv_to_FSC(bundle%mx_iprim),mapchunk(bundle%mx_iprim),STAT=stat)
    if(stat/=0) call cq_abort('buildK: Error allocating send_fsc, recv_to_FSC and mapchunk',stat)
    send_fsc = 0
    recv_to_FSC = 0
    mapchunk = 0
    send_proc = myid
    recv_proc = myid
    ! Two messages per process; use tag and tag + 1 below
    sendtag = 1
    recvtag = 1
    do i=1,numprocs
       send_size = len*norb_send(send_proc+1)!num_send(send_proc+1)*nsf
       recv_size = len*recv_info(recv_proc+1)%orbs!current_loc_atoms(recv_proc+1)*nsf
       if(iprint_DM>=4.AND.myid==0) write(io_lun,*) 'Send and recv sizes: ',send_size, recv_size
       ! Fill SendBuffer
       allocate(SendBuffer(len,norb_send(send_proc+1)),STAT=stat)
       if(stat/=0) call cq_abort('buildK: Unable to allocate SendBuffer !',stat)
       if(iprint_DM>=4.AND.myid==0) write(io_lun,*) 'Filling SendBuffer'
       orb_count = 0
       do j=1,num_send(send_proc+1)
          do nsf1=1,send_orbs(send_proc+1,j)
             orb_count = orb_count+1
             SendBuffer(1,orb_count) = localEig(iband,send_off(send_proc+1,j)+nsf1)
             SendBuffer(2,orb_count) = localEig(jband,send_off(send_proc+1,j)+nsf1)
          end do
          ! We also need to send a list of what FSC each primary atom sent corresponds to - use bundle%ig_prim
          send_FSC(j) = bundle%ig_prim(send_info(send_proc+1,j))
          if(iprint_DM>=4.AND.myid==0) write(io_lun,*) 'Building send_FSC: ',send_info(send_proc+1,j), &
               bundle%ig_prim(send_info(send_proc+1,j)),send_FSC(j)
       end do
       if(orb_count/=norb_send(send_proc+1)) call cq_abort("Orbital mismatch in buildK: ",orb_count,norb_send(send_proc+1))
       if(iprint_DM>=4.AND.myid==0) write(io_lun,*) 'Sending'
       ! Now send
       if(send_size>0) then
          if(send_proc/=myid) then
             call MPI_issend(send_FSC,num_send(send_proc+1),MPI_INTEGER,send_proc,sendtag,MPI_COMM_WORLD,req1,ierr)
             call MPI_issend(SendBuffer,send_size,MPI_DOUBLE_COMPLEX,send_proc,sendtag+1,MPI_COMM_WORLD,req2,ierr)
          end if
       end if
       ! Now receive data
       if(iprint_DM>=4.AND.myid==0) write(io_lun,*) 'Alloc RecvBuffer ',len,recv_info(recv_proc+1)%orbs
       !allocate(RecvBuffer(len,current_loc_atoms(recv_proc+1)*nsf),STAT=stat)
       allocate(RecvBuffer(len,recv_info(recv_proc+1)%orbs),STAT=stat)
       if(stat/=0) call cq_abort('buildK: Unable to allocate RecvBuffer !',stat)
       if(iprint_DM>=4.AND.myid==0) write(io_lun,*) 'Recving'
       if(recv_size>0) then
          if(recv_proc/=myid) then
             call MPI_recv(recv_to_FSC,current_loc_atoms(recv_proc+1),MPI_INTEGER,recv_proc,recvtag,MPI_COMM_WORLD,mpi_stat,ierr)
             if(iprint_DM>=4.AND.myid==0) write(io_lun,*) 'Got recv_to_FSC'
             call MPI_recv(RecvBuffer,recv_size,MPI_DOUBLE_COMPLEX,&
                  recv_proc,recvtag+1,MPI_COMM_WORLD,mpi_stat,ierr)
             if(iprint_DM>=4.AND.myid==0) write(io_lun,*) 'Got RecvBuffer'
          else
             if(iprint_DM>=4.AND.myid==0) write(io_lun,*) 'On-proc: getting recv_to_FSC'
             recv_to_FSC(1:current_loc_atoms(recv_proc+1)) = send_FSC(1:current_loc_atoms(recv_proc+1))
             if(iprint_DM>=4.AND.myid==0) write(io_lun,*) 'On-proc: getting RecvBuffer'
             RecvBuffer(1:len,1:recv_info(recv_proc+1)%orbs) = SendBuffer(1:len,1:recv_info(recv_proc+1)%orbs)
          end if
          if(iprint_DM>=4.AND.myid==0) write(io_lun,*) 'Doing the mapchunk', recv_to_FSC
          do j=1,current_loc_atoms(recv_proc+1)
             mapchunk(j) = LocalAtom(recv_to_FSC(j))
          end do
          if(iprint_DM>=4.AND.myid==0) write(io_lun,*) 'filling buffer'
          orb_count = 0
          do atom = 1,current_loc_atoms(recv_proc+1)
             locatom = mapchunk(atom)
             if(iprint_DM>=4.AND.myid==0) write(io_lun,*) 'Atom, loc: ',atom,locatom,recv_info(recv_proc+1)%ints(locatom)
             ! Scale the eigenvector coefficients we've received
             ! The factor of 0.5 is because the occupation numbers are from 0->2 (we expect 0->1 in K)
             ! The occupation numbers contain the k-point weight
             ! When we're doing M12 not K then occ also contains the eigenvalue
             !do col_sup = 1,nsf
             !   do j=1,len
             !      RecvBuffer(j,(atom-1)*nsf+col_sup) = RecvBuffer(j,(atom-1)*nsf+col_sup)*0.5_double*occs(j)
             !   end do
             !end do
             ! N.B. the routine used for dot is zdotc which takes the complex conjugate of the first vector
             do inter = 1,recv_info(recv_proc+1)%ints(locatom)
                prim = recv_info(recv_proc+1)%prim_atom(inter,locatom)
                if(iprint_DM>=4.AND.myid==0) write(io_lun,*) 'Inter: ',inter,prim
                phase = kps(1)*recv_info(recv_proc+1)%dx(inter,locatom) + kps(2)*recv_info(recv_proc+1)%dy(inter,locatom) + &
                     kps(3)*recv_info(recv_proc+1)%dz(inter,locatom)
                if(iprint_DM>=5.AND.myid==0) write(io_lun,*) 'Prim, where, phase: ',prim, whereMat, phase
                rfac = cos(phase)
                ifac = sin(phase)
                zphase=zexp(dcmplx(0d0,phase))
                do row_sup = 1,recv_info(recv_proc+1)%ndimj(locatom)
                   do col_sup = 1,nsf_species(bundle%species(prim))
                      whereMat = matrix_pos(matBand(1),recv_info(recv_proc+1)%prim_atom(inter,locatom), &
                           recv_info(recv_proc+1)%locj(inter,locatom),col_sup,row_sup)                                            
                            
                            zsum = conjg(localEig(jband,prim_orbs(prim)+col_sup))*RecvBuffer(1,orb_count+row_sup)*zexp(dcmplx(0d0,phase))
!~                             write(0,fmt='(3i5,9e16.8)') orb_count+row_sup,prim_orbs(prim)+col_sup,wheremat,conjg(localEig(jband,prim_orbs(prim)+col_sup)),RecvBuffer(iband,orb_count+row_sup),zexp(dcmplx(0d0,phase)),recv_info(recv_proc+1)%dz(inter,locatom)
                            call store_matrix_value_pos(matBand(1),whereMat,real(zsum))
                            call store_matrix_value_pos(matBand(2),whereMat,aimag(zsum))                                               
                   end do ! col_sup=nsf
                end do ! row_sup=nsf
             end do ! inter=recv_info%ints
             ! Careful - we only want to increment after ALL interactions done
             orb_count = orb_count + recv_info(recv_proc+1)%ndimj(locatom)
          end do ! atom=current_loc_atoms
       end if ! recv_size>0
       if(iprint_DM>=4.AND.myid==0) write(io_lun,*) 'Calling MPI_Wait'
       if(send_size>0.AND.myid/=send_proc) then
          call MPI_Wait(req1,mpi_stat,ierr)
          call MPI_Wait(req2,mpi_stat,ierr)
       end if
       if(iprint_DM>=4.AND.myid==0) write(io_lun,*) 'Calling dealloc'
       deallocate(RecvBuffer,STAT=stat)
       if(stat/=0) call cq_abort("buildK: Failed to dealloc buffer",stat)
       if(iprint_DM>=4.AND.myid==0) write(io_lun,*) 'Calling dealloc'
       deallocate(SendBuffer,STAT=stat)
       if(stat/=0) call cq_abort("buildK: Failed to dealloc buffer",stat)
       ! Increment/decrement recv and send, and wrap
       ! Remember that we go from 0->numprocs-1
       if(iprint_DM>=4.AND.myid==0) write(io_lun,*) 'Doing proc thang'
       send_proc = send_proc +1
       if(send_proc.GT.numprocs-1) send_proc = 0
       recv_proc = recv_proc -1
       if(recv_proc.LT.0) recv_proc = numprocs-1
    end do ! do i=numprocs
    ! Now deallocate all arrays
    deallocate(send_fsc,recv_to_FSC,mapchunk,STAT=stat)
    if(stat/=0) call cq_abort('buildK: Error deallocating send_fsc, recv_to_FSC and mapchunk !',stat)
    do i=numprocs,1,-1
       deallocate(recv_info(i)%ints,recv_info(i)%prim_atom,recv_info(i)%locj,&
            recv_info(i)%dx,recv_info(i)%dy,recv_info(i)%dz,recv_info(i)%ndimj,STAT=stat)
       if(stat/=0) call cq_abort('buildK: Error deallocating recvinfo !',i,stat)
    end do
    deallocate(prim_orbs,send_off,send_orbs,send_info,STAT=stat)
    if(stat/=0) call cq_abort('buildK: Error allocating send_info !',stat)
    deallocate(recv_info,STAT=stat)
    if(stat/=0) call cq_abort('buildK: Error allocating recv_info !',stat)
    deallocate(ints,current_loc_atoms,atom_list,LocalAtom,send_prim,num_send,norb_send,STAT=stat)
    if(stat/=0) call cq_abort('buildK: Error allocating ints etc !',stat)
    
    return
1   format(10x,'Processor: ',i5,' Partition: ',i5)
2   format(10x,'Processor: ',i5,' Atom: ',i5)
3   format(10x,'Processor: ',i5,' Neighbour: ',i5)
4   format(10x,'Proc: ',i5,' Prim, send_prim, num_send: ',3i5)
  end subroutine build_nm  
    
    
    subroutine get_n_grad_m(ngradm, atom_fns_n, atom_fns_grad_m)
  
      use datatypes
      use numbers
      use global_module, only: io_lun
      use dimens,                      only: n_my_grid_points, grid_point_volume
      use block_module,                only: n_pts_in_block
      use primary_module,              only: domain
      use set_blipgrid_module,         only: naba_atoms_of_blocks
      use GenComms,                    only: gsum, inode, ionode
      use global_module,               only: iprint_SC, atomf, ni_in_cell
      use functions_on_grid,           only: gridfunctions, fn_on_grid    
      use maxima_module,   ONLY: maxngrid   
      
      implicit none
  
      ! output
      
      complex(double_cplx) :: ngradm
  
      ! Passed variables
  
      integer, intent(in) :: atom_fns_n(2), atom_fns_grad_m

      ! Local variables
      integer :: blk, i_count_alpha, n, n_i, n_point

      real(double) :: a1,b1,a2,b2
      complex(double_cplx), allocatable :: tmp1(:),tmp2(:)
      complex(double_cplx) :: ii=cmplx(0d0,1d0)
      
      complex(double_cplx), external :: zdotc
  
      if (inode == ionode .and. iprint_SC >= 2) &
          write (io_lun,fmt='(2x,"Entering get_n_grad_m")')
  
      ngradm = 0d0
      allocate(tmp1(maxngrid),tmp2(maxngrid))
      tmp1=0d0
      tmp2=0d0
      i_count_alpha = 0

      do blk = 1, domain%groups_on_node
        n_point = (blk - 1) * n_pts_in_block
        i_count_alpha = (naba_atoms_of_blocks(atomf)%ibegin_blk_orb(blk)-1) * n_pts_in_block
        if (naba_atoms_of_blocks(atomf)%no_of_atom(blk) > 0) then   !TM 30/Jun/2003
            do n_i = 1, naba_atoms_of_blocks(atomf)%no_of_orb(blk)*n_pts_in_block, n_pts_in_block
              do n=1, n_pts_in_block
                  
                  a1 = gridfunctions(atom_fns_n(1))%griddata(i_count_alpha+n_i+n-1)
                  b1 = gridfunctions(atom_fns_n(2))%griddata(i_count_alpha+n_i+n-1)
                  a2 = gridfunctions(atom_fns_grad_m)%griddata(i_count_alpha+n_i+n-1)                
!~                   b2 = gridfunctions(atom_fns_grad_m(2))%griddata(i_count_alpha+n_i+n-1)
                  tmp1(n_point+n)=tmp1(n_point+n)+(a1+ii*b1)*a2
!~                   tmp2(n_point+n)=tmp2(n_point+n)+a2
                            
              end do
            end do
        end if !(naba_atoms_of_blocks(atomf)%no_of_atom(blk) > 0)   !TM 30/Jun/2003
      end do ! blk
      
!~       ngradm = zdotc(n_my_grid_points,tmp1,1,tmp2,1)
!~       write(io_lun,*) n_my_grid_points
!~       write(0,fmt='(2e24.12)') conjg(tmp1(1:n_my_grid_points))*tmp2(1:n_my_grid_points)
!~       ngradm = sum(conjg(tmp1(1:n_my_grid_points))*tmp2(1:n_my_grid_points))
!~       ngradm = sum((tmp1(1:n_my_grid_points))*tmp2(1:n_my_grid_points))
      ngradm = sum((tmp1(1:n_my_grid_points)))
      call gsum(ngradm)
      ngradm = grid_point_volume * ngradm

    end subroutine get_n_grad_m

    function dfdE(x, Ef, T)
      implicit none
      
      real(double) :: dfdE
      real(double), intent(in) :: x, Ef, T
      
      real(double) :: beta
      
      beta = 1d0/(T * kb)      
      dfdE = -beta / ( exp((x - Ef) * beta * 0.5d0) + exp(-(x - Ef) * beta * 0.5d0))**2
            
    end function dfdE
    
    subroutine get_n_min_max(w,Ef,T,eps,n_min,n_max)
      
      use global_module, only: io_lun
      
      implicit none
    
      real(double) :: w(:),Ef,T,eps
      integer :: n_min,n_max
      
      integer :: n,i
      
      n=size(w)
      
      do n_min=1,n
        if (abs(dfdE(w(n_min), Ef, T)).ge.eps) exit
      end do
      
      do n_max=n,1,-1
        if (abs(dfdE(w(n_max), Ef, T)).ge.eps) exit
      end do
      
      write(io_lun,fmt='(A,e16.8,A,2i8)') "For T =", T, " K consider states between ",n_min, n_max
    
    end subroutine get_n_min_max
    
!~ to get weight we need integral_-inf,inf dE df/dE*(eta/(E-E_n)^2+eta^2)*(eta/(E-E_m)^2+eta^2)
!~ df/dE -> delta(E-Ef) for T->0 but lets keep the finite temperature broading for now.
    function get_weight_nm(En, Em, Ef, T, eta)
      use GenComms, only: io_lun
      implicit none
      
      real(double) :: get_weight_nm
      real(double), intent(in) :: En, Em, Ef, T, eta
      
      real(double) :: f1, dx, eps, x, f, xmin , xmax
      integer :: ix, nmax, ii
      
      
    
      
      dx = min(abs(En-Ef),abs(Em-Ef))*1d-3   
      xmin = min(Ef,min(En,Em))-max(kb*T,eta)*100
      xmax = max(Ef,min(En,Em))+max(kb*T,eta)*100      
      nmax = nint((xmax-xmin)/dx) + 1
      
      f1 = 0d0      
      do ix = 2, nmax-1
        x = xmin + dx * (ix -1)
        f = dfdE(x, Ef, T) * lorentz_function(x, En, eta)* lorentz_function(x, Em, eta)
        f1 = f1 + f
      end do
      
      f1 = f1 * 2d0
      x = xmin
      f1 = f1 + dfdE(x, Ef, T) * lorentz_function(x, En, eta)* lorentz_function(x, Em, eta)
      x = xmax
      f1 = f1 + dfdE(x, Ef, T) * lorentz_function(x, En, eta)* lorentz_function(x, Em, eta)
      
      get_weight_nm = f1 * dx
      
      
!~       write(io_lun, fmt='(2e24.12,i8)') get_weight_nm,( eta*eta / ( ( (Ef-En)**2 + eta**2 )  * ( (Ef-Em)**2 + eta**2 ) ) ),nmax
      
    end function get_weight_nm
    
    function lorentz_function(x,xn,eta)
      implicit none
      
        real(double) :: lorentz_function
        real(double), intent(in) :: x, xn, eta
      
        lorentz_function = eta /  ( (x-xn)**2 + eta**2 )
    
    end function lorentz_function

end module kubo_greenwood
