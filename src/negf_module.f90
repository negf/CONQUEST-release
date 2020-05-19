module negf_module

  use datatypes
  implicit none

  contains


!~ no electron phonon coupling constant calculations for now
#ifdef EP_COUPLING  
!~   dump one sided Cartesian derivative of S, i.e. <di/dr|j> to file.
!~   they can be converted to Atomic coordinate displacements using phi(r,R)=phi(r-R) -> d/dr=-d/dR
    subroutine get_dS(mat_dS,di_or_dj,iatom)
    
      use GenComms, ONLY: inode,io_lun
      use global_module,               only: iprint_ops,atomf,numprocs
      use matrix_data,                 only: aSa_range,Hrange
      use mult_module,                 only: matSatomf,allocate_temp_matrix,free_temp_matrix
      use io_module,       only: get_file_name,dump_matrix,dump_charge
      use input_module,    only: io_assign, io_close    
      use PAO_grid_transform_module,   only: single_PAO_to_grid,single_PAO_to_grad
      use functions_on_grid,           only: atomfns,gridfunctions, fn_on_grid, &
                                            allocate_temp_fn_on_grid, &
                                            free_temp_fn_on_grid    
      use set_bucket_module,           only: rem_bucket, atomf_H_atomf_rem ,atomf_atomf_rem
      use functions_on_grid,          only: atomfns
      use calc_matrix_elements_module, only: get_matrix_elements_new
        
      implicit none
      
      integer :: mat_dS(3)
      integer :: di_or_dj,iatom
      
  
      integer :: mat_tmp,tmp_fn,maxijk(3),matrixsize,iunit,icount,direction
      character(256) :: outfile,outfile_petsc,tmpstr,atstr
      
!~       mat_tmp = allocate_temp_matrix(Hrange, 0, atomfns, atomfns)
      
      do direction=1,3
        tmp_fn = allocate_temp_fn_on_grid(atomf)
        gridfunctions(tmp_fn)%griddata = 0d0
        call single_PAO_to_grad(direction, tmp_fn,iatom)    
        if (di_or_dj.eq.1) then
          call get_matrix_elements_new(inode-1, rem_bucket(atomf_atomf_rem), &
                                        mat_dS(direction),tmp_fn , &
                                        atomfns)
        else if (di_or_dj.eq.2) then
          call get_matrix_elements_new(inode-1, rem_bucket(atomf_atomf_rem), &
                                        mat_dS(direction),atomfns, &
                                        tmp_fn)
        else
          write(0,*) "di_or_dj must be 1 or 2, i.e. 1=<i'|j>, 2=<i|j'> ",di_or_dj
        end if
!~         write(outfile_petsc,fmt='(i8)') direction
!~         write(atstr,fmt='(i8)') iatom
!~         write(tmpstr,fmt='(i8)') di_or_dj
        
!~         outfile_petsc="d"//trim(adjustl(tmpstr))//"S_"//trim(adjustl(atstr))//"_"//trim(adjustl(outfile_petsc))
!~         write(io_lun,fmt='(A)') trim(outfile_petsc)
!~         call get_file_name(trim(outfile_petsc),numprocs,inode,outfile)
!~         call io_assign(iunit)
!~         open(unit=iunit,file=outfile)
!~         maxijk=1
!~         matrixsize=0      
!~         call write_mat(mat_dS(direction),Hrange,iunit,matrixsize,icount,maxijk,outfile_petsc)
!~         call io_close(iunit)                                             
        call free_temp_fn_on_grid(tmp_fn)    
      end do
!~       call free_temp_matrix(mat_tmp)
    
    end subroutine get_dS
    
    subroutine get_dH(mat_dS,di_or_dj,iatom,outfile_iat)
    
      use GenComms, ONLY: inode,io_lun
      use global_module,               only: iprint_ops,atomf,numprocs,ni_in_cell
      use matrix_data,                 only: aSa_range,Hrange,aHa_range,Srange,Trange,TSrange,THrange
      use mult_module                
      use io_module,       only: get_file_name,dump_matrix,dump_charge
      use input_module,    only: io_assign, io_close    
      use PAO_grid_transform_module,   only: single_PAO_to_grid,single_PAO_to_grad
      use functions_on_grid,           only: atomfns,gridfunctions, fn_on_grid, &
                                            allocate_temp_fn_on_grid, &
                                            free_temp_fn_on_grid    
      use set_bucket_module,           only: rem_bucket, atomf_H_atomf_rem    ,atomf_atomf_rem
      use functions_on_grid,          only: atomfns
      use calc_matrix_elements_module, only: get_matrix_elements_new,act_on_vectors_new
      use species_module, ONLY : nsf_species, species_label,species
      
      
        
      implicit none
      
      integer :: mat_dS(3)
      integer :: di_or_dj,iatom
      character(256) :: outfile_iat
  
      integer :: tmp_fn,maxijk(3),matrixsize,iunit,icount,mat_BC,mat_ABC,ixyz,Mrange,mat_tmp,i
      character(256) :: outfile,outfile_petsc,strxyz,strat,what
      
      
      mat_BC = allocate_temp_matrix(THrange, 0, atomfns, atomfns)
      call matrix_scale(0d0, mat_BC)
      mat_ABC = allocate_temp_matrix(Hrange, 0, atomfns, atomfns)  
      call matrix_scale(0d0, mat_ABC)
!~       mat_tmp = allocate_temp_matrix(Hrange, 0, atomfns, atomfns)  
!~       call matrix_scale(0d0, mat_tmp)
  
    
      if (di_or_dj.eq.1)  then
        call matrix_product(matT(1),matH(1),mat_BC,mult(T_H_TH))
        what="d1H"          
        outfile_petsc="SinvH"
      else if (di_or_dj.eq.2)  then  
        call matrix_product(matH(1),matT(1),mat_BC,mult(H_T_TH))
        what="d2H"
        outfile_petsc="HSinv"
      end if
      
      what="dH"
  
!~       write(strat,fmt='(i8)') iatom
!~       strat=adjustl(strat)
      
!~       write(io_lun,fmt='(A)') trim(outfile_petsc)      
!~       call get_file_name(trim(outfile_petsc),numprocs,inode,outfile)
!~       call io_assign(iunit)
!~       open(unit=iunit,file=outfile)
!~       maxijk=1
!~       matrixsize=0
!~       call write_mat(mat_BC,THrange,iunit,matrixsize,icount,maxijk,outfile_petsc)
!~       call io_close(iunit)      
  
      matrixsize = 0 !ni_in_cell*nsf
      do i = 1, ni_in_cell
        matrixsize = matrixsize + nsf_species(species(i))
      end do
!~       if (inode.eq.0) write(io_lun,*) "matrixsize ",matrixsize
            
      do ixyz=1,3  
        if (di_or_dj.eq.1)  then
          call matrix_product(mat_dS(ixyz),mat_BC,mat_ABC,mult(S_TH_H))        
        else if (di_or_dj.eq.2) then
          call matrix_product(mat_BC,mat_dS(ixyz),mat_ABC,mult(TH_S_H))        
        end if
        write(strxyz,fmt='(i8)') ixyz
        strxyz=adjustl(strxyz)     
        outfile_petsc=trim(adjustl(outfile_iat))//"_"//trim(strxyz)
        write(io_lun,fmt='(A)') trim(outfile_petsc)      
!~         call get_file_name(trim(outfile_petsc),numprocs,inode,outfile)
!~         call io_assign(iunit)
!~         open(unit=iunit,file=outfile)
        maxijk=1
!~         matrixsize=0
        call matrix_scale(-1d0,mat_ABC)
        call write_mat(mat_ABC,Hrange,iunit,matrixsize,icount,maxijk,outfile_petsc,.true.)
!~         call io_close(iunit)      
      end do
  
    
!~       call free_temp_matrix(mat_tmp)
      call free_temp_matrix(mat_ABC)
      call free_temp_matrix(mat_BC)
    
    end subroutine get_dH  
#endif

    subroutine dump_negf()
      use GenComms, ONLY: inode,ionode
      implicit none
      
      integer :: iunit,ic,ierr,maxijk(3)
      
      if (ionode.eq.inode) then
        open(newunit=iunit,file="sys_info.dat",action="write",status="replace")
        close(iunit)
      end if
    
      call get_mat(1,ic,maxijk)
      call get_sysinfo(1,ic,maxijk)
      call get_sysinfo(2)      
      
    end subroutine dump_negf
    
    subroutine fix_densities_negf(density_file_l,density_file_r,density,scal)
      use datatypes
      use grid_solver
      use maxima_module,    only: maxngrid
      use global_module, only: negf_l_elec_dir,negf_r_elec_dir
      implicit none

      character(*) :: density_file_l,density_file_r
      real(double) :: density(:)
      real(double) :: scal

      real(double), allocatable, target :: denstmp(:)
      real(double), pointer :: rho3(:,:,:)
      integer :: ierr
      character(256) :: file_l,file_r

      allocate(denstmp(maxngrid),stat=ierr)
      if (ierr.ne.0) then
        write(0,*) "allocation err ",ierr
        stop
      end if

      file_l=trim(negf_l_elec_dir)//"/"//trim(density_file_l)
      file_r=trim(negf_r_elec_dir)//"/"//trim(density_file_r)

      denstmp=density/scal
      call rearrange(denstmp,maxngrid,gs_packbx,gs_psendbx,gs_prbx,gs_unpackbx)
      rho3(1:nglx,1:ngly,1:nglz)=>denstmp(1:nglx*ngly*nglz)
      call fix_rho(rho3,file_l,file_r)
      call rearrange(denstmp,maxngrid,gs_packxb,gs_psendxb,gs_prxb,gs_unpackxb)
      density=denstmp*scal
      nullify(rho3)
      deallocate(denstmp,stat=ierr)
      if (ierr.ne.0) then
        write(0,*) "allocation err ",ierr
        stop
      end if

    end subroutine

    subroutine load_knegf()
      use numbers
      use datatypes
      use mult_module,      only: matphi, T_trans, L_trans, LS_trans,  &
                                  SFcoeff_trans, matK, matrix_scale
      use SelfCon,          only: new_SC_potl,q0,flag_kerker,A,flag_wdmetric,q1,&
                                  maxpulaySC,get_pulay_optimal_rho
      use global_module,    only: iprint_init,glob2node,&
                                  MDinit_step,ni_in_cell, &
                                  flag_neutral_atom,atom_coord,nspin,&
                                  atom_coord_diff, rcellx, rcelly, rcellz,&
                                  read_gless,ne_in_cell,io_lun, &
                                  flag_pcc_global, spin_factor, &
                                  negf_density_rescale
      use GenComms,         only: my_barrier,end_comms,inode,ionode,&
                                  cq_abort,gcopy,gmax,gsum
      use density_module,   only: get_electronic_density, density,density_atom,&
                                  density_pcc
      use maxima_module,    only: maxngrid
      use functions_on_grid,only: atomfns, H_on_atomfns
      use dimens,           only: n_my_grid_points,grid_point_volume
      use matrix_data,      only: Lrange,Trange,LSrange,SFcoeff_range,Hrange
      use store_matrix,     only: matrix_store_global, grab_InfoMatGlobal, grab_matrix2, &
                                  InfoMatrixFile
      use UpdateInfo,       only: make_glob2node,Matrix_CommRebuild
      use gridout_module,   only: grid_out_plt, print_grid_data
      use hartree_module,   only: kerker,kerker_and_wdmetric,wdmetric
      use GenBlas,          only: dot,rsum
      use grid_solver

      implicit none


      !H_trans is not prepared. If we need to symmetrise K, we need H_trans
      integer        :: H_trans = 1

      real(double), allocatable, target :: denstmp(:)
      real(double), pointer :: rho3(:,:,:)
      real(double), allocatable :: rho_pul(:,:,:),rho1(:,:)
      real(double), allocatable,target ::   k_resid_pul(:,:,:), resid_pul(:,:,:)
      real(double), pointer, dimension(:,:,:) :: cov_resid_pul(:,:,:)
      real(double) :: electrons(2),electrons1,electrons_tot,R0
      logical :: lexists
      integer :: ipulay,iunit,npulay,jpulay,i,niter,kmix,ig,ierr,nfile
      character(256) :: pltfile

      type(matrix_store_global) :: InfoGlob
      type(InfoMatrixFile),pointer :: Info(:)


      if (inode.eq.ionode) write (io_lun,*) "Get global info to load matrices for Knegf"
      if (inode.eq.ionode) call make_glob2node
      call gcopy(glob2node, ni_in_cell)

      call grab_InfoMatGlobal(InfoGlob,MDinit_step,lnegf=.true.)
      MDinit_step = InfoGlob%MDstep

      do ig = 1, ni_in_cell
        atom_coord_diff(1:3,ig) = atom_coord(1:3,ig) - InfoGlob%atom_coord(1:3,ig)
        if((atom_coord_diff(1,ig)) > half*rcellx) atom_coord_diff(1,ig)=atom_coord_diff(1,ig)-rcellx
        if((atom_coord_diff(1,ig)) < -half*rcellx) atom_coord_diff(1,ig)=atom_coord_diff(1,ig)+rcellx
        if((atom_coord_diff(2,ig)) > half*rcelly) atom_coord_diff(2,ig)=atom_coord_diff(2,ig)-rcelly
        if((atom_coord_diff(2,ig)) < -half*rcelly) atom_coord_diff(2,ig)=atom_coord_diff(2,ig)+rcelly
        if((atom_coord_diff(3,ig)) > half*rcellz) atom_coord_diff(3,ig)=atom_coord_diff(3,ig)-rcellz
        if((atom_coord_diff(3,ig)) < -half*rcellz) atom_coord_diff(3,ig)=atom_coord_diff(3,ig)+rcellz
      end do



      call matrix_scale(zero,matK(1))

      allocate(rho1(maxngrid,nspin),stat=ierr)
      if (ierr.ne.0) then
        write(0,*) "allocation error ",ierr
        stop
      end if
      density=0d0
      rho1=0d0

      if (inode.eq.ionode) write(io_lun,*) "read K_negf matrix..."

      allocate(denstmp(maxngrid),rho_pul(maxngrid,maxpulaySC,nspin), &
              resid_pul(maxngrid,maxpulaySC,nspin), &
              k_resid_pul(maxngrid,maxpulaySC,nspin),stat=ierr)
      if (ierr.ne.0) then
        write(0,*) "allocation error ",ierr
        stop
      end if


      if (flag_pcc_global) then
        call fix_densities_negf("density_pcc_l.dat","density_pcc_r.dat",density_pcc,1d0)
      end if

      if (flag_neutral_atom) then
        if (inode.eq.ionode) write(io_lun,*) "negf_density_rescale ",negf_density_rescale
        call fix_densities_negf("density_atom_l.dat","density_atom_r.dat",density_atom,negf_density_rescale)
      end if



      ipulay=0
      npulay=0
      if (inode.eq.ionode) then
        inquire(file="pulay.negf",exist=lexists)
        if (.not.lexists) then
          open(newunit=iunit,file="pulay.negf",status="replace")
          ipulay=1
          npulay=0
          niter=0
        else if (lexists) then
          open(newunit=iunit,file="pulay.negf",status="old")
          read(iunit,fmt=*,iostat=ierr) niter
          read(iunit,fmt=*,iostat=ierr) ipulay
          read(iunit,fmt=*,iostat=ierr) npulay
          if (ierr.ne.0) then
            write (io_lun, fmt='(A)') 'error reading pulay.negf -> reset '
            ipulay=1
            npulay=0
            niter=0
          end if
        end if
      end if
      call gcopy(ipulay)
      call gcopy(npulay)
      call gcopy(niter)



      do i=1,npulay
        write(pltfile,fmt='(i8)') i
        pltfile=adjustl(pltfile)
        pltfile="rho_pul_"//trim(pltfile)//".dat"
        if (inode.eq.ionode) write(io_lun,*) "pulay ",trim(pltfile),i,npulay
        rho3(1:nglx,1:ngly,1:nglz)=>denstmp(1:nglx*ngly*nglz)
        call read_rho(rho3,pltfile)
        call rearrange(denstmp,maxngrid,gs_packxb,gs_psendxb,gs_prxb,gs_unpackxb)
        rho_pul(:,i,1)=denstmp
        nullify(rho3)

        write(pltfile,fmt='(i8)') i
        pltfile=adjustl(pltfile)
        pltfile="resid_pul_"//trim(pltfile)//".dat"
        if (inode.eq.ionode) write(io_lun,*) "pulay ",trim(pltfile),i,npulay
        rho3(1:nglx,1:ngly,1:nglz)=>denstmp(1:nglx*ngly*nglz)
        call read_rho(rho3,pltfile)
        call rearrange(denstmp,maxngrid,gs_packxb,gs_psendxb,gs_prxb,gs_unpackxb)
        resid_pul(:,i,1)=denstmp
        nullify(rho3)
        k_resid_pul(:,i,1)=resid_pul(:,i,1)
        if (flag_Kerker) then
          if (flag_wdmetric) then
              if (inode.eq.ionode) write(io_lun,*) "wavekerker"
              call kerker_and_wdmetric(k_resid_pul(:,i,1), &
                                      denstmp, &
                                      maxngrid, q0, q1)
          else
              if (inode.eq.ionode) write(io_lun,*) "kerker"
              call kerker(k_resid_pul(:,i,1), maxngrid, q0)
          end if
          ! replace the k_resid_pul history at iPulay to new value
        else
          if (flag_wdmetric) then
              if (inode.eq.ionode) write(io_lun,*) "wave"
              denstmp=k_resid_pul(:,ipulay,1)
              call wdmetric(denstmp, k_resid_pul(:,i,1), &
                            maxngrid, q1)
          end if
        end if ! flag_Kerker

      end do


      pltfile="density_0.dat"
      rho3(1:nglx,1:ngly,1:nglz)=>denstmp(1:nglx*ngly*nglz)
      call read_rho(rho3,pltfile)
      call my_barrier()
      call rearrange(denstmp,maxngrid,gs_packxb,gs_psendxb,gs_prxb,gs_unpackxb)
      density(:,1)=denstmp
      nullify(rho3)


      npulay=npulay+1

      kmix=3
      if (npulay.gt.maxpulaySC) npulay=maxpulaySC


      call grab_matrix2("K_negf_",inode,nfile,Info,InfoGlob,lnegf=.true.)
      call my_barrier()
      call Matrix_CommRebuild(InfoGlob,info,Hrange,H_trans,matK,nfile,n_matrix=nspin)
      electrons=0d0
      call get_electronic_density(rho1, electrons, atomfns,     &
          H_on_atomfns(1), inode, ionode,  &
          maxngrid)
      electrons_tot = spin_factor * sum(electrons)
      if (inode == ionode) &
          write (io_lun, fmt='(A,3f12.6)') 'In Hnegf, electrons: ', electrons_tot,electrons

      pltfile="rho_negf.dat"
      call grid_out_plt(pltfile,"current density",250,rho1(:,1))

!~        pltfile="rhodiff_negf.dat"
!~        call grid_out_plt(pltfile,"current density",250,rho1(:,1)*2d0-density_atom)




      call fix_densities_negf("density_l.dat","density_r.dat",rho1(:,1),1d0)


      electrons1=grid_point_volume * &
                      rsum(n_my_grid_points, rho1(:,1), 1)
      call gsum(electrons1)
      electrons_tot = spin_factor * electrons1
      if (inode == ionode) &
        write (io_lun, fmt='(A,3f12.6)') 'In Hnegf, electrons after fix: ', electrons_tot,ne_in_cell


      call fix_densities_negf("density_l.dat","density_r.dat",density(:,1),1d0)





      if (flag_kerker) then
        cov_resid_pul=>k_resid_pul
      else
        cov_resid_pul=>resid_pul
      end if

      if (1.eq.1) then !((mod(niter+1,kmix).ne.0).and.(niter.ne.0)) then
        if (inode.eq.ionode) write(io_lun,*) "do pulay mixing ",mod(niter+1,2)
        ! in out
        call  update_pulay_history_negf(iPulay,density, rho1, &
                                  rho_pul, resid_pul, k_resid_pul,     &
                                  cov_resid_pul)
        density=0d0
        call get_pulay_optimal_rho(iPulay, density, npulay, A, rho_pul, &
                                resid_pul, k_resid_pul,     &
                                cov_resid_pul)
      else
        if (inode.eq.ionode) write(io_lun,*) "do linear mixing ",mod(niter+1,2)
        rho_pul(1:n_my_grid_points,iPulay,1) = density(1:n_my_grid_points,1)
        resid_pul(1:n_my_grid_points,ipulay,1)= rho1(1:n_my_grid_points,1)-density(1:n_my_grid_points,1)
        k_resid_pul(1:n_my_grid_points,ipulay,1)=resid_pul(1:n_my_grid_points,ipulay,1)

        if (flag_Kerker) then
          if (flag_wdmetric) then
              if (inode.eq.ionode) write(io_lun,*) "wavekerker"
              call kerker_and_wdmetric(k_resid_pul(:,ipulay,1), &
                                      denstmp, &
                                      maxngrid, q0, q1)
          else
              if (inode.eq.ionode) write(io_lun,*) "kerker"
              call kerker(k_resid_pul(:,ipulay,1), maxngrid, q0)
          end if
          ! replace the k_resid_pul history at iPulay to new value
        else
          if (flag_wdmetric) then
              if (inode.eq.ionode) write(io_lun,*) "wave"
              denstmp=k_resid_pul(:,ipulay,1)
              call wdmetric(denstmp, k_resid_pul(:,ipulay,1), &
                            maxngrid, q1)
          end if
        end if ! flag_Kerker
          density(:,1)=density(:,1)+A(1)*cov_resid_pul(:,ipulay,1)
      end if


      R0 = dot(n_my_grid_points, cov_resid_pul(1:n_my_grid_points,ipulay,1), 1,cov_resid_pul(1:n_my_grid_points,ipulay,1), 1)
      call gsum(R0)
      R0 = sqrt(grid_point_volume * R0) / ne_in_cell


      if (inode.eq.ionode) write(io_lun,fmt='(A,4e24.12)') "maxres ",0d0,R0,A(1),spin_factor

!      pltfile="density.info"
!      call fix_densities_negf(pltfile,density(:,1),1d0)

      write(pltfile,fmt='(i8)') ipulay
      pltfile=adjustl(pltfile)
      pltfile="rho_pul_"//trim(pltfile)//".dat"
      call grid_out_plt(pltfile,"rho_pul: ",250,rho_pul(:,ipulay,1))

      write(pltfile,fmt='(i8)') ipulay
      pltfile=adjustl(pltfile)
      pltfile="resid_pul_"//trim(pltfile)//".dat"
      call grid_out_plt(pltfile,"resid_pul: ",250,resid_pul(:,ipulay,1))


      if (inode.eq.ionode) then
        niter=niter+1
        ipulay=ipulay+1
        if (ipulay.gt.maxpulaySC) ipulay=1
        rewind(iunit)
        write(iunit,*) niter
        write(iunit,*) ipulay
        write(iunit,*) npulay
        close(iunit)
        inquire(file="convergence.negf",exist=lexists)
        if (.not.lexists) then
          open(newunit=iunit,file="convergence.negf",action="write",status="new",position="append")
          close(iunit)
        end if
        open(newunit=iunit,file="convergence.negf",action="write",status="old",position="append")
        write(iunit,fmt='(i8,e24.12)') niter,R0
        close(iunit)
      end if








!~         else if ((restart_Knegf).and.(read_gless)) then

!~           write(io_lun,*) "read current_density matrix..."
!~           call grab_matrix2("Gless_",inode,nfile,Info)

!~           call my_barrier()
!~           call Matrix_CommRebuild(Info,Hrange,H_trans,matK(1),nfile)
!~           write(io_lun,*) "Rebuild K mat",H_trans

!~           do xyz=1,3
!~             density=0d0
!~             call get_phi_m_grad_phi(density, electrons, atomfns,     &
!~                                      H_on_atomfns(1), inode, ionode, &
!~                                      maxngrid,xyz)
!~             write(pltfile,*) xyz
!~             pltfile=adjustl(pltfile)
!~             pltfile="current_density_"//trim(pltfile)//".plt"
!~             write(io_lun,*) "output charge density ",trim(pltfile)
!~             call grid_out_plt(pltfile,"current density",250,density(:,1))
!~           end do
!~           call my_barrier ()
!~           call end_comms ()
!~           stop
!~         end if




        call print_grid_data(4)
        call print_grid_data(3)




    end subroutine

    subroutine set_atomic_density_set_dxyz(flag_set_density,lr,level)

      use datatypes
      use numbers
      use global_module,       only: rcellx, rcelly, rcellz, id_glob, &
                                    ni_in_cell, iprint_SC,           &
                                    species_glob, dens, ne_in_cell,  &
                                    ne_spin_in_cell,                 &
                                    IPRINT_TIME_THRES3, nspin,       &
                                    spin_factor,                     &
                                    flag_fix_spin_population,        &
                                    flag_neutral_atom, io_lun,       &
                                    dx_negf,dy_negf,dz_negf,io_lun
      use block_module,        only: nx_in_block, ny_in_block,        &
                                    nz_in_block, n_pts_in_block
      use group_module,        only: blocks, parts
      use primary_module,      only: domain
      use cover_module,        only: DCS_parts
      use set_blipgrid_module, only: naba_atoms_of_blocks
      use GenComms,            only: my_barrier, cq_abort, inode, ionode, gsum
      use atomic_density,      only: atomic_density_table      
      use dimens,              only: n_my_grid_points, grid_point_volume
      use GenBlas,             only: rsum, scal
      use timer_module,        only: WITH_LEVEL
      use maxima_module,  only: maxngrid
      use species_module, only: charge, charge_up, charge_dn
      use density_module, only: density_atom,density_scale,flag_InitialAtomicSpin

      implicit none

      ! Passed Variables
      integer, optional :: level
      logical :: flag_set_density
      integer :: lr

      ! Local Variables
      integer :: ipart, jpart, ind_part, ia, ii, icover, ig_atom
      integer :: the_species
      integer :: ix, iy, iz, j, iblock, ipoint, igrid
      integer :: no_of_ib_ia, offset_position
      integer :: position, iatom, icheck, spin

      real(double) :: dcellx_block, dcelly_block, dcellz_block
      real(double) :: dcellx_grid, dcelly_grid, dcellz_grid
      real(double) :: dcellx_grid2, dcelly_grid2, dcellz_grid2
      real(double) :: xatom, yatom, zatom, step, loc_cutoff
      real(double) :: xblock, yblock, zblock, alpha
      real(double) :: xblock2, yblock2, zblock2
      real(double) :: dx, dy, dz, rx, ry, rz, rrr2, r_from_i,rr2,rx2, ry2, rz2
      real(double) :: a, b, c, d, r1, r2, r3, r4, rr
      ! local charge density returned from splint routine
      real(double)   :: local_density,local_density2, scale, scale_spin_up, scale_spin_dn,dz_off

      ! logical flag to warn if splint routine called out of the
      ! tabulated range.
      logical :: range_flag
      integer :: stat,nnn,mmm
      character(256) :: outfile
      integer :: iunit
!*  ***lat<$




      if (inode == ionode .and. iprint_SC >= 2) then
        write (io_lun, fmt='(2x,"Entering set_density")')
        if (flag_InitialAtomicSpin) write (io_lun, fmt='(2x,"Initial atomic spins are read from input file")')
      endif


      !call sub_enter_output('set_density',2,'5')

      ! initialize density

      density_atom = zero

      ! determine the block and grid spacing
      dcellx_block = rcellx / blocks%ngcellx
      dcelly_block = rcelly / blocks%ngcelly
      dcellz_block = rcellz / blocks%ngcellz
      dcellx_grid = dcellx_block / nx_in_block
      dcelly_grid = dcelly_block / ny_in_block
      dcellz_grid = dcellz_block / nz_in_block

      dcellx_grid2 = dx_negf
      dcelly_grid2 = dy_negf
      dcellz_grid2 = dz_negf
!~       dcellx_block2=dx_negf*nx_in_block
!~       dcelly_block2=dy_negf*ny_in_block
!~       dcellz_block2=dz_negf*nz_in_block

      if (lr.eq.1) then
        dz_off=0d0
      else if (lr.eq.2) then
        dz_off=rcellz-nz_in_block*blocks%ngcellz*dz_negf
      end if


      ! loop around grid points in this domain, and for each
      ! point, get contributions to the charge density from atoms which are
      ! within the cutoff distance to that grid point


      nnn=0
      mmm=0


      do iblock = 1, domain%groups_on_node ! loop over blocks of grid points
        !write(io_lun,*) 'Block ',iblock
        ! determine the position of this block
        xblock = (domain%idisp_primx(iblock) + domain%nx_origin - 1) * dcellx_block
        yblock = (domain%idisp_primy(iblock) + domain%ny_origin - 1) * dcelly_block
        zblock = (domain%idisp_primz(iblock) + domain%nz_origin - 1) * dcellz_block
        ! if there are neighbour partitions
        if (naba_atoms_of_blocks(dens)%no_of_part(iblock) > 0) then
            iatom = 0
            ! loop over neighbour partitions of this block
            do ipart = 1, naba_atoms_of_blocks(dens)%no_of_part(iblock)
              jpart = naba_atoms_of_blocks(dens)%list_part(ipart,iblock)
              if (jpart > DCS_parts%mx_gcover) then
                  call cq_abort('set_ps: JPART ERROR ', ipart, jpart)
              endif
              ind_part = DCS_parts%lab_cell(jpart)
              ! .... then atoms in this partition.
              do ia = 1, naba_atoms_of_blocks(dens)%no_atom_on_part(ipart, iblock)
                  iatom = iatom + 1
                  ii = naba_atoms_of_blocks(dens)%list_atom(iatom, iblock)
                  icover = DCS_parts%icover_ibeg(jpart) + ii - 1
                  ig_atom = id_glob(parts%icell_beg(ind_part) + ii - 1)

                  if (parts%icell_beg(ind_part) + ii - 1 > ni_in_cell) then
                    call cq_abort('set_ps: globID ERROR ', &
                                  ii, parts%icell_beg(ind_part))
                  endif
                  if (icover > DCS_parts%mx_mcover) then
                    call cq_abort('set_ps: icover ERROR ', &
                                  icover, DCS_parts%mx_mcover)
                  endif

                  ! determine the position of the current atom
                  xatom = DCS_parts%xcover(icover)
                  yatom = DCS_parts%ycover(icover)
                  zatom = DCS_parts%zcover(icover)

                  ! determine the type of the current atom
                  the_species = species_glob(ig_atom)
                  loc_cutoff = atomic_density_table(the_species)%cutoff
                  ! step in the density table
                  step = loc_cutoff / &
                        real(atomic_density_table(the_species)%length - 1, &
                              double)

                  ! determine the spin-dependent scaling factor
                  if (flag_InitialAtomicSpin) then
                    if (charge_up(the_species).eq.zero .and. charge_dn(the_species).eq.zero) then
                        scale_spin_up = half
                        scale_spin_dn = half
                    else
                        scale_spin_up = charge_up(the_species)/charge(the_species)
                        scale_spin_dn = charge_dn(the_species)/charge(the_species)
                    endif
                  endif

                  icheck = 0
                  ipoint = 0
                  ! loop over the grid points in the current block
                  do iz = 1, nz_in_block
                    do iy = 1, ny_in_block
                        do ix = 1, nx_in_block
                          ipoint = ipoint + 1
                          igrid = n_pts_in_block * (iblock - 1) + ipoint
                          ! position= offset_position + ipoint
                          if (igrid > n_my_grid_points) &
                                call cq_abort('set_density: igrid error ', &
                                              igrid, n_my_grid_points)
                          dx = dcellx_grid * (ix - 1)
                          dy = dcelly_grid * (iy - 1)
                          dz = dcellz_grid * (iz - 1)
                          ! determine separation between the current
                          ! grid point and atom
                          rx = xblock + dx - xatom
                          ry = yblock + dy - yatom
                          rz = zblock + dz - zatom
                          rrr2 = rx * rx + ry * ry + rz * rz
                          rx = (xblock + dx)/dcellx_grid*dcellx_grid2 - xatom
                          ry = (yblock + dy)/dcelly_grid*dcelly_grid2 - yatom
                          rz = (zblock + dz)/dcellz_grid*dcellz_grid2+dz_off - zatom
                          rr2 = rx * rx + ry * ry + rz * rz
                          if (rrr2< loc_cutoff*loc_cutoff) nnn=nnn+1

                          if (rr2 < loc_cutoff * loc_cutoff) then
                            r_from_i = sqrt(rr2)
                            j = floor(r_from_i/step)+1
                            if(j+1<=atomic_density_table(the_species)%length) then
                               rr = real(j,double) * step
                               a = ( rr - r_from_i ) / step
                               b = one - a
                               c = a * ( a * a - one ) * step * step / six
                               d = b * ( b * b - one ) * step * step / six
                               r1 = atomic_density_table(the_species)%table(j)
                               r2 = atomic_density_table(the_species)%table(j+1)
                               r3 = atomic_density_table(the_species)%d2_table(j)
                               r4 = atomic_density_table(the_species)%d2_table(j+1)
                               ! Store the density for this grid point
                               local_density = a * r1 + b * r2 + c * r3 + d * r4
                               density_atom(igrid) = density_atom(igrid) + local_density ! both up and down
                            end if
                          end if ! if this point is within cutoff
                        end do !ix gridpoints
                    end do  !iy gridpoints
                  end do   !iz gridpoints
              end do ! end loop over atoms in this partition
            end do ! end loop over partitions
        end if ! end if there are neighbour atoms
      end do ! end loop over blocks

      ! renormalise density
      local_density = zero
      do iz = 1, n_my_grid_points
        local_density = local_density + density_atom(iz)
      end do

      local_density = local_density * (dcellx_grid2*dcelly_grid2*dcellz_grid2)!grid_point_volume
      call gsum(local_density)
      call gsum(nnn)
      call gsum(mmm)


      scale = ne_in_cell/local_density
!~       write(io_lun,fmt='(A,2i8,e32.16)') "nnn,mmm ",nnn,mmm,scale
!~       write(io_lun,*) "rcut ", loc_cutoff, loc_cutoff * loc_cutoff
!~       store_density_atom(:) = store_density_atom(:) * scale
!~       write(0,fmt='(3e24.12)') dcellx_grid2,dcelly_grid2,dcellz_grid2
!~       write(0,fmt='(3e24.12)') dcellx_grid,dcelly_grid,dcellz_grid
!      write(io_lun,fmt='(A,4e24.12)') "grid ratio",dcellx_grid2*dcelly_grid2*dcellz_grid2,dcellx_grid*dcelly_grid*dcellz_grid,grid_point_volume,ne_in_cell



      call my_barrier()


    end subroutine set_atomic_density_set_dxyz

    subroutine set_density_pcc_set_dxyz(lr)

      use datatypes
      use numbers,             only: zero, one, six
      use global_module,       only: rcellx, rcelly, rcellz, id_glob, &
                                    ni_in_cell, iprint_SC,           &
                                    species_glob, dens, ne_in_cell,  &
                                    IPRINT_TIME_THRES3,&
                                    dx_negf,dy_negf,dz_negf,io_lun
      use block_module,        only: nx_in_block, ny_in_block,        &
                                    nz_in_block, n_pts_in_block
      use group_module,        only: blocks, parts
      use primary_module,      only: domain
      use cover_module,        only: DCS_parts
      use set_blipgrid_module, only: naba_atoms_of_blocks
      use GenComms,            only: my_barrier, cq_abort, inode, ionode, gsum
      use pseudo_tm_info,      only: pseudo
      use dimens,              only: n_my_grid_points, grid_point_volume
      use GenBlas,             only: rsum, scal
      use timer_module
      use density_module, only: density_pcc

      implicit none

      integer :: lr

      ! Local Variables
      integer :: ipart, jpart, ind_part, ia, ii, icover, ig_atom
      integer :: the_species
      integer :: ix, iy, iz, j, iblock, ipoint, igrid
      integer :: no_of_ib_ia, offset_position
      integer :: position, iatom, icheck

      real(double) :: dcellx_block, dcelly_block, dcellz_block
      real(double) :: dcellx_grid, dcelly_grid, dcellz_grid
      real(double) :: dcellx_grid2, dcelly_grid2, dcellz_grid2
      real(double) :: xatom, yatom, zatom
      real(double) :: pcc_cutoff, pcc_step
      real(double) :: xblock, yblock, zblock,  alpha
      real(double) :: dx, dy, dz, rx, ry, rz, rr2, r_from_i, dz_off
      real(double) :: a, b, c, d, r1, r2, r3, r4, rr
      real(double) :: pcc_density !P.C.C. charge density returned from splint routine

      ! logical flag to warn if splint routine called out of the tabulated range.
      logical :: range_flag

      if (inode == ionode .and. iprint_SC >= 2) &
          write (io_lun, fmt='(2x,"Entering set_density_pcc")')


      density_pcc = zero  ! initialize density
      ! write(io_lun,*) 'Size of density: ',size(density)
      ! call scal(n_my_grid_points,zero,density,1)

      ! determine the block and grid spacing
      dcellx_block=rcellx/blocks%ngcellx; dcellx_grid=dcellx_block/nx_in_block
      dcelly_block=rcelly/blocks%ngcelly; dcelly_grid=dcelly_block/ny_in_block
      dcellz_block=rcellz/blocks%ngcellz; dcellz_grid=dcellz_block/nz_in_block

      dcellx_grid2 = dx_negf
      dcelly_grid2 = dy_negf
      dcellz_grid2 = dz_negf
!~       dcellx_block2=dx_negf*nx_in_block
!~       dcelly_block2=dy_negf*ny_in_block
!~       dcellz_block2=dz_negf*nz_in_block

      if (lr.eq.1) then
        dz_off=0d0
      else if (lr.eq.2) then
        dz_off=rcellz-nz_in_block*blocks%ngcellz*dz_negf
      end if


      ! loop around grid points in this domain, and for each
      ! point, get contributions to the charge density from atoms which are
      ! within the cutoff distance to that grid point

      call my_barrier()

      do iblock = 1, domain%groups_on_node ! loop over blocks of grid points
        !write(io_lun,*) 'Block ',iblock
        ! determine the position of this block
        xblock=(domain%idisp_primx(iblock)+domain%nx_origin-1)*dcellx_block
        yblock=(domain%idisp_primy(iblock)+domain%ny_origin-1)*dcelly_block
        zblock=(domain%idisp_primz(iblock)+domain%nz_origin-1)*dcellz_block
        if(naba_atoms_of_blocks(dens)%no_of_part(iblock) > 0) then ! if there are neighbour partitions
            iatom=0
            ! loop over neighbour partitions of this block
            do ipart=1,naba_atoms_of_blocks(dens)%no_of_part(iblock)   !for now, dens is used even for P.C.C.
              jpart=naba_atoms_of_blocks(dens)%list_part(ipart,iblock)
              if(jpart > DCS_parts%mx_gcover) then
                  call cq_abort('set_ps: JPART ERROR ',ipart,jpart)
              endif
              ind_part=DCS_parts%lab_cell(jpart)
              ! .... then atoms in this partition.
              do ia=1,naba_atoms_of_blocks(dens)%no_atom_on_part(ipart,iblock)
                  iatom=iatom+1
                  ii = naba_atoms_of_blocks(dens)%list_atom(iatom,iblock)
                  icover= DCS_parts%icover_ibeg(jpart)+ii-1
                  ig_atom= id_glob(parts%icell_beg(ind_part)+ii-1)
                  ! determine the type of the current atom
                  the_species=species_glob(ig_atom)
                  ! for P.C.C.
                  if (.NOT. pseudo(the_species)%flag_pcc) then
                    cycle
                  endif

                  if(parts%icell_beg(ind_part) + ii-1 > ni_in_cell) then
                    call cq_abort('set_ps: globID ERROR ', &
                          ii,parts%icell_beg(ind_part))
                  endif
                  if(icover > DCS_parts%mx_mcover) then
                    call cq_abort('set_ps: icover ERROR ', &
                          icover,DCS_parts%mx_mcover)
                  endif

                  ! determine the position of the current atom
                  xatom=DCS_parts%xcover(icover)
                  yatom=DCS_parts%ycover(icover)
                  zatom=DCS_parts%zcover(icover)
!                  write(6,fmt='(i8,3e24.12)') ia,xatom,yatom,zatom
                  pcc_cutoff = pseudo(the_species)%chpcc%cutoff
                  ! step in the density table
                  pcc_step = pseudo(the_species)%chpcc%delta
                  icheck=0
                  ipoint=0
                  ! loop over the grid points in the current block
                  do iz=1,nz_in_block
                    do iy=1,ny_in_block
                        do ix=1,nx_in_block
                          ipoint=ipoint+1
                          igrid=n_pts_in_block*(iblock-1)+ipoint
                          ! position= offset_position + ipoint
                          if(igrid > n_my_grid_points) &
                                call cq_abort('set_density: igrid error ', &
                                igrid, n_my_grid_points)
                          dx=dcellx_grid*(ix-1)
                          dy=dcelly_grid*(iy-1)
                          dz=dcellz_grid*(iz-1)
                          ! determine separation between the current
                          ! grid point and atom
!~                            rx=xblock+dx-xatom
!~                            ry=yblock+dy-yatom
!~                            rz=zblock+dz-zatom
                            rx = (xblock + dx)/dcellx_grid*dcellx_grid2 - xatom
                            ry = (yblock + dy)/dcelly_grid*dcelly_grid2 - yatom
                            rz = (zblock + dz)/dcellz_grid*dcellz_grid2+dz_off - zatom
                          rr2 = rx * rx + ry * ry + rz * rz
                          if(rr2 < pcc_cutoff * pcc_cutoff) then
                            r_from_i = sqrt(rr2)
                            j = floor(r_from_i/pcc_step) + 1
                            if(j+1<=pseudo(the_species)%chpcc%n) then
                               rr = real(j,double) * pcc_step
                               a = ( rr - r_from_i ) / pcc_step
                               b = one - a
                               c = a * ( a * a - one ) * pcc_step * pcc_step / six
                               d = b * ( b * b - one ) * pcc_step * pcc_step / six
                               r1 = pseudo(the_species)%chpcc%f(j)
                               r2 = pseudo(the_species)%chpcc%f(j+1)
                               r3 = pseudo(the_species)%chpcc%d2(j)
                               r4 = pseudo(the_species)%chpcc%d2(j+1)
                               pcc_density = a * r1 + b * r2 + c * r3 + d * r4
                               ! recalculate the density for this grid point
                               density_pcc(igrid) = density_pcc(igrid) + pcc_density
                            end if
                          endif ! if this point is within cutoff
                          ! test output of densities
                          ! print '(4(f10.6,a))', xblock+dx, " *", &
                          !      yblock+dy, " *", zblock+dz, " *", &
                          !      density(igrid), " *"
                        enddo !ix gridpoints
                    enddo  !iy gridpoints
                  enddo   !iz gridpoints
              enddo ! end loop over atoms in this partition
            enddo ! end loop over partitions
        endif ! end if there are neighbour atoms
      enddo ! end loop over blocks
      ! local_density = zero
      ! do iz=1,n_my_grid_points
      !   local_density = local_density + density(iz)
      ! end do
      ! local_density = local_density*grid_point_volume
      ! local_density = grid_point_volume * rsum( n_my_grid_points, density, 1 )

      ! call gsum(local_density)
      ! Correct electron density
      ! density_scale = ne_in_cell/local_density
      ! density = density_scale*density
      ! if(inode.eq.ionode.AND.iprint_SC>0) &
      !      write(io_lun,fmt='(10x,"In set_density, electrons: ",f20.12)') &
      !      density_scale*local_density
      call my_barrier()

      return
    end subroutine set_density_pcc_set_dxyz

    subroutine single_PAO_to_grid_set_dxyz(pao_fns,lr)

      use datatypes
      use primary_module, ONLY: bundle
      use GenComms, ONLY: my_barrier, cq_abort, mtime,gsum
      use dimens, ONLY: r_h
      use GenComms, ONLY: inode, ionode,io_lun
      use numbers
      use global_module, ONLY: rcellx,rcelly,rcellz,id_glob,ni_in_cell, iprint_basis, species_glob, atomf,dx_negf,dy_negf,dz_negf
      use species_module, ONLY: species, npao_species
      !  At present, these arrays are dummy arguments.
      use block_module, ONLY : nx_in_block,ny_in_block,nz_in_block, n_pts_in_block
      use group_module, ONLY : blocks, parts
      use primary_module, ONLY: domain
      use cover_module, ONLY: DCS_parts
      use set_blipgrid_module, ONLY : naba_atoms_of_blocks
      use angular_coeff_routines, ONLY : evaluate_pao
      use functions_on_grid, ONLY: gridfunctions, fn_on_grid
      use pao_format
      use PAO_grid_transform_module, only : check_block

      implicit none
      integer,intent(in) :: pao_fns,lr

      !local
      real(double):: dcellx_block,dcelly_block,dcellz_block
      real(double):: dcellx_block2,dcelly_block2,dcellz_block2
      real(double) :: dcellx_grid, dcelly_grid, dcellz_grid
      real(double) :: dcellx_grid2, dcelly_grid2, dcellz_grid2
      integer :: ipart,jpart,ind_part,ia,ii,icover,ig_atom
      real(double):: xatom,yatom,zatom,alpha,step
      real(double):: xblock,yblock,zblock
      integer :: the_species
      integer :: j,iblock,the_l,ipoint, igrid
      real(double) :: r_from_i
      real(double) :: rr,a,b,c,d,x,y,z,nl_potential
      integer :: no_of_ib_ia, offset_position
      integer :: position,iatom
      integer :: stat, nl, npoint, ip, this_nsf
      integer :: i,m, m1min, m1max,acz,m1,l1,count1
      integer     , allocatable :: ip_store(:)
      real(double), allocatable :: x_store(:)
      real(double), allocatable :: y_store(:)
      real(double), allocatable :: z_store(:)
      real(double), allocatable :: r_store(:)
      real(double) :: coulomb_energy
      real(double) :: rcut
      real(double) :: r1, r2, r3, r4, core_charge, gauss_charge,dz_off
      real(double) :: val, theta, phi, r_tmp
      integer :: nnn,mmm

      allocate(ip_store(n_pts_in_block),x_store(n_pts_in_block),y_store(n_pts_in_block),z_store(n_pts_in_block), &
          r_store(n_pts_in_block))

      ! --  Start of subroutine  ---

      dcellx_block=rcellx/blocks%ngcellx
      dcelly_block=rcelly/blocks%ngcelly
      dcellz_block=rcellz/blocks%ngcellz

      dcellx_grid=dcellx_block/nx_in_block
      dcelly_grid=dcelly_block/ny_in_block
      dcellz_grid=dcellz_block/nz_in_block

      dcellx_grid2 = dx_negf
      dcelly_grid2 = dy_negf
      dcellz_grid2 = dz_negf

      dcellx_block2=dx_negf*nx_in_block!rcellx/blocks%ngcellx
      dcelly_block2=dy_negf*ny_in_block!rcelly/blocks%ngcelly
      dcellz_block2=dz_negf*nz_in_block!rcellz/blocks%ngcellz


      if (lr.eq.1) then
        dz_off=0d0
      else if (lr.eq.2) then
        dz_off=rcellz-nz_in_block*blocks%ngcellz*dz_negf
      end if


      call my_barrier()

      no_of_ib_ia = 0
      gridfunctions(pao_fns)%griddata = zero

      ! loop arround grid points in the domain, and for each
      ! point, get the PAO values
      nnn=0
      do iblock = 1, domain%groups_on_node ! primary set of blocks
        xblock=(domain%idisp_primx(iblock)+domain%nx_origin-1)*dcellx_block2
        yblock=(domain%idisp_primy(iblock)+domain%ny_origin-1)*dcelly_block2
        zblock=(domain%idisp_primz(iblock)+domain%nz_origin-1)*dcellz_block2+dz_off
        if(naba_atoms_of_blocks(atomf)%no_of_part(iblock) > 0) then ! if there are naba atoms
            iatom=0
            do ipart=1,naba_atoms_of_blocks(atomf)%no_of_part(iblock)
              jpart=naba_atoms_of_blocks(atomf)%list_part(ipart,iblock)
              if(jpart > DCS_parts%mx_gcover) then
                  call cq_abort('single_PAO_to_grid: JPART ERROR ',ipart,jpart)
              endif
              ind_part=DCS_parts%lab_cell(jpart)
              do ia=1,naba_atoms_of_blocks(atomf)%no_atom_on_part(ipart,iblock)
                  iatom=iatom+1
                  ii = naba_atoms_of_blocks(atomf)%list_atom(iatom,iblock)
                  icover= DCS_parts%icover_ibeg(jpart)+ii-1
                  ig_atom= id_glob(parts%icell_beg(ind_part)+ii-1)

                  if(parts%icell_beg(ind_part) + ii-1 > ni_in_cell) then
                    call cq_abort('single_PAO_to_grid: globID ERROR ', &
                          ii,parts%icell_beg(ind_part))
                  endif
                  if(icover > DCS_parts%mx_mcover) then
                    call cq_abort('single_PAO_to_grid: icover ERROR ', &
                          icover,DCS_parts%mx_mcover)
                  endif

                  xatom=DCS_parts%xcover(icover)
                  yatom=DCS_parts%ycover(icover)
                  zatom=DCS_parts%zcover(icover)
                  the_species=species_glob(ig_atom)

                  !calculates distances between the atom and integration grid points
                  !in the block and stores which integration grids are neighbours.
                  rcut = r_h + RD_ERR
!~                   call check_block (xblock,yblock,zblock,xatom,yatom,zatom, rcut, &  ! in
!~                        npoint,ip_store,r_store,x_store,y_store,z_store,n_pts_in_block) !out
                  call check_block_set_dxyz (xblock,yblock,zblock,xatom,yatom,zatom, rcut, &  ! in
                      npoint,ip_store,r_store,x_store,y_store,z_store,n_pts_in_block,dcellx_grid2,dcelly_grid2,dcellz_grid2) !out

                  r_from_i = sqrt((xatom-xblock)**2+(yatom-yblock)**2+ &
                      (zatom-zblock)**2 )

                  if(npoint > 0) then

                    !offset_position = (no_of_ib_ia-1) * npao * n_pts_in_block
                    offset_position = no_of_ib_ia
                    do ip=1,npoint
                      nnn=nnn+1
                        ipoint=ip_store(ip)
                        position= offset_position + ipoint
                        if(position > gridfunctions(pao_fns)%size) call cq_abort &
                            ('single_pao_to_grid: position error ', position, gridfunctions(pao_fns)%size)

                        r_from_i = r_store(ip)
!~                         x = (x_store(ip)+xatom)/dcellx_grid*dcellx_grid2-xatom
!~                         y = (y_store(ip)+yatom)/dcelly_grid*dcelly_grid2-yatom
!~                         z = (z_store(ip)+zatom)/dcellz_grid*dcellz_grid2+dz_off-zatom
                        x = x_store(ip)
                        y = y_store(ip)
                        z = z_store(ip)

                        ! For this point-atom offset, we accumulate the PAO on the grid
                        count1 = 1
                        do l1 = 0,pao(the_species)%greatest_angmom
                          do acz = 1,pao(the_species)%angmom(l1)%n_zeta_in_angmom
                              do m1=-l1,l1
                                call evaluate_pao(the_species,l1,acz,m1,x,y,z,val)
!~                                  write(6,fmt='(4i8,4e24.12)') ip,l1,acz,m1,x,y,z,val
                                if(position+(count1-1)*n_pts_in_block > gridfunctions(pao_fns)%size) &
                                      call cq_abort('single_pao_to_grid: position error ', &
                                      position, gridfunctions(pao_fns)%size)
                                gridfunctions(pao_fns)%griddata(position+(count1-1)*n_pts_in_block) = val
                                count1 = count1+1
                              end do ! m1
                          end do ! acz
                        end do ! l1
                    enddo ! ip=1,npoint
                  endif! (npoint > 0) then
                  no_of_ib_ia = no_of_ib_ia + npao_species(the_species)*n_pts_in_block
              enddo ! naba_atoms
            enddo ! naba_part
        endif !(naba_atoms_of_blocks(atomf)%no_of_part(iblock) > 0) !naba atoms?
      enddo ! iblock : primary set of blocks
      call my_barrier()

      deallocate(ip_store,x_store,y_store,z_store,r_store)
      call gsum(nnn)
      if (inode.eq.ionode) write(io_lun,*) "nnn for density ",nnn
      return
    end subroutine Single_PAO_to_grid_set_dxyz

    subroutine check_block_set_dxyz &
        (xblock,yblock,zblock,xatom,yatom,zatom,rcut, &
        npoint, ip_store, r_store, x_store, y_store, z_store,blocksize,&
        dcellx_grid,dcelly_grid,dcellz_grid)

      use numbers
      use global_module, ONLY: rcellx,rcelly,rcellz,io_lun
      use group_module,  ONLY: blocks
      use block_module,  ONLY: nx_in_block,ny_in_block,nz_in_block!, &
!           n_pts_in_block


      implicit none
      !Passed
      integer :: blocksize
      real(double), intent(in):: xblock, yblock, zblock
      real(double), intent(in):: xatom, yatom, zatom, rcut
      integer, intent(out) :: npoint, ip_store(blocksize)
      real(double),intent(out) :: r_store(blocksize)
      real(double),intent(out) :: x_store(blocksize)
      real(double),intent(out) :: y_store(blocksize)
      real(double),intent(out) :: z_store(blocksize)
      !Local
      real(double):: dcellx_block,dcelly_block,dcellz_block
      real(double):: dcellx_grid, dcelly_grid, dcellz_grid
      real(double):: dx, dy, dz
      integer :: ipoint, iz, iy, ix
      real(double) ::  r2, r_from_i, rx, ry, rz, x, y, z, rcut2


      rcut2 = rcut* rcut
!~       dcellx_block=rcellx/blocks%ngcellx; dcellx_grid=dcellx_block/nx_in_block
!~       dcelly_block=rcelly/blocks%ngcelly; dcelly_grid=dcelly_block/ny_in_block
!~       dcellz_block=rcellz/blocks%ngcellz; dcellz_grid=dcellz_block/nz_in_block
!      write(6,*) rcut2,rcut
      ipoint=0
      npoint=0
      do iz=1,nz_in_block
        do iy=1,ny_in_block
            do ix=1,nx_in_block
              ipoint=ipoint+1


              dx=dcellx_grid*(ix-1)
              dy=dcelly_grid*(iy-1)
              dz=dcellz_grid*(iz-1)

              rx=xblock+dx-xatom
              ry=yblock+dy-yatom
              rz=zblock+dz-zatom
              r2 = rx * rx + ry * ry + rz * rz

              if(r2 < rcut2) then
                  npoint=npoint+1
                  if(npoint>blocksize) write(io_lun,*) 'ALERT: Too many points ! ',npoint,blocksize
                  r_from_i = sqrt( r2 )

!                  ! direction cosines needed for spher harms
!                  if ( r_from_i > RD_ERR ) then
!                     x = rx / r_from_i
!                     y = ry / r_from_i
!                     z = rz / r_from_i
!                  else
!                     x = zero ; y = zero ; z = zero
!                     !!   r_from_i = RD_ERR !!   04/07/2002 TM
!                  end if
                  x = rx
                  y = ry
                  z = rz
                  ip_store(npoint)=ipoint
                  r_store(npoint)=r_from_i
                  x_store(npoint)=x
                  y_store(npoint)=y
                  z_store(npoint)=z
              endif  ! (r2 < rcut2) then

            enddo ! ix=1,nx_in_block
        enddo ! iy=1,ny_in_block
      enddo ! iz=1,nz_in_block
      return
    end subroutine check_block_set_dxyz



!   optain current density on the grid, i.e., j_{x,y,z}=sum_nm phi_n*(G^<-G^<')_nm*grad_{x,y,z}*phi_m
!   essentially slightly modified get_electronic_density, can maybe merged with this subroutine
!   turned off timers
      subroutine get_phi_m_grad_phi(denout, electrons, atom_fns, &
                                      atom_fns_K, inode, ionode, size,xyz,level)

        use datatypes
        use numbers
        use GenBlas,                     only: scal, rsum
        use mult_module,                 only: matK, matKatomf, SF_to_AtomF_transform, &
                                              free_temp_matrix
        use matrix_data,                 only: Hrange
        use dimens,                      only: n_my_grid_points, grid_point_volume
        use block_module,                only: n_pts_in_block
        use set_bucket_module,           only: rem_bucket, atomf_H_atomf_rem
        use calc_matrix_elements_module, only: act_on_vectors_new
        use primary_module,              only: domain
        use set_blipgrid_module,         only: naba_atoms_of_blocks
        use GenComms,                    only: gsum
        use global_module,               only: iprint_SC, ni_in_cell, &
                                              flag_Becke_weights, nspin, &
                                              spin_factor, sf, paof, atomf
        use functions_on_grid,           only: gridfunctions, fn_on_grid, &
                                              allocate_temp_fn_on_grid, &
                                              free_temp_fn_on_grid
        use PAO_grid_transform_module,   only: single_PAO_to_grad
        use global_module,          only: io_lun, iprint_SC

        implicit none

        ! Passed variables
        real(double), dimension(:)   :: electrons
        real(double), dimension(:,:) :: denout
        integer :: inode, ionode, size, xyz
        integer :: atom_fns, atom_fns_K
        integer, optional :: level

        ! Local variables
!~         type(cq_timer) :: backtrace_timer
!~         integer        :: backtrace_level
        integer        :: blk, i_count_alpha, n, n_i, n_point, spin, tmp_fn, direction

    !****lat<$
!~         if (       present(level) ) backtrace_level = level+1
!~         if ( .not. present(level) ) backtrace_level = -10
!~         call start_backtrace(t=backtrace_timer,who='get_electronic_density',&
!~              where=area,level=backtrace_level,echo=.true.)
    !****lat>$

!~         call start_timer(tmr_std_chargescf)

        if (inode == ionode .and. iprint_SC >= 2) &
            write (io_lun,fmt='(2x,"Entering get_electronic_density")')
        if (inode.eq.ionode) write(io_lun,*) "domain%groups_on_node",domain%groups_on_node
        spin=1
        direction=xyz
          ! matK->matKatomf backtransformation for contracted SFs
          if (atomf.ne.sf) call SF_to_AtomF_transform(matK(spin), matKatomf(spin), spin, Hrange)

          tmp_fn = allocate_temp_fn_on_grid(atomf)
          gridfunctions(tmp_fn)%griddata = zero
          call single_PAO_to_grad(direction, tmp_fn)



          gridfunctions(atom_fns_K)%griddata = zero

          call act_on_vectors_new(inode-1, rem_bucket(atomf_H_atomf_rem), &
                                  matKatomf(spin), atom_fns_K,atom_fns) !tmp_fn)
          denout(:,spin) = zero
          i_count_alpha = 0
          do blk = 1, domain%groups_on_node
              n_point = (blk - 1) * n_pts_in_block
              i_count_alpha = (naba_atoms_of_blocks(atomf)%ibegin_blk_orb(blk)-1) * n_pts_in_block
              !if(blk == 1) then
              !   i_count_alpha = 0
              !else
              !   i_count_alpha= i_count_alpha+ &
              !        naba_atoms_of_blocks(atomf)%no_of_atom(blk-1)*nsf*n_pts_in_block
              !endif
              if (naba_atoms_of_blocks(atomf)%no_of_atom(blk) > 0) then   !TM 30/Jun/2003
                !do n_i=1, naba_atoms_of_blocks(atomf)%no_of_atom(blk)*NSF*n_pts_in_block, &
                do n_i = 1, naba_atoms_of_blocks(atomf)%no_of_orb(blk)*n_pts_in_block, n_pts_in_block
                    do n=1, n_pts_in_block
                        denout(n_point+n, spin) =                                       &
                            denout(n_point+n, spin) +                                  &
                            gridfunctions(atom_fns_K)%griddata(i_count_alpha+n_i+n-1) * &
                            gridfunctions(tmp_fn)%griddata(i_count_alpha+n_i+n-1)
                    end do
                end do
              end if !(naba_atoms_of_blocks(atomf)%no_of_atom(blk) > 0)   !TM 30/Jun/2003
          end do ! blk

    !~     end do ! direction


        call free_temp_fn_on_grid(tmp_fn)

        ! atom_fns_K is using the same memory as h_on_atomfns, (in other
        ! words we are using h_on_atomfns as temorary storage), so lets be
        ! safe and set it back to zero
        gridfunctions(atom_fns_K)%griddata = zero

!~         call stop_timer(tmr_std_chargescf)

    !****lat<$
!~         call stop_backtrace(t=backtrace_timer,who='get_electronic_density',echo=.true.)
    !****lat>$

        return
      end subroutine get_phi_m_grad_phi

      subroutine get_sysinfo(iwhat,icount,maxijk)

        use GenComms, ONLY : my_barrier
        use global_module, ONLY: species_glob, atom_coord,ni_in_cell,rcellx,rcelly,rcellz,&
      &                         efermi_out
        use primary_module, ONLY: bundle
        use input_module, ONLY: io_assign, io_close
        use species_module, ONLY : nsf_species, species_label
        use support_spec_format, ONLY : flag_one_to_one,flag_paos_atoms_in_cell
        use GenComms,        only: inode, ionode
        use units, only : BohrToAng
        implicit none

!i  nput:
        integer :: iwhat
        integer, optional :: icount,maxijk(3)
!l  ocal
        integer :: i,species_i,n_sup,iosysinfo


        if (inode.eq.ionode) then
          call io_assign(iosysinfo)
          if (iwhat.eq.1) then
            open(unit=iosysinfo,file="sys_info.dat",action="write",status="old")
          else if (iwhat.eq.2) then
            open(unit=iosysinfo,file="sys_coord.xyz",action="write",status="replace")
          end if
          if (iwhat.eq.1) then
            write(unit=iosysinfo,fmt='(3e24.12)') rcellx,0d0,0d0
            write(unit=iosysinfo,fmt='(3e24.12)') 0d0,rcelly,0d0
            write(unit=iosysinfo,fmt='(3e24.12)') 0d0,0d0,rcellz
            if (.not.present(icount).or..not.present(maxijk)) stop ("icount/maxijk not present for iwhat=1")
            write(unit=iosysinfo,fmt='(3i22)') maxijk
            write(unit=iosysinfo,fmt='(2i22)') icount
            write(unit=iosysinfo,fmt='(2e24.12)') efermi_out
            write(unit=iosysinfo,fmt='(e24.12,A)') 0d0,"  this should be removed in the future"
            write(unit=iosysinfo,fmt='(e24.12,A)') 0d0,"  this should be removed in the future"
          end if
          write(unit=iosysinfo,fmt='(i8)') ni_in_cell
          if (iwhat.eq.2) write(unit=iosysinfo,fmt='(A)') "#  system xyz"
          do i = 1, ni_in_cell
            species_i = species_glob(i)
            n_sup = nsf_species(species_i)
            if (iwhat.eq.1) then
              write(unit=iosysinfo,fmt='(3i8,3e24.12,A64)') i,species_i,n_sup,atom_coord(1:3,i),trim(species_label(species_i))
            else if (iwhat.eq.2) then
              write(unit=iosysinfo,fmt='(A4,3e24.12)') trim(species_label(species_i)),atom_coord(1:3,i)*BohrToAng
            end if
          end do
          close(iosysinfo)
        end if



        call my_barrier()

      end subroutine get_sysinfo

      subroutine get_mat(imat,icount,maxijk)

        use global_module,   only: ni_in_cell, numprocs, nspin,io_lun,&
      &    dx_negf,dy_negf,dz_negf,spin_factor,flag_neutral_atom,ne_in_cell,&
      &    gridsolver_select,gridsolver_use,gridsolver_bc,flag_pcc_global

        use GenComms,        only: inode, ionode,gsum,gmax
        use species_module,  only: species, nsf_species
        use io_module,       only: get_file_name,dump_matrix,dump_charge
        use store_matrix,    only: dump_matrix2,dump_InfoMatGlobal,dump_matrix_update
        use input_module,    only: io_assign, io_close
        use mult_module,     only: matH, matS, matT,matK,mat_p
        use matrix_data,     only: Hrange, Srange,lrange
        use gridout_module
        use functions_on_grid,   only: atomfns, H_on_atomfns,gridfunctions, fn_on_grid
        use density_module,      only: get_electronic_density, density,density_atom,density_pcc
        use maxima_module,       only: maxngrid
        use dimens, only: n_my_grid_points
        use hartree_module, only: hartree
        use grid_solver

        implicit none
!   input :
        integer :: imat,icount,maxijk(3),ic,imd

        integer :: matrixsize,ispin,I,ierr
        integer :: lunH1,lunH2,lunS1,lunS2,lunT1,lunT2,lunH1B,lunH2B,lunS1B,lunS2B,lunT1B,lunT2B
        character(len=15) :: filenameH1, filenameH2, filenameS1, filenameS2, filenameT1, filenameT2, &
                            filenameH1B,filenameH2B,filenameS1B,filenameS2B,filenameT1B,filenameT2B
        character(256) :: pltfile
        real(double) :: electrons(nspin),elec_tot,henergy

        matrixsize = 0 !ni_in_cell*nsf
        do i = 1, ni_in_cell
          matrixsize = matrixsize + nsf_species(species(i))
        end do
        if (inode.eq.ionode) write(io_lun,*) "matrixsize ",matrixsize
        imd=0
        call write_matinfo('H',matH(1),hrange)

! write matricies in PETSc format
        icount=0
        do ispin=1,nspin

          if (ispin.eq.1) then
            call write_mat(matS(ispin),Srange,lunS1,matrixsize,ic,maxijk,"S")
            icount=max(icount,ic)
            if (inode.eq.ionode) write(io_lun,*) "S icount",icount
          end if


          if (ispin.eq.1) then
            call write_mat(matK(ispin),Hrange,lunT1,matrixsize,ic,maxijk,"K_u")
            icount=max(icount,ic)
            if (inode.eq.ionode) write(io_lun,*) "Kup icount",icount
          end if


          if (ispin.eq.2) then
            call write_mat(matK(ispin),Hrange,lunT2,matrixsize,ic,maxijk,"K_d")
            icount=max(icount,ic)
            if (inode.eq.ionode) write(io_lun,*) "Kdn icount",icount
          end if


          if (ispin.eq.1) then
            call write_mat(matH(ispin),Hrange,lunH1,matrixsize,ic,maxijk,"H_u")
            icount=max(icount,ic)
            if (inode.eq.ionode) write(io_lun,*) "Hup icount",icount
          end if

          if (ispin.eq.2) then
            call write_mat(matH(ispin),Hrange,lunH2,matrixsize,ic,maxijk,"H_d")
            icount=max(icount,ic)
            if (inode.eq.ionode) write(io_lun,*) "Hdn icount",icount
          end if

        end do

!~         write(io_lun,*) "max(icount)",icount
        call dump_InfoMatGlobal()
! ----  ------------------------------------------------------------------
! dump   grid data density and hartree potential ! assuming neutral atoms for now
        call print_grid_data(1)
        call print_grid_data(2)

        pltfile="density_0.dat"
        call grid_out_plt(pltfile,"density_0",250,density(:,1))


        pltfile="density_atom.dat"
        call grid_out_plt(pltfile,"density_atom",250,density_atom)

        if (flag_pcc_global) then
          pltfile="density_pcc.dat"
          call grid_out_plt(pltfile,"density_pcc",250,density_pcc)
        end if


        if (all( (/ dx_negf,dy_negf,dz_negf /).gt.0d0)) then

          gridfunctions(atomfns)%griddata=0d0
          call single_PAO_to_grid_set_dxyz(atomfns,1)
          call get_electronic_density(density, electrons, atomfns,     &
              H_on_atomfns(1), inode, ionode,  &
              maxngrid)
          elec_tot=sum(density(1:n_my_grid_points,1))
          call gsum(elec_tot)
          elec_tot=elec_tot*dx_negf*dy_negf*dz_negf
          if (inode.eq.ionode) write(io_lun,*) "electrons l",elec_tot
!~           density=density/elec_tot*ne_in_cell
          pltfile="density_l.dat"
          call grid_out_plt(pltfile,"density_l",250,density(:,1))
          if(flag_neutral_atom) then
            call set_atomic_density_set_dxyz(.false.,1)
            pltfile="density_atom_l.dat"
            call grid_out_plt(pltfile,"density_atom_l",250,density_atom(:))
          end if


          call single_PAO_to_grid_set_dxyz(atomfns,2)
          call get_electronic_density(density, electrons, atomfns,     &
              H_on_atomfns(1), inode, ionode,  &
              maxngrid)

          pltfile="density_r.dat"
          call grid_out_plt(pltfile,"density_r",250,density(:,1))
          if(flag_neutral_atom) then
            call set_atomic_density_set_dxyz(.false.,2)
            pltfile="density_atom_r.dat"
            call grid_out_plt(pltfile,"density_atom_r",250,density_atom(:))
          end if


          if (flag_pcc_global) then
            call set_density_pcc_set_dxyz(1)
            pltfile="density_pcc_l.dat"
            call grid_out_plt(pltfile,"density_pcc_l",250,density_pcc(:))

            call set_density_pcc_set_dxyz(2)
            pltfile="density_pcc_r.dat"
            call grid_out_plt(pltfile,"density_pcc_r",250,density_pcc(:))
          end if

        end if

      end subroutine get_mat

      subroutine write_matinfo(stub,matA,mat_range)

        use mpi
        use GenComms, ONLY: inode, ionode, cq_abort,numprocs, my_barrier
        use global_module, ONLY: numprocs, id_glob, nspin
        use io_module, ONLY: get_file_name, get_file_name_2rank
        use store_matrix, only: matrix_store, set_matrix_store

        implicit none
        character(len=*),intent(in) :: stub
        integer,intent(in) :: matA(1)
        integer,intent(in) :: mat_range
        type(matrix_store):: tmp_matrix_store

        integer :: lun, iprim, nprim, jmax, jj, ibeg, jbeta_alpha, len
        integer :: ierr,iunit,i,buflen,j
        character(256) :: outstr

        logical :: binout

        binout=.false.
        

        call set_matrix_store(stub,matA,mat_range,1,tmp_matrix_store)
        nprim=tmp_matrix_store%n_prim

        call MPI_File_delete("mat_info.dat",MPI_INFO_NULL,ierr)
        call MPI_FILE_OPEN(MPI_COMM_WORLD, "mat_info.dat",MPI_MODE_WRONLY + MPI_MODE_CREATE, MPI_INFO_NULL, iunit, ierr)

        if (binout) then
          if (inode.eq.ionode) then

            buflen=1
            call MPI_FILE_WRITE_SHARED(iunit,numprocs,1,MPI_CHARACTER,MPI_STATUS_IGNORE,ierr)
          end if

          call my_barrier()

          do i=1,numprocs

            if (i.eq.inode) then
              buflen=1
              call MPI_FILE_WRITE_SHARED(iunit,inode,buflen,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
              call MPI_FILE_WRITE_SHARED(iunit,nprim,buflen,MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
              buflen=size(tmp_matrix_store%nsf_spec_i(1:nprim))
              call MPI_FILE_WRITE_SHARED(iunit,tmp_matrix_store%nsf_spec_i(1:nprim),buflen,&
            &MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
              buflen=size(tmp_matrix_store%jmax_i(1:nprim))
              call MPI_FILE_WRITE_SHARED(iunit,tmp_matrix_store%jmax_i(1:nprim),buflen,&
            &MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
              buflen=size(tmp_matrix_store%jbeta_max_i(1:nprim))
              call MPI_FILE_WRITE_SHARED(iunit,tmp_matrix_store%jbeta_max_i(1:nprim),buflen,&
            &MPI_INTEGER,MPI_STATUS_IGNORE,ierr)

              if(nprim.gt.0) then
                do iprim=1,nprim
                    jmax = tmp_matrix_store%jmax_i(iprim)
                    ibeg = tmp_matrix_store%ibeg_Rij(iprim)
                    buflen=1
                    call MPI_FILE_WRITE_SHARED(iunit,tmp_matrix_store%idglob_i(iprim),buflen,MPI_INTEGER,&
                  &MPI_STATUS_IGNORE,ierr)
                    buflen=size(tmp_matrix_store%beta_j(ibeg:ibeg+jmax-1))
                    call MPI_FILE_WRITE_SHARED(iunit,tmp_matrix_store%beta_j(ibeg:ibeg+jmax-1),buflen,&
                  &MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
                    buflen=size(tmp_matrix_store%idglob_j(ibeg:ibeg+jmax-1))
                    call MPI_FILE_WRITE_SHARED(iunit,tmp_matrix_store%idglob_j(ibeg:ibeg+jmax-1),buflen,&
                  &MPI_INTEGER,MPI_STATUS_IGNORE,ierr)
                  do jj=1,jmax
                    buflen=size(tmp_matrix_store%vec_Rij(1:3,ibeg+jj-1))
                    call MPI_FILE_WRITE_SHARED(iunit,tmp_matrix_store%vec_Rij(1:3,ibeg+jj-1),buflen,&
                  &MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,ierr)
                  enddo !jj=1,jmax

                enddo !iprim=1,nprim
              endif  ! (nprim .GT. 0)


            end if

            call my_barrier()

          end do
        else

          if (inode.eq.ionode) then
            write(outstr,fmt='(i8,A)') numprocs,new_line('A')
            buflen=len(trim(outstr))
            call MPI_FILE_WRITE_SHARED(iunit,trim(outstr),buflen,MPI_CHARACTER,MPI_STATUS_IGNORE,ierr)
          end if

          call my_barrier()

          ! set_matrix_store : build tmp_matrix_store

  !   I don't know the actual offset per file, probably this could be worked out but as the
  !   ammount of actual data to be written is quite small, lets do this process by process
  !   with a shared file pointer.
          do i=1,numprocs

            if (i.eq.inode) then

              write(outstr,fmt='(2i8,A)') inode, nprim,new_line('A')
              buflen=len(trim(outstr))
              call MPI_FILE_WRITE_SHARED(iunit,trim(outstr),buflen,MPI_CHARACTER,MPI_STATUS_IGNORE,ierr)

              do j=1,nprim
                if (mod(j,7).eq.0) then
                  write(outstr,fmt='(A)') new_line('A')
                  buflen=len(trim(outstr))
                  call MPI_FILE_WRITE_SHARED(iunit,trim(outstr),buflen,MPI_CHARACTER,MPI_STATUS_IGNORE,ierr)
                end if
                write(outstr,fmt='(i8)') tmp_matrix_store%nsf_spec_i(j)
                buflen=len(trim(outstr))
                call MPI_FILE_WRITE_SHARED(iunit,trim(outstr),buflen,MPI_CHARACTER,MPI_STATUS_IGNORE,ierr)
              end do
              write(outstr,fmt='(A)') new_line('A')
              buflen=len(trim(outstr))
              call MPI_FILE_WRITE_SHARED(iunit,trim(outstr),buflen,MPI_CHARACTER,MPI_STATUS_IGNORE,ierr)

              do j=1,nprim
                if (mod(j,7).eq.0) then
                  write(outstr,fmt='(A)') new_line('A')
                  buflen=len(trim(outstr))
                  call MPI_FILE_WRITE_SHARED(iunit,trim(outstr),buflen,MPI_CHARACTER,MPI_STATUS_IGNORE,ierr)
                end if
                write(outstr,fmt='(i8)') tmp_matrix_store%jmax_i(j)
                buflen=len(trim(outstr))
                call MPI_FILE_WRITE_SHARED(iunit,trim(outstr),buflen,MPI_CHARACTER,MPI_STATUS_IGNORE,ierr)
              end do
              write(outstr,fmt='(A)') new_line('A')
              buflen=len(trim(outstr))
              call MPI_FILE_WRITE_SHARED(iunit,trim(outstr),buflen,MPI_CHARACTER,MPI_STATUS_IGNORE,ierr)

              do j=1,nprim
                if (mod(j,7).eq.0) then
                  write(outstr,fmt='(A)') new_line('A')
                  buflen=len(trim(outstr))
                  call MPI_FILE_WRITE_SHARED(iunit,trim(outstr),buflen,MPI_CHARACTER,MPI_STATUS_IGNORE,ierr)
                end if
                write(outstr,fmt='(i8)') tmp_matrix_store%jbeta_max_i(j)
                buflen=len(trim(outstr))
                call MPI_FILE_WRITE_SHARED(iunit,trim(outstr),buflen,MPI_CHARACTER,MPI_STATUS_IGNORE,ierr)
              end do
              write(outstr,fmt='(A)') new_line('A')
              buflen=len(trim(outstr))
              call MPI_FILE_WRITE_SHARED(iunit,trim(outstr),buflen,MPI_CHARACTER,MPI_STATUS_IGNORE,ierr)
              if(nprim .gt. 0) then
                do iprim=1,nprim
                  jmax = tmp_matrix_store%jmax_i(iprim)
                  ibeg = tmp_matrix_store%ibeg_Rij(iprim)

                  write(outstr,fmt='(i8,A)') tmp_matrix_store%idglob_i(iprim),new_line('A')
                  buflen=len(trim(outstr))
                  call MPI_FILE_WRITE_SHARED(iunit,trim(outstr),buflen,MPI_CHARACTER,MPI_STATUS_IGNORE,ierr)

                  do j=1,jmax
                    if (mod(j,7).eq.0) then
                      write(outstr,fmt='(A)') new_line('A')
                      buflen=len(trim(outstr))
                      call MPI_FILE_WRITE_SHARED(iunit,trim(outstr),buflen,MPI_CHARACTER,MPI_STATUS_IGNORE,ierr)
                    end if
                    write(outstr,fmt='(i8)') tmp_matrix_store%beta_j(ibeg+j-1)
                    buflen=len(trim(outstr))
                    call MPI_FILE_WRITE_SHARED(iunit,trim(outstr),buflen,MPI_CHARACTER,MPI_STATUS_IGNORE,ierr)
                  end do

                  write(outstr,fmt='(A)') new_line('A')
                  buflen=len(trim(outstr))
                  call MPI_FILE_WRITE_SHARED(iunit,trim(outstr),buflen,MPI_CHARACTER,MPI_STATUS_IGNORE,ierr)

                  do j=1,jmax
                    if (mod(j,7).eq.0) then
                      write(outstr,fmt='(A)') new_line('A')
                      buflen=len(trim(outstr))
                      call MPI_FILE_WRITE_SHARED(iunit,trim(outstr),buflen,MPI_CHARACTER,MPI_STATUS_IGNORE,ierr)
                    end if
                    write(outstr,fmt='(i8)') tmp_matrix_store%idglob_j(ibeg+j-1)
                    buflen=len(trim(outstr))
                    call MPI_FILE_WRITE_SHARED(iunit,trim(outstr),buflen,MPI_CHARACTER,MPI_STATUS_IGNORE,ierr)
                  end do

                  write(outstr,fmt='(A)') new_line('A')
                  buflen=len(trim(outstr))
                  call MPI_FILE_WRITE_SHARED(iunit,trim(outstr),buflen,MPI_CHARACTER,MPI_STATUS_IGNORE,ierr)

                  do jj=1,jmax
                    write(outstr,fmt='(3e24.12,A)') tmp_matrix_store%vec_Rij(1:3,ibeg+jj-1),new_line('A')
                    buflen=len(trim(outstr))
                    call MPI_FILE_WRITE_SHARED(iunit,trim(outstr),buflen,MPI_CHARACTER,MPI_STATUS_IGNORE,ierr)
                  enddo !jj=1,jmax

                end do !iprim=1,nprim
              end if  ! (nprim .GT. 0)


            end if
            call my_barrier()
          end do
        end if

        call MPI_File_Close(iunit,ierr)
      end subroutine write_matinfo


      subroutine write_mat(matA,Arange,lun,matrix_size,icount,maxijk,fname,ladd)
  ! ----------------------------------------------------------------------
  !
  !   Write out matrix values with indeces
  !   (primary atom, neighbour atom, sf1 in primary atom, sf2 in neighbour atom,
  !    which are all in global numbers).
  !

  ! ----------------------------------------------------------------------
  ! (modfied version based on Nakata-san's source)
  ! Dump matricies into PETSc format (AIJ storage) to load them with MatLoad
  ! http://www.mcs.anl.gov/petsc/petsc-current/docs/manualpages/Mat/MatLoad.html
  !
  ! 1. First run: get non-zero structure
  ! 2. Second run: Write matrix to disk, seperate file for different cells
  ! ----------------------------------------------------------------------


      use datatypes
      use numbers,             ONLY: zero
      use global_module,       ONLY: id_glob, species_glob, numprocs, ni_in_cell, &
                                     x_atom_cell,y_atom_cell,z_atom_cell,io_lun,sf
      use species_module,      ONLY: species, nsf_species, npao_species
      use GenComms,            ONLY: myid, cq_abort, gsum, my_barrier, inode, ionode,numprocs,gmax,gsum,int_copy
  !    use ScalapackFormat,     ONLY: pgid, pg_kpoints
      use matrix_module,       ONLY: matrix, matrix_halo
      use group_module,        ONLY: parts
      use primary_module,      ONLY: bundle
      use cover_module,        ONLY: BCS_parts
      use support_spec_format, ONLY: flag_paos_atoms_in_cell, flag_one_to_one
      use matrix_data,         ONLY: mat, halo,mat_name
      use mult_module,         ONLY: matrix_pos, mat_p
      use input_module,        ONLY: io_assign, io_close
      use timer_module
      use dimens,         only: r_super_x, r_super_y, r_super_z
      use functions_on_grid,           only: atomfns


      implicit none

      ! Passned variables
      integer, intent(in) :: matA, Arange, matrix_size
      integer, intent(out) :: icount,maxijk(3)
      integer :: lun, lunb
      character(*) :: fname
      logical, optional :: ladd




      ! Local variables
      integer :: iprim, np, i, j, ist, atom_num, atom_i, atom_spec, nprim
      integer :: neigh_global_part, atom_j, neigh_species, j_in_halo
      integer :: nsf1, nsf2, nsf, nsf0, max_nab, last_m, &
                 k, l, kl, kg, lg, kstart, lstart, jj,ip,itmp
      integer :: gcspart
      integer :: wheremat
      integer, parameter :: Ilast=9999999
      real(double) :: dx, dy, dz, r2
      real(double) :: rval, phase, rfac, ifac
      integer, allocatable, dimension(:) :: label_nsf_g
      integer, allocatable, dimension(:,:) :: label, start_m
      complex(double_cplx) :: cval
      integer, allocatable :: Kmat(:), Lmat(:),iunit(:,:,:),isort(:)
      integer, allocatable :: nzrow(:,:,:,:),zcount(:,:,:),nzoff(:,:,:,:),icoloff(:,:,:,:)
      integer :: i1,i2,i3,irun,irow,icol,nbatch,ii,ierr
      integer :: lunLBL
      integer :: FSCpart
      integer :: maxi,iheaderoff,ipos1,ipos2
      character(256),allocatable :: petscfile(:,:,:)
      character(256) :: si1,si2,si3
      complex(double_cplx) :: zout,zin
      integer :: ikind
      logical :: lreadandadd
      
      lreadandadd=.false.
      if (present(ladd)) lreadandadd=ladd
      
      ikind=kind(ierr)

!      write(0,*) "start write",trim(fname)
!~       if (inode.eq.ionode) write(io_lun,*) "sf",sf

  !   --- make label_nsf_g of starting positions for global atoms ---

      allocate(label_nsf_g(ni_in_cell),STAT=ierr)
      label_nsf_g(:) = 0

      max_nab = 0
      iprim = 0; atom_i = 0
      do np = 1,bundle%groups_on_node ! Loop over primary set partitions
         if(bundle%nm_nodgroup(np)>0) then ! If there are atoms in partition
            do i = 1,bundle%nm_nodgroup(np) ! Loop over atoms
               atom_num = bundle%nm_nodbeg(np) + i - 1
               iprim = iprim + 1
               atom_i = bundle%ig_prim(iprim)                 ! global number of i
               atom_spec = bundle%species(atom_num)           ! atomic species of i
               label_nsf_g(atom_i) = nsf_species(atom_spec)   ! number of SFs on i
               max_nab = max(max_nab,mat(np,Arange)%n_nab(i)) ! maximum number of n_nab
            enddo
         endif
      enddo

      call my_barrier()
      call gsum(label_nsf_g,ni_in_cell)
      call my_barrier()

      nsf0 = 1
      nsf  = 1
      do i = 1, ni_in_cell
         nsf0 = nsf
         nsf  = nsf + label_nsf_g(i)
         label_nsf_g(i) = nsf0         ! starting position in global square matrix
      enddo
      nsf = nsf - 1
      if (nsf.ne.matrix_size) then
         write(io_lun,*) "fname ",trim(fname)
         write(io_lun,*) 'nsf=',nsf, '  matrix_size=',matrix_size
         call cq_abort("Error!! Size of square matrix is not correct.")
      endif
      if (matrix_size.gt.1000000000) &
         call cq_abort("Error!! matrix_size is larger than 10^10.")





      nprim = bundle%n_prim
      allocate(label(max_nab,nprim),STAT=ierr)
      allocate(start_m(max_nab,nprim),STAT=ierr)
      label(:,:) = 0
      start_m(:,:) = 0

      last_m = 0
      iprim = 0; atom_i = 0
      maxijk=0
!~       write(lun,fmt='(i22)') sum(bundle%nm_nodgroup(1:bundle%groups_on_node))
      do np = 1,bundle%groups_on_node ! Loop over primary set partitions
         if(bundle%nm_nodgroup(np)>0) then ! If there are atoms in partition
            do i = 1,bundle%nm_nodgroup(np) ! Loop over atoms
               atom_num = bundle%nm_nodbeg(np) + i - 1
               iprim = iprim + 1
               atom_i = bundle%ig_prim(iprim)         ! global number of i
               atom_spec = bundle%species(atom_num)   ! atomic species of i
               nsf1 = nsf_species(atom_spec)          ! number of SFs on i


               do j = 1, mat(np,Arange)%n_nab(i) ! Loop over neighbours
                  ist = mat(np,Arange)%i_acc(i) + j - 1
                  gcspart = BCS_parts%icover_ibeg(mat(np,Arange)%i_part(ist))+mat(np,Arange)%i_seq(ist)-1
                  neigh_global_part = BCS_parts%lab_cell(mat(np,Arange)%i_part(ist))
                  atom_j = id_glob(parts%icell_beg(neigh_global_part)+mat(np,Arange)%i_seq(ist)-1)   ! global number of j
                  neigh_species = species_glob(atom_j)                                               ! atomic species of j
                  nsf2 = nsf_species(neigh_species)                                                  ! number of SFs on j
                  FSCpart = BCS_parts%lab_cell(mat(np,Arange)%i_part(ist))!gcspart)
      ! - this is the original definition used during the SCF. The phase depends on the actual atomic position
      !                    dx = BCS_parts%xcover(gcspart)-bundle%xprim(atom_num)
      !                    dy = BCS_parts%ycover(gcspart)-bundle%yprim(atom_num)
      !                    dz = BCS_parts%zcover(gcspart)-bundle%zprim(atom_num)
      !--------------------------------------------------------------------------------------------------------
      ! - this is the usual definition of real-space vs. k-space, that is in terms of real-space units cells.
                  dx = BCS_parts%xcover(gcspart)-x_atom_cell(parts%icell_beg(FSCpart)+mat(np,Arange)%i_seq(ist)-1)
                  dy = BCS_parts%ycover(gcspart)-y_atom_cell(parts%icell_beg(FSCpart)+mat(np,Arange)%i_seq(ist)-1)
                  dz = BCS_parts%zcover(gcspart)-z_atom_cell(parts%icell_beg(FSCpart)+mat(np,Arange)%i_seq(ist)-1)
      !------------------------------------------------------------------------------------------------------------

                  i1=nint(dx/r_super_x)
                  i2=nint(dy/r_super_y)
                  i3=nint(dz/r_super_z)
                  maxijk(1)=max(maxijk(1),abs(i1))
                  maxijk(2)=max(maxijk(2),abs(i2))
                  maxijk(3)=max(maxijk(3),abs(i3))

                  do jj = 1,j
                     if (label(jj,iprim).eq.atom_j) then
                        exit
                     else if (label(jj,iprim).eq.0) then
                        label(jj,iprim) = atom_j               ! neighbor list in global number
                        start_m(jj,iprim) = last_m + 1         ! starting position in local square matrix
                        last_m = last_m + nsf1*nsf2
                        exit
                     endif
                  enddo ! jj

               end do ! j

            end do ! i
         end if ! End if nm_nodgroup > 0
      end do ! np


      call gmax(maxijk(1))
      call gmax(maxijk(2))
      call gmax(maxijk(3))
      write(io_lun,fmt='(A,3i8)') "maxijk new ",maxijk
      write(io_lun,fmt='(A,A)') "matname ",trim(mat_name(Arange))      

      allocate(nzrow(matrix_size,-maxijk(1):maxijk(1),-maxijk(2):maxijk(2),-maxijk(3):maxijk(3)),&
     &  nzoff(matrix_size,-maxijk(1):maxijk(1),-maxijk(2):maxijk(2),-maxijk(3):maxijk(3)),&
     &  icoloff(matrix_size,-maxijk(1):maxijk(1),-maxijk(2):maxijk(2),-maxijk(3):maxijk(3)),&
     &  zcount(-maxijk(1):maxijk(1),-maxijk(2):maxijk(2),-maxijk(3):maxijk(3)),&
     &  iunit(-maxijk(1):maxijk(1),-maxijk(2):maxijk(2),-maxijk(3):maxijk(3)),&
     &  petscfile(-maxijk(1):maxijk(1),-maxijk(2):maxijk(2),-maxijk(3):maxijk(3)),&
     &  stat=ierr)
      if (ierr.ne.0) then
        write(0,*) "allocation error ",ierr
      end if


      zcount=0
      icoloff=0
      nzrow=0
      nzoff=0

      allocate(Kmat(last_m),STAT=ierr)
      allocate(Lmat(last_m),STAT=ierr)
!~       allocate(isort(last_m),STAT=ierr)
!      write(6,*) "last_m ",last_m,nsf1

      Kmat(:) = 0
      Lmat(:) = 0
!      isort=-1

      icount = 0

      iprim = 0; atom_i = 0
      do np = 1,bundle%groups_on_node ! Loop over primary set partitions
         if(bundle%nm_nodgroup(np)>0) then ! If there are atoms in partition
            do i = 1,bundle%nm_nodgroup(np) ! Loop over atoms
               atom_num  = bundle%nm_nodbeg(np) + i - 1
               iprim     = iprim + 1
               atom_i    = bundle%ig_prim(iprim)      ! global number of i
               atom_spec = bundle%species(atom_num)   ! atomic species of i
               nsf1      = nsf_species(atom_spec)     ! number of SFs on i
               kstart    = label_nsf_g(atom_i)        ! starting position in global square matrix for i

               do jj = 1, mat(np,Arange)%n_nab(i) ! Loop over neighbours (global)
                  if (label(jj,iprim).eq.0) exit
                  atom_j = label(jj,iprim)                   ! global number of j
                  neigh_species = species_glob(atom_j)   ! atomic species of j
                  nsf2 = nsf_species(neigh_species)      ! number of SFs on j
                  lstart = label_nsf_g(atom_j)           ! starting position in global square matrix for j

                  kl = start_m(jj,iprim) - 1
                  do k = 1, nsf1
                     kg = kstart + k - 1
                     do l = 1, nsf2
                        lg = lstart + l - 1
                        kl = kl + 1
                        icount = icount + 1
                        Kmat(kl) = kg
                        Lmat(kl) = lg
!~                         isort(kg)=kl
                     enddo ! l
                  enddo ! k

               enddo ! jj
            enddo ! i
         endif
      enddo ! np


      ! for now write out all cells. Later one should avoid the unecessary terms which can be obtained
      ! by translational symmetry X(-R)=X(R)
      do i3=-maxijk(3),maxijk(3)
        write(si3,*) i3
        si3=adjustl(si3)
        do i2=-maxijk(2),maxijk(2)
          write(si2,*) i2
          si2=adjustl(si2)
          do i1=-maxijk(1),maxijk(1)
            write(si1,*) i1
            si1=adjustl(si1)
            petscfile(i1,i2,i3)=trim(fname)//"_"//trim(si1)//"_"//trim(si2)&
           &  //"_"//trim(si3)//"_petsc.dat"
          end do
        end do
      end do

      if (inode.eq.ionode) then
        ! using fortran 2003 standard: access="stream", read/write controlled by pos=x where pos is here given in byte
        ! -> integer(4) -> 4 byte, double precision -> 8 byte, double complex -> 16 byte
        ! PETSc MAT_LOAD need ordered nonzero columns, that is very cumbersome to work out here. -> Use MPI_READs in transport code
        ! to load matrices and work out process partition. No BIG_ENDIAN necessary as for MAT_LOAD
        if (.not.lreadandadd) then
          do i3=-maxijk(3),maxijk(3)
            do i2=-maxijk(2),maxijk(2)
              do i1=-maxijk(1),maxijk(1)
                open(newunit=iunit(i1,i2,i3),file=trim(petscfile(i1,i2,i3)),access="stream",&
              & status="replace") !,CONVERT='BIG_ENDIAN')
                close(iunit(i1,i2,i3))
              end do
            end do
          end do
        end if
        ! -----
      end if

      call mpi_barrier(mpi_comm_world,ierr)

      do i3=-maxijk(3),maxijk(3)
        do i2=-maxijk(2),maxijk(2)
          do i1=-maxijk(1),maxijk(1)
            open(newunit=iunit(i1,i2,i3),file=trim(petscfile(i1,i2,i3)),access="stream",&
           & status="old") !,CONVERT='BIG_ENDIAN')
          end do
        end do
      end do

      do irun=1,2
        icount=0
        do ip=1,numprocs
          if (inode.eq.ip) then
            maxi=0
            iprim = 0; atom_i = 0
            do np = 1,bundle%groups_on_node ! Loop over primary set partitions
               if(bundle%nm_nodgroup(np)>0) then ! If there are atoms in partition
                  do i = 1,bundle%nm_nodgroup(np) ! Loop over atoms
                     atom_num = bundle%nm_nodbeg(np) + i - 1
                     iprim = iprim + 1
                     atom_i = bundle%ig_prim(iprim)         ! global number of i
                     atom_spec = bundle%species(atom_num)   ! atomic species of i
                     nsf1 = nsf_species(atom_spec)          ! number of SFs on i

                     do j = 1, mat(np,Arange)%n_nab(i) ! Loop over neighbours
                        ist = mat(np,Arange)%i_acc(i) + j - 1
                        ! Build the distances between atoms
                        gcspart = BCS_parts%icover_ibeg(mat(np,Arange)%i_part(ist))+mat(np,Arange)%i_seq(ist)-1
                        FSCpart = BCS_parts%lab_cell(mat(np,Arange)%i_part(ist))!gcspart)
                        ! Displacement vector

      ! - this is the original definition used during the SCF. The phase depends on the actual atomic position
      !                    dx = BCS_parts%xcover(gcspart)-bundle%xprim(atom_num)
      !                    dy = BCS_parts%ycover(gcspart)-bundle%yprim(atom_num)
      !                    dz = BCS_parts%zcover(gcspart)-bundle%zprim(atom_num)
      !--------------------------------------------------------------------------------------------------------
      ! - this is the usual definition of real-space vs. k-space, that is in terms of real-space units cells.
                        dx = BCS_parts%xcover(gcspart)-x_atom_cell(parts%icell_beg(FSCpart)+mat(np,Arange)%i_seq(ist)-1)
                        dy = BCS_parts%ycover(gcspart)-y_atom_cell(parts%icell_beg(FSCpart)+mat(np,Arange)%i_seq(ist)-1)
                        dz = BCS_parts%zcover(gcspart)-z_atom_cell(parts%icell_beg(FSCpart)+mat(np,Arange)%i_seq(ist)-1)
      !------------------------------------------------------------------------------------------------------------

                        i1=nint(dx/r_super_x)
                        i2=nint(dy/r_super_y)
                        i3=nint(dz/r_super_z)


                        neigh_global_part = BCS_parts%lab_cell(mat(np,Arange)%i_part(ist))
                        atom_j = id_glob(parts%icell_beg(neigh_global_part)+mat(np,Arange)%i_seq(ist)-1) ! global number of j
                        j_in_halo = halo(Arange)%i_halo(gcspart)                                         ! halo-ID of j
                        neigh_species = species_glob(atom_j)                                             ! atomic species of j
                        nsf2 = nsf_species(neigh_species)                                                ! number of SFs on j
                        ! find starting position
                        kl = 0
                        do jj = 1,j
                           if (label(jj,iprim).eq.atom_j) then
                              kl = start_m(jj,iprim)               ! starting position of elements for j
                              exit
                           endif
                        enddo ! jj
                        if (kl.eq.0) call cq_abort("Fail to find starting pojition for neighbor atom",j)
                        kl = kl - 1
                        do k = 1, nsf1
                           do l = 1, nsf2
                              kl = kl + 1
                              wheremat = matrix_pos(matA,iprim,j_in_halo,k,l)
                              rval = mat_p(matA)%matrix(wheremat)
                              icount=icount+1

!                              if (irun.eq.2) write(unit=lun,fmt='(5I7,es45.24e5,5i8)') i1,i2,i3, Kmat(kl),Lmat(kl),rval,atom_i,&
!                                             atom_j,mat(np,Arange)%n_nab(i),inode,wheremat
!~                               if (rval.eq.0d0) cycle
                              if (irun.eq.1) then
                                nzrow(Kmat(kl),i1,i2,i3)=nzrow(Kmat(kl),i1,i2,i3)+1
                                zcount(i1,i2,i3)=zcount(i1,i2,i3)+1
                              end if
                              if (irun.eq.2) then

                                irow=Kmat(kl)
                                icol=Lmat(kl)
                                icoloff(irow,i1,i2,i3)=icoloff(irow,i1,i2,i3)+1
                                ipos1=(iheaderoff+(nzoff(irow,i1,i2,i3)+icoloff(irow,i1,i2,i3)-1))*ikind+1
                                ipos2=(iheaderoff+zcount(i1,i2,i3))*ikind+(nzoff(irow,i1,i2,i3)+icoloff(irow,i1,i2,i3)-1)*16+1
                                zout=rval
!~
!                                write(6,fmt='(i4,i9,7i9,A,4i12,2e24.12)') inode,zcount(i1,i2,i3),i1,i2,i3,(iheaderoff+zcount(i1,i2,i3)-1)*4,&
!                                &iheaderoff,nzoff(irow,i1,i2,i3),icoloff(irow,i1,i2,i3)," "//trim(fname)//" ",ipos1,ipos2,irow,icol,zout
!                                call flush(6)
                                write(iunit(i1,i2,i3),pos=ipos1) icol-1 !int(icol-1,8)       ! C style arrays, i.e., index starts at 0
                                if (lreadandadd) then
                                  read(iunit(i1,i2,i3),pos=ipos2) zin
                                  zout=zout+zin                                
                                end if
                                write(iunit(i1,i2,i3),pos=ipos2) zout 
                              end if
                           enddo ! l
                        enddo ! k
                     end do ! j
                  end do ! i
               end if ! End if nm_nodgroup > 0
            end do ! np


!~             write(6,fmt='(A,6i8)') "nz",irun,inode,zcount
!~             call flush(6)
          end if


          call mpi_barrier(mpi_comm_world,ierr)



          do i3=-maxijk(3),maxijk(3)
            do i2=-maxijk(2),maxijk(2)
              do i1=-maxijk(1),maxijk(1)
                if (irun.eq.1) call MPI_bcast(zcount(i1,i2,i3), 1, MPI_Integer, ip-1,MPI_comm_world, ierr )
                if (irun.eq.2) call MPI_bcast(icoloff(1:matrix_size,i1,i2,i3), matrix_size, MPI_Integer,ip-1, MPI_comm_world, ierr )
              end do
            end do
          end do

        end do !ip

        do i3=-maxijk(3),maxijk(3)
          do i2=-maxijk(2),maxijk(2)
            do i1=-maxijk(1),maxijk(1)
              do ii=1,matrix_size
                call gsum(nzrow(ii,i1,i2,i3))
              end do
            end do
          end do
        end do


        call my_barrier()

        if ((inode.eq.ionode).and.(irun.eq.1)) then

!~           write(io_lun,*) "matrix_size",matrix_size
!~           write(io_lun,*) "nz",zcount
!~           write(io_lun,*) "nzrow",nzrow
          nzoff(1,:,:,:)=0

          do ii=2,matrix_size
            nzoff(ii,:,:,:)=nzoff(ii-1,:,:,:)+nzrow(ii-1,:,:,:)
          end do

          iheaderoff=4+matrix_size ! 4 is 121126,matrix_size,matrix_size,zcount(i1,i2,i3). nzrow is added to nzoff one by one

!~           write(io_lun,*) "nzoff",nzoff


!~    int    MAT_FILE_CLASSID
!~    int    number of rows
!~    int    number of columns
!~    int    total number of nonzeros
!~    int    *number nonzeros in each row
!~    int    *column indices of all nonzeros (starting index is zero)
!~    PetscScalar *values of all nonzeros

!~ for the time being use integer(4) in transport (and petsc for that matter)

          do i3=-maxijk(3),maxijk(3)
            do i2=-maxijk(2),maxijk(2)
              do i1=-maxijk(1),maxijk(1)
                rewind(iunit(i1,i2,i3))
                write(iunit(i1,i2,i3)) 1211216 !int(1211216,8) 0
                write(iunit(i1,i2,i3)) matrix_size !int(matrix_size,8) 5
                write(iunit(i1,i2,i3)) matrix_size !int(matrix_size,8) 10
                write(iunit(i1,i2,i3)) zcount(i1,i2,i3) !int(zcount(i1,i2,i3),8) 14
                write(iunit(i1,i2,i3)) nzrow(1:matrix_size,i1,i2,i3) !int(nzrow(1:matrix_size,i1,i2,i3),8) 17
              end do
            end do
          end do


        end if

        if (irun.eq.1) then
          call MPI_bcast( iheaderoff, 1, MPI_Integer, ionode-1, &
                MPI_comm_world, ierr )
          nbatch=(maxijk(1)*2+1)*(maxijk(2)*2+1)*(maxijk(3)*2+1)*matrix_size
          call MPI_bcast( nzoff, nbatch, MPI_Integer,ionode-1, &
                MPI_comm_world, ierr )
        end if

      end do !irun


      call mpi_barrier(mpi_comm_world,ierr)

      call gsum(icount)

      do i3=-maxijk(3),maxijk(3)
        do i2=-maxijk(2),maxijk(2)
          do i1=-maxijk(1),maxijk(1)
            close(iunit(i1,i2,i3))
          end do
        end do
      end do

      deallocate(Kmat, Lmat,nzoff,icoloff,nzrow,stat=ierr)
      if (ierr.ne.0) then
        write(0,*) "deallocation error ",ierr
        stop
      end if

    end subroutine write_mat

    subroutine update_pulay_history_negf(iPulay, rho_in, rho_out, &
                                    rho_pul, resid_pul, k_resid_pul,     &
                                    cov_resid_pul)
      use datatypes
      use GenBlas
      use dimens,         only: n_my_grid_points, grid_point_volume
      use GenComms,       only: gsum, cq_abort
      use hartree_module, only: kerker, kerker_and_wdmetric, wdmetric
      use global_module,  only: nspin, flag_fix_spin_population
      use maxima_module,  only: maxngrid
      use memory_module,  only: reg_alloc_mem, reg_dealloc_mem, type_dbl
      use SelfCon,             only:q0,flag_kerker,A,flag_wdmetric,q1,maxpulaySC

      implicit none

      ! Passed parameters
      integer,      intent(in)  :: iPulay

      real(double), dimension(maxngrid,nspin),             intent(in)    :: &
          rho_in,rho_out
      real(double), dimension(maxngrid,maxpulaySC,nspin),  intent(inout) :: &
          rho_pul, resid_pul, k_resid_pul, cov_resid_pul

      ! local variables
      integer :: spin, stat
      real(double), dimension(:,:), allocatable ::  resid



      allocate(resid(maxngrid,nspin), STAT=stat)
      if (stat /= 0) &
          call cq_abort("update_pulay_history: Error alloc mem: ", maxngrid, nspin)


      ! store rho_in to history
      do spin = 1, nspin
        rho_pul(1:n_my_grid_points,iPulay,spin) = &
              rho_in(1:n_my_grid_points,spin)
      end do

      ! calculate residual of rho_in
      resid = zero

      do spin = 1, nspin
        resid(1:n_my_grid_points,spin) = &
              rho_out(1:n_my_grid_points,spin) - &
              rho_in(1:n_my_grid_points,spin)
        ! replace the resid_pul history at iPulay to new value
        resid_pul(1:n_my_grid_points,iPulay,spin) = &
              resid(1:n_my_grid_points,spin)
      end do

      ! calculate the Kerker preconditioned or wave-dependent metric
      ! covarient residuals if requires to
      do spin = 1, nspin
        if (flag_Kerker) then
            if (flag_wdmetric) then
              call kerker_and_wdmetric(resid(:,spin), &
                                        cov_resid_pul(:,iPulay,spin), &
                                        maxngrid, q0, q1)
            else
              call kerker(resid(:,spin), maxngrid, q0)
            end if
            ! replace the k_resid_pul history at iPulay to new value
            k_resid_pul(1:n_my_grid_points,iPulay,spin) = &
                resid(1:n_my_grid_points,spin)
        else
            if (flag_wdmetric) then
              call wdmetric(resid(:,spin), cov_resid_pul(:,iPulay,spin), &
                            maxngrid, q1)
            end if
        end if ! flag_Kerker
      end do ! spin

      deallocate(resid, STAT=stat)
      if (stat /= 0) &
          call cq_abort("update_pulay_history: Error dealloc mem")

    end subroutine update_pulay_history_negf

end module negf_module
