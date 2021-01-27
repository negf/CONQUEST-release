module negf_output_module
  implicit none
  
  contains
  

     subroutine write_mat(matA,Arange,maxijk,fname,ladd)
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
      
      use MPI


      implicit none
      
      type mat_aij
        integer, allocatable :: cols(:), offsets(:)
      end type mat_aij

      ! Passned variables
      integer, intent(in) :: matA, Arange
      integer, intent(out) :: maxijk(3)      
      character(*) :: fname
      logical, optional :: ladd




      ! Local variables
      
      type(mat_aij), allocatable :: mat_info(:,:,:,:)
      
      integer :: iprim, np, i, j, ist, atom_num, atom_i, atom_spec, nprim
      integer :: neigh_global_part, atom_j, neigh_species, j_in_halo
      integer :: nsf1, nsf2, nsf, nsf0, max_nab, last_m, &
                 k, l, kl, kg, lg, kstart, lstart, jj,ip,itmp
      integer :: gcspart, matrix_size
      integer :: wheremat
      integer, parameter :: Ilast=9999999
      real(double) :: dx, dy, dz, r2
      real(double) :: rval, phase, rfac, ifac
      real(double), allocatable :: values_local(:), values_global(:)
            
      integer, allocatable ::  cols_global(:), cols_offset(:), cols_global_idx(:)
      
      integer, allocatable :: label_nsf_g(:)
      integer, allocatable :: label(:,:), start_m(:,:)
      complex(double_cplx) :: cval
      complex(double_cplx), allocatable :: buf_out(:), buf_in(:)
      integer, allocatable :: iunit(:,:,:), IamMat(:,:), displ(:), nz_receive(:), at_npaos(:)
      integer, allocatable :: nzrow(:,:,:,:), nzcount(:,:,:), nzrow_local(:,:,:,:), pos_out(:), &
        pos_col_out(:), buf_col_out(:)
      integer :: i1,i2,i3,irun,irow,icol,nbatch,ii,ierr, nz_l, nz_g, nzoff, max_nsf
      integer :: lun
      integer :: FSCpart,iii(1)
      integer :: maxi,iheaderoff,ipos, kkl, ipos_col, ipos_val, nzrow_i, row_offset, header_offset
      character(256),allocatable :: petscfile(:,:,:)
      character(256) :: si1,si2,si3
      complex(double_cplx) :: zout,zin
      integer :: isize, zsize
      logical :: lreadandadd, ldebug
      integer(8) :: counti, count_rate, countf, i_offset
      
      ldebug = .false.
      lreadandadd=.false.
      if (present(ladd)) lreadandadd=ladd
      
      isize=kind(ierr)
      zsize=kind(cval)*2

     
      matrix_size = 0 !ni_in_cell*nsf
      do i = 1, ni_in_cell
        matrix_size = matrix_size + nsf_species(species(i))
      end do
      max_nsf = maxval(nsf_species)
      
      write(io_lun,fmt='(A,A,A,i8)') "start writing ",trim(fname)," matrix size ",matrix_size

      allocate(buf_out(max_nsf), buf_in(max_nsf), buf_col_out(max_nsf), &
        pos_out(max_nsf), pos_col_out(max_nsf), stat=ierr)
      
      allocate(at_npaos(ni_in_cell), displ(numprocs), nz_receive(numprocs), &
        stat=ierr)
   
      iprim = 0
      maxijk=0
      at_npaos=0d0    
      
      do i = 1, ni_in_cell
        at_npaos(i)=nsf_species(species_glob(i))
      end do

      do irun = 0, 3
        call system_clock(counti, count_rate)
        write(io_lun,fmt='(i3," / 3 ")', advance="no") irun
        iprim = 0
        do np = 1,bundle%groups_on_node ! Loop over primary set partitions
          do i = 1,bundle%nm_nodgroup(np) ! Loop over atoms
            atom_num = bundle%nm_nodbeg(np) + i - 1
            iprim = iprim + 1
            atom_i = bundle%ig_prim(iprim)         ! global number of i
            atom_spec = bundle%species(atom_num)   ! atomic species of i
            nsf1 = nsf_species(atom_spec)          ! number of SFs on i
            
            ii = 0 ! column index 
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
!~                   write(0,fmt='(A,2i8,e24.12)') "atoms ",atom_i,atom_j,dz
              i1=nint(dx/r_super_x)
              i2=nint(dy/r_super_y)
              i3=nint(dz/r_super_z)
              
              if (irun.eq.0) then
                maxijk(1)=max(maxijk(1),abs(i1))
                maxijk(2)=max(maxijk(2),abs(i2))
                maxijk(3)=max(maxijk(3),abs(i3))
                cycle
              end if

              j_in_halo = halo(Arange)%i_halo(gcspart)                                         ! halo-ID of j
              neigh_species = species_glob(atom_j)                                             ! atomic species of j
              nsf2 = nsf_species(neigh_species)                                                ! number of SFs on j

  
              do k = 1, nsf1
                irow = k+sum(at_npaos(1:atom_i-1))
                do l = 1, nsf2
                  icol = l+sum(at_npaos(1:atom_j-1))                      
                  
                  if (irun.eq.1) then
                  
                    nzrow_local(irow,i1,i2,i3) = nzrow_local(irow,i1,i2,i3) + 1
                    
                  else if (irun .eq. 2) then
                  
                    nzrow_local(irow,i1,i2,i3) = nzrow_local(irow,i1,i2,i3) + 1
                    mat_info(irow,i1,i2,i3)%cols(nzrow_local(irow,i1,i2,i3)) = icol
                  
                  else if (irun .eq. 3) then
                  
                    nzrow_local(irow,i1,i2,i3) = nzrow_local(irow,i1,i2,i3) + 1
                    
                    wheremat = matrix_pos(matA,iprim,j_in_halo,k,l)
                    rval = mat_p(matA)%matrix(wheremat)
                    cval = cmplx(rval,0d0,double_cplx)
                                
                    
                    row_offset = sum(nzrow(1:irow-1,i1,i2,i3))
                    ipos = header_offset + row_offset * isize + mat_info(irow,i1,i2,i3)%offsets(nzrow_local(irow,i1,i2,i3)) * isize
                    pos_col_out(l) = ipos
                    buf_col_out(l) = icol -1
!~                     write(unit=iunit(i1,i2,i3),pos=ipos) icol - 1
                    
                    ipos = header_offset + nzcount(i1,i2,i3) * isize + row_offset * zsize + mat_info(irow,i1,i2,i3)%offsets(nzrow_local(irow,i1,i2,i3)) * zsize
                    pos_out(l) = ipos
                    buf_out(l) = cval
!~                     write(unit=iunit(i1,i2,i3),pos=ipos) cval
                    
!~ debug needs offset for all i1,i2,i3 or we also create seperate files                    
!~                     write(si1,fmt='(5i7,es45.24e5,2i8,X,A,X,i8)')  i1,i2,i3,irow,icol,rval,atom_i,atom_j,trim(fname),inode
!~                     si1=adjustl(si1)
!~                     ipos = ( sum(nzrow(1:irow-1,i1,i2,i3)) + mat_info(irow,i1,i2,i3)%offsets(nzrow_local(irow,i1,i2,i3))  ) * (len(trim(si1))+1) +1
!~                     write(lun,fmt='(A)', pos=ipos)  trim(si1)
!~                     write(0,fmt='(5i7,es45.24e5,2i8,X,A,X,i8)')  i1,i2,i3,irow,icol,rval,ipos,len(trim(si1)),trim(fname),inode
                    
                  end if
                end do
                if (irun .eq. 3) then
                  do l = 1, nsf2-1
                    if ((pos_col_out(l+1)-pos_col_out(l)).ne.4) then
                      write(0,*) "buf_col allignment error ",pos_col_out(1:nsf2)
                      stop
                    end if
                    if ((pos_out(l+1)-pos_out(l)).ne.16) then
                      write(0,*) "buf allignment error ",pos_out(1:nsf2)
                      stop
                    end if                    
                  end do
                  if (lreadandadd) then
                    i_offset = pos_out(1) - 1
                    call MPI_FILE_READ_AT(iunit(i1,i2,i3), i_offset, buf_col_out(1:nsf2), &
                     nsf2, MPI_INTEGER, MPI_STATUS_IGNORE, ierr)                  
!~                     read(unit = iunit(i1,i2,i3), pos = pos_out(1)) buf_in(1:nsf2)
                    buf_in(1:nsf2) = swap_endian_z(buf_in(1:nsf2))
                    buf_out(1:nsf2) = buf_in(1:nsf2) + buf_out(1:nsf2)
                  end if
                  i_offset = pos_col_out(1) - 1
                  buf_col_out(1:nsf2) = swap_endian_i4(buf_col_out(1:nsf2))
                  call MPI_FILE_WRITE_AT(iunit(i1,i2,i3), i_offset, buf_col_out(1:nsf2), nsf2, &
                  MPI_INTEGER, MPI_STATUS_IGNORE, ierr)
                  i_offset = pos_out(1) - 1
                  buf_out(1:nsf2) = swap_endian_z(buf_out(1:nsf2))
                  call MPI_FILE_WRITE_AT(iunit(i1,i2,i3), i_offset, buf_out(1:nsf2), nsf2, &
                  MPI_DOUBLE_COMPLEX, MPI_STATUS_IGNORE, ierr)
!~                   write(unit=iunit(i1,i2,i3),pos=pos_col_out(1)) swap_endian_i4(buf_col_out(1:nsf2))
!~                   write(unit=iunit(i1,i2,i3),pos=pos_out(1)) swap_endian_z(buf_out(1:nsf2))
                end if
              end do

            end do ! j
          end do ! i          
        end do ! np
        
        if (irun .eq. 0) then
        
          call gmax(maxijk(1))
          call gmax(maxijk(2))
          call gmax(maxijk(3))
          
          allocate(nzrow(matrix_size,-maxijk(1):maxijk(1),-maxijk(2):maxijk(2),-maxijk(3):maxijk(3)),&
            nzrow_local(matrix_size,-maxijk(1):maxijk(1),-maxijk(2):maxijk(2),-maxijk(3):maxijk(3)),&
            nzcount(-maxijk(1):maxijk(1),-maxijk(2):maxijk(2),-maxijk(3):maxijk(3)),&
            mat_info(matrix_size,-maxijk(1):maxijk(1),-maxijk(2):maxijk(2),-maxijk(3):maxijk(3)),&
            petscfile(-maxijk(1):maxijk(1),-maxijk(2):maxijk(2),-maxijk(3):maxijk(3)),&
            iunit(-maxijk(1):maxijk(1),-maxijk(2):maxijk(2),-maxijk(3):maxijk(3)), stat=ierr)
          nzrow=0
          nzrow_local=0
        else if (irun .eq. 1) then
        
          nzrow = nzrow_local
          do i3=-maxijk(3),maxijk(3)
            do i2=-maxijk(2),maxijk(2)
              do i1=-maxijk(1),maxijk(1)
                do irow=1,matrix_size
                  call gsum(nzrow(irow,i1,i2,i3))
                  nz_l=nzrow_local(irow,i1,i2,i3)
                  allocate(mat_info(irow,i1,i2,i3)%cols(nz_l),&
                    mat_info(irow,i1,i2,i3)%offsets(nz_l), stat=ierr)  
                  if (ierr.ne.0) then
                    write(0,fmt='(i8,A,i8)') inode," allocation error mat_info ",ierr
                  end if
                end do
                nzcount(i1,i2,i3) = sum(nzrow(:,i1,i2,i3))                         
              end do
            end do
          end do
          
          nzrow_local = 0

          
        else if (irun .eq. 2) then
          nz_g = maxval(nzrow)
          allocate(cols_global(nz_g), cols_offset(nz_g), cols_global_idx(nz_g), &
            stat=ierr)        
          do i3=-maxijk(3),maxijk(3)
            do i2=-maxijk(2),maxijk(2)
              do i1=-maxijk(1),maxijk(1)
                do irow=1,matrix_size
                  nz_l=nzrow_local(irow,i1,i2,i3)
                  nz_g=nzrow(irow,i1,i2,i3)   
                  if (nz_g.eq.0) cycle
                  nz_receive=-1
                  call MPI_Gather(nz_l, 1, MPI_INTEGER, nz_receive, 1, MPI_INTEGER, ionode-1, MPI_COMM_WORLD,ierr)      
                  if (ierr.ne.0) then
                    write(0,fmt='(i8,A,i8)') inode," MPI_Gather nz_l failed ",ierr
                    stop
                  end if                  
                  displ=cshift(nz_receive,-1)
                  displ(1)=0
                  
                  cols_global=-1
                  if ((nz_g.ne.sum(nz_receive)).and.(inode.eq.ionode)) then
                    write(io_lun,fmt='(i4,A,6i12,A)') "nz_g .ne. sum nz_receive ",nz_g,sum(nz_receive),i1,i2,i3,irow                    
                    stop
                  end if
                  call MPI_Gatherv(mat_info(irow,i1,i2,i3)%cols(1:nz_l), nz_l, MPI_INTEGER, &
                    cols_global(1:nz_g), nz_receive, displ, MPI_INTEGER, ionode-1, MPI_COMM_WORLD, ierr)
                  if (ierr.ne.0) then
                    write(0,fmt='(i8,A,i8)') inode," MPI_Gatherv cols_local -> cols_global failed ",ierr
                    stop
                  end if                  
!~                   write(io_lun,fmt='(A,6i6," XX ",999i6)') "cols ",i1,i2,i3,irow,matrix_size,nz_g,cols_global(1:nz_g)  
                                  
                  if (ionode .eq. inode) then
                    cols_offset(1:nz_g) = 0
         
                    do icol = 1, nz_g
                      cols_global_idx(icol) = icol
                    end do
                    call sort_ab(cols_global(1:nz_g), cols_global_idx(1:nz_g))    
                    
                    do i = 1, nz_g
                      cols_offset(cols_global_idx(i)) = i - 1
                    end do                
                                        
                    
!~                     do icol = 1, nz_g  
!~                       do k = 1, nz_g
!~                         if (cols_global(k).lt.cols_global(icol)) cols_offset(icol)=cols_offset(icol) + 1
!~                       end do
!~                     end do

                  end if
                  
                  call MPI_SCATTERV(cols_offset(1:nz_g), nz_receive, displ, MPI_INTEGER,&
                   mat_info(irow,i1,i2,i3)%offsets(1:nz_l), nz_l, MPI_INTEGER, ionode-1, MPI_COMM_WORLD, ierr)
                  if (ierr.ne.0) then
                    write(0,fmt='(i8,A,i8)') inode," MPI_Scatterv cols_offset -> local offsets ",ierr
                    stop
                  end if                     
!~                   write(io_lun,fmt='(A,6i6," XX ",999i6)') "offs ",i1,i2,i3,irow,matrix_size,nz_g,cols_offset(1:nz_g)
                  
                end do ! irow
                
! prepare filenames                                
                
                write(si3,*) i3
                si3=adjustl(si3)
                write(si2,*) i2
                si2=adjustl(si2)
                write(si1,*) i1
                si1=adjustl(si1)
                petscfile(i1,i2,i3)=trim(fname)//"_"//trim(si1)//"_"//trim(si2)//&
                  "_"//trim(si3)//"_petsc.dat"
            
! using fortran 2003 standard: access="stream", read/write controlled by pos=x where pos is here given in byte
! -> integer(4) -> 4 byte, double precision -> 8 byte, double complex -> 16 byte
! PETSc MAT_LOAD need ordered nonzero columns, that is very cumbersome to work out here. -> Use MPI_READs in transport code
! to load matrices and work out process partition. No BIG_ENDIAN necessary as for MAT_LOAD
                if (ionode.eq.ionode) then
                  if (.not.lreadandadd) then
                  
                    open(newunit=iunit(i1,i2,i3),file=trim(petscfile(i1,i2,i3)),access="stream",&
                      status="replace", BUFFERED='YES' )
!~                       status="replace",CONVERT='BIG_ENDIAN', BUFFERED='YES' )
! header ---- start    
                      ipos = 1 ! 1
                      write(iunit(i1,i2,i3),pos=ipos) swap_endian_i4(1211216)
                      ipos = ipos + isize ! 5
                      write(iunit(i1,i2,i3),pos=ipos) swap_endian_i4(matrix_size)
                      ipos = ipos + isize ! 9
                      write(iunit(i1,i2,i3),pos=ipos) swap_endian_i4(matrix_size) 
                      ipos = ipos + isize ! 13
                      write(iunit(i1,i2,i3),pos=ipos) swap_endian_i4(nzcount(i1,i2,i3))
                      ipos = ipos + isize ! 17
                      write(iunit(i1,i2,i3),pos=ipos) swap_endian_i4(nzrow(1:matrix_size,i1,i2,i3))
                      ipos = ipos + isize * matrix_size   
! header ---- start                        
                      close(iunit(i1,i2,i3))
                      if (ldebug) then
                        open(newunit=lun, file=trim(fname)//".dat", status="replace", access="stream", &
                          form="formatted")
                        close(lun)
                      end if
                  end if
                  
                end if ! ionode
                
                header_offset = 1 + 4 * isize + matrix_size * isize
                
                call MPI_BARRIER(MPI_COMM_WORLD,ierr)
                                
                
                call MPI_FILE_OPEN(MPI_COMM_WORLD, trim(petscfile(i1,i2,i3)), &
                  MPI_MODE_WRONLY + MPI_MODE_CREATE, MPI_INFO_NULL, iunit(i1,i2,i3), ierr)
                                    
!~                 open(newunit=iunit(i1,i2,i3),file=trim(petscfile(i1,i2,i3)),access="stream",&
!~                   status="old", BUFFERED='YES')
!~                   status="old",CONVERT='BIG_ENDIAN', BUFFERED='YES')
            
              end do ! i1
            end do ! i2
          end do !i3   
          
          if (ldebug) open(newunit=lun, file=trim(fname)//".dat", status="old", &
            access="stream", form="formatted")
          nzrow_local = 0
          
          deallocate(cols_global, cols_offset, cols_global_idx, stat=ierr)        
          
        end if  ! irun = 2
        call system_clock(countf)
        write (io_lun, fmt='(X,e24.12)') real(countf - counti, 8)/real(count_rate, 8)              
      end do ! irun

  
      if (ldebug) close(lun)
      do i3=-maxijk(3),maxijk(3)
        do i2=-maxijk(2),maxijk(2)
          do i1=-maxijk(1),maxijk(1)
            call MPI_File_Close(iunit(i1,i2,i3), ierr)
!~             close(iunit(i1,i2,i3))
          end do
        end do
      end do


    end subroutine write_mat
  
    subroutine sort_ab(A, B)
      implicit none

        integer :: A(:), B(:)
        integer :: i, j, n, l ,ir, a_save, b_save
        
        n = size(A)        
        
        l = n/2+1
        ir = n

        do 
          if(l .gt. 1) then
            L = l - 1
            a_save = A(l)
            b_save = B(l)
          else
            a_save = A(ir)
            b_save = B(ir)
            A(ir) = A(1)
            B(ir) = B(1)
            ir = ir - 1
            if(ir .eq. 1) then
              A(1)=a_save
              B(1)=b_save
              return
            end if
          end if
          i = l
          j = l + l
          do 
            if(j .le. ir) then
              if(j .lt. ir) then
                if(A(j) .le. A(j+1)) j = j + 1
              end if
              if(a_save .lt. A(j)) then
                A(i) = A(j)
                B(i) = B(j)
                i = j
                j = j + j
              else
                j = ir + 1
              end if
            else
              exit
            end if        
          end do           
          A(i) = a_save
          B(i) = b_save
        end do
    
    end subroutine sort_ab 
    
    elemental function swap_endian_z(x) result(res)
      complex(8) :: res
      double complex,intent(in)::x
      integer(1),dimension(8):: bytes_r, bytes_i
      real(8) :: r_r, r_i

      bytes_r = transfer(real(x,8),bytes_r,8)
      bytes_i = transfer(aimag(x),bytes_i,8)
      bytes_r(1:8) = bytes_r(8:1:-1)
      bytes_i(1:8) = bytes_i(8:1:-1)
      r_r = transfer(bytes_r,r_r)
      r_i = transfer(bytes_i,r_i)
      res = cmplx(r_r, r_i, 8)
      
    end function swap_endian_z
    
    elemental function swap_endian_i4(x) result(res)
      integer(4) :: res
      integer(4), intent(in)::x
      integer(1),dimension(4):: bytes

      bytes = transfer(x,bytes,4)
      bytes(1:4) = bytes(4:1:-1)
      res = transfer(bytes,res)
      
    end function swap_endian_i4    
  
end module negf_output_module
