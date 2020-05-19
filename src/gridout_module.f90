module gridout_module
  implicit none

  contains

! quicksort maybe faster but heapsort much easier to implement

  subroutine  heapsort (ra,rb)


    use datatypes

    implicit none

! input,output:
    integer :: ra(:)
    real(double) :: rb(:)
! local :
    real(double) :: rrb
    integer :: l,ir,n,rra,j,i,m

    n=size(ra,1)

    l=n/2+1
    ir=n

  do
    if(l.gt.1)then
      l=l-1
      rra=ra(l)
      rrb=rb(l)
    else
      rra=ra(ir)
      rrb=rb(ir)
      ra(ir)=ra(1)
      rb(ir)=rb(1)
      ir=ir-1
      if(ir.eq.1)then
        ra(1)=rra
        rb(1)=rrb
        return
      end if
    end if
    i=l
    j=l+l
    do while (j.le.ir)
      if(j.lt.ir)then
        if(ra(j).lt.ra(j+1)) j=j+1
      end if
      if(rra.lt.ra(j))then
        ra(i)=ra(j)
        rb(i)=rb(j)
        i=j
        j=j+j
      else
        j=ir+1
      end if

    end do
      ra(i)=rra
      rb(i)=rrb
    end do
  end subroutine heapsort



  subroutine print_grid_data(iwhat)

    use datatypes
    use global_module, only: nspin, spin_factor,io_lun, flag_neutral_atom,gridsolver_use
    use hartree_module, only: hartree
    use density_module, only: density,density_atom
    use maxima_module, only: maxngrid
    use block_module, only: n_blocks, n_pts_in_block, nx_in_block,ny_in_block,&
                            nz_in_block,n_pts_in_block
    use dimens, only: grid_point_volume, n_my_grid_points, n_grid_z, &
                      n_grid_x,n_grid_y
    use gencomms, only : inode,my_barrier,gcopy,ionode
    use grid_solver, only: get_hartree_gridsolver
!~     use io_module,       only: dump_charge

    implicit none

    !input
    integer :: iwhat,ii,ierr

    !local
    real(double),allocatable, target :: xgrid1(:),xgrid2(:)
    real(double) :: rdummy,mmm
    integer :: ix1,ix2,iy1,iy2,iz1,iz2,iix,iiy,iiz,idisp1,idisp2,igrid


    select case(iwhat)
      case(1) ! spin(difference) density and total density
        allocate(xgrid1(maxngrid),stat=ierr)
        if (ierr.ne.0) then
          write(0,*) "allocation error ",ierr
          stop
        end if
        xgrid1=0d0
        do ii=1,nspin
          xgrid1(:)=xgrid1(:)+spin_factor*density(:,ii)
        end do

!~         if (nspin.eq.1) then
!~           call dump_charge(xgrid1, n_my_grid_points, inode, spin=0)
!~         else
!~           call dump_charge(density(:,1), n_my_grid_points, inode, spin=1)
!~           call dump_charge(density(:,2), n_my_grid_points, inode, spin=2)
!~         end if

        if(flag_neutral_atom) then
          call grid_out_plt("td.plt","total density",250,xgrid1)
          xgrid1=xgrid1-density_atom
          call grid_out_plt("td_diff.plt","total density",250,xgrid1)
        else
          call grid_out_plt("td.plt","total density",250,xgrid1)
        end if

        if (nspin.eq.2) then
          xgrid1=density(:,1)-density(:,2)
          call grid_out_plt("sd.plt","spin(difference) density",250,xgrid1)
        end if
        deallocate(xgrid1)
      case(2) ! hartree potential

        allocate(xgrid1(maxngrid),xgrid2(maxngrid),stat=ierr)
        xgrid1=0d0
        do ii=1,nspin
          xgrid1(:)=xgrid1(:)+spin_factor*density(:,ii)
        end do
        xgrid2=0d0

        if(flag_neutral_atom) then

          xgrid1=xgrid1-density_atom

          xgrid2=0d0
          if (gridsolver_use) then
            call get_hartree_gridsolver(xgrid1,xgrid2)
          else
            call hartree(xgrid1,xgrid2,maxngrid,rdummy)
          end if
          call grid_out_plt("hartree_pot_diff.plt","hartree potential",250,xgrid2)
        else
          if (gridsolver_use) then
            call get_hartree_gridsolver(xgrid1,xgrid2)
          else
            call grid_out_plt("dens.plt","dens",250,xgrid1)
            call hartree(xgrid1,xgrid2,maxngrid,rdummy)
          end if
          call grid_out_plt("hartree_pot.plt","hartree potential",250,xgrid2)
        end if

        deallocate(xgrid1)
      case(3)
        allocate(xgrid1(maxngrid),xgrid2(maxngrid),stat=ierr)
        xgrid1=0d0
        do ii=1,nspin
          xgrid1(:)=xgrid1(:)+spin_factor*density(:,ii)
        end do

        if(flag_neutral_atom) xgrid1=xgrid1-density_atom

        if (gridsolver_use) then
          xgrid2=0d0
          call get_hartree_gridsolver(xgrid1,xgrid2)
        else
          call hartree(xgrid1,xgrid2,maxngrid,rdummy)
        end if


        call grid_out_plt("hartree_pot_init.plt","hartree potential",250,xgrid2)
        deallocate(xgrid1)

      case(4) ! spin(difference) density and total density
        allocate(xgrid1(maxngrid),stat=ierr)
        if (ierr.ne.0) then
          write(0,*) "allocation error ",ierr
          stop
        end if
        xgrid1=0d0
        do ii=1,nspin
          xgrid1(:)=xgrid1(:)+spin_factor*density(:,ii)
        end do
        call grid_out_plt("td_init.plt","total density",250,xgrid1)
        if(flag_neutral_atom) xgrid1=xgrid1-density_atom
        call grid_out_plt("td_diff_init.plt","total density",250,xgrid1)
        if (nspin.eq.2) then
          xgrid1=density(:,1)-density(:,2)
          call grid_out_plt("sd_init.plt","spin(difference) density",250,xgrid1)
        end if
        deallocate(xgrid1)
      case default
        if (inode.eq.ionode) write(io_lun,*) "no plt output for ",iwhat
    end select


  end subroutine



  subroutine grid_out_plt(fileout,whatout,itype,xgrid)

    use global_module, only: io_lun,numprocs
    use GenComms, ONLY: inode,ionode,gmax,gmin,gcopy
    use datatypes
    use dimens, only: n_grid_x,n_grid_y,n_grid_z
    use global_module, only: rcellx,rcelly,rcellz
    use block_module, only: nx_in_block,ny_in_block,nz_in_block,n_pts_in_block
    use primary_module, only: domain
    use group_module, only: blocks,parts
    use units, only : BohrToAng
    use mpi
    use grid_index, ONLY: ind_block_x, ind_block_y, ind_block_z, &
         grid_point_x, grid_point_y, grid_point_z, &
         grid_point_block, grid_point_position


    implicit none

!input
    character(*) :: fileout,whatout
    real(double), intent(in) :: xgrid(:)
    integer :: itype
!local
    character(256) :: fileout_z,errorstr
!~     integer(kind=MPI_OFFSET_KIND) :: disp
    integer(8) :: disp
    integer :: filetype
    integer :: ierr,ifh,ifh_z,iblock,iix,iiy,iiz,ix,iy,iz,igrid,ipoint,ipp,io1,io2,ierr2,il
    real(double) :: dx,dy,dz,zmax,ymax,xmax,maxz,minz
!~     real(4) :: buf(n_pts_in_block*domain%groups_on_node )
!~     real(4) :: buf_z(1)
    real(8) :: buf(n_pts_in_block*domain%groups_on_node )
    real(8) :: buf_z(1)
!~     integer(kind=MPI_OFFSET_KIND) :: offsets(n_pts_in_block*domain%groups_on_node )
!~     integer(kind=MPI_OFFSET_KIND) :: offsets_z(1)
    integer :: offsets(n_pts_in_block*domain%groups_on_node )
    integer :: offsets_z(1),iunit
    character(256) :: istr







    if (inode.eq.ionode) then
      call MPI_FILE_DELETE(trim(fileout), MPI_INFO_NULL, ierr)
    end if

    call MPI_BARRIER(MPI_COMM_WORLD,ierr)

    call MPI_FILE_OPEN(MPI_COMM_WORLD, trim(fileout), &
                       MPI_MODE_WRONLY + MPI_MODE_CREATE, &
                       MPI_INFO_NULL, ifh, ierr)

    if (ierr.ne.0) then
      call MPI_Error_string(ierr,errorstr,iblock,iiz)
      write(0,*) trim(errorstr)
      write(0,*) "MPI_FILE_OPEN error ",ierr
      call flush(0)
      stop
    end if

    dx = rcellx / (blocks%ngcellx * nx_in_block)
    dy = rcelly / (blocks%ngcelly * ny_in_block)
    dz = rcellz / (blocks%ngcellz * nz_in_block)

    if (inode.eq.ionode) then
      write(io_lun,*) "output   ",trim(whatout)," to ",trim(fileout),ifh
      call flush(io_lun)
      disp=0  ! file starts at 0 not at 1, C convention !!!!
      call mpi_file_seek(ifh,disp,MPI_SEEK_SET,ierr)
      call mpi_file_write(ifh, int((/ 3,itype,n_grid_z,n_grid_y,n_grid_x /),8), 5, &
                          MPI_INTEGER8, MPI_STATUS_IGNORE, ierr)

      xmax=(rcellx)!*BohrToAng !-dx)
      ymax=(rcelly)!*BohrToAng !-dy)
      zmax=(rcellz)!*BohrToAng !-dz)
      call mpi_file_write(ifh, real((/ 0d0,zmax,0d0,ymax,0d0,xmax /),8), 6, &
                          MPI_DOUBLE, MPI_STATUS_IGNORE, ierr)
    end if

    maxz=0d0
    minz=0d0

    ipp=0
    io1=-1
    io2=-1

    write(istr,fmt='(i8)') inode

    do iblock = 1, domain%groups_on_node ! loop over blocks of grid points
      ipoint = 0
      ! loop over the grid points in the current block
      do iz = 1, nz_in_block
        iiz=(domain%idisp_primz(iblock) + domain%nz_origin - 1)*nz_in_block+iz
        do iy = 1, ny_in_block
          do ix = 1, nx_in_block
            ipp=ipp+1
            ipoint=ipoint+1
            iiy=(domain%idisp_primy(iblock) + domain%ny_origin - 1)*ny_in_block+iy
            igrid = n_pts_in_block * (iblock - 1) + ipoint
           !coordinates of igrid point in cartesian coordiantes R=origin+iix*dx+iiy*dy+iiz*dz
            iix=(domain%idisp_primx(iblock) + domain%nx_origin - 1)*nx_in_block+ix
            disp=(iix+(iiy-1)*n_grid_x+(iiz-1)*n_grid_y*n_grid_x+10)
            offsets(ipp)=disp
            buf(ipp)=real(xgrid(igrid),8)

          end do
        end do
      end do
    end do

    if (domain%groups_on_node.ge.1) call heapsort(offsets,buf)

    call MPI_TYPE_CREATE_INDEXED_BLOCK(domain%groups_on_node*n_pts_in_block, 1, offsets, &
      MPI_DOUBLE, filetype, ierr)
    call MPI_TYPE_COMMIT(filetype, ierr)

    disp = 0

    call MPI_FILE_SET_VIEW(ifh, disp, MPI_DOUBLE, filetype, 'native', MPI_INFO_NULL, ierr)

! For ifort v17 (tested up to update 4) this can end up in a deadlock, version v16 or v18 work.
! Workaround would be to using mpi_file_write which is however noncollective and much slower than
! collective mpi_file_write_all
    call mpi_file_write(ifh,buf,domain%groups_on_node*n_pts_in_block,MPI_DOUBLE,MPI_STATUS_IGNORE,ierr)

    call MPI_FILE_CLOSE(ifh, ierr)

    call gmax(maxz)
    call gmin(minz)
    call gmax(io1)
    call gmax(io2)

    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    if (inode.eq.ionode) then
      write(io_lun,*) "finished plt"
    end if
  end subroutine grid_out_plt


end module gridout_module
