module grid_solver
  use datatypes
  implicit none

!~ poisson grid solver switches
  integer, parameter :: gridsolver_dlmg=1,gridsolver_mud=2

  integer :: gs(3),ge(3)
  integer :: nglx,ngly,nglz,nldx,nldy,nldz,gridsolver_slice,inmax

  real(double), pointer :: pot3(:,:,:)
  integer, allocatable, target :: gs_prbx(:),gs_prxb(:),gs_packbx(:),gs_packxb(:)
  integer, allocatable :: gs_rprbs(:),gs_rprxb(:),slice(:,:)
  integer, pointer :: gs_unpackbx(:), gs_unpackxb(:),gs_psendxb(:), gs_psendbx(:)
  integer, allocatable :: gs_x_columns_node(:), gs_y_columns_node(:), gs_z_columns_node(:)

!~  dl_mg
  integer :: mg_comm,mg_unit

  contains

    subroutine rearrange( cdata, size, pack, psend, pr, unpack )

      use global_module, ONLY: numprocs, area_SC
      use maxima_module, ONLY: maxngrid
      use GenComms, ONLY: exchv, cq_abort,inode
      use memory_module, ONLY: reg_alloc_mem, reg_dealloc_mem, type_dbl, type_int

      ! Passed variables
      integer :: size
      real(double) :: cdata(size)
      real(double), allocatable :: send(:)
      real(double), allocatable :: recv(:)

      integer :: pack(:), psend(:), pr(:), unpack(:)

      ! Local variables
      integer :: I, stat


      allocate(send(maxngrid), recv(maxngrid), STAT=stat)
      if(stat/=0) call cq_abort("Allocation error in rearrange_data: ",size,stat)
      call reg_alloc_mem(area_SC,2*maxngrid,type_dbl)

      ! pack data into sending format
      do I = 1, psend(numprocs+1)-1
         send(pack(I)) = cdata(I)
      end do

      ! call global exchange routine
      call exchv(send, psend, recv, pr, maxngrid)

      ! unpack data
      if(pr(numprocs+1)-1>size) call cq_abort("Error in rearrange_data: ",size,pr(numprocs+1)-1)
      do I = 1, pr(numprocs+1)-1
         cdata(I) = recv(unpack(I))
      end do

      deallocate(send, recv, STAT=stat)
      if(stat/=0) call cq_abort("Deallocation error in rearrange_data: ",size,stat)
      call reg_dealloc_mem(area_SC,2*maxngrid,type_dbl)

      return
    end subroutine rearrange


!~ similar to subroutine set_map but slices the system along x, y, or, z direction and
!~ the slices are distributed to each process such that it can be used (after reordering)
!~ as input for dl_mg gridsolver or mpi fftw
    subroutine set_grid_map ( )
      use fft_module
      use global_module, ONLY: numprocs, area_SC, iprint_SC
      use numbers
      use dimens, ONLY: n_my_grid_points, n_grid_x, n_grid_y, n_grid_z, &
           r_super_x_squared, r_super_y_squared, r_super_z_squared, &
           r_super_x, r_super_y, r_super_z
      use grid_index, ONLY: grid_point_x, grid_point_y, grid_point_z
      use GenComms, ONLY: gcopy_diff, exch, cq_abort, my_barrier, inode, myid,ionode
      use maxima_module, ONLY: maxngrid
      use memory_module, ONLY: reg_alloc_mem, reg_dealloc_mem, type_dbl, type_int

      implicit none

      ! Local variables
      integer :: I, send_count(numprocs), count(numprocs), z, y, node_to, &
           lz, ly, x, lx, JNODE, xcol, get_count(numprocs), site, cnt, &
           col_num, N,ngrid_imax
      integer :: ierr
      integer :: sgs_x_columns_node(numprocs)

      integer :: remote_n_my_grid_points,nlx,nly,nlz,imax,ngrid(3),nlo,nup,nmax,nmin,nlgridx,nlgridy,nlgridz
      integer :: nslice(numprocs),xyz(3),ndown
      integer, allocatable :: remote_grid_point_x(:), remote_grid_point_y(:), remote_grid_point_z(:)


      ! Allocate variables

      allocate(gs_prbx(numprocs+1), gs_prxb(numprocs+1), STAT=ierr)
      if(ierr/=0) call cq_abort("Allocation error for pr in set_fft: ",numprocs+1,ierr)
      call reg_alloc_mem(area_SC,2*(numprocs+1),type_int)
      allocate(gs_rprbs(numprocs), gs_rprxb(numprocs), STAT=ierr)
      if(ierr/=0) call cq_abort("Allocation error for rpr in set_fft: ",numprocs,ierr)
      call reg_alloc_mem(area_SC,2*numprocs,type_int)
      allocate(gs_packbx(maxngrid),gs_packxb(maxngrid), STAT=ierr)
      if(ierr/=0) call cq_abort("Allocation error for pack in set_fft: ",maxngrid,ierr)
      call reg_alloc_mem(area_SC,2*maxngrid,type_int)
      allocate(gs_x_columns_node(numprocs), gs_y_columns_node(numprocs), gs_z_columns_node(numprocs), STAT=ierr)
      if(ierr/=0) call cq_abort("Allocation error for columns in set_fft: ",numprocs,ierr)
      call reg_alloc_mem(area_SC,3*numprocs,type_int)
      if (.not.allocated(remote_grid_point_x)) &
     & allocate(remote_grid_point_x(maxngrid),remote_grid_point_y(maxngrid),remote_grid_point_z(maxngrid),STAT=ierr)
      if(ierr/=0) call cq_abort("Allocation error for remote grids in set_fft: ",maxngrid,ierr)
      call reg_alloc_mem(area_SC,3*maxngrid,type_int)
      if (.not.allocated(slice)) allocate(slice(2,numprocs))
      ! Map the equivalent pointers onto their arrays
      gs_unpackxb => gs_packbx
      gs_unpackbx => gs_packxb
      gs_psendbx => gs_prxb
      gs_psendxb => gs_prbx


      ! first of all, count the number of data points to be SENT to each other
      ! node.
      call my_barrier
      do I = 1, numprocs
         send_count(I) = 0
         count(I) = 0
      end do
      imax=3

      if (n_grid_y.gt.n_grid_z) imax=2
      if (n_grid_x.gt.max(n_grid_y,n_grid_z)) imax=1
      gridsolver_slice=imax


      ngrid=(/ n_grid_x,n_grid_y,n_grid_z /)
      ndown=floor(real(ngrid(imax))/numprocs)
      slice=0
      nslice=0
      if (ndown.lt.2) then
        inmax=floor(real(ngrid(imax))/2d0)
        nlo=inmax-mod(ngrid(imax),2)
        nup=inmax+mod(ngrid(imax)+1,2)
        ndown=2
      else
        inmax=numprocs
        nlo=numprocs*ndown+numprocs-ngrid(imax)
        nup=nlo+1
      end if

      if (inode.le.nlo) then
        ngrid(imax)=ndown
      else if (inode.ge.nup) then
        ngrid(imax)=ndown+1
      end if
      nslice(1:nlo)=ndown
      if (nup.le.inmax) nslice(nup:inmax)=ndown+1
      slice(1,1)=1
      slice(2,1)=ndown
      do n=2,inmax
        slice(1,n)=slice(2,n-1)+1
        slice(2,n)=slice(1,n)+nslice(n)-1
      end do


      if (inode.ge.inmax+1) ngrid(imax)=0

      gs=0
      ge=-1
      if (inode.le.inmax) then
        gs=(/ 1,1,1 /)
        ge=(/ n_grid_x,n_grid_y,n_grid_z /)
        gs(imax)=slice(1,inode)
        ge(imax)=slice(2,inode)
      end if

!--debug
!~        write(io_lun,fmt='(2i8,A,15i8)') inode,numprocs," loc grid ",imax,inmax,nlo,nup,ndown,&
!~        &ngrid(imax),ngrid,ngrid(1)*ngrid(2)*ngrid(3),maxngrid,nslice(inode),slice(1,inode),slice(2,inode)
!~        call flush(io_lun)
!~        call my_barrier()
!-------

      do n = 1, n_my_grid_points
        xyz=(/ grid_point_x(n),grid_point_y(n),grid_point_z(n) /)
         do i=1,numprocs
            if ((xyz(imax).ge.slice(1,i)).and.(xyz(imax).le.slice(2,i))) node_to=i
         end do
         send_count(node_to) = send_count(node_to) + 1
      end do
      gs_psendbx(1) = 1
      do I = 1, numprocs
         gs_psendbx(I+1) = gs_psendbx(I) + send_count(I)
         if(iprint_SC>3.AND.myid==0) write(io_lun,*) inode," proc sending bx to ",I,send_count(I)
      end do
      if(iprint_SC>3.AND.myid==0) write(io_lun,*) inode," bx total: ",gs_psendbx(numprocs)+send_count(numprocs)


      ! Make the gs_packbx
      do n = 1, n_my_grid_points
        xyz=(/ grid_point_x(n),grid_point_y(n),grid_point_z(n) /)
         do i=1,numprocs
            if ((xyz(imax).ge.slice(1,i)).and.(xyz(imax).le.slice(2,i))) node_to=i
         end do
         gs_packbx(n) = count(node_to) + gs_psendbx(node_to)
         count(node_to) = count(node_to) + 1
      end do

      ! We need to work out where on the remote node this
      ! data is going to need to be sent.
      call exch(send_count, get_count, 1)

      ! gs_prbx(I) is where we will unpack data from node I (for PVM)
      gs_prbx(1) = 1
      do I = 1, numprocs
         gs_prbx(I+1) = gs_prbx(I) + get_count(I)
      end do
      gs_x_columns_node(inode) = (gs_prbx(numprocs)+get_count(numprocs))/n_grid_x
      if(iprint_SC>3.AND.myid==0) write(io_lun,*) inode," xcols: ",gs_prbx(numprocs)+get_count(numprocs),gs_x_columns_node(inode)
      if (gs_prbx(numprocs)+get_count(numprocs)>maxngrid+1) &
           call cq_abort("Grid overflow (bx) in Gridsolver: ",gs_prbx(numprocs)+get_count(numprocs), maxngrid)

      ! return this value to be gs_rprbs (for SMAL)
      ! note that we can all return the inverse, gs_prxb, because it is the same
      ! as gs_psendbx
      call exch(gs_prbx, gs_rprbs, 1)
      call exch(gs_prxb, gs_rprxb, 1)
      ! This can probably be done better (?allgather?)
      do i=1,numprocs
         sgs_x_columns_node(i) = gs_x_columns_node(inode)
      enddo
      call exch(sgs_x_columns_node, gs_x_columns_node,1)
      if(iprint_SC>3.AND.myid==0) write(io_lun,*) inode," x columns: ",gs_x_columns_node

      ! now we need to figure out the unpacking.
      ! for each node, JNODE, we loop over all points and identify the ones which
      ! are sending to us, and then identify 'site', where in our x-columns data
      ! that point will be stored.
      cnt = 0
      do jnode = 1, numprocs
         call gcopy_diff(remote_n_my_grid_points,n_my_grid_points,jnode)
         Call gcopy_diff(remote_grid_point_x, grid_point_x, jnode, maxngrid )
         Call gcopy_diff(remote_grid_point_y, grid_point_y, jnode, maxngrid )
         Call gcopy_diff(remote_grid_point_z, grid_point_z, jnode, maxngrid )
         do n = 1, remote_n_my_grid_points

            xyz=(/ remote_grid_point_x(n),remote_grid_point_y(n),remote_grid_point_z(n) /)
            if ((xyz(imax).ge.slice(1,inode)).and.(xyz(imax).le.slice(2,inode))) then
               xyz(imax) = xyz(imax)-slice(1,inode)+1
               cnt = cnt + 1
               site = xyz(1)+(xyz(2)-1)*ngrid(1)+(xyz(3)-1)*ngrid(1)*ngrid(2)
               gs_unpackbx(site) = cnt
            end if
         end do
      end do


      ! now, we need to work out where on the remote node this
      ! data is going to need to be sent. Assume numprocs = 2**NDIMEN
      call exch(send_count, get_count, 1)

      deallocate(remote_grid_point_x,remote_grid_point_y,remote_grid_point_z,STAT=ierr)
      if(ierr/=0) call cq_abort("Deallocation error for remote grids in set_fft: ",maxngrid,ierr)
      call reg_dealloc_mem(area_SC,4*maxngrid,type_int)

      return
    end subroutine set_grid_map

!~ grid points are reordered along 1. x, 2. y, 3. z
    subroutine reorder_to_xyz(a)

      use GenComms, ONLY: inode
      use grid_index, ONLY: ind_block_x, ind_block_y, ind_block_z, &
          grid_point_x, grid_point_y, grid_point_z, &
          grid_point_block, grid_point_position

      use datatypes

      implicit none

      real(double), intent(inout), target, contiguous :: a(:)
      real(double), allocatable :: b(:)


      integer :: ix,iy,iz,ipos,ipts,ierr

      allocate(b(nglx*ngly*nglz),stat=ierr)
      if (ierr.ne.0) then
        write(0,*) "allocation error ",ierr
        stop
      end if


      do ipts=1,nglx*ngly*nglz
        ix=grid_point_x(ipts)-grid_point_x(1)+1
        iy=grid_point_y(ipts)-grid_point_y(1)+1
        iz=grid_point_z(ipts)-grid_point_z(1)+1
        ipos=ix+(iy-1)*nglx+(iz-1)*ngly*nglx
        b(ipos)=a(ipts)
      end do
      a=b
      deallocate(b)


    end subroutine reorder_to_xyz

    subroutine get_dirichlet(a,bc_file_phi_l,bc_file_phi_r)

      use GenComms, ONLY: inode,io_lun
      use grid_index, ONLY: ind_block_x, ind_block_y, ind_block_z, &
          grid_point_x, grid_point_y, grid_point_z, &
          grid_point_block, grid_point_position
      use datatypes
      use global_module, only: gridsolver_bc,negf_mul,negf_mur
      use dimens,           only: n_grid_x,n_grid_y,n_grid_z

      implicit none

      real(double), intent(inout) :: a(:,:,:)
      character(*),intent(in) :: bc_file_phi_l,bc_file_phi_r

      real(double) :: rin_phi_l,rin_phi_r,mul,mur
      integer :: ipos,ipts,ierr,iunit1,jx,jy,jz,jpos,idirichlet,iunit2,iunit3,iunit4
      integer(8) :: mxyz(3),idum
      integer :: mx,my,mz,nxyz(3),jxyz(3),ixyz(3),kxyz(3)


      open(newunit=iunit1,file=trim(bc_file_phi_l),status="old",access="stream",share='DENYNONE',ACTION='READ')
      if (trim(bc_file_phi_l).eq.trim(bc_file_phi_r)) then
        iunit3=iunit1
      else
        open(newunit=iunit3,file=trim(bc_file_phi_r),status="old",access="stream",share='DENYNONE',ACTION='READ')
      end if


      read(iunit1,pos=(1-1)*8+1) idum
      read(iunit1,pos=(2-1)*8+1) idum
      read(iunit1,pos=(3-1)*8+1) mxyz(3)
      read(iunit1,pos=(4-1)*8+1) mxyz(2)
      read(iunit1,pos=(5-1)*8+1) mxyz(1)

      idirichlet=0
      if (gridsolver_bc(1).eq.2) then
        idirichlet=1
      else if (gridsolver_bc(2).eq.2) then
        idirichlet=2
      else if (gridsolver_bc(3).eq.2) then
        idirichlet=3
      end if

      nxyz=(/ n_grid_x,n_grid_y,n_grid_z /)

      do jz=gs(3),ge(3)+nldz
        do jy=gs(2),ge(2)+nldy
          do jx=gs(1),ge(1)+nldx

            ixyz(1)=jx
            ixyz(2)=jy
            ixyz(3)=jz
            jxyz=ixyz
            jpos=jxyz(1)+(jxyz(2)-1)*mxyz(1)+(jxyz(3)-1)*mxyz(2)*mxyz(1)+11

            if (.not.((ixyz(idirichlet).eq.1).or.(ixyz(idirichlet).eq.nxyz(idirichlet)))) cycle

            if (.not.( &
              & ( (ixyz(idirichlet).ge.1).and.(ixyz(idirichlet).le.mxyz(idirichlet)) ).or.&
              & ( (ixyz(idirichlet).ge.nxyz(idirichlet)-mxyz(idirichlet)+1).and.(ixyz(idirichlet).le.nxyz(idirichlet)) ) &
              & )) cycle

            if ( (ixyz(idirichlet).ge.1).and.(ixyz(idirichlet).le.mxyz(idirichlet)) ) then
              jxyz(idirichlet)= ixyz(idirichlet)
            else if ( (ixyz(idirichlet).ge.nxyz(idirichlet)-mxyz(idirichlet)+1).and.(ixyz(idirichlet).le.nxyz(idirichlet)) ) then
              jxyz(idirichlet)=ixyz(idirichlet)-nxyz(idirichlet)+mxyz(idirichlet)
            end if

            jpos=jxyz(1)+(jxyz(2)-1)*mxyz(1)+(jxyz(3)-1)*mxyz(2)*mxyz(1)+11

            kxyz(1)=ixyz(1)-gs(1)+1
            kxyz(2)=ixyz(2)-gs(2)+1
            kxyz(3)=ixyz(3)-gs(3)+1

            if ( (ixyz(idirichlet).ge.1).and.(ixyz(idirichlet).le.mxyz(idirichlet)) ) then
              read(unit=iunit1,pos=(jpos-1)*8+1) rin_phi_l
              a(kxyz(1),kxyz(2),kxyz(3))=(rin_phi_l+negf_mul)
            end if

            if ( (ixyz(idirichlet).ge.nxyz(idirichlet)-mxyz(idirichlet)+1).and.(ixyz(idirichlet).le.nxyz(idirichlet)) ) then
              read(unit=iunit3,pos=(jpos-1)*8+1) rin_phi_r
              a(kxyz(1)+nldx,kxyz(2)+nldy,kxyz(3)+nldz)=(rin_phi_r+negf_mur)
            end if

          end do
        end do
      end do

      close(iunit1)
      close(iunit3)

    end subroutine get_dirichlet


!~ this contains redundant code with get_dirichlet and could be merged if we account for possible shifts (chemical potential,...)
!~ and scaling factors (spin,...)
    subroutine fix_rho(b,bc_file_rho_l,bc_file_rho_r)

      use GenComms, ONLY: inode,io_lun,ionode
      use grid_index, ONLY: ind_block_x, ind_block_y, ind_block_z, &
          grid_point_x, grid_point_y, grid_point_z, &
          grid_point_block, grid_point_position
      use datatypes
      use global_module, only: gridsolver_bc,spin_factor
      use dimens,           only: n_grid_x,n_grid_y,n_grid_z

      implicit none

      real(double), parameter :: pi=4d0*datan(1d0)

      character(*),intent(in) :: bc_file_rho_l,bc_file_rho_r
      real(double), intent(inout) :: b(:,:,:)

      real(double) :: rin_phi_l,rin_rho_l,rin_phi_r,rin_rho_r,mul,mur,z1,z2

      integer :: ipos,ipts,ierr,iunit1,jx,jy,jz,jpos,idirichlet,iunit2,iunit3,iunit4
      integer(8) :: mxyz(3),idum,mcut
      integer :: mx,my,mz,nxyz(3),jxyz(3),ixyz(3),kxyz(3)
      logical :: lexist


      inquire(file=trim(bc_file_rho_l),exist=lexist)
      if (inode.eq.ionode) write(io_lun,fmt='(A,i8,4A,l)') "bc_file",inode,trim(bc_file_rho_l)," ",trim(bc_file_rho_r)," ",lexist
      if (.not.lexist) return


      open(newunit=iunit2,file=trim(bc_file_rho_l),status="old",access="stream")
      if (trim(bc_file_rho_l).eq.trim(bc_file_rho_r)) then
        iunit4=iunit2
      else
        open(newunit=iunit4,file=trim(bc_file_rho_r),status="old",access="stream")
      end if


      read(iunit2,pos=8*(1-1)+1) idum
      read(iunit2,pos=8*(2-1)+1) idum
      read(iunit2,pos=8*(3-1)+1) mxyz(3)
      read(iunit2,pos=8*(4-1)+1) mxyz(2)
      read(iunit2,pos=8*(5-1)+1) mxyz(1)

      idirichlet=0
      if (gridsolver_bc(1).eq.2) then
        idirichlet=1
      else if (gridsolver_bc(2).eq.2) then
        idirichlet=2
      else if (gridsolver_bc(3).eq.2) then
        idirichlet=3
      end if

      nxyz=(/ n_grid_x,n_grid_y,n_grid_z /)

      mcut=mxyz(idirichlet)

      do jz=gs(3),ge(3)
        do jy=gs(2),ge(2)
          do jx=gs(1),ge(1)

            ixyz(1)=jx
            ixyz(2)=jy
            ixyz(3)=jz
            jxyz=ixyz
            jpos=jxyz(1)+(jxyz(2)-1)*mxyz(1)+(jxyz(3)-1)*mxyz(2)*mxyz(1)+11

            if (.not.( &
              & ( (ixyz(idirichlet).ge.1).and.(ixyz(idirichlet).le.mcut) ).or.&
              & ( (ixyz(idirichlet).ge.nxyz(idirichlet)-mcut+1).and.(ixyz(idirichlet).le.nxyz(idirichlet)) ) &
              & )) cycle

            if ( (ixyz(idirichlet).ge.1).and.(ixyz(idirichlet).le.mcut) ) then
              jxyz(idirichlet)= ixyz(idirichlet)
            else if ( (ixyz(idirichlet).ge.nxyz(idirichlet)-mcut+1).and.(ixyz(idirichlet).le.nxyz(idirichlet)) ) then
              jxyz(idirichlet)=ixyz(idirichlet)-nxyz(idirichlet)+mcut
            end if

            jpos=jxyz(1)+(jxyz(2)-1)*mxyz(1)+(jxyz(3)-1)*mxyz(2)*mxyz(1)+11

            kxyz(1)=ixyz(1)-gs(1)+1
            kxyz(2)=ixyz(2)-gs(2)+1
            kxyz(3)=ixyz(3)-gs(3)+1

            if ( (ixyz(idirichlet).ge.1).and.(ixyz(idirichlet).le.mcut) ) then
              read(unit=iunit2,pos=8*(jpos-1)+1) rin_rho_l
              b(kxyz(1),kxyz(2),kxyz(3))=rin_rho_l !/spin_factor
            end if

            if ( (ixyz(idirichlet).ge.nxyz(idirichlet)-mcut+1).and.(ixyz(idirichlet).le.nxyz(idirichlet)) ) then
              read(unit=iunit4,pos=8*(jpos-1)+1) rin_rho_r
              b(kxyz(1)+nldx,kxyz(2)+nldy,kxyz(3)+nldz)=rin_rho_r !/spin_factor
            end if

          end do
        end do
      end do

      close(iunit2)
      close(iunit4)

    end subroutine fix_rho

    subroutine read_rho(b,data_file)

      use GenComms, ONLY: inode,io_lun,my_barrier,ionode
      use grid_index, ONLY: ind_block_x, ind_block_y, ind_block_z, &
          grid_point_x, grid_point_y, grid_point_z, &
          grid_point_block, grid_point_position
      use datatypes
      use global_module, only: gridsolver_bc,spin_factor
      use dimens,           only: n_grid_x,n_grid_y,n_grid_z

      implicit none

      real(double), intent(inout) :: b(:,:,:)
      character(*), intent(in) :: data_file

      character(256) :: instr
      integer :: ipos,ipts,ierr,iunit1,jx,jy,jz,jpos,idirichlet,iunit2,iunit3,iunit4
      integer(8) :: mxyz(3),idum
      integer :: nxyz(3),jxyz(3),ixyz(3),kxyz(3)
      logical :: lexist



      inquire(file=trim(data_file),exist=lexist)
      if (inode.eq.ionode) write(io_lun,*) "open rho: ",trim(data_file),lexist
      if (.not.lexist) then
        write(io_lun,*) "rho file not found ",trim(data_file)
        stop
      end if
      open(newunit=iunit2,file=trim(data_file),status="old",access="stream")

      read(iunit2,pos=8*(1-1)+1) idum
      read(iunit2,pos=8*(2-1)+1) idum
      read(iunit2,pos=8*(3-1)+1) mxyz(3)
      read(iunit2,pos=8*(4-1)+1) mxyz(2)
      read(iunit2,pos=8*(5-1)+1) mxyz(1)

      do jz=gs(3),ge(3)
        do jy=gs(2),ge(2)
          do jx=gs(1),ge(1)
            ixyz(1)=jx
            ixyz(2)=jy
            ixyz(3)=jz
            jxyz=ixyz
            jpos=jxyz(1)+(jxyz(2)-1)*mxyz(1)+(jxyz(3)-1)*mxyz(2)*mxyz(1)+11

            kxyz(1)=ixyz(1)-gs(1)+1
            kxyz(2)=ixyz(2)-gs(2)+1
            kxyz(3)=ixyz(3)-gs(3)+1
            read(unit=iunit2,pos=8*(jpos-1)+1,iostat=ierr) b(kxyz(1),kxyz(2),kxyz(3))
            if (ierr.ne.0) then
              write(0,*) inode," read error ",trim(data_file),ierr,jpos
              write(6,*) inode," read error ",trim(data_file),ierr,jpos
              stop
            end if

          end do
        end do
      end do

      close(iunit2)
      if (inode.eq.ionode) write(io_lun,*) "close rho: ",trim(data_file),lexist

    end subroutine read_rho

    subroutine gridsolver_init(isolver,bc,dxx,dyy,dzz)

      use global_module,    only: numprocs, io_lun,gridsolver_select,&
                                  gridsolver_bc
      use dimens,           only: r_super_x, r_super_y, r_super_z, &
     &                      n_grid_x,n_grid_y,n_grid_z
      use block_module,     only: nx_in_block,ny_in_block,nz_in_block
      use primary_module,   only: domain
      use datatypes
      use units,            only: dist_conv
      use input_module, only: io_assign
      use GenComms, ONLY: inode,gmax,ionode
!~  dl_mg grid solver ----
      use dl_mg_mpi_header
      use dl_mg


      implicit none

      integer,intent(in) :: bc(3),isolver
      real(double), optional :: dxx,dyy,dzz

!~ local
      integer :: ndims,pdims(3),ngx,ngy,ngz,&
     &           maxcycles,splane,ierr
      logical :: pperiods(3)
      character(DL_MG_MAX_ERROR_STRING) :: msg
      real(double) :: dx,dy,dz


      ndims=3
      maxcycles=100

      splane=gridsolver_slice
      pdims=1
      pdims(splane)=inmax !numprocs

      pperiods=(/ DL_MG_BC_PERIODIC.eq.bc(1) , DL_MG_BC_PERIODIC.eq.bc(2) ,&
     & DL_MG_BC_PERIODIC.eq.bc(3)   /)

      if (inode.eq.ionode) then
        write(io_lun,fmt='(A,3i8)') "gridsolver: slice",splane
        write(io_lun,fmt='(A,3i8)') "gridsolver: pdims",pdims
        write(io_lun,fmt='(A)') "gridsolver: boundary conditions (0=periodic,2=dirichlet)"
        write(io_lun,fmt='(3(A,i1))') "gridsolver: x=",bc(1)," y=",bc(2)," z=",bc(3)
        write(io_lun,fmt='(A,3l)') "gridsolver: periods",pperiods
      end if

      call mpi_cart_create(MPI_COMM_WORLD,3,pdims,pperiods,.false.,mg_comm,ierr)

      nglx=ge(1)-gs(1)+1
      ngly=ge(2)-gs(2)+1
      nglz=ge(3)-gs(3)+1


      if (present(dxx)) then
        dx=dxx
      else
        dx=dist_conv*(r_super_x/n_grid_x)
      end if
      if (present(dyy)) then
        dy=dyy
      else
        dy=dist_conv*(r_super_y/n_grid_y)
      end if
      if (present(dzz)) then
        dz=dzz
      else
        dz=dist_conv*(r_super_z/n_grid_z)
      end if
	    if (inode.eq.ionode) write(io_lun,fmt='(A,3e24.12)') "init dl_mg with dxyz",dx,dy,dz

      call io_assign(mg_unit)

      if (inode.eq.ionode) then        
        open(unit = mg_unit, file = "dl_mg_report.out", action = "write", status = "replace")
        close(mg_unit, status = "delete")
      end if

      if (inode.le.inmax) call dl_mg_init	(n_grid_x,n_grid_y,n_grid_z,dx,dy,dz,bc,gs,ge,mg_comm,&
      & mg_unit,"dl_mg_report.out",errors_return = .true.,ierror=ierr)
      if (ierr.ne.0) then
        call dl_mg_error_string(ierr,msg)
        write(io_lun,*) "dl_mg_init",trim(msg)
        stop
      end if
            
      allocate(pot3(1:nglx,1:ngly,1:nglz),stat=ierr)
      if (ierr.ne.0) then
        write(0,*) "allocation error ",ierr
        stop
      end if

      pot3=0d0

    end subroutine gridsolver_init

    subroutine gridsolver_free()

!~  dl_mg grid solver ----
      use dl_mg_mpi_header
      use dl_mg      
      use mpi
!~  ----------------------

      implicit none

      call dl_mg_free()

    end subroutine gridsolver_free

    subroutine get_hartree_dlmg(rho,pot)


      use global_module, only: io_lun,gridsolver_select, &
                               gridsolver_bc, restart_Knegf, &
                               negf_l_elec_dir,negf_r_elec_dir, &
                               dump_negf_data, flag_self_consistent, &
                               flag_neutral_atom
      use datatypes
      use units
      use GenComms, ONLY : my_barrier, gmax,inode,ionode
      use maxima_module,  only : maxngrid
      use input_module, only: io_close

!~  dl_mg grid solver ----
      use dl_mg_mpi_header
      use dl_mg

      implicit none

      real(double), intent(inout), target, contiguous :: rho(:)
      real(double), intent(inout), target, contiguous :: pot(:)

!~ local
      real(double), pointer :: rho3(:,:,:)
      real(double), allocatable :: eps_rel(:,:,:,:),res(:,:,:),eps(:,:,:)
      real(double) :: tol,hmax
      integer :: ierr,q0
      integer :: counti,count_rate,countf
      character(DL_MG_MAX_ERROR_STRING) :: msg
      character(256) :: file_hartree
      logical :: l_exist



      call rearrange(rho,maxngrid,gs_packbx,gs_psendbx,gs_prbx,gs_unpackbx)
      call rearrange(pot,maxngrid,gs_packbx,gs_psendbx,gs_prbx,gs_unpackbx)

      rho3(1:nglx,1:ngly,1:nglz)=>rho(1:nglx*ngly*nglz)
      pot3(1:nglx,1:ngly,1:nglz)=>pot(1:nglx*ngly*nglz)


      allocate(eps_rel(1:nglx,1:ngly,1:nglz,3),&
     &res(1:nglx,1:ngly,1:nglz),&
     &eps(1:nglx,1:ngly,1:nglz), stat=ierr)

!~ eps_rel is the relative permitivity which should here be equal to 1
!~ for now use dummy array as dl_mg requires it. This should be changed
!~ in the dl_mg source code
      eps_rel=1d0
      eps=1d0
      

      if (restart_Knegf) then
        if (flag_neutral_atom) then
          file_hartree = "./hartree_pot_diff.plt"            
        else
          file_hartree = "./hartree_pot.plt"            
        end if
        inquire(file=trim(file_hartree), exist = l_exist)          
        if (l_exist) call read_rho(pot3,trim(file_hartree))
      end if
              

      if (any(gridsolver_bc.eq.2)) then
        if (flag_neutral_atom) then
          call get_dirichlet(pot3(1:nglx,1:ngly,1:nglz),&
            trim(negf_l_elec_dir)//"/hartree_pot_diff.plt",&
            trim(negf_r_elec_dir)//"/hartree_pot_diff.plt")
        else
          call get_dirichlet(pot3(1:nglx,1:ngly,1:nglz),&
            trim(negf_l_elec_dir)//"/hartree_pot.plt",&
            trim(negf_r_elec_dir)//"/hartree_pot.plt")        
        end if
      end if
! call dl_mg_solver twice, first with low convergence tolerance to avoid convergence issues. then restart with a higher convergence threshold.
! this is more stable. if the convergence tolerance is not reached, do one addional run of dl_mg_solver. (this requieres a slightly modified version of dl_mg)
      if (inode.le.inmax) then

!~         pot3=-(pot3)/(4d0*pi)
        rho3=-rho3*4d0*pi
        
        res=0d0
        tol=1e-8

        call system_clock(counti,count_rate)

!~         call dl_mg_solver(eps=eps, eps_mid=eps_rel, alpha=1d0, rho=rho3, &
!~           pot=pot3, fd_order=4,	use_damping=.false., &
!~           tol_res_rel=tol,&
!~           tol_res_abs=tol,&
!~           tol_pot_rel=tol,&
!~           tol_pot_abs=tol,&          
!~           mg_max_conv_rate = 0.95d0, &
!~           max_iters_defco=150,&
!~           max_iters_vcycle=50,&
!~           v_iterations=(/ 2,1 /)*8,&
!~           omega_sor=1.9d0,&
!~           use_pot_in=.true.,res=res,ierror=ierr)
        call dl_mg_solver(eps=eps, eps_mid=eps_rel, alpha=1d0, rho=rho3, &
          pot=pot3, fd_order=4,	use_damping=.false., &  
          tol_pot_rel=tol,&
          tol_pot_abs=tol,&                 
          mg_max_conv_rate = 0.95d0, &
          max_iters_defco=150,&
          max_iters_vcycle=50,&
          v_iterations=(/ 2,1 /)*8,&
          omega_sor=1.9d0,&
          use_pot_in=.true.,res=res,ierror=ierr)
       
       
        if (ierr.ne.0) then
          call dl_mg_error_string	(ierr,msg)
          if (inode.eq.ionode)  write(io_lun,*) "dl_mg_solver: ",ierr,trim(msg)
        end if
      end if

      rho3=-rho3/(4d0*pi)

      hmax=maxval(abs(res))
      call gmax(hmax)
      call system_clock(countf)

      call my_barrier()

      call rearrange(pot,maxngrid,gs_packxb,gs_psendxb,gs_prxb,gs_unpackxb)
      call rearrange(rho,maxngrid,gs_packxb,gs_psendxb,gs_prxb,gs_unpackxb)
  
!~       pot=-pot*4d0*pi
      
      call io_close(mg_unit)

    end subroutine get_hartree_dlmg

    subroutine get_hartree_gridsolver(rho,pot)
      use global_module, only : gridsolver_select
      use datatypes

      implicit none

      real(double), intent(inout), target, contiguous :: rho(:)
      real(double), intent(inout), target, contiguous :: pot(:)

      if (gridsolver_select.eq.1) then
        call get_hartree_dlmg(rho,pot)
      else if (gridsolver_select.eq.2) then
!~       not yet
      end if

    end subroutine

end module grid_solver
