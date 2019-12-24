! -*- mode: F90; mode: font-lock; column-number-mode: true; vc-back-end: CVS -*-
! ------------------------------------------------------------------------------
! $Id$
! ------------------------------------------------------------------------------
! Module ionic_data
! ------------------------------------------------------------------------------
! Code area 1: initialisation
! ------------------------------------------------------------------------------

!!****h* Conquest/ionic_data *
!!  NAME
!!   ionic_data
!!  PURPOSE
!!   Collects the reading and processing of ionic data into one module
!!  USES
!!
!!  AUTHOR
!!   M.J.Gillan
!!  CREATION DATE
!!   June 2002
!!  MODIFICATION HISTORY
!!   14:29, 17/09/2002 drb 
!!    Added module header and a check for gauss/pao before calling read_pao
!!   13:56, 2003/09/22 dave
!!    Added read for basis set: defaults to blips
!!   2008/02/06 08:08 dave
!!    Changed for output to file not stdout
!!   2015/06/08 lat
!!    Added experimental backtrace
!!  SOURCE
!!
module ionic_data

  use global_module, only: io_lun
  use timer_module,  only: start_timer, stop_timer, cq_timer 
  use timer_module,  only: start_backtrace, stop_backtrace

  implicit none

  ! Area identification
  integer, parameter, private :: area = 1

  ! RCS tag for object file identification
  character(len=80), private, save :: &
       RCSid = "$Id$"
!!***

contains

  !!****f* ionic_data/get_ionic_data *
  !!
  !!  NAME 
  !!   get_ionic_data
  !!  USAGE
  !! 
  !!  PURPOSE
  !!   Gets ionic data - will read PAOs, charge density and 
  !!   pseudopotentials
  !!  INPUTS
  !! 
  !! 
  !!  USES
  !! 
  !!  AUTHOR
  !!   M.J.Gillan
  !!  CREATION DATE
  !!   June 2002
  !!  MODIFICATION HISTORY
  !!   14:31, 17/09/2002 drb 
  !!    Added headers and a check for PAOs before reading them
  !!   11:16, 24/09/2002 mjg & drb 
  !!    Added lines to check on whether PAOS are needed
  !!   12:01, 24/09/2002 mjg & drb 
  !!    Added calls to make or read atomic densities (and used
  !!    atomic_density module)
  !!   15:33, 25/09/2002 mjg & drb 
  !!    Added flag_no_atomic_densities to specify that atomic
  !!    densities have not been made
  !!   07:56, 2003/09/22 dave
  !!    Added flag to choose basis set
  !!   08:54, 26/05/2005 dave 
  !!    Added various bits to stop PAO reading if SIESTA pseudos are
  !!    used - the PAOs are now read during .ion file reading.
  !!    **NOTE** we assume that the FIRST call to init_pseudo_siesta
  !!    is during read_and_write and that the tables are read in there
  !!   10:55, 13/02/2006 drb 
  !!    Line length change
  !!   16:00, 2017/02/21 nakata
  !!    Commented out get_support_pao_rep which is no longer used
  !!   15:00, 2017/03/08 nakata
  !!    Removed get_support_pao_rep which is no longer used
  !!  SOURCE
  !!
  subroutine get_ionic_data(inode,ionode,flag_no_atomic_densities,level)

    use datatypes
    use read_pao_info,          only: read_pao
    use atomic_density,         only: make_atomic_density_from_paos, &
                                      spline_atomic_density,         &
                                      atomic_density_method
    use species_module,         only: n_species
    use GenComms,               only: cq_abort
    use global_module,          only: flag_basis_set, blips, PAOs,   &
                                      iprint_init
    use pseudopotential_common, only: pseudo_type, SIESTA, ABINIT
    use blip,                   only: init_blip_flag
    use input_module,           only: leqi

    implicit none

    ! Passed variables
    integer, optional    :: level
    integer, intent(in)  :: inode, ionode
    logical, intent(out) :: flag_no_atomic_densities

    ! Local variables
    character(len=10) :: init_blip_method
    logical           :: flag_blips_from_pao
    type(cq_timer)    :: backtrace_timer
    integer           :: backtrace_level

!****lat<$    
    if (       present(level) ) backtrace_level = level+1
    if ( .not. present(level) ) backtrace_level = -10
    call start_backtrace(t=backtrace_timer,who='get_ionic_data',&
         where=area,level=backtrace_level)
!****lat>$

! 2019/Dec/24 tsuyoshi
!   Now we assume that we always have *.ion files, thus
!    flag_atomic_density_from_pao = .true.
!        flag_no_atomic_densities = .false.
!   (pseudo_type) = SIESTA or ABINIT 
!   (pao)
!   (blips) 
!    flag_blips_from_pao = .true. or .false.
!         init_blip_flag = 'pao' or 'gauss'
!      -> removed flag_blips_from_pao  (we don't need it anymore)
! 2019/Dec/24 tsuyoshi

! We will remove these two flags in the near future 
!       flag_atomic_density_from_pao = .true.
            flag_no_atomic_densities = .false.
!
!! NEW GET_IONIC_DATA : start
       call make_atomic_density_from_paos(inode, ionode, n_species)
       call spline_atomic_density(n_species)
!! NEW GET_IONIC_DATA : end

!    ! Decide whether blips are to be initialised from PAOs
!    flag_blips_from_pao = .false.
!    if (leqi(init_blip_flag,'pao')) flag_blips_from_pao = .true.
!    ! Decide whether atomic densities are to be initialised from PAOs
!    flag_atomic_density_from_pao = .false.
!    if (leqi(atomic_density_method, 'pao')) &
!         flag_atomic_density_from_pao = .true.
!    ! If PAOs are needed, then read them
!    ! Note that we now read PAOs from the .ion file if we're using
!    ! that pseudopotential
!    if ((pseudo_type /= SIESTA .and. pseudo_type /= ABINIT) .and. &
!         (flag_blips_from_pao .or. flag_atomic_density_from_pao .or. &
!          flag_basis_set == PAOs)) then
!
!       if (inode == ionode .and. iprint_init > 1) &
!            write (unit=io_lun, fmt='(//" **** get_ionic_data: about &
!                                      &to call read_pao_info")')
!       call read_pao(inode, ionode, n_species)
!
!    else if ((pseudo_type == SIESTA .or. pseudo_type == ABINIT) .and. &
!             (flag_blips_from_pao .or. flag_atomic_density_from_pao .or. &
!              flag_basis_set == PAOs)) then
!
!       if (inode == ionode .and. iprint_init > 1) &
!            write (unit=io_lun, fmt='(//" **** get_ionic_data: PAO &
!                                      &input from init_pseudo_tm already done")')
!
!    end if
!
!    ! Get the atomic densities
!    if (flag_atomic_density_from_pao) then ! Use PAOs for atomic densities
!       call make_atomic_density_from_paos(inode, ionode, n_species)
!    else if (.not. flag_atomic_density_from_pao) then
!       if (leqi(atomic_density_method, 'read')) then ! User has specified a file
!          call read_atomic_density(inode, ionode, n_species)
!       else 
!          ! Signal to the calling routine that no atomic densities
!          ! have been created
!          flag_no_atomic_densities = .true.
!       end if
!    end if
!
!    if (.not. flag_no_atomic_densities) then ! Spline tables
!       call spline_atomic_density(n_species)
!    end if

!****lat<$
    call stop_backtrace(t=backtrace_timer,who='get_ionic_data')
!****lat>$

    return
  end subroutine get_ionic_data
  !!***

end module ionic_data





