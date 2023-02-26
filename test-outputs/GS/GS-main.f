PROGRAM GS

  use sirisconstants
  use sirismath
  use sirisgeometry
  use sirismesh
  use sirisgaussiansphere
  
  integer :: lmin, lmax, nnod, ntri, ntr, j1, j2, j3, infu, nsample, seed
  integer, dimension(:,:), pointer :: IT
  real(kind=dp) :: nuc, sig, beta, rmax, vol, ell, cs2d, cs4d
  real(kind=dp), dimension(0:256) :: CSCF
  real(kind=dp), dimension(:), pointer :: MUN, PHIN
  real(kind=dp), dimension(:,:), allocatable :: XN0
  real(kind=dp), dimension(:,:), allocatable :: NT0
  real(kind=dp), dimension(0:256,0:256) :: ACF, BCF, SCFSTD
  character(len=file_name_length) :: fbn, cfn, infn, tstr
  logical :: tlog, outml=.false., outidl=.false., outvtk=.false., outoff=.false., &
    outobj=.false.

  ! Input file name as command argument
  if (command_argument_count() < 1) then
    write(error_unit, *) "Give input file name as the first command line argument, please."
    stop
  end if
  call get_command_argument(1,infn)
  inquire(file=infn, exist=tlog)
  if(.not. tlog) then
    write(error_unit, '(A,A,A)') "Input file '", trim(infn), " does not exist."
    stop
  end if

  ! Input file should be OK, open unit and read
  open(newunit=infu, file=infn, action='read')
  ! Three first lines are discarded
  read(infu,*)
  read(infu,*)
  read(infu,*)
  read(infu,*) nsample
  read(infu,*) seed
  read(infu,*) sig
  read(infu,*) nuc
  read(infu,*) lmin
  read(infu,*) lmax
  read(infu,*) ntr
  read(infu,*)
  read(infu,'(A)') fbn
  read(infu,*) j1
  if(j1==1) outml = .true.
  read(infu,*) j1
  if(j1==1) outidl = .true.
  read(infu,*) j1
  if(j1==1) outvtk = .true.
  read(infu,*) j1
  if(j1==1) outoff = .true.
  read(infu,*) j1
  if(j1==1) outobj = .true.
  
  write(output_unit,*) ""
  write(output_unit,'(A,A,A)') "Input read from file '", trim(infn), "':"
  write(output_unit,'(2X,A,I0)') "number of samples: ", nsample
  if(seed>0) then
    write(output_unit,'(2X,A,I0)') "random generator seed: ", seed
  else
    write(output_unit,'(2X,A)') "random generator seed: from system clock"
  end if
  write(output_unit,'(2X,A,G12.6)') "radius relative standard deviation: ", sig
  write(output_unit,'(2X,A,G9.4)') "power law index in correlation function: ", nuc
  write(output_unit,'(2X,A,I0,A,I0)') "min and max orders in correlation series: ", lmin, ", ", lmax
  write(output_unit,'(2X,A,I0)') "number of triangle rows per octant: ", ntr
  tstr = ""
  if(outml) write(tstr,'(A,A)') trim(tstr), " matlab "
  if(outidl) write(tstr,'(A,A)') trim(tstr), " IDF "
  if(outvtk) write(tstr,'(A,A)') trim(tstr), " VTK "
  if(outoff) write(tstr,'(A,A)') trim(tstr), " OFF "
  if(outobj) write(tstr,'(A,A)') trim(tstr), " OBJ "
  write(output_unit,'(2X,A,A)') "output formats selected: ", trim(tstr)
  write(output_unit,*) ""

  ! Input ready, initialize
  CALL init_random(seed)
  beta=sqrt(log(sig**2+1.0_dp))
  
  ! Correlation function
  call cs1cf(CSCF,nuc,lmin,lmax)
  call csini(CSCF,ell,cs2d,cs4d,lmax)

  ! Spherical harmonics
  call sgscfstd(SCFSTD,CSCF,beta,lmax)
  ! Spherical discretization
  call trids(MUN,PHIN,IT,nnod,ntri,ntr)

  ! Allocate vertex tables
  allocate(XN0(nnod,3), NT0(ntri,3))
  
  do j1=1,nsample

    ! Generates the sample spherical harmonics coefficients
    call sgscf(ACF,BCF,SCFSTD,lmax)
    ! Discretized triangle representation
    call rgstd(XN0,NT0,MUN,PHIN,ACF,BCF,rmax,beta,IT,nnod,ntri,lmax)

    ! Compute volume
    vol = compute_volume(XN0,IT,ntri)
    write(output_unit,'(A,I0,A,G12.6)') "Sample ", j1, ", volume ", vol
        
    ! Save geometry
    if(nsample == 1) then
      write(cfn,'(A)') trim(fbn)
    else
      write(cfn,'(A,I0)') trim(fbn), j1
    end if
    ! Clean XN0
    do j2=1,nnod
      do j3=1,3
        if(abs(XN0(j2,j3)) < 1e-12) XN0(j2,j3) = 0.0_dp
      end do
    end do
    if(outml) call save_matlab(cfn,XN0,nnod)
    if(outidl) call save_idl(cfn,XN0,IT,nnod,ntri)
    if(outvtk) call save_vtk(cfn,XN0,IT,nnod,ntri)
    if(outoff) call save_off(cfn,XN0,IT,nnod,ntri)
    if(outobj) call save_obj(cfn,XN0,IT,nnod,ntri)

  end do

END PROGRAM GS
