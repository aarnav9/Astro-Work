PROGRAM singleparticle

  use sirisconstants
  use sirismath
  use sirismaterial
  use sirisgeometry
  use sirismesh
  use sirisgaussiansphere
  use sirisradtrans
  use sirisray
  
  integer  :: nbin, j, jray, jsub, jbox, pout, nray, seed1, &
    ntri,j1,j2,totref,lmax, npar,pin,pinmax,poutmax,ntr,nis, savenum, &
    l1, flgin, np, nrn, pflg,ncm, nnod, nsub, fu, run_num, &
    extmeshtype = -1, cstat
  integer, dimension(50) :: PBOX, NISBOX
  integer, dimension(:,:), pointer :: IT
  logical :: ignabs = .false., extmesh = .false., trans_vert = .false., scale_vert = .false.
  real(kind=dp):: wrayh, wraym, norm, ran2, mui, nui, phii, sphii, cphii, &
    rmax, r0, phi0, qsca,qabs,qin,qout,qis,qbox,abscf,len,lenabs,kekf,beta,nie, &
    nk,Fstop,bin,bin0,norml,wray, wavelen, xa, m2real, &
    m2imag, qscaG, qabsG, qextG, renorm, lin, otin, lenrn, ghg1, ghg11, ghg21, &
    whg1, pmx1, dthe, pnorm1, radius, mmed, sig, nuc, cs2d, cs4d, &
    ell, fintmp
  real(kind=dp), dimension(2):: MA0, MAOUT, MAIN
  real(kind=dp), dimension(3) :: KEOUT, KFOUT, X, ELOUT, EROUT, Y, N, KEIN, KFIN, &
    ELIN, ERIN
  real(kind=dp), dimension(0:256) :: CSCF
  real(kind=dp), dimension(361) :: XP1
  real(kind=dp), dimension(:), allocatable :: CSRN1
  real(kind=dp), dimension(:), pointer :: MUN, PHIN
  real(kind=dp), dimension(2,50) :: MABOX  
  real(kind=dp), dimension(3,50) :: KEBOX, KFBOX, XBOX
  real(kind=dp), dimension(4,4) :: FOUT, FIN  
  real(kind=dp), dimension(0:256,0:256) :: ACF, BCF, SCFSTD
  real(kind=dp), dimension(:,:), allocatable :: XN, NT
  real(kind=dp), dimension(:,:), pointer :: XN0, NT0
  real(kind=dp), dimension(4,4,50) :: FBOX
  real(kind=dp), dimension(361,4,4) :: S, YP1, YP21
  real(kind=dp), dimension(0:360,4,4) :: P1
  complex(kind=dp):: m1,m2
  complex(kind=dp), dimension(3) :: HLOUT, HROUT, HLIN, HRIN
  complex(kind=dp), dimension(3,50) :: HLBOX, HRBOX    
  character(len=1) :: inhmg
  character(len=file_name_length) :: fname, fname2, infile, meshfile, tempfn

  ! Read inputs from file:
  if(command_argument_count() == 0) then
    infile = 'Input.in'
  else
    call get_command_argument(1,infile)
  endif
  open(newunit=fu, file=infile, status='old', action='read')

  read(fu, *) nray        ! Number of rays
  read(fu, *) npar        ! Number of sample particles.
  read(fu, *) Fstop       ! Minimum relative flux
  read(fu, *) pinmax      ! Max internal chord 
  read(fu, *) poutmax     ! Max external chord
  read(fu, *) nbin        ! Number of scattering angle bins
  read(fu, *) mmed        ! Refractive index of the medium
  read(fu, *) m2real      ! Particle refractive index (real)
  read(fu, *) m2imag      ! Particle refractive index (imag)
  read(fu, *) seed1       ! Seed for random number generation 
  read(fu, *) wavelen     ! Wavelength
  read(fu, *) radius      ! Mean radius (microns)
  read(fu, *) run_num     ! Identification number for run -> output file is ice_go_'run'.out 
  !read(fu, *) inhmg       ! use inhomogeneous waves (y/n)  
  inhmg="y"
  read(fu, *) flgin       ! Internal medium: 1=yes
  read(fu, *) pflg        ! Scattering matrix: input(1), Rayleigh(2)
  read(fu, *) otin        ! Internal medium: single-scattering albedo
  read(fu, *) lin         ! Internal medium: mean free path (in micrometers)
  read(fu, *) ghg1        ! Diffuse-medium asymmetry parameter (internal)
  read(fu, *) ghg11       ! Diffuse-medium forward asymmetry (internal)
  read(fu, *) ghg21       ! Diffuse-medium backward parameter (internal)
  read(fu, *) pmx1        ! Diffuse-medium single-scattering max. polarization
  read(fu, *) np          ! Number of phase matrix angular points
  read(fu, *) nrn         ! Number of points in random number array
  read(fu, *) ncm         ! Number of G-L integration points for cosine map
  read(fu, *) sig         ! Relative standard deviation of radius.
  read(fu, *) nuc         ! Power law index for C_3 correlation.
  read(fu, *) lmax        ! Maximum degree in C_1, C_2, C_3.
  read(fu, *) ntr         ! Number of triangle rows in an octant.
  read(fu, *) fname       ! Prefix for the output file names
  ! Optional, read external geometry mesh file
  read(fu, *, IOSTAT=cstat) extmeshtype
  if(cstat==0 .and. extmeshtype > 0) then
    extmesh = .true.
    read(fu, *) j1   ! Translate vertex mean to origin
    if(j1==1) trans_vert = .true.
    read(fu, *) j1    ! Scale vertex mean radius to one
    if(j1==1) scale_vert = .true.
    read(fu, '(A)') tempfn  ! External mesh geometry file name
    call strip_string(tempfn,meshfile,j1)
    npar=1
  end if
  close(fu)
  
  if(extmesh) then
    write(output_unit,'(A)') 'Geometric optics approximation for triangulated mesh particle...'
    write(output_unit,'(A,A,A)') "Particle given in file '", trim(meshfile), "'"
  else
    write(output_unit,'(A)') 'Geometric optics approximation for a Gaussian random sphere...'
    write(output_unit,'(A,F7.1,A)') 'Particle mean radius ', radius, ' micrometers'
  end if
  write(output_unit,'(A,I0)') 'Number of rays: ', nray
    

  call init_random(seed1)

  nsub=nray/npar
  if(nsub*npar /= nray) stop 'Trouble in GEO: ray/particle number mismatch.'

  ! Refractive indices:
  m2 = cmplx(m2real, m2imag, dp)
  m1 = cmplx(mmed, 0.0_dp, dp)   ! ref index of the medium
  call refindapp(MA0,m1,1.0d0)

  ! Initialization of the Gaussian random sphere:
  if(.not. extmesh) then
    beta=sqrt(log(sig**2+1.0_dp))
    call cs1cf(CSCF,nuc,2,lmax)  !power law
    call csini(CSCF,ell,cs2d,cs4d,lmax)

    call sgscfstd(SCFSTD,CSCF,beta,lmax)

    whg1 = (ghg1-ghg21)/(ghg11-ghg21)

    ! Discretization:
    call TRIDS(MUN,PHIN,IT,nnod,ntri,ntr)
    allocate(XN(nnod,3), XN0(nnod,3), NT(ntri,3), NT0(ntri,3))
  
    call sgscf(ACF,BCF,SCFSTD,lmax)
    call rgstd(XN0,NT0,MUN,PHIN,ACF,BCF,rmax,beta,IT,nnod,ntri,lmax)
    
  ! Initialization of the external mesh geometry:
  else
  
    if(extmeshtype == 1) then
      call read_off(meshfile,XN0,IT,NT0,nnod,ntri,cstat)
    else if(extmeshtype == 2) then
      call read_obj(meshfile,XN0,IT,NT0,nnod,ntri,cstat)
    else
      write(error_unit,'(A,I0)') 'Unknown geometry type: ', extmeshtype
      stop
    end if
    if(cstat /= 0) then
      write(error_unit,*) 'Error in reading the mesh geometry from file'
      stop
    end if
    
    allocate(XN(nnod,3), NT(ntri,3))
    
    if(trans_vert) then
      ! Translate vertex mean to origin
      X(:) = 0.0_dp
      do j1=1,nnod
        X(1) = X(1) + XN0(j1,1)
        X(2) = X(2) + XN0(j1,2)
        X(3) = X(3) + XN0(j1,3)
      end do
      X(1) = X(1)/nnod
      X(2) = X(2)/nnod
      X(3) = X(3)/nnod
      do j1=1,nnod
        XN0(j1,1) = XN0(j1,1)-X(1)
        XN0(j1,2) = XN0(j1,2)-X(2)
        XN0(j1,3) = XN0(j1,3)-X(3)
      end do
    end if
    
    if(scale_vert) then
      ! Scale vertices to mean radius of one
      r0 = 0.0_dp
      do j1=1,nnod
        r0 = r0 + sqrt(XN0(j1,1)**2+XN0(j1,2)**2+XN0(j1,3)**2)
      end do
      r0 = r0/nnod
      do j1=1,nnod
        XN0(j1,1) = XN0(j1,1)/r0
        XN0(j1,2) = XN0(j1,2)/r0
        XN0(j1,3) = XN0(j1,3)/r0
      end do
    end if

    ! Compute maximum radius
    rmax = 0.0_dp
    do j1=1,nnod
      r0 = XN0(j1,1)**2+XN0(j1,2)**2+XN0(j1,3)**2
      if(rmax < r0) rmax = r0
    end do
    rmax = sqrt(rmax)

  end if
  
  ! Size parameter, absorption parameter, scaling of free path:
  xa = 2.0_dp*pi*radius/wavelen
  abscf = -2.0_dp*xa*m2imag
  lin=lin/radius

  ! Initialize:
  qsca = 0.0_dp
  qabs = 0.0_dp
  qin = 0.0_dp
  qout = 0.0_dp
  qis = 0.0_dp
  qbox = 0.0_dp
  S = 0.0_dp

  ! Initialize the diffuse scattering phase matrices:
  if(flgin > 0) then
    dthe=pi/np
    allocate(CSRN1(0:nrn))
    if (pflg == 1) then
      call pmatrix1_single(P1,np)
      call psplivi(P1,CSRN1,XP1,YP1,YP21,pnorm1,dthe,nrn,np,ncm)
    end if
  end if

  ! Storage grid parameters:
  bin = pi/nbin
  bin0 = cos(0.5_dp*bin)

  ! Externally propagating rays. Initialization:
  jray = 0
  jsub = 0
  jbox = 0
  wrayh = 0.0_dp ! rays that hit
  wraym = 0.0_dp ! rays that miss


100 if(jbox == 0) then
    pout = 0
    jray = jray + 1

    ! Print status
    if(mod(jray,5000) == 0) then
      write(output_unit,'(A19,I0,A1,I0)') " computing at ray ", jray, "/", nray
    end if

    ! Normalizations:
    if(jray > nray) then ! number of rays done --> quit program
      norm = nray/wrayh
      qsca = qsca*norm
      qabs = qabs*norm
      qin = qin*norm
      qout = qout*norm
      qis = qis*norm
      qbox = qbox*norm
      S = S*norm
      qscaG = qsca/nray           
      qabsG = qabs/nray
      qextG = qscaG + qabsG 

      renorm = 1.0_dp/qscaG
      norm = renorm*2.0_dp/(nray*(1.0_dp-bin0))
       
      do j1 = 1, 4
        do j2 = 1, 4
          S(1,j1,j2) = norm*S(1,j1,j2) ! forward scattering
          S(nbin+1,j1,j2) = norm*S(nbin+1,j1,j2) ! backscattering
        end do
      end do
       
      do j = 1, nbin-1           
        norm = renorm*2.0_dp/(nray*(cos((j-0.5_dp)*bin)-cos((j+0.5_dp)*bin)))      
        do j1 = 1, 4
          do j2 = 1, 4
            S(j+1,j1,j2) = norm*S(j+1,j1,j2)
          end do
        end do
      end do    
       
      write(output_unit,'(A)')           '----------------'
      write(output_unit,'(A)')           'Results:'
      write(output_unit,'(A,E12.4)') 'Qsca = ', qscaG
      write(output_unit,'(A,E12.4)') 'Qabs = ', qabsG
      write(output_unit,'(A,F8.5)')  'SSA = ', qscaG/qextG
      write(output_unit,'(A,F20.1)') 'Rays hit = ', wrayh
      write(output_unit,'(A,F20.1)') 'Rays miss = ', wraym

      write(fname2,'(A,I0,A)') trim(fname), run_num, '-S.out' ! output
      write(output_unit,'(A,A,A)') "Scattering matrix elements written in file '", trim(fname2), "'"
   
      open(newunit=fu, file=trim(fname2))
      savenum = 0
      do l1 = 1, nbin+1
        write (fu,'((F5.1),1X,16(E16.7, 1X))') (l1-1)*180.0/nbin, S(l1,1,1), S(l1,1,2), S(l1,1,3), S(l1,1,4), &
          S(l1,2,1), S(l1,2,2), S(l1,2,3), S(l1,2,4), &
          S(l1,3,1), S(l1,3,2), S(l1,3,3), S(l1,3,4), &
          S(l1,4,1), S(l1,4,2), S(l1,4,3), S(l1,4,4)     
        savenum = savenum + 1
      end do
      close (fu)

      write(fname2,'(A,I0,A)') trim(fname), run_num, '-pmat.out' ! output
      write(output_unit,'(A,A,A)') "...and to file '", trim(fname2), "'"

      open(newunit=fu, file=trim(fname2))
      savenum = 0
      do l1 = 1, nbin+1
        write (fu,'((F5.1),1X,16(E16.7, 1X))') (l1-1)*180.0/nbin, S(l1,1,1), S(l1,1,2), &
          S(l1,2,2), S(l1,3,3), S(l1,3,4), S(l1,4,4)     
        savenum = savenum + 1
      end do
      close(fu)

      write(fname2,'(A,I0,A)') trim(fname), run_num, '-Q.out' ! output
      write(output_unit,'(A,A,A)') "Scattering efficiencies written in file '", trim(fname2), "'"

      open(newunit=fu, file=trim(fname2))
      write(fu,'(A,F15.8)') 'qsca = ',qscaG
      write(fu,'(A,F15.8)') 'qabs = ',qabsG
      write(fu,'(A,F15.8)') 'qext = ',qextG
      write(fu,'(A,F15.8)') 'rhit = ',sqrt(wrayh/nray)
      write(fu,'(A,F15.8)') 'qsca/qext = ',qscaG/qextG
      write(fu,'(A,F15.8)') 'wavelen = ',wavelen
      write(fu,'(A,F15.8)') 'meanRadius = ',radius
      write(fu,'(A,F15.8)') 'mmed = ',mmed
      write(fu,'(A,F15.8)') 'mreal = ',m2real
      write(fu,'(A,F15.8)') 'mimg = ',m2imag
      write(fu,'(A,F15.8)') 'meanFreePath = ',lin*radius
      write(fu,'(A,F15.8)') 'qin = ',qin/nray
      write(fu,'(A,F15.8)') 'qout = ',qout/nray
      write(fu,'(A,F15.8)') 'qis = ',qis/nray
      write(fu,'(A,F15.8)') 'qbox = ',qbox/nray
      close(fu)
      stop
    end if

    ! rayinic initializes the Mueller matrix, wave vector, and 
    ! parallel/perpendicular coordinate axes.
    MAOUT = MA0 
    call rayinic(FOUT,KEOUT,KFOUT,HLOUT,HROUT)

    call random_number(ran2)
    mui = 1.0-2.0*ran2 ! -1...1 ! tama on cos(theta), tasaisesti jakautunut -1...1
    nui = sqrt(1.0 - mui**2) ! taman on sin(theta)
    call random_number(ran2)          
    phii = 2.0*pi*ran2
    cphii = cos(phii)
    sphii = sin(phii)

    ! Rotation of KE, KF, HL, and HR to the particle coordinate system:
    call raypart0(KEOUT,mui,nui,cphii,sphii)
    call raypart0(KFOUT,mui,nui,cphii,sphii)
    call raypart0(HLOUT,mui,nui,cphii,sphii)
    call raypart0(HROUT,mui,nui,cphii,sphii)

    ! Generation of the spherical harmonics coefficients for the
    ! logradius. Discretization and random orientation of the Gaussian
    ! sample sphere. Not with external mesh:
    if(.not. extmesh) then
      jsub=jsub+1
      if(jsub > nsub) then
        call sgscf(ACF,BCF,SCFSTD,lmax)
        call rgstd(XN0,NT0,MUN,PHIN,ACF,BCF,rmax,beta,IT,nnod,ntri,lmax)
        jsub=1
      endif
    end if

    call pranor(XN,NT,XN0,NT0,nnod,ntri)

    ! Weight of incident ray.
    wray = rmax**2

    ! Offset of incident ray:
    call random_number(ran2)
    r0 = rmax*sqrt(ran2)
    call random_number(ran2)
    phi0 = 2.0*pi*ran2
    X = (/ r0*cos(phi0), r0*sin(phi0), -sqrt(rmax**2-r0**2) /) ! - vai ei z-komponenttiin...?

    call raypart0(X,mui,nui,cphii,sphii) ! X is the location where the incident ray starts (particle coordinates)

    ! istri determines the possible interaction point on the particle surface.
    nie = 1.0_dp
    call istri(X,KEOUT,N,XN,NT,IT,nk,len,ntri,nis,nie)
    if(nis == 0) then
      wraym = wraym+wray
      goto 100
    endif
    wrayh = wrayh+wray
    do j1 = 1, 4
      do j2 = 1, 4
        FOUT(j1,j2) = wray*FOUT(j1,j2)
      end do
    end do
   
  else

    ! RAYGETD extracts one ray from the pile of external rays. ISTRI 
    ! determines whether that ray further interacts with the particle 
    ! surface. If it doesn't, SCATTER stores the ray among other scattered 
    ! rays.
    call raygetd1c(FBOX,KEBOX,KFBOX,HLBOX,HRBOX,XBOX,MABOX,PBOX,NISBOX,FOUT, &
      KEOUT,KFOUT,HLOUT,HROUT,X,MAOUT,pout,nis,jbox)
    jbox = jbox - 1

    if(FOUT(1,1) <= wray*Fstop .or. pout > poutmax) then
      qout = qout+FOUT(1,1)
      goto 100
    end if
    N = NT(nis, 1:3)
    call prosca(nk,N,KEOUT)

  end if

  ! Scattering or external incidence:
  if (nk >= 0.0_dp) then ! jos kyseessa sisaltapain ulostulo
    norml = 0.0_dp
    do j1 = 1, 3
      KFOUT(j1)=KEOUT(j1)
      norml=norml+real(HLOUT(j1),dp)**2
    end do

    norml = sqrt(norml)

    do j1 = 1, 3
      ELOUT = real(HLOUT(1:3),dp)/norml
    end do
    call provec(EROUT,KEOUT,ELOUT)

    do j1=1,3
      HLOUT(j1)=cmplx(ELOUT(j1),0.0_dp,dp)
      HROUT(j1)=cmplx(EROUT(j1),0.0_dp,dp)
    end do
    MAOUT(1)=MA0(1)
    MAOUT(2)=MA0(2)
    nie = 1.0
    call istri(X,KEOUT,N,XN,NT,IT,nk,len,ntri,nis,nie) ! osuuko enaa uuteen kolmioon?
    if(nis == 0) then ! ei osu. tallennetaan sade sirontamatriisiin.
      call partray0(KEOUT,mui,nui,cphii,sphii)
      call partray0(ELOUT,mui,nui,cphii,sphii)
      call partray0(EROUT,mui,nui,cphii,sphii)
      call scatter(S,FOUT,KEOUT,ELOUT,EROUT,qsca,bin,bin0,nbin)
      goto 100
    end if

  end if

  ! INCIDE determines the reflected and refracted Mueller matrices, 
  ! wave vectors, and parallel/perpendicular coordinate axes.
  pout = pout+1
  if(inhmg == 'y') then
    call incide(FOUT,KEOUT,KFOUT,HLOUT,HROUT,MAOUT, &
      FIN,KEIN,KFIN,HLIN,HRIN,MAIN, &
      N,m1,m2,totref)
  end if

  ! Internally propagating rays. RAYPUTD stores the reflected ray into 
  ! the pile of external rays.
  pin = pout  
300 jbox = jbox+1
  if(jbox > 49) then
    qbox = qbox+FIN(1,1)
    goto 100
  end if
  call rayputd1c(FBOX,KEBOX,KFBOX,HLBOX,HRBOX,XBOX, &
    MABOX,PBOX,NISBOX,FOUT,KEOUT,KFOUT, &
    HLOUT,HROUT,X,MAOUT,pout,nis,jbox)
  if(pinmax == 0) then
    qabs = qabs+FIN(1,1)
    goto 100
  end if

310 if(FIN(1,1) <= wray*Fstop .or. pin >= pinmax) then
    qin = qin+FIN(1,1)
    goto 100
  end if

  ! ISTRI determines the next interaction point for the internal ray, 
  ! computes the inner normal at the interaction point, ABSORB determines the 
  ! amount of absorbed flux, and INCIDE determines the reflected and refracted 
  ! Mueller matrices, wave vectors, and parallel/perpendicular coordinate axes.
  Y = X

  nie = -1.0_dp
  call istri(Y,KEIN,N,XN,NT,IT,nk,len,ntri,nis,nie)
  if(nis == 0) then ! tahan ei kai pitais paatya mitaan? jos paatyy, niin kolmioinnissa vika. (tjs.)
    qis = qis+FIN(1,1)
    goto 100
  end if

  !-------------------
  !Internal diffuse medium
  if(flgin == 1) then

    call random_number(ran2)
    lenrn=-lin*log(ran2)
    if(lenrn < len) then

      pin=pin+1
      pout=pout+1

      do j1 = 1, 3
        X(j1)=X(j1)+lenrn*KEIN(j1)
      end do

      call prosca(kekf,KEIN,KFIN)
      lenabs = lenrn*abs(kekf)

      if(.not. ignabs) call absorb(FIN,qabs,lenabs,abscf)

      fintmp=FIN(1,1)
      if(FIN(1,1) <= wray*Fstop) then
        qin=qin+fintmp
        goto 100
      endif

      call absorbrt(FIN,qabs,otin)

      norml = 0.0
      do j1 = 1, 3
        norml=norml+real(HLIN(j1),dp)**2
      end do

      norml = sqrt(norml)

      ELIN = real(HLIN(1:3),dp)/norml

      call provec(ERIN,KEIN,ELIN)

      call scart(FIN,KEIN,ERIN,ELIN,CSRN1,XP1,YP1,YP21,np,nrn,pflg)

      do j1 = 1, 3
        KFIN(j1)=KEIN(j1)
        HLIN(j1)=cmplx(ELIN(j1),0.0_dp,dp)
        HRIN(j1)=cmplx(ERIN(j1),0.0_dp,dp)
      end do

      goto 310
    endif

  endif

  !------------------
  X = Y

  call prosca(kekf,KEIN,KFIN)

  lenabs = len*abs(kekf)

  if(.not. ignabs) call absorb(FIN,qabs,lenabs,abscf)

  call incide(FIN,KEIN,KFIN,HLIN,HRIN,MAIN, &
    FOUT,KEOUT,KFOUT,HLOUT,HROUT,MAOUT, &
    N,m2,m1,totref)
  pin = pin+1
  pout = pout+1
  if(totref == 1) goto 310 ! total internal reflection
  if(totref == -1) goto 300

END PROGRAM singleparticle
