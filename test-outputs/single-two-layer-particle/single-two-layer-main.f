PROGRAM singletwolayer

  use sirisconstants
  use sirismath
  use sirismaterial
  use sirisgeometry
  use sirisgaussiansphere
  use sirisradtrans
  use sirisray
  
  integer :: nbin, j, jray, jsub, jbox, jboxin, pout, nray, seed1, ntri, &
    j0, j1, j2, totref, lmax, npar, pin, pinmax, poutmax, ntr, &
    nis, nis2, savenum, l1, rtflg1, rtflg2, np, nrn, ncm, nnod, nsub, &
    run_num, fu
  integer, dimension(500) :: PBOX, NISBOX, PINBOX, NISINBOX
  integer, dimension(:,:), pointer :: IT
  real(kind=dp):: wrayh, wraym, norm, ran2, mui, nui, phii, sphii, &
    cphii, rmax, r0, phi0, qsca, qabs, qin, qout, qis, qbox, abscf1, abscf2, &
    len, len1, len2, len3, lenabs, kekf, beta, nk1, nk2, Fstop, bin, bin0, &
    wray, wavelen, xa, m1real, m1imag, m2real, m2imag, qscaG, qabsG,qextG, &
    renorm, lenrt1, lenrt2, omg1, omg2, lenrn, dthe, pnorm1, pnorm2, radius1, &
    radius2, mmed, sig, nuc, rho, cs2d, cs4d, ell, rsc, nie, nie2
  real(kind=dp), dimension(2):: MA0, MAOUT, MAIN, MACORE
  real(kind=dp), dimension(3) :: KEOUT, KFOUT, X, Y1, Y2, ELOUT, EROUT, N1, N2, &
    KEIN, KFIN, KECORE, KFCORE
  real(kind=dp), dimension(361) :: XP1, XP2
  real(kind=dp), dimension(0:256) :: CSCF
  real(kind=dp), dimension(:), allocatable :: CSRN1, CSRN2
  real(kind=dp), dimension(:), pointer :: MUN, PHIN
  real(kind=dp), dimension(2,500) :: MABOX, MAINBOX
  real(kind=dp), dimension(3,500) :: KEBOX, KFBOX, XBOX, KEINBOX, KFINBOX, XINBOX
  real(kind=dp), dimension(4,4) :: FOUT, FIN, FCORE  
  real(kind=dp), dimension(0:256,0:256) :: ACF, BCF, SCFSTD
  real(kind=dp), dimension(:,:), allocatable :: XN1, XN2
  real(kind=dp), dimension(:,:), allocatable :: NT
  real(kind=dp), dimension(4,4,500) :: FBOX, FINBOX
  real(kind=dp), dimension(361,4,4) :: S, YP1, YP2, YP21, YP22
  real(kind=dp), dimension(0:360,4,4) :: P1, P2
  complex(kind=dp) :: m0, m1, m2
  complex(kind=dp), dimension(3) :: HLOUT, HROUT, HLIN, HRIN, HLCORE, HRCORE
  complex(kind=dp), dimension(3,500) :: HLBOX, HRBOX, HLINBOX, HRINBOX
  character(len=file_name_length) :: fname, fname2, infile

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
  read(fu, *) m1real      ! Mantle: Particle refractive index (real)
  read(fu, *) m1imag      ! Mantle: Particle refractive index (imag)
  read(fu, *) m2real      ! Core: Particle refractive index (real)
  read(fu, *) m2imag      ! Core: Particle refractive index (imag)
  read(fu, *) seed1       ! Seed for random number generation 
  read(fu, *) wavelen     ! Wavelength
  read(fu, *) radius1     ! Outer radius (microns)
  read(fu, *) radius2     ! Inner radius (microns)
  read(fu, *) run_num     ! Identification number for run -> output file is    outputS_'run'.out 
  read(fu, *) rtflg1      ! External medium: 1=yes
  read(fu, *) rtflg2      ! Internal medium: 1=yes
  read(fu, *) omg1        ! External medium: single-scattering albedo
  read(fu, *) omg2        ! Internal medium: single-scattering albedo
  read(fu, *) lenrt1      ! External medium: mean free path (in micrometers)
  read(fu, *) lenrt2      ! Internal medium: mean free path (in micrometers)
  read(fu, *) np          ! Number of phase matrix angular points
  read(fu, *) nrn         ! Number of points in random number array
  read(fu, *) ncm         ! Number of G-L integration points for cosine map
  read(fu, *) sig         ! Relative standard deviation of radius.
  read(fu, *) nuc         ! Power law index for C_3 correlation.
  read(fu, *) lmax        ! Maximum degree in C_1, C_2, C_3.
  read(fu, *) ntr         ! Number of triangle rows in an octant.
  read(fu, *) fname       ! Prefix for the output file names
  close(fu)

  write(output_unit,'(A)') 'Geometric optics approximation for a Gaussian random sphere...'
  write(output_unit,'(A,F7.1,A)') 'Particle outer radius ', radius1, ' micrometers'
  write(output_unit,'(A,F7.1,A)') 'Particle inner radius ', radius2, ' micrometers'
  write(output_unit,'(A,I0)') 'Number of rays: ', nray

  call init_random(seed1)

  nsub=nray/npar
  if(nsub*npar /= nray) stop 'Trouble in GEO: ray/particle number mismatch.'

  !   Initialization of the Gaussian random sphere:
  beta=sqrt(log(sig**2+1.0_dp))
  call cs1cf(CSCF,nuc,2,lmax)  !power law
  call csini(CSCF,ell,cs2d,cs4d,lmax)

  call sgscfstd(SCFSTD,CSCF,beta,lmax)

  rho=beta/ell

  !   Refractive indices:
  m0 = cmplx(mmed, 0.0_dp, dp)    ! ref index of the medium
  m1 = cmplx(m1real, m1imag, dp) ! Mantle
  m2 = cmplx(m2real, m2imag, dp)  ! Core

  ! Size parameter, absorption parameter, refractive indices:
  rsc = radius2/radius1
  xa = 2.0_dp*pi*radius1/(wavelen/real(m0,dp))
  abscf1 = -2.0_dp*xa*aimag(m1/m0)
  abscf2 = -2.0_dp*xa*aimag(m2/m0)

  lenrt1 = lenrt1/radius1
  lenrt2 = lenrt2/radius1

  ! Discretization:
  call TRIDS(MUN,PHIN,IT,nnod,ntri,ntr)
  allocate(XN1(nnod,3), XN2(nnod,3), NT(ntri,3))
  
  ! Initialize:
  qsca = 0.0_dp
  qabs = 0.0_dp
  qin = 0.0_dp
  qout = 0.0_dp
  qis = 0.0_dp
  qbox = 0.0_dp
  S = 0.0_dp

  ! Initialize the scattering phase matrices:

  ! Mantle:
  dthe = pi/np
  allocate(CSRN1(0:nrn))
  if(rtflg1 == 1) then
    call pmatrix1_single(P1,np)
    call PSPLIVI(P1,CSRN1,XP1,YP1,YP21,pnorm1,dthe,nrn,np,ncm)
  end if

  ! Core:
  dthe = pi/np
  allocate(CSRN2(0:nrn))
  if(rtflg2 == 1) then
    call pmatrix1_single(P2,np,2)
    call PSPLIVI(P2,CSRN2,XP2,YP2,YP22,pnorm2,dthe,nrn,np,ncm)
   end if

  ! Storage grid parameters:
  bin = pi/nbin
  bin0 = cos(0.5_dp*bin)

  ! Refractive indices:
  call refindapp(MA0,m0,1.0_dp)

  ! Externally propagating rays. Initialization:
  jray = 0
  jsub = 0
  jbox = 0
  jboxin = 0
  wrayh = 0.0_dp ! rays that hit
  wraym = 0.0_dp ! rays that miss

  call sgscf(ACF,BCF,SCFSTD,lmax)
  call rgstd(XN1,NT,MUN,PHIN,ACF,BCF,rmax,beta,IT,nnod,ntri,lmax)

  do j0 = 1, nnod
    do j1 = 1, 3
      XN2(j0,j1) = rsc*XN1(j0,j1)
    end do
  end do

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
      write(fu,'(A,F15.8)') 'OuterRadius = ',radius1
      write(fu,'(A,F15.8)') 'InnerRadius = ',radius2
      write(fu,'(A,F15.8)') 'mmed = ',mmed
      write(fu,'(A,F15.8)') 'm1real = ',m1real
      write(fu,'(A,F15.8)') 'm1imag = ',m1imag
      write(fu,'(A,F15.8)') 'm2real = ',m2real
      write(fu,'(A,F15.8)') 'm2img = ',m2imag
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
    ! sample sphere:
    jsub=jsub+1
    if(jsub > nsub) then
      call sgscf(ACF,BCF,SCFSTD,lmax)
      call rgstd(XN1,NT,MUN,PHIN,ACF,BCF,rmax,beta,IT,nnod,ntri,lmax)
      do j0 = 1, nnod
        do j1 = 1, 3
          XN2(j0,j1) = rsc*XN1(j0,j1)
        end do
      end do
      jsub=1
    end if

    ! Weight of incident ray.
    wray = rmax**2

    ! Offset of incident ray:
    call random_number(ran2)
    r0 = rmax*sqrt(ran2)
    call random_number(ran2)
    phi0 = 2.0*pi*ran2
    X = (/ r0*cos(phi0), r0*sin(phi0), -sqrt(rmax**2-r0**2) /) 
    call raypart0(X,mui,nui,cphii,sphii) ! X is the location where the incident ray starts (particle coordinates)

    ! istri determines the possible interaction point on the particle surface.
    nie = 1.0_dp
    call istri(X,KEOUT,N1,XN1,NT,IT,nk1,len,ntri,nis,nie)
    if(nis == 0) then
      wraym = wraym+wray
      goto 100
    end if
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
      
    call raygetd1c(FBOX,KEBOX,KFBOX,HLBOX,HRBOX,XBOX,MABOX,PBOX,&
                   NISBOX,FOUT,KEOUT,KFOUT,HLOUT,HROUT,X,MAOUT,pout,nis,jbox)
    jbox = jbox - 1
      
    if(FOUT(1,1) <= wray*Fstop .or. pout > poutmax) then
      qout = qout+FOUT(1,1)
      goto 100
    end if
    N1 = NT(nis, 1:3)
    call prosca(nk1,N1,KEOUT)

  end if
  
  ! Scattering or external incidence:
  if(nk1 >= 0.0_dp) then ! jos kyseessa sisaltapain ulostulo
    call istri(X,KEOUT,N1,XN1,NT,IT,nk1,len,ntri,nis,1.0_dp)
    if(nis == 0) then
      call ehk(ELOUT,EROUT,HLOUT,KEOUT)
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
  pin = pout
  nis2 = nis
  call incide2l(FOUT,KEOUT,KFOUT,HLOUT,HROUT,MAOUT, &
              FIN,KEIN,KFIN,HLIN,HRIN,MAIN,&
              N1,m0,m1,totref)

  if(totref == 0) then

    jbox = jbox+1
    if(jbox > 499) then
      qbox=qbox+FOUT(1,1)
      goto 100
    end if
    call rayputd1c(FBOX,KEBOX,KFBOX,HLBOX,HRBOX,XBOX, &
                   MABOX,PBOX,NISBOX,FOUT,KEOUT,KFOUT, &
                   HLOUT,HROUT,X,MAOUT,pout,nis,jbox)
    jboxin = jboxin+1
    if(jboxin > 499) then
      qbox=qbox+FIN(1,1)
      goto 320
    end if
    call rayputd1c(FINBOX,KEINBOX,KFINBOX,HLINBOX,HRINBOX,XINBOX, &
                   MAINBOX,PINBOX,NISINBOX,FIN,KEIN,KFIN, &
                   HLIN,HRIN,X,MAIN,pin,nis,jboxin)
    goto 320

  else if(totref == 1) then

    jbox = jbox+1
    if(jbox > 499) then
      qbox=qbox+FOUT(1,1)
      goto 100
    endif
    call rayputd1c(FBOX,KEBOX,KFBOX,HLBOX,HRBOX,XBOX, &
                   MABOX,PBOX,NISBOX,FOUT,KEOUT,KFOUT, &
                   HLOUT,HROUT,X,MAOUT,pout,nis,jbox)
    goto 100

  else

    jboxin = jboxin+1
    if(jboxin > 499) then
      qbox=qbox+FOUT(1,1)
      goto 100
    end if
    call rayputd1c(FINBOX,KEINBOX,KFINBOX,HLINBOX,HRINBOX,XINBOX, &
                   MAINBOX,PINBOX,NISINBOX,FOUT,KEOUT,KFOUT, &
                   HLOUT,HROUT,X,MAOUT,pout,nis,jboxin)
    goto 320

  end if

  ! Internally propagating rays. Ray tracing in the internal medium:
320 if(jboxin == 0) then
    goto 100
  end if

  call raygetd1c(FINBOX,KEINBOX,KFINBOX,HLINBOX,HRINBOX,XINBOX,&
                 MAINBOX,PINBOX,NISINBOX,FIN,KEIN,KFIN,&
                 HLIN,HRIN,X,MAIN,pin,nis2,jboxin)
  jboxin=jboxin-1

  if(FIN(1,1) <= wray*Fstop .or. pin > pinmax) then
    qin=qin+FIN(1,1)
    goto 320
  end if

  ! ISTRI determines the next interaction point for the internal ray, 
  ! computes the inner normal at the interaction point, ABSORB determines the 
  ! amount of absorbed flux, and INCIDE determines the reflected and refracted 
  ! Mueller matrices, wave vectors, and parallel/perpendicular coordinate axes.
  Y1 = X
  Y2 = X

  nie = -1.0_dp
  nie2 = 1.0_dp
  call istri(Y1,KEIN,N1,XN1,NT,IT,nk1,len1,ntri,nis,nie)
  call istri(Y2,KEIN,N2,XN2,NT,IT,nk2,len2,ntri,nis2,nie2)
  if (nis == 0) then ! tahan ei kai pitais paatya mitaan? jos paatyy, niin kolmioinnissa vika. (tjs.)
    qis = qis+FIN(1,1)
    goto 320
  end if

  ! MANTLE INTERSECTION:

  ! Radiative transfer in the mantle:
  if (rtflg1 == 1) then
    call random_number(ran2)
    lenrn = -lenrt1*log(ran2)
    if (lenrn < min(len1,len2)) then
      pin = pin + 1
      pout = pout + 1

      call rtmc(FIN,KEIN,KFIN,HLIN,HRIN,X,XP1,YP1,YP21,CSRN1,&
        qabs,lenrn,abscf1,omg1,np,nrn)

      jboxin=jboxin+1
      if(jboxin > 499) then
        qbox=qbox+FIN(1,1)
        goto 320
      end if

      call rayputd1c(FINBOX,KEINBOX,KFINBOX,HLINBOX,HRINBOX,XINBOX, &
                     MAINBOX,PINBOX,NISINBOX,FIN,KEIN,KFIN, &
                     HLIN,HRIN,X,MAIN,pin,nis2,jboxin)

      goto 320
    end if
  end if

  ! Ray tracing in the mantle:
  if(len1 < len2) then

    X = Y1

    call prosca(kekf,KEIN,KFIN)
    lenabs=len1*abs(kekf)
    call absorb(FIN,qabs,lenabs,abscf1)

    if(FIN(1,1) <= wray*Fstop) then
      qin=qin+FIN(1,1)
      goto 320
    endif

    call incide2l(FIN,KEIN,KFIN,HLIN,HRIN,MAIN, &
                FOUT,KEOUT,KFOUT,HLOUT,HROUT,MAOUT, &
                N1,m1,m0,totref)

    pin=pin+1
    pout=pout+1

    if(totref == 0) then

      jbox=jbox+1
      if(jbox > 499) then
        qbox=qbox+FOUT(1,1)
        goto 320
      end if
      call rayputd1c(FBOX,KEBOX,KFBOX,HLBOX,HRBOX,XBOX, &
                     MABOX,PBOX,NISBOX,FOUT,KEOUT,KFOUT, &
                     HLOUT,HROUT,X,MAOUT,pout,nis,jbox)
      jboxin=jboxin+1
      if (jboxin > 499) then
        qbox=qbox+FIN(1,1)
        goto 320
      end if
      call rayputd1c(FINBOX,KEINBOX,KFINBOX,HLINBOX,HRINBOX,XINBOX, &
                     MAINBOX,PINBOX,NISINBOX,FIN,KEIN,KFIN, &
                     HLIN,HRIN,X,MAIN,pin,nis,jboxin)

    else if(totref == 1) then

      jboxin=jboxin+1
      if(jboxin > 499) then
        qbox=qbox+FIN(1,1)
        goto 320
      endif
      call rayputd1c(FINBOX,KEINBOX,KFINBOX,HLINBOX,HRINBOX,XINBOX, &
                     MAINBOX,PINBOX,NISINBOX,FIN,KEIN,KFIN, &
                     HLIN,HRIN,X,MAIN,pin,nis,jboxin)

    else

      jbox=jbox+1
      if(jbox > 499) then
        qbox=qbox+FIN(1,1)
        goto 320
      endif
      call rayputd1c(FBOX,KEBOX,KFBOX,HLBOX,HRBOX,XBOX, &
                     MABOX,PBOX,NISBOX,FIN,KEIN,KFIN, &
                     HLIN,HRIN,X,MAIN,pin,nis,jbox)

    end if

    goto 320

  else
    ! CORE INTERSECTION:

    ! Ray tracing in the mantle:
    X = Y2

    call prosca(kekf,KEIN,KFIN)
    lenabs=len2*abs(kekf)
    call absorb(FIN,qabs,lenabs,abscf1)

    if(FIN(1,1) <= wray*Fstop) then
      qin=qin+FIN(1,1)
      goto 320
    end if

    call incide2l(FIN,KEIN,KFIN,HLIN,HRIN,MAIN, &
                FCORE,KECORE,KFCORE,HLCORE,HRCORE,MACORE, &
                N2,m1,m2,totref)

    pin=pin+1
    pout=pout+1

    if(totref == 0) then

      jboxin=jboxin+1
      if(jboxin > 499) then
        qbox=qbox+FIN(1,1)
        goto 320
      end if
      call rayputd1c(FINBOX,KEINBOX,KFINBOX,HLINBOX,HRINBOX,XINBOX, &
                     MAINBOX,PINBOX,NISINBOX,FIN,KEIN,KFIN, &
                     HLIN,HRIN,X,MAIN,pin,nis2,jboxin)

    else if(totref == 1) then

      jboxin=jboxin+1
      if(jboxin > 499) then
        qbox=qbox+FIN(1,1)
        goto 320
      end if
      call rayputd1c(FINBOX,KEINBOX,KFINBOX,HLINBOX,HRINBOX,XINBOX, &
                     MAINBOX,PINBOX,NISINBOX,FIN,KEIN,KFIN, &
                     HLIN,HRIN,X,MAIN,pin,nis2,jboxin)
      goto 320

    else

      goto 410

    end if

    ! Ray tracing in the core:
410 if(FCORE(1,1) <= wray*Fstop .or. pin >= pinmax) then
      qin=qin+FCORE(1,1)
      goto 320
    end if

    Y2 = X
    call istri(Y2,KECORE,N2,XN2,NT,IT,nk2,len3,ntri,nis2,-1.0_dp)

    if(nis2 == 0) then
      qis=qis+FCORE(1,1)
      goto 320
    end if

    ! Radiative transfer in the core:
    if(rtflg2 == 1) then

      call random_number(ran2)
      lenrn=-lenrt2*log(ran2)
      if(lenrn < len3) then
        pin=pin+1
        pout=pout+1
        call rtmc(FCORE,KECORE,KFCORE,HLCORE,HRCORE, &
          X,XP2,YP2,YP22,CSRN2,qabs,lenrn,abscf2,omg2,np,nrn)
        goto 410
      end if
    end if

    ! Ray tracing in the core continued.
    X = Y2

    call prosca(kekf,KECORE,KFCORE)
    lenabs=len3*abs(kekf)
    call absorb(FCORE,qabs,lenabs,abscf2)

    if(FCORE(1,1) <= wray*Fstop) then
      qin=qin+FCORE(1,1)
      goto 320
    end if

    call incide2l(FCORE,KECORE,KFCORE,HLCORE,HRCORE,MACORE, &
                FIN,KEIN,KFIN,HLIN,HRIN,MAIN, &
                N2,m2,m1,totref)

    pin=pin+1
    pout=pout+1

    if(totref == 0) then

      jboxin=jboxin+1
      if(jboxin > 499) then
        qbox=qbox+FCORE(1,1)
        goto 320
      end if
      call rayputd1c(FINBOX,KEINBOX,KFINBOX,HLINBOX,HRINBOX,XINBOX, &
                     MAINBOX,PINBOX,NISINBOX,FIN,KEIN,KFIN, &
                     HLIN,HRIN,X,MAIN,pin,nis2,jboxin)
      goto 410

    else if(totref == 1) then

      goto 410

    else

      jboxin=jboxin+1
      if (jboxin > 499) then
        qbox=qbox+FCORE(1,1)
        goto 320
      end if
      call rayputd1c(FINBOX,KEINBOX,KFINBOX,HLINBOX,HRINBOX,XINBOX, &
                     MAINBOX,PINBOX,NISINBOX,FCORE,KECORE,KFCORE, &
                     HLCORE,HRCORE,X,MAIN,pin,nis2,jboxin)
      goto 320

    end if
  end if


END PROGRAM singletwolayer
