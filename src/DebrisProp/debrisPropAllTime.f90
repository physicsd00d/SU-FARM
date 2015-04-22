subroutine debrisPropagation(finalConditions,numFinalSteps, initialState,DebrisVel,&
     mass,Sref,&
     nCD,MinfCD,CD,& !Drag inputs
     CLOption,nCL,MinfCL,CL,& !Lift inputs
     LoverD,& !LoverD inputs
     atmosOption,altitudeList, densityList,UList,VList,WList,& !! atmospheric parameter inputs
     GEOptions,filename,nList,PlanetModel,dtInterval,ndtInterval,thetag0)


  use planetEarthCharacteristics
  use constants
  !use newpos_e10_I
  implicit none

  ! Debris propagation routine
  ! Created by Francisco Capristan for the FAA COE-CST effort
  ! July 2011

  ! Modified by Francisco Capristan
  ! June 2012
  ! Notes : Added CD and CL as a function of Minf

  !#### Inputs#################################################################
  ! InitialState     : [x,y,z,Vx,Vy,Vz] : inputs in ECI coordinate frame SI UNITS

  ! DebrisVel          : [Vx_debris,Vy_debris,Vz_debris] : impulse debris velocities ECI coordinate frame
  ! mass                : debris piece mass [kg]
  ! Sref                : debris Sref [m2]
  ! nCD                 : length of CD array. if nCd = 1 then Constant Cd
  ! MinfCD              : array containing Mach numbers for CD interpolation. Must be length nCD 
  ! CD                  : debris piece drag coefficient function of Minf. Must be length nCD
  ! CLOption            : if 1 then CL is used. if 0 then constant LoverD used instead
  ! nCL                 : length of CL array. if nCL = 1, then constant CL
  ! MinfCL              : array containing Mach numbers for CL interpolation. Must be length nCL 
  ! CL                  : debris piece lift coefficient as a function of Minf.
  ! LoverD              : used if ClOption =0. Specifies constant L over D
  ! atmosOption         : if 0 Exponential density used NO WIND. if 1 given density profile used NO WIND. If 2 given density and wind profile used. -1 Vacuum
  ! altitudeList        : altitude vector in decreasing order (to be used in atmospheric profile)
  ! densityList         : density profile corresponding to altitudeList
  ! U/V/W  list         : wind profile corresponding to altitudeList. Usually obtained from GRAM
  ! GEOptions           : if 1, then the trajectory will be writen to filename 
  ! filename            : filename to write trajectory solution suitable for google earth visualization. Only used if GEOptions =1
  ! nList               : number of data points in atmospheric profile    
  ! PlanetModel         : 0 => spherical , 1=> oblate. J2 effects ignored (not calculating orbits, Drag and Lift are the dominating forces)         
  ! dtInterval          : approx value to do rk45 stepping...will be used to set upper limit for dt
  ! ndtInterval         : size of desired returning vector...must be ndtinterval*dt > TOF

  !### Outputs:

  ! FinalConditions :
  !                     finalConditions = (/altitudeFinal,latitudeFinal,longitudeFinal,time,VrelMag,flag/) 

  !                      altitudeFinal       : Final altitude [m], if all successful, then it must be <0
  !                      latitudeFinal       : final latitude [degrees]
  !                      longitudeFinal      : final longitude [degrees]
  !                      VrelMag             : velocity (relative to rotating planet) magnitude
  !                      flag                : returns potential issues (0 -> OK RUN , 1 -> debris landing, but time step appears too big,
  !                                            2 -> debris appears to be in orbit, or more iterations are needed...If not suppose to be in orbit =>try increasing number of max RK4 iterations...or decreasing dt
  !                                            3-> debris is in orbit
  !#############################################################################


  ! Input Parameters
  !##############################################
  double precision , dimension(3) , intent(in) :: DebrisVel
  double precision , dimension(6) , intent(in) :: initialState
  double precision , intent(in) :: mass,Sref, dtInterval
  character*16 ,intent(in) :: filename
  ! Aerodynamic inputs
  integer , intent(in) :: nCd,nCL
  double precision , dimension(nCD),intent(in) :: MinfCD,CD
  double precision , dimension(nCL),intent(in) :: MinfCL,CL
  double precision , intent(in) :: LoverD,thetag0



  ! Input Options
  integer , intent(in) :: GEOptions,CLOption,atmosOption,nList,PlanetModel
  double precision , dimension(nList),intent(in)::altitudeList,densityList,UList,VList,WList
  !##############################################

  !Output Parameters
  !##############################################
  double precision :: latitudeFinal, longitudeFinal, altitudeFinal, VrelMag !TJC, I'm actually outputting VrelativeVector
  double precision :: flag
  integer, parameter :: numConds = 7               !TJC, The number of elements being returned in finalConditions
  integer, intent(in) :: ndtInterval
  double precision , dimension(ndtInterval*numConds),intent(out) :: finalConditions
  integer, intent(out) :: numFinalSteps     !TJC Addition

  !##############################################



  !##############################################

  !For Calculations
  !##############################################
  double precision :: Vrelative , flightPathAngle, HeadingAngle!, angRate
  double precision :: Vmain,mainFlightPathAngle,mainHeadingAngle,altitude, longitude,latitude

  double precision :: maxLat,minLat,maxLon,minLon,MainLat,MainLon, HeightRef,dtmax,dtmin

  double precision, dimension(3) :: V,r,Vrk4,r_np1,Vrk4_np1, rlocal, Vc, VmainC !r and V in cartesian coordinates centered at the center of the earth
  double precision :: t,distCenter, normr, oldAltitude,localR,latLocal,dt,tref,tref_np1! for calculations
  double precision , dimension(3,3) :: localRotation , Rzgamma, Rxbeta, Rwindx, Rwindz,MainLocalRotation    !TJC says only localRotation is in this file
  double precision , dimension(3,3) :: RzMaingamma, RxMainbeta                                      !TJC says only none is in this file
  double precision , dimension(7) :: state
  double precision , dimension(3) :: VrelativeVector, VrelativeInInertialFrame, Vrotfinal           !TJC says only Vrotfinal is in this file, i'm going to use VrelativeVector
  double precision , dimension(3) :: VmainRelativeVector, VmainInInertialFrame, debRange, Rorig     !TJC says only Rorig is in this file
  double precision , dimension(3) :: latlonalt, Rmain
  integer :: totalCounter
  !##############################################



  ! RK suite parameters
  integer ( kind = 4 ), parameter :: neqn = 6

  double precision abserr
  !external f01
  integer  index1
  integer ( kind = 4 ) iflag
  integer ( kind = 4 ) iwork(5)
  !real ( kind = 8 ) :: pi = 3.141592653589793D+00
  double precision relerr
  integer ( kind = 4 ), parameter :: step_num = 12
  !  double precision t
  double precision tout
  double precision work(100+21*neqn)
  double precision y(neqn)

  integer printIX

  finalConditions = -9999 ! initializing to all negative values



  maxLat = -99999.9 ! initialing min and max vals
  maxLon = -99999.9
  minLon =  99999.9
  minLat =  99999.9






  t=0 ! initializing time 
  tref = 0!
  r(1) = initialState(1)
  r(2) = initialState(2)
  r(3) = initialState(3)

  V(1) = initialState(4) + DebrisVel(1) ! debris velocity due to explosion added as an impulse
  V(2) = initialState(5) + DebrisVel(2)
  V(3) = initialState(6) + DebrisVel(3)
  Vrk4 = V
  Y(1) = r(1)
  Y(2) = r(2)
  Y(3) = r(3)
  Y(4) = V(1)
  Y(5) = V(2)
  Y(6) = V(3)
  tout = 0.0
  latlonalt = latLonAltCalculation(y(1:3),tout,thetag0,PlanetModel)

  latitude = latlonalt(1)
  longitude = latlonalt(2)
  altitudeFinal = latlonalt(3)


  if(GEOptions.EQ.1)then 

     open(unit=7,file=filename,status='unknown')
     write(7,*) 'Longitude(deg) Latitude(deg) Altitude(m) Time(sec)'
     write(7,'(4F25.13)')  longitude,latitude,altitudeFinal,tref

  endif

  abserr = 0.0001D+00
  relerr = 0.0001D+00

  iflag = 1
  index1 = 1

  Vrotfinal(1) = -Y(2)*omega
  Vrotfinal(2) = Y(1)*omega
  Vrotfinal(3) = 0
  !VrelMag = norm(Y(4:6)-Vrotfinal)

  VrelativeVector = Y(4:6)-Vrotfinal    !TJC
  VrelMag = norm(VrelativeVector)       !TJC

  !write(*,*) 'index1 = ', index1


  finalConditions(index1:numConds) = (/latlonalt(1),latlonalt(2),latlonalt(3),&
    VrelativeVector(1),VrelativeVector(2),VrelativeVector(3),flag/)
!  write(*,*) 'ndtInterval = ', ndtInterval, '   dtInterval = ', dtInterval, '   tout = ', tout

  do while ((latlonalt(3)>=0.0).AND.(tout <=5.0*3600).AND.(index1 < ndtInterval)) ! stop when it reaches the ground, or 5 hours and still not on the ground
     tout = tout + dtinterval

     !write(*,*) 'ndtInterval = ', ndtInterval, '   dtInterval = ', dtInterval, '   tout = ', tout

     call ode (F, neqn, y, t, tout, relerr, abserr, iflag, work, iwork )



     ! TJC: I think this is crashing my parallel python stuff.  Comment it out because I don't care if the piece hasn't landed
     if ( iflag /= 2 ) then
        !write ( *, '(a)' ) ' '
        !write ( *, '(a)' ) 'TEST01 - Fatal error!'
        !write ( *, '(a,i8)' ) '  ODE returned IFLAG = ', iflag
        !exit
     end if
     Vrotfinal(1) = -y(2)*omega
     Vrotfinal(2) = y(1)*omega
     Vrotfinal(3) = 0
     !VrelMag = norm(Y(4:6)-Vrotfinal)

     VrelativeVector = Y(4:6)-Vrotfinal    !TJC
     VrelMag = norm(VrelativeVector)       !TJC


     latlonalt = latLonAltCalculation(y(1:3),tout,thetag0,PlanetModel)
     finalConditions( (numConds*index1 +1):(numConds*index1 +numConds) ) &
        = (/latlonalt(1),latlonalt(2),latlonalt(3),&
        VrelativeVector(1),VrelativeVector(2),VrelativeVector(3),flag/) !TJC changed order
     index1 = index1 + 1

     if(GEOptions.EQ.1)then
        open(unit=7,file=filename,status='unknown')
        write(7,'(4F25.13)') latlonalt(2),latlonalt(1),latlonalt(3),tout
     endif

  end do

  !write(*,*) 'index1 = ', index1
  numFinalSteps = index1

  printIX = 1
  do while (printIX < numConds*index1)
    !write(*, *) 'fC[', printIX, '] = ', finalConditions(printIX)
    printIX = printIX +1
  end do


  flag = 0

  if (latlonalt(3)>0)then
     flag = 2
  end if

  if (ndtInterval*dtinterval<=tout) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'ERROR in debrisPropagation. Increase ndtInterval'
        write ( *, * ) ' current ndtInterval = ', ndtInterval
        write ( *, * ) ' try ndtInterval >= ', tout/dtinterval
        stop
     end if



  !finalConditions = (/altitudeFinal,latitudeFinal,longitudeFinal,VrelMag,flag/) 
  !finalConditions = (/latlonalt(3),latlonalt(1),latlonalt(2),VrelMag,flag/) 

  close(unit=7)


  !*********

contains 


  subroutine F(ccT,Y,YP)
    implicit none
    double precision :: ccT
    double precision :: Y(6),YP(6)
    double precision , dimension(3) :: resdmdvdt

    resdmdvdt = dvdt(Y(1:3),Y(4:6),ccT)
    YP(1) = Y(4)
    YP(2) = Y(5)
    YP(3) = Y(6)
    YP(4) = resdmdvdt(1)
    YP(5) = resdmdvdt(2)
    YP(6) = resdmdvdt(3)

    return
  END subroutine F



  function dvdt(r,V,ct)
    use constants
    !use GRAMAtmos
    double precision , intent(in) :: ct
    double precision , dimension(3) , intent(in) :: r, V
    double precision , dimension(3) :: weight, drag, dvdt, Vinf, Vrot, interVec, liftDir,lift, r2, LiftRot, latlonalt,dragDir
    double precision , dimension(3) :: VwindLocal,Vwind
    double precision :: altitude , rho, g, newden,newU,newV,newW!,rotAngle
    double precision :: clat,clon,Minf,speedOfSound,CDLocal,CLLocal,VinfMag,tempdebug
    double precision , dimension(3,3) :: localRotation
    ! print *, 'Solving dvdt'   



    latlonalt = latLonAltCalculation(r,ct,thetag0,PlanetModel)!
    clat = latlonalt(1)
    clon = latlonalt(2)
    altitude = latlonalt(3)


    !rotAngle = ct*angRate
    !altitude = norm(r) - planetRadius
    Vrot(1) = -r(2)*omega
    Vrot(2) = r(1)*omega
    Vrot(3) = 0
    !rho = density(altitude)
    Vwind(1)=0
    Vwind(2)=0
    Vwind(3)=0


    if((atmosOption==2).OR.(atmosOption==1)) then

       call atmosInterpolation(newden,newU,newV,newW,altitude,altitudeList,densityList,UList,VList,WList,nList)


       if (altitude>altitudeList(1)) then
          if (altitude<200000) then !above 150 km set density to zero
             rho = (densityList(1)/density(altitudeList(1)))*density(altitude)
             !rho = 0
             !altitude = altitudeList(1)
             !   print *, 'Using exponential model for density'
          else
             rho = 0.0
          end if
       else
          rho = newden !using density given by user
       end if
       !      rho = density(altitude)

       if (altitude>200000.) then
          rho = 0.0
       end if




       if (atmosOption==2) then

          VwindLocal = (/newW,newU,newV/)
          localRotation = transpose(matmul(rRotation(2,-clat*Pi/180),rRotation(3,thetag0+omega*ct+clon*Pi/180)))
          Vwind = matmul(localRotation,VwindLocal)

       end if
    elseif (atmosOption==0) then

       rho = density(altitude)

    elseif (atmosOption ==-1) then
       rho = 0

    else
       print *, 'Error in debrisPropagation: Unknown atmosOption Flag'
       STOP

    end if



    !print *,rho-newden
    g = gravity(altitude)
    weight = mass*(-g*r/norm(r))
    Vinf = V- Vwind-Vrot
    VinfMag = norm(Vinf)
    call getspeedofSound(altitude,speedOfSound)
    Minf = VinfMag/speedOfSound
    !print *,Minf

    if (atmosOption.NE.-1) then

       ! Drag Calculation
       if (nCD>1) then !using CD as a function of Minf

          call LinearCDInterpolation(MinfCD,CD,nCD,Minf,CDLocal)
          ! print *,CDlocal
       else
          CDlocal = CD(1)
       end if

       dragdir = -Vinf !not a unit vector
       drag = 0.5*CDLocal*rho*Sref*VinfMag*dragdir


       !  drag = drag*(1-0.3*sin(rotAngle))

       !! Lift calculation
       !! getting lift direction
       interVec = cross(r,dragdir)
       if (norm(interVec)> 0.00001) then
          interVec = intervec/norm(intervec)
          liftDir = cross(dragdir,interVec)
          liftDir = liftDir/norm(liftDir)
       else
          r2(1) = r(1) - Rorig(1)
          r2(2) = r(2) - Rorig(2)
          r2(3) = r(3) - Rorig(3)
          interVec = cross(r2,dragdir)
          liftDir = cross(dragdir,interVec)
          liftDir = liftDir/norm(liftDir)
          !print *, 'Switching to offset r for LIFT \n'

       end if

       if (CLOption==1)then

          if (nCL>1)then
             call LinearCLInterpolation(MinfCL,CL,nCL,Minf,CLLocal)



          else
             CLLocal = CL(1)
          end if
          lift = 0.5D0*CLLocal*rho*Sref*VinfMag*VinfMag*liftDir
       else ! using constant L over D
          !lift = (LoverD*norm(drag))*liftDir
          lift = (LoverD*CDLocal*0.5D0*rho*VinfMag**2*Sref)*liftDir
          !lift = -norm(drag)*LoverD*liftDir
       end if

    else
       drag = (/0,0,0/)
       lift = (/0,0,0/)
    end if


    ! liftRot = lift*cos(rotAngle) - interVec*sin(rotAngle)*norm(lift) ! this line if a rotating lift vector is desired

    !liftRot = sin(rotAngle)*lift ! this line if a sinusoidal lift vector wrt Time is desired. Set rotAngle to zero for regular lift
    !dvdt = (drag + liftRot + weight)/mass

    dvdt = (drag + lift + weight)/mass
    ! print *,norm(weight)
    !print *,Cdlocal,speedofSound,Minf

  end function dvdt

  function density(altitude)
    use constants
    implicit none

    double precision, intent(in) :: altitude
    double precision :: density

    !  g = gravity(altitude)
    density = rho0*exp(-altitude/Hscale)


  end function density


  function gravity(altitude)
    use constants
    implicit none

    double precision , intent(in) :: altitude
    double precision :: gravity

    gravity = muPlanet/((altitude+planetRadius)**2.D0)
    ! print *, gravity
  end function gravity


  function norm(x)
    use constants
    double precision, dimension(3), intent(in) :: x
    double precision :: norm

    norm = sqrt(x(1)**2.D0 + x(2)**2.D0 + x(3)**2.D0)

  end function norm

  function rk45(r,V,t,dtmax,dtmin)
    use constants
    implicit none

    integer :: counter
    double precision , intent(in) :: t
    double precision , dimension(3), intent(in) :: r,V
    double precision , dimension(7) :: rk45
    double precision ,dimension(3) :: resrvec,resVvec
    double precision :: resr,resV, residual,tol,delta
    double precision , dimension(3) :: r1, r2, r3 ,r4,r5,r6, v1, v2, v3, v4,v5,v6, a1, a2, a3, a4,a5,a6, rf, vf
    double precision , intent(in):: dtmax ,dtmin
    !double precision ,save:: dt = dtmax
    double precision :: dt
    tol = 2.0D-5
    !dtmin = 1.D-5
    dt = dtmax

    !dt = dtmax
    counter = 0

    do while(dtmax>dtmin)
       counter = counter +1
       r1 = r
       v1 = V
       a1 = dvdt(r1,V1,t)

       r2 = r1 + k2_1*dt*v1
       v2 = v1 + k2_1*dt*a1
       a2 = dvdt(r2,v2,t+k2_1*dt)


       r3 = r1 + dt*(k3_2*v1 + k3_3*v2)
       v3 = v1 + dt*(k3_2*a1 + k3_3*a2)
       a3 = dvdt(r3,v3,t + k3_1*dt)

       r4 = r1 + dt*(k4_2*v1 + k4_3*v2 + k4_4*v3)
       v4 = v1 + dt*(k4_2*a1 + k4_3*a2 + k4_4*a3)
       a4 = dvdt(r4,v4,t + k4_1*dt)

       r5 = r1 + dt*(k5_1*v1 + k5_2*v2 + k5_3*v3 + k5_4*v4)
       v5 = v1 + dt*(k5_1*a1 + k5_2*a2 + k5_3*a3 + k5_4*a4)
       a5 = dvdt(r5,v5,t+dt)

       r6 = r1 + dt*(k6_2*v1 + k6_3*v2 + k6_4*v3 + k6_5*v4 + k6_6*v5)
       v6 = v1 + dt*(k6_2*a1 + k6_3*a2 + k6_4*a3 + k6_5*a4 + k6_6*a5)
       a6 = dvdt(r6,v6,t + k6_1*dt)

       !resrvec = 1/360*r1 - 128/4275*r3 - 2197/75240*r4 + .02*r5 + 2/55*r6
       ! resr = norm(resrvec)/dt

       resVvec = res_1*v1 + res_2*v3 + res_3*v4 + res_4*v5 + res_5*v6
       resV = norm(resVvec)

       residual = resV




       ! print *, residual,dt
       if (residual<=tol) then
          rf = r1 + dt*(ap_1*v1 + ap_2*v3 + ap_3*v4 + ap_4*v5)
          vf = v1 + dt*(ap_1*a1 + ap_2*a3 + ap_3*a4 + ap_4*a5)
          !rf = r1 + dt*(16.0/135.0*v1 + 6656.0/12825.0*v3 + 28561.0/56430.0*v4 - 9.0/50.0*v5 + 2.0/55.0*v6)
          !vf = v1 + dt*(16.0/135.0*a1 + 6656.0/12825.0*a3 + 28561.0/56430.0*a4 - 9.0/50.0*a5 + 2.0/55.0*a6)

          ! print *,' exiting loop'
          !print *, norm(rf),residual
          exit

       end if
       delta = 0.84D0*(tol/residual)**(0.25D0)
       if (delta <= 0.1) then
          dt = .1D0*dt
       elseif (delta>=4.0) then
          dt = 4.*dt
       else
          dt = delta*dt
       end if



       if (dt >dtmax) then
          dt = dtmax
       else if (dt<=dtmin) then
          !  print *, 'Min dt reached'
          ! print *, residual
       end if

       if ((residual+1)==residual) then
          dt = .5*dt
       elseif (dt< dtmin) then
          ! print *,norm(rf),dt,norm(r1)

          !dt = dtmin
          print *, 'ERROR in propagator/n Check adaptive time step'
          STOP

       end if

    end do

    !rf = r1 + (dt/6)*(v1+ 2*v2 +2*v3 +v4)
    !vf = v1 + (dt/6)*(a1+ 2*a2 +2*a3 +a4)

    rk45 = (/rf(1),rf(2),rf(3),vf(1),vf(2),vf(3),dt/)
    ! print *,norm(rf),dt
    !rk4(2,:)=vf
    !print *, 'Done with one step'
    !if (dt>.00001) then
    !print *,norm(rf),dt,norm(r1),residual
    !end if


  end function rk45

  function latlonaltCalculation(r,t,thetag0,model)
    use constants
    use planetEarthCharacteristics
    implicit none
    ! calculates spherical lat lon 
    ! uses ellipse for altitude if model = 1

    double precision , dimension(3) , intent(in) :: r
    integer , intent(in) :: model
    double precision , intent(in) :: t,thetag0
    double precision :: normr, rdelta,phigd,Nphi,phigd_np1,tol,height
    double precision , dimension(3):: rlocal, latlonaltCalculation,RR

    tol = 999.
    rlocal =  matmul(Rrotation(3,thetag0+omega*t),r) ! accounting for Rotating planet
    normr = norm(rlocal)


    latlonaltCalculation(2) = atan2(rlocal(2),rlocal(1))*180.0D0/Pi  !Longitude in degrees

    phigd = asin(rlocal(3)/normr) !Latitude in radians

    if (model==1)then
       rdelta = sqrt(rlocal(1)**2+rlocal(2)**2)
       do while (tol>1e-7)
          Nphi = Requator/sqrt(1-ecc**2*(sin(phigd))**2)
          phigd_np1 = atan((rlocal(3) + Nphi*ecc**2*sin(phigd))/rdelta)
          tol = abs(phigd - phigd_np1)
          phigd = phigd_np1
       end do
       height = rdelta/cos(phigd) - Nphi
    elseif (model==0)then
       height = normr - planetRadius
    end if
    latlonaltCalculation(1) = phigd*180.0/pi ! latitude in degrees
    latlonaltCalculation(3) = height

  end function latlonaltCalculation










  function Rrotation(axis,angle) ! rotation from N frame to axis desired cRn
    use constants
    implicit none
    double precision , intent(in) :: angle
    integer , intent(in) :: axis
    double precision , dimension(3,3) :: Rrotation

    if (axis==1) then
       Rrotation(1,1) = 1
       Rrotation(1,2) = 0
       Rrotation(1,3) = 0
       Rrotation(2,1) = 0
       Rrotation(2,2) = cos(angle)
       Rrotation(2,3) = sin(angle)
       Rrotation(3,1) = 0
       Rrotation(3,2) = -sin(angle)
       Rrotation(3,3) = cos(angle)
    elseif (axis ==2) then
       Rrotation(1,1) = cos(angle)
       Rrotation(1,2) = 0
       Rrotation(1,3) = -sin(angle)
       Rrotation(2,1) = 0
       Rrotation(2,2) = 1
       Rrotation(2,3) = 0
       Rrotation(3,1) = sin(angle)
       Rrotation(3,2) = 0
       Rrotation(3,3) = cos(angle)
    elseif (axis ==3) then
       Rrotation(1,1) = cos(angle)
       Rrotation(1,2) = sin(angle)
       Rrotation(1,3) = 0
       Rrotation(2,1) = -sin(angle)
       Rrotation(2,2) = cos(angle)
       Rrotation(2,3) = 0
       Rrotation(3,1) = 0
       Rrotation(3,2) = 0
       Rrotation(3,3) = 1
    end if
  end function Rrotation

  function cross(a,b) 
    use constants
    implicit none

    double precision , dimension (3), intent(in) :: a,b
    double precision, dimension (3)::cross

    cross(1) = a(2) * b(3) - a(3) * b(2)
    cross(2) = a(3) * b(1) - a(1) * b(3)
    cross(3) = a(1) * b(2) - a(2) * b(1)
  end function cross





end subroutine debrisPropagation


