
subroutine propagate(finalConditions,finalDerivs,initialState,&
     mass,mass_time_alt_final,Sref,&
     nCD,MinfCD,CD,& !Drag inputs
     CLOption,nCL,MinfCL,CL,& !Lift inputs
     LoverD,& !LoverD inputs
     atmosOption,altitudeList, densityList,UList,VList,WList,& !! atmospheric parameter inputs
     Tlist,timeList,ntime,Isp,&
     GEOptions,filename,nList,PlanetModel,dtInterval,ndtInterval,thetag,mass_time_alt_opt,thrustOffAngleDeg)




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
  ! InitialState     : [x,y,z,Vx,Vy,Vz] : inputs in ECEF coordinate frame SI UNITS

  ! mass                : initial mass [kg]
  ! mass_time_alt_final : mass, time or altitude final (to stop RK propagation) depending on mass_time_alt_opt
  ! alt_final           : final altitude to propagate[m]
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
  ! Tlist               : Thrust profile as a function of time [N]
  ! timeList            : list with time values that corresponds to Tlist
  ! ntime               : length(timeList)
  ! Isp                 ; Isp value to be used (constant) [sec]
  ! GEOptions           : if 1, then the trajectory will be writen to filename 
  ! filename            : filename to write trajectory solution suitable for google earth visualization. Only used if GEOptions =1
  ! nList               : number of data points in atmospheric profile    
  ! PlanetModel         : 0 => spherical , 1=> oblate. J2 effects ignored (not calculating orbits, Drag and Lift are the dominating forces)         
  ! dtInterval          : approx value to do rk45 stepping...will be used to set upper limit for dt
  ! thetag              : angle in radians for conversion between ECI to ECEF
  ! mass_time_alt_opt   : 1 for mass option, 2 for time, 3 for altitude
  ! thrustOffAngleDeg   : a list with 2 elements that specify an small offset to the thrust vector (degrees)
   !### Outputs:

  ! FinalConditions :
  !                     finalConditions = (x,y,z,Vx,Vy,Vz,mass,t) statevectors + mass + time at different steps
  !#############################################################################



!    Input/output, integer ( kind = 4 ) IFLAG, indicates the status of integration.
!    On input, IFLAG is normally 1 (or -1 in the special case where TOUT is
!    not to be exceeded.)  On normal output, IFLAG is 2.  Other output values
!    are:
!    * 3, integration did not reach TOUT because the error tolerances
!      were too small.  But RELERR and ABSERR were increased appropriately
!      for continuing;
!    * 4, integration did not reach TOUT because more than 500 steps were taken;
!    * 5, integration did not reach TOUT because the equations appear to
!      be stiff;
!    * 6, invalid input parameters (fatal error).
!    The value of IFLAG is returned negative when the input value is



  !external F
  ! Input Parameters
  !##############################################
  double precision , dimension(6) , intent(in) :: initialState
  double precision , intent(in) :: mass,Sref, dtInterval,thetag,mass_time_alt_final
  character*16 ,intent(in) :: filename
  double precision , dimension(2) , intent(in) :: thrustOffAngleDeg
  ! Aerodynamic inputs
  integer , intent(in) :: nCd,nCL,ntime,mass_time_alt_opt
  double precision , dimension(nCD),intent(in) :: MinfCD,CD
  double precision , dimension(nCL),intent(in) :: MinfCL,CL
  double precision , intent(in) :: LoverD,Isp
  double precision , dimension(ntime,3) , intent(in) :: Tlist
  double precision , dimension(ntime) , intent(in) :: timeList



  ! Input Options
  integer , intent(in) :: GEOptions,CLOption,atmosOption,nList,PlanetModel
  double precision , dimension(nList),intent(in)::altitudeList,densityList,UList,VList,WList
  !##############################################

  !Output Parameters
  !##############################################
  double precision ::   altitudeFinal
  !double precision :: flag
  !double precision , allocatable,dimension(:,:) :: tempFinalConditions
  !double precision , dimension(1,1),intent(out) :: finalConditions
  !##############################################
  !##############################################

  !For Calculations
  !##############################################
  double precision ::  longitude,latitude
  double precision, dimension(3) :: V,r,Vrk4 !r and V in cartesian coordinates centered at the center of the earth
  double precision :: t, tref! for calculations
  !double precision , dimension(3,3) ::   Rwindx
  !double precision , dimension(3,3) ::  RxMainbeta
  !double precision , dimension(8) :: state
  double precision , dimension(3) ::   Rorig
  double precision , dimension(3) :: latlonalt
  double precision , dimension(4) :: resdmdvdt
  double precision , dimension(2) :: thrustOffAngleRad
  ! integer :: totalCounter

  !##############################################

  ! RK suite parameters
  integer ( kind = 4 ), parameter :: neqn = 8

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
  integer, intent(in) :: ndtInterval
  !nfinal =  = ceiling(2*timelist(ntime)/dtinterval)
  double precision , dimension(ndtInterval,12),intent(out) :: finalConditions
  double precision , dimension(6),intent(out) :: finalDerivs
  finalConditions = -9999



  thrustOffAngleRad = thrustOffAngleDeg*Pi/180.0


  t=0 ! initializing time 
  tref = 0!
  r(1) = initialState(1)
  r(2) = initialState(2)
  r(3) = initialState(3)

  V(1) = initialState(4) 
  V(2) = initialState(5) 
  V(3) = initialState(6)
  Vrk4 = V
  latlonalt = latLonAltCalculation(r,t,thetag,PlanetModel)

  latitude = latlonalt(1)
  longitude = latlonalt(2)
  altitudeFinal = latlonalt(3)


  y(1) = r(1)
  y(2) = r(2)
  y(3) = r(3)
  y(4) = mass
  y(5) = V(1)
  y(6) = V(2)
  y(7) = V(3)
  y(8) = -norm(tlist(1,:))/(9.8*Isp)


  !abserr = 0.1D+00
  !relerr = 0.1D+00

  !print *,'OISHDFOSDIHFSODIFHSDOFIHSDOFIHDSOFISHDFOIHSDOFIHSDOFIHSDOFISHDFOISDHFOISHDFOSDIHFOISDHFOSIDHFOSIHFODIHF'
  abserr = 0.00001D+00
  relerr = 0.00001D+00

  iflag = 1
  tout = 0.0
  !  101  format(F,F,F,F,F,F)
  !101  format('(5F20.13)')

  if(GEOptions.EQ.1)then 
     open(unit=7,file=filename,status='unknown')

     write(7,*) 'Longitude(deg) Latitude(deg) Altitude(m) Time(sec)'
     latlonalt = latLonAltCalculation(y(1:3),tout,thetag,PlanetModel)
     write(7,'(5F25.13)') latlonalt(2),latlonalt(1),latlonalt(3),tout,y(4)
  endif

  index1 = 1
  finalConditions(index1,:) = (/y(1),y(2),y(3),y(5),y(6),y(7),latlonalt(2),latlonalt(1),latlonalt(3),y(4),tout,thetag + tout*omega/)
  !                            3 positions, 3 velocities,lat,lon,alt, mass, time

  do while (.TRUE.) ! (latlonalt(3)<=alt_final)

     if (mass_time_alt_opt==1) then
        if(mass_time_alt_final>=y(4)) then
           exit
        end if
     elseif (mass_time_alt_opt==2) then
        if(mass_time_alt_final<=tout) then
           exit
        end if
     elseif(mass_time_alt_opt==3) then  
        if(mass_time_alt_final<=latlonalt(3)) then  !Go above a threshold (good for launch)
           exit
        end if
     elseif(mass_time_alt_opt==4) then
        if(mass_time_alt_final>=latlonalt(3)) then  !Drop below a threshold (good for reentry)
           exit
        end if
     else
        print *,'mass_time_alt_opt must be an integer between 1-3. ERROR in trajPop.f90'
        STOP
     end if
     if (latlonalt(3)<=0) then
        exit
     end if

     tout = tout + dtinterval
     index1 = index1 + 1
     call ode (F, neqn, y, t, tout, relerr, abserr, iflag, work, iwork )

     if ( iflag /= 2 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TEST01 - Fatal error!'
        write ( *, '(a,i8)' ) '  ODE returned IFLAG = ', iflag
        exit
     end if
     latlonalt = latLonAltCalculation(y(1:3),tout,thetag,PlanetModel)

     finalConditions(index1,:) =  (/y(1),y(2),y(3),y(5),y(6),y(7),latlonalt(2), &
     latlonalt(1),latlonalt(3),y(4),tout,thetag + tout*omega/)
     ! print *,'alt',latlonalt(3)
     !print *,'time',tout
     if(GEOptions.EQ.1)then 
        open(unit=7,file=filename,status='unknown')
        write(7,'(5F25.13)') latlonalt(2),latlonalt(1),latlonalt(3),tout,y(4)
     endif

     ! now let's store the final derivatives
     
     finalDerivs(1) = y(5)
     finalDerivs(2) = y(6)
     finalDerivs(3) = y(7)
     resdmdvdt = dmdvdt(Y(1:3),Y(5:8),Y(4),tout)
     finalDerivs(4) = resdmdvdt(2)
     finalDerivs(5) = resdmdvdt(3)
     finalDerivs(6) = resdmdvdt(4)
     !dmdvdt = (/dm,dvdt(1),dvdt(2),dvdt(3)/)


  end do








  close(unit=7)


  !*********

contains 


  subroutine F(ccT,Y,YP)
    implicit none
    double precision :: ccT
    double precision :: Y(8),YP(8)
    double precision , dimension(4) :: resdmdvdt

    resdmdvdt = dmdvdt(Y(1:3),Y(5:8),Y(4),ccT)
    YP(1) = Y(5)
    YP(2) = Y(6)
    YP(3) = Y(7)
    YP(4) = resdmdvdt(1)
    YP(5) = resdmdvdt(2)
    YP(6) = resdmdvdt(3)
    YP(7) = resdmdvdt(4)
    YP(8) = 0.0

    return
  END subroutine F

  function dmdvdt(r,V,mass_c,ct)
    use constants
    !use GRAMAtmos
    double precision , intent(in) :: ct,mass_c
    double precision , dimension(3) , intent(in) :: r, V
    double precision , dimension(3) :: weight, drag,  Vinf, Vrot, interVec, liftDir,lift, r2, latlonalt,dragDir
    double precision , dimension(4) :: dmdvdt
    double precision , dimension(3) :: VwindLocal,Vwind,dvdt
    double precision :: altitude , rho, g, newden,newU,newV,newW,dm!,rotAngle
    double precision :: clat,clon,Minf,speedOfSound,CDLocal,CLLocal,VinfMag
    double precision , dimension(3,3) :: localRotation
    double precision , dimension(3) :: Tvector,uvector,rxuvector,u1vector,c1unit,c2unit,rxu1vector,d2unit,u2vector
    ! print *, 'Solving dvdt'   

    call LinearThrustInterpolation(timeList,Tlist,ntime,ct,Tvector)!getting thrust vector
    !print *,'Length ntime ',ntime
    dm = -norm(Tvector)/(Isp*9.8)
    latlonalt = latLonAltCalculation(r,ct,thetag,PlanetModel)!
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
          localRotation = transpose(matmul(rRotation(2,-clat*Pi/180),rRotation(3,thetag+omega*ct+clon*Pi/180)))
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
    weight = mass_c*(-g*r/norm(r))
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

   !!!! this section to handle offsets in thrust angle
    if (norm(Tvector)>0)then
    uvector = Tvector/norm(Tvector)
    rxuvector = cross(r,uvector)
   !uvector,rxuvector,u1vector,c1unit,c2unit,rxu1vector,d2unit
    c1unit = rxuvector/norm(rxuvector)
    u1vector = cos(thrustOffAngleRad(1))*uvector + sin(thrustOffAngleRad(1))*c1unit
    rxu1vector = cross(r,u1vector)
    c2unit = rxu1vector/norm(rxu1vector)
    d2unit = cross(u1vector,c2unit)/norm(cross(u1vector,c2unit))
    u2vector = cos(thrustOffAngleRad(2))*u1vector + sin(thrustOffAngleRad(2))*d2unit
    Tvector = norm(Tvector)*u2vector ! recalculating thrust vector
    end if
   

    dvdt = (drag + lift + weight + Tvector)/mass_c
    dmdvdt = (/dm,dvdt(1),dvdt(2),dvdt(3)/)
    ! print *,norm(weight)
    !print *,Cdlocal,speedofSound,Minf
  end function dmdvdt

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


  function latlonaltCalculation(r,t,thetag0,model)
    use constants
    use planetEarthCharacteristics
    implicit none
    ! calculates spherical lat lon 
    ! uses ellipse for altitude if model = 1

    double precision , dimension(3) , intent(in) :: r
    integer , intent(in) :: model
    double precision , intent(in) :: t
    double precision , intent(in) :: thetag0
    double precision :: normr, rdelta,phigd,Nphi,phigd_np1,tol,height
    double precision , dimension(3):: rlocal, latlonaltCalculation

    tol = 999.
    rlocal =  matmul(Rrotation(3,thetag0+omega*t),r) ! accounting for Rotating planet
    normr = norm(rlocal)
    
    height = -999.
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



  function rk45(mass,r,V,t,dtmax,dtmin)
    use constants
    implicit none
    ! my own rk45 function...currently not used, but seems to work fine

    integer :: counter
    double precision , intent(in) :: t
    double precision , dimension(3), intent(in) :: r,V
    double precision , dimension(8) :: rk45
    double precision ,dimension(3) :: resrvec,resVvec
    double precision :: resr,resV, residual,tol,delta,m1,m2,m3,m4,m5,m6
    double precision :: dm1,dm2,dm3,dm4,dm5,dm6,mf
    double precision , dimension(3) :: r1, r2, r3 ,r4,r5,r6, v1, v2, v3, v4,v5,v6, a1, a2, a3, a4,a5,a6, rf, vf
    double precision , dimension(4) :: ma1,ma2,ma3,ma4,ma5,ma6
    double precision , intent(in):: dtmax ,dtmin,mass
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
       m1 = mass
       ma1 = dmdvdt(r1,V1,m1,t)
       a1 = ma1(2:4)
       dm1 = ma1(1)

       m2 = m1 + k2_1*dm1
       r2 = r1 + k2_1*dt*v1
       v2 = v1 + k2_1*dt*a1
       ma2 = dmdvdt(r2,v2,m2,t+k2_1*dt)
       a2 = ma2(2:4)
       dm2 = ma2(1) 

       m3 = m1 + dt*(k3_2*dm1 + k3_3*dm2)
       r3 = r1 + dt*(k3_2*v1 + k3_3*v2)
       v3 = v1 + dt*(k3_2*a1 + k3_3*a2)
       ma3 = dmdvdt(r3,v3,m3,t + k3_1*dt)
       a3 = ma3(2:4)
       dm3 = ma3(1)

       m4 = m1 + dt*(k4_2*dm1 + k4_3*dm2 + k4_4*dm3)
       r4 = r1 + dt*(k4_2*v1 + k4_3*v2 + k4_4*v3)
       v4 = v1 + dt*(k4_2*a1 + k4_3*a2 + k4_4*a3)
       ma4 = dmdvdt(r4,v4,m4,t + k4_1*dt)
       a4 = ma4(2:4)
       dm4 = ma4(1)

       m5 = m1 + dt*(k5_1*dm1 + k5_2*dm2 + k5_3*dm3 + k5_4*dm4)
       r5 = r1 + dt*(k5_1*v1 + k5_2*v2 + k5_3*v3 + k5_4*v4)
       v5 = v1 + dt*(k5_1*a1 + k5_2*a2 + k5_3*a3 + k5_4*a4)
       ma5 = dmdvdt(r5,v5,m5,t+dt)
       a5 = ma5(2:4)
       dm5 = ma5(1)

       m6 = m1 + dt*(k6_2*dm1 + k6_3*dm2 + k6_4*dm3 + k6_5*dm4 + k6_6*dm5)
       r6 = r1 + dt*(k6_2*v1 + k6_3*v2 + k6_4*v3 + k6_5*v4 + k6_6*v5)
       v6 = v1 + dt*(k6_2*a1 + k6_3*a2 + k6_4*a3 + k6_5*a4 + k6_6*a5)
       ma6 = dmdvdt(r6,v6,m6,t + k6_1*dt)
       a6 = ma6(2:4)
       dm6 = ma6(1)

       !resrvec = 1/360*r1 - 128/4275*r3 - 2197/75240*r4 + .02*r5 + 2/55*r6
       ! resr = norm(resrvec)/dt
       resrvec = res_1*r1 + res_2*r3 + res_3*r4 + res_4*r5 + res_5*r6
       resVvec = res_1*v1 + res_2*v3 + res_3*v4 + res_4*v5 + res_5*v6
       resV = norm(resVvec)
       resr = norm(resrvec)
       residual = resV! + resr




       ! print *, residual,dt
       if (residual<=tol) then
          rf = r1 + dt*(ap_1*v1 + ap_2*v3 + ap_3*v4 + ap_4*v5)
          vf = v1 + dt*(ap_1*a1 + ap_2*a3 + ap_3*a4 + ap_4*a5)
          mf = m1 + dt*(ap_1*dm1 + ap_2*dm3 + ap_3*dm4 + ap_4*dm5)
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

    rk45 = (/rf(1),rf(2),rf(3),vf(1),vf(2),vf(3),dt,mf/)
    ! print *,norm(rf),dt
    !rk4(2,:)=vf
    !print *, 'Done with one step'
    !if (dt>.00001) then
    !print *,norm(rf),dt,norm(r1),residual
    !end if


  end function rk45

end subroutine propagate




