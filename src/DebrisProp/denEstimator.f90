subroutine atmosInterpolation(newden,newU,newV,newW,newalt,altitude,density,U,V,W,n)
use constants  
! newden = interpolated density
  ! newU/newV/newW = interpolated U/V/W velocities
  ! newalt = altitude[m] to be interpolated
  
  ! altitude = vector with altitude values [m] !!!MUST BE IN DECREASING ORDER!!!!! ********
  
  ! density/U/V/W = list of density/U/V/W of size n
  ! n = length of the vectors to do linear interpolation from
  
  
  implicit none
  integer , intent(in) :: n
  integer, save :: nloc = 1
  double precision ,dimension(n), intent(in):: altitude, density, u, v, w
  double precision , intent(in) :: newalt
  double precision, intent(out) :: newden, newU, newV, newW
  double precision :: ratio
  integer :: counter,limitn


 
  if (nloc>=n) then
     nloc = 1
  end if
! print *, nloc
!print *,altitude(1),altitude(n)
!print *, density(1),density(n)
limitn = 100000
counter = 0

if (newalt>altitude(1)) then
     newden = 0
     newU = 0
     newV = 0
     newW = 0
     return
end if




  do while(counter<limitn)
     if (newalt+10==newalt) then
        exit
     end if

     if ((newalt<=altitude(nloc)).AND.(newalt>=altitude(nloc+1))) then
           exit

     elseif (newalt >= altitude(1)) then
        nloc = 1
        exit
     elseif (newalt <=0) then
        nloc = n
        exit
     elseif (newalt <= altitude(nloc)) then
        nloc = nloc+1
     elseif (newalt>=altitude(nloc)) then
        nloc = nloc-1
     end if
     counter = counter+1
! if (counter==limitn) then
!    print*, newalt
!  end if

  end do


  if (nloc >= n) then
     nloc = n-1
  elseif(nloc==0) then
     nloc =1
  elseif(nloc<0) then
     nloc = 1
     print *, 'Warning in denEstimator.f90: nloc < 0 '
  end if
 
  ratio = (newalt-altitude(nloc))/(altitude(nloc+1)-altitude(nloc))
  newden = density(nloc) + ratio*(density(nloc+1)-density(nloc))
  newU = U(nloc) + ratio*(U(nloc+1)-U(nloc))
  newV = V(nloc) + ratio*(V(nloc+1)-V(nloc))
  newW = W(nloc) + ratio*(W(nloc+1)-W(nloc))

  if (newalt+10==newalt) then
     print *, 'Careful atmospheric interpolation unsuccessful. Altitude already NAN'
     newden = 0
     newU = 0
     newV = 0
     newW = 0

  end if


 ! print *, nloc,n
end
