subroutine LinearCDInterpolation(x,y,n,xval,yval)
use constants  
! x => array of values in x axis
! y => array of values in y axis
! xval => value used for interpolation
! yval => value to be returned
! n => length if x or y array
  
!! NOTES: x must be in increasing order...else THIS WILL FAIL
  
  implicit none
  integer , intent(in) :: n
  integer, save :: nloc = 1
  double precision ,dimension(n), intent(in):: x,y
  double precision , intent(in) :: xval
  double precision, intent(out) :: yval
  double precision :: m
  integer :: counter,limitn


 
  if (nloc>=n) then ! to avoid segmentation faults 
     nloc = 1
  end if


   if(xval>=x(n))then ! this IF statment avoids extrapolation, instead use closest known value
      yval = y(n)
      return
    elseif(xval<=x(1))then
      yval = y(1)
      return
    endif

limitn = 100000
counter = 0
  do while(counter<limitn)
  

     if ((xval>=x(nloc)).AND.(xval<=x(nloc+1))) then
           exit
  
     elseif (xval <= x(nloc)) then
        nloc = nloc-1
     elseif (xval >= x(nloc)) then
        nloc = nloc+1
     end if
     counter = counter+1


  end do

if (counter>=limitn) then
   print *, 'ERROR in CD linear Inerpolation'
   print *,'Minf = ',xval
  ! STOP
end if

!! this section just for error handling
!! It just specifies reference points for linear extrapolation

  if (nloc >= n) then
     nloc = n-1
  elseif(nloc<1) then
     nloc =1
  end if
 
  m = (y(nloc+1)-y(nloc))/(x(nloc+1)-x(nloc))

  yval = m*(xval-x(nloc)) + y(nloc)



 ! print *, nloc,n
end


subroutine LinearCLInterpolation(x,y,n,xval,yval)
use constants  
! x => array of values in x axis
! y => array of values in y axis
! xval => value used for interpolation
! yval => value to be returned
! n => length if x or y array
  
!! NOTES: x must be in increasing order...else THIS WILL FAIL
  
  implicit none
  integer , intent(in) :: n
  integer, save :: nloc2 = 1
  double precision ,dimension(n), intent(in):: x,y
  double precision , intent(in) :: xval
  double precision, intent(out) :: yval
  double precision :: m
  integer :: counter,limitn


 
  if (nloc2>=n) then ! to avoid segmentation faults 
     nloc2 = 1
  end if

   if(xval>=x(n))then ! this IF statment avoids extrapolation, instead use closest known value
      yval = y(n)
      return
    elseif(xval<=x(1))then
      yval = y(1)
      return
    endif



limitn = 100000
counter = 0
  do while(counter<limitn)
  

     if ((xval>=x(nloc2)).AND.(xval<=x(nloc2+1))) then
           exit

     elseif (xval <= x(nloc2)) then
        nloc2 = nloc2-1
     elseif (xval >= x(nloc2)) then
        nloc2 = nloc2+1
     end if
     counter = counter+1


  end do


!! this section just for error handling
!! It just specifies reference points for linear extrapolation

  if (nloc2 >= n) then
     nloc2 = n-1
  elseif(nloc2<1) then
     nloc2 =1
  end if
 
  m = (y(nloc2+1)-y(nloc2))/(x(nloc2+1)-x(nloc2))

  yval = m*(xval-x(nloc2)) + y(nloc2)


if (counter>=limitn) then
   print *, 'ERROR in CL linear Inerpolation'
   print *,'Minf = ',xval
  ! STOP
end if





 ! print *, nloc,n
end subroutine LinearCLInterpolation



subroutine LinearThrustInterpolation(x,y,n,xval,yval)
use constants  
! x => array of values in x axis
! y => array of values in y axis
! xval => value used for interpolation
! yval => value to be returned
! n => length if x or y array
  
!! NOTES: x must be in increasing order...else THIS WILL FAIL
  
  implicit none
  integer , intent(in) :: n
  integer, save :: nloc3 = 1
  double precision ,dimension(n), intent(in):: x
  double precision ,dimension(n,3), intent(in):: y !n x 3 matrix
  double precision , intent(in) :: xval
  double precision, dimension(3),intent(out) :: yval
  double precision :: m1,m2,m3
  integer :: counter,limitn


 
  if (nloc3>=n) then ! to avoid segmentation faults 
     nloc3 = 1
  end if

   if(xval>=x(n))then ! this IF statment avoids extrapolation, instead use closest known value
      yval(1) = y(n,1)
      yval(2) = y(n,2)
      yval(3) = y(n,3)

      return
    elseif(xval<=x(1))then
      yval(1) = y(1,1)
      yval(2) = y(1,2)
      yval(3) = y(1,3)   
   return
    endif



limitn = 100000
counter = 0
  do while(counter<limitn)
  

     if ((xval>=x(nloc3)).AND.(xval<=x(nloc3+1))) then
           exit

     elseif (xval <= x(nloc3)) then
        nloc3 = nloc3-1
     elseif (xval >= x(nloc3)) then
        nloc3 = nloc3+1
     end if
     counter = counter+1


  end do


!! this section just for error handling
!! It just specifies reference points for linear extrapolation
 ! PRINT *, 'Nloc3 ',nloc3,n
  if (nloc3 >= n) then
     nloc3 = n-1
  elseif(nloc3<1) then
     nloc3 =1
  end if

  m1 = (y(nloc3+1,1)-y(nloc3,1))/(x(nloc3+1)-x(nloc3))
  m2 = (y(nloc3+1,2)-y(nloc3,2))/(x(nloc3+1)-x(nloc3))
  m3 = (y(nloc3+1,3)-y(nloc3,3))/(x(nloc3+1)-x(nloc3))
  yval(1) = m1*(xval-x(nloc3)) + y(nloc3,1)
  yval(2) = m2*(xval-x(nloc3)) + y(nloc3,2)
  yval(3) = m3*(xval-x(nloc3)) + y(nloc3,3)



if (counter>=limitn) then
   print *, 'ERROR in Thrust linear Inerpolation'
   print *,'time = ',xval
   STOP
end if





 ! print *, nloc,n
end subroutine LinearThrustInterpolation
