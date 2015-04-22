subroutine getspeedofsound(h,c)
use constants  
! h is the altitude in meters
! c is im m/sec  
  
  implicit none
  double precision , intent(in):: h
  double precision, intent(out) :: c
  double precision :: m,c2,c1,h1,h2,hkm

  hkm = h*.001
  
if (hkm<0)then
   c = 340.3
   return
elseif(hkm>86) then
   c = 274.1
   return
end if


if (hkm<=11)then
    c2 = 295.2
    c1 = 340.3
    h2 = 11
    h1 = 0
elseif (hkm<=20)then
    c2 = 295.11
    c1 = 295.1
    h2 = 20
    h1 = 11
elseif (hkm<=32)then
    c2 = 303.0
    c1 = 295.1
    h2 = 32
    h1 = 20
elseif (hkm<=47)then
    c2 = 329.2
    c1 = 303
    h2 = 47
    h1 = 32
elseif (hkm<=48)then
    c2 = 329.8
    c1 = 329.2
    h2 = 48
    h1 = 47
elseif (hkm<=51)then
    c2 = 329.81
    c1 = 329.8
    h2 = 51
    h1 = 48
elseif (hkm<=52)then
    c2 = 328.8
    c1 = 329.82
    h2 = 52
    h1 = 51
elseif (hkm<=72)then
    c2 = 293.4
    c1 = 328.8
    h2 = 72
    h1 = 52
elseif (hkm<=86)then
    c2 = 274.1
    c1 = 293.4
    h2 = 86
    h1 = 72
end if

m = (c2-c1)/(h2-h1)
c = m*(hkm-h1) + c1

end
