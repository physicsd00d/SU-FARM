module planetEarthCharacteristics
use constants

! Planet, atmospheric parameters and other constants
!############################################## 

 !double precision , parameter :: planetRadius = 6371.009*1000 !planet mean Radius [m]
 double precision , parameter :: planetRadius = 6371000.8D0 !
 double precision , parameter :: muPlanet = 398601.2D9 !m3/sec^2
 double precision , parameter :: omega = 2*Pi/(86164.0906) !angular velocity 
 double precision , parameter :: Rgas = 287.0D0
 double precision , parameter :: T0 = 288.15D0
 double precision , parameter :: gamma = 1.4D0
 double precision , parameter :: rho0 = 1.225D0
 double precision , parameter :: g0 = 9.8D0
 double precision , parameter :: Hscale = Rgas*T0/g0
 double precision , parameter :: Requator = 6378137.D0! [m] planet Equatorial Radius
 double precision , parameter :: ecc = 0.081819221456D0 ! Earth eccentricity. Used for oblateness calculations

end
