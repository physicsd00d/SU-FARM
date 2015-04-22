module constants

! Constants and precision
!##############################################
 !integer, parameter :: b8 = selected_real_kind(15)
 double precision , parameter :: Pi = 3.14159265358979323846
! Definitions for RK45
 double precision , parameter :: k2_1 = 0.25
 double precision , parameter :: k3_1 = 3.0/8.0
 double precision , parameter :: k3_2 = 3.0/32.0
 double precision , parameter :: k3_3 = 9.0/32.0
 double precision , parameter :: k4_1 = 12.0/13.0
 double precision , parameter :: k4_2 = 1932.0/2197.0
 double precision , parameter :: k4_3 = -7200.0/2197.0
 double precision , parameter :: k4_4 = 7296.0/2197.0
 double precision , parameter :: k5_1 = 439.0/216.0
 double precision , parameter :: k5_2 = -8.0
 double precision , parameter :: k5_3 = 3680.0/513.0
 double precision , parameter :: k5_4 = -845.0/4104.0
 double precision , parameter :: k6_1 = 0.5
 double precision , parameter :: k6_2 = -8.0/27.0
 double precision , parameter :: k6_3 = 2.0
 double precision , parameter :: k6_4 = -3544.0/2565.0
 double precision , parameter :: k6_5 = 1859.0/4104.0
 double precision , parameter :: k6_6 = -11.0/40.0
 double precision , parameter :: res_1 = 1.0/360.0
 double precision , parameter :: res_2 = -128.0/4275.0
 double precision , parameter :: res_3 = -2197.0/75240
 double precision , parameter :: res_4 = 1.0/50.0
 double precision , parameter :: res_5 = 2.0/55.0
 double precision , parameter :: ap_1 = 25.0/216.0
 double precision , parameter :: ap_2 = 1408.0/2565.0
 double precision , parameter :: ap_3 = 2197.0/4104.0
 double precision , parameter :: ap_4 = -1.0/5.0

 integer gene_size
 end
