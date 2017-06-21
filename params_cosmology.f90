module params_cosmo
  implicit none
!!cosmological parameters
  real,parameter :: pi=3.1415926
  integer,parameter :: Nc=256 ! number of grid cells per dimension
  real,parameter :: H_o=100. !(km/s*h/Mpc)
  real,parameter :: Omega_m=0.32   !Plank
  real,parameter :: Omega_l=1-Omega_m
  real,parameter :: zi=0.5,zf=0. ! initial/final redshift
  real,parameter ::ai=1./(1.+zi), af=1./(1.+zf) ! initial/final scale factor
  real,parameter :: gamma=5./3.

!!physical values
  real,parameter :: Gr=4.301e-9 !(km^2*Mpc/M_sun/s^2) ! gravitational constant
  real,parameter :: C_s=500 ! (km/s)
  real,parameter :: L=5. !(Mpc/h)
  real,parameter :: Ma=1.0e14 !(M_sun/h)  ! mass of source
  real,parameter :: Rho_c=3*Omega_m*H_o**2/(8*pi*Gr)
!  real,parameter :: Amp=1000

!!computational values
  real,parameter :: Ampli=100 
  real,parameter :: den=1.
!  real,parameter :: le=L/N
!  real mass
!  real gra
!  real t
!  real cs

!!units
  real,parameter :: U_den=Rho_c/den
  real,parameter :: U_l=L/Nc
  real,parameter :: U_mass=U_den*U_l**3
  real U_g
  real U_v
  real U_time
!  real p


  real,parameter :: mass=Ma/U_mass
  real,parameter :: gra=Ampli/mass
!  U_mass=U_den*U_l**3
!  mass=M/U_mass
!  gra=Amp/mass
!  U_g=Gr/gra
!  U_time=1./sqrt(U_den*U_g)
!  U_v=U_l/U_time
!  cs=C_s/U_v
!  p=cs**2*den/gamma











end module params_cosmo
