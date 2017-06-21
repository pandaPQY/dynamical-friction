program Unity
  implicit none
!!cosmological parameters
  real,parameter :: pi=3.1415926
  integer,parameter :: N=256 ! number of grid cells per dimension
!  integer,parameter :: nt=100 ! number of time steps
  real,parameter :: H_o=100. !(km/s*h/Mpc)
  real,parameter :: Omega_m=0.32   !Plank
  real,parameter :: Omega_l=1-Omega_m
 ! real z(nt) ! redshift
 ! real,parameter :: zi=0.5,zf=0. ! initial/final redshift
  real,parameter :: a=1 ! scale factor
 ! real,parameter ::ai=1./(1.+zi), af=1./(1.+zf) ! initial/final scale factor
  real,parameter :: gamma=5./3.

!!physical values
  real,parameter :: Gr=4.301e-9 !(km^2*Mpc/M_sun/s^2) ! gravitational constant
 ! real,parameter :: gra=1/(160*pi) ! 1/(20*pi)
  real,parameter :: C_s=500 ! (km/s)
  real,parameter :: L=5. !(Mpc/h)
  real,parameter :: M=1.0e14 !(M_sun/h)
  real,parameter :: Rho_c=3*Omega_m*H_o**2/(8*pi*Gr) !rho_c=3*Omega_m*H_o**2/(8*pi*Gr*a**3) !critical density
  real,parameter :: Amp=0.0001
!!computational values
  real,parameter :: den=1 
!  real,parameter :: le=L/N
  real mass
  real gra
  real t
  real cs
  
!!units
  real,parameter :: U_den=Rho_c/den
  real,parameter :: U_l=L/N
  real,parameter :: U_mass=U_den*U_l**3
  real U_g
  real U_v
  real U_time
  real p


!  real,parameter :: mass=M/U_mass
!  U_mass=U_den*U_l**3
  mass=M/U_mass
  gra=Amp/mass
  U_g=Gr/gra
  U_time=1./sqrt(U_den*U_g)
  U_v=U_l/U_time
  cs=C_s/U_v
  p=cs**2*den/gamma
  print*,'cs=',cs,'p=',p,'mass=',mass,'gra=',gra,' Prho_crit=',U_den
  

end
