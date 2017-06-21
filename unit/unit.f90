program units
  use parameters
  implicit none
!!physical value
  real,parameter :: pi=3.1415926
!  integer,parameter :: N=n ! number of grid cells per dimension
  integer,parameter :: nt=100 ! number of time steps
  real,parameter :: H_o=100. !(km/s*h/Mpc)
  real,parameter :: Omega_m=0.32   !Plank
  real,parameter :: Omega_l=1-Omega_m
  real z(nt) ! redshift
  real,parameter :: zi=1.5,zf=0. ! initial/final redshift
  real a ! scale factor
  real,parameter ::ai=1./(1.+zi), af=1./(1.+zf) ! initial/final scale factor
  real,parameter :: Gr=4.301e-9 !(km^2*Mpc/M_sun/s^2) ! gravitational constant
  real,parameter :: gra=1/(6*pi)!5e-8!1/(160*pi) ! 1/(20*pi)
  real,parameter :: C_s=500 ! (km/s)
  real,parameter :: Rho_c=3*Omega_m*H_o**2/(8*pi*Gr)!rho_c=3*Omega_m*H_o**2/(8*pi*Gr*a**3) !critical density
  real,parameter :: cs=1!/sqrt(2.)
  real,parameter :: den=1.
  real,parameter :: mass=Amp/gra !(M_sun/h) ! mass of source
  real Mass_s !(M_sun/h)
  real d ! grid matter density
  real L !(Mpc) ! physical box siz 
!  real,parameter :: le=L/n !(Mpc)! grid length unit
!  real,parameter :: rs=ratio*le !(Mpc)
!!others
  real moment(5,nt)
  real drag(nt)
  real force(nt)
  real gridtime(nt)
  real time(nt)
  real da1,da2
  real(4), parameter      :: wde = -1.0
  real(4), parameter      :: w_a = 0.0
!!Units:
!  real U_force !!unit of force
  real,parameter :: U_g=Gr/gra
  real,parameter :: U_den=Rho_c/den
  real,parameter :: U_time=1./sqrt(U_g*U_den)!(s/km*Mpc/h)   !!unit of time
!  real U_Mass  !!unit of mass
  real,parameter :: U_v=C_s/cs
  real,parameter :: U_l=U_v*U_time
  real,parameter :: U_mass=U_den*U_l**3
  real,parameter :: U_phi=U_g*U_mass/U_l
  real,parameter :: U_force=U_phi/U_l


!  real V_Gr=Gr/gra
!  real,parameter :: U_den=Rho_c/den
!  U_time=2*a**2/(3*H_o*sqrt(Omega_m))!(s/km*Mpc/h)       ! When  U_Gr=6*pi*Gr

 ! U_mass=U_den*U_l**3 !(M_sun/h)  ! When U_Gr=6*pi*Gr
 ! U_force=U_phi/U_l!3*Omega_m*H_o**4*le/(8*pi*Gr)*(a**3*(Omega_m+(a**3*(1-Omega_m))))
!  call printparams 


!  cs=1/sqrt(2.)
  
!  V_Gr=3*Omega_m*H_o**2*le/(8*pi*C_s/cs)

!
!  Mass_s=mass*U_mass
!  L=N*1.*U_l
!  print*,'Mass=',Mass_s,'L=',L,'cs=',cs,'v=',cs*M
!
!  write(npd,'(i3)') n
!  write(mn,'(f4.1)') M
!  write(am,'(f4.1)') Amp
!  write(rat,'(i3)') ratio
!  fn='../../olddata/me_n'//trim(adjustl(npd))//'M'//trim(adjustl(mn))//'g'//&
!     trim(adjustl(am))//'rs'//trim(adjustl(rat))//'_enc.dat'
!  print*,fn
!  open(100,file=fn,access='stream')
!  read (100) moment
!  close(100)
!  fn='../../olddata/df_n'//trim(adjustl(npd))//'M'//trim(adjustl(mn))//'g'//&
!     trim(adjustl(am))//'rs'//trim(adjustl(rat))//'_enc.dat'
!  print*,fn
!  open(15,file=fn,access='stream') 
!  read(15) drag
!  close(15)
!  gridtime=moment(1,:)-cshift(moment(1,:),-1)
!  gridtime(1)=0
!  a=ai
!  do i=1,nt
!     gridtime(i)=gridtime(i)!2*a**2/3/sqrt(Omega_m*H_o**2)
!     call expansion(a,gridtime(i),da1,da2)
!     a=a+da1+da2
!     z(i)=1./a-1.
!!     U_force=3*Omega_m*H_o**4*U_l/(8*pi*Gr)*(a**3*(Omega_m+(a**3*(1-Omega_m))))
!     force(i)=drag(i)*U_force*a**3
!     print*,moment(1,i),'distance=',moment(1,i)*cs*M,' gridtime = z',z(i),' drag force =',force(i)
!  enddo
contains
  subroutine expansion(a0,dt0,da1,da2)
    implicit none

!#    include "cubepm.par"

    real :: a0,dt0,dt_x,da1,da2
    real(8) :: a_x,adot,addot,atdot,arkm,a3rlm,omHsq
    real(8), parameter :: e = 2.718281828459046


!#ifdef Chaplygin
!    call Chaplygin(a0,dt0,da1,da2)
!    return
!    write(*,*) '***** IF I SEE THIS, SOMETHING IS WRONG!! *****'
!#endif

    !! Expand Friedman equation to third order and integrate
    dt_x=dt0/2
    a_x=a0
!    omHsq=4.0/9.0
    omHsq=8.*pi*gra/3.
!    a3rlm=a_x**3*Omega_l/Omega_m
    a3rlm=a_x**(-3*wde)*Omega_l/Omega_m
!    a3rlm=a_x**(-3*wde - 3*w_a)*(Omega_l/Omega_m)*e**(3*w_a*(a_x - 1))
    arkm=a_x*(1.0-Omega_m-Omega_l)/Omega_m

    adot=sqrt(omHsq*a_x**3*(1.0+arkm+a3rlm))
!    addot=a_x**2*omHsq*(1.5+2.0*arkm+3.0*a3rlm)
    addot=a_x**2*omHsq*(1.5+2.0*arkm+1.5*(1.0-wde)*a3rlm)
!    addot=a_x**2*omHsq*(1.5+2.0*arkm+1.5*(1.0-wde + w_a*(a_x - 1))*a3rlm)
!    atdot=a_x*adot*omHsq*(3.0+6.0*arkm+15.0*a3rlm)
    atdot=a_x*adot*omHsq*(3.0+6.0*arkm+1.5*(2.0-3.0*wde)*(1.0-wde)*a3rlm)
!    atdot=a_x*adot*omHsq*(3.0+6.0*arkm+1.5*(3*w_a**2*a_x**2 +
!    6*w_a*(1-wde-w_a)+ (2.0-3.0*(wde+w_a))*(1.0-(wde+w_a)))*a3rlm)

    da1=adot*dt_x+(addot*dt_x**2)/2.0+(atdot*dt_x**3)/6.0

    a_x=a0+da1
    omHsq=4.0/9.0
!    a3rlm=a_x**3*Omega_l/Omega_m
    a3rlm=a_x**(-3*wde)*Omega_l/Omega_m
!    a3rlm=a_x**(-3*wde - 3*w_a)*(Omega_l/Omega_m)*e**(3*w_a*(a_x - 1))
    arkm=a_x*(1.0-Omega_m-Omega_l)/Omega_m

    adot=sqrt(omHsq*a_x**3*(1.0+arkm+a3rlm))
!    addot=a_x**2*omHsq*(1.5+2.0*arkm+3.0*a3rlm)
    addot=a_x**2*omHsq*(1.5+2.0*arkm+1.5*(1.0-wde)*a3rlm)
!    addot=a_x**2*omHsq*(1.5+2.0*arkm+1.5*(1.0-wde + w_a*(a_x - 1))*a3rlm)
!    atdot=a_x*adot*omHsq*(3.0+6.0*arkm+15.0*a3rlm)
    atdot=a_x*adot*omHsq*(3.0+6.0*arkm+1.5*(2.0-3.0*wde)*(1.0-wde)*a3rlm)
!    atdot=a_x*adot*omHsq*(3.0+6.0*arkm+1.5*(3*w_a**2*a_x**2 +
!    6*w_a*(1-wde-w_a)+ (2.0-3.0*(wde+w_a))*(1.0-(wde+w_a)))*a3rlm)

    da2=adot*dt_x+(addot*dt_x**2)/2.0+(atdot*dt_x**3)/6.0

  end subroutine expansion
end 
