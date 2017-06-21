! TVD split MHD code
! copyright (C) 2001,2003, Ue-Li Pen
! written November 2001 by Ue-Li Pen, pen@cita.utoronto.ca
!This program is free software; you can redistribute it and/or
!modify it under the terms of the GNU General Public License
!as published by the Free Software Foundation; either version 2
!of the License, or (at your option) any later version.
!
!This program is distributed in the hope that it will be useful,
!but WITHOUT ANY WARRANTY; without even the implied warranty of
!MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!GNU General Public License for more details.
!
!You should have received a copy of the GNU General Public License
!along with this program; if not, write to the Free Software
!Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
!
! debugged to pass alven test March 29, 2002
!
! release 0.1 May 2003
! release 0.2 Sept 2009 : fixed b2 interpolation in cfl
!
! General Notes:
! restrictions: requires nx=nz or nz=1
! see http://arxiv.org/astroph/abs/astro-ph/0305088
! or http://www.cita.utoronto.ca/~pen/MHD
! for questions contact pen@cita.utoronto.ca
use parameters
use params_cosmo
use workfile
use output
call init(u,nx,ny,nz,M,c,rho0)
call initg(g,phi,u,nx,ny,nz,Amp,ratio)
call printparams
phi=0
a=ai
do
   fpp=fp
   fp=f
!if (a<=af) then
   call me(u,nx,ny,nz,t)
!endif
!if (a>af .and. t>n/vx) exit
!if (t<=n/vx) then
   call writeu(u,n)
   call writev(u,n,vx)
   call writep(u,n,vx,rho0)
   call vcrl(u,n)
   call dragforce(u,g,n,drag)
   call entropy(u,n)
!endif
   call force(u,g,n,df)
   iter=iter+1
   f=df(1)
!   if (iter>=3) exit
   !if (p==pp) exit
!   if (iter==400) exit
!    if (t>n/vx) exit
!   if(fpp==fp .and. fp==f .and. iter>10) exit
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!test!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   if (iter==50) call slice(u,n)    !test!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   if (iter==51) exit               !test!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!test!!!!!!!!!!!!!!!!!!!!!!!!!!!
   dt=0.9*cfl(u)  
   dt=min(dt,minval(1/sqrt(3*abs(g))))
!   print*,'cfl,dt_cosmo,dt_grav=',0.9*cfl(u),dt_cosmo(a),minval(1/sqrt(3*abs(g)))
if (iter>2000 .or. t>n/vx) exit
!   if(a>1) exit
   call expansion(a,dt,da1,da2)
   a=a+da1+da2
!   dt=min(dt,(tf-t)/2)
   t=t+2*dt
   call fluidx(u,nx,ny,nz,dt)
! the y sweep
   call transpose12(u,u,nx,ny,nz)
   call fluidx(u,ny,nx,nz,dt)
! z sweep
   call transpose13(u,u,ny,nx,nz)
   call fluidx(u,nz,nx,ny,dt)
   call gravity(g,u,n,2*dt)
!   call gravity(g,u,n,dt)
   call fluidx(u,nz,nx,ny,dt)
! back
   call transpose13(u,u,nz,nx,ny)
   call fluidx(u,ny,nx,nz,dt)
! x again
   call transpose12(u,u,ny,nx,nz)
   call fluidx(u,nx,ny,nz,dt)
   if (mod(iter,10) .eq. 1) then
      write(*,*) 't=',t,iter,a,maxval(u(1,:,:,:)),minval(u(1,:,:,:))
      print*,'cfl,dt_cosmo,dt_grav=',0.9*cfl(u),dt_cosmo(a),minval(1/sqrt(3*abs(g)))
!     print*, 'a=',a
   endif
print*, 'a=',a
enddo
call closefile

contains
  function cfl(u)
    real, dimension(:,:,:,:)::u
    real cfl
! locals
    real v,ps,p,c,gamma
    integer i,j,k,ip,jp,kp

    gamma=5./3.
    c=0
    !$omp parallel do private(j,i,kp,jp,ip,v,ps,p) reduction(max:c)
    do k=1,nz
       do j=1,ny
          do i=1,nx
             kp=mod(k,nz)+1
             jp=mod(j,ny)+1
             ip=mod(i,nx)+1
             v=maxval(abs(u(2:4,i,j,k)/u(1,i,j,k)))
             ps=(u(5,i,j,k)-sum(u(2:4,i,j,k)**2,1)/u(1,i,j,k)/2)*(gamma-1)
             p=ps
!   c=max(c,v+sqrt(abs(  (sum(b(:,i,j,k)**2)*(gamma-1)+gamma*p)/u(1,i,j,k))))
             c=max(c,v+sqrt(abs(gamma*p/u(1,i,j,k))))
          end do
       end do
    end do
    cfl=1./c
  end function cfl

  function dt_cosmo(a_x)
    real,parameter :: ratio = 0.02
    real,parameter :: wde=-1.
    real da,a_x
    real dt_cosmo
    real(8) :: omHsq,arkm,a3rlm
    omHsq=8.*pi*gra/3.
    a3rlm=a_x**(-3*wde)*Omega_l/Omega_m
    arkm=a_x*(1.0-Omega_m-Omega_l)/Omega_m
    da=ratio*a_x
    dt_cosmo=da/sqrt(omHsq*a_x**3*(1.0+arkm+a3rlm))
  end function dt_cosmo

  subroutine expansion(a0,dt0,da1,da2)
    implicit none

!#    include "cubepm.par"

    real :: a0,dt0,dt_x,da1,da2
    real(8) :: a_x,adot,addot,atdot,arkm,a3rlm,omHsq
    real(8), parameter :: e = 2.718281828459046
    real,parameter :: wde=-1.
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

    da2=adot*dt_x+(addot*dt_x**2)/2.0+(atdot*dt_x**3)/6.0

  end subroutine expansion

end
