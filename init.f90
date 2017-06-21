!!!!!!!!!!!!change potential model: have to change  3 places !!!!!!!!!!!!!!!!!!
module parameters                            
implicit none
integer, parameter :: n=128,nx=n,ny=n,nz=n
real, dimension(5,nx,ny,nz) :: u
real, dimension(3,nx,ny,nz) :: g
real, dimension(nx,ny,nz) :: phi
real t,dt,tf
integer iter,i
real, parameter :: M=2.0
!integer, parameter :: Amp=100
real, parameter :: r200=50.
real(8),parameter :: Amp=1
integer, parameter :: ratio=0   !rm:5, rs:20, r:20  rh: 5                                             !!!!!!!here!!!!!!!!!!!!!!!!!!!
real, parameter :: rho0=1.
real vx,c,df(3),drag(3,n,n,n)
real a,da1,da2
real p,pp,f,fp,fpp
character(10) npd,mn,am,rat
character(100) fn
logical check(18)

contains

subroutine init(u,nx,ny,nz,M,c,rho0)
  implicit none
  integer nx,ny,nz
  real u(5,nx,ny,nz)
  real p0,c,M,vx,rho0
  integer i
  real, parameter :: epsilon=.1
  real, parameter :: gamma=5./3.
!  p0=3./5./2.
!  p0=1.362328
  p0=0.6811641/4 !cosmo galaxy
!  p0=0.006811641 
  u=0
  u(1,:,:,:)=rho0
!  u(1,n/2,n/2,n/2)=100*rho0
  c=sqrt(gamma*p0/rho0)
  vx=M*c
  u(2,:,:,:)=vx*rho0
  u(5,:,:,:)=p0/(gamma-1.) !*1.5
  u(5,:,:,:)=u(5,:,:,:)+sum(u(2:4,:,:,:)**2,1)/u(1,:,:,:)/2
!  u(5,:,:,:)=u(5,:,:,:)+p0*u(2,:,:,:)*5./3./(2./3.)
  return
end subroutine init

subroutine initg(g,phi,u,nx,ny,nz,Amp,ratio)
  implicit none
  integer nx,ny,nz
  real g(3,nx,ny,nz)
  real phi(nx,ny,nz)
  real u(5,nx,ny,nz)
  integer i,j,k
  integer ratio
  real(8) Amp
!  integer Amp
  real radi
  real sigma
  real C1,C2
  sigma=real(nx)/real(ratio)
  g=0
  phi=0
  C1=Amp*(log(1+r200/ratio)-r200/(r200+ratio))
  C2=C1/r200-Amp/r200*log(1+r200/ratio)
  !!$omp parallel do default(none) shared(g) private(i,j,k)
  do k=1,nz
     do j=1,ny
        do i=1,nx
           radi=sqrt((i-0.5-nx/2)**2+(j-0.5-ny/2)**2+(k-0.5-nz/2)**2)                              !here!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! phi(i,j,k)=-Amp*exp(-1.*radi**2/2/sigma**2)   ! gaussian   r           20
          ! phi(i,j,k)=Amp*-1./radi                       ! 1/r        nothing

!           if (radi<=r200) then
!              phi(i,j,k)=Amp*-1./radi*log(1.+radi/ratio)    ! galaxy     rs          20
!           else
!              phi(i,j,k)=-1*C1/radi+C2
!           endif
           phi(i,j,k)=Amp*-1./(sqrt(ratio**2+radi**2))   ! plummer     rm           5
          ! phi(i,j,k)=Amp*-1./(radi+ratio)               ! Hernquist   rh           5
          ! phi(i,j,k)=Amp*-1./radi**3                 ! square root rq
!           radi=sqrt((k-0.5-nx/2)**2+(j-0.5-ny/2)**2)
!           phi(i,j,k)=Amp*-1./sqrt((i-0.5-nx/2)**2+(ratio*radi)**2)  !reh
!           phi(i,j,k)=Amp*-1./sqrt((ratio*(i-0.5-nx/2))**2+radi**2)  !rep
        enddo
     enddo
  enddo
  g(1,:,:,:)=-1./2*(cshift(phi,1,1)-cshift(phi,-1,1))
  g(2,:,:,:)=-1./2*(cshift(phi,1,2)-cshift(phi,-1,2))
  g(3,:,:,:)=-1./2*(cshift(phi,1,3)-cshift(phi,-1,3))
  do k=1,nz
     do j=1,ny
        do i=1,nx
           radi=sqrt((i-0.5-nx/2)**2+(j-0.5-ny/2)**2+(k-0.5-nz/2)**2)
           if (radi>20) g(:,i,j,k)=0
        enddo
     enddo
  enddo
!   g(1,:,:,:)=g(1,:,:,:)+(phi-cshift(phi,-1,1))
!   g(1,:,:,:)=g(1,:,:,:)+(cshift(phi,-1,2)-cshift(cshift(phi,-1,2),-1,1))
!   g(1,:,:,:)=g(1,:,:,:)+(cshift(phi,-1,3)-cshift(cshift(phi,-1,3),-1,1))
!   g(1,:,:,:)=g(1,:,:,:)+(cshift(cshift(phi,-1,2),-1,3)-cshift(cshift(cshift(phi,-1,2),-1,3),-1,1))
!   g(1,:,:,:)=-1.*g(1,:,:,:)/4.
!   g(2,:,:,:)=g(2,:,:,:)+(phi-cshift(phi,-1,2))
!   g(2,:,:,:)=g(2,:,:,:)+(cshift(phi,-1,1)-cshift(cshift(phi,-1,1),-1,2))
!   g(2,:,:,:)=g(2,:,:,:)+(cshift(phi,-1,3)-cshift(cshift(phi,-1,3),-1,2))
!   g(2,:,:,:)=g(2,:,:,:)+(cshift(cshift(phi,-1,1),-1,3)-cshift(cshift(cshift(phi,-1,1),-1,3),-1,2))
!   g(2,:,:,:)=-1.*g(2,:,:,:)/4.
!   g(3,:,:,:)=g(3,:,:,:)+(phi-cshift(phi,-1,3))
!   g(3,:,:,:)=g(3,:,:,:)+(cshift(phi,-1,2)-cshift(cshift(phi,-1,2),-1,3))
!   g(3,:,:,:)=g(3,:,:,:)+(cshift(phi,-1,1)-cshift(cshift(phi,-1,1),-1,3))
!   g(3,:,:,:)=g(3,:,:,:)+(cshift(cshift(phi,-1,2),-1,1)-cshift(cshift(cshift(phi,-1,2),-1,1),-1,3))
!   g(3,:,:,:)=-1.*g(3,:,:,:)/4.
!  g(5,:,:,:)=sum(g(2:4,:,:,:)*u(2:4,:,:,:),1)
  write(*,*) 'sigma=',sigma
!  phi=0
!  g=0
  !!$omp end parallel do
end subroutine initg

subroutine printparams
check= (/   1,    1,     0,    0,     0,    0,    0,     0,    0,      0,      0,       0,       1,   1,     1,     1,      1,           1/)
!          rho    me     vx    px     vy    vz    py     pz   vdiv    vcrlz   vcrly     vcrlx    df   dragx dragy  dragz  entropyy  entropy
!           1     2      3     4      5     6     7      8     9       10      11       12       13   14     15     16      17         18
check=1
!check(2)=1
!check(13)=1
!check(1)=1
!check(3:8)=1
!check(9:12)=1
!check(1)=1
!check(2)=1
! B field is stored on the left side of each cell
write(npd,'(i3)') n
write(mn,'(f4.1)') M
write(am,'(f10.4)') Amp
!write(am,'(i4)') Amp
write(rat,'(i3)') ratio
!call init(u,nx,ny,nz,M,c,rho0)
vx=M*c
print*,'**************************************************************'
print*,'parameters:                          *************************'
print*,'**************************************************************'
print*,'n    =',n
print*,'Amp  =',Amp
print*,'r200 =',r200
print*,'ratio=',ratio
print*,'Mach =',M
print*,'c    =',c
print*,'vx   =',vx
print*,'rho0 =',rho0
if (M>1) write(*,*) 'supersonic'
if (M==1) write(*,*) 'critical'
if (M<1) write(*,*) 'subsonic'
!call initg(g,phi,u,nx,ny,nz,Amp,ratio)
write(*,*) 'gravity information:'
print*,'max potential=',maxval(phi),maxval(g)
print*,'min potential=',minval(phi),minval(g)
print*,'*************************************************************'
tf=nx*4
t=0
iter=0
print*,'writing in :                      **************************'                               !!!!!!!!here!!!!!!!!!!!!!!!
if (check(1)==1) then
fn='lgrho_n'//trim(adjustl(npd))//'M'//trim(adjustl(mn))//'g'//&
     trim(adjustl(am))//'rm'//trim(adjustl(rat))//'_cut4.dat'
print*,fn
open(99,file=fn,status='replace',access='stream')
endif
if (check(2)==1) then
fn='me_n'//trim(adjustl(npd))//'M'//trim(adjustl(mn))//'g'//&
     trim(adjustl(am))//'rm'//trim(adjustl(rat))//'_cut4.dat'
print*,fn
open(100,file=fn,status='replace',access='stream')
endif
if (check(3)==1) then
fn='vx_n'//trim(adjustl(npd))//'M'//trim(adjustl(mn))//'g'//&
     trim(adjustl(am))//'rm'//trim(adjustl(rat))//'_cut4.dat'
print*,fn
open(101,file=fn,status='replace',access='stream')
endif
if (check(4)==1) then
fn='px_n'//trim(adjustl(npd))//'M'//trim(adjustl(mn))//'g'//&
     trim(adjustl(am))//'rm'//trim(adjustl(rat))//'_cut4.dat'
print*,fn
open(102,file=fn,status='replace',access='stream')
endif
if (check(5)==1) then
fn='vy_n'//trim(adjustl(npd))//'M'//trim(adjustl(mn))//'g'//&
     trim(adjustl(am))//'rm'//trim(adjustl(rat))//'_cut4.dat'
print*,fn
open(111,file=fn,status='replace',access='stream')
endif
if (check(6)==1) then
fn='vz_n'//trim(adjustl(npd))//'M'//trim(adjustl(mn))//'g'//&
     trim(adjustl(am))//'rm'//trim(adjustl(rat))//'_cut4.dat'
print*,fn
open(112,file=fn,status='replace',access='stream')
endif
if (check(7)==1) then
fn='py_n'//trim(adjustl(npd))//'M'//trim(adjustl(mn))//'g'//&
     trim(adjustl(am))//'rm'//trim(adjustl(rat))//'_cut4.dat'
print*,fn
open(113,file=fn,status='replace',access='stream')
endif
if (check(8)==1) then
fn='pz_n'//trim(adjustl(npd))//'M'//trim(adjustl(mn))//'g'//&
     trim(adjustl(am))//'rm'//trim(adjustl(rat))//'_cut4.dat'
print*,fn
open(114,file=fn,status='replace',access='stream')
endif
if (check(9)==1) then

fn='vdiv_n'//trim(adjustl(npd))//'M'//trim(adjustl(mn))//'g'//&
     trim(adjustl(am))//'rm'//trim(adjustl(rat))//'_cut4.dat'
print*,fn
open(11,file=fn,status='replace',access='stream')
endif
if (check(10)==1) then
fn='vcrlz_n'//trim(adjustl(npd))//'M'//trim(adjustl(mn))//'g'//&
     trim(adjustl(am))//'rm'//trim(adjustl(rat))//'_cut4.dat'
print*,fn
open(12,file=fn,status='replace',access='stream')
endif
if (check(11)==1) then
fn='vcrly_n'//trim(adjustl(npd))//'M'//trim(adjustl(mn))//'g'//&
     trim(adjustl(am))//'rm'//trim(adjustl(rat))//'_cut4.dat'
print*,fn
open(13,file=fn,status='replace',access='stream')
endif
if (check(12)==1) then
fn='vcrlx_n'//trim(adjustl(npd))//'M'//trim(adjustl(mn))//'g'//&
     trim(adjustl(am))//'rm'//trim(adjustl(rat))//'_cut4.dat'
print*,fn
open(14,file=fn,status='replace',access='stream')
endif
if (check(13)==1) then
fn='df_n'//trim(adjustl(npd))//'M'//trim(adjustl(mn))//'g'//&
     trim(adjustl(am))//'rm'//trim(adjustl(rat))//'_cut4.dat'
print*,fn
open(15,file=fn,status='replace',access='stream')
endif
if (check(14)==1) then
fn='dragx_n'//trim(adjustl(npd))//'M'//trim(adjustl(mn))//'g'//&
     trim(adjustl(am))//'rm'//trim(adjustl(rat))//'_cut4.dat'
print*,fn
open(16,file=fn,status='replace',access='stream')
endif
if (check(15)==1) then
fn='dragy_n'//trim(adjustl(npd))//'M'//trim(adjustl(mn))//'g'//&
     trim(adjustl(am))//'rm'//trim(adjustl(rat))//'_cut4.dat'
print*,fn
open(17,file=fn,status='replace',access='stream')
endif
if (check(16)==1) then
fn='dragz_n'//trim(adjustl(npd))//'M'//trim(adjustl(mn))//'g'//&
     trim(adjustl(am))//'rm'//trim(adjustl(rat))//'_cut4.dat'
print*,fn
open(18,file=fn,status='replace',access='stream')
endif
if (check(17)==1) then
fn='entropy_n'//trim(adjustl(npd))//'M'//trim(adjustl(mn))//'g'//&
     trim(adjustl(am))//'rm'//trim(adjustl(rat))//'_cut4.dat'
print*,fn
open(19,file=fn,status='replace',access='stream')
endif
if (check(18)==1) then
fn='entropyall_n'//trim(adjustl(npd))//'M'//trim(adjustl(mn))//'g'//&
     trim(adjustl(am))//'rm'//trim(adjustl(rat))//'_cut4.dat'
print*,fn
open(20,file=fn,status='replace',access='stream')
endif
print*,'*************************************************************'
end subroutine printparams

subroutine closefile
close(11)
close(12)
close(13)
close(14)
close(15)
close(16)
close(17)
close(18)
close(19)
close(20)
close(99)
close(100)
close(101)
close(102)
close(111)
close(112)
close(113)
close(114)
end subroutine closefile
end module parameters
