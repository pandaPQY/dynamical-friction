!!!!!!!!!!!!change potential model: have to change  3 places !!!!!!!!!!!!!!!!!!
module parameters                            
implicit none
integer, parameter :: n=256,nx=n,ny=n,nz=n
real, dimension(5,nx,ny,nz) :: u
real, dimension(3,nx,ny,nz) :: g
real, dimension(nx,ny,nz) :: phi
real t,dt,tf
integer iter,i
real, parameter :: M=10.0
real, parameter :: Amp=1.0
integer, parameter :: ratio=20   !rm:5, rs:20, r:20                                              !!!!!!!here!!!!!!!!!!!!!!!!!!!
real, parameter :: rho0=1.
real vx,c,df(3)
real p,pp
character(5) npd,mn,am,rat
character(100) fn
logical check(13)

contains

subroutine init(u,nx,ny,nz,M,c,rho0)
  implicit none
  integer nx,ny,nz
  real u(5,nx,ny,nz)
  real p0,c,M,vx,rho0
  integer i
  real, parameter :: epsilon=.1
  real, parameter :: gamma=5./3.
  p0=3./5./2.
  u=0
  u(1,:,:,:)=rho0
  c=sqrt(gamma*p0/sum(u(1,:,:,:))*nx*ny*nz)
  vx=M*c
  u(2,:,:,:)=vx*rho0
  u(5,:,:,:)=p0*1.5
  u(5,:,:,:)=u(5,:,:,:)+sum(u(2:4,:,:,:)**2,1)/u(1,:,:,:)/2
  u(5,:,:,:)=u(5,:,:,:)+p0*u(2,:,:,:)*5./3./(2./3.)
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
  real Amp,radi
  real sigma
  sigma=real(nx)/real(ratio)
  phi=0
  !!$omp parallel do default(none) shared(g) private(i,j,k)
  do k=1,nz
     do j=1,ny
        do i=1,nx
           radi=sqrt((i-0.5-nx/2)**2+(j-0.5-ny/2)**2+(k-0.5-nz/2)**2)                              !here!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! phi(i,j,k)=-Amp*exp(-1.*radi**2/2/sigma**2)   ! gaussian   r           20
          ! phi(i,j,k)=Amp*-1./radi                       ! 1/r        nothing
           phi(i,j,k)=Amp*-1./radi*log(1.+radi/ratio)    ! galaxy     rs          20
          ! phi(i,j,k)=Amp*-1./(sqrt(ratio**2+radi**2))   ! plummer     rm           5
        enddo
     enddo
  enddo
  g(1,:,:,:)=-1./2*(cshift(phi,1,1)-cshift(phi,-1,1))
  g(2,:,:,:)=-1./2*(cshift(phi,1,2)-cshift(phi,-1,2))
  g(3,:,:,:)=-1./2*(cshift(phi,1,3)-cshift(phi,-1,3))
!  g(5,:,:,:)=sum(g(2:4,:,:,:)*u(2:4,:,:,:),1)
  write(*,*) 'sigma=',sigma
!  write(*,*),minval(g),maxval(g),minval(abs(g))
!  write(*,*) 1./sqrt(maxval(g)*sqrt(3.)),1./sqrt(minval(abs(g))*sqrt(3.))
!  phi=0
!  g=0
  !!$omp end parallel do
end subroutine initg

subroutine printparams
check= (/   1,    1,     1,    1,     1,    1,    1,     1,    1,      1,      1,       1,       1    /)
!          rho    me     vx    px     vy    vz    py     pz   vdiv    vcrl    vcrly     vcrlx    df 
!           1     2      3     4      5     6     7      8     9       10      11       12       13
check=0
! B field is stored on the left side of each cell
write(npd,'(i3)') n
write(mn,'(f4.1)') M
write(am,'(f4.1)') Amp
write(rat,'(i3)') ratio
!call init(u,nx,ny,nz,M,c,rho0)
vx=M*c
print*,'**************************************************************'
print*,'parameters:                          *************************'
print*,'**************************************************************'
print*,'n    =',n
print*,'Amp  =',Amp
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
print*,'max potential=',maxval(phi)
print*,'min potential=',minval(phi)
print*,'*************************************************************'
tf=nx*4
t=0
iter=0
print*,'writing in :                      **************************'                               !!!!!!!!here!!!!!!!!!!!!!!!
if (check(1)==1) then
fn='lgrho_n'//trim(adjustl(npd))//'M'//trim(adjustl(mn))//'g'//&
     trim(adjustl(am))//'rs'//trim(adjustl(rat))//'.dat'
print*,fn
open(99,file=fn,status='replace',access='stream')
endif
if (check(2)==1) then
fn='me_n'//trim(adjustl(npd))//'M'//trim(adjustl(mn))//'g'//&
     trim(adjustl(am))//'rs'//trim(adjustl(rat))//'.dat'
print*,fn
open(100,file=fn,status='replace',access='stream')
endif
if (check(3)==1) then
fn='vx_n'//trim(adjustl(npd))//'M'//trim(adjustl(mn))//'g'//&
     trim(adjustl(am))//'rs'//trim(adjustl(rat))//'.dat'
print*,fn
open(101,file=fn,status='replace',access='stream')
endif
if (check(4)==1) then
fn='px_n'//trim(adjustl(npd))//'M'//trim(adjustl(mn))//'g'//&
     trim(adjustl(am))//'rs'//trim(adjustl(rat))//'.dat'
print*,fn
open(102,file=fn,status='replace',access='stream')
endif
if (check(5)==1) then
fn='vy_n'//trim(adjustl(npd))//'M'//trim(adjustl(mn))//'g'//&
     trim(adjustl(am))//'rs'//trim(adjustl(rat))//'.dat'
print*,fn
open(111,file=fn,status='replace',access='stream')
endif
if (check(6)==1) then
fn='vz_n'//trim(adjustl(npd))//'M'//trim(adjustl(mn))//'g'//&
     trim(adjustl(am))//'rs'//trim(adjustl(rat))//'.dat'
print*,fn
open(112,file=fn,status='replace',access='stream')
endif
if (check(7)==1) then
fn='py_n'//trim(adjustl(npd))//'M'//trim(adjustl(mn))//'g'//&
     trim(adjustl(am))//'rs'//trim(adjustl(rat))//'.dat'
print*,fn
open(113,file=fn,status='replace',access='stream')
endif
if (check(8)==1) then
fn='pz_n'//trim(adjustl(npd))//'M'//trim(adjustl(mn))//'g'//&
     trim(adjustl(am))//'rs'//trim(adjustl(rat))//'.dat'
print*,fn
open(114,file=fn,status='replace',access='stream')
endif
if (check(9)==1) then

fn='vdiv_n'//trim(adjustl(npd))//'M'//trim(adjustl(mn))//'g'//&
     trim(adjustl(am))//'rs'//trim(adjustl(rat))//'.dat'
print*,fn
open(11,file=fn,status='replace',access='stream')
endif
if (check(10)==1) then
fn='vcrl_n'//trim(adjustl(npd))//'M'//trim(adjustl(mn))//'g'//&
     trim(adjustl(am))//'rs'//trim(adjustl(rat))//'.dat'
print*,fn
open(12,file=fn,status='replace',access='stream')
endif
if (check(11)==1) then
fn='vcrly_n'//trim(adjustl(npd))//'M'//trim(adjustl(mn))//'g'//&
     trim(adjustl(am))//'rs'//trim(adjustl(rat))//'.dat'
print*,fn
open(13,file=fn,status='replace',access='stream')
endif
if (check(12)==1) then
fn='vcrlx_n'//trim(adjustl(npd))//'M'//trim(adjustl(mn))//'g'//&
     trim(adjustl(am))//'rs'//trim(adjustl(rat))//'.dat'
print*,fn
open(14,file=fn,status='replace',access='stream')
endif
if (check(13)==1) then
fn='df_n'//trim(adjustl(npd))//'M'//trim(adjustl(mn))//'g'//&
     trim(adjustl(am))//'rs'//trim(adjustl(rat))//'.dat'
print*,fn
open(15,file=fn,status='replace',access='stream')
endif
print*,'*************************************************************'
end subroutine printparams

subroutine closefile
close(11)
close(12)
close(13)
close(14)
close(15)
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
