module workfile
contains
subroutine fluidx(u,nx,ny,nz,dt)
  implicit none
  integer nx,ny,nz
  real u(5,nx,ny,nz),dt
  real, dimension(5,nx) :: u1x
  integer j,k

  !$omp parallel do private(j,u1x)
  do k=1,nz
     do j=1,ny
!        jp=mod(j,ny)+1
!        kp=mod(k,nz)+1
!        u1x=u(:,:,j,k)
        call tvd1(u(:,:,j,k),nx,dt)
     end do
  end do
end subroutine fluidx

subroutine gravity(g,u,n,dt)
  implicit none
  integer n
  real,dimension(5,n,n,n) :: u,g5
  real,dimension(3,n,n,n) :: g
  real dt
  g5=0
  g5(2,:,:,:)=g(1,:,:,:)*u(1,:,:,:)
  g5(3,:,:,:)=g(2,:,:,:)*u(1,:,:,:)
  g5(4,:,:,:)=g(3,:,:,:)*u(1,:,:,:)
  g5(5,:,:,:)=sum(g*u(2:4,:,:,:),1)
  u=u+g5*dt
end subroutine gravity

recursive subroutine tvd1(u,n,dt)
  implicit none
  integer n
  real dt
  real u(5,n)
!locals    
  real, dimension(5,n) :: v, u1,wr,wl,fr,fl,flux,dfrp,dfrm&
       ,dflm,dflp,dfl,dfr
  real c
!  print*, maxval(u(1,:)),minval(u(1,:))
  call mhdflux(v,c,u,n)
  wr=u+v
  wl=u-v
  fr=c*wr
  fl=cshift(c*wl,1,2)
  flux=(fr-fl)/2
  u1=u-(flux-cshift(flux,-1,2))*dt/2
  call mhdflux(v,c,u1,n)
  wr=u1+v
  wl=u1-v
  fr=c*wr
  dfrp=(cshift(fr,1,2)-fr)/2
  dfrm=(fr-cshift(fr,-1,2))/2
  dfr=0
  where(dfrp*dfrm>0)
     dfr=2*dfrp*dfrm/(dfrp+dfrm)  !vanleer
  end where
  fl=cshift(c*wl,1,2)
  dflp=(fl-cshift(fl,1,2))/2
  dflm=(cshift(fl,-1,2)-fl)/2
  dfl=0
  where(dflp*dflm>0)
     dfl=2*dflp*dflm/(dflp+dflm)
  end where
  flux=(fr-fl+(dfr-dfl))/2
  u=u-(flux-cshift(flux,-1,2))*dt
!  print*, maxval(u(1,:)),minval(u(1,:))
  return
end subroutine tvd1

recursive subroutine mhdflux(v,c,u,n)
  implicit none
  integer n
  real, dimension(5,n)::v,u
  real, intent(out) :: c
  ! locals
  real, dimension(n) :: vx,ps,p
  real gamma
  
  gamma=5./3.
!  print*, maxval(u(1,:)),minval(u(1,:))
  vx=u(2,:)/u(1,:)
  ps=(u(5,:)-sum(u(2:4,:)**2,1)/u(1,:)/2)*(gamma-1)
  v(1,:)=u(2,:)
  v(2,:)=u(2,:)*vx+ps
  v(3,:)=u(3,:)*vx
  v(4,:)=u(4,:)*vx
  v(5,:)=(u(5,:)+ps)*vx
  p=ps
  c=maxval(abs(vx)+sqrt(abs(gamma*p)/u(1,:)))
  if (c>0) v=v/c
!  print*, maxval(u(1,:)),minval(u(1,:))
end subroutine mhdflux

subroutine transpose12(ut,u,nx,ny,nz)
  implicit none
  integer nx,ny,nz
  real u(5,nx,ny,nz)
  real ut(5,ny,nx,nz)

  integer i,j,k
  real u2(5,ny,nx)

!$omp parallel do default(none) shared(u,ut,ny,nx,nz) private(i,j,u2)
  do k=1,nz
     do j=1,ny
        do i=1,nx
           u2(:,j,i)=u((/1,3,2,4,5/),i,j,k)
        end do
     end do
     ut(:,:,:,k)=u2
  end do
end subroutine transpose12


subroutine transpose13(ut,u,nx,ny,nz)
  implicit none
  integer nx,ny,nz
  real u(5,nx,ny,nz)
  real ut(5,nz,ny,nx)

  integer i,j,k
  real u2(5,nz,nx)

  if (nx .eq. nz) then
!$omp parallel do default(none) shared(u,ut,ny,nx,nz) private(i,k,u2)
  do j=1,ny
     do k=1,nz
        do i=1,nx
           u2(:,k,i)=u((/1,4,3,2,5/),i,j,k)
        end do
     end do
     ut(:,:,j,:)=u2
  end do
  else if (nz .eq. 1) then
     call transpose12(ut,u,nx,ny,nz)
     do i=1,nx
        do j=1,ny
           ut(:,1,j,i)=ut((/1,4,2,3,5/),1,j,i)
        end do
     end do
  else if (nx .eq. 1) then
     call transpose12(ut,u,ny,nz,nx)
     do i=1,ny
        do j=1,nz
           ut(:,j,i,1)=ut((/1,4,2,3,5/),j,i,1)
        end do
     end do
  else
     write(*,*) 'nz<>nx not supported'
     stop
  endif
end subroutine transpose13
end module workfile
