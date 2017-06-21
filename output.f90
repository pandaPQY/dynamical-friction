module output
contains
subroutine me(u,nx,ny,nz,t)
   implicit none
   integer nx,ny,nz
   real u(5,nx,ny,nz)
   real mom(3), ene,t
   mom(1)=sum(u(2,:,:,:))
   mom(2)=sum(u(3,:,:,:))
   mom(3)=sum(u(4,:,:,:))
   ene=sum(u(2:4,:,:,:)**2/spread(u(1,:,:,:),1,3))
   print*, 'mom:',mom(1),'ene:',ene
   write(100) t,mom,ene
end subroutine me

subroutine writeu(u,n)
  implicit none
  integer n
  real u(5,n,n,n)
  write(99) log(u(1,:,n/2,:))
end subroutine writeu

subroutine writev(u,n,vx)
  implicit none
  integer n
  real u(5,n,n,n)
  real dv(n,n),dv2(n,n),vx,dv1(n,n)
   dv1=u(2,:,n/2,:)/u(1,:,n/2,:)-vx !x
   dv=u(3,:,n/2,:)/u(1,:,n/2,:)     !y
   dv2=u(4,:,n/2,:)/u(1,:,n/2,:)    !z
!   dv1=1./2*(u(2,:,n/2,:)/u(1,:,n/2,:)+u(2,:,n/2+1,:)/u(1,:,n/2+1,:))-vx !x
!   dv=1./2*(u(3,:,n/2,:)/u(1,:,n/2,:)+u(3,:,n/2+1,:)/u(1,:,n/2+1,:))     !y
!   dv2=1./2*(u(4,:,n/2,:)/u(1,:,n/2,:)+u(4,:,n/2+1,:)/u(1,:,n/2+1,:))    !z
  write(101) dv1
  write(111) dv
  write(112) dv2
end subroutine writev

subroutine writep(u,n,vx,rho0)
  implicit none
  integer n
  real u(5,n,n,n)
  real dp(n,n),dp2(n,n),vx,rho0,dp1(n,n)
  dp1=u(2,:,n/2,:)-vx*rho0
  write(102) dp1
   dp=u(3,:,n/2,:)
   dp2=u(4,:,n/2,:)
   write(113) dp
   write(114) dp2
end subroutine writep

subroutine vcrl(u1,n)
integer n,k,kp,j,jp,i,ip
real u1(5,n,n,n),u2(5,n,n,n),u(5,n,2,n),curlx,curly,curlz,vzx,vyx,vxy,vzy,vxz,vyz,divx,divy,divz
j=n/2
jp=j+1
u2(2,:,:,:)=u1(2,:,:,:)/u1(1,:,:,:)
u2(3,:,:,:)=u1(3,:,:,:)/u1(1,:,:,:)
u2(4,:,:,:)=u1(4,:,:,:)/u1(1,:,:,:)
u=u2(:,:,j:jp,:)
j=1
jp=j+1
do k=1,n
   kp=mod(k,n)+1
   do i=1,n
      ip=mod(i,n)+1
      divx=(u(2,ip,j,k)-u(2,i,j,k)+u(2,ip,jp,k)-u(2,i,jp,k) + &
          u(2,ip,j,kp)-u(2,i,j,kp)+u(2,ip,jp,kp)-u(2,i,jp,kp) )/4
      divy=(u(3,i,jp,k)-u(3,i,j,k)+u(3,ip,jp,k)-u(3,ip,j,k) + &
          u(3,i,jp,kp)-u(3,i,j,kp)+u(3,ip,jp,kp)-u(3,ip,j,kp) )/4
      divz=(u(4,i,j,kp)-u(4,i,j,k)+u(4,ip,j,kp)-u(4,ip,j,k) + &
          u(4,i,jp,kp)-u(4,i,jp,k)+u(4,ip,jp,kp)-u(4,ip,jp,k) )/4
      write(11) divx+divy+divz
      vzx=(u(4,ip,j,k)-u(4,i,j,k)+u(4,ip,jp,k)-u(4,i,jp,k) + &
          u(4,ip,j,kp)-u(4,i,j,kp)+u(4,ip,jp,kp)-u(4,i,jp,kp) )/4
      vyx=(u(3,ip,j,k)-u(3,i,j,k)+u(3,ip,jp,k)-u(3,i,jp,k) + &
          u(3,ip,j,kp)-u(3,i,j,kp)+u(3,ip,jp,kp)-u(3,i,jp,kp) )/4
      vxy=(u(2,i,jp,k)-u(2,i,j,k)+u(2,ip,jp,k)-u(2,ip,j,k) + &
          u(2,i,jp,kp)-u(2,i,j,kp)+u(2,ip,jp,kp)-u(2,ip,j,kp) )/4
      vzy=(u(4,i,jp,k)-u(4,i,j,k)+u(4,ip,jp,k)-u(4,ip,j,k) + &
          u(4,i,jp,kp)-u(4,i,j,kp)+u(4,ip,jp,kp)-u(4,ip,j,kp) )/4
      vxz=(u(2,i,j,kp)-u(2,i,j,k)+u(2,ip,j,kp)-u(2,ip,j,k) + &
          u(2,i,jp,kp)-u(2,i,jp,k)+u(2,ip,jp,kp)-u(2,ip,jp,k) )/4
      vyz=(u(3,i,j,kp)-u(3,i,j,k)+u(3,ip,j,kp)-u(3,ip,j,k) + &
          u(3,i,jp,kp)-u(3,i,jp,k)+u(3,ip,jp,kp)-u(3,ip,jp,k) )/4
      curlx=vzy-vyz
      curly=vxz-vzx
      curlz=vyx-vxy
      write(12) curlz!**2+curly**2+curlz**2
      write(13) curly
      write(14) curlx
   enddo
  enddo

end subroutine vcrl

subroutine slice(u,n)
  implicit none
  integer n
  real u(5,n,n,n)
  open(77,file='u_n256M0.3g0.1r20.dat',status='replace',access='stream')
  write(77) u(1:4,:,:,:)
  close(77)
end subroutine slice

subroutine force(u,g,n,df)
  implicit none
  integer n
  real u(5,n,n,n),g(3,n,n,n),df(3)
  df(1)=sum(g(1,:,:,:)*u(1,:,:,:))
  df(2)=sum(g(2,:,:,:)*u(1,:,:,:))
  df(3)=sum(g(3,:,:,:)*u(1,:,:,:))
  write(15) df
end subroutine force

subroutine dragforce(u,g,n,drag)
  implicit none
  integer n
  real u(5,n,n,n),g(3,n,n,n),drag(3,n,n,n)
  drag(1,:,:,:)=g(1,:,:,:)*u(1,:,:,:)
  drag(2,:,:,:)=g(2,:,:,:)*u(1,:,:,:)
  drag(3,:,:,:)=g(3,:,:,:)*u(1,:,:,:)
  write(16) drag(1,:,n/2,:)
  write(17) drag(2,:,n/2,:)
  write(18) drag(3,:,n/2,:)
end subroutine dragforce

subroutine entropy(u,n)
  implicit none
  integer n
  real u(5,n,n,n)
  real entro(n,n,n)
  real p(n,n,n)
  real,parameter :: gamma=5./3.
  p=(gamma-1.)*(u(5,:,:,:)-1./2*sum(u(2:4,:,:,:),1)/u(1,:,:,:))
  entro=p/u(1,:,:,:)**gamma
  write(19) entro(:,n/2,:)
  write(20) sum(entro)
end subroutine entropy
end module output
