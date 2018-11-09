program lens
  implicit none
  real*8, dimension(3,500000) :: in,as
  real*8 :: a,rmax,r,theta,psi,pi,y,perm,gammae,gammain,gammaas
  integer :: maxin,maxas,n
  real*8, allocatable,dimension(:) :: hypin,hypas
  a=0.60583d+0
  rmax=7.0d+0
  pi = 4.0d+0*atan(1.0d+0)
  perm=pi*4.0d-7
  gammae=-1.760859644d+11
  gammain=5.8908d+7
  gammaas=4.595d+7
  call blende(a,in,maxin,as,maxas)
  allocate(hypin(maxin))
  allocate(hypas(maxas))
  open(2,file="hyperfine")
  do n=1,maxin
    r=sqrt(in(1,n)*in(1,n)+in(2,n)*in(2,n)+in(3,n)*in(3,n))
    theta = acos(in(3,n)/r)
    y = psi(r,theta,rmax)
    hypin(n) = -(2.0d+0/3.0d+0)*perm*gammae*gammain*y*y*1.0d+27

  enddo
  do n=1,maxas
    r=dsqrt(as(1,n)*as(1,n)+as(2,n)*as(2,n)+as(3,n)*as(3,n))
    theta = dacos(as(3,n)/r)
    y = psi(r,theta,rmax)
    hypas(n) = -(2.0d+0/3.0d+0)*perm*gammae*gammaas*y*y*1.0d+27
    write(2,*)hypas(n)
  enddo
end program lens

subroutine blende(a,in,k,as,s)
  implicit none
  real*8, dimension(3,500000) :: p,h,in,as
  real*8 :: x,y,z,a,r
  integer :: n,j,s,k
  n=1
  z=0.0d+0
  do while (z<=14.0d+0)
    x=0.0d+0
    do while (x<=14.0d+0)
      y=0.0d+0
      do while (y<=14.0d+0)
        p(:,n)=[x,y,z]
        h(:,n)=[x + a/4.0d+0, y + a/4.0d+0, z + a/4.0d+0]
        n=n+1
        p(:,n)=[x, y + a/2.0d+0, z + a/2.0d+0]
        h(:,n)=[x + a/4.0d+0, y + a*0.75d+0, z + a*0.75d+0]
        n=n+1
        p(:,n)=[x + a/2.0d+0, y, z + a/2.0d+0]
        h(:,n)=[x + a*0.75d+0, y + a/4.0d+0, z + a*0.75d+0]
        n=n+1
        p(:,n)=[x + a/2.0d+0, y + a/2.0d+0, z]
        h(:,n)=[x + a*0.75d+0, y + a*0.75d+0, z + a/4.0d+0]
        n=n+1
        y=y+a
      enddo
      x=x+a
    enddo
    z=z+a
  enddo
  s=1
  k=1
  open(1,file="hemisphere")
  do j= 1,n-1
    p(:,j) = p(:,j) - 7.0d+0
    h(:,j) = h(:,j) - 7.0d+0
    r=sqrt(p(1,j)*p(1,j)+p(2,j)*p(2,j)+p(3,j)*p(3,j))
    if ( r<=7.0d+0 .and. p(3,j)>=0.0d+0 ) then
      in(1,k)=p(1,j)
      in(2,k)=p(2,j)
      in(3,k)=p(3,j)
      k=k+1
    end if
    r=sqrt(h(1,j)*h(1,j)+h(2,j)*h(2,j)+h(3,j)*h(3,j))
    if ( r<=7.0d+0 .and. h(3,j)>=0.0d+0 ) then
      as(1,s)=h(1,j)
      as(2,s)=h(2,j)
      as(3,s)=h(3,j)
      write(1,*)as(:,s)
      s=s+1
    end if
  enddo
  k=k-1
  s=s-1
end subroutine blende


function psi(r,theta,rmax)
  implicit none
  real*8 :: r,theta,rmax,pi,bigpsi,bigtheta,bigr,psi
  pi = 4.0d+0*atan(1.0d+0)
  bigpsi = 1.0d+0/sqrt(2*pi)
  bigtheta = sqrt(3.0d+0/2.0d+0)*cos(theta)
  bigr = sqrt(2.0d+0/rmax)*sin(r*pi/rmax)/r
  psi = bigpsi*bigtheta*bigr
end function psi
