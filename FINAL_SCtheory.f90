program improvedSW
  implicit none
  real*8, dimension(3) :: wi,s,si
  real*8, allocatable,dimension(:) :: a
  real*8, allocatable,dimension(:,:) :: i
  real*8, allocatable,dimension(:) :: th,phi
  real*8, dimension(10000) :: mrxx,rxx,mrxy,rxy,mrzz,rzz
  real*8 :: dt,maxtime,ms,mi,sinth
  integer :: nucspins,n,k,f,kmax,l,fl
  nucspins=16
  allocate(a(nucspins))
  allocate(i(3,nucspins))
  allocate(phi(nucspins+1))
  allocate(th(nucspins+1))
  call random_seed()
  !Hyperfine constant values taken from paper in mT
  a(1)=-0.999985d+0
  a(2)=-0.7369246d+0
  a(3)=0.511210d+0
  a(4)=-0.0826998d+0
  a(5)=0.0655341d+0
  a(6)=-0.562082d+0
  a(7)=-0.905911d+0
  a(8)=0.357729d+0
  a(9)=0.358593d+0
  a(10)=0.869386d+0
  a(11)=-0.232996d+0
  a(12)=0.0388327d+0
  a(13)=0.661931d+0
  a(14)=-0.930856d+0
  a(15)=-0.893077d+0
  a(16)=0.0594001d+0
  !Applied magnetic field in mT. Where w=-yB
  wi=[0.0d+0,0.0d+0,0.5d+0]
  !Time steps in mT-1
  dt=0.1d+0
  !Total time interval in mT-1
  maxtime=30.0d+0
  kmax=maxtime/dt
  !Number of Monte Carlo integration steps
  fl=65536*64
  mrxx=0.0d+0
  mrxy=0.0d+0
  mrzz=0.0d+0
  ms=dsqrt(0.5d+0*(0.5d+0+1.0d+0))
  mi=dsqrt(0.5d+0*(0.5d+0+1.0d+0))
  do l=1,fl
    call randpointsphere(th,phi,nucspins)
    !Assign random angles to S(0)
    sinth = dsin(th(nucspins+1))
    s(1)=ms*sinth*dcos(phi(nucspins+1))
    s(2)=ms*sinth*dsin(phi(nucspins+1))
    s(3)=ms*dcos(th(nucspins+1))
    si=s
    do n=1,nucspins
      !Assign random angles to Ik(0)
      sinth = dsin(th(n))
      i(1,n)= mi*sinth*dcos(phi(n))
      i(2,n)= mi*sinth*dsin(phi(n))
      i(3,n)= mi*dcos(th(n))
    enddo
    do k=1,kmax+1
      !Summation for Monte Carlo integration
      mrxx(k)=mrxx(k) + si(1)*s(1)
      mrxy(k)=mrxy(k) + si(1)*s(2)
      mrzz(k)=mrzz(k) + si(3)*s(3)
      !Propagation routine
      call spropagate(nucspins,wi,dt,s,i,a)
      call ipropagate(nucspins,wi,dt,s,i,a)
      call spropagate(nucspins,wi,dt,s,i,a)
    enddo
  enddo
  !Divide by number of Monte Carlo steps to find integral solutions
  !Multiplication by 2 due to normalisation issue
  rxx=2.0d+0*mrxx/fl
  rxy=2.0d+0*mrxy/fl
  rzz=2.0d+0*mrzz/fl
  !Saves results to file
  call output (rxx,rxy,rzz,dt,maxtime)
  deallocate(a,i,th,phi)
end program improvedSW

!Outputs the tensor results into files for plotting
subroutine output(rxx,rxy,rzz,dt,maxtime)
  implicit none
  integer :: tmax,n
  real*8, dimension(10000) :: rxx,rxy,rzz
  real*8 :: dt,maxtime,t
  open(1,file='Rxx, n=16')
  open(2,file='Rxy, n=16')
  open(3,file='Rzz, n=16')
  t=0.0d+0
  tmax=maxtime/dt
  do n=1,tmax+1
    write(1,*)t,rxx(n)
    write(2,*)t,rxy(n)
    write(3,*)t,rzz(n)
    t=t+dt
  enddo
end subroutine output

!Propagates the electron spin angular momentum vectors
subroutine spropagate(nucspins,wi,dt,s,i,a)
  implicit none
  integer :: nucspins,b
  real*8, dimension(3) :: final,wi,s,s1,s2,s3,w
  real*8, dimension(3,nucspins) :: i
  real*8, dimension(nucspins) :: a
  real*8 :: dt,mw
  final=0.0d+0
  !w = wi + sum over nuclear spins of ak*Ik
  do b=1,nucspins
    final=final+a(b)*i(:,b)
  enddo
  w=wi+final
  !Normalise resultant w vector
  mw=norm2(w)
  w=w/mw
  s1=w*dot_product(w,s)
  s2=s-s1
  call cross(w,s,s3)
  s=s1+s2*cos(mw*dt/2.0d+0)+s3*sin(mw*dt/2.0d+0)
end subroutine spropagate

!Propagates the nuclear spin angular momentum vectors
subroutine ipropagate(nucspins,wi,dt,s,i,a)
  implicit none
  integer :: nucspins,k
  real*8, dimension(3) :: final,wi,s,i1,i2,i3,w
  real*8, dimension(3,nucspins) :: i
  real*8, dimension(nucspins) :: a
  real*8 :: dt,mw
  final=0.0d+0
  do k=1,nucspins
    w=a(k)*s
    mw=norm2(w)
    w=w/mw
    i1=w*dot_product(w,i(:,k))
    i2=i(:,k)-i1
    call cross(w,i(:,k),i3)
    i(:,k)=i1+i2*dcos(mw*dt)+i3*dsin(mw*dt)
  enddo
end subroutine ipropagate

!Calculates the cross product of a and b and returns c
subroutine cross(a,b,c)
  real*8, dimension(3) :: c,a,b
  c(1) = a(2)*b(3)-a(3)*b(2)
  c(2) = a(3)*b(1)-a(1)*b(3)
  c(3) = a(1)*b(2)-a(2)*b(1)
end subroutine cross

!Picks out random angles on a sphere - note that theta and phi are defined
!as in the paper but opposite to on wolfram alpha article on picking points on sphere
subroutine randpointsphere(th,phi,nucspins)
  implicit none
  real*8, dimension(nucspins+1) :: th,phi
  integer :: n,nucspins
  real*8 :: pi,u,v
  pi=4.D0*DATAN(1.D0)
  do n=1,nucspins+1
    call random_number(u)
    call random_number(v)
    phi(n)=2.0d+0*pi*u
    th(n)=dacos(2.0d+0*v-1.0d+0)
  enddo
end subroutine randpointsphere
