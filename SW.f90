program SW
  implicit none
  real*8, allocatable,dimension(:) :: a
  real*8::ti,ti1,w,ws,t,ts,Rxx,Rxy,Rzz,Rmax,dt,i
  integer :: j,k,ni,tmax,n
  ni = 16
  tmax = 30
  dt = 0.01
  !Increasing n leads to better integration via trapezium rule
  n = 13
  allocate(a(ni))
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
  ti1=ti(a,ni)
  w=0.5d+0
  ws=w*ti1
  open(1,file='Rxx,n=16')
  open(2,file='Rxy,n=16')
  open(3,file='Rzz,n=16')
  t = 0
  do while (t <= tmax)
    ts=t/ti1
    call trapz (ts,i,n,ws)
    write(1,*) t,Rxx(ws,ts,i)
    write(2,*) t,Rxy(ws,ts)
    write(3,*) t,Rzz(ws,ts,i)
    t = t + dt
  enddo
end program SW


function ti(a,ni)
  implicit none
  integer :: j,ni
  real*8, dimension(ni) :: a
  real*8 :: ti,sum,x,i
  sum = 0.0d+0
  i=0.5d+0
  do j = 1,ni
    x=(a(j)**2.0d+0)*i*(i+1.0d+0)
    sum=sum+x
  enddo
  ti = sqrt(6.0d+0/sum)
end function ti

function Rxx(w,t,i)
  implicit none
  real*8 :: Rxx,w,t,i
  Rxx=(w*(2.0d+0+exp(-t**2.0d+0)*((w**2.0d+0-2.0d+0)*cos(w*t)-2.0d+0*w*t*sin(w*t)))-4.0d+0*i)/(2.0d+0*w**3.0d+0)
end function Rxx

function Rxy(w,t)
	implicit none
  real*8 :: Rxy,w,t
  Rxy=exp(-t**2.0d+0)*(2.0d+0*w*t*cos(w*t)+(w**2.0d+0-2.0d+0)*sin(w*t))/(2.0d+0*w**2.0d+0)
end function Rxy


function Rzz(w,t,i)
  implicit none
  real*8 :: Rzz,w,t,i
  Rzz=(w*(w**2.0d+0+4.0d+0*exp(-t**2.0d+0)*cos(w*t)-4.0d+0)+8.0d+0*i)/(2.0d+0*w**3.0d+0)
end function Rzz



subroutine trapz(b,i,n,w)
  implicit none
  real*8 :: a,b,i,spac,z,x,sum,w,g
  integer:: it,n,j
  a = 0.0d+0
  it=2**(n-2)
  z=it
  spac=(b-a)/z
  x=a+0.5d+0*spac
  do j=1,it
    sum=sum+g(w,x)
  	x=x+spac
  enddo
  i=(b-a)*sum/z
end subroutine trapz


function g(ws,s)
	implicit none
  real*8 :: g,ws,s
	g=exp(-s**2.0d+0)*sin(ws*s)
end function g
