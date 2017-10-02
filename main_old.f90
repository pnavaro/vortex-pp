program biot 
!Pierre Navaro
!navaro@unistra.fr
!IRMA http://www-irma.u-strasbg.fr

use sorties, only: gnuplot_output
use initial, only: lecture

implicit none

!$hmpp <myGPU> group, target=OPENCL
!$hmpp <myGPU> mapbyname, xp
!$hmpp <myGPU> mapbyname, yp
!$hmpp <myGPU> mapbyname, up
!$hmpp <myGPU> mapbyname, vp
integer, parameter :: idm = 100000
integer :: nbpart, i, j, k, n, iplot, istep, imov, nstep
real :: xp(idm), yp(idm), op(idm), up(idm), vp(idm)
real, dimension(64,64) :: w
real :: time

real(8) :: t0, t1, tcpu

!$hmpp <myGPU> resident
real :: delta
!$hmpp <myGPU> resident
real :: dt

call cpu_time(tcpu)
t0 = getRealTimer()
call lecture( nstep, imov, xp, yp, op, delta, idm, dt, nbpart )
do k = 1, nbpart
   write(9,*) xp(k), yp(k), op(k)
end do
call calcul_w(xp, yp, op, w, nbpart)  
do i = 1, 64
do j = 1, 64
   write(10,*) i, j, w(i,j)
end do
write(10,*) 
end do

iplot   = 0
time = 0.0
!$hmpp <myGPU> allocate
!$hmpp <myGPU> vitesse advancedload, args[up,vp,::delta]
!$hmpp <myGPU> deplace advancedload, args[::dt]

do istep = 1, nstep       !loop over time
   
   !up & vp not reloaded before the callsite and retrieved after
   !$hmpp <myGPU> vitesse callsite, args[up,vp].noupdate=true
   call vitesse(nbpart, xp, yp, op, up, vp, delta)
   !$hmpp <myGPU> deplace callsite, args[up,vp].noupdate=true
   call deplace(nbpart, xp, yp, up, vp, dt)
!   call calcul_w(xp, yp, op, w, nbpart)  
!   do i = 1, 64
!   do j = 1, 64
!      write(10,*) i, j, w(i,j)
!   end do
!      write(10,*) 
!   end do
  
   time = time + dt
   iplot = iplot + 1

   !call silo_output(istep, xp, yp, up, vp, op, w, nbpart)
!   call gnuplot_output(istep, iplot, time, xp, yp, op, nbpart, nstep)
   write(*,"(i5,1x,1g12.3)")istep, time

end do      !next time step

call cpu_time(tcpu)
t1 = getRealTimer()
write(*,"(5x,' CPU time = ', G15.3)") tcpu
write(*,"(5x,' CALLSITE time = ', G15.3)") t1-t0

!$hmpp <myGPU> vitesse delegatedstore, args[up,vp]
!$hmpp <myGPU> release

contains

!$hmpp <myGPU> vitesse codelet, args[nbpart,delta,op].const=true
subroutine vitesse (nbpart, xp, yp, op, up, vp, delta)
implicit none
integer, intent(in)  :: nbpart
real, intent(in)     :: xp(nbpart), yp(nbpart), op(nbpart)
real, intent(inout)  :: up(nbpart), vp(nbpart)
real, intent(in)     :: delta
real :: xo, yo, dx, dy, usum, vsum, dpi, d2
integer :: j, k

real(8) :: r2, r22, r2a1, r2a13, xm
real(8) :: a1, a12, a122

dpi   = 8.0 * atan( 1.0 )

a1    = 0.01
a12   = a1*a1
a122  = a12*a12
   
up = 0.
vp = 0.

do k = 1, nbpart
   
   usum = 0.0; vsum = 0.0
   xo = xp(k); yo = yp(k)
  
   do j = 1 , nbpart
      if( j .ne. k ) then
         dx    = xp( j ) - xo
         dy    = yp( j ) - yo
         r2    = dx * dx + dy * dy
         r22   = r2 * r2
         r2a1  = r2 + a12
         r2a13 = r2a1 * r2a1 * r2a1
         xm    = (r22+3.0*a12*r2+4.0*a122) / r2a13
         usum  = usum + dy * op(j) * xm
         vsum  = vsum - dx * op(j) * xm
      end if
   end do
   up(k) = usum / dpi
   vp(k) = vsum / dpi

   
end do
   
end subroutine vitesse

!$hmpp <myGPU> deplace codelet, args[nbpart,dt].const=true
subroutine deplace (nbpart, xp, yp, up, vp, dt)
implicit none
integer, intent(in)  :: nbpart
real, intent(inout)  :: xp(nbpart), yp(nbpart)
real, intent(in)     :: up(nbpart), vp(nbpart)
real, intent(in)     :: dt
integer :: k

do k = 1, nbpart
   xp(k) = xp(k) + dt * up(k)
   yp(k) = yp(k) + dt * vp(k)
end do
   
end subroutine deplace

subroutine calcul_w(xp, yp, op, w, nbpart)  
implicit none
integer, intent(in) :: nbpart
real, dimension(nbpart), intent(in) :: xp, yp, op
integer :: ip, jp, kp, i, j
real, parameter :: xmin    = - 2., xmax = 2.
real, parameter :: ymin    = - 2., ymax = 2.
real :: a1, a2, a3, a4, xt, yt
logical :: lflag
integer, parameter :: nxw=64, nyw=64
real, dimension(nxw,nyw), intent(inout) :: w

lflag   = .false.

do j = 1, nyw
do i = 1, nxw
w(i,j) = 0.0
end do
end do

do kp = 1,nbpart

   xt = (xp(kp)-xmin)/(xmax-xmin)*nxw
   yt = (yp(kp)-ymin)/(ymax-ymin)*nyw

   ip = floor(xt)
   jp = floor(yt)

   if (      (ip <= 0 .or. ip >= nxw) &
        .or. (jp <= 0 .or. jp >= nyw)) then

      lflag = .true.

   else

      a1 = (ip+1 - xt) * (jp+1 - yt)
      a2 = (xt   - ip) * (jp+1 - yt)
      a3 = (ip+1 - xt) * (yt   - jp)
      a4 = (xt   - ip) * (yt   - jp)

      w(ip  ,jp  )=w(ip  ,jp  )+a1*op(kp)
      w(ip+1,jp  )=w(ip+1,jp  )+a2*op(kp)
      w(ip  ,jp+1)=w(ip  ,jp+1)+a3*op(kp)
      w(ip+1,jp+1)=w(ip+1,jp+1)+a4*op(kp)

   end if

end do

!if (lflag) then
!
!   write(*,*) " xmin, xmax = ", xmin, xmax
!   write(*,*) " ymin, ymax = ", ymin, ymax
!   write(*,*) " xpmin, xpmax = ", minval(xp(:)), maxval(yp(:))
!   write(*,*) " ypmin, ypmax = ", minval(xp(:)), maxval(yp(:))
!   stop 'Calcul de la w. Revoir les parametres'
!
!end if

end subroutine calcul_w

function getRealTimer()
implicit none
real(8) :: out, getRealTimer
integer(8) :: count, count_rate, max
call system_clock(count, count_rate)
count = count - 1254348000*count_rate
out = count
out = out / count_rate
getRealTimer = out
end function getRealTimer

end program biot
