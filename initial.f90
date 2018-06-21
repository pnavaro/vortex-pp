module initial

implicit none

real :: pi, gam0

integer :: nr, nproc, nray=6
integer, private :: i, j, k

logical :: gauss = .true.

contains

subroutine lecture( nstep, imov, xp, yp, op, delta, idm, dt, nbpart ) 

namelist/donnees/nstep, dt, imov, amach, nray, delta

integer, parameter :: ix=0, jx=200, iy=0, jy=200
integer :: nstep, imov, idm, nbpart

real :: xp( idm ), yp( idm ), op( idm )
real :: rf( idm ), zf( idm ), gam( idm )

real :: circ, al, ur, tau, aom, u0, r0
real :: amach, gomeg, delta, dt

nstep = 1
dt    = 0.01
imov  = 1
amach = 0.1
nproc = 1
nray  = 6
delta = 0.01

open(10, file = "input" )
read(10,donnees)
close(10)

pi = 4. * atan(1.)
r0 = 0.5 
u0 = amach 

gam0   = u0 * 2.0 * pi / 0.7 * r0	!gaussienne
!gam0   = 2. * pi * r0 * u0		!constant
!gam0   = 2. * pi / 10.0

aom    = gam0 / ( pi * r0**2 )	!Amplitude du vortex
tau    = 8.0 * pi**2 / gam0	!Periode de co-rotation
gomeg  = gam0/ (4.0*pi)		!Vitesse angulaire
ur     = gomeg 			!Vitesse tangentielle du vortex
al     = 0.5 * tau 	 	!Longeur d'onde

write(*,*) " iterations : ", nstep
write(*,*) " pas de temps : ", dt
write(*,*) " animation : ", imov, " steps "
write(*,*) " aom = ", aom
write(*,*) " r0 = ", r0
write(*,*) " Circulation = ", gam0
write(*,*) " Vitesse de rotation gomeg = ", gomeg
write(*,*) " ur = ", ur
write(*,*) " periode de corotation = ", tau
write(*,*) " --------------------------------------------- "

call distrib( rf, zf, gam, r0, idm )

do k = 1, nr
   xp( k    ) = rf( k ) 
   yp( k    ) = zf( k ) + 1.
   op( k    ) = gam( k )
   xp( k+nr ) = rf( k ) 
   yp( k+nr ) = zf( k ) - 1.
   op( k+nr ) = gam( k )
end do
nbpart = 2 * nr

!Calcul de la vitesse de propagation du systeme

circ = sum(op(1:nr))

write(*,*) ' Nombre total de particules =',nbpart

end subroutine lecture

!---------------------------------------------------------------

subroutine distrib(rf, zf, cir, ray, idm )

integer :: kd, nsec, nsec0, idm

real :: rf( * ), zf( * ), cir( * ), ds( idm )
real :: ssurf, q, sigma, teta, dss, r1, r2, s1, s2, eps
real :: gamt, sgam, dteta, surf, ray, dray, r, dr

pi    = 4.0 * atan( 1.0 )
nsec  = 6

!     rf,zf : position de la particule
!     ds    : taille de l'element de surface
!     cir   : circulation de la particule
!     dr    : pas d'espace dans la direction radiale 
!     nray  : nb de pts ds la direction radiale.
!     dray  : rayon de la particule placee au centre
!     ray   : rayon de la section
!     gam0  : circulation totale 
!     surf  : surface de la section
!     nsec  : nombre de points dans la premiere couronne

dr      = ray / ( nray + 0.5 )
dray    = 0.5 * dr	!rayon de la section centrale
surf    = pi * ray * ray 
dteta   = 2.0 * pi / float( nsec)
gamt    = gam0 / ( 1. - exp( -1.0 ) ) ! total strength of gaussian vortex

k       = 1
rf(  1) = 0.0
zf(  1) = 0.0
ds(  1) = pi * dray * dray

if ( gauss ) then
   cir( 1) = gamt * ( 1.-exp(-(dray/ray)**2)) !gauss
else
   cir( 1) = gam0 * ds( 1 ) / surf   !uniforme
end if
sgam    = cir( 1 )

r1    = dray
s1    = pi * r1**2
nsec0 = nsec
nsec  = 0

!cpn   *** parametre de l'ellipse ***
!c      eps = 0.01		! 0.0 --> disque
      eps = 0.0
!cpn   ******************************

do i = 1, nray

   nsec  = nsec + nsec0
   dteta = 2.0 * pi / float(nsec)
   r     = float( i ) * dr 

   r2  = r + 0.5 * dr
   s2  = pi * r2**2 
   dss = s2 - s1
   s1  = s2

   do j = 1, nsec

      k = k + 1
      if( k .gt. idm ) stop ' !! idm < Nr !! '

      teta    = float( j ) * dteta 
      sigma   = r * ( 1.0 + eps * cos( 2.0*teta ) )
      rf( k ) = sigma * cos( teta )
      zf( k ) = sigma * sin( teta )

      ds(  k ) = dss / float( nsec )

      if ( gauss ) then
         q       = 0.5 * (exp(-(r1/ray)**2)-exp(-(r2/ray)**2) )
         cir( k ) = gamt * dteta / pi * q	! gauss
      else
         cir( k ) = gam0 * ds( k ) / surf	! uniforme
      end if

      sgam     = sgam + cir( k )

    end do

    r1  = r2 

    kd = k - nsec + 1 

end do

nr = k

ssurf = 0.0
do i = 1, nr
   ssurf = ssurf + ds(i)
end do

write(*,*) 'Nb de pts sur la section :', nr,'(',idm,' max )'
write(*,*) 'surface theorique - pratique :', (surf),' - ',(ssurf)
if (gauss) then
   write(*,*) 'circulation theorique - pratique :',(gam0),' ; ',(sgam)
else
   write(*,*) 'circulation theorique - pratique :',(gam0),' ; ',(sgam)
end if

1000  format( F11.5, 2X, F11.5, 2X, F11.5, 1X, F11.5 )

end subroutine distrib

end module initial
