module sorties

implicit none
real, parameter :: pi = 3.141592654
integer, private :: i, j, k
integer, private :: iplot1, iplot2
integer, parameter :: ndims=2
integer, dimension(ndims) :: dims
integer :: dbfile


contains

subroutine gnuplot_output( istep, iplot, time, xp, yp, op, nbpart, nstep )

integer :: nbpart, nstep, iplot, istep
real :: xp(*), yp(*), op(*)
character(len=40):: nom
integer :: ifich, ifich1, ifich2
real :: time
integer :: kk0, kk1, kk2, kk3, kk4
character(len=04):: fin
character(len=01):: aa,bb,cc,dd

ifich = 80

nom ='/tmp/part'

kk0 = iplot
kk1 = kk0/1000
aa  = char(kk1 + 48)
kk2 = (kk0 - kk1*1000)/100
bb  = char(kk2 + 48)
kk3 = (kk0 - (kk1*1000) - (kk2*100))/10
cc  = char(kk3 + 48)
kk4 = (kk0 - (kk1*1000) - (kk2*100) - (kk3*10))/1
dd  = char(kk4 + 48)
fin = aa//bb//cc//dd

open(  ifich, file = 'part.gnu', position="append" )
if ( iplot == 1 ) then
   rewind(ifich)
   write( ifich, "('set xr[-2:2]')" )
   write( ifich, "('set yr[-2:2]')" )
end if
write( ifich, 1000 )trim(nom)//fin
close(ifich)

open(ifich, file = trim(nom)//fin )
rewind(ifich)
do k = 1, nbpart
   write(ifich,"(3f10.3)") xp(k), yp(k), op(k)
end do
close(ifich )

1000 format( "plo '",a,"' w p ")

end subroutine gnuplot_output

end module sorties
