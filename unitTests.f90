program unitTests
use utilities, only: get_dist, getLocalIndices
use SphericalHarmonics
implicit none

!for get_dist
real, dimension(3) :: pos1, pos2
real :: dist = 0.0

!for getLocalIndices
real, dimension(3) :: centre
integer, dimension(3) :: localList
real, dimension(3,3) :: posns
real :: r_c = 3.0

!for legendre
real(8) :: P_l, x = 0.5
integer :: l

!get_dist test
if (.False.) then

	pos1 = (/ 1.0, 0.0, 0.0 /)
	pos2 = (/ 0.0, 1.0, 0.0 /)

	dist = get_dist(pos1,pos2)

	print*,dist

end if

!getLocalIndices test
if (.False.) then
	centre = (/ 0.0,0.0,0.0/)
	localList = (/0,0,0/)
	posns(1,:) = (/ 1.0, 0.0, 0.0 /)
	posns(2,:) = (/ 4.0, 0.0, 0.0 /)
	posns(3,:) = (/ 0.0, 2.0, 0.0 /)
	call getLocalIndices(posns,r_c, centre,localList)
	
	print*,localList
	

end if

!legendre test
if (.True.) then
	l = 5
	P_l = legendre(l,x)
	!print *,P_l
	
end if

end program unitTests
