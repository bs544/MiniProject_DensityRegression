program unitTests
use utilities
use BasisFunctions
implicit none

!for get_dist
real(8), dimension(3) :: pos1, pos2
real(8) :: dist = 0.0d0

!for getLocalIndices and checkInCell
real(8), dimension(3) :: centre
integer, dimension(3) :: localList
real(8), dimension(3,3) :: posns, polarPosns
real(8) :: r_c = 3.0
real(8), dimension(3,3) :: cell
logical :: check
integer :: ii
type(pointType) :: pointInfo

!for legendre
real(8) :: P_l, x = 0.5
integer :: l

!get_dist test
if (.True.) then

	pos1 = (/ 1.0, 0.0, 0.0 /)
	pos2 = (/ 0.0, 1.0, 0.0 /)

	dist = getDist(pos1,pos2)

	print*,dist

end if


!getLocalIndices test
if (.False.) then
	centre = (/ 0.0,0.0,0.0/)
	localList = (/0,0,0/)
	posns(:,1) = (/ 1.0, 0.0, 0.0 /)
	posns(:,2) = (/ 4.0, 0.0, 0.0 /)
	posns(:,3) = (/ 0.0, 2.0, 0.0 /)
	call getLocalIndices(posns,r_c, centre,localList,cell)
	
	print*,localList

end if

!checkInCell test
if(.False.) then
	cell(:,1) = (/ 1.0,0.0,0.0 /)
	cell(:,2) = (/ 0.0,1.5,1.5 /)
	cell(:,3) = (/ 0.0,1.5,-1.5 /)
	
	posns(:,1) = (/0.0,0.0,0.0/)
	posns(:,2) = (/2.5,0.0,0.0/)
	posns(:,3) = (/0.0,0.75,0.9/)
	
	do ii = 1,3
		check = checkInCell(posns(:,ii),cell)
		print*,check
	end do
end if

!dot test
if (.False.) then
	posns(:,1) = (/0.0,0.0,0.1/)
	posns(:,2) = (/2.5,0.0,0.0/)
	posns(:,3) = (/0.0,0.75,0.9/)
	
	print*,dot(posns(:,1),posns(:,2))
	print*,dot(posns(:,1),posns(:,3))
end if
	
!localCart2Polar test
if (.False.) then
	posns(:,1) = (/0.0,0.0,0.1/)
	posns(:,2) = (/2.5,0.0,0.0/)
	posns(:,3) = (/0.0,0.75,0.9/)
	
	centre = (/ 0.0,0.0,0.0 /)
	localList = (/ 1,2,0 /)
	
	polarPosns = 0
	
	
	
	call localCart2Polar(polarPosns, point)
	
	print*, polarPosns(:,1)
	print*, posns(:,1)
	
end if	

!factorial test
if (.False.) then
	do ii = 1, 7

		print *, factorial(ii)
	end do
end if


!minAbsFloor test
if (.True.) then
	dist = 1.8
	print*, minAbsFloor(dist)
	dist = -1.6
	print*, minAbsFloor(dist)
	dist = 102.4
	print*, minAbsFloor(dist)
end if


!legendre test
if (.False.) then
	l = 5
	P_l = legendre(l,x)
	!print *,P_l
	
end if

end program unitTests
