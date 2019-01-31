module utilities
	implicit none
	
	public getDist
	public getLocalIndices
	public checkInCell
	public localCart2Polar
	public factorial
	public minAbsFloor

	contains
	
	
	
	!gets distance between two cartesian coordinates
	function getDist(pos1, pos2)
		implicit none
		real(8), intent(in) :: pos1(3)
		real(8), intent(in) :: pos2(3)
		real(8) :: getDist
		
		getDist = sqrt(sum((pos1-pos2)**2))
		
	end function getDist



	subroutine getLocalIndices(posns,r_c, centre,localList,cell)
	!gets an array of indices for the positions in a cutoff sphere around the centre
		implicit none
		real(8), intent(in) :: r_c
		real(8), intent(in) :: posns(:,:)
		real(8), intent(in) :: centre(3)
		real(8), intent(in) :: cell(3,3)!vector elements are given by the first index, the vector labels by the second (first index varies faster in fortran)
		integer :: localList(:)
		logical :: checkNear
		integer :: array_shape(2)
		integer :: idx, natoms, counter
		integer, allocatable :: inRange(:)
		!inRange gets the indices of the atoms in range, with trailing zeros
		!this becomes localList, and any do loop can be terminated for an index of 0

		array_shape = shape(posns)!should be 3,N since the first index varies more rapidly
		natoms = array_shape(2)
		counter = 0
		
		allocate(inRange(natoms))
		inRange = 0
		
		do idx=1,natoms
			checkNear = isNear(posns(:,idx),centre,cell,r_c)
			if(checkNear) then
				counter = counter + 1
				inRange(counter) = idx		
			end if
		end do
		
		localList = inRange
		
		deallocate(inRange)
		
	end subroutine getLocalIndices
	
	
	
	logical function isNear(pos1,pos2,cell,r_c)
	!checks if two positions are within r_c of each other in the periodic cell
	!moves pos1 forward and backward by the lattice vector in each of the 3 directions while checking
		implicit none
		
		real(8), intent(in) :: pos1(3), pos2(3)
		real(8), intent(in) :: cell(3,3) !vector elements are given by the first index, the vector labels by the second (first index varies faster in fortran)
		real(8), intent(in) :: r_c
		real(8) :: dist, smallestDisp(3)
		real(8), dimension(3) :: zeroVect = (/ 0,0,0 /)
		integer :: ii
		
		smallestDisp = smallestDisplacement(pos1,pos2,cell)
		
		dist = getDist(smallestDisp,zeroVect)
		
		if (dist.le.r_c) then
			isNear = .True.
		else
			isNear = .False.
		end if
	end function isNear



	logical function checkInCell(posn,cell)
	!checks if the function exists inside the bounds of the cell
	!project the position onto each lattice vector and see if it's longer than the vector
		implicit none
		real(8), intent(in) :: posn(:)
		real(8), intent(in) :: cell(3,3)!vector elements are given by the first index, the vector labels by the second (first index varies faster in fortran)
		real(8) :: latticeVectLengths(3)
		real(8), dimension(3) :: zeroVect = (/0,0,0/)
		real(8) :: projectLength, sideLength
		integer :: ii
		
		checkInCell = .True.
		
		do ii = 1,3
			sideLength = sqrt(dot(cell(:,ii),cell(:,ii)))
			projectLength = dot(cell(:,ii),posn)/sideLength
			
			print*, projectLength, sideLength
			
			if (projectLength.gt.sideLength) then
				checkInCell = .False.
			end if
		end do
	end function checkInCell
		
		
		
		
		
		
	real(8) function dot(vect1,vect2)
	!gets the dot product between two cartesian vectors
		implicit none
		real(8), intent(in) :: vect1(3)
		real(8), intent(in) :: vect2(3)
		integer :: ii
		
		dot = 0
		do ii = 1,3
			dot = dot + vect1(ii)*vect2(ii)
		end do
	end function dot


	function smallestDisplacement(pos1,pos2,cell)
	!figures out the smallest relative displacements between two positions in the periodic cell
		implicit none
		
		real(8), intent(in) :: pos1(3), pos2(3)
		real(8), intent(in) :: cell(3,3)
		real(8) :: dist(26)
		real(8) :: latt_vect(3,26)!vectors to get to the 26 cells neighbouring the central cell
		real(8) :: currentMin, currentMinVec(3)
		logical :: inCell
		integer :: ii
		real(8), dimension(3) :: smallestDisplacement
		
		currentMin = 0
		currentMinVec = 0	
		
		latt_vect(:,1:3) = cell
			
		latt_vect(:,4) = cell(:,1) + cell(:,2)
		latt_vect(:,5) = cell(:,1) + cell(:,3)
		latt_vect(:,6) = cell(:,1) - cell(:,2)
		latt_vect(:,7) = cell(:,1) - cell(:,3)
		latt_vect(:,8) = cell(:,2) + cell(:,3)
		latt_vect(:,9) = cell(:,2) - cell(:,3)
		latt_vect(:,10) = cell(:,1) + cell(:,2) + cell(:,3)
		latt_vect(:,11) = cell(:,1) + cell(:,2) - cell(:,3)
		latt_vect(:,12) = cell(:,1) - cell(:,2) + cell(:,3)
		latt_vect(:,13) = cell(:,1) - cell(:,2) - cell(:,3)
		
		latt_vect(:,14:26) = -1*latt_vect(:,1:13)	
		
		do ii = 1,26
			dist(ii) = getDist(pos1+latt_vect(:,ii),pos2)
			if (dist(ii).lt.currentMin) then
				currentMin = dist(ii)
				currentMinVec = pos1+latt_vect(:,ii)-pos2
			end if
		end do

		smallestDisplacement = currentMinVec
	end function smallestDisplacement
		
	
	
	
	subroutine localCart2Polar(posns,polarPosns, centre, localList)
		implicit none
		real(8), intent(in) :: posns(:,:)!elements of the vector on the left, vector number on the right
		real(8), intent(in) :: centre(3)
		integer, intent(in) :: localList(:)
		real(8) :: polarposns(:,:)
		real(8), allocatable :: r(:), theta(:), phi(:), posnsShifted(:,:)
		integer, dimension(2) :: natoms_dim
		integer :: idx, ii
		
		natoms_dim = shape(posns)
		
		allocate(r(natoms_dim(1)))
		allocate(theta(natoms_dim(1)))
		allocate(phi(natoms_dim(1)))
		
		allocate(posnsShifted(natoms_dim(1),natoms_dim(2)))
		
		do ii=1,natoms_dim(1)
			!make this local with localList, it contains indices of local atoms, with trailing zeros
			idx = localList(ii)
			if (idx.eq.0)then
				exit
			end if
			posnsShifted(:,idx) = posns(:,idx)-centre
			r(idx) = getDist(posns(:,idx),centre)
			phi(idx) = acos(posnsShifted(3,idx)/r(idx))
			if (abs(posnsShifted(1,idx)).gt.1e-15 .or. abs(posnsShifted(2,idx)).gt.1e-15) then
				theta(idx) = atan(posnsShifted(2,idx)/posnsShifted(1,idx))
			else
				theta(idx) = 0
			end if
			polarposns(1, idx) = r(idx)
			polarposns(2, idx) = theta(idx)
			polarposns(3, idx) = phi(idx)
		end do

		
	end subroutine localCart2Polar
	
	
	
	integer function factorial(N)
		implicit none
		integer, intent(in) :: N
		integer :: i
		
		if (N.eq.0) then
		factorial = 1
		else if (N.lt.0) then
		write(*,*) "Integer input to factorial needs to be greater than 0"
		else if (N.gt.0) then
		factorial = 1
		
		do i = 1,N
			factorial = factorial*i
		end do
		end if
	end function factorial
	
	
	
	subroutine factorial_array(array)
	!if the maximum number of the factorial is known, then it may be more efficient to write factorials up to that number in an array.
	!this subroutine assumes that the array has already been allocated to the size of the largest factorial needed
	
		implicit none
		integer, intent(inout) :: array(:)
		integer :: ii, array_size
		
		array_size = size(array)
		
		if (array_size.eq.0) then
			write(*,*) "factorial_array appears to have size 0"
		else if (array_size .eq.1) then
			array = 1
		else if (array_size.gt.1) then
			array(1) = 1
			do ii = 2, array_size
				array(ii) = array(ii-1)*ii
			end do
		end if
	end subroutine factorial_array
		
		
	
	integer function minAbsFloor(x)
	!floor function that returns the integer value closest to zero
	!semi copied from Andrew Fowler's code (I looked at it before I wrote this)
	!see his github: https://github.com/andrew31416/densityregression and look for python_floor in spherical_harmonics.f90
		implicit none
		real(8), intent(in) :: x
		
		if (x.ge.0.0d0) then
			minAbsFloor = int(floor(x))
		else
			minAbsFloor = int(ceiling(x))
		end if
	end function minAbsFloor
		
		
	
end module utilities
