module utilities
	implicit none
	
	public get_dist

	contains
	
	!gets distance between two cartesian coordinates
	function get_dist(pos1, pos2)
		implicit none
		real, intent(in) :: pos1(3)
		real, intent(in) :: pos2(3)
		real :: get_dist
		
		get_dist = sqrt(sum((pos1-pos2)**2))
		
	end function get_dist


	!gets an array of indices for the positions in a cutoff sphere around the centre
	subroutine getLocalIndices(posns,r_c, centre,localList,cell)
		implicit none
		real, intent(in) :: r_c
		real, intent(in) :: posns(:,:)
		real, intent(in) :: centre(3)
		real, intent(in) :: cell(3,3)
		integer :: localList(:)
		real :: dist
		integer :: array_shape(2)
		integer :: idx, natoms, counter
		integer, allocatable :: inRange(:)
		!inRange gets the indices of the atoms in range, with trailing zeros
		!this becomes localList, and any do loop can be terminated for an index of 0

		array_shape = shape(posns)
		natoms = array_shape(1)
		counter = 0
		
		allocate(inRange(natoms))
		inRange = 0
		
		do idx=1,natoms
			dist = get_dist(posns(idx,:),centre)
			if(dist < r_c) then
				counter = counter + 1
				inRange(counter) = idx		
			end if
		end do
		
		localList = inRange
		
		deallocate(inRange)
		
	end subroutine getLocalIndices
	
	
	
	logical function isNear(pos1,pos1,cell,r_c)
		implicit none
		
		real, intent(in) :: pos1(3), pos2(3)
		real, intent(in) :: cell(3,3) !vector elements are given by the first index, the vector labels by the second (first index varies faster in fortran)
		real, intent(in) :: r_c
		real :: dist(6), latt_vect(3,6)
		integer :: ii
		
		latt_vect(:,1:3) = cell
		latt_vect(:,4:6) = -1*cell
		
		isNear = .False.
		
		do ii = 1,6
		
		
		
		
	end function isNear
	
	
	
	subroutine localCart2Polar(posns,polarPosns, centre)
		implicit none
		real, intent(in) :: posns(:,:)
		real, intent(in) :: centre(3)
		real :: polarposns(:,:)
		real, allocatable :: r(:), theta(:), phi(:), posnsShifted(:,:)
		integer, dimension(2) :: natoms_dim
		integer :: idx
		
		natoms_dim = shape(posns)
		
		allocate(r(natoms_dim(1)))
		allocate(theta(natoms_dim(1)))
		allocate(phi(natoms_dim(1)))
		
		allocate(posnsShifted(natoms_dim(1),natoms_dim(2)))
		
		do idx=1,natoms_dim(1)
			posnsShifted(idx,:) = posns(idx,:)-centre
			r(idx) = get_dist(posns(idx,:),centre)
		end do
		
		theta = acos(posnsShifted(:,3)/r)
		phi = atan(posnsShifted(:,2)/posnsShifted(:,1))
		
		polarposns(:,1) = r
		polarposns(:,2) = theta
		polarposns(:,3) = phi
		
	end subroutine localCart2Polar
	
	
	
end module utilities
