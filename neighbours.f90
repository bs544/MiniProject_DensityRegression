module neighbours
	!module for creating a list of vectors of atoms within a cutoff radius
	use utilities
	implicit none
	
	
	integer, private :: buffer_size_default = 100
	
	contains
	
	subroutine getNeighbourList(neighbourList,centre,systemState,buffer_size)
	!find the vectors of atoms within r_c of the centre
		implicit none
		type(systemStateType), intent(in) :: systemState
		real(8), intent(in) :: centre(3,3)
		real(8), allocatable, intent(inout) :: neighbourList(:,:)
		integer, intent(in), optional :: buffer_size
		
		if (present(buffer_size))then
			buffer_size_default = buffer_size
		end if
		
		real(8) :: temporaryList(3,buffer_size_default)
		integer :: counter = 0
		
		

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
	
	
end module neighbours
