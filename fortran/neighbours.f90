module neighbours
	!module for creating a list of vectors of atoms within a cutoff radius
	use util
	implicit none
	
	!public getLocalIndices	
	public getNeighbours
	
	integer, private :: buffer_size_default = 100
	
	contains
			
		

!	subroutine getLocalIndices(posns,r_c, centre,localList,cell)
!	!gets an array of indices for the positions in a cutoff sphere around the centre
!		implicit none
!		real(8), intent(in) :: r_c
!		real(8), intent(in) :: posns(:,:)
!		real(8), intent(in) :: centre(3)
!		real(8), intent(in) :: cell(3,3)!vector elements are given by the first index, the vector labels by the second (first index varies faster in fortran)
!		integer :: localList(:)
!		logical :: checkNear
!		integer :: array_shape(2)
!		integer :: idx, natoms, counter
!		integer, allocatable :: inRange(:)
!		!inRange gets the indices of the atoms in range, with trailing zeros
!		!this becomes localList, and any do loop can be terminated for an index of 0

!		array_shape = shape(posns)!should be 3,N since the first index varies more rapidly
!		natoms = array_shape(2)
!		counter = 0
		
!		allocate(inRange(natoms))
!		inRange = 0
		
!		do idx=1,natoms
!			checkNear = isNear(posns(:,idx),centre,cell,r_c)
!			if(checkNear) then
!				counter = counter + 1
!				inRange(counter) = idx		
!			end if
!		end do
		
!		localList = inRange
		
!		deallocate(inRange)
		
!	end subroutine getLocalIndices
	


	subroutine getNeighbours(systemState,point,r_c, buffer_size)
	!gets the number of atoms within r_c, accounting for the periodic cell
	!also get the unshifted vectors of all atoms within r_c of the centre
	!these two values are stored in point
		type(systemStateType), intent(in) :: systemState
		type(pointType), intent(inout) :: point
		real(8), intent(in) :: r_c 
		integer, intent(in), optional :: buffer_size
		
		integer :: maxLattVect(3)!given the ratio of r_c to cell sides, figure out how many lattice vectors to try
		integer :: ii, jj, kk, idx, counter
		real(8) :: dist, shiftedPosn(3), sideLengths(3)
		real(8), allocatable :: temporaryList(:,:)
		
		if (present(buffer_size)) then
			buffer_size_default = buffer_size
		end if
		
		allocate(temporaryList(3,buffer_size_default))
		
		
		counter = 0
		
		do ii = 1,3
			sideLengths(ii) = sqrt(dot(systemState%cell(:,ii),systemState%cell(:,ii)))
		end do
		
		maxLattVect = ceiling(sideLengths/r_c)
		
		do idx = 1, systemState%nAtoms
			do ii = -maxLattVect(1),maxLattVect(1)
				do jj = -maxLattVect(2),maxLattVect(2)
					do kk = -maxLattVect(3),maxLattVect(3)
						shiftedPosn = systemState%atomPositions(:,idx)+ii*systemState%cell(:,1)+jj*systemState%cell(:,2)+kk*systemState%cell(:,3)
						dist = getDist(shiftedPosn,point%pointPosition)
						if (dist.lt.r_c) then
							counter = counter + 1
							temporaryList(:,counter) = shiftedPosn
						end if
					end do
				end do
			end do
		end do
		
		allocate(point%neighbourList(3,counter))
		point%neighbourList = temporaryList(:,1:counter)
		point%numNeighbours = counter
		deallocate(temporaryList)

	
	end subroutine getNeighbours
	
	
end module neighbours
