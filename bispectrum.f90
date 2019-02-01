module bispectrum
	use utilities
	use BasisFunctions, only: sphericalHarm, radialBasis, CG
	use neighbours
	implicit none
	!bispectrum computed as in PHYSICAL REVIEW B 87, 184115 (2013)
	!assumption is made that the cutoff radius is less than the cell length
	
	
	contains
	

	
	
	subroutine localProjector(atomPositions,centre,cell,coeffs,n_max,l_max,r_c)
	!around a general point in space, project atomic positions onto the basis functions
	
		implicit none
		real(8), intent(in) :: atomPositions(:,:)
		real(8), intent(in) :: centre(3), cell(3,3)
		real(8), intent(in) :: coeffs(1:n_max,0:l_max,-l_max:l_max)
		real(8), intent(in) :: r_c
		integer, intent(in) :: n_max,l_max
		
		
		real(8), allocatable :: shiftedPositions(:,:), polarPositions(:,:), localList(:)
		integer :: arrayShape(2)
		integer :: ii, n, l, m
		logical :: allGood !if this isn't true, then the atomic positions are out of the cell bounds
		
		coeffs = 0
		allGood = .True.
		
		arrayShape = shape(atomPositions)!should be (/ #dimensions,#atoms /)
		allocate(shiftedPositions(arrayShape))
		allocate(polarPositions(arrayShape))
		allocate(localList(arrayShape(2)))
		
		do ii = 1, arrayShape(2)
			allGood = checkInCell(atomPositions(:,ii),cell)
		end do
		
		
		call getLocalIndices(atomPositions,r_c,centre,localList,cell)
		call localCart2Polar(posns,polarPosns, centre, localList)
		
			
		
	
	
	
	
	end subroutine localProjector
	
	
	!As above, but this is done around an atom, and summed over every atom
	subroutine globalProjector()
	
	
	end subroutine globalProjector
	
	
	

		
		
		
