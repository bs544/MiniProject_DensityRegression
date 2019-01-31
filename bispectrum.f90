module bispectrum
	use utilities
	use BasisFunctions, only: sphericalHarm, radialBasis, CG
	implicit none
	!bispectrum computed as in PHYSICAL REVIEW B 87, 184115 (2013)
	
	
	type bispectParams
	!derived type for bispectrum parameters
		real(8) :: r_c
		integer :: l_max, n_max
	end type bispectParams
	
	type systemState
	!derived type containing information about the system
		real(8) :: cell(3,3)!vector elements are given by the first index, the vector labels by the second (first index varies faster in fortran)
		real(8) :: atom_positions(:,:)!first index gives vector element, second give atom number
		integer :: nAtoms ! number of atoms
	end type systemState
	
	
	
	
	
	
	
	
	subroutine localProjector(atomPositions,centre,C_nlm,n,l,m,r_c)
	!around a general point in space, project atomic positions onto the basis functions
	
		implicit none
		real(8), intent(in) :: atomPositions(:,:)
		real(8), intent(in) :: centre(:)
		real(8), intent(in) :: C_nlm, r_c
		integer, intent(in) :: n,l,m
		
		real(8), allocatable :: shiftedPositions(:,:), polarPositions(:,:)
		integer :: arrayShape(2)
		integer :: ii
		
		arrayShape = shape(atomPositions)!should be dimensions,#atoms
		allocate(shiftedPositions(arrayShape))
		allocate(polarPositions(arrayShape))
		
		
		do ii = 1, arrayShape(2)
			shiftedPositions(:,ii) = atomPositions(:,ii) - centre
			
		
	
	
	
	
	end subroutine localProjector
	
	
	!As above, but this is done around an atom, and summed over every atom
	subroutine globalProjector()
	
	
	end subroutine globalProjector
	
	
	

		
		
		
