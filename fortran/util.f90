module util
	implicit none
	
	real(8), external :: dnrm2
	real(8), external :: ddot
	real(8), external :: dysev
	real(8), external :: dgemm
	
	real(8) :: diff = dble(1e-8)!tolerance in difference between two values that should be the same
	!I made this one up, check to see what the proper difference should be
	
	type bispectParamType
	!derived type for bispectrum parameters
		real(8) :: r_c
		integer :: l_max, n_max
	end type bispectParamType
	
	type systemStateType
	!derived type containing information about the system
		real(8) :: cell(3,3)!vector elements are given by the first index, the vector labels by the second (first index varies faster in fortran)
		integer :: nAtoms ! number of atoms
		real(8), allocatable :: atomPositions(:,:)!first index gives vector element, second give atom number
	end type systemStateType	
	
	type pointType
	!properties that will vary for each point in the space
		real(8) :: pointPosition(3)
		real(8), allocatable :: neighbourList(:,:)!list of vectors within r_c of the point position (unshifted)
		integer :: numNeighbours
	end type pointType
	

	public getDist
	public checkInCell
	public localCart2Polar
	public factorial
	public minAbsFloor
	public assignBispectParams
	public assignSystemState
	
	
	
	type(bispectParamType), public :: bispectParams
	type(systemStateType), public :: systemState
	type(pointType), public :: point
	
	

	contains
	
	
	
	subroutine assignPoint(point, centre)
		implicit none
		
		type(pointType), intent(inout) :: point
		real(8), intent(in) :: centre
		point%pointPosition = centre
	end subroutine assignPoint
	
	
	
	subroutine assignBispectParams(bispectParams,r_c,l_max,n_max)
		implicit none
		
		type(bispectParamType), intent(inout) :: bispectParams
		integer, intent(in) :: l_max, n_max
		real(8), intent(in) :: r_c
		
		bispectParams%r_c = r_c
		bispectParams%l_max = l_max
		bispectParams%n_max = n_max	
	
	end subroutine assignBispectParams
	
	
	
	subroutine assignSystemState(systemState, cell, atomPositions)
		implicit none
		
		type(systemStateType), intent(inout) :: systemState
		real(8), intent(in) :: cell(3,3)
		real(8), intent(in) :: atomPositions(:,:)
		
		integer :: natoms, arrayShape(2)
		
		arrayShape = shape(atomPositions)
		natoms = arrayShape(2)
		
		systemState%cell = cell
		systemState%atomPositions = atomPositions
		systemState%natoms = natoms
	
	end subroutine assignSystemState
	
	
	
	function getDist(pos1, pos2)
	!gets distance between two cartesian coordinates	
		implicit none
		real(8), intent(in) :: pos1(3)
		real(8), intent(in) :: pos2(3)
		
		real(8), external :: dnrm2
		
		real(8) :: getDist, diff(3)
		
		diff = pos1-pos2
		
		getDist = dnrm2(3,diff,1)!sqrt(sum((pos1-pos2)**2))
		
	end function getDist


	
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
		real(8), intent(in) :: vect1(:)
		real(8), intent(in) :: vect2(:)
		integer :: ii, dim1, dim2
		
		real(8), external :: ddot
		
		dim1 = size(vect1)
		dim2 = size(vect2)
		
		if (dim1.ne.dim2) then
			print*,"The two vectors don't have the same dimensions"
		end if
		
		dot = ddot(dim1,vect1,1,vect2,1)
!		dot = 0
!		do ii = 1,3
!			dot = dot + vect1(ii)*vect2(ii)
!		end do

		
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
		
		call assignLatticeVectors(latt_vect,cell)
		
		do ii = 1,26
			dist(ii) = getDist(pos1+latt_vect(:,ii),pos2)
			if (dist(ii).lt.currentMin) then
				currentMin = dist(ii)
				currentMinVec = pos1+latt_vect(:,ii)-pos2
			end if
		end do

		smallestDisplacement = currentMinVec
	end function smallestDisplacement
	
	
	
	subroutine assignLatticeVectors(latt_vect,cell)
		implicit none
		
		real(8), intent(in) :: cell(3,3)
		real(8), intent(inout) :: latt_vect(3,26)
		
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
	end subroutine assignLatticeVectors
		
	
	
	subroutine localCart2Polar(polarPosns, point)
		implicit none
		type(pointType), intent(in) :: point
		real(8), intent(inout) :: polarposns(:,:)
		
		real(8), allocatable :: r(:), theta(:), phi(:), posnsShifted(:,:)
		
		integer :: idx, ii
		
		allocate(r(point%numNeighbours))
		allocate(theta(point%numNeighbours))
		allocate(phi(point%numNeighbours))
		allocate(posnsShifted(3,point%numNeighbours))
		
		do idx=1,point%numNeighbours
			!make this local with localList, it contains indices of local atoms, with trailing zeros

			posnsShifted(:,idx) = point%neighbourList(:,idx)-point%pointPosition
			r(idx) = getDist(point%neighbourList(:,idx),point%pointPosition)
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

		deallocate(r)
		deallocate(theta)
		deallocate(phi)
		deallocate(posnsShifted)
		
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
	
	integer function odd_factorial(N)
		!gets the factorial as a product of all odd numbers up to N
		implicit none
		integer, intent(in) :: N
		integer :: N_odd, max_idx, i
		
		if ( (N.gt.0).AND.(MOD(N,2).ne.0)) then
			!N is odd and positive
			N_odd = N
		else if (N.gt.0) then
			!N is even and positive (greater than 1)
			N_odd = N-1
		else if (N.eq.0) then
			N_odd = 0
		else
			write(*,*) "Integer input to odd factorial needs to be positive"
		end if
		
		if (N_odd.eq.0) then
			odd_factorial = 1
		else
			odd_factorial = 1
			max_idx = (N_odd+1)/2
			do i=1,max_idx
				odd_factorial = odd_factorial*((2*i)-1)
			end do
		end if
	end function odd_factorial
	
	
	
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
		
		

	subroutine tidyPoint(point)
	!deallocates everything in the point type once it's finished its use
		implicit none
		type(pointType), intent(inout) :: point
		
		if (allocated(point%neighbourList)) then
			deallocate(point%neighbourList)
		end if
	end subroutine tidyPoint
	
	
	
	subroutine tidySystem(systemState)
	!deallocated anything still allocated in the systemState type
		implicit none
		type(systemStateType), intent(inout) :: systemState
		
		if (allocated(systemState%atomPositions)) then
			deallocate(systemState%atomPositions)
		end if
	end subroutine tidySystem
	
	

	subroutine sqrtInvSymmMatrix(symmMat,sqrtInv)
	!gets the square root of the inverse of a symmetric matrix
	!First do an eigenvalue decomposition S = R L R^-1
	!Where R is a unitary matrix and L is the matrix of eigenvalues
	!S^0.5 = R L^0.5 R^-1 since the square of a matrix has the same eigenvectors, with a squared eigenvalues
	!S^-1 = R L^-1 R^-1
	!Thus S^-0.5 = R L^-0.5 R^-1
		implicit none
		real(8), intent(in) :: symmMat(:,:)
		real(8), intent(inout) :: sqrtInv(:,:)
		
		real(8), allocatable :: R(:,:)
		real(8), allocatable :: invR(:,:)
		real(8), allocatable :: Lambda(:,:)
		real(8), allocatable :: invSqrtLambda(:,:)
		real(8), allocatable :: intermediate(:,:)
		integer :: arrayShape(2), ii
		logical :: goodDecomp
		
		
		arrayShape = shape(symmMat)
		allocate(R(arrayShape(1),arrayShape(2)))
		allocate(invR(arrayShape(1),arrayShape(2)))
		allocate(Lambda(arrayShape(1),arrayShape(2)))
		allocate(invSqrtLambda(arrayShape(1),arrayShape(2)))
		allocate(intermediate(arrayShape(1),arrayShape(2)))
		
		call eigenDecomp(symmMat,R,invR,Lambda)
		
		goodDecomp = checkDecomp(symmMat,R,invR,Lambda)
		
		if (goodDecomp) then
			do ii = 1,arrayShape(1)
				invSqrtLambda(ii,ii) = 1/sqrt(Lambda(ii,ii))
			end do
			call MatrixMult(invSqrtLambda,invR,intermediate)
			call MatrixMult(R,intermediate,sqrtInv)
			
			call checkSqrtInv(sqrtInv,symmMat)
		else
			print *, "something went wrong with the eigenvector decomposition"
		end if
		
		deallocate(R)
		deallocate(invR)
		deallocate(Lambda)
		deallocate(invSqrtLambda)
		deallocate(intermediate)
	
	end subroutine sqrtInvSymmMatrix
	
	
	
	subroutine eigenDecomp(symmMat,R,invR,Lambda)
	!gets the eigenvalue decomposition of a symmetric matrix
	!so symmMat = R * Lambda * invR
		implicit none
		real(8), intent(in) :: symmMat(:,:)
		real(8), intent(inout) :: R(:,:)
		real(8), intent(inout) :: invR(:,:)
		real(8), intent(inout) :: Lambda(:,:)
		
		external :: dsyev
		
		logical :: isSymm
		integer :: arrayShape(2), info, workSize, ii, jj
		real(8), allocatable :: eigenValues(:), work(:)
		
		isSymm = checkSquareSymm(symmMat)
		
		if (isSymm) then
			info = 0
			Lambda = 0.0d0
			arrayShape = shape(symmMat)
			R = symmMat
			workSize = 3*arrayShape(1)-1
			allocate(eigenValues(arrayShape(1)))
			allocate(work(workSize))
			call dsyev('V','U',arrayShape(1),R,arrayShape(1),eigenValues,work, workSize, info)
			
			if (info.ne.0) then
				print *, "eigenvector decomposition failed, info value: ", info
			end if
			
			do ii = 1, arrayShape(1)
				do jj = 1, arrayShape(2)
					!can get the inverse of R from the transpose since R should be unitary
					invR(jj,ii) = R(ii,jj)
				end do
				Lambda(ii,ii) = eigenValues(ii)
			end do
			
			deallocate(work)
			deallocate(eigenValues)
			
		else
			print *, "Not doing eigen decomp on non symmetric matrix"
		end if
		
	
	
	end subroutine eigenDecomp
	
	
	
	logical function checkSquareSymm(matrix)
	!checks the matrix is square and symmetric
		implicit none
		real(8), intent(in) :: matrix(:,:)
		
		integer :: dim(2)
		integer :: ii, jj
		
		
		checkSquareSymm = .True.
		
		dim = shape(matrix)
		
		if (dim(1).ne.dim(2)) then
			print *, "Matrix isn't square when it's supposed to be"
			checkSquareSymm = .False.
		end if
		
		do ii = 1, dim(1)
			do jj = ii, dim(2)
				if (abs(matrix(jj,ii)-matrix(ii,jj)).gt.diff) then
					checkSquareSymm = .False.		
					print *, "matrix not symmetric"		
				end if
			end do		
		end do
	
	end function checkSquareSymm
	
	
	
	logical function checkDecomp(matrix, R, invR, Lambda)
	!checks that the eigenvector decomposition worked as expected
	!ie that matrix = R * lambda * invR
		implicit none
		real(8), intent(in) :: matrix(:,:), R(:,:), invR(:,:), Lambda(:,:)
		
		external :: dgemm
		
		real(8), allocatable :: matProduct(:,:), interim(:,:)
		integer :: arrayShape(2), length, ii ,jj
		
		checkDecomp = .True.
		
		arrayShape = shape(matrix)
		length = arrayShape(1)
		
		allocate(matProduct(arrayShape(1),arrayShape(2)))
		allocate(interim(arrayShape(1),arrayShape(2)))
		
		call MatrixMult(Lambda,invR,interim)
		call MatrixMult(R,interim,matProduct)
		
		do ii = 1, length
			do jj = 1, length
				if (abs(matProduct(ii,jj)-matrix(ii,jj)).gt.diff) then
					checkDecomp = .False.
				end if
			end do	
		end do
		
		deallocate(matProduct)
		deallocate(interim)		
	
	end function checkDecomp
	
	
	
	subroutine MatrixMult(mat1,mat2,matProduct)
	!multiplies two matrices: matProduct = mat1*mat2
	!matProducts needs to have dimensions consistent with mat1 and mat2
		implicit none
		real(8), intent(in) :: mat1(:,:)
		real(8), intent(in) :: mat2(:,:)
		real(8), intent(inout) :: matProduct(:,:)
		
		external :: dgemm
		
		integer :: shape1(2), shape2(2), shape3(2)
		
		shape1 = shape(mat1)
		shape2 = shape(mat2)
		shape3 = shape(matProduct)
		
		if (shape3(1).ne.shape1(1).or.shape3(2).ne.shape2(2).or.shape1(2).ne.shape2(1)) then
			print *, "the shape of the matrices is incompatible"
			print *, "Matrix 1 shape: ", shape1
			print *, "Matrix 2 shape: ", shape2
			print *, "Product shape: ", shape3
		else
			call dgemm('n','n',shape3(1),shape3(2),shape1(2),1.0,mat1,shape1(1),mat2,shape2(1),0.0,matProduct,shape3(1))
		end if
	end subroutine MatrixMult
		
	
	
	subroutine checkSqrtInv(sqrtInv,symmMat)
	!checks that you indeed have the square root of the inverse
	!does so by checking that symmMat*(sqrtInv)^2 = I
		implicit none
		real(8), intent(in) :: sqrtInv(:,:)
		real(8), intent(in) :: symmMat(:,:)
		real(8), allocatable :: inv(:,:)
		real(8), allocatable :: I(:,:)
		integer :: arrayShape(2), ii, jj
		
		arrayShape = shape(sqrtInv)
		allocate(inv(arrayShape(1),arrayShape(2)))
		allocate(I(arrayShape(1),arrayShape(2)))
		
		call MatrixMult(sqrtInv,sqrtInv,inv)
		call MatrixMult(symmMat,inv,I)
		
		!check that I is the identity
		do ii = 1, arrayShape(1)
			do jj = 1, arrayShape(1)
				if (ii.ne.jj) then
					if(abs(I(ii,jj)-0.0d0).gt.diff) then
						print *, "The calculation for the square root of the inverse has failed"
					end if
				else if (ii.eq.jj) then
					if (abs(I(ii,jj)-dble(1)).gt.diff) then
						print *, "The calculation for the square root of the inverse has failed"
					end if
				end if
			end do		
		end do
	
	end subroutine checkSqrtInv
	
	
	
end module util
