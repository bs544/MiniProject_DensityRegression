module util
	implicit none
	
	!To Do: either fix or get rid of MatMult, and next time look if something's built in before you write it yourself
	
	real(8), external :: dnrm2
	real(8), external :: ddot
	real(8), external :: dysev
	real(8), external :: dgemm
	
	real(8) :: diff = dble(1e-5)!tolerance in difference between two values that should be the same
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
		real(8), intent(in) :: centre(3)
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
		!integer :: ii
		
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
		real(8), intent(in) :: posn(3)!:
		real(8), intent(in) :: cell(3,3)!vector elements are given by the first index, the vector labels by the second (first index varies faster in fortran)
		!real(8) :: latticeVectLengths(3)
		!real(8), dimension(3) :: zeroVect = (/0,0,0/)
		real(8) :: projectLength, sideLength
		integer :: ii
		
		checkInCell = .True.
		
		do ii = 1,3
			sideLength = sqrt(dot(cell(:,ii),cell(:,ii)))
			projectLength = dot(cell(:,ii),posn)/sideLength
			
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
		integer :: dim1, dim2!, ii
		
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
		!logical :: inCell
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
		
	
	
	subroutine localCart2Polar(polarPosns, point, systemState)
		implicit none
		type(pointType), intent(in) :: point
		type(systemStateType), intent(in) :: systemState
		real(8), intent(inout) :: polarposns(3,systemState%natoms)
		
		real(8), allocatable :: r(:), theta(:), phi(:), posnsShifted(:,:)
		
		integer :: idx!, ii
		real(8), parameter :: pi = 4*atan(1.0d0)
		
		allocate(r(point%numNeighbours))
		allocate(theta(point%numNeighbours))
		allocate(phi(point%numNeighbours))
		allocate(posnsShifted(3,point%numNeighbours))
		
		do idx=1,point%numNeighbours
			!make this local with localList, it contains indices of local atoms, with trailing zeros

			posnsShifted(:,idx) = point%neighbourList(:,idx)-point%pointPosition
			r(idx) = getDist(point%neighbourList(:,idx),point%pointPosition)
			if (r(idx).gt.dble(1e-15)) then
				phi(idx) = acos(posnsShifted(3,idx)/r(idx))
			else
				phi = 0.0d0
			end if
			if (abs(posnsShifted(1,idx)).gt.dble(1e-15)) then
				theta(idx) = atan(posnsShifted(2,idx)/posnsShifted(1,idx))
				!get theta in two negative x quadrants that are neglected
				if(posnsShifted(1,idx).lt.0.0d0.and.posnsShifted(2,idx).lt.0.0d0) then
					!in the bottom left quadrant
					theta(idx) = theta(idx) - pi
				else if(posnsShifted(1,idx).lt.0.0d0.and.posnsShifted(2,idx).ge.0.0d0) then
					!in the top left quadrant
					theta(idx) = theta(idx) + pi
				end if
			else if (posnsShifted(2,idx).ge.0.0d0) then
				theta(idx) = pi/2.0d0!0.0d0
			else 
				theta(idx) = -pi/2.0d0
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
		integer :: ii
		
		if (N.eq.0) then
			factorial = 1
		else if (N.lt.0) then
			write(*,*) "Integer input to factorial needs to be greater than 0"
			factorial = 0
		else if (N.gt.0) then
			factorial = 1
		
			do ii = 1,N
				factorial = factorial*ii
			end do
		end if
	end function factorial
	
	
	
	integer function odd_factorial(N)
		!gets the factorial as a product of all odd numbers up to N
		implicit none
		integer, intent(in) :: N
		integer :: N_odd, max_idx, ii
		
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
			N_odd = 0
		end if
		
		if (N_odd.eq.0) then
			odd_factorial = 1
		else
			odd_factorial = 1
			max_idx = (N_odd+1)/2
			do ii=1,max_idx
				odd_factorial = odd_factorial*((2*ii)-1)
			end do
		end if
	end function odd_factorial
	
	
	
	subroutine factorial_array(array,array_size)
	!if the maximum number of the factorial is known, then it may be more efficient to write factorials up to that number in an array.
	!this subroutine assumes that the array has already been allocated to the size of the largest factorial needed
	
		implicit none
		integer, intent(in) :: array_size
		integer(8), intent(inout) :: array(0:array_size)
		integer :: ii
		
		!array_size = size(array)
		
		if (array_size.lt.0) then
			write(*,*) "bad array size"
		else if (array_size.eq.0) then
			array(0) = 1
		else if (array_size .eq.1) then
			array(0) = 1
			array(1) = 1
		else if (array_size.gt.1) then
			array(0) = 1
			do ii = 1, array_size
				array(ii) = array(ii-1)*ii
			end do
		end if
	end subroutine factorial_array
		
		
	
	integer function minAbsFloor(x)
	!floor function that returns the integer value closest to zero
	!copied from Andrew Fowler's code
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
	
	

	subroutine sqrtInvSymmMatrix(symmMat,sqrtInv,array_size)
	!gets the square root of the inverse of a symmetric matrix
	!First do an eigenvalue decomposition S = R L R^-1
	!Where R is a unitary matrix and L is the matrix of eigenvalues
	!S^0.5 = R L^0.5 R^-1 since the square of a matrix has the same eigenvectors, with a squared eigenvalues
	!S^-1 = R L^-1 R^-1
	!Thus S^-0.5 = R L^-0.5 R^-1
		implicit none
		integer, intent(in) :: array_size
		real(8), intent(in) :: symmMat(array_size,array_size)
		real(8), intent(inout) :: sqrtInv(array_size,array_size)
		
		real(8), dimension(array_size,array_size) :: R
		real(8), dimension(array_size,array_size) :: invR
		real(8), dimension(array_size,array_size) :: Lambda
		real(8), dimension(array_size,array_size) :: invSqrtLambda
		real(8), dimension(array_size,array_size) :: intermediate
		integer :: ii
		logical :: goodDecomp
		
		
		call eigenDecomp(symmMat,R,invR,Lambda,array_size)
		
		goodDecomp = checkDecomp(symmMat,R,invR,Lambda,array_size)
		
		
		if (goodDecomp) then
			invSqrtLambda = 0
			do ii = 1,array_size
				invSqrtLambda(ii,ii) = 1/sqrt(Lambda(ii,ii))
			end do
			!call MatrixMult(invSqrtLambda,invR,intermediate)
			!call MatrixMult(R,intermediate,sqrtInv)
			
			
			intermediate = matmul(invSqrtLambda,invR)
			sqrtInv = matmul(R,intermediate)
			
			call checkSqrtInv(sqrtInv,symmMat,array_size)
		else
			print *, "something went wrong with the eigenvector decomposition"
		end if
	
	end subroutine sqrtInvSymmMatrix
	
	
	
	subroutine eigenDecomp(symmMat,R,invR,Lambda,array_size)
	!gets the eigenvalue decomposition of a symmetric matrix
	!so symmMat = R * Lambda * invR
		implicit none
		integer, intent(in) :: array_size
		real(8), intent(in) :: symmMat(array_size,array_size)
		real(8), intent(inout) :: R(array_size,array_size)
		real(8), intent(inout) :: invR(array_size,array_size)
		real(8), intent(inout) :: Lambda(array_size,array_size)
		
		external :: dsyev
		
		logical :: isSymm
		integer :: arrayShape(2), info, workSize, ii, jj
		real(8), dimension(array_size) :: eigenValues
		real(8), allocatable :: work(:)
		

		
		isSymm = checkSquareSymm(symmMat)
		
		if (isSymm) then
			info = 0
			Lambda = 0.0d0
			arrayShape = shape(symmMat)
			R = symmMat
			workSize = 3*arrayShape(1)-1
			allocate(work(workSize))
!~ 			allocate(eigenValues(arrayShape(1)))

			call dsyev('V','U',arrayShape(1),R,arrayShape(1),eigenValues,work, workSize, info)
			
			if (info.ne.0) then
				print *, "eigenvector decomposition failed, info value: ", info
			end if
			
			do ii = 1, arrayShape(1)
				do jj = 1, arrayShape(2)
					!can get the inverse of R from the transpose since R should be unitary
					invR(jj,ii) = R(ii,jj)
					if (ii.eq.jj) then
						Lambda(ii,ii) = eigenValues(ii)
					else
						Lambda(ii,jj) = 0.0d0
					end if
				end do
				!Lambda(ii,ii) = eigenValues(ii)
			end do
			
			deallocate(work)
!~ 			deallocate(eigenValues)
			
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
	
	
	
	logical function checkDecomp(matrix, R, invR, Lambda,array_size)
	!checks that the eigenvector decomposition worked as expected
	!ie that matrix = R * lambda * invR
		implicit none
		integer, intent(in) :: array_size
		real(8), intent(in) :: matrix(array_size,array_size)
		real(8), intent(in) :: R(array_size,array_size)
		real(8), intent(in) :: invR(array_size,array_size)
		real(8), intent(in) :: Lambda(array_size,array_size)
		
		external :: dgemm
		
		real(8), dimension(array_size,array_size) :: matProduct, interim
		integer :: length, ii ,jj
		
		checkDecomp = .True.
		
		!arrayShape = shape(matrix)
		length = array_size
		
!~ 		allocate(matProduct(array_size,array_size))
!~ 		allocate(interim(array_size,array_size))
		
		!call MatrixMult(Lambda,invR,interim)
		!call MatrixMult(R,interim,matProduct)
		interim = matmul(Lambda,invR)
		matProduct = matmul(R,interim)
		
		
		do ii = 1, length
			do jj = 1, length
				if (abs(matProduct(ii,jj)-matrix(ii,jj)).gt.diff) then
					checkDecomp = .False.
				end if
			end do	
		end do
		
		
!~ 		deallocate(matProduct)
!~ 		deallocate(interim)		
	
	end function checkDecomp
	
	

	subroutine checkSqrtInv(sqrtInv,symmMat,array_size)
	!checks that you indeed have the square root of the inverse
	!does so by checking that symmMat*(sqrtInv)^2 = I
		implicit none
		integer, intent(in) :: array_size
		real(8), intent(in) :: sqrtInv(array_size,array_size)
		real(8), intent(in) :: symmMat(array_size,array_size)
		real(8), dimension(array_size,array_size) :: inv
		real(8), dimension(array_size,array_size) :: I
 		integer :: ii, jj
 		logical :: failure = .False.

!~ 		allocate(inv(array_size,array_size))
!~ 		allocate(I(array_size,array_size))
		
		!inv = 0
		!I = 0
		
		inv = matmul(sqrtInv,sqrtInv)
		I = matmul(symmMat,inv)
		

		
		
		!check that I is the identity
		do ii = 1, array_size
			do jj = 1, array_size
				if (ii.ne.jj) then
					if(abs(I(ii,jj)-0.0d0).gt.diff) then
						failure = .True.
					end if
				else if (ii.eq.jj) then
					if (abs(I(ii,jj)-dble(1)).gt.diff) then
						failure = .True.
					end if
				end if
			end do		
		end do
		
		if (failure) then
			print *, "The calculation for the square root of the inverse has failed"	
			print *, "I"
			do ii = 1, array_size
				print *, I(:,ii)
			end do	
		end if
!~ 		deallocate(inv)
!~ 		deallocate(I)
	
	end subroutine checkSqrtInv
	
	
	subroutine inverseChecker(mat,invmat,array_size)
	!checks that one square array is the inverse of another square array
		implicit none
		integer, intent(in) :: array_size
		real(8), intent(in) :: mat(array_size,array_size)
		real(8), intent(in) :: invmat(array_size,array_size)
		
		
		real(8),dimension(array_size,array_size) :: I
		integer :: ii,jj
		logical :: failure = .False.
		
		
		I = matmul(mat,invmat)
		 
		!check that I is the identity
		do ii = 1, array_size
			do jj = 1, array_size
				if (ii.ne.jj) then
					if(abs(I(ii,jj)-0.0d0).gt.diff) then
						failure = .True.
					end if
				else if (ii.eq.jj) then
					if (abs(I(ii,jj)-dble(1)).gt.diff) then
						failure = .True.
					end if
				end if
			end do		
		end do
		
		if (failure) then
			print *, "The calculation for the inverse has failed"
			print *, "I"
			do ii = 1, array_size
				print *, I(:,ii)
			end do
		end if
		
	end subroutine inverseChecker
	
	
	
	real(8) function CG(l_1,m_1,l_2,m_2,l,m,fact_array,l_max)
	!Pretty much all copied from Andrew Fowler's code, should go back and write myself
	!See his github here: https://github.com/andrew31416 and look for spherical_harmonics.f90
	!This is the Clebsch-Gordan coefficient C_{l_1,m_1,l_2,m_2}^{l,m}
	!The formula to calculate this can be found on page 238 of 'Quantum Theory of Angular Momentum' by Varshalovich

		implicit none
	
		real(8),intent(in) :: l_1,m_1,l_2,m_2,l,m
		integer, intent(in) :: l_max
		integer(8), intent(in) :: fact_array(0:l_max*3)
	
		!* scratch
		real(8) :: minimum,min_array(1:7),sqrtres
		real(8) :: imin,imax,valu,sumres,sqrtarg
		real(8) :: dble_ii
		integer :: ii
	
		if (abs(m_1 + m_2 - m).gt.1e-15) then
			CG = 0.0d0
		else
			min_array(1) = l_1 + l_2 - l
			min_array(2) = l_1 - l_2 + l
			min_array(3) = -l_1 + l_2 + l
			min_array(4) = l_1 + l_2 + l + 1.0d0
			min_array(5) = l_1 - abs(m_1)
			min_array(6) = l_2 - abs(m_2)
			min_array(7) = l - abs(m)
	
			minimum = minval(min_array)
	
			if (minimum.lt.0.0d0) then
				CG = 0.0d0
			else
				sqrtarg = 1.0d0
				sqrtarg = sqrtarg * fact_array(minAbsFloor(l_1+m_1))
				sqrtarg = sqrtarg * fact_array(minAbsFloor(l_1-m_1))
				sqrtarg = sqrtarg * fact_array(minAbsFloor(l_2+m_2))
				sqrtarg = sqrtarg * fact_array(minAbsFloor(l_2-m_2))
				sqrtarg = sqrtarg * fact_array(minAbsFloor(l+m))
				sqrtarg = sqrtarg * fact_array(minAbsFloor(l-m))
				sqrtarg = sqrtarg * dble((int(2.0d0*l) + 1))
				sqrtarg = sqrtarg * fact_array(minAbsFloor(min_array(1)))
				sqrtarg = sqrtarg * fact_array(minAbsFloor(min_array(2)))
				sqrtarg = sqrtarg * fact_array(minAbsFloor(min_array(3)))
	
				! sqrtarg is int so need to divide after casting to double
				sqrtres = sqrt(sqrtarg / fact_array(minAbsFloor(min_array(4))))
				
				min_array(1) = l_1 + m_2 - l
				min_array(2) = l_2 - m_1 - l
				min_array(3) = 0.0d0
				min_array(4) = l_2 + m_2
				min_array(5) = l_1 - m_1
				min_array(6) = l_1 + l_2 - l
	
				imin = maxval(min_array(1:3))
				imax = minval(min_array(4:6))
				sumres = 0.0d0
				do ii=minAbsFloor(imin),minAbsFloor(imax)
					dble_ii = dble(ii)
					valu = 1.0d0
					valu = valu * fact_array(ii)
					valu = valu * fact_array(minAbsFloor(l_1 + l_2 - l - dble_ii ))
					valu = valu * fact_array(minAbsFloor(l_1 - m_1 - dble_ii ))
					valu = valu * fact_array(minAbsFloor(l_2 + m_2 - dble_ii ))
					valu = valu * fact_array(minAbsFloor(l - l_2 + m_1 + dble_ii ))
					valu = valu * fact_array(minAbsFloor(l - l_1 - m_2 + dble_ii ))
					sumres = sumres + (-1.0d0)**ii / valu
				end do
				CG = sqrtres * sumres
			end if
		end if
		
	
!~ 		implicit none
		
!~ 		real(8), intent(in) :: l_1,m_1,l_2,m_2,l,m
	
!~ 		real(8) :: minimum,min_array(1:7),sqrtres
!~ 		real(8) :: imin,imax,val,sumres,sqrtarg
!~ 		real(8) :: dble_ii
!~ 		integer :: ii
	
!~ 		if (abs(m_1 + m_2 - m).gt.1e-15) then
!~ 			CG = 0.0d0
!~ 		else
!~ 			min_array(1) = l_1 + l_2 - l + 0.0d0
!~ 			min_array(2) = l_1 - l_2 + l + 0.0d0
!~ 			min_array(3) = -l_1 + l_2 + l + 0.0d0
!~ 			min_array(4) = l_1 + l_2 + l + 1.0d0
!~ 			min_array(5) = l_1 - abs(m_1) + 0.0d0
!~ 			min_array(6) = l_2 - abs(m_2) + 0.0d0
!~ 			min_array(7) = l - abs(m) + 0.0d0
	
!~ 			minimum = minval(min_array)
	
!~ 			if (minimum.lt.0.0d0) then
!~ 				CG = 0.0d0
!~ 			else
				
!~ 				sqrtarg = 1.0d0
!~ 				sqrtarg = sqrtarg * factorial(minAbsFloor(l_1+m_1))!+ 0.0d0))
!~ 				sqrtarg = sqrtarg * factorial(minAbsFloor(l_1-m_1))!+ 0.0d0))
!~ 				sqrtarg = sqrtarg * factorial(minAbsFloor(l_2+m_2))!+ 0.0d0))
!~ 				sqrtarg = sqrtarg * factorial(minAbsFloor(l_2-m_2))!+ 0.0d0))
!~ 				sqrtarg = sqrtarg * factorial(minAbsFloor(l+m))!+ 0.0d0))
!~ 				sqrtarg = sqrtarg * factorial(minAbsFloor(l-m))!+ 0.0d0))
!~ 				sqrtarg = sqrtarg * dble((int(2.0d0*l) + 1))
!~ 				sqrtarg = sqrtarg * factorial(minAbsFloor(min_array(1)))
!~ 				sqrtarg = sqrtarg * factorial(minAbsFloor(min_array(2)))
!~ 				sqrtarg = sqrtarg * factorial(minAbsFloor(min_array(3)))
				
!~ 				! sqrtarg is int so need to divide after casting to double
!~ 				sqrtres = sqrt(sqrtarg / factorial(minAbsFloor(min_array(4))))
				
!~ 				print *, sqrtarg
!~ 				print *, factorial(minAbsFloor(min_array(4)))
				
!~ 				min_array(1) = l_1 + m_2 - l!+ 0.0d0
!~ 				min_array(2) = l_2 - m_1 - l!+ 0.0d0
!~ 				min_array(3) = 0.0d0
!~ 				min_array(4) = l_2 + m_2!+ 0.0d0
!~ 				min_array(5) = l_1 - m_1!+ 0.0d0
!~ 				min_array(6) = l_1 + l_2 - l!+ 0.0d0
				
!~ 				imin = maxval(min_array(1:3))
!~ 				imax = minval(min_array(4:6))
!~ 				sumres = 0.0d0
!~ 				do ii=minAbsFloor(imin),minAbsFloor(imax)
!~ 					dble_ii = dble(ii)
!~ 					val = 1.0d0
!~ 					val = val * factorial(ii)
!~ 					val = val * factorial(minAbsFloor(l_1 + l_2 - l - dble_ii ))
!~ 					val = val * factorial(minAbsFloor(l_1 - m_1 - dble_ii ))
!~ 					val = val * factorial(minAbsFloor(l_2 + m_2 - dble_ii ))
!~ 					val = val * factorial(minAbsFloor(l - l_2 + m_1 + dble_ii ))
!~ 					val = val * factorial(minAbsFloor(l - l_1 - m_2 + dble_ii ))
!~ 					sumres = sumres + (-1.0d0)**ii / val
!~ 				end do
!~ 				print *, sumres
!~ 				CG = sqrtres * sumres
!~ 			end if
!~ 		end if
	
	end function CG	
	
	
end module util
