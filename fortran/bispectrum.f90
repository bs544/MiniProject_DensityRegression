module bispectrum
	use util
	use basis, only: sphericalHarm, radialBasis, overlapMatrix
	use neighbours
	implicit none
	!bispectrum computed as in PHYSICAL REVIEW B 87, 184115 (2013)

	public getLocalBispectrum
	public get_CG_tensor
	
	
	contains
	
	
	
	function getBispectrum(n_max,l_max,r_c, atomPositions,cell,nAtoms,pointPosition,CG_tensor_in,array_length,local)
	!gets the local or global bispectrum, depending on the logical value local
	!returns a vector array, but you'll need to remove the zero values elsewhere
		implicit none
		integer, intent(in) :: n_max
		integer, intent(in) :: l_max
		real(8), intent(in) :: r_c
		real(8), intent(in) :: atomPositions(:,:)
		real(8), intent(in) :: cell(3,3)
		integer, intent(in) :: nAtoms
		real(8), intent(in) :: pointPosition(3)
		real(8), intent(in) :: CG_tensor_in(0:l_max,0:l_max,0:l_max,0:2*l_max,0:2*l_max,0:2*l_max)
		integer, intent(in) :: array_length
		logical, intent(in) :: local
		
		type(systemStateType) :: systemState
		type(pointType) :: point
		type(bispectParamType) :: bispectParams
		
		real(8), dimension(0:l_max,0:l_max,0:l_max,-l_max:l_max,-l_max:l_max,-l_max:l_max) :: CG_tensor
		real(8), dimension(1:n_max,0:l_max,0:l_max,0:l_max) :: bispect_tensor
		real(8), dimension(array_length) :: getBispectrum
		logical, dimension(1:n_max,0:l_max,0:l_max,0:l_max) :: mask
		
		mask = get_output_mask(n_max,l_max)	
			
		
		call assignBispectParams(bispectParams,r_c,l_max,n_max)
		call format_CG(CG_tensor_in,CG_tensor,l_max)
		call assignSystemState(systemState,cell,atomPositions)
		call assignPoint(point,pointPosition)
		
		
		if (local) then
			call getLocalBispectrum(systemState,point,bispectParams,bispect_tensor,CG_tensor)
		else
			call getGlobalBispectrum(systemState,bispectParams,bispect_tensor,CG_tensor)
		end if
		
		!I might add a non zero mask later but I need to add a program to check the number of non zero bispectrum values
		!write(*,*) bispect_tensor
		getBispectrum = 0
		getBispectrum = pack(bispect_tensor,mask)
		
		call tidySystem(systemState)
		call tidyPoint(point)
		
	
	end function getBispectrum
	
	
	
	subroutine getLocalBispectrum(systemState,point,biPrms,bispect_local,CG_tensor)
		implicit none
		type(systemStateType), intent(in) :: systemState
		type(pointType), intent(inout) :: point
		type(bispectParamType), intent(in) :: biPrms
		real(8), intent(inout) :: bispect_local(1:biPrms%n_max,0:biPrms%l_max,0:biPrms%l_max,0:biPrms%l_max)
		real(8), intent(in) :: CG_tensor(0:biPrms%l_max,0:biPrms%l_max,0:biPrms%l_max,&
		-biPrms%l_max:biPrms%l_max,-biPrms%l_max:biPrms%l_max,-biPrms%l_max:biPrms%l_max)	
		
		complex(8), allocatable :: coeffs(:,:,:)
		integer :: n, l, l_1, l_2, m, m_1, m_2
		integer :: n_max, l_max
		real(8) :: r_c
	
		real(8) :: tmp
		
!		logical :: allGood !if this isn't true, then the atomic positions are out of the cell bounds
		
!		do ii = 1, arrayShape(2)
!			allGood = checkInCell(systemState%atomPositions(:,ii),cell)
!		end do
		
		n_max = biPrms%n_max
		l_max = biPrms%l_max
		r_c = biPrms%r_c
		
		allocate(coeffs(1:n_max,0:l_max,-l_max:l_max))
		coeffs = complex(0.0d0,0.0d0)
				
		call localProjector(systemState,point,coeffs,r_c,n_max,l_max)
		
		if (all(coeffs.eq.complex(0.0d0,0.0d0))) then
			bispect_local = 0
			!print *, "zero coefficients"
		else
		
		bispect_local = 0

		do n=1,n_max
			do l=0,l_max
				do l_1=0,l_max
					do l_2=0,l_max
						tmp = 0
						do m=-l,l
							do m_1=-l_1,l_1
								do m_2=-l_2,l_2
									tmp = tmp + RealPart(conjg(coeffs(n,l,m))*CG_tensor(l,l_1,l_2,m,m_1,m_2)*coeffs(n,l_1,m_1)*coeffs(n,l_2,m_2))
								end do
							end do
						end do
						bispect_local(n,l,l_1,l_2) = tmp
					end do
				end do
			end do
		end do	
		end if
				
		deallocate(coeffs)

	
	end subroutine getLocalBispectrum
	
	
	
	subroutine getGlobalBispectrum(systemState,biPrms,bispect_global,CG_tensor)
		implicit none
		type(systemStateType), intent(in) :: systemState
		type(bispectParamType), intent(in) :: biPrms
		real(8), intent(inout) :: bispect_global(1:biPrms%n_max,0:biPrms%l_max,0:biPrms%l_max,0:biPrms%l_max)
		real(8), intent(in) :: CG_tensor(0:biPrms%l_max,0:biPrms%l_max,0:biPrms%l_max,&
		-biPrms%l_max:biPrms%l_max,-biPrms%l_max:biPrms%l_max,-biPrms%l_max:biPrms%l_max)	
		
		complex(8), allocatable :: coeffs(:,:,:)
		integer :: n, l, l_1, l_2, m, m_1, m_2
		integer :: n_max, l_max
	
		real(8) :: tmp
		
		n_max = biPrms%n_max
		l_max = biPrms%l_max
		
		allocate(coeffs(1:n_max,0:l_max,-l_max:l_max))
		
		call globalProjector(systemState,coeffs,biPrms)
		
		if (all(coeffs.eq.0)) then
			bispect_global = 0
		else
		
		do n=1,n_max
			do l=0,l_max
				do l_1=0,l_max
					do l_2=0,l_max
						tmp = 0
						do m=-l,l
							do m_1=-l_1,l_1
								do m_2=-l_2,l_2
									tmp = tmp + RealPart(conjg(coeffs(n,l,m))*CG_tensor(l,l_1,l_2,m,m_1,m_2)*coeffs(n,l_1,m_1)*coeffs(n,l_2,m_2))
								end do
							end do
						end do
						bispect_global(n,l,l_1,l_2) = tmp
					end do
				end do
			end do
		end do	
		
		end if
		
		deallocate(coeffs)
	
	end subroutine getGlobalBispectrum
	
	
	
	subroutine localProjector(systemState,point,coeffs,r_c,n_max,l_max)
	!around a general point in space, project atomic positions onto the basis functions
	
		implicit none
		type(systemStateType), intent(in) :: systemState
		type(pointType), intent(inout) :: point
		integer, intent(in) :: n_max
		integer, intent(in) :: l_max
		real(8), intent(in) :: r_c
		complex(8), intent(inout) :: coeffs(1:n_max,0:l_max,-l_max:l_max)

		
		complex(8), dimension(1:n_max,0:l_max,-l_max:l_max) :: coeffs_
		real(8), allocatable :: polarPositions(:,:)
		complex(8) :: temp
		integer :: ii, n, n_1, n_2, l, m, buffer_size = 100
		real(8), dimension(n_max,n_max) :: W
		real(8), dimension(n_max,n_max) :: inv_S
		

		coeffs = complex(0.0d0,0.0d0)
		temp = complex(0.0d0,0.0d0)

		

		call getNeighbours(systemState, point,r_c,buffer_size)
		if (point%numNeighbours.gt.0) then
			allocate(polarPositions(3,point%numNeighbours))

			
			call getW(W,n_max)
			inv_S = invOverlap(n_max)
			
			
			call localCart2Polar(polarPositions, point,systemState)
			!print *, polarPositions
			!polar positions elements are order r, phi, theta
			
			!create W here, it's independent of position, so shouldn't be recomputed all the time
			coeffs = complex(0.0d0,0.0d0)
			coeffs_ = complex(0.0d0,0.0d0)

			do l = 0,l_max
				do m = -l,l
					do n = 1, n_max
						do ii = 1,point%numNeighbours
							temp = temp + sphericalHarm(polarPositions(3,ii),polarPositions(2,ii),m,l)*radialBasis(polarPositions(1,ii),r_c,n,n_max,W)
						end do
						coeffs_(n, l, m) = temp
						temp = complex(0.0d0,0.0d0)
					end do
					!convert c' to c
					do n_1=1,n_max
						do n_2=1,n_max
							temp = temp + inv_S(n_2,n_1)*coeffs_(n_1,l,m)
						end do
						coeffs(n_1,l,m) = temp
						temp = complex(0.0d0,0.0d0)
					end do
				end do
			end do
			

			deallocate(polarPositions)
!~ 			deallocate(W)
!~ 			deallocate(coeffs_)
!~ 			deallocate(inv_S)
			


		else
			coeffs = complex(0.0d0,0.0d0)
			!print *, "no neighbours, projected coefficients set to zero"
		end if
		call tidyPoint(point)
			
		
	end subroutine localProjector
	
	
	
	subroutine globalProjector(systemState,coeffs,biPrms)
	!sets the point position to an atom position and removes the atom position from the systemState, then calculates the local bispectrum like that
		implicit none
		
		type(systemStateType), intent(in) :: systemState
		type(bispectParamType), intent(in) :: biPrms
		complex(8), intent(inout) :: coeffs(1:biPrms%n_max,0:biPrms%l_max,-biPrms%l_max:biPrms%l_max)
		
		type(systemStateType) :: tmpState
		type(pointType) :: point
		integer :: ii, jj, n_max, l_max, n_atoms
		real(8) :: r_c
		real(8), allocatable :: tmpAtomPosns(:,:)
		complex(8), allocatable :: tmpCoeffs(:,:,:)
		
		n_max = biPrms%n_max
		l_max = biPrms%l_max
		r_c = biPrms%r_c
		n_atoms = systemState%natoms
		
		allocate(tmpCoeffs(1:n_max,0:l_max,-l_max:l_max))
		allocate(tmpAtomPosns(3,n_atoms-1))
		
		do ii=1, n_atoms
			!first create a system state without the central atom, and a point type centered around the central atom
			do jj=1,n_atoms
				if (ii.ne.jj.and.jj.gt.ii) then
					tmpAtomPosns(:,jj-1) = systemState%atomPositions(:,jj)
				else if (ii.ne.jj.and.jj.lt.ii) then
					tmpAtomPosns(:,jj) = systemState%atomPositions(:,jj)
				end if
			end do
			call assignSystemState(tmpState,systemState%cell,tmpAtomPosns)
			call assignPoint(point,systemState%atomPositions(:,ii))
			
			call localProjector(tmpState,point,tmpCoeffs,r_c,n_max,l_max)
			
			coeffs = coeffs + tmpCoeffs
		end do
		
		coeffs = coeffs/n_atoms
		
		deallocate(tmpCoeffs)
		deallocate(tmpAtomPosns)
		call tidySystem(tmpState)
		
		
	end subroutine globalProjector
			

	
	subroutine format_CG(CG_tensor_in,CG_tensor,l_max)
	!The CG_tensor is calculated once for a given l_max in python and passed in each time to fortran (should really check if there's a better way)
	!Since fortran indices start from 0, the CG_tensor has to be formatted to get passed into python
	!Since I don't want to add all those little bits to the program, I'm changing the CG_tensor to the way I originally planned it
		implicit none
		real(8), intent(in) :: CG_tensor_in(0:l_max,0:l_max,0:l_max,0:2*l_max,0:2*l_max,0:2*l_max)
		real(8), intent(inout) :: CG_tensor(0:l_max,0:l_max,0:l_max,-l_max:l_max,-l_max:l_max,-l_max:l_max)
		
		integer, intent(in) :: l_max
		real(8) :: diff = dble(1e-5)
		
		integer :: l_1,m_1,l_2,m_2,l_,m_

		CG_tensor = 0
		do l_=0,l_max
			do l_1=0,l_max
				do l_2=0,l_max
					do m_=-l_,l_
						do m_1 = -l_1,l_1
							do m_2 = -l_2,l_2
								if (abs(CG_tensor_in(l_,l_1,l_2,m_+l_max,m_1+l_max,m_2+l_max)).lt.diff)then
									CG_tensor(l_,l_1,l_2,m_,m_1,m_2) = 0.0
								else
									CG_tensor(l_,l_1,l_2,m_,m_1,m_2) = CG_tensor_in(l_,l_1,l_2,m_+l_max,m_1+l_max,m_2+l_max)
								end if
							end do
						end do
					end do
				end do
			end do		
		end do		
	end subroutine format_CG
		
	
	
	function get_CG_tensor(l_max)
	!calculate all of the Clebsch Gordan coefficients and store them in an array
		implicit none
		integer, intent(in) :: l_max
		real(8), dimension(0:l_max,0:l_max,0:l_max,0:2*l_max,0:2*l_max,0:2*l_max) :: get_CG_tensor !note that m value doesn't correspond with index, this is to get to python
		!should have dimensions (0:l_max,0:l_max,0:lmax,-l_max:l_max,-l_max:l_max,-l_max:l_max)
		integer :: l_,l_1,l_2,m_,m_1,m_2
		
		do l_=0,l_max
			do l_1=0,l_max
				do l_2=0,l_max
					do m_=-l_,l_
						do m_1 = -l_1,l_1
							do m_2 = -l_2,l_2
								!CG_tensor(l_,l_1,l_2,m_,m_1,m_2) = 0
								get_CG_tensor(l_,l_1,l_2,m_+l_max,m_1+l_max,m_2+l_max) = CG(l_1,m_1,l_2,m_2,l_,m_)
							end do
						end do
					end do
				end do
			end do		
		end do
	end function get_CG_tensor
	
	
	
	real(8) function CG(l_1,m_1,l_2,m_2,l,m)
	!Pretty much all copied from Andrew Fowler's code, should go back and write myself
	!See his github here: https://github.com/andrew31416 and look for spherical_harmonics.f90
	!This is the Clebsch-Gordan coefficient C_{l_1,m_1,l_2,m_2}^{l,m}
	!The formula to calculate this can be found on page 238 of 'Quantum Theory of Angular Momentum' by Varshalovich
	
		implicit none
		
		integer, intent(in) :: l_1,m_1,l_2,m_2,l,m
	
		real(8) :: minimum,min_array(1:7),sqrtres
		real(8) :: imin,imax,val,sumres,sqrtarg
		real(8) :: dble_ii
		integer :: ii
	
		if (abs(m_1 + m_2 - m).gt.1e-15) then
			CG = 0.0d0
		else
			min_array(1) = l_1 + l_2 - l + 0.0d0
			min_array(2) = l_1 - l_2 + l + 0.0d0
			min_array(3) = -l_1 + l_2 + l + 0.0d0
			min_array(4) = l_1 + l_2 + l + 1.0d0
			min_array(5) = l_1 - abs(m_1) + 0.0d0
			min_array(6) = l_2 - abs(m_2) + 0.0d0
			min_array(7) = l - abs(m) + 0.0d0
	
			minimum = minval(min_array)
	
			if (minimum.lt.0.0d0) then
				CG = 0.0d0
			else
				
				sqrtarg = 1.0d0
				sqrtarg = sqrtarg * factorial(minAbsFloor(l_1+m_1+ 0.0d0))
				sqrtarg = sqrtarg * factorial(minAbsFloor(l_1-m_1+ 0.0d0))
				sqrtarg = sqrtarg * factorial(minAbsFloor(l_2+m_2+ 0.0d0))
				sqrtarg = sqrtarg * factorial(minAbsFloor(l_2-m_2+ 0.0d0))
				sqrtarg = sqrtarg * factorial(minAbsFloor(l+m+ 0.0d0))
				sqrtarg = sqrtarg * factorial(minAbsFloor(l-m+ 0.0d0))
				sqrtarg = sqrtarg * dble((int(2.0d0*l) + 1))
				sqrtarg = sqrtarg * factorial(minAbsFloor(min_array(1)))
				sqrtarg = sqrtarg * factorial(minAbsFloor(min_array(2)))
				sqrtarg = sqrtarg * factorial(minAbsFloor(min_array(3)))
				
				! sqrtarg is int so need to divide after casting to double
				sqrtres = sqrt(sqrtarg / factorial(minAbsFloor(min_array(4))))
				
				min_array(1) = l_1 + m_2 - l+ 0.0d0
				min_array(2) = l_2 - m_1 - l+ 0.0d0
				min_array(3) = 0.0d0
				min_array(4) = l_2 + m_2+ 0.0d0
				min_array(5) = l_1 - m_1+ 0.0d0
				min_array(6) = l_1 + l_2 - l+ 0.0d0
				
				imin = maxval(min_array(1:3))
				imax = minval(min_array(4:6))
				sumres = 0.0d0
				do ii=minAbsFloor(imin),minAbsFloor(imax)
					dble_ii = dble(ii)
					val = 1.0d0
					val = val * factorial(ii)
					val = val * factorial(minAbsFloor(l_1 + l_2 - l - dble_ii ))
					val = val * factorial(minAbsFloor(l_1 - m_1 - dble_ii ))
					val = val * factorial(minAbsFloor(l_2 + m_2 - dble_ii ))
					val = val * factorial(minAbsFloor(l - l_2 + m_1 + dble_ii ))
					val = val * factorial(minAbsFloor(l - l_1 - m_2 + dble_ii ))
					sumres = sumres + (-1.0d0)**ii / val
				end do
				CG = sqrtres * sumres
			end if
		end if
	
	end function CG
	
	
	
	function get_output_mask(n_max,l_max)
	!returns an array of logical values so that the redundant values of the bispectrum tensor can be removed as it's flattened
	!conditions as in https://academic.oup.com/mnras/article/318/2/584/1027450
	!section 2 conditions (i) to (iii)
		implicit none
		integer, intent(in) :: n_max
		integer, intent(in) :: l_max
		logical, dimension(1:n_max,0:l_max,0:l_max,0:l_max) :: get_output_mask
		integer :: n,l,l_1,l_2
		
		get_output_mask = .False.
		
		do n=1,n_max
			do l=0,l_max
				do l_1=0,l_max
					do l_2=l_1,l_max !bispectrum is the same if l_1 and l_2 are swapped, so only need one
						!condition (i): |l_1-l_2| <= l <= l_1+l_2
						if (abs(l_1-l_2).le.l.and.l.le.l_1+l_2) then
							if (mod(l+l_1+l_2,2).eq.0) then
							get_output_mask(n,l,l_1,l_2) = .True.
							end if
						end if
					end do
				end do
			end do
		end do
		
		 
	end function get_output_mask
	
	
	
	integer function bispect_length(n_max,l_max)
	!uses the conditions in section 2 of https://academic.oup.com/mnras/article/318/2/584/1027450 
	!figures out the length of the bispectrum array that will be returned
		implicit none
		integer, intent(in) :: n_max
		integer, intent(in) :: l_max
		integer :: counter, n,l,l_1,l_2
		counter = 0
	
		do n=1,n_max
			do l=0,l_max
				do l_1=0,l_max
					do l_2=l_1,l_max !bispectrum is the same if l_1 and l_2 are swapped, so only need one
						!condition (i): |l_1-l_2| <= l <= l_1+l_2
						if (abs(l_1-l_2).le.l.and.l.le.l_1+l_2) then
							if (mod(l+l_1+l_2,2).eq.0) then
								counter = counter+1
							end if
						end if
					end do
				end do
			end do
		end do
		
		bispect_length = counter
	end function bispect_length
	
	
	
	function invOverlap(n_max)
	!gets the inverse of the overlap matrix so that the proper coefficients can be calculated
		implicit none
		integer, intent(in) :: n_max
		real(8), dimension(n_max,n_max) :: invOverlap
		real(8), dimension(n_max,n_max) :: overlap
		integer, dimension(n_max) :: ipiv, work
		integer :: info1, info2
		
		overlap = overlapMatrix(n_max)
		
		
		invOverlap = overlap
		
		call dgetrf(n_max,n_max,invOverlap,n_max,ipiv,info1)
		call dgetri(n_max, invOverlap, n_max, ipiv, work, n_max, info2)
!~ 		call dsytrf('U',n_max,invOverlap,n_max,ipiv,overlap,n_max,info1)
!~ 		call dsytri('U',n_max,invOverlap,n_max,ipiv,overlap,info2)
			
		
		call inverseChecker(overlap,invOverlap,n_max)
	
	end function invOverlap
	
		
		
	subroutine getW(W,n_max)
	!gets the matrix used to combine the phi values to get the radial basis
	!First get the overlap matrix, then do an eigenvalue decomposition S = R L R^-1
	!Where R is a unitary matrix and L is the matrix of eigenvalues
	!S^0.5 = R L^0.5 R^-1 since the square of a matrix has the same eigenvectors, with a squared eigenvalues
	!S^-1 = R L^-1 R^-1
	!Thus S^-0.5 = R L^-0.5 R^-1
		implicit none
		real(8), intent(inout) :: W(:,:)
		integer, intent(in) :: n_max
		
		real(8) :: overlap(n_max,n_max)!, eigenValues(n_max,n_max)
		overlap = overlapMatrix(n_max)
		
		call sqrtInvSymmMatrix(overlap,W,n_max)
	
	end subroutine getW			
	
	
	
end module bispectrum

		
		
		
