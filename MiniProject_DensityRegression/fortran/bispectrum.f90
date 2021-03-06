module bispectrum
	use util
	use basis, only: sphericalHarm, radialBasis, overlapMatrix, getW
	use neighbours
	implicit none
	!bispectrum computed as in PHYSICAL REVIEW B 87, 184115 (2013)
	!To do: check that the thetas and phis are in the right place for the spherical harmonic

	public getLocalBispectrum
	public get_CG_tensor
	
	
	contains
	
	function getPowerSpectrum(n_max,l_max,r_c,atomPositions,cell,nAtoms,pointPosition,W,array_length,local)
	!gets the local or global power spectrum, depending on the logical value local
	!returns a vector array
		implicit none
		integer, intent(in) :: n_max
		integer, intent(in) :: l_max
		real(8), intent(in) :: r_c
		real(8), intent(in) :: atomPositions(:,:)
		real(8), intent(in) :: cell(3,3)
		integer, intent(in) :: nAtoms
		real(8), intent(in) :: pointPosition(3)
		real(8), intent(in) :: W(n_max,n_max)
		integer, intent(in) :: array_length
		logical, intent(in) :: local
		
		type(systemStateType) :: systemState
		type(pointType) :: point
		type(bispectParamType) :: bispectParams
		

		real(8), dimension(1:n_max,0:l_max) :: power_matrix
		real(8), dimension(array_length) :: getPowerSpectrum
		integer :: n,l,m, ii

		
		call assignBispectParams(bispectParams,r_c,l_max,n_max)
		call assignSystemState(systemState,cell,atomPositions)
		call assignPoint(point,pointPosition)
		
		
		if (local) then
			call getLocalPowerSpectrum(systemState,point,bispectParams,power_matrix,W)
		else
			call getGlobalPowerSpectrum(systemState,bispectParams,power_matrix,W)
		end if
		

		getPowerSpectrum = 0.0d0
		getPowerSpectrum = pack(power_matrix,.True.)
		
		do ii = 1, array_length
			if (getPowerSpectrum(ii).lt.dble(1e-18)) then
				getPowerSpectrum(ii) = 0.0d0
			end if
		end do


		
		call tidySystem(systemState)
		call tidyPoint(point)
		
	
	end function 
	
	function getBispectrum(n_max,l_max,r_c, atomPositions,cell,nAtoms,pointPosition,W,CG_tensor_in,array_length,local)
	!gets the local or global bispectrum, depending on the logical value local
	!returns a vector array, with mask selecting the non zero bispectrum tensor elements (see get_output_mask for conditions)
		implicit none
		integer, intent(in) :: n_max
		integer, intent(in) :: l_max
		real(8), intent(in) :: r_c
		real(8), intent(in) :: atomPositions(:,:)
		real(8), intent(in) :: cell(3,3)
		integer, intent(in) :: nAtoms
		real(8), intent(in) :: pointPosition(3)
		real(8), intent(in) :: W(n_max,n_max)
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
		integer :: n,l,l_1,l_2, ii
		!real(8) :: check
		
		mask = get_output_mask(n_max,l_max)	
		
		!check = CG(4+0.0d0,4+0.0d0,4+0.0d0,0+0.0d0,4+0.0d0,4+0.0d0)
		!print *, check
		
		call assignBispectParams(bispectParams,r_c,l_max,n_max)
		call format_CG(CG_tensor_in,CG_tensor,l_max)
		call assignSystemState(systemState,cell,atomPositions)
		call assignPoint(point,pointPosition)
		
		
		if (local) then
			call getLocalBispectrum(systemState,point,bispectParams,bispect_tensor,W,CG_tensor)
		else
			call getGlobalBispectrum(systemState,bispectParams,bispect_tensor,W,CG_tensor)
		end if
		
		!flatten tensor, keeping only the non zero tensor elements
		getBispectrum = 0
		getBispectrum = pack(bispect_tensor,mask)
		
		do ii = 1, array_length
			if (getBispectrum(ii).lt.dble(1e-18)) then
				getBispectrum(ii) = 0.0d0
			end if
		end do
		
		
		call tidySystem(systemState)
		call tidyPoint(point)
		
	
	end function getBispectrum
	
	
	
	subroutine getLocalBispectrum(systemState,point,biPrms,bispect_local,W,CG_tensor)
		implicit none
		type(systemStateType), intent(in) :: systemState
		type(pointType), intent(inout) :: point
		type(bispectParamType), intent(in) :: biPrms
		real(8), intent(in) :: W(biPrms%n_max,biPrms%n_max)
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
		
		call localProjector(systemState,point,coeffs,W,r_c,n_max,l_max)

		if (all(coeffs.eq.complex(0.0d0,0.0d0))) then
			bispect_local = 0
			!print *, "zero coefficients"
		else
		
		bispect_local = 0.0d0

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
	
	
	
	subroutine getGlobalBispectrum(systemState,biPrms,bispect_global,W,CG_tensor)
		implicit none
		type(systemStateType), intent(in) :: systemState
		type(bispectParamType), intent(in) :: biPrms
		real(8), intent(in) :: W(biPrms%n_max,biPrms%n_max)
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
		
		call globalProjector(systemState,coeffs,W,biPrms)
		
		if (all(coeffs.eq.0)) then
			bispect_global = 0.0d0
		else
		bispect_global = 0.0d0
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
	
	subroutine getLocalPowerSpectrum(systemState,point,biPrms,power_local,W)
		implicit none
		type(systemStateType), intent(in) :: systemState
		type(pointType), intent(inout) :: point
		type(bispectParamType), intent(in) :: biPrms
		real(8), intent(in) :: W(biPrms%n_max,biPrms%n_max)
		real(8), intent(inout) :: power_local(1:biPrms%n_max,0:biPrms%l_max)
		
		complex(8), allocatable :: coeffs(:,:,:)
		integer :: n, l, m
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
				
		call localProjector(systemState,point,coeffs,W,r_c,n_max,l_max)
		
		if (all(coeffs.eq.complex(0.0d0,0.0d0))) then
			power_local = 0.0d0
			!print *, "zero coefficients"
		else
		
		power_local = 0.0d0

		do n=1,n_max
			do l=0,l_max
				tmp = 0.0d0
				do m=-l,l_max
					tmp = tmp + RealPart(conjg(coeffs(n,l,m))*coeffs(n,l,m))
				end do
				power_local(n,l) = tmp
			end do
		end do	
		end if
		

				
		deallocate(coeffs)
	
	
	end subroutine getLocalPowerSpectrum
	
	subroutine getGlobalPowerSpectrum(systemState,biPrms,power_global,W)
		implicit none
		type(systemStateType), intent(in) :: systemState
		type(bispectParamType), intent(in) :: biPrms
		real(8), intent(in) :: W(biPrms%n_max,biPrms%n_max)
		real(8), intent(inout) :: power_global(1:biPrms%n_max,0:biPrms%l_max)
		
		complex(8), allocatable :: coeffs(:,:,:)
		integer :: n, l, m
		integer :: n_max, l_max
	
		real(8) :: tmp
		
		n_max = biPrms%n_max
		l_max = biPrms%l_max
		
		allocate(coeffs(1:n_max,0:l_max,-l_max:l_max))
		
		call globalProjector(systemState,coeffs,W,biPrms)
		
		if (all(coeffs.eq.0)) then
			power_global = 0.0d0
		else
		power_global = 0.0d0
		do n=1,n_max
			do l=0,l_max
				tmp = 0.0d0
				do m=-l,l
					tmp = tmp + RealPart(conjg(coeffs(n,l,m))*coeffs(n,l,m))
				end do
				power_global(n,l) = tmp
			end do
		end do	
		
		end if
		
		deallocate(coeffs)
	
	end subroutine getGlobalPowerSpectrum
	
	
	
	subroutine localProjector(systemState,point,coeffs,W,r_c,n_max,l_max)
	!around a general point in space, project atomic positions onto the basis functions
	
		implicit none
		type(systemStateType), intent(in) :: systemState
		type(pointType), intent(inout) :: point
		integer, intent(in) :: n_max
		integer, intent(in) :: l_max
		real(8), intent(in) :: r_c
		real(8), intent(in) :: W(n_max,n_max)
		complex(8), intent(inout) :: coeffs(1:n_max,0:l_max,-l_max:l_max)

		
		complex(8), dimension(1:n_max,0:l_max,-l_max:l_max) :: coeffs_
		real(8), allocatable :: polarPositions(:,:)
		complex(8) :: temp
		integer :: ii, n, n_1, n_2, l, m, buffer_size = 100
!~ 		real(8), dimension(n_max,n_max) :: W
!~ 		real(8), dimension(n_max,n_max) :: inv_S

		coeffs = complex(0.0d0,0.0d0)
		temp = complex(0.0d0,0.0d0)

		call getNeighbours(systemState, point,r_c,buffer_size)
		
		if (point%numNeighbours.gt.0) then
			allocate(polarPositions(3,point%numNeighbours))

!~ 			call getW(W,n_max)
!~ 			inv_S = invOverlap(n_max)
			
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
						coeffs(n, l, m) = temp
						temp = complex(0.0d0,0.0d0)
					end do
					!convert c' to c
!~ 					do n_1=1,n_max
!~ 						do n_2=1,n_max
!~ 							temp = temp + inv_S(n_2,n_1)*coeffs_(n_2,l,m)
!~ 						end do
!~ 						coeffs(n_1,l,m) = temp
!~ 						temp = complex(0.0d0,0.0d0)
!~ 					end do
				end do
			end do

			deallocate(polarPositions)


		else
			coeffs = complex(0.0d0,0.0d0)
			!print *, "no neighbours, projected coefficients set to zero"
		end if
		call tidyPoint(point)
			
		
	end subroutine localProjector
	
	
	
	subroutine globalProjector(systemState,coeffs,W,biPrms)
	!sets the point position to an atom position and removes the atom position from the systemState, then calculates the local bispectrum like that
		implicit none
		
		type(systemStateType), intent(in) :: systemState
		type(bispectParamType), intent(in) :: biPrms
		real(8), intent(in) :: W(biPrms%n_max,biPrms%n_max)
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
		
		tmpAtomPosns = 0.0d0
		tmpCoeffs = complex(0.0d0,0.0d0)
		coeffs = complex(0.0d0,0.0d0)
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
			
			tmpCoeffs = complex(0.0d0,0.0d0)			
			call localProjector(tmpState,point,tmpCoeffs,W,r_c,n_max,l_max)
			
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
		real(8), dimension(0:l_max,0:l_max,0:l_max,0:2*l_max,0:2*l_max,0:2*l_max) :: get_CG_tensor!note that m value doesn't correspond with index, this is to get to python
		!should have dimensions (0:l_max,0:l_max,0:lmax,-l_max:l_max,-l_max:l_max,-l_max:l_max)
		integer :: l_,l_1,l_2,m_,m_1,m_2
		integer(8) :: fact_array(0:l_max*4)
		integer :: cntr = 0
		integer :: tot_cntr = 0
		
		get_CG_tensor = 0.0d0
		call factorial_array(fact_array,l_max*4)
		
		do l_=0,l_max
			do l_1=0,l_max
				do l_2=0,l_max
					do m_=-l_,l_
						do m_1 = -l_1,l_1
							do m_2 = -l_2,l_2
								tot_cntr = tot_cntr + 1
								!CG_tensor(l_,l_1,l_2,m_,m_1,m_2) = 0
								get_CG_tensor(l_,l_1,l_2,m_+l_max,m_1+l_max,m_2+l_max) = CG(l_1+0.0d0,m_1+0.0d0,l_2+0.0d0,&
								m_2+0.0d0,l_+0.0d0,m_+0.0d0,fact_array,l_max)
								if (abs(CG(l_1+0.0d0,m_1+0.0d0,l_2+0.0d0,m_2+0.0d0,l_+0.0d0,m_+0.0d0,fact_array,l_max)).gt.1.0) then
									!print *, CG(l_1+0.0d0,m_1+0.0d0,l_2+0.0d0,m_2+0.0d0,l_+0.0d0,m_+0.0d0,fact_array,l_max)
									print *, "CG error"
								end if
							end do
						end do
					end do
				end do
			end do		
		end do
		
	end function get_CG_tensor
	
	
	

	
	
	
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
		
		
		invOverlap = invSqrtOverlap(n_max)
		invOverlap = matmul(invOverlap,invOverlap)
		
!~ 		call dgetrf(n_max,n_max,invOverlap,n_max,ipiv,info1)
!~ 		call dgetri(n_max, invOverlap, n_max, ipiv, work, n_max, info2)
!~ 		call dsytrf('U',n_max,invOverlap,n_max,ipiv,overlap,n_max,info1)
!~ 		call dsytri('U',n_max,invOverlap,n_max,ipiv,overlap,info2)
			
		
		call inverseChecker(overlap,invOverlap,n_max)
	
	end function invOverlap
	
		
	function invSqrtOverlap(n_max)
	!just something to present getW as a function without changing anything too much
		implicit none
		integer, intent(in) :: n_max
		real(8), dimension(n_max,n_max) :: invSqrtOverlap
		
		call getW(invSqrtOverlap,n_max)
		
	end function invSqrtOverlap
	
	
end module bispectrum

		
		
		
