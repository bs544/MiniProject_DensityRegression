module bispectrum
	use util
	use basis, only: sphericalHarm, radialBasis, invOverlap, getW
	use neighbours
	implicit none
	!bispectrum computed as in PHYSICAL REVIEW B 87, 184115 (2013)

	public getLocalBispectrum
	public get_CG_tensor
	
	
	contains
	
	subroutine getLocalBispectrum(systemState,point,bispectParams,bispect_local,CG_tensor)
		implicit none
		type(systemStateType), intent(in) :: systemState
		type(pointType), intent(inout) :: point
		type(bispectParamType), intent(in) :: bispectParams
		real(8), intent(inout) :: bispect_local(:,:,:,:)
		real(8), intent(in) :: CG_tensor(:,:,:,:,:,:)	
		
		complex(8), allocatable :: coeffs(:,:,:)
		integer :: n, l, l_1, l_2, m, m_1, m_2
		integer :: n_max, l_max
	
		real(8) :: tmp
		
!		logical :: allGood !if this isn't true, then the atomic positions are out of the cell bounds
		
!		do ii = 1, arrayShape(2)
!			allGood = checkInCell(systemState%atomPositions(:,ii),cell)
!		end do
		
		n_max = bispectParams%n_max
		l_max = bispectParams%l_max
		
		allocate(coeffs(1:n_max,0:l_max,-l_max:l_max))
		!allocate(CG_tensor(0:l_max,0:l_max,0:lmax,-l_max:l_max,-l_max:l_max,-l_max:l_max))
		
		!CG_tensor = 0
		!call get_CG_tensor(CG_tensor,l_max)
		
		call localProjector(systemState,point,coeffs,bispectParams)
		
		if (all(coeffs.eq.0)) then
			bispect_local = 0
		else
		
		do n=1,n_max
			do l=0,l_max
				do l_1=0,l_max
					do l_2=0,l_max
						tmp = 0
						do m=-l,l
							do m_1=-l_1,l_1
								do m_2=-l_2,l_2
									tmp = tmp + conjg(coeffs(n,l,m))*CG_tensor(l,l_1,l_2,m,m_1,m_2)*coeffs(n,l_1,m_1)*coeffs(n,l_2,m_2)
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
	
	
	
	subroutine localProjector(systemState,point,coeffs,bispectParams)
	!around a general point in space, project atomic positions onto the basis functions
	
		implicit none
		type(systemStateType), intent(in) :: systemState
		type(pointType), intent(inout) :: point
		type(bispectParamType), intent(in) :: bispectParams
		complex(8), intent(inout) :: coeffs(:,:,:)


		
		complex(8), allocatable :: coeffs_(:,:,:)		
		real(8), allocatable :: polarPositions(:,:)
		complex(8) :: temp
		integer :: ii, n, n_1, n_2, l, m, buffer_size = 100
		integer :: n_max, l_max
		real(8) :: r_c
		real(8), allocatable :: W(:,:)
		real(8), allocatable :: inv_S(:,:)


		coeffs = complex(0.0d0,0.0d0)
		temp = complex(0.0d0,0.0d0)
		r_c = bispectParams%r_c
		n_max = bispectParams%n_max
		l_max = bispectParams%l_max
		

		
		
		call getNeighbours(systemState, point,r_c,buffer_size)
		if (point%numNeighbours.gt.0) then
			allocate(polarPositions(3,point%numNeighbours))
			allocate(coeffs_(1:n_max,0:l_max,-l_max:l_max))
			allocate(W(n_max,n_max))
			allocate(inv_S(n_max,n_max))
			
			call getW(W,n_max)
			inv_S = invOverlap(n_max)
			
			call localCart2Polar(polarPositions, point)
			!polar positions elements are order r, theta, phi
			
			!create W here, it's independent of position, so shouldn't be recomputed all the time
			
			do l = 0,l_max
				do m = -l,l
					do n = 1, n_max
						do ii = 1,point%numNeighbours
							temp = temp + sphericalHarm(polarPositions(2,ii),polarPositions(3,ii),m,l)*radialBasis(polarPositions(1,ii),r_c,n,n_max,W)
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
			deallocate(W)
			deallocate(coeffs_)

		else
			coeffs = 0
		end if
		call tidyPoint(point)
		
	end subroutine localProjector
	
	
	
	subroutine globalProjector(systemState,bispectParams,coeffs)
	!As above, but this is done around an atom, and summed over every atom
		implicit none
		
		type(systemStateType), intent(in) :: systemState
		type(bispectParamType), intent(in) :: bispectParams
		complex(8), intent(inout) :: coeffs(:,:,:)
		
		integer :: ii, n_max, l_max
		real(8) :: atomPosition(3)
		type(pointType) :: point
		complex(8), allocatable :: tempCoeffs(:,:,:)
		
		n_max = bispectParams%n_max
		l_max = bispectParams%l_max
		
		allocate(tempCoeffs(1:n_max,0:l_max,-l_max:l_max))
		
		do ii = 1, systemState%nAtoms
			point%pointPosition = systemState%atomPositions(:,ii)
			call localProjector(systemState,point,tempCoeffs,bispectParams)
			coeffs = coeffs + tempCoeffs
		end do
		
		coeffs = coeffs/systemState%nAtoms
			
	
	end subroutine globalProjector
	
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
	
end module bispectrum

		
		
		
