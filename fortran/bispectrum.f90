module bispectrum
	use utilities
	use BasisFunctions, only: sphericalHarm, radialBasis, CG
	use neighbours
	implicit none
	!bispectrum computed as in PHYSICAL REVIEW B 87, 184115 (2013)

	
	
	contains
	
	subroutine getLocalBispectrum(systemState,point,bispectParams)
		implicit none
		type(systemStateType), intent(in) :: systemState
		type(pointType), intent(inout) :: point
		type(bispectParamType), intent(in) :: bispectParams
		
		complex(8), allocatable :: coeffs(:,:,:)
		integer :: n, n_1, n_2, l, l_1, l_2, m, m_1, m_2
		integer :: n_max, l_max
		real(8), allocatable :: W(:,:)
		
!		logical :: allGood !if this isn't true, then the atomic positions are out of the cell bounds
		
!		do ii = 1, arrayShape(2)
!			allGood = checkInCell(systemState%atomPositions(:,ii),cell)
!		end do
		
		n_max = bispectParams%n_max
		l_max = bispectParams%l_max
		
		allocate(W(1:n_max,1:n_max))
		call getW(W,n_max)
		
		allocate(coeffs(1:n_max,0:l_max,-l_max:l_max))
		
		call localProjector(systemState,point,coeffs,bispectParams)
		
		
		
		
	
	end subroutine getLocalBispectrum
	
	
	
	subroutine localProjector(systemState,point,coeffs,bispectParams)
	!around a general point in space, project atomic positions onto the basis functions
	
		implicit none
		type(systemStateType), intent(in) :: systemState
		type(pointType), intent(inout) :: point
		type(bispectParamType), intent(in) :: bispectParams
		complex(8), intent(inout) :: coeffs(:,:,:)

		
		
		real(8), allocatable :: polarPositions(:,:)
		complex(8) :: temp
		integer :: ii, n, l, m, buffer_size = 100
		integer :: n_max, l_max
		real(8) :: r_c


		coeffs = 0
		temp = 0
		r_c = bispectParams%r_c
		n_max = bispectParams%n_max
		l_max = bispectParams%l_max
		

		
		
		call getNeighbours(systemState, point,r_c,buffer_size)
		
		allocate(polarPositions(3,point%numNeighbours))
		
		call localCart2Polar(polarPositions, point)
		!polar positions elements are order r, theta, phi
		
		!create W here, it's independent of position, so shouldn't be recomputed all the time
		
		do n = 1, n_max
			do l = 0,l_max
				do m = -l,l
					do ii = 1,point%numNeighbours
						temp = temp + sphericalHarm(polarPositions(2,ii),polarPositions(3,ii),m,l)*radialBasis(polarPositions(1,ii),r_c,n)
					end do
					coeffs(n, l, m) = temp
					temp = 0
				end do
			end do
		end do
		
		deallocate(polarPositions)
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
	
	
end module bispectrum

		
		
		
