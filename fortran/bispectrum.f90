module bispectrum
	use utilities
	use BasisFunctions, only: sphericalHarm, radialBasis, CG, invOverlap
	use neighbours
	implicit none
	!bispectrum computed as in PHYSICAL REVIEW B 87, 184115 (2013)

	
	
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
	
	subroutine get_CG_tensor(CG_tensor,l_max)
	!calculate all of the Clebsch Gordan coefficients and store them in an array
		implicit none
		real(8), intent(inout) :: CG_tensor(:,:,:,:,:,:)
		integer, intent(in) :: l_max
		!should have dimensions (0:l_max,0:l_max,0:lmax,-l_max:l_max,-l_max:l_max,-l_max:l_max)
		integer :: l_,l_1,l_2,m_,m_1,m_2
		
		do l_=0,l_max
			do l_1=0,l_max
				do l_2=0,l_max
					do m_=-l_,l_
						do m_1 = -l_1,l_1
							do m_2 = -l_2,l_2
								CG_tensor(l_,l_1,l_2,m_,m_1,m_2) = CG(l_1,m_1,l_2,m_2,l_,m_)
							end do
						end do
					end do
				end do
			end do		
		end do
	end subroutine get_CG_tensor
	
	
end module bispectrum

		
		
		
