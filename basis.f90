module BasisFunctions
	use utilities
	
	implicit none
	
	public sphericalHarm
	public radialBasis
	public CG

	contains
	
	real(8) function legendre(l,x)
	! get the legendre polynomial from the iterative scheme:
	! (n+1)P_(n+1) (x) = (2n+1)xP_n(x) - nP_(n-1)(x)
		implicit none
		
		integer, intent(in) :: l
		real(8), intent(in) ::x
		real(8) :: tmp1, tmp2 !legendre polynomial, and temporary storage for the series
		integer :: ii
		
		if (x.lt. -1 .or. x.gt.1) then
			write(*,*) "x is out of the range of Legendre Polynomial"
			call exit(0)
		end if
		
		if (l.lt.0) then
			write(*,*) "l needs to be a non-negative integer"
			call exit(0)
		end if
		
		!start with 0 order legendre
		legendre = 1
		
		if (l == 1) then
			legendre = x
		
		elseif (l.gt.1) then
			
			tmp1 = 1! set up temporary variable
			legendre = x
			do ii = 2, l
				tmp2 = tmp1 !update temporary variables and legendre
				tmp1 = legendre
				legendre = ((2*ii-1)*x*tmp1 - (ii-1)*tmp2)/ii
				!print *, legendre
			end do
		end if
		
	end function legendre
	
	complex(8) function sphericalHarm(theta,phi,m,l)
		implicit none
		
		integer, intent(in) :: l
		integer, intent(in) :: m
		real(8), intent(in) :: phi
		real(8), intent(in) :: theta
		real(8) :: P_l
		complex(8) :: comp_exp
		
		P_l = legendre(l,cos(theta))
		comp_exp = exp(complex(0,m*phi))
		
		sphericalHarm = P_l*comp_exp
		
	end function sphericalHarm
	
	real(8) function radialPhi(r, r_c, alpha)
	!gives the value of the nth radial basis function at r
		implicit none
		
		real(8), intent(in) :: r, r_c
		integer, intent(in) :: alpha
		real(8) :: denom
		
		denom = sqrt(r_c**(2*alpha+5)/(2*alpha+5))
		
		radialPhi = (r_c-r)**(alpha+2)/denom
	end function radialPhi
	
	
	
	real(8) function radialBasis(r, r_c, n, n_max,W)
	!get the g_n(r) basis
		implicit none
		real(8), intent(in) :: r, r_c
		integer, intent(in) :: n, n_max
		real(8), intent(in) :: W(:,:)
		
		real(8) :: phi(n_max)
		integer :: ii
		

		
		do ii = 1, n_max
			phi(ii) = radialPhi(r,r_c,ii)
		end do	
		
		radialBasis = dot(W(n,:),phi)
	
	end function radialBasis
	
	
	
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
		
		real(8) :: overlap(n_max,n_max), eivenValues(n_max,n_max)
		overlap = overlapMatrix(n_max)
		
		call sqrtInvSymmMatrix(overlap,W)
	
	end subroutine getW
	
	
	
	function overlapMatrix(n_max)
	!outputs the overlap matrix of the phi functions
	!the inverse square root of this is useful to get the g_n(r) basis
		implicit none
		integer, intent(in) :: n_max
		
		real(8), dimension(n_max,n_max) :: overlapMatrix
		integer :: ii, jj
		
		do ii = 1, n_max
			do jj = 1, n_max
				overlapMatrix(ii,jj) = sqrt(dble((5+2*ii)*(5+2*jj)))/dble((5+ii+jj))
			end do
		end do
		
	
	end function overlapMatrix
	
	
	
	real(8) function CG(l_1,m_1,l_2,m_2,l,m)
	!Pretty much all copied from Andrew Fowler's code, should go back and write myself
	!See his github here: https://github.com/andrew31416 and look for spherical_harmonics.f90
	!This is the Clebsch-Gordan coefficient C_{l_1,m_1,l_2,m_2}^{l,m}
	!The formula to calculate this can be found on page 238 of 'Quantum Theory of Angular Momentum' by Varshalovich
	
	implicit none
	
	real(8), intent(in) :: l_1,m_1,l_2,m_2,l,m

	real(8) :: minimum,min_array(1:7),sqrtres
	real(8) :: imin,imax,val,sumres,sqrtarg
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
			sqrtarg = sqrtarg * factorial(minAbsFloor(l_1+m_1))
			sqrtarg = sqrtarg * factorial(minAbsFloor(l_1-m_1))
			sqrtarg = sqrtarg * factorial(minAbsFloor(l_2+m_2))
			sqrtarg = sqrtarg * factorial(minAbsFloor(l_2-m_2))
			sqrtarg = sqrtarg * factorial(minAbsFloor(l+m))
			sqrtarg = sqrtarg * factorial(minAbsFloor(l-m))
			sqrtarg = sqrtarg * dble((int(2.0d0*l) + 1))
			sqrtarg = sqrtarg * factorial(minAbsFloor(min_array(1)))
			sqrtarg = sqrtarg * factorial(minAbsFloor(min_array(2)))
			sqrtarg = sqrtarg * factorial(minAbsFloor(min_array(3)))
			
			! sqrtarg is int so need to divide after casting to double
			sqrtres = sqrt(sqrtarg / factorial(minAbsFloor(min_array(4))))
			
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
	
	
end module BasisFunctions	
		
