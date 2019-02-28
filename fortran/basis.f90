module BasisFunctions
	use utilities
	
	implicit none
	
	public sphericalHarm
	public radialBasis
	public CG
	public invOverlap

	contains
	
	real(8) function legendre_0(l,x)
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
		legendre_0 = 1
		
		if (l == 1) then
			legendre_0 = x
		
		elseif (l.gt.1) then
			
			tmp1 = 1! set up temporary variable
			legendre_0 = x
			do ii = 2, l
				tmp2 = tmp1 !update temporary variables and legendre
				tmp1 = legendre_0
				legendre_0 = ((2*ii-1)*x*tmp1 - (ii-1)*tmp2)/ii
				!print *, legendre
			end do
		end if
		
	end function legendre_0

	real(8) function legendre(m,l,x)
	!gets the associated legendre polynomials
	!use the following:
	!P^{m}_{m} = (-1)^{m}(2m-1)!!(1-x^{2})^{m/2}
	!P^{m}_{m+1}  =  x(2m+1)p^{m}_{m}
	!(l-m)P^{m}_{l} = x(2l-1)P^{m}_{l-1} - (l+m-1)p^{m}_{l-2}
		implicit none
		integer, intent(in) :: m
		integer, intent(in) :: l
		real(8), intent(in) :: x
		integer :: i, m_
		real(8) :: tmp1, tmp2, diff, coeff1, coeff2
		real(8) :: factor
		

		if (m.lt.0) then
			m_ = -1*m
			factor = (-1)**m_ * factorial(l-m_)/(factorial(l+m_)+0.0d0)
		else 
			m_ = m
			factor = 1.0d0
		end if
		diff = l-m_
		
		!deal with the cases that don't require the recursion formula before starting
		!l=1 and l=0 fall under this umbrella anyway, which is nice
		if (diff.lt.0) then
			write(*,*) "l must be greater than or equal to |m|"
		else if (diff.eq.0) then
			legendre = p_mm(m_,x)
		else if (diff.eq.1) then
			legendre = p_mm1(m_,x)
		else if (diff.gt.1) then
			legendre = p_mm1(m_,x)
			tmp1 = p_mm(m_,x)
			do i=m_+2,l
				coeff1 = 2*i-1
				coeff2 = i+m_-1
				diff = i-m_
				tmp2 = tmp1
				tmp1 = legendre
				legendre = (x*coeff1*tmp1-coeff2*tmp2)/diff
			end do
		end if
		legendre = factor*legendre
	end function legendre
	
	real(8) function p_mm(m,x)
	!gets the associated legendre polynomial p^{m}_{m}(x)
	!formula: p_mm = (-1)^{m}(2m-1)!!(1-x^{2})^{m/2}
		implicit none
		integer, intent(in) :: m
		real(8), intent(in) :: x
		integer :: oddfact
		
		if (m.gt.0) then
			oddfact = odd_factorial(2*m-1)
			p_mm = (-1)**m * oddfact * (1-x**2)**(m/2.0)
		else
			p_mm=1
		end if
		
	end function p_mm
	
	real(8) function p_mm1(m,x)
	!gets the associated legendre polynomial p^{m}_{m+1}(x)
	!formula: p_mm1 = x(2m+1)p^{m}_{m+1}
		implicit none
		integer, intent(in) :: m
		real(8), intent(in) :: x
		real(8) :: pmm
		
		pmm = p_mm(m,x)
		
		p_mm1 = x*(2*m+1)*pmm
		
	end function p_mm1
		
	
	complex(8) function sphericalHarm(theta,phi,m,l)
		implicit none
		
		integer, intent(in) :: l
		integer, intent(in) :: m
		real(8), intent(in) :: phi
		real(8), intent(in) :: theta
		real(8) :: P_lm, factor
		complex(8) :: comp_exp
		real(8) :: pi = 4*atan(1.0d0)
		
		P_lm = legendre(m,l,cos(theta))
		comp_exp = exp(complex(0,m*phi))
		factor = sqrt(((2*l+1)*factorial(l-m))/((4*pi)*factorial(l+m)))
		
		sphericalHarm = factor*P_lm*comp_exp
		
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
		
		real(8) :: overlap(n_max,n_max)!, eigenValues(n_max,n_max)
		overlap = overlapMatrix(n_max)
		
		call sqrtInvSymmMatrix(overlap,W)
	
	end subroutine getW
	
	function invOverlap(n_max)
	!gets the inverse of the overlap matrix so that the proper coefficients can be calculated
		implicit none
		integer, intent(in) :: n_max
		real(8), dimension(n_max,n_max) :: invOverlap
		real(8), dimension(n_max,n_max) :: overlap
		integer, dimension(n_max) :: ipiv
		integer :: info1, info2
		
		overlap = overlapMatrix(n_max)
		invOverlap = overlap
		overlap = 0
		call dsytrf('U',n_max,invOverlap,n_max,ipiv,overlap,n_max,info1)
		call dsytri('U',n_max,invOverlap,n_max,ipiv,overlap,info2)
	
	end function invOverlap
	
	
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
	
	
end module BasisFunctions	
		
