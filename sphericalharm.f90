module SphericalHarmonics
	use utilities
	
	implicit none

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
	
	real(8) function sphericalHarm(theta,phi,m,l)
		implicit none
		
		integer, intent(in) :: l
		integer, intent(in) :: m
		real(8), intent(in) :: phi
		real(8), intent(in) :: theta
		real(8) :: P_l
		complex(8) :: comp_exp
		
		P_l = legendre(l,cos(theta))
		comp_exp = exp(complex(0,m*phi))
		
		
		
	end function sphericalHarm
	
end module SphericalHarmonics	
		
