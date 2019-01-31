module bispectrum
	use utilities
	use SphericalHarmonics, only: sphericalHarm
	implicit none
	!bispectrum computed as in PHYSICAL REVIEW B 87, 184115 (2013)
	
	!derived type for parameters
	type param
		real(dp) :: r_c
		integer :: l_max, n_max
	
	
	real, parameter :: r_c !cutoff radius
	
	
	!around a general point in space, project atomic positions onto the basis functions
	subroutine localProjector()
	!isolate local space
	!go to polar coordinates from central point
	!get spherical harmonic values for each atom and multiply by radial basis and sum
	
	
	end subroutine localProjector
	
	
	!As above, but this is done around an atom, and summed over every atom
	subroutine globalProjector()
	
	
	end subroutine globalProjector
	
	
	
	function cutoff(r, r_c)
		implicit none
		real, parameter :: r_c
		real :: r
		
		
		
