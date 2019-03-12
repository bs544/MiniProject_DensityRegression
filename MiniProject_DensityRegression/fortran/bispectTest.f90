program unitTests

use bispectrum
use util
implicit none

integer, parameter :: n_max=5
integer, parameter :: l_max=4
real(8), parameter :: r_c=1.5
integer, parameter :: nAtoms = 2
integer, parameter :: nIterations = 10
real(8) :: atomPositions(3,nAtoms)
real(8), dimension(3,3) :: cell
real(8), dimension(3) :: pointPosition, shift
real(8) :: W(n_max,n_max)
real(8) :: CG_tensor_in(0:l_max,0:l_max,0:l_max,0:2*l_max,0:2*l_max,0:2*l_max)
integer :: array_length
logical :: local = .False.
real(8), allocatable :: bispect(:,:)
real(8), allocatable :: bi(:)
real(8) :: maxDiff(nIterations)
integer :: ii,jj
complex(8) :: coeffs(1:n_max,0:l_max,-l_max:l_max)
complex(8) :: tot_coeffs(nIterations,1:n_max,0:l_max,-l_max:l_max)
real(8) :: powerSpect(n_max*(l_max+1))
real(8) :: tot_powerSpect(nIterations,n_max*(l_max+1))
type(systemStateType) :: systemvar
type(pointType) :: pointvar

atomPositions = 0.0d0
do ii=1,nAtoms
	do jj=1,3
		atomPositions(jj,ii) = dble(ii-1)/5.0d0 + dble(jj-1)*dble(1e-2)
	end do
end do

shift(1) = 0.01d0
shift(2) = 0.02d0
shift(3) = 0.03d0

!initialise everything
array_length = 0
CG_tensor_in = 0.0d0
W = 0.0d0
maxDiff = 0.0d0

array_length = bispect_length(n_max,l_max)
allocate(bispect(nIterations,array_length))
allocate(bi(array_length))
bispect = 0.0d0
bi = 0.0d0

cell(1,1) = 4.0d0
cell(2,1) = 1.0d0
cell(3,1) = 0.0d0
cell(1,2) = 0.2d0
cell(2,2) = 3.5d0
cell(1,3) = 0.0d0
cell(2,3) = 0.0d0
cell(3,3) = 3.0d0

pointPosition(1) = 0.1d0
pointPosition(2) = 0.2d0
pointPosition(3) = 1.2d0

W = invSqrtOverlap(n_max)
CG_tensor_in = get_CG_tensor(l_max)

!~ do jj = 1,nIterations
!~ 	print *, "Atom Positions"
!~ 	print *, atomPositions
!~ 	print *, "point position"
!~ 	print *, pointPosition
!~ 	print *, "bispectrum"
!~ 	bi = getBispectrum(n_max,l_max,r_c, atomPositions,cell,nAtoms,pointPosition,W,CG_tensor_in,array_length,local)
!~ 	bispect(jj,:) = bi
!~ 	print *, bi
!~ 	call translate(atomPositions,shift,nAtoms)
!~ 	call translate(pointPosition,shift,1)	
!~ end do

!~ do jj = 2,nIterations
!~ 	maxDiff(jj-1) = maxval(bispect(jj,:)-bispect(jj-1,:))
!~ end do


!~ do ii = 1,nIterations
!~ 	print *, maxDiff(ii)
!~ 	print *, bispect(ii,1:9)
!~ end do
!~ deallocate(bispect)
!~ deallocate(bi)

!try power spectrum
array_length = n_max*(l_max+1)
do jj = 1,nIterations
	call assignSystemState(systemvar,cell,atomPositions)
	call assignPoint(pointvar,pointPosition)
	powerSpect = getPowerSpectrum(n_max,l_max,r_c,atomPositions,cell,nAtoms,pointPosition,W,array_length,local)
	tot_powerSpect(jj,:) = powerSpect
	call tidyPoint(pointvar)
	call tidySystem(systemvar)
	call translate(atomPositions,shift,nAtoms)
	call translate(pointPosition,shift,1)
end do

maxDiff = 0.0d0
do ii=1,nIterations-1
	maxDiff(ii) = maxval(tot_powerSpect(ii+1,:)-tot_powerSpect(ii,:))
end do

do ii = 1,nIterations-1
	print *, maxDiff(ii)
	print *, tot_powerSpect(ii,1:9)
end do


contains
	subroutine translate(atomPositions,shift,nAtoms)
	!moves atoms in atomPositions by amout shift
	!expect first index to be atom index
		implicit none
		integer, intent(in) :: nAtoms
		real(8), intent(inout) :: atomPositions(3,nAtoms)
		real(8), intent(in) :: shift(3)
		integer:: ii
		do ii = 1,nAtoms
			atomPositions(:,ii) = atomPositions(:,ii) + shift
		end do
		
	end subroutine translate
	
	subroutine rotate(atomPositions,nAtoms)
	!rotates atoms about z axis
		implicit none
		integer, intent(in) :: nAtoms
		real(8), intent(inout) :: atomPositions(3,nAtoms)
		call apply_rotation_matrix(atomPositions)
	end subroutine rotate
	
   subroutine apply_rotation_matrix(matrix)
		real(8),intent(inout) :: matrix(:,:)

		real(8) :: R(1:3,1:3),alpha
		integer :: ii,jj,kk,dim(1:2)        
		real(8),allocatable :: new_matrix(:,:)

		alpha = 1.0d0

		R = 0.0d0
		R(3,1) = 1.0d0
		R(2,2) = cos(alpha)
		R(1,2) = sin(alpha)
		R(2,3) = -sin(alpha)
		R(1,3) = cos(alpha)
		
		! tranpose
		R(3,1) = 1.0d0
		R(2,2) = cos(alpha)
		R(2,3) = sin(alpha)
		R(1,2) = -sin(alpha)
		R(1,3) = cos(alpha)
		

		dim=shape(matrix)
		allocate(new_matrix(dim(1),dim(2)))
		new_matrix = 0.0d0

		do ii=1,dim(1)
			do jj=1,dim(2)
				do kk=1,3
					new_matrix(ii,jj) = new_matrix(ii,jj) + R(ii,kk)*matrix(kk,jj)
				end do
			end do
		end do
		matrix = new_matrix
	end subroutine 

end program unitTests

