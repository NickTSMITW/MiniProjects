program LL
implicit none

integer, parameter :: dp=selected_real_kind(15,30) !double precision definiton
integer :: num_pw
real(kind=dp) :: a0,fermi_wavevector,LL_G,LL_prefactor,delta
double precision, dimension(1000) :: LL_matrix


LL_G = 0.0_dp
fermi_wavevector = 1.0_dp
a0 = 1.0_dp
delta = 0.4_dp

LL_prefactor = 2.0_dp / (fermi_wavevector * a0 * 3.14159_dp)

print*, LL_prefactor

do num_pw=1,1000

LL_G = LL_G + 0.007

LL_matrix(num_pw) = 1 + LL_prefactor * ( 1/(LL_G**2) - (delta/(2*(LL_G**3)))*( ATAN( (2*LL_G - (LL_G)**2) / delta ) &
                                   + ATAN( (2.0_dp*LL_G + (LL_G)**2) / delta )  ) + ( delta**2/(8*LL_G**5) +        &
                                   1 /(2.0_dp*LL_G**3) - 1/(8*LL_G) ) * ( log( delta**2 + (2*LL_G + LL_G**2)**2 )   &
                                   - log( delta**2 + (2*LL_G - LL_G**2)**2 ) ) )
end do


open(unit=20, file="LL_matrix.dat")

LL_G = 0

do num_pw=1,1000

LL_G = LL_G + 0.007

!print*, LL_matrix(num_pw), LL_G
write(20,*) LL_matrix(num_pw), LL_G

end do

end program
