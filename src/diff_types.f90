module diff_types
    integer, parameter :: i4b = selected_int_kind(9)
    integer, parameter :: sp = kind(1.0)
    integer, parameter :: dp = kind(1.0d0)
    real(dp), parameter :: pi=3.141592653589793238462643383279502884197_dp
    real(dp), parameter :: pio2=1.57079632679489661923132169163975144209858_dp
    real(dp), parameter :: twopi=6.283185307179586476925286766559005768394_dp
    real(dp), parameter :: sqrt2=1.41421356237309504880168872420969807856967_dp
    real(dp), parameter :: euler=0.5772156649015328606065120900824024310422_dp
    type diff_specie
        real(dp) :: ff_np, phi_np, sigma_np
        real(dp), dimension(:,:), allocatable :: ff, phi, sigma
        real(dp), dimension(:,:), allocatable :: muLat, muLon 
        real(dp) :: ff_sp, phi_sp, sigma_sp
    end type diff_specie
end module diff_types 
