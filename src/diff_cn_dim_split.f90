module diff_cn_dim_split
use diff_types

contains
subroutine solve_lon_problem_cn(phi, mu, R, sinlat, tau)
    use diff_global_data, only : dlon_

    use diff_linalg, only : cyclic
    implicit none 
    real(dp), dimension(:), intent(inout) :: phi
    real(dp), dimension(:), intent(in) :: mu!, vel
    real(dp), intent(in) :: R, sinlat, tau

    integer(i4b) :: n, k
    real(dp) :: dlon, sf, fact
    real(dp), dimension(size(phi)) :: a, b, c, phi_b    

    n = size(phi)
    dlon = twopi/n


    sf = R * dlon * sinlat
    fact = -tau / (4 * sf * sf)
    k=1
    a(k) =  2 * mu(k) * fact
    b(k) = (-2.0 * ( mu(k) + mu(k + 1) ) ) * fact
    c(k) =  2 * mu(k + 1) * fact

    phi_b(k) = -phi(n) * a(k) + (1.0 - b(k)) * phi(k) - phi(k + 1) * c(k)
    b(k) = b(k) + 1.0
    do k=2, n-1
        a(k) =  2 * mu(k) * fact
        b(k) = (-2.0 * ( mu(k) + mu(k + 1) ) ) * fact
        c(k) =  2 * mu(k + 1) * fact

        phi_b(k) = -phi(k - 1) * a(k) + (1.0 - b(k)) * phi(k) - phi(k + 1) * c(k)
        b(k) = b(k) + 1.0
    enddo
    k=n
    a(k) =  2 * mu(k) * fact
    b(k) = (-2.0 * ( mu(k) + mu(1) ) ) * fact
    c(k) =  2 * mu(1) * fact
  
    phi_b(k) = -phi(k - 1) * a(k) + (1.0 - b(k)) * phi(k) - phi(1) * c(k)
    b(k) = b(k) + 1.0

    call cyclic(a, b, c, c(n), a(1), phi_b, phi)
end subroutine solve_lon_problem_cn

subroutine solve_lat_problem_cn1(north,phi,south,mu,sinC,sinB,R,fact,dlat,g,h)
    use diff_linalg, only : tridag3
    implicit none 
    real(dp), dimension(:), intent(inout) :: phi
    real(dp), dimension(:), intent(in) :: mu
    real(dp), dimension(:), intent(in) :: sinC, sinB
    real(dp), dimension(:), intent(inout) :: g, h
    real(dp), intent(in) :: fact,north,south,dlat,R 

    integer(i4b) :: k, n    
    real(dp), dimension(size(phi)) :: a, b, c, phi_b    

    n = size(phi)
    k  = 1
    a(k) = fact * (- 2 * mu(k) * sinB(k)) / sinC(k)
    b(k) = 2*fact * (mu(k) * sinB(k) + mu(k + 1) * sinB(k + 1)) / sinC(k);
    c(k) = fact * (- 2 * mu(k + 1) * sinB(k + 1)) /sinC(k)    
    phi_b(k) = -a(k) * north + (1.0 - b(k)) * phi(k) - c(k) * phi(k + 1)
    b(k) = b(k) + 1.0
    g(k) = 0.0
    do k=2,n-1
        a(k) = fact * (- 2 * mu(k) * sinB(k)) / sinC(k);
        b(k) = 2 * fact * (mu(k) * sinB(k) + mu(k+1) * sinB(k + 1)) / sinC(k);
        c(k) = fact * (- 2 * mu(k+1) * sinB(k + 1)) / sinC(k)
        phi_b(k) = -a(k) * phi(k-1) + (1.0 - b(k)) * phi(k) - c(k) * phi(k+1)
        b(k) = b(k) + 1.0
        h(k) = 0.0
        g(k) = 0.0        
    end do
    k=n
    h(k) = 0.0
    a(k) = fact * ( - 2 * mu(k) * sinB(k)) / sinC(k)
    b(k) = 2 * fact * (mu(k) * sinB(k) + mu(k+1) * sinB(k+1)) / sinC(k)
    c(k) = fact * ( - 2 * mu(k + 1) * sinB(k + 1)) / sinC(k)
    phi_b(k) = -a(k) * phi(k-1) + (1.0 - b(k)) * phi(k) - c(k) * south
    b(k) = b(k)+1.0

    h(1) = a(1)
    g(n) = c(n)
    call tridag3(a(2:), b, c(:n-1), h, g, phi_b, h, g, phi)
   
end subroutine solve_lat_problem_cn1    

subroutine solve_lat_problem_cn(north,phi,south,mu,sinC,sinB,R,tau)
    implicit none 
    real(dp), dimension(:,:), intent(inout) :: phi
    real(dp), intent(inout) :: north, south
    real(dp), dimension(:,:), intent(in) :: mu
    real(dp), dimension(:), intent(in) :: sinC, sinB
    real(dp), intent(in) :: R, tau

    real(dp) :: dlat, fact, dem, dem1, b1, b2, F1, F2
    real(dp) :: Xp, Xh, Xg, Up, Uh, Ug, dd, dd1, dd2 
    integer(i4b) :: nlats, nlons, k, i1, i2
    real(dp), dimension(size(phi)) :: g, h
    real(dp), dimension(size(phi,2)) :: X, U
    logical :: nonlinear

    nlats = size(phi,1)
    nlons = size(phi,2)
    dlat = pi / (nlats + 1)

    fact = 0.5 * tau / (2 * R * R * dlat * dlat)
    dem = nlons * (R * R * dlat * dlat)
    dem1 = 0.5 * tau / dem

    if (.false.) then 
        b1 = 4.0 * sum(mu(1,:)) * dem1
        b2 = 4.0 * sum(mu(nlats+1,:)) * dem1
        X = ( - 4 * mu(1,:)      ) * dem1
        U = ( - 4 * mu(nlats+1,:)) * dem1
    else 
        b1 = 8.0 * sinB(1) * sum(mu(1,:)) * dem1 / dlat
        b2 = 8.0 * sinB(1) * sum(mu(nlats+1,:)) * dem1 / dlat
        X = sinB(1) * ( - 8 * mu(1,:)) * dem1 / dlat
        U = sinB(1) * ( - 8 * mu(nlats+1,:)) * dem1 / dlat
    endif 

    F1 = (1.0 - b1) * north - dot_product(X, phi(1,:))
    b1 = b1 + 1.0
    F2 = (1.0 - b2) * south - dot_product(U, phi(nlats,:))
    b2 = b2 + 1.0
   
    Xp = 0.0
    Xh = 0.0
    Xg = 0.0
    Up = 0.0
    Uh = 0.0
    Ug = 0.0

    !$omp parallel do private(k,i1,i2),reduction(+:xg,xh,xp,ug,uh,up)     
    do k=1, nlons
        i1 = 1+(k-1)*nlats
        i2 = k*nlats

        call solve_lat_problem_cn1(north,phi(:,k),south,mu(:,k),sinC,sinB,R,fact,dlat,g(i1:i2),h(i1:i2))

        xg = xg + g(i1)*X(k)
        xh = xh + h(i1)*X(k)
        xp = xp + phi(1, k)*X(k)

        ug = ug + g(i2)*U(k)
        uh = uh + h(i2)*U(k)
        up = up + phi(nlats, k)*U(k)
    end do
    !$omp end parallel do

    dd = (Ug - b2) * (Xh - b1) - Uh * Xg;
    dd1 = Xg * (F2 - Up) + (F1 - Xp) * (b2 - Ug)
    dd2 = Uh * (F1 - Xp) + (Up - F2) * (Xh - b1)

    north = dd1 / dd  ! ec. (1.86)
    south = dd2 / dd  ! ec. (1.87)

    !$omp parallel do private(k, i1, i2)
    do k=1, nlons
        i1 = 1+(k-1)*nlats
        i2 = k*nlats        
        phi(:,k) = phi(:,k) - north * h(i1:i2) - south * g(i1:i2)
    end do
    !$omp end parallel do   
end subroutine solve_lat_problem_cn


end module diff_cn_dim_split 
