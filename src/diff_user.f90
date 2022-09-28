module diff_user
use diff_types

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!   UPDATE MU                   !!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine update_mu(sp, iter, tt)
    implicit none 
    type (diff_specie), dimension(:), intent(inout) :: sp 
    real(dp), intent(in) :: tt 
    integer(i4b), intent(in) :: iter 

    call update_mulon(sp, iter, tt)
    call update_mulat(sp, iter, tt)
        
end subroutine update_mu

subroutine update_mulon(sp, iter, tt)
    use diff_global_data, only : R_,lonsB_,lats_
    implicit none 
    type (diff_specie), dimension(:), intent(inout) :: sp 
    real(dp), intent(in) :: tt 
    integer(i4b), intent(in) :: iter 
    integer(i4b) :: k, nlats, nlons, nspecies 

    nspecies = size(sp)
    nlats = size(sp(1)%phi,1)
    nlons = size(sp(1)%phi,2)

    !$omp parallel do  
    do k=1, nlats
        !call phi2mulongrid(sp(1)%phi(k,:), sp(1)%mulon(k,:))
        !sp(1)%mulon(k,:) =  kk*(sp(1)%mulon(k,:))**alpha 
    end do        
    !$omp end parallel do    
    !write(*,*) "MULON", maxval(sp(1)%mulon)
end subroutine update_mulon

subroutine update_mulat(sp, iter, tt)
    use diff_global_data, only : R_,lons_,latsB_
    implicit none 
    type (diff_specie), dimension(:), intent(inout) :: sp 
    real(dp), intent(in) :: tt 
    integer(i4b), intent(in) :: iter 

    integer(i4b) :: k, nlats, nlons, nspecies 

    nspecies = size(sp)
    nlats = size(sp(1)%phi,1)
    nlons = size(sp(1)%phi,2)

    !$omp parallel do  
    do k=1, nlons
        !call phi2mulatgrid(sp(1)%phi_np, sp(1)%phi(:,k), sp(1)%phi_sp, sp(1)%mulat(:,k))
        !sp(1)%mulat(:,k) = kk*(sp(1)%mulat(:,k))**alpha
    end do        
    !$omp end parallel do    

end subroutine update_mulat



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!   FF and SIGMA                !!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine update_ffsigma(sp, iter, tt)
    use diff_global_data, only : R_,lons_,lats_, dlat_
    implicit none 
    type (diff_specie), dimension(:), intent(inout) :: sp 
    real(dp), intent(in) :: tt 
    integer(i4b), intent(in) :: iter 

    integer(i4b) :: k, nlats, nlons, nspecies 

    nspecies = size(sp)

    !! update north pole 
    !sp(1)%sigma_np = 
    !sp(1)%ff_np =   
    
    !$omp parallel do  
    do k=1, nlons    
        !sp(1)%sigma(:,k) =  
        !sp(1)%ff(:,k) = 
    end do        
    !$omp end parallel do    

    !! update south pole       
    !sp(1)%sigma_sp =  
    !sp(1)%ff_sp = 
   

end subroutine update_ffsigma

subroutine dpdt(tt, lat, lon, ff, sigma, sp_in, sp_out)
    implicit none 
    real(dp), dimension(:), intent(in) :: ff, sigma, sp_in  
    real(dp), dimension(:), intent(out) :: sp_out 
    real(dp), intent(in) :: tt, lat, lon 
   
    integer(i4b) :: nspecies 
    nspecies = size(sp_in,1)

    sp_out = ff-sigma*sp_in    
end subroutine dpdt



subroutine update_all(sp, iter, ttime)
    implicit none    
    integer(i4b), intent(in) :: iter    
    type (diff_specie), dimension(:), intent(inout) :: sp 
    real(dp), intent(in) :: ttime

    !call update_mu(sp, iter, ttime)
    call update_ffsigma(sp, iter, ttime)

end subroutine update_all    



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!! AUX ROUTINES !!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine phi2mulongrid(phi, mu)
    implicit none    
    real(dp), dimension(:), intent(in) :: phi
    real(dp), dimension(:), intent(inout) :: mu
    integer(i4b) :: n
    n = size(phi)
    mu(1) = 0.5*(phi(1)+phi(n))
    mu(2:) = 0.5*( phi(1:n-1)+phi(2:) )
end subroutine phi2mulongrid

subroutine phi2mulatgrid(np, phi, sp, mu)
    implicit none    
    real(dp), dimension(:), intent(in) :: phi
    real(dp), dimension(:), intent(inout) :: mu
    real(dp), intent(in) :: np, sp
    integer(i4b) :: n
    n = size(phi)
    mu(1) = 0.5*(phi(1)+np)
    mu(n+1) = 0.5*(phi(n)+sp)
    mu(2:n) = 0.5*( phi(1:n-1)+phi(2:) )
end subroutine phi2mulatgrid


end module diff_user 
