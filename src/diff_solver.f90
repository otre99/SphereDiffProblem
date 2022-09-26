module diff_solver
use diff_types

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine reduce_by_kfactor(phi, kfactor)
    implicit none 
    real(dp), dimension(:), intent(inout) :: phi 
    integer(i4b), intent(in) :: kfactor
    integer(i4b) :: n, i, i1
    n = size(phi)
    i1=1
    do i=1, n/kfactor
        phi(i) = sum(phi(i1:i1+kfactor-1))/kfactor 
        i1=i1+kfactor
    enddo
end subroutine reduce_by_kfactor

subroutine expand_by_kfactor(phi, kfactor)
    implicit none 
    real(dp), dimension(:), intent(inout) :: phi 
    integer(i4b), intent(in) :: kfactor
    integer(i4b) :: n, i, i1
    n = size(phi)
    i1=n-kfactor+1
    do i=n/kfactor, 1, -1 
        phi(i1:i1+kfactor-1) = phi(i)  
        i1=i1-kfactor
    enddo
end subroutine expand_by_kfactor
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine cn_dim_split_linear(sp, tau)
    use diff_global_data, only : I_,J_,current_iter_, current_time_,R_,&
    sinlatsC_,sinlatsB_,kfactor_by_lats_,lats_,lons_
    use diff_cn_dim_split, only : solve_lon_problem_cn,  solve_lat_problem_cn
    use diff_user, only : update_all, dpdt
    implicit none 
    type (diff_specie), dimension(:) :: sp 

    real(dp), dimension(size(sp)) :: sol, k2, sigma, ff          
    integer(i4b) :: i1,i2  

    real(dp), intent(in) :: tau 
    integer(i4b) :: n_species, i, k  


    n_species = size(sp)   

    !!!!!!!!!!!!! Update (u,v) and (mulon, mulat) on t=tn+tau
    call update_all(sp, current_iter_, current_time_+tau)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!! AdvDiff Lat direction [tn-1, tn] !!!!!!!!!!!!!! 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do i=1, n_species        
        call solve_lat_problem_cn(sp(i)%phi_np, sp(i)%phi, sp(i)%phi_sp, sp(i)%mulat,&
        sinlatsC_, sinlatsB_, R_, tau)
    end do

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!! AdvDiff Lon direction [tn-1, tn] !!!!!!!!!!!!!! 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do i=1, n_species        
        !$omp parallel do  
        do k=1, J_
            if (kfactor_by_lats_(k) .ne. 1) then 
                call reduce_by_kfactor(sp(i)%phi(k,:), kfactor_by_lats_(k))
                call solve_lon_problem_cn( sp(i)%phi(k,:I_/kfactor_by_lats_(k)),&
                sp(i)%mulon(k,::kfactor_by_lats_(k)),R_, sinlatsC_(k), tau)
                call expand_by_kfactor(sp(i)%phi(k,:), kfactor_by_lats_(k))                
            else                 
                call solve_lon_problem_cn( sp(i)%phi(k,:), sp(i)%mulon(k,:), R_, sinlatsC_(k), tau)
            endif 
        end do        
        !$omp end parallel do    
    end do    

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!! SIGMA and FORCING [tn-1, tn+1]   !!!!!!!!!!!!!! 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do i=1, n_species        
        sp(i)%phi_np = (2*tau*sp(i)%ff_np+(1.0-tau*sp(i)%sigma_np)*sp(i)%phi_np)/(1+tau*sp(i)%sigma_np)        
        !$omp parallel do  
        do k=1, I_
            sp(i)%phi(:,k) = (2*tau*sp(i)%ff(:,k)+(1.0-tau*sp(i)%sigma(:,k))*sp(i)%phi(:,k))&
            / (1.0+tau*sp(i)%sigma(:,k))        
        end do 
        !$omp end parallel do    
        sp(i)%phi_sp = (2*tau*sp(i)%ff_sp+(1.0-tau*sp(i)%sigma_sp)*sp(i)%phi_sp)/(1+tau*sp(i)%sigma_sp)                
    end do

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!! AdvDiff Lon direction [tn, tn+1] !!!!!!!!!!!!!! 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do i=1, n_species        
        !$omp parallel do  
        do k=1, J_
            if (kfactor_by_lats_(k) .ne. 1) then 
                call reduce_by_kfactor(sp(i)%phi(k,:), kfactor_by_lats_(k))
                call solve_lon_problem_cn( sp(i)%phi(k,:I_/kfactor_by_lats_(k)),&
                sp(i)%mulon(k,::kfactor_by_lats_(k)),R_, sinlatsC_(k), tau)
                call expand_by_kfactor(sp(i)%phi(k,:), kfactor_by_lats_(k))                
            else                 
                call solve_lon_problem_cn( sp(i)%phi(k,:), sp(i)%mulon(k,:), R_, sinlatsC_(k), tau)
            endif 
        end do        
        !$omp end parallel do    
    end do
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!! AdvDiff Lat direction [tn, tn+1] !!!!!!!!!!!!!! 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do i=1, n_species        
        call solve_lat_problem_cn(sp(i)%phi_np, sp(i)%phi, sp(i)%phi_sp, sp(i)%mulat,&
        sinlatsC_, sinlatsB_, R_, tau)
    end do

    do i=1, n_species        
        !$omp parallel do  
        do k=1, J_
            if (kfactor_by_lats_(k) .ne. 1) then 
                call reduce_by_kfactor(sp(i)%phi(k,:), kfactor_by_lats_(k))
                call expand_by_kfactor(sp(i)%phi(k,:), kfactor_by_lats_(k))                
            endif 
        end do        
        !$omp end parallel do    
     end do    
     
end subroutine cn_dim_split_linear


subroutine cn_dim_split_nonlinear(sp, tau)
    use diff_global_data, only : I_,J_,current_iter_, current_time_,R_,&
    sinlatsC_,sinlatsB_,kfactor_by_lats_,lats_, lons_
    use diff_cn_dim_split, only : solve_lon_problem_cn,  solve_lat_problem_cn
    use diff_user, only : update_mulat, update_mulon, dpdt,update_all
    implicit none 
    type (diff_specie), dimension(:) :: sp 
    real(dp), intent(in) :: tau 

    real(dp), dimension(size(sp)) :: sol, k2, sigma, ff      
    integer(i4b) :: n_species, i, i1, i2, k  
    n_species = size(sp)   

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!! AdvDiff Lon direction [tn-1, tn] !!!!!!!!!!!!!! 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    call update_mulon(sp, current_iter_, current_time_)  
    do i=1, n_species        
        !$omp parallel do  
        do k=1, J_
            if (kfactor_by_lats_(k) .ne. 1) then 
                call reduce_by_kfactor(sp(i)%phi(k,:), kfactor_by_lats_(k))
                call solve_lon_problem_cn( sp(i)%phi(k,:I_/kfactor_by_lats_(k)),&
                sp(i)%mulon(k,::kfactor_by_lats_(k)), R_, sinlatsC_(k), tau)
                call expand_by_kfactor(sp(i)%phi(k,:), kfactor_by_lats_(k))                
            else                 
                call solve_lon_problem_cn( sp(i)%phi(k,:), sp(i)%mulon(k,:), R_, sinlatsC_(k), tau)
            endif 
        end do        
        !$omp end parallel do    
    end do 

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!! AdvDiff Lat direction [tn-1, tn] !!!!!!!!!!!!!! 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    call update_mulat(sp, current_iter_, current_time_)    
    do i=1, n_species        
        call solve_lat_problem_cn(sp(i)%phi_np, sp(i)%phi, sp(i)%phi_sp, sp(i)%mulat,&
        sinlatsC_, sinlatsB_, R_, tau)
    end do
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!! SIGMA and FORCING [tn-1,tn+1]    !!!!!!!!!!!!!!! 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do i=1, n_species
        sol(i) = sp(i)%phi_np
        sigma(i) = sp(i)%sigma_np
        ff(i) = sp(i)%ff_np
    enddo
    call dpdt(current_time_, 0.0_dp, 0.0_dp, ff, sigma, sol, k2)
    sol = sol + tau*k2
    call dpdt(current_time_+tau, 0.0_dp, 0.0_dp, ff, sigma, sol, k2)
    do i=1, n_species
        sp(i)%phi_np = sp(i)%phi_np + 2*tau*k2(i)
    enddo
    !$omp parallel do private(k, sol, k2, sigma, ff, i, i1, i2)
    do k=1, I_*J_
        i1 = (k-1)/I_ + 1
        i2 = mod(k-1,I_)+1
        do i=1, n_species
            sol(i) = sp(i)%phi(i1,i2)
            sigma(i) = sp(i)%sigma(i1,i2)
            ff(i) = sp(i)%ff(i1,i2)
        enddo
        call dpdt(current_time_, lats_(i1), lons_(i2), ff, sigma, sol, k2)
        sol = sol + tau*k2
        call dpdt(current_time_+tau, lats_(i1), lons_(i2), ff, sigma, sol, k2)
        do i=1, n_species
            sp(i)%phi(i1,i2) = sp(i)%phi(i1,i2) + 2*tau*k2(i)
        enddo
    end do 
    !$omp end parallel do
    do i=1, n_species
        sol(i) = sp(i)%phi_sp
        sigma(i) = sp(i)%sigma_sp
        ff(i) = sp(i)%ff_sp
    enddo
    call dpdt(current_time_, pi, 0.0_dp, ff, sigma, sol, k2)
    sol = sol + tau*k2
    call dpdt(current_time_+tau, pi, 0.0_dp, ff, sigma, sol, k2)
    do i=1, n_species
        sp(i)%phi_sp = sp(i)%phi_sp + 2*tau*k2(i)
    enddo

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!! AdvDiff Lat direction [tn, tn+1] !!!!!!!!!!!!!! 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    call update_mulat(sp, current_iter_, current_time_+tau)    
    do i=1, n_species        
        call solve_lat_problem_cn(sp(i)%phi_np, sp(i)%phi, sp(i)%phi_sp, sp(i)%mulat,&
        sinlatsC_, sinlatsB_, R_, tau)
    end do

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!! AdvDiff Lon direction [tn, tn+1] !!!!!!!!!!!!!! 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    call update_mulon(sp, current_iter_, current_time_+tau)  
    do i=1, n_species        
        !$omp parallel do  
        do k=1, J_
            if (kfactor_by_lats_(k) .ne. 1) then 
                call reduce_by_kfactor(sp(i)%phi(k,:), kfactor_by_lats_(k))
                call solve_lon_problem_cn( sp(i)%phi(k,:I_/kfactor_by_lats_(k)),&
                sp(i)%mulon(k,::kfactor_by_lats_(k)),R_, sinlatsC_(k), tau)
                call expand_by_kfactor(sp(i)%phi(k,:), kfactor_by_lats_(k))                
            else                 
                call solve_lon_problem_cn( sp(i)%phi(k,:), sp(i)%mulon(k,:), R_, sinlatsC_(k), tau)
            endif 
        end do        
        !$omp end parallel do    
    end do      
end subroutine cn_dim_split_nonlinear


subroutine show_mat(A)
    USE diff_types 
    implicit none
    REAL(DP), DIMENSION(:,:), INTENT(IN) :: A
    INTEGER :: row, col
    DO row=1,size(A, 1)
        write(*,"(100F14.9 )") ( A(row,col), col=1,size(A, 2) )
    ENDDO
end subroutine show_mat

subroutine show_vector(v)
    USE diff_types 
    implicit none
    REAL(DP), DIMENSION(:), INTENT(IN) :: v
    INTEGER :: i, n
    n = size(v)
    write(*,"(100F14.9 )") (v(i), i=1, n)

end subroutine show_vector



end module diff_solver 
