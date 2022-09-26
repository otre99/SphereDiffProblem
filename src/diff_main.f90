PROGRAM AdvectionDiffusionProblem
    USE diff_global_data
    USE diff_io
    use diff_numeric_utils, only : throw_exception
    implicit none 
    integer(i4b) :: nlats, nlons, nsave, npar, method, n_species
    real(dp) :: dt, tf, radius 
    logical  :: standard_grid 
    character*256 :: phi0_path,ff0_path,sigma0_path,mulon0_path,&
    mulat0_path,output_folder

    integer(i4b) :: rows, cols, ntotal
    real(dp) :: ttime 

    !!!!!!!!!!!!!!!!!!!!! PARAMETERS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
    namelist /input/ nlats,nlons,n_species,nsave,dt,tf,method,standard_grid,radius,&
    phi0_path,ff0_path,sigma0_path,mulon0_path,mulat0_path,output_folder,npar
    open(10,file='namelist')
    read(10,nml=input)
    close(10)
    !!!!!!!!!!!!!!!!!!!!! END PARAMETERS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

    call init_problem(nlats, nlons, radius, standard_grid, n_species)

    !!! read phi, mu, sigma, forcing
    call read_initial_data(species_(1),phi0_path,ff0_path,sigma0_path,mulon0_path,mulat0_path)

    call print_info_base(method, 0.5*dt, standard_grid)
    call omp_set_num_threads(npar)
    if ( (method.eq.1) .or. (method.eq.2) .or. (method.eq.3)) then 
        call integrate(output_folder, tf, nsave, dt, method)
    else
        write(*,*) "Method must be 1,2, or 3"
        stop
    endif  

    contains
    subroutine read_initial_data(specie, phi_path, ff_path, sigma_path, mulon_path,mulat_path)
        implicit none 
        type (diff_specie), intent(inout) :: specie 
        character(len=*), intent(in) :: phi_path,ff_path,sigma_path,mulon_path,mulat_path 
        if (phi_path .ne. "") then      
            write(*,*) "Reading ", phi_path  
            call read_specie_from_file(phi_path, ttime, specie%phi_np, specie%phi, specie%phi_sp)
        else 
            write(*,*) "Not phi0  file found, so `phi0  = 0`"
        endif 
        
        if (ff_path .ne. "") then 
            write(*,*) "Reading ", ff_path  
            call read_specie_from_file(ff_path, ttime, specie%ff_np, specie%ff, specie%ff_sp)
        else 
            write(*,*) "Not ff0   file found, so `ff0 = 0`"
        endif 
    
        if (sigma_path .ne. "") then 
            write(*,*) "Reading ", sigma_path  
            call read_specie_from_file(sigma_path, ttime, specie%sigma_np, specie%sigma, specie%sigma_sp)
        else 
            write(*,*) "Not sigma0   file found, so `sigma0 = 0`"
        endif 

        if (mulon_path .ne. "") then 
            write(*,*) "Reading ", mulon_path  
            call read_matrix_from_file(mulon_path, ttime, specie%mulon)
        else 
            write(*,*) "Not mulon file found, so `mulon = 0`"
        endif 
    
        if (mulat_path .ne. "") then 
            write(*,*) "Reading ", mulat_path  
            call read_matrix_from_file(mulat_path, ttime, specie%mulat)
        else 
            write(*,*) "Not mulat file found, so `mulat = 0`"
        endif 
    end subroutine
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!! CN-DIM-SPLIT !!!!!!!!!!!!!!!!!!!!!!!!!! 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine integrate(output_folder, tf, nsave, dt, mt)
    use diff_global_data, only : current_iter_, current_time_
    use diff_solver, only: cn_dim_split_linear, cn_dim_split_nonlinear  
    implicit none         
    character(len=*), intent(in) :: output_folder
    integer(i4b), intent(in) :: nsave, mt
    real(dp), intent(in) :: dt, tf     

    real(dp) :: tau    

    current_time_ = 0.0
    tau = 0.5*dt    

    write(*,*) "********************************************************************************************"
    write(*,*) "Current time    :", 0.0, "(initial conditions)"     
    call print_info(species_, tau) 
    call save_sol(output_folder, 0.0_dp, 0, species_)
    write(*,*) "********************************************************************************************"

    do while ( current_time_ .lt. tf) 
        if (mt .eq. 1) then 
            call cn_dim_split_linear(species_, tau)      
        else if (mt .eq. 2) then
            call cn_dim_split_nonlinear(species_, tau)    
        endif                 
        current_iter_ = current_iter_ + 1
        current_time_ = current_iter_*dt

        if ( mod(current_iter_, nsave) .eq. 0 ) then
            write(*,*) "********************************************************************************************"
            write(*,*) "Current time    :", current_time_    
            call print_info(species_, tau) 
            call save_sol(output_folder, current_time_, current_iter_, species_)
            write(*,*) "********************************************************************************************"
        endif
    end do 
    end subroutine integrate
END PROGRAM AdvectionDiffusionProblem


    
