module diff_global_data
use diff_types


integer(i4b) :: I_, J_, total_cells_, vcells_
integer(i4b) :: nonpole_cells_  
real(dp) :: dlat_, dlon_, R_, current_time_
integer(i4b) :: current_iter_ 
real(dp), dimension(:), allocatable :: lons_, lats_, latsB_, lonsB_
real(dp), dimension(:), allocatable :: sinlatsC_, sinlatsB_ 

integer(i4b), dimension(:), allocatable :: kfactor_by_lats_
type (diff_specie), dimension(:), allocatable :: species_

contains
recursive function next_factor(n, k) result(r)
    implicit none 
    integer(i4b), intent(in) :: n
    integer(i4b) :: r, k
    k=k+1
    if ( mod(n,k) .eq. 0 ) then 
        r = k
    else
        r = next_factor(n,k)
    endif 
end function next_factor

subroutine init_problem(nlats, nlons, radius, standard_grid, n_species)
    use diff_types 
    use diff_numeric_utils, only : arange    
    implicit none
    real(dp), intent(in) :: radius
    integer(i4b), intent(in) :: nlats, nlons, n_species 
    logical, intent(in) :: standard_grid
    integer(i4b) :: i, f

    I_ = nlons
    J_ = nlats
    R_ = radius

    dlat_ = pi / (J_ + 1)
    dlon_ = twopi / I_
    
    nonpole_cells_ = I_ * J_
    total_cells_ = nonpole_cells_ + 2;
    vcells_ = I_ * (J_ + 1)

    allocate( lats_(J_), sinlatsC_(J_), sinlatsB_(J_+1))
    allocate( latsB_(J_+1), lonsB_(I_) )    
    allocate( kfactor_by_lats_(J_))
    allocate( species_(n_species) )

    !!! coordinates
    lats_ = arange(dlat_, dlat_, J_)    
    do i=1, J_
        sinlatsC_(i) = sin( lats_(i) )
        sinlatsB_(i) = sin( lats_(i) - 0.5*dlat_ )
        f = 1
        if (standard_grid .eqv. .FALSE.) then 
            do while ( f*sinlatsC_(i) < 0.5 ) 
                f = next_factor(I_, f)    
            end do 
        endif
        kfactor_by_lats_(i) = f
    end do 
    latsB_(1:J_) = lats_ - 0.5*dlat_
    latsB_(J_+1) = lats_(J_)+0.5*dlat_

    sinlatsB_(J_+1) = sin( lats_(J_) + 0.5*dlat_)     
    lons_ = arange(dlon_/2, dlon_, I_)
    lonsB_ = lons_-0.5*dlon_ 

    !!! init each specie
    do i=1, n_species
        call init_specie(species_(i))
    enddo
end subroutine init_problem

subroutine init_specie(sp)
    implicit none 
    type (diff_specie), intent(inout) :: sp 
    allocate( sp%ff(J_, I_), sp%phi(J_, I_), sp%sigma(J_, I_) )
    allocate( sp%muLat(J_+1, I_), sp%muLon(J_, I_)  )
    sp%ff_np = 0.0
    sp%ff(:,:) = 0.0
    sp%ff_sp = 0.0

    sp%phi_np = 0.0
    sp%phi(:,:) = 0.0
    sp%phi_sp = 0.0

    sp%sigma_np = 0.0
    sp%sigma(:,:) = 0.0
    sp%sigma_sp = 0.0

    sp%muLat(:,:) = 0.0
    sp%muLon(:,:) = 0.0
end subroutine init_specie

end module diff_global_data 
