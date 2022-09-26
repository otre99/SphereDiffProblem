module diff_io
use diff_types

contains
subroutine save_matrix_to_file(fname, ttime, data)
    implicit none    
    character(len=*), intent(in) :: fname
    real(dp), intent(in) :: ttime 
    real(dp), dimension(:,:), intent(in) :: data     

    real(dp) :: rows, cols, n 
    n = size(data)
    rows = size(data,1)
    cols = size(data,2)
    open(unit=1, file=fname, action="write", status='replace', form='unformatted', access='stream')
    write(1) ttime, rows, cols, n, data     
    close(1)
end subroutine save_matrix_to_file

subroutine save_specie_to_file(fname, ttime, north, data, south)
    implicit none    
    character(len=*), intent(in) :: fname
    real(dp), intent(in) :: ttime, north, south 
    real(dp), dimension(:,:), intent(in) :: data     
    
    real(dp) :: rows, cols, n 
    n = size(data) + 2
    rows = size(data,1)
    cols = size(data,2)

    write(*,*) "Saving file: ", fname
    open(unit=1, file=fname, action="write", status='replace', form='unformatted', access='stream')
    write(1) ttime, rows, cols, n, north, data, south     
    close(1)
end subroutine save_specie_to_file

subroutine get_header(funit, ttime, rows, cols, ntotal)
    implicit none    
    integer(i4b), intent(in):: funit 
    real(dp), intent(inout) :: ttime
    integer(i4b), intent(inout):: rows, cols, ntotal 
    real(dp) :: d_rows, d_cols, d_n 
    read(funit) ttime, d_rows, d_cols, d_n
    rows = int(d_rows)
    cols = int(d_cols) 
    ntotal = int(d_n)     
end subroutine get_header

subroutine read_header_from_file(fname, ttime, rows, cols, ntotal)
    implicit none    
    character(len=*), intent(in) :: fname
    real(dp), intent(inout) :: ttime
    integer(i4b), intent(inout):: rows, cols, ntotal 

    open(unit=1, file=fname, action="read", status='old', form='unformatted', access='stream')
    call get_header(1, ttime, rows, cols, ntotal)
    close(1)
end subroutine read_header_from_file

subroutine read_matrix_from_file(fname, ttime, data)
    use diff_numeric_utils, only : throw_exception
    implicit none    
    character(len=*), intent(in) :: fname
    real(dp), intent(inout) :: ttime
    real(dp), dimension(:,:), intent(inout) :: data     

    integer(i4b) :: rr, cc, nn

    open(unit=1, file=fname, action="read", status='old', form='unformatted', access='stream')
    call get_header(1, ttime, rr, cc, nn)

    if (size(data,1) .ne. rr .OR. size(data,2) .ne. cc) then
        call throw_exception('read_matrix_from_file(): unexpected input dimensions in file : '//fname)
    endif 

    if (size(data) .ne. nn) then 
        call throw_exception('read_matrix_from_file(): total != rows x cols')
    endif 
    read(1) data 

    close(1)
end subroutine read_matrix_from_file

subroutine read_specie_from_file(fname, ttime, north, data, south)
    use diff_numeric_utils, only : throw_exception
    implicit none    
    character(len=*), intent(in) :: fname
    real(dp), intent(inout) :: ttime
    real(dp), dimension(:,:), intent(inout) :: data     
    real(dp), intent(inout) :: north, south     

    integer(i4b) :: rr, cc, nn

    open(unit=1, file=fname, action="read", status='old', form='unformatted', access='stream')
    call get_header(1, ttime, rr, cc, nn)

    if (size(data,1) .ne. rr .OR. size(data,2) .ne. cc) then
        call throw_exception('read_specie_from_file(): Wrong input dimensions in file: '//fname)
    endif 

    if (size(data) .ne. nn-2) then 
        call throw_exception('read_specie_from_file(): total != rows x cols + 2')
    endif 
    read(1) north, data, south 
    close(1)
end subroutine read_specie_from_file

function gen_outname(k, iter)
    implicit none 
    integer(i4b), intent(in) :: k, iter    
    character(len=25) :: gen_outname
    write(gen_outname, '(a,I0.2,a,I0.16,a)') 'sp',k,'_',iter,'.out'       
end function gen_outname

subroutine save_sol(outf, ttime, iter, sp)
    use diff_numeric_utils, only : throw_exception
    implicit none    
    character(len=*), intent(in) :: outf
    integer(i4b), intent(in) :: iter    
    type (diff_specie), dimension(:) :: sp 
    real(dp), intent(in) :: ttime

    integer(i4b) :: k 

    do k=1, size(sp)
        call save_specie_to_file(trim(outf)//gen_outname(k-1,iter),&
             ttime, sp(k)%phi_np, sp(k)%phi, sp(k)%phi_sp) 
    enddo
end subroutine save_sol



subroutine print_info_base(method, tau, standard_grid)
    use diff_global_data, only : dlat_, dlon_, I_, J_, R_, species_
    implicit none 
    real(dp), intent(in) :: tau 
    logical, intent(in) :: standard_grid
    integer(i4b), intent(in) :: method

    integer(i4b) :: nspecies
    nspecies = size(species_)
    write(*,"(a,I14)") "Species       : ", nspecies
    write(*,"(a,F14.8)") "Radius        : ", R_
    write(*,"(a,F14.8,a,F14.8)") "Time step (dt): ", 2*tau, "      Tau : ", tau 
    write(*,"(a,F14.8,a,F14.8,a32,I6)") "Lat. res      : ", 180.0*dlat_/pi, " (    len = ", R_*dlat_, " ) J = ", J_ 
    write(*,"(a,F14.8,a,F14.8,a,F14.8,a,I6)") "Lon. res      : ", 180.0*dlon_/pi, " ( minlen = ", dlon_ * R_ * sin(dlat_),&
    " | maxlen = ", dlon_ * R_ , ") I = ", I_
    write(*,"(a,l14)") "Standard Grid : ", standard_grid

    if (method .eq. 1)  then 
        write(*,"(a,a)") "Int. method   : CN_DIM_SPLIT (LINEAR)" 
    else if (method .eq. 2) then 
        write(*,"(a,a)") "Int. method   : CN_DIM_SPLIT (NONLINEAR)" 
    else if (method .eq. 3) then 
        write(*,*) "Int. method   : RK2_TVD_DIM_SPLIT" 
    else 
        write(*,*) "Int. method   : UNKNOW" 
    endif 
end subroutine print_info_base

subroutine print_specie_info(sp, tau)
    use diff_global_data, only : dlat_, dlon_, I_, J_, R_,&
    sinlatsC_,kfactor_by_lats_
    implicit none 
    real(dp), intent(in) :: tau 
    type (diff_specie), intent(in) :: sp 

    integer(i4b) :: k, kf
    real(dp) :: l2norm, mass, area_pole, area  
    real(dp) :: lenlat, lenlat2, cfldiff_lat,  cfldiff_lon
    real(dp) :: eps, lenlon, lenlon2  
    real(dp) :: vmin, vmax 
                   
    cfldiff_lat = -1
    cfldiff_lon = -1
    eps = 1e-31
    lenlat = R_ * dlat_
    lenlat2 = lenlat * lenlat     

    area_pole = 0.25 * pi * R_ * R_ * dlat_ * dlat_
    l2norm = area_pole * sp%phi_np**2
    mass   = area_pole * sp%phi_np
    do k=1,J_
        area = R_**2*dlat_*dlon_*sinlatsC_(k)
        l2norm = l2norm + area*dot_product(sp%phi(k,:), sp%phi(k,:))
        mass = mass + area*sum(sp%phi(k,:))

        kf = kfactor_by_lats_(k)
        lenlon =  R_ * sinlatsC_(k) * dlon_ * kf
        lenlon2 = lenlon**2      
        cfldiff_lon = max( cfldiff_lon, maxval( abs(sp%mulon(k, ::kf)) )/lenlon2)
        cfldiff_lat = max( cfldiff_lat, maxval( abs(sp%mulat(k, :)) )/lenlat2)
    enddo
    l2norm = l2norm + area_pole * sp%phi_sp**2;
    mass = mass + area_pole * sp%phi_sp

    k = J_+1
    vmin = min(sp%phi_np, sp%phi_sp)
    vmin = min(minval(sp%phi), vmin)

    vmax = max(sp%phi_np, sp%phi_sp)
    vmax = max(maxval(sp%phi), vmax)

    write(*,"(a,F24.8)") "    l2-norm      : ", l2norm 
    write(*,"(a,F24.8)") "    mass         : ", mass
    write(*,"(a,F16.10,a,F16.10)") "    D (LON) :", cfldiff_lon * tau
    write(*,"(a,F16.10,a,F16.10)") "    D (LAT) :", cfldiff_lat * tau
    write(*,"(a,F24.12)") "    Min. val.    : ", vmin 
    write(*,"(a,F24.12)") "    Max. val.    : ", vmax 
    write(*,"(a,F24.12)") "    N. Pole.     : ", sp%phi_np 
    write(*,"(a,F24.12)") "    S. Pole.     : ", sp%phi_sp

end subroutine print_specie_info

subroutine print_info(sp, tau)
    implicit none 
    real(dp), intent(in) :: tau 
    type (diff_specie), dimension(:), intent(in) :: sp 
    integer(i4b) :: i

    do i=1,size(sp)
        write(*,*) "Specie: ", i 
        call print_specie_info(sp(i), tau)
    enddo
end subroutine print_info

end module diff_io
