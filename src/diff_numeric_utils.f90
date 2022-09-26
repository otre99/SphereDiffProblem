module diff_numeric_utils
    use diff_types
    implicit none

    interface swap
        module procedure swap_i, swap_r, swap_rv, swap_iv
    end interface
    interface imaxloc
        module procedure imaxloc_r, imaxloc_i
    end interface
    interface arange
        module procedure arange_r, arange_i
    end interface
    
contains
    !!!!!!!!!!!!!!!!! arange !!!!!!!!!!!!!!!!!1
    function arange_r(first, increment,n)
    real(dp), intent(in) :: first, increment
    integer(i4b), intent(in) :: n
    real(dp), dimension(n) :: arange_r
    integer(i4b) :: k, k2
    real(dp) :: temp
    
    if (n > 0) arange_r(1) = first
    if (n <= 16) then
        do k=2,n
            arange_r(k)=arange_r(k-1)+increment
        end do
    else
        do k=2,8
            arange_r(k)=arange_r(k-1)+increment
        end do
        temp=increment*8
        k=8
        do
            if (k >= n) exit
            k2 = k + k
            arange_r(k+1:min(k2,n) ) = temp + arange_r(1:min(k,n-k))
            temp = temp + temp
            k = k2
        end do
    end if
    end function arange_r

    function arange_i(first,increment,n)
    integer(i4b), intent(in) :: first,increment,n
    integer(i4b), dimension(n) :: arange_i
    integer(i4b) :: k,k2,temp
    if (n > 0) arange_i(1)=first
    if (n <= 16) then
        do k=2,n
            arange_i(k)=arange_i(k-1)+increment
        end do
    else
        do k=2,8
            arange_i(k)=arange_i(k-1)+increment
        end do
        temp=increment*8
        k=8
        do
            if (k >= n) exit
            k2=k+k
            arange_i(k+1:min(k2,n))=temp+arange_i(1:min(k,n-k))
            temp=temp+temp
            k=k2
        end do
    end if
    end function arange_i
    !!!!!!!!!!!!!!!!! end_arange !!!!!!!!!!!!!!!!!1    
    
    
    
    !!! linspace() 
    function linspace(first, eend, n)
    real(dp), intent(in) :: first, eend
    integer(i4b), intent(in) :: n
    real(dp), dimension(n) :: linspace
    linspace = arange(first, (eend-first)/(n-1), n)
    end function linspace
    
    subroutine throw_exception(msg)
        character(len=*), intent(in) :: msg
        write (*,*) msg
        stop 'program terminated by runtime error'
    end subroutine throw_exception
    
    
    !!!!!!!!!!!!!!!!! swap !!!!!!!!!!!!!!!!!1
    subroutine swap_i(a,b)
        integer(i4b), intent(inout) :: a,b
        integer(i4b) :: dum
       dum = a
        a = b
        b = dum
    end subroutine swap_i

    subroutine swap_r(a,b)
        real(dp), intent(inout) :: a,b
        real(dp) :: dum
        dum=a
        a=b
        b=dum
    end subroutine swap_r
    
    subroutine swap_rv(a,b)
        real(dp), dimension(:), intent(inout) :: a,b
        real(dp), dimension(size(a)) :: dum
        dum=a
        a=b
        b=dum
    end subroutine swap_rv

    subroutine swap_iv(a,b)
        integer(i4b), dimension(:), intent(inout) :: a,b
        integer(i4b), dimension(size(a)) :: dum
        dum=a
        a=b
        b=dum
    end subroutine swap_iv
    !!!!!!!!!!!!!!!!! end_swap !!!!!!!!!!!!!!!!!

    function outerand(a,b)
        logical, dimension(:), intent(in) :: a,b
        logical, dimension(size(a),size(b)) :: outerand
        outerand = spread(a,dim=2,ncopies=size(b)) .and. &
            spread(b,dim=1,ncopies=size(a))
    end function outerand
    
    function outerprod(a,b)
        real(dp), dimension(:), intent(in) :: a,b
        real(dp), dimension(size(a),size(b)) :: outerprod
        outerprod = spread(a,dim=2,ncopies=size(b)) * &
            spread(b,dim=1,ncopies=size(a))
    end function outerprod

    !!!!!!!!!!!!!!!!! imaxloc !!!!!!!!!!!!!!!!!1
    function imaxloc_r(arr)
        real(dp), dimension(:), intent(in) :: arr
        integer(i4b) :: imaxloc_r
        integer(i4b), dimension(1) :: imax
        imax=maxloc(arr(:))
        imaxloc_r=imax(1)
    end function imaxloc_r

    function imaxloc_i(iarr)
        integer(i4b), dimension(:), intent(in) :: iarr
        integer(i4b), dimension(1) :: imax
        integer(i4b) :: imaxloc_i
        imax=maxloc(iarr(:))
        imaxloc_i=imax(1)
    end function imaxloc_i    
    !!!!!!!!!!!!!!!!! end_imaxloc !!!!!!!!!!!!!!!!!1
        
end module diff_numeric_utils
