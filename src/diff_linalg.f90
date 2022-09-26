module diff_linalg
use diff_types 
    
contains
subroutine tridag(a,b,c,r,u)
    use diff_numeric_utils, only : throw_exception
    implicit none
    real(dp), dimension(:), intent(in) :: a,b,c,r
    real(dp), dimension(:), intent(out) :: u
    real(dp), dimension(size(b)) :: gam 
    integer(i4b) :: n,j
    real(dp) :: bet

    n = size(a)+1
    if ( (n .ne. size(b) ) .OR. (n .ne. size(c)+1) .or. (n .ne. size(r) ) .OR. ( n .ne. size(u) ) ) then 
        call throw_exception('tridag(): Error with matrix dimensions')
    end if 

    bet=b(1)
    if (bet == 0.0) call throw_exception('tridag: Error at code stage 1')

    u(1)=r(1)/bet
    do j=2,n
        gam(j)=c(j-1)/bet
        bet=b(j)-a(j-1)*gam(j)
        if (bet == 0.0) call throw_exception('tridag: Error at code stage 2')
        u(j)=(r(j)-a(j-1)*u(j-1))/bet
    end do

    do j=n-1,1,-1
        u(j)=u(j)-gam(j+1)*u(j+1)
    end do
END subroutine  tridag

subroutine tridag3(a,b,c,r1,r2,r3,u1,u2,u3)
    use diff_numeric_utils, only : throw_exception
    implicit none
    real(dp), dimension(:), intent(in) :: a,b,c,r1,r2, r3
    real(dp), dimension(:), intent(out) :: u1, u2, u3
    real(dp), dimension(size(b)) :: gam 
    integer(i4b) :: n,j
    real(dp) :: bet

    n = size(a)+1
    if ( (n .ne. size(b) ) .OR. (n .ne. size(c)+1) .or. (n .ne. size(r1) ) .OR. ( n .ne. size(u1) ) ) then 
        call throw_exception('tridag3(): Error with matrix dimensions')
    end if 

    bet=b(1)
    if (bet == 0.0) call throw_exception('tridag3(): Error at code stage 1')

    u1(1)=r1(1)/bet
    u2(1)=r2(1)/bet
    u3(1)=r3(1)/bet
    do j=2,n
        gam(j)=c(j-1)/bet
        bet=b(j)-a(j-1)*gam(j)
        if (bet == 0.0) call throw_exception('tridag3(): Error at code stage 2')
        u1(j)=(r1(j)-a(j-1)*u1(j-1))/bet
        u2(j)=(r2(j)-a(j-1)*u2(j-1))/bet
        u3(j)=(r3(j)-a(j-1)*u3(j-1))/bet
    end do

    do j=n-1,1,-1
        u1(j)=u1(j)-gam(j+1)*u1(j+1)
        u2(j)=u2(j)-gam(j+1)*u2(j+1)
        u3(j)=u3(j)-gam(j+1)*u3(j+1)
    end do
END subroutine  tridag3

subroutine tridag2(a,b,c,r1,r2,u1,u2)
    use diff_numeric_utils, only : throw_exception
    implicit none
    real(dp), dimension(:), intent(in) :: a,b,c,r1,r2
    real(dp), dimension(:), intent(out) :: u1, u2
    real(dp), dimension(size(b)) :: gam 
    integer(i4b) :: n,j
    real(dp) :: bet

    n = size(a)+1
    if ( (n .ne. size(b) ) .OR. (n .ne. size(c)+1) .or. (n .ne. size(r1) ) .OR. ( n .ne. size(u1) ) ) then 
        call throw_exception('tridag2(): Error with matrix dimensions')
    end if 

    bet=b(1)
    if (bet == 0.0) call throw_exception('tridag2(): Error at code stage 1')

    u1(1)=r1(1)/bet
    u2(1)=r2(1)/bet
    do j=2,n
        gam(j)=c(j-1)/bet
        bet=b(j)-a(j-1)*gam(j)
        if (bet == 0.0) call throw_exception('tridag2(): Error at code stage 2')
        u1(j)=(r1(j)-a(j-1)*u1(j-1))/bet
        u2(j)=(r2(j)-a(j-1)*u2(j-1))/bet
    end do

    do j=n-1,1,-1
        u1(j)=u1(j)-gam(j+1)*u1(j+1)
        u2(j)=u2(j)-gam(j+1)*u2(j+1)
    end do
END subroutine  tridag2

subroutine cyclic(a, b, c, alpha, beta, r, x)
    use diff_numeric_utils, only :  throw_exception
    implicit none
    real(dp), dimension(:), intent(in):: a, b, c, r
    real(dp), intent(in) :: alpha, beta
    real(dp), dimension(:), intent(out):: x

    integer(i4b) :: n
    real(dp) :: fact,gamma
    real(dp), dimension(size(x)) :: bb, u, z
    n = size(a)
    if ( (n .ne. size(b)) .or. (n .ne. size(c)) .or. (n .ne. size(r)) .or. (n .ne. size(x))) then
        call throw_exception('cyclic(): Error with matrix dimensions')
    endif

    if (n<=2) call throw_exception('cyclic(): System dimension must be greater than 2')

    gamma = -b(1) !Avoid subtraction error in forming bb(1).
    bb(1)=b(1) - gamma
    bb(n)=b(n) - alpha*beta/gamma
    bb(2:n-1)=b(2:n-1)
    !call tridag(a(2:n),bb,c(1:n-1),r,x)
    u(1)=gamma
    u(n)=alpha
    u(2:n-1)=0.0
    !call tridag(a(2:n),bb,c(1:n-1),u,z)
    call tridag2(a(2:n),bb,c(1:n-1),r,u,x,z)
    fact=(x(1)+beta*x(n)/gamma)/(1.0_dp+z(1)+beta*z(n)/gamma)
    x=x-fact*z
END subroutine cyclic

END module diff_linalg
