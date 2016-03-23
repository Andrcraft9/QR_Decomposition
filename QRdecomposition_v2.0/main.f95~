program main
    implicit none
    ! Variables
    integer, parameter :: n = 512
    integer*4 :: k, siz
    real*8 :: A(n, n), Q(n, n), R(n, n), P(n, n), vT(n, 1), vrT(n, 1), vq(n, 1)
    real*8 :: w1, w2
    real*8 :: error_of_qr
    real*8 :: newA(n, n)
    !real*8 :: QQT(n, n)
    real :: start, finish
    ! Функции
    real*8 :: dnrm2, ddot

    ! Program
    call init_random_seed

    call RANDOM_NUMBER(A)
    !call create_bad_matrix(n, 20, A)

    call create_I(n, Q)
    R = A
    P = 0

    call cpu_time(start)
    ! A = Q*R
	! P = I - (2/vT*v) * v*vT - householder matrix, v - householder
    ! Pn*Pn-1*...*P1*A = R
    ! P1*P2*...*Pn = Q
    do k = 1, n
        ! Compute Rn = Pn * Rn-1 and Qn = Qn-1 * PnT

        ! Set size our vectors
        siz = n - k + 1

        ! Compute Householder vector vT
        vT(1 : siz, 1) = R(k : k + siz - 1, k)

        w1 = dot_product(vT(1 : siz, 1), vT(1 : siz, 1))
        w1 = sqrt(w1)

        if (R(k, k) > 0) then
            vT(1, 1) = vT(1, 1) + w1
        else
            vT(1, 1) = vT(1, 1) - w1
        endif

		! Compute vrT = vT * R and compute w2 = (v, v)        
		vrT(1 : siz, 1) = matmul(vT(1 : siz, 1), R(k : k + siz - 1, k : k + siz - 1))
        w2 = dot_product(vT(1 : siz, 1), vT(1 : siz, 1))

		!Compute vq = Q * v
        vq(1 : n, 1) = matmul(Q(1 : n, k : k + siz - 1), vT(1 : siz, 1))
      
        ! Compute w2 = 2/(v, v)
        w2 = 2.0 / w2

		! Compute new R = R - w2(v * vrT)
        P(k : k + siz - 1, k : k + siz - 1) = w2 * matmul(vT(1 : siz, 1 : 1), transpose(vrT(1 : siz, 1 : 1)))
        R(k : k + siz - 1, k : k + siz - 1) = R(k : k + siz - 1, k : k + siz - 1) - P(k : k + siz - 1, k : k + siz - 1)

    	! Compute new Q = Q - w2(vq * vT)
        P(1 : n, k : k + siz - 1) = w2 * MATMUL(vq(1 : n, 1 : 1), transpose(vT(1 : siz, 1 : 1)))
        Q(1 : n, k : k + siz - 1) = Q(1 : n, k : k + siz - 1) - P(1 : n, k : k + siz - 1)
     
    enddo
    call cpu_time(finish)
    print *, finish - start

    !QQT = MATMUL(Q, TRANSPOSE(Q))
    newA = MATMUL(Q, R)

    print *, error_of_qr(A, newA, n)
end

subroutine create_I(n, A)
    implicit none
    ! Variables
    integer*4 n
    real*8 A(n, n)
    integer*4 i,j

    ! Function
    do j = 1, n
        do i = 1, n
            if (i == j) then
                A(i, i) = 1
            else
                A(i, j) = 0
            endif
        enddo
    enddo

    return
end

! Generate bad matrix Ak = H1 * diag(yk) * H2,
! where y1 / yn = 10**k and y1 > y2 ... yn-1 > yn,
! H1, H2 - Householder matrix
subroutine create_bad_matrix(n, k, A)
    implicit none
    ! Variables
    integer*4 :: n, k
    real*8 :: A(n, n)

    integer*4 :: i
    real*8 :: w1, w2
    real*8 :: v1(n, 1), v2(n, 1)
    real*8 :: v1T(1, n), v2T(1, n)

    real*8 :: H1(n, n), H2(n, n), Y(n, n)
    real*8 :: h, deg
    ! Function
    call init_random_seed

    ! Generate random Householder vector v1 and v2 and compute w = (v,v)
    call RANDOM_NUMBER(v1)
    call RANDOM_NUMBER(v2)


    w1 = dot_product(v1(1 : n, 1), v1(1 : n, 1))
    w2 = dot_product(v2(1 : n, 1), v2(1 : n, 1))
    !w1 = 0
    !w2 = 0
    !do i = 1, n
    !    w1 = w1 + v1(i, 1)**2
    !    w2 = w2 + v2(i, 1)**2
    !enddo

    w1 = 2 / w1
    w2 = 2 / w2

    ! Create Householder matrix H1 and H2
    call create_I(n, H1)
    call create_I(n, H2)
    v1T = TRANSPOSE(v1) ! !!!
    v2T = TRANSPOSE(v2) ! !!!
    H1 = H1 - w1*MATMUL(v1, v1T)
    H2 = H2 - w2*MATMUL(v2, v2T)

    ! Create daig(yk)
    call create_I(n, Y)
    h = REAL(k, 8) / (n - 1)
    deg = 0
    do i = 1, n
        Y(n - i + 1, n - i + 1) = 10**deg
        deg = deg + h
    enddo

    ! Create Ak = H1 * Y * H2
    A = MATMUL(H1, MATMUL(Y, H2))

    return
end

real*8 function error_of_qr(A, newA, n)
    implicit none
    ! Variables
    integer, intent(in) :: n
    real*8 A(n, n), newA(n, n)
    real*8 norm
    integer*4 i, j
    ! Function
    error_of_qr = 0
    norm = 0
    do j = 1, n
        do i = 1, n
            error_of_qr = error_of_qr + (A(i, j) - newA(i, j))**2
            norm = norm + A(i, j)**2
        enddo
    enddo
    norm = sqrt(norm)
    error_of_qr = sqrt(error_of_qr) / norm
    return
end function

subroutine init_random_seed()
    use iso_fortran_env, only: int64
    implicit none
    integer, allocatable :: seed(:)
    integer :: i, n, un, istat, dt(8), pid
    integer(int64) :: t

    call random_seed(size = n)
    allocate(seed(n))
    ! First try if the OS provides a random number generator
    open(newunit=un, file="/dev/urandom", access="stream", &
         form="unformatted", action="read", status="old", iostat=istat)
    if (istat == 0) then
       read(un) seed
       close(un)
    else
       ! Fallback to XOR:ing the current time and pid. The PID is
       ! useful in case one launches multiple instances of the same
       ! program in parallel.
       call system_clock(t)
       if (t == 0) then
          call date_and_time(values=dt)
          t = (dt(1) - 1970) * 365_int64 * 24 * 60 * 60 * 1000 &
               + dt(2) * 31_int64 * 24 * 60 * 60 * 1000 &
               + dt(3) * 24_int64 * 60 * 60 * 1000 &
               + dt(5) * 60 * 60 * 1000 &
               + dt(6) * 60 * 1000 + dt(7) * 1000 &
               + dt(8)
       end if
       pid = getpid()
       t = ieor(t, int(pid, kind(t)))
       do i = 1, n
          seed(i) = lcg(t)
       end do
    end if
    call random_seed(put=seed)
  contains
    ! This simple PRNG might not be good enough for real work, but is
    ! sufficient for seeding a better PRNG.
    function lcg(s)
      integer :: lcg
      integer(int64) :: s
      if (s == 0) then
         s = 104729
      else
         s = mod(s, 4294967296_int64)
      end if
      s = mod(s * 279470273_int64, 4294967291_int64)
      lcg = int(mod(s, int(huge(0), int64)), kind(0))
    end function lcg
  end subroutine init_random_seed
