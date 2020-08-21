!
! Solve linear equation for x (Ax = b => x = bA^-1) using LU decomposition and back substitution.
! Input:
!   X  = An m by n array.
!   TX = Precalculated transpose array of X, size n by m
!   Y  = An m-element vector containing the right-hand side of the linear system Ax = b.
! Output:
!   B  = An n-element vector.

subroutine least_squares (x, y, tx, b)
  use type
  implicit none
 
  interface
    subroutine ludcmp (a, indx, d)
      use type
      real (dp), dimension (:, :), intent (inout) :: a
      integer (i4b), dimension (:), intent (out) :: indx
      real (dp), intent (out) :: d
    end subroutine ludcmp
 
    subroutine lubksb (a, indx, b)
      use type
      real (dp), dimension (:, :), intent (in) :: a
      integer (i4b), dimension (:), intent (in) :: indx
      real (dp), dimension (:), intent (inout) :: b
    end subroutine lubksb
  end interface
 
  real (dp), intent (in) :: x (:, :)
  real (dp), intent (in) :: y (:)
  real (dp), intent (in) :: tx (:, :)
  real (dp), allocatable, intent (out) :: b (:)
 
  real (dp), allocatable :: a (:, :)
  integer (i4b), allocatable :: indx (:)
  integer (i4b) :: nvars, ntimes
  real (dp) :: d
 
  nvars = size (x, 2) - 1
  ntimes = size (y)
 
!  print *,'least'
  allocate (b(nvars+1))
  allocate (a(nvars+1, nvars+1))
  allocate (indx(nvars+1))
!  print *,'allocated'
  !print *, "Y = ", Y
  !print *, "X = ", X
  b = matmul (tx, y)
  a = matmul (tx, x)
  !print *, "A = ", A
  !print *, "B = ", B
!print *,'lu start'
  call ludcmp (a, indx, d)
!print *,'ludcmp'
  if (any(abs(a) < 9.99999968E-15)) then
    b (:) = 0.0d0
!AJN
!     print *, "Warning, LUdcmp produced a zero."
    return
  end if
!print *,'lubksb'
  !print *, "LU A = ", A
  call lubksb (a, indx, b)
 
  !print *, matmul(matmul(TX, X), B)
!print *,'deallocate'
  deallocate (a)
  deallocate (indx)
!print *,'done'
end subroutine least_squares
 
 
! -----------------------------
subroutine logistic_regression (x, y, tx, yp, b)
  use type
  implicit none
 
  interface
    subroutine least_squares (x, y, tx, b)
      use type
      real (dp), intent (in) :: x (:, :)
      real (dp), intent (in) :: y (:)
      real (dp), intent (in) :: tx (:, :)
      real (dp), allocatable, intent (out) :: b (:)
    end subroutine least_squares
  end interface
 
  real (dp), intent (in) :: x (:, :)
  real (dp), intent (in) :: y (:)
  real (dp), intent (in) :: tx (:, :)
  integer (i4b), intent (in) :: yp (:)
  real (dp), allocatable, intent (out) :: b (:)
 
  real (dp), allocatable :: ypd (:), p (:), yn (:), bn (:), v (:, :), xv (:, :)
  !  real(DP), allocatable :: P(:), YN(:), BN(:), V(:,:), XV(:,:)
  integer (i4b) :: nvars, ntimes, i, t, f, it
  !real (dp) :: d
 
  nvars = size (x, 2) - 1
  ntimes = size (y)
 
  allocate (b(nvars+1))
  allocate (ypd(ntimes))
  allocate (yn(ntimes))
  allocate (p(ntimes))
  allocate (v(ntimes, ntimes))
  allocate (xv(ntimes, nvars+1))
 
  do t = 1, ntimes, 1
    if (yp(t) .gt. 0.0) then
      ypd (t) = 1.0d0
    else
      ypd (t) = 0.0d0
    end if
  end do
 
  b = 0.0d0
  i = 0
  it = 0
  f = 0 ! AWW added initialization
  do while (f /=  1)
    !print*, 'matmul(x,b):'
    !print*, matmul(x, b)
    p = 1.0d0 / (1.0d0+exp(-matmul(x, b)))
    if (any(p > 0.97)) then
      ! print *, "WARNING: logistic regression diverging"
      f = 1
    else
 
      yn = ypd - p
      v = 0.0d0
      do t = 1, ntimes, 1
        v (t, t) = p (t) * (1.0d0-p(t))
      end do
      xv = matmul (v, x)

      call least_squares (xv, yn, tx, bn)
 
      f = 1
      do i = 1, nvars + 1, 1
        if (bn(i) .gt. 1.0E-04 .or. bn(i) .lt.-1.0E-04) then
          f = 0
        end if
      end do
      if (it > 8) then
        ! print *, "WARNING: logistic regression failed to converge"
        f = 1
      end if
 
      b = b + bn
      deallocate (bn)
 
    end if
    it = it + 1
  end do
  ! print *, "Iterations = ", it
  ! print *, "Final B = ", B
 
end subroutine logistic_regression
 
 
subroutine logistic_regressionrf (x, y, tx, b)
  use type
  implicit none
 
  interface
    subroutine least_squares (x, y, tx, b)
      use type
      real (dp), intent (in) :: x (:, :)
      real (dp), intent (in) :: y (:)
      real (dp), intent (in) :: tx (:, :)
      real (dp), allocatable, intent (out) :: b (:)
    end subroutine least_squares
  end interface
 
  real (dp), intent (in) :: x (:, :)
  real (dp), intent (in) :: y (:)
  real (dp), intent (in) :: tx (:, :)
  real (dp), allocatable, intent (out) :: b (:)
 
  real (dp), allocatable :: ypd (:), p (:), yn (:), bn (:), v (:, :), xv (:, :)
  !  real(DP), allocatable :: P(:), YN(:), BN(:), V(:,:), XV(:,:)
  integer (i4b) :: nvars, ntimes, i, t, f, it
  !real (dp) :: d
 
  nvars = size (x, 2) - 1
  ntimes = size (y)
 
  !print *, 'logistic regression',ntimes,nvars,Yp
 
  allocate (b(nvars+1))
  allocate (ypd(ntimes))
  allocate (yn(ntimes))
  allocate (p(ntimes))
  allocate (v(ntimes, ntimes))
  allocate (xv(ntimes, nvars+1))
 
  do t = 1, ntimes, 1
    if (y(t) .ne. 0.0) then
      ypd (t) = 1.0d0
    else
      ypd (t) = 0.0d0
    end if
  end do
 
  b = 0.0d0
  i = 0
  it = 0
  f = 1 ! AWW added initialization
  !print *, "B = ", B
  do while (f /=  1)
    p = 1.0d0 / (1.0d0+exp(-matmul(x, b)))
    if (any(p > 0.97)) then
        !print *, "WARNING: logistic regression diverging"
      f = 1
    else
 
      yn = ypd - p
      v = 0.0d0
      do t = 1, ntimes, 1
        v (t, t) = p (t) * (1.0d0-p(t))
      end do
      xv = matmul (v, x)
      call least_squares (xv, yn, tx, bn)
 
      f = 1
      do i = 1, nvars + 1, 1
        if (bn(i) .gt. 1.0E-04 .or. bn(i) .lt.-1.0E-04) then
          f = 0
        end if
      end do
      if (it > 8) then
        ! print *, "WARNING: logistic regression failed to converge"
        f = 1
      end if
 
      b = b + bn
      deallocate (bn)
 
    end if
    it = it + 1
  end do
  ! print *, "Iterations = ", it
  ! print *, "Final B = ", B
 
end subroutine logistic_regressionrf


! Hongli add k-fold cross-validation subroutine
! -----------------------------
subroutine cross_validation (x, y, w, err)
  use type
  implicit none
 
  interface
    subroutine least_squares (x, y, tx, b)
      use type
      real (dp), intent (in) :: x (:, :)
      real (dp), intent (in) :: y (:)
      real (dp), intent (in) :: tx (:, :)
      real (dp), allocatable, intent (out) :: b (:)
    end subroutine least_squares
  end interface
 
  real (dp), intent (in) :: x (:, :)                                        ! spatoial attributes of sta_limit stations
  real (dp), intent (in) :: y (:)                                           ! observations of sta_limit stations
  real (dp), intent (in) :: w (:, :)                                        ! weights of sta_limit stations
  real (dp), intent (out) :: err                                            ! mean error of kfold tests 
 
  integer (i4b), parameter :: kfold = 5                                     ! k-fold
  real (dp), allocatable :: x_train (:, :), y_train (:), w_train (:, :)     ! divide x,y,w into training and test groups
  real (dp), allocatable :: x_test (:, :), y_test (:), w_test (:, :)        
  real (dp), allocatable :: tx_train (:, :), twx_train (:, :), b_train (:)  ! intermediate variables of regression
  real (dp) :: errsum, wgtsum                                               ! error sum and weight sum for one of kfold tests
  integer (i4b) :: nstns, nvars, ntests, ntrains                            ! number of stations, attributes, test and train samples
  integer (i4b) :: i, j, k, ctrain, ctest                                   ! variables for loop and count
    
  nstns = size (x, 1)                                                       ! total number of stations 
  nvars = size (x, 2)                                                       ! total number of spatial attributes 
  ntests = floor(nstns/real (kfold, kind(dp)))                              ! number of stations in the test group
  ntrains = nstns-ntests                                                    ! number of stations in the training group  
  err = 0.0                    
  
  do k = 1, kfold, 1
      allocate (x_train(ntrains, nvars))
      allocate (x_test(ntests, nvars))
      allocate (y_train(ntrains))
      allocate (y_test(ntests))
      allocate (w_train(ntrains, ntrains))
      allocate (w_test(ntests, ntests))
      allocate (tx_train(nvars, ntrains))
      allocate (twx_train(nvars, ntrains))
      
      ctrain = 0
      ctest = 0
      x_train = 0.0
      x_test = 0.0
      y_train = 0.0
      y_test = 0.0
      w_train = 0.0
      w_test = 0.0
      errsum = 0.0
      wgtsum = 0.0
 
      ! divide arrays into training and test groups
      do i = 1, nstns, 1
          if ( i >= ((k-1)*ntests+1) .AND. i <= (k*ntests) ) then
              ctest = ctest + 1
              x_test(ctest, :) = x(i,:)
              y_test(ctest) = y(i)
              w_test(ctest, ctest) = w(i, i)              
          else
              ctrain = ctrain + 1
              x_train(ctrain, :) = x(i,:)
              y_train(ctrain) = y(i)
              w_train(ctrain, ctrain) = w(i, i)               
          end if
      end do ! end i loop 
                
      ! estimate regression coefficients with training data
      tx_train = transpose(x_train)
      twx_train = matmul(tx_train, w_train)
      call least_squares(x_train, y_train, twx_train, b_train)
              
      ! estimate regression error with test data
      do j = 1, ntests, 1
          wgtsum = wgtsum + w_test (j, j)
          errsum = errsum + (w_test (j, j)*(real (dot_product(x_test(j, :), b_train), kind(sp))-y_test(j))**2)    
      end do ! end j loop
      
      err = err + real ((errsum/wgtsum)**(1.0/2.0), kind(sp))
       
      deallocate (x_train)  
      deallocate (x_test)  
      deallocate (y_train)  
      deallocate (y_test)  
      deallocate (w_train)  
      deallocate (w_test)  
      deallocate(tx_train)  
      deallocate(twx_train)
      deallocate (b_train)        
  
  end do ! end k loop
  
  err = err/kfold  
 
end subroutine cross_validation

