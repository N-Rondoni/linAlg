module LinAl
  implicit none
  
  integer, save :: msize, nsize
  real(kind(0.d0)), dimension(:,:), allocatable, save :: mat
  

contains

  !********************************************************

  subroutine readMat(filename)

    implicit none
    character(len=*) :: filename

    integer :: i,j

    ! Reads a file containing the matrix A 
    ! Sample file:
    !
    ! 4 4 
    ! 2.0 1.0 1.0 0.0
    ! 4.0 3.0 3.0 1.0
    ! 8.0 7.0 9.0 5.0
    ! 6.0 7.0 9.0 8.0
    !
    ! Note that the first 2 numbers in the first line are the matrix dimensions, i.e., 4x4,
    ! then the next msize lines are the matrix entries. This matrix is found in Eq. 2.18 of the lecture note.
    ! Note that entries must be separated by a tab.


    open(10,file=filename)

    ! Read the matrix dimensions
    !read(10,*) i,j 

    ! Read matrix
    do i=1,msize
       read(10,*) ( mat(i,j), j=1,nsize )
    enddo

    close(10)
    
  end subroutine readMat

  !*********************************************************

  subroutine computeTrace(A, m, trace) !needs to be a funciton of trace?
    
    implicit none
    real (kind(0.d0)), intent(in)              :: A(:,:)
    integer, intent(in)                        :: m
    real (kind(0.d0)), intent(out)             :: trace
    ! Loop var
    integer :: i


    trace = 0.d0
    do i=1,m
       trace = trace + A(i,i)
    enddo


  end subroutine computeTrace

  !*********************************************************
  
  subroutine euclidNorm(v, m, norm)

    implicit none
    real (kind(0.d0)), intent(in)                :: v(:) 
    integer, intent(in)                          :: m 
    real (kind(0.d0)), intent(out)               :: norm
    ! Loop var
    integer                                      :: i
    
    
    norm = 0.d0
    do i=1,m
       norm = norm + v(i)**2
    enddo
    norm =SQRT(norm)
  
  end subroutine euclidNorm
  
  !*********************************************************

  subroutine matPrint(A, m, n)

    implicit none
    real(kind(0.d0)), intent(in)                :: A(:,:)
    integer, intent(in)                         :: m, n
    ! Loop vars
    integer :: i

    
    print *, m,n 
    do i = 1,m
       print *, A(i,:)
    enddo
   
  end subroutine 

!*********************************************************

 subroutine outerProd(u, v, A)
    ! need and outer product in Householder, writing it here
    ! takes in two vectors u, v
    ! returns their outer product (a matrix) A
    
    implicit none
    real(kind(0.d0)), intent(in)         :: u(:)
    real(kind(0.d0)), intent(in)         :: v(:)
    real(kind(0.d0)), intent(out)        :: A(:,:)    
    ! loop and temp vars
    integer                              :: m,n 
    integer                              :: i,j
 
    m = size(u)
    n = size(v)
    
    A = 0.d0
    do i = 1, m
       do j = 1, n
       A(i,j) = u(i)*v(j) 
       enddo
    enddo

 end subroutine   
!*********************************************************

 subroutine reconstruct(U, s, Vt, k, Appx)
 !takes in U (a matrix) s (a vector of singular values), and Vt (a matrix), k to denote the subset size
 ! returns U S VT, reconstructing the matrix A_sigma_k = Appx
    
    implicit none
    real(kind(0.d0)), intent(in)       :: U(:,:), Vt(:,:)
    real(kind(0.d0)), intent(in)       :: s(:)
    integer                            :: k
    real(kind(0.d0)), intent(out)      :: Appx(:,:)
    ! loop and temp vars
    integer                            :: m, n, rank
    real(kind(0.d0)), allocatable      :: tempMat(:,:)
    integer                            :: i

    m = size(U, 1)
    n = size(Vt, 1)
    rank = min(m, n) !big ol assumption
    ! in full decomp, U is m x m
    ! Vt is n x n 


    allocate(tempMat(m,n))
    tempMat = 0.d0
    Appx = 0.d0
    !print *, size(s(1)*U(:,1))
    !print *, size(Vt(1, :))
    do i = 1, k
       call outerProd(s(k)*U(:, k), Vt(k,:), tempMat)
       Appx = Appx + tempMat
    enddo

    deallocate(tempMat)

 end subroutine
!*********************************************************
 
 subroutine reconstruct2(U, s, Vt, Appx, r1, r2) !note arguemnts have changed 
    ! similiar to last reconstruct, but trying to get rid of call to outerProd
    ! include start and end (r1, r2) to call for higher values as needed in assignment
    implicit none
    real(kind(0.d0)), intent(in)       :: U(:,:), Vt(:,:)
    real(kind(0.d0)), intent(in)       :: s(:)
    integer                            :: r1, r2
    real(kind(0.d0)), intent(out)      :: Appx(:,:)
    ! temp and loop vars
    integer                            :: m, n, i, j, k
    m = size(U, 1)
    n = size(Vt, 1)

    ! There has to be a way to do this with slicing but I'm too simple to figure it out
    do j = 1, n
       do i = 1, m
          do k = r1, r2
             Appx(i, j) = Appx(i, j) + s(k)*U(i, k)*Vt(k, j)
          enddo
       enddo
    enddo
    
 end subroutine 
!*********************************************************

 subroutine write_output(A, fileName)
  
   implicit none
   real(kind(0.d0)), intent(in) :: A(:,:)   
   character(len=*) :: fileName !must be passed in as a string, e.g., 'fileName'

   open(20, file=fileName, status="replace")  
     write(20,*) Transpose(A) 
   close(20) 

  end subroutine 

!********************************************************

 subroutine inputMatMaker(D, A, ldA, b)
 ! returns a matrix full of ones, aside from main diagonal filled with D
 ! also returns a bector b filled with integers 1, .., D
 ! note choice of ldA is choice of dimension for sq matrix (allocate matrix in driver to be size A(ldA,ldA))
   implicit none
   integer, intent(in)              :: D
   real(kind(0.d0)), intent(out)    :: A(:,:)
   real(kind(0.d0)), intent(out)    :: b(:)
   integer, intent(in)              :: ldA
   integer                          :: i

   A = 1.d0
   ! fill A and b with appropriate things
   do i = 1, ldA
      A(i,i) = D
      b(i) = i
   enddo

 end subroutine

!***********************************************************
 
 subroutine Jacobi(A, x, b, accuracy, counter)
 ! get slighted gauss (lowercase, extra shade) you have enough things in your name
 ! Arguments:
 ! Solves system Ax = b  to desired accuracy (returns x given A, b)
 ! returns counter, number of iterations it took to achieve desired accuracy.
   implicit none
   real(kind(0.d0)), intent(in)  :: A(:,:)
   real(kind(0.d0)), intent(out) :: x(:)
   real(kind(0.d0)), intent(in)  :: b(:)
   real(kind(0.d0)), intent(in)  :: accuracy
   integer, intent(out)          :: counter
   ! loop and temp vars
   real(kind(0.d0)), allocatable :: Dmat(:), R(:,:)
   integer                       :: ldA, i
   real(kind(0.d0)), allocatable :: res(:)
   ! create D and R from A
   ldA = size(A, 1) !A square
   allocate(Dmat(ldA))
   allocate(R(ldA, ldA))
   allocate(res(ldA))
   Dmat = 0.d0
   R = A
   do i = 1, ldA
     Dmat(i) = A(i,i)
     R(i,i) = 0.d0 
   enddo
   !call matPrint(Dmat, ldA, ldA)
   !call matPrint(R, ldA, ldA)
   ! set x equal to guess
   x(1) = 1.d0  
   ! set res to something big to get loop cookin
   res = 20.d0
  
   counter = 0

   do while (norm2(res) > accuracy)
      ! update x
      x = (b - matmul(R,x))/Dmat !A(i,i) could be dmat, since they share the same diag
      !print *, "made it to here
      ! update res to get out
      res = b - matmul(A, x)
      counter = counter + 1
      !print *, counter, norm2(res)
   enddo


   deallocate(res)
   deallocate(Dmat)
   deallocate(R)
 end subroutine

!***********************************************************
 
 subroutine Seidel(A, x, b, accuracy, counter)
 ! get slighted gauss (lowercase, extra shade) you have enough things in your name
 ! Solves system Ax = b  to desired accuracy
   implicit none
   real(kind(0.d0)), intent(in)  :: A(:,:)
   real(kind(0.d0)), intent(out) :: x(:)
   real(kind(0.d0)), intent(in)  :: b(:)
   real(kind(0.d0)), intent(in)  :: accuracy
   integer, intent(out)          :: counter
   ! loop and temp vars
   real(kind(0.d0)), allocatable :: Dmat(:,:), R(:,:)
   integer                       :: ldA, i, j
   real(kind(0.d0)), allocatable :: res(:)
   real(kind(0.d0))              :: tempSum1, tempSum2
 
   ! create D and R from A
   ldA = size(A, 1) !A square
   allocate(Dmat(ldA, ldA))
   allocate(R(ldA, ldA))
   allocate(res(ldA))
   Dmat = 0.d0
   R = A
   do i = 1, ldA
     Dmat(i,i) = A(i,i) !I never use this but I could have
     R(i,i) = 0.d0 
   enddo
   !call matPrint(Dmat, ldA, ldA)
   !call matPrint(R, ldA, ldA)
   ! set x equal to guess
   x(1) = 1.d0  
   res = 20.d0
   !tempSum = 0.d0
   counter = 0

   do while (norm2(res) > accuracy)
      ! update x
      do i = 1, ldA
         tempSum1 = 0.d0   
         tempSum2 = 0.d0
         do j=i+1, ldA !this is rh sum in eqn 6.15
            tempSum1 = tempSum1 + R(i,j)*x(j)  
         enddo
         do j=1, i-1
            tempSum2 = tempSum2 + R(i,j)*x(j)
         enddo
         x(i) = (b(i) - tempSum2 - tempSum1)/A(i,i) !A(i,i) could be dmat, since they share the same diag
         !print *, "made it to here"
      enddo 
      ! update res to get out
      res = b - matmul(A, x)
      counter = counter + 1
      !print *, counter, norm2(res)
   enddo


   deallocate(res)
   deallocate(Dmat)
   deallocate(R)
 end subroutine

!***********************************************************

 subroutine conjGrad(A, x, b, accuracy, counter)
   implicit none
   real(kind(0.d0)), intent(in)  :: A(:,:)
   real(kind(0.d0)), intent(out) :: x(:)
   real(kind(0.d0)), intent(in)  :: b(:)
   real(kind(0.d0)), intent(in)  :: accuracy
   integer, intent(out)          :: counter 
   real(kind(0.d0)), allocatable :: y(:), res(:), p(:)
   integer                       :: ldA, i, j
   real(kind(0.d0))              :: error, errorN, alpha, beta

   ldA = size(A, 1)
   allocate(y(ldA))
   allocate(res(ldA))
   allocate(p(ldA))

   ! check if matrix is symmetric
   ! assume matrix is symmetric unless proven guilty (presumption of innocence applies here)
   do i = 1, Lda
      do j = 1, ldA
         if (A(i,j) .NE. A(j,i)) then
            print *, "CG needs a symmetric matrix, ya donkey"
            stop
         endif
      enddo
   enddo
   ! set x to guess, for this method x = 0 vector
   x = 0.d0
   res = b - matmul(A, x)
   p = res
   error = norm2(res)**2
   counter = 0

   do while(sqrt(error) > accuracy)
      y = matmul(A, p)
      alpha = dot_product(p, res)/(dot_product(p, y))
      x = x + alpha*p
      res = res - alpha*y
      errorN = norm2(res)**2
      beta = -dot_product(res, y)/dot_product(p, y) !was using errorN in numerator, converged but real slow?
      p = res + beta*p
      counter = counter + 1 
      error = errorN
   enddo

   deallocate(y)
   deallocate(res)
   deallocate(p) 
 end subroutine
!***********************************************************

end module LinAl
