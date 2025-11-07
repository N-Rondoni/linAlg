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
    read(10,*) i,j 

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
    integer :: i, j

    
    print *, m,n 
    do i = 1,m
       print *, A(i,:)
    enddo
   
  end subroutine 

  !*********************************************************

  subroutine gaussElim(A, B, Am, An, Bm, Bn, logicFlag)
      ! A,B are matrices. A sqaure, B non square
      ! Am, An the dimensions m x n of A
      ! Logic flag to determine if problem singular or not
      ! Should return an uppter triangular mx A
      ! B contains modified set of rhs vectors  

     implicit none 
     real (kind(0.d0)), intent(in out)    :: A(:,:), B(:,:)
     integer, intent(in)                  :: Am, An, Bm, Bn
     logical, intent(out)                 :: logicFlag
     
     ! Loop and temp vars
     integer                              :: i,j
     integer                              :: pos
     real (kind(0.d0))                    :: pivot
     real (kind(0.d0)), allocatable       :: temp(:)
     

     allocate(temp(Am))
     logicFlag = .FALSE.
     temp = 0.0

     !looping over rows
     do j=1, Am-1
       ! the below loop finds pivot, position of pivot (max in column, index of max)
       pivot = 0.0
       ! loop over columns in this row
       do i=j, An
         if (ABS(A(i,j)) > pivot) then 
            pivot =  ABS(A(i,j)) 
            pos = i
         endif
       enddo
       ! next we interchange rows if pivot position (pos) != row position (j)  
       ! this ensures max value is found on the diagonal
       if (pos .NE. j) then
         temp(:) = A(pos,:) 
         A(pos,:) = A(j,:)
         A(j,:) = temp(:)
         ! don't forget to make the same changes to B
         temp(:) = B(pos,:)
         B(pos,:) = B(j,:)
         B(j,:) = temp(:)
       endif 
       ! next we check if the matrix is singular, e.g., zero on diagonal
       if (A(j,j) .EQ. 0) then
         logicFlag = .TRUE. !meaning the matrix is singular
         stop
       endif
       ! now loop over rows below row j 
       do i=j+1, Am
          B(i,:) = B(i,:) - A(i,j)*B(j,:)/A(j,j)
          A(i,:) = A(i,:) - A(i,j)*A(j,:)/A(j,j)   
       enddo
     enddo
      
     deallocate(temp)

  end subroutine gaussElim    

  !*********************************************************
  
  subroutine backSub(U, B, X)
  ! Takes an uppertriangular U
  ! A mx of n output vectors B
  ! returns A matrix containing n solution vectors X
  ! Solves UX = B

  implicit none
  real (kind(0.d0)), intent(in out)    :: U(:,:), B(:,:)
  real (kind(0.d0)), intent(out)       :: X(:,:)
  integer                              :: Am, An, Bm, Bn 
  !loop vars
  integer                              :: i,j,k
  real (kind(0.d0))                    :: tempSum   

  !pull dimensions from the varying matrices, m rows n columns
  !Am corresponds to the number of rows of U
  Am = size(U(:,1))
  An = size(U(1,:)) 
  Bm = size(B(:,1))
  Bn = size(B(1,:))
  !print *, Am, An, Bn, Bm
  ! loop through columns of B, solve UX = B one column at a time
  ! so column one of X contains the (vector) solution, call it x, to Ax = b, where b is the first column of B 
  X = 0.0
  do j=1, Bn 
     X(Bm,j) = B(Bm,j)/U(Am,Am) !since A and B have same # of rows, Am=Bm. !pretty sure this is right
     do i = Am-1, 1, -1 
       if (U(i,i) .EQ. 0) then
         stop
       endif
       tempSum = 0.0
       do k=i+1,Am
          tempSum = tempSum + U(i,k)*X(k,j)
       enddo
     X(i,j) = (B(i,j)-tempSum)/U(i,i)  
     enddo
  enddo
  end subroutine backSub
  
  !*********************************************************

  subroutine LU_decomp(A, m, logicFlag, s)
  
  implicit none
  real (kind(0.d0)), intent(in out)    :: A(:,:) ! matrix to be decomposed
  integer, intent(out)                 :: s(:)   ! permutation vector
  logical, intent(out)                 :: logicFlag
  integer                              :: m      ! m rows in A
  ! Loop/temp vars
  real (kind(0.d0))                    :: pivot
  real (kind(0.d0)), allocatable       :: temp(:)
  integer                              :: i,j,k,pos, tempC

    allocate(temp(m))
    logicFlag = .FALSE.
    s = 0.0

    ! initialize permutation vector
    do j=1, m
       s(j) = j
    enddo

    do j=1, m
    ! the below loop finds pivot, position of pivot (max in column, index of max)
       pivot = 0
       ! loop over columns in this row
       do i=j, m
         if (ABS(A(i,j)) > pivot) then 
            pivot =  ABS(A(i,j)) 
            pos = i
         endif
       enddo
       ! next we interchange rows if pivot position (pos) != row position (j)  
       ! this ensures max value is found on the diagonal
       if (pos .NE. j) then
         temp(:) = A(pos,:) 
         A(pos,:) = A(j,:)
         A(j,:) = temp(:)
         ! don't forget to make the same changes to S
         tempC = s(pos)
         s(pos) = s(j)
         s(j) = tempC
       endif
       if (A(j,j) .EQ. 0.d0) then
         logicFlag = .TRUE. !meaning the matrix is singular
         stop
       endif
       ! differs from gaussian  Elim here really; couple changes to indicies as well
       do i = j+1, m
          A(i,j) = A(i,j)/A(j,j)
          do k = j+1, m
             A(i,k) = A(i,k) - A(i,j)*A(j,k) 
          enddo
       enddo  
      enddo

       
      deallocate(temp)     
  end subroutine LU_decomp

  !*********************************************************

  subroutine backsubLU(A, B, m, s, X)
  implicit none !this was added and not tested, but feel like it should be here lol.
  real (kind(0.d0)), intent (in out)          :: A(:,:), B(:,:)
  integer, intent(in)                         :: m
  integer, intent(in)                         :: s(:)
  real (kind(0.d0)), intent(out)              :: X(:,:)
  !temp and loop vars
  integer                                     :: k, i, j
  !integer                                    :: tempS
  real (kind(0.d0)), allocatable              :: tempSum(:)
 
  
  allocate(tempSum(size(B,2)))
  !allocate(X(size(B,1),size(B,2)))


    ! permutes B in the same way s (and thus A) were permuted in LU decomp
    do j=1, m
       X(j,:) =  B(s(j),:)
    enddo 
    
    ! Now do forward substitution
    do j= 1, m-1
      do i= j+1, m
         X(i,:) = X(i,:) - X(j,:)*A(i,j) 
      enddo
    enddo

    ! Now do backwards substitution, looping over columns of B to do it for each one. 
    ! Outer loop does this
    !X = 0.0
       
     do i=m, 1, -1
       if (A(i,i) .EQ. 0) then
          stop !matrix is singular  
       endif
     tempSum = 0.0
       do k=i+1, m
         tempSum = tempSum + A(i,k)*X(k,:) 
       enddo
       X(i,:) = (X(i,:) -tempSum)/A(i,i)
     enddo
    

  deallocate(tempSum)  
  !deallocate(X)  
  end subroutine backsubLU

 !*********************************************************
 !
 ! Homework 3 subroutines below here
 !
 !*********************************************************


  subroutine cholDecomp(A, logicFlag)
  ! takes in a Hermitian positive definite matrix A
  
  implicit none  
  real (kind(0.d0)), intent(in out)        :: A(:,:)
  logical, intent(out)                     :: logicFlag
  ! loop and temp vars
  integer                                  :: m, n
  integer                                  :: i, j, k
    
    logicFlag = .FALSE.
    m = size(A,1) ! number of rows
    n = size(A,2) ! number of columns
    
    !note A is square, so m = n
    do j = 1, n !looping over columns
      ! computing new diagonal elements
      do k = 1, j-1
         A(j,j) = A(j,j) - A(j,k)*A(j,k) !want A^*(j,k).     
         ! note that A hermitian so could just do A(j,k)**2?
         ! leaving like this to be more general
      enddo
      A(j,j) = SQRT(A(j,j))
      ! now calculate elements below diagonal
      do i = j+1, m
         do k = 1, j-1
            A(i,j) = A(i,j) - A(i,k)*A(j,k) !again this conjugate may be off. Try this for now.
         enddo
         A(i,j) = A(i,j)/A(j,j)
      enddo
    enddo
  end subroutine cholDecomp

 !*********************************************************
 
 subroutine backsubChol(A, B, X)
 ! takes in a cholesky decomposed matrix A
 ! and an output mx B. 
 ! Finds the X in AX = B. 
    
    implicit none
    real(kind(0.d0)), intent(in out)  :: A(:,:)
    real(kind(0.d0)), intent(in out)  :: B(:)
    real(kind(0.d0)), intent(out)     :: X(:)
    ! loop and temp vars
    integer                               :: m, n
    real (kind(0.d0))                     :: tempSum
    integer                               :: i, j, k, L, d !d is degree of poly fitting
    real (kind(0.d0)), allocatable        :: Y(:)

    d = 6
    m = size(A,1)
    n = size(B)

   
    allocate(Y(d)) 

    ! need to intialize Y and X in an intelligent way
    ! in LUbacksub did X(k,:) = B(s(j),:), loop through j. 
    do i = 1, m
       Y(i) = B(i)
    enddo

    ! Foward Substitution: 
    ! add an outer loop to do all columns of X perhaps? (NO) 
    ! changed to handle only one column vector B. 
    do i = 1, m
       tempSum = B(i)
       do j = 1, i-1
          tempSum = tempSum - Y(j)*A(i,j)
       enddo
       Y(i) = tempSum/A(i,i)
    enddo

    ! Back Sub:
    do i = m, 1, -1
       if (A(i,i) .EQ. 0) then
          print *, 'Matrix is singular'
          stop
       endif        
       do k = i+1, m
          Y(i) = Y(i) - A(k,i)*X(k) !psuedo code uses X instead of Y here. 
       enddo                                    
       X(i) = Y(i)/A(i,i)
    enddo

   
    deallocate(Y)
    
 end subroutine backsubChol

 !*********************************************************

 subroutine vanderMaker(A, Vandermonde, y)
 ! takes in a matrix A, should be two columns corresponding to input / output
 ! returns Vandermonde matrix. This is plugged into cholesky decomp. 
 ! After backsub, the least squares problem is solved.

    implicit none
    real(kind(0.d0)), intent(in)               :: A(:,:)
    real(kind(0.d0)), intent(out)              :: Vandermonde(:,:)
    real(kind(0.d0)), intent(out)              :: y(:)
    integer                                    :: d   ! degree of fitting
    ! loop and temp vars
    integer                                    :: i,j ! loops
    integer                                    :: m,n ! dims   
    
    
    !for higher degree of polynomial fitting, increase d
    d = 6 
    
    ! store number of input/output pairs
    m = size(A(:,1))
    n = size(A(1,:)) !this better be 2
    if (n .NE. 2) then
       print *, "you dun goofed, I can't Vandermonde under these conditions."
    endif
    
    
    ! allocate/populate Vandermonde mx, number of columns is degree of fitting
    !allocate(Vandermonde(m, d))
    Vandermonde = 0.d0

    do i = 1, m
        do j = 1, d
           Vandermonde(i,j) = A(i,1)**(j-1) 
        enddo
    enddo

    
    ! time to populate the output vector, y
    !allocate(y(m))
    y = 0.d0
    
    do i = 1, m
        y(i) = A(i,2)
    enddo

    !deallocate(y)
    !deallocate(Vandermonde)

 end subroutine vanderMaker

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
    
    A = 0.e0
    do i = 1, m
       do j = 1, n
       A(i,j) = u(i)*v(j) 
       enddo
    enddo

 end subroutine   


 !*********************************************************
 
 subroutine householderQR(A, Q)
    ! takes in A, slowly transforms it into R, returns it
    ! also returns Q 
    
    implicit none
    real(kind(0.d0)), intent(in out)     :: A(:,:)
    real(kind(0.d0)), intent(out)        :: Q(:,:)
    ! loop and temp vars
    real(kind(0.d0)), allocatable        :: Vmat(:,:) 
    real(kind(0.d0)), allocatable        :: v(:)
    real(kind(0.d0)), allocatable        :: s(:)
    integer                              :: m, n !dims
    integer                              :: i, j, k
    real(kind(0.d0))                     :: norm
    real(kind(0.d0))                     :: normTwo 
    real(kind(0.d0)), allocatable        :: tempMat(:,:)
    real(kind(0.d0)), allocatable        :: vt(:) 
    ! the following declarations are used in creation of Q
    real(kind(0.d0)), allocatable        :: Id(:,:)
    real(kind(0.d0)), allocatable        :: H(:,:)


    m = size(A,1) ! m rows  
    n = size(A,2) ! n columns - should be degree of polynomial fitting for vanderMx

    allocate(Vmat(m,n))
    allocate(v(m))   
    allocate(s(n))
    allocate(tempMat(m, m))

    ! loop over columns of A
    do j=1, n 
       norm = 0.d0 
       !compute signed norm for each column 
       do i = j, m !psuedocode loops from j to m. Shouldn't the norm be 1 to m no matter what? 
          norm = norm + A(i,j)**2
       enddo
       norm = SQRT(norm)   
       s(j) = SIGN(norm, A(j,j))
       
       ! time to start populating v
       
       v = 0.d0 ! fill v full of zeros
       ! a do while loop to fill the remainder of v with the entries of A
       i = j
       do while (i<m+1)  ! this implies each v_j should be length m. We have n v_j
          v(i) = A(i,j) 
          ! hopefully this if statement assigns exactly the jth entry of v to be different
          if (i .EQ. j) then
             v(i) = A(j,j) + s(j)
          endif
          i = i+1
       enddo
   
       ! lmao just found the norm2 function so excited
       ! it requires double precision and we're required to do single :(
       normTwo = 0.d0
       do i = 1, m
          normTwo = normTwo + v(i)**2 
       enddo
       normTwo = SQRT(normTwo)
      
       v = v / normTwo
       ! I'ma store these bad boys I got the memory why not
       ! This aint the moon landing
       Vmat(:,j) = v
      
    
       call outerProd(v, v, tempMat) 
       ! initialize Q to be the identity
       ! apply tempMat to it each pass

       A = A - 2.d0*matmul(tempMat, A)
    enddo
    
    ! attempt to create H, giving Q
    ! this can be done using v
    ! first create ID mx, used in many steps of computation
    allocate(Id(m,m))
    allocate(H(m,m))

    Id = 0.d0
    do i = 1, m
       Id(i,i) = 1
    enddo

    ! Initialize Q to be identity
    Q = Id
    H = 0.d0

    ! compute each H, apply it to Q
    do j = 1, n
       call outerProd(Vmat(:,j), Vmat(:,j), tempMat)
       H = Id - 2.d0*tempMat
       Q = matmul(H, Q)    
    enddo


    Q = transpose(Q) ! this was added in assignment 4, for previous assignments Q is transposed.

    !deallocate(vt)
    deallocate(H)
    deallocate(Id)
    deallocate(tempMat) 
    deallocate(v)
    deallocate(Vmat)
 end subroutine

 !*********************************************************

 subroutine backsubQR(R, b, x)
    !takes in upper tri R made from HH
    ! b made as QT y
    ! x a solution vector
    implicit none
    real(kind(0.d0)), intent(in)   :: R(:,:)
    real(kind(0.d0)), intent(in)   :: b(:)
    real(kind(0.d0)), intent(out)  :: x(:) 
    ! temp and loop vars
    real(kind(0.d0))               :: tempSum
    integer                        :: i, k, m

    m = size(R,1)

    do i=m, 1, -1
       if (R(i,i) .EQ. 0) then
          stop !matrix is singular  
       endif
       tempSum = 0.d0
       do k=i+1, m
         tempSum = tempSum + R(i,k)*x(k) 
       enddo
       x(i) = (b(i) -tempSum)/R(i,i)
     enddo
  
 end subroutine

 !*********************************************************
 !
 ! Homework 4 subroutines below here
 !
 !*********************************************************

 ! this is the householder routine edited to give tri diag
 ! also changed to be double precision
 subroutine diagHH(A)
    ! takes in A, slowly transforms it into R, returns it
    ! also returns Q (hopefully if I can get it to work)
    
    implicit none
    real(kind(0.d0)), intent(in out)     :: A(:,:)
    ! loop and temp vars
    real(kind(0.d0)), allocatable        :: Vmat(:,:) 
    real(kind(0.d0)), allocatable        :: v(:)
    real(kind(0.d0)), allocatable        :: s(:)
    integer                              :: m, n !dims
    integer                              :: i, j, k
    real(kind(0.d0))                     :: norm
    real(kind(0.d0))                     :: normTwo 
    real(kind(0.d0)), allocatable        :: tempMat(:,:)
    real(kind(0.d0)), allocatable        :: vt(:) 
    ! the following declarations are used in creation of Q
    real(kind(0.d0)), allocatable        :: Id(:,:)
    real(kind(0.d0)), allocatable        :: H(:,:)


    m = size(A,1) ! m rows  
    n = size(A,2) ! n columns - should be degree of polynomial fitting for vanderMx

    allocate(Vmat(m,n))
    allocate(v(m))   
    allocate(s(n))
    allocate(tempMat(m, m))

    
    do j=1, m-1
       norm = 0.d0 
       !compute signed norm for each column 
       do i = j+1, m !psuedocode loops from j to m. Shouldn't the norm be 1 to m no matter what? 
          norm = norm + A(i,j)**2
       enddo
       norm = SQRT(norm)   
       s(j) = SIGN(norm, A(j+1,j))
       
       ! time to start populating v
       
       !v = 0.d0 ! fill v full of zeros
       ! a do while loop to fill the remainder of v with the entries of A

       do i = 1, m
          if (i<j+1) then
             v(i) = 0.d0
          elseif (i .EQ. j+1) then
             v(i) = A(i, i-1) + s(j)
          else
              v(i) = A(i,j)
          endif
       enddo    
       ! lmao just found the norm2 function so excited
       ! it requires double precision and we're required to do single :(
       normTwo = 0.d0
       do i = 1, m
          normTwo = normTwo + v(i)**2 
       enddo
       normTwo = SQRT(normTwo)
      
       v = v / normTwo
       ! I'ma store these bad boys I got the memory why not
       ! This aint the moon landing
       Vmat(:,j) = v
         
       call outerProd(v, v, tempMat) 
     
       ! apply tempMat to it each pass

       A = A - 2.d0*matmul(tempMat, A)
       A = A - 2.d0*matmul(A, tempMat)
    enddo
    
    deallocate(tempMat) 
    deallocate(v)
    deallocate(Vmat)
 end subroutine

 !*********************************************************

 subroutine eigenvalQR(A, V)
    ! takes square, tri diag matrix and finds eigenvalues
    ! replaces A with eigenvalues along main diag, other entries trend towards zero
    implicit none
    real(kind(0.d0)), intent(in out)   :: A(:,:)
    real(kind(0.d0)), intent(out)      :: V(:,:)  
    real(kind(0.d0)), allocatable      :: Id(:,:)
    real(kind(0.d0)), allocatable      :: Q(:,:)
    integer                            :: m, n
    ! temp and loop vars
    integer                            :: i, j, k
    real(kind(0.d0))                   :: thresh
    real(kind(0.d0)), allocatable      :: errorVec(:)
    real(kind(0.d0))                   :: error
    integer                            :: counter

    m = size(A, 1) ! m rows
    n = size(A, 2) ! n columns
    thresh = 10.d0**(-10.d0) !set thresh to some smol value, like 1E-8
 
    allocate(Id(m,m))
    allocate(Q(m,m))
    allocate(errorVec(m-1))

    Id = 0.d0
    do i=1, m
      Id(i,i) = 1
    enddo

    V = Id
    error = 1
    counter = 0

    do while (error > thresh)  ! do while error is too large.  
       call householderQR(A, Q)
       !call matprint(A, m, n) confrims A has been turned into upper triangular R
       A = matmul(A, Q)
       V = matmul(V, Q)  ! returns eigenvectors, if you're into that kinda thing
       
       ! compute error to exit loop 
       ! this takes norm of vector comprised of entries immediately under main diagonal
       ! when the above nerd is small enough, stop looping.
       errorVec = 0.d0 
       do j = 1, m-1
             errorVec(j) = A(j+1, j)
       enddo
       error = norm2(errorVec)
       counter = counter + 1 ! curious how many times it passes through the do while
    enddo
    
    !print *, "number of iterations required to reach a threshold of ", thresh
    !print *, counter

    deallocate(errorVec)
    deallocate(Q)
    deallocate(Id)
 end subroutine

 !*********************************************************

 subroutine shiftQR(A)
    ! takes square, tri diag matrix and finds eigenvalues
    ! replaces A with eigenvalues along main diag, other entries trend towards zero
    implicit none
    real(kind(0.d0)), intent(in out)   :: A(:,:)
    real(kind(0.d0)), allocatable      :: Q(:,:), Id(:,:)
    ! temp and loop vars
    integer                            :: i, j, k, m, n
    real(kind(0.d0)), allocatable      :: tempMat(:,:)
    real(kind(0.d0))                   :: mu
    real(kind(0.d0))                   :: thresh
    real(kind(0.d0)), allocatable      :: errorVec(:)
    real(kind(0.d0))                   :: subEntry
    integer                            :: counter

    m = size(A,1)
    n = size(A,2) !supposed to be square but why not

    allocate(Id(m,m))
    allocate(Q(m,m))
    allocate(tempMat(m,m))
    
    Id = 0.d0
    do i=1, m
      Id(i,i) = 1
    enddo
    
    thresh = 10.d0**(-10.d0) !set thresh to some smol value, like 1E-10   
    subEntry = 1.d0
    tempMat = 0.d0


    do while (subEntry > thresh)    ! if you set this bound to 10, you will get hella NaN
       ! set shift
       mu = A(m,m)
       ! compute shifted mat
       tempMat = A - mu*Id
       !call matPrint(tempMat, m, m) 

       ! QR factor shifted mat, turns tempMat into R
       
       call householderQR(tempMat, Q)
       !call matPrint(tempMat, m, m)
       !call matPrint(Q, m, n)

       ! Update A with RQ + mu*I
       A = matmul(tempMat, Q) + mu*Id

       ! compute error to get out of this trap
       ! looks at final sub diagonal entry
       subEntry = ABS(A(m, m-1))
    enddo

    deallocate(tempMat)
    deallocate(Q)
    deallocate(Id)
 end subroutine 

!*********************************************************
 subroutine eigvecIter(A, mu, x)
 ! takes in a matrix A, eigenvalue guess (NOT exact eig val).
 ! returns eigenvector
    implicit none 
    real(kind(0.d0)), intent(in out) :: A(:,:)
    real(kind(0.d0)), intent(out)    :: x(:,:)
    real(kind(0.d0)), intent(in)     :: mu
    real(kind(0.d0)), allocatable    :: Id(:,:), B(:,:)
    ! temp and loop vars
    integer                          :: m, n, i
    real(kind(0.d0)), allocatable    :: y(:,:)  
    logical                          :: logicFlag
    integer, allocatable             :: s(:)

    m = size(A,1)
    n = size(A,2) !supposed to be square but why not

    allocate(Id(m,m))
    allocate(B(m,m))
    allocate(y(m,1))
    allocate(s(m))

    Id = 0.d0
    do i=1, m
      Id(i,i) = 1
    enddo
    
    ! populate initial vector, must have norm 1
    x = 0.d0
    x(m,1) = 1 
    !print *, x !confirms x is populated correctly 
    
    B = A - mu*Id
    
    call LU_decomp(B, m, logicFlag, s)

    !print *, "A = "
    !call matprint(A, m, m)
    !print *, "B = "
    !call matPrint(B, m, m)   !confirms B is populated correctly for initial pass
   
    do i=1, 10
       x = x/norm2(x(:,1)) 
       !solves for y in By = x
       call backsubLU(B, x, m, s, y) 
       y(:,1) = y(:,1) / norm2(y(:,1))    
       x(:,1) = y(:,1)
    enddo

   deallocate(s)
   deallocate(y)
   deallocate(Id)
   deallocate(B) 

 end subroutine


end module LinAl
