module LinAl
  implicit none
  
  integer, save :: msize, nsize
  real, dimension(:,:), allocatable, save :: mat
  

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
    real (kind(0.d0)), intent(in)               :: A(:,:)
    integer, intent(in)                         :: m, n
    ! Loop vars
    integer :: i, j

    
    print *, m,n 
    do i = 1,m
       print *, A(i,:)
       ! if the above doesn't work, try A(:,i)
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
  ! A matrix containing n solution vectors X
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

  subroutine backsubLU(A, B, m, s, X)
  real (kind(0.d0)), intent (in out)          :: A(:,:), B(:,:)
  integer, intent(in)                         :: m
  integer, intent(in)                         :: s(:)
  real (kind(0.d0)), intent(out)              :: X(:,:)
  !temp and loop vars
  integer                                     :: k, i, j
  !integer                                     :: tempS
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
end module LinAl
