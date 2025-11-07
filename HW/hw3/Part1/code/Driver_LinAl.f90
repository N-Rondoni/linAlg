Program Driver_LinAl

  use LinAl

  implicit none
  
  character(len=100)                 :: myFileName
  character(len=100)                 :: myFileName2
  integer                            :: i,j
  ! section 1 type declarations below
  real(kind(0.e0)), allocatable   :: lsData(:,:)
  real(kind(0.e0)), allocatable   :: vanderMX(:,:)
  real(kind(0.e0)), allocatable   :: y(:)
  integer                            :: d !degree of poly fitting in vanderMX
  real(kind(0.e0)), allocatable   :: A(:,:)
  real(kind(0.e0)), allocatable   :: b(:)
  logical                         :: logicFlag
  real(kind(0.e0)), allocatable   :: X(:)
  real(kind(0.e0))                :: temp
  real(kind(0.e0)), allocatable   :: solcheck(:)
  real(kind(0.e0))                :: norm
  real(kind(0.e0)), allocatable   :: classicError(:)
  ! section 2 type declarations below
  real(kind(0.e0)), allocatable   :: Q(:,:)
  real(kind(0.e0)), allocatable   :: vanderMXtwo(:,:) 
  real(kind(0.e0)), allocatable   :: Id(:,:)
  real(kind(0.e0)), allocatable   :: tempDif(:,:) 
  real(kind(0.e0)), allocatable   :: tempDifQ(:,:)
  real(kind(0.e0)), allocatable   :: xSol(:)
  real(kind(0.e0)), allocatable   :: jankBmat(:,:)
  real(kind(0.e0)), allocatable   :: Rhat(:,:)
  real(kind(0.e0)), allocatable   :: Qhat(:,:)
  real(kind(0.e0)), allocatable   :: outVec(:)
  !This portion brings in the least squares data, to be used to create Vandermonde matrix.
  


  myFileName = 'least_squares_data.dat'

  open(10,file=myFileName)
  read(10,*) msize,nsize
  close(10)

  allocate(mat(msize,nsize))
  ! Always initialize with zeros
  mat = 0.0

  call readMat(myFileName)
  
  allocate(lsData(msize, nsize))
  
  lsData = mat

  deallocate(mat)
  !call matPrint(lsData, 21, 2) confirms lsData contains what was read

  
  ! with least squares data here, time to make Vandermonde matrix
  ! this hardcoded d can be changed to desired degree of poly fitting
  ! also must be changed in module, where it is named d still.
  ! d is found in the vanderMaker and backsubChol  modules. Change them there too.

  !d = 4 !uncomment this and comment the below for a degree 3 poly
  d = 6
  allocate(vanderMX(msize,d))  
  allocate(vanderMXtwo(msize,d))
  allocate(y(msize))

  call vanderMaker(lsData, vanderMX, y)
  vanderMXtwo = vanderMX ! this is done so I have an extra copy to compar QR against
  ! needs to happen since my QR changed the input matrix (in this case vandermx)

  !call matPrint(vanderMX, msize, d)
  ! print *, y
  ! the above prints confirm vanderMx, y are populated correctly. 
  ! matPrint won't take y cause it is nerd.
  
  ! create the normal equations to be solved
  allocate(A(d,d))
  allocate(b(d))
  

  A = matmul(transpose(vanderMX),vanderMX)
  b = matmul(transpose(vanderMX),y)
  
 
 
  call cholDecomp(A, logicFlag)
  !call matPrint(A, d, d) 
  
  allocate(X(d))
 
  call backsubChol(A, b, X)
  !with this call, the solution vector X has been created.
 


  print *, "The coefficients of the line of best fit are (in ascending order of power)"
  print *, X
  print *, 'This polynomial is degree'
  print *, (d-1)
  print *, 'These coefficients have also been written out to coeff.dat'
  
  myFileName2 = 'coeff.dat'
  open(20, file=myFileName2, status='replace')
     write(20,*) X
  close(20)
 
  allocate(classicError(msize))
  ! though the below looks like a mess, it is really just matmul(A,X) - b
  ! doing it this way since A gets updated in cholDecomp
  ! as such the original definition of A is reused here.
  ! this is error for the normal equations
  classicError = matmul(matmul(transpose(vanderMX),vanderMX), X) - b
  call euclidNorm(classicError, msize, norm)
  print *, 'The 2-norm of Ax - b, representing our error, is'
  print *, norm 
  

  ! this is error for the original set up, not normal equations!
  allocate(solcheck(msize))
  solcheck = matmul(vanderMX, X) - y
  call euclidNorm(solcheck, msize, norm)
  print *, 'The 2-norm error between the fitted curve and data is'
  print *, norm
  
  
  ! --------------------------------------------------------------
  !                                                              -
  !    Section 2 Driver routine follows:                         -
  !                                                              -
  ! --------------------------------------------------------------



  allocate(Q(msize,msize)) 
  print *, '---------------------------------------------'
  print *, 'Section 2 results:'
  print *, ' '
  print *, 'The matrix A - QR has dimensions and representation ' 
  call householderQR(vanderMX, Q)
  allocate(tempDif(msize, d))
  tempDif = vanderMXtwo - matmul(transpose(Q), vanderMx) !This Q is transposed. A = Qt vanderMX
  call matPrint(tempDif, msize, d) 
  
  
  print *, ' '
  
  do j=1, d
     print *, 'The Euclidean norm of column', j, 'is',  norm2(tempDif(:,j))
  enddo
  print *, ' '


  allocate(tempDifQ(msize,msize))
  allocate(Id(msize, msize))
  
  Id = 0.e0
  do j = 1, msize
     Id(j,j) = 1
  enddo

  tempDifQ = matmul(Q,transpose(Q)) - Id
  print *, "The matrix Q^TQ - I has the following dimensions and representation "
  call matPrint(tempDifQ, msize, msize) ! confirms Q is created accurately.
  

  print *, ' '
  do j = 1, msize
     print *, 'The Euclidean norm of column', j, 'is', norm2(tempDifQ(:,j))
  enddo
  


  ! solve the least square problem by looking at the first n lines of Rx = Q^t b
  ! backsub on these n lines
  allocate(xSol(d)) !was d by d
  allocate(jankBmat(d,d)) !I'd rather turn b into a jank matrix than rewrite backsub don't @ me
    
 
  allocate(Qhat(msize,d))
  allocate(Rhat(d,d))
  ! too hard to continue thinking of Q as already transposed.
  !print *, 'the size of Q is'
  !print *, size(Q,1)
  !print *, size(Q,2)
  Q = transpose(Q)


  Qhat = (Q(:, 1:d)) !pulls the first through dth columns
                     
  !print *, 'the size of Qhat is'
  !print *, size(Qhat,1)
  !print *, size(Qhat,2)

  Rhat = vanderMX(1:d,:) 

  allocate(outVec(msize))
  outVec = matmul(transpose(Qhat), y)
  !print *, outVec
  


  !print *, size(outVec,1)
  !print *, size(Rhat, 1)
  !print *, size(Rhat, 2)

  jankBmat = 0.e0
  do i = 1, d
      jankBmat(i,1) = outVec(i)
  enddo
  


  !call backsub(Rhat, jankBmat(1:d,1:d), xSol)
  call backsubQR(Rhat, outVec, xSol)
  print *, 'The coefficients of the least squares equations projected onto the span of A are '
  print *, xSol



  print *, 'The 2norm error Rx - Qtb is'
  print *, norm2(matmul(Rhat, xsol) - outVec)
  
  deallocate(outVec)
  deallocate(Qhat)
  deallocate(Rhat)
  deallocate(jankBmat)
  deallocate(xSol) 
  deallocate(Id)
  deallocate(tempDifQ)
  deallocate(tempDif)
  deallocate(Q)
  deallocate(classicError)
  deallocate(solcheck)
  deallocate(X)
  deallocate(b)
  deallocate(A)
  deallocate(vanderMX)
  deallocate(vanderMXtwo)
  deallocate(y)  
  deallocate(lsData)
End Program Driver_LinAl
