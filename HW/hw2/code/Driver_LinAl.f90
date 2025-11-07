Program Driver_LinAl

  use LinAl !only: mat, msize, nsize, readMat

  implicit none
  
  character(len=100) :: myFileName, myFileName2
  integer                                 :: i,j
  real                                    :: traceA, norm
  ! Declarations of type below here are for section 3
  real (kind(0.d0)), allocatable, save    :: As(:,:), Bs(:,:)
  integer, save                           :: Am, An, Bm, Bn
  logical                                 :: isSingular
  real (kind(0.d0)), allocatable, save    :: X(:,:)
  real (kind(0.d0)), allocatable, save    :: E(:,:) 
  ! Declaration of type below here are used in section 4
  real (kind(0.d0)), allocatable          :: Afour(:,:), Bfour(:,:) 
  logical                                 :: SingularOrNah
  integer, allocatable                    :: s(:)        
  real (kind(0.d0)), allocatable          :: Xfour(:,:)
  real (kind(0.d0)), allocatable, save    :: Efour(:,:)
  ! Declaration of type below here are used in section 5
  real (kind(0.d0))                       :: planeMat(3,3)
  real (kind(0.d0))                       :: planeVec(3,3)


  !Begin actually driving, Section 2
  myFileName = 'Amat.dat'

  open(10,file=myFileName)
  read(10,*) msize,nsize
  close(10)

  allocate(mat(msize,nsize))
  

  ! Always initialize with zeros
  mat = 0.0
  
  print *, ' '
  print *, 'Results from section 2: Warming-up ' 

  call readMat(myFileName)
   
  do i = 1, msize
     write(*,*) (mat(i,j) , j = 1, nsize )
  end do
  
  
  !call matPrint(mat, msize, nsize)
  ! apparently the above is not needed, the write statements provided print out as desired
  ! Leaving in place to serve as an example of a way to print without dimensions

  call computeTrace(mat, msize, traceA)
  print *,'Trace of above matrix is:', traceA

  do j=1,nsize
    call euclidNorm(mat(:,j), msize, norm)    
    print *, 'Euclidean norm of column', j, 'is', norm
  enddo
  
  
  
  
  
  allocate(As(msize, nsize)) 
  As = mat
  An = nsize
  Am = msize
  
  !print *, As
  !print *, An 
  !print *, Am 

  deallocate(mat) 
   
  !print *, As : print statements confirm As is saved despite mat being deallocated


  ! #########################################################################
  !                                                                        ##
  ! The above answers section 2                                            ##
  !                                                                        ##
  ! Driver routine for section 3 below:                                    ##
  !                                                                        ##
  ! #########################################################################

  print *, "______________________________"
  print *, ' ' 
  print *, 'Section 3 results: '
    
  myFileName2 = 'Bmat.dat'
        
  open(10,file=myFileName2)
  read(10,*) msize,nsize
  close(10)

  allocate(mat(msize,nsize))
  

  !Always initialize with zeros
  mat = 0.0
  
  call readMat(myFileName2)
  
  allocate(Bs(msize,nsize))
  Bs = mat
  Bn = nsize
  Bm = msize
  
  

  deallocate(mat)
  !print *, Bs

  print *, 'Matrix A has the following dimensions and representation before gaussian elimination'
  call matPrint(As, Am, An)
  
  print *, 'Matrix B has the following dimensions and representation before gaussian elimination'
  call matPrint(Bs, Bm, Bn)
  ! if you don't want the dimensions to be shown, use the following instead
  !do i = 1, Bmsize
     !write(*,*) (Bmat(i,j) , j = 1, nsize )
  !end do
  
  !sneak another assignment in here to use for the routine in section 4
  ! doing this since gauss elim writes over A, B
  allocate(Afour(Am, An))
  allocate(Bfour(Bm, Bn))
  Afour = As
  Bfour = Bs


  call gaussElim(As, Bs, Am, An, Bm, Bn, isSingular)  
  
  print *, 'Matrix A has the following dimensions and representation after gaussian elimination'
  call matPrint(As, Am, An)

  print *, 'Matrix B has the following dimensions and representation after gaussian elimination'
  call matPrint(Bs, Bm, Bn)
 
  
  ! Finally compute solution via backsub
  allocate(X(Am, Bn)) 
  call backSub(As, Bs, X)
  ! with the call, the solution matrix X has been written.
  
  print *, 'The solution matrix X to AX=B has the following dimensions and representation: (note Ax = b, column vectors x,b)'
  call matPrint(X, Am, Bn) 
  
  ! define the error matrix, print it to screen. E = As*X - Bs
  allocate(E(Am, Bn))
  print *, 'The Error Matrix E = AX - B and its dimensions are'
  E(:,:) = matmul(As(:,:), X(:,:)) - Bs(:,:)
  call matPrint(E, Am, Bn)
  ! now compute/print the norms of each column of error matrix  
  
  do j=1,nsize
    call euclidNorm(E(:,j), Bn, norm)    
    print *, 'Euclidean norm of column', j, 'is', norm
  enddo

  ! #########################################################################
  !                                                                        ##
  ! The above answers section 3                                            ##
  !                                                                        ##
  ! Driver routine for section 4 below:                                    ##
  !                                                                        ##
  ! #########################################################################
  print *, "______________________________"
  print *, ' ' 
  print *, 'Section 4 results: '
  print *, ' '
  print *, 'Matrix A has the following dimensions and representation before LU decomposition '
  call matPrint(Afour, Am, An)


  allocate(s(Am))
  call LU_decomp(Afour, Am, SingularOrNah, s)
  print *, 'Matrix A has the following dimensions and representation after LU decomposition'
  print *, 'Note: the lower tri is taken to be the lower part of A with ones on the diagonal'
  print *, 'Thus U is stored in the upper part of the below mx (everything on and above the diagonal).  '
  call matPrint(Afour, Am, An)
  !print *, s  from this you can reverse intuit the permutation matrix
  ! for instance, s in the 4x4 case always starts as [1,2,3,4]. If post algorithm s -> [3,4,2,1]
  ! simply see what matrix permutes these values in that way
  ! then it follows that PA = LU 
  
  allocate(Xfour(Bm, Bn)) 


  print *, 'After backsubstitution we have the following solution matrix'
  call backsubLU(Afour, Bfour, Am, s, Xfour)
  call matPrint(Xfour, Bm, Bn)

  
  ! define the error matrix, print it to screen. E = As*X - Bs
  allocate(Efour(Am, Bn))
  print *, 'The Error Matrix E = AX - B and its dimensions are'
  Efour(:,:) = matmul(As(:,:), Xfour(:,:)) - Bs(:,:)
  call matPrint(Efour, Am, Bn)
  ! now compute/print the norms of each column of error matrix  
  
  do j=1,nsize
    call euclidNorm(Efour(:,j), Bn, norm)    
    print *, 'Euclidean norm of column', j, 'is', norm
  enddo








  ! #########################################################################
  !                                                                        ##
  ! The above answers section 4                                            ##
  !                                                                        ##
  ! Driver routine for section 5 below:                                    ##
  !                                                                        ##
  ! #########################################################################
  
  !don't judge me for entering these directly Ian I couldn't get the reader to work on another file and it's midnight
  
  planeMat(1,1) = 1
  planeMat(1,2) = 2
  planeMat(1,3) = 3
  planeMat(2,1) = -3
  planeMat(2,2) = 2
  planeMat(2,3) = 5
  planeMat(3,1) = 2*ACOS(0.d0) 
  planeMat(3,2) = EXP(1.d0)
  planeMat(3,3) = -SQRT(2.d0)
  !call matPrint(planeMat,3,3)
  
  planeVec(1,1) = 1
  planeVec(2,1) = 1
  planeVec(3,1) = 1
  planeVec(1,2) = 0
  planeVec(1,3) = 0
  planeVec(2,2) = 0
  planeVec(2,3) = 0
  planeVec(3,2) = 0
  planeVec(3,3) = 0

  !print *, planeVec
  print *, "______________________________"
  print *, ' ' 
  print *, 'Section 5 results: '
  print *, ' '
  !print *, "the first column of the below matrix represents the coefficients for the plane equation containing the three given points"
  call gaussElim(planeMat, planeVec, 3, 3, 3, 3, isSingular)
  call backSub(planeMat, planeVec, X)
 
  print *, 'Coefficients of plane equation: (ignore the zero, not sure why that is printing)'
  print *, X(:,1)






  deallocate(Xfour)
  deallocate(Bfour)
  deallocate(s)
  deallocate(Afour)
  deallocate(E)
  deallocate(As)
  deallocate(Bs)
  deallocate(X)
End Program Driver_LinAl
