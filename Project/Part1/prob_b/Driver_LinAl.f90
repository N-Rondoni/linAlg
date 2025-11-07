! A driver routine to answer the prob_b.1 questions
Program Driver_LinAl

  use LinAl
  implicit none
  real(kind(0.d0)), allocatable  :: A(:,:)
  real(kind(0.d0)), allocatable  :: b(:)            ! counter is number of iterations to desired accuracy
  integer                        :: i, D, ldA, counter ! respectively diagonal entry, dimension of mx (lda x lda)
  real(kind(0.d0)), allocatable  :: x(:)
  real(kind(0.d0))               :: accuracy

  ! Leading dim of A, 
  ldA = 10
  accuracy = 5.d0**(-10.d0)

  allocate(A(ldA, ldA))
  allocate(b(ldA))   
  allocate(x(ldA))
  
  ! This should have been a subroutine. I printed in a weird way so I could work with the outputs easier.
  print *, "Required outputs for b.1:  "
  print *, " "

  print *, "Jacobi outputs: iterations, error for D = 2, 5, 10, 100, 1000"

  D = 2
  !call inputMatMaker(D, A, ldA, b) 
  !!call Jacobi(A, x, b, accuracy, counter)
  !print *, counter, norm2(b - matmul(A, x))
  print *, "Jacobi doesn't converge for D=2." 
  
  D = 5 
  !call inputMatMaker(D, A, ldA, b)
  !call Jacobi(A, x, b, accuracy, counter)
  !print *, counter, norm2(b - matmul(A, x))
  print *, "Jacobi doesn't converge for D=5."

  D = 10 
  call inputMatMaker(D, A, ldA, b)
  call Jacobi(A, x, b, accuracy, counter)
  print *, counter, norm2(b - matmul(A, x))


  D = 100
  call inputMatMaker(D, A, ldA, b)
  call Jacobi(A, x, b, accuracy, counter)
  print *, counter, norm2(b - matmul(A, x))
  
  D = 1000
  call inputMatMaker(D, A, ldA, b)
  call Jacobi(A, x, b, accuracy, counter)
  print *, counter, norm2(b - matmul(A, x))
 
  print *, "Seidel outputs: iterations, error for D = 2, 5, 10, 100"
  D = 2
  call inputMatMaker(D, A, ldA, b) 
  call Seidel(A, x, b, accuracy, counter)
  print *, counter, norm2(b - matmul(A, x))
   
  D = 5 
  call inputMatMaker(D, A, ldA, b)
  call Seidel(A, x, b, accuracy, counter)
  print *, counter, norm2(b - matmul(A, x))
  
  D = 10 
  call inputMatMaker(D, A, ldA, b)
  call Seidel(A, x, b, accuracy, counter)
  print *, counter, norm2(b - matmul(A, x))
  
  D = 100
  call inputMatMaker(D, A, ldA, b)
  call Seidel(A, x, b, accuracy, counter)
  print *, counter, norm2(b - matmul(A, x))
 
  D = 1000
  call inputMatMaker(D, A, ldA, b)
  call Seidel(A, x, b, accuracy, counter)
  print *, counter, norm2(b - matmul(A, x))


  !Modify the program s.t. A(i,i) = i, 1 else.  
  print *, "A(i,i) modified to have i on main diagonal, 1 else:"
  
  do i =1, ldA
     A(i,i) = i
  enddo
  
  !call Jacobi(A, x, b, accuracy, counter)
  !print *, counter, norm2(b - matmul(A, x))
  print *, "Jacobi doesn't converge in this case "

  print *, "Seidel however, converges: "
  call Seidel(A, x, b, accuracy, counter)
  print *, counter, norm2(b - matmul(A, x))
   
  !###############################################
  !                                              #
  !   Driver routine for b.2 to follow.          #
  !                                              #
  !                                              #
  !################################################

  print *, "_______________________________________ "
  print *, " "
  print *, "Problem b.2 outputs:  "
  print *, "Conjugate Gradient method: iterations, error for D = 2, 5, 10, 100, 1000"
  
  LdA = 10
  accuracy = 5.d0**(-10.d0)
  
  D = 2
  call inputMatMaker(D, A, ldA, b)
  call conjGrad(A, x, b, accuracy, counter)
  print *, counter, norm2(b - matmul(A, x)) 

  D = 5
  call inputMatMaker(D, A, ldA, b)
  call conjGrad(A, x, b, accuracy, counter)
  print *, counter, norm2(b - matmul(A, x)) 

  D = 10
  call inputMatMaker(D, A, ldA, b)
  call conjGrad(A, x, b, accuracy, counter)
  print *, counter, norm2(b - matmul(A, x)) 

  D = 100
  call inputMatMaker(D, A, ldA, b)
  call conjGrad(A, x, b, accuracy, counter)
  print *, counter, norm2(b - matmul(A, x)) 

  D = 1000 
  call inputMatMaker(D, A, ldA, b)
  call conjGrad(A, x, b, accuracy, counter)
  print *, counter, norm2(b - matmul(A, x)) 

  !Modify the program s.t. A(i,i) = i, 1 else.  
  print *, "A(i,i) modified to have i on main diagonal, 1 else:"
  print *, "if A is a 10x10 matrix, we have iteration, error "
  do i =1, ldA
     A(i,i) = i
  enddo
  
  call conjGrad(A, x, b, accuracy, counter) 
  print *, counter, norm2(b - matmul(A, x))
  
  deallocate(A)
  deallocate(b)
  deallocate(x)
  
  print *, "A(i,i) still modified to have i on main diagonal, 1 else: "
  print *, "Now A is a 100x100 matrix, yielding iteration, error "
  ldA = 100

  allocate(A(ldA, ldA))
  allocate(b(ldA))   
  allocate(x(ldA))
  
  call inputMatMaker(D, A, ldA, b)

  do i = 1, ldA
     A(i,i) = i
  enddo 
  
  call conjGrad(A, x, b, accuracy, counter)
  print *, counter, norm2(b -  matmul(A,x))
  

  deallocate(x)
  deallocate(A)
  deallocate(b)
end Program Driver_LinAl
