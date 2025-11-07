Program Driver_LinAl

  use LinAl

  implicit none  
  character(len=100)                 :: myFileName
  character(len=100)                 :: myFileName2
  integer                            :: i,j
  ! problem 1 type declarations below
  real(kind(0.d0)), dimension(4,4)      :: A
  real(kind(0.d0)), dimension(1:16)     :: b = (/5, 4, 1, 1, 4, 5, 1, 1, 1, 1, 4, 2, 1, 1, 2, 4/) 
  ! problem 2 type declarations below
  real(kind(0.d0)), dimension(3,3)      :: A2 
  real(kind(0.d0)), dimension(3,3)      :: Atwo
  real(kind(0.d0)), dimension(1:9)      :: b2 = (/3, 1, 0, 1, 2, 1, 0, 1, 1/)
  real(kind(0.d0)), allocatable         :: V(:,:)
  ! problem 3 type declarations below
  real(kind(0.d0))                      :: mu
  real(kind(0.d0)), dimension(1:16)     :: b3 = (/2, 1, 3, 4, 1, -3, 1, 5, 3, 1, 6, -2, 4, 5, -2, -1/)
  real(kind(0.d0)), dimension(4,4)      :: C
  real(kind(0.d0)), allocatable         :: x(:,:)
  A = reshape(b, (/4, 4 /))
  A2 = reshape(b2, (/3, 3 /)) 
  C = reshape(b3, (/4, 4/))

  print *, "The given matrix has the following representation and dimensions before reduction:"
  call matPrint(A, 4, 4)
  call diagHH(A)

  print *, "The above matrix and its dimensions after reduction to tridiagonal form are"
  call matPrint(A, 4, 4)
  

  !-------------------------------------
  !
  ! Problem 2 type driver routine below:
  !
  !-------------------------------------
  print *, "---------------------------- "
  print *, " "
  print *, "Problem 2 outputs to follow: "

  allocate(V(size(A2,1), size(A2,2)))
 
  print *, "Consider the following matrix"
  call matPrint(A2, 3, 3)
  
  !call diagHH(A2) not needed, matrix is already tridiagonal

  print *, "The eigenvalues of the above matrix, found without shift, are located along the main diagonal of: "
  call eigenvalQR(A2, V) 
  
  call matprint(A2, 3, 3)
  
  ! Now do it with shift
  print *, "Part (ii):"

  Atwo = reshape(b2, (/3, 3/))
  
  call shiftQR(Atwo)

  call matPrint(Atwo, 3, 3)
 
 

  !-------------------------------------
  !
  ! Problem 3 type driver routine below:
  !
  !-------------------------------------
  print *, "---------------------------- "
  print *, " "
  print *, "Problem 3 outputs to follow: "
  allocate(x(size(C, 1), 1))
   

  mu = 5.669
  print *, "The (incorrect) eigenvector assosciated with eigenvalue ", mu, "is"
  call eigvecIter(C, mu, x)
  print *, x
  print *, " "

  mu = -8.028
  print *, "The (incorrect) eigenvector assosciated with eigenvalue ", mu, "is"
  call eigvecIter(C, mu, x)
  print *, x
  print *, " "

  mu = 7.932
  print *, "The (incorrect) eigenvector assosciated with eigenvalue ", mu, "is"
  call eigvecIter(C, mu, x)
  print *, x
  print *, " "

  mu = -1.573
  print *, "The (incorrect) eigenvector assosciated with eigenvalue ", mu, "is"
  call eigvecIter(C, mu, x)
  print *, x
  print *, " "




  deallocate(x)
  deallocate(V)
end Program Driver_LinAl
