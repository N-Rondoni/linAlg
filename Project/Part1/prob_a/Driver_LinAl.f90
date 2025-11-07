! A driver routine to answer the prob_a questions
Program Driver_LinAl

  use LinAl

  implicit none
  
  character(len=100)               :: myFileName, fileName2!, fileNameTest
  integer                          :: i
  !prob 1_a  type declarations below
  real(kind(0.d0)), allocatable   :: dogData(:,:), dogData2(:,:)
  real(kind(0.d0)), allocatable   :: s(:) !singular values in a vector length m
  integer                         :: ldu, ldvt, lwork, info
  real(kind(0.d0)), allocatable   :: U(:,:)
  real(kind(0.d0)), allocatable   :: Vt(:,:)
  real(kind(0.d0)), allocatable   :: Work(:) !they call this the work array. Should be rank 2? 
  real(kind(0.d0)), allocatable   :: Appx(:,:)


  ! pull in dog (data)
  myFileName = 'dog_bw_data.dat'
  !myFileName = 'dog_bw_data_lowRes_811x1280.dat'

  open(10,file=myFileName)
  read(10,*) msize,nsize
  close(10)

  allocate(mat(msize,nsize))
  ! Always initialize with zeros
  mat = 0.d0

  call readMat(myFileName)
  
  allocate(dogData(msize, nsize))
  allocate(dogData2(msize, nsize))

  dogData = mat
  dogData2 = dogData

  deallocate(mat)

  !print *, size(dogData, 1)
  !print *, size(dogData, 2)
  !fileNameTest = 'testDog.dat'
  !call write_output(dogData, fileNameTest) confirms write_output isn't scrambling my life

  ! Matrix is loaded in, time to start driving.
  ! desired leading dimension of U
  ldu = msize
  ! desired leading dimension of Vt 
  ldvt = nsize
  ! set lwork = -1 to query workspace size
  lwork = -1
  
  allocate(U(ldu, ldu)) 
  allocate(Vt(ldvt, ldvt))
  allocate(work(1)) !or just work(1)? 
  allocate(s(msize))

  call dgesvd('A', 'A', msize, nsize, dogData, msize, s, U, ldu, Vt, ldvt, work, lwork, info)  
      
  lwork = nint(work(1)) 
  deallocate(work)

  allocate(work(lwork))
  
  call dgesvd('A', 'A', msize, nsize, dogData, msize, s, U, ldu, Vt, ldvt, work, lwork, info)  
  
  deallocate(work)
  
   
  Allocate(Appx(msize,nsize))
  ! print statements will be sent to output.txt, which can be read by python
  ! if you plan to work with the errors, make sure to comment out the below print routine for singular values.
  ! The errors produced from these commented out print commands are stored in errors.dat

  call reconstruct2(U, s, Vt, Appx, 1, 20)
  call write_output(Appx, 'dog_20.dat')
  !print *, 20,  (norm2(DogData2 - Appx)/(msize*nsize)) 

  call reconstruct2(U, s, Vt, Appx, 21, 40)
  call write_output(Appx, 'dog_40.dat')
  !print *, 40,  (norm2(DogData2 - Appx)/(msize*nsize)) 

  call reconstruct2(U, s, Vt, Appx, 41, 80)
  call write_output(Appx, 'dog_80.dat')
  !print *, 80,  (norm2(DogData2 - Appx)/(msize*nsize)) 

  call reconstruct2(U, s, Vt, Appx, 81, 160)
  call write_output(Appx, 'dog_160.dat')
  !print *, 160,  (norm2(DogData2 - Appx)/(msize*nsize)) 

  call reconstruct2(U, s, Vt, Appx, 161, 320)
  call write_output(Appx, 'dog_320.dat')
  !print *, 320,  (norm2(DogData2 - Appx)/(msize*nsize)) 

  call reconstruct2(U, s, Vt, Appx, 321, 640)
  call write_output(Appx, 'dog_640.dat')
  !print *, 640,  (norm2(DogData2 - Appx)/(msize*nsize)) 

  call reconstruct2(U, s, Vt, Appx, 641, 1280)
  call write_output(Appx, 'dog_1280.dat')
  !print *, 1280,  (norm2(DogData2 - Appx)/(msize*nsize)) 

  call reconstruct2(U, s, Vt, Appx, 1281, 2560)
  call write_output(Appx, 'dog_2560.dat')
  !print *, 2560,  (norm2(DogData2 - Appx)/(msize*nsize)) 

  call reconstruct2(U, s, Vt, Appx, 2561, 3355)
  call write_output(Appx, 'dog_3355.dat')
  !print *, 3355,  (norm2(DogData2 - Appx)/(msize*nsize)) 



  ! finally write out singular values, to be reported in report
  print *, "The singular values asked for are: "
  do i = 1, 10
     print *, s(i)
  enddo
  print *, s(20)
  print *, s(40)
  print *, s(80)
  print *, s(160)
  print *, s(320)
  print *, s(640)
  print *, s(1280)
  print *, s(2560)
  print *, s(3355)


  deallocate(Appx)
  deallocate(U)
  deallocate(Vt)
  deallocate(s)
  deallocate(dogData)
End Program Driver_LinAl
