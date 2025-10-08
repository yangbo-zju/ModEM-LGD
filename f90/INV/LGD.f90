module LGD

! inherits SensComp, DataIO and all modules they use. Also Main_MPI and Sub_MPI

use invcore


implicit none

public  :: LGDsolver

! iteration control for the LGD solver is initialized once
! and saved in the module to be used by most subroutines

  type  :: LGDiterControl_t
     ! maximum number of iterations in one call to iterative solver
     integer            :: maxIter
     ! convergence criteria: return from solver if rms < rmsTol
     real (kind=prec)   :: rmsTol
     ! exit if lambda < lambdaTol approx. 1e-4
     real (kind=prec)   :: lambdaTol
     ! set lambda_i = lambda_{i-1}/k when the inversion stalls
     real (kind=prec)   :: k
     ! initial value of lambda
     real (kind=prec)   :: lambda
     ! scale factor for the step length
     real (kind=prec)   :: alpha
     ! exponential factor for the Levy destribution
     real (kind=prec)   :: beta
     ! threshold for the stagnation test
     real (kind=prec)   :: epsilon
     ! lower bounds of parameters m
     real (kind=prec)   :: m_l
     ! upper bounds of parameters m
     real (kind=prec)   :: m_u
     ! model and data output file name
     character(80)              :: fname
  end type LGDiterControl_t

  type(LGDiterControl_t), private, save :: iterControl

Contains

!**********************************************************************
   subroutine set_LGDiterControl(iterControl)
   !  at present this is a simple routine to setup for
   !  iteration control; allowing parameters to be read from
   !   a file would be desireable.

   type(LGDiterControl_t), intent(inout)	:: iterControl

     ! maximum number of iterations in one call to iterative solver
     iterControl%maxIter = 1000
     ! convergence criteria: return from solver if rms < rmsTol
     iterControl%rmsTol  = 1.05
     ! stop updating lambda if lambda < lambdaTol approx. 1e-4
     iterControl%lambdaTol = 1.0e-4 ! makes no sense to use too small a value
     ! set lambda_i = lambda_{i-1}/k when the inversion stalls
     iterControl%k = 10.
     ! initial value of lambda
     iterControl%lambda = 1.
     ! scale factor for the step length
     iterControl%alpha = 0.01
     ! exponential factor for the Levy destribution
     iterControl%beta = 1.5
     ! threshold for the stagnation test
     iterControl%epsilon = 1.0e-12
     ! lower bounds of parameters m
     iterControl%m_l = 1.
     ! upper bounds of parameters m
     iterControl%m_u = 1000.
     ! model and data output file name
     iterControl%fname = 'Modular'

   end subroutine set_LGDiterControl


   ! **************************************************************************
   ! * read_LGDiterControl reads the inverse solver configuration from file

   subroutine read_LGDiterControl(iterControl,rFile,fileExists)

    type(LGDiterControl_t), intent(inout)  :: iterControl
    character(*), intent(in)                :: rFile
    logical, intent(out), optional          :: fileExists
    integer                                 :: ios
    logical                                 :: exists
    character(80)                           :: string

    ! Initialize inverse solver configuration

    call set_LGDiterControl(iterControl)

    inquire(FILE=rFile,EXIST=exists)
    if (present(fileExists)) then
       fileExists = exists
    end if

    if (.not. exists) then
       return
    else
       write(*,*) 'Reading inverse configuration from file ',trim(rFile)
    end if

    open (unit=ioInvCtrl,file=rFile,status='old',iostat=ios)

    if(ios/=0) then
       write(0,*) 'Error opening file: ', rFile
    end if

    ! This is the list of options specified in the startup file

    read (ioInvCtrl,'(a36,a80)') string,iterControl%fname
    if (output_level > 2) then
       write (*,*)
       write (*,'(a36,a80)') string,iterControl%fname
    end if
    iterControl%fname = adjustl(iterControl%fname)
    read (ioInvCtrl,'(a36,g15.7)') string,iterControl%lambda
    if (output_level > 2) then
       write (*,'(a36,g15.7)') string,iterControl%lambda
    end if
    read (ioInvCtrl,'(a36,g15.7)') string,iterControl%k
    if (output_level > 2) then
       write (*,'(a36,g15.7)') string,iterControl%k
    end if
    read (ioInvCtrl,'(a36,g15.7)') string,iterControl%alpha
    if (output_level > 2) then
       write (*,'(a36,g15.7)') string,iterControl%alpha
    end if
    read (ioInvCtrl,'(a36,g15.7)') string,iterControl%beta
    if (output_level > 2) then
       write (*,'(a36,g15.7)') string,iterControl%beta
    end if
    read (ioInvCtrl,'(a36,g15.7)') string,iterControl%epsilon
    if (output_level > 2) then
       write (*,'(a36,g15.7)') string,iterControl%epsilon
    end if
    read (ioInvCtrl,'(a36,g15.7)') string,iterControl%m_l
    if (output_level > 2) then
       write (*,'(a36,g15.7)') string,iterControl%m_l
    end if
    read (ioInvCtrl,'(a36,g15.7)') string,iterControl%m_u
    if (output_level > 2) then
       write (*,'(a36,g15.7)') string,iterControl%m_u
    end if
    read (ioInvCtrl,'(a36,g15.7)') string,iterControl%rmsTol
    if (output_level > 2) then
       write (*,'(a36,g15.7)') string,iterControl%rmsTol
    end if
    read (ioInvCtrl,'(a36,g15.7)') string,iterControl%lambdaTol
    if (output_level > 2) then
       write (*,'(a36,g15.7)') string,iterControl%lambdaTol
    end if
    read (ioInvCtrl,'(a36,i5)') string,iterControl%maxIter
    if (output_level > 2) then
       write (*,'(a36,i5)') string,iterControl%maxIter
       write (*,*)
    end if

    close(ioInvCtrl)

   end subroutine read_LGDiterControl

!**********************************************************************
   subroutine update_damping_parameter(lambda,mHat,F,grad)

   real(kind=prec), intent(inout)  :: lambda
   type(modelParam_t), intent(in)              :: mHat
   real(kind=prec), intent(inout)  :: F
   type(modelParam_t), intent(inout)             :: grad

   real(kind=prec) :: SS, mNorm, Nmodel
	 type(modelParam_t)          :: dSS

   ! compute the model norm
   mNorm = dotProd(mHat,mHat)
   Nmodel = countModelParam(mHat)

   ! (scaled) sum of squares = penalty functional - scaled model norm
   SS = F - (lambda * mNorm/Nmodel)

   ! initialize
   dSS = mHat

   ! subtract the model norm derivative from the gradient of the penalty functional
   call linComb(ONE,grad,MinusTWO*lambda/Nmodel,mHat,dSS)

   ! update the damping parameter lambda
   lambda = lambda/iterControl%k

   ! penalty functional = (scaled) sum of squares + scaled model norm
   F = SS + (lambda * mNorm/Nmodel)
   ! write(ioLog, *)'Roughness of model: ', mNorm/Nmodel

	 ! add the model norm derivative to the gradient of the penalty functional
   call linComb(ONE,dSS,TWO*lambda/Nmodel,mHat,grad)

   call deall_modelParam(dSS)

   end subroutine update_damping_parameter

!**********************************************************************
   subroutine LGDsolver(d,lambda,m0,m,fname)

   !  d is data; on output it contains the responses for the inverse model
   type(dataVectorMTX_t), intent(inout)         :: d
   !  lambda is regularization parameter
   real(kind=prec), intent(inout)               :: lambda
   !  m0 is prior model parameter
   type(modelParam_t), intent(in)               :: m0
   !  m is solution parameter ... on input m contains starting guess
   type(modelParam_t), intent(inout)            :: m
   !  fname is a string that specifies the control file
   character(*), intent(in), optional           :: fname


   !  local variables
   type(dataVectorMTX_t)                        :: dHat, dBest, res
   type(modelParam_t)                           :: mHat, m_minus_m0, mBest, mMean, mStdDev, mMeanOld, mStdDevOld
   type(modelParam_t)                           :: grad, g, h, gPrev, Normalized_g
   real(kind=prec)                              :: value, valueBest, rms, rmsBest, rmsBestPrev
   real(kind=prec)                              :: alpha, beta, epsilon
   real(kind=prec)                              :: gnorm, mNorm, mNormBest, Nmodel
   real(kind=prec)                              :: g_dot_g, g_dot_gPrev, gPrev_dot_gPrev
   real(kind=prec)                              :: LevyStep
   real(kind=prec)                              :: m_u, m_l, p
   real(kind=prec)                              :: coefficientPR
   integer                                      :: iter, iterBest, iter_no_update, ios
   integer                                      :: lambda_update_level
   logical                                      :: ok
   character(3)                                 :: iterChar
   character(100)                               :: mFile, mHatFile, gradFile
   character(100)                               :: dataFile, resFile, logFile
   type(solnVectorMTX_t)                        :: eAll


   if (present(fname)) then
      call read_LGDiterControl(iterControl,fname,ok)
      if (ok) then
         lambda = iterControl%lambda
      end if
   else
      call set_LGDiterControl(iterControl)
   end if

   ! initialize the output to log file
   logFile = trim(iterControl%fname)//'_LGD.log'
   open (unit=ioLog,file=logFile,status='unknown',position='append',iostat=ios)

   write(*,'(a33,es8.1)') 'The damping parameter lambda is ',lambda
   write(ioLog,'(a33,es8.1)') 'The damping parameter lambda is ',lambda

   ! Initialize parameters
   iter = 0
   lambda_update_level = 0
   alpha = iterControl%alpha
   beta = iterControl%beta
   epsilon = iterControl%epsilon

   ! Initialize boundarys and scale factors for individual parameters
   m_u = iterControl%m_u
   m_l = iterControl%m_l
   p = m_u - m_l

   ! starting model contains the rough deviations from the prior
   mHat = m



   !  compute the penalty functional and predicted data
   call func(lambda,d,m0,mHat,value,mNorm,dHat,eAll,rms)
   call printfLevy('START',lambda,value,mNorm,rms)
   call printfLevy('START',lambda,value,mNorm,rms,logFile)
   rmsBest = rms
   iter_no_update = 0

   write(iterChar,'(i3.3)') 0
   ! output (smoothed) initial model and responses for later reference
   call CmSqrtMult(mHat,m_minus_m0)
   call linComb(ONE,m_minus_m0,ONE,m0,m)
   if (output_level > 1) then
     mFile = trim(iterControl%fname)//'_LGD_'//iterChar//'.rho'
     call write_modelParam(m,trim(mFile))
   end if
   if (output_level > 2) then
     dataFile = trim(iterControl%fname)//'_LGD_'//iterChar//'.dat'
     call write_dataVectorMTX(dHat,trim(dataFile))
   end if

   ! compute gradient of the full penalty functional
   call gradient(lambda,d,m0,mHat,grad,dHat,eAll)
   if (output_level > 3) then
     gradFile = trim(iterControl%fname)//'_LGD_'//iterChar//'.grt'
     call write_modelParam(grad,trim(gradFile))
   end if

   ! check if gradient is ok
   gnorm = sqrt(dotProd(grad,grad))
   write(*,'(a42,es12.5)') '    GRAD: initial norm of the gradient is',gnorm
   write(ioLog,'(a42,es12.5)') '     GRAD: initial norm of the gradient is',gnorm
   if (gnorm < epsilon) then
      call errStop('Problem with your gradient computations: first gradient is zero')
   end if

	g = grad
   call linComb(R_ZERO, g, ONE/gnorm, g, Normalized_g)
   call linComb(R_ZERO, g, MinusONE, Normalized_g, h)


!---------------------------------loooooops!--------------------------------
   do
      !  test for convergence ...
      if((rmsBest .lt. iterControl%rmsTol) .or. (iter .ge. iterControl%maxIter)) then
         exit
      end if
      iter = iter + 1

	   ! save the values of the functional and the directional derivative
	   gPrev = g

      ! initialize LevyStep
      call LevyDistribution(beta, LevyStep)

      ! LGD
      call linComb(ONE, mHat, alpha * LevyStep * p, h, mHat)
      write(ioLog, '(a16, es12.6)') 'StepLength = ', alpha * LevyStep * p
      ! write(ioLog, *) 'LevyStep = ', LevyStep

      ! call UpperLowerBounds(m_u, m_l, mHat)
      call func(lambda,d,m0,mHat,value,mNorm,dHat,eAll,rms)
      call gradient(lambda,d,m0,mHat,grad,dHat,eAll)
      g = grad
      gnorm = sqrt(dotProd(g, g))
      call linComb(R_ZERO, g, ONE/gnorm, g, Normalized_g)
      call linComb(R_ZERO, h, MinusONE, Normalized_g, h)
	   write(*,'(a25,i5)') 'Completed LGD iteration ',iter
	   write(ioLog,'(a25,i5)') 'Completed LGD iteration ',iter
      write(ioLog, *) 'rms =', rms, 'rmsBest =', rmsBest
	   Nmodel = countModelParam(mHat)
	   mNorm = dotProd(mHat,mHat)/Nmodel
      write(ioLog, '(a30, es12.6)') 'Roughness of model: ', mNorm/Nmodel


      ! check if new model is better
      if (rms .le. rmsBest) then
         rmsBestPrev = rmsBest
         iter_no_update = 0
         iterBest = iter
         mBest = mHat
         dBest = dHat
         rmsBest = rms
         mNormBest = mNorm
         valueBest = value
         call printfLevy('updated',lambda,value,mNorm,rms)
         call printfLevy('updated',lambda,value,mNorm,rms, logFile)
      else
         iter_no_update = iter_no_update + 1
         call printfLevy('not updated',lambda,value,mNorm,rms)
         call printfLevy('not updated',lambda,value,mNorm,rms, logFile)
      end if

      ! write out the intermediate model solution and responses
      call CmSqrtMult(mHat,m_minus_m0)
      call linComb(ONE,m_minus_m0,ONE,m0,m)

      write(iterChar,'(i3.3)') iter
      if ((output_level > 1) .and. (iter_no_update .eq. 0))then
        mFile = trim(iterControl%fname)//'_LGD_'//iterChar//'.rho'
        call write_modelParam(m,trim(mFile))
      end if
      if ((output_level > 2) .and. (iter_no_update .eq. 0))then
        mHatFile = trim(iterControl%fname)//'_LGD_'//iterChar//'.prm'
        call write_modelParam(mHat,trim(mHatFile))
      end if
      if ((output_level > 2) .and. (iter_no_update .eq. 0)) then
        dataFile = trim(iterControl%fname)//'_LGD_'//iterChar//'.dat'
        call write_dataVectorMTX(dHat,trim(dataFile))
      end if
      ! compute residual for output: res = d-dHat; do not normalize by errors
      if ((output_level > 2) .and. (iter_no_update .eq. 0)) then
        res = d
        call linComb(ONE,d,MinusONE,dHat,res)
        resFile = trim(iterControl%fname)//'_LGD_'//iterChar//'.res'
        call write_dataVectorMTX(res,trim(resFile))
      end if

      ! compute mean model and StdDev model
      ! write(*,'(a25,i5)')'UPDATING MEAN AND STDDEV MODEL'
      if (iter .eq. 1) then
         mMeanOld = m
         mMean = m
         ! call linComb(ONE, m, MinusONE, m, mStdDevOld)
         ! call linComb(ONE, m, MinusONE, m, mStdDev)
      else if (iter .eq. 2) then
         call linComb(ONE/TWO, mMeanOld, ONE/TWO, m, mMean)
         call linComb(ONE/TWO, mMeanOld, MinusONE/TWO, m, mStdDev)
         call ElementWiseMult(mStdDev)
         call SqrtModel(mStdDev)
      else
         mMeanOld = mMean
         mStdDevOld = mStdDev
         call RunningStdDev(m, iter, mMeanOld, mMean, mStdDevOld, mStdDev)
      end if

      ! if ((abs(rmsBestPrev - rmsBest) .le. 2.0e-3) .or. (iter_no_update .eq. (iterControl%maxIter / 5))) then
      ! ! if (iter_no_update .eq. 10) then
      !    call update_damping_parameter(lambda,mHat,value,grad)
      !    lambda_update_level = lambda_update_level + 1
      !    if (lambda .le. iterControl%lambdaTol) then
      !       lambda = iterControl%lambdaTol
      !       write(*,'(a42)') 'lambda too small, stop updating lambda...'
      !       write(ioLog,'(a42)') 'lambda too small, stop updating lambda...'
      !    else
      !       write(*,'(a28, i5, a8, es12.6)') 'Lambda updated in iteration:', iter, ' lambda=',lambda
      !       write(ioLog,'(a28, i5, a8, es12.6)') 'Lambda updated in iteration:', iter, ' lambda=',lambda
      !    end if
      ! end if

      if (lambda_update_level .lt. 4) then
         if (iter_no_update .ge. ((iterControl%maxIter - iter)/(4 - lambda_update_level))) then
            lambda_update_level = lambda_update_level + 1
            call update_damping_parameter(lambda,mHat,value,grad)
            write(*,'(a28, i5, a8, es12.6)') 'Lambda updated in iteration:', iter, ' lambda=',lambda
            write(ioLog,'(a28, i5, a8, es12.6)') 'Lambda updated in iteration:', iter, ' lambda=',lambda
         end if
      else if (lambda_update_level .eq. 4) then
         lambda_update_level = lambda_update_level + 1
         write(*,'(a42)') 'lambda too small, stop updating lambda...'
         write(ioLog,'(a42)') 'lambda too small, stop updating lambda...'
      end if

      ! g_dot_g = dotProd(g,g)
      ! g_dot_gPrev = dotProd(g,gPrev)
      ! gPrev_dot_gPrev = dotProd(gPrev,gPrev)

      ! ! Polak-Ribiere variant
      ! coefficientPR = ( g_dot_g - g_dot_gPrev ) / gPrev_dot_gPrev

      ! if (gnorm .le. epsilon) then
      !     ! Stagnation test failed
      !     write(*,'(a45)') 'Stagnation test failed, updating the gradients orthogonally'
      !     write(ioLog,'(a45)') 'Stagnation test failed, updating the gradients orthogonally'
      !     call linComb(ONE, g, coefficientPR, h, h)
      ! end if

   end do


   ! write out the mean and StdDev model solution
   write(iterChar,'(i3.3)') iter
   if (output_level > 1)then
      mFile = trim(iterControl%fname)//'_LGD_MEAN.rho'
      call write_modelParam(mMean,trim(mFile))
   end if

   write(iterChar,'(i3.3)') iter
   if (output_level > 1)then
      mFile = trim(iterControl%fname)//'_LGD_StdDev.rho'
      call write_modelParam(mStdDev,trim(mFile))
   end if


   ! multiply by C^{1/2} and add m_0
   call CmSqrtMult(mBest,m_minus_m0)
   call linComb(ONE,m_minus_m0,ONE,m0,m)
   d = dBest
   write(*,'(a20,i5,a25,i5)') 'LGD iterations:',iter, 'Best iteration:', iterBest
   write(ioLog,'(a20,i5,a25,i5)') 'LGD iterations:',iter, 'Best iteration:', iterBest

   call printfLevy('Final Best',lambda,valueBest,mNormBest,rmsBest)
   call printfLevy('Final Best',lambda,valueBest,mNormBest,rmsBest, logFile)




   close(ioLog,iostat=ios)

   ! cleaning up
   call deall_dataVectorMTX(dHat)
   call deall_dataVectorMTX(res)
   call deall_dataVectorMTX(dBest)
   call deall_modelParam(mHat)
   call deall_modelParam(m_minus_m0)
   call deall_modelParam(mBest)
   call deall_modelParam(mMean)
   call deall_modelParam(mStdDev)
   call deall_modelParam(mMeanOld)
   call deall_modelParam(mStdDevOld)
   call deall_modelParam(grad)
   call deall_modelParam(g)
   call deall_modelParam(h)
   call deall_modelParam(gPrev)
   call deall_modelParam(Normalized_g)
   call deall_solnVectorMTX(eAll)

   end subroutine LGDsolver


  ! The step length of LGD is given by the Levy distribution, so we don't
  ! need code for line search any more. But relatively, we need a new subroutine
  ! to generate step lengths that follow the Levy distribution.
  !**********************************************************************
   subroutine gamma_lanczos(x, gamma_value)

      implicit none
      real (kind=prec), intent(in) :: x
      real (kind=prec), intent(out) :: gamma_value
      real (kind=prec), parameter :: g = 7.0_prec
      real (kind=prec), parameter :: p(0:8) = [ &
            0.99999999999980993_prec, 676.5203681218851_prec, -1259.1392167224028_prec, &
            771.32342877765313_prec, -176.61502916214059_prec, 12.507343278686905_prec, &
            -0.13857109526572012_prec, 9.9843695780195716e-6_prec, 1.5056327351493116e-7_prec]
      real (kind=prec) :: y, z, sum, t
      integer :: i, flag

      if (x .lt. 0.5_prec) then
         flag = 1
         y = ONE - x
      else
         flag = 0
         y = x
      end if

      z = y - ONE
      sum = p(R_ZERO)
      do i = 1, 8
         sum = sum + p(i) / (z + i)
      end do
      t = z + g + 0.5_prec
      gamma_value = sqrt(TWO * PI) * t**(z + 0.5_prec) * exp(-t) * sum

      if (flag .eq. 1) then
         gamma_value = PI / (sin(PI * x) * gamma_value)
      end if
   end subroutine gamma_lanczos

  !**********************************************************************
   subroutine NormalDistribution1(mean, std, x)
  ! Generate random number follows given normal distribution

      real(kind=prec), intent(in) :: mean, std
      real(kind=prec), intent(out) :: x
      real(kind=prec) :: u, v
      real(kind=prec) :: r, theta
      ! integer :: clock, n, i
      ! integer, allocatable :: seed(:)

      ! call system_clock(count=clock)
      ! call random_seed(size=n)
      ! allocate(seed(n))
      ! seed = clock + 42 * [(i - 1, i = 1, n)]
      ! call random_seed(put=seed)

      call random_seed()
      call random_number(u)
      call random_number(v)
      ! deallocate(seed)
      r = dsqrt(MinusTWO * log(u))
      theta = TWO * PI * v
      x = mean + std * r * sin(theta)

  end subroutine NormalDistribution1

  !**********************************************************************
   subroutine NormalDistribution2(mean, std, x)
  ! Generate random number follows given normal distribution

      real(kind=prec), intent(in) :: mean, std
      real(kind=prec), intent(out) :: x
      real(kind=prec) :: u, v
      real(kind=prec) :: r, theta
      integer :: clock, n, i
      integer, allocatable :: seed(:)

      call system_clock(count=clock)
      call random_seed(size=n)
      allocate(seed(n))
      seed = clock + 42 * [(i - 1, i = 1, n)]
      call random_seed(put=seed)

      ! call random_seed()
      call random_number(u)
      call random_number(v)
      deallocate(seed)
      r = dsqrt(MinusTWO * log(u))
      theta = TWO * PI * v
      x = mean + std * r * cos(theta)

  end subroutine NormalDistribution2

!***************************************************************************
   subroutine LevyDistribution(beta, LevyStep)
      ! Generate Levy Step

      real(kind=prec), intent(in) :: beta
      real(kind=prec), intent(out) :: LevyStep
      real(kind=prec) :: U, V
      real(kind=prec) :: sigma_u, sigma_v
      real(kind=prec) :: gamma_u1, gamma_u2

      call gamma_lanczos(ONE+beta, gamma_u1)
      call gamma_lanczos((ONE+beta)/TWO, gamma_u2)
      sigma_u = (gamma_u1*sin(pi*beta/TWO))/(gamma_u2*beta*(TWO**((beta-ONE)/TWO)))**(ONE/beta)
      sigma_v = ONE

      call NormalDistribution1(R_ZERO, sigma_u, U)
      write(ioLog, *) 'U = ', U
      call NormalDistribution2(R_ZERO, sigma_v, V)
      write(ioLog, *) 'V = ', V
      LevyStep = abs(U / (abs(V) ** (ONE / beta)))

  end subroutine LevyDistribution


!    subroutine LevyDistribution(beta, LevyStep)
!    implicit none
    
!    real (kind=prec), intent(out) :: LevyStep
!    real (kind=prec), intent(in) :: beta
!    real (kind=prec) :: U, V
!    real (kind=prec) :: sigma_u, sigma_v

!     call random_number(U)
!     call random_number(V)

!     sigma_u = (gamma(1+beta)*sin(pi*beta/2))/(gamma((1+beta)/2)*beta*(2**((beta-1)/2)))**(1/beta)
!     sigma_v = 1

!     LevyStep = abs(U/sigma_u)**(1/beta)
    
! contains

!     function gamma(z)
!         real (kind=prec), intent(in) :: z
!         real (kind=prec) :: gamma
!         gamma = exp(gammln(z))
!     contains
!         function gammln(z)
!             real (kind=prec), intent(in) :: z
!             real (kind=prec) :: gammln
!             gammln = log(sqrt(2 * PI)) + (z - 0.5)*log(z) - z
!         end function gammln
!     end function gamma

! end subroutine LevyDistribution

   subroutine RunningMean(m_iter, iter, mMeanOld, mMean)
      ! compute mean model dynamically

      type(modelParam_t), intent(in) :: m_iter
      type(modelParam_t), intent(in) :: mMeanOld
      type(modelParam_t), intent(inout) :: mMean
      integer, intent(in) :: iter

      ! if (iter .eq. 1) then
      !    mMeanOld = m_iter
      !    mMean = m_iter
      !    ! write(*,*),m_iter%Ny, m_iter%Nx, m_iter%NzEarth
      !    write(*,'(a25,i5)')'MeanIterOne'
      ! else
      !    write(*,'(a25,i5)')'Mean update'
      !    ! mMean = (mMeanOld * (iter - 1) + m_iter) / iter
      !    call linComb((iter - ONE)/iter, mMeanOld, (ONE/iter), m_iter, mMean)
      !    write(*,'(a25,i5)')'Mean update end'
      ! end if

      ! write(*,'(a25,i5)')'Mean update'
      ! mMean = (mMeanOld * (iter - 1) + m_iter) / iter
      call linComb((iter - ONE)/iter, mMeanOld, (ONE/iter), m_iter, mMean)
      ! write(*,'(a25,i5)')'Mean update end'

  end subroutine RunningMean

!***************************************************************************
  subroutine RunningStdDev(m_iter, iter, mMeanOld, mMean, mStdDevOld, mStdDev)
      ! compute popular StdDev model dynamically

      type(modelParam_t), intent(in) :: m_iter
      type(modelParam_t), intent(in) :: mMeanOld
      type(modelParam_t), intent(inout) :: mMean
      type(modelParam_t), intent(inout) :: mStdDevOld
      type(modelParam_t), intent(inout) :: mStdDev
      integer, intent(in) :: iter
      character(80) :: paramType

      ! if (iter .eq. 1) then
      !    ! after the first iteration, StdDev = 0
      !    call linComb(ONE, m_iter, MinusONE, m_iter, mStdDev)
      !    write(*,'(a25,i5)')'StdDevIterOne'
      ! else
      !    ! get new mean model
      !    write(*,'(a25,i5)')'StdDev update'
      !    call RunningMean(m_iter, iter, mMeanOld, mMean)
      !    write(*,'(a25,i5)')'RunningMean'
      !    call linComb(ONE, mMean, MinusONE, m_iter, mStdDev)
      !    write(*,'(a25,i5)')'linComb 1'
      !    call ElementWiseMult(mStdDev)
      !    write(*,'(a25,i5)')'ElementWiseMult 1'
      !    call ElementWiseMult(mStdDevOld)
      !    write(*,'(a25,i5)')'ElementWiseMult 2'
      !    call linComb(ONE/(iter - ONE), mStdDev, (iter - ONE)/iter, mStdDevOld, mStdDev)
      !    write(*,'(a25,i5)')'linComb 2'
      !    call SqrtModel(mStdDev)
      !    write(*,'(a25,i5)')'StdDev update end'
      ! end if

      ! get new mean model
      ! write(*,'(a25,i5)')'StdDev update'
      call RunningMean(m_iter, iter, mMeanOld, mMean)
      ! write(*,'(a25,i5)')'RunningMean'
      call linComb(ONE, mMean, MinusONE, m_iter, mStdDev)
      ! write(*,'(a25,i5)')'linComb 1'
      call ElementWiseMult(mStdDev)
      ! write(*,'(a25,i5)')'ElementWiseMult 1'
      call ElementWiseMult(mStdDevOld)
      ! write(*,'(a25,i5)')'ElementWiseMult 2'
      call linComb(ONE/(iter - ONE), mStdDev, (iter - ONE)/iter, mStdDevOld, mStdDev)
      write(*,'(a25,i5)')'linComb 2'
      call SqrtModel(mStdDev)


      call getType_modelParam(mStdDev,paramType)
      write(*,*)'ParamType of StdDev model is ', paramType
      ! write(*,'(a25,i5)')'StdDev update end'

  end subroutine RunningStdDev

end module LGD





