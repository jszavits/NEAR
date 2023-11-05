!==============================================================================
!
!   MODEL
!   -----
!   Totally asymmetric simple exclusion process with open boundary conditions,
!   elongated particles and site-dependent hopping rates.
!
!   PROGRAM
!   -------
!   The program computes elongation rates relative to the initiation rate from 
!   the Ribo-seq A-site density profile using NLOPT package and the third order 
!   of the power series method. 
!
!==============================================================================

module arrays
  real(kind=8),allocatable :: f1(:),f2(:),f3(:),f22(:,:),f21(:),&
                  	      & f33a(:,:),f33b(:,:),f33L(:,:),f32(:,:),f31(:)
end module

module objfun
  real(kind=8),allocatable :: rhoexp(:),omega(:)
  integer(kind=4),allocatable :: included(:),positions(:)
  integer(kind=4) :: L,ll
  integer(kind=8) :: Niter
end module

program dTASEP
  
  use mt19937
  use arrays
  use objfun
  implicit none

  interface
    subroutine mc(L,ll,omega,tau,iter,t,rho,current)
      integer(kind=4),intent(in) :: L,ll
      real(kind=8),intent(in) :: omega(2:L)
      integer(kind=4),intent(inout) :: tau(2:L)
      integer(kind=8),intent(in) :: iter     
      real(kind=8),intent(inout) :: t
      real(kind=8),intent(out) :: rho(2:L),current(1:L)
    end subroutine mc

    subroutine psm3(n,s,x,y)
      integer(kind=4),intent(in) :: n,s
      real(kind=8),intent(in) :: x(2:n)
      real(kind=8),intent(out) :: y(2:n)
    end subroutine psm3
  end interface

  ! nlopt constants
  include 'nlopt.f'

  ! program parameters
  integer(kind=4) :: err11,err22,err33,err44,err55,stat22,stat33,i,j,k,&
		     & filenamelen,Nequations,filenamelen2,ires
  integer(kind=8) :: opt,iter
  integer(kind=4),allocatable :: tau(:),mfok(:),psmok(:)
  integer(kind=8),allocatable :: reads(:)
  character(len=80) :: inputfile
  character(len=30) :: str
  character(len=10) :: genename
  character(len=3) :: codoni
  character(len=3),allocatable :: codon(:)
  logical :: eof,eof2,found,under1
  real(kind=8),allocatable :: x(:),lb(:),ub(:),grad(:),rhomc(:),current(:),&
                              & rhomc0(:),current0(:),omega0(:),rhopsm(:)
  real(kind=8) :: t1,t2,ftol,mftol,psmtol,readsi,h1,hi,density,avreads,minf,time,&
                  & tden1,tden2,rel,fracerr,maxtime,ttime,val0,val
  external fnc3

!======================================================================
! INPUT
!======================================================================

  call cpu_time(t1)

  open(unit=11,file='NEAR.dat',status='OLD',iostat=err11)
  read(11,*) ll			! ribosome footprint length
  read(11,*) mftol		! tolerance for the relative error between rhoexp and rhomc0
  read(11,*) psmtol		! tolerance for the relative error between rhomc0 and rhopsm
  read(11,*) ftol		! stop iterating when the relative change in the objective 
						! function between two iterations is ftol
  read(11,*) maxtime    ! stop iterating after reaching maxtime seconds
  read(11,*) iter		! number of iterations for MC (Gillespie algorithm)
  read(11,*) inputfile	! list with genes to optimise	
  close (11)

!======================================================================
! MAIN PROGRAM
!======================================================================

  ! initial seed for random number generator
  call init_genrand(4357)

  ! list of genes to analyse
  filenamelen=len(trim(adjustl(inputfile)))
  open(22,file=inputfile(1:filenamelen),status='OLD',iostat=err22)

  ! main loop over a list of genes to analyse
  eof=.FALSE.
  do while (eof.EQV..FALSE.)
    read(22,*,iostat=stat22) genename,density
    if (stat22.EQ.0) then

	  ! begin measuring time
      call cpu_time(t1)

      ! write to log file
      filenamelen2=len(trim(adjustl(genename)))
      open(44,file='output/'//genename(1:filenamelen2)//'-NEAR.log',&
        & status='REPLACE',iostat=err44)
      write(44,'(A)') 'INFO: Program started.'
      write(44,'(A)') 'INFO: Reading A-site_'//genename(1:filenamelen2)//'.dat.'
      close(44)

      ! read file with reads      
      open(33,file='ribo-seq/A-site_'//genename(1:filenamelen2)//'.dat',&
        & status='OLD',iostat=err33)
		
	  if (err33.eq.0) then
        k=0
        eof2=.FALSE.
        do while (eof2.EQV..FALSE.)
          read(33,*,iostat=stat33) codoni,readsi
          if (stat33.EQ.0) then		
            k=k+1
          else
      	    eof2=.TRUE.
          endif
        enddo
        close(33)

        ! number of codons
        L=k+1

        ! array to store absolute ribosome density
        allocate(rhoexp(2:L))

        ! read file with reads again
        allocate(codon(2:L),reads(2:L))
        open(33,file='ribo-seq/A-site_'//genename(1:filenamelen2)//'.dat',&
		  & status='OLD',iostat=err33)
        j=0
        avreads=0.0d0
        do i=2,L
          read(33,*) codon(i),reads(i)
          ! normalization 
          if (reads(i).gt.0) then
            rhoexp(i)=dble(reads(i))*density
            avreads=avreads+dble(reads(i))
          else
            j=j+1
            rhoexp(i)=density
            avreads=avreads+1.0d0
          endif
        enddo
        avreads=avreads/dble(L-1)
        rhoexp=rhoexp/avreads
        close(33)

        ! check if all local densities are <1
        under1=.true.
        i=2
        do while ((under1.eqv..true.).and.(i.LE.L))
          if (rhoexp(i).gt.1.0d0) then
            under1=.false.
          endif
          i=i+1
        enddo

        ! proceed if all local densities are <1
        if (under1.eqv..true.) then
        
          ! write to log file
          open(44,file='output/'//genename(1:filenamelen2)//'-NEAR.log',&
            & status='OLD',position='APPEND',iostat=err44)
          write(str,'(I4)') L
          write(44,'(A)') 'INFO: Number of codons including START = ' // trim(adjustl(str))
          write(str,'(I4)') j
          write(44,'(A)') 'INFO: Number of codons with zero reads = ' // trim(adjustl(str))
          write(44,'(A)') 'INFO: Getting elongation-to-initation rates from the mean-field solution.' 

          ! initial estimate of the rates from the mean-field solution
          allocate(omega0(2:L),omega(2:L))
          h1=1.0d0-sum(rhoexp(2:ll+1))
          do i=2,L-ll
            hi=rhoexp(i)*(1.0d0-sum(rhoexp(i+1:i+ll)))/(1.0d0-sum(rhoexp(i+1:i+ll))+rhoexp(i+ll))
            omega0(i)=h1/hi
            if (omega0(i).lt.0.0d0) then
              write(44,'(A60,I7)') 'WARNING: Negative MF elongation rate encountered at codon = ',i
              omega0(i)=1.0d0/rhoexp(i)
            endif
          enddo
          do i=L-ll+1,L
            hi=rhoexp(i)
            omega0(i)=h1/hi
            if (omega0(i).lt.0.0d0) then
              write(44,'(A60,I7)') 'WARNING: Negative MF elongation rate encountered at codon = ',i
              omega0(i)=1.0d0/rhoexp(i)
            endif
          enddo
          close(44)

          ! write to log file
          open(44,file='output/'//genename(1:filenamelen2)//'-NEAR.log',&
            & status='OLD',position='APPEND',iostat=err44)
          write(44,'(A)') 'INFO: Computing local densities using MC simulation &
		    & with the mean-field rates.'
          close(44)

          ! run Monte Carlo simulation with MF rates (to reach the steady state)
          allocate(tau(2:L))
          allocate(rhomc0(2:L),current0(1:L))
          ttime=0.0d0
          time=0.0d0
          tau=0
          rhomc0=0.0d0
          current0=0.0d0
          found=.FALSE.
          do while (found.EQV..FALSE.)
            tden1=sum(rhomc0)/dble(L-1)
            call mc(L,ll,omega0,tau,100_8,time,rhomc0,current0)
            ttime=ttime+time
            tden2=sum(rhomc0)/dble(L-1)
            if (tden2.gt.epsilon(1.0d0)) then
              rel=abs(tden2-tden1)/tden2
              if (rel.lt.1.d-3) then
                found=.TRUE.
              endif
            endif
          enddo

          ! run Monte Carlo simulation with MF rates (in the steady state)
		  ! and compute density profile
          time=0.0d0
          rhomc0=0.0d0
          current0=0.0d0
          call mc(L,ll,omega0,tau,iter,time,rhomc0,current0)

          ! compute the initial value of the objective function
          val0=0.0d0
          do i=2,L
            val0=val0+(rhomc0(i)-rhoexp(i))*(rhomc0(i)-rhoexp(i))
          enddo

          ! write to log file
          open(44,file='output/'//genename(1:filenamelen2)//'-NEAR.log',&
            & status='OLD',position='APPEND',iostat=err44)
          write(44,'(A42,ES23.15E3)') 'INFO: Sum of squares (mean-field rates) = ',val0
          close(44)

          ! find codons for which the relative error between rhomc0 and rhoexp is < mftol
          allocate(mfok(2:L))
          mfok=0
          do i=2,L
            fracerr=abs(rhoexp(i)-rhomc0(i))/rhoexp(i) 
            if (fracerr.lt.mftol) then
              mfok(i)=1  ! mf is a good approximation
            endif
          enddo

          ! auxillary arrays for the power series method
          allocate(f22(2:L-ll,2+ll:L),f21(2:L))
          allocate(f33a(2+ll:L-ll,2*ll+2:L),f33b(2+ll:L-ll,2*ll+2:L),f33L(2:L-2*ll,2+ll:L-ll))
          allocate(f32(2:L-ll,2+ll:L),f31(2:L))
          allocate(f1(2:L),f2(2:L),f3(2:L))
        
          ! find codons for which the relative error between rhopsm and rhomc0 is < psmtol
          allocate(rhopsm(2:L))
          call psm3(L,ll,omega0,rhopsm)
          allocate(psmok(2:L))
          do i=2,L
            fracerr=abs(rhomc0(i)-rhopsm(i))/rhomc0(i) 
            if (fracerr.lt.psmtol) then
              psmok(i)=1  ! power series method is a good approximation
            endif 
          enddo

          ! write to log file
          open(44,file='output/'//genename(1:filenamelen2)//'-NEAR.log',&
            & status='OLD',position='APPEND',iostat=err44)
 
          allocate(included(2:L))
          j=0
          k=0
          do i=2,L
            if (mfok(i).eq.1) then
              included(i)=0
            else
              ! count codons with with mfok=0
              j=j+1
              if (psmok(i).eq.1) then
                included(i)=1  
              else
                included(i)=-1
                k=k+1
                write(44,'(A71,I7)') 'WARNING: Power series method gives error larger than &
			      & PSMTOL at codon = ',i
              endif
            endif
          enddo
          if (k.gt.0) then
            write(44,'(A47,I7)') 'INFO: Power series method may not be realiable.'
          endif
          close(44)

          ! number of equations for local optimisation
          Nequations=j

          ! mfok and psmok are not needed anymore
          deallocate(mfok)
          deallocate(psmok)

          ! nonlinear optimisation
          if (Nequations.gt.0) then

            ! select codons for nonlinear optimisation
            allocate(positions(1:Nequations))
            j=0
            do i=2,L
              if ((included(i).eq.1).or.(included(i).eq.-1)) then
                j=j+1
                positions(j)=i
              endif
            enddo

            ! write to log file
            open(44,file='output/'//genename(1:filenamelen2)//'-NEAR.log',&
              & status='OLD',position='APPEND',iostat=err44)
            write(44,'(A52,I7)') 'INFO: Number of codons for nonlinear optimisation = ', Nequations
            write(44,'(A)') 'INFO: Started derivative-free nonlinear least-squares &
			  & optimization (BOBYQA algoritm).'
            close(44)

            allocate(x(1:Nequations))	      

            ! select algorithm
            opt=0
            call nlo_create(opt,NLOPT_LN_BOBYQA,Nequations)

            ! box constraints
            allocate(lb(1:Nequations),ub(1:Nequations))
            lb=1.d-3
            ub=1.d6
            call nlo_set_lower_bounds(ires,opt,lb)
            call nlo_set_upper_bounds(ires,opt,ub)

            ! initial rates
            omega=omega0
            do j=1,Nequations
              i=positions(j)
              if (omega(i).LE.ub(j)) then
                x(j)=omega(i)
              else
                x(j)=ub(j)-epsilon(1.0d0)
              endif
            enddo

            ! objective function
            call nlo_set_min_objective(ires,opt,fnc3,0)

            ! stopping criteria
            call nlo_set_ftol_rel(ires,opt,ftol)
            call nlo_set_maxtime(ires,opt,maxtime)

            ! optimisation step
            Niter=0
            call nlo_optimize(ires,opt,x,minf)

            deallocate(x)
            deallocate(lb,ub)
	        deallocate(positions)

            if (ires.lt.0) then  ! optimisation unsuccessful
              ! write to log file
              open(44,file='output/'//genename(1:filenamelen2)//'-NEAR.log',&
                & status='OLD',position='APPEND',iostat=err44)
              write(44,'(A29,I7)') 'INFO: Number of iterations = ', Niter
              write(44,'(A37,I4)') 'INFO: NLopt FAILED with IRES value = ',ires
			  if (ires.eq.-1) then
			    write(44,'(A49)') 'INFO: NLopt stopped because of a generic failure.'
			  elseif (ires.eq.-2) then
			    write(44,'(A49)') 'INFO: NLopt stopped because of invalid arguments.'
			  elseif (ires.eq.-3) then
			    write(44,'(A30)') 'INFO: NLopt ran out of memory.'
			  else
			    write(44,'(A69)') 'INFO: NLopt stopped because rounding errors limited &
			      & further progress.' 
			  endif
              close(44)

            else  ! optimisation successful
              ! write to log file
              open(44,file='output/'//genename(1:filenamelen2)//'-NEAR.log',&
                & status='OLD',position='APPEND',iostat=err44)
              write(44,'(A29,I7)') 'INFO: Number of iterations = ', Niter
              write(44,'(A45,I4)') 'INFO: NLopt was SUCCESSFUL with IRES value = ',ires
              if (ires.eq.3) then
                write(44,'(A50)') 'INFO: NLopt finished because FTOL_REL was reached.'
              endif
              if (ires.eq.6) then
                write(44,'(A49)') 'INFO: NLopt finished because MAXTIME was reached.'
              endif
              write(44,'(A)') 'INFO: Computing local densities using MC simulation with the optimised rates.'
              close(44)

              ! run Monte Carlo simulation with MF rates (to reach the steady state)
              allocate(rhomc(2:L),current(1:L))
              ttime=0.0d0
              time=0.0d0
              tau=0
              rhomc=0.0d0
              current=0.0d0
              found=.FALSE.
              do while (found.EQV..FALSE.)
                tden1=sum(rhomc)/dble(L-1)
                call mc(L,ll,omega,tau,100_8,time,rhomc,current)
                ttime=ttime+time
                tden2=sum(rhomc)/dble(L-1)
                if (tden2.gt.epsilon(1.0d0)) then
        	      rel=abs(tden2-tden1)/tden2
        	      if (rel.lt.1.d-3) then
        	        found=.TRUE.
        	      endif
                endif
              enddo

              ! run Monte Carlo simulation with MF rates (in the steady state) and compute density profile
              time=0.0d0
              rhomc=0.0d0
              current=0.0d0
              call mc(L,ll,omega,tau,iter,time,rhomc,current)

              ! compute optimised value of the objective function
              val=0.0d0
              do i=2,L
                val=val+(rhomc(i)-rhoexp(i))*(rhomc(i)-rhoexp(i))
              enddo

              ! write to log file
              open(44,file='output/'//genename(1:filenamelen2)//'-NEAR.log',&
                & status='OLD',position='APPEND',iostat=err44)
              write(44,'(A39,ES23.15E3)') 'INFO: Sum of squares (optimised rates) = ',val
              close(44)

              ! write output
              open(55,file='output/'//genename(1:filenamelen2)//'-NEAR-results.dat',status='REPLACE',iostat=err55)
              do i=2,L
                write(55,'(A3,A1,I1,A1,ES23.15E3,A1,ES23.15E3,A1,ES23.15E3,A1,ES23.15E3,A1,ES23.15E3,A1,ES23.15E3,A1,& 
                  & ES23.15E3,A1,ES23.15E3,A1,ES23.15E3,A1,ES23.15E3,A1,ES23.15E3)') codon(i),' ',included(i),' ',&
                  & rhoexp(i),' ',rhomc0(i),' ',rhomc(i),' ',f1(i)+f2(i)+f3(i),' ',omega0(i),' ',omega(i),' ',&
                  & current0(i),' ',current(i),' ',f1(i),' ',f2(i),' ',f3(i)
              enddo
              close(55)
              deallocate(rhomc,current)
            endif
            call nlo_destroy(opt)
          else

            ! write to log file
            open(44,file='output/'//genename(1:filenamelen2)//'-NEAR.log',&
              & status='OLD',position='APPEND',iostat=err44)
            write(44,'(A)') 'INFO: Nothing to optimise. Rates are estimated from the mean-field approximation.'
            close(44)

            ! write output
            open(55,file='output/'//genename(1:filenamelen2)//'-NEAR-results.dat',status='REPLACE',iostat=err55)
            do i=2,L
              write(55,'(A3,A1,I1,A1,ES23.15E3,A1,ES23.15E3,A1,ES23.15E3,A1,ES23.15E3,A1,ES23.15E3,A1,ES23.15E3,A1,& 
                & ES23.15E3,A1,ES23.15E3,A1,ES23.15E3,A1,ES23.15E3,A1,ES23.15E3)') codon(i),' ',included(i),' ',&
                & rhoexp(i),' ',rhomc0(i),' ',rhomc0(i),' ',f1(i)+f2(i)+f3(i),' ',omega0(i),' ',omega0(i),' ',&
                & current0(i),' ',current0(i),' ',f1(i),' ',f2(i),' ',f3(i)
            enddo
            close(55)
          endif
 
		  ! stop measuring time
          call cpu_time(t2)

          ! write to log file
          open(44,file='output/'//genename(1:filenamelen2)//'-NEAR.log',&
            & status='OLD',position='APPEND',iostat=err44)
          write(44,'(A39,F9.3,A9)') 'INFO: Program finished successfully in ',t2-t1,' seconds.'
          close(44)

          deallocate(included) 
          deallocate(tau)
          deallocate(rhomc0,current0)
          deallocate(omega0,omega)
          deallocate(rhopsm)
          deallocate(f22,f21)
          deallocate(f33a,f33b,f33L)
          deallocate(f32,f31)
          deallocate(f1,f2,f3)
        else

          ! write to log file
          open(44,file='output/'//genename(1:filenamelen2)//'-NEAR.log',&
	        & status='OLD',position='APPEND',iostat=err44)
	      write(44,'(A)') 'ERROR: Unable to normalise Ribo-seq. Gene skipped.'
          close(44)
        endif

        deallocate(rhoexp)
        deallocate(codon,reads)
	  else
	    ! write to log file
          open(44,file='output/'//genename(1:filenamelen2)//'-NEAR.log',&
	        & status='OLD',position='APPEND',iostat=err44)
	      write(44,'(A)') 'ERROR: Unable to find a file with A-site reads. Gene skipped.'
          close(44)
	  endif
    else
      ! end of file
      eof=.TRUE.
    endif
  enddo
  close(22)

end program dTASEP

!======================================================================
! SUBROUTINES
!======================================================================

! computes local densities using power series method (third-order)
subroutine psm3(n,s,x,y)

  use arrays
  implicit none

  integer(kind=4),intent(in) :: n,s
  real(kind=8),intent(in) :: x(2:n)
  real(kind=8),intent(out) :: y(2:n)
  
  ! local variables

  integer(kind=4) :: i,j,k
  real (kind=8) :: b1,b2,e0

  ! FIRST ORDER
  !=============

  ! one particle
  b1=0.0d0
  do i=2,n
    f1(i)=1.0d0/x(i)
    b1=b1+f1(i)
  enddo

  ! SECOND-ORDER 
  !==============

  f2=0.0d0
  b2=0.0d0
  
  ! two-particle configurations
  f22=0.0d0
  do i=2,n-s
    do j=i+s,n
      e0=x(j)
      if (i.eq.2) then
        f22(i,j)=f1(j)
      else
        f22(i,j)=x(i-1)*f22(i-1,j)
      endif
      if (j-i.gt.s) then
        e0=e0+x(i)
        f22(i,j)=f22(i,j)+x(j-1)*f22(i,j-1)
      endif
      f22(i,j)=f22(i,j)/e0
      f2(i)=f2(i)+f22(i,j)
      f2(j)=f2(j)+f22(i,j)
      b2=b2+f22(i,j)
    enddo
  enddo

  ! one-particle configurations
  f21=0.0d0
  do i=2,n
    e0=x(i)
    if (i.eq.2) then
      f21(i)=x(n)*f22(i,n)
    elseif ((i.gt.2).and.(i.le.s+1)) then
      f21(i)=x(i-1)*f21(i-1)+x(n)*f22(i,n)
    elseif ((i.gt.s+1).and.(i.le.n-s)) then
      f21(i)=x(i-1)*f21(i-1)+x(n)*f22(i,n)-f1(i)
    else
      f21(i)=x(i-1)*f21(i-1)-f1(i)
    endif
    f21(i)=f21(i)/e0
    f2(i)=f2(i)+f21(i)-b1*f1(i)
    b2=b2+f21(i)
  enddo

  ! THIRD-ORDER 
  !=============

  f3=0.0d0

  ! three-particle configurations
  f33a=0.0d0
  f33b=0.0d0
  f33L=0.0d0
  do i=2,n-2*s
    do j=i+s,n-s
      do k=j+s,n
        e0=x(k)
        if (i.EQ.2) then
          f33b(j,k)=f22(j,k)
        else
          f33b(j,k)=x(i-1)*f33a(j,k)
        endif
        if (j-i.gt.s) then
          e0=e0+x(i)
          f33b(j,k)=f33b(j,k)+x(j-1)*f33b(j-1,k)
        endif
        if (k-j.gt.s) then
          e0=e0+x(j)
          f33b(j,k)=f33b(j,k)+x(k-1)*f33b(j,k-1)
        endif
        f33b(j,k)=f33b(j,k)/e0
        f3(i)=f3(i)+f33b(j,k)
        f3(j)=f3(j)+f33b(j,k)
        f3(k)=f3(k)+f33b(j,k)
      enddo
      f33L(i,j)=f33b(j,n)
    enddo
    f33a(i+s:n-s,i+2*s:n)=f33b(i+s:n-s,i+2*s:n)
  enddo

  ! two-particle configurations
  f32=0.0d0
  do i=2,n-s
    do j=i+s,n
      e0=x(j)
      if (i.eq.2) then
        f32(i,j)=f21(j)
      elseif ((i.gt.2).and.(i.le.s+1)) then
        f32(i,j)=x(i-1)*f32(i-1,j)
      else
        f32(i,j)=x(i-1)*f32(i-1,j)-f22(i,j)
      endif
      if (j-i.gt.s) then
        e0=e0+x(i)
        f32(i,j)=f32(i,j)+x(j-1)*f32(i,j-1)
      endif
      if (j.le.n-s) then
        f32(i,j)=f32(i,j)+x(n)*f33L(i,j)
      endif
      f32(i,j)=f32(i,j)/e0
      f3(i)=f3(i)+f32(i,j)
      f3(j)=f3(j)+f32(i,j)   
    enddo
  enddo

  ! one-particle configurations
  f31=0.0d0
  do i=2,n
    e0=x(i)
    if (i.eq.2) then
      f31(i)=x(n)*f32(i,n)
    elseif ((i.gt.2).and.(i.le.s+1)) then
      f31(i)=x(i-1)*f31(i-1)+x(n)*f32(i,n)
    elseif ((i.gt.s+1).and.(i.le.n-s)) then
      f31(i)=x(i-1)*f31(i-1)+x(n)*f32(i,n)-f21(i)
    else
      f31(i)=x(i-1)*f31(i-1)-f21(i)
    endif
    f31(i)=f31(i)/e0
    f3(i)=f3(i)+f31(i)-b1*f2(i)-b2*f1(i)
  enddo

  do i=2,n
    y(i)=f1(i)+f2(i)+f3(i)
  enddo
  return
end subroutine psm3

! computes the value of objective function
subroutine fnc3(val,n,x,grad,need_gradient,f_data)

  use arrays
  use objfun
  implicit none

  real(kind=8) :: val,x(n),grad(n)
  integer(kind=4) :: n,need_gradient,f_data

  ! local variables

  integer(kind=4) :: i,j,k
  real (kind=8) :: b1,b2,e0,rhoi

  ! FIRST ORDER
  !=============

  ! one particle
  b1=0.0d0
  j=0
  do i=2,L
    if ((included(i).eq.1).or.(included(i).eq.-1)) then
      j=j+1
      omega(i)=x(j)
    endif
    f1(i)=1.0d0/omega(i)
    b1=b1+f1(i)
  enddo

  ! SECOND-ORDER 
  !==============

  f2=0.0d0
  b2=0.0d0
  
  ! two-particle configurations
  f22=0.0d0
  do i=2,L-ll
    do j=i+ll,L
      e0=omega(j)
      if (i.eq.2) then
        f22(i,j)=f1(j)
      else
        f22(i,j)=omega(i-1)*f22(i-1,j)
      endif
      if (j-i.gt.ll) then
        e0=e0+omega(i)
        f22(i,j)=f22(i,j)+omega(j-1)*f22(i,j-1)
      endif
      f22(i,j)=f22(i,j)/e0
      f2(i)=f2(i)+f22(i,j)
      f2(j)=f2(j)+f22(i,j)
      b2=b2+f22(i,j)
    enddo
  enddo

  ! one-particle configurations
  f21=0.0d0
  do i=2,L
    e0=omega(i)
    if (i.eq.2) then
      f21(i)=omega(L)*f22(i,L)
    elseif ((i.gt.2).and.(i.le.ll+1)) then
      f21(i)=omega(i-1)*f21(i-1)+omega(L)*f22(i,L)
    elseif ((i.gt.ll+1).and.(i.le.L-ll)) then
      f21(i)=omega(i-1)*f21(i-1)+omega(L)*f22(i,L)-f1(i)
    else
      f21(i)=omega(i-1)*f21(i-1)-f1(i)
    endif
    f21(i)=f21(i)/e0
    f2(i)=f2(i)+f21(i)-b1*f1(i)
    b2=b2+f21(i)
  enddo

  ! THIRD-ORDER 
  !=============

  f3=0.0d0

  ! three-particle configurations
  f33a=0.0d0
  f33b=0.0d0
  f33L=0.0d0
  do i=2,L-2*ll
    do j=i+ll,L-ll
      do k=j+ll,L
        e0=omega(k)
        if (i.EQ.2) then
          f33b(j,k)=f22(j,k)
        else
          f33b(j,k)=omega(i-1)*f33a(j,k)
        endif
        if (j-i.gt.ll) then
          e0=e0+omega(i)
          f33b(j,k)=f33b(j,k)+omega(j-1)*f33b(j-1,k)
        endif
        if (k-j.gt.ll) then
          e0=e0+omega(j)
          f33b(j,k)=f33b(j,k)+omega(k-1)*f33b(j,k-1)
        endif
        f33b(j,k)=f33b(j,k)/e0
        f3(i)=f3(i)+f33b(j,k)
        f3(j)=f3(j)+f33b(j,k)
        f3(k)=f3(k)+f33b(j,k)
      enddo
      f33L(i,j)=f33b(j,L)
    enddo
    f33a(i+ll:L-ll,i+2*ll:L)=f33b(i+ll:L-ll,i+2*ll:L)
  enddo

  ! two-particle configurations
  f32=0.0d0
  do i=2,L-ll
    do j=i+ll,L
      e0=omega(j)
      if (i.eq.2) then
        f32(i,j)=f21(j)
      elseif ((i.gt.2).and.(i.le.ll+1)) then
        f32(i,j)=omega(i-1)*f32(i-1,j)
      else
        f32(i,j)=omega(i-1)*f32(i-1,j)-f22(i,j)
      endif
      if (j-i.gt.ll) then
        e0=e0+omega(i)
        f32(i,j)=f32(i,j)+omega(j-1)*f32(i,j-1)
      endif
      if (j.le.L-ll) then
        f32(i,j)=f32(i,j)+omega(L)*f33L(i,j)
      endif
      f32(i,j)=f32(i,j)/e0
      f3(i)=f3(i)+f32(i,j)
      f3(j)=f3(j)+f32(i,j)   
    enddo
  enddo

  ! one-particle configurations
  f31=0.0d0
  do i=2,L
    e0=omega(i)
    if (i.eq.2) then
      f31(i)=omega(L)*f32(i,L)
    elseif ((i.gt.2).and.(i.le.ll+1)) then
      f31(i)=omega(i-1)*f31(i-1)+omega(L)*f32(i,L)
    elseif ((i.gt.ll+1).and.(i.le.L-ll)) then
      f31(i)=omega(i-1)*f31(i-1)+omega(L)*f32(i,L)-f21(i)
    else
      f31(i)=omega(i-1)*f31(i-1)-f21(i)
    endif
    f31(i)=f31(i)/e0
    f3(i)=f3(i)+f31(i)-b1*f2(i)-b2*f1(i)
  enddo

  ! objective function value
  val=0.0d0
  do i=2,L
    rhoi=f1(i)+f2(i)+f3(i)
    val=val+(rhoi-rhoexp(i))*(rhoi-rhoexp(i))
  enddo
  Niter=Niter+1
  return
end subroutine fnc3

! runs MC simulations
subroutine mc(L,ll,omega,tau,iter,t,rho,current)

  use mt19937
  implicit none

  integer(kind=4),intent(in) :: L,ll
  real(kind=8),intent(in) :: omega(2:L)
  integer(kind=4),intent(inout) :: tau(2:L)
  integer(kind=8),intent(in) :: iter
  real(kind=8),intent(inout) :: t
  real(kind=8),intent(out) :: rho(2:L),current(1:L)

  ! local variables
  integer(kind=8) :: i,k
  real(kind=8) :: r1,r2,a0,atest,dt
  real(kind=8),allocatable :: a(:)

  ! initial propensity function
  allocate(a(1:L))
  if (sum(tau(2:ll+1)).EQ.0) then
    a(1)=1.0d0
  else
    a(1)=0.0d0
  endif
  do i=2,L-ll
    a(i)=omega(i)*dble(tau(i))*dble(1-tau(i+ll))
  enddo
  do i=L-ll+1,L
    a(i)=omega(i)*dble(tau(i))
  enddo  

  t=0.0d0
  rho=0.0d0
  current=0.0d0

  do k=1,iter*L

    a0=sum(a)
   
    r1=grnd()
    r2=grnd()

    ! chooses time until the next event
    dt=-dlog(dble(r1))/a0
    t=t+dt

    ! density profile
    rho=rho+dble(tau)*dt

    ! chooses which event happens next
    i=1
    atest=a(1)
    do while (r2.GT.atest/a0)
      i=i+1
      atest=atest+a(i)
    enddo

    ! updates the propensity function and moves the particle
    a(i)=0.0d0
    if (i.EQ.1) then
      a(2)=omega(2)*dble(1-tau(2+ll))
      tau(2)=1
    elseif ((i.GE.2).AND.(i.LE.ll)) then
      a(i+1)=omega(i+1)*dble(1-tau(i+1+ll))
      tau(i)=0
      tau(i+1)=1
    elseif (i.EQ.ll+1) then
      a(1)=1.0d0
      a(i+1)=omega(i+1)*dble(1-tau(i+1+ll))
      tau(i)=0
      tau(i+1)=1 
    elseif ((i.GE.ll+2).AND.(i.LE.L-ll-1)) then
      a(i-ll)=omega(i-ll)*dble(tau(i-ll))
      a(i+1)=omega(i+1)*dble(1-tau(i+1+ll))
      tau(i)=0
      tau(i+1)=1
    elseif ((i.GE.L-ll).AND.(i.LE.L-1)) then
      a(i-ll)=omega(i-ll)*dble(tau(i-ll))
      a(i+1)=omega(i+1)
      tau(i)=0
      tau(i+1)=1
    else
      a(i-ll)=omega(i-ll)*dble(tau(i-ll))
      tau(L)=0
    endif
    current(i)=current(i)+1.0d0
  enddo
  deallocate(a)
  rho=rho/t
  current=current/t
  return
end subroutine mc
