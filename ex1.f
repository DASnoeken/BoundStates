c     Program to calculate the harmonic oscillator eigenvalues.
c     Computer exercise 1 for the quantum dynamics course adapted for FORTRAN
c     compile with gfortran -g -c -O ex1.f
c     then         gfortran -o ex1 ex1.o -llapack
c     We need the lapack library for the zheev function
c     I also put the compilation commands in a bash script called "C"
c     DAS, 11-4-2019

      program ex1
      implicit none
      integer i, j, info, a,b,c,d
      integer maxX,n,m
      parameter(maxX=5,n=499,m=n+1)
      double precision delta 
      double precision, dimension(m) :: x
      double precision, dimension(m,m) :: V, T
      complex*16, dimension(m,m) :: H
      double precision, dimension(m) :: E, Eexact, error
      integer lwork
      parameter(lwork=2*(m)-1)
      real, dimension(3*(m)-2) :: rwork
      complex*16, dimension(lwork) :: work
      character*4 eval
      
      info=0
      do a=1,3*m-2           					!I try to initialize as much as possible
        rwork(a)=0
      enddo
      delta=2*dfloat(maxX)/dfloat(n)				!difference between each point
      do i=1,n+1
        x(i)=-dfloat(maxX)+dfloat(i-1)*delta			!fill x vector
      enddo			

      do c=1,m
        do d=1,m
          T(c,d)=0						!Initialize T with all zeros
          H(c,d)=0						!Initialize H with all zeros
        enddo
      enddo
      do j=1,m							!Kinetic energy operator
        T(j,j)=-2
        if (j.ne.1.and.j.ne.m) then				!Cases for when T entries are 1
          T(j,j+1)=1
          T(j,j-1)=1
        elseif(j.eq.1) then					!edges
          T(j,j+1)=1
        elseif(j.eq.m) then
          T(j,j-1)=1
        endif
      enddo
      T=T*(-1.d0/(2.d0*delta*delta))				!adjust with constant
      do c=1,m
        do d=1,m
          if (c.eq.d) then 
               V(c,d)=0.5*x(c)*x(d)				!Potential (diagonal case)
          else 
               V(c,d)=0						!Potential (off diagonals are all zero)
          endif
        enddo
      enddo
      do a=1,m
        E(a)=0
        Eexact(a)=(a-1)+0.5					!Exact solutions to the harmonic oscillator energies
        do b=1,m
          H(a,b)=T(a,b)+V(a,b)		        		!hamiltonian
        enddo
      enddo
      call zheev('N','U',m,H,m,E,work,lwork,rwork,info)		!diagonalize, to get eigenvectors you should change 'N' to 'V' and they will be stored in H. Doesn't seem to work yet
      if(info.ne.0) then					!Also uncomment below lines related to file unit=12 for eigenvectors
        write(*,*)'Error!', info
        stop 111
      endif
      open(unit=10,file='eigenvals.dat',form='formatted')
      open(unit=11,file='errors.dat',form='formatted')
c      open(unit=12,file='eigenvec.dat',form='formatted')
      do a=1,m
        write(10,*) E(a)					!Write eigenvalues to eigenvals.dat
        error(a)=Eexact(a)-E(a)					!Calculate difference between exact solutions and second order finite difference solutions in absolute numbers
        if(abs(1.d2*error(a)) .le. 1.d0) then			!Set a threshold value to 1% error from the actual value
          eval='Good'						!Prints Good if error is below 1%
        else
          eval='Bad'						!Prints Bad if error is over 1%
        endif
        write(11,*) error(a),'   ',eval
c        do b=1,m
c          write(12,*) H(a,b)
c        enddo
c        write(12,*)'----------------------------------------'
      enddo
      close(10)
      close(11)
c      close(12)
      end
