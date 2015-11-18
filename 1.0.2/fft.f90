! one dimensional fourier transform, to use this subroutine 
! 1.define {emin,de} and make them common in the main program
! 2.define an output array for {out}
	program main 

	implicit real*8(a-h,o-z) 

!        integer*4, parameter :: nt = 2001 

        real*8, allocatable, dimension(:) :: t
	complex*16, allocatable, dimension(:) :: out,data
        complex*16 :: im
        
	im = (0d0,1d0)

        open(103,file='fft.dat')

        open(11,file='INFT') 
        read(11,*) nt 
        close(11) 

        allocate(t(nt),out(nt),data(nt))

        open(10,file='corr')
        do i=1,nt 
          read(10,*) t(i),a,b,c,d 
          data(i) = a + b*im 
        enddo 
        close(10)
 
        dt = t(2) - t(1) 
        Emax = -0.01d0

	Pi = 4.0d0*atan(1d0)

!        de = 1d0/(Nt*dt)
        de = 0.001
        emin = -0.18d0
        
        out = (0d0,0d0)
        alfa = -log(1d-4/abs(data(nt)))/t(nt)**2

        write(*,*) 'damping coefficient', afa 

	i = 1
	do while (Emin+(i-1)*de < Emax)
        
          en = Emin + (i-1)*de
          do j=1,Nt
            out(i) = data(j)*exp(im*en*t(j))*exp(-alfa*t(j)**2)*dt + out(i)
          enddo

          write(103,1000) en,out(i)
          
          i = i+1 

        enddo
        
1000    format(1000(e14.7,1x))

        return 
        end program 


        
        
        
