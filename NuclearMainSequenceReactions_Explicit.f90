!*****************************************************************************
!                                                                            *
!  PROGRAM: Nuclear Main Sequence Reactions                                  *
!                                                                            *
!  PURPOSE:  This Program simulate the evolution in time of the percentual   *
!  abbundances of light elements involved in nuclear reactions in main       *
!  Sequence(Sun Like) stars, the simulation will be computed using a explicit*
!  Euler method and it will be possible to select a fixed precision with a   *
!  fixed time lapse or a adapting precision per step using a variable time   *
!  lapse automatically fixed based on the percentual precision requested.     *
!                                                                            *
!*****************************************************************************

! Main Program
    program NuclearMainSequenceReactions_Explicit

    implicit none

    ! Variables !

  
    integer (kind=4) :: i,o,op
    
    real (kind=8) :: t,dt,tend,Temp,Temp9,Temp9A,rho,cs1,cs3,cs4,TSTART,TSTOP
    real (kind=8) :: Z
    ! abbundances are defined as vector so that during simulation second slot refer to step n and first slot to step n+1
    real (kind=8,len=2) y1,y3,y4
    
    character (len=1) :: c
    
    logical (kind=1) opt
    

    ! Body of NuclearMainSequenceReactions !

    ! Defining initial values
     
    opt= .TRUE.

    do while(opt)

        t=0.

1   	write (*,*) '-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~'
        write (*,*) 'Please insert the structure central temperature in kelvins(>0)'
    	read (*,*) Temp
    	if(Temp<=0.) then
		write (*,*) 'Warning! Temperature must be positive!'
        	goto 1
   	endif
        
        Temp9 = Temp*1.D-9

2   	write (*,*) '-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~'
        write (*,*) 'Please insert the structure central density in g/cm^3(>0)'
    	read (*,*) rho
    	if(rho<=0.) then
		write (*,*) 'Warning! Density must be positive!'
        	goto 1
   	endif

3
        write (*,*) '-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~'
        write (*,*) 'Please insert the initial abundance for 3He(0<=y3<=1)'
    	read (*,*) y3(2)
    	if((y3(2)<0.).or.y3(2)>1.) then
		write (*,*) 'Warning! abundances must be in [0,1] interval!'
        	goto 3
   	endif
        y3(2) = y3(2)/3.

4
        write (*,*) '-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~'
        write (*,*) 'Please insert the initial abundance for 4He(0<=y4<=1)'
    	read (*,*) y4(2)
    	if((y4(2)<0.).or.y4(2)>1.) then
		write (*,*) 'Warning! abundances must be in [0,1] interval!'
        	goto 4
   	endif
        if((y4(2)+ (3.*y3(2)))>1. ) then
                write (*,*) 'Warning! total abbundance must be in [0,1] interval!'
        	goto 3
        endif
        y4(2) = y4(2)/4.

5
        write (*,*) '-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~'
        write (*,*) 'Please insert the metallicity of the structure(0<=Z<=1)'
    	read (*,*) Z
    	if((Z<0.).or.Z>1.) then
		write (*,*) 'Warning! metallicity must be in [0,1] interval!'
        	goto 5
   	endif
        if(((4.*y4(2))+ (3.*y3(2))+Z)>1. ) then
                write (*,*) 'Warning! total abbundance must be in [0,1] interval!'
        	goto 3
        endif

        write (*,*) '-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~'
        write (*,*) 'The abbundance of H will be calculated from other values'
        y1 = 1. - Z - 4.*y4(2) - 3.*y3(2)

    ! Calculation of reactions cross sections in function of selected temperature

 
       call crosssection1 (Temp9,cs1)
       call crosssection3 (Temp9,cs3)
       call crosssection4 (Temp9,cs4)

    ! Opening support files for the memorization of the results

    	open (unit = 1, file = 'RisultatiEsercitazione2.txt')
    	write (1,*) '################################################################'
   	write (1,*) '#           Risultati numerici dell`esercitazione 2            #'
    	write (1,*) '################################################################'
    	write (1,*) ''
    	write (1,*) '#    t(years)      y1       y3      y4                         #'
    	write (1,*) ''
101    	format (4f20.15)
    
      ! Beginning of the menu
11        
        write (*,*) '-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~'
        write (*,*) 'Please select simulation method :'
        write (*,*) '1 -  Simulation with fixed time lapse (faster)'
        write (*,*) '2 -  Simulation with variable time lapse (greater accuracy)'
        write (*,*) '0 - Stop Program Execution'
   	write (*,*) '-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~'
        write (*,*) 'Make your choice :'
    	read (*,*) o
    	if(o<0.or.o>2) then
		write (*,*) 'Warning ! Selected simulation isn`t available'
        	goto 11
   	endif
 
        write (*,*) '-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~'

        SELECT CASE (o) 
              
        ! Fixed time lapse simulation
          case (1) 
     
               ! Setting the time lapse of the simulation
12             write (*,*) '-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~'
               write (*,*) 'Please insert the time lapse used for the simulation'
               read (*,*) dt
    	       if(dt<=0) then
		   write (*,*) 'Warning! Time lapse must be positive!'
        	   goto 12
   	       endif
               
    	       
               ! Defining the ending condition
13             write (*,*) '-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~'
               write (*,*) 'Select the stop condition :'
               write (*,*) '1 -  Simulate untill time Tend as been reached'
               write (*,*) '2 -  Simulate untill hydrogen abbundance has reached a value y1min'
               write (*,*) '0 - Stop Program Execution'
   	       write (*,*) '-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~'
               write (*,*) 'Make your choice :'
    	       read (*,*) op
    	       if(op<0.or.op>2) then
		    write (*,*) 'Warning ! Selected simulation isn`t available'
        	    goto 13
   	       endif
              
               Select case(op)
               
               case (1)
                    ! Setting the time lapse of the simulation
14                  write (*,*) '-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~'
                    write (*,*) 'Please insert the time of ending for the simulation'
                    read (*,*) tend
    	            if(tend<=0) then
		       write (*,*) 'Warning! Tend must be positive!'
        	       goto 14
   	            endif

                    write (*,*) 'The numerical simulation has started !'
                    ! Inizialization of system tyme for the determination of total time elapsed
    	            call CPU_TIME(TSTART)
    	            write (*,*) '-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~'
    	            do while (t.lt.tend)
        	       write (1,101) (t/31536000.),y1(2),y3(2),y4(2)
		       t = t + dt
                       if(t>tend) t=tend
        	       y1(1) = y1(2) + dt*rho*(-(3./2.)*(y1(2)^2.)*cs1 - y3(2)*y4(2)*cs4 + (y3(2)**2.)*cs3)
                       y3(1) = y3(2) + dt*rho*(-(y3(2)**2.)*cs3 - y3(2)*y4(2)*cs4 + 0.5*(y1(2)**2.)*cs1)
                          
    	            end do
               
               !Ending program
               case (0) 
         
	           goto 999

               end select


          !
          case (2) 
            
           
                
          
          !Ending program
          case (0) 
         
	       goto 999

          end select

    ! Closing files
      close (unit = 1)
      
        

    ! creation of the graphics
      CALL SYSTEM ('gnuplot -persist RisultatiEsercitazione1.plt')     
    
    ! Tempo di calcolo finale della CPU
    	call CPU_TIME (TSTOP)
    
    ! Execution Completed!
        write (*,*) '-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~'
        write (*,*) 'Simulation completed in ',TSTOP-TSTART,' seconds'

999   	write (*,*) '-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~'
        write (*,*) 'Do you wish to relaunch the program (y/n)?'
    	read (*,*) c
    	if(c/='y'.and.c/='n'.and.c/='Y'.and.c/='N') then
		write (*,*) 'Wrong Immission!'
        	goto 999
   	endif
        
        write (*,*) '-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~'
        if(c=='n'.or.c=='N') then
		opt = .FALSE.
   	endif
    end do  

!

    ! End Program !
    stop
    end program NuclearMainSequenceReactions_Explicit

! Subroutine Section
    
   !Subroutine for the calculation of p-p reaction's cross section 
    subroutine crosssection1 (Temp9,cs1)
        real (kind=8) :: Temp9,cs1
        cs1 = (3.82E-15/(Temp9**(2./3.)*DEXP(-3.38/(Temp9**(1./3.)))*(1.+0.123* &
        (Temp9**(1./3.))+1.09*((Temp9**(2./3.))+0.938*Temp9))
        return
        end
   
    !Subroutine for the calculation of he3 reaction's cross section 
    subroutine crosssection3 (Temp9,cs3)
        real (kind=8) :: Temp9,cs3
        cs3 = (5.96E10/(Temp9**(2./3.)*DEXP(-12.276/(Temp9**(1./3.))*(1.+0.034* &
        (Temp9**(1./3.)-0.199*(Temp9**(2./3.)-.047*Temp9+.032*(Temp9**(4./3.)+0.019*(Temp9**(5./3.)))
        return
        end

    !Subroutine for the calculation of he4 reaction's cross section 
    subroutine crosssection4 (Temp9,cs4)
        real (kind=8) :: Temp9,Temp9A,T9A13,T9A56,cs4
        Temp9A=Temp9/(1.D0+4.95D-2*Temp9)
        WT9A=DLOG(Temp9A)
        T9A13=DEXP(1.D0/3.D0*WT9A)
        T9A56=DEXP(5.D0/6.D0*WT9A)
        cs4 = (5.97E6*T9A56/Temp9**(3./2.)*DEXP(-12.826/T9A13))
        return
        end
