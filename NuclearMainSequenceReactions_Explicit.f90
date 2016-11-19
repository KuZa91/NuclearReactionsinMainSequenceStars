!*****************************************************************************
!                                                                            *
!  PROGRAM: Nuclear Main Sequence Reactions                                  *
!                                                                            *
!  PURPOSE:  This Program simulate the evolution in time of the percentual   *
!  abbundances of light elements involved in nuclear reactions in main       *
!  Sequence(Sun Like) stars, the simulation will be computed using a explicit*
!  Euler method and it will be possible to select a fixed precision with a   *
!  fixed time lapse or a adapting precision per step using a variable time   *
!  lapse automatically fixed based on the percentual precision requested.    *
!                                                                            *
!*****************************************************************************

! Main Program
    program NuclearMainSequenceReactions_Explicit

    implicit none

    ! Variables !

  
    integer (kind=4) :: i,wp,o,op
    
    real (kind=8) :: t,dt,tend,Temp,Temp9,Temp9A,rho,cs1,cs3,cs4,cs12,cs14c,cs14o,cs16,TSTART,TSTOP,Z
    ! abbundances are defined as vector so that during simulation second slot refer to step n and first slot to step n+1
    real (kind=8), dimension (2) :: y1,y3,y4,y12,y14,y16
    
    character (len=1) :: c
    
    logical (kind=1) :: opt
    

    ! Body of NuclearMainSequenceReactions !

    ! Defining initial values
     
    opt = .TRUE.

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
        	goto 2
   	endif

3       write (*,*) '-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~'
        write (*,*) 'Please insert the initial abundance for 3He(0<=y3<=1)'
    	read (*,*) y3(2)
    	if((y3(2)<0.).or.y3(2)>1.) then
		write (*,*) 'Warning! abundances must be in [0,1] interval!'
        	goto 3
   	endif
        y3(2) = (y3(2)/3.)

4       write (*,*) '-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~'
        write (*,*) 'Please insert the initial abundance for 4He(0<=y4<=1)'
    	read (*,*) y4(2)
    	if((y4(2)<0.).or.y4(2)>1.) then
		write (*,*) 'Warning! abundances must be in [0,1] interval!'
        	goto 4
   	endif
        if((y4(2)+ (3.*y3(2)))>1. ) then
                write (*,*) 'Warning! total abundance must be in [0,1] interval!'
        	goto 3
        endif
        y4(2) = (y4(2)/4.)

5       write (*,*) '-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~'
        write (*,*) 'Please insert the metallicity of the structure(0<=Z<=1)'
    	read (*,*) Z
    	if((Z<0.).or.Z>1.) then
		write (*,*) 'Warning! metallicity must be in [0,1] interval!'
        	goto 5
   	endif
        if(((4.*y4(2))+ (3.*y3(2))+Z)>1. ) then
                write (*,*) 'Warning! total abundance must be in [0,1] interval!'
        	goto 3
        endif

6       write (*,*) '-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~'
        write (*,*) 'Please insert the initial abundance for 12C(0<=y12<=Z)'
    	read (*,*) y12(2)
    	if((y12(2)<0.).or.y12(2)>Z) then
		write (*,*) 'Warning! abundances must be in [0,Z] interval!'
        	goto 6
   	endif
        y12(2) = (y12(2)/12.)

7       write (*,*) '-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~'
        write (*,*) 'Please insert the initial abundance for 14N(0<=y14<=Z)'
    	read (*,*) y14(2)
    	if((y14(2)<0.).or.y14(2)>Z) then
		write (*,*) 'Warning! abundances must be in [0,Z] interval!'
        	goto 7
   	endif
        if(((12.)*y12(2)+ y14(2)) > Z ) then
                write (*,*) 'Warning! total abundance must be in [0,Z] interval!'
        	goto 6
        endif
        y14(2) = (y14(2)/14.)

8       write (*,*) '-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~'
        write (*,*) 'Please insert the initial abundance for 16O(0<=y16<=Z)'
    	read (*,*) y16(2)
    	if((y16(2)<0.).or.y16(2)>Z) then
		write (*,*) 'Warning! abundances must be in [0,Z] interval!'
        	goto 8
   	endif
        if(((12.)*y12(2) + (14.)*y14(2) + y16(2))> Z ) then
                write (*,*) 'Warning! total abundance must be in [0,Z] interval!'
        	goto 6
        endif
        y16(2) = (y16(2)/16.)

        write (*,*) '-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~'
        write (*,*) 'The abundance of H will be calculated from other values'
        y1 = 1. - Z - 4.*y4(2) - 3.*y3(2)

9      	write (*,*) '-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~'
        write (*,*) 'On each many steps do you wish to write on file ?(suggested 1000)'
    	read (*,*) wp 
    	if(wp<=0.) then
		write (*,*) 'Warning! Writing Precision must be positive!'
        	goto 9
   	endif

    ! Calculation of reactions cross sections in function of selected temperature

 
       call crosssection1 (Temp9,cs1)
       call crosssection3 (Temp9,cs3)
       call crosssection4 (Temp9,cs4)
       call crosssection12 (Temp9,cs12)
       call crosssection14O (Temp9,cs14o)
       call crosssection14C (Temp9,cs14c)
       call crosssection16 (Temp9,cs16)

    ! Opening support files for the memorization of the results

    	open (unit = 1, file = 'NuclearAbundances.txt')
    	write (1,*) '################################################################'
   	write (1,*) '#           Risultati numerici dell`esercitazione 2            #'
    	write (1,*) '################################################################'
    	write (1,*) ''
    	write (1,*) '#    t(years)      y1       y3      y4    y12    y14    y16    #'
    	write (1,*) ''
101    	format (f25.5,6f20.15)
    
      ! Beginning of the menu

11      write (*,*) '-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~'
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

        SELECT CASE (o) 
              
        ! Fixed time lapse simulation
          case (1) 
     
               ! Setting the time lapse of the simulation
12             write (*,*) '-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~'
               write (*,*) 'Please insert the time lapse used for the simulation in seconds'
               read (*,*) dt
    	       if(dt<=0) then
		   write (*,*) 'Warning! Time lapse must be positive!'
        	   goto 12
   	       endif
               
    	       
               ! Defining the ending condition
13             write (*,*) '-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~'
               write (*,*) 'Select the stop condition :'
               write (*,*) '1 -  Simulate untill time Tend as been reached'
               write (*,*) '2 -  Simulate untill hydrogen abundance has reached a value y1min'
               write (*,*) '0 - Stop Program Execution'
   	       write (*,*) '-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~'
               write (*,*) 'Make your choice :'
    	       read (*,*) op
    	       if(op<0.or.op>2) then
		    write (*,*) 'Warning ! Selected simulation isn`t available'
        	    goto 13
   	       endif
              
               Select case(op)
               
               ! Run untill Tend is reached

               case (1)
                    ! Setting the ending time of the simulation
14                  write (*,*) '-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~'
                    write (*,*) 'Please insert the time of ending for the simulation in years'
                    read (*,*) tend
    	            if(tend<=0) then
		       write (*,*) 'Warning! Tend must be positive!'
        	       goto 14
   	            endif
                    tend = tend * 31536000. 
                   ! Beginning of the numerical simulation
                    write (*,*) 'The numerical simulation has started !'
                    ! Inizialization of system tyme for the determination of total time elapsed
    	            call CPU_TIME(TSTART)
    	            write (*,*) '-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~'
                    i=0
    	            do while (t.lt.tend)
        	       ! The program will write on file on each 1000 steps
                       if((i.eq.0).or.(i.eq.wp)) then
                           write (1,101) (t/(31536000.)),y1(2),y3(2),y4(2),y12(2),y14(2),y16(2)
                           i=0
                       endif
		       t = t + dt
                       i = i + 1
                       if(t>tend) t=tend
                       
        	       y1(1) = y1(2) + dt*rho*(-(3./2.)*(y1(2)**2.)*cs1 - y3(2)*y4(2)*cs4 + (y3(2)**2.)*cs3 &
                       - (2.)*y1(2)*y12(2)*cs12 - (2.)*y1(2)*y14(2)*cs14c - (2.)*y1(2)*y14(2)*cs14o &
                       -(2.)*y1(2)*y16(2)*cs16)
                       
                       y3(1) = y3(2) + dt*rho*(-(y3(2)**2.)*cs3 - y3(2)*y4(2)*cs4 + 0.5*(y1(2)**2.)*cs1)

                       y4(1) = y4(2) + dt*rho*((0.5)*(y3(2)**2)*cs3 + y3(2)*y4(2)*cs4 + y1(2)*y14(2)*cs14c&
                       + y1(2)*y16(2)*cs16)
                       
                       y12(1) = y12(2) + dt*rho*(-y1(2)*y12(2)*cs12 + y1(2)*y14(2)*cs14c)

                       y14(1) = y14(2) + dt*rho*(y1(2)*y12(2)*cs12 - y1(2)*y14(2)*cs14c - y1(2)*y14(2)*cs14o &
                       + y1(2)*y16(2)*cs16)
                       
                       y16(1) = y16(2) + dt*rho*(- y1(2)*y16(2)*cs16 + y1(2)*y14(2)*cs14o)

                       y1(2) = y1(1)
                       y3(2) = y3(1)
                       y4(2) = y4(1)
                       y12(2) = y12(1)
                       y14(2) = y14(1)
                       y16(2) = y16(1)
    	            end do


               ! Run untill the hydrogen abundance has reached a value y1min 

               case (2) 
                    ! Setting the minimum value of y1
15                  write (*,*) '-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~'
                    write (*,*) 'Please insert the minimum abbundance y1 that has to be reached in the simulation'
                    read (*,*) tend
                    ! The program will write on file on each 1000 steps
    	            if((tend<=0.).or.tend>1.) then
		       write (*,*) 'Warning! Tend must be in [0,1] range!'
        	       goto 15
   	            endif
                    ! Beginning of the numerical simulation
                    write (*,*) 'The numerical simulation has started !'
                    ! Inizialization of system tyme for the determination of total time elapsed
    	            call CPU_TIME(TSTART)
    	            write (*,*) '-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~'
                    i=0
    	            do while (y1(2).gt.tend)
                     
        	       if((i.eq.0).or.(i.eq.wp)) then
                           write (1,101) (t/31536000.),y1(2),y3(2),y4(2),y12(2),y14(2),y16(2)
                           i=0
                       endif
		       t = t + dt
                       i = i + 1

        	       y1(1) = y1(2) + dt*rho*(-(3./2.)*(y1(2)**2.)*cs1 - y3(2)*y4(2)*cs4 + (y3(2)**2.)*cs3 &
                       - (2.)*y1(2)*y12(2)*cs12 - (2.)*y1(2)*y14(2)*cs14c - (2.)*y1(2)*y14(2)*cs14o &
                       -(2.)*y1(2)*y16(2)*cs16)
                       
                       y3(1) = y3(2) + dt*rho*(-(y3(2)**2.)*cs3 - y3(2)*y4(2)*cs4 + 0.5*(y1(2)**2.)*cs1)

                       y4(1) = y4(2) + dt*rho*((0.5)*(y3(2)**2)*cs3 + y3(2)*y4(2)*cs4 + y1(2)*y14(2)*cs14c&
                       + y1(2)*y16(2)*cs16)
                       
                       y12(1) = y12(2) + dt*rho*(-y1(2)*y12(2)*cs12 + y1(2)*y14(2)*cs14c)

                       y14(1) = y14(2) + dt*rho*(y1(2)*y12(2)*cs12 - y1(2)*y14(2)*cs14c - y1(2)*y14(2)*cs14o &
                       + y1(2)*y16(2)*cs16)
                       
                       y16(1) = y16(2) + dt*rho*(- y1(2)*y16(2)*cs16 + y1(2)*y14(2)*cs14o)

                       y1(2) = y1(1)
                       y3(2) = y3(1)
                       y4(2) = y4(1)
                       y12(2) = y12(1)
                       y14(2) = y14(1)
                       y16(2) = y16(1)
                       
    	            end do
                    write (*,*) 'y1 has reached ',y1(2),' abundance after ',(t/31536000.),' years!'
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
      CALL SYSTEM ('gnuplot -persist NuclearAbundancesPP.plt')
      CALL SYSTEM ('gnuplot -persist NuclearAbundancesCNO.plt')     
    
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
        cs1 = (3.82E-15/(Temp9**(2./3.))*DEXP(-3.38/(Temp9**(1./3.)))*(1.+0.123* &
        (Temp9**(1./3.))+1.09*(Temp9**(2./3.))+0.938*Temp9))
        return
        end
   
    !Subroutine for the calculation of he3 reaction's cross section 
    subroutine crosssection3 (Temp9,cs3)
        real (kind=8) :: Temp9,cs3
        cs3 = (5.96E10/(Temp9**(2./3.))*DEXP(-12.276/(Temp9**(1./3.)))*(1.+0.034*(Temp9**(1./3.)) &
        - 0.199*(Temp9**(2./3.))- 0.047*Temp9 + 0.032*(Temp9**(4./3.))+0.019*(Temp9**(5./3.))))
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

    !Subroutine for the calculation of C12 reaction's cross section 
    subroutine crosssection12 (Temp9,cs12)
        real (kind=8) :: Temp9,cs12
        cs12 = (2.04D7/(Temp9**(2./3.))*DEXP(-13.690/(Temp9**(1./3.))-(Temp9/1.500)**2) &
        *(1.+.030*(Temp9**(1./3.))+1.19*(Temp9**(2./3.))+.254*Temp9+2.06*(Temp9**(4./3.))&
        +1.12*(Temp9**(5./3.))) +1.08D5/(Temp9**(3./2.))*DEXP(-4.925/Temp9) & 
        +2.15D5/(Temp9**(3./2.))*DEXP(-18.179/Temp9))
        return
        end

    !Subroutine for the calculation of N14 to C12 reaction's cross section 
    subroutine crosssection14C (Temp9,cs14c)
        real (kind=8) :: Temp9,cs14c,S5,S32,S7

        S5=(5.08D7/(Temp9**(2./3.))*DEXP(-15.228/(Temp9**(1./3.))-(Temp9/3.090)**2) &
        *(1.+.027*(Temp9**(1./3.))-.778*(Temp9**(2./3.))-.149*Temp9+.261*(Temp9**(4./3.)) &
        +.127*(Temp9**(5./3.)))+2.28D3/(Temp9**(3./2.))*DEXP(-3.011/Temp9)&
        +1.65D4*(Temp9**(1./3.))*DEXP(-12.007/Temp9))
        
        S32=(9.78D8/(Temp9**(2./3.))*DEXP(-15.251/(Temp9**(1./3.))-(Temp9/0.450)**2) &
        *(1.+.027*(Temp9**(1./3.))+.219*(Temp9**(2./3.))+0.042*Temp9+6.83*(Temp9**(4./3.))&
        +3.32*(Temp9**(5./3.)))+1.11D4/(Temp9**(3./2.))*DEXP(-3.328/Temp9)&
        +1.49D4/(Temp9**(3./2.))*DEXP(-4.665/Temp9)+3.80D6/(Temp9**(3./2.))*DEXP(-11.048/Temp9))

        S7=(1.08D12/(Temp9**(2./3.))*DEXP(-15.251/(Temp9**(1./3.))-(Temp9/0.522)**2) &
        *(1.+.027*(Temp9**(1./3.))+2.62*(Temp9**(2./3.))+0.501*Temp9+5.36*(Temp9**(4./3.)) &
        +2.60*(Temp9**(5./3.)))+1.19D8/(Temp9**(3./2.))*DEXP(-3.676/Temp9)+5.41D8/(Temp9**(1./2.))&
        *DEXP(-8.926/Temp9)+0.1d0*4.72D8/(Temp9**(3./2.))*DEXP(-7.721/Temp9)+2.20D9/(Temp9**(3./2.))&
        *DEXP(-11.418/Temp9))
     
        cs14c=S5*S7/(S7+S32)

        return
        end

    !Subroutine for the calculation of N14 to O16 reaction's cross section 
    subroutine crosssection14O (Temp9,cs14o)
        real (kind=8) :: Temp9,cs14o,S5,S32,S7

        S5=(5.08D7/(Temp9**(2./3.))*DEXP(-15.228/(Temp9**(1./3.))-(Temp9/3.090)**2) &
        *(1.+.027*(Temp9**(1./3.))-.778*(Temp9**(2./3.))-.149*Temp9+.261*(Temp9**(4./3.)) &
        +.127*(Temp9**(5./3.)))+2.28D3/(Temp9**(3./2.))*DEXP(-3.011/Temp9)&
        +1.65D4*(Temp9**(1./3.))*DEXP(-12.007/Temp9))
        
        S32=(9.78D8/(Temp9**(2./3.))*DEXP(-15.251/(Temp9**(1./3.))-(Temp9/0.450)**2) &
        *(1.+.027*(Temp9**(1./3.))+.219*(Temp9**(2./3.))+0.042*Temp9+6.83*(Temp9**(4./3.))&
        +3.32*(Temp9**(5./3.)))+1.11D4/(Temp9**(3./2.))*DEXP(-3.328/Temp9)&
        +1.49D4/(Temp9**(3./2.))*DEXP(-4.665/Temp9)+3.80D6/(Temp9**(3./2.))*DEXP(-11.048/Temp9))

        S7=(1.08D12/(Temp9**(2./3.))*DEXP(-15.251/(Temp9**(1./3.))-(Temp9/0.522)**2) &
        *(1.+.027*(Temp9**(1./3.))+2.62*(Temp9**(2./3.))+0.501*Temp9+5.36*(Temp9**(4./3.)) &
        +2.60*(Temp9**(5./3.)))+1.19D8/(Temp9**(3./2.))*DEXP(-3.676/Temp9)+5.41D8/(Temp9**(1./2.))&
        *DEXP(-8.926/Temp9)+0.1d0*4.72D8/(Temp9**(3./2.))*DEXP(-7.721/Temp9)+2.20D9/(Temp9**(3./2.))&
        *DEXP(-11.418/Temp9))

        cs14o = S5*(1.D0-S7/(S7+S32))

        return
        end
    
    !Subroutine for the calculation of C12 reaction's cross section 
    subroutine crosssection16 (Temp9,cs16)
        real (kind=8) :: Temp9,cs16
        cs16 = (1.50D8/((Temp9**(2./3.))*(1.+2.13*(1.-DEXP(-.728*(Temp9**(2./3.))))))* &
        DEXP(-16.692/(Temp9**(1./3.))))
        return
        end


