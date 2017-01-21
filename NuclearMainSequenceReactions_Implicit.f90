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
    program NuclearMainSequenceReactions_Implicit

    implicit none

    ! Variables !

  
    integer (kind=4) :: i,j,wp,o,op,Nreac,check
    
    real (kind=8) :: t,dt,tend,Temp,Temp9,Temp9A,rho,TSTART,TSTOP,Z,tollE,tollDY,Emax,DYmax
    ! abbundances are defined as vector so that during simulation second slot refer to step n and first slot to step n+1
    
    real (kind=8), dimension (2) :: y1,y3,y4,y12,y14,y16

    real (kind=8), dimension (6) :: E

    real (kind=8), dimension (7) :: cs
    
    real (kind=8), dimension (6,6) :: A,MXapp
    
    character (len=1) :: c
    
    logical (kind=1) :: opt
    

    ! Body of NuclearMainSequenceReactions !

    ! Defining initial values
     
    opt = .TRUE.
    Nreac = 6

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

10      write (*,*) '-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~'
        write (*,*) 'Please insert the time lapse used for the simulation in seconds'
        read (*,*) dt
    	if(dt<=0) then
	        write (*,*) 'Warning! Time lapse must be positive!'
        	goto 10
   		endif

    ! Calculation of reactions cross sections in function of selected temperature

 
        call crosssection1 (Temp9,cs(1))
        call crosssection3 (Temp9,cs(2))
        call crosssection4 (Temp9,cs(3))
        call crosssection12 (Temp9,cs(4))
        call crosssection14O (Temp9,cs(5))
        call crosssection14C (Temp9,cs(6))
        call crosssection16 (Temp9,cs(7))

    ! Opening support files for the memorization of the results

    	open (unit = 1, file = 'NuclearAbundances.txt')
    	write (1,*) '################################################################'
   	    write (1,*) '#           Risultati numerici dell`esercitazione 2            #'
    	write (1,*) '################################################################'
    	write (1,*) ''
    	write (1,*) '#    t(years)      y1       y3      y4    y12    y14    y16    #'
    	write (1,*) ''
101    	format (f25.5,6f20.15)
    
        ! Defining the ending condition
11      write (*,*) '-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~'
        write (*,*) 'Select the stop condition :'
        write (*,*) '1 -  Simulate untill time Tend as been reached'
        write (*,*) '2 -  Simulate untill hydrogen abundance has reached a value y1min'
        write (*,*) '0 - Stop Program Execution'
   	    write (*,*) '-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~'
        write (*,*) 'Make your choice :'
        read (*,*) o
    	if(o<0.or.o>2) then
		    write (*,*) 'Warning ! Selected simulation isn`t available'
        	goto 11
   	    endif

        SELECT CASE (o) 
              
            ! Run untill Tend is reached
            case (1) 
     
                ! Setting the ending time of the simulation

12              write (*,*) '-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~'
                write (*,*) 'Please insert the time of ending for the simulation in years'
                read (*,*) tend
    	        if(tend<=0) then
		        	write (*,*) 'Warning! Tend must be positive!'
        	        goto 12
   	          	endif
                tend = tend * 31536000. 

13              write (*,*) '-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~'
                write (*,*) 'Please insert the tollerance on the E value for the simulation'
                read (*,*) tollE
    	        if(tollE<=0) then
		        	write (*,*) 'Warning! TollE must be positive!'
        	        goto 13
   	          	endif

14              write (*,*) '-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~'
                write (*,*) 'Please insert the tollerance on DY for the simulation'
                read (*,*) tollDY
    	        if(tollDY<=0) then
		       	  	write (*,*) 'Warning! TollDY must be positive!'
        	        goto 14
   	            endif

                ! Beginning of the numerical simulation
                write (*,*) 'The numerical simulation has started !'
                ! Inizialization of system tyme for the determination of total time elapsed
    	        call CPU_TIME(TSTART)
    	        write (*,*) '-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~'
                i=0
    	        do while (t.lt.tend)
        	    		! The program will write on file on each 1000 steps
15                      if((i.eq.0).or.(i.eq.wp)) then
                       		write (1,101) (t/(31536000.)),y1(2),y3(2),y4(2),y12(2),y14(2),y16(2)
                            i=0
                        endif
                        t = t + dt
                        i = i + 1
                        if(t>tend) t=tend
                       
                        ! Defining vectors of n+1 steps as the vectors of n step for the calculation of discrepancy vector

                        y1(1) = y1(2)
                        y3(1) = y3(2)
                        y4(1) = y4(2)
                        y12(1) = y12(2)
                        y14(1) = y14(2)
                        y16(1) = y16(2)
                       
                        ! Calculating the discrepancy vector
 
        	            E(1) = (y1(1)- y1(2))/(dt*rho) + ((3./2.)*(y1(1)**2.)*cs(1) + y3(1)*y4(1)*cs(3) - (y3(1)**2.)*cs(2) &
                        + (2.)*y1(1)*y12(1)*cs(4) + (2.)*y1(1)*y14(1)*cs(6) + (2.)*y1(1)*y14(1)*cs(5) &
                        +(2.)*y1(1)*y16(1)*cs(7))
                        Emax = E(1)
                       
                        E(2) = (y3(1) - y3(2))/(dt*rho) + ((y3(1)**2.)*cs(2) + y3(1)*y4(1)*cs(3) - 0.5*(y1(1)**2.)*cs(1))
                        if(E(2)>Emax) Emax = E(2)

                        E(3) = (y4(1) - y4(2))/(dt*rho) + (-(0.5)*(y3(1)**2)*cs(2) - y3(1)*y4(1)*cs(3) - y1(1)*y14(1)*cs(6)&
                        - y1(1)*y16(1)*cs(7))
                        if(E(3)>Emax) Emax = E(3)
                       
                        E(4) = (y12(1) - y12(2))/(dt*rho) + (y1(1)*y12(1)*cs(4) - y1(1)*y14(1)*cs(6))
                        if(E(4)>Emax) Emax = E(4)       
      
                        E(5) = (y14(1) - y14(2))/(dt*rho) - (y1(1)*y12(1)*cs(4) - y1(1)*y14(1)*cs(6) - y1(1)*y14(1)*cs(5) &
                        + y1(1)*y16(1)*cs(7)) 
                        if(E(5)>Emax) Emax = E(5)
                       
                        E(6) = (y16(1) - y16(2))/(dt*rho) - (- y1(1)*y16(1)*cs(7) + y1(1)*y14(1)*cs(5))
                        if(E(6)>Emax) Emax = E(6)

                        DYmax = 10*tollDY
                       
                        ! Calculating the coefficients of the discrepancy derivatives matrix
                        do while(Emax>tollE .and. DYmax>tollDY)
                            
                            call DDM (dt,rho,y1,y3,y4,y12,y14,y16,cs,A)                                                     

                            ! Calling the subroutine for the calculation of the increment value for the y vectors
              
                            call KERNEL (Nreac,A,MXapp,E,Nreac,check)
                            if(check == 1) then
                                write (*,*) 'Warning the application runned into a critical error!'
                                goto 999
                            endif

                            ! Creating the new vectors
                                
                            y1(1) = y1(1) + E(1)
                            if(y1(1) > 0.) DYmax = ABS(E(1)/y1(1))

                            y3(1) = y3(1) + E(2)
                            if(ABS(E(2)/y3(1))>DYmax .and. y3(1)>0.) DYmax = ABS(E(2)/y3(1))
 
                            y4(1) = y4(1) + E(3)
                            if(ABS(E(3)/y4(1))>DYmax .and. y4(1)>0.) DYmax = ABS(E(3)/y4(1))

                            y12(1) = y12(1) + E(4)
                            if(ABS(E(4)/y12(1))>DYmax .and. y12(1)>0.) DYmax = ABS(E(4)/y12(1))
               
                            y14(1) = y14(1) + E(5)
                            if(ABS(E(5)/y14(1))>DYmax .and. y14(1)>0.) DYmax = ABS(E(5)/y14(1))           

                            y16(1) = y16(1) + E(6)
                            if(ABS(E(6)/y16(1))>DYmax .and. y16(1)>0.) DYmax = ABS(E(6)/y16(1))
      
                            ! Calculating the new Error vector with the corrected abundances

                            E(1) = (y1(1)- y1(2))/(dt*rho) + ((3./2.)*(y1(1)**2.)*cs(1) + y3(1)*y4(1)*cs(3) - (y3(1)**2.)*cs(2) &
                            + (2.)*y1(1)*y12(1)*cs(4) + (2.)*y1(1)*y14(1)*cs(6) + (2.)*y1(1)*y14(1)*cs(5) &
                            +(2.)*y1(1)*y16(1)*cs(7))
                            Emax = E(1)
                       
                        	E(2) = (y3(1) - y3(2))/(dt*rho) + ((y3(1)**2.)*cs(2) + y3(1)*y4(1)*cs(3) - 0.5*(y1(1)**2.)*cs(1))
                        	if(E(2)>Emax) Emax = E(2)

                        	E(3) = (y4(1) - y4(2))/(dt*rho) + (-(0.5)*(y3(1)**2)*cs(2) - y3(1)*y4(1)*cs(3) - y1(1)*y14(1)*cs(6)&
                        	- y1(1)*y16(1)*cs(7))
                        	if(E(3)>Emax) Emax = E(3)
                       
                        	E(4) = (y12(1) - y12(2))/(dt*rho) + (y1(1)*y12(1)*cs(4) - y1(1)*y14(1)*cs(6))
                        	if(E(4)>Emax) Emax = E(4)       
      
                        	E(5) = (y14(1) - y14(2))/(dt*rho) - (y1(1)*y12(1)*cs(4) - y1(1)*y14(1)*cs(6) - y1(1)*y14(1)*cs(5) &
                        	+ y1(1)*y16(1)*cs(7)) 
                        	if(E(5)>Emax) Emax = E(5)
                       
                        	E(6) = (y16(1) - y16(2))/(dt*rho) - (- y1(1)*y16(1)*cs(7) + y1(1)*y14(1)*cs(5))
                        	if(E(6)>Emax) Emax = E(6)
                       end do          
    	         end do

            ! Run untill the hydrogen abundance has reached a value y1min 
            case (2) 

                ! Setting the minimum value of y1
16              write (*,*) '-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~'
                write (*,*) 'Please insert the minimum abbundance y1 that has to be reached in the simulation'
                read (*,*) tend
    	        if((tend<=0.).or.tend>1.) then
		       		write (*,*) 'Warning! y1 ending abundance must be in [0,1] range!'
        	        goto 16
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

        	        ! Defining vectors of n+1 steps as the vectors of n step for the calculation of discrepancy vector

                    y1(1) = y1(2)
                    y3(1) = y3(2)
                    y4(1) = y4(2)
                    y12(1) = y12(2)
                    y14(1) = y14(2)
                    y16(1) = y16(2)
                       
                    ! Calculating the discrepancy vector
 
        	        E(1) = (y1(1)- y1(2))/(dt*rho) + ((3./2.)*(y1(1)**2.)*cs(1) + y3(1)*y4(1)*cs(3) - (y3(1)**2.)*cs(2) &
                    + (2.)*y1(1)*y12(1)*cs(4) + (2.)*y1(1)*y14(1)*cs(6) + (2.)*y1(1)*y14(1)*cs(5) &
                    +(2.)*y1(1)*y16(1)*cs(7))
                    Emax = E(1)
                       
                    E(2) = (y3(1) - y3(2))/(dt*rho) + ((y3(1)**2.)*cs(2) + y3(1)*y4(1)*cs(3) - 0.5*(y1(1)**2.)*cs(1))
                    if(E(2)>Emax) Emax = E(2)

                    E(3) = (y4(1) - y4(2))/(dt*rho) + (-(0.5)*(y3(1)**2)*cs(2) - y3(1)*y4(1)*cs(3) - y1(1)*y14(1)*cs(6)&
                    - y1(1)*y16(1)*cs(7))
                    if(E(3)>Emax) Emax = E(3)
                       
                    E(4) = (y12(1) - y12(2))/(dt*rho) + (y1(1)*y12(1)*cs(4) - y1(1)*y14(1)*cs(6))
                    if(E(4)>Emax) Emax = E(4)       
      
                    E(5) = (y14(1) - y14(2))/(dt*rho) - (y1(1)*y12(1)*cs(4) - y1(1)*y14(1)*cs(6) - y1(1)*y14(1)*cs(5) &
                    + y1(1)*y16(1)*cs(7)) 
                    if(E(5)>Emax) Emax = E(5)
                       
                    E(6) = (y16(1) - y16(2))/(dt*rho) - (- y1(1)*y16(1)*cs(7) + y1(1)*y14(1)*cs(5))
                    if(E(6)>Emax) Emax = E(6)

                    DYmax = 10*tollDY
                       
                    ! Calculating the coefficients of the discrepancy derivatives matrix
                    do while(Emax>tollE .and. DYmax>tollDY)
                                
                            call DDM (dt,rho,y1,y3,y4,y12,y14,y16,cs,A)

                            ! Calling the subroutine for the calculation of the increment value for the y vectors
              
                            call KERNEL (Nreac,A,MXapp,E,Nreac,check)
                            if(check == 1) then
                                write (*,*) 'Warning the application runned into a critical error!'
                                goto 999
                            endif

                            ! Creating the new vectors
                                
                            y1(1) = y1(1) + E(1)
                            if(y1(1) > 0.) DYmax = ABS(E(1)/y1(1))

                            y3(1) = y3(1) + E(2)
                            if(ABS(E(2)/y3(1))>DYmax .and. y3(1)>0.) DYmax = ABS(E(2)/y3(1))
 
                            y4(1) = y4(1) + E(3)
                            if(ABS(E(3)/y4(1))>DYmax .and. y4(1)>0.) DYmax = ABS(E(3)/y4(1))

                            y12(1) = y12(1) + E(4)
                            if(ABS(E(4)/y12(1))>DYmax .and. y12(1)>0.) DYmax = ABS(E(4)/y12(1))
               
                            y14(1) = y14(1) + E(5)
                            if(ABS(E(5)/y14(1))>DYmax .and. y14(1)>0.) DYmax = ABS(E(5)/y14(1))           

                            y16(1) = y16(1) + E(6)
                            if(ABS(E(6)/y16(1))>DYmax .and. y16(1)>0.) DYmax = ABS(E(6)/y16(1))
      
                            ! Calculating the new Error vector with the corrected abundances

                            E(1) = (y1(1)- y1(2))/(dt*rho) + ((3./2.)*(y1(1)**2.)*cs(1) + y3(1)*y4(1)*cs(3) - (y3(1)**2.)*cs(2) &
                        	+ (2.)*y1(1)*y12(1)*cs(4) + (2.)*y1(1)*y14(1)*cs(6) + (2.)*y1(1)*y14(1)*cs(5) &
                        	+(2.)*y1(1)*y16(1)*cs(7))
                        	Emax = E(1)
                       
                        	E(2) = (y3(1) - y3(2))/(dt*rho) + ((y3(1)**2.)*cs(2) + y3(1)*y4(1)*cs(3) - 0.5*(y1(1)**2.)*cs(1))
                        	if(E(2)>Emax) Emax = E(2)

                        	E(3) = (y4(1) - y4(2))/(dt*rho) + (-(0.5)*(y3(1)**2)*cs(2) - y3(1)*y4(1)*cs(3) - y1(1)*y14(1)*cs(6)&
                        	- y1(1)*y16(1)*cs(7))
                        	if(E(3)>Emax) Emax = E(3)
                       
                        	E(4) = (y12(1) - y12(2))/(dt*rho) + (y1(1)*y12(1)*cs(4) - y1(1)*y14(1)*cs(6))
                       		if(E(4)>Emax) Emax = E(4)       
      
                        	E(5) = (y14(1) - y14(2))/(dt*rho) - (y1(1)*y12(1)*cs(4) - y1(1)*y14(1)*cs(6) - y1(1)*y14(1)*cs(5) &
                        	+ y1(1)*y16(1)*cs(7)) 
                        	if(E(5)>Emax) Emax = E(5)
                       
                        	E(6) = (y16(1) - y16(2))/(dt*rho) - (- y1(1)*y16(1)*cs(7) + y1(1)*y14(1)*cs(5))
                        	if(E(6)>Emax) Emax = E(6)
                       end do          
                       
    	        end do
                write (*,*) 'y1 has reached ',y1(2),' abundance after ',(t/31536000.),' years!'          
          
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
    end program NuclearMainSequenceReactions_Implicit

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
    
    !Subroutine for the calculation of the Discrepancy Derivative Matrix (DDM) given the abundance vectors for each element
  
    subroutine DDM (dt,rho,y1,y3,y4,y12,y14,y16,cs,A)
        real (kind=8) :: dt,rho

        real (kind=8), dimension (2) :: y1,y3,y4,y12,y14,y16

        real (kind=8), dimension (7) :: cs  

        real (kind=8), dimension (6,6) :: A

                                                 !DE1/Dy1!
                                 
                       		A(1,1) = (1. - y1(2))/(dt*rho) + (3.)*y1(1)*cs(1) + (2.)*(y12(1)*cs(4) + y14(1)*cs(6) + y14(1)*cs(5) &
                            + y16(1)*cs(7))  
                      
                                                 !DE1/Dy3!

                      		A(1,2) = y4(1)*cs(3) - (2.)*y3(1)*cs(2)

                                                 !DE1/Dy4!

                       		A(1,3) = y3(1)*cs(3)

                                                 !DE1/Dy12!
 
                       		A(1,4) = (2.)*y1(1)*cs(4)

                                                 !DE1/Dy14!

                       		A(1,5) = (2.)*(y1(1)*cs(6) + y1(1)*cs(5))

                                                 !DE1/Dy16!

                       		A(1,6) = (2.)*y1(1)*cs(7)
 
                                                 !DE2/Dy1!

                       		A(2,1) = - y1(1)*cs(1)
      
                                                 !DE2/Dy3!

                       		A(2,2) = (1. - y3(2))/(dt*rho) + (2.)*y3(1)*cs(2) + y4(1)*cs(3)

                                                 !DE2/Dy4!

                       		A(2,3) = y3(1)*cs(3)
                     
                                                 !DE2/Dy12!

                      		A(2,4) = 0.

                                                 !DE2/Dy14!

                       		A(2,5) = 0.

                                                 !DE2/Dy16!

                       		A(2,6) = 0.

                                                 !DE3/Dy1!

                       		A(3,1) = - y14(1)*cs(6) - y16(1)*cs(7)

                                                 !DE3/Dy3!

                       		A(3,2) = - y3(1)*cs(2) - y4(1)*cs(3)

                                                 !DE3/Dy4!

                       		A(3,3) = (1. - y4(2))/(dt*rho) - y3(1)*cs(3)

                                                 !DE3/Dy12!

                       		A(3,4) = 0.

                                                 !DE3/Dy14!

                       		A(3,5) = - y1(1)*cs(6)

                                                 !DE3/Dy16!

                       		A(3,6) = - y1(1)*cs(7)
 
                                                 !DE4/Dy1!

                       		A(4,1) = y12(1)*cs(4) - y14(1)*cs(6)

                                                 !DE4/Dy3!

                       		A(4,2) = 0.
 
                                                 !DE4/Dy4!

                       		A(4,3) = 0.

                                                 !DE4/Dy12!

                       		A(4,4) = (1. - y12(2))/(dt*rho) + y1(1)*cs(4)

                                                 !DE4/Dy14!

                       		A(4,5) = -y1(1)*cs(6)
 
                                                 !DE4/Dy16!

                       		A(4,6) = 0.

                                                 !DE5/Dy1!

                       		A(5,1) = -y12(1)*cs(4) + y14(1)*cs(6) + y14(1)*cs(5) -y16(1)*cs(7)
      
                                                 !DE5/Dy3!

                       		A(5,2) = 0.

                                                 !DE5/Dy4!

                       		A(5,3) = 0.

                                                 !DE5/Dy12!

                       		A(5,4) = -y1(1)*cs(4)

                                                 !DE5/Dy14!

                       		A(5,5) = (1. - y14(1))/(dt*rho) + y1(1)*cs(6) + y1(1)*cs(5)

                                                 !DE5/Dy16!

                       		A(5,6) = -y1(1)*cs(7)
 
                                                 !DE6/Dy1!

                       		A(6,1) = y16(1)*cs(7) - y14(1)*cs(5)

                                                 !DE6/Dy3!

                       		A(6,2) = 0.
  
                                                 !DE6/Dy4!

                       		A(6,3) = 0.

                                                 !DE6/Dy12!

                       		A(6,4) = 0.

                                                 !DE6/Dy14!

                       		A(6,5) = - y1(1)*cs(5)

                                                 !DE6/Dy16!

                       		A(6,6) = (1. - y16(2))/(dt*rho) + y1(1)*cs(7)       
    	return
    	end

    ! Subroutine for the calculation of increment vector given a DDM
    ! MAXNE = numero totale di isotopi
    ! DXX = matrice dei coefficienti del sistema, dimensionata a (maxne,maxne)
    ! A = matrice allocata nel programma chiamante dimensionata a (maxne*maxne)
    ! B = vettore dei termini noti dimensionato a (maxne) in entrata - soluzioni del sistema in uscita
    ! N = numero di isotopi considerati nel network nel vostro caso N=MAXNE
    ! KS = flag di controllo KS=0 se operazione andata abuon fine, KS=1 in caso di errori

    SUBROUTINE KERNEL (MAXNE,DXX,A,B,N,KS)
        
        implicit double precision (A-H,O-Z)
        DIMENSION DXX(maxne,maxne),B(maxne),A(maxne*maxne)
        DO IE=1,N
           DO II=1,N
               IA=(IE-1)*N+II
               A(IA)=DXX(II,IE)
           END DO
        END DO
        TOL=0.D0
        KS=0
        JJ=-N
      DO 65 J=1,N
      JY=J+1
      JJ=JJ+N+1
      BIGA=0.D0
      IT=JJ-J
      DO 30 I=J,N
      IJ=IT+I
      IF(DABS(BIGA)-DABS(A(IJ))) 20,30,30
   20 BIGA=A(IJ)
      IMAX=I
   30 CONTINUE
      IF(DABS(BIGA)-TOL) 35,35,40
   35 KS=1
      RETURN
   40 I1=J+N*(J-2)
      IT=IMAX-J
      DO 50 K=J,N
      I1=I1+N
      I2=I1+IT
      SAVE=A(I1)
      A(I1)=A(I2)
      A(I2)=SAVE
   50 A(I1)=A(I1)/BIGA
      SAVE=B(IMAX)
      B(IMAX)=B(J)
      B(J)=SAVE/BIGA
      IF(J-N) 55,70,55
   55 IQS=N*(J-1)
      DO 65 IX=JY,N
      IXJ=IQS+IX
      IT=J-IX
      DO 60 JX=JY,N
      IXJX=N*(JX-1)+IX
      JJX=IXJX+IT
   60 A(IXJX)=A(IXJX)-(A(IXJ)*A(JJX))
   65 B(IX)=B(IX)-(B(J)*A(IXJ))
   70 NY=N-1
      IT=N*N
      DO 80 J=1,NY
      IA=IT-J
      IB=N-J
      IC=N
      DO 80 K=1,J
      B(IB)=B(IB)-A(IA)*B(IC)
      IA=IA-N
   80 IC=IC-1
      RETURN
      END