include 'testnum.fi'
program CurWork2
    use testN
    implicit none
    integer :: N = 4, nn = 3
    real :: R0 = 1.0, R1 = 3.0, L = 2.0

    test = 2
    capa = 5
    open (1, file="output.txt" )
    do while (N <= 128)
        write (1,*) "N= ", N
        call Simulation(N,N,R0,R1,L,test)
        N = 2**nn
        nn= nn+1
    enddo
end program CurWork2


!-----------------------------------------------------------!
!                       Subroutines                         !
!-----------------------------------------------------------!
    Subroutine Simulation(Nr,Nz,R0,R1,L,test)
        integer Nr, Nz
        real R0,R1,L
        real Hr,Hz
        real r(0:Nr), z(0:Nz)
        real r12(1:Nr), z12(1:Nz)   

        real C_A(0:Nr-1)   !  main diag matrix C
        real C_B(0:Nr-2)   !  naddiag matrix C
        real C_C(1:Nr-1)     !  poddiag matrix C

        real F(0:Nr-1,0:Nz)   ! right side of equation
        real V(0:Nr-1,0:Nz)   ! solutions
        real U(0:Nr-1,0:Nz)   ! exact solution
        real p(0:Nr-1,0:Nz)   ! pj = Fj in method of complete reduction
        real :: EpsMax = 0.0 
        integer test
!----------------------STEP---------------------------
        real NNr,NNz

        NNr = Nr
        NNz = Nz
        Hr = (R1 - R0)/NNr
        Hz = L/NNz
!----------------------GRID---------------------------
        r(0) = R0
        do i = 1, Nr
            r(i) = r(0) + i*Hr
            r12(i) = r(i-1) + Hr/2
        enddo 
        

        z(0) = 0 
        do j = 1, Nz
            z(j) = z(0)+ j*Hz
            z12(j) = z(j-1) + Hz/2
        enddo

	print *, 'start simulation'

        1 format (100es17.5)  
        2 format (7a7)
        open (1, file="output.txt" )
        
        call FindU(Nr,Nz,r,z,U)
	print *, 'U has founded'
        call FillMatrixC(Nr,C_A,C_B,C_C,Hr,Hz,r,r12)
	print *, 'matrix C intialized'
        call InitStartV(Nr,Nz,V,r,L)
	print *, 'V initialized'
        call InitF(R0,R1,L,Nr,Nz,Hr,Hz,r,r12,z,F)
        
	print *, 'F initialized'

        call FindNev(Nr,Nz,U,F,C_A,C_B,C_C)

	print *, 'nev has founded'	

        call ForwardWay(Nr,Nz,C_A,C_B,C_C,F,p) 

	print *, 'fw is finished'
 
        call ReverseWay(Nr,Nz,C_A,C_B,C_C,V,F,p)

	print *, 'rw is finished'

        !write (1,2)"V="
        !call PrintMatrix(V,Nr,Nz)
        call CalcEps(Nr,Nz,r,z,U,V,EpsMax)
  	
	print *, 'eps is counte'

        write (1,*)"EpsMax="
        write (1,1) EpsMax
    End subroutine

    Subroutine FillMatrixC(Nr,A,B,C,Hr,Hz,r,r12)
        integer Nr
        real A(0:Nr-1)    
        real B(0:Nr-2)
        real C(1:Nr-1)
        real Hr,Hz
        real r(0:Nr),r12(1:Nr)
        A = 0
        B = 0
        C = 0

	A(0) = (2*(Hz**2)/((Hr**2)*r(i))) * (r12(1)*funcK(r12(1))+capa) +2
	B(0) = -(2*(Hz**2)/((Hr**2)*r(i))) * r12(1)*funcK(r12(1))
        do i = 1 , Nr-1
            A(i)= ((Hz**2)/(Hr**2))*( (r12(i+1)/r(i))*funcK(r12(i+1)) + (r12(i)/r(i))*funcK(r12(i)) ) + 2.0 
        enddo
        do i = 1, Nr-2
            B(i) = -((Hz**2)/(Hr**2))* ( (r12(i+1)/r(i))*funcK(r12(i+1)) )                    
	enddo
        do i = 1, Nr-1
            C(i) = -((Hz**2)/((Hr**2)*r(i)))* (r12(i)*funcK(r12(i)) ) 
        enddo
    End subroutine
    
    Subroutine InitStartV(Nr,Nz,V,r,L)
        integer Nr, Nz
        real L
        real V(Nr-1,0:Nz)
        real r(0:Nr)
        
        V = 0
        do i = 0, Nr-1
            V(i,0)= fi3(r(i))
            V(i,Nz) = fi4(r(i),L)
        enddo
    End subroutine

    Subroutine InitF(R0,R1,L,Nr,Nz,Hr,Hz,r,r12,z,F)
        real R0,R1,L
        integer Nr, Nz
        real Hr,Hz
        real r(0:Nr),r12(1:Nr),z(0:Nz)
        real F(Nr-1,0:Nz)
        F = 0
        do j = 1, Nz-1      
            F(Nr-1,j) = (Hz**2)*funcF(r(Nr-1),z(j)) + ( (Hz**2)/(Hr**2) )*(r12(Nr)/r(Nr-1))*funcK(r12(Nr))*fi2(R1,z(j))
        enddo
        
        do i = 0, Nr-2
            do j = 1, Nz-1
                F(i,j) = (Hz**2)*funcF(r(i),z(j))
            enddo
        enddo

        do i = 0, Nr-1
            F(i,0) = fi3(r(i))
            F(i,Nz) = fi4(r(i),L)
        enddo
    End subroutine

    Subroutine SweepMethod(N,A,B,C,V,F)
        integer N
        real A(N-1) !   main diag 
        real B(N-2) !   naddiag 
        real C(2:N-1)!  poddiag 
        real V(N-1), F(N-1)

        real Alpha(2:N-1), Betta(2:N-1)
        ! forward process
        Alpha(2) = - B(1)/A(1)
        Betta(2) = F(1)/A(1)
        do i = 2,N-2
            if(abs(Alpha(i))>=1) then
                print *, "Alpha >=1  (badly)"
            else
                Alpha(i+1) = -B(i)/(C(i)*Alpha(i)+A(i))
                Betta(i+1) = (F(i) - C(i)*Betta(i))/(C(i)*Alpha(i)+A(i))
            endif
        enddo
        V(N-1) = (F(N-1) - C(N-1)*Betta(N-1))/(C(N-1)*Alpha(N-1)+A(N-1))
        ! reverce process
        do i = N-2,1,-1
            V(i) = Alpha(i+1)*V(i+1)+Betta(i+1)
        enddo
    End subroutine
    
    Subroutine ForwardWay(Nr,Nz,a,b,c,F,p)
        integer Nz,Nr
        real a(Nr-1)     !  main diag matrix C
        real b(Nr-2)     !  naddiag matrix C
        real c(2:Nr-1)   !  poddiag matrix C
        real F(Nr-1,0:Nz)
        real p(Nr-1,0:Nz)
        real n_real
        integer n_count
        real S(Nr-1,0:Nz)
        real alpha
        real clk(Nr-1)
        integer :: index = 0, j = 0
        real buf_v(Nr-1), buf_f(Nr-1)
        p = 0
        S = 0
        n_real = ( log(float(Nz))/log(float(2)) )
        n_count = int(n_real)
        p = F
        do k = 1, n_count-1
            index = 1
            j = index*(2**k)
            do while(j<=(Nz-2**k))
                buf_v = 0
                buf_f = 0
                S = 0
                do l = 1, 2**(k-1)
                    alpha = FAlpha(l,k-1)
                    call FClk(l,k-1,a,clk,Nr)
                    buf_f = alpha*(  p( :,j-(2**(k-1)) ) + p( :,j+(2**(k-1)) ))
                    call SweepMethod(Nr,clk,b,c,buf_v,buf_f)
                    S(:,j) = S(:,j) + buf_v
                enddo
                p(:,j) = 0.5*(p(:,j) + S(:,j))
                index = index + 1
                j = index*(2**k)
            enddo !end while
        enddo

    End subroutine 

    Subroutine ReverseWay(Nr,Nz,a,b,c,V,F,p)
        integer Nz,Nr
        real a(Nr-1)     !  main diag matrix C
        real b(Nr-2)     !  naddiag matrix C
        real c(2:Nr-1)   !  poddiag matrix C
        real V(Nr-1,0:Nz)
        real F(Nr-1,0:Nz)
        real p(Nr-1,0:Nz)
        real n_real
        integer n_count
        real alpha
        real clk(Nr-1)
        integer index , j ,k
        real buf_v(Nr-1), buf_f(Nr-1)

        n_real = ( log(float(Nz))/log(float(2)) )
        n_count = int(n_real)
        do k = n_count,1,-1
            index = 1 
            j = index*(2**(k-1))
            do while(j<=(Nz-2**(k-1)))
                buf_v = 0
                buf_f = 0
                do l=1,(2**(k-1))
                    alpha = FAlpha(l,k-1)
                    call FClk(l,k-1,a,clk,Nr)
                    buf_f = p(:,j) + alpha*( V(:,j-2**(k-1)) + V(:,j+2**(k-1)) )
                    call SweepMethod(Nr,clk,b,c,buf_v,buf_f)
                    V(:,j) = V(:,j) + buf_v
                enddo
                index = index + 2
                j = index*(2**(k-1))
            enddo ! end while
        enddo
    End subroutine 

    Subroutine FClk(l,k,A,Clk,Nr)
        integer l,k,Nr
        real A(Nr-1)
        real Clk(Nr-1)
        real pi
            pi = 3.141592653
            Clk = A - 2*COS( ((2*l-1)*pi)/(2**(k+1)) )
    End subroutine

    Subroutine CalcEps(Nr,Nz,r,z,U,V,EpsMax)
        integer Nr, Nz
        real r(0:Nr), z(0:Nz)
        real V(Nr-1,0:Nz)
        real U(Nr-1,0:Nz)
        real Eps(Nr-1,0:Nz)
        EpsMax = 0.0
        do i = 0, Nr-1
            do j = 0, Nz
                Eps(i,j) = abs(U(i,j) - V(i,j))
                EpsMax = MAX(EpsMax,Eps(i,j))
            enddo
        enddo
    End subroutine

    Subroutine PrintMatrix(matrix,Nr,Nz)
        integer Nr,Nz
        real matrix(Nr-1,0:Nz)

        1 format (100es17.5)  
        2 format (7a7)
        open (1, file="output.txt" )

        do i=1,Nr-1
            write (1,1), (matrix(i,j), j=0,Nz)
        enddo
    End subroutine

    Subroutine PrintVectorNr(vector,Nr)
        integer Nr
        real vector(Nr-1)

        1 format (100f17.5)  
        2 format (7a7)
        open (1, file="output.txt" )

        do i=1,Nr-1
            write (1,1000), vector(i)
        enddo

        1000 format (10(F15.5))  
    End subroutine

    Subroutine FindU(Nr,Nz,r,z,U)
        integer Nr, Nz
        real r(0:Nr), z(0:Nz)
        real U(Nr-1,0:Nz)

        do i = 0,Nr-1
            do j = 0,Nz
                U(i,j) = funcU(r(i),z(j))
            enddo
        enddo
    End subroutine

    Subroutine FindNev(Nr,Nz,U,F,a,b,c)
        integer Nr,Nz
        real U(Nr-1,0:Nz)
        real F(Nr-1,0:Nz)
        real Nev(Nr-1,0:Nz)
        real a(Nr-1)     !  main diag matrix C
        real b(Nr-2)     !  naddiag matrix C
        real c(2:Nr-1)   !  poddiag matrix C
        real dopF(Nr-1)

        Nev(:,0) = F(:,0) - U(:,0)
        Nev(:,Nz) = F(:,Nz) - U(:,Nz)
        do j = 1, Nz-1
            call MatrixMultVector(Nr,Nz,a,b,c,U(:,j),dopF)
            Nev(:,j) = F(:,j)-( - U(:,j-1) - U(:,j+1) + dopF)
        enddo 
    End Subroutine
    
    Subroutine MatrixMultVector(Nr,Nz,a,b,c,V,F)
        integer Nr,Nz
        real a(Nr-1)     !  main diag matrix C
        real b(Nr-2)     !  naddiag matrix C
        real c(2:Nr-1)   !  poddiag matrix C
        real V(Nr-1)
        real F(Nr-1)

        F(1) = a(1)*V(1) + b(1)*V(2)
        F(Nr-1) = c(Nr-1)*V(Nr-2) + a(Nr-1)*V(Nr-1)
        
        do i = 2,Nr-2
            F(i) = c(i)*V(i-1) + a(i)*V(i)+b(i)*V(i+1)
        enddo 

    End subroutine
!-----------------------------------------------------------!
!                       Functions                           !
!-----------------------------------------------------------!
    Function funcU(r,z)
    use testN
        real r,z
        real funcU
        
        select case (test)
            case(0)
                funcU = 1.0
            case(1)
                funcU = 3*r+z
            case(2)
                funcU = exp(r) + z
            
            case default
                funcU = 1.0
        end select
    End function 
        
    Function funcF(r,z)
    use testN
        real r,z
        real funcF
         
        select case (test)
            case(0)
                funcF = 0.0
            case(1)
                funcF = -6.0
            case(2)
                funcF = -exp(r)*( 2 + r )
            case(3) 
                funcF = -14*(r**2) + 16*(z**2) - 24*(r**4) + 8*(r**2)*(z**2)
            case default
                funcF = 0.0
        end select
    End function

    Function funcK(r)
    use testN
        real r
        real funcK 

        select case (test)
            case(0)
                funcK = 1.0
            case(1)
                funcK = r
            case(2)
                funcK = r 
            case(3)
                funcK = 1 + r**2
            case default
                funcK = 1.0
        end select
    End function
    
    Function fi1(R0,z)
    use testN
        real R0,z
        real fi1
        
        select case (test)
            case(3)
                fi1 = 1 - z**2 - z**4
            case default
                fi1 = capa*funcU(R0,z) - 5; !-5 is fi1(z_i)
        end select
    End function

    Function fi2(R1,z)
    use testN
        real R1, z 
        real fi2
        select case (test)
            case(3)
                fi2 = 81 - 9*(z**2) - z**4
            case default
                fi2 = funcU(R1,z)
        end select
    End function
    
    Function fi3(r)
    use testN
        real r 
        real fi3
        select case (test)
            case(3)
                fi3 = r**4 
            case default
                fi3 = funcU(r,0.0)
        end select
    End function

    Function fi4(r,L)
    use testN
        real r ,L
        real fi4
        select case (test)
            case(3)
                fi4 = r**4 - 4*(r**2) -16 
            case default
                fi4 = funcU(r,L)
        end select
    End function

    Function FAlpha(ll,kk)
        integer ll,kk
        real FAlpha
        real pi
        real l,k
            pi = 3.141592653
            l=ll
            k=kk
            FAlpha = ( (-1)**(l+1) )*  SIN( ((2*l - 1 )*pi)/( 2**(k+1)) )/(2**k)
    End function
