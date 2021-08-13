program mc_main
!Rotina principal do código -> Simulação de Monte Carlo em Ensemble NPT utilizando potenciais anisotrópicos (SW)
        use global_variables              !Variaveis do programa
        use init_confani2                 !configuração inicial
        use gr                            !calcula do histograma para g(r)
        implicit none
        integer*8 :: i, j, k              !Contadores para os fors
        integer*8 :: ciclo                !contador de ciclos
        real*8    :: a                    !parametro de largura da celula unitaria
	    real*8    :: beta                 !1/kBT   
        real*8    :: v, vn, vm, dv, dvb   !potencial para n m, diferença de potencial e potencial * beta
	    real*8    :: dv_new               !delta de volume
	    real*8    :: deltah               !entalpia do sistema
        real*8    :: pot                  !energia potencial
        real*8    :: pressure             !pressao
        real*8    :: drmax                !variação maxima no deslocamento
        real*8    :: max_v, deltav        !variação maxima de volume
	    real*8    :: ratio_trans          !razao de aceite para translacao 
	    real*8    :: ratio_rot            !razao de aceite para rotacao
    	real*8    :: ratio_vol            !razão de aceite para volume
        real*8    :: max_ang_dummy        !variavel que armazena o valor de angulo maximo anterior
        real*8    :: random_n             !numero aleatorio gerado
        real*8    :: rxn, ryn, rzn        !posições do estado atual
        real*8    :: rxm, rym, rzm        !posições do novo estado
	    real*8    :: rnew(3,100000)       !posicoes apos novo volume
	    real*8    :: boxlnew(3)           !novo tamanho dos eixos da caixa
	    real*8    :: volume_new           !novo volume
        real*8    :: volume_dummy         !estoca a soma dos valores dos volumes
        real*8    :: volume_dummy2
        real*8    :: h_dummy              !estoca a soma dos valores de entalpia
        real*8    :: h_dummy2
        real*8    :: vh_dummy 
        real*8    :: rijsq                !distancia entre 2 moleculas ao quadrdo
        real*8    :: rij                  !distancia entre duas moleculas]
        real*8    :: g_r(100000)          !vetor para armezenar g(r)
	    real*8    :: quat(100000,4)       !Armazena os quaternions atuais
	    real*8    :: q_old(4)             !armazena a orientaçao antiga da particula
        real*8    :: max_ang              !angulo maximo para rotação
	    real*8    :: e1(3), e2(3)         !armazena as orientações
	    real*8    :: initang              !efixed for initial orientation and initital angle
	    real*8    :: e_old(3)             !armazena a orientação antiga da particula
	    real*8    :: e_new(3)             !armazena a orientação nova da particula
	    real*8    :: quat_new(4)          !novo quaternion
	    real*8    :: rxij, ryij, rzij     !armazena a diferença de duas posiçoes
	    real*8    :: urij(3)              !vetor unitario da distancia de 2 particulas
	    real*8    :: new_pot, potij, dpot !variaveis para o calculo da nova energia com mudança de volume
	    real*8    :: rxi, rxj
	    real*8    :: ryi, ryj
	    real*8    :: rzi, rzj
        real*8    :: cp 
        integer*8 :: max_ciclo            !numero mximo de ciclos por step
        integer*8 :: coord                !coordenada inicial de rotacao
        integer*8 :: seed                 !semente para numeros aletorios
        integer*8 :: ntrans, nrot, nvol   !contador para o numero de movimentos rotaçao translacao e volume
        integer*8 :: nacc_trans, nacc_rot, nacc_vol !contador para o numero de movimentos ACEITOS  de rotacao translacao e volume
        logical   :: overlap              !booleano para checar overlap
        character :: filename*40          !file name for outpout files
        character :: get*100              !variable for reading input file
        !write(*,*) "Use reduced units"
        !write(*,*) "Enter with number of particles, packing factor, reduced pressure, reduced temperature and elongation (K) "
        !read(*,*) n_particles, eta0, pressure, temp, ke
        !write(*,*) "Enter the initial rotation axis (1-x, 2-y, 3-z) and angle (radians)"
        !read(*,*) coord, initang
        !write(*,*) "Enter the number of steps"
        !read(*,*) max_steps
       
        !Input file --> see example input1.sci 
        open(101, File='input1.sci')
            read(101,*) get, n_particles
            read(101,*) get, eta0
            read(101,*) get, pressure
            read(101,*) get, temp 
            read(101,*) get, ke
            read(101,*) get, coord
            read(101,*) get, initang
            read(101,*) get, max_steps
        close(101)
        !Parametros gerais da simulação          								
        max_ciclo   = 2*n_particles                     !numero max de ciclos
        n_equil     = 0.7d0*max_steps                   !Steps equilibração
        n_prod      = max_steps-n_equil                 !Steps produçao
        n_save      = 1000                              !Bloco de steps que mostra rodando na tela
        n_adjust    = 1000                              !Bloco de steps para arrumar o drmax
        seed        = 202419                            !seed do numero aleatorio, meu RA
        mtoang      = 10.d0**(-10.d0)
        !Parametros para o SW e Monte Carlo
        lambda      = 1.626d0                           !tamanho do poço
        lambdasq    = lambda*lambda
        overlap     = .FALSE.                           !booleano de overlap, inicialmente falso
        g_r(:)      = 0.d0                              !vetor para calcular a g(r)
               
        !Parametros anisotrópicos
        sigs        = 2.144d0                                !sigma side by side --> reduzido
        sige        = ke*sigs                                !k = sigs / sige
        epslon      = 275.d0 
        !sigma check --> 21-07
        sigma       = ((sigs*sigs*sige)**(1.d0/3.d0))/sigs   !sigma sphere on a anisotropic potential (Joyce)	
        chi         = (ke*ke - 1.d0)/(ke*ke + 1.d0)          !anisotropia da particula    
        v_particle  = (pi*sigs*sigs*sige)/6.d0               !volume da particula, vai ser reduzido na init_conf
  
        !reference system
        efixed = (/0.d0, 0.d0, 1.d0/) !Z Axis

        !calcula densidade reduzida a partir do empacotamento
        rho_red = eta0/(v_particle/sigs**3.0d0)

        !Gera configuração inicial cubico face centrada
        call fcc()

        !Aloca o vetor de orientação  para todas as particulas baseado no angulo e eixo de rotação escolhidos
        call initial_orientation(coord, initang, quat(:,:))
      
        !printa alguns parametros iniciais
        write(*,*) "-----------------------Simulation Parameters--------------------------------"
        write(*,*) "Total steps", max_steps
        write(*,*) "Number Particles:", n_particles
        write(*,*) "Particle volume:", v_particle
        write(*,*) "Anisotropy:", ke
        write(*,*) "Initial Density:", rho_red
        write(*,*) "Packing Factor:", eta0
        write(*,*) "Pressure:", pressure
        write(*,*) "Initial box length:", box_length(1), box_length(2), box_length(3)
        write(*,*) "Initial Volume:", volume          
        write(*,*) "----------------------------------------------------------------------------"
          
        !Convert real variables to string, useful for file handling. 
        write(press_c, '(I3)') int((pressure*epslon*boltzmann)/((sigs**3.d0*(10.d0**(-10.d0))**3.d0))/1000000)
        write(temp_c, '(I3)')  int(temp*epslon)

        !Maximum displacement for initial moves
        drmax   = 0.05d0
        max_ang = 0.05d0
        max_v   = 0.001d0

        !Registro da configuração inicial 
        filename = 'conf_' // press_c // 'MPa_' // temp_c // 'K' // '.xyz'
        open(1, file = filename)
        write(1,*) n_particles ! "positions x, y, z", "orientation x, y, z"
        write(1,*) ''
        do i = 1, n_particles
                write(1,'(a1,10e15.7)') "C", rx(i), ry(i), rz(i), quat(i,:), sigs*0.5d0/sigs, sigs*0.5d0/sigs, sige*0.5d0/sigs
        end do

        !Find the longest axis
        sig_max = max(sigs/sigs,sige/sigs)

        !Calculo da energia inicial
        call compute_total_energy(v)

        write(*,*) "Initial state energy", v
        
        !Stochastic sampling using Metropolis

        !File for accept ratio 
        filename = 'ratio_' // press_c // '_MPa_' // temp_c // '_K' // '.dat'
        open(2,file= filename)
        write(2,*) "steps  ratio_tran | max_dr | ratio_rot | max_ang | ratio_vol | max_vol" 
        
        !File for potencial, volume and density storage   
        filename = 'potencial_' // press_c // '_MPa_' // temp_c // '_K' // '.dat'
        open(3, file = filename)
         
        !File for histogram and used for g(r)
        filename = 'hist_' // press_c // '_MPa_' // temp_c // '_K' // '.dat'
        open(4, file  = filename)

        write(*,*) "Initializing Monte Carlo"
       
        ntrans     = 0                  !total de movimentos de translacao
        nrot       = 0                  !total de movimentos de rotacao
        nvol       = 0                  !total de movimentos de volume

        nacc_trans = 0                  !contador de movimentos translacionais aceitos
        nacc_rot   = 0                  !contador de movimentos de rotação aceitos
        nacc_vol   = 0                  !contador de movimentos de volume aceitos

        volume_dummy = 0.d0             !armazena o volume 
        h_dummy      = 0.d0             !armazena a entalpia
        h_dummy2     = 0.d0
        vh_dummy     = 0.d0
        do steps = 1 , max_steps
                !contador na tela para acompanhar os passos
                if (mod(steps, n_save) == 0) then
                        write(*,*) steps, v, volume           
                end if
    
        do ciclo = 1, max_ciclo  
                
                call ranf(seed, random_n)
                random_n = idint(random_n*dble((n_particles + 1))) + 1

                if (random_n .le. dble(n_particles)) then
!----------------------------------------------------------------------------------------------------
!                           	        Translational Move                 		    	            !										
!----------------------------------------------------------------------------------------------------
                                !incremento no numero total de translacoes
                                ntrans = ntrans + 1                         
                                
                                !seleciona uma particula aleatotriamente
                                
                                !usar o random_n de cima 
                                call ranf(seed,random_n)
                                i = nint(random_n * dble(n_particles-1)) + 1
              
                                rxn      = rx(i)  
                                ryn      = ry(i)
                                rzn      = rz(i)
                                e_old(:) = e(i,:)
                                q_old(:) = quat(i,:)
                                !calcula energia do estado n
                                call compute_energy(i, rxn, ryn, rzn, vn, overlap, e_old)

                                !gera novo estado
                                call ranf(seed, random_n)
                                rxm = rxn+(2.d0*random_n-1.d0)*drmax
                                call ranf(seed, random_n)
                                rym = ryn+(2.d0*random_n-1.d0)*drmax
                                call ranf(seed, random_n)
                                rzm = rzn+(2.d0*random_n-1.d0)*drmax
            
                                !Imagem mínima
                                rxm = rxm - box_length(1)*dnint(rxm/box_length(1)) 
                                rym = rym - box_length(2)*dnint(rym/box_length(2))
                                rzm = rzm - box_length(3)*dnint(rzm/box_length(3))
               
                                !Mudança na orientação
                           
                                call rotation_quat(seed, max_ang, q_old, quat_new)
                                call rotation(efixed, e_new, quat_new) 
                                
                                call compute_energy(i, rxm, rym, rzm, vm, overlap, e_new)

                                if (.not. overlap) then
                                        dv         = vm - vn
                                        dvb        = dv/temp
                                        if (dvb .le. 0.d0) then
                                                rx(i)      = rxm
                                                ry(i)      = rym
                                                rz(i)      = rzm
                                                v          = v + dv
                                                e(i,:)     = e_new(:)
                                                quat(i,:)  = quat_new(:)
                                                nacc_trans = nacc_trans + 1
                                        else if (dvb .lt. 75.d0) then
                                                call ranf(seed, random_n)
                                                if (dexp(-dvb) .ge. random_n) then
                                                        rx(i)      = rxm
                                                        ry(i)      = rym
                                                        rz(i)      = rzm
                                                        v          = v + dv
                                                        e(i,:)     = e_new(:)
                                                        quat(i,:)  = quat_new(:)
                                                        nacc_trans = nacc_trans + 1
                                                end if 
                                        end if
                                end if 
!----------------------------------------------------------------------------------------------
! 				                        Volume Move                                       
!----------------------------------------------------------------------------------------------
                else
                        !Incrementa o numero de movimentos no volume
                        nvol  = nvol + 1
                
                        !zerar newpot
                        new_pot = 0.d0

                        !gera um novo volume
                        call new_volume(seed, rnew, boxlnew, volume_new, max_v)
                
                        deltav = volume_new - volume

                        !Novo potencial 
                        call new_potential (rnew, boxlnew, overlap, new_pot)

                        if (.not. overlap) then
                                dpot   = new_pot - v
                                deltah = pressure*deltav + dpot   !pdv/epslon + potencial/epslon
                 
                                deltah = deltah/temp - (dble(n_particles)+1.d0)*dlog(volume_new/volume)

                                if (deltah .le. 0.0d0) then
                                        rx(:)           = rnew(1,:)
                                        ry(:)           = rnew(2,:)
                                        rz(:)           = rnew(3,:)
                                        box_length(:)   = boxlnew(:)
                                        volume          = box_length(1) * box_length(2) * box_length(3)
                                        rho_red         = dble(n_particles)/volume
                                        v               = new_pot
                                        nacc_vol        = nacc_vol + 1

                                elseif (deltah .lt. 75.d0) then
                                        call ranf(seed, random_n)
                                        if (dexp(-deltah) .gt. random_n) then
                                                rx(:)               = rnew(1,:)
                                                ry(:)               = rnew(2,:)
                                                rz(:)               = rnew(3,:)
                                                box_length(:)       = boxlnew(:)
                                                volume              = box_length(1) * box_length(2) * box_length(3)
                                                rho_red             = dble(n_particles)/volume
                                                v                   = new_pot
                                                nacc_vol            = nacc_vol + 1
                                       end if 
                                end if
                        end if
                end if  !fim do movimento aleatorio
          end do        !fim do ciclo   
    
!----------------------------------------------------------------------------------------------
!                          Inicio do ajuste de propriedades de deslocamento                   !
!----------------------------------------------------------------------------------------------
    
                !Ajuste do deslocamento e angulos maximos 
                if (steps .le. n_equil) then

                        if (mod(steps, n_adjust) == 0) then                          !durante a equilibração apenas
                                ratio_trans = dble(nacc_trans)/dble(ntrans)
                                if (ratio_trans .lt. 0.5d0) then 
                                        drmax   = 0.95d0*drmax
                                        max_ang = 0.95d0*max_ang
                                else 
                                        drmax   = 1.05d0*drmax
                                        max_ang = 1.05d0*max_ang
                                end if
                                ntrans     = 0 
                                nacc_trans = 0
                        end if                         
                               
                        if (mod(steps, n_adjust) == 0) then
                                ratio_vol = dble(nacc_vol)/dble(nvol)
                                if (ratio_vol .lt. 0.5d0) then
                                        max_v = 0.95d0*max_v
                                else 
                                        max_v = 1.05d0*max_v
                                end if
                                nvol       = 0
                                nacc_vol   = 0
                        end if
           
    
                        if (mod(steps, n_adjust) == 0) then
                            write(2,*) steps, ratio_trans, drmax, ratio_rot, max_ang, ratio_vol,  max_v   !salva os dados de razao de aceite, deslocamento maximo, angulo de rotacao maximo e variacao de volume maxima	
                        end if
                end if
!-------------------------------------------------------------------------------------------------------------!
!                                          Preparação de Arquivos e g(r)  
!-------------------------------------------------------------------------------------------------------------!
                 if (mod(steps,n_save) == 0) then 
                        !arquivo para salvar as posições das particulas
                        write(1,*) n_particles
                        write(1,*) ''   
                        do i = 1 , n_particles
                                write(1, '(a1,10e15.7)') "C",  rx(i), ry(i), rz(i), quat(i,:), sigs*0.5d0/sigs, sigs*0.5d0/sigs&
                                &,sige*0.5d0/sigs
                        end do

                        !calcular a g(r)
                        if (steps .ge. n_equil) then
                                call rdf()
                                g_r(:) = g_r(:) + hist(:)
                        end if                   
                end if
        
                if (steps .ge. n_equil) then
                        volume_dummy  = volume_dummy  + volume
                        volume_dummy2 = volume_dummy2 + volume**2.d0
                        h_dummy       = h_dummy  + v
                        h_dummy2      = h_dummy2 + v**2.d0 
                        vh_dummy      = vh_dummy + volume*v
                end if
                
                if (mod(steps,n_save) == 0) then 
                    write(3,*) steps, v, volume, rho_red 
                end if
        end do !fim dos  steps
        
!----------------------------------------------------------------------------------------
! 		              Fim do metodo de Monte Carlo, normalização da g(r)                
!----------------------------------------------------------------------------------------
        !normalização da g(r)
        g_r(:) = g_r(:)*dble(n_save)/dble(n_prod)
         
        do i = 1, nbin
                write(4,*) i*deltar, g_r(i)
        end do 
!---------------------------------------------------------------------------------------
!                           Properties calculation and file closing
!---------------------------------------------------------------------------------------
        call calculate_properties(h_dummy, h_dummy2, volume_dummy, volume_dummy2, vh_dummy)

        !Fecha todos os arquivos com dados 
        close(1)  !conf.xyz
        close(2)  !ratio.dat
        close(3)  !potencial.dat
        close(4)  !histograma.dat
!--------------------------------------------------------------------------------------
!                       Monte Carlo End --> Printing final parameters
!--------------------------------------------------------------------------------------

        write(*,*) "-------------- Parâmetros finais ------------------"
        write(*,*) "Volume Final Reduzido:", volume_dummy/n_prod
        write(*,*) "Entalpia Final Reduzida:", h_dummy/n_prod
        write(*,*) "Box Length", box_length(:)
        write(*,*) "Pressão:", pressure
        write(*,*) "Temperatura:", temp 
        write(*,*) "Sigma_s e Sigma_e", sigs, sige 
        write(*,*) "---------------------------------------------------"
end program mc_main

subroutine initial_orientation(coord, initang, quat)
!Essa subrotina gera as orientações iniciais e as armazena no vetor e(:,3) e quat(:,4)
        use global_variables
        implicit none
        real*8                :: random_n
        integer*8             :: n
        real*8                :: initang
	    real*8                :: initaxis(3)
	    real*8                :: e_new(3)
	    real*8                :: quat(10000, 4)
	    real*8                :: q_old(4)
        integer*8             :: coord 
        
        quat(:,:)  = 0.d0 
        if (coord == 1) then
                initaxis = (/1.d0, 0.d0, 0.d0/)
        else if (coord == 2) then
                initaxis = (/0.d0, 1.d0, 0.d0/)
        else if (coord == 3) then 
                initaxis = (/0.d0, 0.d0, 1.d0/)
        else
                write(*,*) "Invalid axis index please choose between 1-x, 2-y or 3-z"
                STOP
        end if 

        do n = 1, n_particles 
                quat(n,1) = dcos(initang*0.5)             !real (w)
                quat(n,2) = initaxis(1)*dsin(initang*0.5) !imaginario (x)
                quat(n,3) = initaxis(2)*dsin(initang*0.5) !imaginario (y)
                quat(n,4) = initaxis(3)*dsin(initang*0.5) !imaginario (z)
                q_old(:)  = quat(n,:)

                !generate orientation based on initial quat
                call rotation(efixed, e_new, q_old)
                e(n,:)    = e_new(:) 
       end do
end subroutine initial_orientation

subroutine compute_total_energy (v)
!This subroutine iterates all particles on a double loop to compute the initial state energy and check if the initial configuration
!is possibles, i.e no overlaps detected before the stochastic sampling begin. If an overlap is detect the WHOLE program is stopped
!and it must be repeated
        use global_variables
        implicit none
        integer*8 :: i                                      !Index for counter and particle orientations                           
        integer*8 :: j                                      !"
        real*8    :: rxi,rxj                                !Particle position holder
        real*8    :: ryi,ryj                                !"
        real*8    :: rzi,rzj                                !"
        real*8    :: rxij                                   !Molecules distances for each axis
        real*8    :: ryij                                   !"
        real*8    :: rzij                                   !"
        real*8    :: rijsq, rij                             !Distance module and distance module **2
        real*8    :: vij                                    !Pairwise interacton potential
	    real*8    :: e1(3), e2(3)                           !Orientation vector for each molecule
        real*8    :: v                                      !Initial total energy 
	    real*8    :: urij(3)                                !Distance between molecules in a unit vector
        logical   :: overlap                                !Logical variable for overlap
       
        
        !Initialization
        v       = 0.d0    
        overlap = .FALSE.

        do i= 1, n_particles-1
                rxi    = rx(i)
                ryi    = ry(i)
                rzi    = rz(i)
                e1(:)  = e(i,:)
                do j=i+1,n_particles
                        e2(:) = e(j,:)
                        rxj   = rx(j)
                        ryj   = ry(j)
                        rzj   = rz(j)

                        !distance
                        rxij  = rxi-rxj
                        ryij  = ryi-ryj
                        rzij  = rzi-rzj

                        ! Minimum image convention
                        rxij  = rxij-box_length(1)*anint(rxij/box_length(1))
                        ryij  = ryij-box_length(2)*anint(ryij/box_length(2))
                        rzij  = rzij-box_length(3)*anint(rzij/box_length(3))
                       
                        rijsq = rxij*rxij+ryij*ryij+rzij*rzij !distance between two atoms
                        rij   = dsqrt(rijsq)

                        !versor for HGO
                        urij(1)  = rxij/rij
                        urij(2)  = ryij/rij
                        urij(3)  = rzij/rij
    
                        ! Cut-off radius
                        call check_overlap(e1, e2, rij, overlap, urij)
                        if (overlap) then
                                print*, "Overlap detected on initial configuration"
                                STOP
                        else
                                call compute_potential(rij, vij)
                                v =  v + vij
                        end if
                       
                end do
        end do

end subroutine compute_total_energy

subroutine compute_energy(i, rxi, ryi, rzi, v_part, overlap, e1)
!This subroutine computes the energy when a trial move is made. The loop iterates the moved particle all over the others to see if
!the new configuration is possible, i.e no overlaps between molecules (check overlap), if a overlap is detected the subroutine
!returns to the main program with overlap flag = TRUE. If no overlap is detected the subroutine continues to compute the potential
!energy between these two molecules. Then, the total energy of the trial system is returned to main program.
        use global_variables
        implicit none
        integer*8 :: i                                     !Index for counter and particle orientations
        integer*8 :: j                                     ! "
        real*8    :: rxi,rxj                               !Particle position holder
        real*8    :: ryi,ryj                               ! "
        real*8    :: rzi,rzj                               ! "
        real*8    :: rxij                                  !Molecules distances for each axis
        real*8    :: ryij                                  ! "
        real*8    :: rzij                                  ! "
        real*8    :: rijsq, rij                            !Distance module and distance module **2
        real*8    :: vij                                   !Pairwise interaction potential
        real*8    :: v_part                                !Total potential for new system
	    real*8    :: e1(3), e2(3)                          !Orientation vector for each molecule
	    real*8    :: urij(3)                               !Distance between molecules in a unit vector
        logical   :: overlap                               !Logical variable for overlap
        !Initialization
        v_part  = 0.d0
        overlap = .FALSE.
        
        !Iterates the "i"-particle (moved and rotate) over all the others
        do j = 1, n_particles
            if (j .ne. i) then
                        !Holds the orientation and position for j-particle
                        e2(:) = e(j,:)
                        rxj   = rx(j)
                        ryj   = ry(j)
                        rzj   = rz(j)
                
                        !Distance between i-j
                        rxij  = rxi-rxj
                        ryij  = ryi-ryj
                        rzij  = rzi-rzj

                        !Minimum image convention
                        rxij  = rxij-box_length(1)*anint(rxij/box_length(1))
                        ryij  = ryij-box_length(2)*anint(ryij/box_length(2))
                        rzij  = rzij-box_length(3)*anint(rzij/box_length(3))
  
                        !Distance between i-j
                        rijsq = rxij*rxij+ryij*ryij+rzij*rzij
                        rij   = dsqrt(rijsq)

                        !Unit versor
                        urij(1)   = rxij/rij
                        urij(2)   = ryij/rij
                        urij(3)   = rzij/rij
                        
                        !Checks if the molecules  overlap, if true -> return immediately to main and the move is denied, if false
                        !--> the potential for them is calculated
                        call check_overlap(e1, e2, rij, overlap,  urij)
                        
                        if (overlap) then
                                return 
                        else
                                call compute_potential(rij, vij)
                                v_part = v_part + vij
                        end if
            end if
        end do 
end subroutine compute_energy   

subroutine ranf(seed, random_n)
!This subroutine is a pseudo-random number generator, direct from Allen and Tildesley Appendix
        implicit none
        integer*8, parameter  :: l = 1029
        integer*8, parameter  :: c = 221591
        integer*8, parameter  :: m = 1048576
        integer*8             :: seed
        real*8                :: random_n
        !Generates a random number based on a seed 
        seed     = mod((seed*l+c),m)
        random_n = real(seed)/real(m)

end subroutine ranf

subroutine check_overlap(e1, e2, rij, overlap, urij)
!This subroutine computes if there is overlap between two molecules based on their contact distance (sig_hgo), computed on a
!subroutine below, and on their maximum overlap distance (sig_max) and returns the logical variable OVERLAP to compute_energy,
!compute_total_energy and new_potential subroutines.
        use global_variables
        implicit none
	    real*8    :: sig_hgo                   !HGO Contact distance between two molecules
	    real*8    :: e1(3), e2(3)              !Orientation vector for each molecule
	    real*8    :: rij                       !Distance between two molecules
	    real*8    :: urij(3)                   !unit versor between two centers of mass
        logical   :: overlap                   !Flag variables, TRUE if overlap is detected and immediately exits the subroutine
        
        !Intialization of contact distance
        sig_hgo = 0.d0  
        
        !Compares the distancee to the maximum overlap possibility (the biggest axis), if it's bigger than sig_max, then it's
        !impossible to have any overlap.
        if (rij .lt. sig_max) then
                call sigma_hgo(e1, e2, sig_hgo, urij)
                if (rij .le. sig_hgo) then
                        overlap = .TRUE.
                        return 
                end if       
        end if
        
        !If no overlap detected, returns false to outer subroutines
        overlap = .FALSE.
end subroutine check_overlap
        
subroutine compute_potential(rij, vij)
!This subroutine simply computes the potention on a pair contribution between 2 molecules based on the SQUARE WELL potential
!the depth width is defined by the parameter SIGMA SPHERE and the range (lambda)
        use global_variables
        implicit none
        real*8    :: rij                !Distance between 2 molecules
        real*8    :: vij                !Potential to each pair of molecules
              
        vij = 0.d0
         
        if (rij .lt. lambda*sigma) then
                vij = -1.d0
        else 
                vij = 0.d0
        end if

end subroutine compute_potential

subroutine rotation_quat(seed, max_ang, qold, qnew)
!This subroutine generates a random ROTATION QUATERNION that, multiplied by the current quaternion,
!generates the new quaterion (quat_new).

!This rotation quaternion is generated by a random axis (subroutine random_vector)
        use global_variables
        implicit none
        integer*8           :: seed         !seed for random number generator
        real*8              :: max_ang      !maximum angle displacement
        real*8              :: random_n     !random number
        real*8              :: axis(3)      !axis generated on random_vector subroutine (axis for rotation move)
        real*8              :: ang          !angle of rotation trial move
	    real*8              :: quat_rot(4)  !rotation quaternion generated
	    real*8              :: qold(4)      !old quaternion for molecule
        real*8,intent(out)  :: qnew(4)      !new quaternion for molecule
        
        !Initialization 
        quat_rot(:) = 0.d0    
        qnew(:)     = 0.d0 

        !Random number generator       
        call ranf(seed, random_n)
        random_n  = 2.d0*random_n - 1.d0    !random number in range [-1,1]
        ang       = max_ang*random_n        !random angle in (radian)
        
        !Create random axis -> ver rand_vec (joyce) -> unitario
        call random_vector(seed, axis)
       
        !Rotation quaternion generated
        quat_rot(1) = cos(ang*0.5)
        quat_rot(2) = sin(ang*0.5)*axis(1)
        quat_rot(3) = sin(ang*0.5)*axis(2)
        quat_rot(4) = sin(ang*0.5)*axis(3)
       
        !Multiply old and rotation quaternion, generating a new one
        qnew(1) = quat_rot(1)*qold(1) - quat_rot(2)*qold(2) - quat_rot(3)*qold(3) - quat_rot(4)*qold(4)       
        qnew(2) = quat_rot(2)*qold(1) + quat_rot(1)*qold(2) - quat_rot(4)*qold(3) + quat_rot(3)*qold(4)     
        qnew(3) = quat_rot(3)*qold(1) + quat_rot(4)*qold(2) + quat_rot(1)*qold(3) - quat_rot(2)*qold(4)     
        qnew(4) = quat_rot(4)*qold(1) - quat_rot(3)*qold(2) + quat_rot(2)*qold(3) + quat_rot(1)*qold(4)
        !q_new**2 = 1 --> Check 21-07
end subroutine rotation_quat

subroutine random_vector(seed, v)
!Generates a vector on the surface of a unit sphere.
!Allen and Tildesley 2th edition, page 514 (routine:Marsaglia 1972).
        implicit none
    	Real*8,    Intent(out) :: v(3) 
    	Real*8                 :: n1,n2,nsq 
        Integer*8, Intent(in)  :: seed  
        
        nsq = 2.d0
 
        Do while (nsq .gt. 1.d0)
                Call ranf(seed,n1) 
                Call ranf(seed,n2) 
                n1      = 2d0*n1 - 1.d0
                n2      = 2d0*n2 - 1.d0
                nsq     = n1*n1 + n2*n2
        End Do 
        v(1)    = 2.d0*n1*sqrt(1.d0-nsq)
        v(2)    = 2.d0*n2*sqrt(1.d0-nsq)
        v(3)    = 1.d0-2.d0*nsq

end subroutine random_vector

subroutine rotation(e_old, e_new, q)
!This subroutine will generate the rotation matrix, ALREADY TRANSPOSED, based on the new generated quaternion and the fixed
!orientation (e_old). It results on the new trial orientation to the particle (e_new)
        use global_variables
        implicit none
        real*8,dimension(3,3)  :: rotMT         !rotation matrix TRASNPOSED
        real*8,intent(in)      :: e_old(3)      !old orientation
        real*8,intent(in)      :: q(4)          !Quaternion
	    real*8,intent(out)     :: e_new(3)      !resulting new orientantion
               
        !initialization 
        e_new(:) = 0.d0
        
        !Rotation matrix transposed --> REFERENCE: Frenkel 2002 pg. 49 or Allen 2ed pg.110
        rotMT(1,1)  = q(1)*q(1) + q(2)*q(2) - q(3)*q(3) - q(4)*q(4)
        rotMT(1,2)  = 2.d0*(q(2)*q(3) - q(1)*q(4))
        rotMT(1,3)  = 2.d0*(q(2)*q(4) + q(1)*q(3))
        rotMT(2,1)  = 2.d0*(q(2)*q(3) + q(1)*q(4))
        rotMT(2,2)  = q(1)*q(1) - q(2)*q(2) + q(3)*q(3) - q(4)*q(4)
        rotMT(2,3)  = 2.d0*(q(3)*q(4) - q(1)*q(2))
        rotMT(3,1)  = 2.d0*(q(2)*q(4) - q(1)*q(3))
        rotMT(3,2)  = 2.d0*(q(3)*q(4) + q(1)*q(2))
        rotMT(3,3)  = q(1)*q(1) - q(2)*q(2) - q(3)*q(3) + q(4)*q(4)


        !e_new = AT * e -> AT = Rotation Matrix, e_old = fixed orientation 
        e_new(1) = rotMT(1,1)*e_old(1) + rotMT(1,2)*e_old(2) + rotMT(1,3)*e_old(3)
        e_new(2) = rotMT(2,1)*e_old(1) + rotMT(2,2)*e_old(2) + rotMT(2,3)*e_old(3)   
        e_new(3) = rotMT(3,1)*e_old(1) + rotMT(3,2)*e_old(2) + rotMT(3,3)*e_old(3)  
        !e_new**2 = 1 --> Check 21-07
end subroutine rotation

subroutine sigma_hgo(e1, e2, sig_hgo, r_versor)
!This subroutine simply computes the contact distance between 2 molecules (sigma HGO) based on their orientation and the unit vector
!between them, all inputs from check_overlap. Then, it returns the distance sig_hgo.
        use global_variables
        implicit none
        integer*8           :: i,j                                                               !index for particles
	    real*8              :: rxi, ryi, rzi, rxj, ryj, rzj, rxij, ryij, rzij                    !Positions for partciles
	    real*8              :: re1, re2, e1e2                                                    !Scalar product between orientation and distance
	    real*8              :: e1(3), e2(3)                                                      !Orientation vector for particles
	    real*8              :: r_versor(3)                                                       !Unit vector between 2 molecules
	    real*8, intent(out) :: sig_hgo                                                           !HGO contact distance --> output

        !Scalar products for each particle orientation 
        re1   = e1(1)*r_versor(1) + e1(2)*r_versor(2) + e1(3)*r_versor(3)
        re2   = e2(1)*r_versor(1) + e2(2)*r_versor(2) + e2(3)*r_versor(3)
        e1e2  = e1(1)*e2(1) + e1(2)*e2(2) + e1(3)*e2(3) 

        !Sigma_HGO calculation
        sig_hgo = (1.d0 - 0.5d0*chi*((re1 + re2)**2.d0/(1.d0+chi*e1e2) + (re1 - re2)**2.d0&
        /(1.d0-chi*e1e2)))**(-0.5d0)  
        !sig_hgo < sig_max  --> check 21-07 
end subroutine sigma_hgo

subroutine new_volume(seed, rnew, boxlnew, vnew, max_v)
!This subroutine generates the new trial volume move, it also escalate all positions on the new box volume
!the result of this subroutine will serve as input for new_potential below.
        use global_variables
        implicit none
        integer*8           :: i                                            !index for particles
        integer*8           :: seed                                         !seed for random number generator
	    real*8, intent(out) :: boxlnew(3), rnew(3, n_particles), vnew       !new box length, new positions and new potential 
	    real*8              :: random_n                                     !random number 
	    real*8              :: lnvnew, boxr                                 !log of new trial volume, box scale
	    real*8              :: max_v                                        !max volume displacement
        
        !Initialization of variables 
        lnvnew    = 0.d0
        vnew      = 0.d0
        boxr      = 0.d0
        rnew(:,:) = 0.d0
        
        !generate random number to volume trial move
        call ranf(seed, random_n)
        
        !Then, a new log volume is sampled --> REFERENCE: Allen & Tildesley 2ed pg. 167
        lnvnew = dlog(volume) + (2.d0*random_n - 1.d0)*max_v
        vnew   = dexp(lnvnew)

        !Based on the new volume, you can rescale box length based on a ratio of new and old volumes
        boxr   = (vnew/volume)**(1.d0/3.d0)
        boxlnew(:) = boxr * box_length(:) 
       
        !And also for all particles positions, all of them shall be rescaled based on the new volume
        do i = 1, n_particles
                rnew(1, i) = boxr * rx(i)
                rnew(2, i) = boxr * ry(i)
                rnew(3, i) = boxr * rz(i)
        end do
end subroutine new_volume

subroutine new_potential(rnew, boxlnew, overlap, new_pot) 
!This subsroutine recalculates the system total energy after having it's volume shifted, it has the
!same purpose of compute_total_energy used above, but with new box_lengths 
        use global_variables
        implicit none
        integer*8           :: i, j                                                     !Particle indexes and counters
	    real*8              :: rnew(3, n_particles)                                     !Array that holds the new scaled positions
	    real*8              :: rxi, ryi, rzi, rxj, ryj, rzj, rxij, ryij, rzij           !particle positions on each axis (x, y, z)
	    real*8              :: urij(3)                                                  !unit vector between 2 molecules
	    real*8              :: rij, rijsq                                               !module of the distance, module **2
	    real*8              :: boxlnew(3)                                               !box length after volume shift
	    real*8              :: e1(3), e2(3)                                             !orientations for each molecule
	    real*8, intent(out) :: new_pot                                                  !New total energy of the system
	    real*8              :: potij                                                    !potential for each 2 atoms being compared
        logical             :: overlap                                                  !flag for overlap detection of molecules
        
        !Initialization of pair potential and total system potential
        potij   = 0.d0
        new_pot = 0.d0
        
        !Loop for every particle, excluding the last which will be already counted
        do i = 1, n_particles-1
                e1  = e(i,:)
                rxi = rnew(1,i)
                ryi = rnew(2,i)
                rzi = rnew(3,i)
                do j = i+1, n_particles
                        e2      = e(j,:)
                        rxj     = rnew(1,j)
                        ryj     = rnew(2,j) 
                        rzj     = rnew(3,j)
                        
                        !Distance between two atoms 
                        rxij    = rxi - rxj
                        ryij    = ryi - ryj
                        rzij    = rzi - rzj

                        !Minimum Image Convention with new box lengths
                        rxij  = rxij - boxlnew(1)*anint(rxij/boxlnew(1))
                        ryij  = ryij - boxlnew(2)*anint(ryij/boxlnew(2))
                        rzij  = rzij - boxlnew(3)*anint(rzij/boxlnew(3))
                        
                        !Module of the distance between 2 molecules
                        rijsq   = rxij*rxij + ryij*ryij + rzij*rzij
                        rij     = dsqrt(rijsq)
                        
                        !Unit vector 
                        urij(1) = rxij/rij
                        urij(2) = ryij/rij
                        urij(3) = rzij/rij
                        
                        !Check if there is any overlap in new configuration
                        call check_overlap(e1, e2, rij, overlap, urij) 
                        
                        !If an overlap is detect, exit immediately to main 
                        if (overlap) then
                                return
                        !if no overlap detected, compute the potential 
                        else 
                                call compute_potential(rij, potij)
                                new_pot = new_pot + potij
                        end if
                end do
        end do
end subroutine new_potential

subroutine calculate_properties (h_dummy, h_dummy2, volume_dummy, volume_dummy2, vh_dummy)
!This subroutine is the final part for calculating fluid properties. Here you can find Cp, alfa, beta, and u_jt calculations.
!Main references: Allen & Tildesley 2nd ed, 2017, p.68 and 
!Aimoli, Maggin and Abreu, Force field comparison and thermodynamic property calculation of supercritical CO2 and CH4 using
!molecular dynamics simulations, 2014

        use global_variables
        implicit none
        real*8  :: h_dummy, h_dummy2, volume_dummy, volume_dummy2, vh_dummy
        real*8  :: b, c, e2, f    !constants for Passut and Danner correlation
        real*8  :: cp_id          !ideal part
        real*8  :: cp_r           !residual part 
        real*8  :: conversion     !conversion from BTU/lb*R -> J/kg K
        real*8  :: h, h2          !enthalpy, enthalpy**2 with units
        real*8  :: v, v2          !volume, volume**2 with units
        real*8  :: vh             !volume*enthalpy with units
        real*8  :: temp_unit      !temperature in K
        real*8  :: temp_unit_r    !temperature in R (used for Cp ideal contribution)
        real*8  :: cv, cp         !specific heat
        real*8  :: alfa, kt       !thermal expansion and isothermal compressibility 
        real*8  :: rho            !final density with units
        real*8  :: u_jt, cs       !Joule-Thomsom Coef and sound speed
        real*8  :: mmolar         !molar mass of co2
        logical :: file_exists    !logical for file opening
        character:: filename*100  !file name 
!-----------------------------------------------------------------------------------------
!                                    Returning units
!-----------------------------------------------------------------------------------------
        mtoang    = 10.d0**(-10.d0)
        mmolar    = 44.0095d0                                                                   !g/mol  
        temp_unit = temp*epslon                                                                 !K
        v         = volume_dummy*(sigs*mtoang)**3.d0/n_prod                                     !m^3
        v2        = volume_dummy2*((sigs*mtoang)**3.d0)**2.d0/n_prod                            !m^6                        
        h         = h_dummy*epslon*boltzmann*avogadro/n_particles/n_prod                        !J/mol
        h2        = h_dummy*(epslon*boltzmann*avogadro/n_particles)**2.d0/n_prod                !(J/mol)**2.d0
        vh        = vh_dummy*epslon*boltzmann*avogadro*((sigs*mtoang)**3.d0)/n_prod             !J/(mol*m^3)
!-----------------------------------------------------------------------------------------
!                                   Density (kg/m3)
!----------------------------------------------------------------------------------------- 
        !rho = n_particles/(volume_dummy/n_prod) 
        rho = dble(n_particles)/42168.4d0
        rho = rho*(mmolar/1000.d0)/(avogadro*(sigs*mtoang)**3.d0)
        write(*,*) rho
        stop
!-----------------------------------------------------------------------------------------
!                                   Specific Heat (Cp)
!-----------------------------------------------------------------------------------------
        
        !Passut and Danner coeficients for ideal gas cp calculation
        b = 0.114433d0
        c = 0.101132*10**(-3.d0)
        e2 = 0.034706*10**(-10.d0)
        f = -0.013140*10**(-14.d0)
        conversion = 4186.8d0 

        !For Cp calculation, first it's necessary to compute the ideal part
        !Main Reference: Passut and Danner 1972
        !Cp_id = B + 2CT + 3CT^2 + 4ET^3 + 5FT^4

        temp_unit_r = (temp*epslon)*1.8d0 !Conversion to Rankine 
        cp_id       = b + 2.d0*c*temp_unit_r + 3.d0*c*temp_unit_r**2.d0 + 4.d0*e2*temp_unit_r**3.d0 + 5.d0*f*temp_unit_r**4.d0
        cp_id       = cp_id*conversion    !converted to J/kg K

        !Cp_residual calculatiion             
        cp_r     = (h2 - h**2.d0)/(boltzmann*temp_unit**2.d0)
        cp       = cp_r + cp_id

!---------------------------------------------------------------------------------------
!                             Isothermal compressibility (Kt) 
!---------------------------------------------------------------------------------------
        kt = (v2 - v**2.d0)/(v*boltzmann*temp_unit)
       
!-------------------------------------------------------------------------------------
!                                      Alfa 
!-------------------------------------------------------------------------------------
        alfa = (vh - v*h)/(v*boltzmann*temp_unit**2.d0)

!-------------------------------------------------------------------------------------
!                               Specficic Heat (Cv)
!-------------------------------------------------------------------------------------
        cv = cp - temp_unit*v*alfa**2.d0/kt      

!-------------------------------------------------------------------------------------   
!                               Joule-Thomsom Coef (u_jt)
!-------------------------------------------------------------------------------------
        u_jt = v*(temp_unit*alfa-1)/cp

!------------------------------------------------------------------------------------
!                                 Sound Speed (cs)
!------------------------------------------------------------------------------------
        cs = (cp/(cv*kt*rho))**0.5d0

!------------------------------------------------------------------------------------
!                               Prop. writing file
!------------------------------------------------------------------------------------
        
!File opening --> if it exists, simply append the new calculations on the end of the file
!if it doesn't --> create the isothermal file
        filename = temp_c // "_CO2.out"
        inquire(FILE=filename, exist=file_exists) 
        if (file_exists) then
            open(101, FILE = filename, status="old", position="append", action="write")
                write(101,'(2a4, 7E15.7)')  temp_c, press_c, rho, cv, cp, kt, alfa, u_jt, cs
            close(101)
        else 
            open(101, FILE=filename)
                write(101, '(9a17)') "T (K)", "P (MPa)", "Density(kg/m^3)", "Cv (J/molK)", "Cp (J/molK)", "K_T (1/Pa)",&
                & "Alfa (1/K) ", "u_jt (K/MPa) ", "Sound_Speed (m/s) "
                write(101, '(2a4, 7E15.7)') temp_c, press_c, rho, cv, cp, kt, alfa, u_jt, cs
            close(101)
        end if  

!This file holds the reduces properties for Volume, Enthalpy and Volume*Enthalpy for posterior computation         
        filename = "Dummy_CO2.out"
        inquire(FILE=filename, exist=file_exists)
        if (file_exists) then
            open(101, FILE=filename, status="old", position="append",action="write")
                write(101, '(2a4, 5E15.7)') temp_c, press_c, volume_dummy, volume_dummy2, h_dummy, h_dummy2, vh_dummy
            close(101)
        else
            open(101, FILE=filename)
                write(101, '(7a17)') "T (K)", "P (MPa)", "V_reduced", "V^2_reduced", "H_reduced", "H^2_reduced", "V*H_reduced"
                write(101, '(2a4, 5E15.7)') temp_c, press_c, volume_dummy/n_prod, volume_dummy2/n_prod, h_dummy/n_prod,&
                &h_dummy2/n_prod, vh_dummy/n_prod
            close(101)
        end if 

end subroutine calculate_properties
