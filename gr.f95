module gr
        use global_variables
        contains
subroutine rdf ()
                implicit none
                integer*8           :: i, j      !contadores
                integer*8           :: bin !numero de bins total, bin atuak
                real*8              :: ru, rl  !variação em r, raio inferior, raio superior
                real*8              :: rij, rijsq, rxi, ryi, rzi, rxj, ryj, rzj, rxij, ryij, rzij !coordenads de posição, distancia entre
        !particulas
                real*8              :: ngi, cte !n particulas pro gas ideal, constante para g(r)
                

                nbin    = 400 
                deltar  = 0.5d0 * min(box_length(1), box_length(2), box_length(3))/dble(nbin)
                cte     = 4d0 * pi * rho_red/3d0
                hist(:) = 0d0
        
                do i = 1, n_particles - 1
                        do j = i+1, n_particles
                                rxi  = rx(i)
                                ryi  = ry(i)
                                rzi  = rz(i)
                                rxj  = rx(j)
                                ryj  = ry(j)
                                rzj  = rz(j)
                                rxij = rxi - rxj
                                ryij = ryi - ryj
                                rzij = rzi - rzj
                                !Minimum image convention
                                rxij = rxij - box_length(1)*anint(rxij/box_length(1))
                                ryij = ryij - box_length(2)*anint(ryij/box_length(2))
                                rzij = rzij - box_length(3)*anint(rzij/box_length(3))
                                !distance between moleules 
                                rijsq     = rxij*rxij + ryij*ryij + rzij*rzij        
                                rij       = dble(sqrt(rijsq))
                                bin       = floor(rij/deltar) + 1
                                if (bin <= nbin) then
                                        hist(bin) = hist(bin) + 2.d0
                                end if  
                        end do
                end do

                do i = 1, nbin
                        rl      = deltar*dble(i-1)
                        ru      = rl + deltar
                        ngi     = cte*(ru**3.d0-rl**3.d0)
                        hist(i) = hist(i)/ngi/dble(n_particles)
                end do
        
    
end subroutine rdf
end module gr
