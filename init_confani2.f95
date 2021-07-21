 module init_confani2
        use global_variables

        contains
              subroutine fcc()
                implicit none 
                integer*8  :: i, j, k   !contadores
                integer*8  :: cont      !contadores
                real*8     :: lx,ly,lz  !Parametro da aresta
		real*8     :: lbox, ucl
		real*8     :: l, d
		real*8     :: a
                integer*8  :: uc 
                l           = sige
                d           = sigs  

                volume     = n_particles*v_particle/eta0/(sigs**3.d0)
                lbox       = volume**(1.d0/3.d0)

                !Numero de celulas por eixo para CFC
                uc      = nint((dble(n_particles)*0.25d0)**(1.d0/3.d0)) 
                ucl     = lbox/uc

                !tamanho da caixa em x -> rho = n_particles/volume -> volume = lx*ly*lz -> lx*lx*sige*lx   
                lx = (ucl**3.d0*d/l)**(1.d0/3.d0)           
                lz = lx*(sige/sigs)/(sigs/sigs)
                ly = lx

                !ajustando tamanho da caixa
                box_length(1) = lx*uc
                box_length(2) = ly*uc
                box_length(3) = lz*uc
               
                !Aqui está ok !				
 
                !Start the count
                cont = 1
                do i=1, uc
                        do j=1, uc
                                do k=1, uc
                                        !Vértices
                                        rx(cont) = (dble(i)-1.d0)*lx - 0.5d0*box_length(1)
                                        ry(cont) = (dble(j)-1.d0)*ly - 0.5d0*box_length(2)
                                        rz(cont) = (dble(k)-1.d0)*lz - 0.5d0*box_length(3)
                                        !Face frontal
                                        rx(cont+1) = (dble(i)-1.d0)*lx + 0.5d0*lx - 0.5d0*box_length(1)
                                        ry(cont+1) = (dble(j)-1.d0)*ly - 0.5d0*box_length(2)
                                        rz(cont+1) = (dble(k)-1.d0)*lz + 0.5d0*lz - 0.5d0*box_length(3)
                                        !Face lateral esquerda
                                        rx(cont+2) = (dble(i)-1.d0)*lx - 0.5d0*box_length(1)
                                        ry(cont+2) = (dble(j)-1.d0)*ly + 0.5d0*ly - 0.5d0*box_length(2)
                                        rz(cont+2) = (dble(k)-1.d0)*lz + 0.5d0*lz - 0.5d0*box_length(3)
                                        !Face lateral direita
                                        rx(cont+3) = (dble(i)-1.d0)*lx + 0.5d0*lx - 0.5d0*box_length(1)
                                        ry(cont+3) = (dble(j)-1.d0)*ly + 0.5d0*ly - 0.5d0*box_length(2)
                                        rz(cont+3) = (dble(k)-1.d0)*lz - 0.5d0*box_length(3)
                                        !passa para  as proximas particulas 
                                        cont = cont + 4
                                end do
                        end do
                 end do

    end subroutine fcc
end module init_confani2
