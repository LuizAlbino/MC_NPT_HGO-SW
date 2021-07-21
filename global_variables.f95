module global_variables
        implicit none
        
        integer*8         :: n_particles                        !numero de particulas
        real*8, parameter :: pi        = 3.14159265359          !numero pi
        real*8, parameter :: avogadro  = 6.0221409d23           !numero de avogadro
        real*8, parameter :: boltzmann = 1.38064852d-23         !cte de boltzmann
        real*8, parameter :: r         = avogadro*boltzmann     !Coeficiente do gás ideal
        real*8            :: rho, molar_mass                    !kg/m3, g/mol e Angstrons
        real*8            :: box_length(3)                      !tamanho das arestas da caixa
        real*8            :: eta0, rho_red                      !densidade numerica e reduzida  
        real*8            :: mtoang                             !conversão de metros para angstrons
        real*8            :: rx(100000)                         !Posição das particulas em X
        real*8            :: ry(100000)                         !Posição das particulas em Y
        real*8            :: rz(100000)                         !Posição das particulas em Z
        real*8            :: e(100000,3)                        !vetor de orientações X, Y, Z
        real*8            :: efixed(3)              
        real*8            :: volume, v_particle                 !Volume da caixa de simulação
        real*8            :: rcut                               !rcutoff
        real*8            :: rcutsq                             !r cut off ao quadrado
        real*8            :: temp                               !temeperatura
        real*8            :: sigma,  sig_max                    !sigma e sigma2
        real*8            :: lambda, epslon, lambdasq           !tamanho do poço, valor dod poço
        real*8            :: sigs, sige, ke, chi                !parametros anisotrópicos Joyce --> p.44
        real*8            :: hist(100000)                       !histograma pro calculo da g(r)
        real*8            :: deltar                             !variação de r pro calculo da g(r)
        integer*8         :: steps                              !numero de passos
        integer*8         :: nbin                               !numero de bins para g(r)
        integer*8         :: max_steps                          !max de passos
        integer*8         :: n_save                             !intervalo de salvar os valores
        integer*8         :: n_equil                            !passos para equilibração
        integer*8         :: n_prod                             !passos para produção
        integer*8         :: n_adjust                           !intervalo para ajustar drmax
        integer*8         :: nacc                               !numero de transições aceitas
        character         :: atom*1, dummy*1                    !letra identificadora dos atomos

end module global_variables
