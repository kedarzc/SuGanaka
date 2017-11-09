program CST_main

    use MainModule

    implicit none

    ! PART I
    ! DECLARE THE VARIABLES
    real, allocatable :: geom(:,:)          ! Geom Stores the nodal co-ordinates
    real, allocatable :: node_disp(:,:)     ! Stores the nodal displacements
    real, allocatable :: dee(:,:)
    real, allocatable :: nf(:,:)            ! Stores the nodal degrees of freedom
    real, allocatable :: nodal_loads(:,:)   ! Stores the loads on the nodes
    real, allocatable :: fg(:,:)
    real, allocatable :: bee(:,:)   
    real, allocatable :: ke(:,:)            ! Element stiffness matrix
    real, allocatable :: kk(:,:)            ! Global stiffness matrix
    real, allocatable :: kk_inv(:,:)        ! Inverse of the global stiffness matrix
    real, allocatable :: delta(:,:)         ! Stores calculated displacements
    real, allocatable :: ndnums(:)
    real, allocatable :: eld(:,:)           ! Stores the element displacement vector 

    real, allocatable :: EPS(:,:)           ! Store the elemental strains
    real, allocatable :: SIGMA(:,:)         ! Store the elemental stresses
    
    integer, allocatable :: g(:,:)
    integer, allocatable :: connec(:,:)        ! Stores element Connectivity

    integer :: nnd              ! Total number of nodes in the model
    integer :: nodof            ! Number of degrees of freedom per element
    integer :: nne              ! Number of nodes per element
    integer :: eldof            ! Element degree of freedom
    integer :: nel              ! Number of elements
    integer :: active_nf        ! Active degrees of freedom
    integer :: elem_num         ! Stores the current element number
    integer :: i, j, k, l       ! Counters 

    real :: E       ! Elastic modulus
    real :: nu      ! This is the poissons ratio
    real :: thick   ! Thickness of the element
    real :: c
    real :: Area    ! Area of the element
    
    character(len=50)    :: fileName_Nodes   ! File with nodal info
     character(len=50)   :: fileName_Connec  ! Nodal Connectivity

    ! Variables for reading the input file
    integer :: status

    ! Declaring variables to be user by the visualizer
    integer :: scaling_factor = 100


    ! INPUTS
    ! =====================================================

    ! PART II
    ! INITIALIZE THE VARIABLES

    ! Please enter the inputs
    ! ----------------------

    fileName_Nodes  = 'MeshData_ABQ.csv'
    fileName_Connec = 'Connectivity_ABQ.csv'

    ! Count the number of nodes and elements
    call countLines(fileName_Nodes,nnd)
    call countLines(fileName_Connec,nel)

    nodof = 2   ! Number of degrees of freedom per element
    nne   = 3   ! Number of nodes per Element
    eldof = nne*nodof

    ! Materials
    ! ---------
    E  = 200000 ! MPa
    nu = 0.2    ! Poisson's ratio 

    ! Section thickness
    thick = 5! Beam thickness in mm

    ! Initialize the nodal freedom matrix to 1
    allocate(nf(nnd,nodof))
    call ones_matrix(nnd,nodof,nf)

    ! Boundary Conditions
    ! -------------------
    nf(25,1) = 0.0   ; nf(25,2) = 0.0;
    nf(50,1) = 0.0   ; nf(50,2) = 0.0;
    nf(75,1) = 0.0   ; nf(75,2) = 0.0;
    nf(100,1) = 0.0   ; nf(100,2) = 0.0;
    nf(125,1) = 0.0   ; nf(125,2) = 0.0;   
    nf(150,1) = 0.0   ; nf(150,2) = 0.0;
    nf(175,1) = 0.0   ; nf(175,2) = 0.0;
    nf(200,1) = 0.0   ; nf(200,2) = 0.0;
    nf(225,1) = 0.0   ; nf(225,2) = 0.0;
    nf(250,1) = 0.0   ; nf(250,2) = 0.0;
    nf(275,1) = 0.0   ; nf(275,2) = 0.0;

    ! Loading
    ! -------
    allocate(nodal_loads(nnd,nodof))
    call zeros_matrix(nnd,nodof,nodal_loads)

    !  Node on load 2
    nodal_loads(1,1) = 0.0
    nodal_loads(1,2) = -1000.0

    ! END OF INPUTS
    ! =====================================================



    ! PART III
    ! ALLOCATE SPACE TO THE ARRAYS

    allocate(geom(nnd,nodof))
    allocate(node_disp(nnd,nodof))
    allocate(connec(nel,nne))

    ! Read the Mesh Data
    open(unit=15, file=fileName_Nodes, status='old', action='read', iostat=status)
    read(15,*) ((geom(i,j),j=1,nodof),i=1,nnd)

    ! Read the Element Connectivity
    open(unit=16, file=fileName_Connec, status='old', action='read', iostat=status)
    read(16,*) ((connec(i,j),j=1,nne),i=1,nel)

    allocate(dee(nodof+1,nodof+1))
    allocate(bee(nne,eldof))
    allocate(g(eldof,1))
    allocate(ke(eldof,eldof))
    allocate(delta(active_nf,1))
    allocate(eld(eldof,1))

    allocate(EPS(nel,nne))
    allocate(SIGMA(nel,nne))

    call formdsig(E,nu,dee)

    ! Find the active degrees of freedom and update the 
    ! nf matrix
    call find_active_dof(nnd,nodof,active_nf,nf)

    ! Global force vector
    ! ------------------
    ! This force vector will have rows = active_nf 
    !                      and columns = 1

    ! Allocate memory to the global force vector
    allocate(fg(active_nf,1))
    allocate(kk(active_nf,active_nf))
    allocate(kk_inv(active_nf,active_nf))

    
    ! Form the gobal force vector
    call form_global_force_vector(nnd,nodof,nf,fg, active_nf, nodal_loads)

    do elem_num=1,nel

        call elem_CST(geom,connec,elem_num,Area,bee,nnd,nodof,nel,nne,eldof,nf,g)

        ! Form the element stiffness matrix
        ke = thick*Area*matmul(matmul(transpose(bee),dee),bee)

        ! Form the global stiffness matrix
        call form_kk(kk,ke,g,active_nf,eldof)

    enddo

    ! Invert the stiffness matrix
    call inverse(kk,kk_inv,active_nf)

    ! Compute the required dedlections
    delta = matmul(kk_inv,fg)

    ! Compute the stresses and strains in the elements
    do elem_num=1,nel

        call elem_CST(geom,connec,elem_num,Area,bee,nnd,nodof,nel,nne,eldof,nf,g)

        ! Retrieve the displacements for the element at each node

        do j=1,eldof
            if (g(j,1)==0) then
                eld(j,1) = 0
            else
                eld(j,1) = delta(g(j,1),1)
            endif
        end do

        ! Compute the strains for the element
        EPS(elem_num:elem_num,1:3)   = transpose(matmul(bee,eld))
        SIGMA(elem_num:elem_num,1:3) = transpose(matmul(dee,matmul(bee,eld)))


    end do 

    ! Compute the final co-ordinates(co-ordinate + delta)
    do i =1,nnd

        if (nf(i,1) /= 0.0) then
            node_disp(i,1) = delta(int(nf(i,1)),1)
        else
            node_disp(i,1) = 0.0
        endif

        if (nf(i,2) /= 0.0) then
            node_disp(i,2) = delta(int(nf(i,2)),1)
        else
            node_disp(i,2) = 0.0
        endif

    end do

    open(unit = 2789, file = "Results.sdb")

    write(2789,*) 'Nodal displacements'
    write(2789,*) ''

    do i =1,nnd
        write(2789,*) node_disp(i,1), node_disp(i,2)
    end do        

    write(2789,*) '' 
    write(2789,*) 'Strains'                  

    do i =1,nel
        write(2789,*) EPS(i,1),EPS(i,2),EPS(i,3)
    end do 
    
    write(2789,*) '' 
    write(2789,*) 'Streses'                  

    do i =1,nel
        write(2789,*) SIGMA(i,1), SIGMA(i,2), SIGMA(i,3)
    end do 


    close (2789)        





























    ! This is a major addition. The next part of code prints out the 
    ! VTK file for visualization in Paraview


    ! List of variables needed
    ! 1. Name of results file
    ! 2. 










    open(unit = 2788, file = "Results.vtk")
    write(2788,"(A)") '# vtk DataFile Version 1.0'
    write(2788,"(A)") '2D Unstructured Grid of Linear Triangles'
    write(2788,"(A)") 'ASCII'

    ! Nodes
    write(2788,"(A)") ''
    write(2788,"(A)") 'DATASET UNSTRUCTURED_GRID'
    write(2788,"(A,I10,A)") 'POINTS ', nnd, ' float'

    do i=1,nnd
        write(2788,"(F15.4,F15.4,F15.4)") geom(i,1) + scaling_factor*node_disp(i,1), &
                                          geom(i,2) + scaling_factor*node_disp(i,2), &
                                          0.0000    + scaling_factor*0.0000
    end do

    ! Elements
    write(2788,"(A)") ''
    write(2788,"(A,I10,I10)") 'CELLS ', nel, nel*(1+nne)

    ! Subtract one from the Connectivity since the index starts from 0 
    ! VERY IMPORTANT
    do i=1,nel
        write(2788,*) nne, int(connec(i,1))-1, int(connec(i,2))-1, int(connec(i,3))-1
    end do

    ! Types of cells in vtk
    write(2788,"(A)") ''
    write(2788,"(A,I10)") 'CELL_TYPES ', nel

    do i=1,nel
        write(2788,"(I3)") 5
    end do
    
    ! Write the nodal displacements
    write(2788,"(A)") ''
    write(2788,"(A,I10)") 'POINT_DATA ', nnd
    write(2788,"(A)") 'VECTORS displacement float'
    do i=1,nnd
        write(2788,"(E20.4,E20.4,E20.4)") node_disp(i,1), node_disp(i,2), 0.0000
    end do

    close (2788)


end program CST_main
