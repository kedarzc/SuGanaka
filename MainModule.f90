module MainModule

	contains 

    ! This subroutine forms the elasticity matrix for a plane stress problem
    ! ----------------------------------------------------------------------
    subroutine formdsig(E,nu,dee)

        implicit none

        real :: c
        real, intent (in) :: E,nu
        real, intent(out) :: dee(3,3)

        c = E/(1-nu**2)

        dee = c*reshape((/ 1.0, nu, 0.0, nu, 1.0, 0.0, 0.0, 0.0, (1.0-nu)/2.0 /), (/3,3/))

    end subroutine formdsig



    ! This subroutine sets all elenents of the matrix to 1
    ! ----------------------------------------------------
    subroutine ones_matrix(rows,cols,a)

        implicit none

        integer :: i,j, rows, cols
        real, dimension(rows,cols) :: a

        do i=1,rows
            do j=1,cols
                a(i,j) = 1.0
            enddo
        enddo

    end subroutine ones_matrix


    ! This subroutine sets all elenents of the matrix to 0
    ! ----------------------------------------------------
    subroutine zeros_matrix(rows,cols,a)

        implicit none

        integer :: i,j, rows, cols
        real, dimension(rows,cols) :: a

        ! Initialize the global stiffness matrix
        do i=1,rows
            do j=1,cols
                a(i,j) = 0.0
            enddo
        enddo

    end subroutine zeros_matrix




    
    ! Find the active degrees of freedom and update the nf matrix
    subroutine find_active_dof(nnd,nodof,active_nf,nf)

        implicit none
        
        integer :: active_nf, i, j
        integer, intent(in) :: nnd, nodof
        real, intent(out) :: nf(nnd,nodof)

        active_nf = 0
        do i=1,nnd
            do j=1,nodof
                if (nf(i,j) == 1) then
                    active_nf = active_nf + 1
                    nf(i,j) = active_nf
                end if
            enddo
        enddo

    end subroutine find_active_dof




    ! Form the global force vector (load vector)
    subroutine form_global_force_vector(nnd,nodof,nf,fg,active_nf,nodal_loads)

        integer, intent(in) :: nnd,nodof,active_nf
        real, intent (in)   :: nodal_loads(nnd,nodof), nf(nnd,nodof)
        real, intent(out)   :: fg(active_nf,1)

        do i=1,nnd
            do j=1,nodof
                if (nf(i,j) > 0) then
                    fg(int(nf(i,j)),1) = nodal_loads(i,j)
                end if
            enddo
        enddo 

    end subroutine form_global_force_vector


    ! A subroutine that calculates the determinant of a 3 by 3 matrix
    subroutine determinant_3by3(myMatrix,det)

        implicit none

        real :: det, myMatrix(3,3)

        det = myMatrix(1,1)*myMatrix(2,2)*myMatrix(3,3)  &
            - myMatrix(1,1)*myMatrix(2,3)*myMatrix(3,2)  &
            - myMatrix(1,2)*myMatrix(2,1)*myMatrix(3,3)  &
            + myMatrix(1,2)*myMatrix(2,3)*myMatrix(3,1)  &
            + myMatrix(1,3)*myMatrix(2,1)*myMatrix(3,2)  &
            - myMatrix(1,3)*myMatrix(2,2)*myMatrix(3,1)

    end subroutine determinant_3by3




    subroutine form_kk(kk,ke,g,active_nf,eldof)

        integer, intent(in) :: eldof, active_nf, g(eldof,1)
        real, intent(in) :: ke(eldof,eldof)
        real, intent(out) :: kk(active_nf,active_nf)

        ! update the global stiffness matrix
        do l=1,eldof
            if ((g(l,1) == 0) .eqv. .False.) then
                do m=1,eldof
                    if ((g(m,1) == 0) .eqv. .False.) then
                        kk(g(l,1), g(m,1)) = &
                            kk(g(l,1), g(m,1)) + ke(l,m)
                    end if
                enddo
            end if
        enddo


    end subroutine form_kk



subroutine inverse(KK,KK_inv,active_nf)

    !============================================================
    ! Inverse matrix
    ! Method: Based on Doolittle LU factorization for Ax=b
    ! Alex G. December 2009
    !-----------------------------------------------------------
    ! input ...
    ! KK(active_nf,active_nf) - array of coefficients for matrix A
    ! active_nf      - dimension
    ! output ...
    ! KK_inv(active_nf,active_nf) - inverse matrix of A
    ! comments ...
    ! the original matrix KK(active_nf,active_nf) will be destroyed
    ! during the calculation
    !===========================================================
    
    implicit none
    integer active_nf
    real KK(active_nf,active_nf), KK_inv(active_nf,active_nf)
    real L(active_nf,active_nf), U(active_nf,active_nf), b(active_nf),&
         d(active_nf), x(active_nf)
    real coeff
    integer i, j, k

    ! step 0: initialization for matrices L and U and b
    ! Fortran 90/95 aloows such operations on matrices
    L=0.0
    U=0.0
    b=0.0

    ! step 1: forward elimination
    do k=1, active_nf-1
       do i=k+1,active_nf
          coeff=KK(i,k)/KK(k,k)
          L(i,k) = coeff
          do j=k+1,active_nf
             KK(i,j) = KK(i,j)-coeff*KK(k,j)
          end do
       end do
    end do

    ! Step 2: prepare L and U matrices
    ! L matrix is KK matrix of the elimination coefficient
    ! + the diagonal elements are 1.0
    do i=1,active_nf
      L(i,i) = 1.0
    end do
    ! U matrix is the upper triangular part of A
    do j=1,active_nf
      do i=1,j
        U(i,j) = KK(i,j)
      end do
    end do

    ! Step 3: compute columns of the inverse matrix C
    do k=1,active_nf
      b(k)=1.0
      d(1) = b(1)
    ! Step 3a: Solve Ld=b using the forward substitution
      do i=2,active_nf
        d(i)=b(i)
        do j=1,i-1
          d(i) = d(i) - L(i,j)*d(j)
        end do
      end do
    ! Step 3b: Solve Ux=d using the back substitution
      x(active_nf)=d(active_nf)/U(active_nf,active_nf)
      do i = active_nf-1,1,-1
        x(i) = d(i)
        do j=active_nf,i+1,-1
          x(i)=x(i)-U(i,j)*x(j)
        end do
        x(i) = x(i)/u(i,i)
      end do
    ! Step 3c: fill the solutions x(active_nf) into column k of C
      do i=1,active_nf
        KK_inv(i,k) = x(i)
      end do
      b(k)=0.0
    end do
    
end subroutine inverse




    ! Calculate the Area of CST element
    ! ---------------------------------
    subroutine elem_CST(geom,connec,elem_num,Area,bee,nnd,nodof,nel,nne,eldof,nf,g)

        implicit none

        integer, intent(in) :: elem_num, connec(nel,nne),nnd,nodof,nel,nne,eldof
        real, intent(in) :: geom(nnd,nodof), nf(nnd,nodof)

        real, intent(out) :: Area, bee(nne,eldof)
        integer, intent(out) :: g(eldof,1)

        real :: det, x1,y1,x2,y2,x3,y3, elem_coords(3,3), m(3,3)
        integer :: j,k,l

        ! Extract the x,y-coordinates 
        ! --------------------------

        ! Node 1
        x1 = geom(connec(elem_num,1),1)
        y1 = geom(connec(elem_num,1),2)

        ! Node 2
        x2 = geom(connec(elem_num,2),1)
        y2 = geom(connec(elem_num,2),2)

        ! Node 3
        x3 = geom(connec(elem_num,3),1)
        y3 = geom(connec(elem_num,3),2)

        ! Find the area of the element
        elem_coords = reshape((/ 1.0, 1.0, 1.0, &
                                x1, x2, x3, &
                                y1, y2, y3 /), (/3,3/))
        
        call determinant_3by3(transpose(elem_coords),det)

        Area = 0.5*det

        ! Compute B matrix

        m(1,1) = (x2*y3 - x3*y2)/(2.0*Area)
        m(1,2) = (y2-y3)/(2.0*Area)
        m(1,3) = (x3-x2)/(2.0*Area)

        m(2,1) = (x3*y1-x1*y3)/(2.0*Area)
        m(2,2) = (y3-y1)/(2.0*Area)
        m(2,3) = (x1-x3)/(2.0*Area)

        m(3,1) = (x1*y2-x2*y1)/(2.0*Area)
        m(3,2) = (y1-y2)/(2.0*Area)
        m(3,3) = (x2-x1)/(2.0*Area)

        bee = reshape((/ m(1,2), 0.0, m(1,3), & 
                         0.0, m(1,3), m(1,2), &
                         m(2,2), 0.0, m(2,3), &
                         0.0, m(2,3), m(2,2), &
                         m(3,2), 0.0, m(3,3), &
                         0.0, m(3,3), m(3,2) /), (/ 3,6 /))

        ! Form the steering vector
        l = 0

        do k=1,nne
            do j=1,nodof
                l = l+1
                g(l,1) = nf(connec(elem_num,k),j)
            enddo
        enddo


    end subroutine elem_CST


    ! This subroutine counts the number of lines in a file
    subroutine countLines(fileName,nlines)

        implicit none

        character(len=*), intent(in) :: fileName
        integer, intent(out) :: nlines

        integer :: io, status

        nlines = 0

        open(10,file=fileName,status='old', action='read', iostat=status)

        do 
            read(10,*,iostat=io) 

            if (io/=0) then
                exit 
            else
                nlines = nlines + 1
            end if
        enddo


        close(10)

    end subroutine countLines


end module MainModule