program chgsum
    ! Read two CHGCAR file and compute a linear combination. Three arguments 
    ! are mandatory:
    !     * args 1: CHGCAR file 1
    !     * args 2: CHGCAR file 2
    !     * args 3: scale_factor
    !
    ! Linear combination:
    !    CHGCAR1 + scale_factor * CHGCAR2
    !
    ! Output is written into the file CHGCAR_sum
    !

    implicit none

    ! arguments:
    character(len=100)                            :: chgcar1, chgcar2, fact

    ! local:
    integer                                       :: i, j, k, nline
    integer                                       :: ngxf1, ngyf1, ngzf1
    integer                                       :: ngxf2, ngyf2, ngzf2
    character(len=100)                            :: line
    character(len=100), dimension(:), allocatable :: header
    real(kind=8)                                  :: scale_factor
    real(kind=8), dimension(:,:,:), allocatable   :: rho1, rho2

    ! get arguments
    call get_command_argument(1, chgcar1)
    call get_command_argument(2, chgcar2)
    call get_command_argument(3, fact)
    read(fact, *) scale_factor

    !write(*,*) trim(chgcar1)
    !write(*,*) trim(chgcar2)
    !write(*,*) scale_factor

    ! read file 1
    write(*, "('Reading file 1  : ', a)") trim(chgcar1)
    open(unit=10, file=chgcar1, action="read")
    nline = 1
    read(10, "(a)") line
    do while (.not. line == "")
        read(10, "(a)") line
        nline = nline + 1
    end do
    !write(*,*) nline
    
    read(10, *) ngxf1, ngyf1, ngzf1
    !write(*,"(3i6)") ngxf1, ngyf1, ngzf1
    allocate(rho1(ngxf1, ngyf1, ngzf1))
    allocate(header(nline))
    read(10, *) (((rho1(i, j, k), i = 1, ngxf1), j = 1, ngyf1), k = 1, ngzf1)
    close(10)

    ! read file 2
    write(*, "('Reading file 2  : ', a)") trim(chgcar2)
    open(unit=10, file=chgcar2, action="read")
    do i = 1, nline
        read(10, "(a)") header(i)
    end do
    read(10, *) ngxf2, ngyf2, ngzf2
    !write(*,"(3i6)") ngxf2, ngyf2, ngzf2
    if( ngxf1 /= ngxf2 .or. ngyf1 /= ngyf2 .or. ngzf1 /= ngzf2) then
        stop "unconsistant grid"
    end if
    allocate(rho2(ngxf2, ngyf2, ngzf2))
    read(10, *) (((rho2(i, j, k), i = 1, ngxf2), j = 1, ngyf2), k = 1, ngzf2)
    close(10)

    ! do the sum
    write(*, "('Compute ', a, ' + ', F5.2, a)") trim(chgcar1), scale_factor, &
        " " // trim(chgcar2)
    rho1(:,:,:) = rho1(:,:,:) + scale_factor * rho2(:,:,:)

    ! write sum file
    write(*, "('Writting file 2 : ', a)") "CHGCAR_sum"
    open(unit=10, file="CHGCAR_sum", action="write")
    do i = 1, nline
        write(10, "(a)") trim(header(i))
    end do
    write(10, "(3i6)") ngxf1, ngyf1, ngzf1
    write(10, "(5E19.11)") (((rho1(i, j, k), i = 1, ngxf1), j = 1, ngyf1), & 
        k = 1, ngzf1)
    close(10)

end program chgsum
