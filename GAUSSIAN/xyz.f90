! ******************************************************************************
!
! programme : xyz
!
! description : genere un fichier au format xyz a partir d'une geometrie lu
!               sur un fichier .log ou .com de gaussian
!
! format xyz :
! ligne 1 = nbre d'atome
! ligne 2 = titre
! ligne 3 -> nbre d'atome = nomatome x y z
!
! syntaxe : xyz mon_calcul.log
! fichier de sortie mon_calcul.xyz
!
! dans un .log en entree
! Extrait la derniere geometrie Input orientation trouvee
!
! dans un .com contenant des coordonnées cartesienne
! Extrait la geometrie et l'ecrit au format xyz
!
! ******************************************************************************

program xyz

  implicit none

  integer, parameter                 :: ndim = 1000
  integer                            :: nat, iat, k, io
  integer                            :: klog, kcom
  character(len=100)                 :: fichier, ligne
  character(len=4), dimension(ndim)  :: nomcom
  character(len=4)                   :: nom
  integer, dimension(ndim)           :: mass
  double precision,dimension(ndim,3) :: x
  logical                            :: existe, fin, trouve

  write(*,*)

  call get_command_argument(1,fichier)
  write(*,*)"lecture du fichier "//trim(fichier)

! test d'existence du fichier et reconnaissance
  inquire( file=trim(fichier) , exist=existe)
  if( .not. existe ) then
     write(*,*)" fichier non trouve: "//trim(fichier)
     stop
  end if

  ! teste .com ou .log
  klog = 0
  kcom = 0
  if( index(fichier,".log") /= 0 .or. index(fichier,".out") /= 0 ) klog = 1
  if( index(fichier,".com") /= 0 ) kcom = 1
  if( klog == 0 .and. kcom == 0) then
    write(*,*)" type .log ou .com non identifie "//trim(fichier)
    write(*,*)"le nom du fichier doit finir par .log ou .com"
    stop
  end if
  if( klog == 1 .and. kcom == 1) then
    write(*,*)" type .log ou .com non identifie "//trim(fichier)
    write(*,*)" le nom du fichier contient .com et .log"
    stop
  end if

  nom(:) = "XX"

! LECTURE DU FICHIER

  if( klog == 1 ) then

! CAS DU FICHIER LOG

    write(*,*)"c'est un fichier .log"
    open(10,file=trim(fichier),action="read")
    fin = .false.
    do while ( .not. fin )
      trouve = .false.
      do while( .not. trouve )
        read( 10, "(a)", iostat=io)ligne
        if( io < 0 ) then
          fin = .true.
          trouve = .true.
        end if
        if( index(ligne,"Input orientation") /= 0 ) trouve = .true.
        if( index(ligne,"Z-Matrix orientation") /= 0 ) trouve = .true.
      end do

      if( .not. fin ) then
        read(10,*)
        read(10,*)
        read(10,*)
        read(10,*)
        iat = 0
        trouve = .false.
        do while ( .not. trouve )
          ligne = ""
          read(10,"(a)") ligne
          if(index(ligne,"-------") /= 0 ) then
            trouve = .true.
          else
            iat = iat + 1
            read(ligne,*) k, mass(iat), k, (x(iat,k), k=1,3)
          end if
        end do
        nat = iat
      end if

    end do
    close(10)

  else if( kcom == 1 ) then

! CAS DU FICHIER COM

    write(*,*)"c'est un fichier .com"
    open(10,file=trim(fichier),action="read")
    ! calcul
    fin = .false.
    do while ( .not. fin )
      ligne=""
      read(10,"(a)")ligne
      if( len_trim(ligne) == 0) fin = .true. ! ligne blanche
    end do
    ! titre
    fin = .false.
    do while ( .not. fin )
      ligne=""
      read(10,"(a)")ligne
      if( len_trim(ligne) == 0) fin = .true. ! ligne blanche
    end do
    ! ligne charge multiplicite
    read(10,"(a)")ligne
    ! atomes
    iat=0
    fin = .false.
    do while ( .not. fin )
      ligne=""
      read(10,"(a)")ligne
      if( ligne == "" ) then
        fin = .true.
      else
        iat = iat  + 1
        read(ligne,*) nomcom(iat), (x(iat,k), k=1,3)
      end if
    end do
    nat=iat
    close(10)
  else

    write(*,*)"type de fichier inconnu"
    stop

  end if

  write(*,*)"nombre d'atome lu ",nat

! ecriture du fichier xyz

  if( klog == 1 ) then
    if( index(fichier,".log") /= 0 ) k=index(fichier,".log")
    if( index(fichier,".out") /= 0 ) k=index(fichier,".out")

    open(11,file=trim(fichier(:k-1) )//".xyz" )
    write(*,*)"fichier xyz => "//trim(fichier(:k-1) )//".xyz"
  else if(kcom == 1 ) then
    k=index(fichier,".com")
    open(11,file=trim(fichier(:k-1) )//".xyz" )
    write(*,*)"fichier xyz => "//trim(fichier(:k-1) )//".xyz"
  else
    write(*,*)"type de fichier inconnu"
    stop
  end if

  write(ligne,*)nat

  write(11,"(a)")trim(adjustl(ligne))
  write(11,"(a)")"geometrie issue du fichier "//trim(fichier)

  do iat=1,nat

    if( klog == 1 ) then

      if( mass(iat) == 1 ) nom = "H"
      if( mass(iat) == 2 ) nom = "He"

      if( mass(iat) == 3 ) nom = "Li"
      if( mass(iat) == 4 ) nom = "Be"
      if( mass(iat) == 5 ) nom = "B"
      if( mass(iat) == 6 ) nom = "C"
      if( mass(iat) == 7 ) nom = "N"
      if( mass(iat) == 8 ) nom = "O"
      if( mass(iat) == 9 ) nom = "F"
      if( mass(iat) ==10 ) nom = "Ne"

      if( mass(iat) ==11 ) nom = "Na"
      if( mass(iat) ==12 ) nom = "Mg"
      if( mass(iat) ==13 ) nom = "Al"
      if( mass(iat) ==14 ) nom = "Si"
      if( mass(iat) ==15 ) nom = "P"
      if( mass(iat) ==16 ) nom = "S"
      if( mass(iat) ==17 ) nom = "Cl"
      if( mass(iat) ==18 ) nom = "Ar"

      if( mass(iat) ==19 ) nom = "K"
      if( mass(iat) ==20 ) nom = "Ca"
      if( mass(iat) ==21 ) nom = "Sc"
      if( mass(iat) ==22 ) nom = "Ti"
      if( mass(iat) ==23 ) nom = "V"
      if( mass(iat) ==24 ) nom = "Cr"
      if( mass(iat) ==25 ) nom = "Mn"
      if( mass(iat) ==26 ) nom = "Fe"
      if( mass(iat) ==27 ) nom = "Co"
      if( mass(iat) ==28 ) nom = "Ni"
      if( mass(iat) ==29 ) nom = "Cu"
      if( mass(iat) ==30 ) nom = "Zn"
      if( mass(iat) ==31 ) nom = "Ga"
      if( mass(iat) ==32 ) nom = "Ge"
      if( mass(iat) ==33 ) nom = "As"
      if( mass(iat) ==34 ) nom = "Se"
      if( mass(iat) ==35 ) nom = "Br"
      if( mass(iat) ==36 ) nom = "Kr"

      if( mass(iat) ==37 ) nom = "Rb"
      if( mass(iat) ==38 ) nom = "Sr"
      if( mass(iat) ==39 ) nom = "Y"
      if( mass(iat) ==40 ) nom = "Zr"
      if( mass(iat) ==41 ) nom = "Nb"
      if( mass(iat) ==42 ) nom = "Mo"
      if( mass(iat) ==43 ) nom = "Tc"
      if( mass(iat) ==44 ) nom = "RU"
      if( mass(iat) ==45 ) nom = "Rh"
      if( mass(iat) ==46 ) nom = "Pd"
      if( mass(iat) ==47 ) nom = "Ag"
      if( mass(iat) ==48 ) nom = "Cd"
      if( mass(iat) ==49 ) nom = "In"
      if( mass(iat) ==50 ) nom = "Sn"
      if( mass(iat) ==51 ) nom = "Sb"
      if( mass(iat) ==52 ) nom = "Te"
      if( mass(iat) ==53 ) nom = "I"
      if( mass(iat) ==54 ) nom = "Xe"

      if( mass(iat) ==55 ) nom = "Cs"
      if( mass(iat) ==56 ) nom = "Ba"

      if( mass(iat) ==72 ) nom = "Hf"
      if( mass(iat) ==73 ) nom = "Ta"
      if( mass(iat) ==74 ) nom = "W"
      if( mass(iat) ==75 ) nom = "Re"
      if( mass(iat) ==76 ) nom = "Os"
      if( mass(iat) ==77 ) nom = "Ir"
      if( mass(iat) ==78 ) nom = "Pt"
      if( mass(iat) ==79 ) nom = "Au"
      if( mass(iat) ==80 ) nom = "Hg"
      if( mass(iat) ==81 ) nom = "Tl"
      if( mass(iat) ==82 ) nom = "Pb"
      if( mass(iat) ==83 ) nom = "Bi"
      if( mass(iat) ==84 ) nom = "Po"
      if( mass(iat) ==85 ) nom = "At"
      if( mass(iat) ==86 ) nom = "Rn"

      write(11,"(a4,3F12.6)")adjustl(nom),(x(iat,k),k=1,3)
    else if( kcom == 1 ) then
      write(11,"(a4,3F12.6)")adjustl(nomcom(iat)),(x(iat,k),k=1,3)
    else
      write(*,*)"type de fichier inconnu"
      stop
    end if

  end do
  write(11,*)
  close(11)

  write(*,*)

end program xyz
