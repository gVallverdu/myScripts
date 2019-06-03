! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
! Germain Vallverdu - germain.vallverdu@gmail.com
!
! 26/11/2008
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

MODULE topology
  implicit none
  ! FICHIER TOPOLOGY : prmtop 
  character(len=100)::top
  
  ! nombre d'atome et de residue
  integer::nat,nres
  ! ires(iat) = numero du residue de l'atome iat
  integer,dimension(:),allocatable::ires 
  ! pointres(jres) = numero du premier atome du residue jres
  integer,dimension(:),allocatable::pointres
  ! mass(iat) = masse de l'atome iat
  ! q(iat) = charge de l'atome iat
  double precision,dimension(:),allocatable::mass,q
  ! at_name(iat) = nom de l'atome iat
  ! res_name(iat) = nom du residue de l'atome iat
  character(len=4),dimension(:),allocatable::at_name,res_name
  ! reslist(jres) = nom du residue jres
  character(len=4),dimension(:),allocatable::reslist
  
  CONTAINS
    SUBROUTINE lectop
      implicit none
      integer::i,iat,jres
      integer,dimension(20)::ipt
      character(len=100)::ligne
      
      open(10, file=trim(top), action="read")
      
      ! lecture des pointeurs
      Do
        read(10,"(a)")ligne
        If( index(ligne,"%FLAG POINTERS") /= 0 ) exit
      Enddo
      read(10,*)
      read(10,"(10i8)")(ipt(i),i=1,20)
      
      nat  = ipt(1)
      nres = ipt(12)
      ! dimensions des tableaux
      allocate( mass(nat), q(nat), at_name(nat), res_name(nat) )
      allocate( ires(nat) )     
      allocate( reslist(nres), pointres(nres) )
      
      ! lecture des noms d'atomes
      Do
        read(10,"(a)")ligne
        If( index(ligne,"%FLAG ATOM_NAME") /= 0 ) exit
      Enddo
      read(10,*)
      read(10,"(20a4)")(at_name(iat),iat=1,nat)
      
      ! lecture des charges
      Do
        read(10,"(a)") ligne
        If (index(ligne,"%FLAG CHARGE") /= 0 ) exit
      Enddo
      read(10,*)
      read(10,"(5E16.8)") (q(iat),iat=1,nat)
      
      ! lecture des masses
      Do
        read(10,"(a)") ligne
        If (index(ligne,"%FLAG MASS") /= 0 ) exit
      Enddo
      read(10,*)
      read(10,"(5E16.8)") (mass(iat),iat=1,nat)
      
      ! lecture de la listes d'acides aminés
      Do
        read(10,"(a)")ligne
        If( index(ligne,"%FLAG RESIDUE_LABEL") /= 0 ) exit
      Enddo
      read(10,*)
      read(10,"(20a4)")(reslist(jres),jres=1,nres)
      
      ! lecture des pointeurs pour les résidues
      Do
        read(10,"(a)")ligne
        If( index(ligne,"%FLAG RESIDUE_POINTER") /= 0 ) exit
      Enddo
      read(10,*)
      read(10,"(10I8)")(pointres(jres),jres=1,nres)
      
      close(10)
      
      ! nom des résidus et numéro en fonction du numero d'atome
      jres = 1
      Do iat=1,nat
        If( iat == pointres(jres+1) ) jres = jres + 1 
        If(jres == nres )exit
        ires(iat) = jres
        res_name(iat) = reslist(jres)
      Enddo
      ! le derniers residus
      Do iat=pointres(jres),nat
        ires(iat)=jres
        res_name(iat)=reslist(jres)
      Enddo
      
      write(*,"('nbre d atome     : ',i7)")nat
      write(*,"('nbre de residues : ',i7)")nres  
      write(*,"('top lu           : ',a)")trim(top)
      write(*,*)
    
    END SUBROUTINE lectop

END MODULE topology

! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

