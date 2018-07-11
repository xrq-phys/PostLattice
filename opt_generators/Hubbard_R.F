C=======================================================================
C Hubbard Parameter FCIDUMP Generator - R-Space
C Author: RuQing Xu, Dept. of Mod. Phys., USTC.
C +-----------+---------------------------------+
C | Stop Code | Explanation                     |
C | 192 - 255 | Wrong Input                     |
C | 128 - 191 | Not Implemented, Not On The Way |
C |   1 - 127 | Not Implemented, On The Way     |
C +-----------+---------------------------------+
C=======================================================================

      PROGRAM HubbardR

        IMPLICIT NONE
        CHARACTER(32) :: ifnm, ofnm
        INTEGER :: lattice(3), norb, ndim, ir, iw, ip, i1, j1, k1
        LOGICAL :: exist_t, uhf
        REAL(8) :: t, U

        NAMELIST /FCI/ norb, uhf
        NAMELIST /HUB/ lattice, t, U

        CALL GET_COMMAND_ARGUMENT(1, ifnm)
        CALL GET_COMMAND_ARGUMENT(2, ofnm)

        IF (LEN_TRIM(ifnm) == 0) ifnm = 'HUB'
        IF (LEN_TRIM(ofnm) == 0) ofnm = 'FCIDUMP'

        INQUIRE(FILE=ifnm, EXIST=exist_t)
        IF (.NOT.exist_t) STOP 255
        OPEN(NEWUNIT=iw, FILE=ofnm)
        OPEN(NEWUNIT=ir, FILE=ifnm)

        t = 1
        U = 1
        lattice = 0
        READ(ir, HUB)

        ndim = COUNT(lattice /= 0)
        IF (ANY(lattice(1:ndim) < 1)) STOP 254
        IF (ndim < 1) STOP 253

        uhf = .TRUE.
        norb = PRODUCT(lattice(1:ndim))*2
        IF (ndim /= 3) lattice(ndim+1:) = 1
        WRITE(iw, FCI)

        SELECT CASE(ndim)
        CASE(1)
          STOP 127
        CASE(2)
          DO i1 = 1, lattice(1)
          DO j1 = 1, lattice(2)
            ! Eigenvalue.
            ! The tiny negative value is to avoid initial distribution 
            ! having doublons.
            WRITE(iw, 6500) 1D-12*(-1)**(i1+j1), 
     <        cood_fold(ndim, lattice, (/i1, j1, 0/))*2, 0, 0, 0
            WRITE(iw, 6500) 1D-12*(-1)**(i1+j1+1), 
     <        cood_fold(ndim, lattice, (/i1, j1, 0/))*2-1, 0, 0, 0
            ! Kinetic
            WRITE(iw, 6500) -t, 
     <        cood_fold(ndim, lattice, (/i1, j1, 0/))*2,
     <        cood_fold(ndim, lattice, MOD((/i1+1, j1, 1/)+lattice-1, 
     <                                     lattice)+1)*2, 0, 0
            WRITE(iw, 6500) -t, 
     <        cood_fold(ndim, lattice, (/i1, j1, 0/))*2-1,
     <        cood_fold(ndim, lattice, MOD((/i1+1, j1, 1/)+lattice-1, 
     <                                     lattice)+1)*2-1, 0, 0
            WRITE(iw, 6500) -t, 
     <        cood_fold(ndim, lattice, (/i1, j1, 0/))*2,
     <        cood_fold(ndim, lattice, MOD((/i1-1, j1, 1/)+lattice-1, 
     <                                     lattice)+1)*2, 0, 0
            WRITE(iw, 6500) -t, 
     <        cood_fold(ndim, lattice, (/i1, j1, 0/))*2-1,
     <        cood_fold(ndim, lattice, MOD((/i1-1, j1, 1/)+lattice-1, 
     <                                     lattice)+1)*2-1, 0, 0
            WRITE(iw, 6500) -t, 
     <        cood_fold(ndim, lattice, (/i1, j1, 0/))*2,
     <        cood_fold(ndim, lattice, MOD((/i1, j1+1, 1/)+lattice-1, 
     <                                     lattice)+1)*2, 0, 0
            WRITE(iw, 6500) -t, 
     <        cood_fold(ndim, lattice, (/i1, j1, 0/))*2-1,
     <        cood_fold(ndim, lattice, MOD((/i1, j1+1, 1/)+lattice-1, 
     <                                     lattice)+1)*2-1, 0, 0
            WRITE(iw, 6500) -t, 
     <        cood_fold(ndim, lattice, (/i1, j1, 0/))*2,
     <        cood_fold(ndim, lattice, MOD((/i1, j1-1, 1/)+lattice-1, 
     <                                     lattice)+1)*2, 0, 0
            WRITE(iw, 6500) -t, 
     <        cood_fold(ndim, lattice, (/i1, j1, 0/))*2-1,
     <        cood_fold(ndim, lattice, MOD((/i1, j1-1, 1/)+lattice-1, 
     <                                     lattice)+1)*2-1, 0, 0
            ! Potential
            WRITE(iw, 6500) U, 
     <        cood_fold(ndim, lattice, (/i1, j1, 0/))*2-1,
     <        cood_fold(ndim, lattice, (/i1, j1, 0/))*2-1,
     <        cood_fold(ndim, lattice, (/i1, j1, 0/))*2,
     <        cood_fold(ndim, lattice, (/i1, j1, 0/))*2
          END DO
          END DO
        CASE(3)
          STOP 126
        END SELECT

        CLOSE(ir)
        CLOSE(iw)
 6500   FORMAT(ES20.13, 4I4)

      CONTAINS

        PURE FUNCTION cood_fold(ndim, lattice, indice) RESULT(idx_fold)

          INTEGER, INTENT(IN) :: ndim, lattice(3), indice(3)
          INTEGER :: idx_fold

          idx_fold = indice(1)
          IF (ndim >= 2) idx_fold = idx_fold + lattice(1)*(indice(2)-1)
          IF (ndim >= 3) idx_fold = idx_fold + 
     <                              PRODUCT(lattice(1:2))*(indice(3)-1)

        END FUNCTION cood_fold

      END PROGRAM HubbardR

