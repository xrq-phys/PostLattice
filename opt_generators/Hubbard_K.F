C=======================================================================
C Hubbard Parameter FCIDUMP Generator - K-Space
C Author: RuQing Xu, Dept. of Mod. Phys., USTC.
C +-----------+---------------------------------+
C | Stop Code | Explanation                     |
C | 192 - 255 | Wrong Input                     |
C | 128 - 191 | Not Implemented, Not On The Way |
C |   1 - 127 | Not Implemented, On The Way     |
C +-----------+---------------------------------+
C=======================================================================

      PROGRAM HubbardK

        IMPLICIT NONE
        CHARACTER(32) :: ifnm, ofnm, pfnm
        INTEGER :: lattice(3), norb, ndim, ir, iw, ip
        INTEGER :: i1, j1, k1, i2, j2, k2, i3, j3, k3, i4, j4, k4
        LOGICAL :: exist_t, uhf
        REAL(8) :: ukv(3), t, U
        REAL(8), PARAMETER :: pi = 3.1415926535897932384626433832795D+00

        NAMELIST /FCI/ norb, uhf
        NAMELIST /OPR/ norb, uhf
        NAMELIST /HUB/ lattice, t, U

        CALL GET_COMMAND_ARGUMENT(1, ifnm)
        CALL GET_COMMAND_ARGUMENT(2, ofnm)
        CALL GET_COMMAND_ARGUMENT(3, pfnm)

        IF (LEN_TRIM(ifnm) == 0) ifnm = 'HUB'
        IF (LEN_TRIM(ofnm) == 0) ofnm = 'FCIDUMP'
        IF (LEN_TRIM(pfnm) == 0) pfnm = 'SCDUMP'

        INQUIRE(FILE=ifnm, EXIST=exist_t)
        IF (.NOT.exist_t) STOP 255
        OPEN(NEWUNIT=iw, FILE=ofnm)
        OPEN(NEWUNIT=ir, FILE=ifnm)
        OPEN(NEWUNIT=ip, FILE=pfnm)

        t = 1
        U = 1
        lattice = 0
        READ(ir, HUB)

        ndim = COUNT(lattice /= 0)
        IF (ANY(lattice(1:ndim) < 1)) STOP 254
        IF (ndim < 1) STOP 253

        uhf = .TRUE.
        norb = PRODUCT(lattice(1:ndim))*2
        WRITE(iw, FCI)
        WRITE(ip, OPR)

        ukv = 0
        ukv(1:ndim) = 2.D+00*pi/lattice(1:ndim)

        SELECT CASE(ndim)
        CASE(1)
          STOP 127
        CASE(2)
          DO i1 = 1, lattice(1)
          DO j1 = 1, lattice(2)
          DO i2 = 1, lattice(1)
          DO j2 = 1, lattice(2)
          DO i3 = 1, lattice(1)
          DO j3 = 1, lattice(2)
            ! <34|12>
            i4 = MOD(i1 + i2 - i3 - 1 + 2*lattice(1), lattice(1)) + 1
            j4 = MOD(j1 + j2 - j3 - 1 + 2*lattice(2), lattice(2)) + 1
            WRITE(iw, 6500) U/PRODUCT(lattice(1:ndim)),
     <                      cood_fold(ndim, lattice, (/i1, j1, 0/))*2,
     <                      cood_fold(ndim, lattice, (/i3, j3, 0/))*2,
     <                      cood_fold(ndim, lattice, (/i2, j2, 0/))*2-1,
     <                      cood_fold(ndim, lattice, (/i4, j4, 0/))*2-1
            WRITE(iw, 6500) U/PRODUCT(lattice(1:ndim)),
     <                      cood_fold(ndim, lattice, (/i1, j1, 0/))*2-1,
     <                      cood_fold(ndim, lattice, (/i3, j3, 0/))*2-1,
     <                      cood_fold(ndim, lattice, (/i2, j2, 0/))*2,
     <                      cood_fold(ndim, lattice, (/i4, j4, 0/))*2
          IF (MOD(i1 + i2 - 2, lattice(1)) == 0) THEN
            WRITE(ip, 6500) 1.,
     <                      cood_fold(ndim, lattice, (/i1, j1, 0/))*2,
     <                      cood_fold(ndim, lattice, (/i3, j3, 0/))*2,
     <                      cood_fold(ndim, lattice, (/i2, j2, 0/))*2-1,
     <                      cood_fold(ndim, lattice, (/i4, j4, 0/))*2-1
            WRITE(ip, 6500) 1.,
     <                      cood_fold(ndim, lattice, (/i1, j1, 0/))*2-1,
     <                      cood_fold(ndim, lattice, (/i3, j3, 0/))*2-1,
     <                      cood_fold(ndim, lattice, (/i2, j2, 0/))*2,
     <                      cood_fold(ndim, lattice, (/i4, j4, 0/))*2
          END IF
          END DO
          END DO
          END DO
          END DO
          END DO
          END DO
          DO i1 = 1, lattice(1)
          DO j1 = 1, lattice(2)
            ! Eigenvalues and <3|1>
            WRITE(iw, 6500) -2*t*SUM(COS(ukv(1:2)*(/i1-1, j1-1/))), 
     <                      cood_fold(ndim, lattice, (/i1, j1, 0/))*2,
     <                      0, 0, 0
            WRITE(iw, 6500) -2*t*SUM(COS(ukv(1:2)*(/i1-1, j1-1/))), 
     <                      cood_fold(ndim, lattice, (/i1, j1, 0/))*2-1,
     <                      0, 0, 0
            WRITE(iw, 6500) -2*t*SUM(COS(ukv(1:2)*(/i1-1, j1-1/))), 
     <                      cood_fold(ndim, lattice, (/i1, j1, 0/))*2,
     <                      cood_fold(ndim, lattice, (/i1, j1, 0/))*2, 
     <                      0, 0
            WRITE(iw, 6500) -2*t*SUM(COS(ukv(1:2)*(/i1-1, j1-1/))), 
     <                      cood_fold(ndim, lattice, (/i1, j1, 0/))*2-1,
     <                      cood_fold(ndim, lattice, (/i1, j1, 0/))*2-1,
     <                      0, 0
          END DO
          END DO
        CASE(3)
          STOP 126
        END SELECT

        CLOSE(ir)
        CLOSE(iw)
        CLOSE(ip)
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

      END PROGRAM HubbardK

