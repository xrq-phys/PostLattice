C=======================================================================
C Compute 2-Body Expectation Values With RDM
C Author: RuQing Xu, Dept. of Mod. Phys., USTC.
C NB: This program doesn't take any input containing information about
C     number of electrons. It is up to the user to multiply output by
C     a factor of N(N-1)/2.
C +-----------+---------------------------------+
C | Stop Code | Explanation                     |
C | 192 - 255 | Wrong Input                     |
C | 128 - 191 | Not Implemented, Not On The Way |
C |   1 - 127 | Not Implemented, On The Way     |
C +-----------+---------------------------------+
C=======================================================================

      PROGRAM RDMMean

        USE ISO_FORTRAN_ENV, ONLY: IOSTAT_END

        IMPLICIT NONE
        CHARACTER(32) :: hmfnm, dmfnm, opfnm
        REAL(8), ALLOCATABLE :: rdm(:, :, :, :), opcoeff(:, :, :, :)
        REAL(8) :: x, rdm_trace, opr_est
        INTEGER :: i, a, j, b
        INTEGER :: ihm, idm, iop, ios
        INTEGER :: norb, norb_l, nelec, ms2, 
     <             orbsym(100), isym(100), syml(100), symlz(100), 
     <             propbitlen(100), nprop(3)
        LOGICAL :: uhf, uhf_l, exist_t

        NAMELIST /FCI/ norb, nelec, ms2, orbsym, isym, syml, symlz,
     <                 propbitlen, nprop, uhf
        NAMELIST /OPR/ norb, uhf

        CALL GET_COMMAND_ARGUMENT(1, dmfnm)
        CALL GET_COMMAND_ARGUMENT(2, hmfnm)
        CALL GET_COMMAND_ARGUMENT(3, opfnm)

        IF (LEN_TRIM(dmfnm) == 0) dmfnm = 'RDM'
        IF (LEN_TRIM(hmfnm) == 0) hmfnm = 'FCIDUMP'
        IF (LEN_TRIM(opfnm) == 0) opfnm = 'SCDUMP'

        INQUIRE(FILE=dmfnm, EXIST=exist_t)
        IF (.NOT.exist_t) STOP 255
        INQUIRE(FILE=hmfnm, EXIST=exist_t)
        IF (.NOT.exist_t) STOP 255
        INQUIRE(FILE=opfnm, EXIST=exist_t)
        IF (.NOT.exist_t) STOP 255
        OPEN(NEWUNIT=idm, FILE=dmfnm)
        OPEN(NEWUNIT=ihm, FILE=hmfnm)
        OPEN(NEWUNIT=iop, FILE=opfnm)

        READ(ihm, FCI)
        IF (.NOT. uhf) STOP 191
        ALLOCATE(opcoeff(norb, norb, norb, norb))

        uhf_l = uhf
        norb_l = norb
        READ(iop, OPR)
        IF (uhf .NEQV. uhf_l .OR. norb /= norb_l) STOP 255

        rdm_trace = 0
        opr_est = 0
        opcoeff = 0
        DO
          READ(iop, *, IOSTAT=ios) x, i, a, j, b
          IF (ios == IOSTAT_END) EXIT
          IF (ANY((/i, a, j, b/) > norb) .OR. ANY((/i, a, j, b/) <= 0))
     <        CYCLE
          opcoeff(i, a, j, b) = x
        END DO

        DO
          READ(idm, *, IOSTAT=ios) a, b, i, j, x
          IF (ios == IOSTAT_END) EXIT
          IF (ANY((/i, a, j, b/) > norb) .OR. ANY((/i, a, j, b/) <= 0))
     <        CYCLE
          IF (i == a .AND. j == b) rdm_trace = rdm_trace + x
          opr_est = opr_est + opcoeff(i, a, j, b)*x
        END DO

        CLOSE(idm)
        CLOSE(ihm)
        CLOSE(iop)
        DEALLOCATE(opcoeff)
        WRITE(*, *) opr_est/rdm_trace

      END PROGRAM RDMMean

