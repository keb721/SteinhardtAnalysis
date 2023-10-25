PROGRAM MAIN

  USE CLUSTERINF

  IMPLICIT NONE

  CHARACTER(LEN=8)                               :: dcd_file = "traj.dcd"
  INTEGER                                        :: N, snapshots, snap, solid_threshold, l, atom, &
                                                    cluster, i, bg_clus, curr_clus, outfile, ierr
  LOGICAL                                        :: uc_present
  REAL                                           :: q6_threshold
  REAL, DIMENSION(:,:), ALLOCATABLE              :: positions
  REAL(KIND=REAL64), DIMENSION(:,:), ALLOCATABLE :: vec_ql_mat_rl !, vec_ql_mat_im
  REAL(KIND=REAL64), DIMENSION(3,3)              :: uc_mat
  REAL(KIND=REAL64)                              :: uc_len, uc_inv, distance_cutoff, mnq, theta, phi
  REAL(KIND=REAL64), DIMENSION(:), ALLOCATABLE   :: rl, im, angles
  INTEGER, ALLOCATABLE                           :: cluster_mat(:), neighbours(:) 
                                                 !! Dealing with vectors to be assigned later
  INTEGER, DIMENSION(:), ALLOCATABLE             :: con_mat
  INTEGER, DIMENSION(1)                          :: nshape



  l               = 6
  q6_threshold    = 0.5
  distance_cutoff = 1.43_REAL64 ! LAMMPS is working in LJ units not lattice units
  solid_threshold = 8           ! Need 8 neighbours to be solid


  outfile=99

  CALL READ_DCD_HEADER(dcd_file, N, snapshots, uc_present)

  
  ALLOCATE(positions(3, N))
  ALLOCATE(vec_ql_mat_rl(N, -l:l))
!  ALLOCATE(vec_ql_mat_im(N, -l:l))
  ALLOCATE(con_mat(N))
  ALLOCATE(rl(-l:l))
!  ALLOCATE(im(-l:l))

  nshape = SHAPE(rl)

  WRITE(*,*) "# INFO: Currently relies on a cubic simulation box"
  
  OPEN(UNIT=outfile, FILE='q6.txt', STATUS='REPLACE', FORM='FORMATTED', IOSTAT=ierr)  
  
  DO snap = 1,101

     WRITE(outfile, *) ' '
     WRITE(outfile,*) "TIME:", snap
     WRITE(outfile, *) ' '
     
     vec_ql_mat_rl  = 0.0_REAL64
!     vec_ql_mat_im  = 0.0_REAL64
     con_mat = 0

     
     CALL READ_DCD_SNAPSHOT(dcd_file, N, snap, uc_present, positions, uc_mat)
  
     uc_len = uc_mat(1,1)
     
     uc_inv = 1/uc_len

     ! Get positions etc.

     ! Unsophisticated measure of getting neighbours relies on using each atom individually

     DO atom = 1, N       
        neighbours             = GET_NEIGHS(positions, uc_len, uc_inv, distance_cutoff, atom)       
!        CALL ALL_VEC_Ql(l, positions, atom, uc_len, uc_inv, neighbours, rl, im)
        rl = VEC_Ql(l, positions, atom, uc_len, uc_inv, neighbours)
        vec_ql_mat_rl(atom, :) = rl
!        vec_ql_mat_im(atom, :) = im
     END DO
     
     DO atom = 1, N
        neighbours             = GET_NEIGHS(positions, uc_len, uc_inv, distance_cutoff, atom)       
        con_mat(atom) = CONNECTIONS(atom, neighbours, vec_ql_mat_rl, q6_threshold)
     END DO

     cluster_mat = CLUSTERNUMS(con_mat, solid_threshold, q6_threshold, vec_ql_mat_rl, positions, uc_len, uc_inv, distance_cutoff)
        


     DO atom = 1, N
!        mnq        = LOCAL_Ql(atom, neighbours, vec_ql_mat_rl, vec_ql_mat_im) 
        WRITE(outfile,*) atom, cluster_mat(atom)
    END DO
    
 END DO

CLOSE(outfile)

 
END PROGRAM MAIN
