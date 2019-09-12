
C#################################################################

      subroutine ss_initialize()
      implicit none
      integer nprow, npcol, nprocs, k, myrank
      integer context
      common /cont/context
      integer dummy, m,n, nrhs
      double precision a
      character*1 type

      external blacs_pinfo, BLACS_GET, BLACS_GRIDINIT


      call BLACS_PINFO(myrank, nprocs)

C     Factor a number in two
      do k = int(sqrt(dble(nprocs))), 1, -1
        if (mod(nprocs, k) .eq. 0) goto 101
      enddo
 101  nprow = nprocs / k
      npcol = k

      CALL BLACS_GET(-1, 0, context)
      CALL BLACS_GRIDINIT(context, 'Row-major', nprow, npcol)

      if (myrank .ne. 0) then
         call ss_linsys(m, n, nrhs, a, dummy, a, dummy, type, dummy)
         call ss_finalize
         stop
      endif

      end

C#################################################################

      subroutine ss_finalize()
      implicit none
      integer myrank, nprocs
      integer context
      external blacs_gridexit, blacs_exit, blacs_pinfo, igebs2d
      common /cont/context

      call BLACS_PINFO(myrank, nprocs)
      if (myrank. eq. 0) call igebs2d(context, 'a', 'i' , 1, 1, 0, 1)

      CALL BLACS_GRIDEXIT(context)
      CALL BLACS_EXIT(0)
      end


C#################################################################

      subroutine ss_linsys(m, n, nrhs, matrix, gld_mat, rhs, gld_rhs, 
     &        type, info)
      implicit none
      integer m, n, nrhs, gld_mat, gld_rhs
      double precision matrix(gld_mat, *), rhs(gld_rhs, *)
      character*1 type
      integer info

      integer block_size_row, block_size_col
      parameter (block_size_row = 16)
      parameter (block_size_col = block_size_row)
      integer lld_mat, fld_mat, lld_rhs, fld_rhs, numroc
      integer nprow, npcol, myrow, mycol, myrank, nprocs
      integer cont, context
      integer desc_mat(9), desc_rhs(9)
      logical first
      save first, desc_mat, desc_rhs, lld_mat, fld_mat, lld_rhs, fld_rhs
      data first/.true./
      common /cont/context
      external igebs2d, igebr2d, numroc
      external BLACS_PINFO, BLACS_GRIDINFO, DESCINIT
      intrinsic max

      call BLACS_PINFO(myrank, nprocs)

C     broadcast 1 to signal continuation. Finalize sends 0.
 100  if (myrank.eq.0) then
         cont = 1
         call igebs2d(context, 'a', 'i' , 1, 1, cont, 1)
      else
         call igebr2d(context, 'a', 'i' , 1, 1, cont, 1, 0, 0)
         if (cont .ne. 1) return
      endif

C#
        if (first) then
C#         first = .false.
C#
         if (myrank .eq. 0) then
            call broadcast(context, m, n, nrhs, type) 
         else
            call receive(context, m, n, nrhs, type)
         endif

      call BLACS_GRIDINFO(context, nprow, npcol, myrow, mycol)
      lld_mat = max(1,NUMROC(m, block_size_row, myrow, 0, nprow))
      lld_rhs = lld_mat

      fld_mat = max(1,NUMROC(n, block_size_col, mycol, 0, npcol))
      fld_rhs = max(1,NUMROC(n, nrhs, mycol, 0, npcol))

C     construct descs for matrix and rhs
         call DESCINIT(desc_mat, m, n, block_size_row, block_size_col, 
     &           0, 0, context, lld_mat, info)
         if (info .ne. 0) return

         call DESCINIT(desc_rhs, m, nrhs, block_size_row, nrhs, 
     &           0, 0, context, lld_rhs, info)
         if (info .ne. 0) return
      endif

      call solveA(block_size_row, m, nrhs, desc_mat, matrix, gld_mat, 
     &        desc_rhs, rhs, gld_rhs, type, lld_mat, fld_mat,
     &        lld_rhs, fld_rhs, info)

      if (myrank .ne. 0) goto 100
      end


C#################################################################

      subroutine solveA(bsr, m, nrhs, desc_mat, matrix, gld_mat, 
     &        desc_rhs, rhs, gld_rhs, type, lld_mat, fld_mat,
     &        lld_rhs, fld_rhs, info)
      implicit none
      integer bsr
      integer desc_mat(9), desc_rhs(9)
      integer m, nrhs, gld_mat, gld_rhs
      integer lld_mat, lld_rhs, fld_mat, fld_rhs
      double precision matrix(gld_mat, *), rhs(gld_rhs, *)
      character*1 type
      integer info
      integer ipiv(lld_mat+bsr)
      double precision local(lld_mat, fld_mat), locrhs(lld_rhs, fld_rhs)
      external PDGESV, PDPOSV

      call scatter(desc_mat, matrix, gld_mat, local, lld_mat)
      call scatter(desc_rhs, rhs, gld_rhs, locrhs, lld_rhs)

      if (type .eq. 'g') CALL PDGESV(m, nrhs, local, 1, 1, 
     &        desc_mat, ipiv, locrhs, 1, 1, desc_rhs, info)
      
      if (type .eq. 's') CALL PDPOSV('u', m, nrhs, local, 1, 1, 
     &        desc_mat, locrhs, 1, 1, desc_rhs, info)
      if (info .ne. 0) return

      call gather(desc_rhs, rhs, gld_rhs, locrhs, lld_rhs)

      end
      
      
C#################################################################

      subroutine broadcast(context, m, n, nrhs, type)
      implicit none
      integer context, m,n,nrhs
      character*1 type
      external igebs2d

      call igebs2d(context, 'a', 'i' , 1, 1, m, 1)
      call igebs2d(context, 'a', 'i' , 1, 1, n, 1)
      call igebs2d(context, 'a', 'i' , 1, 1, nrhs, 1)
      if (type .eq. 'g') call igebs2d(context, 'a', 'i' , 1, 1, 1, 1)
      if (type .eq. 's') call igebs2d(context, 'a', 'i' , 1, 1, 2, 1)

      end

C#################################################################
      subroutine receive(context, m, n, nrhs, type)
      implicit none
      integer context, m,n,nrhs
      character*1 type
      integer typ
      external igebr2d

      call igebr2d(context, 'a', 'i' , 1, 1, m, 1, 0, 0)
      call igebr2d(context, 'a', 'i' , 1, 1, n, 1, 0, 0)
      call igebr2d(context, 'a', 'i' , 1, 1, nrhs, 1, 0, 0)

      call igebr2d(context, 'a', 'i' , 1, 1, typ, 1, 0, 0)
      if (typ .eq. 1) type = 'g'
      if (typ .eq. 2) type = 's'
      end

C#################################################################

      subroutine scatter(desc, global, gld, local, lld)
      implicit none
      
      integer desc(*), gld, lld
      double precision global(gld,*), local(lld,*)

      integer mrow, mcol, block_size_row, block_size_col
      integer i, j, context, row, col, stride_row, stride_col
      integer nprow, npcol, elem_r, elem_c, myrow, mycol
      external BLACS_GRIDINFO, dgesd2d, dgerv2d
      intrinsic min

      context        = desc(2)
      mrow           = desc(3)
      mcol           = desc(4)
      block_size_row = desc(5)
      block_size_col = desc(6)

      call BLACS_GRIDINFO(context, nprow, npcol, myrow, mycol)

      stride_row = block_size_row * nprow
      stride_col = block_size_col * npcol

C$OMP PARALLEL DO
      do i = 0, mrow-1, block_size_row
        row = mod(i, stride_row) / block_size_row
        elem_r = min(block_size_row, mrow - i)

C$OMP PARALLEL DO
        do j = 0, mcol-1, block_size_col
          col = mod(j, stride_col) / block_size_col
          elem_c = min(block_size_col, mcol - j)

          if ((myrow .eq. 0) .and. (mycol .eq. 0)) then
             call DGESD2D(context, elem_r, elem_c, 
     &               global(i+1, j+1), gld, row, col)
          endif
          
          if ((myrow .eq. row) .and. (mycol .eq. col)) then
             call DGERV2D(context, elem_r, elem_c, 
     &               local((i / stride_row) * block_size_row + 1, 
     &               (j / stride_col) * block_size_col + 1), lld, 0, 0)
          endif
        enddo
      enddo

      end


C#################################################################

      subroutine gather(desc, global, gld, local, lld)
      implicit none
      
      integer desc(*), gld, lld
      double precision global(gld,*), local(lld,*)

      integer mrow, mcol, block_size_row, block_size_col
      integer i, j, context, row, col, stride_row, stride_col
      integer nprow, npcol, elem_r, elem_c, myrow, mycol
      external BLACS_GRIDINFO, dgesd2d, dgerv2d
      intrinsic min

      context        = desc(2)
      mrow           = desc(3)
      mcol           = desc(4)
      block_size_row = desc(5)
      block_size_col = desc(6)

      call BLACS_GRIDINFO(context, nprow, npcol, myrow, mycol)

      stride_row = block_size_row * nprow
      stride_col = block_size_col * npcol

C$OMP PARALLEL DO
      do i = 0, mrow-1, block_size_row
        row = mod(i, stride_row) / block_size_row
        elem_r = min(block_size_row, mrow - i)

C$OMP PARALLEL DO
        do j = 0, mcol-1, block_size_col
          col = mod(j, stride_col) / block_size_col
          elem_c = min(block_size_col, mcol - j)

          if ((myrow .eq. row) .and. (mycol .eq. col)) then
             call DGESD2D(context, elem_r, elem_c, 
     &               local((i / stride_row) * block_size_row + 1, 
     &               (j / stride_col) * block_size_col + 1), lld, 0, 0)
          endif

          if ((myrow .eq. 0) .and. (mycol .eq. 0)) then
             call DGERV2D(context, elem_r, elem_c, 
     &               global(i+1, j+1), gld, row, col)
          endif          
        enddo
      enddo

      end
C########################################################################
