program main

!#define MKL

use funmod2
#ifndef MKL
use mpi
#endif

implicit none
intrinsic count
intrinsic len
intrinsic product
intrinsic cmplx

include 'amica15_header.f90'



#ifdef MKL
include 'mpif.h'
include 'mkl_vml.f90'
#endif


!-------------- INITIALIZE MPI ----------------

!--- init mpi and get myrank
call MPI_INIT(ierr)
call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)
call MPI_COMM_SIZE(MPI_COMM_WORLD,tot_procs,ierr)

!--- get number of procs on this host
call MPI_GET_PROCESSOR_NAME(host_name,name_len,ierr)
print *, myrank+1, 'processor name = ', trim(host_name); call flush(6)
!--- turn hostname into a number (hash)
host_num = 0
do i = 1,name_len
   host_num = host_num + iachar(host_name(i:i)) * 31**(name_len-i)
end do
host_num = abs(host_num)
print *, myrank+1, 'host_num = ', host_num; call flush(6)

call MPI_COMM_SPLIT(MPI_COMM_WORLD,host_num,0,node_comm,ierr)
call MPI_COMM_RANK(node_comm,node_rank,ierr)
call MPI_COMM_SIZE(node_comm,node_procs,ierr)

print *, 'This is MPI process', myrank+1, 'of', tot_procs, &
     '; I am process', node_rank+1, 'of', node_procs, 'on node: ', trim(host_name)
call flush(6)
call MPI_BARRIER(MPI_COMM_WORLD,ierr)

!--- set up MPI communicators
call MPI_COMM_SPLIT(MPI_COMM_WORLD,node_rank,0,seg_comm,ierr)

!--- make unused processes wait
if (node_rank > 0) then
   call MPI_BCAST(tot_procs,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
   call MPI_FINALIZE(ierr)
   stop
end if
call MPI_COMM_RANK(seg_comm,seg_rank,ierr)
call MPI_COMM_SIZE(seg_comm,seg_nodes,ierr)
print *, myrank+1, ' : node root process', seg_rank+1, 'of', seg_nodes; call flush(6)
call MPI_BARRIER(seg_comm,ierr)


!--------------- GET ARGS -----------------

if (myrank == 0) then
   print "(a)", 'Processing arguments ...'
   call flush(6)

   call get_cmd_args ! process input arguments from command line and file

   call date_and_time(date,time) ! tag output with (hopefully) unique month/day/hr/min/sec of execution
   call system_clock(c0)
   if (outdirparam == '') then
      call system('mkdir '//outdir//date(5:8)//'_'//time(1:6))
      open(unit=20,file=outdir//date(5:8)//'_'//time(1:6)//'/'//printoutfile,status='replace')
   else
      call system('mkdir '//trim(outdirparam)) ! Make the output directory specified in the parameter file
      open(unit=20,file=trim(outdirparam)//'/'//printoutfile,status='replace')
   end if
   print *, 'output directory = ', trim(outdirparam); call flush(6)

end if


!--------------- INITIALIZE OPENMP ----------------

call MPI_BCAST(max_thrds,1,MPI_INTEGER,0,seg_comm,ierr)

!if ((node_procs == 1) .and. (seg_nodes == 1)) then
   num_thrds = max_thrds
!else
!   num_thrds = min(max_thrds,node_procs)
!end if
print *, myrank+1, ': setting num_thrds to ', num_thrds, ' ...'; call flush(6)
call omp_set_num_threads(num_thrds)

!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP & PRIVATE(thrdnum) SHARED(num_thrds_used)
thrdnum = omp_get_thread_num()
if (thrdnum == 0) then
   num_thrds_used = omp_get_num_threads()
end if
!$OMP END PARALLEL

num_thrds = num_thrds_used
print *, myrank+1, ': using', num_thrds, 'threads.'
call flush(6)
call MPI_BARRIER(seg_comm,ierr)

!--- gather num thrds used by all procs on root
if (myrank == 0) then
   allocate(node_thrds(seg_nodes),stat=ierr); call tststat(ierr)
end if
call MPI_GATHER(num_thrds,1,MPI_INTEGER,node_thrds,1,MPI_INTEGER,0,seg_comm,ierr)
if (myrank == 0) then
   print *, myrank+1, ': node_thrds = ', node_thrds
end if

!---------------- GET RECL MULTIPLIER ----------------
if (myrank == 0) then
   call bytes_in_rec(nbyte)
   print *, myrank+1, ': REAL nbyte = ', nbyte; call flush(6)
end if

!--------------- BROADCAST VARIABLES ----------------

call MPI_BCAST(num_files,1,MPI_INTEGER,0,seg_comm,ierr)
if (seg_rank > 0) then
   allocate(infile(num_files))
   allocate(num_samples(num_files))
   allocate(field_dim(num_files))
end if
do k = 1,num_files
   call MPI_BCAST(infile(k),500,MPI_CHARACTER,0,seg_comm,ierr)
end do
call MPI_BCAST(num_samples,num_files,MPI_INTEGER,0,seg_comm,ierr)
call MPI_BCAST(data_dim,1,MPI_INTEGER,0,seg_comm,ierr)
call MPI_BCAST(field_dim,num_files,MPI_INTEGER,0,seg_comm,ierr)

call MPI_BCAST(filter_length,1,MPI_INTEGER,0,seg_comm,ierr)

call MPI_BCAST(do_opt_block,1,MPI_LOGICAL,0,seg_comm,ierr)

call MPI_BCAST(do_sphere,1,MPI_LOGICAL,0,seg_comm,ierr)
call MPI_BCAST(do_mean,1,MPI_LOGICAL,0,seg_comm,ierr)
call MPI_BCAST(load_sphere,1,MPI_LOGICAL,0,seg_comm,ierr)
call MPI_BCAST(load_mean,1,MPI_LOGICAL,0,seg_comm,ierr)
call MPI_BCAST(update_A,1,MPI_LOGICAL,0,seg_comm,ierr)
call MPI_BCAST(update_c,1,MPI_LOGICAL,0,seg_comm,ierr)
call MPI_BCAST(update_mu,1,MPI_LOGICAL,0,seg_comm,ierr)
call MPI_BCAST(update_alpha,1,MPI_LOGICAL,0,seg_comm,ierr)
call MPI_BCAST(update_gm,1,MPI_LOGICAL,0,seg_comm,ierr)
call MPI_BCAST(update_beta,1,MPI_LOGICAL,0,seg_comm,ierr)
call MPI_BCAST(dorho,1,MPI_LOGICAL,0,seg_comm,ierr)
call MPI_BCAST(load_A,1,MPI_LOGICAL,0,seg_comm,ierr)
call MPI_BCAST(load_c,1,MPI_LOGICAL,0,seg_comm,ierr)
call MPI_BCAST(load_mu,1,MPI_LOGICAL,0,seg_comm,ierr)
call MPI_BCAST(load_alpha,1,MPI_LOGICAL,0,seg_comm,ierr)
call MPI_BCAST(load_gm,1,MPI_LOGICAL,0,seg_comm,ierr)
call MPI_BCAST(load_beta,1,MPI_LOGICAL,0,seg_comm,ierr)
call MPI_BCAST(load_rho,1,MPI_LOGICAL,0,seg_comm,ierr)

call MPI_BCAST(outdirparam,500,MPI_CHARACTER,0,seg_comm,ierr)
call MPI_BCAST(indirparam,500,MPI_CHARACTER,0,seg_comm,ierr)

call MPI_BCAST(write_LLt,1,MPI_LOGICAL,0,seg_comm,ierr)
call MPI_BCAST(write_nd,1,MPI_LOGICAL,0,seg_comm,ierr)

call MPI_BCAST(block_size,1,MPI_INTEGER,0,seg_comm,ierr)

call MPI_BCAST(blk_min,1,MPI_INTEGER,0,seg_comm,ierr)
call MPI_BCAST(blk_max,1,MPI_INTEGER,0,seg_comm,ierr)
call MPI_BCAST(blk_step,1,MPI_INTEGER,0,seg_comm,ierr)

call MPI_BCAST(do_reject,1,MPI_LOGICAL,0,seg_comm,ierr)
call MPI_BCAST(maxrej,1,MPI_INTEGER,0,seg_comm,ierr)
call MPI_BCAST(rejstart,1,MPI_INTEGER,0,seg_comm,ierr)
call MPI_BCAST(rejint,1,MPI_INTEGER,0,seg_comm,ierr)
call MPI_BCAST(rejsig,1,MPI_DOUBLE_PRECISION,0,seg_comm,ierr)
call MPI_BCAST(invsigmax,1,MPI_DOUBLE_PRECISION,0,seg_comm,ierr)
call MPI_BCAST(load_rej,1,MPI_LOGICAL,0,seg_comm,ierr)

call MPI_BCAST(chpdfstart,1,MPI_INTEGER,0,seg_comm,ierr)
call MPI_BCAST(chpdfint,1,MPI_INTEGER,0,seg_comm,ierr)
call MPI_BCAST(numchpdf,1,MPI_INTEGER,0,seg_comm,ierr)

!call MPI_BCAST(decwindow,1,MPI_INTEGER,0,seg_comm,ierr)

call MPI_BCAST(do_newton,1,MPI_LOGICAL,0,seg_comm,ierr)
call MPI_BCAST(newt_start,1,MPI_INTEGER,0,seg_comm,ierr)
call MPI_BCAST(newt_ramp,1,MPI_INTEGER,0,seg_comm,ierr)
call MPI_BCAST(newtrate,1,MPI_DOUBLE_PRECISION,0,seg_comm,ierr)

call MPI_BCAST(dft_length,1,MPI_INTEGER,0,seg_comm,ierr)
call MPI_BCAST(doscaling,1,MPI_LOGICAL,0,seg_comm,ierr)
call MPI_BCAST(scalestep,1,MPI_INTEGER,0,seg_comm,ierr)
call MPI_BCAST(writestep,1,MPI_INTEGER,0,seg_comm,ierr)
call MPI_BCAST(max_iter,1,MPI_INTEGER,0,seg_comm,ierr)
call MPI_BCAST(num_models,1,MPI_INTEGER,0,seg_comm,ierr)
call MPI_BCAST(nbyte,1,MPI_INTEGER,0,seg_comm,ierr)
call MPI_BCAST(pdftype,1,MPI_INTEGER,0,seg_comm,ierr)

call MPI_BCAST(lrate,1,MPI_DOUBLE_PRECISION,0,seg_comm,ierr)
call MPI_BCAST(lrate0,1,MPI_DOUBLE_PRECISION,0,seg_comm,ierr)
call MPI_BCAST(rholrate0,1,MPI_DOUBLE_PRECISION,0,seg_comm,ierr)
call MPI_BCAST(minrho,1,MPI_DOUBLE_PRECISION,0,seg_comm,ierr)
call MPI_BCAST(maxrho,1,MPI_DOUBLE_PRECISION,0,seg_comm,ierr)
call MPI_BCAST(rho0,1,MPI_DOUBLE_PRECISION,0,seg_comm,ierr)
call MPI_BCAST(minlrate,1,MPI_DOUBLE_PRECISION,0,seg_comm,ierr)
call MPI_BCAST(lratefact,1,MPI_DOUBLE_PRECISION,0,seg_comm,ierr)
call MPI_BCAST(rholratefact,1,MPI_DOUBLE_PRECISION,0,seg_comm,ierr)

!--- set up the random number generator for this node
call system_clock(c1)
call random_seed(PUT = c1 * (myrank+1) * (seed+myrank+1))
!call DRANDINITIALIZE(1,1,(myrank+1)*(c1/tot_procs + 1),lseed,state,lstate,info)


!-------------------- GET THE DATA ------------------------

!--- calculate file/offsets for the segment nodes
allocate(seg_list(seg_nodes,2,3)) ! start/stop, file/sample/coord
if (myrank == 0) then
   print *, 'getting segment list ...'; call flush(6)
   call get_seg_list
end if
call MPI_BCAST(seg_list,seg_nodes*2*3,MPI_INTEGER,0,seg_comm,ierr)
call MPI_BCAST(all_blks,1,MPI_INTEGER,0,seg_comm,ierr)
numgoodsum = all_blks

!-------------- load the data and build samplist for this node -------------
call get_data
print *, myrank+1, ': data = ', dataseg(1)%data(1:2,1); call flush(6)
!print *, myrank+1, ': numsegs = ', numsegs; call flush(6)

allocate(blk_size(numsegs),stat=ierr); call tststat(ierr)
do seg = 1,numsegs
   fieldsize = dataseg(seg)%lastdim
   !print *, myrank+1, 'seg ', seg, ': fieldsize = ', fieldsize; call flush(6)
   if (fieldsize < num_thrds) then
      print *, 'Error: minimum segment size too small'; call flush(6)
      quittmp = 1
      exit
   else
      quittmp = 0
   end if
   blk_size(seg) = min(fieldsize,block_size)
end do
call MPI_ALLREDUCE(quittmp,quit,1,MPI_INTEGER,MPI_MAX,seg_comm,ierr)
if (quit == 1) then
   stop
end if

!---------------------------- get the mean --------------------------------
nx = data_dim
allocate( work(5*nx*nx),stat=ierr); call tststat(ierr)
allocate(mean(nx))
mean = dble(0.0)

if (load_mean) then
   if (seg_rank == 0) then
      print *, myrank+1, ': reading mean from: ', trim(adjustl(indirparam))//'/'//'mean'; call flush(6)
      open(unit=15,file=trim(adjustl(indirparam))//'/'//'mean',access='direct',recl=2*nbyte*nx)
      read(15,rec=1) mean
      close(15)
      !print *, 'mean = ', mean(1:min(5,nw)); call flush(6)
   end if
elseif (do_mean) then

   if (myrank == 0) then
      print *, 'getting the mean ...'; call flush(6)
   end if
   allocate(meantmp(nx))
   meantmp = dble(0.0)
   call DSCAL(nx,dble(0.0),meantmp,1)
   cnt1 = 0
   do seg = 1,numsegs
      fnum = dataseg(seg)%filenum
      fieldsize = dataseg(seg)%lastdim
      ldim = dataseg(seg)%lastdim

      call DAXPY(nx,dble(1.0),sum(dataseg(seg)%data(:,1:ldim),2),1,meantmp,1)
      cnt1 = cnt1 + fieldsize
   end do

   call DSCAL(nx,dble(0.0),mean,1)
   cnt = 0
   call MPI_REDUCE(meantmp,mean,nx,MPI_DOUBLE_PRECISION,MPI_SUM,0,seg_comm,ierr)
   call MPI_REDUCE(cnt1,cnt,1,MPI_INTEGER,MPI_SUM,0,seg_comm,ierr)

   if (myrank == 0) then
      call DSCAL(nx,dble(1.0)/dble(cnt),mean,1)
      print *, ' mean = ', mean(1:3); call flush(6)
   end if
end if


call MPI_BCAST(mean,nx,MPI_DOUBLE_PRECISION,0,seg_comm,ierr)

!--- subtract the mean
if (myrank == 0) then
   print *, 'subtracting the mean ...'; call flush(6)
end if

do seg = 1,numsegs
   fnum = dataseg(seg)%filenum         
   fieldsize = dataseg(seg)%lastdim

   num_blocks = fieldsize /  blk_size(seg)
   lastblocksize = mod(fieldsize,blk_size(seg))

   !print *, 'num_blocks = ', num_blocks, ' lastblocksize = ', lastblocksize, ' fieldsize = ', fieldsize; call flush(6);

   do k = 1,num_blocks
      bstrt = (k-1)*blk_size(seg) + 1
      bstp = bstrt + blk_size(seg) - 1
      call DAXPY(nx*blk_size(seg),dble(-1.0),spread(mean,2,blk_size(seg)),1,dataseg(seg)%data(:,bstrt:bstp),1)
   end do
   if (lastblocksize .ne. 0) then
      bstrt = num_blocks*blk_size(seg) + 1
      bstp = fieldsize
      call DAXPY(nx*lastblocksize,dble(-1.0),spread(mean,2,lastblocksize),1,dataseg(seg)%data(:,bstrt:bstp),1)
   end if
end do



!------------------------ sphere the data -------------------------------
allocate(S(nx,nx)); S = dble(0.0)
allocate(Stmp(nx,nx)); Stmp = dble(0.0)
!--- get the covariance matrix
if (myrank == 0) then
   print *, 'getting the covariance matrix ...'; call flush(6)
end if

call DSCAL(nx*nx,dble(0.0),Stmp,1)

cnt1 = 0
do seg = 1,numsegs
   fnum = dataseg(seg)%filenum
   fieldsize = dataseg(seg)%lastdim
   ldim = dataseg(seg)%lastdim
   !print *, 'seg ', seg, ' fnum = ', fnum, ' fieldsize = ', fieldsize, ' ldim = ', ldim; call flush(6)
   num_blocks = fieldsize / blk_size(seg)
   lastblocksize = mod(fieldsize,blk_size(seg))
   !print *, 'seg ', seg, ' num_blocks = ', num_blocks, ' lastblocksize = ', lastblocksize; call flush(6)
   do k = 1,num_blocks
      bstrt = (k-1)*blk_size(seg) + 1
      bstp = bstrt + blk_size(seg) - 1
      call DSYRK('L','N',nx,blk_size(seg),dble(1.0),dataseg(seg)%data(:,bstrt:bstp),nx,dble(1.0),Stmp,nx)
   end do
   if (lastblocksize .ne. 0) then
      bstrt = num_blocks*blk_size(seg) + 1
      bstp = ldim
      call DSYRK('L','N',nx,lastblocksize,dble(1.0),dataseg(seg)%data(:,bstrt:bstp),nx,dble(1.0),Stmp,nx)
   end if
   cnt1 = cnt1 + fieldsize

end do
call DSCAL(nx*nx,dble(0.0),S,1)
cnt = 0
call MPI_REDUCE(Stmp,S,nx*nx,MPI_DOUBLE_PRECISION,MPI_SUM,0,seg_comm,ierr)
call MPI_REDUCE(cnt1,cnt,1,MPI_INTEGER,MPI_SUM,0,seg_comm,ierr)
if (myrank == 0) then
   print *, 'cnt = ', cnt; call flush(6)
   call DSCAL(nx*nx,dble(1.0)/dble(cnt),S,1)
end if

call MPI_BCAST(S,nx*nx,MPI_DOUBLE_PRECISION,0,seg_comm,ierr)
allocate(eigs(nx))
allocate(eigv(nx))
lwork = 10*nx*nx
print *, 'doing eig nx = ', nx, ' lwork = ', lwork; call flush(6)
call DCOPY(nx*nx,S,1,Stmp,1)

call DSYEV('V','L',nx,Stmp,nx,eigs,work,lwork,info)

if (info == 0) then
   print *, 'minimum eigenvalues = ', eigs(1:min(nx/2,3))
   write(20,*) 'minimum eigenvalues = ', eigs(1:min(nx/2,3))
   print *, 'maximum eigenvalues = ', eigs(nx:(nx-min(nx/2,3)+1):-1)
   write(20,*) 'maximum eigenvalues = ', eigs(nx:(nx-min(nx/2,3)+1):-1)
   call flush(6)
else
   print "(a)", 'Error doing eigenvalue decomposition!!!'
   write(20,"(a)") 'Error doing eigenvalue decomposition!!!'
   call flush(6)
end if

numeigs = min(pcakeep,count(eigs > mineig))
print *, 'num eigs kept = ', numeigs; call flush(6)
write(20,*) 'num eigs kept = ', numeigs; call flush(6)


if (load_sphere) then

   if (seg_rank == 0) then         
      print *, myrank+1, ': reading S from: ', trim(adjustl(indirparam))//'/'//'S'; call flush(6)
      open(unit=15,file=trim(adjustl(indirparam))//'/'//'S',access='direct',recl=2*nbyte*nx*nx)
      read(15,rec=1) S
      close(15)
      !print *, 'initial S = ', S(1:min(5,nx*nx)); call flush(6)

      call DCOPY(nx*nx,S,1,Stmp,1)
      lwork = 5*nx*nx
      allocate( wr(nx),stat=ierr); call tststat(ierr)
      call DGEQRF(nx,nx,Stmp,nx,wr,work,lwork,info)
      if (info < 0) then
         print *, 'error getting QR!!!'
      end if
      sldet = dble(0.0)
      do i = 1,numeigs
         sldet = sldet + log(Stmp(i,i))
      end do
      deallocate(wr)
   end if
else
   !--- get the sphering matrix
   if (myrank == 0) then
      print *, 'getting the sphering matrix ...'; call flush(6)
   end if

   if (do_sphere) then
      if (seg_rank == 0) then
         lwork = 10*nx*nx
         !print *, 'doing eig nx = ', nx, ' lwork = ', lwork; call flush(6)
         call DCOPY(nx*nx,S,1,Stmp,1)

         call DSYEV('V','L',nx,Stmp,nx,eigs,work,lwork,info)

         if (info == 0) then
            print *, 'minimum eigenvalues = ', eigs(1:min(nx/2,3))
            write(20,*) 'minimum eigenvalues = ', eigs(1:min(nx/2,3))
            print *, 'maximum eigenvalues = ', eigs(nx:(nx-min(nx/2,3)+1):-1)
            write(20,*) 'maximum eigenvalues = ', eigs(nx:(nx-min(nx/2,3)+1):-1)
            call flush(6)
         else
            print "(a)", 'Error doing eigenvalue decomposition!!!'
            write(20,"(a)") 'Error doing eigenvalue decomposition!!!'
            call flush(6)
         end if

         numeigs = min(pcakeep,count(eigs > mineig))
         print *, 'num eigs kept = ', numeigs; call flush(6)
         write(20,*) 'num eigs kept = ', numeigs; call flush(6)

         allocate(Stmp2(nx,nx)); Stmp2 = dble(0.0)
         !call DSCAL(nx*nx,dble(0.0),Stmp2,1)
         !---reverse the order of eigenvectors
         do i = 1,nx
            eigv(i) = eigs(nx-i+1)
            do j = 1,nx
               Stmp2(i,j) = Stmp(j,nx-i+1)
            end do
            !call DCOPY(nx,Stmp(:,nx-i+1),1,Stmp2(:,i),1)
         end do
         call DCOPY(nx*nx,Stmp2,1,Stmp,1)
         sldet = dble(0.0)
         do i = 1,numeigs
            do j = 1,nx
               Stmp2(i,j) = Stmp2(i,j) / sqrt(eigv(i))
               !if (isNaN(Stmp2(i,j))) then
               !   print *, 'NaN! i,j = ', i, ',', j, ' eigv = ', eigv(i); call flush(6)
               !end if
            end do
            sldet = sldet - dble(0.5)*log(eigv(i))
            !call DSCAL(nx,dble(1.0)/sqrt(eigs(nx-i+1)),Stmp(:,i),1)
         end do

         if (numeigs == nx) then
            call DSCAL(nx*nx,dble(0.0),S,1)
            if (do_approx_sphere) then
               call DGEMM('T','N',nx,nx,nx,dble(1.0),Stmp,nx,Stmp2,nx,dble(1.0),S,nx)
            else
               call DCOPY(nx*nx,Stmp2,1,S,1)       
            endif
         else
            if (do_approx_sphere) then
               allocate(Stmp3(numeigs,numeigs)); Stmp3 = dble(0.0);
               call DCOPY(nx*nx,Stmp,1,S,1)
               call DGESVD( 'O', 'A', numeigs, numeigs, S, nx,  eigs,  S,  nx, Stmp,  nx, work, lwork, info )
               call DSCAL(numeigs*numeigs,dble(0.0),Stmp3,1)
               call DGEMM('T','T',numeigs,numeigs,numeigs,dble(1.0),Stmp,nx,S,nx,dble(1.0),Stmp3,numeigs)
               call DSCAL(nx*nx,dble(0.0),S,1)
               call DGEMM('N','N',numeigs,nx,numeigs,dble(1.0),Stmp3,numeigs,Stmp2,nx,dble(1.0),S,nx)
            else
               call DCOPY(nx*nx,Stmp2,1,S,1)       
            end if
         end if
      end if

   else
      !--- just normalize by the channel variances (don't sphere)
      call DCOPY(nx*nx,S,1,Stmp,1)
      call DSCAL(nx*nx,dble(0.0),S,1)

      sldet = dble(0.0)
      do i = 1,nx
         if (Stmp(i,i) .gt. dble(0.0)) then
            S(i,i) = dble(1.0) / sqrt(Stmp(i,i))
            sldet = sldet + dble(0.5)*log(S(i,i))
         end if
      end do
      numeigs = nx
   end if
end if

if (myrank == 0) then
   print *, 'sphering the data ...'
   call flush(6)
end if

call MPI_BCAST(S,nx*nx,MPI_DOUBLE_PRECISION,0,seg_comm,ierr)
call MPI_BCAST(sldet,1,MPI_DOUBLE_PRECISION,0,seg_comm,ierr)
call MPI_BCAST(numeigs,1,MPI_INTEGER,0,seg_comm,ierr)

allocate(xtmp(nx,maxval(blk_size))); xtmp = dble(0.0);
do seg = 1,numsegs
   fnum = dataseg(seg)%filenum
   fieldsize = dataseg(seg)%lastdim

   num_blocks = fieldsize / blk_size(seg)
   lastblocksize = mod(fieldsize,blk_size(seg))
   do k = 1,num_blocks
      bstrt = (k-1)*blk_size(seg) + 1
      bstp = bstrt + blk_size(seg) - 1
      call DSCAL(nx*blk_size(seg),dble(0.0),xtmp(:,1:blk_size(seg)),1)
      call DGEMM('N','N',nx,blk_size(seg),nx,dble(1.0),S,nx,dataseg(seg)%data(:,bstrt:bstp),nx,dble(1.0),xtmp(:,1:blk_size(seg)),nx)
      call DCOPY(nx*blk_size(seg),xtmp(:,1:blk_size(seg)),1,dataseg(seg)%data(:,bstrt:bstp),1)
   end do
   if (lastblocksize .ne. 0) then
      bstrt = num_blocks * blk_size(seg) + 1
      bstp = bstrt + lastblocksize - 1
      call DSCAL(nx*lastblocksize,dble(0.0),xtmp(:,1:lastblocksize),1)
      call DGEMM('N','N',nx,lastblocksize,nx,dble(1.0),S,nx,dataseg(seg)%data(:,bstrt:bstp),nx,dble(1.0),xtmp(:,1:lastblocksize),nx)
      call DCOPY(nx*lastblocksize,xtmp(:,1:lastblocksize),1,dataseg(seg)%data(:,bstrt:bstp),1)
   end if
end do
deallocate(xtmp)
nw = numeigs


if (seg_rank == 0) then
   ! get the pseudoinverse of S
   call DCOPY(nx*nx,S,1,Stmp2,1)
   allocate(Spinv(nx,numeigs)); Spinv = dble(0.0)
   allocate(sUtmp(numeigs,numeigs))
   allocate(sVtmp(numeigs,nx))
   print *, 'numeigs = ', numeigs; call flush(6)
   call DGESVD( 'A', 'S', numeigs, nx, Stmp2, nx, eigs, sUtmp, numeigs, sVtmp, numeigs, work, lwork, info )
   do i = 1,numeigs
      sVtmp(i,:) = sVtmp(i,:) / eigs(i)
   end do
   call DGEMM('T','T',nx,numeigs,numeigs,dble(1.0),sVtmp,numeigs,sUtmp,numeigs,dble(0.0),Spinv,nx)

end if

if (seg_rank == 0 .and. print_debug) then
   print *, 'S = '; call flush(6)
   call matout(S(1:2,1:2),2,2)
   print *, 'Sphered data = '; call flush(6)
   call matout(dataseg(1)%data(1:2,1:2),2,2)
end if

call flush(6)
call MPI_BARRIER(seg_comm,ierr)   

!-------------------- ALLOCATE VARIABLES ---------------------

if (num_comps .eq. -1) then
   num_comps = nw * num_models
end if

print *, myrank + 1, ': Allocating variables ...'; call flush(6) 


allocate( Dtemp (num_models),stat=ierr); call tststat(ierr); Dtemp = dble(0.0)
allocate( Dsum (num_models),stat=ierr); call tststat(ierr); Dsum = dble(0.0)
if (myrank == 0) then
   allocate( LL (max(1,max_iter)),stat=ierr); call tststat(ierr); LL = dble(0.0)
end if
allocate( pdtype(nw,num_models),stat=ierr); call tststat(ierr)
allocate( m2(nw,num_models),stat=ierr); call tststat(ierr); m2 = dble(0.0)
allocate( m4(nw,num_models),stat=ierr); call tststat(ierr); m4 = dble(0.0)
allocate( m2sum(nw,num_models),stat=ierr); call tststat(ierr); m2sum = dble(0.0)
allocate( m4sum(nw,num_models),stat=ierr); call tststat(ierr); m4sum = dble(0.0)
allocate( nsub(num_models),stat=ierr); call tststat(ierr);
pdtype = pdftype
if (pdftype == 1) then
   do_choose_pdfs = .true.
   numchpdf = 0
end if

allocate( comp_list (nw,num_models),stat=ierr); call tststat(ierr)
allocate( comp_used (num_comps),stat=ierr); call tststat(ierr)

allocate( c (nw,num_models),stat=ierr); call tststat(ierr); c = dble(0.0)
allocate( dc_numer_tmp (nw,num_models),stat=ierr); call tststat(ierr); dc_numer_tmp = dble(0.0)
allocate( dc_numer (nw,num_models),stat=ierr); call tststat(ierr); dc_numer = dble(0.0)
allocate( dc_denom_tmp (nw,num_models),stat=ierr); call tststat(ierr); dc_denom_tmp = dble(0.0)
allocate( dc_denom (nw,num_models),stat=ierr); call tststat(ierr); dc_denom = dble(0.0)
allocate( cx (nx,num_models),stat=ierr); call tststat(ierr); cx = dble(0.0)

allocate( A (nw,num_comps),stat=ierr); call tststat(ierr); A = dble(0.0)
allocate( Ax (nx,num_comps),stat=ierr); call tststat(ierr); Ax = dble(0.0)
allocate( W (nw,nw,num_models),stat=ierr); call tststat(ierr); W = dble(0.0)
allocate( wc (nw,num_models),stat=ierr); call tststat(ierr); wc = dble(0.0)


if (do_newton) then
   allocate( dbaralpha_numer_tmp (num_mix,nw,num_models),stat=ierr); call tststat(ierr); dbaralpha_numer_tmp = dble(0.0)
   allocate( dbaralpha_numer (num_mix,nw,num_models),stat=ierr); call tststat(ierr); dbaralpha_numer = dble(0.0)
   allocate( dbaralpha_denom_tmp (num_mix,nw,num_models),stat=ierr); call tststat(ierr); dbaralpha_denom_tmp = dble(0.0)
   allocate( dbaralpha_denom (num_mix,nw,num_models),stat=ierr); call tststat(ierr); dbaralpha_denom = dble(0.0)
   allocate( dlambda_numer_tmp (num_mix,nw,num_models),stat=ierr); call tststat(ierr); dlambda_numer_tmp = dble(0.0)
   allocate( dlambda_numer (num_mix,nw,num_models),stat=ierr); call tststat(ierr); dlambda_numer = dble(0.0)
   allocate( dlambda_denom_tmp (num_mix,nw,num_models),stat=ierr); call tststat(ierr); dlambda_denom_tmp = dble(0.0)
   allocate( dlambda_denom (num_mix,nw,num_models),stat=ierr); call tststat(ierr); dlambda_denom = dble(0.0)
   allocate( dkappa_numer_tmp (num_mix,nw,num_models),stat=ierr); call tststat(ierr); dkappa_numer_tmp = dble(0.0)
   allocate( dkappa_numer (num_mix,nw,num_models),stat=ierr); call tststat(ierr); dkappa_numer = dble(0.0)
   allocate( dkappa_denom_tmp (num_mix,nw,num_models),stat=ierr); call tststat(ierr); dkappa_denom_tmp = dble(0.0)
   allocate( dkappa_denom (num_mix,nw,num_models),stat=ierr); call tststat(ierr); dkappa_denom = dble(0.0)
   allocate( dsigma2_numer_tmp (nw,num_models),stat=ierr); call tststat(ierr); dsigma2_numer_tmp = dble(0.0)
   allocate( dsigma2_numer (nw,num_models),stat=ierr); call tststat(ierr); dsigma2_numer = dble(0.0)
   allocate( dsigma2_denom_tmp (nw,num_models),stat=ierr); call tststat(ierr); dsigma2_denom_tmp = dble(0.0)
   allocate( dsigma2_denom (nw,num_models),stat=ierr); call tststat(ierr); dsigma2_denom  = dble(0.0)
end if

allocate( Wtmp(nw,nw),stat=ierr); call tststat(ierr); Wtmp = dble(0.0)
allocate( Wtmp2(nw,nw,num_thrds),stat=ierr); call tststat(ierr); Wtmp2 = dble(0.0)
allocate( dAk (nw,num_comps),stat=ierr); call tststat(ierr); dAk = dble(0.0)
allocate( dA (nw,nw,num_models),stat=ierr); call tststat(ierr); dA = dble(0.0)
allocate( dWtmp (nw,nw,num_models),stat=ierr); call tststat(ierr); dWtmp = dble(0.0)
allocate( wr(nw),stat=ierr); call tststat(ierr); wr = dble(0.0)
allocate(nd(max(1,max_iter),num_comps),stat=ierr); call tststat(ierr); nd = dble(0.0)

allocate( zeta (num_comps),stat=ierr); call tststat(ierr); zeta = dble(0.0)
allocate( baralpha (num_mix,nw,num_models),stat=ierr); call tststat(ierr); baralpha = dble(0.0)
allocate( kappa (nw,num_models),stat=ierr); call tststat(ierr); kappa = dble(0.0)
allocate( lambda (nw,num_models),stat=ierr); call tststat(ierr); lambda = dble(0.0)
allocate( sigma2 (nw,num_models),stat=ierr); call tststat(ierr); sigma2 = dble(0.0)

allocate( gm (num_models),stat=ierr); call tststat(ierr); gm = dble(0.0)
!if (update_gm) then
   allocate( dgm_numer_tmp (num_models),stat=ierr); call tststat(ierr); dgm_numer_tmp = dble(0.0)
   allocate( dgm_numer (num_models),stat=ierr); call tststat(ierr); dgm_numer = dble(0.0)
!end if
allocate( alpha (num_mix, num_comps),stat=ierr); call tststat(ierr); alpha = dble(0.0)
if (update_alpha) then
   allocate( dalpha_numer (num_mix, num_comps),stat=ierr); call tststat(ierr); dalpha_numer = dble(0.0)
   allocate( dalpha_numer_tmp (num_mix, num_comps),stat=ierr); call tststat(ierr); dalpha_numer_tmp = dble(0.0)

   allocate( dalpha_denom (num_mix, num_comps),stat=ierr); call tststat(ierr); dalpha_denom = dble(0.0)
   allocate( dalpha_denom_tmp (num_mix, num_comps),stat=ierr); call tststat(ierr); dalpha_denom_tmp = dble(0.0)
end if

allocate( mu (num_mix, num_comps),stat=ierr); call tststat(ierr); mu = dble(0.0)
allocate( mutmp (num_mix, num_comps),stat=ierr); call tststat(ierr); mutmp = dble(0.0)
if (update_mu) then
   allocate( dmu_numer (num_mix, num_comps),stat=ierr); call tststat(ierr); dmu_numer = dble(0.0)
   allocate( dmu_numer_tmp (num_mix, num_comps),stat=ierr); call tststat(ierr); dmu_numer_tmp = dble(0.0)

   allocate( dmu_denom (num_mix, num_comps),stat=ierr); call tststat(ierr); dmu_denom = dble(0.0)
   allocate( dmu_denom_tmp (num_mix, num_comps),stat=ierr); call tststat(ierr); dmu_denom_tmp = dble(0.0)
end if

allocate( sbeta (num_mix, num_comps),stat=ierr); call tststat(ierr); sbeta = dble(0.0);
allocate( sbetatmp (num_mix, num_comps),stat=ierr); call tststat(ierr); sbetatmp = dble(0.0);
if (update_beta) then
   allocate( dbeta_numer (num_mix, num_comps),stat=ierr); call tststat(ierr); dbeta_numer = dble(0.0);
   allocate( dbeta_denom (num_mix, num_comps),stat=ierr); call tststat(ierr); dbeta_denom = dble(0.0);
   allocate( dbeta_numer_tmp (num_mix, num_comps),stat=ierr); call tststat(ierr); dbeta_numer_tmp = dble(0.0);
   allocate( dbeta_denom_tmp (num_mix, num_comps),stat=ierr); call tststat(ierr); dbeta_denom_tmp = dble(0.0);
end if

allocate( rho (num_mix, num_comps),stat=ierr); call tststat(ierr); rho = dble(0.0)
if (dorho) then
   allocate( rhotmp (num_mix, num_comps),stat=ierr); call tststat(ierr); rhotmp = dble(0.0)
   allocate( drho_numer (num_mix, num_comps),stat=ierr); call tststat(ierr); drho_numer = dble(0.0)
   allocate( drho_numer_tmp (num_mix, num_comps),stat=ierr); call tststat(ierr); drho_numer_tmp = dble(0.0)
   allocate( drho_denom (num_mix, num_comps),stat=ierr); call tststat(ierr); drho_denom = dble(0.0)
   allocate( drho_denom_tmp (num_mix, num_comps),stat=ierr); call tststat(ierr); drho_denom_tmp = dble(0.0)
end if
allocate(ipivnw(nw)); ipivnw = dble(0.0)
allocate(ipivnx(nx)); ipivnx = dble(0.0)

allocate(tmpvec3(nw)); tmpvec3 = dble(0.0)
allocate(tmpvec4(nw)); tmpvec4 = dble(0.0)

allocate(xtmp(100,100)); xtmp = dble(0.0)
call random_number(xtmp)
!call DRANDUNIFORM(100*100,dble(0.0),dble(1.0),state,xtmp,info)
deallocate(xtmp)


!------------------- INITIALIZE VARIABLES ----------------------
print *, myrank+1, ': Initializing variables ...'; call flush(6);


if (seg_rank == 0) then
   if (load_gm) then
      print *, myrank+1, ': reading gm init from: ', trim(adjustl(indirparam))//'/'//'gm'; call flush(6)
      open(unit=15,file=trim(adjustl(indirparam))//'/'//'gm',access='direct',recl=2*nbyte*num_models)
      read(15,rec=1) gm
      close(15)
   else
      gm = dble(1.0) / dble(num_models)
   end if

   if (load_alpha) then
      print *, myrank+1, ': reading alpha init from: ', trim(adjustl(indirparam))//'/'//'alpha'; call flush(6)
      open(unit=15,file=trim(adjustl(indirparam))//'/'//'alpha',access='direct',recl=2*nbyte*num_mix*num_comps)
      read(15,rec=1) alpha
      close(15)
   else
      alpha(1:num_mix,:) = dble(1.0) / dble(num_mix)
   end if

   if (load_mu) then
      print *, myrank+1, ': reading mu init from: ', trim(adjustl(indirparam))//'/'//'mu'; call flush(6)
      open(unit=15,file=trim(adjustl(indirparam))//'/'//'mu',access='direct',recl=2*nbyte*num_mix*num_comps)
      read(15,rec=1) mu
      close(15)
   else
      do k = 1,num_comps
         do j = 1,num_mix
            mu(j,k) = dble(j) - dble(1.0) - dble(num_mix-1)/dble(2.0)
         end do
      end do
      if (.not. fix_init) then
         call random_number(mutmp)
         !call DRANDUNIFORM(num_mix*num_comps,dble(0.0),dble(1.0),state,mutmp,info)
         mu(1:num_mix,:) = mu(1:num_mix,:) + dble(0.05)*(dble(1.0) - dble(2.0)*mutmp)
      end if
   end if

   if (load_beta) then
      print *, myrank+1, ': reading sbeta init from: ', trim(adjustl(indirparam))//'/'//'sbeta'; call flush(6)
      open(unit=15,file=trim(adjustl(indirparam))//'/'//'sbeta',access='direct',recl=2*nbyte*num_mix*num_comps)
      read(15,rec=1) sbeta
      close(15)
   else
      if (fix_init) then
         sbeta(1:num_mix,:) = dble(1.0)
      else
         call random_number(sbetatmp)
         !call DRANDUNIFORM(num_mix*num_comps,dble(0.0),dble(1.0),state,sbetatmp,info)
         sbeta(1:num_mix,:) = dble(1.0) + dble(0.1)*(dble(0.5) - sbetatmp)
      end if
   end if

   if (load_rho) then
      print *, myrank+1, ': reading rho init from: ', trim(adjustl(indirparam))//'/'//'rho'; call flush(6)
      open(unit=15,file=trim(adjustl(indirparam))//'/'//'rho',access='direct',recl=2*nbyte*num_mix*num_comps)
      read(15,rec=1) rho
      close(15)
   else
      rho(1:num_mix,:) = rho0
   end if

   if (load_c) then
      print *, myrank+1, ': reading c init from: ', trim(adjustl(indirparam))//'/'//'c'; call flush(6)
      open(unit=15,file=trim(adjustl(indirparam))//'/'//'c',access='direct',recl=2*nbyte*nw*num_models)
      read(15,rec=1) c
      close(15)
   else
      c = dble(0.0)
   end if

   if (load_A) then
      print *, myrank+1, ': reading A init from: ', trim(adjustl(indirparam))//'/'//'A'; call flush(6)
      open(unit=15,file=trim(adjustl(indirparam))//'/'//'A',access='direct',recl=2*nbyte*nw*num_comps)
      read(15,rec=1) A
      close(15)
      do h = 1,num_models
         do i = 1,nw
            comp_list(i,h) = (h-1) * nw + i
         end do            
      end do
   else
      !print *, myrank+1, ': Ainit is identity'; call flush(6)
      do h = 1,num_models
         if (fix_init) then
            A(:,(h-1)*nw+1:h*nw) = dble(0.0)
            do i = 1,nw
               A(i,(h-1)*nw+i) = dble(1.0)
            end do
         else
            call random_number(Wtmp)
            !CALL DRANDUNIFORM(nw*nw,dble(0.0),dble(1.0),state,Wtmp,info)
            A(:,(h-1)*nw+1:h*nw) = dble(0.01)*(dble(0.5)-Wtmp)
            do i = 1,nw
               A(i,(h-1)*nw+i) = dble(1.0)
               !wnrm = DNRM2(nw,W(i,:,h),1)
               Anrmk = sqrt(sum(A(:,(h-1)*nw+i)*A(:,(h-1)*nw+i)))
               A(:,(h-1)*nw+i) = A(:,(h-1)*nw+i) / Anrmk
               comp_list(i,h) = (h-1) * nw + i
            end do
         end if
      end do
   end if

   call get_unmixing_matrices      !print *, myrank+1, ': Ainit is identity'; call flush(6)


   if (load_comp_list) then
      print *, myrank+1, ': reading comp_list init from: ', trim(adjustl(indirparam))//'comp_list'; call flush(6)
      open(unit=15,file=trim(adjustl(indirparam))//'/'//'c',access='direct',recl=2*nbyte*nw*num_models)
      read(15,rec=1) comp_list
      close(15)
      comp_used = .false.
      do h = 1,num_models
         do i = 1,nw
            comp_used(comp_list(i,h)) = .true.
         end do
      end do
   else
      do h = 1,num_models
         do i = 1,nw
            comp_list(i,h) = (h-1) * nw + i
         end do
      end do
      comp_used = .true.
   end if



   if (print_debug) then
      print *, 'data ='; call flush(6)
      call matout(dataseg(1)%data(1:2,1:6),2,6)
      print *, 'gm = ', gm; call flush(6)
      print *, 'alpha = '; call flush(6)
      call matout(alpha(1:2,1:min(6,num_comps)),2,min(6,num_comps))
      print *, 'mu = '; call flush(6)
      call matout(mu(1:2,1:min(6,num_comps)),2,min(6,num_comps))
      print *, 'sbeta = '; call flush(6)
      call matout(sbeta(1:2,1:min(6,num_comps)),2,min(6,num_comps))
      print *, 'rho = '; call flush(6)
      call matout(rho(1:2,1:min(6,num_comps)),2,min(6,num_comps))
      print *, 'c = '; call flush(6)
      call matout(c(1:2,1:num_models),2,num_models)
      print *, 'A = '; call flush(6)
      call matout(A(1:2,1:min(6,num_comps)),2,min(6,num_comps))
      do h = 1,num_models
         print *, 'W ', h, ' = '; call flush(6)
         call matout(W(1:2,1:2,h),2,2)
      end do
   end if

end if


call MPI_BCAST(W,nw*nw*num_models,MPI_DOUBLE_PRECISION,0,seg_comm,ierr)
call MPI_BCAST(wc,nw*num_models,MPI_DOUBLE_PRECISION,0,seg_comm,ierr)
call MPI_BCAST(comp_list,nw*num_models,MPI_INTEGER,0,seg_comm,ierr)
call MPI_BCAST(gm,num_models,MPI_DOUBLE_PRECISION,0,seg_comm,ierr)
call MPI_BCAST(alpha,num_mix*num_comps,MPI_DOUBLE_PRECISION,0,seg_comm,ierr)
call MPI_BCAST(mu,num_mix*num_comps,MPI_DOUBLE_PRECISION,0,seg_comm,ierr)      
call MPI_BCAST(sbeta,num_mix*num_comps,MPI_DOUBLE_PRECISION,0,seg_comm,ierr)      
call MPI_BCAST(rho,num_mix*num_comps,MPI_DOUBLE_PRECISION,0,seg_comm,ierr)      
call MPI_BCAST(num_mix,num_comps,MPI_INTEGER,0,seg_comm,ierr)      


if (load_rej) then
   if (seg_rank == 0) then
      print *, myrank+1, ': reading LLt from: ', trim(adjustl(indirparam))//'/'//outfile_LLt; call flush(6)
   end if
   k = 0
   do j = 1,seg_nodes
      if (j .eq. seg_rank+1) then
         ngood = 0
         do seg = 1,numsegs
            open(unit=19,file=trim(adjustl(indirparam))//'/'//outfile_LLt,access='direct',recl=2*nbyte)
            do i = 1,dataseg(seg)%lastdim
               do jj = 1,num_models
                  read(19,rec=(k+i-1)*(num_models+1)+jj) dataseg(seg)%modloglik(jj,i)
               end do
               read(19,rec=(k+i)*(num_models+1)) dataseg(seg)%loglik(i)
            end do
            close(19)

            jj = 1
            do i = 1,dataseg(seg)%lastdim
               if (sum(dataseg(seg)%modloglik(:,i)) == dble(0.0)) then
                  dataseg(seg)%gooddata(i) = .false.
               else
                  dataseg(seg)%gooddata(i) = .true.
                  dataseg(seg)%goodinds(jj) = i
                  jj = jj + 1
               end if
            end do

            dataseg(seg)%numgood = count(dataseg(seg)%gooddata)
            ngood = ngood + dataseg(seg)%numgood

            k = k + dataseg(seg)%lastdim
         end do
      end if
      call MPI_BCAST(k,1,MPI_INTEGER,j-1,seg_comm,ierr)
   end do
   call MPI_ALLREDUCE(ngood,numgoodsum,1,MPI_INTEGER,MPI_SUM,seg_comm,ierr)
end if

!-------------------- Determine optimal block size -------------------
if (do_opt_block) then
    call determine_block_size
else
    call allocate_blocks
endif
print *, myrank + 1, ': block size = ', block_size
do seg = 1,numsegs
   blk_size(seg) = min(dataseg(seg)%lastdim,block_size)
end do
call flush(6)
call MPI_BARRIER(seg_comm,ierr)   


v = dble(1.0)/dble(num_models)

leave = .false.

print *, myrank+1, ': entering the main loop ...'; call flush(6)


!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX main loop XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
iter = 1
numrej = 0


call system_clock(c1)
do
   if (iter > max_iter) then
      exit
   end if

   !----- get determinants
   if (seg_rank == 0) then ! .and. (update_A .or. (iter == 1))) then
      do h = 1,num_models
         call DCOPY(nw*nw,W(:,:,h),1,Wtmp,1)
         lwork = 5*nx*nx
 
         call DGEQRF(nw,nw,Wtmp,nw,wr,work,lwork,info)
         if (info < 0) then
            print *, 'error getting QR!!!'
         end if
         Dtemp(h) = dble(0.0)
         do i = 1,nw
            if (Wtmp(i,i) .ne. dble(0.0)) then
               Dtemp(h) = Dtemp(h) + log(abs(Wtmp(i,i)))
            else
               print *, 'model ', h, ' determinant zero!'; call flush(6)
               Dtemp(h) = Dtemp(h) + minlog
            end if
         end do
      end do
      Dsum = Dtemp
      if (print_debug) then
         print *, 'Dsum = ', Dsum; call flush(6)
      end if
   end if
   call MPI_BCAST(Dsum,num_models,MPI_DOUBLE_PRECISION,0,seg_comm,ierr)

   !----- get updates: gm, alpha, mu, sbeta, rho, W
   !print *, myrank+1, ' calling get_updates ....'; call flush(6)
   !print *, myrank, ': gm = ', gm(1:min(5,num_models)); call flush(6)
   !print *, myrank, ': alpha = ', alpha(1:min(5,2*num_mix)); call flush(6)
   !print *, myrank, ': mu = ', mu(1:min(5,2*num_mix)); call flush(6)
   !print *, myrank, ': sbeta = ', sbeta(1:min(5,2*num_mix)); call flush(6)
   !print *, myrank, ': rho = ', rho(1:min(5,2*num_mix)); call flush(6)
   !print *, myrank, ': c = ', c(1:min(5,nw)); call flush(6)
   !print *, myrank, ': W = ', W(1:min(5,nw)); call flush(6)
   !print *, myrank+1, ': getting likelihood ....'; call flush(6)
   call get_updates_and_likelihood
   call accum_updates_and_likelihood

   !----- display log likelihood of data
   if (seg_rank == 0) then

      call system_clock(c2,counts_per_sec)
      t0 = dble(c2-c1)/dble(counts_per_sec)

      if (mod(iter,outstep) == 0) then
         print "(a,i6,a,f13.10,a,f14.10,a,f13.10,a,2(e13.5),a,f6.2,a,f5.1,a)",&
              ' iter', iter,' lrate = ', lrate, ' LL = ', LL(iter), &
              ' nd = ', ndtmpsum,', D = ', maxval(Dsum), minval(Dsum), &
              '  (', t0, ' s, ', dble(max_iter-iter)*t0/dble(3600), ' h)'
         write(20,"(a,i6,a,f13.10,a,f14.10,a,f13.10,a,2(e13.5),a,f6.2,a,f5.1,a)") ' iter', iter,' lrate = ', lrate,  ' LL = ', LL(iter), &
              ' nd = ', ndtmpsum,', D = ', maxval(Dsum), minval(Dsum), &
              '  (', t0, ' s, ', dble(max_iter-iter)*t0/dble(3600), ' h)'
         call flush(6)
      end if

      call system_clock(c1)
   end if

   !----- check whether likelihood is increasing
   if (seg_rank == 0) then
      ! if we get a NaN early, try to reinitialize and startover a few times 
      if ((iter .le. restartiter) .and. isNaN(LL(iter))) then
         if (numrestarts > maxrestarts) then
            leave = .true.
         else
            do h = 1,num_models
               if (fix_init) then
                  A(:,(h-1)*nw+1:h*nw) = dble(0.0)
                  do i = 1,nw
                     A(i,(h-1)*nw+i) = dble(1.0)
                  end do
               else
                  call random_number(Wtmp)
                  !CALL DRANDUNIFORM(nw*nw,dble(0.0),dble(1.0),state,Wtmp,info)
                  A(:,(h-1)*nw+1:h*nw) = dble(0.01)*(dble(0.5)-Wtmp)
                  do i = 1,nw
                     A(i,(h-1)*nw+i) = dble(1.0)
                     !wnrm = DNRM2(nw,W(i,:,h),1)
                     Anrmk = sqrt(sum(A(:,(h-1)*nw+i)*A(:,(h-1)*nw+i)))
                     A(:,(h-1)*nw+i) = A(:,(h-1)*nw+i) / Anrmk
                     comp_list(i,h) = (h-1) * nw + i
                  end do
               end if
            end do
            call get_unmixing_matrices        
            startover = .true.
            print *, 'Reinitializaing and starting over ...'; call flush(6)
            numrestarts = numrestarts + 1
         end if
      end if
      if (iter > 1) then
         if (isNaN(LL(iter)) .and. (iter > restartiter)) then
            leave = .true.;
            print *, 'Got NaN! Exiting ...'; call flush(6)
         end if
         if (LL(iter) < LL(iter-1)) then
            print *, 'Likelihood decreasing!';  call flush(6)
            if ((lrate .le. minlrate) .or. (ndtmpsum .le. min_nd)) then
               leave = .true.
               print *, 'minimum change threshold met, exiting loop ...';  call flush(6)
            else
               lrate = lrate * lratefact
               rholrate = rholrate * rholratefact
               numdecs = numdecs + 1
               if (numdecs .ge. maxdecs) then
                  lrate0 = lrate0 * lratefact
                  if (iter > newt_start) then
                     rholrate0 = rholrate0 * rholratefact
                  end if
                  if (do_newton .and. (iter > newt_start)) then
                     print *, 'Reducing maximum Newton lrate'; call flush(6)
                     newtrate = newtrate * lratefact
                  end if
                  numdecs = 0
               end if
            end if
         end if
         if (use_min_dll) then
             if ((LL(iter) - LL(iter-1)) < min_dll) then
                 numincs = numincs + 1
                 if (numincs > maxincs) then
                    leave = .true.
                    print *, 'Exiting because likelihood increasing by less than ', min_dll, ' for more than ', maxincs, ' iterations ...'
                    write(20,"(a,e13.5,a,i6,a)") 'Exiting because likelihood increasing by less than ', min_dll, ' for more than ', maxincs, &
                         ' iterations ...'
                end if
             else
                numincs = 0
             end if 
         end if
         if (use_grad_norm) then
             if (ndtmpsum .le. min_nd) then
                 leave = .true.
                 print *, 'Exiting because norm of weight gradient less than ', min_nd, ' ...'
                 write(20,"(a,e13.5,a)") 'Exiting because norm of weight gradient less than ', min_nd, ' ...'
             end if
         end if 
      end if
      if (do_newton .and. (iter == newt_start)) then
         print *, 'Starting Newton ... setting numdecs to 0'; call flush(6)
         numdecs = 0
      end if
  
      call flush(20)
   end if

   call MPI_BCAST(leave,1,MPI_LOGICAL,0,seg_comm,ierr)
   call MPI_BCAST(startover,1,MPI_LOGICAL,0,seg_comm,ierr)
   !call MPI_BCAST(redoiter,1,MPI_LOGICAL,0,seg_comm,ierr)

   if (leave) then
      exit
   end if

   if (startover) then
      call MPI_BCAST(W,nw*nw*num_models,MPI_DOUBLE_PRECISION,0,seg_comm,ierr)
      call MPI_BCAST(wc,nw*num_models,MPI_DOUBLE_PRECISION,0,seg_comm,ierr)
   else

      !----- do updates: gm, alpha, mu, sbeta, rho, W
      !print *, myrank+1, ': calling do_updates ....'; call flush(6)
      call update_params
      
      if ((writestep .ge. 0) .and. mod(iter,writestep) == 0) then
         !print *, myrank+1, ': calling write_output ...'; call flush(6)
         call write_output
      end if
      
      !----- reject data
      if ((do_reject .and. (maxrej > 0)) .and. ((iter==rejstart) .or. ((mod(max(1,iter-rejstart),rejint) == 0) .and. (numrej < maxrej)))) then
         !print *, myrank+1, ': calling reject data ....'; call flush(6)
         call reject_data
         numrej = numrej + 1
      end if
      
      iter = iter + 1
   end if
end do !iter

call write_output

call MPI_BCAST(tot_procs,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
call MPI_FINALIZE(ierr)

if (myrank == 0) then
   call system_clock(c2,counts_per_sec)
   t0 = dble(c2-c0)/dble(counts_per_sec)

   print "(a,f6.2,a)", '... done. Execution time: ', t0/dble(3600), ' h '; call flush(6)
   write(20,"(a,f6.2,a)") '... done. Execution time: ', t0/dble(3600), ' h '
   print *, 'output directory = ', trim(outdirparam); call flush(6)
   write(20,"(a,a)") 'output directory = ', trim(outdirparam)
   close(20)
end if


!----------------------------------------------------------------------
contains
!----------------------------------------------------------------------


subroutine get_updates_and_likelihood


  !if (update_gm) then
     dgm_numer_tmp = dble(0.0)
  !end if
  if (update_alpha) then
     dalpha_numer_tmp = dble(0.0)
     dalpha_denom_tmp = dble(0.0)
  end if
  if (update_mu) then
     dmu_numer_tmp = dble(0.0)
     dmu_denom_tmp = dble(0.0)
  end if
  if (update_beta) then
     dbeta_numer_tmp = dble(0.0)
     dbeta_denom_tmp = dble(0.0)
  end if
  if (dorho) then
     drho_numer_tmp = dble(0.0)
     drho_denom_tmp = dble(0.0)
  end if
  
  if (update_A .and. do_newton) then
     dbaralpha_numer_tmp = dble(0.0)
     dbaralpha_denom_tmp = dble(0.0)
     dkappa_numer_tmp = dble(0.0)
     dkappa_denom_tmp = dble(0.0)
     dlambda_numer_tmp = dble(0.0)
     dlambda_denom_tmp = dble(0.0)
     dsigma2_numer_tmp = dble(0.0)
     dsigma2_denom_tmp = dble(0.0)
  end if
  if (update_c) then
     dc_numer_tmp = dble(0.0)
     dc_denom_tmp = dble(0.0)
  end if

  dWtmp = dble(0.0)
  LLtmp = dble(0.0)

  !--------- loop over the segments ----------  
  !print *, 'Looping over the segments ...'; call flush(6)
  do seg = 1,numsegs
     ldim = dataseg(seg)%lastdim
     ngood = dataseg(seg)%numgood          
     if (do_reject) then
        num_blocks = ngood / (num_thrds*block_size)
     else
        num_blocks = ldim / (num_thrds*block_size)
     end if
     
     !--------- loop over the blocks ----------
     !print *, 'Looping over the blocks in segment ', seg; call flush(6)
     do blk = 1,num_blocks
        !print *, 'Doing block ', blk; call flush(6)        
        x0strt = (blk-1)*num_thrds*block_size + 1
        if (blk < num_blocks) then
           x0stp = blk*num_thrds*block_size
        else
           if (do_reject) then
              x0stp = ngood
           else
              x0stp = ldim
           end if
        end if

        !print *, myrank+1, ': Setting bsize ... '; call flush(6)
        
        bsize = x0stp - x0strt + 1
        if (bsize .le. 0) then
           exit
        end if
        
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP & PRIVATE (thrdnum,tblksize,t,h,i,j,k,xstrt,xstp,bstrt,bstp,LLinc,tmpsum,usum,vsum)
        
        thrdnum = omp_get_thread_num()
        tblksize = bsize / num_thrds

        !print *, myrank+1, thrdnum+1, ': Inside openmp code ... '; call flush(6)
        
        xstrt = x0strt + thrdnum*tblksize
        bstrt = thrdnum*tblksize + 1
        if (thrdnum+1 < num_thrds) then
           xstp = xstrt + tblksize - 1
           bstp = bstrt + tblksize - 1
        else
           xstp = x0stp
           bstp = bsize
        end if
        
        tblksize = bstp - bstrt + 1
        
        do h = 1,num_models

           !print *, myrank+1,':', thrdnum+1,': initializing Ptmp ... D = ', Dsum(h), ' lgm = ', log(gm(h)), ' sldet = ', sldet; call flush(6)
           Ptmp(bstrt:bstp,h) = Dsum(h) + log(gm(h)) + sldet
           !--- get b
           !print *, myrank+1,':', thrdnum+1,': getting b ...'; call flush(6)
           if (update_c .and. update_A) then
              do i = 1,nw
                 !b(bstrt:bstp,i,h) = spread(dble(-1.0)*wc(i,h),1,tblksize)
                 !do t = bstrt,bstp
                    b(bstrt:bstp,i,h) = dble(-1.0)*wc(i,h)
                 !end do
              end do
           else
              call DSCAL(nw*tblksize,dble(0.0),b(bstrt:bstp,:,h),1)
           end if
           !print *, myrank+1,':', thrdnum+1,': calling DGEMM ...'; call flush(6)
           if (do_reject) then
              call DGEMM('T','T',tblksize,nw,nw,dble(1.0),dataseg(seg)%data(:,dataseg(seg)%goodinds(xstrt:xstp)),nx, &
                   W(:,:,h),nw,dble(1.0),b(bstrt:bstp,:,h),tblksize)
           else
              call DGEMM('T','T',tblksize,nw,nw,dble(1.0),dataseg(seg)%data(:,xstrt:xstp),nx,W(:,:,h),nw,dble(1.0), &
                   b(bstrt:bstp,:,h),tblksize)
           end if
           if (print_debug .and. (blk == 1) .and. (thrdnum == 0)) then
              print *, myrank+1, ': b ', h, ' = '; call flush(6)
              call matout(b(1:6,1:2,h),6,2)
           end if
           
           !--- get y z
           !print *, myrank+1,':', thrdnum+1,': getting y ...'; call flush(6)
           do i = 1,nw
              !--- get probability
              select case (pdtype(i,h))
              case (0)
                 do j = 1,num_mix
                    y(bstrt:bstp,i,j,h) = sbeta(j,comp_list(i,h)) * ( b(bstrt:bstp,i,h) - mu(j,comp_list(i,h)) )
                    if (rho(j,comp_list(i,h)) == dble(1.0)) then
#ifdef MKL
                       !call vdExp(tblksize,-abs(y(bstrt:bstp,i,j,h)),tmpvec(bstrt:bstp))
                       !z0(bstrt:bstp,j) = alpha(j,comp_list(i,h)) * sbeta(j,comp_list(i,h)) * tmpvec(bstrt:bstp) / dble(2.0)
#else
                       !call vrda_exp(tblksize,-abs(y(j,bstrt:bstp)),tmpvec(bstrt:bstp))
                       !z0(j,bstrt:bstp) = alpha(j,i,h) * sbeta(j,i,h) * tmpvec(bstrt:bstp) / dble(2.0)
#endif
                       z0(bstrt:bstp,j) = log(alpha(j,comp_list(i,h))) + log(sbeta(j,comp_list(i,h)))  &
                            - abs(y(bstrt:bstp,i,j,h)) - log(dble(2.0))
                    else if (rho(j,comp_list(i,h)) == dble(2.0)) then
                       !call vrda_exp(tblksize,-y(j,bstrt:bstp)*y(j,bstrt:bstp),tmpvec(bstrt:bstp))
                       !z0(j,bstrt:bstp) = alpha(j,i,h) * sbeta(j,i,h) * tmpvec(bstrt:bstp) / dble(1.772453851)
                       z0(bstrt:bstp,j) = log(alpha(j,comp_list(i,h))) + log(sbeta(j,comp_list(i,h)))  &
                            - y(bstrt:bstp,i,j,h)*y(bstrt:bstp,i,j,h) - log(dble(1.772453851))
                    else
#ifdef MKL
                       call vdLn(tblksize,abs(y(bstrt:bstp,i,j,h)),tmpvec(bstrt:bstp))
                       call vdExp(tblksize,rho(j,comp_list(i,h))*tmpvec(bstrt:bstp),tmpvec2(bstrt:bstp))                       
#else
                       call vrda_log(tblksize,abs(y(bstrt:bstp,i,j,h)),tmpvec(bstrt:bstp))
                       call vrda_exp(tblksize,rho(j,comp_list(i,h))*tmpvec(bstrt:bstp),tmpvec2(bstrt:bstp))
#endif
                       !call vrda_exp(tblksize,-tmpvec2(bstrt:bstp),tmpvec(bstrt:bstp))
                       z0(bstrt:bstp,j) = log(alpha(j,comp_list(i,h))) + log(sbeta(j,comp_list(i,h))) - tmpvec2(bstrt:bstp) &
                            - gamln(dble(1.0)+dble(1.0)/rho(j,comp_list(i,h))) - log(dble(2.0))
                    end if
                 end do
              case (2)
                 do j = 1,num_mix
                    y(bstrt:bstp,i,j,h) = sbeta(j,comp_list(i,h)) * ( b(bstrt:bstp,i,h) - mu(j,comp_list(i,h)) )
                    !call vrda_exp(tblksize,-dble(0.5)*y(j,bstrt:bstp)*y(j,bstrt:bstp),tmpvec(bstrt:bstp))
                    !z0(j,bstrt:bstp) = alpha(j,i,h) * sbeta(j,i,h) * tmpvec(bstrt:bstp) / dble(2.506628274)
                    z0(bstrt:bstp,j) = log(alpha(j,comp_list(i,h))) + log(sbeta(j,comp_list(i,h))) - &
                         dble(0.5)*y(bstrt:bstp,i,j,h)*y(bstrt:bstp,i,j,h) - log(dble(2.506628274))
                 end do
              case (3)
                 do j = 1,num_mix
                    y(bstrt:bstp,i,j,h) = sbeta(j,comp_list(i,h)) * ( b(bstrt:bstp,i,h) - mu(j,comp_list(i,h)) )
                    tmpvec(bstrt:bstp) = cosh(dble(0.5)*y(bstrt:bstp,i,j,h))
#ifdef MKL
                    call vdLn(tblksize,tmpvec(bstrt:bstp),tmpvec2(bstrt:bstp))
#else
                    call vrda_log(tblksize,tmpvec(bstrt:bstp),tmpvec2(bstrt:bstp))
#endif
                    !z0(j,bstrt:bstp) = alpha(j,i,h) * sbeta(j,i,h) * dble(0.25) / (tmpvec(bstrt:bstp)*tmpvec(bstrt:bstp))
                    z0(bstrt:bstp,j) = log(alpha(j,comp_list(i,h))) + log(sbeta(j,comp_list(i,h))) - &
                         dble(2.0)*tmpvec2(bstrt:bstp) - log(dble(4.0))
                 end do
              case (4)
                 y(bstrt:bstp,i,1,h) = sbeta(1,comp_list(i,h)) * ( b(bstrt:bstp,i,h) - mu(1,comp_list(i,h)) )              
                    !call vrda_exp(tblksize,-y(1,bstrt:bstp)*y(1,bstrt:bstp)/dble(2.0),tmpvec(bstrt:bstp))
                    !z0(1,bstrt:bstp) = alpha(1,i,h) * sbeta(1,i,h) * tmpvec(bstrt:bstp) * cosh(y(1,bstrt:bstp)) / dble(4.132731354)
                    tmpvec(bstrt:bstp) = cosh(y(bstrt:bstp,i,1,h))
#ifdef MKL
                    call vdLn(tblksize,tmpvec(bstrt:bstp),tmpvec2(bstrt:bstp))
#else
                    call vrda_log(tblksize,tmpvec(bstrt:bstp),tmpvec2(bstrt:bstp))
#endif
                    z0(bstrt:bstp,1) = log(sbeta(1,comp_list(i,h))) - &
                         dble(0.5)*y(bstrt:bstp,i,1,h)*y(bstrt:bstp,i,1,h) + tmpvec2(bstrt:bstp) - log(dble(4.132731354))
              case (1)
                 y(bstrt:bstp,i,1,h) = sbeta(1,comp_list(i,h)) * ( b(bstrt:bstp,i,h) - mu(1,comp_list(i,h)) )
                    !call vrda_exp(tblksize,-y(1,bstrt:bstp)*y(1,bstrt:bstp)/dble(2.0),tmpvec(bstrt:bstp))
                    !z0(1,bstrt:bstp) = alpha(1,i,h) * sbeta(1,i,h) * tmpvec(bstrt:bstp) / cosh(y(1,bstrt:bstp)) / dble(1.858073988)
                    tmpvec(bstrt:bstp) = cosh(y(bstrt:bstp,i,1,h))
#ifdef MKL
                    call vdLn(tblksize,tmpvec(bstrt:bstp),tmpvec2(bstrt:bstp))
#else
                    call vrda_log(tblksize,tmpvec(bstrt:bstp),tmpvec2(bstrt:bstp))
#endif
                    z0(bstrt:bstp,1) = log(sbeta(1,comp_list(i,h))) - &
                         dble(0.5)*y(bstrt:bstp,i,1,h)*y(bstrt:bstp,i,1,h) - tmpvec2(bstrt:bstp) - log(dble(1.858073988))
              end select
              
              !--- add the log likelihood of this component to the likelihood of this time point
              Pmax(bstrt:bstp) = maxval(z0(bstrt:bstp,:),2)

              ztmp(bstrt:bstp) = dble(0.0)
              do j = 1,num_mix
                 ztmp(bstrt:bstp) = ztmp(bstrt:bstp) + exp(z0(bstrt:bstp,j) - Pmax(bstrt:bstp))
              end do
              tmpvec(bstrt:bstp) = Pmax(bstrt:bstp) + log(ztmp(bstrt:bstp))
              Ptmp(bstrt:bstp,h) = Ptmp(bstrt:bstp,h) + tmpvec(bstrt:bstp)
              !--- get normalized z
              do j = 1,num_mix
                 z(bstrt:bstp,i,j,h) = dble(1.0) / exp(tmpvec(bstrt:bstp) - z0(bstrt:bstp,j))
              end do              
           end do
        end do

        !print *, myrank+1,':', thrdnum+1,': getting Pmax and v ...'; call flush(6)

        !--- get LL, v
        Pmax(bstrt:bstp) = maxval(Ptmp(bstrt:bstp,:),2)
        !print *, myrank+1,':', thrdnum+1,': getting Pmax blk = ', blk, '.sum(Pmax)=', sum(Pmax(bstrt:bstp)); call flush(6)

        vtmp(bstrt:bstp) = dble(0.0)
        do h = 1,num_models
           vtmp(bstrt:bstp) = vtmp(bstrt:bstp) + exp(Ptmp(bstrt:bstp,h) - Pmax(bstrt:bstp))
        end do
        P(bstrt:bstp) = Pmax(bstrt:bstp) + log(vtmp(bstrt:bstp))
        LLinc = sum( P(bstrt:bstp) )
!$OMP CRITICAL
        LLtmp = LLtmp + LLinc
!$OMP END CRITICAL
        do h = 1,num_models        
           if (do_reject) then
              dataseg(seg)%modloglik(h,dataseg(seg)%goodinds(xstrt:xstp)) = Ptmp(bstrt:bstp,h)
              dataseg(seg)%loglik(dataseg(seg)%goodinds(xstrt:xstp)) = P(bstrt:bstp)
           else
              dataseg(seg)%modloglik(h,xstrt:xstp) = Ptmp(bstrt:bstp,h)
              dataseg(seg)%loglik(xstrt:xstp) = P(bstrt:bstp)
           end if
           v(bstrt:bstp,h) = dble(1.0) / exp(P(bstrt:bstp) - Ptmp(bstrt:bstp,h))
        end do

        if (print_debug .and. (blk == 1) .and. (thrdnum == 0)) then
           print *, 'v = '; call flush(6)
           call matout(v(1:3,:),3,num_models)
        end if

        !--- get g, u, ufp
        !print *, myrank+1,':', thrdnum+1,': getting g ...'; call flush(6)
        do h = 1,num_models

           vsum = sum( v(bstrt:bstp,h) )

           !if (update_gm) then
!$OMP CRITICAL
              dgm_numer_tmp(h) = dgm_numer_tmp(h) + vsum 
!$OMP END CRITICAL
           !end if
        
           if (update_A) then
              call DSCAL(nw*tblksize,dble(0.0),g(bstrt:bstp,:),1)
           end if
           
           do i = 1,nw

              if (update_A .and. do_newton) then
                 !print *, myrank+1,':', thrdnum+1,': getting dsigma2 ...'; call flush(6)
                 tmpsum = sum( v(bstrt:bstp,h) * b(bstrt:bstp,i,h) * b(bstrt:bstp,i,h) )
!$OMP CRITICAL
                 dsigma2_numer_tmp(i,h) = dsigma2_numer_tmp(i,h) + tmpsum
                 dsigma2_denom_tmp(i,h) = dsigma2_denom_tmp(i,h) + vsum
!$OMP END CRITICAL
              end if
              !print *, myrank+1,':', thrdnum+1,': getting dc ...'; call flush(6)
              if (update_c) then
                 if (do_reject) then
                    tmpsum = sum( v(bstrt:bstp,h) * dataseg(seg)%data(i,dataseg(seg)%goodinds(xstrt:xstp)) )
                 else
                    tmpsum = sum( v(bstrt:bstp,h) * dataseg(seg)%data(i,xstrt:xstp) )
                 end if
!$OMP CRITICAL
                 dc_numer_tmp(i,h) = dc_numer_tmp(i,h) + tmpsum
                 dc_denom_tmp(i,h) = dc_denom_tmp(i,h) + vsum
!$OMP END CRITICAL
              end if

              !print *, myrank+1,':', thrdnum+1,': getting u ...'; call flush(6)
              do j = 1,num_mix

                 u(bstrt:bstp) = v(bstrt:bstp,h) * z(bstrt:bstp,i,j,h)
                 usum = sum( u(bstrt:bstp) )

                 !--- get fp, zfp
                 select case (pdtype(i,h))
                 case (0)
                    if (rho(j,comp_list(i,h)) == dble(1.0)) then
                       fp(bstrt:bstp) = sign(dble(1.0),y(bstrt:bstp,i,j,h))
                    else if (rho(j,comp_list(i,h)) == dble(2.0)) then
                       fp(bstrt:bstp) = y(bstrt:bstp,i,j,h) * dble(2.0)
                    else
#ifdef MKL
                       call vdLn(tblksize,abs(y(bstrt:bstp,i,j,h)),tmpvec(bstrt:bstp))
                       call vdExp(tblksize,(rho(j,comp_list(i,h))-dble(1.0))*tmpvec(bstrt:bstp),tmpvec2(bstrt:bstp))
#else
                       call vrda_log(tblksize,abs(y(bstrt:bstp,i,j,h)),tmpvec(bstrt:bstp))
                       call vrda_exp(tblksize,(rho(j,comp_list(i,h))-dble(1.0))*tmpvec(bstrt:bstp),tmpvec2(bstrt:bstp))
#endif
                       fp(bstrt:bstp) = rho(j,comp_list(i,h)) * sign(dble(1.0),y(bstrt:bstp,i,j,h)) * tmpvec2(bstrt:bstp)
                    end if
                 case (2)
                    fp(bstrt:bstp) = y(bstrt:bstp,i,j,h)
                 case (3)
                    fp(bstrt:bstp) = tanh(y(bstrt:bstp,i,j,h)/dble(2.0))
                 case (4)
                    fp(bstrt:bstp) = y(bstrt:bstp,i,1,h) - tanh(y(bstrt:bstp,i,1,h))
                 case (1)
                    fp(bstrt:bstp) = y(bstrt:bstp,i,1,h) + tanh(y(bstrt:bstp,i,1,h))
                 end select

                 ufp(bstrt:bstp) = u(bstrt:bstp) * fp(bstrt:bstp)
              
                 !--- get g
                 if (update_A) then
                    g(bstrt:bstp,i) = g(bstrt:bstp,i) + sbeta(j,comp_list(i,h)) * ufp(bstrt:bstp)
                    if (do_newton .and. (iter .ge. newt_start)) then
                       tmpsum = sum( ufp(bstrt:bstp) * fp(bstrt:bstp) ) * sbeta(j,comp_list(i,h))**2
!$OMP CRITICAL
                       dkappa_numer_tmp(j,i,h) = dkappa_numer_tmp(j,i,h) + tmpsum
                       dkappa_denom_tmp(j,i,h) = dkappa_denom_tmp(j,i,h) + usum
!$OMP END CRITICAL
                       tmpvec(bstrt:bstp) = fp(bstrt:bstp) * y(bstrt:bstp,i,j,h) - dble(1.0)
                       tmpsum = sum( u(bstrt:bstp) * tmpvec(bstrt:bstp) * tmpvec(bstrt:bstp) )
!$OMP CRITICAL
                       dlambda_numer_tmp(j,i,h) = dlambda_numer_tmp(j,i,h) + tmpsum
                       dlambda_denom_tmp(j,i,h) = dlambda_denom_tmp(j,i,h) + usum
                       dbaralpha_numer_tmp(j,i,h) = dbaralpha_numer_tmp(j,i,h) + usum
                       dbaralpha_denom_tmp(j,i,h) = dbaralpha_denom_tmp(j,i,h) + vsum
!$OMP END CRITICAL
                    end if
                 end if
                 !print *, myrank+1, ': g(',bstrt,':',bstrt+4,',',h,') = '; call flush(6)
                 !call matout(g(1:2,bstrt:bstrt+4,h),2,5)               
                 if (update_alpha) then
!$OMP CRITICAL
                    dalpha_numer_tmp(j,comp_list(i,h)) = dalpha_numer_tmp(j,comp_list(i,h)) + usum
                    dalpha_denom_tmp(j,comp_list(i,h)) = dalpha_denom_tmp(j,comp_list(i,h)) + vsum
!$OMP END CRITICAL
                 end if
                 if (update_mu) then
                    tmpsum = sum( ufp(bstrt:bstp) )
!$OMP CRITICAL
                    dmu_numer_tmp(j,comp_list(i,h)) = dmu_numer_tmp(j,comp_list(i,h)) + tmpsum
!$OMP END CRITICAL
                    if (rho(j,comp_list(i,h)) .le. dble(2.0)) then
                       tmpsum = sbeta(j,comp_list(i,h)) * sum( ufp(bstrt:bstp) / y(bstrt:bstp,i,j,h) )
                    else
                       tmpsum = sbeta(j,comp_list(i,h)) * sum( ufp(bstrt:bstp) * fp(bstrt:bstp) )
                    end if
!$OMP CRITICAL
                    dmu_denom_tmp(j,comp_list(i,h)) = dmu_denom_tmp(j,comp_list(i,h)) + tmpsum 
!$OMP END CRITICAL
                 end if
                 if (update_beta) then
!$OMP CRITICAL
                    dbeta_numer_tmp(j,comp_list(i,h)) = dbeta_numer_tmp(j,comp_list(i,h)) + usum
!$OMP END CRITICAL
                    if (rho(j,comp_list(i,h)) .le. dble(2.0)) then
                       tmpsum = sum( ufp(bstrt:bstp) * y(bstrt:bstp,i,j,h) )
!$OMP CRITICAL
                       dbeta_denom_tmp(j,comp_list(i,h)) =  dbeta_denom_tmp(j,comp_list(i,h)) + tmpsum
!$OMP END CRITICAL
                    end if
                 end if
                 if (dorho) then
                    tmpy(bstrt:bstp) = abs(y(bstrt:bstp,i,j,h))
#ifdef MKL
                    call vdLn(tblksize,tmpy(bstrt:bstp),logab(bstrt:bstp))
                    call vdExp(tblksize,rho(j,comp_list(i,h))*logab(bstrt:bstp),tmpy(bstrt:bstp))
                    call vdLn(tblksize,tmpy(bstrt:bstp),logab(bstrt:bstp))
#else
                    call vrda_log(tblksize,tmpy(bstrt:bstp),logab(bstrt:bstp))
                    call vrda_exp(tblksize,rho(j,comp_list(i,h))*logab(bstrt:bstp),tmpy(bstrt:bstp))
                    call vrda_log(tblksize,tmpy(bstrt:bstp),logab(bstrt:bstp))
#endif
                    where (tmpy(bstrt:bstp) < epsdble)
                       logab(bstrt:bstp) = dble(0.0)
                    end where
                    tmpsum = sum( u(bstrt:bstp) * tmpy(bstrt:bstp) * logab(bstrt:bstp) )
!$OMP CRITICAL
                    drho_numer_tmp(j,comp_list(i,h)) =  drho_numer_tmp(j,comp_list(i,h)) + tmpsum
                    drho_denom_tmp(j,comp_list(i,h)) =  drho_denom_tmp(j,comp_list(i,h)) + usum
!$OMP END CRITICAL
                    if (update_beta .and. (rho(j,comp_list(i,h)) > dble(2.0))) then
                       tmpsum = sum( u(bstrt:bstp) * tmpy(bstrt:bstp) )
!$OMP CRITICAL
                       dbeta_denom_tmp(j,comp_list(i,h)) = dbeta_denom_tmp(j,comp_list(i,h)) + tmpsum
!$OMP END CRITICAL
                    end if
                 else
                    if (update_beta .and. (rho(j,comp_list(i,h)) > dble(2.0))) then
                       tmpy(bstrt:bstp) = abs(y(bstrt:bstp,i,j,h))
#ifdef MKL
                       call vdLn(tblksize,tmpy(bstrt:bstp),logab(bstrt:bstp))
                       call vdExp(tblksize,rho(j,comp_list(i,h))*logab(bstrt:bstp),tmpy(bstrt:bstp))
#else
                       call vrda_log(tblksize,tmpy(bstrt:bstp),logab(bstrt:bstp))
                       call vrda_exp(tblksize,rho(j,comp_list(i,h))*logab(bstrt:bstp),tmpy(bstrt:bstp))
#endif
                       tmpsum = sum( u(bstrt:bstp) * tmpy(bstrt:bstp) )
!$OMP CRITICAL
                       dbeta_denom_tmp(j,comp_list(i,h)) = dbeta_denom_tmp(j,comp_list(i,h)) + tmpsum
!$OMP END CRITICAL
                    end if
                 end if

              end do
           end do

           if (print_debug .and. (blk == 1) .and. (thrdnum == 0)) then
              print *, 'g ', h, ' = '; call flush(6)
              call matout(g(1:6,1:2),6,2)
           end if
           
           if (update_A) then              
              call DSCAL(nw*nw,dble(0.0),Wtmp2(:,:,thrdnum+1),1)              
              call DGEMM('T','N',nw,nw,tblksize,dble(1.0),g(bstrt:bstp,:),tblksize,b(bstrt:bstp,:,h),tblksize, &
                   dble(1.0),Wtmp2(:,:,thrdnum+1),nw)              
!$OMP CRITICAL
              call DAXPY(nw*nw,dble(1.0),Wtmp2(:,:,thrdnum+1),1,dWtmp(:,:,h),1)
!$OMP END CRITICAL
           end if

        end do

!$OMP END PARALLEL

     end do

     !print *, myrank+1,': finished segment ', seg; call flush(6)

  end do

end subroutine get_updates_and_likelihood

!----------------------------------------------------------------------
!----------------------------------------------------------------------

subroutine accum_updates_and_likelihood

  !--- add to the cumulative dtmps
  !print *, myrank+1,': calling reduce ...', seg; call flush(6)
  !if (update_gm) then
     call MPI_REDUCE(dgm_numer_tmp,dgm_numer,num_models,MPI_DOUBLE_PRECISION,MPI_SUM,0,seg_comm,ierr)
  !end if
  if (update_alpha) then
     call MPI_REDUCE(dalpha_numer_tmp,dalpha_numer,num_mix*num_comps,MPI_DOUBLE_PRECISION,MPI_SUM,0,seg_comm,ierr)
     call MPI_REDUCE(dalpha_denom_tmp,dalpha_denom,num_mix*num_comps,MPI_DOUBLE_PRECISION,MPI_SUM,0,seg_comm,ierr)
  end if
  if (update_mu) then
     call MPI_REDUCE(dmu_numer_tmp,dmu_numer,num_mix*num_comps,MPI_DOUBLE_PRECISION,MPI_SUM,0,seg_comm,ierr)
     call MPI_REDUCE(dmu_denom_tmp,dmu_denom,num_mix*num_comps,MPI_DOUBLE_PRECISION,MPI_SUM,0,seg_comm,ierr)
  end if
  if (update_beta) then
     call MPI_REDUCE(dbeta_numer_tmp,dbeta_numer,num_mix*num_comps,MPI_DOUBLE_PRECISION,MPI_SUM,0,seg_comm,ierr)
     call MPI_REDUCE(dbeta_denom_tmp,dbeta_denom,num_mix*num_comps,MPI_DOUBLE_PRECISION,MPI_SUM,0,seg_comm,ierr)
  end if
  if (dorho) then
     call MPI_REDUCE(drho_numer_tmp,drho_numer,num_mix*num_comps,MPI_DOUBLE_PRECISION,MPI_SUM,0,seg_comm,ierr)
     call MPI_REDUCE(drho_denom_tmp,drho_denom,num_mix*num_comps,MPI_DOUBLE_PRECISION,MPI_SUM,0,seg_comm,ierr)
  end if
  if (update_c) then
     call MPI_REDUCE(dc_numer_tmp,dc_numer,nw*num_models,MPI_DOUBLE_PRECISION,MPI_SUM,0,seg_comm,ierr)
     call MPI_REDUCE(dc_denom_tmp,dc_denom,nw*num_models,MPI_DOUBLE_PRECISION,MPI_SUM,0,seg_comm,ierr)
  end if  
  if (update_A) then
     call MPI_REDUCE(dWtmp,dA,nw*nw*num_models,MPI_DOUBLE_PRECISION,MPI_SUM,0,seg_comm,ierr)
     
     if (do_newton .and. (iter .ge. newt_start)) then
        call MPI_REDUCE(dbaralpha_numer_tmp,dbaralpha_numer,num_mix*nw*num_models,MPI_DOUBLE_PRECISION,MPI_SUM,0,seg_comm,ierr)
        call MPI_REDUCE(dbaralpha_denom_tmp,dbaralpha_denom,num_mix*nw*num_models,MPI_DOUBLE_PRECISION,MPI_SUM,0,seg_comm,ierr)
        call MPI_REDUCE(dkappa_numer_tmp,dkappa_numer,num_mix*nw*num_models,MPI_DOUBLE_PRECISION,MPI_SUM,0,seg_comm,ierr)
        call MPI_REDUCE(dkappa_denom_tmp,dkappa_denom,num_mix*nw*num_models,MPI_DOUBLE_PRECISION,MPI_SUM,0,seg_comm,ierr)
        call MPI_REDUCE(dlambda_numer_tmp,dlambda_numer,num_mix*nw*num_models,MPI_DOUBLE_PRECISION,MPI_SUM,0,seg_comm,ierr)
        call MPI_REDUCE(dlambda_denom_tmp,dlambda_denom,num_mix*nw*num_models,MPI_DOUBLE_PRECISION,MPI_SUM,0,seg_comm,ierr)
        call MPI_REDUCE(dsigma2_numer_tmp,dsigma2_numer,nw*num_models,MPI_DOUBLE_PRECISION,MPI_SUM,0,seg_comm,ierr)
        call MPI_REDUCE(dsigma2_denom_tmp,dsigma2_denom,nw*num_models,MPI_DOUBLE_PRECISION,MPI_SUM,0,seg_comm,ierr)
     end if
  end if

  if (seg_rank == 0) then
     if (update_A) then

        if (do_newton .and. (iter .ge. newt_start)) then
           baralpha = dbaralpha_numer / dbaralpha_denom
           sigma2 = dsigma2_numer / dsigma2_denom
           kappa = dble(0.0)
           lambda = dble(0.0)
           do h = 1,num_models
              do i = 1,nw
                 do j = 1,num_mix
                    dkap = dkappa_numer(j,i,h) / dkappa_denom(j,i,h)
                    kappa(i,h) = kappa(i,h) + baralpha(j,i,h) * dkap
                    lambda(i,h) = lambda(i,h) + &
                         baralpha(j,i,h) * ( dlambda_numer(j,i,h)/dlambda_denom(j,i,h) + dkap * mu(j,comp_list(i,h))**2 )
                 end do
              end do
           end do

           if (print_debug) then
              print *, 'sigma2 = '; call flush(6)
              call matout(sigma2(1:2,:),2,num_models)
              print *, 'lambda = '; call flush(6)
              call matout(lambda(1:2,:),2,num_models)
              print *, 'kappa = '; call flush(6)
              call matout(kappa(1:2,:),2,num_models)
              print *, 'baralpha = '; call flush(6)
              call matout(baralpha(1,1:2,:),2,num_models)
           end if
        end if

        nd(iter,:) = dble(0.0)
        no_newt = .false.

        do h = 1,num_models
           if (print_debug) then
              print *, 'dA ', h, ' = '; call flush(6)
              call matout(dA(1:2,1:2,h),2,2)
           end if
           
           if (do_reject) then
              call DSCAL(nw*nw,dble(-1.0)/dgm_numer(h),dA(:,:,h),1)
           else
              call DSCAL(nw*nw,dble(-1.0)/dgm_numer(h),dA(:,:,h),1)
           end if
           
           do i = 1,nw
              dA(i,i,h) = dA(i,i,h) + dble(1.0)
           end do

           if (print_debug) then
              print *, 'dA ', h, ' = '; call flush(6)
              call matout(dA(1:2,1:2,h),2,2)
           end if

           posdef = .true.                      
           if (do_newton .and. (iter .ge. newt_start)) then
                            
              do i = 1,nw
                 do k = 1,nw
                    if (i == k) then
                       Wtmp(i,i) = dA(i,i,h) / lambda(i,h)
                    else
                       sk1 = sigma2(i,h) * kappa(k,h)
                       sk2 = sigma2(k,h) * kappa(i,h)
                       if (sk1*sk2 > dble(1.0)) then
                          Wtmp(i,k) = (sk1*dA(i,k,h) - dA(k,i,h)) / (sk1*sk2 - dble(1.0))
                       else
                          posdef = .false.
                          no_newt = .true.
                       end if
                    end if
                 end do
              end do
           end if

           if ((.not. do_newton) .or. (.not. posdef) .or. (iter < newt_start)) then
              Wtmp = dA(:,:,h)
           end if
           
           call DSCAL(nw*nw,dble(0.0),dA(:,:,h),1)
           call DGEMM('N','N',nw,nw,nw,dble(1.0),A(:,comp_list(:,h)),nw,Wtmp,nw,dble(1.0),dA(:,:,h),nw)         
           
        end do


        dAk = dble(0.0)
        zeta = dble(0.0)
        do h = 1,num_models
           do i = 1,nw
              dAk(:,comp_list(i,h)) = dAk(:,comp_list(i,h)) + gm(h)*dA(:,i,h)
              zeta(comp_list(i,h)) = zeta(comp_list(i,h)) + gm(h)
           end do
        end do
        do k = 1,num_comps
           dAk(:,k) = dAk(:,k) / zeta(k)
        end do
        nd(iter,:) = sum(dAk*dAk,1)
        ndtmpsum = sqrt(sum(nd(iter,:),mask=comp_used) / (nw*count(comp_used)))

     end if
  end if

  
  call MPI_REDUCE(LLtmp,LLtmp2,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,seg_comm,ierr)
  if (seg_rank == 0) then
     if (do_reject) then
        LL(iter) = LLtmp2 / dble(numgoodsum*nw)
     else
        LL(iter) = LLtmp2 / dble(all_blks*nw)
     end if
  end if


end subroutine accum_updates_and_likelihood

!----------------------------------------------------------------------
!----------------------------------------------------------------------

subroutine update_params

  if (seg_rank == 0) then
     
     if (update_gm) then
        if (do_reject) then
           gm = dgm_numer / dble(numgoodsum)
        else
           gm = dgm_numer / dble(all_blks)
        end if
     end if
    
     if (update_alpha) then
        alpha = dalpha_numer / dalpha_denom
     end if

     if (update_c) then
        c = dc_numer / dc_denom
     end if

     !print *, 'updating A ...'; call flush(6)
     if (update_A .and. ((iter < share_start) .or. (mod(iter,share_iter) > 5))) then
        if (do_newton .and. (.not. no_newt) .and. (iter .ge. newt_start)) then
           lrate = min( newtrate, lrate + min(dble(1.0)/dble(newt_ramp),lrate) )
           rholrate = rholrate0
           call DAXPY(nw*num_comps,dble(-1.0)*lrate,dAk,1,A,1)
        else
           if (.not. posdef) then
              print *, 'Hessian not positive definite, using natural gradient'; call flush(6)
           end if
           lrate = min( lrate0, lrate + min(dble(1.0)/dble(newt_ramp),lrate) )
           rholrate = rholrate0
           call DAXPY(nw*num_comps,dble(-1.0)*lrate,dAk,1,A,1)
        end if
     end if

     
     if (update_mu) then
        mu = mu + dmu_numer / dmu_denom
     end if
     
     if (update_beta) then
        sbeta = sbeta * sqrt( dbeta_numer / dbeta_denom )
        sbetatmp = min(invsigmax,sbeta)
        sbeta = max(invsigmin,sbetatmp)
     end if
     
     if (dorho) then
        do k = 1,num_comps
           do j = 1,num_mix
              rho(j,k) = rho(j,k) + rholrate * ( dble(1.0) - &
                   (rho(j,k) / psifun(dble(1.0)+dble(1.0)/rho(j,k))) * drho_numer(j,k) / drho_denom(j,k) )
           end do
        end do
        rhotmp = min(maxrho,rho)
        rho = max(minrho,rhotmp)
     end if


     !--- rescale
     !print *, 'rescaling A ...'; call flush(6)
     if (doscaling) then
        do k = 1,num_comps
           Anrmk = sqrt(sum(A(:,k)*A(:,k)))
           if (Anrmk > dble(0.0)) then
              A(:,k) = A(:,k) / Anrmk
              mu(:,k) = mu(:,k) * Anrmk
              sbeta(:,k) = sbeta(:,k) / Anrmk
           end if
           !print *, myrank+1, 'after rescale: A = '; call flush(6)
           !call matout(A(1:2,1:2),2,2)
        end do
     end if

     if (share_comps .and. (iter .ge. share_start) .and. (mod(iter-share_start,share_iter) == 0)) then
        free_pass = .false.
        call identify_shared_comps
     else
        free_pass = .false.
     end if

     call get_unmixing_matrices        

     if (print_debug) then
        print *, 'gm = ', gm; call flush(6)
        
        print *, 'alpha = '; call flush(6)
        call matout(alpha(1:2,1:min(6,num_comps)),2,min(6,num_comps))
        
        print *, 'mu = '; call flush(6)
        call matout(mu(1:2,1:min(6,num_comps)),2,min(6,num_comps))

        print *, 'sbeta = '; call flush(6)
        call matout(sbeta(1:2,1:min(6,num_comps)),2,min(6,num_comps))

        print *, 'rho = '; call flush(6)
        call matout(rho(1:2,1:min(6,num_comps)),2,min(6,num_comps))

        print *, 'c = '; call flush(6)
        call matout(c(1:2,:),2,num_models)

        print *, 'A = '; call flush(6)
        call matout(A(1:2,1:min(6,num_comps)),2,min(6,num_comps))

        do h = 1,num_models
           print *, 'W ', h, ' = '; call flush(6)
           call matout(W(1:2,1:2,h),2,2)
        end do
     end if


  end if

  call MPI_BCAST(gm,num_models,MPI_DOUBLE_PRECISION,0,seg_comm,ierr)
  call MPI_BCAST(alpha,num_mix*num_comps,MPI_DOUBLE_PRECISION,0,seg_comm,ierr)
  call MPI_BCAST(mu,num_mix*num_comps,MPI_DOUBLE_PRECISION,0,seg_comm,ierr)      
  call MPI_BCAST(sbeta,num_mix*num_comps,MPI_DOUBLE_PRECISION,0,seg_comm,ierr)      
  call MPI_BCAST(rho,num_mix*num_comps,MPI_DOUBLE_PRECISION,0,seg_comm,ierr)      
  call MPI_BCAST(W,nw*nw*num_models,MPI_DOUBLE_PRECISION,0,seg_comm,ierr)
  call MPI_BCAST(wc,nw*num_models,MPI_DOUBLE_PRECISION,0,seg_comm,ierr)

end subroutine update_params

!----------------------------------------------------------------------
!----------------------------------------------------------------------

subroutine get_pmi

    
end subroutine get_pmi

!----------------------------------------------------------------------
!----------------------------------------------------------------------

subroutine identify_shared_comps


do h = 1,num_models
   do hh = h+1,num_models
      do i = 1,nw
         do ii = 1,nw

            call DGEMV('N',nw,nw,dble(1.0),Spinv2,nw,A(:,comp_list(i,h)),1,dble(0.0),tmpvec3,1)
            call DGEMV('N',nw,nw,dble(1.0),Spinv2,nw,A(:,comp_list(ii,hh)),1,dble(0.0),tmpvec4,1)
            t0 = abs(sum(tmpvec3*A(:,comp_list(ii,hh)))) / (sqrt(sum(tmpvec3*A(:,comp_list(i,h))))*sqrt(sum(tmpvec4*A(:,comp_list(ii,hh)))))
            !print *, 't0 = ', t0; call flush(6)
            if (t0 .ge. comp_thresh) then
               do_update = .true.
               if (comp_list(ii,hh) .eq. comp_list(i,h)) then
                  !print *, 'Not identifying already identified components ', hh, '-', comp_list(ii,hh), ' and ', h, '-', comp_list(i,h);
                  call flush(6)
                  do_update = .false.
               else
                  do j = 1,num_models
                     if (any(comp_list(:,j) .eq. comp_list(i,h)) .and. any(comp_list(:,j) .eq. comp_list(ii,hh))) then
                        print *, 'Not identifying correlated components due to common presence in model ', j; call flush(6)
                        do_update = .false.
                     end if
                  end do
               end if
               if (do_update) then
                  print *, 'Identifying component ', hh, '-', ii, ' with component ', h, '-', i; call flush(6)
                  comp_used(comp_list(ii,hh)) = .false.
                  do j = 1,nw
                     do k = 1,num_models
                        if (comp_list(j,k) .eq. comp_list(ii,hh)) then
                           comp_list(j,k) = comp_list(i,h)
                        end if
                     end do
                  end do
                  free_pass = .true.
               end if
            end if
         end do
      end do
   end do
end do

if (share_comps) then
    print *, 'Number of unique components = ', count(comp_used); call flush(6)
end if
end subroutine identify_shared_comps

!----------------------------------------------------------------------
!----------------------------------------------------------------------

subroutine get_unmixing_matrices

  lwork = 5*nx*nx

  do h = 1,num_models

     call DCOPY(nw*nw,A(:,comp_list(:,h)),1,W(:,:,h),1)

     call DGETRF(nw,nw,W(:,:,h),nw,ipivnw,info)
     call DGETRI(nw,W(:,:,h),nw,ipivnw,work,lwork,info)
     if (info .gt. 0) then
        print *, 'Matrix A(:,:,', h, ') is singular!'; call flush(6)
        ! Get the pseudoinverse
        !call DGEMM('N','N',nw,nx,nw,dble(1.0),W(:,:,h),nw,S,nw,dble(1.0),Stmp,nw)
        !call DGESVD( 'O', 'A', nw, nw, Stmp2, nw,  eigs,  Stmp2,  nw, Stmp2,  nw, work, lwork, info )
        !do i = 1,nw
        !   if (eigs(i) > mineigv) then
        !      eigs(i) = dble(1.0) / eigs(i)
        !   else
        !      eigs(i) = dble(0.0)
        !   end if
        !   call DSCAL(nw,eigs(i),Stmp(:,i),1)
        !end do
        !call DGEMM('T','T',nx,nw,nw,dble(1.0),Stmp2,nw,Stmp,nw,dble(1.0),A(:,:,h),nx)
     end if
     if (info .lt. 0) then
        print *, 'Input ', -info, ' to DGETRI in get_unmixing_matrices subroutine had illegal value!'; call flush(6)
     end if

     call DGEMV('N',nw,nw,dble(1.0),W(:,:,h),nw,c(:,h),1,dble(0.0),wc(:,h),1)

  end do


end subroutine get_unmixing_matrices

!----------------------------------------------------------------------
!----------------------------------------------------------------------

subroutine allocate_blocks 

N1 = 2*block_size*num_thrds

allocate( b (N1,nw,num_models),stat=ierr); call tststat(ierr); b = dble(0.0)
allocate( g (N1,nw),stat=ierr); call tststat(ierr); g = dble(0.0)

allocate( v (N1,num_models),stat=ierr); call tststat(ierr); v = dble(0.0)

allocate( y (N1,nw,num_mix,num_models),stat=ierr); call tststat(ierr); y = dble(0.0)
allocate( z (N1,nw,num_mix,num_models),stat=ierr); call tststat(ierr); z = dble(0.0)
allocate( z0 (N1,num_mix),stat=ierr); call tststat(ierr); z0 = dble(0.0)
allocate( fp (N1),stat=ierr); call tststat(ierr); fp = dble(0.0)
allocate( ufp (N1),stat=ierr); call tststat(ierr); ufp = dble(0.0)
allocate( u (N1),stat=ierr); call tststat(ierr); u = dble(0.0)
allocate( utmp (N1),stat=ierr); call tststat(ierr); utmp = dble(0.0)
allocate( ztmp (N1),stat=ierr); call tststat(ierr); ztmp = dble(0.0)
allocate( vtmp (N1),stat=ierr); call tststat(ierr); vtmp = dble(0.0)

allocate( logab (N1),stat=ierr); call tststat(ierr); logab = dble(0.0)
allocate( tmpy(N1),stat=ierr); call tststat(ierr); tmpy = dble(0.0)

allocate( Ptmp (N1,num_models),stat=ierr); call tststat(ierr); Ptmp = dble(0.0)
!allocate( Ptmp2 (N1,num_models),stat=ierr); call tststat(ierr)

allocate( P (N1),stat=ierr); call tststat(ierr); P = dble(0.0)
allocate( Pmax (N1),stat=ierr); call tststat(ierr); Pmax = dble(0.0)

allocate( tmpvec(N1),stat=ierr); call tststat(ierr); tmpvec = dble(0.0)
allocate( tmpvec2(N1),stat=ierr); call tststat(ierr); tmpvec2 = dble(0.0)

end subroutine allocate_blocks 

!----------------------------------------------------------------------
!----------------------------------------------------------------------

subroutine deallocate_blocks 


deallocate( b  )
deallocate( g )

deallocate( v )

deallocate( y )
deallocate( z )
deallocate( z0 )
deallocate( fp )
deallocate( ufp )
deallocate( u )
deallocate( utmp )
deallocate( ztmp )
deallocate( vtmp )

deallocate( logab )
deallocate( tmpy )

deallocate( Ptmp )
!deallocate( Ptmp2 )

deallocate( P )
deallocate( Pmax ) 

deallocate( tmpvec )
deallocate( tmpvec2 )

end subroutine deallocate_blocks 

!----------------------------------------------------------------------
!----------------------------------------------------------------------

subroutine determine_block_size

  print *, myrank+1, ': Determining optimal block size ....'; call flush(6)

  allocate(blk_time(1+floor(dble(blk_max-blk_min)/dble(blk_step))))
  k = 1
  do block_size = blk_min,blk_max,blk_step
     call allocate_blocks
     call system_clock(c1)
     call get_updates_and_likelihood
     call system_clock(c2,counts_per_sec)
     blk_time(k) = dble(c2-c1)/dble(counts_per_sec)
     k = k + 1
     call deallocate_blocks
  end do

  print *, myrank+1, ': blk_times = ', blk_time; call flush(6)

  j = minloc(blk_time,1)
  block_size = blk_min + (j-1)*blk_step
  call allocate_blocks
  print *, myrank + 1, ': Minimum block time = ', blk_time(j); call flush(6)
end subroutine determine_block_size


subroutine determine_block_size2

  print *, myrank+1, ': Determining optimal block size ....'; call flush(6)

  allocate(xtmp(2*blk_max*num_thrds,nw)); xtmp = dble(0.0)
  allocate(tmpvec(2*blk_max*num_thrds)); tmpvec = dble(0.0)
  allocate(tmpvec2(2*blk_max*num_thrds)); tmpvec2 = dble(0.0)
  allocate(Wtmp(nw,nw)); Wtmp = dble(0.0)
  !print *, 'num times = ', 1+floor(dble(blk_max-blk_min)/dble(blk_step)); call flush(6)
  allocate(blk_time(1+floor(dble(blk_max-blk_min)/dble(blk_step))))

  call DCOPY(nw*nw,S(1:nw,1:nw),1,Wtmp,1)

  k = 1
  do block_size = blk_min,blk_max,blk_step


     !print *, myrank+1, ': Doing block_size = ', block_size, ' ...'; call flush(6)
     call system_clock(c1)

     do seg = 1,numsegs
        ldim = dataseg(seg)%lastdim
        num_blocks = ldim / (num_thrds*block_size)
     
        !--------- loop over the blocks ----------

        do blk = 1,num_blocks
           !print *, 'Doing block ', blk; call flush(6)        
           x0strt = (blk-1)*num_thrds*block_size + 1
           if (blk < num_blocks) then
              x0stp = blk*num_thrds*block_size
           else
              x0stp = ldim
           end if

           !print *, myrank+1, ': Setting bsize ... '; call flush(6)        
           bsize = x0stp - x0strt + 1
        
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP & PRIVATE (thrdnum,tblksize,t,i,xstrt,xstp,bstrt,bstp)
        
           thrdnum = omp_get_thread_num()
           tblksize = bsize / num_thrds
           
           !print *, myrank+1, thrdnum+1, ': Inside openmp code ... '; call flush(6)
           
           xstrt = x0strt + thrdnum*tblksize
           bstrt = thrdnum*tblksize + 1
           if (thrdnum+1 < num_thrds) then
              xstp = xstrt + tblksize - 1
              bstp = bstrt + tblksize - 1
           else
              xstp = x0stp
              bstp = bsize
           end if
           
           tblksize = bstp - bstrt + 1
        
           call DGEMM('T','T',tblksize,nw,nw,dble(1.0),dataseg(seg)%data(:,xstrt:xstp),nx,Wtmp,nw,dble(1.0),xtmp(bstrt:bstp,:),tblksize)

           do i = 1,nw
#ifdef MKL
              call vdLn(tblksize,abs(xtmp(bstrt:bstp,i)),tmpvec(bstrt:bstp))
              call vdExp(tblksize,dble(2.0)*tmpvec(bstrt:bstp),tmpvec2(bstrt:bstp))
#else
              call vrda_log(tblksize,abs(xtmp(bstrt:bstp,i)),tmpvec(bstrt:bstp))
              call vrda_exp(tblksize,dble(2.0)*tmpvec(bstrt:bstp),tmpvec2(bstrt:bstp))
#endif
              tmpvec(bstrt:bstp) = tmpvec2(bstrt:bstp)- gamln(dble(1.0)+dble(1.0)/dble(1.5)) - log(dble(2.0))
           end do
!$OMP END PARALLEL

        end do

     end do

     call system_clock(c2,counts_per_sec)
     blk_time(k) = dble(c2-c1)/dble(counts_per_sec)
     k = k + 1

  end do

  !print *, myrank+1, ': blk_time = ', blk_time; call flush(6)

  j = minloc(blk_time,1)
  block_size = blk_min + (j-1)*blk_step
  print *, myrank + 1, ': Minimum block time = ', blk_time(j); call flush(6)
  deallocate(xtmp)
  deallocate(Wtmp)
  deallocate(tmpvec)
  deallocate(tmpvec2)  

end subroutine determine_block_size2

!----------------------------------------------------------------------

subroutine reject_data

  if (myrank == 0) then
     print "(a)", 'Doing rejection ....'
     write(20,"(a)") 'Doing rejection ....'
     call flush(6)
  end if

  !--- get likelihood mean and standard deviation
  llmean = dble(0.0)
  llvar = dble(0.0)
  llmax = dble(-999999999.0)
  llmin = dble(999999999.0)
  do seg = 1,numsegs
     llmean = llmean + sum(dataseg(seg)%loglik,MASK=dataseg(seg)%gooddata)
     llvar = llvar + sum(dataseg(seg)%loglik*dataseg(seg)%loglik,MASK=dataseg(seg)%gooddata)
     llmax = max(llmax,maxval(dataseg(seg)%loglik,MASK=dataseg(seg)%gooddata))
     llmin = min(llmin,minval(dataseg(seg)%loglik,MASK=dataseg(seg)%gooddata))
  end do

  call MPI_ALLREDUCE(llmean,llmeansum,1,MPI_DOUBLE_PRECISION,MPI_SUM,seg_comm,ierr)
  call MPI_ALLREDUCE(llvar,llvarsum,1,MPI_DOUBLE_PRECISION,MPI_SUM,seg_comm,ierr)
  call MPI_ALLREDUCE(llmax,llmaxmax,1,MPI_DOUBLE_PRECISION,MPI_MAX,seg_comm,ierr)
  call MPI_ALLREDUCE(llmin,llminmin,1,MPI_DOUBLE_PRECISION,MPI_MIN,seg_comm,ierr)

  llmeansum = llmeansum / dble(numgoodsum)
  llvarsum = llvarsum / dble(numgoodsum)
  llsig = sqrt(llvarsum - llmeansum**2)

  ngood = 0
  do seg = 1,numsegs
     do j = 1,dataseg(seg)%numgood
        if (dataseg(seg)%loglik(dataseg(seg)%goodinds(j)) < llmeansum - llsig*rejsig) then
           dataseg(seg)%gooddata(dataseg(seg)%goodinds(j)) = .false.
           dataseg(seg)%loglik(dataseg(seg)%goodinds(j)) = dble(0.0)
           dataseg(seg)%modloglik(:,dataseg(seg)%goodinds(j)) = dble(0.0)
        end if
     end do
     !where (dataseg(seg)%loglik < llmeansum - rejsig*llsig)
     !   dataseg(seg)%loglik = dble(0.0)
     !   dataseg(seg)%gooddata = .false.
     !end where
     dataseg(seg)%numgood = count(dataseg(seg)%gooddata)
     ngood = ngood + dataseg(seg)%numgood
     k = 1
     do j = 1,dataseg(seg)%lastdim
        if (dataseg(seg)%gooddata(j)) then
           dataseg(seg)%goodinds(k) = j
           k = k + 1
        end if
     end do
  end do

  call MPI_ALLREDUCE(ngood,numgoodsum,1,MPI_INTEGER,MPI_SUM,seg_comm,ierr)

  if (myrank == 0) then
     print *, 'maximum likelihood value = ', llmaxmax/nx
     write(20,*) 'maximum likelihood value = ', llmaxmax/nx
     print *, 'minimum likelihood value = ', llminmin/nx
     write(20,*) 'minimum likelihood value = ', llminmin/nx
     print *, 'average likelihood value = ', llmeansum/nx
     write(20,*) 'average likelihood value = ', llmeansum/nx
     print *, 'standard deviation       = ', llsig/nx
     write(20,*) 'standard deviation       = ', llsig/nx
     print *, 'rejecting data with likelihood less than ', (llmeansum - rejsig*llsig)/nx
     write(20,*) 'rejecting data with likelihood less than ', (llmeansum - rejsig*llsig)/nx
     call flush(6)

     if (maxrej-numrej-1 > 1) then
        print *, 'rejected ', all_blks - numgoodsum, ' data points so far. Will perform rejection ', &
             maxrej-numrej-1, ' more times at intervals of ',rejint, ' iterations.'
        write(20,*) 'rejected ', all_blks - numgoodsum, ' data points so far. Will perform rejection ', &
             maxrej-numrej-1, ' more times at intervals of ',rejint, ' iterations.'
     elseif (maxrej-numrej-1 == 1) then
        print *, 'rejected ', all_blks - numgoodsum, ' data points so far. Will perform rejection one more time after ', &
             rejint, ' iterations.'
        write(20,*) 'rejected ', all_blks - numgoodsum, ' data points so far. Will perform rejection one more time after ', &
             rejint, ' iterations.'
     else
        print *, 'rejected ', all_blks - numgoodsum, ' data points. No further rejections will be performed.'
        write(20,*) 'rejected ', all_blks - numgoodsum, ' data points. No further rejections will be performed.'
     end if
  end if

end subroutine reject_data

!----------------------------------------------------------------------

subroutine write_output
  k = 0
  if (write_LLt) then
     do j = 1,seg_nodes
        if (j .eq. seg_rank+1) then
           do seg = 1,numsegs
              if (outdirparam == '') then
                 open(unit=19,file=outdir//date(5:8)//'_'//time(1:6)//'/'//outfile_LLt,access='direct',recl=2*nbyte)
              else
                 open(unit=19,file=trim(outdirparam)//'/'//outfile_LLt,access='direct',recl=2*nbyte)
              end if
              
              do i = 1,dataseg(seg)%lastdim
                 do jj = 1,num_models
                    write(19,rec=(k+i-1)*(num_models+1)+jj) dataseg(seg)%modloglik(jj,i)
                 end do
                 write(19,rec=(k+i)*(num_models+1)) dataseg(seg)%loglik(i)
              end do
              call flush(19)
              close(19)
              k = k + dataseg(seg)%lastdim
           end do
        end if
        call MPI_BCAST(k,1,MPI_INTEGER,j-1,seg_comm,ierr)
     end do
  end if



  if (myrank == 0) then
     !print *, 'Writing data ....'
     !print *, 'gm = ', gm(1:min(5,num_models)); call flush(6)
     !print *, 'alpha = ', alpha(1:min(5,2*num_mix)); call flush(6)
     !print *, 'mu = ', mu(1:min(5,2*num_mix)); call flush(6)
     !print *, 'sbeta = ', sbeta(1:min(5,2*num_mix)); call flush(6)
     !print *, 'rho = ', rho(1:min(5,2*num_mix)); call flush(6)
     !print *, 'c = ', c(1:min(5,nw)); call flush(6)
     !print *, 'W = ', W(1:min(5,nw)); call flush(6)

     if (do_history) then
        if (mod(iter,histstep) == 0) then
           if (iter == histstep) then
              call system('mkdir '//trim(outdirparam)//'/history') 
              !call system('mkdir '//trim(outdirparam)//'/history'//' &>/dev/null') 
           end if
           write(tmpstring,"(i15)") iter
           call system('mkdir '//trim(outdirparam)//'/history/'//trim(adjustl(tmpstring)))
           !call system('mkdir '//trim(outdirparam)//'/history/'//trim(adjustl(tmpstring))//' &>/dev/null') 
        end if
     end if

     call DGEMM('N','N',nx,num_models,nw,dble(1.0),Spinv,nx,c,nw,dble(0.0),cx,nx)
     call DGEMM('N','N',nx,num_comps,nw,dble(1.0),Spinv,nx,A,nw,dble(0.0),Ax,nx)

     !do i = 1,nx
     !   do j = 1,nx
     !      if (isNaN(S(i,j))) then
     !         print *, 'NaN! i,j = ', i, ',', j; call flush(6)
     !      end if
     !   end do
     !end do


     if (outdirparam == '') then
        open(unit=29,file=outdir//date(5:8)//'_'//time(1:6)//'/'//outfile_c,access='direct',status='replace',recl=2*nbyte*nw*num_models)
        write(29,rec=1) c; call flush(29); close(29)
        open(unit=9,file=outdir//date(5:8)//'_'//time(1:6)//'/'//outfile_w,access='direct',status='replace',recl=2*nbyte*nw*nw*num_models)
        write(9,rec=1) W; call flush(9); close(9)
        open(unit=10,file=outdir//date(5:8)//'_'//time(1:6)//'/'//outfile_gamma,access='direct',status='replace',recl=2*nbyte*num_models)
        write(10,rec=1) gm; call flush(10); close(10)
        open(unit=11,file=outdir//date(5:8)//'_'//time(1:6)//'/'//outfile_alpha,access='direct',status='replace',recl=2*nbyte*num_mix*num_comps)
        write(11,rec=1) alpha; call flush(11); close(11)
        open(unit=12,file=outdir//date(5:8)//'_'//time(1:6)//'/'//outfile_mu,access='direct',status='replace',recl=2*nbyte*num_mix*num_comps)
        write(12,rec=1) mu; call flush(12); close(12)
        open(unit=13,file=outdir//date(5:8)//'_'//time(1:6)//'/'//outfile_sbeta,access='direct',status='replace',recl=2*nbyte*num_mix*num_comps)
        write(13,rec=1) sbeta; call flush(13); close(13)
        open(unit=14,file=outdir//date(5:8)//'_'//time(1:6)//'/'//outfile_rho,access='direct',status='replace',recl=2*nbyte*num_mix*num_comps)
        write(14,rec=1) rho; call flush(14); close(14)
        open(unit=15,file=outdir//date(5:8)//'_'//time(1:6)//'/'//outfile_mean,access='direct',status='replace',recl=2*nbyte*nx)
        write(15,rec=1) mean; call flush(15); close(15)
        open(unit=16,file=outdir//date(5:8)//'_'//time(1:6)//'/'//outfile_sphere,access='direct',status='replace',recl=2*nbyte*nx*nx)
        write(16,rec=1) S; call flush(16); close(16)
        open(unit=17,file=outdir//date(5:8)//'_'//time(1:6)//'/'//outfile_LL,access='direct',status='replace',recl=2*nbyte*max(1,max_iter))
        write(17,rec=1) LL; call flush(17); close(17)
        if (write_nd) then
           open(unit=18,file=outdir//date(5:8)//'_'//time(1:6)//'/'//outfile_nd,access='direct',status='replace',recl=2*nbyte*max(1,max_iter)*num_comps)
           write(18,rec=1) nd; call flush(18); close(18)            
        end if
        open(unit=21,file=outdir//date(5:8)//'_'//time(1:6)//'/'//outfile_A,access='direct',status='replace',recl=2*nbyte*nw*num_comps)
        write(21,rec=1) A; call flush(21); close(21)
        open(unit=22,file=outdir//date(5:8)//'_'//time(1:6)//'/'//outfile_comp_list,access='direct',status='replace',recl=2*nbyte*nw*num_models)
        write(22,rec=1) comp_list; call flush(21); close(22)
     else
        open(unit=29,file=trim(outdirparam)//'/'//outfile_c,access='direct',status='replace',recl=2*nbyte*nw*num_models)
        write(29,rec=1) c; call flush(29); close(29)
        open(unit=9,file=trim(outdirparam)//'/'//outfile_w,access='direct',status='replace',recl=2*nbyte*nw*nw*num_models)
        write(9,rec=1) W; call flush(9); close(9)
        open(unit=10,file=trim(outdirparam)//'/'//outfile_gamma,access='direct',status='replace',recl=2*nbyte*num_models)
        write(10,rec=1) gm; call flush(10); close(10)
        open(unit=11,file=trim(outdirparam)//'/'//outfile_alpha,access='direct',status='replace',recl=2*nbyte*num_mix*num_comps)
        write(11,rec=1) alpha; call flush(11); close(11)
        open(unit=12,file=trim(outdirparam)//'/'//outfile_mu,access='direct',status='replace',recl=2*nbyte*num_mix*num_comps)
        write(12,rec=1) mu; call flush(12); close(12)
        open(unit=13,file=trim(outdirparam)//'/'//outfile_sbeta,access='direct',status='replace',recl=2*nbyte*num_mix*num_comps)
        write(13,rec=1) sbeta; call flush(13); close(13)
        open(unit=14,file=trim(outdirparam)//'/'//outfile_rho,access='direct',status='replace',recl=2*nbyte*num_mix*num_comps)
        write(14,rec=1) rho; call flush(14); close(14)
        open(unit=15,file=trim(outdirparam)//'/'//outfile_mean,access='direct',status='replace',recl=2*nbyte*nx)
        write(15,rec=1) mean; call flush(15); close(15)
        open(unit=16,file=trim(outdirparam)//'/'//outfile_sphere,access='direct',status='replace',recl=2*nbyte*nx*nx)
        write(16,rec=1) S; call flush(16); close(16)
        open(unit=17,file=trim(outdirparam)//'/'//outfile_LL,access='direct',status='replace',recl=2*nbyte*max(1,max_iter))
        write(17,rec=1) LL; call flush(17); close(17)            
        if (write_nd) then
           open(unit=18,file=trim(outdirparam)//'/'//outfile_nd,access='direct',status='replace',recl=2*nbyte*max(1,max_iter)*num_comps)
           write(18,rec=1) nd; call flush(18); close(18)    
        end if
        open(unit=21,file=trim(outdirparam)//'/'//outfile_A,access='direct',status='replace',recl=2*nbyte*nw*num_comps)
        write(21,rec=1) A; call flush(21); close(21)
        open(unit=22,file=trim(outdirparam)//'/'//outfile_comp_list,access='direct',status='replace',recl=2*nbyte*nw*num_models)
        write(22,rec=1) comp_list; call flush(21); close(22)
     end if
     if (do_history .and. mod(iter,histstep) == 0) then
        open(unit=29,file=trim(outdirparam)//'/history/'//trim(adjustl(tmpstring))//'/'//outfile_c,access='direct',status='replace',recl=2*nbyte*nw*num_models)
        write(29,rec=1) c; call flush(29); close(29)
        open(unit=9,file=trim(outdirparam)//'/history/'//trim(adjustl(tmpstring))//'/'//outfile_w,access='direct',status='replace',recl=2*nbyte*nw*nw*num_models)
        write(9,rec=1) W; call flush(9); close(9)
        open(unit=10,file=trim(outdirparam)//'/history/'//trim(adjustl(tmpstring))//'/'//outfile_gamma,access='direct',status='replace',recl=2*nbyte*num_models)
        write(10,rec=1) gm; call flush(10); close(10)
        open(unit=11,file=trim(outdirparam)//'/history/'//trim(adjustl(tmpstring))//'/'//outfile_alpha,access='direct',status='replace',recl=2*nbyte*num_mix*num_comps)
        write(11,rec=1) alpha; call flush(11); close(11)
        open(unit=12,file=trim(outdirparam)//'/history/'//trim(adjustl(tmpstring))//'/'//outfile_mu,access='direct',status='replace',recl=2*nbyte*num_mix*num_comps)
        write(12,rec=1) mu; call flush(12); close(12)
        open(unit=13,file=trim(outdirparam)//'/history/'//trim(adjustl(tmpstring))//'/'//outfile_sbeta,access='direct',status='replace',recl=2*nbyte*num_mix*num_comps)
        write(13,rec=1) sbeta; call flush(13); close(13)
        open(unit=14,file=trim(outdirparam)//'/history/'//trim(adjustl(tmpstring))//'/'//outfile_rho,access='direct',status='replace',recl=2*nbyte*num_mix*num_comps)
        write(14,rec=1) rho; call flush(14); close(14)
        open(unit=15,file=trim(outdirparam)//'/history/'//trim(adjustl(tmpstring))//'/'//outfile_mean,access='direct',status='replace',recl=2*nbyte*nx)
        write(15,rec=1) mean; call flush(15); close(15)
        open(unit=16,file=trim(outdirparam)//'/history/'//trim(adjustl(tmpstring))//'/'//outfile_sphere,access='direct',status='replace',recl=2*nbyte*nx*nx)
        write(16,rec=1) S; call flush(16); close(16)
        open(unit=21,file=trim(outdirparam)//'/history/'//trim(adjustl(tmpstring))//'/'//outfile_A,access='direct',status='replace',recl=2*nbyte*nw*num_comps)
        write(21,rec=1) A; call flush(21); close(21)
        open(unit=22,file=trim(outdirparam)//'/history/'//trim(adjustl(tmpstring))//'/'//outfile_comp_list,access='direct',status='replace',recl=2*nbyte*nw*num_models)
        write(22,rec=1) comp_list; call flush(21); close(22)
     end if
  end if
end subroutine write_output

!----------------------------------------------------------------------

subroutine get_seg_list
   integer :: i, j, k, gotfile, cnt0

   ! get number of blocks in each segment sample
   allocate(blocks_in_sample(num_files))
   do k = 1,num_files
      blocks_in_sample(k) = field_dim(k)
   end do
   print *, 'blocks in sample = ', blocks_in_sample
   if (minval(blocks_in_sample) < 1) then
      print *, 'Error: minimum segment dim must be greater than '; call flush(6)
      stop
   end if
   all_blks = sum(blocks_in_sample*num_samples)
   print *, 'total blocks = ', all_blks

   allocate(node_blocks(seg_nodes),stat=ierr); call tststat(ierr)
   do i = 1,seg_nodes
      node_blocks(i) = node_thrds(i) * all_blks / sum(node_thrds)
   end do
   print *, 'node blocks = ', node_blocks
   if (maxval(node_blocks) < 1) then
      print *, 'Error: max node_blocks less than 1'; call flush(6)
      stop
   end if

   seg_list(1,1,1) = 1   ! start at first file
   seg_list(1,1,2) = 1   ! first sample
   seg_list(1,1,3) = 1   ! first block

   do j = 1,seg_nodes-1
      bcnt = sum(node_blocks(1:j))
      cnt0 = 0
      gotfile = 0
      do k = 1,num_files
         do i = 1,num_samples(k)
            ! find the file for this endpoint
            if (cnt0 + blocks_in_sample(k) < bcnt) then
               cnt0 = cnt0 + blocks_in_sample(k)
               cycle
            else
               ! got the file and sample
               gotfile = 1

               if (bcnt - cnt0 < num_thrds) then
                  if (i > 1) then
                     seg_list(j,2,1) = k
                     seg_list(j,2,2) = i-1
                     seg_list(j,2,3) = field_dim(k)

                     seg_list(j+1,1,1) = k
                     seg_list(j+1,1,2) = i
                     seg_list(j+1,1,3) = 1
                  else
                     if (k > 1) then ! first sample of later file
                        seg_list(j,2,1) = k-1
                        seg_list(j,2,2) = num_samples(k-1)
                        seg_list(j,2,3) = field_dim(k-1)

                        seg_list(j+1,1,1) = k
                        seg_list(j+1,1,2) = 1
                        seg_list(j+1,1,3) = 1
                     else ! first sample of first file
                        ! avoided by error above
                        print *, 'error: in supposedly avoided seg_list condition!'; call flush(6)
                     end if
                  end if
               else
                  if (field_dim(k) - (bcnt - cnt0) < num_thrds) then
                     seg_list(j,2,1) = k
                     seg_list(j,2,2) = i
                     seg_list(j,2,3) = field_dim(k)
                     if (i < num_samples(k)) then
                        seg_list(j+1,1,1) = k
                        seg_list(j+1,1,2) = i + 1
                        seg_list(j+1,1,3) = 1
                     else
                        if (k < num_files) then
                           seg_list(j+1,1,1) = k + 1
                           seg_list(j+1,1,2) = 1
                           seg_list(j+1,1,3) = 1
                        else
                           print *, 'Error: in get_seg_list at the last file last sample too early'
                           stop
                        end if
                     end if
                  else
                     seg_list(j,2,1) = k
                     seg_list(j,2,2) = i
                     seg_list(j,2,3) = bcnt - cnt0
                     
                     seg_list(j+1,1,1) = k
                     seg_list(j+1,1,2) = i
                     seg_list(j+1,1,3) = bcnt - cnt0 + 1
                  end if
                  !print *, myrank+1, ': bcnt = ', bcnt, ' cnt0 = ', cnt0, ' seg_list = ', seg_list(j+1,1,3); call flush(6)
               end if
               exit
            end if
         end do
         if (gotfile == 1) then
            exit
         end if
      end do
   end do
   ! last endpoint is last block of last file
   seg_list(seg_nodes,2,1) = num_files
   seg_list(seg_nodes,2,2) = num_samples(num_files)
   seg_list(seg_nodes,2,3) = field_dim(num_files)
 
   do j = 1,seg_nodes
      print *, 'node ', j, ' start: file ', seg_list(j,1,1), ' sample ', seg_list(j,1,2), ' index ', seg_list(j,1,3)
      print *, 'node ', j, ' stop : file ', seg_list(j,2,1), ' sample ', seg_list(j,2,2), ' index ', seg_list(j,2,3)
   end do
end subroutine get_seg_list

!----------------------------------------------------------------------

subroutine get_data
  integer :: i, j, k

  filestart = seg_list(seg_rank+1,1,1)
  sampstart = seg_list(seg_rank+1,1,2)
  filestop = seg_list(seg_rank+1,2,1)
  sampstop = seg_list(seg_rank+1,2,2)
  
  if (filestart == filestop) then
     numsegs = sampstop - sampstart + 1
  else
     numsegs = num_samples(filestart) - sampstart + 1 + sampstop
     do k = filestart+1,filestop-1
        numsegs = numsegs + num_samples(k)
     end do
  end if
  
  !print *, myrank+1, ': filestart = ', filestart, ' filestop = ', filestop, ': sampstart = ', sampstart, ' sampstop = ', sampstop
  !CALL flush(6)
  
  allocate(dataseg(numsegs))
  
  ! get the first segment
  if (filestart == filestop .and. sampstop == sampstart) then

     lastd = field_dim(filestart)
     sampsize = data_dim * field_dim(filestart)
     datsize = data_dim

     ! load the data in this segment
     dataseg(1)%filenum = filestart     
     dataseg(1)%lastdim = seg_list(seg_rank+1,2,3) - seg_list(seg_rank+1,1,3) + 1
     if (.not. dble_data) then
        !print *, myrank+1, 'lastdim = ', dataseg(1)%lastdim, ': datsize = ', datsize; call flush(6)
        allocate(real_array(datsize),stat=ierr); call tststat(ierr)
        !print *, myrank+1, 'lastdim 2 = ', dataseg(1)%lastdim, ': datsize = ', datsize; call flush(6)
        allocate(dataseg(1)%data(datsize,dataseg(1)%lastdim),stat=ierr); call tststat(ierr)
        !print *, myrank+1, 'lastdim 3 = ', dataseg(1)%lastdim, ': datsize = ', datsize; call flush(6)
        open(unit=8,file=infile(filestart),access="direct",recl=nbyte*datsize)
        !print *, myrank+1, 'infile = ', trim(infile(filestart)), ' recl = ', nbyte*datsize; call flush(6)
        !print *, myrank+1, ' going in .....'; call flush(6)
        do i = 1,dataseg(1)%lastdim
            read(8,rec=(sampstart-1)*lastd + i + seg_list(seg_rank+1,1,3) - 1) real_array
            dataseg(1)%data(:,i) = dble(real_array)
        end do
        deallocate(real_array)
     else
        allocate(dble_array(datsize),stat=ierr); call tststat(ierr)
        allocate(dataseg(1)%data(datsize,dataseg(1)%lastdim),stat=ierr); call tststat(ierr)
        open(unit=8,file=infile(filestart),access="direct",recl=2*nbyte*datsize)
        do i = 1,dataseg(1)%lastdim
            read(8,rec=(sampstart-1)*lastd + i + seg_list(seg_rank+1,1,3) - 1) dble_array
            dataseg(1)%data(:,i) = dble_array
        end do
        deallocate(dble_array)
     end if 
     close(8)
     allocate(dataseg(1)%loglik(dataseg(1)%lastdim),stat=ierr); call tststat(ierr)
     dataseg(1)%loglik = dble(0.0)
     allocate(dataseg(1)%modloglik(num_models,dataseg(1)%lastdim),stat=ierr); call tststat(ierr)
     dataseg(1)%modloglik = dble(0.0)

     if (do_reject) then
        allocate(dataseg(1)%gooddata(dataseg(1)%lastdim),stat=ierr); call tststat(ierr)
        allocate(dataseg(1)%goodinds(dataseg(1)%lastdim),stat=ierr); call tststat(ierr)
        dataseg(1)%gooddata = .true.
        dataseg(1)%numgood = dataseg(1)%lastdim
        do j = 1,dataseg(1)%numgood
           dataseg(1)%goodinds(j) = j
        end do
     end if

  else if (filestart == filestop .and. sampstop > sampstart) then
     !print *, 'data 1 ...'; call flush(6)

     lastd = field_dim(filestart)
     sampsize = data_dim * field_dim(filestart)
     datsize = data_dim
     
     !print *, 'data 2 ...'; call flush(6)
     ! load the remaining data in the first sample
     dataseg(1)%filenum = filestart
     dataseg(1)%lastdim = lastd - seg_list(seg_rank+1,1,3) + 1

     if (.not. dble_data) then
        !print *, myrank+1, 'lastdim = ', dataseg(1)%lastdim, ': datsize = ', datsize; call flush(6)
        allocate(real_array(datsize),stat=ierr); call tststat(ierr)
        allocate(dataseg(1)%data(datsize,dataseg(1)%lastdim)); call tststat(ierr)
        open(unit=8,file=infile(filestart),access="direct",recl=nbyte*datsize)
        do i = 1,dataseg(1)%lastdim
            read(8,rec=(sampstart-1)*lastd + i + seg_list(seg_rank+1,1,3) - 1) real_array
            dataseg(1)%data(:,i) = dble(real_array)
        end do
        close(8)
        deallocate(real_array)
     else
        allocate(dble_array(datsize),stat=ierr); call tststat(ierr)
        allocate(dataseg(1)%data(datsize,dataseg(1)%lastdim)); call tststat(ierr)
        open(unit=8,file=infile(filestart),access="direct",recl=2*nbyte*datsize)
        do i = 1,dataseg(1)%lastdim
            read(8,rec=(sampstart-1)*lastd + i + seg_list(seg_rank+1,1,3) - 1) dble_array
            dataseg(1)%data(:,i) = dble_array
        end do
        close(8)
        deallocate(dble_array)
     end if 
     allocate(dataseg(1)%loglik(dataseg(1)%lastdim),stat=ierr); call tststat(ierr)
     dataseg(1)%loglik = dble(0.0)
     allocate(dataseg(1)%modloglik(num_models,dataseg(1)%lastdim),stat=ierr); call tststat(ierr)
     dataseg(1)%modloglik = dble(0.0)

     if (do_reject) then
        allocate(dataseg(1)%gooddata(dataseg(1)%lastdim),stat=ierr); call tststat(ierr)
        allocate(dataseg(1)%goodinds(dataseg(1)%lastdim),stat=ierr); call tststat(ierr)
        dataseg(1)%gooddata = .true.
        dataseg(1)%numgood = dataseg(1)%lastdim
        do j = 1,dataseg(1)%numgood
           dataseg(1)%goodinds(j) = j
        end do
     end if

     ! get the last sample segment
     dataseg(numsegs)%filenum = filestop
     dataseg(numsegs)%lastdim = seg_list(seg_rank+1,2,3)
     if (.not. dble_data) then
        allocate(real_array(datsize))
        allocate(dataseg(numsegs)%data(datsize,dataseg(numsegs)%lastdim))
        open(unit=8,file=infile(filestop),access="direct",recl=nbyte*datsize)
        do i = 1,dataseg(numsegs)%lastdim
            read(8,rec=(sampstop-1)*lastd + i) real_array
            dataseg(numsegs)%data(:,i) = dble(real_array)
        end do
        close(8)
        deallocate(real_array)
     else
        allocate(dble_array(datsize))
        allocate(dataseg(numsegs)%data(datsize,dataseg(numsegs)%lastdim))
        open(unit=8,file=infile(filestop),access="direct",recl=2*nbyte*datsize)
        do i = 1,dataseg(numsegs)%lastdim
            read(8,rec=(sampstop-1)*lastd + i) dble_array
            dataseg(numsegs)%data(:,i) = dble_array
        end do
        close(8)
        deallocate(dble_array)
     end if
     allocate(dataseg(numsegs)%loglik(dataseg(numsegs)%lastdim),stat=ierr); call tststat(ierr)
     dataseg(numsegs)%loglik = dble(0.0)
     allocate(dataseg(numsegs)%modloglik(num_models,dataseg(numsegs)%lastdim),stat=ierr); call tststat(ierr)
     dataseg(numsegs)%modloglik = dble(0.0)

     if (do_reject) then
        allocate(dataseg(numsegs)%gooddata(dataseg(numsegs)%lastdim),stat=ierr); call tststat(ierr)
        allocate(dataseg(numsegs)%goodinds(dataseg(numsegs)%lastdim),stat=ierr); call tststat(ierr)
        dataseg(numsegs)%gooddata = .true.
        dataseg(numsegs)%numgood = dataseg(numsegs)%lastdim
        do j = 1,dataseg(numsegs)%numgood
           dataseg(numsegs)%goodinds(j) = j
        end do
     end if

     ! get the middle sample segments
     if (.not. dble_data) then
        allocate(real_array(sampsize))
        open(unit=8,file=infile(filestart),access="direct",recl=nbyte*sampsize)
        do i = 1,sampstop-sampstart-1
            allocate(dataseg(i+1)%data(datsize,lastd))
            read(8,rec=sampstart+i) real_array
            dataseg(i+1)%data = reshape(dble(real_array),(/datsize,lastd/))
            dataseg(i+1)%filenum = filestart
            dataseg(i+1)%lastdim = lastd

            
            allocate(dataseg(i+1)%loglik(dataseg(i+1)%lastdim),stat=ierr); call tststat(ierr)
            dataseg(i+1)%loglik = dble(0.0)
            allocate(dataseg(i+1)%modloglik(num_models,dataseg(i+1)%lastdim),stat=ierr); call tststat(ierr)
            dataseg(i+1)%modloglik = dble(0.0)

            if (do_reject) then
                allocate(dataseg(i+1)%gooddata(dataseg(i+1)%lastdim),stat=ierr); call tststat(ierr)
                allocate(dataseg(i+1)%goodinds(dataseg(i+1)%lastdim),stat=ierr); call tststat(ierr)
                dataseg(i+1)%gooddata = .true.
                dataseg(i+1)%numgood = dataseg(i+1)%lastdim
                do j = 1,dataseg(i+1)%numgood
                    dataseg(i+1)%goodinds(j) = j
                end do
            end if 
         end do
         close(8)
         deallocate(real_array)
     else
        allocate(dble_array(sampsize))
        open(unit=8,file=infile(filestart),access="direct",recl=2*nbyte*sampsize)
        do i = 1,sampstop-sampstart-1
            allocate(dataseg(i+1)%data(datsize,lastd))
            read(8,rec=sampstart+i) dble_array
            dataseg(i+1)%data = reshape(dble_array,(/datsize,lastd/))
            dataseg(i+1)%filenum = filestart
            dataseg(i+1)%lastdim = lastd
            
            allocate(dataseg(i+1)%loglik(dataseg(i+1)%lastdim),stat=ierr); call tststat(ierr)
            dataseg(i+1)%loglik = dble(0.0)
            allocate(dataseg(i+1)%modloglik(num_models,dataseg(i+1)%lastdim),stat=ierr); call tststat(ierr)
            dataseg(i+1)%modloglik = dble(0.0)

            if (do_reject) then
                allocate(dataseg(i+1)%gooddata(dataseg(i+1)%lastdim),stat=ierr); call tststat(ierr)
                allocate(dataseg(i+1)%goodinds(dataseg(i+1)%lastdim),stat=ierr); call tststat(ierr)
                dataseg(i+1)%gooddata = .true.
                dataseg(i+1)%numgood = dataseg(i+1)%lastdim
                do j = 1,dataseg(i+1)%numgood
                    dataseg(i+1)%goodinds(j) = j
                end do
            end if 
         end do
         close(8)
         deallocate(dble_array)
     end if     
  else
     ! there is more than one file
     
     ! get the first segment
     lastd = field_dim(filestart)
     sampsize = data_dim * field_dim(filestart)
     datsize = data_dim
     
     dataseg(1)%filenum = filestart
     dataseg(1)%lastdim = lastd - seg_list(seg_rank+1,1,3) + 1

     if (.not. dble_data) then
        allocate(real_array(datsize))
        allocate(dataseg(1)%data(datsize,dataseg(1)%lastdim))
        open(unit=8,file=infile(filestart),access="direct",recl=nbyte*datsize)
        do i = 1,dataseg(1)%lastdim 
            read(8,rec=(sampstart-1)*lastd + i + seg_list(seg_rank+1,1,3) - 1) real_array
            dataseg(1)%data(:,i) = dble(real_array)
        end do
        close(8)
        deallocate(real_array)
     else
        allocate(dble_array(datsize))
        allocate(dataseg(1)%data(datsize,dataseg(1)%lastdim))
        open(unit=8,file=infile(filestart),access="direct",recl=2*nbyte*datsize)
        do i = 1,dataseg(1)%lastdim 
            read(8,rec=(sampstart-1)*lastd + i + seg_list(seg_rank+1,1,3) - 1) dble_array
            dataseg(1)%data(:,i) = dble_array
        end do
        close(8)
        deallocate(dble_array)
     end if
     allocate(dataseg(1)%loglik(dataseg(1)%lastdim),stat=ierr); call tststat(ierr)
     dataseg(1)%loglik = dble(0.0)
     allocate(dataseg(1)%modloglik(num_models,dataseg(1)%lastdim),stat=ierr); call tststat(ierr)
     dataseg(1)%modloglik = dble(0.0)

     if (do_reject) then
        allocate(dataseg(1)%gooddata(dataseg(1)%lastdim),stat=ierr); call tststat(ierr)
        allocate(dataseg(1)%goodinds(dataseg(1)%lastdim),stat=ierr); call tststat(ierr)
        dataseg(1)%gooddata = .true.
        dataseg(1)%numgood = dataseg(1)%lastdim
        do j = 1,dataseg(1)%numgood
           dataseg(1)%goodinds(j) = j
        end do
     end if

     ! get the remaining segments in the starting file
     if (.not. dble_data) then
        allocate(real_array(sampsize))
        open(unit=8,file=infile(filestart),access="direct",recl=nbyte*sampsize)
        snum = 2
        do i = 1,num_samples(filestart)-sampstart
            allocate(dataseg(i+1)%data(datsize,lastd))
            read(8,rec=sampstart+i) real_array
            dataseg(i+1)%data = reshape(dble(real_array),(/datsize,lastd/))
            dataseg(i+1)%filenum = filestart
            dataseg(i+1)%lastdim = lastd
            snum = snum + 1
            allocate(dataseg(i+1)%loglik(dataseg(i+1)%lastdim),stat=ierr); call tststat(ierr)
            dataseg(i+1)%loglik = dble(0.0)
            allocate(dataseg(i+1)%modloglik(num_models,dataseg(i+1)%lastdim),stat=ierr); call tststat(ierr)
            dataseg(i+1)%modloglik = dble(0.0)

            if (do_reject) then
                allocate(dataseg(i+1)%gooddata(dataseg(i+1)%lastdim),stat=ierr); call tststat(ierr)
                allocate(dataseg(i+1)%goodinds(dataseg(i+1)%lastdim),stat=ierr); call tststat(ierr)
                dataseg(i+1)%gooddata = .true.
                dataseg(i+1)%numgood = dataseg(i+1)%lastdim
                do j = 1,dataseg(i+1)%numgood
                    dataseg(i+1)%goodinds(j) = j
                end do
            end if
        end do
        close(8)
        deallocate(real_array)     
     else
        allocate(dble_array(sampsize))
        open(unit=8,file=infile(filestart),access="direct",recl=2*nbyte*sampsize)
        snum = 2
        do i = 1,num_samples(filestart)-sampstart
            allocate(dataseg(i+1)%data(datsize,lastd))
            read(8,rec=sampstart+i) dble_array
            dataseg(i+1)%data = reshape(dble_array,(/datsize,lastd/))
            dataseg(i+1)%filenum = filestart
            dataseg(i+1)%lastdim = lastd
            snum = snum + 1
            allocate(dataseg(i+1)%loglik(dataseg(i+1)%lastdim),stat=ierr); call tststat(ierr)
            dataseg(i+1)%loglik = dble(0.0)
            allocate(dataseg(i+1)%modloglik(num_models,dataseg(i+1)%lastdim),stat=ierr); call tststat(ierr)
            dataseg(i+1)%modloglik = dble(0.0)

            if (do_reject) then
                allocate(dataseg(i+1)%gooddata(dataseg(i+1)%lastdim),stat=ierr); call tststat(ierr)
                allocate(dataseg(i+1)%goodinds(dataseg(i+1)%lastdim),stat=ierr); call tststat(ierr)
                dataseg(i+1)%gooddata = .true.
                dataseg(i+1)%numgood = dataseg(i+1)%lastdim
                do j = 1,dataseg(i+1)%numgood
                    dataseg(i+1)%goodinds(j) = j
                end do
            end if
        end do
        close(8)
        deallocate(dble_array)     
     end if
     
     ! get the last segment
     lastd = field_dim(filestop)
     sampsize = data_dim * field_dim(filestop)
     datsize = data_dim

     dataseg(numsegs)%filenum = filestop
     dataseg(numsegs)%lastdim = seg_list(seg_rank+1,2,3)
     if (.not. dble_data) then
        allocate(real_array(datsize))
        allocate(dataseg(numsegs)%data(datsize,dataseg(numsegs)%lastdim))
        open(unit=8,file=infile(filestop),access="direct",recl=nbyte*datsize)
        do i = 1,dataseg(numsegs)%lastdim
            read(8,rec=(sampstop-1)*lastd + i) real_array
            dataseg(numsegs)%data(:,i) = dble(real_array)
        end do
        close(8)
        deallocate(real_array)
     else
        allocate(dble_array(datsize))
        allocate(dataseg(numsegs)%data(datsize,dataseg(numsegs)%lastdim))
        open(unit=8,file=infile(filestop),access="direct",recl=2*nbyte*datsize)
        do i = 1,dataseg(numsegs)%lastdim
            read(8,rec=(sampstop-1)*lastd + i) dble_array
            dataseg(numsegs)%data(:,i) = dble_array
        end do
        close(8)
        deallocate(dble_array)
     end if
     allocate(dataseg(numsegs)%loglik(dataseg(numsegs)%lastdim),stat=ierr); call tststat(ierr)
     dataseg(numsegs)%loglik = dble(0.0)
     allocate(dataseg(numsegs)%modloglik(num_models,dataseg(numsegs)%lastdim),stat=ierr); call tststat(ierr)
     dataseg(numsegs)%modloglik = dble(0.0)

     if (do_reject) then
        allocate(dataseg(numsegs)%gooddata(dataseg(numsegs)%lastdim),stat=ierr); call tststat(ierr)
        allocate(dataseg(numsegs)%goodinds(dataseg(numsegs)%lastdim),stat=ierr); call tststat(ierr)
        dataseg(numsegs)%gooddata = .true.
        dataseg(numsegs)%numgood = dataseg(numsegs)%lastdim
        do j = 1,dataseg(numsegs)%numgood
           dataseg(numsegs)%goodinds(j) = j
        end do
     end if
     
     ! get the remaining segments in stop file
     if (.not. dble_data) then
        allocate(real_array(sampsize))
        open(unit=8,file=infile(filestop),access="direct",recl=nbyte*sampsize)
        do i = 1,sampstop-1
            allocate(dataseg(numsegs-i)%data(datsize,lastd))
            read(8,rec=sampstop-i) real_array
            dataseg(numsegs-i)%data = reshape(dble(real_array),(/datsize,lastd/))
            dataseg(numsegs-i)%filenum = filestop
            dataseg(numsegs-i)%lastdim = lastd

            allocate(dataseg(numsegs-i)%loglik(dataseg(numsegs-i)%lastdim),stat=ierr); call tststat(ierr)
            dataseg(numsegs-i)%loglik = dble(0.0)
            allocate(dataseg(numsegs-i)%modloglik(num_models,dataseg(numsegs-i)%lastdim),stat=ierr); call tststat(ierr)
            dataseg(numsegs-i)%modloglik = dble(0.0)

            if (do_reject) then
                allocate(dataseg(numsegs-i)%gooddata(dataseg(numsegs-i)%lastdim),stat=ierr); call tststat(ierr)
                allocate(dataseg(numsegs-i)%goodinds(dataseg(numsegs-i)%lastdim),stat=ierr); call tststat(ierr)
                dataseg(numsegs-i)%gooddata = .true.
                dataseg(numsegs-i)%numgood = dataseg(numsegs-i)%lastdim
                do j = 1,dataseg(numsegs-i)%numgood
                    dataseg(numsegs-i)%goodinds(j) = j
                end do
            end if
            !snum = snum + 1
        end do
        close(8)
        deallocate(real_array)
     else
        allocate(dble_array(sampsize))
        open(unit=8,file=infile(filestop),access="direct",recl=2*nbyte*sampsize)
        do i = 1,sampstop-1
            allocate(dataseg(numsegs-i)%data(datsize,lastd))
            read(8,rec=sampstop-i) dble_array
            dataseg(numsegs-i)%data = reshape(dble_array,(/datsize,lastd/))
            dataseg(numsegs-i)%filenum = filestop
            dataseg(numsegs-i)%lastdim = lastd

            allocate(dataseg(numsegs-i)%loglik(dataseg(numsegs-i)%lastdim),stat=ierr); call tststat(ierr)
            dataseg(numsegs-i)%loglik = dble(0.0)
            allocate(dataseg(numsegs-i)%modloglik(num_models,dataseg(numsegs-i)%lastdim),stat=ierr); call tststat(ierr)
            dataseg(numsegs-i)%modloglik = dble(0.0)

            if (do_reject) then
                allocate(dataseg(numsegs-i)%gooddata(dataseg(numsegs-i)%lastdim),stat=ierr); call tststat(ierr)
                allocate(dataseg(numsegs-i)%goodinds(dataseg(numsegs-i)%lastdim),stat=ierr); call tststat(ierr)
                dataseg(numsegs-i)%gooddata = .true.
                dataseg(numsegs-i)%numgood = dataseg(numsegs-i)%lastdim
                do j = 1,dataseg(numsegs-i)%numgood
                    dataseg(numsegs-i)%goodinds(j) = j
                end do
            end if
            !snum = snum + 1
        end do
        close(8)
        deallocate(dble_array)
     end if
     
     ! get segments from  middle files
     do k = filestart+1,filestop-1
        lastd = field_dim(k)
        sampsize = data_dim * field_dim(k)
        datsize = data_dim
        print *, 'file ', k, 'snum = ', snum, ' lastd = ', lastd, ' numsamp = ', num_samples(k); call flush(6)
        if (.not. dble_data) then
            allocate(real_array(sampsize))
            open(unit=8,file=infile(k),access="direct",recl=nbyte*sampsize)
            do i = 1,num_samples(k)
                dataseg(snum)%filenum = k
                dataseg(snum)%lastdim = lastd
                allocate(dataseg(snum)%data(datsize,lastd))
                read(8,rec=i) real_array
                dataseg(snum)%data = reshape(dble(real_array),(/datsize,lastd/))

                allocate(dataseg(snum)%loglik(dataseg(snum)%lastdim),stat=ierr); call tststat(ierr)
                dataseg(snum)%loglik = dble(0.0)
                allocate(dataseg(snum)%modloglik(num_models,dataseg(snum)%lastdim),stat=ierr); call tststat(ierr)
                dataseg(snum)%modloglik = dble(0.0)

                if (do_reject) then
                    allocate(dataseg(snum)%gooddata(dataseg(snum)%lastdim),stat=ierr); call tststat(ierr)
                    allocate(dataseg(snum)%goodinds(dataseg(snum)%lastdim),stat=ierr); call tststat(ierr)
                    dataseg(snum)%gooddata = .true.
                    dataseg(snum)%numgood = dataseg(snum)%lastdim
                    do j = 1,dataseg(snum)%numgood
                        dataseg(snum)%goodinds(j) = j
                    end do
                end if
                snum = snum + 1
            end do
            close(8)
            deallocate(real_array)
        else
            allocate(dble_array(sampsize))
            open(unit=8,file=infile(k),access="direct",recl=2*nbyte*sampsize)
            do i = 1,num_samples(k)
                dataseg(snum)%filenum = k
                dataseg(snum)%lastdim = lastd
                allocate(dataseg(snum)%data(datsize,lastd))
                read(8,rec=i) dble_array
                dataseg(snum)%data = reshape(dble_array,(/datsize,lastd/))

                allocate(dataseg(snum)%loglik(dataseg(snum)%lastdim),stat=ierr); call tststat(ierr)
                dataseg(snum)%loglik = dble(0.0)
                allocate(dataseg(snum)%modloglik(num_models,dataseg(snum)%lastdim),stat=ierr); call tststat(ierr)
                dataseg(snum)%modloglik = dble(0.0)

                if (do_reject) then
                    allocate(dataseg(snum)%gooddata(dataseg(snum)%lastdim),stat=ierr); call tststat(ierr)
                    allocate(dataseg(snum)%goodinds(dataseg(snum)%lastdim),stat=ierr); call tststat(ierr)
                    dataseg(snum)%gooddata = .true.
                    dataseg(snum)%numgood = dataseg(snum)%lastdim
                    do j = 1,dataseg(snum)%numgood
                        dataseg(snum)%goodinds(j) = j
                    end do
                end if
                snum = snum + 1
            end do
            close(8)
            deallocate(dble_array)
        end if
     end do
     
  end if

end subroutine get_data

!----------------------------------------------------------------------

subroutine get_cmd_args
  integer :: i, j, k
  
  !numargs = iargc() - 4
  !numargs = iargc()
  !print *, 'numargs = ', numargs; call flush(6)
  !if (numargs < 1) then
  !   print "(a)", 'Usage: amica <paramfile>'
  !   stop
  !end if
  numargs = 1
  allocate(arglist(numargs))
  do i = 1,numargs
     call getarg(i,arglist(i))
  end do
  
  read(arglist(1),'(a)') infilename
  !print *, 'infilename = ', infilename(1:len(infilename))
  open(unit=28,file=infilename)
  !print *, 'numargs = ', numargs, '| first arg: ', trim(arglist(1)), ' | all args: ', arglist
  if (numargs >= 2) then
     print *, 'Error: must have either 1 arg (paramfile) or 3 (paramfile num_model_nodes num_segment_nodes)'
     stop
  end if
  
  do
     read(28,'(a)',iostat=ierr) tmpstring
     !print *, 'tmpstring 1 is <begin>', tmpstring, '<end>'
     if (ierr < 0) then
        !print *, 'leaving get_cmd_args ....'
        close(28)
        exit
     end if
     tmpstring = adjustl(tmpstring)
     !print *, 'tmpstring 2 is <begin>', tmpstring, '<end>'
     if (tmpstring(1:1) == '#' .or. tmpstring(1:1) == ' ') then
        cycle
     end if
     i = scan(tmpstring,' ')
     if (i .le. 1) then
        cycle
     end if
     argname = tmpstring(1:i-1)
     !print *, 'argname is <begin>', argname, '<end>'
     tmparg = adjustl(tmpstring(i:))
     !print *, 'tmparg is <begin>', tmparg, '<end>'

     select case(trim(argname))
     case('files')
        ! count the number of files
        call flush(6)
        num_files = 0
        num_files_tokens = 0
        do
           i = scan(tmparg,' ')
           if (i .le. 1) then
              exit
           end if
           !print *, 'tmprarg(i-1) = ', tmparg(i-1:i-1); call flush(6)
           if (tmparg(i-1:i-1) == '/') then
              call system('ls -F1 '//tmparg(1:i-1)//' > tmpflist')
              open(unit=29,file='tmpflist')
              !open(unit=29,file='tmpflist',status='scratch')
              do   
                 read(29,'(a)',iostat=ierr) tmpdirarg
                 if (ierr < 0) then
                    exit
                 end if
                 num_files = num_files + 1
              end do
              close(29)
           elseif (tmparg(i-1:i-1) == '\') then
              call system('dir /B '//tmparg(1:i-1)//' > tmpflist')
              open(unit=29,file='tmpflist')
              !open(unit=29,file='tmpflist',status='scratch')
              do   
                 read(29,'(a)',iostat=ierr) tmpdirarg
                 if (ierr < 0) then
                    exit
                 end if
                 num_files = num_files + 1
              end do
              close(29)
           else
              num_files = num_files + 1
           end if
           num_files_tokens = num_files_tokens + 1
           tmparg = adjustl(tmparg(i:))
        end do
        allocate(infile(num_files))
        allocate(num_samples(num_files))
        allocate(num_dir_files(num_files_tokens))
        ! read in the file names
        i = scan(tmpstring,' ')
        tmpstring = adjustl(tmpstring(i:))
        k = 1
        j = 1
        do
           if (j > num_files) then
              exit
           end if
           i = scan(tmpstring,' ')
           if (i .le. 1) then
              exit
           end if
           if (tmpstring(i-1:i-1) == '/') then
              call system('ls -F1 '//tmpstring(1:i-1)//' > tmpflist')
              open(unit=29,file='tmpflist')
              !open(unit=29,file='tmpflist',status='scratch')
              num_dir_files(k) = 0
              do
                 read(29,'(a)',iostat=ierr) tmpdirarg
                 if (ierr < 0) then
                    exit
                 end if
                 write(tmpstring2,'(a)') tmpstring(1:i-1)//tmpdirarg
                 read(tmpstring2,'(a)') infile(j)
                 !read(tmpstring(1:i-1)//tmpdirarg,'(a)',iostat=ierr) infile(j)
                 num_dir_files(k) = num_dir_files(k) + 1
                 j = j + 1
              end do
              close(29)
           else if (tmpstring(i-1:i-1) == '\') then
              call system('dir /B '//tmpstring(1:i-1)//' > tmpflist')
              open(unit=29,file='tmpflist')
              !open(unit=29,file='tmpflist',status='scratch')
              num_dir_files(k) = 0
              do
                 read(29,'(a)',iostat=ierr) tmpdirarg
                 if (ierr < 0) then
                    exit
                 end if
                 infile(j) = tmpstring(1:i-1)//tmpdirarg
          !       write(infile(j),'(a)') tmpstring(1:i-1)//tmpdirarg
                 !!!write(tmpstring2,'(a)') tmpstring(1:i-1)//tmpdirarg
                 !!!read(tmpstring2,'(a)') infile(j)
                 !read(tmpstring(1:i-1)//tmpdirarg,'(a)',iostat=ierr) infile(j)
                 num_dir_files(k) = num_dir_files(k) + 1
                 j = j + 1
              end do
              close(29)
           else
              read(tmpstring(1:i-1),'(a)') infile(j)
              num_dir_files(k) = 1
              j = j + 1
           end if
           k = k + 1
           tmpstring = adjustl(tmpstring(i:))
        end do
        print *, 'num_files = ', num_files; call flush(6)
        print *, 'FILES: '; call flush(6)
        do k = 1,num_files
           print *, trim(infile(k))
        end do
        print *, 'num_dir_files = ', num_dir_files; call flush(6)
     case('num_samples')
        if (num_files .le. 0) then
           print *, 'Error: files must precede num_samples in paramfile'
           stop
        end if
        j = 1
        do k=1,num_files_tokens
           i = scan(tmparg,' ')
           if (i .le. 1) then
              print *, 'Error: num_samples entries is less than num_files in paramfile'
              stop
           end if
           read(tmparg(1:i-1),'(i12)') tmp_read
           num_samples(j) = tmp_read
           do ii = 1,num_dir_files(k)-1
              num_samples(j+ii) = num_samples(j)                 
           end do
           j = j + num_dir_files(k)
           tmparg = adjustl(tmparg(i:))
        end do
        print *, 'num_samples = ', num_samples; call flush(6)
        call flush(6)
     case('data_dim')
        read(tmparg,'(i12)') data_dim
        print *, 'data_dim = ', data_dim; call flush(6)
     case('field_dim')
        !print *, 'num_files_tokens = ', num_files_tokens; call flush(6)
        if (num_files .le. 0) then
           print *, 'Error: files must precede field_dim in paramfile'
        end if
        allocate(field_dim(num_files))
        j = 1
        do k=1,num_files_tokens
           !print *, 'j = ', j, ' k = ', k; call flush(6)
           i = scan(tmparg,' ')
           !print *, 'i = ', i; call flush(6)
           if (i .le. 1) then
              print *, 'Error: in field_dim'
              stop
           end if
           read(tmparg(1:i-1),'(i12)') tmp_read
           !print *, 'tmp_read = ', tmp_read; call flush(6)
           field_dim(j) = tmp_read
           tmparg = adjustl(tmparg(i:))
           !print *, 'tmparg = ', tmparg; call flush(6)
       
           do kk = 1,num_dir_files(k)-1
              field_dim(j+kk) = field_dim(j)
           end do
           j = j + num_dir_files(k)
        end do
        print *, 'field_dim = ', field_dim; call flush(6)
     case('filter_length')
        read(tmparg,'(i12)') filter_length
        print *, 'filter length = ', filter_length; call flush(6)
     case('block_size')
        read(tmparg,'(i12)') block_size
        print *, 'initial matrix block_size = ', block_size; call flush(6)
     case('blk_min')
        read(tmparg,'(i12)') blk_min 
        print *, 'blk_min = ', blk_min; call flush(6)
     case('blk_max')
        read(tmparg,'(i12)') blk_max 
        print *, 'blk_max = ', blk_max; call flush(6)
     case('blk_step')
        read(tmparg,'(i12)') blk_step 
        print *, 'blk_step = ', blk_step; call flush(6)
     case('dft_length')
        read(tmparg,'(i12)') dft_length
        print *, 'dft_length = ', dft_length; call flush(6)
     case('num_comps')
        read(tmparg,'(i12)') num_comps
        print *, 'number of intial comps = ', num_comps; call flush(6)
     case('num_models')
        read(tmparg,'(i12)') num_models
        print *, 'number of models = ', num_models; call flush(6)
     case('num_mix_comps')
        read(tmparg,'(i12)') num_mix
        print *, 'number of density mixture components = ', num_mix; call flush(6)
     case('num_mix')
        read(tmparg,'(i12)') num_mix
        print *, 'max number of density mixture components = ', num_mix; call flush(6)
     case('max_iter')
        read(tmparg,'(i12)') max_iter
        print *, 'max_iter = ', max_iter; call flush(6)
     case('use_grad_norm')
        read(tmparg,'(i12)') k
        if (k == 1) then
           use_grad_norm = .true.
        else
           use_grad_norm = .false.
        end if
        print *, 'use_grad_norm = ', k; call flush(6)
     case('use_min_dll')
        read(tmparg,'(i12)') k
        if (k == 1) then
           use_min_dll = .true.
        else
           use_min_dll = .false.
        end if
        print *, 'use_min_dll = ', k; call flush(6)
     case('min_grad_norm')
        read(tmparg,'(f15.12)') min_nd
        print *, 'min grad norm = ', min_nd; call flush(6)
     case('min_dll')
        read(tmparg,'(f15.12)') min_dll
        print *, 'min dll = ', min_dll; call flush(6)
     case('do_newton')
        read(tmparg,'(i12)') k
        if (k == 1) then
           do_newton = .true.
        else
           do_newton = .false.
        end if
        print *, 'do_newton = ', k; call flush(6)
     case('newt_start')
        read(tmparg,'(i12)') newt_start
        print *, 'newt_start = ', newt_start; call flush(6)
     case('newt_ramp')
        read(tmparg,'(i12)') newt_ramp
        print *, 'newt_ramp = ', newt_ramp; call flush(6)
     case('do_rho')
        read(tmparg,'(i12)') k
        if (k == 1) then
           dorho = .true.
        else
           dorho = .false.
        end if
        print *, 'do_rho = ', k; call flush(6)
     case('load_rho')
        read(tmparg,'(i12)') k
        if (k == 1) then
           load_rho = .true.
        else
           load_rho = .false.
        end if
        print *, 'load_rho = ', k; call flush(6)
     case('dble_data')
        read(tmparg,'(i12)') k
        if (k == 1) then
           dble_data = .true.
        else
           dble_data = .false.
        end if
        print *, 'dble_data = ', k; call flush(6)
     case('byte_size')
        read(tmparg,'(i12)') nbyte
        print *, 'byte_size = ', nbyte; call flush(6)
     case('load_rej')
        read(tmparg,'(i12)') k
        if (k == 1) then
           load_rej = .true.
        else
           load_rej = .false.
        end if
        print *, 'load_rej = ', k; call flush(6)
     case('print_debug')
        read(tmparg,'(i12)') k
        if (k == 1) then
           print_debug = .true.
        else
           print_debug = .false.
        end if
        print *, 'print_debug = ', k; call flush(6)
     case('update_A')
        read(tmparg,'(i12)') k
        if (k == 1) then
           update_A = .true.
        else
           update_A = .false.
        end if
        print *, 'update_A = ', k; call flush(6)
     case('load_A')
        read(tmparg,'(i12)') k
        if (k == 1) then
           load_A = .true.
        else
           load_A = .false.
        end if
        print *, 'load_A = ', k; call flush(6)
     case('update_c')
        read(tmparg,'(i12)') k
        if (k == 1) then
           update_c = .true.
        else
           update_c = .false.
        end if
        print *, 'update_c = ', k; call flush(6)
     case('load_c')
        read(tmparg,'(i12)') k
        if (k == 1) then
           load_c = .true.
        else
           load_c = .false.
        end if
        print *, 'load_c = ', k; call flush(6)
     case('share_comps')
        read(tmparg,'(i12)') k
        if (k == 1) then
           share_comps = .true.
        else
           share_comps = .false.
        end if
        print *, 'share_comps = ', k; call flush(6)
     case('comp_thresh')
        read(tmparg,'(f15.12)') comp_thresh
        print *, 'comp_thresh = ', comp_thresh; call flush(6)
     case('share_start')
        read(tmparg,'(i12)') share_start
        print *, 'share_start = ', share_start; call flush(6)
     case('share_iter')
        read(tmparg,'(i12)') share_iter
        print *, 'share_int = ', share_iter; call flush(6)
     case('load_comp_list')
        read(tmparg,'(i12)') k
        if (k == 1) then
           load_comp_list = .true.
        else
           load_comp_list = .false.
        end if
        print *, 'load_comp_list = ', k; call flush(6)
     case('do_mean')
        read(tmparg,'(i12)') k
        if (k == 1) then
           do_mean = .true.
        else
           do_mean = .false.
        end if
        print *, 'do_mean = ', k; call flush(6)
     case('do_sphere')
        read(tmparg,'(i12)') k
        if (k == 1) then
           do_sphere = .true.
        else
           do_sphere = .false.
        end if
        print *, 'do_sphere = ', k; call flush(6)
     case('do_approx_sphere')
        read(tmparg,'(i12)') k
        if (k == 1) then
           do_approx_sphere = .true.
        else
           do_approx_sphere = .false.
        end if
        print *, 'do_approx_sphere = ', k; call flush(6)
     case('pcakeep')
        read(tmparg,'(i12)') pcakeep
        print *, 'pcakeep = ', pcakeep; call flush(6)
     case('pcadb')
        read(tmparg,'(f15.12)') pcadb
        print *, 'pcadb = ', pcadb; call flush(6)
     case('invsigmax')
        read(tmparg,'(f15.12)') invsigmax
        print *, 'invsigmax = ', invsigmax; call flush(6)
     case('invsigmin')
        read(tmparg,'(f15.12)') invsigmin
        print *, 'invsigmin = ', invsigmin; call flush(6)
     case('load_all_param')
        read(tmparg,'(i12)') k
        if (k == 1) then
           load_mean = .true.
           load_sphere = .true.
           load_gm = .true.
           load_alpha = .true.
           load_mu = .true.
           load_beta = .true.
           load_rho = .true.
           load_W = .true.
           load_c = .true.
        end if
        print *, 'load_all_param = ', k; call flush(6)
     case('load_sphere')
        read(tmparg,'(i12)') k
        if (k == 1) then
           load_sphere = .true.
        else
           load_sphere = .false.
        end if
        print *, 'load_sphere = ', k; call flush(6)
     case('load_mean')
        read(tmparg,'(i12)') k
        if (k == 1) then
           load_mean = .true.
        else
           load_mean = .false.
        end if
        print *, 'load_mean = ', k; call flush(6)
     case('update_mu')
        read(tmparg,'(i12)') k
        if (k == 1) then
           update_mu = .true.
        else
           update_mu = .false.
        end if
        print *, 'update_mu = ', k; call flush(6)
     case('load_mu')
        read(tmparg,'(i12)') k
        if (k == 1) then
           load_mu = .true.
        else
           load_mu = .false.
        end if
        print *, 'load_mu = ', k; call flush(6)
     case('update_beta')
        read(tmparg,'(i12)') k
        if (k == 1) then
           update_beta = .true.
        else
           update_beta = .false.
        end if
        print *, 'update_beta = ', k; call flush(6)
     case('load_beta')
        read(tmparg,'(i12)') k
        if (k == 1) then
           load_beta = .true.
        else
           load_beta = .false.
        end if
        print *, 'load_beta = ', k; call flush(6)
     case('update_alpha')
        read(tmparg,'(i12)') k
        if (k == 1) then
           update_alpha = .true.
        else
           update_alpha = .false.
        end if
        print *, 'update_alpha = ', k; call flush(6)
     case('load_alpha')
        read(tmparg,'(i12)') k
        if (k == 1) then
           load_alpha = .true.
        else
           load_alpha = .false.
        end if
        print *, 'load_alpha = ', k; call flush(6)
     case('update_gm')
        read(tmparg,'(i12)') k
        if (k == 1) then
           update_gm = .true.
        else
           update_gm = .false.
        end if
        print *, 'update_gm = ', k; call flush(6)
     case('load_gm')
        read(tmparg,'(i12)') k
        if (k == 1) then
           load_gm = .true.
        else
           load_gm = .false.
        end if
        print *, 'load_gm = ', k; call flush(6)
     case('write_nd')
        read(tmparg,'(i12)') k
        if (k == 1) then
           write_nd = .true.
        else
           write_nd = .false.
        end if
        print *, 'write_nd = ', k; call flush(6)
     case('write_LLt')
        read(tmparg,'(i12)') k
        if (k == 1) then
           write_LLt = .true.
        else
           write_LLt = .false.
        end if
        print *, 'write_LLt = ', k; call flush(6)
     case('do_reject')
        read(tmparg,'(i12)') k
        if (k == 1) then
           do_reject = .true.
        else
           do_reject = .false.
        end if
        print *, 'do_reject = ', k; call flush(6)
     case('do_history')
        read(tmparg,'(i12)') k
        if (k == 1) then
           do_history = .true.
        else
           do_history = .false.
        end if
        print *, 'do_history = ', k; call flush(6)
     case('histstep')
        read(tmparg,'(i12)') histstep
        print *, 'histstep = ', histstep; call flush(6)
     case('do_opt_block')
        read(tmparg,'(i12)') k
        if (k == 1) then
           do_opt_block = .true.
        else
           do_opt_block = .false.
        end if
        print *, 'do_opt_block = ', k; call flush(6)
     case('fix_init')
        read(tmparg,'(i12)') k
        if (k == 1) then
           fix_init = .true.
        else
           fix_init = .false.
        end if
        print *, 'fix_init = ', k; call flush(6)
     case('mineig')
        read(tmparg,'(e15.3)') mineig
        print *, 'minimum data covariance eigenvalue = ', mineig; call flush(6)
     case('lrate')
        read(tmparg,'(f15.12)') lrate0
        print *, 'initial lrate = ', lrate0; call flush(6)
        lrate = lrate0
     case('minlrate')
        read(tmparg,'(e15.3)') minlrate
        print *, 'minimum lrate = ', minlrate; call flush(6)
     case('lratefact')
        read(tmparg,'(f15.12)') lratefact
        print *, 'lrate factor = ', lratefact; call flush(6)
     case('rholrate')
        read(tmparg,'(f15.12)') rholrate0
        print *, 'initial rholrate = ', rholrate0; call flush(6)
        rholrate = rholrate0
     case('newtrate')
        read(tmparg,'(f15.12)') newtrate
        print *, 'initial newton lrate = ', newtrate; call flush(6)
     case('rholratefact')
        read(tmparg,'(f15.12)') rholratefact
        print *, 'rho lrate factor = ', rholratefact; call flush(6)
     case('rho0')
        read(tmparg,'(f15.12)') rho0
        print *, 'rho0 = ', rho0; call flush(6)
     case('minrho')
        read(tmparg,'(f15.12)') minrho
        print *, 'min rho = ', minrho; call flush(6)
     case('maxrho')
        read(tmparg,'(f15.12)') maxrho
        print *, 'max rho = ', maxrho; call flush(6)
     case('pdftype')
        read(tmparg,'(i12)') pdftype
        print *, 'pdf type = ', pdftype; call flush(6)
     case('decwindow')
        read(tmparg,'(i12)') decwindow
        print *, 'dec window = ', decwindow; call flush(6)
     case('max_decs')
        read(tmparg,'(i12)') maxdecs
        print *, 'max_decs = ', maxdecs; call flush(6)
     case('writestep')
        read(tmparg,'(i12)') writestep
        print *, 'write step = ', writestep; call flush(6)
     case('kurt_start')
        read(tmparg,'(i12)') chpdfstart
        print *, 'kurt_start = ', chpdfstart; call flush(6)
     case('num_kurt')
        read(tmparg,'(i12)') maxchpdf
        print *, 'num kurt = ', maxchpdf; call flush(6)
     case('kurt_int')
        read(tmparg,'(i12)') chpdfint
        print *, 'kurt interval = ', chpdfint; call flush(6)
     case('rejsig')
        read(tmparg,'(f15.12)') rejsig
        print *, 'reject sigma = ', rejsig; call flush(6)
     case('numrej')
        read(tmparg,'(i12)') maxrej
        print *, 'num reject = ', maxrej; call flush(6)
     case('rejstart')
        read(tmparg,'(i12)') rejstart
        print *, 'reject start = ', rejstart; call flush(6)
     case('rejint')
        read(tmparg,'(i12)') rejint
        print *, 'reject interval = ', rejint; call flush(6)
     case('doscaling')
        read(tmparg,'(i12)') k
        if (k == 1) then
           doscaling = .true.
        else
           doscaling = .false.
        end if
        print *, 'doscaling = ', k; call flush(6)
     case('scalestep')
        read(tmparg,'(i12)') scalestep
        print *, 'scalestep = ', scalestep; call flush(6)
     case('max_threads')
        read(tmparg,'(i12)') max_thrds
        print *, 'max_thrds = ', max_thrds; call flush(6)
     case('outdir')
        read(tmparg,'(a)') outdirparam
     case('indir')
        read(tmparg,'(a)') indirparam
     end select

  end do

  if (pdftype .ne. 0) then
     dorho = .false.
  end if


end subroutine get_cmd_args


!----------------------------------------------------------------------

subroutine tststat(stat)
integer :: stat
if (stat > 0) then
   print *, 'allocation error'
   stop
end if
end subroutine


!----------------------------------------------------------------------
end program main
!----------------------------------------------------------------------

