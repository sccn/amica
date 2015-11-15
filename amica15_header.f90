
type data_segment
   integer :: filenum, lastdim, numgood
   double precision, dimension(:,:), pointer :: data
   logical, dimension(:), pointer :: gooddata(:)
   integer, dimension(:), pointer :: goodinds(:)
   double precision, dimension(:), pointer :: loglik(:)
   double precision, dimension(:,:), pointer :: modloglik(:,:)
end type data_segment

type(data_segment), allocatable :: dataseg(:)

logical :: do_approx_sphere = .true.
logical :: do_sphere = .true., do_mean = .true., dorho = .true., doscaling = .true., print_debug = .false.
logical :: leave = .false., update_mu = .true., update_beta = .true., update_A = .true., fix_init = .false.
logical :: doPCA = .false., load_rho = .false., load_A = .false., load_mu = .false., load_W = .false.
logical :: load_beta = .false., load_alpha = .false., load_gm = .false., load_comp_list = .false.
logical :: update_alpha = .true., update_gm = .true., startover = .false., share_comps = .false.
logical :: write_nd = .false., do_newton = .false., free_pass = .false.
logical :: do_reject = .false., do_choose_pdfs = .false., do_opt_block = .true.
logical :: load_c = .false., update_c = .true., load_rej = .false., redoiter = .false.
logical :: write_LLt = .true., load_mean = .false., load_sphere = .false.
logical :: posdef = .true., declrate = .false., no_newt = .false., do_update, do_history = .false., dble_data = .false.
logical :: use_grad_norm = .true., use_min_dll = .true.
logical, allocatable :: update_component(:,:), comp_used(:)

integer, allocatable, target :: int_array(:)
real, allocatable, target :: real_array(:)
double precision, allocatable, target :: dble_array(:)

DOUBLE PRECISION, ALLOCATABLE :: A(:,:), Ax(:,:), dAk(:,:), dA(:,:,:), W(:,:,:), Wtmp(:,:), Wtmp2(:,:,:)
DOUBLE PRECISION, ALLOCATABLE :: c(:,:), cx(:,:), wc(:,:), dc_numer_tmp(:,:), dc_numer(:,:), dc_denom_tmp(:,:), dc_denom(:,:)
double precision, allocatable :: S(:,:), dCov(:,:), Spinv(:,:), Spinv2(:,:), Stmp(:,:), Stmp2(:,:), Stmp3(:,:), dWtmp(:,:,:), ipivnx(:), ipivnw(:)
DOUBLE PRECISION, ALLOCATABLE :: y(:,:,:,:), logab(:), g(:,:), b(:,:,:), sUtmp(:,:), sVtmp(:,:)
double precision, allocatable :: z(:,:,:,:), z0(:,:), u(:), ufp(:), fp(:)
DOUBLE PRECISION, ALLOCATABLE :: v(:,:), nd(:,:), tmpy(:)
double precision, allocatable :: xtmp(:,:), blk_time(:)
DOUBLE PRECISION, ALLOCATABLE :: P(:), Ptmp(:,:), Ptmp2(:,:), Pmax(:), Dsum(:), Dmax(:), Dmin(:)
double precision, allocatable :: wr(:), Dtemp(:), LL(:), m2(:,:), m4(:,:), m2sum(:,:), m4sum(:,:)

DOUBLE PRECISION, ALLOCATABLE :: gm(:), alpha(:,:), mu(:,:), sbeta(:,:), rho(:,:), zeta(:)
DOUBLE PRECISION, ALLOCATABLE :: mutmp(:,:), sbetatmp(:,:)
double precision, allocatable :: rhotmp(:,:), tmpvec(:), tmpvec2(:), tmpvec3(:), tmpvec4(:), ztmp(:), utmp(:), vtmp(:)

DOUBLE PRECISION, ALLOCATABLE :: lambda(:,:), kappa(:,:), sigma2(:,:), baralpha(:,:,:)

DOUBLE PRECISION, ALLOCATABLE :: dgm_numer_tmp(:), dgm_numer(:)
DOUBLE PRECISION, ALLOCATABLE :: dalpha_numer_tmp(:,:), dalpha_denom_tmp(:,:), dalpha_numer(:,:), dalpha_denom(:,:)
DOUBLE PRECISION, ALLOCATABLE :: dmu_numer_tmp(:,:), dmu_denom_tmp(:,:), dmu_numer(:,:), dmu_denom(:,:)
DOUBLE PRECISION, ALLOCATABLE :: dbeta_numer_tmp(:,:), dbeta_denom_tmp(:,:), dbeta_numer(:,:), dbeta_denom(:,:)
double precision, allocatable :: drho_numer_tmp(:,:), drho_denom_tmp(:,:), drho_numer(:,:), drho_denom(:,:)

DOUBLE PRECISION, ALLOCATABLE :: dbaralpha_numer_tmp(:,:,:), dbaralpha_denom_tmp(:,:,:), dbaralpha_numer(:,:,:), dbaralpha_denom(:,:,:)
DOUBLE PRECISION, ALLOCATABLE :: dlambda_numer_tmp(:,:,:), dlambda_denom_tmp(:,:,:), dlambda_numer(:,:,:), dlambda_denom(:,:,:)
DOUBLE PRECISION, ALLOCATABLE :: dkappa_numer_tmp(:,:,:), dkappa_denom_tmp(:,:,:), dkappa_numer(:,:,:), dkappa_denom(:,:,:)
DOUBLE PRECISION, ALLOCATABLE :: dsigma2_numer_tmp(:,:), dsigma2_denom_tmp(:,:), dsigma2_numer(:,:), dsigma2_denom(:,:)

DOUBLE PRECISION, ALLOCATABLE :: mean(:), meantmp(:), eigs(:)
double precision, allocatable :: eigv(:)

integer, allocatable :: comp_list(:,:)

INTEGER, ALLOCATABLE :: work(:), num_samples(:), field_dim(:), node_thrds(:), node_blocks(:), pdtype(:,:)
INTEGER, ALLOCATABLE :: num_dir_files(:), seg_list(:,:,:), blk_size(:), blocks_in_sample(:), nsub(:)

double precision :: mincond = 1.0e-15, minrho = 1.0, maxrho = 2.0, minlog = -1500.0, LLinc, maxdble = 1.0e32, mineigv = 1.0e-15, mineig = 1.0e-15;
double precision :: ndtmpsum, ndtmpsumsum, rejsig = 3.0, LLtmp, LLtmp2, lrate0, pcadb
double precision :: lrate = 0.1, minlrate = 1.0e-12, lratefact = 0.5, rho0 = 1.5, rholrate = 0.05, rholrate0 = 0.05
double precision :: llmaxmax, llminmin, m2tmp, m4tmp, Ndiv, invsigmax = 1000.0
double precision :: bex, datafrac, t0, rholratefact = 0.1, sldet, Anrmk
double precision :: llvar = 0.0, llsig = 0.0, llvarsum = 0.0, llmean = 0.0
double precision :: llmeansum = 0.0, llmax, llmin, invsigmin = 0.0001
double precision :: newtrate = 0.5, epsdble = 1.0e-16, sk1, sk2, natrate, minhess = 1.0e-5
double precision :: usum, tmpsum, vsum, dkap, comp_thresh = 0.99, min_dll = 1.0e-9, min_nd = 1.0e-7

integer :: num_comps = -1, num_mix = 3, num_mix_init = 3, share_iter = 100, share_start = 100, blk_min = 128, blk_max = 1024
integer :: data_dim, N1, h, hh, t, flen, numgood, numgoodsum, ngood, maxchpdf, blk_step = 128
integer :: numchpdf, chpdfint, chpdfstart, pcakeep
integer :: dft_length = 256, numsegs, bsize, decwindow = 1, numrej = 0
integer :: rejstart = 2, rejint = 3, maxrej = 1, maxrestart = 10
integer :: block_size = 128, fieldsize, fnum, lastblocksize, numeigs
integer :: seg_comm, node_procs, node_comm, sampsize, datsize, lastd, snum
integer :: nodefiles, firstfile, seg, tblksize
integer :: mw, nw, mx, nx, cnt, ntmp, all_blks, ldim
integer :: xstrt, xstp, x0strt, x0stp, bstrt, bstp, quittmp, quit
integer :: field_blocksize, data_blocksize, filter_length, bcnt
integer :: spherescale = 1, scalestep = 1, newt_start = 20, newt_ramp = 10

INTEGER :: num_files = 0, num_files_tokens, tmp_read, tmpind, gotfile
integer :: mstrt, mstp, fixcoords = 0, iargc, outstep = 1
INTEGER :: myrank, node_rank, seg_rank, seg_nodes, thrdnum, tot_procs
integer :: num_procs, ierr, tot_thrds, num_thrds, max_thrds = 24, num_thrds_used
INTEGER :: num_models = 1, max_iter, pdftype = 1
integer :: num_tot, tot_dim, offset, blk, samp
INTEGER :: random_init_W = 0, nbyte = 4, read_Winit = 0
integer :: coststep = 1, writestep = 100, tot_blks, nblks
INTEGER :: seed(2) = (/ 123456, 654321 /)
INTEGER :: numargs, argnum, filenum, filestart, filestop, sampnum, sampstart
integer :: sampstop, blknum, blocknum, fld1, fld2, num_blocks
INTEGER :: i, j, k
integer :: ii, jj, kk, c0, c1, c2
integer :: counts_per_sec=0, cnt1=0, cnt2=0
INTEGER :: iter, len, fh, info, lwork, lstate=16, state(16), lseed = 1
integer :: host_num, ip(4), name_len
integer :: status(5)
integer :: req1, req2

integer :: numdecs = 0, maxdecs = 5, numrestarts = 0, maxrestarts = 3, restartiter = 10, histstep = 10
integer :: numincs = 0, maxincs = 5
integer, parameter :: MAX_CHARS = 500
CHARACTER(len=500), ALLOCATABLE :: arglist(:), infile(:)
CHARACTER(len=500) :: W_infile='', infilename='', argname='', tmparg=''
character(len=500) :: tmpdirarg = '', outdirparam = '', indirparam = ''
CHARACTER(len=500) :: tmpstring='', tmpstring2='', host_name=''
CHARACTER(len=6) :: outdir = 'output'
character(len=1) :: tmpchar
character(len=8) :: ipstr
CHARACTER(len=500) :: tmpdir = '/tmp/'

CHARACTER(len=1), PARAMETER :: outfile_w = 'W'
CHARACTER(len=1), PARAMETER :: outfile_c = 'c'
CHARACTER(len=9), PARAMETER :: outfile_comp_list = 'comp_list'
CHARACTER(len=1), PARAMETER :: outfile_sphere = 'S'
CHARACTER(len=1), PARAMETER :: outfile_A = 'A'
CHARACTER(len=2), PARAMETER :: outfile_LL = 'LL'
CHARACTER(len=3), PARAMETER :: outfile_LLt = 'LLt'
CHARACTER(len=2), PARAMETER :: outfile_nd = 'nd'
CHARACTER(len=2), PARAMETER :: outfile_gamma = 'gm'
CHARACTER(len=2), PARAMETER :: outfile_mu = 'mu'
CHARACTER(len=5), PARAMETER :: outfile_sbeta = 'sbeta'
CHARACTER(len=5), PARAMETER :: outfile_alpha = 'alpha'
CHARACTER(len=3), PARAMETER :: outfile_rho = 'rho'
CHARACTER(len=4), PARAMETER :: outfile_mean = 'mean'
CHARACTER(len=7), PARAMETER :: printoutfile = 'out.txt'
CHARACTER(len=8) date
CHARACTER(len=10) time
CHARACTER(len=5) iterstr
data date /'12345678'/, time /'1234567890'/, iterstr /'_____'/

interface
   function DNRM2(n,x,inc)
      integer :: n, inc
      double precision :: DNRM2, x(n)
   end function DNRM2
end interface
interface
   function DDOT(n,x,incx,y,incy)
      integer :: n, incx, incy
      double precision :: DDOT, x(n), y(n)
   end function DDOT
end interface
interface
   function fastpow(x,y)
      double precision :: fastpow,x,y
   end function fastpow
end interface
interface
   function fastlog(x)
      double precision :: fastlog,x
   end function fastlog
end interface
interface
   function fastexp(x)
      double precision :: fastexp,x
   end function fastexp
end interface
interface
   subroutine omp_set_num_threads( num_thr )
      integer :: num_thr
   end subroutine omp_set_num_threads
end interface
interface
   function omp_get_num_procs()
      integer :: omp_get_num_procs
   end function omp_get_num_procs
end interface
interface
   function omp_get_num_threads()
      integer :: omp_get_num_threads
   end function omp_get_num_threads
end interface
interface
   function omp_get_thread_num()
      integer :: omp_get_thread_num
   end function omp_get_thread_num
end interface
interface
   subroutine omp_set_num_threads__( num_thr )
      integer :: num_thr
   end subroutine omp_set_num_threads__
end interface
interface
   function omp_get_num_procs__()
      integer :: omp_get_num_procs__
   end function omp_get_num_procs__
end interface
interface
   function omp_get_num_threads__()
      integer :: omp_get_num_threads__
   end function omp_get_num_threads__
end interface
interface
   function omp_get_thread_num__()
      integer :: omp_get_thread_num__
   end function omp_get_thread_num__
end interface
interface
   subroutine omp_set_num_threads_( num_thr )
      integer :: num_thr
   end subroutine omp_set_num_threads_
end interface
interface
   function omp_get_num_procs_()
      integer :: omp_get_num_procs_
   end function omp_get_num_procs_
end interface
interface
   function omp_get_num_threads_()
      integer :: omp_get_num_threads_
   end function omp_get_num_threads_
end interface
interface
   function omp_get_thread_num_()
      integer :: omp_get_thread_num_
   end function omp_get_thread_num_
end interface
interface
   function nag_gamma(x)
      double precision :: x, nag_gamma
   end function nag_gamma
end interface
