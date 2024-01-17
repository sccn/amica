% runamica17_nsg()
%   Perform AMICA -- adaptive ICA using multiple models with shared components
%   on Neuroscience Gateway (NSG)
%
% Usage:
%         >> runamica17_nsg(dat,'Key1',Value1,...);
%
% Inputs:
%
%   dat            EEG structure, or a file name (including path) of file containing
%                  floating point data, or a .set filename, or a data
%                  matrix (chans x frames)
%   'Key1'         Keyword string for argument 1
%   Value1         Value of argument 1
%   ...            ...
%
% Outputs:
%
%   Stored in outdir, default to 'amicaouttmp' folder under the current
%   wokring directory. Use EEG = eeg_loadamica(outdir) to read the output files.
%
% Optional keyword inputs:
%
%   outdir              name of directory to write output (does not have to exist), def=pwd/amicaouttmp/
%   indir               optional input directory from which to load init
%   numprocs            number of nodes to use
%   max_threads         maximum number of threads to use if run locally, def=4
%   batch               use a batch script, 0 or 1, def=0
%   num_chans           number of channels in data (only needed if dat input is a filename)
%   block_size          matrix multiplication block size
%   do_opt_block        optimize block size, def=0
%   blk_min             min of block optimization range
%   blk_max             max of block optimization range
%   blk_step            step for block optimization range
%   num_models          number of models to learn, def=1
%   num_mix_comps       number of mixture components in source model, def=3
%   max_iter            maximum number of iterations to perform, def=2500
%   lrate               initial learning rate for natural gradient, def=0.1
%   minlrate            lrate after which to stop, def=1e-8
%   use_grad_norm       stop when average parameter update is les than min_grad_norm, def=1
%   min_grad_norm       min_grad_norm value, def=1e-7
%   use_min_dll         stop when log likelihood change is les than min_dll, def=1
%   min_dll             min_dll value, def=1e-10
%   do_history          save a record of the paramter values over iterations def=0
%   histstep            iteration interval at which to save for history, def=10
%   lratefact           multiplicative factor by which to decrease lrate, def=0.5
%   rholrate            initial lrate for shape parameters, def=0.05
%   rho0                initial shape parameter value, def=1.5
%   minrho              minimum shape parameter value, def=1.0
%   maxrho              maximum shape parameter value, def=2.0
%   rholratefact        multiplicative factor by which to dec rholrate, def=0.5
%   do_newton           flag for newton method, default = 1 (do newton)
%   newt_start          for newton method, iter at which to start newton, def=50
%   newtrate            for newton method, lrate for newton iterations, def=1.0
%   newt_ramp           for newton method, number of iter to ramp up to newtrate, def=10
%   writestep           iteration interval between output writes, def=10
%   write_nd            flag to write history of component update norms, def=1
%   write_llt           flag to write model log likelihoods of time points, def=1
%   do_reject           flag for doing rejection of time points, def=0
%   numrej              number of rejections to perform, def=3
%   rejsig              number of std dev of likelihood below which to reject data, def=3
%   rejstart            for rejection, iteration at which to start reject, def=1
%   rejint              for rejection, iteration interval between reject, def=1
%   decwindow           moving average window to detect likelihood decrease, def=1
%   update_A            flag to update mixing matrices, def=1
%   update_c            flag to update model centers, def=1
%   update_gamma        flag to update model probabilities, def=1
%   update_alpha        flag to update source mixture proportions, def=1
%   update_mu           flag to update source mixture mixture locations, def=1
%   update_sbeta        flag to update source mixture scales, def=1
%   invsigmax           maximum value of inverse scale parameters, def=100.0
%   invsigmin           minimum value of inverse scale parameters, def=0.00001
%   do_rho              flag to update shape parameters, def=1
%   load_comp_list      flag to load component assignment list, def=0
%   load_rej            flag to load LLt to get rejections from, def=0
%   load_param          flag to load parameters, def=0
%   do_mean             flag to remove mean from data, def=1
%   do_sphere           flag to sphere data before ica, def=1
%   doPCA               flag to to PCA dimensionalit reduction, def=0
%   pcakeep             for PCA reduction, number of components to keep, def=chans
%   doscaling           flag to rescale unmixing matrix rows to unit norm, def=1
%   scalestep           iteration interval at which to rescale unmixing rows, def=1
%
% Outputs:
%   EEG - if there is an input EEG, it is loaded with the ICA output, first
%   model in case the
%   To load output use the function loadmodout() after job ends:
%                         
%       mods = eeg_loadamica(outdir);
%
%   mods is a structure containing the output components and density models.
%       mods.A(:,:,h)       the components for model h.
%       mods.varord(:,h)    the index order of the components in variance order
%       mods.Lht            the likelihood of time points for each model (if set)
%       mods.LL             the history of the log likelihood over iterations
%       mods.c(:,h)         the center for model h
%       mods.W(:,:,h)       the unmixing matrix for model h
%       mods.S              the sphering matrix, first num_pcs rows used
%                       
% See also: eeg_loadamica()
%
%

function EEG = runamica17_nsg(dat,varargin)

if nargin < 1
    help runamica17_nsg;
end

% define output directory - create a folder under the current working directory
outdir = [ pwd filesep 'amicaouttmp' filesep ];

% file path of the AMICA binary file
binfile = 'amica17nsg';  % AMICA binary file for Linux
basepath = [fileparts(mfilename('fullpath')) filesep];  % directory where this .m file and the binary file are located
% basepath_nsg = '/home/nsguser/applications/eeglab_new/plugins/amica1.5.1/';   % filepath of amica plugin on NSG 

if exist([basepath binfile],'file')  % use bin file in amica plugin on NSG
    AMBIN = [basepath binfile];
else  % cannot find binary file
    disp('No binary file found. Exit runamica17_nsg.m')
    return
end

% run locally (single node) or globally (multiple nodes)
RUN_LOCALLY = 0;


%%%%%%%%%% Define AMICA parameters %%%%%%%%%%%%%%%

numprocs = 1;           % # nodes. When set to > 1, RUN_LOCALLY is set to 0
max_threads = 4;       % # threads. For NSG, default to max-number (24) threads
bytesize = 4;           % convert data to floating point format

batch = 0;

fix_init = 0;

block_size = 128; % initial block size
do_opt_block = 0; % optimize block size
blk_min = 16;
blk_step = 16;
blk_max = 128;

num_models = 1;         % Number of ICA models
num_mix_comps = 3;      % Number of generalized gaussians mixture for p(s)
pdftype = 0;
max_iter = 2500;

mineig = 1e-9;

do_history = 0;
histstep = 10;

lrate = 0.01;
minlrate = 1e-8;
lratefact = 0.5;

rholrate = 0.05;
rho0 = 1.5;
minrho = 1.0;
maxrho = 2.0;
rholratefact = 0.5;

do_newton = 1;
newt_start = 50;
newtrate = 1.0;
newt_ramp = 10;

do_reject = 0;
numrej = 10;
rejsig = 3.0;
rejstart = 1;
rejint = 3;

use_min_dll = 1;
min_dll = 1e-10;
use_grad_norm = 1;
min_grad_norm = 1e-7;

writestep = 1000;
write_nd = 0;
write_LLt = 1;
decwindow = 1;
maxdecs = 3;

update_A = 1;
update_c = 1;
update_gamma = 1;
update_alpha = 1;
update_mu = 1;
update_sbeta = 1;
invsigmax = 100.0;
invsigmin = 1e-8;
do_rho = 1;

load_comp_list = 0;
load_rej = 0;
load_W = 0;
load_c = 0;
load_gm = 0;
load_alpha = 0;
load_mu = 0;
load_sbeta = 0;
load_rho = 0;
load_param = 0;

do_mean = 1;
do_sphere = 1;
doPCA = 1;
pcadb = 30.0;
doscaling = 1;
scalestep = 1;


%%%%%%%%%% Collect keywords and values from argument list %%%%%%%%%%%%%%%

chans = -1;
pcakeep = -1;

   for i = 1:2:length(varargin) % for each Keyword
      Keyword = varargin{i};
      Value = varargin{i+1};
      if ~isstr(Keyword)
         fprintf('runamica(): keywords must be strings')
         return
      end
      Keyword = lower(Keyword); % convert upper or mixed case to lower

      if strcmp(Keyword,'outdir')
         if ~isstr(Value)
            fprintf('runamica(): outdir must be a string');
            return
         else
            outdir = Value;
         end
      elseif strcmp(Keyword,'numprocs')
         if isstr(Value)
            fprintf('runamica(): numprocs must be a positive integer');
            return
         else
            numprocs = Value;
         end
      elseif strcmp(Keyword,'num_chans')
         if isstr(Value)
            fprintf('runamica(): num_chans must be a positive integer');
            return
         else
            chans = Value;
         end
      elseif strcmp(Keyword,'num_models')
         if isstr(Value)
            fprintf('runamica(): num_models must be a positive integer');
            return
         else
            num_models = Value;
         end
      elseif strcmp(Keyword,'indir')
         if ~isstr(Value)
            fprintf('runamica(): indir must be a string');
            return
         else
            indir = Value;
            load_param = 1;
            load_W = 1;
            load_c = 1;
            load_gm = 1;
            load_alpha = 1;
            load_mu = 1;
            load_sbeta = 1;
            load_rho = 1;
            load_comp_list = 1;
         end
      elseif strcmp(Keyword,'num_mix_comps')
         if isstr(Value)
            fprintf('runamica(): num_mix_comps must be a positive integer');
            return
         else
            num_mix_comps = Value;
         end
      elseif strcmp(Keyword,'pdftype') 
         if isstr(Value)
            fprintf('runamica(): pdftype should be either 0, 1, 2, or 3');
            return
         else
            pdftype = Value;
         end
      elseif strcmp(Keyword,'max_iter')
         if isstr(Value)
            fprintf('runamica(): max_iter must be a positive integer');
            return
         else
            max_iter = Value;
         end
      elseif strcmp(Keyword,'batch')
         if isstr(Value)
            fprintf('runamica(): batch must be 0 or 1');
            return
         else
            batch = Value;
         end
      elseif strcmp(Keyword,'do_opt_block')
         if isstr(Value)
            fprintf('runamica(): do_opt_block must be 0 or 1');
            return
         else
            do_opt_block = Value;
         end
      elseif strcmp(Keyword,'blk_min')
         if isstr(Value)
            fprintf('runamica(): blk_min must be a positive integer');
            return
         else
            blk_min = Value;
         end
      elseif strcmp(Keyword,'blk_step')
         if isstr(Value)
            fprintf('runamica(): blk_step must be a positive integer');
            return
         else
            blk_step = Value;
         end
      elseif strcmp(Keyword,'blk_max')
         if isstr(Value)
            fprintf('runamica(): blk_max must be a positive integer');
            return
         else
            blk_max = Value;
         end
      elseif strcmp(Keyword,'block_size')
         if isstr(Value)
            fprintf('runamica(): block_size must be a positive integer');
            return
         else
            block_size = Value;
         end
      elseif strcmp(Keyword,'lrate')
         if isstr(Value)
            fprintf('runamica(): lrate must be a positive number');
            return
         else
            lrate = Value;
         end
      elseif strcmp(Keyword,'minlrate')
         if isstr(Value)
            fprintf('runamica(): minlrate should be lrate to stop at');
            return
         else
            minlrate = Value;
         end
      elseif strcmp(Keyword,'use_min_dll')
         if isstr(Value)
            fprintf('runamica(): use_min_dll should be 0 or 1');
            return
         else
            use_min_dll = Value;
         end
      elseif strcmp(Keyword,'min_dll')
         if isstr(Value)
            fprintf('runamica(): min_dll should be LL change to stop at');
            return
         else
            min_dll = Value;
         end
      elseif strcmp(Keyword,'use_grad_norm')
         if isstr(Value)
            fprintf('runamica(): use_grad_norm should be 0 or 1');
            return
         else
            use_grad_norm = Value;
         end
      elseif strcmp(Keyword,'min_grad_norm')
         if isstr(Value)
            fprintf('runamica(): min_grad_norm should be weight change norm to stop at');
            return
         else
            min_grad_norm = Value;
         end
      elseif strcmp(Keyword,'mineig')
         if isstr(Value)
            fprintf('runamica(): mineig should be a positive number');
            return
         else
            mineig = Value;
         end
      elseif strcmp(Keyword,'lratefact')
         if isstr(Value)
            fprintf('runamica(): lratefact must be a positive number');
            return
         else
            lratefact = Value;
         end
      elseif strcmp(Keyword,'rholrate')
         if isstr(Value)             
            fprintf('runamica(): rholrate must be a number');
            return
         else
            rholrate = Value;
         end
      elseif strcmp(Keyword,'rho0')
         if isstr(Value)
            fprintf('runamica(): rho0 must be a number');
            return
         else
            rho0 = Value;
         end
      elseif strcmp(Keyword,'minrho')
         if isstr(Value)
            fprintf('runamica(): minrho must be a number');
            return
         else
            minrho = Value;
         end
      elseif strcmp(Keyword,'maxrho')
         if isstr(Value)
            fprintf('runamica(): maxrho must be a number');
            return
         else
            maxrho = Value;
         end
      elseif strcmp(Keyword,'rholratefact')
         if isstr(Value)
            fprintf('runamica(): rholratefact must be a number');
         else 
            rholratefact = Value;
         end
      elseif strcmp(Keyword,'invsigmax')
         if isstr(Value)
            fprintf('runamica(): invsigmax must be a number');
         else 
            invsigmax = Value;
         end
      elseif strcmp(Keyword,'invsigmin')
         if isstr(Value)
            fprintf('runamica(): invsigmin must be a number');
         else 
            invsigmin = Value;
         end
      elseif strcmp(Keyword,'fix_init')
         if isstr(Value)
            fprintf('runamica(): fix_init should be 0 or 1');
            return
         else             
            fix_init = Value;
         end
      elseif strcmp(Keyword,'do_newton')
         if isstr(Value)
            fprintf('runamica(): do_newton should be 0 or 1');
            return
         else             
            do_newton = Value;
         end
      elseif strcmp(Keyword,'update_A') %#ok<STCUL> 
         if isstr(Value)
            fprintf('runamica(): update_A should be 0 or 1');
            return
         else             
            update_A = Value;
         end
      elseif strcmp(Keyword,'update_c')
         if isstr(Value)
            fprintf('runamica(): update_c should be 0 or 1');
            return
         else             
            update_c = Value;
         end
      elseif strcmp(Keyword,'update_gamma')
         if isstr(Value)
            fprintf('runamica(): update_gamma should be 0 or 1');
            return
         else             
            update_gamma = Value;
         end
      elseif strcmp(Keyword,'update_alpha')
         if isstr(Value)
            fprintf('runamica(): update_alpha should be 0 or 1');
            return
         else             
            update_alpha = Value;
         end
      elseif strcmp(Keyword,'update_mu')
         if isstr(Value)
            fprintf('runamica(): update_mu should be 0 or 1');
            return
         else             
            update_mu = Value;
         end
      elseif strcmp(Keyword,'update_sbeta')
         if isstr(Value)
            fprintf('runamica(): update_sbeta should be 0 or 1');
            return
         else             
            update_sbeta = Value;
         end
      elseif strcmp(Keyword,'do_rho')
         if isstr(Value)
            fprintf('runamica(): do_rho should be 0 or 1');
            return
         else             
            do_rho = Value;
         end
      elseif strcmp(Keyword,'load_param')
         if isstr(Value)
            fprintf('runamica(): load_param should be 0 or 1');
            return
         else
            load_param = Value;
            if load_param == 1
                load_W = 1;
                load_c = 1;
                load_gm = 1;
                load_alpha = 1;
                load_mu = 1;
                load_sbeta = 1;
                load_rho = 1;  
                load_comp_list = 1;
            end
         end         
      elseif strcmp(Keyword,'load_rej')
         if isstr(Value)
            fprintf('runamica(): load_rej should be 0 or 1');
            return
         else
            load_rej = Value;
         end         
      elseif strcmp(Keyword,'write_nd')
         if isstr(Value)
            fprintf('runamica(): write_nd should be 0 or 1');
            return
         else             
            write_nd = Value;
         end
      elseif strcmp(Keyword,'write_llt')
         if isstr(Value)
            fprintf('runamica(): write_llt should be 0 or 1');
            return
         else             
            write_LLt = Value;
         end
      elseif strcmp(Keyword,'newtrate')
         if isstr(Value)
            fprintf('runamica(): newtrate must be a number');
            return
         else
            newtrate = Value;
         end
      elseif strcmp(Keyword,'newt_start')
         if isstr(Value)
            fprintf('runamica(): newt_start must be a number');
            return
         else
            newt_start = Value;
         end
      elseif strcmp(Keyword,'newt_ramp')
         if isstr(Value)
            fprintf('runamica(): newt_ramp must be a number');
            return
         else
            newt_ramp = Value;
         end
      elseif strcmp(Keyword,'do_reject')
         if isstr(Value)
            fprintf('runamica(): do_reject must be 0 or 1');
            return
         else
            do_reject = Value;
         end
      elseif strcmp(Keyword,'numrej')
         if isstr(Value)
            fprintf('runamica(): numrej must be a number');
            return
         else
            numrej = Value;
         end
      elseif strcmp(Keyword,'rejstart')
         if isstr(Value)
            fprintf('runamica(): rejstart must be a number');
            return
         else
            rejstart = Value;
         end
      elseif strcmp(Keyword,'rejsig')
         if isstr(Value)
            fprintf('runamica(): rejsig must be a number');
            return
         else
            rejsig = Value;
         end
      elseif strcmp(Keyword,'rejint')
         if isstr(Value)
            fprintf('runamica(): rejint must be a number');
            return
         else
            rejint = Value;
         end                  
      elseif strcmp(Keyword,'do_history')
         if isstr(Value)
            fprintf('runamica(): do_history must be 0 or 1');
            return
         else
            do_history = Value;
         end
      elseif strcmp(Keyword,'histstep')
         if isstr(Value)
            fprintf('runamica(): histstep must be a number');
            return
         else
            histstep = Value;
         end
      elseif strcmp(Keyword,'max_threads')
         if isstr(Value)
            fprintf('runamica(): max_threads must be a number');
            return
         else             
             max_threads = Value;             
         end
      elseif strcmp(Keyword,'writestep')
         if isstr(Value)
            fprintf('runamica(): writestep must be a number');
            return
         else
            writestep = Value;
         end         
      elseif strcmp(Keyword,'decwindow')
         if isstr(Value)
            fprintf('runamica(): decwindow must be a number');
            return
         else
            decwindow = Value;
         end         
      elseif strcmp(Keyword,'maxdecs')
         if isstr(Value)
            fprintf('runamica(): maxdecs must be a number');
            return
         else
            maxdecs = Value;
         end         
      elseif strcmp(Keyword,'pcakeep')
         if isstr(Value)
            fprintf('runamica(): pcakeep must be a number');
            return
         else
            pcakeep = Value;
         end         
      elseif strcmp(Keyword,'bytesize')
         if isstr(Value)
            fprintf('runamica(): bytesize must be a number');
            return
         else
            bytesize = Value;
         end         
      elseif strcmp(Keyword,'doscaling')
         if isstr(Value)
            fprintf('runamica(): doscaling must be a number');
            return
         else
            doscaling = Value;
         end      
      elseif strcmp(Keyword,'do_sphere'), %added by JRI 8/21/13
         if isstr(Value)
            fprintf('runamica(): do_sphere must be a number');
            return
         else
            do_sphere = Value;
         end    
      else
         fprintf(['runamica(): unknown flag: ' Keyword])
         return
      end
   end
   

% get the number of channels and frames
if isnumeric(dat)
    chans = size(dat,1);
    frames = size(dat,2);
elseif isstruct(dat)
    chans = dat.nbchan;
    frames = dat.pnts * dat.trials;
    file = [dat.filepath filesep dat.datfile];
else
    if dat(1) == '.'
        dat = [pwd dat(2:end)];        
    end
    if strcmp(dat(end-3:end),'.set')
        EEG = pop_loadset(dat);
        chans = EEG.nbchan;
        frames = EEG.pnts * EEG.trials;
        file = [EEG.filepath filesep EEG.datfile];
    elseif chans == -1
        error('Need num_chans keyword argument with number of channels, unknown file type');
    else
        fid = fopen(dat,'r');
        dattmp = fread(fid,[chans inf],'float');
        frames = size(dattmp,2);
        fclose(fid);
        clear dattmp;
        file = dat;
    end
end

if isnumeric(dat)
    filename = ['tmpdata' num2str(round(rand(1)*100000)) '.fdt' ];
    disp(['Writing data file: ' fullfile(pwd,filename)]);
    fid = fopen(filename, 'w');
    fwrite(fid, dat, 'float');
    fclose(fid)
    chans = size(dat,1);
    frames = size(dat,2);
    file = fullfile(pwd,filename);
end            
   

%%%%%%%%%%%%%%%%%%%% define pars for multi-nodes run %%%%%%%%%%%%%%%%%%%%%

%if numprocs > 1
%    RUN_LOCALLY = 0;
%    if npernode > max_threads
%        disp('Setting npernode to max_threads')
%        npernode = max_threads;
%    end
%end


%%%%%%%%%%%%%%%%%%%%% create the param file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if outdir(1) == '.'
    outdir = [pwd outdir(2:end)];
end
if outdir(end) ~= filesep
    outdir(end+1) = filesep;
end
system(['mkdir ' outdir]);

fid = fopen([outdir 'input.param'],'w');
if fid < 1
    errordlg(['cannot create file in ' outdir],'Bad Input','modal');
    return;
end

fprintf(fid,'files %s\n',file);
fprintf(fid,'outdir %s\n',outdir);
fprintf(fid,'block_size %d\n',block_size);
fprintf(fid,'do_opt_block %d\n',do_opt_block);
fprintf(fid,'blk_min %d\n',blk_min);
fprintf(fid,'blk_step %d\n',blk_step);
fprintf(fid,'blk_max %d\n',blk_max);
fprintf(fid,'num_models %d\n',num_models);
fprintf(fid,'max_threads %d\n',max_threads);
fprintf(fid,'use_min_dll %d\n',use_min_dll);
fprintf(fid,'min_dll %d\n',min_dll);
fprintf(fid,'use_grad_norm %d\n',use_grad_norm);
fprintf(fid,'min_grad_norm %d\n',min_grad_norm);
fprintf(fid,'num_mix_comps %d\n',num_mix_comps);
fprintf(fid,'pdftype %d\n',pdftype);
fprintf(fid,'max_iter %d\n',max_iter);
fprintf(fid,'num_samples 1\n');
fprintf(fid,'data_dim %d\n',chans);
fprintf(fid,'field_dim %d\n',frames);
fprintf(fid,'field_blocksize 1\n');
fprintf(fid,'do_history %d\n',do_history);
fprintf(fid,'histstep %d\n',histstep);
fprintf(fid,'lrate %f\n', lrate);
fprintf(fid,'minlrate %e\n', minlrate);
fprintf(fid,'mineig %e\n', mineig);
fprintf(fid,'lratefact %f\n', lratefact);
fprintf(fid,'rholrate %f\n', rholrate);
fprintf(fid,'rho0 %f\n', rho0);
fprintf(fid,'minrho %f\n', minrho);
fprintf(fid,'maxrho %f\n', maxrho);
fprintf(fid,'rholratefact %f\n',rholratefact);
fprintf(fid,'do_newton %d\n',do_newton);
fprintf(fid,'newt_start %d\n',newt_start);
fprintf(fid,'newt_ramp %d\n',newt_ramp);
fprintf(fid,'newtrate %f\n', newtrate);
fprintf(fid,'do_reject %d\n',do_reject);
fprintf(fid,'numrej %d\n',numrej);
fprintf(fid,'rejsig %f\n',rejsig);
fprintf(fid,'rejstart %d\n',rejstart);
fprintf(fid,'rejint %d\n',rejint);
fprintf(fid,'writestep %d\n',writestep);
fprintf(fid,'write_nd %d\n',write_nd);
fprintf(fid,'write_LLt %d\n',write_LLt);
fprintf(fid,'decwindow %d\n',decwindow);
fprintf(fid,'max_decs %d\n',maxdecs);
fprintf(fid,'fix_init %d\n',fix_init);
fprintf(fid,'update_A %d\n',update_A);
fprintf(fid,'update_c %d\n',update_c);
fprintf(fid,'update_gm %d\n',update_gamma);
fprintf(fid,'update_alpha %d\n',update_alpha);
fprintf(fid,'update_mu %d\n',update_mu);
fprintf(fid,'update_beta %d\n',update_sbeta);
fprintf(fid,'invsigmax %f\n',invsigmax);
fprintf(fid,'invsigmin %f\n',invsigmin);
fprintf(fid,'do_rho %d\n',do_rho);
fprintf(fid,'load_rej %d\n',load_rej);
fprintf(fid,'load_W %d\n',load_W);
fprintf(fid,'load_c %d\n',load_c);
fprintf(fid,'load_gm %d\n',load_gm);
fprintf(fid,'load_alpha %d\n',load_alpha);
fprintf(fid,'load_mu %d\n',load_mu);
fprintf(fid,'load_beta %d\n',load_sbeta);
fprintf(fid,'load_rho %d\n',load_rho);
fprintf(fid,'load_comp_list %d\n',load_comp_list);
if load_param == 1
    fprintf(fid,'indir %s\n',indir);
end
fprintf(fid,'do_mean %d\n',do_mean);
fprintf(fid,'do_sphere %d\n',do_sphere);
fprintf(fid,'doPCA %d\n',doPCA);
if pcakeep == -1
    pcakeep = chans;
end
fprintf(fid,'pcakeep %d\n',pcakeep);
fprintf(fid,'pcadb %f\n',pcadb);
fprintf(fid,'byte_size %d\n',bytesize);
fprintf(fid,'doscaling %d\n',doscaling);
fprintf(fid,'scalestep %d\n',scalestep);
fclose(fid);


%%%%%%%%%%%%%%%%%%%%% execute AMICA binary %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try
    if RUN_LOCALLY == 1 % single node
        
        disp('Running AMICA on single node ...');
        
        % execute the binary file with the generated parameter file
        tic;
        system([AMBIN ' ' outdir 'input.param']);
        t = toc;
        
        if t < 2, disp('Oops! AMICA did not run succesfully...'); end
        
    else % multiple nodes
        
        disp('Running AMICA across multiple nodes ...')
        
        if batch
            fid = fopen('tmp.sh','w');
            fprintf(fid,'#!/bin/sh\n');
            fprintf(fid,'#SBATCH --job-name="amica"\n');
            strtmp = '#SBATCH --output="amicatst.%j.%N.out"';
            fprintf(fid,'%s\n',strtmp);
            fprintf(fid,'#SBATCH --partition=compute\n');
            fprintf(fid,['#SBATCH --nodes=' int2str(numprocs) '\n']);
            fprintf(fid,'#SBATCH --ntasks-per-node=32\n');
            fprintf(fid,['#SBATCH --cpus-per-task=' int2str(max_threads) '\n']);
            fprintf(fid,'#SBATCH --mem=249208M\n');
            fprintf(fid,'#SBATCH --account=csd403\n');
            fprintf(fid,'#SBATCH --export=ALL\n');
            fprintf(fid,'#SBATCH -t 04:00:00\n');
            fprintf(fid,'module load cpu/0.15.4 slurm intel intel-mkl mvapich2\n');
            fprintf(fid,['export OMP_NUM_THREADS=' int2str(max_threads) ...
                '; export MV2_ENABLE_AFFINITY=0; export SRUN_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK}\n']);
            fprintf(fid,[' srun --export=ALL --mpi=pmi2 ' AMBIN ' ' outdir 'input.param\n']);
            fclose(fid);

            system('sbatch ./tmp.sh');
        else
            str = ['export OMP_NUM_THREADS=' int2str(max_threads) ' ' ...
                '; export MV2_ENABLE_AFFINITY=0 ' ...
                '; export SRUN_CPUS_PER_TASK=' int2str(max_threads) ' ' ...
                '; module purge; module load cpu/0.15.4 slurm intel intel-mkl mvapich2 ; ' ...
                'salloc --partition=compute --nodes=' int2str(numprocs) ' --ntasks-per-node=32 ' ...
                ' --cpus-per-task=' int2str(max_threads) ...
                ' --mem=249208M --account=csd403 -t 04:00:00 ' ...
                'srun --partition=compute --nodes=' int2str(numprocs) ' --ntasks-per-node=32 ' ...
                ' --cpus-per-task=' int2str(max_threads) ...
                ' --mem=249208M --account=csd403 --export=ALL -t 04:00:00 ' ...
                AMBIN ' ' outdir 'input.param'];

            system(str);

            if nargout > 0
                EEG = eeg_loadamica(dat,outdir,1);
                eeglab redraw;
            end
        end

        %system(['ibrun --npernode ' int2str(npernode) ' ' AMBIN ' ' outdir 'input.param']); % added by SHH
        %system(['ssh `qconf -ss` -n ' MPI_BIN 'mpirun -np ' int2str(numprocs) ' -machinefile ' basepath 'machines ' ...
        %             AMBIN ' ' outdir 'input.param']);

    end
    
    
catch    
    
    disp('Cannot execute AMICA binary.')
    
end

% remove the tmp data
%if exist(filename,'file')
%    delete(filename);
%end

