% pop_runamica() - Perform AMICA -- adaptive ICA using multiple models through GUI.

% Optional keywords:
%
%   indir               optional input directory from which to load init
%   num_models          number of models to learn, default = 1
%   num_mix_comps       number of mixture components in source model, def=3
%   numprocs            number or processors (slots) to use, def=8
%   max_iter            maximum number of iterations to perform, def=2000
%   lrate               initial learning rate for natural gradient, def=0.1
%   lratefact           multiplicative factor by which to decrease lrate, def=0.5
%   minlrate            lrate after which to stop, def=1e-8
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
%   numrej              for rejection, number of rejections to perform, def=3
%   rejsig              for rejection, number of standard dev of likelihood
%                           below which to reject data
%   rejstart            for rejection, iteration at which to start reject, def=3
%   rejint              for rejection, iteration interval between reject, def=3
%   max_threads         maximum number of threads to use on a node, def=4
%   decwindow           moving average window to detect likelihood decrease, def=1
%   update_A            flag to update mixing matrices, def=1
%   update_c            flag to update model centers, def=1
%   update_gamma        flag to update model probabilities, def=1
%   update_alpha        flag to update source mixture proportions, def=1
%   update_mu           flag to update source mixture mixture locations, def=1
%   update_sbeta        flag to update source mixture scales, def=1
%   invsigmax           maximum value of inverse scale parameters, def=100.0
%   invsigmin           minimum value of inverse scale parameters, def=1e-8
%   do_rho              flag to update shape parameters, def=1
%   load_rej            flag to load LLt to get rejections from, def=0
%   load_param          flag to load parameters, def=0
%   pcakeep             for PCA reduction, number of components to keep, def=chans
%   doscaling           flag to rescale unmixing matrix rows to unit norm, def=1
%   block_size          matrix block size (for block multiplication), def=128
%
% Disabled:
%   do_mean             flag to remove mean from data, def=1
%   do_sphere           flag to sphere data before ica, def=1
%   doPCA               flag to to PCA dimensionalit reduction, def=0
%   kurt_int            for ext. infomax, iteration interval between calc, def=1
%   kurt_start          for ext. infomax, iter to start kurtosis calc, def=3
%   load_comp_list      flag to load component assignment list, def=0
%   num_kurt            for ext. infomax, number of kurtosis calc, def=5
%   scalestep           iteration interval at which to rescale unmixing rows, def=1
%


function EEG = pop_runamica(EEG)
%global EEG
global blocksize do_opt_block blk_min blk_step blk_max setblk;
global do_sphere mineig setsph;
global do_reject rejsig numrej rejstart rejint setrej;
global num_mix_comps setmod;
global maxiter use_min_dll min_dll use_grad_norm min_grad_norm setstop;
global lrate do_newton newt_start do_history histstep writestep setlrat;
global numprocs machinefile setmpi dompi

basepath = [fileparts(mfilename('fullpath')) filesep]; % directory where this .m file is located

defaultOutputDirectory = [EEG.filepath 'amicaout'];
defaultInFile = [EEG.filepath EEG.datfile];
setsph = 0;
setrej = 0;
setmod = 0;
setstop = 0;
setblk = 0;
setlrat = 0;
setmpi = 0;
dompi = 0;

spheredlg = 'global do_sphere mineig setsph;'; 
spheredlg = [spheredlg 'sphres = inputgui(''geometry'',{ [0.5 1.5], [0.5 0.5 1] },''uilist'',{'];
spheredlg = [spheredlg '{''style'' ''text'' ''string'' ''Pre-processing:''},'];
spheredlg = [spheredlg '{''style'' ''popupmenu'' ''string'' ''Symm. sphering (approx. if reduced rank)|Princ. Comps (Eigenvectors)|No sphering trans''},'];
spheredlg = [spheredlg '{''style'' ''text'' ''string'' ''min. eigv:''},'];
spheredlg = [spheredlg '{''style'' ''edit'' ''string'' ''1e-15''},{} });'];
spheredlg = [spheredlg 'do_sphere = sphres{1}; mineig = str2num(sphres{2}); setsph = 1'];

rejectdlg = 'global do_reject rejsig numrej rejstart rejint setrej;';
rejectdlg = [rejectdlg 'rejres = inputgui(''geometry'',{ 1 [1 1], [1 1], [1 1], [1 1] },''uilist'',{'];
rejectdlg = [rejectdlg '{''style'' ''checkbox'' ''string'' ''Reject Un-Likely Data Samples'' ''value'' 1},'];
rejectdlg = [rejectdlg '{''style'' ''text'' ''string'' ''Num Reject:''},'];
rejectdlg = [rejectdlg '{''style'' ''edit'' ''string'' ''5''},'];
rejectdlg = [rejectdlg '{''style'' ''text'' ''string'' ''Num LL stdev:''},'];
rejectdlg = [rejectdlg '{''style'' ''edit'' ''string'' ''3''},'];
rejectdlg = [rejectdlg '{''style'' ''text'' ''string'' ''Start iter:''},'];
rejectdlg = [rejectdlg '{''style'' ''edit'' ''string'' ''1''},'];
rejectdlg = [rejectdlg '{''style'' ''text'' ''string'' ''Iter btwn rej:''},'];
rejectdlg = [rejectdlg '{''style'' ''edit'' ''string'' ''0''} });'];
rejectdlg = [rejectdlg 'do_reject = rejres{1}; numrej = str2num(rejres{2}); rejsig = str2num(rejres{3});'];
rejectdlg = [rejectdlg 'rejstart = str2num(rejres{4}); rejint = str2num(rejres{5})+1; setrej = 1;'];

modeldlg = 'global num_mix_comps setmod;';
modeldlg = [modeldlg 'modres = inputgui(''geometry'',{ [1 1] },''uilist'',{'];
modeldlg = [modeldlg '{''style'' ''text'' ''string'' ''Num mix pdf:''},'];
modeldlg = [modeldlg '{''style'' ''edit'' ''string'' ''3''} });'];
modeldlg = [modeldlg 'num_mix_comps = str2num(modres{1}); setmod = 1;'];

stoppingdlg = 'global maxiter use_min_dll min_dll use_grad_norm min_grad_norm setstop;';
stoppingdlg = [stoppingdlg 'stopres = inputgui(''geometry'',{ [1 1], [1 1], [1 1] },''uilist'',{'];
stoppingdlg = [stoppingdlg '{''style'' ''text'' ''string'' ''Maximum iterartions:''},'];
stoppingdlg = [stoppingdlg '{''style'' ''edit'' ''string'' ''2000''},'];
stoppingdlg = [stoppingdlg '{''style'' ''checkbox'' ''string'' ''Min gradient norm:'' ''value'' 1},'];
stoppingdlg = [stoppingdlg '{''style'' ''edit'' ''string'' ''1e-7''},'];
stoppingdlg = [stoppingdlg '{''style'' ''checkbox'' ''string'' ''Min LL change:'' ''value'' 1},'];
stoppingdlg = [stoppingdlg '{''style'' ''edit'' ''string'' ''1e-9''} });'];
stoppingdlg = [stoppingdlg 'maxiter = str2num(stopres{1}); use_grad_norm = stopres{2}; min_grad_norm = str2num(stopres{3});'];
stoppingdlg = [stoppingdlg 'use_min_dll = stopres{4}; min_dll = str2num(stopres{5}); setstop = 1;'];

lratedlg = 'global lrate do_newton newt_start do_history histstep writestep setlrat;';
lratedlg = [lratedlg 'lrateres = inputgui(''geometry'',{ [1.5 1], [1.5 1], [1.5 1], [1.5 1]},''uilist'',{'];
lratedlg = [lratedlg '{''style'' ''text'' ''string'' ''Initial learning rate:''},'];
lratedlg = [lratedlg '{''style'' ''edit'' ''string'' ''0.1''},'];
lratedlg = [lratedlg '{''style'' ''checkbox'' ''string'' ''Start Newt. (iter):'' ''value'' 1},'];
lratedlg = [lratedlg '{''style'' ''edit'' ''string'' ''50''},'];
lratedlg = [lratedlg '{''style'' ''checkbox'' ''string'' ''Write out. (iter):'' ''value'' 1},'];
lratedlg = [lratedlg '{''style'' ''edit'' ''string'' ''50''},'];
lratedlg = [lratedlg '{''style'' ''checkbox'' ''string'' ''Save Hist. (iter):'' ''value'' 0},'];
lratedlg = [lratedlg '{''style'' ''edit'' ''string'' ''10''} });'];
lratedlg = [lratedlg 'lrate = str2num(lrateres{1}); do_newton = lrateres{2}; newt_start = str2num(lrateres{3});'];
lratedlg = [lratedlg 'if lrateres{4}>0,  writestep = str2num(lrateres{5}); else writestep = 0, end;'];
lratedlg = [lratedlg 'do_history = lrateres{6}; histstep = str2num(lrateres{7}); setlrat = 1;'];

blockdlg = 'global blocksize do_opt_block blk_min blk_step blk_max setblk;';
blockdlg = [blockdlg 'blockres = inputgui(''geometry'',{ [1 1], 1, [1 1], [1 1], [1 1]},''uilist'',{'];
blockdlg = [blockdlg '{''style'' ''text'' ''string'' ''Matrix Block Size:''},'];
blockdlg = [blockdlg '{''style'' ''edit'' ''string'' ''256''},'];
blockdlg = [blockdlg '{''style'' ''checkbox'' ''string'' ''Optimize Block Size:'' ''value'' 1},'];
blockdlg = [blockdlg '{''style'' ''text'' ''string'' ''Min block size''},'];
blockdlg = [blockdlg '{''style'' ''edit'' ''string'' ''256''},'];
blockdlg = [blockdlg '{''style'' ''text'' ''string'' ''Blk test step:''},'];
blockdlg = [blockdlg '{''style'' ''edit'' ''string'' ''256''},'];
blockdlg = [blockdlg '{''style'' ''text'' ''string'' ''Max block size:''},'];
blockdlg = [blockdlg '{''style'' ''edit'' ''string'' ''1024''} });'];
blockdlg = [blockdlg 'blocksize = str2num(blockres{1}); do_opt_block = blockres{2}; blk_min = str2num(blockres{3});'];
blockdlg = [blockdlg 'blk_step = str2num(blockres{4}); blk_max = str2num(blockres{5}); setblk = 1;'];

mpidlg = 'global numprocs machinefile setmpi dompi;';
mpidlg = [mpidlg 'mpires = inputgui(''geometry'',{ [1 3], [1 3], [1 1 2]},''uilist'',{'];
mpidlg = [mpidlg '{''style'' ''checkbox'' ''string'' ''Use MPI -- Location:'' ''value'' 0},'];
mpidlg = [mpidlg '{''style'' ''edit'' ''string'' ''$MPI_BIN''},'];
mpidlg = [mpidlg '{''style'' ''text'' ''string'' ''Host file:''},'];
mpidlg = [mpidlg '{''style'' ''edit'' ''string'' [fileparts(which(''pop_runamica'')) filesep ''machines'']},'];
mpidlg = [mpidlg '{''style'' ''text'' ''string'' ''Number of Nodes:''},'];
mpidlg = [mpidlg '{''style'' ''edit'' ''string'' ''4''} {}});'];
mpidlg = [mpidlg 'dompi = mpires{1}; if mpires{1} MPI_BIN = mpires{2}; RUN_LOCALLY = 0; else RUN_LOCALLY = 1, end; machinefile = mpires{3}, numprocs = str2num(mpires{4}),'];
mpidlg = [mpidlg 'fid = fopen(machinefile,''r''); clear mchs; k=1; while 1 str = fgetl(fid); if ~ischar(str) break; else mchs{k} = str; k=k+1;end; end; fclose(fid);'];
mpidlg = [mpidlg 'mchs{1:numprocs}; setmpi = 1;'];
%mpidlg = [mpidlg 'ld = getloads(mchs); [mk,ord]=sort(ld); mchs{1:numprocs}'];

uigeom = {...
    [0.2 0.8] ...
    [0.2 0.8] ...
    [0.2 0.2 0.2 0.2 0.2] ...
    [0.2 0.2 0.4 0.2] ...    
    [0.2 0.2 0.6] ...
    [0.2 0.2 0.2 0.2 0.2] };

uilist = {...
      {'style' 'text' 'string' 'In File:'} ...
      {'style' 'edit' 'string' defaultInFile} ...
    {'style' 'text' 'string' 'Out Directory:'} ...
    {'style' 'edit' 'string' defaultOutputDirectory} ...
      {'style' 'text' 'string' '# of Channels:'} ...
      {'style' 'edit' 'string' int2str(EEG.nbchan)}   ...
      {'style' 'text' 'string' '# PCA Dims.:'} ...
      {'style' 'edit' 'string' int2str(EEG.nbchan)}   ...
      {'style' 'pushbutton' 'string' 'Sphering ...' 'callback' spheredlg } ...
    {'style' 'text' 'string' '# of Models:'} ...
    {'style' 'edit' 'string' '1'} ...
    {'style' 'popupmenu' 'string' 'Generalized Gaussian|Gaussian|Extended Infomax'} ...    
    {'style' 'pushbutton' 'string' 'Models ...' 'callback' modeldlg} ...
      {'style' 'pushbutton' 'string' 'Reject ...' 'callback' rejectdlg} ...
      {'style' 'pushbutton' 'string' 'Stopping ...' 'callback' stoppingdlg} ...
      {} ...
    {'style' 'text' 'string' '# of Threads:'} ...
    {'style' 'edit' 'string' '2'} ...
    {'style' 'pushbutton' 'string' 'mpi ...' 'callback' mpidlg} ...
    {'style' 'pushbutton' 'string' 'Block Sizes ...' 'callback' blockdlg} ...
    {'style' 'pushbutton' 'string' 'Lrates, Etc. ...' 'callback' lratedlg} ...
};


guititle = 'AMICA -- pop_runamica()';
result = inputgui(uigeom, uilist, 'pophelp(''pop_runamica'')', guititle, [], 'normal');
datfile = result{1}; outdir = result{2}; nchan = str2num(result{3}); pcakeep = str2num(result{4});
num_mod = str2num(result{5}); pdftype = result{6}(:)-1; max_threads = str2num(result{7});

arglist = {'outdir',outdir,'num_chans',nchan,'pcakeep',pcakeep,'max_threads',max_threads};
if num_mod > 1
    arglist = {arglist{:},'num_models',num_mod};
end
if pdftype ~= 0
    arglist = {arglist{:},'pdftype',pdftype};
end
if setrej
    arglist = {arglist{:},'do_reject',do_reject,'rejsig',rejsig,'numrej',numrej,'rejstart',rejstart,'rejint',rejint};
end
if setsph
    if do_sphere == 1
        do_approx_sphere = 1;
    elseif do_sphere == 2
        do_approx_sphere = 0;
    end
    arglist = {arglist{:},'do_approx_sphere',do_approx_sphere,'mineig',mineig};
end
if setmod
    arglist = {arglist{:},'num_mix_comps',num_mix_comps};
end
if setstop
    arglist = {arglist{:},'max_iter',maxiter,'use_grad_norm',use_grad_norm,'min_grad_norm',min_grad_norm,...
        'use_min_dll',use_min_dll,'min_dll',min_dll};
end
if setlrat
    arglist = {arglist{:},'lrate',lrate,'do_newton',do_newton,'newt_start',newt_start,'writestep',writestep,...
        'do_history',do_history,'histstep',histstep};
end
if setblk
    arglist = {arglist{:},'block_size',blocksize,'do_opt_block',do_opt_block,'blk_min',blk_min,'blk_step',blk_step,'blk_max',blk_max};
end
if setmpi && dompi
    arglist = {arglist{:},'numprocs',numprocs}; %,'machinefile',machinefile};
end

if isfield(EEG,'datfile') && length(EEG.datfile) > 0
    disp('Found datfile');
    [W,S,mods] = runamica15(datfile,arglist{:});
else
    disp('No datfile field found in EEG structure. Will write temp file.');
    [W,S,mods] = runamica15(EEG.data(:,:),arglist{:});
end

EEG.icaweights = W;
EEG.icasphere = S(1:size(W,1),:);
EEG.icawinv = mods.A(:,:,1);
EEG.mods = mods;
