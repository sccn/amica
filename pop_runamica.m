% pop_runamica() - Perform AMICA -- adaptive ICA using multiple models 
%                  through GUI or command line call.
%
% Usage:
% EEG = pop_runamica(EEG) % pop up interface
% EEG = pop_runamica(EEG, 'key', val)
%
% Optional keywords:
%   'num_models'     - number of models to learn, default = 1              
%   'num_mix_comps'  - number of mixture components in source model, def=3 
%   'numprocs'       - number or processors (slots) to use, def=8 
%   'maxiter'        - maximum number of iterations to perform, def=2000   
%   'lrate'          - initial learning rate for natural gradient, def=0.1
%   'do_newton'      - flag for newton method, default = 1 (do newton)
%   'newt_start'     - for newton method, iter at which to start newton, def=50
%   'writestep'      - iteration interval between output writes, def=10
%   'do_reject'      - flag for doing rejection of time points, def=0
%   'numrej'         - for rejection, number of rejections to perform, def=3
%   'rejsig'         - for rejection, number of standard dev of likelihood
%                      below which to reject data {DEFAULT??}
%   'rejstart'       - for rejection, iteration at which to start reject, def=3
%   'rejint'         - for rejection, iteration interval between reject, def=3
%   'max_threads'    - maximum number of threads to use on a node, def=4
%   'pcakeep'        - for PCA reduction, number of components to keep, def=chans
%   'block_size'     - matrix block size (for block multiplication), def=128
%
% Authors: Jason Palmer, Arnaud Delorme and Ramon Martinez-Cancino, SCCN 2017

%   indir               optional input directory from which to load init (MISSING)
%   lratefact           multiplicative factor by which to decrease lrate, def=0.5 (MISSING) 
%   minlrate            lrate after which to stop, def=1e-8                       (MISSING)
%   rholrate            initial lrate for shape parameters, def=0.05              (MISSING)
%   rho0                initial shape parameter value, def=1.5                    (MISSING)
%   minrho              minimum shape parameter value, def=1.0                    (MISSING)
%   maxrho              maximum shape parameter value, def=2.0                    (MISSING)
%   rholratefact        multiplicative factor by which to dec rholrate, def=0.5   (MISSING)

%   decwindow           moving average window to detect likelihood decrease, def=1 (MISSING)
%   doscaling           flag to rescale unmixing matrix rows to unit norm, def=1  (MISSING)

%   newtrate            for newton method, lrate for newton iterations, def=1.0          (MISSING)
%   newt_ramp           for newton method, number of iter to ramp up to newtrate, def=10 (MISSING)

%   write_nd            flag to write history of component update norms, def=1    (MISSING)
%   write_llt           flag to write model log likelihoods of time points, def=1 (MISSING)

% Disabled:
%   do_mean             flag to remove mean from data, def=1
%   do_sphere           flag to sphere data before ica, def=1
%   doPCA               flag to to PCA dimensionalit reduction, def=0
%   kurt_int            for ext. infomax, iteration interval between calc, def=1
%   kurt_start          for ext. infomax, iter to start kurtosis calc, def=3
%   load_comp_list      flag to load component assignment list, def=0
%   num_kurt            for ext. infomax, number of kurtosis calc, def=5
%   scalestep           iteration interval at which to rescale unmixing rows, def=1

% (MISSING) ---------------------------------------------------------------
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
% -------------------------------------------------------------------------


function [EEG, com] = pop_runamica(EEG,varargin)

com = [];

% Options
machinefile = [fileparts(which('pop_runamica')) filesep 'machines'];
try
    options = varargin;
    if ~isempty( varargin )
        for i = 1:2:numel(options)
            g.(options{i}) = options{i+1};
        end
    else g= []; end
catch
    disp('pop_runamica() error: calling convention {''key'', value, ... } error'); return;
end
% NOte : detect if CAR and change pcakeep to nchan-1
try g.datfile;               catch, g.datfile       = fullfile(EEG.filepath, EEG.datfile);   end
try g.outdir;                catch, g.outdir        = fullfile(EEG.filepath, 'amicaout');    end
try g.num_mod;               catch, g.num_mod       = 1;                            end
try g.maxiter;               catch, g.maxiter       = 2000;                         end
try g.max_threads;           catch, g.max_threads   = 4;                            end
try g.pcakeep;               catch, g.pcakeep       = EEG.nbchan;                   end
try g.nchan;                 catch, g.nchan         = EEG.nbchan;                   end
try g.pdftype;               catch, g.pdftype       = 0;                            end

% Block of former global vars
% Flags --------------------------------------------------------
try g.setsph;               catch, g.setsph  = 0;           end
try g.setrej;               catch, g.setrej  = 0;           end
try g.setmod;               catch, g.setmod  = 0;           end
try g.setstop;              catch, g.setstop = 0;           end
try g.setblk;               catch, g.setblk  = 0;           end
try g.setlrat;              catch, g.setlrat = 0;           end
try g.setmpi;               catch, g.setmpi  = 0;           end
%---------------------------------------------------------------
 
try g.do_sphere;        g.setsph = 1;   catch, g.do_sphere     = 1;             end
try g.mineig;           g.setsph = 1;   catch, g.mineig        = 1e-12;         end 

try g.do_reject;        g.setrej = 1;   catch, g.do_reject     = 0;             end 
try g.numrej;           g.setrej = 1;   catch, g.numrej        = 3;             end
try g.rejsig;           g.setrej = 1;   catch, g.rejsig        = 3;             end 
try g.rejstart;         g.setrej = 1;   catch, g.rejstart      = 3;             end 
try g.rejint;           g.setrej = 1;   catch, g.rejint        = 3;             end 

try g.num_mix_comps;    g.setmod = 1;   catch, g.num_mix_comps = 3;             end

try g.maxiter;          g.setstop = 1;  catch, g.maxiter       = 2000;          end
try g.use_grad_norm;    g.setstop = 1;  catch, g.use_grad_norm = [];            end
try g.min_grad_norm;    g.setstop = 1;  catch, g.min_grad_norm = [];            end
try g.use_min_dll;      g.setstop = 1;  catch, g.use_min_dll   = 1;             end
try g.min_dll;          g.setstop = 1;  catch, g.min_dll       = 1e-9;          end

try g.lrate;            g.setlrat = 1;  catch, g.lrate         = 0.1;           end
try g.do_newton;        g.setlrat = 1;  catch, g.do_newton     = 1;             end
try g.newt_start;       g.setlrat = 1;  catch, g.newt_start    = 50;            end
try g.writestep;        g.setlrat = 1;  catch, g.writestep     = 10;            end
try g.do_history;       g.setlrat = 1;  catch, g.do_history    = 0;             end
try g.histstep;         g.setlrat = 1;  catch, g.histstep      = 10;            end

try g.block_size;       g.setblk = 1;   catch, g.block_size    = 128;           end
try g.do_opt_block;     g.setblk = 1;   catch, g.do_opt_block  = 1;             end
try g.blk_min;          g.setblk = 1;   catch, g.blk_min       = 256;           end
try g.blk_step;         g.setblk = 1;   catch, g.blk_step      = 256;           end
try g.blk_max;          g.setblk = 1;   catch, g.blk_max       = 1024;          end

try g.dompi;            g.setmpi = 1;   catch, g.dompi         = 0;             end
try g.numprocs;         g.setmpi = 1;   catch, g.numprocs      = 1;             end
try g.machinefile;      g.setmpi = 1;   catch, g.machinefile   = machinefile;   end

if nargin == 1
    % Formers globals initialized for GUI
    g.setsph  = 0;
    g.setrej  = 0;
    g.setmod  = 0;
    g.setstop = 0;
    g.setblk  = 0;
    g.setlrat = 0;
    g.setmpi  = 0;
    g.dompi   = 0;
    
    % ---
    spheredlg = ['sphres = inputgui(''geometry'',{ [0.5 1.5], [0.5 0.5 1] },''uilist'',{'];
    spheredlg = [spheredlg '{''style'' ''text'' ''string'' ''Pre-processing:''},'];
    spheredlg = [spheredlg '{''style'' ''popupmenu'' ''string'' ''Symm. sphering (approx. if reduced rank)|Princ. Comps (Eigenvectors)|No sphering trans''},'];
    spheredlg = [spheredlg '{''style'' ''text'' ''string'' ''min. eigv:''},'];
    spheredlg = [spheredlg '{''style'' ''edit'' ''string'' ''1e-15''},{} });'];
    spheredlg = [spheredlg 'if ~isempty(sphres)'];
    spheredlg = [spheredlg 'g = get(get(gcbo,''parent''),''UserData'');'];
    spheredlg = [spheredlg 'g.do_sphere = sphres{1}; g.mineig = str2num(sphres{2}); g.setsph = 1;'];
    spheredlg = [spheredlg 'set(get(gcbo,''parent''),''UserData'', g);'];
    spheredlg = [spheredlg 'end;'];
    
    rejectdlg = ['rejres = inputgui(''geometry'',{ 1 [1 1], [1 1], [1 1], [1 1] },''uilist'',{'];
    rejectdlg = [rejectdlg '{''style'' ''checkbox'' ''string'' ''Reject Un-Likely Data Samples'' ''value'' 1},'];
    rejectdlg = [rejectdlg '{''style'' ''text'' ''string'' ''Num Reject:''},'];
    rejectdlg = [rejectdlg '{''style'' ''edit'' ''string'' ''5''},'];
    rejectdlg = [rejectdlg '{''style'' ''text'' ''string'' ''Num LL stdev:''},'];
    rejectdlg = [rejectdlg '{''style'' ''edit'' ''string'' ''3''},'];
    rejectdlg = [rejectdlg '{''style'' ''text'' ''string'' ''Start iter:''},'];
    rejectdlg = [rejectdlg '{''style'' ''edit'' ''string'' ''1''},'];
    rejectdlg = [rejectdlg '{''style'' ''text'' ''string'' ''Iter btwn rej:''},'];
    rejectdlg = [rejectdlg '{''style'' ''edit'' ''string'' ''0''} });'];
    rejectdlg = [rejectdlg 'if ~isempty(rejres)'];
    rejectdlg = [rejectdlg 'g = get(get(gcbo,''parent''),''UserData'');'];
    rejectdlg = [rejectdlg 'g.do_reject = rejres{1}; g.numrej = str2num(rejres{2}); g.rejsig = str2num(rejres{3});'];
    rejectdlg = [rejectdlg 'g.rejstart = str2num(rejres{4}); g.rejint = str2num(rejres{5})+1; g.setrej = 1;'];
    rejectdlg = [rejectdlg 'set(get(gcbo,''parent''),''UserData'', g);'];
    rejectdlg = [rejectdlg 'end;'];
    
    modeldlg = ['modres = inputgui(''geometry'',{ [1 1] },''uilist'',{'];
    modeldlg = [modeldlg '{''style'' ''text'' ''string'' ''Num mix pdf:''},'];
    modeldlg = [modeldlg '{''style'' ''edit'' ''string'' ''3''} });'];
    modeldlg = [modeldlg 'if ~isempty(modres)'];
    modeldlg = [modeldlg 'g = get(get(gcbo,''parent''),''UserData'');'];
    modeldlg = [modeldlg 'g.num_mix_comps = str2num(modres{1}); g.setmod = 1;'];
    modeldlg = [modeldlg 'set(get(gcbo,''parent''),''UserData'', g);'];
    modeldlg = [modeldlg 'end;'];
    
    stoppingdlg = ['stopres = inputgui(''geometry'',{ [1 1], [1 1], [1 1] },''uilist'',{'];
    stoppingdlg = [stoppingdlg '{''style'' ''text'' ''string'' ''Maximum iterations:''},'];
    stoppingdlg = [stoppingdlg '{''style'' ''edit'' ''string'' ''2000''},'];
    stoppingdlg = [stoppingdlg '{''style'' ''checkbox'' ''string'' ''Min gradient norm:'' ''value'' 1},'];
    stoppingdlg = [stoppingdlg '{''style'' ''edit'' ''string'' ''1e-7''},'];
    stoppingdlg = [stoppingdlg '{''style'' ''checkbox'' ''string'' ''Min LL change:'' ''value'' 1},'];
    stoppingdlg = [stoppingdlg '{''style'' ''edit'' ''string'' ''1e-9''} });'];
    stoppingdlg = [stoppingdlg 'if ~isempty(stopres)'];
    stoppingdlg = [stoppingdlg 'g = get(get(gcbo,''parent''),''UserData'');'];
    stoppingdlg = [stoppingdlg 'g.maxiter = str2num(stopres{1}); g.use_grad_norm = stopres{2}; g.min_grad_norm = str2num(stopres{3});'];
    stoppingdlg = [stoppingdlg 'g.use_min_dll = stopres{4}; g.min_dll = str2num(stopres{5}); g.setstop = 1;'];
    stoppingdlg = [stoppingdlg 'set(get(gcbo,''parent''),''UserData'', g);'];
    stoppingdlg = [stoppingdlg 'end;'];
    
    lratedlg = ['lrateres = inputgui(''geometry'',{ [1.5 1], [1.5 1], [1.5 1], [1.5 1]},''uilist'',{'];
    lratedlg = [lratedlg '{''style'' ''text'' ''string'' ''Initial learning rate:''},'];
    lratedlg = [lratedlg '{''style'' ''edit'' ''string'' ''0.1''},'];
    lratedlg = [lratedlg '{''style'' ''checkbox'' ''string'' ''Start Newt. (iter):'' ''value'' 1},'];
    lratedlg = [lratedlg '{''style'' ''edit'' ''string'' ''50''},'];
    lratedlg = [lratedlg '{''style'' ''checkbox'' ''string'' ''Write out. (iter):'' ''value'' 1},'];
    lratedlg = [lratedlg '{''style'' ''edit'' ''string'' ''50''},'];
    lratedlg = [lratedlg '{''style'' ''checkbox'' ''string'' ''Save Hist. (iter):'' ''value'' 0},'];
    lratedlg = [lratedlg '{''style'' ''edit'' ''string'' ''10''} });'];
    lratedlg = [lratedlg 'if ~isempty(lrateres)'];
    lratedlg = [lratedlg 'g = get(get(gcbo,''parent''),''UserData'');'];
    lratedlg = [lratedlg 'g.lrate = str2num(lrateres{1}); g.do_newton = lrateres{2}; g.newt_start = str2num(lrateres{3});'];
    lratedlg = [lratedlg 'if lrateres{4}>0,  g.writestep = str2num(lrateres{5}); else g.writestep = 0, end;'];
    lratedlg = [lratedlg 'g.do_history = lrateres{6}; g.histstep = str2num(lrateres{7}); g.setlrat = 1;'];
    lratedlg = [lratedlg 'set(get(gcbo,''parent''),''UserData'', g);'];
    lratedlg = [lratedlg 'end;'];
    
    blockdlg = ['blockres = inputgui(''geometry'',{ [1 1], 1, [1 1], [1 1], [1 1]},''uilist'',{'];
    blockdlg = [blockdlg '{''style'' ''text'' ''string'' ''Matrix Block Size:''},'];
    blockdlg = [blockdlg '{''style'' ''edit'' ''string'' ''256''},'];
    blockdlg = [blockdlg '{''style'' ''checkbox'' ''string'' ''Optimize Block Size:'' ''value'' 1},'];
    blockdlg = [blockdlg '{''style'' ''text'' ''string'' ''Min block size''},'];
    blockdlg = [blockdlg '{''style'' ''edit'' ''string'' ''256''},'];
    blockdlg = [blockdlg '{''style'' ''text'' ''string'' ''Blk test step:''},'];
    blockdlg = [blockdlg '{''style'' ''edit'' ''string'' ''256''},'];
    blockdlg = [blockdlg '{''style'' ''text'' ''string'' ''Max block size:''},'];
    blockdlg = [blockdlg '{''style'' ''edit'' ''string'' ''1024''} });'];
    blockdlg = [blockdlg 'if ~isempty(blockres)'];
    blockdlg = [blockdlg 'g = get(get(gcbo,''parent''),''UserData'');'];
    blockdlg = [blockdlg 'g.block_size = str2num(blockres{1}); g.do_opt_block = blockres{2}; g.blk_min = str2num(blockres{3});'];
    blockdlg = [blockdlg 'g.blk_step = str2num(blockres{4}); g.blk_max = str2num(blockres{5}); g.setblk = 1;'];
    blockdlg = [blockdlg 'set(get(gcbo,''parent''),''UserData'', g);'];
    blockdlg = [blockdlg 'end;'];
    
    mpidlg = ['mpires = inputgui(''geometry'',{ [1 3], [1 3], [1 1 2]},''uilist'',{'];
    mpidlg = [mpidlg '{''style'' ''checkbox'' ''string'' ''Use MPI -- Location:'' ''value'' 0},'];
    mpidlg = [mpidlg '{''style'' ''edit'' ''string'' ''$MPI_BIN''},'];
    mpidlg = [mpidlg '{''style'' ''text'' ''string'' ''Host file:''},'];
    mpidlg = [mpidlg '{''style'' ''edit'' ''string'' [fileparts(which(''pop_runamica'')) filesep ''machines'']},'];
    mpidlg = [mpidlg '{''style'' ''text'' ''string'' ''Number of Nodes:''},'];
    mpidlg = [mpidlg '{''style'' ''edit'' ''string'' ''4''} {}});'];
    mpidlg = [mpidlg 'if ~isempty(mpires)'];
    mpidlg = [mpidlg 'g = get(get(gcbo,''parent''),''UserData'');'];
    mpidlg = [mpidlg 'g.dompi = mpires{1}; if mpires{1} MPI_BIN = mpires{2}; RUN_LOCALLY = 0; else RUN_LOCALLY = 1, end; g.machinefile = mpires{3}, g.numprocs = str2num(mpires{4}),'];
    mpidlg = [mpidlg 'fid = fopen(g.machinefile,''r''); clear mchs; k=1; while 1 str = fgetl(fid); if ~ischar(str) break; else mchs{k} = str; k=k+1;end; end; fclose(fid);'];
    mpidlg = [mpidlg 'mchs{1:numprocs}; g.setmpi = 1;'];
    mpidlg = [mpidlg 'set(get(gcbo,''parent''),''UserData'', g);'];
    mpidlg = [mpidlg 'end;'];
    %mpidlg = [mpidlg 'ld = getloads(mchs); [mk,ord]=sort(ld); mchs{1:numprocs}'];
    
    uigeom = {...
        [0.2 0.8] ...
        [0.2 0.8] ...
        [0.2 0.2 0.2 0.2 0.2] ...
        [0.2 0.2 0.4 0.2] ...
        [0.2 0.2 0.6] ...
        [0.2 0.2 0.2 0.2 0.2] };
    
    uilist = {...
        {'style' 'text'        'string' 'In File:'} ...
        {'style' 'edit'        'string' g.datfile} ...
        {'style' 'text'        'string' 'Out Directory:'} ...
        {'style' 'edit'        'string' g.outdir} ...
        {'style' 'text'        'string' '# of Channels:'} ...
        {'style' 'edit'        'string' int2str(g.nchan)}   ...
        {'style' 'text'        'string' '# PCA Dims.:'} ...
        {'style' 'edit'        'string' int2str(g.pcakeep)}   ...
        {'style' 'pushbutton'  'string' 'Sphering ...'         'callback' spheredlg } ...
        {'style' 'text'        'string' '# of Models:'} ...
        {'style' 'edit'        'string' int2str(g.num_mod)} ...
        {'style' 'popupmenu'   'string' 'Generalized Gaussian|Gaussian|Extended Infomax'} ...
        {'style' 'pushbutton'  'string' 'Models ...'           'callback' modeldlg} ...
        {'style' 'pushbutton'  'string' 'Reject ...'           'callback' rejectdlg} ...
        {'style' 'pushbutton'  'string' 'Stopping ...'         'callback' stoppingdlg} ...
        {} ...
        {'style' 'text'        'string' '# of Threads:'} ...
        {'style' 'edit'        'string' '2'} ...
        {'style' 'pushbutton'  'string' 'mpi ...'              'callback' mpidlg} ...
        {'style' 'pushbutton'  'string' 'Block Sizes ...'      'callback' blockdlg} ...
        {'style' 'pushbutton'  'string' 'Lrates, Etc. ...'     'callback' lratedlg} ...
        };
    
    
    guititle = 'AMICA -- pop_runamica()';
    [result, g] = inputgui(uigeom, uilist, 'pophelp(''pop_runamica'')', guititle, g, 'normal');
    if ~isempty(result)
        g.datfile = result{1}; g.outdir = result{2}; g.nchan = str2num(result{3}); g.pcakeep = str2num(result{4});
        g.num_mod = str2num(result{5}); g.pdftype = result{6}(:)-1; g.max_threads = str2num(result{7});
    else
        return;
    end
end

% Either call from GUI or command line
arglist = {'outdir',g.outdir,'num_chans',g.nchan,'pcakeep',g.pcakeep,'max_threads',g.max_threads};
if g.num_mod > 1
    arglist = {arglist{:},'num_models',g.num_mod};
end
if g.pdftype ~= 0
    arglist = {arglist{:},'pdftype',g.pdftype};
end
if g.setrej
    arglist = {arglist{:},'do_reject',g.do_reject,'rejsig',g.rejsig,'numrej',g.numrej,'rejstart',g.rejstart,'rejint',g.rejint};
end
if g.setsph
    if g.do_sphere == 1
        do_approx_sphere = 1;
    elseif g.do_sphere == 2
        do_approx_sphere = 0;
    end
    arglist = {arglist{:},'do_approx_sphere',do_approx_sphere,'mineig',g.mineig};
end
if g.setmod
    arglist = {arglist{:},'num_mix_comps',g.num_mix_comps};
end
if g.setstop
    arglist = {arglist{:},'max_iter',g.maxiter,'use_grad_norm',g.use_grad_norm,'min_grad_norm',g.min_grad_norm,...
        'use_min_dll',g.use_min_dll,'min_dll',g.min_dll};
end
if g.setlrat
    arglist = {arglist{:},'lrate',g.lrate,'do_newton',g.do_newton,'newt_start',g.newt_start,'writestep',g.writestep,...
        'do_history',g.do_history,'histstep',g.histstep};
end
if g.setblk
    arglist = {arglist{:},'block_size',g.block_size,'do_opt_block',g.do_opt_block,'blk_min',g.blk_min,'blk_step',g.blk_step,'blk_max',g.blk_max};
end
if all([g.setmpi, g.dompi]) || ~isempty(g.numprocs)
    arglist = {arglist{:},'numprocs',g.numprocs}; %,'machinefile',machinefile};
end

if isfield(EEG,'datfile') && length(EEG.datfile) > 0
    disp('Found datfile');
    [W,S,mods] = runamica15(g.datfile,arglist{:});
else
    disp('No datfile field found in EEG structure. Will write temp file.');
    [W,S,mods] = runamica15(EEG.data(:,:),arglist{:});
end

EEG.icaweights = W;
EEG.icasphere  = S(1:size(W,1),:);
EEG.icawinv    = mods.A(:,:,1);
EEG.mods       = mods;

% Building com for history-------------------------------------------------
if nargout > 1
    varargstring = '';
    gtmp = g;
    
    if ~g.setsph
        gtmp = rmfield(gtmp,'do_sphere');
        gtmp = rmfield(gtmp,'mineig');
        gtmp = rmfield(gtmp,'setsph');
    end
    
    if ~g.setrej
        gtmp = rmfield(gtmp,'do_reject');
        gtmp = rmfield(gtmp,'numrej');
        gtmp = rmfield(gtmp,'rejsig');
        gtmp = rmfield(gtmp,'rejstart');
        gtmp = rmfield(gtmp,'rejint');
        gtmp = rmfield(gtmp,'setrej');
    end
    
    if ~g.setmod
        gtmp = rmfield(gtmp,'num_mix_comps');
        gtmp = rmfield(gtmp,'setmod');
    end
    
    if ~g.setstop
        gtmp = rmfield(gtmp,'maxiter');
        gtmp = rmfield(gtmp,'use_grad_norm');
        gtmp = rmfield(gtmp,'min_grad_norm');
        gtmp = rmfield(gtmp,'use_min_dll');
        gtmp = rmfield(gtmp,'min_dll');
        gtmp = rmfield(gtmp,'setstop');
    end
    
    if ~g.setlrat
        gtmp = rmfield(gtmp,'lrate');
        gtmp = rmfield(gtmp,'do_newton');
        gtmp = rmfield(gtmp,'newt_start');
        gtmp = rmfield(gtmp,'writestep');
        gtmp = rmfield(gtmp,'do_history');
        gtmp = rmfield(gtmp,'histstep');
        gtmp = rmfield(gtmp,'setlrat');
    end
    
    if ~g.setblk
        gtmp = rmfield(gtmp,'block_size');
        gtmp = rmfield(gtmp,'do_opt_block');
        gtmp = rmfield(gtmp,'blk_min');
        gtmp = rmfield(gtmp,'blk_step');
        gtmp = rmfield(gtmp,'blk_max');
        gtmp = rmfield(gtmp,'setblk');
    end
    
    if ~g.setmpi
        gtmp = rmfield(gtmp,'dompi');
        gtmp = rmfield(gtmp,'numprocs');
        gtmp = rmfield(gtmp,'machinefile');
        gtmp = rmfield(gtmp,'setmpi');
    end
    
    gnames =  fieldnames(gtmp); % getting names
    
    for i = 1:length(gnames)
        gval = gtmp.(gnames{i});
        if ~strcmp(gval,'do_approx_sphere')
            if isempty(gval)
                gval = '[]';
            elseif isstr(gval)
                gval= ['''' gval ''''];
            elseif isnumeric(gval)
                gval = num2str(gval);
            end
            varargstring = [varargstring ',''' gnames{i} ''',' gval];
        end
    end
    
    com = ['[EEG,com] = pop_runamica(EEG' varargstring ')'];
end
