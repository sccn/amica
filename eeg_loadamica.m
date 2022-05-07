function EEG = eeg_loadamica( EEG, outdir, modnum )
if nargin < 3
    modnum = 1;
end
mods = loadmodout15(outdir);

EEG.icasphere = mods.S(1:mods.num_pcs,:);
EEG.icaweights = mods.W(:,:,modnum);
EEG.icawinv = mods.A(:,:,modnum);
EEG.icachansind = 1:EEG.nbchan;
EEG = eeg_checkset(EEG);
%EEG.icaact = EEG.icaweights*EEG.icasphere*EEG.data;
%eeglab redraw
%pop_saveset(EEG);
%EEG = pop_loadset(EEG.filename,EEG.filepath);