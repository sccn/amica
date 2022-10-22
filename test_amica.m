eeglab;
EEG = pop_loadset('Memorize.set');
EEG = pop_runamica(EEG, 'maxiter', 10, 'max_threads', 4);