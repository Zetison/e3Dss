
noCoresToUse = 321;
waringMessageId = 'MATLAB:mir_warning_maybe_uninitialized_temporary';
warning('off',waringMessageId)
if verLessThan('matlab', '8.3 (R2014a)')
    isOpen = matlabpool('size') > 0;
    if ~isOpen
        matlabpool open 12
        warning('Change the above number to match the number of workers available!')
    end
else
    noCoresToUse = feature('numCores');
    if noCoresToUse > 12
        noCoresToUse = 12;
    end
    poolobj = gcp('nocreate');
    if isempty(poolobj)
        parpool(noCoresToUse, 'IdleTimeout', Inf)
    end
end
varCol.runInParallell = noCoresToUse; % If set to 0, parfor = for. Else, the number specify the number of workers to use. 