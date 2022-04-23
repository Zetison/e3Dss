waringMessageId = 'MATLAB:mir_warning_maybe_uninitialized_temporary';
warning('off',waringMessageId)

if ~exist('noCoresToUse','var')
    noCoresToUse = feature('numCores');
end
poolobj = gcp('nocreate');
if isempty(poolobj)
    noCoresToUse = min(feature('numCores'),noCoresToUse);
    parpool(noCoresToUse, 'IdleTimeout', Inf)
elseif poolobj.NumWorkers ~= min(feature('numCores'),noCoresToUse)
    delete(poolobj)
    noCoresToUse = min(feature('numCores'),noCoresToUse);
    parpool(noCoresToUse, 'IdleTimeout', Inf)
end
maxNumCompThreads(noCoresToUse);