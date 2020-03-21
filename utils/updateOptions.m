function options = updateOptions(options,newOptions)
% read the acceptable names

if isstruct(newOptions)
    newOptionNames = fieldnames(newOptions);   
    for i = 1:numel(newOptionNames) % pair is {propName;propValue}
        inpName = newOptionNames{i}; % make case insensitive

        options.(inpName) = newOptions.(inpName);
    end
else
    if isstruct(newOptions{1})
        if numel(options) < numel(newOptions)
            options(numel(options)+1:numel(newOptions)) = {struct()};
        end
        for i = 1:numel(newOptions)
            options{i} = updateOptions(options{i},newOptions{i});
        end
    else
        optionNames = fieldnames(options);
        % count arguments
        nArgs = length(newOptions);
        if round(nArgs/2) ~= nArgs/2
            error('Must have propertyName/propertyValue pairs')
        end

        for pair = reshape(newOptions,2,[]) % pair is {propName;propValue}
            inpName = pair{1}; % make case insensitive

            if any(strcmp(inpName,optionNames))
                options.(inpName) = pair{2};
            else
                error('%s is not a recognized parameter name',inpName)
            end
        end
    end
end