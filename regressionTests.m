%% Regression tests
% A successful run of this script should be performed prior to any push to git
close all
clear all
if 1
    startup
    Eps = 1e-5;
    testFolder = 'examples/tests/';
    studyName = {'Chang1994soa_Fig1617','Fender1972sfa_Fig23','Ayres1987ars_Fig1','Ihlenburg1998fea_Fig52','Skelton1997tao_Fig10567','Hetmaniuk2012raa_Fig81217','Sage1979mri_Fig14','Venas2019e3s_Fig161819'};
    studyName = {'Venas2019e3s_Fig161819'};
    stringShift = 60;
    noFailedTests = 0;
    for i_study = 1:numel(studyName)
        fprintf(['\n%-' num2str(stringShift) 's'], ['Running test ''' studyName{i_study} ''' ...'])
        testFailed = false;
        try
            eval(['tasks = ' studyName{i_study} ';']);
            tasks_ref = load([testFolder studyName{i_study} '.mat'],'tasks');
            for i_task = 1:numel(tasks)
                if iscell(tasks)
                    task = tasks{i_task};
                else
                    task = tasks(i_task);
                end
                if iscell(tasks_ref.tasks)
                    task_ref = tasks_ref.tasks{i_task};
                else
                    task_ref = tasks_ref.tasks(i_task);
                end
                
                fieldNames = fieldnames(task_ref);
                for field = fieldNames.'
                    entry = task.(field{1});
                    entry_ref = task_ref.(field{1});
                    
                    if any(isnan(entry_ref(:)))
                        continue
                    end
                    reg_error = norm(entry(:)-entry_ref(:))/norm(entry_ref(:));
                    if reg_error > Eps && ~strcmp(field{1},'visElements') % the delaunay routine has been changed in R2022 compared to when the ref was made and so visElements has also changed
                        testFailed = true;
                        break
                    end
                end
                if testFailed
                    break
                end
            end
            if testFailed
                fprintf('test failed due to incorrect results (i_study=%d, i_task=%d)! (relative error = %g)', i_study, i_task, reg_error)
                noFailedTests = noFailedTests + 1;
            else
                fprintf('successfully!')
            end
        catch ME
            rethrow(ME)
            fprintf('Test failed due to runtime error!')
            noFailedTests = noFailedTests + 1;
        end
    end
    fprintf(['\n\nNumber of failed tests: ' num2str(noFailedTests) '\n'])
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create and store tests
if 0
    testFolder = 'examples/tests/';

    studyName = {'Chang1994soa_Fig1617','Fender1972sfa_Fig23','Ayres1987ars_Fig1','Ihlenburg1998fea_Fig52','Skelton1997tao_Fig10567','Hetmaniuk2012raa_Fig81217','Sage1979mri_Fig14','Venas2019e3s_Fig161819'};
    studyName = {'Venas2019e3s_Fig161819'};

    for i_study = 1:numel(studyName)
        close all
        eval(['tasks = ' studyName{i_study} '(true);']);
        save([testFolder studyName{i_study} '.mat'],'tasks')
    %     keyboard
    end
end