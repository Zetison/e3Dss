function x_extremas = findExtremas(f, a, b, N,options)
newEpsilon = eps;
x = linspace(a, b, N);
if nargin < 5
    options = optimset('TolX',newEpsilon, 'TolFun', newEpsilon, 'MaxIter', 1000,'MaxFunEvals',1000); %,'Display','iter');
end
tic
F = f(x);
fprintf('Completed initial search in %f seconds.\n', toc)

candidates = or(F(2:end-1) > max(F(1:end-2),F(3:end)), F(2:end-1) < min(F(1:end-2),F(3:end)));
candidates = find([0 candidates 0]);

II = [x(candidates-1); x(candidates+1)];
x_extremas2 = zeros(size(II));
% for i = 1:size(II,2)
parfor i = 1:size(II,2)
    tic
    I = II(:,i);
    x0 = mean(I);
    x_min = fminsearchbnd(f, x0, I(1), I(2), options);
    x_max = fminsearchbnd(@(x)-f(x), x0, I(1), I(2), options);
    x_res = zeros(2,1);
    if min(abs(x_min-I)) > 1e6*options.TolX
        x_res(1) = x_min;
    end
    if min(abs(x_max-I)) > 1e6*options.TolX
        x_res(2) = x_max;
    end
    x_extremas2(:,i) = x_res;
    fprintf('Completed %d out of %d candidate intervals. Ellapsed time = %f\n', i, size(II,2), toc);
end
x_extremas = [x_extremas2(1,:) x_extremas2(2,:)];
x_extremas(x_extremas == 0) = [];