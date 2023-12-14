function showprogress(j)
%
if mod(j,10)==0
    fprintf('*')
else
    fprintf('.')
end
if mod(j,50)==0
    fprintf('\n')
end