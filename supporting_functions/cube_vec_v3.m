function res=cube_vec_v3(X, CUBE_RES, collision_alt_limit)
%Only consider RSO below collision_alt_limit for collision

idx_invalid = any(abs(X)'>collision_alt_limit);
X(idx_invalid,:) = NaN;
X_dis = floor(X(:,1:3)/CUBE_RES);

% Shift origin such that X_dis is always positive
shift_lim = max(abs(X_dis(:)))+10;
X_dis = X_dis+shift_lim;
shift_lim2 = 2*shift_lim;
X_idx = X_dis(:,1)*(shift_lim2*shift_lim2)+X_dis(:,2)*shift_lim2+X_dis(:,3);

[~, uniqueIdx] = unique(X_idx);
dupelicateIdx = ismember( X_idx, X_idx( setdiff( 1:numel(X_idx), uniqueIdx ) ) );
duplicates = X_idx(dupelicateIdx);
[~,~,ic] = unique(duplicates);
duplicate_idx = find(dupelicateIdx);

% func = @(x) {cat(1,x)};
func = @(x) {nchoosek(x,2)};
if ~isempty(ic)
    res = splitapply(func,duplicate_idx,ic);
else
    res = {};
end

% %% Backup (method)
% duplicate_group = findgroups(duplicates);
% X_unique = unique(X_idx);
% X_unique = [X_unique;X_unique(end)+1];
% Y = discretize(X_idx,X_unique);

% %% Old method v2
% % [counts, edges] = histcounts(X_idx,unique(X_idx));
% X_unique = unique(X_idx);
% X_unique = [X_unique;X_unique(end)+1];
% B = find(histcounts(X_idx,X_unique)>1);
% 
% res = cell(length(B),1);
% for ik=1:length(B)
%     res{ik} = find(X_idx==X_unique(B(ik)));
% end
% 
% %% Check if output matches
% % Reformat new method to old format
% res2 = cell(max(duplicate_group),1);
% for ik=1:max(duplicate_group)
%     res2{ik} = duplicate_idx(duplicate_group==ik);
% end
% 
% r_valid = zeros(length(res),1);
% for ik=1:length(res)
%     r_valid(ik) = any(cellfun(@(m)isequal(m,res{ik}),res2));
% end
% 
% if ~all(r_valid)
%     pause
% end
% 
% if length(res)~=length(res2)
%     pause
% end


% %% Old method v1
% M = containers.Map('KeyType','char','ValueType','any');
% 
% for k=1:length(X(1,:))
%     out=(floor(X(:,k)/CUBE_RES));
%     key=[num2str(out(1)) ',' num2str(out(2)) ',' num2str(out(3))];
%     if isKey(M,key)
%         M(key)=[M(key) k];
%     else 
%         M(key)=k;
%     end
% end
% K=keys(M);
% res2=cell(1);
% for k=1:M.Count
%     if length(M(K{1,k}))>1
%         res2{end+1} = M(K{1,k});
%     end
% end
% % 
% %% Check if output matches for old v1 and old v2
% r_valid = zeros(length(res),1);
% for ik=1:length(res)
%     r_valid(ik) = any(cellfun(@(m)isequal(m,res{ik}),res2));
% end
% 
% if ~all(r_valid)
%     pause
% end
% 
% if length(res)~=length(res2)-1
%     pause
% end