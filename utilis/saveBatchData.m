%% Save training data for different batches

function saveBatchData(X, Y, param)
warning off

mkdir(param.PATH);
rmdir(param.PATH, 's');
mkdir(param.PATH);

fprintf('Saving Batch\n');
for i=1:param.bnum
    if i==1 || mod(i,50)~=1
        fprintf('.');
    else
        fprintf('\n.');
    end
    
    %%
    idx = randi([1 size(X,2)],param.bsize,1);
    xl = X(:,idx);
    save(sprintf([param.PATH,'XL_%d'],i),'xl');
    clear xl
    
    xh = Y(:,idx);
    save(sprintf([param.PATH,'XH_%d'],i),'xh');
    clear xh
end
fprintf('\n');
end