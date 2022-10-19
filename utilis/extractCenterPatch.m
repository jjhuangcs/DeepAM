%% Extract the ennter region of the patches
function subMat = extractCenterPatch(in, psz)

[sx ~] = size(in);
psz0    = sqrt(sx);

subMat  = zeros(psz^2,sx);
offset  = (psz0 - psz)/2;

count1 = 1;
count2 = 1;
for i=1:psz0
    for j=1:psz0
        if i>= offset+1 && i<= psz0-offset
            if j>= offset+1 && j<= psz0-offset
                subMat(count2, count1) = 1;
                
                count2 = count2 + 1;
            end
        end
        count1 = count1 + 1;
    end
end

end