function showDict(Omega)

[~,sz] = size(Omega);
Eff = 1;
for i=1:sz
    if i==1
        Eff = Omega{i}*Eff;
        drawAtoms(Omega{i},i);
        figure(i);title(sprintf('Layer %d Analysis Dictionary', i))
    elseif i < sz
        Eff = Omega{i}*Eff;
        drawAtoms(Omega{i},i);
        figure(i);title(sprintf('Layer %d Analysis Dictionary', i))
    else
        drawAtoms(Eff,i);
        figure(i);title('Effective Analysis Dictionary')
        
        drawAtoms(Omega{i}',i+1);
        figure(i+1);title('Synthesis Dictionary')
    end
end

end