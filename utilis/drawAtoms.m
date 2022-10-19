%% draw learned dictionary
function A = drawAtoms(M,figNum)

A = draw_atoms(M,sqrt(size(M,2)),round(sqrt(size(M,1))),figNum);

end