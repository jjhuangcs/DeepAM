% (c) Simon Hawe, Lehrstuhl fuer Datenverarbeitung Technische Universitaet
% Muenchen, 2012. Contact: simon.hawe@tum.de 
function A_matrix = draw_atoms(M, sz, per_row,figureNum)

sz = [sz,sz];

if iscell(M)
    Mn = M{end};
    for i=numel(M)-1:-1:1
        Mn = kron(Mn,M{i});
    end
    M=Mn;
end

cols = ceil(size(M,1)/per_row);
rand = 1;
if numel(sz)==3
    c = 3;
else
    c = 1;
end
A_matrix = zeros(cols*sz(1)+(cols+1)*rand,per_row*sz(2)+(per_row+1)*rand,c)+255;

i = 1;
for y=1:cols
        if i > size(M,1)
            break;
        end
    for x=1:per_row
        crnt = reshape(M(i,:),sz);
        
        ys = (y-1)*(sz(1)+rand)+1+rand;
        xs = (x-1)*(sz(2)+rand)+1+rand;
        value = (crnt/(2*max(max(abs(crnt(:)))))+0.5)*255;
        A_matrix(ys:ys+sz(1)-1,xs:xs+sz(2)-1,:) = value;

        i=i+1;
        if i > size(M,1)
            break;
        end
    end
end

if nargin == 3
    figure(123456);
    imshow(uint8(A_matrix));
    drawnow;
end

if nargin == 4
    figure(figureNum);
    imshow(uint8(A_matrix));
    drawnow;
end