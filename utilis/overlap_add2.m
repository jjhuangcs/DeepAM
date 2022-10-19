function result = overlap_add2( patches, img_size, grid )
% Image construction from overlapping patches
result = zeros(img_size,'single');
weight = zeros(img_size);

psz     = sqrt(size(patches,1));
psz0    = size(grid, 1);
offset  = (psz0 - psz)/2;

for i = 1:size(grid,3)
  patch = reshape(patches(:, i), psz, psz);
  result(grid(offset+1:psz0-offset, offset+1:psz0-offset, i)) = result(grid(offset+1:psz0-offset, offset+1:psz0-offset, i)) + patch;
  weight(grid(offset+1:psz0-offset, offset+1:psz0-offset, i)) = weight(grid(offset+1:psz0-offset, offset+1:psz0-offset, i)) + 1;
end
I = logical(weight);
result(I) = result(I) ./ weight(I);
end