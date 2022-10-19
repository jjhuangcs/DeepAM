function result = overlap_add( patches, img_size, grid )
% Image construction from overlapping patches
result = zeros(img_size,'single');
weight = zeros(img_size);
for i = 1:size(grid,3)
  patch = reshape(patches(:, i), size(grid, 1), size(grid, 2));
  result(grid(:, :, i)) = result(grid(:, :, i)) + patch;
  weight(grid(:, :, i)) = weight(grid(:, :, i)) + 1;
end
I = logical(weight);
result(I) = result(I) ./ weight(I);
end