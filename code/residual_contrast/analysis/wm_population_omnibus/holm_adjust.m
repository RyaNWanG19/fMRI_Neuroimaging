function pAdjusted = holm_adjust(pValues)
%HOLM_ADJUST Holm step-down adjusted p-values, preserving input shape.

inputSize = size(pValues);
p = double(pValues(:));
if any(~isfinite(p) | p < 0 | p > 1)
    error('holm_adjust:InvalidPValues', ...
        'pValues must be finite and lie in [0, 1].');
end

m = numel(p);
[pSorted, order] = sort(p, 'ascend');
scaled = (m - (1:m)' + 1) .* pSorted;
adjustedSorted = min(1, cummax(scaled));

pAdjusted = zeros(m, 1);
pAdjusted(order) = adjustedSorted;
pAdjusted = reshape(pAdjusted, inputSize);
end
