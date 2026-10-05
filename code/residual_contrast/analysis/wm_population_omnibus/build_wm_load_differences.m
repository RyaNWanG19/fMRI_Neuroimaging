function [D, mapping] = build_wm_load_differences(residuals)
% Four 2-back minus 0-back curves; output [ROI-category x time x subject].
% Source condition order: 0bk body/faces/places/tools, then 2bk same order.
sz = size(residuals); sz(end+1:4) = 1;
assert(isnumeric(residuals) && isreal(residuals) && numel(sz)==4 && sz(4)==8, ...
    'Expected a real numeric ROI x time x subject x 8 array.');
nr = sz(1); names = ["body";"faces";"places";"tools"];
D = zeros(nr*4,sz(2),sz(3));
for c = 1:4
    rows = (c-1)*nr+(1:nr);
    D(rows,:,:) = double(residuals(:,:,:,c+4))-double(residuals(:,:,:,c));
end
roi = repmat((1:nr)',4,1);
category = repelem(names,nr);
condition0 = repelem((1:4)',nr);
condition2 = condition0+4;
testID = category+"_ROI_"+string(roi);
mapping = table(testID,roi,category,condition0,condition2, ...
    'VariableNames',{'test_id','atlas_roi_index','category','condition_0back','condition_2back'});
end
