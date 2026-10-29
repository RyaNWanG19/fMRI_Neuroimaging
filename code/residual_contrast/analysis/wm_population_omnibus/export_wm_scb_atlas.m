function export_wm_scb_atlas(outputDir)
% Export local canlab2018 parcel metadata and deterministic atlas boundaries.
% Uses the historical project's ordinal ROI convention, not a newly verified
% match to the atlas version used to extract the original responses.
if ~isfolder(outputDir), mkdir(outputDir); end
csvFile=fullfile(outputDir,'atlas_regions.csv');
niiFile=fullfile(outputDir,'canlab2018_labels.nii');
assert(~isfile(csvFile) && ~isfile(niiFile),'Use fresh atlas output files.');
a=load_atlas('canlab2018');
assert(numel(a.labels)==489 && size(a.dat,2)==1);
roi=(1:489)'; label=string(a.labels(:)); network=string(a.labels_2(:));
description=string(a.label_descriptions(:));
assert(numel(network)==489 && numel(description)==489);
T=table(roi,label,network,description,'VariableNames', ...
    {'atlas_roi_index','atlas_label','atlas_network','atlas_description'});
writetable(T,csvFile,'Encoding','UTF-8');
write(a,'fname',niiFile);
fid=fopen(fullfile(outputDir,'atlas_provenance.txt'),'w');
cleanup=onCleanup(@() fclose(fid));
fprintf(fid,'Local canlab2018 atlas; 489 parcels.\n');
fprintf(fid,'load_atlas implementation: %s\n',which('load_atlas'));
fprintf(fid,'atlas_name: %s\n',a.atlas_name);
fprintf(fid,'space: %s\n',a.space_description);
fprintf(fid,['Ordinal ROI convention follows project extraction code. ' ...
    'Historical extraction atlas identity cannot be independently established.\n']);
fprintf('Exported atlas metadata and deterministic labels to %s\n',outputDir);
end
