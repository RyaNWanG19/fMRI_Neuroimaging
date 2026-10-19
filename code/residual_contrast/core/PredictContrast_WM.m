function PredictContrast_WM(DATA_FILE, OUTDIR, cfg)

if nargin < 3, cfg = struct(); end

if nargin<2 || isempty(OUTDIR), OUTDIR = fullfile(pwd,'outputs'); end
if ~exist(OUTDIR,'dir'), mkdir(OUTDIR); end

tmp = load(DATA_FILE, 'Results'); 
assert(isfield(tmp,'Results'), 'Results missing in %s', DATA_FILE);
Rfile = tmp.Results; clear tmp;

sz = size(Rfile); sz(end+1:4)=1;

% Canonicalize: [S x R x T x C]
if     sz(1)==489 && sz(2)==41          % [R x T x S x C] -> [S x R x T x C]
    R = permute(Rfile, [3 1 2 4]);
elseif sz(2)==489 && sz(3)==41          % already [S x R x T x C]
    R = Rfile;
else
    error('Unexpected dims for WM file: [%s]. Expect [489 x 41 x S x 8] or [S x 489 x 41 x 8].', num2str(sz));
end

[S,Rg,T,C] = size(R);
assert(Rg==489 && T==41 && C==8, 'WM expected R=489,T=41,C=8; got [%d,%d,%d]', Rg,T,C);

Diff_Body  = R(:,:,:,5) - R(:,:,:,1);
Diff_Face  = R(:,:,:,6) - R(:,:,:,2);
Diff_Place = R(:,:,:,7) - R(:,:,:,3);
Diff_Tool  = R(:,:,:,8) - R(:,:,:,4);

Results = cat(4, Diff_Body, Diff_Face, Diff_Place, Diff_Tool);  % [S x R x T x 4]
cond_labels = {'Body_LoadDiff','Face_LoadDiff','Place_LoadDiff','Tool_LoadDiff'};

COMP2 = nchoosek(1:4, 2);

temp_file = fullfile(OUTDIR, sprintf('temp_WM_Interaction_%s.mat', datestr(now,'HHMMSSFFF')));
AlreadyCanonical = true;
save(temp_file, 'Results','AlreadyCanonical','-v7.3');

try
    run_observed_stage_flexible(temp_file, OUTDIR, COMP2, cond_labels, 'WM_Interaction', cfg);
catch ME
    if exist(temp_file,'file'), delete(temp_file); end
    rethrow(ME);
end

if exist(temp_file,'file'), delete(temp_file); end
end
