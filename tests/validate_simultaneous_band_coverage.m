function report = validate_simultaneous_band_coverage(outputDir, nRep, B, method)
% Reproducible finite-sample joint coverage experiment, not proof of validity.
% Families are four tests x 41 lags, smaller than the full WM family.
% Gaussian and asymmetric innovations; AR(1) temporal covariance plus shared
% subject variation across tests. Nonzero mean tests centering correctly.
if nargin<1, outputDir=''; end
if nargin<2, nRep=200; end
if nargin<3, B=999; end
if nargin<4, method='subject-bootstrap-t'; end
root=fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'code','residual_contrast','analysis','wm_population_omnibus'));
validateattributes(nRep,{'numeric'},{'scalar','integer','positive'});
validateattributes(B,{'numeric'},{'scalar','integer','>=',19});
names=["Gaussian_AR1";"Asymmetric_AR1"];
nt=41; m=4; n=200; rho=0.85;
time=(0:40)'*0.72;
truth=repmat(0.4*sin(time'/6),m,1);
cfg=struct('seed',1729,'alpha',0.05,'batchSize',100,'verbose',false, ...
    'hrfModelName','synthetic','runMode','verification', ...
    'scb',struct('nBootstrap',B,'minSubjects',10,'method',method));
covered=zeros(2,1); pointwiseCovered=zeros(2,1); widths=zeros(2,1);
for scenario=1:2
    rng(4181+scenario,'twister');
    for rep=1:nRep
        e=randn(n,m,nt);
        if scenario==2
            e=(exp(0.7*e)-exp(0.7^2/2))/sqrt((exp(0.7^2)-1)*exp(0.7^2));
        end
        noise=e;
        for t=2:nt
            noise(:,:,t)=rho*noise(:,:,t-1)+sqrt(1-rho^2)*e(:,:,t);
        end
        shared=randn(n,1);
        D=noise+reshape(0.5*shared,n,1,1)+reshape(truth,1,m,nt);
        primary=struct('config',cfg,'mean_D',reshape(mean(D,1),m,nt), ...
            'subjectMasks',true(n,m),'subjectIDs',string((1:n)'), ...
            'roiIDs',string((1:m)'), 'roiNames',string(compose('Synthetic_%d',(1:m)')), ...
            'effectiveN',repmat(n,m,1),'timeSec',time);
        primary.config.scb.seed=9901+1000*scenario+rep;
        bands=simultaneous_confidence_bands(D,primary);
        covered(scenario)=covered(scenario)+all(bands.lower<=truth & bands.upper>=truth,'all');
        error=abs(primary.mean_D-truth);
        pointwiseCovered(scenario)=pointwiseCovered(scenario)+all(error<=1.96*bands.standardError,'all');
        widths(scenario)=widths(scenario)+mean(bands.upper-bands.lower,'all');
        if scenario==1 && rep==1 && ~isempty(outputDir)
            if ~isfolder(outputDir), mkdir(outputDir); end
            primary.simultaneousBands=bands;
            plot_simultaneous_confidence_bands(primary,1:4,fullfile(outputDir,'synthetic_band_example.png'));
        end
        if mod(rep,50)==0, fprintf('Coverage %s: %d/%d experiments\n',names(scenario),rep,nRep); end
    end
end
coverage=covered/nRep;
mcSE=sqrt(coverage.*(1-coverage)/nRep);
report=table(names,repmat(string(method),2,1),repmat(n,2,1),repmat(m*nt,2,1),repmat(nRep,2,1), ...
    repmat(B,2,1),coverage,mcSE,pointwiseCovered/nRep,widths/nRep, ...
    'VariableNames',{'scenario','method','subjects','family_cells','experiments', ...
    'bootstrap_replicates','simultaneous_coverage','coverage_monte_carlo_se', ...
    'pointwise_joint_coverage','mean_band_width'});
disp(report);
if ~isempty(outputDir)
    if ~isfolder(outputDir), mkdir(outputDir); end
    reportFile=fullfile(outputDir,'coverage_validation.csv');
    assert(~isfile(reportFile),'Refusing to overwrite a coverage report.');
    writetable(report,reportFile);
end
end
