function output = prepare_medium_dsge_data(rawFile, outputFile, options)
%PREPARE_MEDIUM_DSGE_DATA 生成第四章中等规模DSGE模型的七条观测序列。
%
% 直接运行：
%   prepare_medium_dsge_data
%
% 自定义输入、输出与样本：
%   opts = struct('sampleStart', datetime(1984,1,1), ...
%                 'sampleEnd', datetime(2007,4,1), ...
%                 'filterStart', datetime(1960,1,1), 'hpLambda', 1600);
%   out = prepare_medium_dsge_data('us_dsge_raw.xlsx', ...
%                                 'dsge_observables.mat', opts);
%
% 默认输入、输出均位于本m文件所在目录，不依赖当前工作目录。
% 输入Excel：
%   Quarterly!A:H: Date GDPC1 GDP PCEC FPI GDPDEF HOANBS COMPNFB
%   Monthly!A:B:   Date FEDFUNDS
%
% MAT文件中的七个列向量可直接供Dynare varobs使用：
%   GDPobs Cobs Iobs Nobs wobs Piobs Robs
% 对应文档的GDP^obs、C^obs、I^obs、N^obs、w^obs、Pi^obs、R^obs。
% 前五条为原始对数水平的单侧HP周期项乘100。
% 后两条为季度对数总通胀率、总利率乘100后减去估计样本均值。
% GDP、C、I、N采用总量口径。真实工资直接对应模型的w。
%
% 仅需基础MATLAB，无需安装Dynare或Econometrics Toolbox。
% 原始发布机构为BEA、BLS、Federal Reserve，数据由FRED提供。
% 单侧HP使用Stock and Watson (1999, JME)的状态空间表述：
% https://doi.org/10.1016/S0304-3932(99)00027-6
% 下载版本与各序列的单位、来源见原始Excel右侧说明。

dataDir = fileparts(mfilename('fullpath'));
if nargin < 1 || isempty(rawFile)
    rawFile = fullfile(dataDir, 'us_dsge_raw.xlsx');
end
if nargin < 2 || isempty(outputFile)
    outputFile = fullfile(dataDir, 'dsge_observables.mat');
end
if nargin < 3 || isempty(options)
    options = struct();
end
cfg = struct('sampleStart', datetime(1984,1,1), ...
    'sampleEnd', datetime(2007,4,1), ...
    'filterStart', datetime(1960,1,1), 'hpLambda', 1600);
names = fieldnames(options);
for j = 1:numel(names)
    if ~isfield(cfg, names{j})
        error('DSGE:Option', 'Unknown option: %s', names{j});
    end
    cfg.(names{j}) = options.(names{j});
end
dateFields = {'sampleStart','sampleEnd','filterStart'};
for j = 1:numel(dateFields)
    value = cfg.(dateFields{j});
    if ~isdatetime(value), value = datetime(value); end
    if ~isscalar(value) || isnat(value) || day(value) ~= 1 || ...
            ~ismember(month(value), [1 4 7 10])
        error('DSGE:Date', '%s must be a quarter-start date.', dateFields{j});
    end
    cfg.(dateFields{j}) = value;
end
validateattributes(cfg.hpLambda, {'numeric'}, {'scalar','real','finite','positive'});
if cfg.filterStart >= cfg.sampleStart || cfg.sampleStart > cfg.sampleEnd
    error('DSGE:Sample', 'Require filterStart < sampleStart <= sampleEnd.');
end
if ~isfile(rawFile)
    error('DSGE:Input', 'Raw workbook not found: %s', rawFile);
end
if strcmpi(char(rawFile), char(outputFile))
    error('DSGE:Output', 'Input workbook and output MAT paths must differ.');
end

%% 1. 读取原始频率数据。右侧元数据不进入数值表。
Q = readtable(rawFile, 'Sheet','Quarterly', 'Range','A:H', ...
    'VariableNamingRule','preserve', 'UseExcel',false);
M = readtable(rawFile, 'Sheet','Monthly', 'Range','A:B', ...
    'VariableNamingRule','preserve', 'UseExcel',false);
qNames = {'GDPC1','GDP','PCEC','FPI','GDPDEF','HOANBS','COMPNFB'};
require_columns(Q, [{'Date'},qNames], 'Quarterly');
require_columns(M, {'Date','FEDFUNDS'}, 'Monthly');
Q.Date = import_dates(Q.Date);
M.Date = import_dates(M.Date);
Q = Q(~isnat(Q.Date), :);
M = M(~isnat(M.Date), :);
Q = sortrows(Q, 'Date');
M = sortrows(M, 'Date');
if numel(unique(Q.Date)) ~= height(Q) || numel(unique(M.Date)) ~= height(M)
    error('DSGE:Duplicate', 'Duplicate dates in raw workbook.');
end
if any(day(Q.Date) ~= 1 | ~ismember(month(Q.Date), [1 4 7 10])) || ...
        any(day(M.Date) ~= 1)
    error('DSGE:Frequency', 'Quarterly/monthly dates must be period start dates.');
end
for j = 1:numel(qNames)
    Q.(qNames{j}) = numeric_column(Q.(qNames{j}));
end
M.FEDFUNDS = numeric_column(M.FEDFUNDS);
qIndex = 4*year(Q.Date) + ceil(month(Q.Date)/3) - 1;
mIndex = 12*year(M.Date) + month(M.Date) - 1;
if any(diff(qIndex) ~= 1) || any(diff(mIndex) ~= 1)
    error('DSGE:Gap', 'Raw date grid has a missing quarter or month.');
end

%% 2. 月度联邦基金利率取季度算术平均，要求完整的3个月。
mQuarter = 4*year(M.Date) + ceil(month(M.Date)/3) - 1;
[rateQuarter, ~, group] = unique(mQuarter);
quarterCount = accumarray(group, 1);
validCount = accumarray(group, double(isfinite(M.FEDFUNDS)));
meanFunds = accumarray(group, M.FEDFUNDS, [], @mean);
meanFunds(quarterCount ~= 3 | validCount ~= 3) = NaN;
[hasRate, position] = ismember(qIndex, rateQuarter);
fundsQuarter = NaN(height(Q),1);
fundsQuarter(hasRate) = meanFunds(position(hasRate));
% 不完整的最新季度留为NaN。需要该季度的估计样本会在下文明确报错。

%% 3. 价格平减与流量单位转换。GDPDEF的基期值为100。
price = Q.GDPDEF / 100;
levels = [Q.GDPC1/4, Q.PCEC./price/4, Q.FPI./price/4, ...
          Q.HOANBS, Q.COMPNFB./price];
levelNames = {'GDP','C','I','N','w'};
% GDP、C、I: 十亿2017年美元/季度。消费、投资均用GDP平减指数。
% N: 总工时指数。w: 用GDP平减指数平减的名义小时报酬指数。
% 不对任何总量除以人口。每小时报酬继续使用小时口径。
inflationLog = [NaN; 100*diff(log(Q.GDPDEF))];
interestLog = 100*log(1 + fundsQuarter/400);

use = Q.Date >= cfg.filterStart & Q.Date <= cfg.sampleEnd;
expectedDates = (cfg.filterStart:calmonths(3):cfg.sampleEnd)';
if ~isequal(Q.Date(use), expectedDates)
    error('DSGE:Coverage', 'Quarterly data do not cover the requested filter window.');
end
if any(~isfinite(levels(use,:)) | levels(use,:) <= 0, 'all')
    error('DSGE:Levels', 'Missing or nonpositive level data within filter window.');
end
estimateRows = Q.Date >= cfg.sampleStart & Q.Date <= cfg.sampleEnd;
if any(~isfinite(inflationLog(estimateRows))) || ...
        any(~isfinite(interestLog(estimateRows))) || ...
        any(~isfinite(fundsQuarter(estimateRows)))
    error('DSGE:Rates', 'Missing inflation or a complete three-month interest rate.');
end

%% 4. 对对数水平进行单侧HP滤波。先滤波完整预样本，再截取估计样本。
filterDates = Q.Date(use);
logLevels = log(levels(use,:));
[hpTrend, hpCycle] = one_sided_hp(logLevels, cfg.hpLambda);
select = filterDates >= cfg.sampleStart;
observations = 100*hpCycle(select,:);
pi0 = mean(inflationLog(estimateRows));
i0 = mean(interestLog(estimateRows));
observations = [observations, inflationLog(estimateRows)-pi0, ...
    interestLog(estimateRows)-i0];
quarter_dates = Q.Date(estimateRows);
expectedSample = (cfg.sampleStart:calmonths(3):cfg.sampleEnd)';
if ~isequal(quarter_dates, expectedSample) || ...
        size(observations,1) ~= numel(expectedSample) || ...
        any(~isfinite(observations),'all')
    error('DSGE:Observations', 'Invalid observation matrix or sample alignment.');
end
% 单侧性检查：删除未来观测后，已估计的趋势保持不变。
prefixLength = max(2, floor(size(logLevels,1)/2));
[prefixTrend, ~] = one_sided_hp(logLevels(1:prefixLength,:), cfg.hpLambda);
prefixError = max(abs(prefixTrend-hpTrend(1:prefixLength,:)), [], 'all');
if prefixError > 1e-10
    error('DSGE:HP', 'One-sided filtering failed the prefix invariance check.');
end

%% 5. 保存七条列向量、日期、原始校准目标和处理信息。
output = struct();
obs_names = {'GDPobs','Cobs','Iobs','Nobs','wobs','Piobs','Robs'};
for j = 1:numel(obs_names)
    output.(obs_names{j}) = observations(:,j);
end
output.obs_names = obs_names;
output.observations = observations;
output.quarter_dates = quarter_dates;
output.year_quarter = [year(quarter_dates), ceil(month(quarter_dates)/3)];
output.quarter_labels = cellstr(string(year(quarter_dates)) + "Q" + ...
    string(ceil(month(quarter_dates)/3)));

% 稳态目标使用未滤波数据，计算区间与估计区间相同。
calibration = struct();
calibration.s_C = mean(Q.PCEC(estimateRows)./Q.GDP(estimateRows));
calibration.s_I = mean(Q.FPI(estimateRows)./Q.GDP(estimateRows));
calibration.s_G = 1-calibration.s_C-calibration.s_I;
calibration.pi0 = pi0;
calibration.i0 = i0;
calibration.Pi = exp(pi0/100);
calibration.R = exp(i0/100);
calibration.beta = calibration.Pi/calibration.R;
if any(~isfinite([calibration.s_C,calibration.s_I,calibration.s_G])) || ...
        calibration.s_C <= 0 || calibration.s_I <= 0 || calibration.s_G <= 0 || ...
        calibration.beta <= 0 || calibration.beta >= 1
    error('DSGE:Calibration', 'Sample moments imply inadmissible steady-state targets.');
end
output.calibration = calibration;
output.processing = struct('dates',filterDates, 'level_names',{levelNames}, ...
    'levels',levels(use,:), 'log_levels',logLevels, 'hp_trend',hpTrend, ...
    'hp_cycle',hpCycle, 'inflation_log',inflationLog(use), ...
    'interest_log',interestLog(use), 'federal_funds_annual_percent',fundsQuarter(use));
output.metadata = struct('settings',cfg, 'raw_workbook',char(rawFile), ...
    'created_at',datestr(now,'yyyy-mm-dd HH:MM:SS'), ...
    'source_ids',{{'GDPC1','GDP','PCEC','FPI','GDPDEF','HOANBS','COMPNFB','FEDFUNDS'}}, ...
    'source_platform','FRED, Federal Reserve Bank of St. Louis', ...
    'data_vintage','Current revised data. Retrieval date and source URLs are recorded in Excel.', ...
    'hp_initial_state','First log observation repeated twice; predicted covariance 1e5*eye(2).', ...
    'hp_method','Causal Kalman filter for an integrated second-order random-walk trend.', ...
    'units','GDPobs Cobs Iobs Nobs wobs: 100 times log cycle. Piobs Robs: quarterly log percentage points.', ...
    'population_adjustment','None', 'wage_mapping','Hourly compensation deflated by GDPDEF maps to w.', ...
    'prefix_invariance_error',prefixError);
[outputDir,~,extension] = fileparts(char(outputFile));
if ~strcmpi(extension,'.mat')
    error('DSGE:Output', 'Output filename must end in .mat.');
end
if ~isempty(outputDir) && ~isfolder(outputDir), mkdir(outputDir); end
save(outputFile, '-struct','output', '-v7');
fprintf('Saved %s\n', char(outputFile));
fprintf('Sample: %s to %s. Observations: %d x 7.\n', ...
    output.quarter_labels{1}, output.quarter_labels{end}, size(observations,1));
fprintf('HP initialization: %s to %s. Lambda = %.0f.\n', ...
    datestr(cfg.filterStart,'yyyy-mm-dd'), ...
    datestr(cfg.sampleStart-calmonths(3),'yyyy-mm-dd'),cfg.hpLambda);
fprintf('sC = %.6f, sI = %.6f, sG = %.6f, beta = %.6f.\n', ...
    calibration.s_C,calibration.s_I,calibration.s_G,calibration.beta);
end

function require_columns(T, required, sheetName)
if ~all(ismember(required, T.Properties.VariableNames))
    error('DSGE:Columns', 'Missing required columns in sheet %s.', sheetName);
end
end

function result = import_dates(input)
if isdatetime(input)
    result = input;
elseif isnumeric(input)
    result = datetime(input,'ConvertFrom','excel');
else
    strings = string(input);
    result = NaT(size(strings));
    valid = ~ismissing(strings) & strlength(strings)>0;
    % Excel exports may store ISO datetime text. The first 10 characters are the date.
    result(valid) = datetime(extractBefore(strings(valid)+"T",11), ...
        'InputFormat','yyyy-MM-dd');
end
result = result(:);
end

function result = numeric_column(input)
if isnumeric(input)
    result = double(input);
else
    result = str2double(string(input));
end
result = result(:);
end

function [trend, cycle] = one_sided_hp(y, lambda)
% y_t = tau_t + e_t, tau_t = 2*tau_{t-1}-tau_{t-2}+v_t.
% Var(e)=1, Var(v)=1/lambda. 每一步只使用当期及过去观测。
% 初始化仅用第一期数据，不使用第二期或估计期末数据。
F = [2 -1; 1 0];
H = [1 0];
Q = diag([1/lambda,0]);
[n,k] = size(y);
trend = zeros(n,k);
state = [y(1,:);y(1,:)];
P = 1e5*eye(2);                 % 第一期更新前的预测协方差
for t = 1:n
    gain = P*H'/(H*P*H'+1);
    filtered = state + gain*(y(t,:)-H*state);
    trend(t,:) = H*filtered;
    A = eye(2)-gain*H;
    covariance = A*P*A' + gain*gain';  % Joseph形式保持数值对称与半正定
    state = F*filtered;
    P = F*covariance*F' + Q;
    P = (P+P')/2;
end
cycle = y-trend;
end

