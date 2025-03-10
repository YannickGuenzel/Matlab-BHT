function [p, s, TestStatDistribution, c] = BootHypTest(xSample, z, y, N_Boot, seed, nComparisons, userTestStat)
% BootHypTest performs a bootstrap-based randomization test for various 
% scenarios: one-sample, two-sample, two-sample-pairs, ranked consistency,
% ks-test, or a ranksum test-like test.
%
% Usage:
%   [p, s, TestStatDistribution, c] = BootstrapHypothesisTesting(xSample, z, y, N_Boot, seed, nComparisons, userTestStat)
%
% Inputs:
%   xSample      - Type of test ('one-sample', 'two-sample', 'two-sample-pairs',
%                  'ranked-consistency', 'ks-test', or 'ranksum').
%   z            - First sample (or original sample/matrix).
%   y            - Second sample or comparison value (defaults provided based on xSample).
%   N_Boot       - Number of bootstrap resamples (default: 5000).
%   seed         - Seed for reproducibility (default: 1234).
%   nComparisons - Number of multiple comparisons (default: 1).
%   userTestStat - (Optional) User-defined test statistic function.
%
% Outputs:
%   p                   - Holm–Bonferroni corrected p-value.
%   s                   - Shannon information: -log2(p).
%   TestStatDistribution- Struct with fields 'sample' (test statistic for the original data)
%                         and 'boot' (vector of bootstrap test statistics).
%   c                   - Cohen's d (effect size) if applicable; otherwise [].
%
% Version: 10-Mar-2025, Matlab R2023a

% Input Validation and Default Assignment -----------------------------
if nargin < 2
    error('Error: At least two inputs (xSample and z) are required.');
end
% Updated validTests list to include 'ranksum'
validTests = {'one-sample', 'two-sample', 'two-sample-pairs', 'ranked-consistency', 'ks-test', 'ranksum'};
if ~(ischar(xSample) || isstring(xSample))
    error('Error: xSample must be a string.');
end
xSample = char(xSample);
if ~ismember(xSample, validTests)
    error('Error: Unknown test type "%s".', xSample);
end
if ~isnumeric(z)
    error('Error: Input z must be numeric.');
end

% Set default seed and set RNG early to ensure reproducibility for random defaults
if nargin < 5 || isempty(seed)
    seed = 1234;
end
if ~isnumeric(seed)
    error('Error: seed must be numeric.');
end
rng(seed);

% Default for y depends on test type
if nargin < 3 || isempty(y)
    switch xSample
        case 'one-sample'
            y = 0;
        case 'two-sample'
            y = randn(size(z));
        case 'two-sample-pairs'
            y = zeros(size(z));
        case 'ranked-consistency'
            y = 1;
        case 'ks-test'
            y = randsample(z, length(z));
        case 'ranksum'
            y = randn(size(z));
    end
end

if nargin < 4 || isempty(N_Boot)
    N_Boot = 5000;
end
if ~isnumeric(N_Boot) || N_Boot <= 0
    error('Error: N_Boot must be a positive integer.');
end

if nargin < 6 || isempty(nComparisons) || nComparisons < 1
    nComparisons = 1;
end
if ~isnumeric(nComparisons)
    error('Error: nComparisons must be numeric.');
end

if nargin < 7
    userTestStat = [];
end

% Set Test Statistic Function -----------------------------------------
if ~isempty(userTestStat)
    TestStat = userTestStat;
else
    switch xSample
        case 'one-sample'
            TestStat = @(x1, L_x1, PredetVal) abs(mean(x1) - PredetVal) ./ (std(x1) / sqrt(L_x1));
        case 'two-sample'
            TestStat = @(x1, x2, L_x1, L_x2) abs(mean(x1) - mean(x2)) ./ sqrt(var(x1)/L_x1 + var(x2)/L_x2);
        case 'two-sample-pairs'
            TestStat = @(x1, L_x1) abs(mean(x1)) ./ (std(x1) / sqrt(L_x1));
        case 'ranked-consistency'
            TestStat = @(n, k, r) ((12*n)/(k*(k+1))) * sum((r - ((k+1)/2)).^2);
        case 'ks-test'
            TestStat = @(h1, h2) max(abs((cumsum(h1,1)./sum(h1,1)) - (cumsum(h2,1)./sum(h2,1))));
        case 'ranksum'
            TestStat = @ranksum_stat;
    end%switch
end%if user test

% Main Bootstrap Resampling -------------------------------------------
switch xSample
    case 'one-sample'
        z = z(:);
        z_tilde = z - mean(z) + y;
        bootIdx = randsample(1:length(z), length(z)*N_Boot, true);
        z_boot = reshape(z_tilde(bootIdx), length(z), N_Boot);
        TestStatDistribution.sample = TestStat(z, length(z), y);
        TestStatDistribution.boot = TestStat(z_boot, length(z), y);
        c = [];

    case 'two-sample'
        z = z(:);
        y = y(:);
        x = [z; y];
        x_boot = reshape(randsample(x, length(x)*N_Boot, true), length(x), N_Boot);
        z_boot = x_boot(1:length(z), :);
        y_boot = x_boot(length(z)+1:end, :);
        TestStatDistribution.sample = TestStat(z, y, length(z), length(y));
        TestStatDistribution.boot = TestStat(z_boot, y_boot, length(z), length(y));
        c = abs(mean(z) - mean(y)) / std(x);

    case 'two-sample-pairs'
        test_distr = [z(:); y(:)];
        z_boot = reshape(randsample(test_distr, length(z)*N_Boot, true), length(z), N_Boot);
        y_boot = reshape(randsample(test_distr, length(y)*N_Boot, true), length(y), N_Boot);
        x_diff = z(:) - y(:);
        x_boot = z_boot - y_boot;
        TestStatDistribution.sample = TestStat(x_diff, length(x_diff));
        TestStatDistribution.boot = TestStat(x_boot, size(x_boot,1));
        c = abs(mean(z) - mean(y)) / std(x_diff);

    case 'ranked-consistency'
        [nRows, nCols] = size(z);
        if y > 1
            if mod(nRows, y) ~= 0
                error('Error: The number of rows must be a multiple of y.');
            end
            numBlocks = nRows / y;
        else
            numBlocks = nRows;
            y = 1;
        end
        % Pre-compute block indices for the original sample
        blockIdx = arrayfun(@(j) ((j-1)*y + (1:y))', (1:numBlocks)', 'UniformOutput', false);
        m_original = cell2mat(cellfun(@(rows) reshape(tiedrank(reshape(z(rows, :), [], 1)), [y, nCols]), ...
            blockIdx, 'UniformOutput', false));
        avgRank_orig = mean(m_original, 1);
        TestStatDistribution.sample = TestStat(nRows, nCols, avgRank_orig);

        % Generate all bootstrap samples at once (3D array: rows x cols x N_Boot)
        numElements = numel(z);
        indices = randi(numElements, numElements, N_Boot);
        z_boot_all = reshape(z(indices), nRows, nCols, N_Boot);
        m_boot = zeros(size(z_boot_all));
        for j = 1:numBlocks
            rows = blockIdx{j};
            blockData = z_boot_all(rows, :, :);
            blockData_reshaped = reshape(blockData, y*nCols, N_Boot);
            blockRanks = tiedrank(blockData_reshaped);
            blockRanks = reshape(blockRanks, y, nCols, N_Boot);
            m_boot(rows, :, :) = blockRanks;
        end
        avgRank_all = squeeze(mean(m_boot, 1));  % size: (nCols x N_Boot)
        const_factor = (12 * nRows) / (nCols * (nCols+1));
        boot_stat = const_factor * sum((avgRank_all - ((nCols+1)/2)).^2, 1);  % 1 x N_Boot
        TestStatDistribution.boot = boot_stat(:);
        c = [];

    case 'ks-test'
        z = z(:);
        y = y(:);
        x = sum([z, y], 2);
        U = rand(length(x), N_Boot);
        z_boot = floor(U .* (x + 1));
        y_boot = repmat(x, 1, N_Boot) - z_boot;
        TestStatDistribution.sample = TestStat(z, y);
        TestStatDistribution.boot = TestStat(z_boot, y_boot);
        c = [];

    case 'ranksum'
        z = z(:);
        y = y(:);
        x = [z; y];
        x_boot = reshape(randsample(x, length(x)*N_Boot, true), length(x), N_Boot);
        z_boot = x_boot(1:length(z), :);
        y_boot = x_boot(length(z)+1:end, :);
        % Compute the ranksum test statistic for the original data
        TestStatDistribution.sample = TestStat(z, y);
        % For bootstrapped samples, compute the statistic for each resample
        TestStatDistribution.boot = zeros(N_Boot, 1);
        for i = 1:N_Boot
            TestStatDistribution.boot(i) = TestStat(z_boot(:,i), y_boot(:,i));
        end
        c = [];
end%switch

% Compute p-value -----------------------------------------------------
p = mean(TestStatDistribution.boot >= TestStatDistribution.sample);
if p == 0
    p = 1/N_Boot;
    warning('Test resulted in p=0. Thus it was set to 1/N_Boot.');
elseif p < 1/N_Boot
    p = 1/N_Boot;
    warning('Test resulted in p<1/N_Boot. Thus it was set to 1/N_Boot.');
end

% Holm–Bonferroni Correction ------------------------------------------
p = min(1, p * nComparisons);

% Compute Shannon Information -----------------------------------------
s = -log2(p);

    % Helper function for the ranksum test statistic ------------------
    function stat = ranksum_stat(x1, x2)
        % ranksum_stat computes the absolute deviation of the rank sum for x1 from its expected value.
        n1 = length(x1);
        n2 = length(x2);
        combined = [x1; x2];
        r = tiedrank(combined);
        expected = n1*(n1+n2+1)/2;
        stat = abs(sum(r(1:n1)) - expected);
    end%FCN:ranksum_stat
end%FCN:BootstrapHypothesisTesting
