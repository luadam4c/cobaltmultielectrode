function filteredData = zof_mark (data, fc, npoles, si, varargin)
%% Uses a Butterworth filter twice to filter data (each column is a vector of samples)
% Usage: filteredData = zof_mark (data, fc, npoles, si, varargin)
% Explanation: 
%       A Butterworth filter is applied twice (once forward, once backward)
%       to remove lag effect of filter.
% Outputs:
%       filteredData - the filtered version of data
% Arguments:    
%       data    - data where each column is a vector of samples
%               must be a nonempty numeric array
%       fc      - the cutoff frequency(ies) (Hz) for the filter
%               must be a numeric and:
%                   a scalar by default or if ftype == 'low' or 'high'
%                   a two-element vector if ftype == 'bandpass' or 'stop'
%               consistent with the documentation for butter()
%       npoles  - order of filter
%                   i.e., the number of poles in the transfer function
%                   i.e., the order of the polynomial in the denominator 
%                           of the transfer function
%                   Note: the higher the order the steeper the cutoff
%               must be a numeric scalar that is a positive integer
%               consistent with the documentation for butter()
%       si      - sampling interval (seconds)
%               must be a numeric scalar
%       ftype   - (opt) filter type
%               must be an unambiguous, case-insensitive match to one of: 
%                   'low'       - lowpass filter with cutoff frequency fc
%                   'high'      - highpass filter with cutoff frequency fc
%                   'bandpass'  - bandpass filter of order 2*npoles 
%                                   with cutoff frequencies fc(1) & fc(2)
%                   'stop'      - bandstop filter of order 2*npoles 
%                                   with cutoff frequencies fc(1) & fc(2)
%               default == 'low' if fc has one element and 
%                       == 'bandpass' if fc has two elements
%
% Used by:    
%        apply this command in a LINUX terminal to find them:
%             grep --include=*.m -rlw '/home/Matlab/' -e "zof_mark"
%
% File History:
% ---------- Created by Mark P Beenhakker
% 2017-04-27 AL - Added inputParser scheme
% 2017-04-27 AL - Renamed cutoff frequency lp -> fc
% 2017-04-27 AL - Added nTraces
% 2017-04-27 AL - Added ftype as an argument for butter()
% 2017-04-27 AL - Renamed blpf -> numeratorCoeff; alpf -> denominatorCoeff
% 2017-04-27 AL - Renamed filttr to filtTr1, filtTr2 & filtTr3, respectively
% 2017-04-27 AL - Added Wn for clarity and consistency with the documentation 
%                   for butter()
% 2017-04-27 AL - Default for ftype is now 'bandpass' whenever fc 
%                   has two elements
%
% 2017-05-23 AL - Changed line width and indentation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Add required inputs to an input Parser
iP = inputParser;
addRequired(iP, 'data', ...     % vector of samples
    @(x) validateattributes(x, {'numeric'}, {'nonempty'}));
addRequired(iP, 'fc', ...       % the cutoff frequency(ies) (Hz) for the filter
    @(x) isnumeric(x) && isvector(x) && numel(x) <= 2);
addRequired(iP, 'npoles', ...   % order of filter
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive', 'integer'}));
addRequired(iP, 'si', ...       % sampling interval (seconds)
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));

%% Add optional inputs to the input Parser
addOptional(iP, 'ftype', 'low', ...        % filter type
    @(x) any(validatestring(x, {'low', 'high', 'bandpass', 'stop'})));

%% Read parameter values from the input Parser
parse(iP, data, fc, npoles, si, varargin{:});
ftype = validatestring(iP.Results.ftype, {'low', 'high', 'bandpass', 'stop'});

%% Change default for ftype if fc has two elements
if numel(fc) > 1 && strcmp(ftype, 'low')
    ftype = 'bandpass';
end

%% Filter each trace one by one
filteredData = zeros(size(data));       % initialize filtered data
nTraces = size(data, 2);                % number of traces to filter
for i = 1:nTraces
    % Extract the ith trace to filter
    thisTrace = data(:, i);             % ith trace to filter

    % Find the normalized cutoff frequency Wn = fc/(fs/2), 
    %   where fs = sampling frequency (Hz) = 1/si 
    %   and fs/2 is the Nyquist frequency
    Wn = fc * 2 * si;       % normalized cutoff frequency (half-cycles/sample)

    % Find the transfer function coefficients of a Butterworth filter
    %   with order npoles and normalized cutoff frequency Wn
    [numeratorCoeff, denominatorCoeff] = butter(npoles, Wn, ftype);    

    % Apply filter to trace
    filtTr1 = filter(numeratorCoeff, denominatorCoeff, thisTrace);
    
    % Filter the reversed trace to remove time lag from filter
    filtTr2 = filter(numeratorCoeff, denominatorCoeff, filtTr1(end:-1:1)); 

    % Reverse the trace to original orientation
    filtTr3 = filtTr2(end:-1:1);

    % Place filtered trace in output matrix
    filteredData(:, i) = filtTr3;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

function filtered_traces = zof_mark(data,lp,poles,si)
% data is a vector of samples
% lp is the low pass cutoff
% poles is number of poles?
% not so clear on meaning of "poles" here
% si is sample interval in sec
% output is filtered trace that was run through a butterworth
%filter twice, once forward, then once backward to remove lag effect of
%filter
filtered_traces = zeros(size(data,1), size(data,2)) ;
filtered_traces(:,1) = data(:,1) ;

for i = 1:size(data,2)
sweep = data(:,i) ;
[blpf,alpf]=butter(poles,lp*2*si);
filttr=filter(blpf,alpf,sweep);
filttr=filter(blpf,alpf,filttr(end:-1:1)); % filter the reversed trace to remove time lag from filter
filttr=filttr(end:-1:1); % reverse the trace to original orientation
filtered_traces(:,i) = filttr ;

end

%}