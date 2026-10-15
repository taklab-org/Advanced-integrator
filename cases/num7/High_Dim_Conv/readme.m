%% High_Dim_Conv/readme.m
% Usage notes for the hybrid convolution routines.
%
% The hybrid functions automatically choose between:
%   1. direct convolution in small 1D/2D cases with p = 2 or 3
%   2. FFT-based convolution otherwise
%
% Available hybrid entry points:
%   [vpshort,vpfull,vpextended,method] = powerconvcos_hybrid(v,p,opts)
%   [vpshort,vpfull,vpextended,method] = powerconvsin_hybrid(v,p,opts)
%   [vpshort,vpfull,method]           = powerconvtaylor_hybrid(v,p,opts)
%
% Input/output conventions are the same as the original functions:
%   powerconvcos / powerconvsin / powerconvtaylor
% The extra output METHOD is either 'direct' or 'fft'.
%
% ---------------------------------------------------------------------
% Basic usage
% ---------------------------------------------------------------------
%
% Cosine series:
%   v = [1, 0.3, 0.1];
%   [vpshort,vpfull,vpextended,method] = powerconvcos_hybrid(v,2);
%
% Sine series:
%   v = [0, 1, 0.5];
%   [vpshort,vpfull,vpextended,method] = powerconvsin_hybrid(v,3);
%
% Taylor series:
%   v = [1, 1];
%   [vpshort,vpfull,method] = powerconvtaylor_hybrid(v,2);
%
% ---------------------------------------------------------------------
% Automatic switching
% ---------------------------------------------------------------------
%
% Direct mode is currently used only when all of the following hold:
%   - the input is 1D or 2D
%   - p is 2 or 3
%   - the problem size is below an internal threshold
%
% Default thresholds:
%   - 1D: p=2 up to length 64, p=3 up to length 32
%   - 2D double: up to 8x8
%   - 2D INTLAB point intervals: up to 4x4
%   - 2D INTLAB intervals with nonzero radius:
%       p=2 up to 8x8, p=3 up to 6x6
%
% These defaults reflect the current tradeoff:
%   - direct convolution can give tighter interval bounds
%   - FFT is usually faster once the problem is not small
%
% ---------------------------------------------------------------------
% Options
% ---------------------------------------------------------------------
%
% The optional struct OPTS can override the decision:
%   opts.mode = 'auto'   % default
%   opts.mode = 'fft'    % always use FFT implementation
%   opts.mode = 'direct' % always use direct implementation
%
% Example:
%   opts.mode = 'direct';
%   [vpshort,vpfull,vpextended,method] = powerconvcos_hybrid(v,2,opts);
%
% You can also override thresholds:
%   opts.max_1d_p2
%   opts.max_1d_p3
%   opts.max_2d_double_p2
%   opts.max_2d_double_p3
%   opts.max_2d_interval_point_p2
%   opts.max_2d_interval_point_p3
%   opts.max_2d_interval_wide_p2
%   opts.max_2d_interval_wide_p3
%
% Example:
%   opts.max_2d_double_p2 = 12;
%   [vpshort,vpfull,vpextended,method] = powerconvcos_hybrid(v,2,opts);
%
% ---------------------------------------------------------------------
% Notes
% ---------------------------------------------------------------------
%
% 1. Direct mode is not implemented for dimensions >= 3.
%    In those cases the hybrid routines fall back to FFT.
%
% 2. Direct mode is not implemented for p > 3.
%    In those cases the hybrid routines fall back to FFT.
%
% 3. If opts.mode = 'direct' is requested for an unsupported case,
%    the hybrid routine throws an error.
