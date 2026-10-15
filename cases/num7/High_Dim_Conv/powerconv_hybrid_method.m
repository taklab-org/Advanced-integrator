function method = powerconv_hybrid_method(v, p, opts)
% Choose FFT or direct for hybrid convolution routines.

opts = powerconv_hybrid_defaults(opts);
mode = lower(char(opts.mode));

if strcmp(mode, 'fft')
    method = 'fft';
    return
elseif strcmp(mode, 'direct')
    assert_direct_supported(v, p);
    method = 'direct';
    return
elseif ~strcmp(mode, 'auto')
    error('powerconv_hybrid_method: opts.mode must be auto, fft, or direct');
end

if ~is_direct_candidate(v, p)
    method = 'fft';
    return
end

sz = size(v);
max_len = max(sz);
is_vector_input = isvector(v);

intervals = exist('intval','file') && isintval(v(1));
has_radius = false;
if intervals
    has_radius = any(rad(v(:)) > 0);
end

if is_vector_input
    if p == 2 && max_len <= opts.max_1d_p2
        method = 'direct';
    elseif p == 3 && max_len <= opts.max_1d_p3
        method = 'direct';
    else
        method = 'fft';
    end
    return
end

if ~intervals
    if p == 2 && max_len <= opts.max_2d_double_p2
        method = 'direct';
    elseif p == 3 && max_len <= opts.max_2d_double_p3
        method = 'direct';
    else
        method = 'fft';
    end
    return
end

if has_radius
    if p == 2 && max_len <= opts.max_2d_interval_wide_p2
        method = 'direct';
    elseif p == 3 && max_len <= opts.max_2d_interval_wide_p3
        method = 'direct';
    else
        method = 'fft';
    end
else
    if p == 2 && max_len <= opts.max_2d_interval_point_p2
        method = 'direct';
    elseif p == 3 && max_len <= opts.max_2d_interval_point_p3
        method = 'direct';
    else
        method = 'fft';
    end
end
end

function tf = is_direct_candidate(v, p)
tf = (ndims(v) <= 2) && any(p == [2, 3]);
end

function assert_direct_supported(v, p)
if ~is_direct_candidate(v, p)
    error('powerconv_hybrid_method: direct mode supports only vectors/matrices with p = 2 or 3');
end
end
