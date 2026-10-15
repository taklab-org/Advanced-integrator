function cutoff = standardChop(coeffs, tol)
% This code is a COPY of standardChop function by Chebfun.
% Written by Jared Aurentz and Nick Trefethen, July 2015.
% All copyright is by The University of Oxford and The Chebfun Developers.

% Step 1: Convert COEFFS to a new monotonically nonincreasing
%         vector ENVELOPE normalized to begin with the value 1.

[n,sizeCheb]  = size(coeffs);
b = abs(coeffs);
m = b(end,:).*ones(n,sizeCheb);
for j = n-1:-1:1
    m(j,:) = max(b(j,:), m(j+1,:));
end
if all( m(1,:) == 0 )
    cutoff = 1;
    return
end
envelope = m./m(1,:);

% Step 2: Scan ENVELOPE for a value PLATEAUPOINT, the first point J-1, if any,
% that is followed by a plateau.  A plateau is a stretch of coefficients
% ENVELOPE(J),...,ENVELOPE(J2), J2 = round(1.25*J+5) <= N, with the property
% that ENVELOPE(J2)/ENVELOPE(J) > R.  The number R ranges from R = 0 if
% ENVELOPE(J) = TOL up to R = 1 if ENVELOPE(J) = TOL^(2/3).  Thus a potential
% plateau whose starting value is ENVELOPE(J) ~ TOL^(2/3) has to be perfectly
% flat to count, whereas with ENVELOPE(J) ~ TOL it doesn't have to be flat at
% all.  If a plateau point is found, then we know we are going to chop the
% vector, but the precise chopping point CUTOFF still remains to be determined
% in Step 3.
plateauPoint = zeros(sizeCheb,1);
cutoff = n*ones(sizeCheb,1);

for ell = 1:sizeCheb
    for j = 2:n
        j2 = round(1.25*j + 5);
        if ( j2 > n )
            % there is no plateau: exit
            cutoff = max(cutoff);
            return
        end
        e1 = envelope(j,ell);
        e2 = envelope(j2,ell);
        r = 3*(1 - log(e1)/log(tol));
        plateau = (e1 == 0) | (e2/e1 > r);
        if ( plateau )
            % a plateau has been found: go to Step 3
            plateauPoint(ell) = j - 1;
            break
        end
    end
% end

% Step 3: fix CUTOFF at a point where ENVELOPE, plus a linear function
% included to bias the result towards the left end, is minimal.
%
% Some explanation is needed here.  One might imagine that if a plateau is
% found, then one should simply set CUTOFF = PLATEAUPOINT and be done, without
% the need for a Step 3. However, sometimes CUTOFF should be smaller or larger
% than PLATEAUPOINT, and that is what Step 3 achieves.
%
% CUTOFF should be smaller than PLATEAUPOINT if the last few coefficients made
% negligible improvement but just managed to bring the vector ENVELOPE below the
% level TOL^(2/3), above which no plateau will ever be detected.  This part of
% the code is important for avoiding situations where a coefficient vector is
% chopped at a point that looks "obviously wrong" with PLOTCOEFFS.
%
% CUTOFF should be larger than PLATEAUPOINT if, although a plateau has been
% found, one can nevertheless reduce the amplitude of the coefficients a good
% deal further by taking more of them.  This will happen most often when a
% plateau is detected at an amplitude close to TOL, because in this case, the
% "plateau" need not be very flat.  This part of the code is important to
% getting an extra digit or two beyond the minimal prescribed accuracy when it
% is easy to do so.



% for ell = 1:sizeCheb
if ( envelope(plateauPoint(ell)) == 0 )
    cutoff(ell) = plateauPoint(ell);
else
    j3 = sum(envelope(:,ell) >= tol^(7/6));
    if ( j3 < j2 )
        j2 = j3 + 1;
        envelope(j2,ell) = tol^(7/6);
    end
    cc = log10(envelope(1:j2,ell));
    cc = cc(:);
    cc = cc + linspace(0, (-1/3)*log10(tol), j2)';
    [~, d] = min(cc);
    cutoff(ell) = max(d - 1, 1);
end

end
cutoff = max(cutoff);
