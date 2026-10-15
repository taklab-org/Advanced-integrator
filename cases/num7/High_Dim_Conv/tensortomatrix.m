function MT = tensortomatrix(T,M)
% function MT = tensortomatrix(T,M)
% reshapes a tensor T representing a matrix into a matrix MT
% here M=[M(1) M(2) M(3)] contains the Fourier-Chebyshev-Fourier type

sizeT=size(T);
MM=sum(M);
if length(sizeT)==2*MM
    MT=reshape(T,[prod(sizeT(1:MM)),prod(sizeT(MM+1:2*MM))]);
elseif length(sizeT)==2*sum(M)+2
    T=permute(T,[1,3:MM+2,2,MM+3:2*MM+2]);
    sizeTT=size(T);
    MT=reshape(T,[prod(sizeTT(1:MM+1)),prod(sizeTT(MM+2:2*MM+2))]);
else
    error('T does not have size a compatible with M');
end
