function avError = relError( varargin )
% Calculates relative error between multiple n-dimensional arrays.
%
% relError( A,B,C,... ) where size(A)=size(B)=size(C) calculates
% the error between all inputs and displays the error for each
% comparison.
%
% avError = relError( A,B,C,... ) outputs a 2D array of size
% nArgs x nArgs comparing each input to each input. Use 'display'
% parameter to suppress the output text.
%
% Example
%    A=[1 2]; B=[1 2]; C=[2 4];
%    err = relError(A,B,C,'display',false)
%
%    err =
%            A    B    C
%    A -->   0    0   100
%    B -->   0    0   100
%    C -->  50   50    0
%
% So for the error between the first argument and the third
% argument, use err(1,3).
%   Author: Daniel Simmons
%   Copyright 2014 - DansPhD.com
% Deal with not enough arguments
nArgs = length(varargin);
if nArgs<2, error('More arguments are needed'); end;
% Allow user to specify whether an output is displayed
show=true;
for a = 1:nArgs
    if ischar(varargin{a})
        if strcmpi(varargin{a},'display')
            show=varargin{a+1};
            nArgs=nArgs-2;
        end
    end
end
% Go though each combination of pairs of input arguments
avError = zeros(nArgs-1,nArgs-1);
for k=1:nArgs
    for j=k+1:nArgs
        avError(k,j) = compareArray(varargin{k}, varargin{j});
        avError(j,k) = compareArray(varargin{j}, varargin{k});
        if show
            fprintf('\nRelative error between arrays %i and %i = %.3g %%',...
                k,j,max(avError(j,k),avError(k,j)))
        end
    end
end
% Buffer the display output
if show
    fprintf('\n\n')
end
    % Private function comparing 1 pair of arguments
    function [ av_error ] = compareArray( A, B )
        
        % Deal with matrices with non-matching dimensions
        nA = ndims(A); nB = ndims(B);
        if nA~=nB,
            error('Matrices must have the same dimensions')
        end
        for i=1:nA
            N = size(A,i); M = size(B,i);
            if N ~= M
                error('First array has %i length, second has %i',N,M)
            end
        end
        
        % Translate arguments to 1 dimension
        A = A(:); B = B(:);
        
        % Deal with 0 matrices
        if all(A==0)
            error('Matrices must have values other than 0')
        end
        
        % Deal with infs
        inf_idx_A = find(isinf(A)); inf_idx_B = find(isinf(B));
        if numel(inf_idx_A)>0
            if numel(inf_idx_A) == numel(inf_idx_B)
                if inf_idx_A == inf_idx_B
                    A(inf_idx_A) = 10;      % Set inf in both arrays to any real number
                    B(inf_idx_B) = 10;      % if both have them at same position
                else
                    error('Array contains anomalous inf.')
                end
            else
                error('Array contains anomalous inf.')
            end
        end
        
        % Deal with NaNs
        nan_idx_A = find(isnan(A)); nan_idx_B = find(isnan(B));
        if numel(nan_idx_A)>0
            if numel(nan_idx_A) == numel(nan_idx_B)
                if nan_idx_A == nan_idx_B
                    A(nan_idx_A) = 10;      % Set NaN in both arrays to any real number
                    B(nan_idx_B) = 10;      % if both have them at same position
                else
                    error('Array contains anomalous NaN.')
                end
            else
                error('Array contains anomalous NaN.')
            end
        end
        
        %% Compute error
        % av_error = 100*mean(abs(A-B))/mean(abs(A));  % Old method
        
        if (all(isreal(A)) && all(isreal(B)))      % If both arrays are real...
            av_error = 100*( sqrt(sum((A-B).^2)) )/( sqrt(sum(A.^2)) );
        else                                       % If any array has imaginary values...
            av_error = 100*(norm(A-B)/norm(A));
        end
        
    end
end
