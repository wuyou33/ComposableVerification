A=reshape(1:108,6,6,3); 
c=mat2cell(A,2*ones(1,3),2*ones(1,3),1*ones(1,3));

csize = size(c);
% Treat 3+ dimension arrays

% Construct the matrix by concatenating each dimension of the cell array into
%   a temporary cell array, CT
% The exterior loop iterates one time less than the number of dimensions,
%   and the final dimension (dimension 1) concatenation occurs after the loops

% Loop through the cell array dimensions in reverse order to perform the
%   sequential concatenations
for cdim=(length(csize)-1):-1:1
    % Pre-calculated outside the next loop for efficiency
    ct = cell([csize(1:cdim) 1]);
    cts = size(ct);
    ctsl = length(cts);
    mref = {};

    % Concatenate the dimension, (CDIM+1), at each element in the temporary cell
    %   array, CT
    for mind=1:prod(cts)
        [mref{1:ctsl}] = ind2sub(cts,mind);
        % Treat a size [N 1] array as size [N], since this is how the indices
        %   are found to calculate CT
        if ctsl==2 && cts(2)==1
            mref = {mref{1}};
        end
        % Perform the concatenation along the (CDIM+1) dimension
        ct{mref{:}} = cat(cdim+1,c{mref{:},:});
    end
    % Replace M with the new temporarily concatenated cell array, CT
    c = ct;
end

% Finally, concatenate the final rows of cells into a matrix
m = cat(1,c{:});


% A_one=A_cell(:,:,1);
% permcell = permute(A_cell,[1 2 3]);
% permnum=[permcell{:}];
% ans=reshape(permnum,2,2,9,3);
% ans=num2cell(ans,[1,2]);
% disp(A_one{2,2});
% ans_one=ans(:,:,1);
% disp(ans_one{:});
% C = reshape(C,[],size(A_cell,2),1);
