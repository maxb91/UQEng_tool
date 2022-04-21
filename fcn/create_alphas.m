function [p_index, p_index_roots] = create_alphas(M, P)
% [p_index, p_index_roots] = CREATE_ALPHAS(M, P): generate the indices for a polynomial chaos basis
% to be used in evaluating the basis. M is the number of variables and P
% is the maximum polynomial order.

%% parsing and checking the arguments for options and consistency
% parse the command line
% keywords to be parsed
% parse_keys = {'q-norm', 'max-interaction'};
% type ('f' = simple variable, 'p' = key, value pair)
% parse_types = {'p', 'p'};
% parsed = uq_simple_parser(varargin, parse_keys, parse_types);

% now use the parsed inputs to set the options of our basis creation
% q-norm
% if ~strcmp(parsed{1}, 'false')
%     truncation = 1;
%     q = parsed{1};
% else
% default to no truncation
truncation = 0;
% end

% max interaction terms
% if ~strcmp(parsed{2}, 'false')
%     J = 1:min(M,parsed{2});
% end

% default to maximum interaction terms if they don't exist
% if ~exist('J', 'var')
J = 1:min(M,max(P));
% end
P = 0:P;

%% Accumulating the indices as a function of J
% we want this function to be more the less vectorized, or at least support vectorize
% synthax for future improvements
NJ = numel(J);
NP = numel(P);

% we will store the root generators in a 2D cell array. It will be used for q-norm
% selection
p_index_roots = cell(NP, NJ);
p_index = cell(NP, NJ);

% now we generate all the unique J-tuples of degree P
for pp = 1:NP
    for jjj = 1:NJ
        if P(pp) < J(jjj)
            continue;
        end
        % getting the jtuple
        cur_jtuple = generate_jtuples_sum_p(P(pp),J(jjj));
        jtuple = zeros(size(cur_jtuple,1), M);
        try
            % padding to the correct directon to the left (Blatman Thesis, Appendix C)
            %jtuple = uq_padarray(jtuple', M - J(jj), 'pre')';
            
            jtuple(:, (M-size(cur_jtuple,2)+1):end) = cur_jtuple;
            
        catch me
            fprintf('Failed to pad jtuple M = %d jj = %d', M, jjj);
            rethrow(me);
        end
        
        %%% THERE IS NO TRUNCTAION!!
        % applying truncation if necessary
        %         if truncation
        %             switch truncation
        %                 case 1 % q-norm selection
        %                     q_norm = uq_q_norm(jtuple,q);
        %                     jtuple = squeeze(jtuple(q_norm <= max(P),:));
        %             end
        %         else
        q_norm = sum(jtuple,2);
        jtuple = squeeze(jtuple(q_norm <= max(P),:));
        %         end
        %%%
        
        % save the outputs
        p_index_roots{pp,jjj} = jtuple;
        
        % now let's get the full coefficients matrix:
        
        % first we need to loop over every line of the jtuple
        tmpidx = cell(size(jtuple,1),1); % store the intermediate results in a temporary index
        % loopuq
        for kk = 1:size(jtuple,1)
            tmpidx{kk} = double(permute_coefficients(jtuple(kk,:)));
        end
        % and now concatenate all of them together in a sparse matrix
        p_index{pp,jjj} = sparse(cat(1, tmpidx{:}));
    end % end of loop over the polynomial degree
end % end of loop ver ntuples

%% this is a reordering routine  that we won't use at this stage
% clear p_index_roots;
% for ii = 1:NP
%   p_index_roots{ii} = cat(1, p_index{ii,:});
%     [~, pidx] = sortrows(fliplr(p_index_roots{ii}));
%   p_index_roots{ii} = p_index_roots{ii}(pidx,:);
% end

%%  this is where we reshape the final output as desired
% ok, now let's put everything in a unique large bidimensional array ordered rowwise by P
% and columnwise by the variable index (up to M). This is very large!

p_index = reshape(p_index, numel(p_index),1);
if ~min(P) % if we have 0, let's add the constant term
    p_index = [zeros(1,M); p_index];
end

p_index = cat(1,p_index{:});
p_index = full(p_index);
end

function PERMS = generate_jtuples_sum_p(P, J)
%  fuction PERMS = GENERATE_JTUPLES_SUM_P(P,J) generates all the J-tuples of integers
%  a(i) that satisfy the following conditions using Knuth's H algorithm:
%    - a(i) < P
%    - sum(a) = P
%    - a(i+1) < a(i)
%  The additional constraints must be met: P >= J

%% integrity checks
if P < J
    error('uq_generate_jtuples_sum_p: Error, specified P is larger than J!!');
end


%% trivial cases: J = 0 and J = 1
if J == 0
    PERMS = 0 ;
elseif J == 1
    PERMS = P ;
else
    %% main loop
    % let's at the moment return a horizontal vector, but it may be better to return it
    % as a column vector instead. The comments in this section reflect those in the
    % decription of the Knuth algorithm
    
    PERMS = zeros(1, J) ;
    
    %  "Initialize"
    Y = ones(1, J+1) ;
    Y(1) = P-J+1 ;
    Y(J+1) = -1 ;
    
    i = 0 ;
    while 1
        % "Visit"
        i = i + 1 ;
        PERMS(i,:) = Y(1:J) ;
        if Y(2) < Y(1)-1
            % "Tweak" Y(1) and Y(2)
            Y(1) = Y(1) - 1 ;
            Y(2) = Y(2) + 1 ;
        else
            % "Find" j
            j = 3 ;
            s = Y(1) + Y(2) - 1 ;
            while Y(j) >= Y(1) - 1
                s = s + Y(j) ;
                j = j + 1 ;
            end
            % "Increase" Y(j)
            if j > J
                break
            else
                z = Y(j) + 1 ;
                Y(j) = z ;
                j = j-1 ;
            end
            % "Tweak" Z(1) ... Z(j)
            while j > 1
                Y(j) = z;
                s = s - z ;
                j = j - 1 ;
            end
            Y(1) = s ;
        end
    end
    
end

% and finally sort the values
PERMS = sortrows(PERMS);
end


function myperms = permute_coefficients(myset)
% calculating unique permutations of a vector. Using chunk allocation
% initializations

myset = int8(myset);
NTYPE = 1; % corresponding to int8

max_size_in_memory = 8192; % max size in MB of the permutations matrix


% preallocating variables to specify their type
j = zeros(1,1,'int32') ;
l = zeros(1,1,'int32') ;
k = zeros(1,1,'int32') ;
MM = zeros(1,1,'int32') ;
N = zeros(1,1,'int32') ;
i = zeros(1,1,'int32') ; %% let's see if this helps


MM = length(myset) ;
N = length(myset(myset>0));
%tmp_set = myset ;

nrows = prod(double(MM-N+1:MM)); % number of rows in case of all different non-zero elements in myset

% now getting the multiplicity of each of the unique elements in myset
% un = unique(myset(myset>0));
% nunique = length(un);
% for i = 1:nunique
%       mult(i) = length(find(myset(myset>0) == un(i)));
% end

un = myset(myset>0);
totn = 0;
curn = numel(un);
mult = zeros(curn, 1);
unvalue = mult;
ii = 1;
while totn < N
    %curn = numel(un);
    curel = min(un);
    un = un(un>curel);
    tmpcurn = numel(un);
    mult(ii) = curn - tmpcurn;
    unvalue(ii) = curel;
    totn = totn + mult(ii);
    curn = tmpcurn;
    ii = ii + 1;
end

mult = mult(1:ii-1);
unvalue = unvalue(1:ii-1);
% get uniques, multiplicity and reorder
multcumul = cumsum(mult);
tmp_set = zeros(1,MM);
idx = MM - multcumul(end) + 1;
multcumul = multcumul + idx - 1;

for jj = 1:ii-1
    tmp_set(idx:multcumul(jj)) = unvalue(jj);
    idx = multcumul(jj) + 1;
end

%tmp_set = padarray(tmp_set, [0 (M-multcumul(end))], 'pre');

% final number of rows, taking multiplicity into consideration
nrows = nrows / prod(factorial(mult));

% check for total memory
mfingerprint = nrows*MM*NTYPE/2^20; % mem fingerprint in MB
if  mfingerprint > max_size_in_memory
    error('number of permutations too high: would require %d MB, while the specified max is %d MB\n', mfingerprint, max_size_in_memory);
end

% % allocate the necessary memory
myperms = zeros(nrows,MM, 'int8');

% works with row vectors only
tmp_set = reshape(tmp_set, 1, MM);


i = 0;
while 1
    %   L1. Visit
    i=i+1 ;
    myperms(i,:) = tmp_set ;
    
    
    %   L2. Find j
    j = MM - 1 ;
    while j && tmp_set(j) >= tmp_set(j+1)
        j = j - 1 ;
    end
    
    if ~j
        break
    end
    
    %   L3. Increase aj
    l = MM ;
    while tmp_set(j) >= tmp_set(l)
        l = l-1 ;
    end
    
    aux = tmp_set(j) ;
    tmp_set(j) = tmp_set(l) ;
    tmp_set(l) = aux ;
    %   L4. Reverse aj+1...aM
    k = j+1 ;
    l = MM ;
    while k<l
        if ~(tmp_set(k) || tmp_set(l))% do not exchange zero entries
            k=k+1 ; l=l-1 ;
            continue;
        end
        aux = tmp_set(k) ;
        tmp_set(k) = tmp_set(l) ;
        tmp_set(l) = aux ;
        k=k+1 ; l=l-1 ;
    end
end
end


