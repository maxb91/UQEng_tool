%% Compare the model responses with the responses of the metamodel for increasing polynomial degrees
neval = 20;
degrees_conv = ones(neval,1)*[1:5];
incdeg_conv = arrayfun(@(x)struct('p',x),degrees_conv);
for nn = 1:neval
    fprintf('\n Evalution %i:\n',nn)
    for pp = degrees_conv(1,:) % polynomial degree
        tic;
        % FIRST CREATE PCE MODEL OF DEGREE p
        incdeg_conv(nn,pp).P = factorial(M+pp)/(factorial(M)*factorial(pp));% Compute cardinality
        incdeg_conv(nn,pp).n = k*incdeg_conv(nn,pp).P; % resulting number of samples

        % Get sampled points in range [0,1]
        incdeg_conv(nn,pp).u01 = lhsdesign(incdeg_conv(nn,pp).n,M);

        % Transform -- Experimental design
        incdeg_conv(nn,pp).x_sigsr = sigsr_lim(1)+incdeg_conv(nn,pp).u01(:,1)*(sigsr_lim(2)-sigsr_lim(1));
        incdeg_conv(nn,pp).x_srm = srm_lim(1)+incdeg_conv(nn,pp).u01(:,2)*(srm_lim(2)-srm_lim(1));
        incdeg_conv(nn,pp).x_taub1 = tau_lim(1)+incdeg_conv(nn,pp).u01(:,3)*(tau_lim(2)-tau_lim(1)); % According to TCM (Alvarez 1998): taub1 = fct
        incdeg_conv(nn,pp).X = [incdeg_conv(nn,pp).x_sigsr, incdeg_conv(nn,pp).x_srm, incdeg_conv(nn,pp).x_taub1];
        incdeg_conv(nn,pp).y_ED  =M_avg_strain(incdeg_conv(nn,pp).x_sigsr,incdeg_conv(nn,pp).x_srm,incdeg_conv(nn,pp).x_taub1);

        % Transform to [-1,1] for Legendre
        incdeg_conv(nn,pp).X_uniform = incdeg_conv(nn,pp).X;
        for i=1:3
            incdeg_conv(nn,pp).X_uniform(:,i) = ((incdeg_conv(nn,pp).X(:,i) - sampling_limits(i,1)) / (sampling_limits(i,2)-sampling_limits(i,1)))*2-1;
        end

        % Get alphas for Legendre polynomials
        [incdeg_conv(nn,pp).p_index, incdeg_conv(nn,pp).p_index_roots] = create_alphas(M, pp);

        % Set up Psi matrix
        incdeg_conv(nn,pp).Psi = ones(incdeg_conv(nn,pp).n,incdeg_conv(nn,pp).P);
        for i=1:incdeg_conv(nn,pp).n
            for j=1:incdeg_conv(nn,pp).P
                % For each element of the matrix, compute the multivariate basis by
                % a product of the univariate ones:
                incdeg_conv(nn,pp).Psi(i,j) = eval_legendre(incdeg_conv(nn,pp).X_uniform(i,1),incdeg_conv(nn,pp).p_index(j,1))*...
                    eval_legendre(incdeg_conv(nn,pp).X_uniform(i,2),incdeg_conv(nn,pp).p_index(j,2))*...
                    eval_legendre(incdeg_conv(nn,pp).X_uniform(i,3),incdeg_conv(nn,pp).p_index(j,3));
            end
        end

        % PCE coefficients can be obtained from the following equation:
        incdeg_conv(nn,pp).c=pinv(transpose(incdeg_conv(nn,pp).Psi)*incdeg_conv(nn,pp).Psi)*transpose(incdeg_conv(nn,pp).Psi)*incdeg_conv(nn,pp).y_ED;

        % Evaluate PCE on ED
        incdeg_conv(nn,pp).y_PCE=0;
        for j=1:incdeg_conv(nn,pp).P
            incdeg_conv(nn,pp).y_PCE=incdeg_conv(nn,pp).y_PCE+incdeg_conv(nn,pp).c(j)*incdeg_conv(nn,pp).Psi(:,j);
        end

        % SECOND: EVALUATE THE PCE MODEL ON VALIDTION SAMPLE

        % Set up Psi matrix based on validation set (always the same validation
        % set)
        incdeg_conv(nn,pp).Psi_val = ones(nv,incdeg_conv(nn,pp).P);
        for i=1:nv
            for j=1:incdeg_conv(nn,pp).P
                incdeg_conv(nn,pp).Psi_val(i,j) = eval_legendre(XV_uniform(i,1),incdeg_conv(nn,pp).p_index(j,1))*...
                    eval_legendre(XV_uniform(i,2),incdeg_conv(nn,pp).p_index(j,2))*...
                    eval_legendre(XV_uniform(i,3),incdeg_conv(nn,pp).p_index(j,3));
            end
        end

        % Value of PCE based on validation set
        incdeg_conv(nn,pp).y_PCE_val=0;
        for j=1:incdeg_conv(nn,pp).P
            incdeg_conv(nn,pp).y_PCE_val=incdeg_conv(nn,pp).y_PCE_val+incdeg_conv(nn,pp).c(j)*incdeg_conv(nn,pp).Psi_val(:,j);
        end
        incdeg_conv(nn,pp).mu_PCE = mean(incdeg_conv(nn,pp).y_PCE_val);
        incdeg_conv(nn,pp).var_PCE=var(incdeg_conv(nn,pp).y_PCE_val);
        incdeg_conv(nn,pp).t=toc;
        fprintf('PCE with degree %i done, mu = %4.4f, var = %4.4e, t = %4.2f \n',pp,incdeg_conv(nn,pp).mu_PCE,incdeg_conv(nn,pp).var_PCE,incdeg_conv(nn,pp).t)
    end

    for pp = degrees_conv(1,:)
        incdeg_conv(nn,pp).muhat_PCE = incdeg_conv(nn,pp).c(1);
        incdeg_conv(nn,pp).varhat_PCE =sum(incdeg_conv(nn,pp).c(2:end).*incdeg_conv(nn,pp).c(2:end));
        incdeg_conv(nn,pp).cvhat_PCE =sqrt(incdeg_conv(nn,pp).varhat_PCE)/incdeg_conv(nn,pp).muhat_PCE;
        incdeg_conv(nn,pp).h=diag(incdeg_conv(nn,pp).Psi*pinv(transpose(incdeg_conv(nn,pp).Psi)*incdeg_conv(nn,pp).Psi)*transpose(incdeg_conv(nn,pp).Psi));
        incdeg_conv(nn,pp).epsLOO = (1./incdeg_conv(nn,pp).n).*...
            sum(((incdeg_conv(nn,pp).y_ED-incdeg_conv(nn,pp).y_PCE)./(1-incdeg_conv(nn,pp).h)).^2)./...
            sum((incdeg_conv(nn,pp).y_ED-mean(incdeg_conv(nn,pp).y_ED)).^2);
        incdeg_conv(nn,pp).epsMSQ = sum((y_ED_val-incdeg_conv(nn,pp).y_PCE_val).^2)./...
            sum((y_ED_val-mean(y_ED_val)).^2);
        incdeg_conv(nn,pp).epsMSQ_ED = sum((incdeg_conv(nn,pp).y_ED-incdeg_conv(nn,pp).y_PCE).^2)./...
            sum((incdeg_conv(nn,pp).y_ED-mean(incdeg_conv(nn,pp).y_ED)).^2);
    end
end

%% OVERSampling
%
% p = 4;
oversampling_conv = ones(neval,1)*[0.5:0.5:4];
incn_conv = arrayfun(@(x)struct('oversampling',x),oversampling_conv);
for nn=1:neval
    fprintf('\n Evalution %i:\n',nn)
    for kk = 1:length(oversampling_conv(1,:))
        tic;
        incn_conv(nn,kk).total_degree = p;
        

        % FIRST CREATE PCE MODEL OF DEGREE p
        incn_conv(nn,kk).P = factorial(M+p)/(factorial(M)*factorial(p));% Compute cardinality
        incn_conv(nn,kk).n = ceil(incn_conv(nn,kk).oversampling*incn_conv(nn,kk).P); % resulting number of samples

        % Get sampled points in range [0,1]
        incn_conv(nn,kk).u01 = lhsdesign(incn_conv(nn,kk).n,M);

        % Transform -- Experimental design
        incn_conv(nn,kk).x_sigsr = sigsr_lim(1)+incn_conv(nn,kk).u01(:,1)*(sigsr_lim(2)-sigsr_lim(1));
        incn_conv(nn,kk).x_srm = srm_lim(1)+incn_conv(nn,kk).u01(:,2)*(srm_lim(2)-srm_lim(1));
        incn_conv(nn,kk).x_taub1 = tau_lim(1)+incn_conv(nn,kk).u01(:,3)*(tau_lim(2)-tau_lim(1)); % According to TCM (Alvarez 1998): taub1 = fct
        incn_conv(nn,kk).X = [incn_conv(nn,kk).x_sigsr, incn_conv(nn,kk).x_srm, incn_conv(nn,kk).x_taub1];
        incn_conv(nn,kk).y_ED  =M_avg_strain(incn_conv(nn,kk).x_sigsr,incn_conv(nn,kk).x_srm,incn_conv(nn,kk).x_taub1);

        % Transform to [-1,1] for Legendre
        incn_conv(nn,kk).X_uniform = incn_conv(nn,kk).X;
        for i=1:3
            incn_conv(nn,kk).X_uniform(:,i) = ((incn_conv(nn,kk).X(:,i) - sampling_limits(i,1)) / (sampling_limits(i,2)-sampling_limits(i,1)))*2-1;
        end

        % Get alphas for Legendre polynomials
        [incn_conv(nn,kk).p_index, incn_conv(nn,kk).p_index_roots] = create_alphas(M, p);

        % Set up Psi matrix
        incn_conv(nn,kk).Psi = ones(incn_conv(nn,kk).n,incn_conv(nn,kk).P);
        for i=1:incn_conv(nn,kk).n
            for j=1:incn_conv(nn,kk).P
                % For each element of the matrix, compute the multivariate basis by
                % a product of the univariate ones:
                incn_conv(nn,kk).Psi(i,j) = eval_legendre(incn_conv(nn,kk).X_uniform(i,1),incn_conv(nn,kk).p_index(j,1))*...
                    eval_legendre(incn_conv(nn,kk).X_uniform(i,2),incn_conv(nn,kk).p_index(j,2))*...
                    eval_legendre(incn_conv(nn,kk).X_uniform(i,3),incn_conv(nn,kk).p_index(j,3));
            end
        end

        % PCE coefficients can be obtained from the following equation:
        incn_conv(nn,kk).c=pinv(transpose(incn_conv(nn,kk).Psi)*incn_conv(nn,kk).Psi)*transpose(incn_conv(nn,kk).Psi)*incn_conv(nn,kk).y_ED;

        % Evaluate PCE on ED
        incn_conv(nn,kk).y_PCE=0;
        for j=1:incn_conv(nn,kk).P
            incn_conv(nn,kk).y_PCE=incn_conv(nn,kk).y_PCE+incn_conv(nn,kk).c(j)*incn_conv(nn,kk).Psi(:,j);
        end


        % SECOND: EVALUATE THE PCE MODEL ON VALIDTION SAMPLE

        % Set up Psi matrix based on validation set (always the same validation
        % set)
        incn_conv(nn,kk).Psi_val = ones(nv,incn_conv(nn,kk).P);
        for i=1:nv
            for j=1:incn_conv(nn,kk).P
                incn_conv(nn,kk).Psi_val(i,j) = eval_legendre(XV_uniform(i,1),incn_conv(nn,kk).p_index(j,1))*...
                    eval_legendre(XV_uniform(i,2),incn_conv(nn,kk).p_index(j,2))*...
                    eval_legendre(XV_uniform(i,3),incn_conv(nn,kk).p_index(j,3));
            end
        end

        % Value of PCE based on validation set
        incn_conv(nn,kk).y_PCE_val=0;
        for j=1:incn_conv(nn,kk).P
            incn_conv(nn,kk).y_PCE_val=incn_conv(nn,kk).y_PCE_val+incn_conv(nn,kk).c(j)*incn_conv(nn,kk).Psi_val(:,j);
        end
        incn_conv(nn,kk).mu_PCE = mean(incn_conv(nn,kk).y_PCE_val);
        incn_conv(nn,kk).var_PCE=var(incn_conv(nn,kk).y_PCE_val);
        incn_conv(nn,kk).t=toc;
        fprintf('PCE with experimental design of %i samples done, mu = %4.4f, var = %4.4e, t = %4.2f \n',incn_conv(nn,kk).n,incn_conv(nn,kk).mu_PCE,incn_conv(nn,kk).var_PCE,incn_conv(nn,kk).t)
    end

    for kk = 1:length(oversampling_conv(1,:))
        incn_conv(nn,kk).muhat_PCE = incn_conv(nn,kk).c(1);
        incn_conv(nn,kk).varhat_PCE =sum(incn_conv(nn,kk).c(2:end).*incn_conv(nn,kk).c(2:end));
        incn_conv(nn,kk).cvhat_PCE =sqrt(incn_conv(nn,kk).varhat_PCE)/incn_conv(nn,kk).muhat_PCE;
        incn_conv(nn,kk).h=diag(incn_conv(nn,kk).Psi*pinv(transpose(incn_conv(nn,kk).Psi)*incn_conv(nn,kk).Psi)*transpose(incn_conv(nn,kk).Psi));
        incn_conv(nn,kk).epsLOO = (1./incn_conv(nn,kk).n).*...
            sum(((incn_conv(nn,kk).y_ED-incn_conv(nn,kk).y_PCE)./(1-incn_conv(nn,kk).h)).^2)./...
            sum((incn_conv(nn,kk).y_ED-mean(incn_conv(nn,kk).y_ED)).^2);
        incn_conv(nn,kk).epsMSQ = sum((y_ED_val-incn_conv(nn,kk).y_PCE_val).^2)./...
            sum((y_ED_val-mean(y_ED_val)).^2);
        incn_conv(nn,kk).epsMSQ_ED = sum((incn_conv(nn,kk).y_ED-incn_conv(nn,kk).y_PCE).^2)./...
            sum((incn_conv(nn,kk).y_ED-mean(incn_conv(nn,kk).y_ED)).^2);
    end
end

