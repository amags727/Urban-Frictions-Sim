function results = solveModel_final(param,data,settings)
    % 
    % Solve baseline model
    %
    struct2vars(param);
    struct2vars(data);
    struct2vars(settings);
    outerdiff = Inf;
    w = w_init;
    rr = rr_init;
    rf = rf_init;
    A = A_init;
    u = u_init;
    U_bar = U_bar_init;
    iter = 0;

    while outerdiff>tol
        w_tr = w';
        v = (((d.^(-1)).*w_tr).*(u.*(rr.^(beta1-1)))).^theta;
        v_welfare = (((d.^(-1*xi)).*w_tr).*(u.*(rr.^(beta1-1)))).^theta;
        W = gamma(1-(1/theta))*(sum(sum(v,2),1))^(1/theta);
        
        % Residents
        v_r = sum(v,2);
        res_share = v_r./sum(v_r,1); 
        L_r = (L_bar*res_share);

        % Workers
        v_f = sum(v,1);
        l_share = v_f./sum(v_f,2);
        L_f = (L_bar*l_share)';

        % Average wages
        w_adjust_theta = ((d.^(-1)).*w_tr).^(theta);
        CMA = sum(w_adjust_theta,2);
        pi_ij = w_adjust_theta./(sum(w_adjust_theta,2));
        w_bar = sum(pi_ij.*w_tr,2);

        %% Update
        w_prime = alpha1*A.*((H_f./L_f).^(1-alpha1));
        rr_prime = (1-beta1)*w_bar.*(L_r./H_r); 
        rf_prime = (1-alpha1)*(A.^(1/(1-alpha1))).*((w./alpha1).^(-alpha1/(1-alpha1))); 
        D_r = u_bar;
        Omega = sum(repmat(D_r',I,1).*(exp(-delta_u*t)),2);
        u_prime = Omega.^mu_u;
        A_prime = A_bar;
       
        Lr_welfcalc = sum(v_welfare(:))^(1/theta);
        if open_city==1
            L_bar_prime = (Lr_welfcalc./U_bar)'.*L_bar;
            U_bar_prime = U_bar;
        else
            L_bar_prime = L_bar;
            U_bar_prime = Lr_welfcalc;
        end

        w  = psi1*w_prime + (1-psi1)*w;
        w  = w./(prod(w)^(1/I));
        w_prime= w_prime./(prod(w_prime)^(1/I));
        rr = psi1*rr_prime + (1-psi1)*rr;
        rf = psi1*rf_prime + (1-psi1)*rf;
        u  = psi1*u_prime + (1-psi1)*u;
        A  = psi1*A_prime + (1-psi1)*A;
        L_bar = psi1*L_bar_prime + (1-psi1)*L_bar;
        U_bar = psi1*U_bar_prime + (1-psi1)*U_bar;

        w_diff = max(abs(w - w_prime));
        rr_diff = max(abs(rr - rr_prime));
        rf_diff = max(abs(rf - rf_prime));
        u_diff = max(abs(u - u_prime));
        A_diff = max(abs(A - A_prime));
        L_bar_diff = abs(L_bar_prime-L_bar);
        U_bar_diff = abs(U_bar_prime-U_bar);
        outerdiff = max([w_diff; rr_diff; rf_diff; u_diff; A_diff; L_bar_diff;U_bar_diff]);
  
       iter = iter + 1;
       if iter>1000
           disp("failure")
            break
       end
    end
    if iter<1000
    end 
    % Save and export
    results = v2struct(w,rr,rf,v,u,A,L_r,L_f,w_bar,pi_ij,res_share,l_share,W,Omega,u_bar,CMA,L_bar,U_bar);    

