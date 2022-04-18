classdef MESP
    properties
        C double
        C_comp double %inv(C)
        size double
        A double
        b double
        F double
        Fsquare double
        F_comp double
        Fsquare_comp double
        ldetC double
        scaleC double % buffer variable
    end

    %% intialization function
    methods  
        function obj = MESP(C, A, b)
            if nargin < 0
                error("Please input valid values.");
            end
            obj.C=C;
            obj.C_comp=inv(C);
            obj.size=length(C);
            obj.A = A;
            obj.b = b;
            [obj.F,obj.Fsquare,obj.ldetC] = gen_data(C,0);
            [obj.F_comp,obj.Fsquare_comp,~] = gen_data(C,1);
        end
    end

    %% methods for obtain lower bound
    methods
        function [lb,info] = obtain_lb(obj,s)
        % obtain the lower bound and corresponding x indices 
        obtain_lb_inline;
        end
    end

    %% methods block for factorization bound
    methods
        function [fval,dx,info] = DDFact_obj(obj,x,s,Gamma)
        % This function calculate the objective value, gradient, and info of DDFact for given data
        %{
        Input:
        x       - current point for the DDFact problem
        s       - the size of subset we want to choose, also equals to the summation of all elements of x
        Gamma   - diagonal scaling paramter newC = Diag(Gamma)*C*Diag(Gamma)

        Output:
        fval    - objective value of DDFact at current point x
        dx      - the gradient of the obejctive function of DDFact at x
        info    - struct containing necesssary information
        %}
        DDFact_obj_inline;
        end
        
        function [fval,dx] = DDFact_obj_knitro(obj,x,s,Gamma)
            % create a callback function for Knitro specifying objective value and gradient 
            [fval,dx,~] = obj.DDFact_obj(x,s,Gamma);
            fval=-fval;
            dx=-dx;
        end

        function [fval,x,info] = Knitro_DDFact(obj,x0,s,Gamma)
        % calling knitro to solve the DDFact problem
        %{
        Input:
        s       - the size of subset we want to choose, also equals to the
                  summation of all elements of x0
        x0      - initial point
        Gamma   - diagonal scaling paramter newC = Diag(Gamma)*C*Diag(Gamma)

        Output:
        fval    - objective value of DDFact at optimal solution x
        x       - optimal solution x
        info    - struct containing necesssary information
        %}
        Knitro_DDFact_inline;
        end

        function [optGamma,info]=BFGS_DDFact_Gamma(obj,s,GammaInit)
        % BFGS method for optimizaing symmetric diagonal scaling parameter
        % of DDFact objective function (factorization bound)
        BFGS_DDFact_Gamma_inline;
        end
    end
    
    %% methods for the complementaty factorization bound
    methods
        function [fval,dx,info] = DDFact_comp_obj(obj,x,s,Gamma)
        % This function calculate the objective value, gradient, and info of comp DDFact for given data
        %{
        Input:
        x       - current point for the DDFact problem
        s       - the size of subset we want to choose, also equals to the summation of all elements of x
        Gamma   - diagonal scaling paramter newC = Diag(Gamma)*C*Diag(Gamma)

        Output:
        fval    - objective value of DDFact at current point x
        dx      - the gradient of the obejctive function of DDFact at x
        info    - struct containing necesssary information
        %}
        DDFact_comp_obj_inline;
        end
        
        function [fval,dx] = DDFact_comp_obj_knitro(obj,x,s,Gamma)
        % create a callback function for Knitro specifying objective value and gradient 
        [fval,dx,~] = obj.DDFact_comp_obj(x,s,Gamma);
        fval=-fval;
        dx=-dx;
        end

        function [fval,x,info] = Knitro_DDFact_comp(obj,x0,s,Gamma)
        % calling knitro to solve the complementary DDFact problem
        %{
        Input:
        s       - the size of subset we want to choose, also equals to the
                  summation of all elements of x0
        x0      - initial point
        Gamma   - diagonal scaling paramter newC = Diag(Gamma)*C*Diag(Gamma)

        Output:
        fval    - objective value of DDFact at optimal solution x
        x       - optimal solution x
        info    - struct containing necesssary information
        %}
        Knitro_DDFact_comp_inline;
        end

        function [optGamma,info]=BFGS_DDFact_comp_Gamma(obj,s,GammaInit)
        % BFGS method for optimizaing symmetric diagonal scaling parameter
        % of DDFact objective function (factorization bound)
        BFGS_DDFact_comp_Gamma_inline;
        end
    end
    

    %% methods for the linx bound with row scaling
    methods
        function [fval,dx,info] = Linx_obj(obj,x,s,Gamma)
        % This function calculate the objective value and gradient of the objective function of linx
        %{
        Input: 
        x       - current point for linx bound
        s       - the size of subset we want to choose, also equals to the summation of all elements of x
        Gamma   - row diagonal scaling parameter
        
        Output:
        fval    - the linx bound at current point x
        dx      - the gradient of obejctive function of Linx at x
        info    - struct containing necesssary information
        %}
        Linx_obj_inline;
        end
        
        function [fval,dx] = Linx_obj_knitro(obj,x,s,Gamma)
            % create a callback function for Knitro specifying objective value and gradient 
            [fval,dx,~] = obj.Linx_obj(x,s,Gamma);
            fval=-fval;
            dx=-dx;
        end

        function [fval,x,info] = Knitro_Linx(obj,x0,s,Gamma)
        % calling knitro to solve the DDFact problem
        %{
        Input:
        s       - the size of subset we want to choose, also equals to the
                  summation of all elements of x0
        x0      - initial point
        Gamma   - row diagonal scaling paramter newC = Diag(Gamma)*C
        comp    - indicate complementary bound, comp=1 means we need to
                  compare with the lower bound associated with n-s

        Output:
        fval    - objective value of DDFact at optimal solution x
        x       - optimal solution x
        info    - struct containing necesssary information
        %}
        Knitro_Linx_inline;
        end

        function [optGamma,info]=BFGS_Linx_Gamma(obj,s,GammaInit)
        % BFGS method for optimizaing row diagonal scaling parameter
        % of Linx objective function (factorization bound)
        BFGS_Linx_Gamma_inline;
        end
    end

    %% mixing DDFact and Linx
    methods
        function [fval,dx] = mix_DDFact_Linx_obj_knitro(obj,x,s,Gamma1,Gamma2,alpha)
        % create a callback function for Knitro specifying objective value and gradient 
        % This function calculate the objective value, gradient, and info
        % of mixing Linx&DDFact bound
        %{
        Input:
        x       - current point for the mix problem
        s       - the size of subset we want to choose, also equals to the summation of all elements of x
        Gamma1  - symmetric diagonal scaling paramter for DDFact
        Gamma2  - row diagonal scaling paramter for Linx
        alpha   - mixing parameter
    
        Output:
        fval    - objective value at current point x
        dx      - the gradient of the obejctive function at x
        %}
        mix_DDFact_Linx_obj_knitro_inline;
        end

        function [fval,dx,info] = mix_DDFact_Linx(obj,s,Gamma1,Gamma2)
        % mixing Linx and DDFact bound (with optimal mixing parameter)
        %{
        Input:
        s       - the size of subset we want to choose, also equals to the summation of all elements of x
        Gamma1  - diagonal scaling parameter for DDFact
        Gamma2  - row scaling parameter for Linx

        Output:
        fval    - objective value at current point x
        dx      - the gradient of the obejctive function at x
        info    - struct containing necesssary information
        %}
        mix_DDFact_Linx_inline;
        end
    end 

    %% mixing DDFactcomp and Linx
    methods
        function [fval,dx] = mix_DDFact_comp_Linx_obj_knitro(obj,x,s,Gamma1,Gamma2,alpha)
        % create a callback function for Knitro specifying objective value and gradient 
        % This function calculate the objective value, gradient, and info
        % of mixing Linx&DDFact comp bound
        %{
        Input:
        x       - current point for the mix problem
        s       - the size of subset we want to choose, also equals to the summation of all elements of x
        Gamma1  - symmetric diagonal scaling paramter for DDFact comp
        Gamma2  - row diagonal scaling paramter for Linx
        alpha   - mixing parameter
    
        Output:
        fval    - objective value at current point x
        dx      - the gradient of the obejctive function at x
        %}
        mix_DDFact_comp_Linx_obj_knitro_inline;
        end

        function [fval,dx,info] = mix_DDFact_comp_Linx(obj,s,Gamma1,Gamma2)
        % mixing Linx and DDFact comp bound (with optimal mixing parameter)
        %{
        Input:
        s       - the size of subset we want to choose, also equals to the summation of all elements of x
        Gamma1  - diagonal scaling parameter for DDFact comp
        Gamma2  - row scaling parameter for Linx

        Output:
        fval    - objective value at current point x
        dx      - the gradient of the obejctive function at x
        info    - struct containing necesssary information
        %}
        mix_DDFact_comp_Linx_inline;
        end
    end

    %% mixing DDFact and DDFact comp
    methods
        function [fval,dx] = mix_DDFact_DDFact_comp_obj_knitro(obj,x,s,Gamma1,Gamma2,alpha)
        % create a callback function for Knitro specifying objective value and gradient 
        % This function calculate the objective value, gradient, and info
        % of mixing DDFact&DDFact comp bound
        %{
        Input:
        x       - current point for the mix problem
        s       - the size of subset we want to choose, also equals to the summation of all elements of x
        Gamma1  - symmtric diagonal scaling paramter for DDFact
        Gamma2  - symmtric diagonal scaling paramter for DDFact comp
        alpha   - mixing parameter
    
        Output:
        fval    - objective value at current point x
        dx      - the gradient of the obejctive function at x
        %}
        mix_DDFact_DDFact_comp_obj_knitro_inline;
        end

        function [fval,dx,info] = mix_DDFact_DDFact_comp(obj,s,Gamma1,Gamma2)
        % mixing DDFact and DDFact comp bound (with optimal mixing parameter)
        %{
        Input:
        s       - the size of subset we want to choose, also equals to the summation of all elements of x
        Gamma1  - symmtric diagonal scaling parameter for DDFact
        Gamma2  - symmtric diagonal scaling parameter for DDFact comp

        Output:
        fval    - objective value at current point x
        dx      - the gradient of the obejctive function at x
        info    - struct containing necesssary information
        %}
        mix_DDFact_DDFact_comp_inline;
        end
    end 

    %% alternating optimizing mixing parameter and scaling parameter for mixing bound
    methods
        function [fval, x, info] = mix_alteritr(s, mix_pattern)
        %{
        Input:
        s       - the size of subset we want to choose, also equals to the summation of all elements of x
        mix_pattern
                - "DDFact_Linx"
                - "DDFact_comp_Linx"
                - "DDFact_DDFact_comp"
        Output:
        fval    - objective value at current point x
        dx      - the gradient of the obejctive function at x
        info    - struct containing necesssary information
        %}
        mix_alteritr_inline;
        end
    end
end