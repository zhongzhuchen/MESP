classdef MESP
    properties
        C double
        size int16
        A double
        b double
        F double
        Fsquare double
    end

    %% intialization function
    methods  
        function obj = MESP(C, A, b)
            if nargin < 0
                error("Please input valid values.");
            end
            obj.C=C;
            obj.size=length(C);
            obj.A = A;
            obj.b = b;
            [obj.F,obj.Fsquare,~] = gen_data(C,0);
        end
    end

    %% methods for obtain lower bound
    methods
        function [lb,info] = obtain_lb(obj,s)
        % obtain the lower bound and corresponding x indices 
        obtain_lb_inline;
        end
    end

    %% methods block for factorization
    methods
        function [fval,dx,info] = DDFact_obj(obj,x,s)
        % This function calculate the objective value, gradient, and info of DDFact for given data
        %{
        Input:
        x       - current point for the DDFact problem
        s       - the size of subset we want to choose, also equals to the summation of all elements of x
        
        Output:
        fval    - objective value of DDFact at current point x
        dx      - the gradient of the obejctive function of DDFact at x
        info    - struct containing necesssary information
        %}
        DDFact_obj_inline;
        end
        
        function [fval,dx] = DDFact_obj_knitro(obj,x,s)
            % create a callback function for Knitro specifying objective value and gradient 
            [fval,dx,~] = obj.DDFact_obj(x,s);
            fval=-fval;
            dx=-dx;
        end

        function [fval,x,info] = Knitro_DDFact(obj,x0,s,comp)
        % calling knitro to solve the DDFact problem
        %{
        Input:
        s       - the size of subset we want to choose, also equals to the
                  summation of all elements of x0
        x0      - initial point
        comp    - indicate complementary bound
        Output:
        fval    - objective value of DDFact at optimal solution x
        x       - optimal solution x
        info    - struct containing necesssary information
        %}
        Knitro_DDFact_inline;
        end
    end
end