function obj = get_distribution(sigma,mu,n) 

obj.delta = sigma.*randn(n,1) + mu ;
obj.pd = makedist('Normal','mu',mu,'sigma',sigma);
obj.PDF = pdf(obj.pd,obj.delta) ;
