function obj = get_distribution_pce(sigma,mu,n,itr) 

obj.delta = sigma.*n(itr) + mu ;
obj.pd = makedist('Normal','mu',mu,'sigma',sigma);
obj.PDF = pdf(obj.pd,obj.delta) ;
