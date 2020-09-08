function obj = get_hermite_poly_mat(a,b)
    u = zeros(length(a),b);
    u(:,1) = ones(size(a));
    u(:,2) = a;
    u(:,3) = a.*a-1;
    u(:,4) = a.^3-3*a;
    u(:,5) = a.^4-6*a.^2+3;
    u(:,6) = a.^5-10*a.^3+15*a;
    u(:,7) = a.^6-15*a.^4+45*a.^2-15;
    obj = u(:,1:b) ; 

