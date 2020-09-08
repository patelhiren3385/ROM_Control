function obj = get_actuator_matrices(EI, rho, e31, b, x1,x2, ra)
K = @(D,x1,x2)[ -(12*D)/(x1 - x2)^3,  (6*D)/(x1 - x2)^2,  (12*D)/(x1 - x2)^3,  (6*D)/(x1 - x2)^2;
    (6*D)/(x1 - x2)^2,   -(4*D)/(x1 - x2),  -(6*D)/(x1 - x2)^2,   -(2*D)/(x1 - x2);
    (12*D)/(x1 - x2)^3, -(6*D)/(x1 - x2)^2, -(12*D)/(x1 - x2)^3, -(6*D)/(x1 - x2)^2;
    (6*D)/(x1 - x2)^2,   -(2*D)/(x1 - x2),  -(6*D)/(x1 - x2)^2,   -(4*D)/(x1 - x2)];
M = @(pA,x1,x2)[-(13*pA*(x1 - x2))/35, (11*pA*(x1 - x2)^2)/210,     -(9*pA*(x1 - x2))/70, -(13*pA*(x1 - x2)^2)/420;
    (11*pA*(x1 - x2)^2)/210,   -(pA*(x1 - x2)^3)/105,  (13*pA*(x1 - x2)^2)/420,     (pA*(x1 - x2)^3)/140;
    -(9*pA*(x1 - x2))/70, (13*pA*(x1 - x2)^2)/420,    -(13*pA*(x1 - x2))/35, -(11*pA*(x1 - x2)^2)/210;
-(13*pA*(x1 - x2)^2)/420,    (pA*(x1 - x2)^3)/140, -(11*pA*(x1 - x2)^2)/210,    -(pA*(x1 - x2)^3)/105];

P = @(e31, b,z)[0.0;e31.*b.*z;0.0;-e31.*b.*z];

obj.stiffness = K(EI,x1,x2);

obj.mass = M(rho,x1,x2);

obj.force = P(e31,b,ra);
end