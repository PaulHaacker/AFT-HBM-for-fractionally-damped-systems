function f = vanderpol(t,x,nu)
% Van der Pol system

f = [x(2) ; 
     nu*(1-x(1)^2)*x(2) - x(1)  ];

