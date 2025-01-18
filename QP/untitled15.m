H = eye(2);
f = [1;1]; 
A = [1 2; 1 -1];
bupper = [1; 2; 3; 4];
blower = [-1; -2; -3; -4];
sense = zeros(4,1,'int32');

[x,fval,exitflag,info] = daqp.quadprog(H,f,A,bupper,blower,sense);

%%
d = daqp();
d.setup(H,f,A,bupper,blower,sense);
[x,fval,exitflag,info] = d.solve();

d.settings('iter_limit',2000)