J = diag([1; 1.5; 2]);

Q0 = eye(3);
w0 = randn(3,1);
x0 = [Q0(:); w0];

[tsim, xsim] = ode45(@(t,x)Qode(t,x,J), [0 60], x0);


Qtraj = zeros(3,3,length(xsim));
orth_err = zeros(length(xsim),1);
for k=1:length(xsim)
    Qtraj(:,:,k) = reshape(xsim(k,1:9),3,3);
    orth_err(k) = trace(Qtraj(:,:,k)'*Qtraj(:,:,k) - eye(3));
end

plot(tsim,orth_err);