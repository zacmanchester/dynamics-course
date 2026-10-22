J = diag([1; 1.5; 2]);

Q0 = eye(3);
q0 = [1; 0; 0; 0];
w0 = randn(3,1);
x0_1 = [Q0(:); w0];
x0_2 = [q0; w0];

[tsim1, xsim1] = ode45(@(t,x)Qode(t,x,J), [0 60], x0_1);
[tsim2, xsim2] = ode45(@(t,x)quatode(t,x,J), [0 60], x0_2);


Qtraj = zeros(3,3,length(tsim1));
orth_err = zeros(length(tsim1),1);
for k=1:length(tsim1)
    Qtraj(:,:,k) = reshape(xsim1(k,1:9),3,3);
    %[U,S,V] = svd(Qtraj(:,:,k));
    %Qtraj(:,:,k) = U*V';
    orth_err(k) = trace(Qtraj(:,:,k)'*Qtraj(:,:,k) - eye(3));
end

figure(1)
plot(tsim1,orth_err);

qtraj = zeros(4,length(tsim2));
norm_err = zeros(length(tsim2),1);
for k=1:length(tsim2)
    qtraj(:,k) = xsim2(k,1:4);
    %qtraj(:,k) = qtraj(:,k)/norm(qtraj(:,k));
    norm_err(k) = qtraj(:,k)'*qtraj(:,k) - 1;
end

figure(2)
plot(tsim2,norm_err);