clear
//create grid
maxX = 5;                      //bounds of the problem
n =500;                         //define number of points
x = linspace(-maxX,maxX,n);
x = x';
delta = x(2)-x(1);

//Kinetic energy
for i=1:n
    T(i,i)= -2;             //T matrix has -2 on diagonals
    if i~=1 & i ~= n
        T(i,i-1)=1;         //of diagonal elements, followed by a few special cases
        T(i,i+1)=1;
    end
    if i==1                 //Far topleft corner
        T(i,i+1)=1;
    end
    if i==n           //Far bottomright corner
        T(i,i-1)=1;
    end
end
T = T*(-0.5/delta^2);           //need to multiply by this constant to get correct T

//potential energy  (V=1/2 x^2)
for i=1:n
    V(i,i) = 1/2*(x(i))^2;
end

H = T+V;            //Hamiltonian
[u,E]=spec(H);      //u = eigenvectors,     E = eigenvalues

//precision test
z=diag(E);
e=100*abs(3.5-z(4));   //error at the fourth eigenvalue, exact solution is 3.5

plot(x,u(:,1:4))        //plot eigenvectors 1 to 4
title("Harmonic oscillator levels 1 to 4")
xlabel("displacement")
ylabel("amplitude")

for(i=1:n)
    exact(i)=(i-1)+1/2;
end
err=exact-z;
relerr=err*100;
