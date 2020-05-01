//Morse oscillator. Molecule: HF
clear
hbar = 1;                                       //Dirac's constant
omegae = 0.513086529995;                           //vibrational constant unit in ev
mu = 0.95;                                      //reduced mass for HF
De=5.8;                                   //Dissociation energy in ev
alfa = (mu*omegae^2/(2*De))^(1/2);
re = 0.91680;                                   //equilibrium distance
//create grid
maxX = 10;                                    //bounds of the problem
n = 600;                                         //define number of points
x = linspace(-maxX,maxX,n);
x = x';
delta = x(2)-x(1);

//Kinetic energy
for i=1:n
    T(i,i)= -2;                                 //T matrix has -2 on diagonals
    if i~=1 & i ~= n
        T(i,i-1)=1;                             //of diagonal elements, followed by a few special cases
        T(i,i+1)=1;
    end     //if
    if i==1                                     //Far topleft corner
        T(i,i+1)=1;
    end     //if
    if i==n                                     //Far bottomright corner
        T(i,i-1)=1;
    end     //if
end     //for
T = T*(-hbar^2*0.5/(mu*delta^2));               //need to multiply by this constant to get correct T

//potential
for i=1:n
    V(i,i) = De*(1-exp(-alfa*(x(i)-re)))^2;
end     //for

H = T+V;                                        //Hamiltonian
[u,E]=spec(H);                                  //u = eigenvectors,     E = eigenvalues
e = diag(E);                                    //vector of eigenvalues

plot(x,-u(:,1:4))                                //plot eigenvectors 1 to 4
title('Morse potential for HF')
xlabel('displacement')
ylabel('amplitude')
transit = e(2)-e(1);                            //transition energy between level 1 and 2, unit: eV

transcm = transit*8065.510204;                      //transition energy in cm-1, 1eV = 8065.51cm-1
