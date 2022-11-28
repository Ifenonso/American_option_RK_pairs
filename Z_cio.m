% 1. Dormand and Prince

clc
clear all
close all

format long
rd = 0.1; di = 0;      % risk free interest for each regime
sigma = 0.3;  % volatility for each regime

K = 100;    % asset price at expiration time
T = 1;    % expiration time
Tol = 10^-4;

h = 0.02; dt = 0.02; x_max = 3; 

nx = x_max/h; nt = T/dt;   % rectangular domain of x and t (mesh)
m_r = dt/(h*h);   % mesh ratio
                     % tolerance
                     
sfm(1) = K; % allocating first column to be the initial
                        % optimal exercise boundary
tic                       
for i = 1:nx+1
    if i==1
        U_old(i) = 0;  
        W_old(i) = 0;  
        
    else
        U_old(i) = 0;  
        W_old(i) = 0;
    end
end
                    
for i=1:nx
    if i == 1
        A(1,1) = -24/h^2-24/h; A(1,2) = 24/h^2;
        M(1,1) = 7; M(1,2) = 6; M(1,3) = -1;
    elseif i == nx
        A(i,nx-1) = 12/h^2; A(i,nx) = -24/h^2;
        M(i,nx-1) = 1; M(i,nx) = 10;
    else
        A(i,i-1) = 12/h^2;
        A(i,i) = -24/h^2;
        A(i,i+1) = 12/h^2;
        M(i,i-1) = 1;
        M(i,i) = 10;
        M(i,i+1) = 1;
    end
end

for i=1:nx-1
    if i == 1
        B(1,1) = -24/h^2; B(1,2) = 12/h^2;
        N(1,1) = 14; N(1,2) = -5; N(1,3) = 4; N(1,4) = -1;
    elseif i == nx-1
        B(i,nx-2) = 12/h^2; B(i,nx-1) = -24/h^2;
        N(i,nx-1) = 14; N(i,nx-2) = -5; N(i,nx-3) = 4; N(i,nx-4) = -1;
    else
        B(i,i-1) = 12/h^2;
        B(i,i) = -24/h^2;
        B(i,i+1) = 12/h^2;
        N(i,i-1) = 1;
        N(i,i) = 10;
        N(i,i+1) = 1;
    end
end

Hu1(1) = K; Hu1(nx) = 0; Hu1 = Hu1';

mb1 = 0.5*sigma^2; mc1 = rd; ah = 2*h; nu = rd-di-mb1;
kk = 1; t = 0; rr = 0; dta = dt;

m1 = 4980/27; m2 = 600/9; m3 = 32/3; 
n1 = 256; n2 = -48; n3 = 256/27; n4 = -1;

tic
while t<T
    
    % begin.................................................
    nn = 0; 
    if t+dt>T
      dt = T-t;  % ensure that the t+dt<=T
    end
    % end.................................................
    
    while(1)
        
        % begin.................................................
        Uh = U_old(3); Uh2 = U_old(5); Uh3 = U_old(7); Uh4 = U_old(9);    
        aa = di*sfm(kk); Qx = sqrt(rd*K-aa)/sigma;   
    
        a1 = (m3*2*Qx*ah^3)/(3*sigma^4);    
        b1 = (-2*m2*Qx*ah^2)/(3*sigma^2)+(4*m3*nu*Qx*ah^3)/(3*sigma^4)-...
           ((aa*m3*(ah^3))/(3*Qx*sigma^2));   
        gg = (ah*m1*Qx)+(m3*ah^3)*(((2*Qx*nu*nu)/(3*sigma^4))+(aa*nu)/...
           (6*Qx*sigma^4)-((aa^2)/(12*Qx*Qx*Qx*sigma^4))+((rd*Qx)/(2*...
           sigma^2))-(aa/(4*Qx*sigma^2)))-(m2*ah^2)*(((2*Qx*nu)/(3*...
           sigma^2))+(aa/(3*Qx*sigma^2)));
    
        ch = -n1*sqrt(Uh-K+exp(ah)*sfm(kk))-n2*sqrt(Uh2-K+exp(2*ah)*...
           sfm(kk))-n3*sqrt(Uh3-K+exp(3*ah)*sfm(kk))-n4*sqrt(Uh4-K+...
           exp(4*ah)*sfm(kk))+gg;
    
        s1 = ((-b1-sqrt(b1^2-4*a1*ch))/(2*a1))+nu; 
     
        Hw1(1) = -sfm(kk);  Hw1(nx-1) = 0; Huu1(1) = K-sfm(kk); 
        Huu1(nx-1) = 0; Hw1 = Hw1'; Huu1 = Huu1';
                          
        Ru1 = (mb1)*(M\(A*U_old(1:nx)'+(24/h)*Hu1))+s1*W_old(1:nx)'-...
           mc1*U_old(1:nx)'; 
        Ua1 = U_old(1:nx)'+(dt/5)*Ru1; 
        sfa = K-Ua1(1);
            
        Rw1 = (mb1)*(N\(B*W_old(2:nx)'+(12/h^2)*Hw1))+s1*(N\(B*U_old(2:nx)'...
            +(12/h^2)*Huu1))-mc1*W_old(2:nx)';
        Wa1 = W_old(2:nx)'+(dt/5)*Rw1; 
        % end.................................................
    
        % begin.................................................
        Uh = Ua1(3); Uh2 = Ua1(5); Uh3 = Ua1(7); Uh4 = Ua1(9);   
        aa = di*sfa; Qx = sqrt(rd*K-aa)/sigma;   
    
        a1 = (m3*2*Qx*ah^3)/(3*sigma^4); 
        b1 = (-2*m2*Qx*ah^2)/(3*sigma^2)+(4*m3*nu*Qx*ah^3)/(3*sigma^4)-...
           ((aa*m3*(ah^3))/(3*Qx*sigma^2));  
        gg = (ah*m1*Qx)+(m3*ah^3)*(((2*Qx*nu*nu)/(3*sigma^4))+(aa*nu)/...
           (6*Qx*sigma^4)-((aa^2)/(12*Qx*Qx*Qx*sigma^4))+((rd*Qx)/(2*...
            sigma^2))-(aa/(4*Qx*sigma^2)))-(m2*ah^2)*(((2*Qx*nu)/(3*...
            sigma^2))+(aa/(3*Qx*sigma^2)));
    
        ch = -n1*sqrt(Uh-K+exp(ah)*sfa)-n2*sqrt(Uh2-K+exp(2*ah)*sfa)-...
            n3*sqrt(Uh3-K+exp(3*ah)*sfa)-n4*sqrt(Uh4-K+exp(4*ah)*sfa)+gg;
    
        s1 = ((-b1-sqrt(b1^2-4*a1*ch))/(2*a1))+nu; 
  
        Iw1(1) = -sfa;  Iw1(nx-1) = 0; Iuu1(1) = K-sfa; Iuu1(nx-1) = 0;
        Iw1 = Iw1'; Iuu1 = Iuu1';
    
        Su1 = (mb1)*(M\(A*Ua1+(24/h)*Hu1))+s1*[-sfa;Wa1]-mc1*Ua1;
        Ub1 = U_old(1:nx)'+dt*((3/40)*Ru1+(9/40)*Su1); 
        sfb = K-Ub1(1);
    
        Sw1 = (mb1)*(N\(B*Wa1+(12/h^2)*Iw1))+s1*(N\(B*Ua1(2:nx)+...
            (12/h^2)*Iuu1))-mc1*Wa1;
        Wb1 = W_old(2:nx)'+dt*((3/40)*Rw1+(9/40)*Sw1);
        % end.................................................
    
        % begin.................................................
        Uh = Ub1(3); Uh2 = Ub1(5); Uh3 = Ub1(7); Uh4 = Ub1(9);    
        aa = di*sfb; Qx = sqrt(rd*K-aa)/sigma;   
    
        a1 = (m3*2*Qx*ah^3)/(3*sigma^4);
        b1 = (-2*m2*Qx*ah^2)/(3*sigma^2)+(4*m3*nu*Qx*ah^3)/(3*sigma^4)-...
            ((aa*m3*(ah^3))/(3*Qx*sigma^2));   
        gg = (ah*m1*Qx)+(m3*ah^3)*(((2*Qx*nu*nu)/(3*sigma^4))+(aa*nu)/...
            (6*Qx*sigma^4)-((aa^2)/(12*Qx*Qx*Qx*sigma^4))+((rd*Qx)/(2*...
             sigma^2))-(aa/(4*Qx*sigma^2)))-(m2*ah^2)*(((2*Qx*nu)/(3*...
             sigma^2))+(aa/(3*Qx*sigma^2)));
    
        ch = -n1*sqrt(Uh-K+exp(ah)*sfb)-n2*sqrt(Uh2-K+exp(2*ah)*sfb)-n3*...
            sqrt(Uh3-K+exp(3*ah)*sfb)-n4*sqrt(Uh4-K+exp(4*ah)*sfb)+gg;
    
        s1 = ((-b1-sqrt(b1^2-4*a1*ch))/(2*a1))+nu; 
    
        Jw1(1) = -sfb;  Jw1(nx-1) = 0; Juu1(1) = K-sfb; Juu1(nx-1) = 0;
        Jw1 = Jw1'; Juu1 = Juu1'; 
                         
        Tu1 = (mb1)*(M\(A*Ub1+(24/h)*Hu1))+s1*[-sfb;Wb1]-mc1*Ub1;    
        Uc1 = U_old(1:nx)'+dt*((44/45)*Ru1-(56/15)*Su1+(32/9)*Tu1);
        sfc = K-Uc1(1);
    
        Tw1 = (mb1)*(N\(B*Wb1+(12/h^2)*Jw1))+s1*(N\(B*Ub1(2:nx)+...
            (12/h^2)*Juu1))-mc1*Wb1;    
        Wc1 = W_old(2:nx)'+dt*((44/45)*Rw1-(56/15)*Sw1+(32/9)*Tw1); 
        % end.................................................
    
        % begin.................................................
        Uh = Uc1(3); Uh2 = Uc1(5); Uh3 = Uc1(7); Uh4 = Uc1(9);   
        aa = di*sfc; Qx = sqrt(rd*K-aa)/sigma;   
        a1 = (m3*2*Qx*ah^3)/(3*sigma^4);
        b1 = (-2*m2*Qx*ah^2)/(3*sigma^2)+(4*m3*nu*Qx*ah^3)/(3*sigma^4)-...
           ((aa*m3*(ah^3))/(3*Qx*sigma^2));    
        gg = (ah*m1*Qx)+(m3*ah^3)*(((2*Qx*nu*nu)/(3*sigma^4))+(aa*nu)/...
          (6*Qx*sigma^4)-((aa^2)/(12*Qx*Qx*Qx*sigma^4))+((rd*Qx)/(2*...
           sigma^2))-(aa/(4*Qx*sigma^2)))-(m2*ah^2)*(((2*Qx*nu)/(3*...
           sigma^2))+(aa/(3*Qx*sigma^2)));
    
        ch = -n1*sqrt(Uh-K+exp(ah)*sfc)-n2*sqrt(Uh2-K+exp(2*ah)*sfc)-...
            n3*sqrt(Uh3-K+exp(3*ah)*sfc)-n4*sqrt(Uh4-K+exp(4*ah)*sfc)+gg;
    
        s1 = ((-b1-sqrt(b1^2-4*a1*ch))/(2*a1))+nu;
    
        Kw1(1) = -sfc;  Kw1(nx-1) = 0; Kuu1(1) = K-sfc; Kuu1(nx-1) = 0;
        Kw1 = Kw1'; Kuu1 = Kuu1';
    
        Tu2 = (mb1)*(M\(A*Uc1+(24/h)*Hu1))+s1*[-sfc;Wc1]-mc1*Uc1;
        Ud1 = U_old(1:nx)'+dt*((19372/6561)*Ru1-(25360/2187)*Su1+...
            (64448/6561)*Tu1-(212/729)*Tu2);
        sfd = K-Ud1(1);
    
        Tw2 = (mb1)*(N\(B*Wc1+(12/h^2)*Kw1))+s1*(N\(B*Uc1(2:nx)+...
            (12/h^2)*Kuu1))-mc1*Wc1;    
        Wd1 = W_old(2:nx)'+dt*((19372/6561)*Rw1-(25360/2187)*Sw1+...
            (64448/6561)*Tw1-(212/729)*Tw2);
        % end.................................................
    
    
        % begin.................................................
        Uh = Ud1(3); Uh2 = Ud1(5); Uh3 = Ud1(7); Uh4 = Ud1(9);   
        aa = di*sfd; Qx = sqrt(rd*K-aa)/sigma;   
        a1 = (m3*2*Qx*ah^3)/(3*sigma^4);
        b1 = (-2*m2*Qx*ah^2)/(3*sigma^2)+(4*m3*nu*Qx*ah^3)/(3*sigma^4)-...
            ((aa*m3*(ah^3))/(3*Qx*sigma^2));    
        gg = (ah*m1*Qx)+(m3*ah^3)*(((2*Qx*nu*nu)/(3*sigma^4))+(aa*nu)/...
            (6*Qx*sigma^4)-((aa^2)/(12*Qx*Qx*Qx*sigma^4))+((rd*Qx)/(2*...
             sigma^2))-(aa/(4*Qx*sigma^2)))-(m2*ah^2)*(((2*Qx*nu)/(3*...
             sigma^2))+(aa/(3*Qx*sigma^2)));
    
        ch = -n1*sqrt(Uh-K+exp(ah)*sfd)-n2*sqrt(Uh2-K+exp(2*ah)*sfd)-n3*...
            sqrt(Uh3-K+exp(3*ah)*sfd)-n4*sqrt(Uh4-K+exp(4*ah)*sfd)+gg;
    
        s1 = ((-b1-sqrt(b1^2-4*a1*ch))/(2*a1))+nu;
    
        Lw1(1) = -sfd;  Lw1(nx-1) = 0; Luu1(1) = K-sfd; Luu1(nx-1) = 0;
        Lw1 = Lw1'; Luu1 = Luu1';
    
        Tu3 = (mb1)*(M\(A*Ud1+(24/h)*Hu1))+s1*[-sfd;Wd1]-mc1*Ud1;
        Ue1 = U_old(1:nx)'+dt*((9017/3168)*Ru1-(355/33)*Su1+...
            (46732/5247)*Tu1+(49/176)*Tu2-(5103/18656)*Tu3);
        sfe = K-Ue1(1);
    
        Tw3 = (mb1)*(N\(B*Wd1+(12/h^2)*Lw1))+s1*(N\(B*Ud1(2:nx)+(12/h^2)...
            *Luu1))-mc1*Wd1;    
        We1 = W_old(2:nx)'+dt*((9017/3168)*Rw1-(355/33)*Sw1+...
            (46732/5247)*Tw1+(49/176)*Tw2-(5103/18656)*Tw3);
        % end.................................................
        
        % begin.................................................
        Uh = Ue1(3); Uh2 = Ue1(5); Uh3 = Ue1(7); Uh4 = Ue1(9);   
        aa = di*sfe; Qx = sqrt(rd*K-aa)/sigma;   
        a1 = (m3*2*Qx*ah^3)/(3*sigma^4);
        b1 = (-2*m2*Qx*ah^2)/(3*sigma^2)+(4*m3*nu*Qx*ah^3)/(3*sigma^4)-...
           ((aa*m3*(ah^3))/(3*Qx*sigma^2));    
        gg = (ah*m1*Qx)+(m3*ah^3)*(((2*Qx*nu*nu)/(3*sigma^4))+(aa*nu)/...
           (6*Qx*sigma^4)-((aa^2)/(12*Qx*Qx*Qx*sigma^4))+((rd*Qx)/(2*...
           sigma^2))-(aa/(4*Qx*sigma^2)))-(m2*ah^2)*(((2*Qx*nu)/(3*...
           sigma^2))+(aa/(3*Qx*sigma^2)));
    
        ch = -n1*sqrt(Uh-K+exp(ah)*sfe)-n2*sqrt(Uh2-K+exp(2*ah)*sfe)-n3*...
           sqrt(Uh3-K+exp(3*ah)*sfe)-n4*sqrt(Uh4-K+exp(4*ah)*sfe)+gg;
    
        s1 = ((-b1-sqrt(b1^2-4*a1*ch))/(2*a1))+nu;
    
        Mw1(1) = -sfe;  Mw1(nx-1) = 0; Muu1(1) = K-sfe; Muu1(nx-1) = 0;
        Mw1 = Mw1'; Muu1 = Muu1';
    
        Tu4 = (mb1)*(M\(A*Ue1+(24/h)*Hu1))+s1*[-sfe;We1]-mc1*Ue1;
        Uf1 = U_old(1:nx)'+dt*((35/384)*Ru1+(500/1113)*Tu1+(125/192)...
            *Tu2-(2187/6784)*Tu3+(11/84)*Tu4);
        sff = K-Uf1(1);
    
        Tw4 = (mb1)*(N\(B*We1+(12/h^2)*Mw1))+s1*(N\(B*Ue1(2:nx)+(12/h^2)...
            *Muu1))-mc1*We1;    
        Wf1 = W_old(2:nx)'+dt*((35/384)*Rw1+(500/1113)*Tw1+(125/192)...
            *Tw2-(2187/6784)*Tw3+(11/84)*Tw4);
               
        ab = isreal(sum(Uf1)); ac = isreal(sum(Wf1));
        if ab == 1 && ac == 1
            break;
        else
            dt = 0.5*dt;
            Huu1 = Huu1';  Hw1 = Hw1'; Iuu1 = Iuu1';  Iw1 = Iw1';
            Juu1 = Juu1';  Jw1 = Jw1'; Kuu1 = Kuu1';  Kw1 = Kw1';
            Luu1 = Luu1';  Lw1 = Lw1'; Muu1 = Muu1';  Mw1 = Mw1';
        end       
    % end.................................................
    end    
    
    % begin.................................................
    Uh = Uf1(3); Uh2 = Uf1(5); Uh3 = Uf1(7); Uh4 = Uf1(9);   
    aa = di*sff; Qx = sqrt(rd*K-aa)/sigma;   
    a1 = (m3*2*Qx*ah^3)/(3*sigma^4);
    b1 = (-2*m2*Qx*ah^2)/(3*sigma^2)+(4*m3*nu*Qx*ah^3)/(3*sigma^4)-...
       ((aa*m3*(ah^3))/(3*Qx*sigma^2));    
    gg = (ah*m1*Qx)+(m3*ah^3)*(((2*Qx*nu*nu)/(3*sigma^4))+(aa*nu)/...
       (6*Qx*sigma^4)-((aa^2)/(12*Qx*Qx*Qx*sigma^4))+((rd*Qx)/(2*...
       sigma^2))-(aa/(4*Qx*sigma^2)))-(m2*ah^2)*(((2*Qx*nu)/(3*...
       sigma^2))+(aa/(3*Qx*sigma^2)));
    
    ch = -n1*sqrt(Uh-K+exp(ah)*sff)-n2*sqrt(Uh2-K+exp(2*ah)*sff)-n3*...
       sqrt(Uh3-K+exp(3*ah)*sff)-n4*sqrt(Uh4-K+exp(4*ah)*sff)+gg;
    
    s1 = ((-b1-sqrt(b1^2-4*a1*ch))/(2*a1))+nu;
    
    Nw1(1) = -sff;  Nw1(nx-1) = 0; Nuu1(1) = K-sff; Nuu1(nx-1) = 0;
    Nw1 = Nw1'; Nuu1 = Nuu1';
    
    Vu1 = (mb1)*(M\(A*Uf1+(24/h)*Hu1))+s1*[-sff;Wf1]-mc1*Uf1;
    Un1 = Uf1;
    
    Un = U_old(1:nx)'+dt*((5179/57600)*Ru1+(7571/16695)*Tu1+...
        (393/640)*Tu2-(92097/339200)*Tu3+(187/2100)*Tu4+(1/40)*Vu1); 
    
    Vw1 = (mb1)*(N\(B*Wf1+(12/h^2)*Nw1))+s1*(N\(B*Uf1(2:nx)+(12/h^2)...
       *Nuu1))-mc1*Wf1;    
    Wn1 = Wf1;  
   
    Wn = W_old(2:nx)'+dt*((5179/57600)*Rw1+(7571/16695)*Tw1+...
        (393/640)*Tw2-(92097/339200)*Tw3+(187/2100)*Tw4+(1/40)*Vw1); 
    
    Reu(1) = max((Un1-Un)); Reu(2) = max((Wn1-Wn));
    aw = abs(max(Reu));
    
   %end: compute the asset option and option Greeks and their error
   
   % begin: transfer values.
   if aw<=Tol
       dsfm(kk) = ((-b1-sqrt(b1^2-4*a1*ch))/(2*a1))*sfm(kk);
       dsf = ((-b1-sqrt(b1^2-4*a1*ch))/(2*a1))*sfm(kk);
       sfnew = K-Un1(1); U_iter1(nx+1) = 0; 
       for i = 1:nx
           U_iter1(i) = Un1(i);
       end
       
       W_iter1(1) = -sfnew; W_iter1(nx+1) = 0;  
       for i = 2:nx
          W_iter1(i) = Wn1(i-1);
       end
       sfm(kk+1) = sfnew; sf=sfnew  % Storing the exercise boundary
       fff(kk) = dt; kk = kk+1;
        
        for i = 1:nx+1
            U_old(i) = U_iter1(i); % assigning the converged asset, delta 
            W_old(i) = W_iter1(i); % delta and gamma option 
        end       
        Huu1 = Huu1';  Hw1 = Hw1'; Iuu1 = Iuu1';  Iw1 = Iw1';
        Juu1 = Juu1';  Jw1 = Jw1'; Kuu1 = Kuu1';  Kw1 = Kw1';
        Luu1 = Luu1';  Lw1 = Lw1'; Muu1 = Muu1';  Mw1 = Mw1';
        Nuu1 = Nuu1';  Nw1 = Nw1';
        t = t+dt;  
        dt = 0.9*dt*((Tol/abs(max(Reu)))^(0.25));
        % end: transfer values
    else
        % begin: select the optimal step size
        dt = 0.9*dt*((Tol/abs(max(Reu)))^(0.2));
        Huu1 = Huu1';  Hw1 = Hw1'; Iuu1 = Iuu1';  Iw1 = Iw1';
        Juu1 = Juu1';  Jw1 = Jw1'; Kuu1 = Kuu1';  Kw1 = Kw1';
        Luu1 = Luu1';  Lw1 = Lw1'; Muu1 = Muu1';  Mw1 = Mw1';
        Nuu1 = Nuu1';  Nw1 = Nw1';
        rr = rr+1;
    end    % change the column
end
toc

ay(1) = 0;
for i = 2:length(fff)+1
ay(i) = ay(i-1)+fff(i-1);
end

ayy(1) = 0;
for i = 2:length(fff)
ayy(i) = ay(i-1)+fff(i-1);
end

for i = 1:nx+1
    S1(i) = exp((i-1)*h)*sfnew; % first regime asset price calculatio
end

for i = 1:nx+2
    if i==1
        Smain1(i) = 0;  
        Smain3(i) = 0;
        Vmain1(i) = K - Smain1(i); 
        Vmain3(i) = K - Smain3(i);
    else
        Smain1(i) = S1(i-1); 
        Smain3(i) = K+(i-2);       % calculation
        Vmain1(i) = U_iter1(i-1);
        Vmain3(i) = 0;
    end
end

Wmain1(1) = -1;
for i = 2:nx+2
    Wmain1(i) = W_iter1(i-1)/S1(i-1); % first regime delta option 
end 

% start: interpolated data
Tableu(1) = interp1(Smain1,Vmain1,80); 
Tableu(2) = interp1(Smain1,Vmain1,90,'spline');  
Tableu(3) = interp1(Smain1,Vmain1,100,'spline');  
Tableu(4) = interp1(Smain1,Vmain1,110,'spline'); 
Tableu(5) = interp1(Smain1,Vmain1,120,'spline');     


Tablew(1) = interp1(Smain1,Wmain1,80); 
Tablew(2) = interp1(Smain1,Wmain1,90,'spline');  
Tablew(3) = interp1(Smain1,Wmain1,100,'spline');  
Tablew(4) = interp1(Smain1,Wmain1,110,'spline'); 
Tablew(5) = interp1(Smain1,Wmain1,120,'spline');     

% % start: visualization
figure(1)
subplot(1,2,1)
plot(ay,sfm,'LineWidth',2)
set(gca,'FontSize',12)
o=ylabel('Optimal Exercise Boundary','FontSize',14);
get(o)
p=xlabel('Time to Maturity','FontSize',14);
get(p)
grid on
axis([0,0.5,75,100])

subplot(1,2,2)
plot(ayy,dsfm,'LineWidth',2)
set(gca,'FontSize',12)
o=ylabel('Deriv. Opt. Exer. Bound.','FontSize',14);
get(o)
p=xlabel('Time to Maturity','FontSize',14);
get(p)
grid on
axis([0,0.5,-10,10])

figure(2)
subplot(1,2,1)
plot(S1,U_iter1,'b','LineWidth',2)
hold on
plot(Smain3,Vmain3,'k','LineWidth',2)
set(gca,'FontSize',12)
o=ylabel('Asset Option','FontSize',14);
get(o)
p=xlabel('Asset Price','FontSize',14);
get(p)
grid on
%axis([0,200,0,100])

subplot(1,2,2)
plot(Smain1,Wmain1,'LineWidth',2)
set(gca,'FontSize',12)
o=ylabel('Delta Option','FontSize',14);
get(o)
p=xlabel('Asset Price','FontSize',14);
get(p)
grid on
axis([0 200 -1 0])
