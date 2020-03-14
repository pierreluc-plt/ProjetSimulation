clc; clearvars;

% Dimensions d'appartement
Lx=2.4; %[m]
Ly=2.4;  %[m]

omega = 2.*pi.*10e3;

d_ar=[];tini_ar=[];tinv_ar=[];mem_ar=[];Tm_ar=[];
fact=1/4;
d=0.1*fact; %Pas de discrétisation en [m]
d_ar=[d_ar d];
Nx=round(Lx/d)+1;
Ny=round(Ly/d)+1;

alpha_b = 3.21e-4;
rho_b = 10e9;
alpha_e = 1.18;
rho_e = 2.15e9;

tic
% Initialisation de la source de chaleur et de la conductivité thermique
S=sparse(Ny,Nx); B=sparse(Ny,Nx); alpha=sparse(Ny,Nx); rho=sparse(Ny,Nx);
for i=1:Ny
    y=(i-1)*d;
    for j=1:Nx
        x=(j-1)*d;
        
                if (x>=1.4)&&(x<=1.4+d)&&(y>=1)&&(y<=1+d)
                    % À l'intérieur de l'élément chauffant
                    S(i,j)=1;
                else
                    % À l'extérieur de l'élément chauffant
                    S(i,j)=0;
                end
        
        
        % Délimitation des espaces
        if (x>=1.5)&&(x<=2)&&(y>=1)&&(y<=1.5)
            % Bois
            B(i,j) = 640.72;
            alpha(i,j) = alpha_b;
            rho(i,j) = rho_b;
        else
            % eau
            B(i,j) = 998.30;
            alpha(i,j) = alpha_e;
            rho(i,j) = rho_e;
        end

    end
end

M=sparse(Nx*Ny,Nx*Ny);
b=sparse(Nx*Ny,1);

ke = (rho_e.*(alpha_e + 1i.*omega));
kb = (rho_b.*(alpha_b + 1i.*omega));

for i=1:Ny
    y=(i-1)*d;
    for j=1:Nx
        x=(j-1)*d;
        % remplir la ligne pl de la matrice M
        pl=(i-1)*Nx+j;
        
        if ((i>1)&(i<Ny))&&((j>1)&(j<Nx)) ...
                &~((x==1.5) & ((y>=1)&(y<=(1.5)))) ...
                &~((x==(2)) & ((y>=1)&(y<=(1.5)))) ...
                &~((y==1) & ((x>=2)&(x<=(1.5)))) ...
                &~((y==(1.5)) & ((x>=2)&(x<=(2))))
            % noeud qui est strictement à l'intérieur de la cellule de simulation
            pc=pl;                      M(pl,pc)=(omega.^2./B(i,j).*d.^2.*rho(i,j) + 1i.*omega.*alpha(i,j).*d.^2.*rho(i,j) - 4); % contribution de noeud (i,j)
            pc=(i-1)*Nx+j-1;            M(pl,pc)=1; % contribution de noeud (i,j-1)
            pc=(i-1)*Nx+j+1;            M(pl,pc)=1; % contribution de noeud (i,j+1)
            pc=(i-2)*Nx+j;              M(pl,pc)=1; % contribution de noeud (i-1,j)
            pc=(i)*Nx+j;                M(pl,pc)=1; % contribution de noeud (i+1,j)
            b(pl)= -S(i,j).*rho(i,j).*d.^2;
        elseif (i==1)
            % noeud sur le plafond y=0
            pc=pl;                      M(pl,pc)=2; % contribution de noeud (1,j)
            pc=(i)*Nx+j;                M(pl,pc)=-5; % contribution de noeud (2,j)
            pc=(i+1)*Nx+j;              M(pl,pc)=4; % contribution de noeud (3,j)
            pc=(i+2)*Nx+j;              M(pl,pc)=-1; % contribution de noeud (4,j)
            b(pl)=0;
        elseif (i==Ny)
            % noeud sur le plancher y=Ly
            pc=pl;                      M(pl,pc)=2; % contribution de noeud (Ny,j)
            pc=(i-2)*Nx+j;              M(pl,pc)=-5; % contribution de noeud (Ny-1,j)
            pc=(i-3)*Nx+j;              M(pl,pc)=4; % contribution de noeud (Ny-2,j)
            pc=(i-4)*Nx+j;              M(pl,pc)=-1; % contribution de noeud (Ny-3,j)
            b(pl)=0;
        elseif (j==1)
            % noeud à la surface externe du mur x=0
            pc=pl;                      M(pl,pc)=2; % contribution de noeud (i,1)
            pc=(i-1)*Nx+j+1;            M(pl,pc)=-5; % contribution de noeud (i,2)
            pc=(i-1)*Nx+j+2;            M(pl,pc)=4; % contribution de noeud (i,3)
            pc=(i-1)*Nx+j+3;            M(pl,pc)=-1; % contribution de noeud (i,4)
            b(pl)=0;
        elseif (j==Nx)
            % noeud à la surface externe du mur x=Nx
            pc=pl;                      M(pl,pc)=2; % contribution de noeud (i,Nx)
            pc=(i-1)*Nx+j-1;            M(pl,pc)=-5; % contribution de noeud (i,Nx-1)
            pc=(i-1)*Nx+j-2;            M(pl,pc)=4; % contribution de noeud (i,Nx-2)
            pc=(i-1)*Nx+j-3;            M(pl,pc)=-1; % contribution de noeud (i,Nx-3)
            b(pl)=0;
            
            % Conditions entre milieu
        elseif (x==1.5) & ((y>=1)&(y<=(1.5)))
            %Le long de la surface entre mur et appartement gauche
            pc=(i-1)*Nx+j-2;            M(pl,pc)=kb; % contribution de noeud (i,M-2)
            pc=(i-1)*Nx+j-1;            M(pl,pc)=-4*kb; % contribution de noeud (i,M-1)
            pc=pl;                      M(pl,pc)=(3*kb+3*ke); % contribution de noeud (i,M)
            pc=(i-1)*Nx+j+1;            M(pl,pc)=-4*ke; % contribution de noeud (i,M+1)
            pc=(i-1)*Nx+j+2;            M(pl,pc)=ke; % contribution de noeud (i,M+2)
            b(pl)=0;
        elseif (x==(2)) & ((y>=1)&(y<=(1.5)))
            %Le long de la surface entre mur et appartement droite
            pc=(i-1)*Nx+j-2;            M(pl,pc)=ke; % contribution de noeud (i,M-2)
            pc=(i-1)*Nx+j-1;            M(pl,pc)=-4*ke; % contribution de noeud (i,M-1)
            pc=pl;                      M(pl,pc)=(3*kb+3*ke); % contribution de noeud (i,M)
            pc=(i-1)*Nx+j+1;            M(pl,pc)=-4*kb; % contribution de noeud (i,M+1)
            pc=(i-1)*Nx+j+2;            M(pl,pc)=kb; % contribution de noeud (i,M+2)
            b(pl)=0;
        elseif (y==1) & ((x>=1.5)&(x<=(2)))
            %Le long du plancher entre mur et appartement (sol)
            pc=(i-3)*Nx+j;              M(pl,pc)=kb; % contribution de noeud (M-2,j)
            pc=(i-2)*Nx+j;              M(pl,pc)=-4*kb; % contribution de noeud (M-1,j)
            pc=pl;                      M(pl,pc)=(3*kb+3*ke); % contribution de noeud (M,j)
            pc=(i)*Nx+j;                M(pl,pc)=-4*ke; % contribution de noeud (M+1,j)
            pc=(i+1)*Nx+j;              M(pl,pc)=ke; % contribution de noeud (M+2,j)
            b(pl)=0;
        elseif (y==(1.5)) & ((x>=1.5)&(x<=(2)))
            %Le long du plafond entre mur et appartement (haut)
            pc=(i-3)*Nx+j;              M(pl,pc)=ke; % contribution de noeud (M-2,j)
            pc=(i-2)*Nx+j;              M(pl,pc)=-4*ke; % contribution de noeud (M-1,j)
            pc=pl;                      M(pl,pc)=(3*kb+3*ke); % contribution de noeud (M,j)
            pc=(i)*Nx+j;                M(pl,pc)=-4*kb; % contribution de noeud (M+1,j)
            pc=(i+1)*Nx+j;              M(pl,pc)=kb; % contribution de noeud (M+2,j)
            b(pl)=0;
        else
            display('Erreur dans la définition de la matrice de coefficients');
        end
        
    end
end
tini_ar=[tini_ar toc]

tic
%T=M\b;
[L,U]=lu(M);T=U\(L\b);
tinv_ar=[tinv_ar toc]

mem_ar=[mem_ar 8*(Nx*Ny)^2]

Tr=reshape(T,Nx,Ny)';

Tm_ar=[Tm_ar Tr(round(Ly/d/2+1),round(Lx/d/2+1))]




figure(1)
h=pcolor((0:d:Lx),(0:d:Ly),abs(S));set(h,'LineStyle','none')
colorbar;
xlabel('x [m]'); ylabel('y [m]'); title('S(x,y) [W/m^3]')
axis equal
axis tight

figure(2)
h=pcolor((0:d:Lx),(0:d:Ly),B);set(h,'LineStyle','none')
colorbar;
xlabel('x [m]'); ylabel('y [m]'); title('k(x,y) [W/(m^2\cdotK)]')
axis equal
axis tight

figure(3)
h=pcolor((0:d:Lx),(0:d:Ly),log(abs(Tr)));set(h,'LineStyle','none')
colormap(hot);
colorbar;
xlabel('x [m]'); ylabel('y [m]'); title('T(x,y) [^oC]')
axis equal
axis tight


