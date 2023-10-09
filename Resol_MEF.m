% Résolution d'un écoulement de Poiseuille.

clear all;
close all;
clc;

%******************************************************
%********* Taille caractéristique de la section  ******
%******    et fichier contenant le maillage    ********
%******************************************************

Echelle = 1; % En m à priori.

% *** mode : 0 carré        ***
% ***        1 triangle     ***
% ***        2 cacaouette   ***

mode = 2;

switch mode

case {0}
% Si "Mailleur_Simple_..."
NomfichMaillage='triangulation_Rec'; % Nom du fichier de sauvegarde 
% du maillage
Longueur=Echelle; % En m à priori.
Largeur=Longueur; % En m à priori.
Nx=30; % Entrée commune aux rectangles et aux triangles équilatéraux.
Ny=15; % Uniquement pour les rectangles
%Mailleur_Simple_Poiseuille_Rect(Longueur,Largeur,Echelle,Nx,Ny,NomfichMaillage);

case{1}
NomfichMaillage='triangulation_Tri';
Longueur=Echelle; % En m à priori.
Largeur=Longueur; % En m à priori.
Nx=30; % Entrée commune aux rectangles et aux triangles équilatéraux.
Ny=15; % Uniquement pour les rectangles
%Mailleur_Simple_Poiseuille_TriEqui(Longueur,Echelle,Nx,NomfichMaillage);

case{2}
% Si Conversion depuis freefem++
% Le fichier ".msh" est réputé être en coordonnées réduites
NomfichMaillage_freefem='Cacaouette_Avec_Trous.msh';
NomfichMaillage='Cacaouette_Avec_Trous';
%Converti_Mailleur_Freefem(NomfichMaillage_freefem,NomfichMaillage,Echelle);

end

% pause; % Pour examiner, sauvegarder ou imprimer le tracé du maillage.
% close all;

%******************************************************
%      Rechargement du fichier de maillage            *
%******************************************************

% Nomfich = fichier contenant les données crées par le mailleur.
NomfichMaillage=strcat(NomfichMaillage,'.mat');
if (~exist(NomfichMaillage,'file'))
    error('Le maillage n"existe pas.');
end

% On charge le maillage et les diverses données utiles.
load(NomfichMaillage);

a = 1; % Chute de pression
mu = 1; % viscosité constante du modèle

%******************************************************
%            Tableau des matrices locales             *
%******************************************************

% Initialisation/déclaration des cellules de "tab".
tab = cell(2,NT); 

for k=1:NT
    tab{1,k}=0.0; % surface du triangle
    tab{2,k}=zeros(3,3); % Matrice local K_loc
end

% Calcul et affectation dans "tab" des matrices locales 
% K_loc et des surfaces de chaque triangle.

for k=1:NT

    n1 = Num(1,k);
    x1 = Som(1,n1);
    y1 = Som(2,n1);

    n2 = Num(2,k);
    x2 = Som(1,n2);
    y2 = Som(2,n2);

    n3 = Num(3,k);
    x3 = Som(1,n3);
    y3 = Som(2,n3);

    tab{1,k} = (1/2)*abs((x2-x1)*(y3-y1)-(x3-x1)*(y2-y1));

    A = tab{1,k};

    tab{2,k}(1,1) = (mu/(4*A))*((x3-x2)^2+(y3-y2)^2);

    tab{2,k}(1,2) = (mu/(4*A))*(-(x3-x1)*(x3-x2)-(y3-y1)*(y3-y2));

    tab{2,k}(1,3) = (mu/(4*A))*((x2-x1)*(x3-x2)+(y2-y1)*(y3-y2));


    tab{2,k}(2,1) = tab{2,k}(1,2);

    tab{2,k}(2,2) = (mu/(4*A))*((x3-x1)^2+(y3-y1)^2);

    tab{2,k}(2,3) = (mu/(4*A))*(-(x2-x1)*(x3-x1)-(y2-y1)*(y3-y1));


    tab{2,k}(3,1) = tab{2,k}(1,3);

    tab{2,k}(3,2) = tab{2,k}(2,3);

    tab{2,k}(3,3) = (mu/(4*A))*((x2-x1)^2+(y2-y1)^2);

end




close all;

% ******************************************************
% Assemblage de la matrice de rigidité globale 
% et du vecteur de sollicitation global 
% ******************************************************

M_global = zeros(NS,NS);

F_global = zeros(NS,1);

for k=1:NT

    % Matrice de rigidité globale

    T = zeros(NS,NS);

    T(Num(1,k),Num(1,k)) = tab{2,k}(1,1);

    T(Num(1,k),Num(2,k)) = tab{2,k}(1,2);

    T(Num(1,k),Num(3,k)) = tab{2,k}(1,3);


    T(Num(2,k),Num(1,k)) = tab{2,k}(2,1);
    
    T(Num(2,k),Num(2,k)) = tab{2,k}(2,2);

    T(Num(2,k),Num(3,k)) = tab{2,k}(2,3);


    T(Num(3,k),Num(1,k)) = tab{2,k}(3,1);

    T(Num(3,k),Num(2,k)) = tab{2,k}(3,2);

    T(Num(3,k),Num(3,k)) = tab{2,k}(3,3);

    M_global = M_global + T;


    % Vecteur de sollicitation globale

    F = zeros(NS,1);

    F(Num(1,k)) = a*(1/3)*tab{1,k};
    F(Num(2,k)) = a*(1/3)*tab{1,k};
    F(Num(3,k)) = a*(1/3)*tab{1,k};

    F_global = F_global + F;
   
end


% ******************************************************
% ************ Conditions aux limites   ****************
% ******************************************************

for i=N+1:NS

    M_global(i,:) = 0;

    M_global(i,i) = 1;

    F_global(i) = 0;

end

% ******************************************************
% **     Calcul de la vitesse en résolvant        ******
% **            le système linéaire                  ***
% ******************************************************


U = inv(M_global)*F_global;

% Uxy = zeros(3,NS);
% 
% for i = 1:NS
% 
%     Uxy(1,i) = Som(1,i);
%     Uxy(2,i) = Som(2,i);
%     Uxy(3,i) = U(i);
% 
% end



% ******************************************************
% ***      Affichage de la solution par "patches"    ***
% ******************************************************

% Lx=Longueur/Echelle;
% Ly=Largeur/Echelle;
% dx=Lx/Nx;
% dy=Ly/Ny;
% 
% xi=[0:dx:1];
% yi=[0:dy:1];
% 
% Umat = reshape(Uxy(3,:),Nx+1,Ny+1);
% surf(xi,yi,Umat');

figure(1)

for k=1:NT

    n1 = Num(1,k);
    x1 = Som(1,n1);
    y1 = Som(2,n1);

    n2 = Num(2,k);
    x2 = Som(1,n2);
    y2 = Som(2,n2);

    n3 = Num(3,k);
    x3 = Som(1,n3);
    y3 = Som(2,n3);

    hold on

    patch([x1 x2 x3 x1],[y1 y2 y3 y1],[U(n1) U(n2) U(n3) U(n1)])

end

% ******************************************************
% ***      Affichage de la solution théorique    ***
% ******************************************************

Uth=zeros(NS,1);

for i=1:N

    x=Som(1,i);
    y=Som(2,i);

    Uth(i,1)=-a*sqrt(3)/(2*mu*Longueur)*y*((x-Longueur/2)-y/sqrt(3)+Longueur/2)*((x-Longueur/2)-y/sqrt(3)-Longueur/2);
end
 
Ugap = abs(U - Uth);

figure(2) 

for k=1:NT

    n1 = Num(1,k);
    x1 = Som(1,n1);
    y1 = Som(2,n1);

    n2 = Num(2,k);
    x2 = Som(1,n2);
    y2 = Som(2,n2);

    n3 = Num(3,k);
    x3 = Som(1,n3);
    y3 = Som(2,n3);

    hold on

    patch([x1 x2 x3 x1],[y1 y2 y3 y1],[Ugap(n1) Ugap(n2) Ugap(n3) Ugap(n1)])

    %patch([x1 x2 x3 x1],[y1 y2 y3 y1],[Uth(n1) Uth(n2) Uth(n3) Uth(n1)])
    
end

