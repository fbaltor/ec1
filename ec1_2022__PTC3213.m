%%%% ///////////////////////////////////////////////////////////////////////////
%%%%        PTC3213 - EC1 - 2021 - DIFERENÇAS FINITAS
%%%%
%%%%    12552562 Felipe Keller Baltor
%%%%    NUSP2 Aluno B
%%%%    NUSP3 Aluno C
%%%%                           Complete os campos com ???????????????????????
%%%%
%%%%  /////////////////////////////////////////////////////////////////////////
clear;
clf;
%%%
%%%  A linha 15 (pkg install...) deve ser executada uma única vez no
%%% GNU Octave e necessita de conxão com a internet. Pode também ser executada fora do programa.
%%%
warning ("off");
%pkg install -local -forge  matgeom; %% descomentar para rodar a 1a vez, somente
%pkg load matgeom; % executar este comando apenas 1 vez, após abrir o Octave
warning ("on");
clc;
%%%  ======================================================================
%%%   Dados de entrada
%%%
NUSP = 12552562; % NUSP do 1o aluno (ordem alfab.)
%
a=10.0e-02; # dimensões em cm
b=6.0e-02;
c=3.0e-02;
d=3.0e-02;
g=3.0e-02;
h=(b-d)/2;
epsr=2;
sigma=3.0e-03;  # S/m
sigma_dual=3.5e-03; # S/m
eps0=8.85418782e-12;  % F/m
Vmin=0.0;     % Volts
Vmax=100.0;   % Volts
%%%
%%%  ======================================================================
%%%
%%%            Definicao do dominio
%%%%
%%%% A variavel dx abaixo é a discretização utilizada. Valores diferentes
%%%% daqueles sugeridos abaixo não funcionarão. Diminua o dx para gerar a
%%%% versão final a ser entregue.
%dx=0.05; % Tempo de execucao MUITO longo!!
%dx=0.1;   % Tempo de execucao longo!!
%dx=0.25;  % recomendado para a versao final
dx=0.5;   %% Mude para dx=0.25 somente quando for gerar os resultados finais!!!
erro=0.0;
start=start_Dual= 50;
iter=0;
dy=dx;
lx=a;
ly=b;
Nx=round(lx/dx)+1;
Ny=round(ly/dx)+1;
ring1= [0 0; lx 0; lx ly; 0 ly; 0 0];
ring2=[g h; g h+d; g+c h+d; g+c h; g h];
polyg={ring1,ring2};
verts = polygonVertices(polyg);
xgv=((1:Nx)-1)*dx;
ygv=((1:Ny)-1)*dx;
[x,y]=meshgrid(xgv,ygv);
verts1 = polygonVertices(ring1);
verts2 = polygonVertices(ring2);
xv1=verts1(:,1);
yv1=verts1(:,2);
xv2=verts2(:,1);
yv2=verts2(:,2);
[in1,on1] = inpolygon(x,y,xv1,yv1);
[in2,on2] = inpolygon(x,y,xv2,yv2);
%%%
%%%    Atribui condições de contorno
%%%
r=find(in1&~in2|on2); % tudo
p=find(in1&~on1&~in2); % só  nos internos
q=find(on1|on2); % só na fronteira
iVmax=find(on2);
iFuro=find(in2&!on2);
Phi_prev=zeros(size(x));
Phi_new=zeros(size(x));
Phi_new(iVmax)= Vmax;
Phi_new(iFuro)= NaN;
Phi_new(p)= start;
%%% ========================================================================
%%%
%%% Contador de iteraçõeste
iter=0;
%%% Erro máximo entre Phi_new e Phi_prev
erro=max(max(abs(Phi_new-Phi_prev)));
%%%
%%%             Laço iterativo - Método das Diferenças Finitas
%%%
while(erro > 1e-4 && iter < 1e4)% Executa até convergir ou atingir o maximo de iterações
    iter=iter+1; % Incrementa iteração
%%% Atualiza o potencial dos nós internos pela média dos 4 vizinhos - Eq. Laplace - M.D.F.
    for k=1:size(p,1);
        [i,j]=ind2sub(size(x),p(k));
            Phi_new(i,j)=(Phi_new(i-1,j)+Phi_new(i+1,j)+Phi_new(i,j-1)+Phi_new(i,j+1))/4;
    end
%%%% Calcula máximo erro entre Phi_atual e Phi_prev de todo o domínio
    erro=max(max(abs(Phi_new-Phi_prev)));
    eps(iter)=erro;
%%%% Atualiza a matriz de potenciais
    Phi_prev=Phi_new;

end
niter1=iter;
if (niter1 == 1e4 && erro > 1e-4)
	disp([' Numero maximo de iteracoes atingido sem convergencia :', num2stg(niter1), '  iteracoes \? Erro: \n', num2str(erro), 'Os resultados podem nao ter significado!\n']);
end