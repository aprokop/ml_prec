function [] = analysis()

global Az
if ~exist('Az', 'var')
    fprintf('Loading Az...');
    Az = mmread('~/Code/Work/spe_prec/project/az.mm');
    fprintf(' ok\n');
end

global nx
global ny
global nz
nx = 60;
ny = 220;
nz = 85;

%%
if 0
    max_x = 15;
    max_y = 1;
    rho = zeros(max_x,max_y);
    indices = zeros(nz,1);
    for j = 1:max_y
        for i = 1:max_x
            for k = 1:nz
                indices(k) = INDEX(i,j,k);
            end
            Az3 = Az(indices,indices);
            D = diag(diag(Az3));
            rho(i,j) = max(abs(eig(eye(85)-D\Az3)));
        end
    end
    plot(rho);
end
%%
if 1
    c = zeros(nz,1);
    max_x = 1;
    max_y = 5;
    ind = 1;
    for j = 1:max_y
        for i = 1:max_x
            for k = 1:nz
                ind = INDEX(i,j,k);
                c(k) = 1 - C(Az,ind)/Az(ind,ind);
            end
            plot(c);
            grid on;
            input('>');
        end
    end

end
%%

    function [id] = INDEX(i,j,k)
        id = (k-1)*ny*nx + (j-1)*nx + (i-1) + 1;
    end

    function [c] = C(A,i)
        c = sum(A(i,:));
    end
end
