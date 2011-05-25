clear all
load data/RGB
load data/GRA31
load data/surfs1995   

S = RGB;

[u,s,v] = svd(surfs1995);
R = u(:,1:3);

[u,s,v] = svd(GRA31);
I(:,1) = mean(GRA31');
I(:,2:3) = u(:,1:2);
I=u(:,1:3);

T = zeros(3,3,3);
for i=1:3,
    for l=1:3,
        for j=1:3,            
            T(i,l,j) = S(:,i)'*(I(:,l).*R(:,j));
        end
    end
end


Q = zeros(3,3);
for a=1:3,
    for b=1:3,
        Q(a,b) = R(:,a)'*R(:,b);
    end
end

K = surfs1995;

C = zeros(3, size(K,2));
for k=1:size(K,2),
    C(:, k)=inv(Q)*(R' * K(:,k));
end

U = mean(C')';


B = zeros(3,3);
for j=1:3
    B = B + T(:,:,j)*U(j);
end


illum = GRA31(:, randi(size(GRA31, 2)));
surfaces = surfs1995(:,randperm(size(surfs1995, 2)));
surfaces = surfaces(:,1:40);
img = ((repmat(illum, 1, size(surfaces, 2)).*surfaces)' * S);
V = mean(img)';

e = inv(B)*V;
illumest = I*e;
plot(illum./max(illum))
hold on
plot(illumest./max(illumest), 'green')
hold off


A = zeros(3,3);
for j=1:3
    A = A + T(:,:,j)*e(j);
end

C=inv(A)*img';

for s = 1:3;
    surf=surfaces(:,s);
    surfest = R*C(:,s);
    figure;
    plot(surf./max(surf))
    hold on
    plot(surfest./max(surfest), 'green')
    hold off
end