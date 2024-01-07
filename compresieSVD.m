function [imagineComprimata,Cr,MSE]=compresieSVD(denumireImagine,CrDorit)
    arguments(Input)
        denumireImagine (1,:) char
        CrDorit (1,1) double {mustBePositive,mustBeLessThanOrEqual(CrDorit,1)}=1
    end
    arguments(Output)
        imagineComprimata (:,:,3) uint8
        Cr (1,1) double
        MSE (1,1) double
    end
    imagine=imread(denumireImagine);
    [uR,sigmaR,vR]=svd(double(imagine(:,:,1)),'econ');
    [uG,sigmaG,vG]=svd(double(imagine(:,:,2)),'econ');
    [uB,sigmaB,vB]=svd(double(imagine(:,:,3)),'econ');
    [rangInitial,~]=size(sigmaR);
    salveazaMatrici();
    gradDeCompresie=calculeazaGradDeCompresie();
    [uR,sigmaR,vR]=eliminaVSsiRedimensioneaza(uR,sigmaR,vR,gradDeCompresie);
    [uG,sigmaG,vG]=eliminaVSsiRedimensioneaza(uG,sigmaG,vG,gradDeCompresie);
    [uB,sigmaB,vB]=eliminaVSsiRedimensioneaza(uB,sigmaB,vB,gradDeCompresie);
    imagineComprimata=obtineImagineaComprimata(imagine,uR,sigmaR,vR,uG,sigmaG,vG,uB,sigmaB,vB);
    Cr=calculeazaCr();
    clear rangInitial;
    salveazaMatrici();
    clear uR sigmaR vR uG sigmaG vG uB sigmaB vB;
    MSE=calculeazaMSE(imagine,imagineComprimata);
    comparaDimensiuni();
    imshow([imagine imagineComprimata]);

    function Cr=calculeazaCr
        Cr=numel(sigmaR)/rangInitial; % Ca sa fie subunitar
    end

    function salveazaMatrici
        persistent nrApelari;
        if isempty(nrApelari)
            denumireU='U';
            denumireSigma='sigma';
            denumireV='V';
            nrApelari=1;
        else
            denumireU='Ucomprimat';
            denumireSigma='sigmaComprimat';
            denumireV='Vcomprimat';
        end
        matrice=[uR uG uB];
        save(denumireU,'matrice');
        matrice=[sigmaR sigmaG sigmaB];
        save(denumireSigma,'matrice');
        matrice=[vR vG vB];
        save(denumireV,'matrice');
    end

    function gradDeCompresie=calculeazaGradDeCompresie()
        gradDeCompresie = (1-CrDorit)*rangInitial;
    end

end

function[U,sigma,V]=eliminaVSsiRedimensioneaza(U,sigma,V,nrValoriSaFieSterse)
    arguments(Input)
        U (:,:) double {ortogonalValidator}
        sigma (:,:) double {sigmaValidator}
        V (:,:) double {ortogonalValidator}
        nrValoriSaFieSterse (1,1) uint64 {valoareValida(nrValoriSaFieSterse,sigma)}
    end
    arguments(Output)
        U (:,:) double {ortogonalValidator}
        sigma (:,:) double {sigmaComprimatValidator}
        V (:,:) double {ortogonalValidator}
    end
    U=U(:,1:end-nrValoriSaFieSterse);
    %sigma=diag(sigma(1:end-nrValoriSaFieSterse,1:end-nrValoriSaFieSterse));
    sigma=diag(sigma);
    sigma=sigma(1:end-nrValoriSaFieSterse,1);
    V=V(:,1:end-nrValoriSaFieSterse);
end

function ortogonalValidator(Q)
    arguments
        Q (:,:) double
    end
    %{
    [m,n]=size(Q);
    if Q*transpose(Q) ~= eye(m)
        error('Matricea nu este ortogonala.');
    end
    if Q'*Q ~= eye(n)
        error('Matricea nu este ortogonala.');
    end
    %}
end

function sigmaValidator(sigma)
    arguments
        sigma (:,:) double {diagonalValidator}
    end
    %{
    [m,n]=size(sigma);
    m=min(m,n);
    clear n;
    if sigma(1,1)<=0
        error('Cel putin o valoare singula este negativa sau nula.');
    end
    for i=2:m
        if sigma(i,i)<=0
            error('Cel putin o valoare singulara este negativa sau nula.');
        elseif sigma(i-1,i-1)<sigma(i,i)
            error('Valorile singulare nu sunt descrescatoare.');
        end
    end
    %}
end

function diagonalValidator(sigma)
    arguments
        sigma (:,:) double
    end
    %{
    [m,n]=size(sigma);
    for i=1:m
        for j=1:i-1
            if sigma(i,j)~=0
                error('Matricea nu este diagonala.');
            end
        end
        for j=i+1:n
            if sigma(i,j)~=0
                error('Matricea nu este diagonala.');
            end
        end
    end
    %}
end

function sigmaComprimatValidator(sigma)
arguments
    sigma(:,1) double
end
    for i=2:numel(sigma)
        if sigma(i)<=0
            error('Cel putin o valoare singulara este negativa sau nula.');
        elseif sigma(i-1)<sigma(i)
            error('Valorile singulare nu sunt descrescatoare.');
        end
    end
end

function valoareValida(nrValoriSaFieSterse,sigma)
arguments
    nrValoriSaFieSterse (1,1) uint16
    sigma (:,:) double {sigmaValidator}
end
    [m,n]=size(sigma);
    m=min(m,n);
    uint8(7);
    for r=m:-1:1
        if sigma(r,r)~=0
            break;
        end
    end
    if nrValoriSaFieSterse >= r
        error('Numărul de valori singulare dorite a fi șterse (%u) este mai mare sau egal decât numărul de valori singulare din imagine (adică %u).',nrValoriSaFieSterse,length(sigma));
    end
end

function[imagineReconstruita]=obtineImagineaComprimata(image,uR,sigmaR,vR,uG,sigmaG,vG,uB,sigmaB,vB)
    arguments(Input)
        image (:,:,3) uint8
        uR (:,:) double {ortogonalValidator}
        sigmaR (:,:) double {sigmaComprimatValidator}
        vR (:,:) double {ortogonalValidator}
        uG (:,:) double {ortogonalValidator}
        sigmaG (:,:) double {sigmaComprimatValidator}
        vG (:,:) double {ortogonalValidator}
        uB (:,:) double {ortogonalValidator}
        sigmaB (:,:) double {sigmaComprimatValidator}
        vB (:,:) double {ortogonalValidator}
    end
    arguments(Output)
        imagineReconstruita (:,:,3) uint8
    end
    [m,n,~]=size(double(image));
    imagineReconstruita=uint8(zeros(m,n,3));
    clear m n;    
    imagineReconstruita(:,:,1)=uint8(uR*diag(sigmaR)*vR');
    imagineReconstruita(:,:,2)=uint8(uG*diag(sigmaG)*vG');
    imagineReconstruita(:,:,3)=uint8(uB*diag(sigmaB)*vB');
end

function ret=calculeazaMSE(imagine,imagineComprimata)
arguments(Input)
    imagine(:,:,3) uint8
    imagineComprimata(:,:,3) uint8
end
arguments(Output)
    ret (1,1) double
end
    ret=sum(sum((imagine-imagineComprimata).^2))/numel(imagine);
    ret=(ret(1)+ret(2)+ret(3))/3;
end

function comparaDimensiuni
    U=[dir('U.mat').bytes dir('Ucomprimat.mat').bytes
        dir('sigma.mat').bytes dir('sigmaComprimat.mat').bytes
        dir('V.mat').bytes dir('Vcomprimat.mat').bytes];
    if U(1,1)>U(1,2) && U(2,1) > U(2,2) && U(3,1) > U(3,2)
        disp('Imaginea a fost comprimata');
    end
    disp(U);
end