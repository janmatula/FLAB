clc; clear all; close all;

%VUT ID 164212

B = 0.3; %magneticka indukce, musime zvolit [T]
E = 40000; %intenzita elektrostatickeho pole, musime zvolit [V/m]
n = 1000; %pocet iontu
poziceDetektoru = 0.65;


delkaLetoveTrubice = 0.2;
z = 1;
q = 1* 1.6021766208e-19;%naboj iontu


m = 164; %hmotnost iontu v Daltonech, prvni trojcisli VUTID
Ek = 3; %kineticka energie iontu v kev, soucet poslednich dvou nenulovych cislic VUTID
rozptyl_Ek = 0.4 * (1+6+4+2+1+2)/6; %rozptyl kineticke enregie iontu
rozptyl_Osa = 0.1 * (1+6+4+2+1+2); %uhlovy rozptyl iontu od centralni osy zarizeni, stupne

%prevody jednotek
q = 0.001*z; %prevod q na keV
q = q/6.2415064799632e15; %q na jouly
m = m*1.6605e-27; %prevod Daltonu na kg
Ek = Ek/6.2415064799632e15; %Ek na Jouly
rozptyl_Ek = rozptyl_Ek/6.2415064799632e15; %rozptyl Ek na Jouly



rozptylIontu = normrnd(0, rozptyl_Osa, [1, n]);
normaly = tan(deg2rad(rozptylIontu));

rozptylenaEk = normrnd(Ek, rozptyl_Ek, [1, n]);


delka = 0:0.001:0.2;
x = ones(length(delka), n);
y = [];

primkaElPole2 = @(x) 1.6116 * x -0.4;

primkaMagnet1 =  @(x) 1.6116 * x - 0.7667;

primkaMagnet2 = @(x) -0.6205 * x + 0.3;


%vyznaceni magnetickych a elektrickych casti
figure
plot(0.2:0.001:0.28,primkaElPole2(0.2:0.001:0.28), 'g')
xlim([0 0.8])
hold on
plot(0.4:0.001:0.478,primkaMagnet1(0.4:0.001:0.478), 'r')
plot(0.478:0.001:0.63,primkaMagnet2(0.478:0.001:0.63), 'r')
line([0.2 0.2], [-0.07768 0.05], 'Color','g');
line([poziceDetektoru poziceDetektoru], [-0.15 0.1], 'Color','blue');

%vypocet polomeru kruznic
rElek = (2*rozptylenaEk)/(q*E);
rMagnet = (sqrt(2*m*rozptylenaEk))/(q*B);

%hodnoty na y ose pro vypocet pozice detektoru
poziceIontuNaDetektoru = [];

for i = 1:n
    y = [];
    delkovyIterator = 0;
    
    %letova trubice
    for j = 0:0.001:0.2
        y = [y; j.*normaly(i)];
        delkovyIterator = delkovyIterator + 0.001;
    end
    
    %elektricke pole
    pocatek = y(end, :);
    
    for j = 0:0.001:rElek(i)
        %vykresleni kruznice
        y = [y; sqrt(rElek(i)*rElek(i) - j*j)-rElek(i)+pocatek];
        delkovyIterator = delkovyIterator + 0.001;
        if y(end)<primkaElPole2(delkovyIterator)
            smerY =y(end) - y(end-1);
            break;
        end
    end
    %vypocet smernice pro cestu mezi elektro a magnetem
    smerY =y(end) - y(end-1);
    
    %cesta z elektrostaticke do magneticke casti
    while y(end)>primkaMagnet1(delkovyIterator)
        y = [y; y(end)+smerY];
        delkovyIterator = delkovyIterator + 0.001;
    end
    
    %osetreni bugu
    if y(end)<-0.1
        continue
    end
    
    %a hura do magneticke casti
    pocatek = y(end, :);
    %pocatekKruznice = rad2deg(atan(abs(smerY/0.001) + -0.6205));
    smerY = 'n';
    for j = (-(rMagnet(i)*31.28/90)+0.001:0.001:rMagnet(i))
        nuloveY = rMagnet(i) - sqrt(rMagnet(i)*rMagnet(i) - ((rMagnet(i)*31.28/90))*((rMagnet(i)*31.28/90)));
        y = [y; -sqrt(rMagnet(i)*rMagnet(i) - j*j)+rMagnet(i)+pocatek-nuloveY];
        delkovyIterator = delkovyIterator + 0.001;
        if y(end)>primkaMagnet2(delkovyIterator)
            smerY =y(end) - y(end-1);
            break;
        end
        
    end
    
    %osetreni bugu
    if smerY == 'n'
        continue
    end
    
    %detektooor
    while delkovyIterator<=1
        y = [y; y(end) + smerY];
        delkovyIterator = delkovyIterator + 0.001;
        if delkovyIterator > poziceDetektoru
            poziceIontuNaDetektoru = [poziceIontuNaDetektoru y(end)];
            break
        end
    end
    
    
    %vykresleni trajektorie daneho iontu
    plot(0:0.001:(length(y)-1)/1000, y, '.')
    hold on
    
end

%vypocet smerodatne odchylky v miste detektoru
disp(['Smerodatna odchylka iontu v miste detektoru je : ' num2str(std(poziceIontuNaDetektoru))])
