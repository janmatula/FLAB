clc; clear all; close all;
B = 6e20; %magneticka indukce, musime zvolit
E = 1e20; %intenzita elektrostatickeho pole, musime zvolit
delkaLetoveTrubice = 0.2;
z = 1;
q = 1* 1.6021766208e-19;%naboj iontu

m = 164; %hmotnost iontu v Daltonech, prvni trojcisli VUTID
Ek = 3; %kineticka energie iontu v kev, soucet poslednich dvou nenulovych cislic VUTID
rozptyl_Ek = 0.4 * (1+6+4+2+1+2)/6; %rozptyl kineticke enregie iontu
rozptyl_Osa = 0.1 * (1+6+4+2+1+2); %uhlovy rozptyl iontu od centralni osy zarizeni, stupne

%prevody jednotek
% m = 1.66053904*10^(-27)*m; %hmotnost na kg
% Ek = Ek/ 6.2415064799632E+15; % kiteticka na jouly
% rozptyl_Ek = rozptyl_Ek/6.2415064799632E+15; %rozptyl na jouly
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = 3; %pocet iontu

rozptylIontu = normrnd(0, rozptyl_Osa, [1, n]);
normaly = tan(deg2rad(rozptylIontu));

rozptylenaEk = normrnd(Ek, rozptyl_Ek, [1, n]);


delka = 0:0.001:0.2;
x = ones(length(delka), n);
y = [];

primkaElPole2 = @(x) 1.6116 * x -0.4;

primkaMagnet1 =  @(x) 1.6116 * x - 0.7667;

primkaMagnet2 = @(x) -0.6205 * x + 0.3;



figure
plot(0.2:0.001:0.28,primkaElPole2(0.2:0.001:0.28), 'g')
xlim([0 0.8])
hold on
plot(0.4:0.001:0.478,primkaMagnet1(0.4:0.001:0.478), 'r')
plot(0.478:0.001:0.63,primkaMagnet2(0.478:0.001:0.63), 'r')
%line([0.6 0.6], [-0.3 -0.06]);
line([0.2 0.2], [-0.07768 0.05], 'Color','green');
%line([0.9 0.9], [-0.1 0], 'Color', 'red');
rElek = (2*rozptylenaEk)/(q*E);
rMagnet = (sqrt(2*m*rozptylenaEk))/(q*B);
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

% while y(end)>primkaMagnet2(delkovyIterator)
%     y = [y; y(end)+smerY];
%     delkovyIterator = delkovyIterator + 0.001;
% end


while y(end)>primkaMagnet1(delkovyIterator)
    y = [y; y(end)+smerY];
    delkovyIterator = delkovyIterator + 0.001;
end

pocatek = y(end, :);
for j = 0:0.001:rMagnet(i)
    %nuloveY = rMagnet(i) - sqrt(rMagnet(i)*rMagnet(i) - rMagnet(i)/2*rMagnet(i)/2);
    y = [y; -sqrt(rMagnet(i)*rMagnet(i) - j*j)+rMagnet(i)+pocatek];
    delkovyIterator = delkovyIterator + 0.001;
    if y(end)>primkaMagnet2(delkovyIterator)
        smerY =y(end) - y(end-1);
        break;
    end
   
end

%detektor
while delkovyIterator<=1.2
    y = [y; y(end) + smerY];
    delkovyIterator = delkovyIterator + 0.001;
end



plot(0:0.001:(length(y)-1)/1000, y, '.')
hold on

end

