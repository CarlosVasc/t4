clear
clc

i = complex(0,1);
f = logspace(1,8, 100);
w = 2*pi*f;

%-----alterar-----
C1 = 1e-06;
C2 = 1e-03;
RB1 = 80000;
RB2 = 20000;
RC1 = 1000;
RE1 = 100;
%-----------------

vin = 0.01;
VT=25e-3;
BFN=178.7;
VAFN=69.7;
VBEON=0.7;
VCC=12;
RS=100;

RB=1/(1/RB1+1/RB2);
VEQ=RB2/(RB1+RB2)*VCC;
IB1=(VEQ-VBEON)/(RB+(1+BFN)*RE1);
IC1=BFN*IB1;
IE1=(1+BFN)*IB1;
VE1=RE1*IE1;
VO1=VCC-RC1*IC1;
VCE=VO1-VE1;
gm1 = IC1/VT;
rpi1 = BFN/gm1;
ro1 = VAFN/IC1;
Zc1 = 1./(i*w*C1);
Zc2 = 1./(i*w*C2);
ZI1 = 1/(1/RB + 1/rpi1);
ZO1 = 1/(1/ro1 + 1/RC1);

for k=1:length(w)
A = [0,          -RB,            0,            0,               RS+RB;
    -RE1,      -Zc2(k),          0,       RE1 + Zc2(k),           0; 
    RE1 + ro1 + RC1, 0,        -ro1,        -RE1,                0;
    0,         gm1*rpi1,         1,            0,                 0;
    0,        RB+Zc1(k)+rpi1+Zc2(k), 0,     -Zc2(k),            -RB];

B = [vin; 0; 0; 0; 0];

X=linsolve(A,B);
  
voE(k) = abs(RC1 * X(1));
AV1(k) = voE(k)/vin;
end

AV1_db=20*log10(AV1);
fid1 = figure ();
plot (log10(f), AV1_db, "b");
legend("v_o(f)/v_i(f)");
xlabel ("Frequency [Hz]");
ylabel ("Gain [db]");
print (fid1, "Gain1.eps", "-depsc");

%ouput stage
BFP = 227.3;
VAFP = 37.2;
RE2 = 100;
VEBON = 0.7;
VI2 = VO1;
IE2 = (VCC-VEBON-VI2)/RE2;
IC2 = BFP/(BFP+1)*IE2;
VO2 = VCC - RE2*IE2;

gm2 = IC2/VT;
go2 = IC2/VAFP;
gpi2 = gm2/BFP;
ge2 = 1/RE2;
ro2 = 1/go2;
rpi2 = 1/gpi2;

for k=1:length(w)
A2 = [rpi2+RE2,   -RE2,     0;
      gm2*rpi2,     0,      1; 
       -RE2,     RE2+ro2, -ro2];

B2 = [voE(k); 0; 0];

X2=linsolve(A2,B2);
vo2(k) = abs((X2(1)-X2(2))*RE2);
AV2(k) = vo2(k)/voE(k);
end


ZI2 = (gm2+gpi2+go2+ge2)/gpi2/(gpi2+go2+ge2);
ZO2 = 1/(gm2+gpi2+go2+ge2);
ZO = 1/(gm2*(rpi2/(rpi2+ZO1))+(1/(rpi2+ZO1))+go2+ge2);
AV = AV1.*AV2;
AV_db=20*log10(AV);

fid2 = figure ();
plot (log10(f), AV_db,log10(f),max(AV_db)-3);
legend("v_o(f)/v_i(f)");
xlabel ("Frequency [Hz]");
ylabel ("Gain [db]");
print (fid2, "GainFinal.eps", "-depsc");


k = 1;
while 3 < ((max(AV_db) - AV_db(k)))
    k = k + 1;
end
lowerCutoffFreq = f(k);
bandwidth = 10e6 - lowerCutoffFreq;

cost = 1e-3*(RE1 + RC1 + RB1 + RB2 + RE2) + 1e6*(C1 + C2) + 2*0.1;
Merit = (max(AV) * bandwidth)/(cost * lowerCutoffFreq)

