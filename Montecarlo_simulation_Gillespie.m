clear all; close all; clc;

initNA1 = 10000;
initNA2 = 0;
initNB = 10000;
A1size = initNA1;
A2size = initNA2;
Bsize = initNB;

InitialVecPopulations = [initNA1, initNA2, initNB];

Wfit = 1.0;
Rfit = 0.1;
EB = 1;
tB = 0;
fA1 = Rfit;
fA2 = Wfit;
fB = Rfit;

AvTimeEvents = 2.5e - 05;
tauB = 5 * AvTimeEvents;
MA = 1;
MaxNumIter = 200000;
NumbChangEnvir = MaxNumIter / 10;

RelPopA = 0;
RelPopB = 0;
ProbCambio = 1 / 10;
A22A1 = 0;
A12A2 = 0;
i = 0;
t = 0;

while (RelPopA < 0.9 && RelPopB < 0.9)
    RelPopA = (A1size + A2size) / (A1size + A2size + Bsize);
    RelPopB = Bsize / (A1size + A2size + Bsize);
    i = i + 1;
    Medio(i) = MA;
    EstadoB(i) = EB;
    tiempoResp(i) = tB;
    tB = tB + AvTimeEvents;
    if (mod(i, NumbChangEnvir) == 0)
        MA = mod(MA, 2) + 1;
        tB = 0;
    end
    if (MA == 2)
        fA1 = Wfit;
        fA2 = Rfit;
    else
        fA1 = Rfit;
        fA2 = Wfit;
    end
    if (EB == MA)
        fB = Rfit;
    else
        fB = Wfit;
    end
    if (tB >= tauB)
        EB = MA;
    end
    RA1 = A1size * fA1;
    RA2 = A2size * fA2;
    RB = Bsize * fB;
    DA1 = A1size;
    DA2 = A2size;
    DB = Bsize;
    T12 = ProbCambio * A1size;
    T21 = ProbCambio * A2size;
    alpha = RA1 + RA2 + RB + DA1 + DA2 + DB + T12 + T21;
    tiempo(i) = - log(rand(1, 1)) / (alpha);
    steptime = tiempo(i);
    t = tiempo(i) + t;
    ProbRepA1 = RA1 / alpha;
    ProbRepA2 = RA2 / alpha;
    ProbRepB = RB / alpha;
    ProbDeathA1 = DA1 / alpha;
    ProbDeathA2 = DA2 / alpha;
    ProbDeathB = DB / alpha;
    ProbTrans12 = T12 / alpha;
    ProbTrans21 = T21 / alpha;
    z = rand(1, 1);
    if (z < ProbRepA1)
        A1size = A1size + 1;
    elseif (z < ProbRepA1 + ProbRepA2)
        A2size = A2size + 1;
    elseif (z < ProbRepA1 + ProbRepA2 + ProbRepB)
        Bsize = Bsize + 1;
    elseif (z < ProbRepA1 + ProbRepA2 + ProbRepB + ProbDeathA1)
        A1size = A1size - 1;
    elseif (z < ProbRepA1 + ProbRepA2 + ProbRepB + ProbDeathA1 + ProbDeathA2)
        A2size = A2size - 1;
    elseif (z < ProbRepA1 + ProbRepA2 + ProbRepB + ProbDeathA1 + ProbDeathA2 + ProbDeathB)
        Bsize = Bsize - 1;
    elseif (z < ProbRepA1 + ProbRepA2 + ProbRepB + ProbDeathA1 + ...
          ProbDeathA2 + ProbDeathB + ProbTrans12)
        A1size = A1size - 1;
        A2size = A2size + 1;
        A12A2 = A12A2 + 1;
    elseif (z < ProbRepA1 + ProbRepA2 + ProbRepB + ProbDeathA1 + ProbDeathA2 + ...
          ProbDeathB + ProbTrans12 + ProbTrans21)
        A2size = A2size - 1;
        A1size = A1size + 1;
        A22A1 = A22A1 + 1;
    end
    time(i) = t;
    populA1(i) = A1size;
    populA2(i) = A2size;
    populB(i) = Bsize;
    if (i == MaxNumIter)
        RelPopA = 1;
        y = 0;
    else
        y = 1;
    end
end

FinalVecPopulations = [A1size, A2size, Bsize];
x = [1, 2, 3];
K = 0.5;
figure;

bar1 = bar(x, InitialVecPopulations, ’FaceColor’, ’r’, ’EdgeColor’, ’r’);

hold on

bar2 = bar(x, FinalVecPopulations, ’FaceColor’, ’b’, ’EdgeColor’, ’b’);

set(bar1, ’BarWidth’, K / 2);
set(bar2, ’BarWidth’, K / 3);

hold off;

legend(’Initial Populations’, ’Final Populations’)
axis([0 4 0 (100 + max([A1size A2size Bsize initNA1 initNA2 initNB]))])

figure(2);
subplot(3, 1, 1)

stairs(time, populA1)
title(’A_1’)
ylabel(’n_{A1}’)
xlabel(’t’)

subplot(3, 1, 2)
stairs(time, populA2)
title(’A_2’)
ylabel(’n_{A2}’)
xlabel(’t’)

subplot(3, 1, 3)
stairs(time, populB)
title(’B’)
ylabel(’n_{B}’)
xlabel(’t’)
figure;

plot(time, Medio)

hold on

plot(time, EstadoB, ’g’)