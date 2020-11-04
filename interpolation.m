clc;
clear variables;
close all force;

a = 0;
b = 1;

Na = 8;
dxa = (b - a) / Na;
x1 = (a : dxa : b);
N1 = length(x1);
fprintf("N1 =\n");
disp(N1);
f1(:) = fun(x1(:));

Nb = 10 * Na;
dxb = (b - a) / Nb;
x0 = (a : dxb : b);
N0 = length(x0);
fprintf("N0 =\n");
disp(N0);
f0(:) = fun(x0(:));

f2 = sosed(f1, x0, x1);
f3 = line(f1, x0, x1);
f4 = Lagrange(f1, x0, x1);

f5 = forwardNewton(f1, x0, x1);
f6 = backNewton(f1, x0, x1);
f7 = kubSplain(f1, x0, x1);

subplot(2, 1, 1);
%модель аналитической функции
plot(x0, f0, 'LineWidth', 4);
hold on;
%модель экспериментальных данных
plot(x1, f1, 'x', 'LineWidth', 8);
hold on;
%метод ближайшего соседа
stairs(x0, f2);
hold on;
%линейна€ интерпол€ци€
plot(x0, f3, 'LineWidth', 2);
hold on;
%полиномом Ћагранжа
plot(x0, f4, 'LineWidth', 2);
hold off;
ylim([-1, 1]);
grid on;
grid minor;
legend('модель аналитической функции', 'модель экспериментальных данных', 'метод ближайшего соседа', 'линейна€ интерпол€ци€', 'полиномом Ћагранжа', 'Location', 'NorthEast');

subplot(2, 1, 2);
%модель аналитической функции
plot(x0, f0, 'LineWidth', 6);
hold on;
%модель экспериментальных данных
plot(x1, f1, 'x', 'LineWidth', 8);
hold on;
%формула Ќьютона (вперед)
plot(x0, f5, 'LineWidth', 4);
hold on;
%формула Ќьютона (назад)
plot(x0, f6, 'LineWidth', 2);
hold on;
%кубические сплайны
plot(x0, f7, 'LineWidth', 1);
hold off;
ylim([-1, 1]);
grid on;
grid minor;
legend('модель аналитической функции', 'модель экспериментальных данных', 'формула Ќьютона (вперед)', 'формула Ќьютона (назад)', 'кубические сплайны', 'Location', 'NorthEast');

function f = sosed(f1, x, y)
    k = 1;
    f = zeros(size(x));
    for i = 1 : 1 : length(x)
        if k ~= length(y)
            del1 = abs(x(i) - y(k));
            del2 = abs(x(i) - y(k + 1));
            if (del2 <= del1)
                k = k + 1;
            end
        end
        f(i) = f1(k);
    end
end

function f = line(f1, x, y)
    k = 1;
    f = zeros(size(x));
    for i = 1 : 1 : length(x)
        f(i) = f1(k) + (f1(k + 1) - f1(k)) / (y(k + 1) - y(k)) * (x(i) - y(k));
        if (x(i) >= y(k + 1) && k < length(y) - 1)
            k = k + 1;
        end
    end
end

function f = Lagrange(f1, x, y)
    f = zeros(size(x));
    for j = 1 : 1 : length(x)
        L = 0;
        for k = 1 : 1 : length(y)
            I = 1;
            for i = 1 : 1 : length(y)
                if i ~= k
                    I = I .* (x(j) - y(i)) / ((y(k) - y(i)));
                end
            end
            L = L + f1(k) * I;
        end
        f(j) = L;
    end
end

function f = forwardNewton(f1, x, y)
    k = length(y);
    h = y(2) - y(1);
    f = zeros(size(x));
    DEL = diff(f1, k);
    
    for j = 1 : 1 : length(x)
        q = (x(j) - y(1)) / h;
        f(j) = 0;
        Pr = 1;
        n = 0;
        
        while (n <= k - 1)
            f(j) = f(j) + DEL(1, n + 1) * Pr;
            n = n + 1;
            Pr = Pr .* (q - n + 1) / n;
        end
    end
end

function f = backNewton(f1, x, y)
    k = length(y);
    h = y(k) - y(k - 1);
    DEL = diff(f1, k);
    f = zeros(size(x));
    
    for j = 1 : 1 : length(x)
        q = (x(j) - y(k)) / h;
        f(j) = 0;
        Pr = 1;
        n = 0;
        
        while (n <= k - 1)
            f(j) = f(j) + DEL(k - n, n + 1) * Pr;
            n = n + 1;
            Pr = Pr .* (q + n - 1) / n;
        end
    end
end

function f = kubSplain(f1, x, y)
    n = length(y) - 1;
    [c(1), c(n + 1), K(1), L(1)] = deal(0);
    f = zeros(size(x));
    
    for k = 2 : 1 : n
        h(k) = y(k + 1) - y(k);
        h(k - 1) = y(k) - y(k - 1);
        F(k) = 3 * ((f1(k + 1) - f1(k)) / h(k) - (f1(k) - f1(k - 1)) / h(k - 1));
        V(k) = 2 * (h(k) + h(k - 1));
        K(k) = (F(k) - h(k - 1) * K(k - 1)) / (V(k) - h(k - 1) * L(k - 1));
        L(k) = h(k) / (V(k) - h(k - 1) * L(k - 1));
    end
    
    for k = n : - 1 : 2
        c(k) = K(k) - L(k) * c(k + 1);
    end
    
    for k = 1 : 1 : n
        d(k) = (c(k + 1) - c(k)) / (3 * h(k));
        b(k) = (f1(k + 1) - f1(k)) / h(k) - c(k) * h(k) - d(k) * h(k) ^ 2;
        a(k) = f1(k);
    end
    
    k = 1;
    for i = 1 : 1 : length(x)
        h(k) = x(i) - y(k);
        f(i) = a(k) + b(k) * h(k) + c(k) * h(k) ^ 2 + d(k) * h(k) ^ 3;
        if (x(i) >= y(k + 1) && k < n)
            k = k + 1;
        end
    end
end

function DEL = diff(f1, p)
    DEL = zeros(p);
    DEL(:, 1) = f1;
    for n = 2 : 1 : p
        for i = 1 : 1 : p - n + 1
            DEL(i, n) = DEL(i + 1, n - 1) - DEL(i, n - 1);
        end 
    end
end

function DEL = diff2(f1, p)
    DEL = zeros(1, p);
    for n = 1 : 1 : p
        del = 0;
        for i = 0 : 1 : n - 1
            C = factorial(n) / (factorial(i) * factorial(n - i));
            del = del + (-1) ^ i * C * f1(n - i);
        end
        DEL(n) = del;
    end
end

function y = fun(x)
    y = sin(2 * pi * x);
end



















