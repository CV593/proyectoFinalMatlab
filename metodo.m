function metodo()
    coeffs = str2num(input('Ingrese coeficientes del polinomio (separados por espacios, del coeficiente principal al independiente): ', 's'));
    n = length(coeffs) - 1;
    if n < 2
        error('El polinomio debe ser de al menos grado 2.');
    end
    u = str2double(input('Ingrese valor inicial de u (por defecto 0): ', 's'));
    if isnan(u) || ~isreal(u)
        u = 0; 
    end
    v = str2double(input('Ingrese valor inicial de v (por defecto 0): ', 's'));
    if isnan(v) || ~isreal(v)
        v = 0;
    end
    errorTol = str2double(input('Ingrese error porcentual deseado (ej. 1 para 1%, por defecto 0.1): ', 's'));
    if isnan(errorTol) || errorTol <= 0 
        errorTol = 0.1; 
    end
    errorTol = errorTol / 100;
    [roots, iterData] = bairstow(coeffs, u, v, errorTol);
    if isempty(roots)
        disp('No se encontraron raíces reales con el método de Bairstow.');
    else
        disp('Raíces encontradas (forma algebraica):');
        for i = 1:length(roots)
            fprintf('x%d = %f + %fi\n', i-1, real(roots(i)), imag(roots(i)));
        end
        disp('Resultados de las iteraciones:');
        fprintf('Iter\t\tu\t\tv\t\tError Relativo\t\tResiduos\n');
        for i = 1:length(iterData)
            fprintf('%d\t\t%f\t%f\t%f\t\t\t%f\n', iterData(i).iter, iterData(i).u, iterData(i).v, iterData(i).error, iterData(i).residual);
            if iterData(i).error < errorTol && iterData(i).residual < errorTol
                break; 
            end
        end
        figure;
        fplot(@(x) polyval(coeffs, x));
        hold on;
        plot(real(roots), imag(roots), 'ro', 'MarkerSize', 8);
        title('Gráfico de f(x) con raíces encontradas (Método de Bairstow)');
        xlabel('Parte Real (x)');
        ylabel('Parte Imaginaria (y)');
        legend('f(x)', 'Raíces');
        grid on;
        axis equal;
        xlim([min(real(roots))-1, max(real(roots))+1]);
        ylim([min(imag(roots))-1, max(imag(roots))+1]);
    end
end

function [roots, iterData] = bairstow(coeffs, u, v, errorTol)
    n = length(coeffs) - 1;
    roots = [];
    iterData = struct('iter', {}, 'u', {}, 'v', {}, 'error', {}, 'residual', {});
    iter = 0;
    while n >= 3 && iter < 100 
        iter = iter + 1;
        [b, c] = syntheticDivision(coeffs, u, v);
        denom = c(end-2)^2 - c(end-1) * c(end-3);
        if denom ~= 0
            du = (b(end-1) * c(end-1) - b(end) * c(end-2)) / denom;
            dv = (b(end) * c(end-3) - b(end-1) * c(end-2)) / denom;
            u = u + du;
            v = v + dv;
        else
            du = 0;
            dv = 0;
        end
        error = max(abs([du/u, dv/v]));
        residual = abs(b(end));
        iterData(iter).iter = iter;
        iterData(iter).u = u;
        iterData(iter).v = v;
        iterData(iter).error = error;
        iterData(iter).residual = residual;
        if error < errorTol && residual < errorTol
            break;
        end
        delta = u^2 - 4 * v;
        if delta >= 0
            roots = [roots; (u + sqrt(delta)) / 2; (u - sqrt(delta)) / 2];
        else
            roots = [roots; complex(u / 2, sqrt(-delta) / 2); complex(u / 2, -sqrt(-delta) / 2)];
        end
        coeffs = b(1:end-2);
        n = length(coeffs) - 1;
    end
    if n == 2
        delta = coeffs(2)^2 - 4 * coeffs(1) * coeffs(3);
        if delta >= 0
            roots = [roots; (-coeffs(2) + sqrt(delta)) / (2 * coeffs(1)); (-coeffs(2) - sqrt(delta)) / (2 * coeffs(1))];
        else
            roots = [roots; complex(-coeffs(2) / (2 * coeffs(1)), sqrt(-delta) / (2 * coeffs(1))); complex(-coeffs(2) / (2 * coeffs(1)), -sqrt(-delta) / (2 * coeffs(1)))];
        end
    elseif n == 1
        roots = [roots; -coeffs(2) / coeffs(1)];
    end
end
function [b, c] = syntheticDivision(coeffs, u, v)
    n = length(coeffs) - 1;
    b = zeros(1, n + 1);
    c = zeros(1, n + 1);
    b(1) = coeffs(1);
    b(2) = coeffs(2) + u * b(1);
    for i = 3:n + 1
        b(i) = coeffs(i) + u * b(i - 1) + v * b(i - 2);
    end
    c(1) = b(1);
    c(2) = b(2) + u * c(1);
    for i = 3:n
        c(i) = b(i) + u * c(i - 1) + v * c(i - 2);
    end
end
