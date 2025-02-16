
using System.Net.NetworkInformation;

enum TypeExtremum {
    None, Min, Max
}
class OptimumOneArg {
    public delegate double DelOneValueFunc(double x);

    // Функция для нахождения минимума
    // Метод Брента-Деккера
    private const double MACHINE_EPSILON = 2.220446049250313e-16; // Machine epsilon for double precision
    public static double Minimize(
        double a, 
        double b, 
        Func<double, double> f
    ){
        double t = MACHINE_EPSILON;
        double tolerance = 0.0;

        double m = 0.0;
        double p = 0.0;
        double q = 0.0;
        double r = 0.0;
        double s = 0.0;

        double a1 = a;
        double b1 = b;
        double fa = f(a1);
        double fb = f(b1);

        double c = a1;
        double fc = fa;
        double e = b1 - a1;
        double d = e;

        while (true) {
            if (Math.Abs(fc) < Math.Abs(fb)) {
                a1 = b1;
                b1 = c;
                c = a1;
                fa = fb;
                fb = fc;
                fc = fa;
            }

            tolerance = 2.0 * MACHINE_EPSILON * Math.Abs(b1) + t;
            m = 0.5 * (c - b1);

            if (Math.Abs(m) <= tolerance || fb == 0.0) {
                break;
            }

            if (Math.Abs(e) < tolerance || Math.Abs(fa) <= Math.Abs(fb)) {
                e = m;
                d = e;
            } else {
                s = fb / fa;

                if (a1 == c) {
                    p = 2.0 * m * s;
                    q = 1.0 - s;
                } else {
                    q = fa / fc;
                    r = fb / fc;
                    p = s * (2.0 * m * q * (q - r) - (b1 - a1) * (r - 1.0));
                    q = (q - 1.0) * (r - 1.0) * (s - 1.0);
                }

                if (p > 0.0) {
                    q = -q;
                } else {
                    p = -p;
                }

                s = e;
                e = d;

                if ((2.0 * p < 3.0 * m * q - Math.Abs(tolerance * q)) && (p < Math.Abs(0.5 * s * q))) {
                    d = p / q;
                } else {
                    e = m;
                    d = e;
                }
            }

            a1 = b1;
            fa = fb;

            if (tolerance < Math.Abs(d)) {
                b1 += d;
            } else if (0.0 < m) {
                b1 += tolerance;
            } else {
                b1 -= tolerance;
            }

            fb = f(b1);

            if ((fb > 0.0 && fc > 0.0) || (fb <= 0.0 && fc <= 0.0)) {
                c = a1;
                fc = fa;
                e = b1 - a1;
                d = e;
            }
        }

        return b1;
    }

    // Функция определения границ интервала для унимодальных функций
    public static (double, double) BracketMinimum(
        Func<double, double> f,
        // DelOneValueFunc f 
        double x = 0, 
        double s = 1e-2, 
        double k = 2.0
    ) {
        // Инициализация начальных точек
        double a = x, ya = f(a);
        double b = a + s, yb = f(b);
        
        // Проверка направления поиска
        if (yb > ya) {
            // Меняем местами a и b
            (a, b) = (b, a);
            (ya, yb) = (yb, ya);
            s = -s; // Инвертируем шаг
        }
        while (true) {
            // Вычисляем новую точку c
            double c = b + s, yc = f(c);

            // Проверяем условие завершения
            if (yc > yb)
                return a < c ? (a, c) : (c, a);

            // Обновляем точки для следующей итерации
            (a, ya) = (b, yb);
            (b, yb) = (c, yc);
            s *= k; // Увеличиваем шаг
        }
    }

    public static (double, double) FibonacciSearch(
        Func<double, double> f,
        double a, 
        double b, 
        int n, 
        double epsilon = 0.001
    ){
        // Фибоначчи-параметры
        double phi = (1 + Math.Sqrt(5)) / 2.0; // Золотое сечение
        double s = (1 - Math.Sqrt(5)) / (1 + Math.Sqrt(5));
        double rho = 1 / (phi * (1 - Math.Pow(s, n + 1)) / (1 - Math.Pow(s, n)));
        // Первая точка d
        double d = rho * b + (1 - rho) * a;
        double yd = f(d);

        for (int i = 1; i < n; i++) {
            double c;
            if (i == n - 1) {
                // Последняя итерация с учётом ε
                c = epsilon * a + (1 - epsilon) * d;
            } else {
                // Основной шаг алгоритма
                c = rho * a + (1 - rho) * b;
            }

            double yc = f(c);
            
            // Обновляем интервал
            if (yc < yd) {
                b = d;
                d = c;
                yd = yc;
            } else {
                a = b;
                b = c;
            }

            // Обновляем rho
            rho = 1 / (phi * (1 - Math.Pow(s, n - i + 1)) / (1 - Math.Pow(s, n - i)));
        }

        return a < b ? (a, b) : (b, a);
    }

    public static (double, double) GoldenSectionSearch(
        Func<double, double> f, 
        double a, 
        double b, 
        int n
    ) {
        // Золотое сечение
        double rho = (Math.Sqrt(5) - 1) / 2.0; // φ - 1

        // Лямбда-функция для вычисления новой точки c или d
        Func<double, double, double> computePoint = (x, y) => rho * x + (1 - rho) * y;

        // Инициализация точки d
        double d = computePoint(b, a);
        double yd = f(d);

        // Лямбда-функция для обновления интервала
        Action updateInterval = () => {
            double c = computePoint(a, b);
            double yc = f(c);

            if (yc < yd) {
                b = d;
                d = c;
                yd = yc;
            } else {
                a = b;
                b = c;
            }
        };

        for (int i = 1; i < n; i++) {
            updateInterval();
        }

        return a < b ? (a, b) : (b, a);
    }



    public static (double, double, double) QuadraticFitSearch(
        Func<double, double> f, 
        double a, 
        double b, 
        double c, 
        int n
    ) {
        // Вычисление значений функции в начальных точках
        double ya = f(a);
        double yb = f(b);
        double yc = f(c);

        for (int i = 0; i < n - 3; i++) {
            // Формула для нахождения вершины параболы
            double numerator = 0.5 * (ya * (Math.Pow(b, 2) - Math.Pow(c, 2)) +
                                      yb * (Math.Pow(c, 2) - Math.Pow(a, 2)) +
                                      yc * (Math.Pow(a, 2) - Math.Pow(b, 2)));

            double denominator = ya * (b - c) + yb * (c - a) + yc * (a - b);

            double x = numerator / denominator; // Координата вершины
            double yx = f(x); // Значение функции в точке x

            // Обновление интервала в зависимости от значения yx
            if (x > b) {
                if (yx > yb) {
                    c = x;
                    yc = yx;
                } else {
                    a = b;
                    ya = yb;
                    b = x;
                    yb = yx;
                }
            } else if (x < b) {
                if (yx > yb) {
                    a = x;
                    ya = yx;
                } else {
                    c = b;
                    yc = yb;
                    b = x;
                    yb = yx;
                }
            }
        }

        return (a, b, c);
    }

    public static (double optX, double optValue, int iterations) StepByStepMethod(
        double initX,
        double h,
        double eps, 
        DelOneValueFunc func,
        int maxIterations = 1000
    ) {
        double currentX = initX;
        double newX, currentFuncValue, newFuncValue;
        int iterations = 0;
        while (Math.Abs(h) > eps && iterations < maxIterations) {
            iterations++;
            currentFuncValue = func(currentX);
            newX = currentX + h;
            newFuncValue = func(newX);

            if (newFuncValue > currentFuncValue) {
            // Если функция в новой точке больше, то мы промахнулись — уменьшаем шаг и меняем направление
                h = -h / 2;
            } else {
            // Если функция в новой точке меньше или равна, увеличиваем шаг
                h *= 1.2;
                currentX = newX;
            }
        }

        return (currentX, func(currentX), iterations);
    }

    public static OptimumResult StepByStepMethodOR(
        double x,
        double h,
        double eps,
        DelOneValueFunc func
    ) {
        double xs, fs, ft;
        int k = 0;
        while (Math.Abs(h) > eps)
        {
            k++;
            ft = func(x);
            xs = x + h;
            fs = func(xs);

            if (fs > ft)
            {
                h = -h / 2;
            }
            else { h = h * 1.2; }

            x = xs;
            //  ft = fs;
        }
        Vector xr = new Vector(1);
        xr[0] = x;
        Vector fr= new Vector(1);
        fr[0] = func(x);
        OptimumResult result = new OptimumResult(xr, fr, null, k);
        return result;
    }

    public static (double,double) StepByStepMethodDR(
        double x,
        double h, 
        double eps, 
        DelOneValueFunc func
    ) {
        double xs, fs, ft;
        int k = 0;
        while (Math.Abs(h) > eps)
        {
            k++;
            ft = func(x);
            xs = x + h;
            fs = func(xs);

            if (fs > ft)
            {
                h = -h / 2;
            }
            else { h = h * 1.2; }

            x = xs;
            //  ft = fs;
        }
        
        return (x,func(x));
    }

    public static (double, double, int) GoldenRatioMethodDR(
        double a,
        double b,
        double eps,
        DelOneValueFunc func,
        int maxIterations = 1000
    ) {
        int iterations = 0;
        double phi = (1 + Math.Sqrt(5)) / 2;
        double x1 = a + (b - a) / (phi + 1);
        double x2 = b - (b - a) / (phi + 1);
        double f1 = func(x1);
        double f2 = func(x2);
        // double resphi = 2 - phi;  
        // double x1 = a + resphi * (b - a);
        // double x2 = b - resphi * (b - a);
        // double f1 = func(x1);
        // double f2 = func(x2);

        do {
            if (f1 < f2) {
                b = x2;
                x2 = x1;
                f2 = f1;
                x1 = a + (b - a) / (phi + 1);
                // x1 = a + resphi * (b - a);
                f1 = func(x1);
            } else {
                a = x1;
                x1 = x2;
                f1 = f2;
                x2 = b - (b - a) / (phi + 1);
                // x2 = b - resphi * (b - a);
                f2 = func(x2);
            }
            iterations++;
        } while (b - a > eps && iterations < maxIterations);

        double x = (x1 + x2) / 2.0;
        return (x, func(x), iterations);
    }

    private static (double, double) ParabolaVertex(
        double x1, 
        double x2, 
        double x3, 
        DelOneValueFunc func
    ) {
        double y1 = func(x1);
        double y2 = func(x2);
        double y3 = func(x3);

        double a = (y3 - (x3*(y2 - y1) + x2*y1 - x1*y2) / (x2 - x1)) / (x3*(x3 - x1 - x2) + x1*x2);
        double b = (y2 - y1) / (x2 - x1) - a * (x1 + x2);
        double c = (x2*y1 - x1*y2) / (x2 - x1) + a*x1*x2;
        
        // Координаты вершины параболы
        double xVertex = -b / (2 * a);
        double yVertex = a*xVertex*xVertex + b*xVertex + c;

        return (xVertex, yVertex);
    }

    public static (double, double, int) SquareApproxMethodDR(
        double x1, 
        double x2, 
        double x3, 
        double eps, 
        DelOneValueFunc func, 
        out TypeExtremum type
    ) {
        double xPrev, xCurr;
        (xPrev, _) = ParabolaVertex(x1, x2, x3, func);
        int iterations = 0;

        do {
            iterations++;
            double y1 = func(x1);
            double y2 = func(x2);
            double y3 = func(x3);

            (xCurr, _) = ParabolaVertex(x1, x2, x3, func);

            if (xCurr > x1 && xCurr < x2) {
                if (y1 > y2 && y1 > y3) {
                    x3 = xCurr;
                } else {
                    x1 = xCurr;
                }
            } else if (xCurr > x2 && xCurr < x3) {
                if (y2 > y1 && y2 > y3) {
                    x1 = xCurr;
                } else {
                    x2 = xCurr;
                }
            }
        } while (Math.Abs(xPrev - xCurr) > eps);

        double a, b, c;
        double y1Final = func(x1);
        double y2Final = func(x2);
        double y3Final = func(x3);

        a = (y3Final - (x3*(y2Final - y1Final) + x2*y1Final - x1*y2Final) / (x2 - x1)) / (x3*(x3 - x1 - x2) + x1*x2);
        type = a > 0 ? TypeExtremum.Min : TypeExtremum.Max;

        return (xCurr, func(xCurr), iterations);    
    }
} 