class OptimumMoreOneArg {
    public delegate double Function(Vector x);

    public static OptimumResult RandomSearch(
        Vector startPoint,
        double eps,
        Function func,
        double h,
        int maxIter = 1000
    ) {
        Random rand = new Random();
        Vector xBest = startPoint.Copy();
        double fBest = func(xBest);
        int iter = 0;
        int pointsPerIter = 3 * startPoint.GetSize(); // 3n points per iteration

        // Функции для генерации случайных точек и вычисления лучшей точки
        Func<int, Vector> generateRandomDirection = (int dimension) => Vector.NormalizeRandom(dimension);
        Func<Vector, Vector, double, Vector> computeRandomPoint = (Vector basePoint, Vector direction, double stepSize) =>
            basePoint + direction * stepSize;
        Func<Vector, double, (Vector, double)> findBestRandomPoint = (Vector currentBest, double stepSize) => {
            Vector localBestPoint = new Vector(startPoint.GetSize());
            double localBestValue = double.MaxValue;

            for (int i = 0; i < pointsPerIter; i++) {
                Vector direction = generateRandomDirection(startPoint.GetSize());
                Vector candidatePoint = computeRandomPoint(currentBest, direction, stepSize);
                double candidateValue = func(candidatePoint);

                if (candidateValue < localBestValue) {
                    localBestValue = candidateValue;
                    localBestPoint = candidatePoint.Copy();
                }
            }

            return (localBestPoint, localBestValue);
        };

        while (Math.Abs(h) > eps && iter < maxIter) {
            // Найти локально лучшую точку среди случайных
            var (xMin, fMin) = findBestRandomPoint(xBest, h);

            // Обновление текущей лучшей точки и адаптация шага
            if (fMin < fBest) {
                xBest = xMin.Copy();
                fBest = fMin;
                h *= 1.2; // Увеличить шаг при улучшении
            } else {
                h *= 0.5; // Уменьшить шаг при отсутствии улучшений
            }

            iter++;
        }

        return new OptimumResult(xBest, null, new Vector(fBest), iter);
    }


    // Модифицированный метод случайного поиска
    public static OptimumResult ModifiedRandomSearch(
        Vector startPoint,
        Function func,
        double initialStep = 1.0,
        double eps = 1e-6,
        int maxIter = 1000
    ) {
        int n = startPoint.GetSize();
        int pointsPerIter = 3 * n;
        Vector xCurrent = startPoint.Copy();
        double fCurrent = func(xCurrent);
        double h = initialStep;
        int iter = 0;

        // Лямбда-функции
        Func<Vector, Vector, Vector> extrapolate = (Vector x, Vector best) => x + (best - x) * 2;
        Func<int, Vector> generateRandomPoint = (int dimension) => Vector.NormalizeRandom(dimension);
        Func<Vector, Vector, double, Vector> randomStep = (Vector x, Vector randomDirection, double stepSize) => x + stepSize * randomDirection;

        while (Math.Abs(h) > eps && iter < maxIter) {
            // Лямбда для поиска лучшей точки среди случайных
            var findBestPoint = () => {
                Vector bestPoint = new Vector(n);
                double bestValue = double.MaxValue;

                for (int i = 0; i < pointsPerIter; i++) {
                    Vector xNew = randomStep(xCurrent, generateRandomPoint(n), h);
                    double fNew = func(xNew);

                    if (fNew < bestValue) {
                        bestValue = fNew;
                        bestPoint = xNew.Copy();
                    }
                }

                return (bestPoint, bestValue);
            };

            // Найти лучшую точку
            var (xBest, fBest) = findBestPoint();

            if (fBest < fCurrent) {
                // Попробовать экстраполяцию
                Vector xExtra = extrapolate(xCurrent, xBest);
                double fExtra = func(xExtra);

                if (fExtra < fBest) {
                    xCurrent = xExtra;
                    fCurrent = fExtra;
                } else {
                    xCurrent = xBest;
                    fCurrent = fBest;
                }

                h *= 1.2; // Увеличить шаг
            } else {
                h *= 0.5; // Уменьшить шаг
            }

            iter++;
        }

        return new OptimumResult(xCurrent, null, new Vector(fCurrent), iter);
    }

    // Градиентные методы
    // x^(k+1) = x^k - lambda^k * grad f(x^k)
    // где lambda^k:
    // - постоянная => метод может расходиться
    // - дроблный шаг, т.е. длина шага в процессе спуска делится на некое число
    // - наискорейший спуск: lambda^k = argmin[lambda] f(x^k - lambda * grad f(x^k))

    private static Vector CalculateGradient(Vector x0, Function f, double eps, Vector xCurrent, double fCurrent)
    {
        Vector grad = new Vector(x0.GetSize());
        for (int i = 0; i < x0.GetSize(); i++) {
            Vector xTemp = xCurrent.Copy();
            xTemp[i] += eps / 2;
            grad[i] = (f(xTemp) - fCurrent) / (eps / 2);
        }

        return grad;
    }

    // простейший вариант градиента
    public static OptimumResult SimpleGradient(
        Vector x0,                     // Starting point
        Function f,                    // Target function
        double lambda = 0.1,           // Step size
        double eps = 1e-6,             
        int maxIter = 20000             
    ) {
        Vector xCurrent = x0.Copy();
        double fCurrent = f(xCurrent);
        int iter = 0;

        while (Math.Abs(lambda) > eps && iter < maxIter) {
            // Calculate gradient numerically
            Vector grad = CalculateGradient(x0, f, eps, xCurrent, fCurrent);

            // Check stopping criterion (расходится)
            if (grad.Norma1() < eps)
                break;

            // Make step
            Vector xNew = xCurrent - lambda * grad;
            double fNew = f(xNew);

            // Update current point
            xCurrent = xNew;
            fCurrent = fNew;
            iter++;
        }

        return new OptimumResult(
            xCurrent,           // Best point
            null,              // No additional info
            new Vector(fCurrent), // Function value
            iter               // Iterations count
        );
    }

    // метод изменения шага
    public static OptimumResult GradientWithStepChanging(
        Vector x0,                     // Starting point
        Function f,        // Target function
        double lambda = 0.1,                // Initial step size
        double eps = 1e-6,            
        int maxIter = 20000
    ) {
        Vector xCurrent = x0.Copy();
        double fCurrent = f(xCurrent);
        int iter = 0;

        while (Math.Abs(lambda) > eps && iter < maxIter) {
            // Calculate gradient numerically
            Vector grad = CalculateGradient(x0, f, eps, xCurrent, fCurrent);

            // Make step
            Vector xNew = xCurrent - lambda * grad;
            double fNew = f(xNew);

            // Update step size based on improvement
            if (fNew < fCurrent) {
                lambda *= 1.2; // Increase step size on improvement
            } else {
                lambda /= 2; // Decrease step size on no improvement
            }

            // Update current point
            xCurrent = xNew.Copy();
            fCurrent = fNew;
            iter++;
        }

        return new OptimumResult(
            xCurrent,           // Best point
            null,              // No additional info
            new Vector(fCurrent), // Function value
            iter               // Iterations count
        );
    }

    // метод наискорейшего спуска
    public static OptimumResult FastestDescent(
        Vector x0,                     // Starting point
        Function f,                    // Target function
        double eps = 1e-6,            
        int maxIter = 20000)          
    {
        Vector xCurrent = x0.Copy();
        double fCurrent = f(xCurrent);
        int iter = 0;

        while (iter < maxIter && Math.Abs(fCurrent) > eps)
        {
            // Calculate gradient numerically
            Vector grad = CalculateGradient(x0, f, eps, xCurrent, fCurrent);

            if (grad.Norma1() < eps)
                break;

            // Find optimal step size
            double lambda = OptimumOneArg.StepByStepMethodDR(0, 0.5, eps, hk => f(xCurrent - hk * grad)).Item1;

            // Make step
            Vector xNew = xCurrent - lambda * grad;
            double fNew = f(xNew);

            // Update current point
            xCurrent = xNew;
            fCurrent = fNew;
            iter++;
        }

        return new OptimumResult(
            xCurrent,           // Best point
            null,              // No additional info
            new Vector(fCurrent), // Function value
            iter               // Iterations count
        );
    }

    // Метод сопряженных градиентов
    public static OptimumResult ConjugateGradient(
        Vector x0,
        double h,
        double eps,
        Function func
    ) {
        Vector x = x0.Copy(); // Начальная точка
        int iter = 0; // Счётчик итераций
        int n = x.GetSize(); // Размерность задачи

        // Вспомогательные лямбда-функции
        Func<Vector, Vector> gradient = (Vector point) => -1 * Gradient(point, eps, func); // Антиградиент
        Func<Vector, double> normSquared = (Vector vec) => vec * vec; // Квадрат нормы вектора

        Vector r = gradient(x); // Начальный антиградиент r(0)
        Vector d = r.Copy(); // Начальное направление d(0)
        double sigmaNew = normSquared(r); // Квадрат нормы антиградиента
        double sigma0 = sigmaNew; // Норма начального градиента для проверки сходимости

        // Основной цикл метода
        while (iter < 10000 && sigmaNew > eps * eps * sigma0) {
            // Поиск оптимального шага α (с использованием лямбда-функции для одномерной оптимизации)
            double alpha = OptimumOneArg.StepByStepMethodDR(0, h, eps, a => func(x + a * d)).Item1;

            // Переход в новую точку
            x = x + alpha * d;

            // Вычисление нового антиградиента
            Vector rNew = gradient(x);
            double sigmaOld = sigmaNew;
            sigmaNew = normSquared(rNew);

            // Рассчёт коэффициента β
            double beta = sigmaNew / sigmaOld;

            // Обновление направления с условием рестарта
            d = (iter % n == 0 || rNew * d <= 0)
                ? rNew // Рестарт: направление равно антиградиенту
                : rNew + beta * d; // Сопряжённое направление

            r = rNew.Copy();
            iter++;
        }

        // Возврат результата
        return new OptimumResult(x, null, new Vector(func(x)), iter);
    }


    // Метод для вычисления градиента функции в точке
    private static Vector Gradient(Vector x, double eps, Function func) {
        Vector grad = new Vector(x.GetSize());

        for (int i = 0; i < x.GetSize(); i++) {
            Vector xForward = x.Copy();
            Vector xBackward = x.Copy();

            xForward[i] += eps / 2.0;
            xBackward[i] -= eps / 2.0;

            grad[i] = (func(xForward) - func(xBackward)) / eps;
        }

        return grad;
    }

    // Метод Нелдера-Мида (Nelder-Mead Simplex Method)
    public static Vector NelderMead(
        Func<Vector, double> f,
        Vector[] simplex, 
        double epsilon = 1e-6, 
        double alpha = 1.0, 
        double beta = 2.0, 
        double gamma = 0.5
    ) {
        int n = simplex.Length - 1; // Размерность задачи
        double[] y = new double[simplex.Length]; // Значения функции в точках симплекса

        // Вычисляем начальные значения функции
        for (int i = 0; i < simplex.Length; i++)
            y[i] = f(simplex[i]);

        double delta = double.MaxValue;

        while (delta > epsilon) {
            // Сортируем симплекс по значениям функции
            var sortedIndices = y.Select((value, index) => 
                            new { value, index })
                                .OrderBy(x => x.value)
                                .Select(x => x.index)
                                .ToArray();

            // Переназначаем точки симплекса согласно отсортированным индексам
            simplex = sortedIndices.Select(i => simplex[i]).ToArray();
            y = sortedIndices.Select(i => y[i]).ToArray();

            // Находим наилучшую, худшую и вторую худшую точки
            Vector xl = simplex[0];  // Лучшая точка
            Vector xh = simplex[simplex.Length - 1];  // Худшая точка
            Vector xs = simplex[simplex.Length - 2];  // Вторая худшая точка
            Vector xm = new Vector(xl.Size);

            // Вычисляем центроид для первых n точек
            for (int i = 0; i < n; i++) {
                xm = xm.Add(simplex[i]);
            }
            xm = xm * (1.0 / n);

            // Шаг отражения
            Vector xr = xm + (xm - xh) * alpha;
            double yr = f(xr);

            if (yr < y[0]) {
                // Шаг расширения
                Vector xe = xm + (xr - xm) * beta;
                double ye = f(xe);

                if (ye < yr) {
                    simplex[simplex.Length - 1] = xe;
                    y[simplex.Length - 1] = ye;
                } else {
                    simplex[simplex.Length - 1] = xr;
                    y[simplex.Length - 1] = yr;
                }
            } else if (yr >= y[simplex.Length - 2]) {
                if (yr < y[simplex.Length - 1]) {
                    simplex[simplex.Length - 1] = xr;
                    y[simplex.Length - 1] = yr;
                }

                // Шаг сжатия
                Vector xc = xm + (xh - xm) * gamma;
                double yc = f(xc);

                if (yc >= y[simplex.Length - 1]) {
                    // Уменьшаем симплекс, сдвигая все точки ближе к лучшей
                    for (int i = 1; i < simplex.Length; i++) {
                        simplex[i] = xl + (simplex[i] - xl) * 0.5;
                        y[i] = f(simplex[i]);
                    }
                } else {
                    simplex[simplex.Length - 1] = xc;
                    y[simplex.Length - 1] = yc;
                }
            } else {
                simplex[simplex.Length - 1] = xr;
                y[simplex.Length - 1] = yr;
            }

            // Вычисляем стандартное отклонение значений функции
            delta = Math.Sqrt(y.Select(v => Math.Pow(v - y.Average(), 2)).Average());
        }

        // Возвращаем точку с минимальным значением функции
        return simplex[0];
    }

    public static OptimumResult NelderMeadeMethod(
        Vector initialGuess,
        double stepSize, 
        double epsilon, 
        Func<Vector, double> f, 
        double alpha = 1.0, 
        double beta = 2.0, 
        double gamma = 0.5
    ) {
        int dimension = initialGuess.Size;
        int maxIterations = 10000; // Ограничение на количество итераций
        int iterationCount = 0;

        // Инициализация начального симплекса
        Vector[] simplex = new Vector[dimension + 1];
        simplex[0] = initialGuess;
        for (int i = 1; i <= dimension; i++) {
            Vector point = initialGuess.Copy();
            point[i - 1] += stepSize;
            simplex[i] = point;
        }

        double[] y = simplex.Select(v => f(v)).ToArray();
        double delta = double.MaxValue;

        while (delta > epsilon && iterationCount < maxIterations) {
            iterationCount++;

            // Сортировка точек симплекса по значениям функции
            var sortedIndices = y.Select((value, index) =>
                            new { value, index })
                                .OrderBy(x => x.value)
                                .Select(x => x.index)
                                .ToArray();

            simplex = sortedIndices.Select(i => simplex[i]).ToArray();
            y = sortedIndices.Select(i => y[i]).ToArray();

            Vector xl = simplex[0]; // Лучшая точка
            Vector xh = simplex[simplex.Length - 1]; // Худшая точка
            Vector xs = simplex[simplex.Length - 2]; // Вторая худшая точка
            Vector xm = new Vector(dimension);

            // Вычисляем центроид для первых n точек
            for (int i = 0; i < dimension; i++)
                xm = xm.Add(simplex[i]);
            xm = xm * (1.0 / dimension);

            // Шаг отражения
            Vector xr = xm + (xm - xh) * alpha;
            double yr = f(xr);

            if (yr < y[0]) {
                // Шаг расширения
                Vector xe = xm + (xr - xm) * beta;
                double ye = f(xe);

                if (ye < yr) {
                    simplex[simplex.Length - 1] = xe;
                    y[simplex.Length - 1] = ye;
                } else {
                    simplex[simplex.Length - 1] = xr;
                    y[simplex.Length - 1] = yr;
                }
            } else if (yr >= y[simplex.Length - 2]) {
                if (yr < y[simplex.Length - 1]) {
                    simplex[simplex.Length - 1] = xr;
                    y[simplex.Length - 1] = yr;
                }

                // Шаг сжатия
                Vector xc = xm + (xh - xm) * gamma;
                double yc = f(xc);

                if (yc >= y[simplex.Length - 1]) {
                    // Уменьшение симплекса
                    for (int i = 1; i < simplex.Length; i++) {
                        simplex[i] = xl + (simplex[i] - xl) * 0.5;
                        y[i] = f(simplex[i]);
                    }
                } else {
                    simplex[simplex.Length - 1] = xc;
                    y[simplex.Length - 1] = yc;
                }
            } else {
                simplex[simplex.Length - 1] = xr;
                y[simplex.Length - 1] = yr;
            }

            // Вычисляем стандартное отклонение значений функции
            delta = Math.Sqrt(y.Select(v => Math.Pow(v - y.Average(), 2)).Average());
        }

        // Возвращаем результат
        Vector optimalPoint = simplex[0];
        double optimalValue = y[0];
        Vector fopt = new Vector(new double[] { optimalValue });

        return new OptimumResult(optimalPoint, null, fopt, iterationCount);
    }

    public static (Vector x, double fValue) Newton(
        Vector x, 
        double eps, 
        Func<Vector, double> f
    ) {
        int k = 0;

        while (true) {
            Vector grad = Gradient(x, f); // Вычисление градиента
            Matrix hess = Hesse(x, f);    // Вычисление гессиана
            Matrix hessInv = Matrix.MatrixInverse(hess); // Обратная матрица Гессиана

            if (hessInv == null) {
                Console.WriteLine("Матрица Гессе не обратима в точке: " + x);
                return (x, f(x)); // Возвращаем текущую точку как результат
            }

            Vector direction = hessInv * grad; // Направление оптимизации
            double hopt = 1.0; // Для простоты: фиксированный шаг (можно заменить на line search)
            Vector x1 = x - hopt * direction;

            Vector dx = x1 - x;
            if (dx.Norma2() < eps) {  // Проверка сходимости
                return (x1, f(x1)); // Возвращаем точку минимума и значение функции
            }

            x = x1; // Обновляем текущую точку
            k++;
        }
    }

    public static Vector Gradient(
        Vector x, 
        Func<Vector, double> f
    ) {
        int n = x.Size;
        Vector grad = new Vector(n);
        double delta = 1e-5; // Малое значение для вычисления численных производных

        for (int i = 0; i < n; i++) {
            Vector x1 = x.Copy();
            Vector x2 = x.Copy();

            x1[i] += delta;
            x2[i] -= delta;

            grad[i] = (f(x1) - f(x2)) / (2 * delta); // Центральная разностная схема
        }

        return grad;
    }


    public static Matrix Hesse(
        Vector x, 
        Func<Vector, double> f
    ) {
        int n = x.Size;
        Matrix hess = new Matrix(n, n);
        double delta = 1e-5; // Малое значение для вычисления численных производных

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                // Создаём копии вектора для сдвига по координатам
                Vector x1 = x.Copy();
                Vector x2 = x.Copy();
                Vector x3 = x.Copy();
                Vector x4 = x.Copy();

                // Сдвиги по координатам i и j
                x1[i] += delta; x1[j] += delta;
                x2[i] += delta; x2[j] -= delta;
                x3[i] -= delta; x3[j] += delta;
                x4[i] -= delta; x4[j] -= delta;

                // Вычисление смешанной частной производной
                hess[i, j] = (f(x1) - f(x2) - f(x3) + f(x4)) / (4 * delta * delta);
            }
        }

        return hess;
    }


    // Линейное изменение шага
    public static double[] LineSearch(
        Func<double[], double> f, 
        double[] x, 
        double[] d
    ) {
        // Определяем целевую функцию как функцию от α
        Func<double, double> objective = alpha => {
            double[] temp = new double[x.Length];
            for (int i = 0; i < x.Length; i++) {
                temp[i] = x[i] + alpha * d[i];
            }
            return f(temp);
        };

        // Нахождение интервала (a, b), содержащего минимум
        var (a, b) = OptimumOneArg.BracketMinimum(objective);

        // Нахождение значения α, минимизирующего objective
        double alpha = OptimumOneArg.Minimize(a, b, objective);

        // Вычисляем новую точку x + α * d
        double[] result = new double[x.Length];
        for (int i = 0; i < x.Length; i++) {
            result[i] = x[i] + alpha * d[i];
        }

        return result;
    }

    // Метод сопряженных градиентов с обновлением Polak-Ribière
    private double[] d; // Направление спуска
    private double[] g; // Текущий градиент

    private void Init(
        Func<double[], double> f,
        Func<double[], double[]> gradF, 
        double[] x
    ) {
        g = gradF(x);           // Вычисляем начальный градиент
        d = Negate(g);          // Направление спуска: -g
    }

    private static double[] Negate(double[] v) {
        double[] result = new double[v.Length];
        for (int i = 0; i < v.Length; i++)
            result[i] = -v[i];
        return result;
    }

    public double[] Step(
        Func<double[], double> f, 
        Func<double[], double[]> gradF, 
        double[] x
    ) {
        var gPrev = g;          // Сохраняем текущий градиент
        g = gradF(x);           // Вычисляем новый градиент

        // Вычисляем коэффициент β
        double beta = Math.Max(0, DotProduct(g, Subtract(g, gPrev)) / DotProduct(gPrev, gPrev));

        // Вычисляем новое направление: d' = -g + β * d
        d = Add(Negate(g), Scale(beta, d));

        // Выполняем линейный поиск для нахождения x'
        double[] xNew = LineSearch(f, x, d);

        return xNew;            // Возвращаем новое значение x
    }

    // Вспомогательные функции для операций с векторами
    private static double[] Add(double[] v1, double[] v2) {
        double[] result = new double[v1.Length];
        for (int i = 0; i < v1.Length; i++)
            result[i] = v1[i] + v2[i];
        return result;
    }

    private static double[] Subtract(double[] v1, double[] v2) {
        double[] result = new double[v1.Length];
        for (int i = 0; i < v1.Length; i++)
            result[i] = v1[i] - v2[i];
        return result;
    }

    private static double[] Scale(double scalar, double[] v) {
        double[] result = new double[v.Length];
        for (int i = 0; i < v.Length; i++)
            result[i] = scalar * v[i];
        return result;
    }

    private static double DotProduct(double[] v1, double[] v2) {
        double result = 0;
        for (int i = 0; i < v1.Length; i++)
            result += v1[i] * v2[i];
        return result;
    }
}