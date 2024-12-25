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
        Vector xBest = startPoint;
        double fBest = func(xBest);
        int iter = 0;
        // 3n points per iteration, where n is the number of dimensions
        int pointsPerIter = 3 * startPoint.GetSize(); 

        while (Math.Abs(h) > eps) {
            Vector xMin = new Vector(startPoint.GetSize());
            double fMin = double.MaxValue;
            
            // Generate and test random points
            for (int i = 0; i < pointsPerIter; i++) {
                // Generate normalized random direction
                Vector direction = Vector.NormalizeRandom(startPoint.GetSize());
                Vector xNew = xBest + direction * h;
                double fNew = func(xNew);
                
                if (fNew < fMin) {
                    fMin = fNew;
                    xMin = xNew.Copy();
                }
            }
            
            // Update best point and adapt step size
            if (fMin < fBest) {
                xBest = xMin.Copy();
                fBest = fMin;
                h *= 1.2; // Increase step size on improvement
            }
            else {
                h *= 0.5; // Decrease step size on no improvement
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
        int maxIter = 1000)
    {
        int n = startPoint.GetSize();
        int pointsPerIter = 3 * n;
        Vector xCurrent = startPoint.Copy();
        double fCurrent = func(xCurrent);
        double h = initialStep;
        int iter = 0;

        while (Math.Abs(h) > eps && iter < maxIter)
        {
            // Generate and evaluate random points
            Vector xBest = new Vector(n);
            double fBest = double.MaxValue;
            
            for (int i = 0; i < pointsPerIter; i++)
            {
                Vector xNew = xCurrent + h * Vector.NormalizeRandom(n);
                double fNew = func(xNew);
                
                if (fNew < fBest)
                {
                    fBest = fNew;
                    xBest = xNew.Copy();
                }
            }

            // If improvement found, try extrapolation
            if (fBest < fCurrent)
            {
                // Extrapolation step
                Vector xExtra = xCurrent + (xBest - xCurrent) * 2;
                double fExtra = func(xExtra);

                if (fExtra < fBest)
                {
                    xCurrent = xExtra;
                    fCurrent = fExtra;
                }
                else
                {
                    xCurrent = xBest;
                    fCurrent = fBest;
                }
                h *= 1.2; // Increase step size
            }
            else
            {
                h *= 0.5; // Decrease step size
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
        for (int i = 0; i < x0.GetSize(); i++)
        {
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
    ){
        Vector xCurrent = x0.Copy();
        double fCurrent = f(xCurrent);
        int iter = 0;

        while (Math.Abs(lambda) > eps && iter < maxIter)
        {
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
        int maxIter = 20000)          
    {
        Vector xCurrent = x0.Copy();
        double fCurrent = f(xCurrent);
        int iter = 0;

        while (Math.Abs(lambda) > eps && iter < maxIter)
        {
            // Calculate gradient numerically
            Vector grad = CalculateGradient(x0, f, eps, xCurrent, fCurrent);

            // Make step
            Vector xNew = xCurrent - lambda * grad;
            double fNew = f(xNew);

            // Update step size based on improvement
            if (fNew < fCurrent)
            {
                lambda *= 1.2; // Increase step size on improvement
            }
            else
            {
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
        Vector r = new Vector(x.GetSize()); // Антиградиент
        Vector d = new Vector(x.GetSize()); // Направление спуска

        int iter = 0; // Счётчик итераций
        int n = x.GetSize(); // Размерность задачи
        double sigmaNew, sigmaOld, sigma0;

        // Вычисление начального антиградиента r(0) и направления d(0)
        r = -1 * Gradient(x, eps, func);
        d = r.Copy();

        sigmaNew = r * r; // Квадрат нормы антиградиента
        sigma0 = sigmaNew; // Норма начального градиента для проверки сходимости

        while (iter < 10000 && sigmaNew > eps * eps * sigma0)
        {
            // Цикл одномерной минимизации (поиск оптимального шага a)
            double alpha = OptimumOneArg.StepByStepMethodDR(0, h, eps, a => func(x + a * d)).Item1;

            // Переход в новую точку
            x = x + alpha * d;

            // Вычисление нового антиградиента
            Vector rNew = -1 * Gradient(x, eps, func);
            sigmaOld = sigmaNew;
            sigmaNew = rNew * rNew;

            // Проверка на рестарт (каждые n шагов или если направление теряет положительность)
            double beta = sigmaNew / sigmaOld;
            if (iter % n == 0 || rNew * d <= 0)
            {
                d = rNew;
            }
            else
            {
                d = rNew + beta * d; // Обновление сопряжённого направления
            }

            r = rNew.Copy();
            iter++;
        }

        // Возврат результата
        Vector f = new Vector(func(x));
        return new OptimumResult(x, null, f, iter);
    }

    // Метод для вычисления градиента функции в точке
    private static Vector Gradient(Vector x, double eps, Function func)
    {
        Vector grad = new Vector(x.GetSize());

        for (int i = 0; i < x.GetSize(); i++)
        {
            Vector xForward = x.Copy();
            Vector xBackward = x.Copy();

            xForward[i] += eps / 2.0;
            xBackward[i] -= eps / 2.0;

            grad[i] = (func(xForward) - func(xBackward)) / eps;
        }

        return grad;
    }

}