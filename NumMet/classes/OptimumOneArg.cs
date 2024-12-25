
enum TypeExtremum {
    None, Min, Max
}
class OptimumOneArg {
    public delegate double DelOneValueFunc(double x);

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
        DelOneValueFunc func)
    {
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
        DelOneValueFunc func)
    {
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

//     private static double Derivative1(double x, double eps, DelOneValueFunc func)
//     {
//         return (func(x + eps) - func(x)) / eps;
//     }

//     private static double Derivative2(double x, double eps, DelOneValueFunc func)
//     {
//         return (Derivative1(x + eps, eps, func) - Derivative1(x, eps, func)) / eps;
//     }

//     // Модифицированный метод Ньютона для нахождения экстремума (может расходится)
//     // Вычисляю вторую производную только в начале и потом ее использую
//     public static (double, double) NewtonMethodDR(double x, double eps, DelOneValueFunc func, out TypeExtr type)
//     {
//         double d1 = Derivative1(x, eps, func);
//         double d2 = Derivative2(x, eps, func);

//         if (d2 > 0)
//         {
//             type = TypeExtr.Min;
//         }
//         else
//         {
//             type = TypeExtr.Max;
//         }

//         double x_prev = x;
//         double x_curr = x_prev - d1 / d2;
//         double diff_prev = x_curr - x_prev;
//         do
//         {
//             d1 = Derivative1(x_curr, eps, func);
//             x_prev = x_curr;
//             x_curr = x_prev - d1 / d2;
//             double diff_curr = x_curr - x_prev;
//             if (Math.Abs(diff_curr) > Math.Abs(diff_prev))
//             {
//                 throw new Exception("NewtonMethodDR: расходится");
//             }
//             diff_prev = diff_curr;
//         } while (Math.Abs(x_prev - x_curr) > eps);

//         return (x_curr, func(x_curr));
//     }
} 