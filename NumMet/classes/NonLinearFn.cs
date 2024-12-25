class NonLinearFn {
    public static double BisectionMethod(double leftBorder, double rightBorder, double eps, Func<double,double> func) {
        int countIteration = 0;
        // leftBorder – левая граница отрезка
        // rightBorder – правая граница отрезка
        // eps – требуемая точность
        // Func<double,double> func – делегат функции func(double x) которая возвращает значение double
        double f_left= func(leftBorder);
        double f_right= func(rightBorder);
        if(f_left * f_right > 0)
            throw new ArgumentException();
            // нет корня на отрезке
        double middle; // середина отрезка
        do {
            middle = (leftBorder + rightBorder)/2;
            double funcFromMiddle = func(middle);
            // значение функции в середине отрезка
            if (f_left * funcFromMiddle < 0) {
                rightBorder = middle;
            } else {
                leftBorder = middle;
                f_left = funcFromMiddle;
            }
            countIteration++;
        } while (Math.Abs(rightBorder - leftBorder) > eps);
        Console.WriteLine($"Count of iteration: {countIteration}");
        return middle;  // результат
    }

    public static double Newton(double x0, double eps, Func<double, double> func) {
        double delta = eps / 2;         
        double xCurrent = x0;           //current value (Xn)
        double dx = double.MaxValue;    //initial distance between xCurrent and xNext
        int countIteration = 0;

        // Вычисление начального значения x1
        double diffForX1 = (func(xCurrent + delta) - func(xCurrent)) / delta;
        double xNext = xCurrent - (func(xCurrent) / diffForX1);
        // Вычисление начального значения x1

        do {
            xCurrent = xNext;
            double diff = (func(xCurrent + delta) - func(xCurrent)) / delta; //differential formula   
            xNext = xCurrent - (func(xCurrent) / diff);     //main problem solve (Newton method)                         

            //Проверка на расходимость ряда
            double newDx = xNext - xCurrent;                //new distance between xCurrent and xNext
            if (Math.Abs(newDx) > Math.Abs(dx)) {
                return double.NaN;                  //return NaN (ряд расходится)
            }
            dx = newDx;
            //Проверка на расходимость ряда 
            countIteration++;

        } while (Math.Abs(xNext - xCurrent) > eps);
        
        Console.WriteLine($"Count of iteration: {countIteration}");
        return xCurrent;
    }

    public static double SimpleIteration(double x0, double eps, Func<double, double> func) {
        double xCurrent = x0;
        double xNext;
        double dx = double.MaxValue;
        double newDx;
        int countIteration = 0;

        do {
            xNext = func(xCurrent);

            //условие сходимости
            newDx = xCurrent - xNext;
            if (Math.Abs(newDx) > Math.Abs(dx)) 
                return double.NaN;
            dx = newDx;
            //условие сходимости
            
            xCurrent = xNext;

            countIteration++;

        } while (Math.Abs(dx) > eps);

        Console.WriteLine($"Count of iteration: {countIteration}");
        return xCurrent;
    }

    //метод хорд
    public static double ChordMethod(double a, double b, double eps, Func<double,double> func) {
        double fa = func(a);
        double fb = func(b);
        if (fa * fb > 0) return double.NaN;
        while (Math.Abs(b - a) > eps) {
            a -= (b - a) * func(a) / (func(b) - func(a));
            b -= (a - b) * func(b) / (func(a) - func(b));            
        }
        return (a + b) / 2;
    }
}