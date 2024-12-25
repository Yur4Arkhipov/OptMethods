class Integral {
    public delegate double DoubleIntegralFunc(double x, double y);
    public static double RectangleMethod(double a, double b, double eps, Func<double, double> f) {
        double previousFunctionValue, currentFunctionValue, sum;
        double diference; // разность между текущим и предыдущим значенем функции
        int n = 1; // количество интервлов в начале
        double h = (b - a) / n; // шаг разбиения
        sum = f(a); //начальное значение функции
        previousFunctionValue = sum * h;

        do {
            n *= 2; 
            h = (b - a) / n; 
            for (int i = 1; i < n; i += 2)
                sum += f(a + i * h);
        
            // вычисление нового приближенного значения интеграла
            currentFunctionValue = sum * h;

            diference = currentFunctionValue - previousFunctionValue;
            previousFunctionValue = currentFunctionValue;
        } while (Math.Abs(diference) > eps);

        return previousFunctionValue;
    }

    /// Метод трапеций для вычисления определенного интеграла.
    /// a - Нижний предел интегрирования.
    /// b - Верхний предел интегрирования.
    /// e - Точность вычисления.
    /// f - Функция, для которой вычисляется интеграл.
    /// Значение определенного интеграла.
    public static double TrapezoidMethod(double a, double b, double eps, Func<double, double> f) {
        int n = 1;
        double h = (b - a) / n;
        double sum = f(a) + f(b);
        double previousFunctionValue = sum * h / 2;
        double difference;

        do {
            n *= 2;
            h = (b - a) / n;
            sum = f(a) + f(b); // включаем в сумму значения функции в точках a и b
            for (int i = 1; i < n; i++) {
                sum += 2 * f(a + i * h);
            }
            difference = sum * (h / 2) - previousFunctionValue;
            previousFunctionValue = sum * (h / 2);
        } while (Math.Abs(difference) > eps);

        return previousFunctionValue;
    }

    public static double SimpsonMethod(double a, double b, double eps, Func<double, double> f) {
        int n = 1;
        double h, previousFunctionValue = 1, y1, y2 = 0;
        double y0 = f(a) + f(b);
        double difference;

        do {
            n *= 2;
            h = (b - a) / n;
            difference = previousFunctionValue;
            y1 = 0;

            for (int i = 1; i < n; i += 2) 
                y1 += f(a + h * i);

            previousFunctionValue = h / 3 * (y0 + 2 * y2 + 4 * y1);
            y2 = y1 + y2;
            difference -= previousFunctionValue;
        } while (Math.Abs(difference) > eps);

        return previousFunctionValue;
    }   

    // Двойной интеграл
    public static double DoubleIntegral(double a, double b, double c, double d, double eps, DoubleIntegralFunc function) {
        double intervalCount = 1;
        double intervalWidthX = (b - a) / intervalCount;
        double intervalWidthY = (d - c) / intervalCount;
        double previousSum;
        double currentSum = 0;

        do {
            previousSum = currentSum;
            intervalCount *= 2;
            intervalWidthX = (b - a) / intervalCount;
            intervalWidthY = (d - c) / intervalCount;
            currentSum = 0;

            for (int i = 0; i <= intervalCount; i++) {
                for (int j = 0; j <= intervalCount; j++) {
                    double x = a + i * intervalWidthX;
                    double y = c + j * intervalWidthY;
                    double weightX = (i == 0 || i == intervalCount) ? 1 : (i % 2 == 0) ? 2 : 4;
                    double weightY = (j == 0 || j == intervalCount) ? 1 : (j % 2 == 0) ? 2 : 4;
                    currentSum += weightX * weightY * function(x, y);
                }
            }

            currentSum = currentSum * intervalWidthX * intervalWidthY / 9;
        } while (Math.Abs(previousSum - currentSum) > eps);

        return currentSum;
    }
}