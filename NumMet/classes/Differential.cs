delegate Vector PravDU(double t,Vector x);
class Differential {
    // вычисление дифференциальных уравнений
    public delegate Vector Derivative(double t, Vector x);
    // 
   /* два параметра: t (время) и x (вектор).
    возвращает новый вектор result, который является результатом вычисления дифференциальных уравнений.
    В данном случае, функция решает простое дифференциальное уравнение второго порядка. 
    Она использует значения вектора x для вычисления новых значений, которые затем записываются в result.
    Конкретно, функция выполняет следующие действия:
    Создает новый вектор result с размерностью 2.
    Присваивает первый элемент result[0] значение x[1].
    Присваивает второй элемент result[1] значение -x[0].
    Возвращает вектор result.
    Таким образом, функция DifferentialEquations вычисляет производные второго порядка для заданного времени t и вектора x.*/
    public static Vector DifferentialEquations(double t, Vector x)
    {
        Vector result = new Vector(2);
        result[0] = x[1];
        result[1] = -x[0];
        return result;
    }

    public static Matrix EulerMethod(double initialTime, double finalTime, Vector initialState, int numSteps, Derivative derivative)
    {
        // определение размерности вектора состояния
        int dimension = initialState.Size;
        // создание матрицы результатов
        Matrix result = new Matrix(dimension + 1, numSteps + 1);
        // вычисление шага
        double stepSize = (finalTime - initialTime) / numSteps;
        // создание вектора-столбца для хранения времени и состояния
        Vector column = new Vector(dimension + 1);
        column[0] = initialTime;

        // заполнение первого столбца матрицы результатов
        for (int i = 0; i < dimension; i++)
            column[i + 1] = initialState[i];

        result.SetColumn(0, column);
        Vector xt = initialState.Copy();
        double t = initialTime;
        Vector derivativeValue;

        // выполнение метода Эйлера для каждого шага
        for (int k = 1; k <= numSteps; k++)
        {
            // вычисление производной в текущей точке
            derivativeValue = derivative(t, xt);
            // вычисление нового состояния
            xt = xt + derivativeValue * stepSize;
            // обновление времени
            t += stepSize;
            column[0] = t;

            // заполнение столбца матрицы результатов
            for (int i = 0; i < dimension; i++)
                column[i + 1] = xt[i];

            result.SetColumn(k, column);
        }

        return result;
    }


    public static Matrix RungeKutta2Method(double initialTime, double finalTime, Vector initialState, int numSteps, Derivative derivative)
    {
        // определение размерности вектора состояния
        int dimension = initialState.Size;
        // создание матрицы результатов
        Matrix result = new Matrix(dimension + 1, numSteps + 1);
        // вычисление шага
        double stepSize = (finalTime - initialTime) / numSteps;
        // создание вектора-столбца для хранения времени и состояния
        Vector column = new Vector(dimension + 1);
        column[0] = initialTime;

        // заполнение первого столбца матрицы результатов
        for (int i = 0; i < dimension; i++)
            column[i + 1] = initialState[i];

        result.SetColumn(0, column);
        Vector xt = initialState.Copy();
        double t = initialTime;

        // выполнение метода Рунге-Кутты второго порядка для каждого шага
        for (int k = 1; k <= numSteps; k++)
        {
            // вычисление k1
            Vector k1 = derivative(t, xt);
            // вычисление k2
            Vector k2 = derivative(t + stepSize, xt + k1 * stepSize);
            // вычисление нового состояния xt
            xt = xt + (k1 + k2) * (stepSize / 2.0);
            // обновление времени t
            t += stepSize;
            column[0] = t;

            // заполнение столбца матрицы результатов
            for (int i = 0; i < dimension; i++)
                column[i + 1] = xt[i];

            result.SetColumn(k, column);
        }

        return result;
    }

    public static Matrix RungeKutta4Method(double initialTime, double finalTime, Vector initialState, int numSteps, Derivative derivative)
    {
        // определение размерности вектора состояния
        int dimension = initialState.Size;
        // создание матрицы результатов
        Matrix result = new Matrix(dimension + 1, numSteps + 1);
        // вычисление шага
        double stepSize = (finalTime - initialTime) / numSteps;
        // создание вектора-столбца для хранения времени и состояния
        Vector column = new Vector(dimension + 1);
        column[0] = initialTime;

        // заполнение первого столбца матрицы результатов
        for (int i = 0; i < dimension; i++)
            column[i + 1] = initialState[i];

        result.SetColumn(0, column);
        Vector xt = initialState.Copy();
        double t = initialTime;

        // выполнение метода Рунге-Кутты четвертого порядка для каждого шага
        for (int k = 1; k <= numSteps; k++)
        {
            // вычисление k1
            Vector k1 = derivative(t, xt);
            // вычисление k2
            Vector k2 = derivative(t + stepSize / 2.0, xt + k1 * (stepSize / 2.0));
            // вычисление k3
            Vector k3 = derivative(t + stepSize / 2.0, xt + k2 * (stepSize / 2.0));
            // вычисление k4
            Vector k4 = derivative(t + stepSize, xt + k3 * stepSize);
            // вычисление нового состояния xt
            xt = xt + (k1 + 2 * k2 + 2 * k3 + k4) * (stepSize / 6.0);
            // обновление времени t
            t += stepSize;
            column[0] = t;

            // заполнение столбца матрицы результатов
            for (int i = 0; i < dimension; i++)
                column[i + 1] = xt[i];

            result.SetColumn(k, column);
        }

        return result;
    }

    public static Matrix AdamsMethod(double initialTime, double finalTime, Vector initialState, int numSteps, Derivative derivative)
    {
        // определение размерности вектора состояния
        int dimension = initialState.Size;
        // создание матрицы результатов
        Matrix result = new Matrix(dimension + 1, numSteps + 1);
        // вычисление шага
        double stepSize = (finalTime - initialTime) / numSteps;
        // создание вектора-столбца для хранения времени и состояния
        Vector column = new Vector(dimension + 1);
        column[0] = initialTime;

        // заполнение первого столбца матрицы результатов
        for (int i = 0; i < dimension; i++)
            column[i + 1] = initialState[i];

        result.SetColumn(0, column);
        Vector[] xt = new Vector[numSteps + 1];
        double[] t = new double[numSteps + 1];
        Vector[] derivatives = new Vector[numSteps + 1];

        t[0] = initialTime;
        xt[0] = initialState.Copy();
        derivatives[0] = derivative(t[0], xt[0]);

        // выполнение метода Адамса для первых 3 шагов
        for (int k = 1; k <= Math.Min(3, numSteps); k++)
        {
            // вычисление времени t[k]
            t[k] = t[k - 1] + stepSize;
            // вычисление производной в текущей точке
            derivatives[k] = derivative(t[k], xt[k - 1]);
            // вычисление нового состояния xt[k]
            xt[k] = xt[k - 1] + derivatives[k] * stepSize;
            // обновление времени в столбце матрицы результатов
            column[0] = t[k];

            // заполнение столбца матрицы результатов
            for (int i = 0; i < dimension; i++)
                column[i + 1] = xt[k][i];

            result.SetColumn(k, column);
        }

        // выполнение метода Адамса для остальных шагов
        for (int k = 3; k < numSteps; k++)
        {
            // вычисление времени t[k+1]
            t[k + 1] = t[k] + stepSize;
            // вычисление нового состояния xt[k+1] с использованием формулы Адамса
            xt[k + 1] = xt[k] + (stepSize / 24) * (55 * derivatives[k] - 59 * derivatives[k - 1] + 37 * derivatives[k - 2] - 9 * derivatives[k - 3]);
            // вычисление производной в новой точке
            derivatives[k + 1] = derivative(t[k + 1], xt[k + 1]);
            // обновление времени в столбце матрицы результатов
            column[0] = t[k + 1];

            // заполнение столбца матрицы результатов
            for (int i = 0; i < dimension; i++)
                column[i + 1] = xt[k + 1][i];

            result.SetColumn(k + 1, column);
        }

        return result;
    }
}