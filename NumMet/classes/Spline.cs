class Spline {
    Vector xVector, yVector;
    int countPoints, countIntervals;
    double[] splineParameters;

    public Spline(Vector x, Vector y) {
        this.xVector = x;
        this.yVector = y;
        if (x.Size == y.Size) {
            this.countPoints = x.Size;
            this.countIntervals = countPoints - 1;
            this.splineParameters = new double[countPoints];
        }
    }

    public void solveParameters() {
        // Инициализация начальных значений параметров сплайна
        splineParameters[0] = 0.0;
        splineParameters[countIntervals] = 0.0;
        
        // Расчет расстояния между точками
        double[] h = new double[countIntervals];
        for (int i = 0; i < countIntervals; i++) {
            h[i] = xVector[i + 1] - xVector[i];
        }
        
        // Инициализация вспомогательных векторов и матрицы
        Vector low = new Vector(countIntervals - 1);
        Vector middle = new Vector(countIntervals - 1);
        Vector up = new Vector(countIntervals - 1);
        Vector res = new Vector(countIntervals - 1);
        Vector b = new Vector(countIntervals - 1);

        // Расчет значений вспомогательных векторов и вектора b
        for (int i = 0; i < countIntervals - 1; i++) {
            if (i == 0)
                low[0] = 0.0;
            if (i > 0)
                low[i] = h[i] / 6.0;
            middle[i] = (h[i] + h[i + 1]) / 3.0;
            if (i < (countIntervals - 2)) 
                up[i] = h[i + 1] / 6.0;
            else
                up[i] = 0.0;
            b[i] = (yVector[i + 2] - yVector[i + 1]) / h[i + 1] - (yVector[i + 1] - yVector[i]) / h[i];
        }
        
        // Создание матрицы для прогонки
        Matrix progonkaMatrix = new Matrix(low, middle, up);
        
        // Решение системы уравнений методом прогонки
        res = LinearFn.ProgonkaMethod(progonkaMatrix, b);
        
        // Заполнение массива параметров сплайна
        for (int i = 1; i < countIntervals; i++) {
            splineParameters[i] = res[i - 1];
        }
    }

    public double getValue(double x) {
        int i;

        // Проверка, находится ли x в диапазоне значений xVector
        if (x < xVector[0] || x > xVector[countIntervals]) return 0;

        // Нахождение интервала, в котором находится x
        for (i = 0; i < countIntervals; i++) {
            if (x >= xVector[i] && x <= xVector[i + 1]) break;
        }

        // Расчет разности между соседними значениями x
        double h = xVector[i + 1] - xVector[i];

        // Инициализация результата
        double result = 0;

        // Расчет значения сплайна в точке x
        result += splineParameters[i + 1] * Math.Pow(x - xVector[i], 3) / (6 * h);
        result += splineParameters[i] * Math.Pow(xVector[i + 1] - x, 3) / (6 * h);
        result += (yVector[i + 1] - splineParameters[i + 1] * h * h / 6) * (x - xVector[i]) / h;
        result += (yVector[i] - splineParameters[i] * h * h / 6) * (xVector[i + 1] - x) / h;

        // Возвращение результата
        return result;
    }
}
