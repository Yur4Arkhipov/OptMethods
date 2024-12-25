class MinRect {
    public delegate double Fn(double x);
    private Fn[] fns; // Array of function delegates

    private Vector x; // Vector of x values
    private Vector y; // Vector of y values

    public Vector param; // вектор параметров

    public MinRect(Vector x, Vector y, Fn[] fns) {
        if (x.Size != y.Size) {
            throw new Exception("Размеры векторов не совпадают");
        }

        this.x = x;
        this.y = y;
        this.fns = fns;
    }

    public void ComputeParameters() {
        int n = x.Size; // количество точек
        int m = fns.Length; // количество функций
        Matrix A = new Matrix(n, m); // создаем матрицу A размером n x m
        for (int i = 0; i < n; i++) {
            A.SetRow(i, GetVector(x[i])); // устанавливаем каждую строку матрицы A в вектор значений функций в x[i]
        }
        param = LinearFn.SquareRootsMethod(A.Trans() * A, A.Trans() * y); // вычисляем параметры методом квадратных корней
    }


    // создает вектор, содержащий значения каждой функции в указанной точке
    // точка, в которой вычисляются функции
    // вектор, содержащий значения каждой функции в указанной точке
    private Vector GetVector(double x) {
        Vector result = new Vector(fns.Length); // создаем вектор такой же длины, как и количество функций
        for (int i = 0; i < fns.Length; i++) {
            result[i] = fns[i](x); // вычисляем значение каждой функции в точке x и записываем в вектор
        }
        return result;
    }

    public double GetError() {
        Vector result = new Vector(x.Size);
        for (int i = 0; i < x.Size; i++) {
            result[i] = y[i] - param * GetVector(x[i]); // Compute the error for each data point and store it in the vector
        }
        return result.Norma1();
    }
}