class LinearFn {

    // решение СЛАУ методом Гаусса
    public static Vector GaussMethod(Matrix matrix, Vector vector) {
        //проверка на корректность входных данных
        if (matrix.Columns != vector.Size) {
            throw new ArgumentException("Несоответствие количества аргументов и вектора");
        };
        for (int i = 0; i < matrix.Rows - 1; i++) {
            int j = i;
            double mainDiagonalElement = matrix[i,j];
            if (mainDiagonalElement == 0) {
                throw new ArgumentException("Главная диагональ имеет 0");
            }
            if (mainDiagonalElement != 1) {
                for (int k = 0; k < matrix.Columns; k++) {
                    matrix[i, k] /= mainDiagonalElement;
                }
                //изменение элемента вектора
                vector[j] /= mainDiagonalElement;
            }
            // прямой ход
            for (int k = i + 1; k < matrix.Rows; k++) {
                Vector currentRow = matrix.GetRow(k);
                for (int p = 0; p < matrix.Columns; p++) {
                    if (matrix[k, p] != 0) {
                        double newElement = matrix[k, p] - matrix[k, j] * matrix[i, p]; 
                        currentRow.SetElement(newElement, p);
                    }
                }
                vector[k] = vector[k] - matrix[k, j] * vector[i];
                matrix.SetRow(k, currentRow);
            }

        }

        //Получение решения после приведения к верхней треугольной матрице
        return Matrix.Upper(matrix, vector);
    }

    //решение СЛАУ методом последовательных приближений
    public static Vector SuccessiveApproximationMethod(Matrix A, Vector B) {
        if (A.Rows != A.Columns) return null;
        if (A.Rows != B.Size) return null;

        int n = A.Rows;
        double eps = 0.00000001;
        int max;
        double tmp;

        for (int j = 0; j < n; j++)
        {
            max = j;
            for (int i = j + 1; i < n; i++)
            {
                if (Math.Abs(A[i, j]) > Math.Abs(A[max, j])) { max = i; };
            }

            if (max != j)
            {
                Vector temp = A.GetRow(max); A.SetRow(max, A.GetRow(j)); A.SetRow(j, temp);
                tmp = B[max]; B[max] = B[j]; B[j] = tmp;

            }
            if (Math.Abs(A[j, j]) < eps) return null;

        }

        Vector beta = new Vector(n);
        for (int i = 0; i < beta.Size; i++)
        {
            beta[i] = B[i] / A[i, i];
        }

        Matrix alpha = new Matrix(n, n);
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                if (i != j) alpha[i, j] = A[i, j] / A[i, i];
                else alpha[i, j] = 0;
            }
            beta[i] = B[i] / A[i, i];
        }

        Vector prev_x = beta;
        Vector current_x;
        Vector delta;

        do
        {
            current_x = beta - alpha * prev_x;
            delta = current_x - prev_x;
            prev_x = current_x;
        }
        while (delta.Norma1() > eps);

        return prev_x;
    }

    // метод квадратных корней
    public static Vector SquareRootsMethod(Matrix matrix, Vector vector) {
        //проверка на корректность входных данных
        if (matrix.Columns != vector.Size) {
            throw new ArgumentException("Несоответствие количества аргументов и вектора");
        };
        int n = matrix.Rows;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (matrix[i,j] != matrix[j,i]) {
                    throw new ArgumentException("Матрица не симметричаная");
                }
            }
        }
        // L одно и то же, что T'
        var L = new Matrix(new double[n, n]);
        var y = new Vector(new double[n]);
        // вектор ответов
        var x = new Vector(new double[n]);

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < i + 1; j++) {
                double sum = 0.0;
                if (j == i) {
                    for (int k = 0; k < j; k++) {
                        sum += Math.Pow(L[i, k], 2);
                    }
                    L[i, j] = Math.Sqrt(matrix[i, i] - sum);
                }
                else {
                    for (int k = 0; k < j; k++) {
                        sum += (L[i, k] * L[j, k]);
                    }
                    L[i, j] = (matrix[i, j] - sum) / L[j, j];
                }
            }
        }

        // Решение Ly = b
        for (int i = 0; i < n; i++) {
            double sum = 0.0;
            for (int j = 0; j < i; j++) {
                sum += L[i, j] * y[j];
            }
            y[i] = (vector[i] - sum) / L[i, i];
        }

        // Решение L^T x = y
        // то есть транспонируем матрицу L(T'), получаем матрицу T и решаем ее
        var LT = L.Trans();
        for (int i = n - 1; i >= 0; i--) {
            double sum = 0.0;
            for (int j = i + 1; j < n; j++) {
                sum += LT[i, j] * x[j];
            }
            x[i] = (y[i] - sum) / LT[i, i];
        }

        return x;
    }

    public static Vector ReverseHod(Matrix L, Vector Y) {
        int n = L.Rows;
        Vector X = new(n);
        for (int i = n - 1; i >= 0; i--) {
            double sum = 0.0;
            for (int j = i + 1; j < n; j++)
                sum += L[i, j] * X[j];
            X[i] = (Y[i] - sum) / L[i, i];
        }

        return X;
    }

    // метод прогонки (для трехдиагональной матрицы)
    public static Vector ProgonkaMethod(Matrix matrix, Vector vector) {
        //проверка на корректность входных данных
        if (matrix.Columns != vector.Size) {
            throw new ArgumentException("Несоответствие количества аргументов и вектора");
        };
        
        int n = matrix.Rows;
        Vector result = new Vector(n);

        double[] alpha = new double[n];
        double[] beta = new double[n];
        double y, c, d;
        /* b1 c1 0  0
           a2 b2 c2 0
           0  a3 b3 c3
           0  0  a4 b4  */
        // y1 = b1 (b - элемент главной диагонали)
        // alpha1 = -c1 / y1 (c1 - диагональ выше главной)
        // beta1 = d1 / y1 (d1 - элемент вектора)
        y = matrix[0,0];
        c = matrix[0,1];
        d = vector[0];
        alpha[0] = -c / y;
        beta[0] = d / y;

        for (int i = 1; i < n; i++) {
            // yi = bi + ai*alpha{i-1}
            // alphai = -ci / yi
            // betai = (di - alphai*beta{i-1}) / yi
            y = matrix[i, i] + matrix[i, i - 1] * alpha[i - 1];
            if (i + 1 < n) {
                alpha[i] = -matrix[i, i + 1] / y;   
            }
            beta[i] = (vector[i] - matrix[i, i - 1] * beta[i - 1]) / y;
        }
        result[n - 1] = beta[n - 1];

        for (int i = n - 2; i >= 0; i--) {
            result[i] = alpha[i] * result[i + 1] + beta[i];
        }

        return result;
    }    

    public static void GramSchmidt(Matrix A, out Matrix Q, out Matrix R) {
        int m = A.Rows;
        int n = A.Columns;

        Q = new Matrix(m, n);
        R = new Matrix(n, n);

        // Итерация по каждому столбцу матрицы A
        for (int j = 0; j < n; j++) {
            // Вычисление j-ого столбца матрицы Q
            Vector temp = A.GetColumn(j);

            // Ортогонализация текущего столбца относительно предыдущих столбцов
            for (int k = 0; k < j; k++) {
                // Скалярное произведение q_k и a_j
                double dotProduct = Q.GetColumn(k) * A.GetColumn(j);

                // Вычитание проекции текущего столбца на уже ортогонализированные столбцы
                temp -= dotProduct * Q.GetColumn(k);
            }

            // Нормализация j-ого столбца матрицы Q
            double norm = Math.Sqrt(temp * temp);
            Vector normalizedTemp = new Vector(temp.Size);
            for (int i = 0; i < temp.Size; i++) {
                normalizedTemp[i] = temp[i] / norm;
            }
            Q.SetColumn(j, normalizedTemp);

            // Вычисление j-ой строки матрицы R
            for (int k = 0; k <= j; k++) {
                // Скалярное произведение q_k и a_j
                double dotProduct = Q.GetColumn(k) * A.GetColumn(j);
                R[k, j] = dotProduct;
            }
        }
    }

    public static Vector Solve_GramSchmidt(Matrix A, Vector B) {
        GramSchmidt(A, out var Q, out var R);
        return ReverseHod(R, Q.Trans() * B);
    }
}