class Matrix {
    protected int rows, columns;
    protected double[,] data;
    public Matrix(int r, int c) {
        this.rows = r; 
        this.columns = c;
        data = new double[rows, columns];
        for (int i = 0; i < rows; i++)
            for (int j = 0; j < columns; j++) data[i, j] = 0;
    }
    public Matrix(double[,] mm) {
        this.rows = mm.GetLength(0);
        this.columns = mm.GetLength(1);
        data = new double[rows, columns];
        for (int i = 0; i < rows; i++)
            for (int j = 0; j < columns; j++)
                data[i, j] = mm[i, j];
    }

    public Matrix(Vector low, Vector sr, Vector up) {
        if (low.Size != sr.Size || sr.Size != up.Size)
            throw new ArgumentException("Несоответствие количества аргументов и вектора");
        int n = low.Size;
        this.rows = n;
        this.columns = n;
        data = new double[n, n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (i == j) {
                    data[i, j] = sr[i];
                }
                else if (i == j - 1) {
                    data[i, j] = low[i];
                }
                else if (i == j + 1) {
                    data[i, j] = up[i];
                }
                else {
                    data[i, j] = 0;
                }
            }
        }
    }
    public int Rows { get { return rows; } }
    public int Columns { get { return columns; } }

    public double this[int i, int j] {
        get {
            if (i < 0 && j < 0 && i >= rows && j >= columns) {
                // Console.WriteLine(" Индексы вышли за пределы матрицы ");
                return Double.NaN;
            } else return data[i, j];
        }
        set {
            if (i < 0 && j < 0 && i >= rows && j >= columns) {
                //Console.WriteLine(" Индексы вышли за пределы матрицы ");
            } else data[i, j] = value;
        }
    }
    public Vector GetRow(int r) {
        if (r >= 0 && r < rows)
        {
            Vector row = new Vector(columns);
            for (int j = 0; j < columns; j++) row[j] = data[r, j];
            return row;
        }
        return null;
    }
    public Vector GetColumn(int c) {
        if (c >= 0 && c < columns)
        {
            Vector column = new Vector(rows);
            for (int i = 0; i < rows; i++) column[i] = data[i, c];
            return column;
        }
        return null;
    }
    public double Norma() {
        double sum = 0.0;
        for (int i = 0; i < rows; i++)
        {
            for (int j = 0; j < columns; j++)
            {
                sum += data[i, j] * data[i, j];
            }
        }
        return Math.Sqrt(sum);
    }
    public bool SetRow(int index, Vector r) {
        if (index < 0 || index > rows) return false;
        if (r.Size != columns) return false;
        for (int k = 0; k < columns; k++) data[index, k] = r[k];
        return true;
    }
    public bool SetColumn(int index, Vector c) {
        if (index < 0 || index > columns) return false;
        if (c.Size != rows) return false;
        for (int k = 0; k < rows; k++) data[k, index] = c[k];
        return true;
    }
    public void SwapRows(int r1, int r2) {
        if (r1 < 0 || r2 < 0 || r1 >= rows || r2 >= rows || (r1 == r2)) return;
        Vector v1 = GetRow(r1);
        Vector v2 = GetRow(r2);
        SetRow(r2, v1);
        SetRow(r1, v2);
    }
    public Matrix Copy() {
        Matrix r = new Matrix(rows, columns);
        for (int i = 0; i < rows; i++)
            for (int j = 0; j < columns; j++) r[i, j] = data[i, j];
        return r;
    }
    public Matrix Trans() {
        Matrix transposeMatrix = new Matrix(columns, rows);
        for (int i = 0; i < rows; i++)
        {
            for (int j = 0; j < columns; j++)
            {
                transposeMatrix.data[j, i] = data[i, j];
            }
        }
        return transposeMatrix;
    }
    //Сложение матриц
    public static Matrix operator +(Matrix m1, Matrix m2) {
        if (m1.rows != m2.rows || m1.columns != m2.columns)
        {
            throw new Exception("Матрицы не совпадают по размерности");
        }
        Matrix result = new Matrix(m1.rows, m1.columns);

        for (int i = 0; i < m1.rows; i++)
        {
            for (int j = 0; j < m2.columns; j++)
            {
                result[i, j] = m1[i, j] + m2[i, j];
            }
        }
        return result;
    }
    
    //Вычитание матриц
    public static Matrix operator -(Matrix m1, Matrix m2) {
        if (m1.rows != m2.rows || m1.columns != m2.columns)
        {
            throw new Exception("Матрицы не совпадают по размерности");
        }
        Matrix result = new Matrix(m1.rows, m2.columns);

        for (int i = 0; i < m1.rows; i++)
        {
            for (int j = 0; j < m2.columns; j++)
            {
                result[i, j] = m1[i, j] - m2[i, j];
            }
        }
        return result;
    }
    
    //Умножение матриц
    public static Matrix operator *(Matrix m1, Matrix m2) {
        if (m1.columns != m2.rows)
        {
            throw new Exception("Количество столбцов первой матрицы не равно количеству строк второй");
        }
        Matrix result = new Matrix(m1.rows, m2.columns);
        for (int i = 0; i < m1.rows; i++)
        {
            for (int j = 0; j < m2.columns; j++)
            {
                for (int k = 0; k < m1.columns; k++)
                {
                    result[i, j] += m1[i, k] * m2[k, j];
                }
            }
        }
        return result;
    }
   
   //умножение матрицы на вектор
    public static Vector operator *(Matrix a, Vector b) {
        if (a.columns != b.Size) return null;
        Vector r = new Vector(a.rows);
        for (int i = 0; i < a.rows; i++)
        {
            r[i] = a.GetRow(i) * b;
        }
        return r;
    }
   
    //Умножение матрицы на число
    public static Matrix operator *(Matrix a, double number) {
        for (int i = 0; i < a.rows; i++)
        {
            for (int j = 0; j < a.columns; j++) a[i, j] *= number;
        }
        return a;
    }

    //Умножение числа на матрицу
    public static Matrix operator *(double number, Matrix a) {   
        for (int i = 0; i < a.rows; i++)
        {
            for (int j = 0; j < a.columns; j++) a[i, j] *= number;
        }
        return a;

    }

    //Решение СЛАУ с нижней треугольной матрицей
    public static Vector Lower(Matrix a, Vector b) {
        if (a.rows != b.Size)
            throw new Exception("Матрицы и вектор не совпадают по размерности");
        for (int i = 0; i < a.rows; i++)
            if (a[i, i] == 0)
                throw new Exception("Матрицы null diaganal");

        Vector x = new Vector(b.Size);
        x[0] = b[0] / a[0, 0];
        for (int i = 1; i < a.rows; i++)
        {
            double c = 0;
            for (int j = 0; j < i; j++) { c += x[j] * a[i, j]; }
            x[i] = (b[i] - c) / (a[i, i]);
        }
        return x;
    }

    //Решение СЛАУ с верхней треугольной матрицей
    public static Vector Upper(Matrix a, Vector b) {
        if (a.rows != b.Size)
            throw new Exception("Матрицы и вектор не совпадают по размерности");
        for (int i = 0; i < a.rows; i++)
            if (a[i, i] == 0)
                throw new Exception("Матрицы null diaganal");

        Vector x = new Vector(b.Size);
        x[b.Size - 1] = b[b.Size - 1] / a[a.rows - 1, a.columns - 1];
        for (int i = a.rows - 2; i >= 0; i--)
        {
            double c = 0;
            for (int j = b.Size - 1; j > i; j--) { c += a[i, j] * x[j]; }
            x[i] = (b[i] - c) / (a[i, i]);
        }
        return x;
    }

    //Выыод матрицы
    public static void PrintMatrix(Matrix matrix, string nameMatrix = " ") {
        Console.WriteLine("Matrix {0}:", nameMatrix);
        for (int i = 0; i < matrix.Rows; i++)
        {
            for (int j = 0; j < matrix.Columns; j++)
            {
                Console.Write("{0}\t", matrix[i, j]);
            }
            Console.WriteLine("");
        }
        Console.Write('\n');
    }
}
