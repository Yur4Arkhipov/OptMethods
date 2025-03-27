// #define NonLinear
// #define LinearAlgEquations
// #define Spline
// #define MinRect
// #define Integral
// #define Differential
// #define OptimumOneArg 
// #define OptimumMoreOneArg
// #define LinearProgramming
// #define Graph
#define CreditOptimizer


#if NonLinear //Non-linear alg
double root1 = NonLinearFn.BisectionMethod(1.0, 2.0, 0.0001, x => x*x*x - x - 2);
Console.WriteLine($"BisectionMethod / Root: {root1}\n");
double root2 = NonLinearFn.Newton(1.5, 0.001, Math.Cos);
Console.WriteLine($"Newton / Root: {root2}\n");
double root3 = NonLinearFn.SimpleIteration(0.5, 0.001, x => x * x);
Console.WriteLine($"SimpleIteration / Root: {root3}\n");
double root4 = NonLinearFn.ChordMethod(8, 3, 0.001, x => x * x * x - 18*x - 83);
Console.WriteLine($"ChordMethod / Root: {root4}\n");

#elif LinearAlgEquations //Systems lianear alg equations
double[,] matrix = {{3,2,-5},
                    {2,-1,3},
                    {1,2,-1}};
double[] vector = {-1,13,9};
            
double[,] matrixDop = {{1,3,-2, 0, -2},
                    {3,4,-5,1,-3},
                    {-2,-5,3,-2,2},
                    {0,1,-2,5,3},
                    {-2,-3,2,3,4}};
double[] vectorDop = {0.5,5.4,5.0,7.5,3.3};

double[,] matrixToSquare = {{2,1,4},
                            {1,1,3},
                            {4,3,14}};
double[] vectorToSquare = {16,12,52};

double[,] matrixToProgonka = {{2,-1,0},
                              {5,4,2},
                              {0,1,-3}};
double[] vectorToProgonka = {3,6,2};

double[,] matrixToSquare2 = {{4,1},
                            {1,0.3472}};
double[] vectorToSquare2 = {20,6.1666};


Matrix A1 = new Matrix(matrix);
Vector B1 = new Vector(vector);
Matrix A2 = new Matrix(matrix);
Vector B2 = new Vector(vector);
Matrix A3 = new Matrix(matrixToSquare);
Vector B3 = new Vector(vectorToSquare);
Vector BDop = new Vector(vectorDop);
Matrix ADOp = new Matrix(matrixDop);
Vector progonkaVector = new Vector(vectorToProgonka);
Matrix progonkaMatrix = new Matrix(matrixToProgonka);
Matrix matrixSquare2 = new Matrix(matrixToSquare2);
Vector vectorSquare2 = new Vector(vectorToSquare2);

Matrix.PrintMatrix(A1, "draft");
System.Console.Write("Vector: ");
foreach (var item in vector) {
    System.Console.Write(item + " ");
}
System.Console.WriteLine();

System.Console.WriteLine("Gauss: ");
Vector Gauss = LinearFn.GaussMethod(A1, B1);
Console.WriteLine(Gauss.ToString());
System.Console.WriteLine("Successive approximation: ");
Vector SuccessiveApproximationMethod = LinearFn.SuccessiveApproximationMethod(A2, B2);
Console.WriteLine(SuccessiveApproximationMethod.ToString());
System.Console.WriteLine("Square roots:");
Vector SquareRootsMethod = LinearFn.SquareRootsMethod(A3, B3);
Console.WriteLine(SquareRootsMethod.ToString());
System.Console.WriteLine("Progonka:");
Vector ProgonkaMethod = LinearFn.ProgonkaMethod(progonkaMatrix, progonkaVector);
Console.WriteLine(ProgonkaMethod.ToString());
System.Console.WriteLine("Square roots method 2");
Vector SquareRootsMethod2 = LinearFn.SquareRootsMethod(matrixSquare2, vectorSquare2);
System.Console.WriteLine(SquareRootsMethod2.ToString());

Matrix A = new(3, 3);
A.SetRow(0, new Vector(new double[] {1, 5, 9}));
A.SetRow(1, new Vector(new double[] {6, 2, 1}));
A.SetRow(2, new Vector(new double[] {7, 7, 4}));
Matrix Q, R;
LinearFn.GramSchmidt(A, out Q, out R);
Matrix.PrintMatrix(Q, "Q:");
Matrix.PrintMatrix(R, "R:");

double[,] matrixToGramShmidt = {{2,-1,0},
                                {5,4,2},
                                {0,1,-3}};
double[] vectorToGramShmidt = {3,6,2};
Matrix gramShmidtMatrix = new Matrix(matrixToGramShmidt);
Vector gramShmidtVector = new Vector(vectorToGramShmidt);

Vector result_gaus = LinearFn.GaussMethod(gramShmidtMatrix, gramShmidtVector);
Vector result_gramSchmidt = LinearFn.Solve_GramSchmidt(gramShmidtMatrix, gramShmidtVector);
Console.WriteLine("Result gauss: ");
Console.WriteLine(result_gaus);
Console.WriteLine("Result gramSchmidt: ");
Console.WriteLine(result_gramSchmidt);

#elif Spline
Vector xx = new Vector(new double[] { 1, 2.2, 3, 5, 7});
Vector yy = new Vector(new double[] { 6, 4.3, 9.8, 1.6, 2.9});
Spline spline = new Spline(xx, yy);
spline.solveParameters();
List<double> arrX = new List<double>(1);
List<double> arrY = new List<double>(1);
for (double i = 1; i <= 7; i += 0.1) {
    arrX.Add(i);
}
for (double valueX = 1; valueX <= 7; valueX += 0.1) {
    double valueY = spline.getValue(valueX);
    arrY.Add(valueY);
    // Console.WriteLine($" Value x: {valueX} Value y: {valueY}");
}
System.Console.WriteLine($"X: {string.Join(",", arrX)}");
System.Console.WriteLine();
System.Console.WriteLine($"Y: {string.Join(",", arrY)}");


System.Console.WriteLine(arrX.Count);
System.Console.WriteLine(arrY.Count);


#elif MinRect //МНК
Vector x1 = new Vector(new double[] { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 });
Vector y1 = new Vector(new double[] { 0, 14, 18, 20, 22, 24, 26, 28, 30, 32, 34 });
MinRect.Fn[] fns1 = new MinRect.Fn[] { x => x, x => x * x };

MinRect kv1 = new MinRect(x1, y1, fns1);
kv1.ComputeParameters();
Console.WriteLine("Параметры: {0}", kv1.param);
Console.WriteLine("Ошибка: {0}", kv1.GetError());

Vector x2 = new Vector(new double[] { 2, 4, 6, 12 });
Vector y2 = new Vector(new double[] { 8, 5.25, 3.5, 3.25});
MinRect.Fn[] fns2 = new MinRect.Fn[] { x => 1, x => 1 / x };

MinRect kv2 = new MinRect(x2, y2, fns2);
kv2.ComputeParameters();
Console.WriteLine("Параметры: {0}", kv2.param);
Console.WriteLine("Ошибка: {0}", kv2.GetError());

#elif Integral //Интегралы
Console.WriteLine("\nИнтегралы");
System.Console.WriteLine("Функция: x^2");
System.Console.WriteLine("Нижний предел: -1");
System.Console.WriteLine("Верхний предел: 1");
Console.WriteLine("Метод прямоугольника");
double resultRect = Integral.RectangleMethod(-1, 1, 0.001, x => x * x * x);
System.Console.WriteLine($"Результат: {resultRect}");

Console.WriteLine("Метод трапеций");
double resultTrap = Integral.TrapezoidMethod(-1, 1, 0.001, x => x * x * x);
System.Console.WriteLine($"Результат: {resultTrap}");

Console.WriteLine("Метод Симпсона");
double resultSimps = Integral.SimpsonMethod(-1, 1, 0.001, x => x * x * x);
System.Console.WriteLine($"Результат: {resultSimps}");

Console.WriteLine("Двойной интеграл");
double resultDoubleIntegral = Integral.DoubleIntegral(0, 1, 0, 2, 0.001, (x, y) => x * y * y);
System.Console.WriteLine($"Результат: {resultDoubleIntegral}");

#elif Differential //Дифференциальные уравнения
Console.WriteLine("\nДифференциальные уравнения");
Console.WriteLine("tn = 0");
Console.WriteLine("tk = 1");
Console.WriteLine("xn[0] = 0, xn[1] = 1");
Console.WriteLine("m = 10");
Vector xn = new Vector(new double[] { 0, 1 });
Console.WriteLine("Аналитическое решение");

for (double t = 0; t < 1.0; t += 0.1) {
    Console.WriteLine("t={0} Sin(t)={1} Cos(t)={2}", Math.Round(t,4), Math.Round(Math.Sin(t), 4), Math.Round(Math.Cos(t), 4));
}

Console.WriteLine("\nМетод Эйлера");
Matrix Euler = Differential.EulerMethod(0, 1, xn, 10, D);
for (int i = 0; i < Euler.Rows; i++) {
    for (int j = 0; j < Euler.Columns; j++) {
        Euler[i, j] = Math.Round(Euler[i, j], 4);
    }
}
Console.WriteLine($"Ответ:");
Matrix.PrintMatrix(Euler, "Эйлер");

Console.WriteLine("\nМетод Рунге - Кутты 2-го порядка");
Matrix RK2 = Differential.RungeKutta2Method(0, 1, xn, 10, D);
RK2 = roundMatrix(RK2);
Console.WriteLine($"Ответ:");
Matrix.PrintMatrix(RK2, "Рунге - Кутты 2");

Console.WriteLine("\nМетод Рунге - Кутты 4-го порядка");
Matrix RK4 = Differential.RungeKutta4Method(0, 1, xn, 10, D);
RK4 = roundMatrix(RK4);
Console.WriteLine($"Ответ:");
Matrix.PrintMatrix(RK4, "Рунге - Кутты 4");

Console.WriteLine("\nМетод Адамса");
Matrix Adams = Differential.AdamsMethod(0, 1, xn, 10, D);
Adams = roundMatrix(Adams);
Console.WriteLine($"Ответ:");
Matrix.PrintMatrix(Adams, "Адамс");

static Vector D(double t, Vector zx) {
    return new Vector([zx[1], -zx[0]]);
}

static Matrix roundMatrix(Matrix m) {
    for (int i = 0; i < m.Rows; i++) {
        for (int j = 0; j < m.Columns; j++) {
            m[i, j] = Math.Round(m[i, j], 4);
        }
    }
    return m;
}

#elif OptimumOneArg

double xopt, fopt;
int iterations;

// Метод bracketMinimum определения границ интервала, содержащего экстремум
System.Console.WriteLine("BracketMinimum:");
{
    var (start1, end1) = OptimumOneArg.BracketMinimum(x => Math.Pow(x - 3, 2));
    Console.WriteLine($"func = (x-3)^2 >>> Интервал минимума: [{start1}, {end1}]");
    var (start2, end2) = OptimumOneArg.BracketMinimum(x => Math.Pow(x + 2, 2) + 6);
    Console.WriteLine($"func = (x+2)^2 + 6 >>> Интервал минимума: [{start2}, {end2}]");
}
System.Console.WriteLine();

// Метод Фидбоначчи
System.Console.WriteLine("Fibonacci method:");
{
    Func<double, double> f1 = x => Math.Pow(x - 3, 2);
    // Интервал поиска
    double a1 = 0;
    double b1 = 6;
    // Число итераций
    int n1 = 100;
    // Поиск минимума
    var (start1, end1) = OptimumOneArg.FibonacciSearch(f1, a1, b1, n1);
    Console.WriteLine($"func = (x-3)^2 >>> Интервал минимума: [{start1}, {end1}]");
    Func<double, double> f2 = x => Math.Pow(x +2, 2) + 6;
    // Интервал поиска
    double a2 = -5;
    double b2 = 10;
    // Число итераций
    int n2 = 100;
    // Поиск минимума
    var (start2, end2) = OptimumOneArg.FibonacciSearch(f2, a2, b2, n2);
    Console.WriteLine($"func = (x+2)^2 + 6 >>> Интервал минимума: [{start2}, {end2}]");
}
System.Console.WriteLine();

// Метод Золотого сечения
System.Console.WriteLine("GoldenSectionSearch: ");
{
    // Функция f(x) = (x - 3)^2
    Func<double, double> f1 = x => Math.Pow(x - 3, 2);
    double a1 = 0;
    double b1 = 6;
    int n1 = 100;

    var (start1, end1) = OptimumOneArg.GoldenSectionSearch(f1, a1, b1, n1);
    Console.WriteLine($"func = (x-3)^2 >>> Интервал минимума: [{start1}, {end1}]");

    // Функция f(x) = sin(x) в интервале [2, 4]
    Func<double, double> f2 = x => Math.Sin(x);
    double a2 = 2;
    double b2 = 4;
    int n2 = 15;

    var (start2, end2) = OptimumOneArg.GoldenSectionSearch(f2, a2, b2, n2);
    Console.WriteLine($"func = sin(x) >>> Интервал минимума: [{start2}, {end2}]");

    // Функция f(x) = exp(x-2)-x в интервале [-2, 6]
    Func<double, double> f3 = x => Math.Exp(x-2)-x;
    double a3 = -2;
    double b3 = 6;
    int n3 = 30;

    var (start3, end3) = OptimumOneArg.GoldenSectionSearch(f3, a3, b3, n3);
    Console.WriteLine($"func = exp(x-2)-x [-2, 6] >>> Интервал минимума: [{start3}, {end3}]");
}
System.Console.WriteLine();

// QuaranticFitSearch (Метод квадратичной аппроксимации)
System.Console.WriteLine("QuadranticFitSearch:");
{ 
    // Функция f(x) = (x - 3)^2
    Func<double, double> f1 = x => Math.Pow(x - 3, 2);
    double a1 = 0;
    double b1 = 2;
    double c1 = 4;
    int n1 = 100;

    var (start1, mid1, end1) = OptimumOneArg.QuadraticFitSearch(f1, a1, b1, c1, n1);
    Console.WriteLine($"func = (x-3)^2 >>> Интервал минимума: ({start1}, {mid1}, {end1})");

    // Функция f(x) = sin(x) в интервале [0, π]
    Func<double, double> f2 = x => Math.Sin(x);
    double a2 = 0;
    double b2 = Math.PI / 2;
    double c2 = Math.PI;
    int n2 = 20;

    var (start2, mid2, end2) = OptimumOneArg.QuadraticFitSearch(f2, a2, b2, c2, n2);
    Console.WriteLine($"func = sin(x) [0, π] >>> Интервал минимума: ({start2}, {mid2}, {end2})");
}
System.Console.WriteLine();

// Метод может искать только min
Console.WriteLine("StepByStepMethod: ");
(xopt, fopt, iterations) = OptimumOneArg.StepByStepMethod(1.0, 0.5, 0.000001, x => x * x + x);
Console.WriteLine("func = x * x + x >>> x = {0}, f(x) = {1}, iterations = {2}", xopt, fopt, iterations);
Console.WriteLine("StepByStepMethodOR: ");
Console.WriteLine("func = x * x + x >>> {0}", OptimumOneArg.StepByStepMethodOR(1.0, 0.5, 0.000001, x => x * x + x)); 
Console.WriteLine("StepByStepMethodDR: ");
(xopt,fopt) = OptimumOneArg.StepByStepMethodDR(1.0, 0.5, 0.000001, x => x * x + x);
Console.WriteLine("func = x * x + x >>> x = {0}  f(x) = {1}", xopt, fopt);
System.Console.WriteLine();

// Метод может искать только min
Console.WriteLine("GoldenRatioMethod: ");
(xopt, fopt, iterations) = OptimumOneArg.GoldenRatioMethodDR(-10.0, 10.0, 0.000001, x => x * x + x);
Console.WriteLine("func = x * x + x >>> x = {0} f(x) = {1}, iterations = {2}", xopt, fopt, iterations);
(xopt, fopt, iterations) = OptimumOneArg.GoldenRatioMethodDR(-10.0, 10.0, 0.000001, x => 2*x * x - 3*x + 1);
Console.WriteLine("func = 2x * x - 3x  + 1 >>> x = {0} f(x) = {1}, iterations = {2}", xopt, fopt, iterations);
System.Console.WriteLine();

System.Console.WriteLine("SquareApproxMethodDR:");
(xopt, fopt, iterations) = OptimumOneArg.SquareApproxMethodDR(-10.0, 0.0, 10.0, 0.000001, x => x * x + x, out TypeExtremum type1);
Console.WriteLine("func = x * x + x >>> x = {0} f(x) = {1}, type = {2}, iterations = {3}", xopt, fopt, type1, iterations);
(xopt, fopt, iterations) = OptimumOneArg.SquareApproxMethodDR(-10.0, 0.0, 10.0, 0.000001, x => 2*x * x - 3*x + 1, out TypeExtremum type2);
Console.WriteLine("func = 2x * x - 3x  + 1 >>> x = {0} f(x) = {1}, type = {2}, iterations = {3}", xopt, fopt, type2, iterations);
(xopt, fopt, iterations) = OptimumOneArg.SquareApproxMethodDR(-10.0, 0.0, 10.0, 0.000001, x => -2*x * x - 3*x + 1, out TypeExtremum type3);
Console.WriteLine("func = -2x * x - 3x  + 1 >>> x = {0} f(x) = {1}, type = {2}, iterations = {3}", xopt, fopt, type3, iterations);
// (xopt, fopt, iterations) = OptimumOneArg.SquareApproxMethodDR(-10.0, 0.0, 10.0, 0.000001, x => -0.5*x + 2, out TypeExtremum type4);
// Console.WriteLine("func = -1/2x + 2 >>> x = {0} f(x) = {1}, type = {2}, iterations = {3}", xopt, fopt, type3, iterations);

#elif OptimumMoreOneArg
static double Rosenbrock(Vector x)
{
    double a = 1 - x[0];
    double b = x[1] - x[0] * x[0];
    return a * a + 10 * b * b;
}
static double funcForGradient(Vector x)
{
    return 100 * x[0] * x[0] + x[1] * x[1];
}
{
    // Метод случайного поиска
    OptimumResult result = OptimumMoreOneArg.RandomSearch(new Vector(new double[] { 0.0, 0.0 }), 0.000001, Rosenbrock, 0.1);
    Console.WriteLine("RandomSearch func=Rosenbrock {0}", result);
    OptimumResult result2 = OptimumMoreOneArg.RandomSearch(new Vector(new double[] { 4.0, 4.0 }), 0.000001, funcForGradient, 0.1);
    Console.WriteLine("RandomSearch func=10 x^2 + y^2 {0}", result2);

    // Модифицированный метод случайного поиска
    OptimumResult result3 = OptimumMoreOneArg.ModifiedRandomSearch(new Vector(new double[] { 0.0, 0.0 }), Rosenbrock, 0.1, 0.000001);
    Console.WriteLine("ModifiedRandomSearch func=Rosenbrock {0}", result3);
    OptimumResult result4 = OptimumMoreOneArg.ModifiedRandomSearch(new Vector(new double[] { 4.0, 4.0 }), funcForGradient,  0.1, 0.000001);
    Console.WriteLine("ModifiedRandomSearch func=10 x^2 + y^2 {0}", result4);
}
System.Console.WriteLine();

// Метод Ньютона
{
    // Пример функции f(x, y) = (x-1)^2 + (y-2)^2
    Func<Vector, double> f = v => Math.Pow(v[0] - 1, 2) + Math.Pow(v[1] - 2, 2);

    Vector x0 = new Vector(new double[] { 0, 0 }); // Начальная точка
    double eps = 1e-6;

    (Vector optimum, double fValue) = OptimumMoreOneArg.Newton(x0, eps, f);

    Console.WriteLine($"Newton's method: func=(x-1)^2 + (y-2)^2 >>> Min: x={optimum}, f(x)={fValue}");
}
System.Console.WriteLine();

// Градиентные методы
{
    // Простой градиентный метод
    OptimumResult result1 = OptimumMoreOneArg.SimpleGradient(new Vector(new double[] { 0.0, 0.0 }), Rosenbrock, 0.1, 0.000001);
    OptimumResult result2 = OptimumMoreOneArg.SimpleGradient(new Vector(new double[] { 0.0, 0.0 }), Rosenbrock, 0.01, 0.000001);
    OptimumResult result3 = OptimumMoreOneArg.SimpleGradient(new Vector(new double[] { 0.0, 0.0 }), Rosenbrock, 0.001, 0.000001);
    OptimumResult result4 = OptimumMoreOneArg.SimpleGradient(new Vector(new double[] { 0.0, 0.0 }), Rosenbrock, 0.0005, 0.000001);
    Console.WriteLine("Gradient func=Rosenbrock {0}, step: 0.1", result1);
    Console.WriteLine("Gradient func=Rosenbrock {0}, step: 0.01", result2);
    Console.WriteLine("Gradient func=Rosenbrock {0}, step: 0.001", result3);
    Console.WriteLine("Gradient func=Rosenbrock {0}, step: 0.0005", result4);
    OptimumResult result5 = OptimumMoreOneArg.SimpleGradient(new Vector(new double[] { 4.0, 4.0 }), funcForGradient, 0.1, 0.000001);
    OptimumResult result6 = OptimumMoreOneArg.SimpleGradient(new Vector(new double[] { 4.0, 4.0 }), funcForGradient, 0.01, 0.000001);
    OptimumResult result7 = OptimumMoreOneArg.SimpleGradient(new Vector(new double[] { 4.0, 4.0 }), funcForGradient, 0.001, 0.000001);
    OptimumResult result8 = OptimumMoreOneArg.SimpleGradient(new Vector(new double[] { 4.0, 4.0 }), funcForGradient, 0.0005, 0.000001);
    Console.WriteLine("Gradient func=10 x^2 + y^2 {0}, step: 0.1", result5);
    Console.WriteLine("Gradient func=10 x^2 + y^2 {0}, step: 0.01", result6);
    Console.WriteLine("Gradient func=10 x^2 + y^2 {0}, step: 0.001", result7);
    Console.WriteLine("Gradient func=10 x^2 + y^2 {0}, step: 0.0005", result8);
    // Из полученных результатов можно сделать вывод, что при слишком большом чаге метод расходится,
    // при слишком малом сходится медленно и точность хуже.
    // Нужно выбирать шаг наибольшим из тех, при которых метод сходится.
    // Вывод: шаг 0.001 - оптимальный
}
System.Console.WriteLine();
{
    // Метод изменения шага
    OptimumResult result1 = OptimumMoreOneArg.GradientWithStepChanging(new Vector(new double[] { 0.0, 0.0 }), Rosenbrock, 0.01, 0.000001);
    Console.WriteLine("Gradient func=Rosenbrock {0}, step: 0.1", result1);
    OptimumResult result2 = OptimumMoreOneArg.GradientWithStepChanging(new Vector(new double[] { 4.0, 4.0 }), funcForGradient, 0.1, 0.000001);
    Console.WriteLine("Gradient func=10 x^2 + y^2 {0}, step: 0.1", result2);
}
System.Console.WriteLine();
{
    // Метод наискорейшего спуска
    OptimumResult result1 = OptimumMoreOneArg.FastestDescent(new Vector(new double[] { 0.0, 0.0 }), Rosenbrock, 0.000001);
    Console.WriteLine("FastestDescent func=Rosenbrock {0}", result1);
    OptimumResult result2 = OptimumMoreOneArg.FastestDescent(new Vector(new double[] { 4.0, 4.0 }), funcForGradient, 0.000001);
    Console.WriteLine("FastestDescent func=10 x^2 + y^2 {0}", result2);
}
System.Console.WriteLine();
{
    // Метод сопряженных градиентов
    OptimumResult result1 = OptimumMoreOneArg.ConjugateGradient(new Vector(new double[] { 0.0, 0.0 }), 0.1,  0.000001, Rosenbrock);
    Console.WriteLine("ConjugateGradients func=Rosenbrock {0}", result1);
    OptimumResult result2 = OptimumMoreOneArg.ConjugateGradient(new Vector(new double[] { 4.0, 4.0 }), 0.1,  0.000001, funcForGradient);
    Console.WriteLine("ConjugateGradients func=10 x^2 + y^2 {0}", result2);
}
System.Console.WriteLine();

// Nealder-Mead (деформированный многогранник)
{
    Func<Vector, double> Rosenbrock1 = (Vector v) => {
        double x = v.GetElement(0);
        double y = v.GetElement(1);
        return 100 * Math.Pow(y - x * x, 2) + Math.Pow(1 - x, 2);
    };

    Vector initialGuess = new Vector(new double[] { 0.0, 0.0 });
    double stepSize = 0.1;
    double epsilon = 1e-6;

    OptimumResult result = OptimumMoreOneArg.NelderMeadeMethod(initialGuess, stepSize, epsilon, Rosenbrock);
    Console.WriteLine($"NelderMeadeMethod func=Rosenbrock {result}");
}
{
    Func<Vector, double> f = (Vector v) => v.GetElement(0) * v.GetElement(0) + v.GetElement(1) * v.GetElement(1); // Пример квадратичной функции
    Vector[] simplex = new Vector[] {
        new Vector(new double[] { 1.0, 2.0 }),
        new Vector(new double[] { 1.5, 2.5 }),
        new Vector(new double[] { 2.0, 1.0 })
    };
    Vector result1 = OptimumMoreOneArg.NelderMead(f, simplex);
    Vector result2 = OptimumMoreOneArg.NelderMead(Rosenbrock, simplex);
    Console.WriteLine($"Optimized point: {result1}");
    Console.WriteLine($"Optimized point: {result2}");
}
System.Console.WriteLine();


// Методы спуска
// Метод линейного поиска (LineSearch)
System.Console.WriteLine("Densent");
{
    // Пример функции f(x) = (x₁ - 1)² + (x₂ - 2)²
    Func<double[], double> f = x => Math.Pow(x[0] - 1, 2) + Math.Pow(x[1] - 2, 2);

    double[] x = { 0.0, 0.0 }; // Начальная точка
    double[] d = { 1.0, 1.0 }; // Направление поиска

    double[] result = OptimumMoreOneArg.LineSearch(f, x, d);

    Console.WriteLine($"Новая точка: [{string.Join(", ", result)}]");
}


#elif LinearProgramming

/* System.Console.WriteLine("TransportationProblem: ");
{
    
} */


#elif Graph

System.Console.WriteLine("FordFulkerson:");
{
    Graphs.Graph g = new Graphs.Graph(isDirected: true);
    
    // Create vertices
    Graphs.Vertex a = new Graphs.Vertex("A"); // source
    Graphs.Vertex f = new Graphs.Vertex("F"); // sink
    Graphs.Vertex b = new Graphs.Vertex("B");
    Graphs.Vertex c = new Graphs.Vertex("C");
    Graphs.Vertex d = new Graphs.Vertex("D");
    Graphs.Vertex e = new Graphs.Vertex("E");
    
    // Add vertices to graph
    g.allVertexs.AddRange(new List<Graphs.Vertex> {a, f, b, c, d, e});
    
    // // Add edges with capacities
    g.allEdges.AddRange(new List<Graphs.Edge> {
        new Graphs.Edge(a, b, 7),
        new Graphs.Edge(a, c, 4),
        new Graphs.Edge(b, c, 4),
        new Graphs.Edge(b, e, 2),
        new Graphs.Edge(c, e, 8),
        new Graphs.Edge(c, d, 4),
        new Graphs.Edge(e, d, 4),
        new Graphs.Edge(e, f, 5),
        new Graphs.Edge(d, f, 12),
    });
    
    // Find maximum flow
    double maxFlow = g.FordFulkerson(a, f);
    Console.WriteLine($"Maximum flow: {maxFlow}");
}
{
    Graphs.Graph g = new Graphs.Graph(isDirected: false);
    
    // Create vertices
    Graphs.Vertex a = new Graphs.Vertex("A"); // source
    Graphs.Vertex f = new Graphs.Vertex("F"); // sink
    Graphs.Vertex b = new Graphs.Vertex("B");
    Graphs.Vertex c = new Graphs.Vertex("C");
    Graphs.Vertex d = new Graphs.Vertex("D");
    Graphs.Vertex e = new Graphs.Vertex("E");
    
    // Add vertices to graph
    g.allVertexs.AddRange(new List<Graphs.Vertex> {a, f, b, c, d, e});
    
    // Add edges with capacities
    g.allEdges.AddRange(new List<Graphs.Edge> {
        new Graphs.Edge(a, b, 10),
        new Graphs.Edge(b, a, 10),
        new Graphs.Edge(a, d, 7),
        new Graphs.Edge(d, a, 7),
        new Graphs.Edge(a, f, 6),
        new Graphs.Edge(f, a, 6),
        new Graphs.Edge(a, e, 8),
        new Graphs.Edge(e, a, 8),

        new Graphs.Edge(b, f, 9),
        new Graphs.Edge(f, b, 9),
        new Graphs.Edge(b, d, 10),
        new Graphs.Edge(d, b, 10),
        new Graphs.Edge(b, c, 8),
        new Graphs.Edge(c, b, 8),
        new Graphs.Edge(b, e, 9),
        new Graphs.Edge(e, b, 9),

        new Graphs.Edge(c, e, 10),
        new Graphs.Edge(e, c, 10),
        new Graphs.Edge(c, d, 10),
        new Graphs.Edge(d, c, 10),

        new Graphs.Edge(d, f, 12),
        new Graphs.Edge(f, d, 12),
        new Graphs.Edge(d, e, 6),
        new Graphs.Edge(e, d, 6),
    });
    
    // Find maximum flow
    double maxFlow = g.FordFulkerson(e, f);
    Console.WriteLine($"Maximum flow: {maxFlow}");
}
{
    Graphs.Graph g = new Graphs.Graph(isDirected: false);
    
    // Create vertices
    Graphs.Vertex a = new Graphs.Vertex("A"); // source
    Graphs.Vertex f = new Graphs.Vertex("F"); // sink
    Graphs.Vertex b = new Graphs.Vertex("B");
    Graphs.Vertex c = new Graphs.Vertex("C");
    Graphs.Vertex d = new Graphs.Vertex("D");
    Graphs.Vertex e = new Graphs.Vertex("E");
    
    // Add vertices to graph
    g.allVertexs.AddRange(new List<Graphs.Vertex> {a, f, b, c, d, e});
    
    // // Add edges with capacities
    g.allEdges.AddRange(new List<Graphs.Edge> {
        new Graphs.Edge(a, b, 7),
        new Graphs.Edge(b, a, 7),
        new Graphs.Edge(a, c, 4),
        new Graphs.Edge(c, a, 4),
        new Graphs.Edge(b, c, 4),
        new Graphs.Edge(c, b, 4),
        new Graphs.Edge(b, e, 2),
        new Graphs.Edge(e, b, 2),
        new Graphs.Edge(c, e, 8),
        new Graphs.Edge(e, c, 8),
        new Graphs.Edge(c, d, 4),
        new Graphs.Edge(d, c, 4),
        new Graphs.Edge(e, d, 4),
        new Graphs.Edge(d, e, 4),
        new Graphs.Edge(e, f, 5),
        new Graphs.Edge(f, e, 5),
        new Graphs.Edge(d, f, 12),
        new Graphs.Edge(f, d, 12),
    });
    
    // Find maximum flow
    double maxFlow = g.FordFulkerson(a, f);
    Console.WriteLine($"Maximum flow: {maxFlow}");
}
{
    Console.WriteLine("Поток:");
    // Создаём граф
    Graphs.Graph graph = new Graphs.Graph(isDirected: true);

    Graphs.Vertex v1 = new Graphs.Vertex("1"); // source
    Graphs.Vertex v2 = new Graphs.Vertex("2"); // sink
    Graphs.Vertex v3 = new Graphs.Vertex("3");
    Graphs.Vertex v4 = new Graphs.Vertex("4");
    Graphs.Vertex v5 = new Graphs.Vertex("5");
    Graphs.Vertex v6 = new Graphs.Vertex("6");

    graph.allVertexs.AddRange(new List<Graphs.Vertex> {v1, v2, v3, v4, v5, v6});

    // Добавляем рёбра между вершинами
    graph.AddDualEdge(v1, v2, 7); 
    graph.AddDualEdge(v1, v3, 4); 
    graph.AddDualEdge(v2, v4, 2); 
    graph.AddDualEdge(v2, v3, 4); 
    graph.AddDualEdge(v3, v4, 8);
    graph.AddDualEdge(v3, v5, 4);
    graph.AddDualEdge(v5, v4, 4);
    graph.AddDualEdge(v6, v4, 5);
    graph.AddDualEdge(v6, v5, 12);

    double maxFlow = graph.FordFulkerson(v1, v6);
    // Выводим результат
    Console.WriteLine("Максимальный поток: " + maxFlow);
}

#elif CreditOptimizer
System.Console.WriteLine("Resource Allocation Problem:");
{
    int[] funds = { 0, 20, 40, 60, 80, 100 };
    int[][] enterprisesData = 
    {
        new int[] { 0, 10, 31, 42, 62, 76 },  // Предприятие 1 (F1)
        new int[] { 0, 12, 24, 36, 52, 74 },  // Предприятие 2
        new int[] { 0, 11, 36, 45, 60, 77 },  // Предприятие 3
        new int[] { 0, 16, 37, 46, 63, 80 }   // Предприятие 4
    };

    var optimizer = new CreditOptimizer();
    var intermediateTables = optimizer.CalculateIntermediateTables(funds, enterprisesData);
    
    optimizer.PrintIntermediateTables(funds, intermediateTables);
}
System.Console.WriteLine();
System.Console.WriteLine("Resource Allocation Problem:");
{
    // Пример входных данных
    int[] funds = { 0, 1, 2, 3, 4 };
    int[][] enterprisesData = new int[][]
    {
        new int[] { 0, 10, 18, 29, 38 },  // Предприятие 1
        new int[] { 0, 11, 18, 27, 39 },    // Предприятие 2
        new int[] { 0, 9, 19, 30, 40 },    // Предприятие 3
        new int[] { 0, 12, 20, 28, 39 }     // Предприятие 4
    };

    var optimizer = new CreditOptimizer();
    var intermediateTables = optimizer.CalculateIntermediateTables(funds, enterprisesData);
    optimizer.PrintIntermediateTables(funds, intermediateTables);
}
System.Console.WriteLine();
System.Console.WriteLine("Resource Allocation Problem:");
{
    // Пример входных данных
    int[] funds = { 0, 1, 2, 3, 4 };
    int[][] enterprisesData = new int[][]
    {
        new int[] { 0, 10, 18, 29, 41 },  // Предприятие 1
        new int[] { 0, 12, 20, 28, 42 },    // Предприятие 2
        new int[] { 0, 13, 19, 30, 40 },    // Предприятие 3
    };

    var optimizer = new CreditOptimizer();
    var intermediateTables = optimizer.CalculateIntermediateTables(funds, enterprisesData);
    optimizer.PrintIntermediateTables(funds, intermediateTables);
}


#endif