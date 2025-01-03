﻿// #define NonLinear
// #define LinearAlgEquations
// #define Spline
// #define MinRect
// #define Integral
// #define Differential
// #define OptimumOneArg 
// #define OptimumMoreOneArg
// #define Simplex



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
    // при слишком малом сходится медленно и точчность хуже.
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
{
    // метод деформированного многогранника
    OptimumResult result = OptimumVector.NelderMeadeMethod(new Vector(new double[] { 0.0, 0.0 }), 0.1, 0.000001, Rosenbrock);
    Console.WriteLine("NedlerMeadOptimization func=Rosenbrock {0}", result);
}


#elif Simplex



#endif


