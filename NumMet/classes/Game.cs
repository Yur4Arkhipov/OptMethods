
public class Game {
    private double[,] payoffMatrix;

    public Game(double[,] matrix) {
        payoffMatrix = matrix;
    }

    public PureStrategyResult FindPureStrategy() {
        int rows = payoffMatrix.GetLength(0);
        int cols = payoffMatrix.GetLength(1);

        // Находим минимумы по строкам
        double[] rowMin = new double[rows];
        int[] rowMinIndices = new int[rows];
        for (int i = 0; i < rows; i++) {
            rowMin[i] = payoffMatrix[i, 0];
            rowMinIndices[i] = 0;
            for (int j = 1; j < cols; j++) {
                if (payoffMatrix[i, j] < rowMin[i]) {
                    rowMin[i] = payoffMatrix[i, j];
                    rowMinIndices[i] = j;
                }
            }
        }

        // Находим максимум из минимумов по строкам (maximin)
        double maximin = rowMin[0];
        int maximinRow = 0;
        for (int i = 1; i < rows; i++) {
            if (rowMin[i] > maximin) {
                maximin = rowMin[i];
                maximinRow = i;
            }
        }

        // Находим максимумы по столбцам
        double[] colMax = new double[cols];
        int[] colMaxIndices = new int[cols];
        for (int j = 0; j < cols; j++) {
            colMax[j] = payoffMatrix[0, j];
            colMaxIndices[j] = 0;
            for (int i = 1; i < rows; i++) {
                if (payoffMatrix[i, j] > colMax[j]) {
                    colMax[j] = payoffMatrix[i, j];
                    colMaxIndices[j] = i;
                }
            }
        }

        // Находим минимум из максимумов по столбцам (minimax)
        double minimax = colMax[0];
        int minimaxCol = 0;
        for (int j = 1; j < cols; j++) {
            if (colMax[j] < minimax) {
                minimax = colMax[j];
                minimaxCol = j;
            }
        }

        // Ищем седловые точки
        List<SaddlePoint> saddlePoints = new List<SaddlePoint>();
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                bool isRowMin = true;
                for (int k = 0; k < cols; k++) {
                    if (payoffMatrix[i, k] < payoffMatrix[i, j]) {
                        isRowMin = false;
                        break;
                    }
                }

                bool isColMax = true;
                for (int k = 0; k < rows; k++) {
                    if (payoffMatrix[k, j] > payoffMatrix[i, j]) {
                        isColMax = false;
                        break;
                    }
                }

                if (isRowMin && isColMax) {
                    saddlePoints.Add(new SaddlePoint(i, j, payoffMatrix[i, j]));
                }
            }
        }

        if (saddlePoints.Count > 0) {
            return new PureStrategyResult(
                true,
                saddlePoints,
                maximin,
                minimax
            );
        }

        return new PureStrategyResult(
            false,
            new List<SaddlePoint>(),
            maximin,
            minimax
        );
    }

    public MixedStrategyResult FindMixedStrategy() {
        int rows = payoffMatrix.GetLength(0);
        int cols = payoffMatrix.GetLength(1);
        
        // Для 2x2 игр используем аналитическое решение
        if (rows == 2 && cols == 2) {
            return SolveFor2x2MixedStrategy();
        }
        
        return SolveMixedStrategy();   
    }
    
    public MixedStrategyResult SolveMixedStrategy() {
        var r = FindPureStrategy();
        if (r.SaddlePoints.Count > 0) {
            throw new Exception("Есть седловые точки");
        }

        int rows = payoffMatrix.GetLength(0);
        int cols = payoffMatrix.GetLength(1);
        double min = double.MaxValue;
        MixedStrategyResult bestStrategy = null;
        int[] bestColumns = null;

        for (int i = 0; i < cols; i++) {
            for (int j = i + 1; j < cols; j++) {
                double[,] B = {
                    {payoffMatrix[0, i], payoffMatrix[0, j]},
                    {payoffMatrix[1, i], payoffMatrix[1, j]},
                };
                try {
                    var strategy = new Game(B).SolveFor2x2MixedStrategy();
                    if (strategy.GameValue < min) {
                        min = strategy.GameValue;
                        
                        // Создаем новый результат с информацией об используемых столбцах
                        double[] secondPlayerFullStrategy = new double[cols];
                        secondPlayerFullStrategy[i] = strategy.SecondPlayerStrategy[0];
                        secondPlayerFullStrategy[j] = strategy.SecondPlayerStrategy[1];
                        
                        bestStrategy = new MixedStrategyResult(
                            strategy.FirstPlayerStrategy,
                            secondPlayerFullStrategy,
                            strategy.GameValue,
                            new int[] { i, j }
                        );
                        
                        bestColumns = new int[] { i, j };
                    }
                } catch (DivideByZeroException) {
                    continue;
                }
            }
        }

        if (bestStrategy == null) {
            throw new Exception("Не удалось найти смешанную стратегию");
        }

        return bestStrategy;
    }

    private MixedStrategyResult SolveFor2x2MixedStrategy() {
        double a = payoffMatrix[0, 0];
        double b = payoffMatrix[0, 1];
        double c = payoffMatrix[1, 0];
        double d = payoffMatrix[1, 1];
        
        double denominator = a - b - c + d;
        
        // Оптимальные стратегии для первого игрока (p1, p2)
        double p1 = (d - c) / denominator;
        double p2 = 1 - p1;
        
        // Оптимальные стратегии для второго игрока (q1, q2)
        double q1 = (d - b) / denominator;
        double q2 = 1 - q1;
        
        // Цена игры
        double gameValue = (a * d - b * c) / denominator;
        
        // Используемые столбцы для 2x2 матрицы - это просто столбцы 0 и 1
        int[] usedColumns = new int[] { 0, 1 };
        
        return new MixedStrategyResult(
            new double[] { p1, p2 },
            new double[] { q1, q2 },
            gameValue,
            usedColumns
        );
    }

    public void PrintStrategy() {
        var pureResult = FindPureStrategy();
        if (pureResult.HasPureStrategy) {
            Console.WriteLine("Найдена чистая стратегия!");
            Console.WriteLine($"Цена игры: {pureResult.GameValue}");
            Console.WriteLine("Седловые точки:");
            foreach (var point in pureResult.SaddlePoints){
                Console.WriteLine($"  Строка: {point.Row + 1}, Столбец: {point.Column + 1}, Значение: {point.Value}");
            }
        } else {
            Console.WriteLine("Чистая стратегия не существует");
            Console.WriteLine($"Maximin = {pureResult.Maximin}, Minimax = {pureResult.Minimax}");
            
            // Добавляем вывод смешанной стратегии
            var mixedResult = FindMixedStrategy();
            Console.WriteLine("\nНайдена смешанная стратегия:");
            Console.WriteLine($"Цена игры: {mixedResult.GameValue}");
            
            Console.WriteLine("Оптимальная стратегия первого игрока (строки):");
            for (int i = 0; i < mixedResult.FirstPlayerStrategy.Length; i++) {
                Console.WriteLine($"  Стратегия {i + 1}: {mixedResult.FirstPlayerStrategy[i]:F4}");
            }
            
            Console.WriteLine("Оптимальная стратегия второго игрока (столбцы):");
            for (int j = 0; j < mixedResult.SecondPlayerStrategy.Length; j++) {
                if (mixedResult.SecondPlayerStrategy[j] > 0) {
                    Console.WriteLine($"  Стратегия {j + 1}: {mixedResult.SecondPlayerStrategy[j]:F4}");
                }
            }
            
            // Выводим информацию об используемых столбцах
            if (mixedResult.UsedColumns != null) {
                Console.WriteLine("\nИспользуемые столбцы в смешанной стратегии:");
                foreach (int col in mixedResult.UsedColumns) {
                    Console.WriteLine($"  Столбец {col + 1}");
                }
            }
        }
    }
}

public class PureStrategyResult {
    public bool HasPureStrategy { get; }
    public List<SaddlePoint> SaddlePoints { get; }
    public double Maximin { get; }
    public double Minimax { get; }
    public double GameValue => HasPureStrategy ? SaddlePoints[0].Value : 0;

    public PureStrategyResult(bool hasPureStrategy, List<SaddlePoint> saddlePoints, double maximin, double minimax) {
        HasPureStrategy = hasPureStrategy;
        SaddlePoints = saddlePoints;
        Maximin = maximin;
        Minimax = minimax;
    }
}

public class SaddlePoint {
    public int Row { get; }
    public int Column { get; }
    public double Value { get; }

    public SaddlePoint(int row, int column, double value) {
        Row = row;
        Column = column;
        Value = value;
    }
}


public class MixedStrategyResult {
    public double[] FirstPlayerStrategy { get; }
    public double[] SecondPlayerStrategy { get; }
    public double GameValue { get; }
    public int[] UsedColumns { get; }
    
    public MixedStrategyResult(double[] firstPlayerStrategy, double[] secondPlayerStrategy, double gameValue) {
        FirstPlayerStrategy = firstPlayerStrategy;
        SecondPlayerStrategy = secondPlayerStrategy;
        GameValue = gameValue;
        UsedColumns = null; // Если столбцы не указаны
    }
    
    public MixedStrategyResult(double[] firstPlayerStrategy, double[] secondPlayerStrategy, double gameValue, int[] usedColumns) {
        FirstPlayerStrategy = firstPlayerStrategy;
        SecondPlayerStrategy = secondPlayerStrategy;
        GameValue = gameValue;
        UsedColumns = usedColumns;
    }
}