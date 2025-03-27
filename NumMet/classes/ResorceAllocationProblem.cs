public class CreditOptimizer {
    public int[][] CalculateIntermediateTables(int[] funds, int[][] enterprisesData) {
        int numEnterprises = enterprisesData.Length;
        int numFundsOptions = funds.Length;

        // Инициализация таблиц F1-F4
        int[][] F = new int[numEnterprises][];
        
        // F1 - просто копируем данные первого предприятия
        F[0] = enterprisesData[0].ToArray();

        // Вычисляем последующие таблицы
        for (int k = 1; k < numEnterprises; k++) {
            F[k] = new int[numFundsOptions];
            
            for (int i = 0; i < numFundsOptions; i++) {
                int maxValue = 0;
                
                for (int j = 0; j <= i; j++) {
                    int currentValue = enterprisesData[k][j] + F[k-1][i-j];
                    if (currentValue > maxValue) {
                        maxValue = currentValue;
                    }
                }
                
                F[k][i] = maxValue;
            }
        }

        return F;
    }

    public void PrintIntermediateTables(int[] funds, int[][] F) {
        Console.WriteLine("Итоговая таблица:");
        
        for (int k = 0; k < F.Length; k++) {
            Console.Write($"F{k+1}: {{");
            for (int i = 0; i < F[k].Length; i++) {
                Console.Write(F[k][i]);
                if (i < F[k].Length - 1) Console.Write(",");
            }
            Console.WriteLine("}");
        }
    }
}