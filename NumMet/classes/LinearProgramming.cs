class LinearProgramming {
    public int[,] SolveTransportationProblem(int[] supply, int[] demand, int[,] cost) {
        int m = supply.Length;
        int n = demand.Length;
        int[,] allocation = new int[m, n];
        
        int i = 0, j = 0;
        while (i < m && j < n) {
            int allocated = Math.Min(supply[i], demand[j]);
            allocation[i, j] = allocated;
            supply[i] -= allocated;
            demand[j] -= allocated;

            if (supply[i] == 0) i++;
            if (demand[j] == 0) j++;
        }

        return allocation;
    }

    public void PrintSolution(int[,] allocation, int[,] cost) {
        int totalCost = 0;
        Console.WriteLine("Transportation Matrix:");
        for (int i = 0; i < allocation.GetLength(0); i++) {
            for (int j = 0; j < allocation.GetLength(1); j++) {
                Console.Write("{0,4}", allocation[i, j]);
                totalCost += allocation[i, j] * cost[i, j];
            }
            Console.WriteLine();
        }
        Console.WriteLine("Total Cost: " + totalCost);
    }
}