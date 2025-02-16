// #include "simplex.h"
#include <iostream>
#include <fstream>
#include <ostream>
#include <string>
#include <vector>
#include <iomanip>
#include <limits>

// Большое число
const double BIG_M = 1000;

void open_files(
  std::ifstream &input_file, 
  std::ofstream &output_file, 
  const std::string &input_fileName, 
  const std::string &output_fileName
) {
  // Открываем входной файл
  input_file.open(input_fileName);
  // Открываем выходной файл
  output_file.open(output_fileName);

  // Проверяем, удалось ли открыть входной файл
  if (!input_file.is_open()) {
      // Если нет, выводим сообщение об ошибке и завершаем программу
      std::cerr << "Error opening input file: " << input_fileName << std::endl;
      exit(1);
  }
  // Проверяем, удалось ли открыть выходной файл
  if (!output_file.is_open()) {
      // Если нет, выводим сообщение об ошибке и завершаем программу
      std::cerr << "Error opening output file: " << output_fileName << std::endl;
      exit(1);
  }
}

void read_input(
   std::ifstream &input_file, int &num_equations, int &num_variables,
   std::vector<std::vector<double>> &restriction_coefficients, 
   std::vector<double> &b, 
   std::vector<int> &signs,
   std::vector<double> &function_coefficients
) {
   // Читаем количество уравнений и переменных из файла
   input_file >> num_equations >> num_variables;

   // Изменяем размеры векторов в соответствии с количеством уравнений и переменных
   restriction_coefficients.resize(num_equations, std::vector<double>(num_variables));
   b.resize(num_equations);
   signs.resize(num_equations);
   function_coefficients.resize(num_variables);

   // Читаем коэффициенты ограничений, значения b и знаки из файла
   for (int i = 0; i < num_equations; ++i) {
       for (int j = 0; j < num_variables; ++j) {
           input_file >> restriction_coefficients[i][j];
       }
       input_file >> b[i];
       input_file >> signs[i];
   }

   // Читаем коэффициенты функции из файла
   for (int j = 0; j < num_variables; ++j) {
       input_file >> function_coefficients[j];
   }
}

void initialize_table(
    const std::vector<std::vector<double>> &restriction_coefficients, 
    const std::vector<double> &b, 
    const std::vector<double> &function_coefficients, 
    std::vector<std::vector<double>> &table,
    std::vector<int>& cure_bases_idx, 
    std::vector<double>& costs
) {
    int num_equations = restriction_coefficients.size();
    int num_variables = restriction_coefficients[0].size();
    int num_cols = table[0].size();

    // Заполнение таблицы (основные коэффициенты)
    for (int i = 0; i < num_equations; ++i) {
        for (int j = 0; j < num_variables; ++j) {
            table[i][j] = restriction_coefficients[i][j];
        }
        table[i][num_cols - 1] = b[i]; // Значение свободного члена
    }

    cure_bases_idx.clear(); // Очистка индексов базисов
    for (int j = 0; j < num_cols - 1; ++j) { // Пробегаем по всем столбцам, кроме свободных членов
        int one_count = 0;
        for (int i = 0; i < num_equations; ++i) {
            if (table[i][j] == 1) {
                one_count++;
            } else if (table[i][j] != 0) {
                one_count = -1; // Если встретился не 0 или 1, столбец не базисный
                break;
            }
        }

        if (one_count == 1) {
            cure_bases_idx.push_back(j); // Добавляем индекс базисного столбца
        }
    }

    // Вычисление значений для последней строки (строка целевой функции)
    for (int j = 0; j < num_cols - 1; ++j) { // Пробегаем по всем столбцам, кроме столбца свободных членов
        double sum = 0.0;
        // Суммируем для всех базисных переменных
        for (int k = 0; k < cure_bases_idx.size(); ++k) {
            int base_idx = cure_bases_idx[k];
            sum += costs[base_idx] * table[k][j]; // c2 * a0j, c4 * a1j и т.д.
        }

        // Пересчитываем коэффициент в строке цели для переменной j
        table[num_equations][j] = sum - costs[j]; // k_i = сумма - Ci
    }

    // Для свободного члена в строке целевой функции
    double sum_b = 0.0;
    for (int k = 0; k < cure_bases_idx.size(); ++k) {
        int base_idx = cure_bases_idx[k];
        sum_b += costs[base_idx] * table[k][num_cols - 1]; // c2 * b0, c4 * b1 и т.д.
    }

    // Записываем пересчитанное значение свободного члена
    table[num_equations][num_cols - 1] = sum_b; 
}

void canonicalize(
    std::vector<std::vector<double>> &restriction_coefficients, 
    std::vector<double> &b, 
    std::vector<int> &signs,
    std::vector<double> &function_coefficients, 
    std::vector<std::vector<double>> &canonical_rest_coefs, 
    std::vector<double> &canonical_fun_coefs, 
    std::vector<double>& costs,
    std::ofstream &output_file
) {
    int num_equations = restriction_coefficients.size();
    int num_variables = function_coefficients.size();

    // Рассчитаем количество новых переменных
    int slack_count = 0, artifical_count = 0;

    for (int sign : signs) {
        if (sign == -1) { // <=
            ++slack_count;
        } else if (sign == 1) { // >=
            ++slack_count;
            ++artifical_count;
        } else if (sign == 0) { // =
            ++artifical_count;
        }
    }

    int total_variables = num_variables + slack_count + artifical_count;

    // Расширяем матрицу коэффициентов
    canonical_rest_coefs.resize(num_equations, std::vector<double>(total_variables, 0.0f));
    canonical_fun_coefs.resize(total_variables, 0.0f);
    costs.resize(total_variables);

    int slack_idx = num_variables;                   // Индекс slack-переменных (базисных не искусственных)
    int artifical_idx = num_variables + slack_count; // Индекс искусственных переменных

    for (int i = 0; i < num_equations; ++i) {
        for (int j = 0; j < num_variables; ++j) {
            canonical_rest_coefs[i][j] = restriction_coefficients[i][j];
        }

        if (signs[i] == -1) { // <=
            canonical_rest_coefs[i][slack_idx++] = 1.0f;
        } else if (signs[i] == 1) { // >=
            canonical_rest_coefs[i][slack_idx++] = -1.0f; 
            canonical_rest_coefs[i][artifical_idx++] = 1.0f;
        } else if (signs[i] == 0) { // =
            canonical_rest_coefs[i][artifical_idx++] = 1.0f;
        }
    }
    for (int j = 0; j < slack_count; ++j) {
        slack_idx--;
    }
    for (int j = 0; j < artifical_count; ++j) {
        artifical_idx--;
    }

    // Целевая функция: переносим только основные переменные
    for (int j = 0; j < num_variables; ++j) {
        canonical_fun_coefs[j] = function_coefficients[j];
    }
    // std::cout << artifical_count << std::endl;
    // std::cout << slack_count << std::endl;
    // std::cout << total_variables << std::endl;
    // std::cout << function_coefficients.size() << std::endl;
    // std::cout << artifical_idx << std::endl;
    // std::cout << slack_idx << std::endl;

    for (int i = 0; i < total_variables; i++) {
        if (i < function_coefficients.size()) {
            costs[i] = function_coefficients[i];
        } else if (i == artifical_idx) {
            costs[i] = -BIG_M;
            artifical_idx++;
        } else {
            costs[i] = 0.0f;
        }
    }
}

void print_table(
    std::ofstream &output_file,
    const std::vector<std::vector<double>> &table,
    int &count_iteration,
    const std::vector<int>& cur_base
) {
    int num_rows = table.size();
    int num_cols = table[0].size();
    count_iteration++;

    // Ширина для каждого столбца
    const int col_width = 12;

    // Заголовок для итерации
    output_file << "\n==========================================\n";
    output_file << "Iteration " << count_iteration << ":\n";
    output_file << "==========================================\n";

    // Заголовки столбцов
    output_file << std::setw(col_width) << "Base";
    for (int j = 0; j < num_cols - 1; ++j) {
        output_file << std::setw(col_width) << ("X" + std::to_string(j + 1));
    }
    output_file << std::setw(col_width) << "b" << "\n";

    // Разделитель
    output_file << std::string(col_width * (num_cols + 1), '-') << "\n";

    // Вывод таблицы
    for (int i = 0; i < num_rows; ++i) {
        // Базисная переменная
        if (i < num_rows - 1) {
            output_file << std::setw(col_width) << ("X" + std::to_string(cur_base[i] + 1));
        } else {
            output_file << std::setw(col_width) << " ";
        }

        // Элементы строки
        for (int j = 0; j < num_cols; ++j) {
            output_file << std::setw(col_width) << std::fixed << std::setprecision(2) << table[i][j];
        }
        output_file << "\n";
    }

    // Разделитель после таблицы
    output_file << std::string(col_width * (num_cols + 1), '-') << "\n";
}

bool simplex_iteration(
    std::vector<std::vector<double>> &table,
    std::ofstream &output_file,
    int num_equations,
    int num_variables,
    std::vector<int>& cur_bases
) {
    int num_rows = table.size();
    int num_cols = table[0].size();

    // Найти опорный столбец (самое отрицательное значение в строке цели)
    int pivot_col = -1;
    double min_value = 0.0f;
    for (int j = 0; j < num_cols - 1; ++j) {
        if (table[num_rows - 1][j] < min_value) {
            min_value = table[num_rows - 1][j];
            pivot_col = j;
        }
    }

    if (pivot_col == -1) return false; // Оптимальное решение найдено

    // Найти опорную строку (минимальное положительное отношение правой части к элементу в опорном столбце)
    int pivot_row = -1;
    double min_ratio = std::numeric_limits<double>::max();
    for (int i = 0; i < num_equations; ++i) {
        if (table[i][pivot_col] > 0) {
            double ratio = table[i][num_cols - 1] / table[i][pivot_col];
            if (ratio < min_ratio) {
                min_ratio = ratio;
                pivot_row = i;
            }
        }
    }

    if (pivot_row == -1) {
        output_file << "\nUnbounded solution\n";
        return false;
    }
    cur_bases[pivot_row] = pivot_col;

    // Пересчет элементов матрицы
    double pivot_element = table[pivot_row][pivot_col];
    for (int j = 0; j < num_cols; ++j) {
        table[pivot_row][j] /= pivot_element;
    }

    for (int i = 0; i < num_rows; ++i) {
        if (i != pivot_row) {
            double factor = table[i][pivot_col];
            for (int j = 0; j < num_cols; ++j) {
                table[i][j] -= factor * table[pivot_row][j];
            }
        }
    }

    return true;
}

void write_solution(
    std::ofstream &output_file, 
    double result,
    const std::vector<std::vector<double>> &table,
    const std::vector<int>& cur_bases
) {
    int num_equations = table.size() - 1; // Количество ограничений
    int num_variables = table[0].size() - 1; // Количество переменных (без b)

    // Запись итогового значения целевой функции
    output_file << "\n==========================================\n";
    output_file << "Final Solution:\n";
    output_file << "==========================================\n";
    output_file << "Objective Function Value: " << std::fixed << std::setprecision(2) << result << "\n";

    // Вектор для хранения значений переменных
    std::vector<double> solution(num_variables, 0.0);

    // Заполнение значений для базисных переменных
    for (int i = 0; i < num_equations; ++i) {
        if (cur_bases[i] < num_variables) { // Проверка, что индекс базиса в пределах переменных
            solution[cur_bases[i]] = table[i][num_variables]; // Значение из столбца b
        }
    }

    // Запись значений переменных
    output_file << "Variable Values:\n";
    for (int i = 0; i < num_variables; ++i) {
        output_file << "X" << (i + 1) << " = " << std::fixed << std::setprecision(2) << solution[i] << "\n";
    }

    // Разделитель
    output_file << "==========================================\n";
}


double calculate_fun(
    const std::vector<int>& cure_bases_idx,
    const std::vector<double>& costs,
    const std::vector<std::vector<double>>& table
) {
    double result = 0.0;
    int num_rows = table.size();
    int num_cols = table[0].size();

    // Проходим по всем базисным индексам
    for (size_t i = 0; i < cure_bases_idx.size(); ++i) {
        int base_idx = cure_bases_idx[i];
        if (base_idx < num_cols - 1) { // Убедимся, что индекс в пределах переменных
            result += costs[base_idx] * table[i][num_cols - 1]; // cost * b
        }
    }

    return result;
}

int main() {
    int num_equations, num_variables;
    std::ifstream input_file;
    std::ofstream output_file;

    open_files(input_file, output_file, "input.txt", "output.txt");
    std::vector<std::vector<double>> restriction_coefficients;
    std::vector<double> b;
    std::vector<int> signs;
    std::vector<double> function_coefficients;
    read_input(input_file, num_equations, num_variables, restriction_coefficients, b, signs, function_coefficients);
    input_file.close();

    // Приведение к канонической форме
    std::vector<std::vector<double>> canonical_rest_coefs;
    std::vector<double> canonical_fun_coefs;
    std::vector<double> costs;
    canonicalize(restriction_coefficients, b, signs, function_coefficients, canonical_rest_coefs, canonical_fun_coefs, costs, output_file);
    
    // Инициализация симплекс таблицы
    int num_rows = num_equations + 1; // +1 для f
    int num_cols = canonical_fun_coefs.size() + 1; // +1 для b
    std::vector<std::vector<double>> table(num_rows, std::vector<double>(num_cols, 0.0f));
    std::vector<int> cure_bases_idx;
    initialize_table(canonical_rest_coefs, b, canonical_fun_coefs, table, cure_bases_idx, costs);

    int count_iteration = 0;
    print_table(output_file, table, count_iteration, cure_bases_idx);

    // пересчет элементов матрицы
    while (simplex_iteration(table, output_file, num_equations, num_variables, cure_bases_idx)) {
        print_table(output_file, table, count_iteration, cure_bases_idx);
    }
    
    double result = calculate_fun(cure_bases_idx, costs, table);
    write_solution(output_file, result, table, cure_bases_idx);
    output_file.close();

    std::cout << "Result: " << result << std::endl;
    std::vector<double> solution(num_variables, 0.0);
    for (int i = 0; i < num_equations; ++i) {
        if (cure_bases_idx[i] < num_variables) { 
            solution[cure_bases_idx[i]] = table[i][num_variables];
        }
    }
    std::cout << "Variable Values:\n";
    for (int i = 0; i < num_variables; ++i) {
        std::cout << "X" << (i + 1) << " = " << std::fixed << std::setprecision(2) << solution[i] << "\n";
    }
    std::cout << "b: \n";
    for (int i = 0; i < cure_bases_idx.size(); ++i) {
        std::cout << table[i][num_cols - 1] << " ";
    }
    std::cout << std::endl;

    return 0;
}