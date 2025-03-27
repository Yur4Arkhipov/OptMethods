using System.Transactions;

public class Graphs {
   public enum COLORS_VERTEX { WHITE, GRAY, BLACK}

    public class Edge {
        public Vertex BeginPoint; // Начальная вершина
        public Vertex EndPoint;  // Конечная вершина

        public double distance; // Длина ребра
        
        // Конструктор
        public Edge(Vertex begin, Vertex end, double distance)
        {
            this.BeginPoint = begin;
            this.EndPoint = end;
            this.distance = distance;
        }
        public override string ToString() {
            string sout = "";
            sout = "{"+BeginPoint.label + "  " + EndPoint.label + " Distance: " + distance.ToString()+"}";
            return sout;
        }
    }

   public class Vertex {
        private static int IDV = 0; 
        private int ID;
        public string label; // Метка (имя вершины)
        private List<Edge> edges; // Список ребер, связанных с вершиной
        public double sumDistance = 0; // Сумма растояний
        public COLORS_VERTEX color; // Цвет вершины
        public Vertex previousVertex; // Ссылка на предшественника
        public bool visited;

        //конструктор
        public Vertex(string label) {
            this.label = label;
            IDV++;
            edges = new List<Edge>();
            sumDistance = Double.MaxValue;
            color = COLORS_VERTEX.WHITE;
            previousVertex = null;
            ID = IDV;
            this.visited = false;
        }

        public int GetID() { return ID; }
        // Получение списка ребер
        public List<Edge> GetEdges() { return edges; }
        public override string ToString() {
            string sout = "";
            sout += label;
            sout = sout + "(ID="+ID.ToString()+")";
            return sout ;
        }
        // Просмотр ребер, связанных с вершиной
        public void ViewEdges() {
            Console.Write("Edges for {0}", this);
            foreach(Edge curedge in edges)
                Console.Write("  {0}", curedge);
            Console.WriteLine();
        }
        // Добавление ребра
        public bool AddEdge(Edge edge) {
            if (edge.BeginPoint != this) return false;
            for (int i = 0; i < edges.Count; i++)
            {
                Edge CurEdge = edges[i];
                if (edge.EndPoint.Equals(CurEdge.EndPoint)) return false;
            }
            edges.Add(edge);
            return true;
        }        
    } 

   public class Graph {
        public List<Vertex> allVertexs; // Список всех вершин
        public List<Edge> allEdges; // Список всех ребер
        private bool isDirected; // флаг для определения типа графа
        //конструктор
        public Graph() {
            allVertexs = new List<Vertex>();
            allEdges = new List<Edge>();
        }
        public Graph(bool isDirected = true) {
            this.isDirected = isDirected;
            allVertexs = new List<Vertex>();
            allEdges = new List<Edge>();
        }

        //добавление ребра
        public bool AddEdge(Vertex v1, Vertex v2, double dist) {
            if (!allVertexs.Contains(v1)) return false;
            if (!allVertexs.Contains(v2)) return false;
            foreach (Edge item in v1.GetEdges()) {
                if (item.EndPoint.GetID() == v2.GetID()) return false;
            }
            
            Edge ev1v2 = new Edge(v1, v2, dist);
            v1.GetEdges().Add(ev1v2); allEdges.Add(ev1v2);
            return true;
        }

        public bool AddDualEdge(Vertex v1, Vertex v2, double dist) {
            if(!allVertexs.Contains(v1)) return false;
            if(!allVertexs.Contains(v2)) return false;
            foreach (Edge item in v1.GetEdges()) {
                if (item.EndPoint.GetID() == v2.GetID()) return false;
            }
            foreach (Edge item in v2.GetEdges()) {
                if (item.EndPoint.GetID() == v1.GetID()) return false;
            }
            Edge ev1v2 = new Edge(v1, v2, dist);
            Edge ev2v1 = new Edge(v2, v1, dist);
            v1.GetEdges().Add(ev1v2); allEdges.Add(ev1v2);
            v2.GetEdges().Add(ev2v1); allEdges.Add(ev2v1);
            
            return true;
        }

        public void PrintAdjacencyList(Graph graph) {
            Console.WriteLine("Adjacency list of Graph: ");
            foreach (var vertex in graph.allVertexs) {
                Console.Write($"Vertex {vertex.label}: ");
                foreach (var edge in vertex.GetEdges()) {
                    Console.Write($"{edge.EndPoint.label}({edge.distance}) ");
                }
                Console.WriteLine();
            }
        }


        // Поиск в ширину
        public void BFS(Vertex source) {
            Queue<Vertex> queue = new Queue<Vertex>(); // Очередь вершин
            List<Vertex> traversedEdges = new List<Vertex>(); // Список пройденных рёбер
            // Инициализация
            foreach (Vertex vertex in allVertexs) {
                vertex.sumDistance = double.MaxValue;
                vertex.previousVertex = null;
                vertex.color = COLORS_VERTEX.WHITE;
            }
            source.previousVertex = null;
            source.color = COLORS_VERTEX.GRAY;
            source.sumDistance = 0;
            queue.Enqueue(source);
            Vertex u, v; // u - следующая вершина
            Edge tr;
            List<Edge> edges_u;
            // Основной цикл
            while (queue.Count > 0) {
                u = queue.Dequeue();
                edges_u = u.GetEdges(); // список рёбер, исходящих из u
                foreach (Edge edge in edges_u) {
                    v = edge.EndPoint; // Получаем конечную вершину рассматриваемого ребра
                    if (u.sumDistance + edge.distance < v.sumDistance) {
                        v.sumDistance = u.sumDistance + edge.distance;
                    }
                    if (v.color == COLORS_VERTEX.WHITE) { // Если вершина v белая, то она не была посещена
                        v.color = COLORS_VERTEX.GRAY; // Помечаем вершину как серую
                        v.sumDistance = u.sumDistance + edge.distance; // Обновляем расстояние до вершины v
                        v.previousVertex = u; // Запоминаем, что вершина v была достигнута из вершины u
                        queue.Enqueue(v); // Добавляем вершину v в очередь для дальнейшего рассмотрения
                        traversedEdges.Add(v); // Добавляем ребро в список пройденных рёбер
                    }
                }
                u.color = COLORS_VERTEX.BLACK;
            }
            /* Console.WriteLine("BFS list edgees:"); // Печатаем список пройденных рёбер
            Console.Write(source + "->");
            foreach(var item in traversedEdges) {
                Console.Write(item + "->");
            } */
        }


        // поиск в глубину
        private static int time = 0;
        public void DFS(Vertex source) {
            List<Vertex> traversedEdges = new List<Vertex>(); // Список пройденных вершин

            foreach (Vertex vertex in allVertexs) {
                vertex.color = COLORS_VERTEX.WHITE;
                vertex.previousVertex = null;
            }
            time = 0;
            DFSVisit(source, traversedEdges);
            /* Console.WriteLine("DFS list edgees:"); // Печатаем список пройденных рёбер
            Console.Write(source + "->");
            foreach(var item in traversedEdges) {
                Console.Write(item + "->");
            } */
        }

        // Рекурсивная функция для обхода в глубину
        private void DFSVisit(Vertex u, List<Vertex> traversedEdges) {
            u.color = COLORS_VERTEX.GRAY;
            time++;
            u.sumDistance = time;

            foreach (Edge edge in u.GetEdges()) {
                Vertex v = edge.EndPoint;
                if (v.color == COLORS_VERTEX.WHITE) {
                    v.previousVertex = u;
                    traversedEdges.Add(v); //добавление вершины в список пройденных вершин
                    DFSVisit(v, traversedEdges); // Рекурсивно вызываем DFS_VISIT для вершины v
                }
            }

            u.color = COLORS_VERTEX.BLACK; // Помечаем вершину как обработанную
            time++;
            u.visited = true;
        }

        public void PrintPath(Vertex s, Vertex v) {
            List<Vertex> listVertex = new List<Vertex>(); // Список пройденных вершин
            if (v == s) {
                Console.Write(s.label + " ");
            } else if (v.previousVertex == null) {
                Console.WriteLine("Path from " + s.label + " to " + v.label + " is unvailable");
            } else {
                PrintPath(s, v.previousVertex);
                Console.Write(v.label + " ");
            }
        }

        public List<Vertex> GetPath(Vertex s, Vertex v) {
            List<Vertex> path = new List<Vertex>();

            // Если вершина v не достижима из вершины s, вернуть пустой путь
            if (v.previousVertex == null) {
                Console.WriteLine("Путь из " + s.label + " в " + v.label + " отсутствует");
                return path;
            }

            // Построение пути от вершины v до вершины s
            Vertex currentVertex = v;
            while (currentVertex != s) {
                path.Add(currentVertex);
                currentVertex = currentVertex.previousVertex;
            }
            path.Add(s); // Добавляем начальную вершину s в конец пути
            path.Reverse(); // Переворачиваем путь, чтобы он шел от вершины s до вершины v

            return path;
        }

        /* Находим эксцентриситеты (длину максимально удаленной вершины) для каждой вершины
        центр графа - вершины с минимальным эксцентриситетом */
        public List<Vertex> GetCenters() {
            List<Vertex> centers = new List<Vertex>();
            List<Vertex> vertices = allVertexs;
            int minEccentricity = int.MaxValue;
            //находим расстояние до максимально удаленной вершлины - эксцентриситет
            // Найдем минимальный эксцентриситет
            //перебираем все вершины
            foreach (Vertex vertex in vertices) {
                //находим все расстояния от текущей вершины до всех других вершин 
                BFS(vertex);
                //обнуляем максимальное значение дистанции для текущей вершины
                int maxDistance = int.MinValue;
                // Найдем максимальное расстояние от текущей вершины до других вершин
                foreach (Vertex otherVertex in vertices) {
                    if (otherVertex.sumDistance > maxDistance) {
                        maxDistance = (int)otherVertex.sumDistance;
                    }
                }
                // Обновим минимальный эксцентриситет
                if (maxDistance < minEccentricity) {
                    minEccentricity = maxDistance;
                }
            }

            // Найдем все вершины с минимальным эксцентриситетом
            foreach (Vertex vertex in vertices) {
                BFS(vertex);
                int maxDistance = int.MinValue;
                foreach (Vertex otherVertex in vertices) {
                    if (otherVertex.sumDistance > maxDistance) {
                        maxDistance = (int)otherVertex.sumDistance;
                    }
                }
                // Если текущая вершина имеет минимальный эксцентриситет, добавляем ее в список центров
                if (maxDistance == minEccentricity) {
                    centers.Add(vertex);
                }
            }

            return centers;
        }

        // Находит путь с помощью BFS в остаточной сети
        private bool BFSPath(Vertex source, Vertex sink, Dictionary<Vertex, Vertex> parent) {
            // Инициализация всех вершин как непосещенные
            foreach (Vertex vertex in allVertexs) {
                vertex.color = COLORS_VERTEX.WHITE;
            }

            Queue<Vertex> queue = new Queue<Vertex>();
            queue.Enqueue(source);
            source.color = COLORS_VERTEX.GRAY;
            parent[source] = null;

            // Поиск в ширину
            while (queue.Count > 0) {
                Vertex u = queue.Dequeue();
                foreach (Edge edge in u.GetEdges()) {
                    Vertex v = edge.EndPoint;
                    // Проверяем, что вершина не посещена и есть доступная пропускная способность
                    if (v.color == COLORS_VERTEX.WHITE && edge.distance > 0) {
                        queue.Enqueue(v);
                        v.color = COLORS_VERTEX.GRAY;
                        parent[v] = u;
                    }
                }
                u.color = COLORS_VERTEX.BLACK;
            }

            // Возвращаем true, если достигли стока
            return sink.color != COLORS_VERTEX.WHITE;
        }

        private Edge GetEdge(Vertex u, Vertex v) {
            return u.GetEdges().Find(e => e.EndPoint == v);
        }

        // Обновляет пропускные способности в остаточной сети
        private void UpdateFlow(Vertex u, Vertex v, double flow) {
            Edge forward = GetEdge(u, v);
            Edge backward = GetEdge(v, u);

            if (forward != null) {
                forward.distance -= flow;
                if (!isDirected && backward != null) {
                    // Для неориентированного графа синхронизируем пропускные способности
                    backward.distance = forward.distance;
                }
            }
            // Для ориентированного графа увеличиваем обратный поток
            if (isDirected && backward != null) {
                backward.distance += flow;
            }
        }

        // Основной алгоритм Форда-Фалкерсона
        public double FordFulkerson(Vertex source, Vertex sink) {
            if (source == sink) return 0;

            double maxFlow = 0;
            Dictionary<Vertex, Vertex> parent = new Dictionary<Vertex, Vertex>();
            
            // Сохраняем начальные пропускные способности
            Dictionary<(Vertex, Vertex), double> initialCapacities = new Dictionary<(Vertex, Vertex), double>();
            Dictionary<(Vertex, Vertex), double> currentFlows = new Dictionary<(Vertex, Vertex), double>();
            foreach (Edge edge in allEdges) {
                initialCapacities.Add((edge.BeginPoint, edge.EndPoint), edge.distance);
                currentFlows.Add((edge.BeginPoint, edge.EndPoint), 0);
            }

            // Создаем остаточную сеть
            Graph residualGraph = new Graph(isDirected);
            foreach (Vertex v in allVertexs) {
                residualGraph.allVertexs.Add(v);
            }

            // Копируем ребра с учетом типа графа
            foreach (Edge edge in allEdges) {
                residualGraph.AddEdge(edge.BeginPoint, edge.EndPoint, edge.distance);
                if (isDirected) {
                    if (GetEdge(edge.EndPoint, edge.BeginPoint) == null) {
                        residualGraph.AddEdge(edge.EndPoint, edge.BeginPoint, 0);
                    }
                }
            }

            Console.WriteLine("\nНачальные пропускные способности ребер:");
            foreach (var edge in allEdges) {
                Console.WriteLine($"Ребро {edge.BeginPoint.label} -> {edge.EndPoint.label}: {edge.distance}");
            }

            while (BFSPath(source, sink, parent)) {
                double pathFlow = double.MaxValue;
                for (Vertex v = sink; v != source; v = parent[v]) {
                    Vertex u = parent[v];
                    Edge edge = GetEdge(u, v);
                    pathFlow = Math.Min(pathFlow, edge.distance);
                }

                // Вывод найденного пути и его пропускной способности
                Console.WriteLine($"\nНайден путь с пропускной способностью {pathFlow}:");
                for (Vertex v = sink; v != source; v = parent[v]) {
                    Console.Write($"{v.label} <- ");
                }
                Console.WriteLine(source.label);

                // Обновление потока
                for (Vertex v = sink; v != source; v = parent[v]) {
                    Vertex u = parent[v];
                    UpdateFlow(u, v, pathFlow);
                    currentFlows[(u, v)] += pathFlow;
                }

                maxFlow += pathFlow;
            }

            Console.WriteLine("\nОстаточные пропускные способности ребер:");
            foreach (var edge in allEdges) {
                double initialCapacity = initialCapacities[(edge.BeginPoint, edge.EndPoint)];
                double currentFlow = currentFlows[(edge.BeginPoint, edge.EndPoint)];
                double remainingCapacity = edge.distance;
                // double usedCapacity = initialCapacity - remainingCapacity;
                Console.WriteLine($"Ребро {edge.BeginPoint.label} -> {edge.EndPoint.label}: " +
                    $"использовано {currentFlow} из {initialCapacity}");
            }

            Console.WriteLine($"\nМаксимальный поток: {maxFlow}");
            return maxFlow;
        }    
    }
}