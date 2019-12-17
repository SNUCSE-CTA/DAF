public class j2k_graph {
    int N;
    int[] degree;
    int[][] adjList;
    int[] label;
    int id;
    public j2k_graph(int n){
        N=n;
        degree = new int[n];
        adjList = new int[n][n];
        label = new int[n];
    }
}
