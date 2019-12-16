import java.util.*;

class main{
    public static void main(String[] args){
        int[][] ad = new int[11][10];
        int[] degree = new int[11];
        String e = "01 12 23 24 27 36 45 49 56 59 69 89";
        Scanner r = new Scanner(e);
        while (r.hasNextInt()){
            int a = r.nextInt();
            int u = a/10;
            int v = a%10;
            ad[u][degree[u]++] = v;
            ad[v][degree[v]++] = u;
        }
    }
}

public class minDup{
    public minDup(){

    }

    public class heap{
        int size;
        int[][] data;
        public heap(int n){
            size = 0;
            data = new int[n+2][2];
        }
        public void pudh(int v, int n){
            size++;
            data[size][0] = v;
            data[size][1] = n;
            for(int i = size; i > 1; i/=2){
                if(data[i][1] > data[i/2][1]){
                    data[0][0] = data[i][0];
                    data[0][1] = data[i][1];
                    data[i][0] = data[i/2][0];
                    data[i][1] = data[i/2][1];
                    data[i/2][0] = data[0][0];
                    data[i/2][1] = data[0][1];
                }
            }
        }
        public int[] pop(){
            int[] rtn = new int[2];
            rtn[0] = data[1][0];
            rtn[1] = data[1][1];
            data[1][0] = data[size][0];
            data[1][1] = data[size][1];
            size--;
            int look;
            for(int i = 1; 2*i <= size; i = look){
                if(2*i == size) look = 2*i;
                else {
                    look = (data[2*i+1][1] > data[2*i][1])?(2*i+1):(2*i);
                }
                data[0][0] = data[i][0];
                data[0][1] = data[i][1];
                data[i][0] = data[look][0];
                data[i][1] = data[look][1];
                data[look][0] = data[0][0];
                data[look][1] = data[0][1];
            }
            return rtn;
        }
    }

    //downOrder는 인접리스트와 차수, 정점의 수를 받아 차수의 내림차순 정렬을 리턴한다.
    //리턴 배열 a[i]는 입력의 i번째 정점이 몇 번 정점이어야 하는지이다.
    //i는 0에서 시작한다.
    public int[] downOrder(int[][] adjList, int degree[], int N){
        int[] rtn = new int[N];
        heap hp = new heap(N);
        for(int i = 0; i < N; i++){
            hp.pudh(i, degree[i]);
        }
        for(int i = 0; i < N; i++){
            rtn[i] = hp.pop()[0];
        }
        return rtn;
    }
    // 아니면 인접리스트와 차수만 넣거나 인접리스트만 넣어도 된다.
    public int[] downOrder(int degree[], int N){
        int[] rtn = new int[N];
        heap hp = new heap(N);
        for(int i = 0; i < N; i++){
            hp.pudh(i, degree[i]);
        }
        for(int i = 0; i < N; i++){
            rtn[i] = hp.pop()[0];
        }
        return rtn;
    }
    public int[] downOrder(int degree[]){
        int N = degree.length;
        int[] rtn = new int[N];
        heap hp = new heap(N);
        for(int i = 0; i < N; i++){
            hp.pudh(i, degree[i]);
        }
        for(int i = 0; i < N; i++){
            rtn[i] = hp.pop()[0];
        }
        return rtn;
    }



}