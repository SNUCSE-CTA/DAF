import java.util.*;
import java.math.*;

class minDupMain{
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

        minDup md = new minDup();
        int[] arr = new int[10];
        arr = md.greedyLog(ad, degree, 10);
        for(int i = 0; i < 10; i++)
            System.out.println(arr[i]);
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
    //리턴 배열 a[i]는 적합하게 재배열된 배열이다.
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

    // 그리디스럽게 정점을 0번부터 보면서 어디에 넣으면 중복 정점이 적을지 봐보자.
    // 패스트리 겹치는 정점 줄이기 이슈의 4번 코멘트 참조
    // 입력은 인접리스트, 차수, 크기
    public int[] greedyLog(int[][] adjList, int[] degree, int N){
        int[] dup = new int[N];
        double[] logdup = new double[N];
        int[] arr = new int[N];
        int[] where = new int[N];
        arr[0] = 0;
        where[0] = 0;
        double sum;

        for(int v = 1; v < N; v++){
            dup = new int[N];
            for(int i = 0; i < v; i++){
                for(int j = 0; j < degree[arr[i]]; j++){
                    if(adjList[arr[i]][j] >= v) continue;
                    else if(where[adjList[arr[i]][j]] < where[arr[i]]){
                        dup[arr[i]] += dup[adjList[arr[i]][j]];
                    }
                }
                if(dup[arr[i]] == 0)
                    dup[arr[i]] = 1;
            }
            for(int i = 0; i < v; i++){
                logdup[i] = Math.log(dup[i])/(where[i]+1);
            }

            sum = 0;
            for(int i = 0; i < degree[v]; i++){
                sum += logdup[adjList[v][i]];
            }
            sum /= 2;
            double nujuckhap = 0;
            int a=0, b=0;
            for(int loc = 0; loc <= v; loc++){
                boolean chk = false;
                for(int i = 0; i < degree[v]; i++){
                    if(arr[loc] == adjList[v][i]){
                        chk = true;
                        break;
                    }
                }
                if(chk){
                    nujuckhap+=dup[arr[loc]];
                    if(nujuckhap<sum) a = loc;
                    else{
                        b = loc;
                        break;
                    }
                }
            }
            a = (a+b+1)/2;
            where[v] = a;
            for(int i = 0; i < v; i++) {
                if (where[i] >= a) where[i]++;
            }
            for(int i = v; i>a; i--){
                arr[i] = arr[i-1];
            }
            arr[a] = v;
        }

        return arr;
    }
}