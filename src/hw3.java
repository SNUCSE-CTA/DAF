import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;

public class hw3{

	public static void main(String[] args){

			if(args.length!=3){
				System.out.println("Usage: ./program dataFile queryFile numQuery");
				return;
			}

			String dataFileName = args[0];
			String queryFileName = args[1];
			int numQuery = Integer.parseInt(args[2]);

			ProcessIO pio = new ProcessIO(dataFileName, queryFileName, numQuery);

			AdjacentList[] querySet = pio.readQuery();


		
	}

}

/*
 * class for processing I/O
 *
 * Input :
 *   Data graph : has 3112 vertex
 *   Query set : 100 graph, non-sparse(degree > 3)
 *
 * Output :
 *	 .dag file : string of "v1 v2 v3 ..." for each query graph
 *   *total number of line of .dag must be 'numQuery'
 *   **just print output on console(spec)
 *
 * execute code in spec :
 *  java hw3 yeast yeast_400n 100 > yeast_400.dag
 *
 */
class ProcessIO{

    String dataFileName;
    String queryFileName;
    int numQuery;

    ProcessIO(String dataFileName, String queryFileName, int numQuery){
        this.dataFileName = dataFileName;
        this.queryFileName = queryFileName;
        this.numQuery = numQuery;
    }

    //proto
    void readDate(){}


    /*
     *
     * data format for query graph :
     * - 100 graph instances
     * - in each graph :
     *      at first line : t [graph id] [number of vertices] [sum of degree]
     *      2~v+1 line : [vertex id] [label of vertex] [degree of vertex] [v_1] [v_2] ... [v_d]
     *      * v_i are the neighbors of the node
     *
     */
    AdjacentList[] readQuery(){
        AdjacentList[] querySet = new AdjacentList[numQuery];

        String line;
        String[] graphTag;
        String[] vertexTag;

        int numOfVertices;

        try{
            BufferedReader br = new BufferedReader(new FileReader(new File(queryFileName)));

            //for each query
            for(int i=0;i<numQuery;i++){
                line = br.readLine();

                graphTag = line.split(" ");
                numOfVertices = Integer.parseInt(graphTag[2]);

                //for each vertices in one query graph
                for(int j =0;j<numOfVertices;j++){
                    line = br.readLine();
                    vertexTag = line.split(" ");

                }

            }




            return querySet;
        }catch (Exception e) {
            System.out.println(e);
            return null;
        }


    }

    //print all DAG string to console
    void printAllDAGs(){}

}


/*
 * class for building Directed Acyclic Graph(DAG)
 *
 *
 *
 */
class DAG{

    int root;       //root vertex of this DAG
    int numOfNode;  //number of total vertex

    //proto
    DAG(){}
    void buildDAG(){}
    //int selectRoot(){}

}
class AdjacentList{

    int numOfVertex;
    Vertex vertices[];

    AdjacentList(int numOfVertex){
        this.numOfVertex = numOfVertex;
        vertices = new Vertex[numOfVertex];
    }

}
class Vertex{

    Label label;
    int degreeData;
    int degreeQuery;
    boolean visited = false;

}

class Label{

    int labelNum;
    int frequency = 0;

}
