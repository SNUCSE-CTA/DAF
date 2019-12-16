import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.HashSet;
import java.util.Set;

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


    /*
     *
     * Data graph format :
     *  t [graph id] [number of vertices]
     *  v [vertex id] [label of node]
     *  e [vertex 1] [vertex 2] [label of edge]
     *
     * this method read only one date graph in this case.
     *
     */
    AdjacentList readDate(){
        AdjacentList dataGraph = new AdjacentList(1);//initialize for compile
        String line;
        String[] lineTag;
        Set<Label> labelSet = new HashSet<>();

        int[] originalLabels;
        int[] renamedLabels;
        int[] dataLabels;
        int largestLabel = -1;

        int vid, src, dst;
        int labelRenameIndex = 0;
        int numOfLabel = 0;

        try{
            BufferedReader br1 = new BufferedReader(new FileReader(new File(queryFileName)));
            line = br1.readLine();

            //for each query
            while(line!=null){

                lineTag = line.split(" ");

                if(lineTag[0].equals("t")){

                    dataGraph = new AdjacentList(Integer.parseInt(lineTag[2]));

                } else if(lineTag[0].equals("v")){
                    vid = Integer.parseInt(lineTag[1]);

                    //add vertex in data graph
                    Vertex v = new Vertex(vid);
                    dataGraph.vertices[vid] = v;

                    //add label in set
                    int labelValue = Integer.parseInt(lineTag[2]);
                    Label tmpLabel = new Label(labelValue);
                    labelSet.add(tmpLabel);
                    if(largestLabel<labelValue){
                        largestLabel = labelValue;
                    }

                } else if(lineTag[0].equals("e")){
                    //ignore label of edge in this case

                    src = Integer.parseInt(lineTag[1]);
                    dst = Integer.parseInt(lineTag[2]);

                    dataGraph.vertices[src].addNeighbor(dataGraph.vertices[dst]);
                    dataGraph.vertices[dst].addNeighbor(dataGraph.vertices[src]);

                }else{
                    //data graph input error
                }
                line = br1.readLine();
            }
            numOfLabel = labelSet.size();
            Label[] labels = new Label[numOfLabel];

            originalLabels = new int[largestLabel+1]; //data label
            renamedLabels = new int[numOfLabel];
            dataLabels = new int[dataGraph.numOfVertex];

            BufferedReader br2 = new BufferedReader(new FileReader(new File(queryFileName)));
            line = br2.readLine();

            while(line!=null){
                lineTag = line.split(" ");

                if(lineTag[0].equals("t")){

                } else if(lineTag[0].equals("v")){

                    vid = Integer.parseInt(lineTag[1]);
                    int labelValue = Integer.parseInt(lineTag[2]);

                    //dataLabels[vid] =



                } else if(lineTag[0].equals("e")){

                }else{
                    //data graph input error
                }

                line = br2.readLine();
            }

            return dataGraph;
        }catch (Exception e) {
            System.out.println(e);
            return null;
        }
    }


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

    private int id;
    Label label=null;
    int degreeData=0;
    int degreeQuery=0;
    boolean visited = false;
    Set<Vertex> neighborVertices;

    Vertex(int id){
        this.id = id;
        neighborVertices = new HashSet<>();
    }
    public void addNeighbor(Vertex v){
        neighborVertices.add(v);
        degreeData++;
    }

    public boolean equals(Vertex v){
        return id==v.id;
    }

}

class Label{

    int renamedLabel;
    int originalLabel;
    int frequency = 0;

    Label(int originalLabel){
        this.originalLabel = originalLabel;
    }

    @Override
    public int hashCode() {
        return originalLabel%524287;
    }

    public boolean equals(Label label){
        return originalLabel==label.originalLabel;
    }

}
