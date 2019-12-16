public class hw3{

	public static void main(String[] args){

			if(args.length!=3){
				System.out.println("Usage: ./program dataFile queryFile numQuery");
				return;
			}

			String dataFileName = args[0];
			String queryFileName = args[1];
			int numQuery = Integer.parseInt(args[2]);



		
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
 *   total number of line of .dag must be 100
 *
 */
class processIO{


}


/*
 * class for building Directed Acyclic Graph(DAG)
 *
 *
 *
 */
class DAG{


}
