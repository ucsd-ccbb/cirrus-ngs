import java.io.BufferedReader;
import java.io.InputStreamReader;

/**
 * The class is to track the PBS queue status. For example, how many jobs are running in the queue.
 * @author guorong
 *
 */
public class PBSTracker {
	   public void trackPBSQueue(int minutes, String shellScript)
	    {
		   int jobs = 0;
	        while (true) {
	            try {
	                Thread.sleep(minutes * 60000);
	                
		        		String[] commandLine = new String[]{"qstat"};
		        		Process process = Runtime.getRuntime().exec(commandLine);

	                String line;
	                BufferedReader in = new BufferedReader(new InputStreamReader(process.getInputStream()));
	                
	                boolean isDone = true;
	                int index = 0;
	                while ((line = in.readLine()) != null) {
	                    if (line.indexOf(shellScript) > -1)
	                    		isDone = false;
	                    index++;
	                }
	                
	                if (index > 2)
	                {
	                		if ((index - 2) != jobs)
	                		{
	                			jobs = index - 2;
	                			System.out.println("");
	                			System.out.println((index - 2) + "  job(s) are running...");
	                		}
	                		else
	                			System.out.print("...");
	                }
	                
	                if (isDone)
	                {
	                		System.out.println("");
	                		System.out.println("No jobs is running...");
	                		break;
	                }
	                
	                in.close();
	            } catch (Exception e) {
	            }
	        }
	    }
	   
		public int checkNodeNum()
		{
	        int num = 0;
	        try {
	    	   		String[] commandLine = new String[]{"qhost"};
				Process process = Runtime.getRuntime().exec(commandLine);
			
				String line;
				BufferedReader in = new BufferedReader(new InputStreamReader(process.getInputStream()));
				    
				while ((line = in.readLine()) != null) {
					System.out.println(line);
				    if (line.startsWith("ip-"))
				    num++;
				}
	            in.close();
	        } catch (Exception e) {
	        		System.out.println(e);
	        }
	       
	       return num;
		}
}
