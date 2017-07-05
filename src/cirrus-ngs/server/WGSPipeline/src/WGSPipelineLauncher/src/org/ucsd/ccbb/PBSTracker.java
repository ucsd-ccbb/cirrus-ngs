package org.ucsd.ccbb;
import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.sql.Timestamp;
import java.util.Date;

/**
 * The class is to track the PBS queue status. 
 * For example, how many jobs are running in the queue.
 * 
 * @author Guorong Xu
 */
public class PBSTracker {
	public void trackPBSQueue(int minutes, String shellScript)
	{
		int jobs = 0;
	    while (true) {
	    	try {
	    			Date date= new Date();
	    			Thread.sleep(minutes * 60000);
	                
		        	String[] commandLine = new String[]{"qstat"};
		        	Process process = Runtime.getRuntime().exec(commandLine);

	                String line;
	                BufferedReader in = new BufferedReader(new InputStreamReader(process.getInputStream()));
	                
	                boolean isDone = true;
	                int count = 0;
	                while ((line = in.readLine()) != null) {
	                    if (line.indexOf(shellScript) > -1)
	                    		isDone = false;
	                    count++;
	                }
	                //ignore the first two lines 
	                if (count > 2)
	                {
	                	if ((count - 2) != jobs)
	                	{
	                		jobs = count - 2;
	                		System.out.println("");
	                		//System.out.println(new Timestamp(date.getTime()) + "  job(s) are running...");
	                		System.out.println(new Timestamp(date.getTime()) + " running job(s): " + (count - 2));
	                	}
	                	//else
	                		//need to be modified here: print out the current time
	                		//how to print out how long the program has run?
	                		//System.out.print("...");
	                }
	                
	                if (isDone)
	                {
	                	System.out.println("");
	                	System.out.println(new Timestamp(date.getTime()) + " No jobs is running...");
	                	break;
	                }
	                
	                in.close();
	    	} catch (Exception e) {
	    	}
	    }
	}
	   
	/**
	 * Checks how many command lines started with the word "node" ?
	 * @return	an integer indicating the number of commands started with the 
	 * 			word "node"
	 */
	public int checkNodeNum()
	{
		int num = 0;
	    try {
	    		String[] commandLine = new String[]{"qhost"};
	    	   	Process process = Runtime.getRuntime().exec(commandLine);
			
	    	   	String line;
	    	   	BufferedReader in = new BufferedReader(new InputStreamReader(process.getInputStream()));
				    
	    	   	while ((line = in.readLine()) != null) {
	    			Date date= new Date();
	    	   		System.out.println(new Timestamp(date.getTime()) + line);
	    	   		if (line.startsWith("node"))
	    	   			num++;
				}
	            in.close();
	    } catch (Exception e) {
			Date date2= new Date();
			System.out.println(new Timestamp(date2.getTime()));
	        System.out.println(e);
	    }
	       
	    return num;
	}
}
