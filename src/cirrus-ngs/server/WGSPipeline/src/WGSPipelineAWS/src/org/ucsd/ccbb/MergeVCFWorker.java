package org.ucsd.ccbb;

import java.util.ArrayList;

/**
 * MergeVCFWorker merges VCF files
 * 
 * @author Guorong Xu
 */
public class MergeVCFWorker implements Runnable {
	/**an arraylist stores a list of VCF files*/
	public ArrayList<String> m_vcfFileList = null;
	/**a string stores the output file name*/
	public String m_outputFileName = "";
	/**a string indicate whether or not merge by region?*/
	public String m_byRegion = "";
	/**a string stores the log path*/
	public String m_logFile = "";
	/**a ConfigLoader stores the configuration loader*/
	public ConfigLoader m_loader = null;
	
	/**
	 * Constructs a MergeVCFWorker instance with specified configuration 
	 * loader, a list of VCF files?. output file name, a string 
	 * indicating whether or not merge by region and log path.
	 * @param loader			load the configuration file
	 * @param vcfFileList		a list of VCF files
	 * @param outputFileName	the string will be stores in m_outputFileName
	 * @param byRegion			"true" indicates merge by region?
	 * @param logFile			the log path
	 */
	public MergeVCFWorker(ConfigLoader loader, ArrayList<String> vcfFileList, String outputFileName, String byRegion, String logFile) {
		this.m_loader = loader;
		this.m_vcfFileList = vcfFileList;
		this.m_outputFileName = outputFileName;
		this.m_byRegion = byRegion;
		this.m_logFile = logFile;
	}

	@Override
	/**
	 * Overrides the function run() from interface Runnable.
	 * The function calls the function processCommand() to run. If an
	 * InterruptedExecption is catched, its throwable and backtrace will be 
	 * print to the standard error stream
	 */
    public void run() {     
        try 
        {
			processCommand();
		} catch (InterruptedException e) {
			e.printStackTrace();
		}
    }

	/**
	 * Applies the MergeVCFExecutor to process commands.
	 * @throws InterruptedException
	 */
	private void processCommand() throws InterruptedException 
	{					
	    MergeVCFExecutor mve= new MergeVCFExecutor();
	    mve.setShellFile(m_loader.get("merge.vcf.shell", "/shared/workspace/WGSPipeline/scripts/merge_vcf.sh"));
	    mve.setInputFileNames(m_vcfFileList);
	    mve.setByRegion(m_byRegion);
	    mve.setOutputFileName(m_outputFileName);
	    mve.setLogFile(m_logFile);
	    mve.execute("", m_outputFileName);
	}
}