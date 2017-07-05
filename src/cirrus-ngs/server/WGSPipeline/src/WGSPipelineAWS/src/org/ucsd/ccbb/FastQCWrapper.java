 package org.ucsd.ccbb;

/**
 * FastQCWrapper runs fastQC
 * 
 * @author Guorong Xu
 */
public class FastQCWrapper implements Runnable {
	/**a string stores the file name*/
	public String m_fileName = "";
	/**a string stores the log path*/
	public String m_logFile = "";
	/**a ConfigLoader stores the configuration loader*/
	public ConfigLoader m_loader = null;
	
	/** Construct a new FastQCWrapper instance with specified configuration 
	 * loader, file name and log path.
	 */
	public FastQCWrapper(ConfigLoader loader, String fileName, String logFile) {
		this.m_loader = loader;
		this.m_fileName = fileName;
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
	 * Applies the PrintReadsExecutor to process commands.
	 * @throws InterruptedException
	 */
	private void processCommand() throws InterruptedException 
	{					
		FastqcExecutor qcExecutor = new FastqcExecutor();
		qcExecutor.setShellFile(m_loader.get("fastqc.shell", "/shared/workspace/WGSPipeline/scripts/fastqc.sh"));			
		qcExecutor.setFastqFile(m_fileName);
		qcExecutor.setLogFile(m_logFile);
		qcExecutor.execute("", "");
	}
}