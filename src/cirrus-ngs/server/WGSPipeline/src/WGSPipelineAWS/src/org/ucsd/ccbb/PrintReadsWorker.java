package org.ucsd.ccbb;

/**
 * PrintReadsWorker prints final BAM for variants calling.
 * 
 * @author Guorong Xu
 */
public class PrintReadsWorker implements Runnable {
	/**a string stores the name of a BAM file*/
	public String m_bamFile = "";
	/**a string stores the log path*/
	public String m_logFile = "";
	/**a ConfigLoader stores the configuration file*/
	public ConfigLoader m_loader = null;
	
	/**
	 * Construct a new PrintReadsWorker instance with specified configuration 
	 * loader, BAM file and log path.
	 * @param loader	load the configuration file
	 * @param bamFile	a BAM file contains information of readsd
	 * @param logFile	the log path will be set to the field m_logFile
	 */
	public PrintReadsWorker(ConfigLoader loader, String bamFile, String logFile) {
		this.m_loader = loader;
		this.m_bamFile = bamFile;
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
        System.out.println(Thread.currentThread().getName()+ " processing file: " + m_bamFile + ".realign.bam");
        
        try 
        {
			processCommand();
		} catch (InterruptedException e) {
			e.printStackTrace();
		}
        System.out.println(Thread.currentThread().getName()+ "End.");
    }

	/**
	 * Applies the PrintReadsExecutor to process commands.
	 * @throws InterruptedException
	 */
	private void processCommand() throws InterruptedException 
	{					
		PrintReadsExecutor pre = new PrintReadsExecutor();
		pre.setShellFile(m_loader.get("post-alignment.print.reads.shell", "/shared/workspace/WGSPipeline/scripts/print_reads.sh"));
		pre.setFileName(m_bamFile);
		pre.setLogFile(m_logFile);
		pre.execute("", m_bamFile);
	}
}