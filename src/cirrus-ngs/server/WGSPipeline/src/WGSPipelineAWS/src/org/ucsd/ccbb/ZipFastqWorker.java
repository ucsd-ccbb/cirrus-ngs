package org.ucsd.ccbb;


/**
 * ZipFastqWorker compresses fastq file.
 * 
 * @author Guorong Xu
 */
public class ZipFastqWorker implements Runnable {
	/**a string stores the name of the fastqFile*/
	public String m_fastqFile = "";
	/**a string stores the log path*/
	public String m_logFile = "";
	/**a ConfigLoader stores the configuration loader*/
	public ConfigLoader m_loader = null;
	
	/**
	 * Constructs a ZipFastqWorker instance with specified configuration loader,
	 * a fastq format file name and log path.  
	 * @param loader		load the configuration file
	 * @param fastqFile		a fastq format file that needs to be compressed
	 * @param logFile		the log path will be set to the field m_logFile
	 */
	public ZipFastqWorker(ConfigLoader loader, String fastqFile, String logFile) {
		this.m_loader = loader;
		this.m_fastqFile = fastqFile;
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
        System.out.println(Thread.currentThread().getName()+ " processing file: " + m_fastqFile + ".fastq");
        
        try 
        {
			processCommand();
		} catch (InterruptedException e) {
			e.printStackTrace();
		}
        System.out.println(Thread.currentThread().getName()+ "End.");
    }

	/**
	 * Applies the ZipFastqExecutor to process commands.
	 * @throws InterruptedException
	 */
	private void processCommand() throws InterruptedException 
	{					
		ZipFastqExecutor zipFastq = new ZipFastqExecutor();
		zipFastq.setShellFile(m_loader.get("compressFastq.shell", "/shared/workspace/WGSPipeline/scripts/compress_fastq.sh"));
		zipFastq.setLogFile(m_logFile);
		zipFastq.execute("", m_fastqFile);
	}
}