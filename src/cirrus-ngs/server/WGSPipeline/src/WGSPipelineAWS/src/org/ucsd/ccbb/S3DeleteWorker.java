package org.ucsd.ccbb;

import java.sql.Timestamp;

/**
 * S3DeleteWorker deletes file from S3
 * 
 * @author Guorong Xu
 */
public class S3DeleteWorker implements Runnable {
	/**a string stores the S3 URL*/
	public String m_S3URL = "";
	/**a string stores the file name*/
	public String m_fileName = "";
	/**a string stores the log path*/
	public String m_logFile = "";
	/**a ConfigLoader stores the configuration loader*/
	public ConfigLoader m_loader = null;
	
	/**
	 * Construct a S3DeleteWorker instance with specified configuration loader, S3 URL, 
	 * a file name and log path.
	 * @param loader		load the configuration file
	 * @param S3URL			the S3 URL
	 * @param fileName		the file name will be set to the field 
	 * 						m_fileName
	 * @param logFile		the log path with be set to the field m_logFile
	 */
	public S3DeleteWorker(ConfigLoader loader, String S3URL, String fileName, String logFile) {
		this.m_loader = loader;
		this.m_S3URL = S3URL;
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
	 * Applies the S3DeleteExecutor to process commands.
	 * @throws InterruptedException
	 */
	private void processCommand() throws InterruptedException 
	{					
	    S3DeleteExecutor s3deleter = new S3DeleteExecutor();
	    s3deleter.setShellFile(m_loader.get("s3delete.shell", "/shared/workspace/WGSPipeline/scripts/s3delete.sh"));
	    s3deleter.setS3URL(m_S3URL);
	    s3deleter.setLogFile(m_logFile);
	    s3deleter.setFileName(m_fileName);
	    s3deleter.execute("", m_fileName);
	}
}