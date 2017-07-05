package org.ucsd.ccbb;

/**
 * S3DownloadWorker downloads files from S3
 * 
 * @author Guorong Xu
 */
public class S3DownloadWorker implements Runnable {
	/**a string stores the S3 download URL*/
	public String m_S3DownloadURL = "";
	/**a string stores the file name*/
	public String m_fileName = "";
	/**a string stores the log path*/
	public String m_logFile = "";
	/**a ConfigLoader stores the configuration loader*/
	public ConfigLoader m_loader = null;
	
	/**
	 * Construct a S3DownloadWorker instance with specified configuration loader, S3 
	 * download URL, a file name and log path.
	 * @param loader			load the configuration file
	 * @param S3DownloadURL		the S3 upload URL
	 * @param fileName			the file name will be set to the field 
	 * 							m_fileName
	 * @param logFile			the log path with be set to the field m_logFile
	 */
	public S3DownloadWorker(ConfigLoader loader, String S3DownloadURL, String fileName, String logFile) {
		this.m_loader = loader;
		this.m_S3DownloadURL = S3DownloadURL;
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
	 * Applies the S3DownloadExecutor to process commands.
	 * @throws InterruptedException
	 */
	private void processCommand() throws InterruptedException 
	{					
	    S3DownloadExecutor s3d= new S3DownloadExecutor();
	    s3d.setShellFile(m_loader.get("s3download.shell", "/shared/workspace/WGSPipeline/scripts/s3download.sh"));
	    s3d.setFileName(m_fileName);
	    s3d.setS3DownloadURL(m_S3DownloadURL);
	    s3d.setLogFile(m_logFile);
	    s3d.execute("", m_fileName);
	}
}