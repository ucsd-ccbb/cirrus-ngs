package org.ucsd.ccbb;

/**
 * S3UploadWorker uploads files to S3
 * 
 * @author Guorong Xu
 */
public class S3UploadWorker implements Runnable {
	/**a string stores the S3 upload URL*/
	public String m_S3UploadURL = "";
	/**a string stores the file name*/
	public String m_fileName = "";
	/**a string stores the log path*/
	public String m_logFile = "";
	/**a ConfigLoader stores the configuration loader*/
	public ConfigLoader m_loader = null;
	
	/**
	 * Construct a S3UploadWorker instance with specified configuration loader, S3 
	 * upload URL, a file name and log path.
	 * @param loader		load the configuration file
	 * @param S3UploadURL	the S3 upload URL
	 * @param fileName		the file name will be set to the field m_fileName
	 * @param logFile		the log path will be set to the field m_logFile
	 */
	public S3UploadWorker(ConfigLoader loader, String S3UploadURL, String fileName, String logFile) {
		this.m_loader = loader;
		this.m_S3UploadURL = S3UploadURL;
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
        //System.out.println(Thread.currentThread().getName()+ " processing file: " + m_bamFile + ".sort.split.bam");       
        try 
        {
			processCommand();
		} catch (InterruptedException e) {
			e.printStackTrace();
		}
        //System.out.println(Thread.currentThread().getName()+ "End.");
    }

	/**
	 * Applies the S3UploadExecutor to process commands.
	 * @throws InterruptedException
	 */
	private void processCommand() throws InterruptedException 
	{					
	    S3UploadExecutor s3u= new S3UploadExecutor();
	    s3u.setShellFile(m_loader.get("s3upload.shell", "/shared/workspace/WGSPipeline/scripts/s3upload.sh"));
	    s3u.setFileName(m_fileName);
	    s3u.setS3UploadURL(m_S3UploadURL);
	    s3u.setLogFile(m_logFile);
	    s3u.execute("", m_fileName);
	}
}